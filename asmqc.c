/* asmqc.c -- quick assembly-vs-reference QC from a `ropebwt3 lift` liftover map.
 *
 * Motivation: PHGv2 assembly QC (SNPs, indel bp, % of the reference covered) is
 * computed from an AnchorWave gVCF, which takes hours and hundreds of GB per
 * assembly. The liftover map already pairs every unique reference k-mer sampled
 * at a fixed stride with its exact occurrence(s) in each assembly of the index,
 * so the anchor-level picture (where the assembly is collinear with the
 * reference, where it is not, and roughly how diverged it is) can be read off it
 * in seconds -- before any base-level alignment.
 *
 * Method (lift-only mode), for every assembly sequence x reference sequence pair:
 *  1. Sort that pair's liftover points by reference position and chain them in
 *     two stages, both orientations. Stage 1 (blocks): minimap2-style DP with
 *     bounded look-back along a near-constant diagonal (|delta_asm - delta_ref|
 *     <= --block-diag), at least -n anchors per block. Stage 2: join same-strand
 *     blocks in order with a log-scaled gap penalty across gaps up to --max-gap
 *     (3 Mb by default for syntenic chains, to match AnchorWave's multi-Mb
 *     deletions; --max-gap-off for other chromosomes / scaffolds). Two stages keep
 *     a single stray anchor from being pulled across a large gap.
 *  2. Greedily extract chains (best score first, points used once). Chains to the
 *     assembly sequence with the same chromosome suffix as the reference sequence
 *     are syntenic; the orientation carrying more reference bp is "forward", the
 *     other "inverted". Chains to other chromosomes / scaffolds are reported apart.
 *  3. Coverage: union of syntenic chain spans on the reference = mapped bp (the
 *     analogue of PHG percentRefMapped); the rest is gap.
 *  4. Indel lower bound: between consecutive chained anchors, d = delta_asm -
 *     delta_ref > 0 is net inserted bp, < 0 net deleted bp. An interval holding
 *     both an insertion and a deletion only shows the net, so these are floors.
 *  5. Divergence estimate: on "tight" intervals (delta_ref <= --tight-gap and
 *     |d| <= --tight-diag, i.e. plainly collinear sequence), count the reference
 *     stride positions that are anchorable in the panel (some assembly carries
 *     that exact k-mer) but missing from this assembly. With f = present/(present
 *     + missing) the fraction of k-mers carried intact, the per-bp variant-event
 *     rate is estimated as 1 - f^(1/k) (independent-sites approximation; variants
 *     cluster, so this is a relative ranking more than an absolute rate).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "rb3priv.h"
#include "io.h"
#include "lift.h"
#include "ketopt.h"
#include "kthread.h"
#define HAVE_KALLOC
#include "ksw2.h"
#include <pthread.h>

static FILE *qc_chain_fp = 0; // --chains output (NULL = off)
static FILE *qc_var_fp = 0;   // --variants output (NULL = off)
static pthread_mutex_t qc_chain_lock = PTHREAD_MUTEX_INITIALIZER;
static const char *qc_cls_name[] = { "syn", "syn_inv", "other_chr", "scaffold" };

#define QC_SYN_FWD 0
#define QC_SYN_INV 1
#define QC_OTHER   2
#define QC_SCAF    3

/* ---- base-level comparison between chained anchors (asmqc-fa) -----------
 *
 * Between two consecutive anchors of a syntenic chain, the reference segment R and the assembly
 * segment A (reverse-complemented on an inverted chain) are compared with AnchorWave's gap model:
 * match 0, mismatch 4, 2-piece affine gaps min(4+2l, 80+l) (ksw2, converted to an equivalent
 * positive-match scoring, see qa_ksw):
 *   - equal length and <= QA_GAPLESS_MM mismatches: gapless, mismatches = SNPs;
 *   - small and near-diagonal: one global alignment, unless deleting R + inserting A scores better;
 *   - otherwise breakpoints first: extend from both flanking anchors with z-drop (crossing SNPs and
 *     small indels, stopping where similarity collapses), re-align the extended pieces, and treat the
 *     middle between the breakpoints as a clean indel (one side empty), or split it on shared unique
 *     17-/13-mers (or place a short piece by 11-mer diagonal voting) and recurse, or take the better
 *     of alignment vs replacement;
 *   - replacement = AnchorWave's representation of unalignable sequence: R as one deletion, A as one
 *     insertion.
 * SNPs are also tallied by local density (>= 10 / >= 40 other SNPs within +-50 bp).
 */
#define QA_MAX_LEN   200000
#define QA_GAPLESS_MM 0.03
#define QA_SK    17
#define QA_SK2   13     // second seeding pass for diverged segments with no shared unique 17-mer
#define QA_VK    11     // diagonal voting k-mer (repeats allowed) for lopsided / repetitive gaps
#define QA_VBIN  256    // diagonal bin width
#define QA_VPAD  256    // window padding around the voted diagonal
#define QA_VOCC  8      // ignore voting k-mers occurring more often than this in the shorter sequence
#define QA_HIST  4096   // exact indel-size bins 1..QA_HIST-1, last bin = larger

typedef struct qa_var_s { int64_t pos, rlen, alen; uint8_t r, a, type; } qa_var_t;

typedef struct {
	uint8_t *tb; int64_t m_tb;
	int32_t *h0, *h1, *f0, *f1; int32_t m_w;
	uint8_t *rc; int64_t m_rc;
	int32_t *ev; int64_t n_ev, m_ev;        // pending events of one DP: >0 insertion, <0 deletion
	uint8_t *ops; int64_t n_ops, m_ops;     // traceback path, end to start: 0 match, 1 mismatch, 2 ins, 3 del
	uint64_t *ka, *kb; int64_t m_ka, m_kb;  // seeding scratch
	int dump;                               // --variants: collect SNP / replacement records
	struct qa_var_s *v, *pv;                // committed / pending (one DP)
	int64_t n_v, m_v, n_pv, m_pv;
} qa_buf_t;
#define QA_VSNP 0
#define QA_VREP 1
#define QA_VINS 2
#define QA_VDEL 3

static inline void qa_var_push(qa_var_t **v, int64_t *n, int64_t *m, int64_t pos, int64_t rlen, int64_t alen, uint8_t r, uint8_t a, uint8_t type)
{
	if (*n == *m) *m = *m? *m * 2 : 1024, *v = (qa_var_t*)realloc(*v, *m * sizeof(qa_var_t));
	(*v)[*n].pos = pos, (*v)[*n].rlen = rlen, (*v)[*n].alen = alen, (*v)[*n].r = r, (*v)[*n].a = a, (*v)[*n].type = type;
	++*n;
}

typedef struct {
	int64_t aln, snp, ins_n, ins_bp, ins_max, del_n, del_bp, del_max, rep_n;
	int64_t snp_d10, snp_d40;              // SNPs with >= 10 / >= 40 other SNPs within +-50 bp (local mismatch >= ~10% / ~40%)
	int64_t *hist;                          // [2*QA_HIST]: insertion sizes, then deletion sizes
} qa_stat_t;

static void qa_buf_free(qa_buf_t *b)
{
	free(b->tb); free(b->h0); free(b->h1); free(b->f0); free(b->f1); free(b->rc); free(b->ev); free(b->ka); free(b->kb); free(b->v); free(b->pv); free(b->ops);
	memset(b, 0, sizeof(*b));
}

static inline void qa_event(qa_stat_t *st, int64_t len, int is_del, qa_buf_t *b, int64_t pos) // pos: reference position (0-based)
{
	if (len > 0 && b && b->dump) qa_var_push(&b->v, &b->n_v, &b->m_v, pos, is_del? len : 0, is_del? 0 : len, 0, 0, is_del? QA_VDEL : QA_VINS);
	int64_t bin = len < QA_HIST? len : QA_HIST - 1;
	if (len <= 0) return;
	if (st->hist == 0) st->hist = (int64_t*)calloc(2 * QA_HIST, sizeof(int64_t));
	st->hist[is_del * QA_HIST + bin]++;
	if (is_del) { st->del_n++, st->del_bp += len; if (len > st->del_max) st->del_max = len; }
	else { st->ins_n++, st->ins_bp += len; if (len > st->ins_max) st->ins_max = len; }
}

static void qa_replace(qa_stat_t *st, int64_t lr, int64_t la, qa_buf_t *b, int64_t r0, int why)
{
	if (b->dump) qa_var_push(&b->v, &b->n_v, &b->m_v, r0, lr, la, (uint8_t)why, 0, QA_VREP);
	if (lr > 0) qa_event(st, lr, 1, 0, r0); // the REP record already describes it
	if (la > 0) qa_event(st, la, 0, 0, r0);
	if (lr > 0 && la > 0) st->rep_n++;
}

static inline void qa_push_ev(qa_buf_t *b, int32_t v)
{
	if (b->n_ev == b->m_ev) b->m_ev = b->m_ev? b->m_ev * 2 : 256, b->ev = (int32_t*)realloc(b->ev, b->m_ev * sizeof(int32_t));
	b->ev[b->n_ev++] = v;
}

static inline int qa_acgt(uint8_t c) { return c >= 1 && c <= 4; }

static inline void qa_push_op(qa_buf_t *b, uint8_t op, int64_t n)
{
	if (b->n_ops + n > b->m_ops) { b->m_ops = (b->n_ops + n) * 2; b->ops = (uint8_t*)realloc(b->ops, b->m_ops); }
	memset(b->ops + b->n_ops, op, n); b->n_ops += n;
}


/* AnchorWave's gap model: match 0, mismatch 4, 2-piece affine gaps min(4+2l, 80+l) (WFA penalties in
 * AnchorWave proali), here with ksw2's banded global alignment. ksw2 needs a positive match score
 * (with match 0 its SIMD scores go wrong), so the penalties are converted exactly: with match +a,
 * mismatch -(x-a), gap extend -(e-a/2) and the same opens, score = a(lr+la)/2 - penalty for every
 * alignment, so both rank alignments identically. a = 2: +2 / -2 / open 4 ext 1 / open 80 ext 0. */
#define QA_KSW_A     2
#define QA_KSW_MIS   4
#define QA_KSW_GO1   4
#define QA_KSW_GE1   2
#define QA_KSW_GO2   80
#define QA_KSW_GE2   1
#define QA_KSW_W     200            // band beyond the length difference
#define QA_KSW_CELLS 400000000LL    // (lr + la) x band cap per alignment (~bytes of traceback)

// local SNP-density classes of the SNPs of one alignment (positions ascending)
static void qa_snp_classes(const int64_t *p, int64_t n, qa_stat_t *st)
{
	int64_t i, lo = 0, hi = 0;
	for (i = 0; i < n; ++i) {
		while (p[lo] < p[i] - 50) ++lo;
		while (hi < n && p[hi] <= p[i] + 50) ++hi;
		if (hi - lo - 1 >= 10) st->snp_d10++;
		if (hi - lo - 1 >= 40) st->snp_d40++;
	}
}

/* Global 2-piece affine alignment of R[0,lr) (target) vs A[0,la) (query). Returns 1 and adds SNPs,
 * indels and aligned bp to st; -1 when the band would be too large; 0 (nothing added) when check_rep
 * is set and deleting R + inserting A scores better under the same model -- the band cannot represent
 * that path itself, AnchorWave's unbanded WFA can. */
static inline int64_t qa_gapcost(int64_t l)
{
	int64_t c1 = QA_KSW_GO1 + QA_KSW_GE1 * l, c2 = QA_KSW_GO2 + QA_KSW_GE2 * l;
	return l <= 0? 0 : c1 < c2? c1 : c2;
}

static int qa_ksw(const uint8_t *R, int64_t lr, const uint8_t *A, int64_t la, qa_buf_t *b, qa_stat_t *st, int64_t r0, int check_rep)
{
#define QA_KSW_B (QA_KSW_MIS - QA_KSW_A)
	static const int8_t mat[25] = { QA_KSW_A,-QA_KSW_B,-QA_KSW_B,-QA_KSW_B,-1, -QA_KSW_B,QA_KSW_A,-QA_KSW_B,-QA_KSW_B,-1,
		-QA_KSW_B,-QA_KSW_B,QA_KSW_A,-QA_KSW_B,-1, -QA_KSW_B,-QA_KSW_B,-QA_KSW_B,QA_KSW_A,-1, -1,-1,-1,-1,-1 };
	int64_t w = (lr > la? lr - la : la - lr) + QA_KSW_W, i, j, k, n_snp = 0;
	uint8_t *t, *qs;
	int64_t *sp;
	ksw_extz_t ez;
	if (lr > INT32_MAX / 4 || la > INT32_MAX / 4 || (lr + la) * (2 * w + 1) > QA_KSW_CELLS) return -1;
	t = (uint8_t*)malloc(lr + la + 2), qs = t + lr + 1;
	for (i = 0; i < lr; ++i) t[i] = qa_acgt(R[i])? R[i] - 1 : 4;
	for (i = 0; i < la; ++i) qs[i] = qa_acgt(A[i])? A[i] - 1 : 4;
	memset(&ez, 0, sizeof(ez));
	ksw_extd2_sse(0, (int)la, qs, (int)lr, t, 5, mat, QA_KSW_GO1, QA_KSW_GE1 - QA_KSW_A / 2, QA_KSW_GO2, QA_KSW_GE2 - QA_KSW_A / 2,
				  (int)w, -1, 0, 0, &ez);
	// back in AnchorWave units: penalty = a(lr+la)/2 - score
	if (check_rep && (ez.score <= KSW_NEG_INF / 2 || (QA_KSW_A / 2) * (lr + la) - ez.score > qa_gapcost(lr) + qa_gapcost(la))) {
		// AnchorWave's own scoring prefers deleting R and inserting A over this alignment
		kfree(0, ez.cigar); free(t);
		return 0;
	}
	sp = (int64_t*)malloc((lr < la? lr : la) * sizeof(int64_t) + 8);
	for (k = 0, i = j = 0; k < ez.n_cigar; ++k) {
		int op = ez.cigar[k] & 0xf, len = ez.cigar[k] >> 4, z;
		if (op == 0) { // M: target i, query j
			for (z = 0; z < len; ++z, ++i, ++j) {
				if (t[i] != qs[j] && t[i] < 4 && qs[j] < 4) {
					sp[n_snp++] = r0 + i;
					if (b->dump) qa_var_push(&b->v, &b->n_v, &b->m_v, r0 + i, 1, 1, R[i], A[j], QA_VSNP);
				}
			}
			st->aln += len;
		} else if (op == 1) { qa_event(st, len, 0, b, r0 + i); j += len; } // I: assembly has extra bases
		else if (op == 2) { qa_event(st, len, 1, b, r0 + i); i += len; }   // D: reference bases missing
	}
	st->snp += n_snp;
	qa_snp_classes(sp, n_snp, st);
	kfree(0, ez.cigar);
	free(sp); free(t);
	return 1;
}

static int qa_cmp_u64(const void *a, const void *b)
{
	uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
	return x < y? -1 : x > y? 1 : 0;
}

// k-mers occurring exactly once in s[0,l): packed (code << 30 | pos), sorted by code; returns count
static int64_t qa_uniq_kmers(const uint8_t *s, int64_t l, int32_t sk, uint64_t **buf, int64_t *m)
{
	int64_t i, n = 0, w = 0, len = 0;
	uint64_t code = 0, mask = (1ULL << (2 * sk)) - 1;
	if (l < sk) return 0;
	if (l > *m) *m = l, *buf = (uint64_t*)realloc(*buf, l * sizeof(uint64_t));
	for (i = 0; i < l; ++i) {
		if (!qa_acgt(s[i])) { len = 0, code = 0; continue; }
		code = (code << 2 | (uint64_t)(s[i] - 1)) & mask;
		if (++len >= sk) (*buf)[n++] = code << 30 | (uint64_t)(i - sk + 1);
	}
	qsort(*buf, n, sizeof(uint64_t), qa_cmp_u64);
	for (i = 0; i < n; ) { // keep codes that occur once
		int64_t j = i + 1;
		while (j < n && (*buf)[j] >> 30 == (*buf)[i] >> 30) ++j;
		if (j - i == 1) (*buf)[w++] = (*buf)[i];
		i = j;
	}
	return w;
}

static void qa_align(const uint8_t *R, int64_t lr, const uint8_t *A, int64_t la, int depth, qa_buf_t *b, qa_stat_t *st, int64_t r0);

/* Seeds for a gap the way an FM-index / minimap2 sees it: every shared QA_CK-mer, repeats included
 * (capped at QA_COCC copies per side), chained by dynamic programming (gain = newly matched bases,
 * gap cost 0.01*k*|d| + 0.5*log2|d|) -- unique-k-mer seeding finds nothing in TE-rich gaps where
 * every k-mer repeats. Returns the best chain's anchors, ascending in both coordinates. */
static int QA_CK = 19, QA_COCC = 32; // tunable for experiments (env QA_CK / QA_COCC / QA_CMIN)
static double QA_CMIN = 0.0;          // minimum chain score as a fraction of the shorter side
#define QA_CMAXP 4000000
#define QA_CH    64
#define QA_CMAXD 5000
static int64_t qa_kmers_all(const uint8_t *s, int64_t l, int32_t k, uint64_t **buf, int64_t *m)
{
	int64_t i, n = 0, len = 0;
	uint64_t code = 0, mask = (1ULL << (2 * k)) - 1;
	if (l < k) return 0;
	if (l > *m) *m = l, *buf = (uint64_t*)realloc(*buf, l * sizeof(uint64_t));
	for (i = 0; i < l; ++i) {
		if (!qa_acgt(s[i])) { len = 0, code = 0; continue; }
		code = (code << 2 | (uint64_t)(s[i] - 1)) & mask;
		if (++len >= k) (*buf)[n++] = code << 30 | (uint64_t)(i - k + 1);
	}
	qsort(*buf, n, sizeof(uint64_t), qa_cmp_u64);
	return n;
}

static int64_t qa_chain_kmers(const uint8_t *R, int64_t lr, const uint8_t *A, int64_t la, qa_buf_t *b, int64_t **sr_, int64_t **sa_)
{
	int64_t nr, na, i, j, np = 0, best = -1, n, *pv;
	uint64_t *pr;
	double *f, bf = 0;
	*sr_ = *sa_ = 0;
	if (lr >= (1LL << 30) || la >= (1LL << 30)) return 0;
	{ static int init = 0; if (!init) { char *e; init = 1;
		if ((e = getenv("QA_CK"))) QA_CK = atoi(e);
		if ((e = getenv("QA_COCC"))) QA_COCC = atoi(e);
		if ((e = getenv("QA_CMIN"))) QA_CMIN = atof(e); } }
	nr = qa_kmers_all(R, lr, QA_CK, &b->ka, &b->m_ka);
	na = qa_kmers_all(A, la, QA_CK, &b->kb, &b->m_kb);
	pr = (uint64_t*)malloc(QA_CMAXP * sizeof(uint64_t));
	for (i = j = 0; i < nr && j < na && np < QA_CMAXP; ) {
		uint64_t ci = b->ka[i] >> 30, cj = b->kb[j] >> 30;
		int64_t i2 = i, j2 = j, x, y;
		if (ci < cj) { while (i < nr && b->ka[i] >> 30 == ci) ++i; continue; }
		if (ci > cj) { while (j < na && b->kb[j] >> 30 == cj) ++j; continue; }
		while (i2 < nr && b->ka[i2] >> 30 == ci) ++i2;
		while (j2 < na && b->kb[j2] >> 30 == cj) ++j2;
		if (i2 - i <= QA_COCC && j2 - j <= QA_COCC)
			for (x = i; x < i2; ++x) for (y = j; y < j2 && np < QA_CMAXP; ++y)
				pr[np++] = (b->ka[x] & 0x3fffffff) << 32 | (b->kb[y] & 0x3fffffff);
		i = i2, j = j2;
	}
	if (np == 0) { free(pr); return 0; }
	qsort(pr, np, sizeof(uint64_t), qa_cmp_u64);
	f = (double*)malloc(np * sizeof(double)), pv = (int64_t*)malloc(np * sizeof(int64_t));
	for (i = 0; i < np; ++i) {
		int64_t ri = (int64_t)(pr[i] >> 32), ai = (int64_t)(pr[i] & 0xffffffff);
		double fi = QA_CK; int64_t pi = -1;
		for (j = i - 1; j >= 0 && i - j <= QA_CH; --j) {
			int64_t dr = ri - (int64_t)(pr[j] >> 32), da = ai - (int64_t)(pr[j] & 0xffffffff), d, g;
			double sc;
			if (dr <= 0 || da <= 0) continue;
			d = dr > da? dr - da : da - dr;
			if (d > QA_CMAXD) continue;
			g = dr < da? dr : da; if (g > QA_CK) g = QA_CK;
			sc = f[j] + g - (d? 0.01 * QA_CK * d + 0.5 * log2((double)d) : 0.0);
			if (sc > fi) fi = sc, pi = j;
		}
		f[i] = fi, pv[i] = pi;
		if (fi > bf) bf = fi, best = i;
	}
	if (best < 0 || bf < QA_CMIN * (lr < la? lr : la)) { free(pr); free(f); free(pv); return 0; }
	for (n = 0, i = best; i >= 0; i = pv[i]) ++n;
	*sr_ = (int64_t*)malloc(n * sizeof(int64_t)), *sa_ = (int64_t*)malloc(n * sizeof(int64_t));
	for (j = n, i = best; i >= 0; i = pv[i]) --j, (*sr_)[j] = (int64_t)(pr[i] >> 32), (*sa_)[j] = (int64_t)(pr[i] & 0xffffffff);
	free(pr); free(f); free(pv);
	return n;
}

/* Lopsided or repetitive gap with no shared unique seeds (e.g. a full TE in the reference, a solo LTR in
 * the assembly): find where the shorter sequence sits in the longer one by voting on diagonals with
 * (possibly repeated) QA_VK-mers, align it to that window, and count the rest of the longer sequence as
 * one deletion (reference longer) or insertion (assembly longer) on each side. Returns 1 if aligned. */
static int qa_diag(const uint8_t *R, int64_t lr, const uint8_t *A, int64_t la, qa_buf_t *b, qa_stat_t *st, int64_t r0)
{
	int ref_long = lr >= la;
	const uint8_t *S = ref_long? A : R, *L = ref_long? R : A;
	int64_t ls = ref_long? la : lr, ll = ref_long? lr : la, i, n = 0, nb, best = -1, bv = 0, tot = 0, ws, we;
	uint64_t code = 0, mask = (1ULL << (2 * QA_VK)) - 1;
	int32_t len = 0, *votes;
	int r;
	if (ls < 200 || ls > QA_MAX_LEN || ll > 50 * (int64_t)QA_MAX_LEN) return 0;
	if (ls > b->m_ka) b->m_ka = ls, b->ka = (uint64_t*)realloc(b->ka, ls * sizeof(uint64_t));
	for (i = 0; i < ls; ++i) { // k-mers of the shorter sequence: (code << 30 | pos), sorted
		if (!qa_acgt(S[i])) { len = 0, code = 0; continue; }
		code = (code << 2 | (uint64_t)(S[i] - 1)) & mask;
		if (++len >= QA_VK) b->ka[n++] = code << 30 | (uint64_t)(i - QA_VK + 1);
	}
	if (n == 0) return 0;
	qsort(b->ka, n, sizeof(uint64_t), qa_cmp_u64);
	nb = (ll + ls) / QA_VBIN + 2; // diagonal = posL - posS + ls, binned
	votes = (int32_t*)calloc(nb, sizeof(int32_t));
	len = 0, code = 0;
	for (i = 0; i < ll; ++i) {
		int64_t lo, hi, j;
		if (!qa_acgt(L[i])) { len = 0, code = 0; continue; }
		code = (code << 2 | (uint64_t)(L[i] - 1)) & mask;
		if (++len < QA_VK) continue;
		{ // range of this code among S's k-mers
			uint64_t key = code << 30;
			int64_t a0 = 0, a1 = n;
			while (a0 < a1) { int64_t m = (a0 + a1) >> 1; if (b->ka[m] < key) a0 = m + 1; else a1 = m; }
			lo = a0, hi = a0;
			while (hi < n && b->ka[hi] >> 30 == code) ++hi;
		}
		if (hi - lo == 0 || hi - lo > QA_VOCC) continue;
		for (j = lo; j < hi; ++j) {
			int64_t d = (i - QA_VK + 1) - (int64_t)(b->ka[j] & 0x3fffffff) + ls;
			votes[d / QA_VBIN]++, ++tot;
		}
	}
	for (i = 0; i < nb; ++i) { // best pair of adjacent bins (an indel can split a diagonal)
		int64_t v = votes[i] + (i + 1 < nb? votes[i + 1] : 0);
		if (v > bv) bv = v, best = i;
	}
	free(votes);
	// require a clear diagonal: >= 10 votes and well above an even spread
	if (best < 0 || bv < 10 || bv * nb < 20 * tot) return 0;
	{
		int64_t d0 = best * QA_VBIN - ls; // posL - posS at the bin start
		ws = d0 - QA_VPAD; we = d0 + QA_VBIN * 2 + ls + QA_VPAD;
		if (ws < 0) ws = 0;
		if (we > ll) we = ll;
		if (we - ws < ls / 2) return 0;
	}
	r = ref_long? qa_ksw(R + ws, we - ws, A, la, b, st, r0 + ws, 1)
	             : qa_ksw(R, lr, A + ws, we - ws, b, st, r0, 1);
	if (r != 1) return 0;
	if (ref_long) { if (ws > 0) qa_event(st, ws, 1, b, r0); if (ll - we > 0) qa_event(st, ll - we, 1, b, r0 + we); }
	else { if (ws > 0) qa_event(st, ws, 0, b, r0); if (ll - we > 0) qa_event(st, ll - we, 0, b, r0 + lr); }
	return 1;
}

static void qa_seeded(const uint8_t *R, int64_t lr, const uint8_t *A, int64_t la, int depth, qa_buf_t *b, qa_stat_t *st, int64_t r0)
{
	int64_t nr, na, i, j, ns = 0, *sr, *sa, *tail, *prev, nl = 0, *keep, nk, pr, pa, rb = r0;
	int32_t sk = QA_CK, pass, chained;
	sr = sa = 0;
	// experimental (off by default, QA_CHAIN=1): all shared k-mers incl. repeats, best chain; on Oh43 chr10
	// it did not recover more AnchorWave SNPs than unique-k-mer seeding (k=11-19 tried)
	ns = getenv("QA_CHAIN")? qa_chain_kmers(R, lr, A, la, b, &sr, &sa) : 0;
	chained = ns > 0;
	for (pass = 0; pass < 2 && ns == 0; ++pass) { // fallback: unique 17-mers; if none are shared, unique 13-mers
		sk = pass == 0? QA_SK : QA_SK2;
		nr = qa_uniq_kmers(R, lr, sk, &b->ka, &b->m_ka);
		na = qa_uniq_kmers(A, la, sk, &b->kb, &b->m_kb);
		free(sr); free(sa);
		sr = (int64_t*)malloc((nr < na? nr : na) * sizeof(int64_t) + 8);
		sa = (int64_t*)malloc((nr < na? nr : na) * sizeof(int64_t) + 8);
		for (i = j = 0; i < nr && j < na; ) { // shared unique k-mers
			uint64_t ci = b->ka[i] >> 30, cj = b->kb[j] >> 30;
			if (ci < cj) ++i;
			else if (ci > cj) ++j;
			else sr[ns] = (int64_t)(b->ka[i] & 0x3fffffff), sa[ns] = (int64_t)(b->kb[j] & 0x3fffffff), ++ns, ++i, ++j;
		}
	}
	if (ns == 0) { free(sr); free(sa); if (!qa_diag(R, lr, A, la, b, st, r0) && qa_ksw(R, lr, A, la, b, st, r0, 1) != 1) qa_replace(st, lr, la, b, r0, 1); return; }
	if (!chained) { // order by reference position, then longest increasing subsequence on assembly position
		uint64_t *o = (uint64_t*)malloc(ns * sizeof(uint64_t));
		int64_t *tr = (int64_t*)malloc(ns * sizeof(int64_t)), *ta = (int64_t*)malloc(ns * sizeof(int64_t));
		for (i = 0; i < ns; ++i) o[i] = (uint64_t)sr[i] << 32 | (uint64_t)i;
		qsort(o, ns, sizeof(uint64_t), qa_cmp_u64);
		for (i = 0; i < ns; ++i) tr[i] = sr[o[i] & 0xffffffff], ta[i] = sa[o[i] & 0xffffffff];
		memcpy(sr, tr, ns * sizeof(int64_t)), memcpy(sa, ta, ns * sizeof(int64_t));
		free(o); free(tr); free(ta);
	}
	tail = (int64_t*)malloc(ns * sizeof(int64_t)), prev = (int64_t*)malloc(ns * sizeof(int64_t));
	if (chained) { // the chain is already collinear: keep every anchor
		for (i = 0; i < ns; ++i) tail[i] = i, prev[i] = i - 1;
		nl = ns;
	} else for (i = 0; i < ns; ++i) {
		int64_t lo = 0, hi = nl;
		while (lo < hi) { int64_t mid = (lo + hi) >> 1; if (sa[tail[mid]] < sa[i]) lo = mid + 1; else hi = mid; }
		prev[i] = lo > 0? tail[lo - 1] : -1;
		tail[lo] = i;
		if (lo == nl) ++nl;
	}
	keep = (int64_t*)malloc(nl * sizeof(int64_t));
	for (nk = nl, i = nl > 0? tail[nl - 1] : -1; i >= 0; i = prev[i]) keep[--nk] = i;
	// walk the seeds; merge overlapping seeds on one diagonal, drop other overlaps, recurse into gaps
	pr = 0, pa = 0;
	for (i = 0; i < nl; ++i) {
		int64_t r0 = sr[keep[i]], a0 = sa[keep[i]], r1 = r0 + sk, a1 = a0 + sk;
		if (r0 < pr || a0 < pa) {
			if (r0 - a0 == pr - pa && r1 > pr) { st->aln += r1 - pr; pr = r1, pa = a1; } // same diagonal: extend
			continue;
		}
		qa_align(R + pr, r0 - pr, A + pa, a0 - pa, depth + 1, b, st, rb + pr);
		st->aln += sk;
		pr = r1, pa = a1;
	}
	qa_align(R + pr, lr - pr, A + pa, la - pa, depth + 1, b, st, rb + pr);
	free(sr); free(sa); free(tail); free(prev); free(keep);
}

/* Extend from the start of R/A (or from the end when rev) with z-drop: SNPs and small indels are
 * crossed, extension stops where similarity collapses (an indel / rearrangement breakpoint).
 * Extension needs a positive match score; the extended pieces are re-aligned with AnchorWave's
 * scoring afterwards. Sets *er, *ea to the extended lengths. */
#define QA_EXT_MATCH 2
#define QA_EXT_ZDROP 200
#define QA_EXT_W     500
static void qa_extend(const uint8_t *R, int64_t lr, const uint8_t *A, int64_t la, int rev, int64_t *er, int64_t *ea)
{
	static const int8_t mat[25] = { QA_EXT_MATCH,-QA_KSW_MIS,-QA_KSW_MIS,-QA_KSW_MIS,-1, -QA_KSW_MIS,QA_EXT_MATCH,-QA_KSW_MIS,-QA_KSW_MIS,-1,
		-QA_KSW_MIS,-QA_KSW_MIS,QA_EXT_MATCH,-QA_KSW_MIS,-1, -QA_KSW_MIS,-QA_KSW_MIS,-QA_KSW_MIS,QA_EXT_MATCH,-1, -1,-1,-1,-1,-1 };
	int64_t i, nr = lr < INT32_MAX / 4? lr : INT32_MAX / 4, na = la < INT32_MAX / 4? la : INT32_MAX / 4;
	uint8_t *t = (uint8_t*)malloc(nr + na + 2), *qs = t + nr + 1;
	ksw_extz_t ez;
	*er = *ea = 0;
	for (i = 0; i < nr; ++i) { uint8_t c = rev? R[lr - 1 - i] : R[i]; t[i] = qa_acgt(c)? c - 1 : 4; }
	for (i = 0; i < na; ++i) { uint8_t c = rev? A[la - 1 - i] : A[i]; qs[i] = qa_acgt(c)? c - 1 : 4; }
	memset(&ez, 0, sizeof(ez));
	ksw_extd2_sse(0, (int)na, qs, (int)nr, t, 5, mat, QA_KSW_GO1, QA_KSW_GE1, QA_KSW_GO2, QA_KSW_GE2, QA_EXT_W, QA_EXT_ZDROP, 0,
				  KSW_EZ_EXTZ_ONLY | KSW_EZ_SCORE_ONLY, &ez);
	if (ez.max > 0) *er = ez.max_t + 1, *ea = ez.max_q + 1;
	kfree(0, ez.cigar); free(t);
}

static void qa_align(const uint8_t *R, int64_t lr, const uint8_t *A, int64_t la, int depth, qa_buf_t *b, qa_stat_t *st, int64_t r0)
{
	int64_t i, Lr, La, Rr, Ra;
	if (lr <= 0 && la <= 0) return;
	if (lr <= 0) { qa_event(st, la, 0, b, r0); return; }
	if (la <= 0) { qa_event(st, lr, 1, b, r0); return; }
	if (lr == la) { // gapless fast path
		int64_t mm = 0, lim = (int64_t)(QA_GAPLESS_MM * lr);
		for (i = 0; i < lr && mm <= lim; ++i)
			if (R[i] != A[i] && qa_acgt(R[i]) && qa_acgt(A[i])) ++mm;
		if (mm <= lim) {
			st->aln += lr, st->snp += mm;
			if (mm) {
				int64_t *sp = (int64_t*)malloc(mm * sizeof(int64_t)), n = 0;
				for (i = 0; i < lr; ++i)
					if (R[i] != A[i] && qa_acgt(R[i]) && qa_acgt(A[i])) {
						sp[n++] = r0 + i;
						if (b->dump) qa_var_push(&b->v, &b->n_v, &b->m_v, r0 + i, 1, 1, R[i], A[i], QA_VSNP);
					}
				qa_snp_classes(sp, n, st);
				free(sp);
			}
			return;
		}
	}
	// small and near-diagonal: one global alignment (or a replacement if that scores better)
	if (lr <= 2000 && la <= 2000 && (lr > la? lr - la : la - lr) <= QA_KSW_W) {
		if (qa_ksw(R, lr, A, la, b, st, r0, 1) == 1) return;
		if (depth > 0) { qa_replace(st, lr, la, b, r0, 3); return; }
	}
	// breakpoints: extend from both flanking anchors
	qa_extend(R, lr, A, la, 0, &Lr, &La);
	if (Lr >= lr && La >= la) { // the extension crossed the whole gap
		if (qa_ksw(R, lr, A, la, b, st, r0, 1) != 1) qa_replace(st, lr, la, b, r0, 3);
		return;
	}
	if (Lr >= lr || La >= la) { // one side used up: aligned piece, then the rest of the other side is one indel
		if (qa_ksw(R, Lr, A, La, b, st, r0, 0) != 1) qa_replace(st, Lr, La, b, r0, 3);
		if (lr > Lr) qa_event(st, lr - Lr, 1, b, r0 + Lr);
		if (la > La) qa_event(st, la - La, 0, b, r0 + Lr);
		return;
	}
	qa_extend(R + Lr, lr - Lr, A + La, la - La, 1, &Rr, &Ra);
	if (Lr > 0 && qa_ksw(R, Lr, A, La, b, st, r0, 0) != 1) qa_replace(st, Lr, La, b, r0, 3);
	if (Rr > 0 && qa_ksw(R + lr - Rr, Rr, A + la - Ra, Ra, b, st, r0 + lr - Rr, 0) != 1) qa_replace(st, Rr, Ra, b, r0 + lr - Rr, 3);
	{ // the middle, between the two breakpoints
		const uint8_t *MR = R + Lr, *MA = A + La;
		int64_t mr = lr - Lr - Rr, ma = la - La - Ra, m0 = r0 + Lr;
		if (mr <= 0 || ma <= 0) { // a clean indel (one side empty)
			if (mr > 0) qa_event(st, mr, 1, b, m0);
			if (ma > 0) qa_event(st, ma, 0, b, m0);
			return;
		}
		if (depth >= 3 || (Lr == 0 && Rr == 0 && depth > 0)) { // nothing more to find: best of alignment vs replacement
			if (qa_ksw(MR, mr, MA, ma, b, st, m0, 1) != 1) qa_replace(st, mr, ma, b, m0, depth >= 3? 2 : 1);
			return;
		}
		qa_seeded(MR, mr, MA, ma, depth, b, st, m0);
	}
}

/* Compare the gap between consecutive chain anchors u -> v (reference ascending). */
static void qa_gap(const uint8_t *rs, const uint8_t *cs, int64_t ur, int64_t uc, int64_t vr, int64_t vc, int strand, int32_t k,
				   qa_buf_t *b, qa_stat_t *st)
{
	int64_t lr = vr - (ur + k), la, i;
	const uint8_t *R = rs + ur + k, *A;
	if (strand == 0) la = vc - (uc + k), A = cs + uc + k;
	else { // reference ascending = assembly descending: assembly segment [vc + k, uc), reverse-complemented
		la = uc - (vc + k);
		if (la > b->m_rc) b->m_rc = la, b->rc = (uint8_t*)realloc(b->rc, la);
		for (i = 0; i < la; ++i) {
			uint8_t c = cs[uc - 1 - i];
			b->rc[i] = qa_acgt(c)? 5 - c : c;
		}
		A = b->rc;
	}
	if (lr < 0 || la < 0) { // anchors overlap on one side: net indel only
		int64_t d = la - lr; // = delta_asm - delta_ref
		if (d > 0) qa_event(st, d, 0, b, ur + k); else if (d < 0) qa_event(st, -d, 1, b, ur + k);
		return;
	}
	qa_align(R, lr, A, la, 0, b, st, ur + k);
}

typedef struct {
	int64_t st, en;  // reference [st,en)
	int32_t cls;
} qc_ivl_t;

typedef struct {
	int64_t n_pts, n_chain, n_chain_inv, n_chained_anchor;
	int64_t ins_n, ins_bp, ins_max, del_n, del_bp, del_max;
	int64_t tight_bp, tight_present, tight_missing;
	int64_t n_ivl, m_ivl;
	qc_ivl_t *ivl;
	qa_stat_t b;            // base-level comparison (asmqc-fa only)
} qc_acc_t;

typedef struct {
	int32_t k, stride, h, min_anchor;
	int64_t max_gap, max_gap_off, min_span, tight_gap, tight_diag, block_diag;
} qc_opt_t;

typedef struct {
	const qc_opt_t *opt;
	const rb3_liftpt_t **sp; // [sid] liftover points of each assembly sequence
	const int64_t *sn;       // [sid] their count
	const rb3_sid_t *sid;
	int32_t n_taxa, n_ref;
	int32_t *seq_taxon;     // [sid] taxon index, -1 for reference sequences
	int32_t *seq_ref;       // [sid] reference index, -1 otherwise
	int64_t *taxon_st;      // sequences of taxon t: those with seq_taxon == t
	int32_t **ref_cum;      // [r] prefix count of panel-anchorable stride positions
	int64_t *ref_nb;        // [r] number of stride bins
	qc_acc_t *acc;          // [t * n_ref + r]
	const uint8_t **seq;    // [sid] nt6 sequences; NULL = anchors only (lift mode)
} qc_ctx_t;

static pthread_mutex_t qc_acc_lock = PTHREAD_MUTEX_INITIALIZER;

typedef struct { int64_t r, c; } qc_pt_t;
typedef struct { double f; int32_t i, s; } qc_cand_t;

static int qc_cmp_lpt(const void *a, const void *b)
{
	const rb3_liftpt_t *x = (const rb3_liftpt_t*)a, *y = (const rb3_liftpt_t*)b;
	if (x->rsid != y->rsid) return x->rsid < y->rsid? -1 : 1;
	if (x->rpos != y->rpos) return x->rpos < y->rpos? -1 : 1;
	if (x->cpos != y->cpos) return x->cpos < y->cpos? -1 : 1;
	return 0;
}

static int qc_cmp_cand(const void *a, const void *b)
{
	const qc_cand_t *x = (const qc_cand_t*)a, *y = (const qc_cand_t*)b;
	if (x->f != y->f) return x->f > y->f? -1 : 1;
	if (x->i != y->i) return x->i < y->i? -1 : 1;
	return x->s - y->s;
}

static int qc_cmp_ivl(const void *a, const void *b)
{
	const qc_ivl_t *x = (const qc_ivl_t*)a, *y = (const qc_ivl_t*)b;
	if (x->st != y->st) return x->st < y->st? -1 : 1;
	return x->en < y->en? -1 : x->en > y->en? 1 : 0;
}

static const char *qc_suffix(const char *name) // part after the first '_' (the chromosome label)
{
	const char *p = strchr(name, '_');
	return p? p + 1 : name;
}

static int qc_is_chr(const char *lab) { return strncasecmp(lab, "chr", 3) == 0; }

// same chromosome label, ignoring a case-insensitive "chr" prefix and leading zeros ("chr01" == "chr1" == "Chr1")
static int qc_same_chrom(const char *a, const char *b)
{
	if (qc_is_chr(a)) a += 3;
	if (qc_is_chr(b)) b += 3;
	while (*a == '0' && a[1] != 0) ++a;
	while (*b == '0' && b[1] != 0) ++b;
	return strcasecmp(a, b) == 0;
}

static inline double qc_gap_pen(int64_t d)
{
	return d == 0? 0.0 : 2.0 + 4.0 * log2(1.0 + (double)d);
}

static void qc_add_ivl(qc_acc_t *a, int64_t st, int64_t en, int32_t cls)
{
	if (a->n_ivl == a->m_ivl) {
		a->m_ivl = a->m_ivl? a->m_ivl * 2 : 64;
		a->ivl = (qc_ivl_t*)realloc(a->ivl, a->m_ivl * sizeof(qc_ivl_t));
	}
	a->ivl[a->n_ivl].st = st, a->ivl[a->n_ivl].en = en, a->ivl[a->n_ivl].cls = cls;
	a->n_ivl++;
}

typedef struct {
	int64_t rst, ren, n, bp;  // reference span, anchor count, reference bp
	int32_t s, i_end, n_idx;  // strand (0 fwd, 1 rev), index of the last anchor
	int32_t *idx;             // anchor indices, ascending reference position
} qc_chain_t;

static __thread const qc_chain_t *g_blk; // qsort context for qc_cmp_blk
static int qc_cmp_blk(const void *a, const void *b)
{
	const qc_chain_t *x = &g_blk[*(const int32_t*)a], *y = &g_blk[*(const int32_t*)b];
	if (x->s != y->s) return x->s - y->s;
	return x->rst < y->rst? -1 : x->rst > y->rst? 1 : 0;
}

/* Chain one (assembly sequence, reference sequence) group; pts sorted by r then c. */
static void qc_chain_group(const qc_ctx_t *q, qc_acc_t *a, int32_t csid, int32_t rsid, int64_t n, const qc_pt_t *pts, qa_buf_t *qb)
{
	const qc_opt_t *o = q->opt;
	int32_t r = q->seq_ref[rsid], s, i, j, k = o->k, n_ch = 0, m_ch = 0;
	double *f[2];
	int32_t *p[2];
	uint8_t *used;
	qc_cand_t *cand;
	qc_chain_t *ch = 0, *blk = 0;
	int32_t n_blk = 0, m_blk = 0;
	int64_t bp_strand[2] = {0, 0}, nc = 0;
	int is_syn = qc_same_chrom(qc_suffix(q->sid->name[csid]), qc_suffix(q->sid->name[rsid]));
	int is_scaf = !qc_is_chr(qc_suffix(q->sid->name[csid]));
	// syntenic chains may bridge AnchorWave-scale deletions; chains to other chromosomes and
	// scaffolds use a stricter cap so sparse stray anchors are not bridged into long false spans
	int64_t max_gap = is_syn? o->max_gap : o->max_gap_off;

	a->n_pts += n;
	if (n < 2) return;
	// stage 1: blocks = anchor chains along a near-constant diagonal (|d| <= --block-diag),
	// so a lone stray anchor can never be pulled across a large gap
	for (s = 0; s < 2; ++s) {
		f[s] = (double*)malloc(n * sizeof(double));
		p[s] = (int32_t*)malloc(n * sizeof(int32_t));
		for (i = 0; i < n; ++i) {
			int64_t ci = s? -pts[i].c : pts[i].c;
			double best = k;
			int32_t bj = -1;
			for (j = i - 1; j >= 0 && i - j <= o->h; --j) {
				int64_t dr = pts[i].r - pts[j].r, dc = ci - (s? -pts[j].c : pts[j].c), d;
				double sc;
				if (dr > max_gap) break;
				if (dr <= 0 || dc <= 0 || dc > max_gap) continue;
				d = dc > dr? dc - dr : dr - dc;
				if (d > o->block_diag) continue;
				sc = f[s][j] + (dr < k || dc < k? (dr < dc? dr : dc) : k) - qc_gap_pen(d);
				if (sc > best) best = sc, bj = j;
			}
			f[s][i] = best, p[s][i] = bj;
		}
	}
	cand = (qc_cand_t*)malloc(2 * n * sizeof(qc_cand_t));
	for (s = 0; s < 2; ++s)
		for (i = 0; i < n; ++i)
			if (p[s][i] >= 0) cand[nc].f = f[s][i], cand[nc].i = i, cand[nc].s = s, ++nc;
	qsort(cand, nc, sizeof(qc_cand_t), qc_cmp_cand);
	used = (uint8_t*)calloc(n, 1);
	for (i = 0; i < nc; ++i) {
		int32_t x = cand[i].i, cnt = 0, y;
		s = cand[i].s;
		if (used[x]) continue;
		for (y = x; y >= 0 && !used[y]; y = p[s][y]) ++cnt;
		if (cnt >= o->min_anchor) {
			qc_chain_t *c;
			int32_t t = cnt;
			if (n_blk == m_blk) m_blk = m_blk? m_blk * 2 : 16, blk = (qc_chain_t*)realloc(blk, m_blk * sizeof(qc_chain_t));
			c = &blk[n_blk++];
			c->idx = (int32_t*)malloc(cnt * sizeof(int32_t));
			c->n_idx = cnt, c->s = s, c->i_end = x, c->n = cnt;
			for (y = x; t > 0; y = p[s][y]) c->idx[--t] = y;
			c->rst = pts[c->idx[0]].r, c->ren = pts[x].r + k, c->bp = c->ren - c->rst;
		}
		for (y = x; y >= 0 && !used[y]; y = p[s][y]) used[y] = 1;
	}
	// stage 2: join same-strand blocks in order across gaps up to max_gap (AnchorWave-scale indels)
	if (n_blk > 0) {
		int32_t *ord = (int32_t*)malloc(n_blk * sizeof(int32_t)), *pb = (int32_t*)malloc(n_blk * sizeof(int32_t));
		double *g = (double*)malloc(n_blk * sizeof(double));
		uint8_t *bused = (uint8_t*)calloc(n_blk, 1);
		qc_cand_t *bc = (qc_cand_t*)malloc(n_blk * sizeof(qc_cand_t));
		g_blk = blk;
		for (i = 0; i < n_blk; ++i) ord[i] = i;
		qsort(ord, n_blk, sizeof(int32_t), qc_cmp_blk); // by (strand, reference start)
		for (i = 0; i < n_blk; ++i) {
			const qc_chain_t *bi = &blk[ord[i]];
			const qc_pt_t *fi = &pts[bi->idx[0]];
			double best = (double)bi->n * k;
			int32_t bj = -1;
			for (j = i - 1; j >= 0 && i - j <= o->h; --j) {
				const qc_chain_t *bjb = &blk[ord[j]];
				const qc_pt_t *lj = &pts[bjb->idx[bjb->n_idx - 1]];
				int64_t dr, dc, d;
				double sc;
				if (bjb->s != bi->s) break;
				dr = fi->r - lj->r, dc = bi->s? lj->c - fi->c : fi->c - lj->c;
				if (dr <= 0 || dc <= 0 || dr > max_gap || dc > max_gap) continue;
				d = dc > dr? dc - dr : dr - dc;
				sc = g[j] + (double)bi->n * k - qc_gap_pen(d);
				if (sc > best) best = sc, bj = j;
			}
			g[i] = best, pb[i] = bj;
			bc[i].f = best, bc[i].i = i, bc[i].s = 0;
		}
		qsort(bc, n_blk, sizeof(qc_cand_t), qc_cmp_cand);
		for (i = 0; i < n_blk; ++i) {
			int32_t x = bc[i].i, y, nb = 0, tot = 0, t;
			qc_chain_t *c;
			if (bused[x]) continue;
			for (y = x; y >= 0 && !bused[y]; y = pb[y]) ++nb, tot += blk[ord[y]].n_idx;
			if (n_ch == m_ch) m_ch = m_ch? m_ch * 2 : 16, ch = (qc_chain_t*)realloc(ch, m_ch * sizeof(qc_chain_t));
			c = &ch[n_ch];
			c->idx = (int32_t*)malloc(tot * sizeof(int32_t));
			c->n_idx = tot, c->n = tot, c->s = blk[ord[x]].s;
			for (t = tot, y = x; nb > 0; y = pb[y], --nb) { // blocks back to front
				const qc_chain_t *b = &blk[ord[y]];
				t -= b->n_idx;
				memcpy(c->idx + t, b->idx, b->n_idx * sizeof(int32_t));
				bused[y] = 1;
			}
			c->i_end = c->idx[tot - 1];
			c->rst = pts[c->idx[0]].r, c->ren = pts[c->i_end].r + k, c->bp = c->ren - c->rst;
			if (c->bp >= o->min_span) {
				bp_strand[c->s] += c->bp;
				++n_ch;
			} else free(c->idx);
		}
		free(ord); free(pb); free(g); free(bused); free(bc);
		for (i = 0; i < n_blk; ++i) free(blk[i].idx);
		free(blk);
	}
	// classify and accumulate
	for (i = 0; i < n_ch; ++i) {
		qc_chain_t *c = &ch[i];
		int32_t cls, dom = bp_strand[1] > bp_strand[0]? 1 : 0;
		cls = !is_syn? (is_scaf? QC_SCAF : QC_OTHER) : (c->s == dom? QC_SYN_FWD : QC_SYN_INV);
		qc_add_ivl(a, c->rst, c->ren, cls);
		if (qc_chain_fp) {
			int64_t cmin = INT64_MAX, cmax = -1;
			for (j = 0; j < c->n_idx; ++j) {
				int64_t cc = pts[c->idx[j]].c;
				if (cc < cmin) cmin = cc;
				if (cc > cmax) cmax = cc;
			}
			pthread_mutex_lock(&qc_chain_lock);
			fprintf(qc_chain_fp, "%s\t%ld\t%ld\t%s\t%ld\t%ld\t%c\t%ld\t%s\n", q->sid->name[rsid], (long)c->rst, (long)c->ren,
					q->sid->name[csid], (long)cmin, (long)cmax + k, "+-"[c->s], (long)c->n, qc_cls_name[cls]);
			pthread_mutex_unlock(&qc_chain_lock);
		}
		if (cls == QC_SYN_FWD || cls == QC_SYN_INV) {
			a->n_chain++;
			if (cls == QC_SYN_INV) a->n_chain_inv++;
			a->n_chained_anchor += c->n;
			for (j = 1; j < c->n_idx; ++j) {
				const qc_pt_t *u = &pts[c->idx[j-1]], *v = &pts[c->idx[j]];
				int64_t dr = v->r - u->r, dc = c->s? u->c - v->c : v->c - u->c, d = dc - dr;
				if (d > 0) { a->ins_n++, a->ins_bp += d; if (d > a->ins_max) a->ins_max = d; }
				else if (d < 0) { a->del_n++, a->del_bp -= d; if (-d > a->del_max) a->del_max = -d; }
				if (cls == QC_SYN_FWD && dr <= o->tight_gap && (d < 0? -d : d) <= o->tight_diag) {
					int64_t b0 = u->r / o->stride + 1, b1 = v->r / o->stride; // bins strictly between
					a->tight_bp += dr;
					a->tight_present++;
					if (b1 > b0 && b1 <= q->ref_nb[r]) a->tight_missing += q->ref_cum[r][b1] - q->ref_cum[r][b0];
				}
			}
		}
	}
	// base level (asmqc-fa): like AnchorWave, every reference base gets at most one alignment.
	// Inverted syntenic chains claim their gaps first (AnchorWave aligns inversions), then forward
	// chains by anchor count; an anchor or gap already >= half claimed (per base) is skipped. Two passes:
	// anchors and balanced gaps first, then lopsided gaps (> 2 kb, one side < half the other) -- so a
	// transposed or rearranged block inside a long gap of the main chain is aligned by its own chain, and
	// the main chain's gap only deletes what is left unclaimed.
	if (q->seq && is_syn && n_ch > 0) {
		int64_t nb = q->sid->len[rsid] + 2, z;
		uint64_t *cl = (uint64_t*)calloc(nb / 64 + 1, 8);
		int32_t *ord = (int32_t*)malloc(n_ch * sizeof(int32_t)), no = 0, dom = bp_strand[1] > bp_strand[0]? 1 : 0, pass;
		for (pass = 0; pass < 2; ++pass)
			for (i = 0; i < n_ch; ++i)
				if ((pass == 0) == (ch[i].s != dom)) ord[no++] = i;
		for (z = 1; z < no; ++z) { // stable by anchor count, descending, within each pass
			int32_t x = ord[z], y = (int32_t)z;
			while (y > 0 && (ch[ord[y-1]].s != dom) == (ch[x].s != dom) && ch[ord[y-1]].n < ch[x].n) ord[y] = ord[y-1], --y;
			ord[y] = x;
		}
#define QC_CLAIMED(x) (cl[(x) >> 6] >> ((x) & 63) & 1)
#define QC_CLAIM(x)   (cl[(x) >> 6] |= 1ULL << ((x) & 63))
		{
			int pass2;
			for (pass2 = 0; pass2 < 2; ++pass2) // 0: anchors + balanced gaps; 1: lopsided gaps
			for (z = 0; z < no; ++z) {
				qc_chain_t *c = &ch[ord[z]];
				for (j = 0; j < c->n_idx; ++j) {
					const qc_pt_t *v = &pts[c->idx[j]];
					int64_t b0, b1, x, tot, got;
					if (j > 0) { // the gap before this anchor
						const qc_pt_t *u = &pts[c->idx[j-1]];
						int64_t g0 = u->r + k, g1 = v->r, la = c->s? u->c - (v->c + k) : v->c - (u->c + k), lr = g1 - g0;
						int64_t mx = lr > la? lr : la, mn = lr < la? lr : la;
						int lopsided = mx > 2000 && mn * 2 < mx;
						if (lopsided == pass2) {
							int skip = 0, partial = 0;
							if (g1 > g0) {
								tot = g1 - g0;
								for (x = g0, got = 0; x < g1; ++x) got += QC_CLAIMED(x);
								if (got * 2 >= tot) skip = 1;
								else if (pass2 && got * 10 >= tot) partial = 1;
							}
							if (partial) { // other chains aligned part of it: delete only what is left
								int64_t run = 0;
								for (x = g0; x <= g1; ++x) {
									if (x < g1 && !QC_CLAIMED(x)) { QC_CLAIM(x); ++run; continue; }
									if (run > 0) qa_event(&a->b, run, 1, qb, x - run), run = 0;
								}
								if (la > 0) qa_event(&a->b, la, 0, qb, g0);
							} else if (!skip) {
								for (x = g0; x < g1; ++x) QC_CLAIM(x);
								qa_gap(q->seq[rsid], q->seq[csid], u->r, u->c, v->r, v->c, c->s, k, qb, &a->b);
							}
						}
					}
					if (pass2 == 0) { // the anchor itself (exact match)
						b0 = v->r, b1 = v->r + k;
						for (x = b0, got = 0; x < b1; ++x) got += QC_CLAIMED(x);
						if (got * 2 < k) {
							a->b.aln += k;
							for (x = b0; x < b1; ++x) QC_CLAIM(x);
						}
					}
				}
			}
		}
#undef QC_CLAIMED
#undef QC_CLAIM
		free(cl); free(ord);
		if (qb->dump && qb->n_v > 0) {
			int64_t z2;
			static const char nt[] = "NACGTN";
			pthread_mutex_lock(&qc_chain_lock);
			for (z2 = 0; z2 < qb->n_v; ++z2) {
				const qa_var_t *w = &qb->v[z2];
				if (w->type == QA_VSNP)
					fprintf(qc_var_fp, "%s\t%s\t%ld\tSNP\t%c\t%c\n", q->sid->name[csid], q->sid->name[rsid], (long)w->pos + 1, nt[w->r < 6? w->r : 5], nt[w->a < 6? w->a : 5]);
				else if (w->type == QA_VREP)
					fprintf(qc_var_fp, "%s\t%s\t%ld\tREP\t%ld\t%ld\t%d\n", q->sid->name[csid], q->sid->name[rsid], (long)w->pos + 1, (long)w->rlen, (long)w->alen, w->r);
				else
					fprintf(qc_var_fp, "%s\t%s\t%ld\t%s\t%ld\t%ld\n", q->sid->name[csid], q->sid->name[rsid], (long)w->pos + 1, w->type == QA_VDEL? "DEL" : "INS", (long)w->rlen, (long)w->alen);
			}
			pthread_mutex_unlock(&qc_chain_lock);
			qb->n_v = 0;
		}
	}
	for (i = 0; i < n_ch; ++i) free(ch[i].idx);
	free(ch); free(used); free(cand);
	for (s = 0; s < 2; ++s) free(f[s]), free(p[s]);
}

static void qc_acc_merge(qc_acc_t *d, const qc_acc_t *x)
{
	int64_t i;
	d->n_pts += x->n_pts, d->n_chain += x->n_chain, d->n_chain_inv += x->n_chain_inv, d->n_chained_anchor += x->n_chained_anchor;
	d->ins_n += x->ins_n, d->ins_bp += x->ins_bp, d->del_n += x->del_n, d->del_bp += x->del_bp;
	if (x->ins_max > d->ins_max) d->ins_max = x->ins_max;
	if (x->del_max > d->del_max) d->del_max = x->del_max;
	d->tight_bp += x->tight_bp, d->tight_present += x->tight_present, d->tight_missing += x->tight_missing;
	for (i = 0; i < x->n_ivl; ++i) qc_add_ivl(d, x->ivl[i].st, x->ivl[i].en, x->ivl[i].cls);
	d->b.aln += x->b.aln, d->b.snp += x->b.snp, d->b.rep_n += x->b.rep_n, d->b.snp_d10 += x->b.snp_d10, d->b.snp_d40 += x->b.snp_d40;
	d->b.ins_n += x->b.ins_n, d->b.ins_bp += x->b.ins_bp, d->b.del_n += x->b.del_n, d->b.del_bp += x->b.del_bp;
	if (x->b.ins_max > d->b.ins_max) d->b.ins_max = x->b.ins_max;
	if (x->b.del_max > d->b.del_max) d->b.del_max = x->b.del_max;
	if (x->b.hist) {
		if (d->b.hist == 0) d->b.hist = (int64_t*)calloc(2 * QA_HIST, sizeof(int64_t));
		for (i = 0; i < 2 * QA_HIST; ++i) d->b.hist[i] += x->b.hist[i];
	}
}

static void qc_worker(void *data, long sx, int tid) // one assembly sequence per job
{
	qc_ctx_t *q = (qc_ctx_t*)data;
	int64_t n = q->sn[sx], i, i0;
	int32_t t = q->seq_taxon[sx];
	rb3_liftpt_t *tmp;
	qc_pt_t *pts;
	qa_buf_t qb;
	(void)tid;
	if (t < 0 || n == 0) return;
	memset(&qb, 0, sizeof(qb));
	qb.dump = (qc_var_fp != 0 && q->seq != 0);
	tmp = (rb3_liftpt_t*)malloc(n * sizeof(rb3_liftpt_t));
	pts = (qc_pt_t*)malloc(n * sizeof(qc_pt_t));
	memcpy(tmp, q->sp[sx], n * sizeof(rb3_liftpt_t));
	qsort(tmp, n, sizeof(rb3_liftpt_t), qc_cmp_lpt);
	for (i0 = 0, i = 1; i <= n; ++i) {
		if (i == n || tmp[i].rsid != tmp[i0].rsid) {
			int64_t j;
			int32_t r = q->seq_ref[tmp[i0].rsid];
			if (r >= 0) {
				qc_acc_t loc;
				memset(&loc, 0, sizeof(loc));
				for (j = i0; j < i; ++j) pts[j - i0].r = tmp[j].rpos, pts[j - i0].c = tmp[j].cpos;
				qc_chain_group(q, &loc, (int32_t)sx, tmp[i0].rsid, i - i0, pts, &qb);
				pthread_mutex_lock(&qc_acc_lock);
				qc_acc_merge(&q->acc[(int64_t)t * q->n_ref + r], &loc);
				pthread_mutex_unlock(&qc_acc_lock);
				free(loc.ivl); free(loc.b.hist);
			}
			i0 = i;
		}
	}
	free(tmp); free(pts);
	qa_buf_free(&qb);
}

static int64_t qc_union(const qc_acc_t *a, uint32_t cls_mask, qc_ivl_t *buf)
{
	int64_t i, n = 0, tot = 0, st = -1, en = -1;
	for (i = 0; i < a->n_ivl; ++i)
		if (cls_mask >> a->ivl[i].cls & 1) buf[n++] = a->ivl[i];
	qsort(buf, n, sizeof(qc_ivl_t), qc_cmp_ivl);
	for (i = 0; i < n; ++i) {
		if (buf[i].st > en) {
			if (en > st) tot += en - st;
			st = buf[i].st, en = buf[i].en;
		} else if (buf[i].en > en) en = buf[i].en;
	}
	if (en > st) tot += en - st;
	return tot;
}

static void qc_print_row(FILE *fp, const char *taxon, const char *chrom, int64_t len, int64_t n_anch,
						 const qc_acc_t *a, int64_t mapped, int64_t inv, int64_t other, int64_t scaf, int32_t k)
{
	double fi = a->tight_present + a->tight_missing > 0? (double)a->tight_present / (a->tight_present + a->tight_missing) : 0.0;
	double div = fi > 0.0? 1.0 - pow(fi, 1.0 / k) : 1.0;
	const qa_stat_t *bs = &a->b;
	int64_t q1 = 0, q2 = 0, q3 = 0, ni = bs->ins_n + bs->del_n, cum = 0, z;
	if (bs->hist && ni > 0) { // indel-size quartiles over insertions and deletions together
		for (z = 1; z < QA_HIST; ++z) {
			int64_t c0 = cum;
			cum += bs->hist[z] + bs->hist[QA_HIST + z];
			if (!q1 && cum * 4 >= ni) q1 = z;
			if (!q2 && cum * 2 >= ni) q2 = z;
			if (!q3 && cum * 4 >= 3 * ni) q3 = z;
			(void)c0;
		}
	}
	fprintf(fp, "%s\t%s\t%ld\t%ld\t%ld\t%ld\t%.6f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.6f\t%ld\t%ld\t%.6f\t%.6f"
			"\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.6f\t%.6f\t%.6f\t%.3f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
			taxon, chrom, (long)len, (long)n_anch, (long)a->n_pts, (long)mapped, len? (double)mapped / len : 0.0,
			(long)inv, (long)other, (long)scaf, (long)(len - mapped), (long)a->n_chain, (long)a->n_chain_inv,
			(long)a->ins_n, (long)a->ins_bp, (long)a->ins_max, (long)a->del_n, (long)a->del_bp, (long)a->del_max,
			(long)a->tight_bp, len? (double)a->tight_bp / len : 0.0, (long)a->tight_present, (long)a->tight_missing, fi, div,
			(long)bs->snp, (long)ni, (long)bs->ins_n, (long)bs->del_n, (long)bs->ins_bp, (long)bs->del_bp,
			len? (double)bs->aln / len : 0.0, len? (double)bs->del_bp / len : 0.0, len? (double)bs->ins_bp / len : 0.0,
			ni? (double)(bs->ins_bp + bs->del_bp) / ni : 0.0, (long)q1, (long)q2, (long)q3,
			(long)(bs->ins_max > bs->del_max? bs->ins_max : bs->del_max), (long)bs->rep_n, (long)bs->snp_d10, (long)bs->snp_d40);
}


/* Shared driver: sid = all sequences (reference ones start with ref_prefix), sp/sn = the
 * liftover points of each assembly sequence. any_in (optional, indexed by reference order)
 * marks the anchorable stride bins used as the divergence denominator; NULL derives it from
 * the points (panel mode: a bin is anchorable if any assembly carries its k-mer). */
static int qc_run(const qc_opt_t *opt, const rb3_sid_t *sid, const char *ref_prefix,
				  const rb3_liftpt_t **sp, const int64_t *sn, uint8_t **any_in, const uint8_t **seqs, int n_threads, FILE *fp, int print_header)
{
	qc_ctx_t q;
	char **taxa = 0;
	int32_t m_taxa = 0, t, r;
	int64_t sx, plen, i, max_ivl = 0;
	qc_ivl_t *buf;
	uint8_t **any;

	memset(&q, 0, sizeof(q));
	q.opt = opt, q.sid = sid, q.sp = sp, q.sn = sn, q.seq = seqs;
	q.seq_taxon = (int32_t*)malloc(sid->n_seq * sizeof(int32_t));
	q.seq_ref = (int32_t*)malloc(sid->n_seq * sizeof(int32_t));
	plen = strlen(ref_prefix);
	for (sx = 0; sx < sid->n_seq; ++sx) {
		const char *nm = sid->name[sx], *u = strchr(nm, '_');
		int32_t l = u? (int32_t)(u - nm) : (int32_t)strlen(nm);
		q.seq_taxon[sx] = q.seq_ref[sx] = -1;
		if (strncmp(nm, ref_prefix, plen) == 0) { q.seq_ref[sx] = q.n_ref++; continue; }
		for (t = q.n_taxa - 1; t >= 0; --t) // taxa are contiguous in practice; scan back
			if ((int32_t)strlen(taxa[t]) == l && strncmp(taxa[t], nm, l) == 0) break;
		if (t < 0) {
			if (q.n_taxa == m_taxa) m_taxa = m_taxa? m_taxa * 2 : 32, taxa = (char**)realloc(taxa, m_taxa * sizeof(char*));
			taxa[q.n_taxa] = (char*)calloc(l + 1, 1);
			memcpy(taxa[q.n_taxa], nm, l);
			t = q.n_taxa++;
		}
		q.seq_taxon[sx] = t;
	}
	if (q.n_ref == 0) { fprintf(stderr, "ERROR: no sequence starts with '%s'\n", ref_prefix); return 1; }

	q.ref_nb = (int64_t*)calloc(q.n_ref, sizeof(int64_t));
	q.ref_cum = (int32_t**)calloc(q.n_ref, sizeof(int32_t*));
	for (sx = 0; sx < sid->n_seq; ++sx)
		if (q.seq_ref[sx] >= 0) q.ref_nb[q.seq_ref[sx]] = sid->len[sx] / opt->stride + 1;
	if (any_in == 0) { // panel-anchorable reference stride bins: any assembly carries the anchor
		int64_t n_off = 0, n_all = 0;
		any = (uint8_t**)calloc(q.n_ref, sizeof(uint8_t*));
		for (r = 0; r < q.n_ref; ++r) any[r] = (uint8_t*)calloc(q.ref_nb[r], 1);
		for (sx = 0; sx < sid->n_seq; ++sx) {
			const rb3_liftpt_t *lp = sp[sx];
			for (i = 0; i < sn[sx]; ++i) {
				int32_t rr = lp[i].rsid >= 0 && lp[i].rsid < sid->n_seq? q.seq_ref[lp[i].rsid] : -1;
				if (rr < 0) continue;
				++n_all;
				if (lp[i].rpos % opt->stride != 0) ++n_off;
				if (lp[i].rpos / opt->stride < q.ref_nb[rr]) any[rr][lp[i].rpos / opt->stride] = 1;
			}
		}
		if (n_off > 0)
			fprintf(stderr, "[W::%s] %ld of %ld points are off the -s %d grid; is -s the stride the map was built with?\n",
					__func__, (long)n_off, (long)n_all, opt->stride);
		fprintf(stderr, "[M::%s] %ld liftover points, %d assemblies, %d reference sequences\n", __func__, (long)n_all, q.n_taxa, q.n_ref);
	} else any = any_in;
	for (r = 0; r < q.n_ref; ++r) {
		q.ref_cum[r] = (int32_t*)malloc((q.ref_nb[r] + 1) * sizeof(int32_t));
		q.ref_cum[r][0] = 0;
		for (i = 0; i < q.ref_nb[r]; ++i) q.ref_cum[r][i + 1] = q.ref_cum[r][i] + (any[r][i] != 0);
	}
	if (any_in == 0) {
		for (r = 0; r < q.n_ref; ++r) free(any[r]);
		free(any);
	}

	q.acc = (qc_acc_t*)calloc((int64_t)q.n_taxa * q.n_ref, sizeof(qc_acc_t));
	kt_for(n_threads, qc_worker, &q, sid->n_seq);

	// report
	for (i = 0; i < (int64_t)q.n_taxa * q.n_ref; ++i)
		if (q.acc[i].n_ivl > max_ivl) max_ivl = q.acc[i].n_ivl;
	buf = (qc_ivl_t*)malloc((max_ivl + 1) * sizeof(qc_ivl_t));
	if (print_header) fprintf(fp, "taxon\tchrom\tref_len\tn_anchorable\tn_pts\tmapped_bp\tfrac_mapped\tinv_bp\tother_chr_bp\tscaf_bp\tgap_bp"
		   "\tn_chain\tn_chain_inv\tins_n\tins_bp_lb\tins_max\tdel_n\tdel_bp_lb\tdel_max"
		   "\ttight_bp\tfrac_tight\ttight_present\ttight_missing\tfrac_intact\tdiv_est"
		   "\tnumSNPs\tnumIndels\tnumIns\tnumDel\tnumBasesInserted\tnumBasesDeleted\tpercentRefAligned\tpercentRefDeleted"
		   "\tpercentRefInserted\tmeanIndelSize\tlowerQuartileIndelSize\tmedianIndelSize\tupperQuartileIndelSize\tlargestIndel\tnumReplacements"
		   "\tnumSNPs_localMismatch10\tnumSNPs_localMismatch40\n");
	for (t = 0; t < q.n_taxa; ++t) {
		qc_acc_t tot;
		int64_t tl = 0, tanch = 0, tm = 0, tinv = 0, toth = 0, tsc = 0;
		memset(&tot, 0, sizeof(tot));
		for (sx = 0; sx < sid->n_seq; ++sx) {
			qc_acc_t *a;
			int64_t m, inv, oth, sc;
			if ((r = q.seq_ref[sx]) < 0) continue;
			a = &q.acc[(int64_t)t * q.n_ref + r];
			m = qc_union(a, 1u<<QC_SYN_FWD | 1u<<QC_SYN_INV, buf);
			inv = qc_union(a, 1u<<QC_SYN_INV, buf);
			oth = qc_union(a, 1u<<QC_OTHER, buf);
			sc = qc_union(a, 1u<<QC_SCAF, buf);
			qc_print_row(fp, taxa[t], sid->name[sx], sid->len[sx], q.ref_cum[r][q.ref_nb[r]], a, m, inv, oth, sc, opt->k);
			tl += sid->len[sx], tanch += q.ref_cum[r][q.ref_nb[r]], tm += m, tinv += inv, toth += oth, tsc += sc;
			{ qc_acc_t tmp_a = *a; tmp_a.n_ivl = 0; qc_acc_merge(&tot, &tmp_a); } // intervals are not needed for ALL
		}
		qc_print_row(fp, taxa[t], "ALL", tl, tanch, &tot, tm, tinv, toth, tsc, opt->k);
		free(tot.ivl); free(tot.b.hist);
	}
	free(buf);

	for (i = 0; i < (int64_t)q.n_taxa * q.n_ref; ++i) free(q.acc[i].ivl), free(q.acc[i].b.hist);
	for (r = 0; r < q.n_ref; ++r) free(q.ref_cum[r]);
	for (t = 0; t < q.n_taxa; ++t) free(taxa[t]);
	free(q.acc); free(q.ref_cum); free(q.ref_nb); free(taxa);
	free(q.seq_taxon); free(q.seq_ref);
	return 0;
}

static ko_longopt_t qc_long_options[] = {
	{ "ref-prefix",  ko_required_argument, 301 },
	{ "max-gap",     ko_required_argument, 302 },
	{ "min-span",    ko_required_argument, 303 },
	{ "tight-gap",   ko_required_argument, 304 },
	{ "tight-diag",  ko_required_argument, 305 },
	{ "max-gap-off", ko_required_argument, 306 },
	{ "ref-name",    ko_required_argument, 307 },
	{ "chains",      ko_required_argument, 308 },
	{ "block-diag",  ko_required_argument, 309 },
	{ "variants",    ko_required_argument, 310 },
	{ 0, 0, 0 }
};

static void qc_opt_init(qc_opt_t *opt)
{
	memset(opt, 0, sizeof(qc_opt_t));
	opt->k = 100, opt->stride = 2000, opt->h = 50, opt->min_anchor = 3;
	opt->max_gap = 3000000, opt->max_gap_off = 1000000, opt->min_span = 1000, opt->tight_gap = 1000, opt->tight_diag = 20, opt->block_diag = 10000;
}

// returns 0 if the option was consumed, -1 otherwise
static int qc_opt_parse1(qc_opt_t *opt, int c, const char *arg, int *n_threads, char **ref_prefix)
{
	if (c == 'k') opt->k = atoi(arg);
	else if (c == 's') opt->stride = (int32_t)rb3_parse_num(arg);
	else if (c == 't') *n_threads = atoi(arg);
	else if (c == 'h') opt->h = atoi(arg);
	else if (c == 'n') opt->min_anchor = atoi(arg);
	else if (c == 301 || c == 307) *ref_prefix = (char*)arg;
	else if (c == 302) opt->max_gap = rb3_parse_num(arg);
	else if (c == 303) opt->min_span = rb3_parse_num(arg);
	else if (c == 304) opt->tight_gap = rb3_parse_num(arg);
	else if (c == 305) opt->tight_diag = rb3_parse_num(arg);
	else if (c == 306) opt->max_gap_off = rb3_parse_num(arg);
	else if (c == 309) opt->block_diag = rb3_parse_num(arg);
	else if (c == 310) {
		qc_var_fp = fopen(arg, "w");
		if (qc_var_fp == 0) { fprintf(stderr, "ERROR: cannot write '%s'\n", arg); exit(1); }
		fprintf(qc_var_fp, "#asm_seq\tref_seq\tpos(1-based)\ttype\tref|ref_len\talt|asm_len\n");
	}
	else if (c == 308) {
		qc_chain_fp = fopen(arg, "w");
		if (qc_chain_fp == 0) { fprintf(stderr, "ERROR: cannot write '%s'\n", arg); exit(1); }
		fprintf(qc_chain_fp, "ref_seq\tref_start\tref_end\tasm_seq\tasm_start\tasm_end\tstrand\tn_anchor\tclass\n");
	}
	else return -1;
	return 0;
}

static void qc_opt_usage(const qc_opt_t *opt, int n_threads)
{
	fprintf(stderr, "  -k INT          anchor k-mer length [%d]\n", opt->k);
	fprintf(stderr, "  -s NUM          anchor stride on the reference [%d]\n", opt->stride);
	fprintf(stderr, "  -t INT          threads [%d]\n", n_threads);
	fprintf(stderr, "  -h INT          chaining look-back in anchors [%d]\n", opt->h);
	fprintf(stderr, "  -n INT          min anchors per chain [%d]\n", opt->min_anchor);
	fprintf(stderr, "  --block-diag NUM  max |asm gap - ref gap| between anchors of one block [%ld]\n", (long)opt->block_diag);
	fprintf(stderr, "  --max-gap NUM   max reference/assembly gap between blocks of a syntenic chain [%ld]\n", (long)opt->max_gap);
	fprintf(stderr, "  --max-gap-off NUM  same, for chains to other chromosomes / scaffolds [%ld]\n", (long)opt->max_gap_off);
	fprintf(stderr, "  --min-span NUM  min reference span of a kept chain [%ld]\n", (long)opt->min_span);
	fprintf(stderr, "  --chains FILE   also write every kept chain (reference span, assembly span, strand, class)\n");
	fprintf(stderr, "  --variants FILE (asmqc-fa) also write every SNP and replaced segment on the reference\n");
	fprintf(stderr, "  --tight-gap NUM / --tight-diag NUM  collinear-interval limits for the divergence estimate [%ld/%ld]\n",
			(long)opt->tight_gap, (long)opt->tight_diag);
}

int main_asmqc(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	qc_opt_t opt;
	char *ref_prefix = 0;
	int32_t c, n_threads = 4, ret;
	int64_t sx;
	rb3_lift_t *lf;
	rb3_sid_t *sid;
	const rb3_liftpt_t **sp;
	int64_t *sn;

	qc_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:s:t:h:n:", qc_long_options)) >= 0)
		qc_opt_parse1(&opt, c, o.arg, &n_threads, &ref_prefix);
	if (argc - o.ind < 2 || ref_prefix == 0) {
		fprintf(stderr, "Usage: ropebwt3 asmqc --ref-prefix=STR [options] <idx.lift> <idx.fmd.len.gz>\n");
		fprintf(stderr, "Options (k and stride are NOT stored in the .lift; pass the values used by `ropebwt3 lift`):\n");
		qc_opt_usage(&opt, n_threads);
		return 1;
	}
	lf = rb3_lift_restore(argv[o.ind]);
	if (lf == 0) { fprintf(stderr, "ERROR: failed to read liftover map '%s'\n", argv[o.ind]); return 1; }
	sid = rb3_sid_read(argv[o.ind + 1]);
	if (sid == 0) { fprintf(stderr, "ERROR: failed to read sequence list '%s'\n", argv[o.ind + 1]); return 1; }
	if (sid->n_seq != rb3_lift_n_seq(lf)) {
		fprintf(stderr, "ERROR: %ld sequences in '%s' but %ld in the liftover map\n", (long)sid->n_seq, argv[o.ind + 1], (long)rb3_lift_n_seq(lf));
		return 1;
	}
	sp = (const rb3_liftpt_t**)malloc(sid->n_seq * sizeof(rb3_liftpt_t*));
	sn = (int64_t*)malloc(sid->n_seq * sizeof(int64_t));
	for (sx = 0; sx < sid->n_seq; ++sx) sp[sx] = rb3_lift_seq_points(lf, sx, &sn[sx]);
	ret = qc_run(&opt, sid, ref_prefix, sp, sn, 0, 0, n_threads, stdout, 1);
	if (qc_chain_fp) fclose(qc_chain_fp), qc_chain_fp = 0;
	free(sp); free(sn);
	rb3_sid_destroy(sid);
	rb3_lift_destroy(lf);
	return ret;
}

/* ---- asmqc-fa: pairwise mode straight from FASTA, no index build ---------
 *
 * Reproduces the anchors `ropebwt3 lift` would produce for an index holding just the
 * reference and one assembly: every reference k-mer at offset i*stride (no N) that occurs
 * exactly once in the reference (both strands) and exactly once in the assembly (both
 * strands) becomes a point (cpos = forward-strand start in the assembly). Occurrences are
 * found with a polynomial rolling hash over both strands and verified base by base. The
 * reference anchor table and its uniqueness are computed once and reused for every
 * assembly given on the command line. The divergence denominator is "unique in the
 * reference" (there is no panel), so div_est is not directly comparable to panel mode.
 */

#define QF_B 0x9E3779B97F4A7C15ULL

typedef struct {
	uint64_t h;           // forward hash of the reference k-mer
	int64_t rpos, cpos;   // reference / assembly forward-strand start
	int32_t rsid, csid;   // sequence index (combined name space)
	uint32_t rocc, cocc;  // occurrences in reference / assembly (both strands)
} qf_anc_t;

typedef struct {
	int32_t n, m;
	char **name;          // "<taxon>_<label>"
	int64_t *len;
	uint8_t **seq;        // nt6 (1..4 = ACGT, else break)
} qf_seqs_t;

typedef struct {
	int32_t k;
	uint64_t bk1, binv;   // B^(k-1), B^-1 (mod 2^64)
	uint64_t mask;
	int32_t *slot;        // open-addressing table of anchor indices, -1 = empty
	qf_anc_t *a;
	int64_t n;
	const qf_seqs_t *ref; // reference sequences (rsid = index in ref)
	const qf_seqs_t *cur; // sequences being scanned
	int32_t csid_off;     // combined-space index of cur->seq[0]
	int is_ref_scan;
} qf_tab_t;

static uint64_t qf_hash_fwd(const uint8_t *s, int32_t k)
{
	uint64_t h = 0;
	int32_t j;
	for (j = 0; j < k; ++j) h = h * QF_B + s[j];
	return h;
}

static inline int32_t qf_find(const qf_tab_t *T, uint64_t h, int64_t from)
{ // next anchor with hash h, starting at probe slot `from`; returns slot or -1
	int64_t x;
	for (x = from & T->mask;; x = (x + 1) & T->mask) {
		int32_t e = T->slot[x];
		if (e < 0) return -1;
		if (T->a[e].h == h) return (int32_t)x;
	}
}

static void qf_hit(qf_tab_t *T, int32_t sidx, const uint8_t *s, int64_t st, uint64_t h, int is_rc)
{
	int64_t x = (int64_t)(h ^ (h >> 29));
	int32_t sl;
	for (sl = qf_find(T, h, x); sl >= 0; sl = qf_find(T, h, sl + 1)) {
		qf_anc_t *e = &T->a[T->slot[sl]];
		const uint8_t *r = T->ref->seq[e->rsid] + e->rpos;
		int32_t j, k = T->k, ok = 1;
		if (!is_rc) { for (j = 0; j < k; ++j) if (s[st + j] != r[j]) { ok = 0; break; } }
		else { for (j = 0; j < k; ++j) if (s[st + j] != 5 - r[k - 1 - j]) { ok = 0; break; } }
		if (!ok) continue;
		if (T->is_ref_scan) __sync_fetch_and_add(&e->rocc, 1);
		else if (__sync_fetch_and_add(&e->cocc, 1) == 0) e->cpos = st, e->csid = T->csid_off + sidx;
	}
}

static void qf_scan_worker(void *data, long sidx, int tid)
{
	qf_tab_t *T = (qf_tab_t*)data;
	const uint8_t *s = T->cur->seq[sidx];
	int64_t i, len = T->cur->len[sidx];
	int32_t k = T->k, l = 0;
	uint64_t h = 0, hr = 0, p = 1;
	(void)tid;
	for (i = 0; i < len; ++i) {
		uint8_t c = s[i];
		if (c < 1 || c > 4) { l = 0, h = hr = 0, p = 1; continue; }
		if (l < k) {
			h = h * QF_B + c;
			hr += (uint64_t)(5 - c) * p, p *= QF_B;
			++l;
		} else {
			uint8_t o = s[i - k];
			h = (h - o * T->bk1) * QF_B + c;
			hr = (hr - (uint64_t)(5 - o)) * T->binv + (uint64_t)(5 - c) * T->bk1;
		}
		if (l == k) {
			qf_hit(T, (int32_t)sidx, s, i - k + 1, h, 0);
			qf_hit(T, (int32_t)sidx, s, i - k + 1, hr, 1);
		}
	}
}

static void qf_scan(qf_tab_t *T, const qf_seqs_t *ss, int32_t csid_off, int is_ref, int n_threads)
{
	T->cur = ss, T->csid_off = csid_off, T->is_ref_scan = is_ref;
	kt_for(n_threads, qf_scan_worker, T, ss->n);
}

static void qf_build_table(qf_tab_t *T, const qf_seqs_t *ref, int32_t k, int32_t stride)
{
	int64_t m = 0, i, cap = 1;
	int32_t r;
	uint64_t x;
	T->k = k, T->ref = ref;
	for (r = 0; r < ref->n; ++r) if (ref->len[r] >= k) m += (ref->len[r] - k) / stride + 1;
	T->a = (qf_anc_t*)calloc(m > 0? m : 1, sizeof(qf_anc_t));
	for (r = 0; r < ref->n; ++r) {
		for (i = 0; i + k <= ref->len[r]; i += stride) {
			const uint8_t *q = ref->seq[r] + i;
			int32_t j;
			for (j = 0; j < k; ++j) if (q[j] < 1 || q[j] > 4) break;
			if (j < k) continue;
			T->a[T->n].h = qf_hash_fwd(q, k), T->a[T->n].rpos = i, T->a[T->n].rsid = r, T->a[T->n].csid = -1;
			T->n++;
		}
	}
	while (cap < T->n * 2) cap <<= 1;
	T->mask = cap - 1;
	T->slot = (int32_t*)malloc(cap * sizeof(int32_t));
	for (i = 0; i < cap; ++i) T->slot[i] = -1;
	for (i = 0; i < T->n; ++i) {
		uint64_t h = T->a[i].h;
		for (x = (h ^ (h >> 29)) & T->mask; T->slot[x] >= 0; x = (x + 1) & T->mask) {}
		T->slot[x] = (int32_t)i;
	}
	for (x = 1, i = 0; i < k - 1; ++i) x *= QF_B;
	T->bk1 = x;
	for (x = QF_B, i = 0; i < 6; ++i) x *= 2 - QF_B * x; // Newton: inverse mod 2^64
	T->binv = x;
}

static char *qf_taxon_from_path(const char *fn)
{
	const char *b = strrchr(fn, '/');
	char *t, *p;
	b = b? b + 1 : fn;
	t = strdup(b);
	if ((p = strchr(t, '.')) != 0) *p = 0;
	if ((p = strstr(t, "_renamed")) != 0 && p[8] == 0) *p = 0;
	for (p = t; *p; ++p) if (*p == '_') *p = '-'; // '_' separates taxon from chromosome label
	return t;
}

static int qf_load(const char *fn, const char *taxon, qf_seqs_t *ss)
{
	rb3_seqio_t *fp = rb3_seq_open(fn, 0);
	const char *nm;
	char *s;
	int64_t len;
	size_t tl = strlen(taxon);
	if (fp == 0) return -1;
	while ((s = rb3_seq_read1(fp, &len, &nm)) != 0) {
		const char *lab = nm;
		char *full;
		// strip an existing "<taxon>_" prefix (case of PHG-renamed FASTAs)
		{
			const char *u = strchr(nm, '_');
			if (u && (size_t)(u - nm) == tl) {
				int32_t j, same = 1;
				for (j = 0; j < (int32_t)tl; ++j)
					if (nm[j] != taxon[j] && !(nm[j] == '_' && taxon[j] == '-')) { same = 0; break; }
				if (same) lab = u + 1;
			}
		}
		if (ss->n == ss->m) {
			ss->m = ss->m? ss->m * 2 : 64;
			ss->name = (char**)realloc(ss->name, ss->m * sizeof(char*));
			ss->len = (int64_t*)realloc(ss->len, ss->m * sizeof(int64_t));
			ss->seq = (uint8_t**)realloc(ss->seq, ss->m * sizeof(uint8_t*));
		}
		full = (char*)malloc(tl + strlen(lab) + 2);
		sprintf(full, "%s_%s", taxon, lab);
		ss->name[ss->n] = full, ss->len[ss->n] = len;
		ss->seq[ss->n] = (uint8_t*)malloc(len + 1);
		memcpy(ss->seq[ss->n], s, len);
		rb3_char2nt6(len, ss->seq[ss->n]);
		ss->n++;
	}
	rb3_seq_close(fp);
	return 0;
}

static void qf_seqs_free(qf_seqs_t *ss)
{
	int32_t i;
	for (i = 0; i < ss->n; ++i) free(ss->name[i]), free(ss->seq[i]);
	free(ss->name); free(ss->len); free(ss->seq);
	memset(ss, 0, sizeof(*ss));
}

int main_asmqc_fa(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	qc_opt_t opt;
	char *ref_name = 0, *ref_taxon, *ref_pfx;
	int32_t c, n_threads = 8, a, r, ret = 0, header_done = 0;
	int64_t i, n_uniq = 0;
	qf_seqs_t ref;
	qf_tab_t T;
	uint8_t **any;
	double t0 = rb3_realtime();

	qc_opt_init(&opt);
	opt.stride = 200;
	while ((c = ketopt(&o, argc, argv, 1, "k:s:t:h:n:", qc_long_options)) >= 0)
		qc_opt_parse1(&opt, c, o.arg, &n_threads, &ref_name);
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: ropebwt3 asmqc-fa [options] <ref.fa> <asm1.fa> [asm2.fa ...]\n");
		fprintf(stderr, "Pairwise QC of each assembly against the reference straight from FASTA (no index build).\n");
		fprintf(stderr, "Taxon names come from the file names (text before the first '.', '_renamed' dropped).\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  --ref-name STR  reference taxon name [from the file name]\n");
		qc_opt_usage(&opt, n_threads);
		return 1;
	}
	if (opt.k > 1000 || opt.k < 16) { fprintf(stderr, "ERROR: -k must be in [16,1000]\n"); return 1; }
	ref_taxon = ref_name? strdup(ref_name) : qf_taxon_from_path(argv[o.ind]);
	ref_pfx = (char*)malloc(strlen(ref_taxon) + 2);
	sprintf(ref_pfx, "%s_", ref_taxon);
	memset(&ref, 0, sizeof(ref));
	if (qf_load(argv[o.ind], ref_taxon, &ref) < 0 || ref.n == 0) { fprintf(stderr, "ERROR: failed to read '%s'\n", argv[o.ind]); return 1; }
	fprintf(stderr, "[M::%s::%.1f] loaded reference %s: %d sequences\n", __func__, rb3_realtime() - t0, ref_taxon, ref.n);

	memset(&T, 0, sizeof(T));
	qf_build_table(&T, &ref, opt.k, opt.stride);
	qf_scan(&T, &ref, 0, 1, n_threads);
	any = (uint8_t**)calloc(ref.n, sizeof(uint8_t*));
	for (r = 0; r < ref.n; ++r) any[r] = (uint8_t*)calloc(ref.len[r] / opt.stride + 1, 1);
	for (i = 0; i < T.n; ++i)
		if (T.a[i].rocc == 1) any[T.a[i].rsid][T.a[i].rpos / opt.stride] = 1, ++n_uniq;
	fprintf(stderr, "[M::%s::%.1f] %ld stride %d-mers, %ld unique in the reference\n", __func__, rb3_realtime() - t0, (long)T.n, opt.k, (long)n_uniq);

	for (a = o.ind + 1; a < argc; ++a) {
		qf_seqs_t as;
		char *taxon = qf_taxon_from_path(argv[a]);
		rb3_sid_t sid;
		const rb3_liftpt_t **sp;
		int64_t *sn, np = 0;
		rb3_liftpt_t *pt;
		int32_t s, n_all;

		memset(&as, 0, sizeof(as));
		if (strcmp(taxon, ref_taxon) == 0) { fprintf(stderr, "ERROR: assembly '%s' has the reference's taxon name\n", argv[a]); ret = 1; free(taxon); continue; }
		if (qf_load(argv[a], taxon, &as) < 0 || as.n == 0) { fprintf(stderr, "ERROR: failed to read '%s'\n", argv[a]); ret = 1; free(taxon); continue; }
		for (i = 0; i < T.n; ++i) T.a[i].cocc = 0, T.a[i].csid = -1;
		qf_scan(&T, &as, ref.n, 0, n_threads);

		// combined name space: reference sequences, then this assembly's
		n_all = ref.n + as.n;
		memset(&sid, 0, sizeof(sid));
		sid.n_seq = n_all;
		sid.name = (char**)malloc(n_all * sizeof(char*));
		sid.len = (int32_t*)malloc(n_all * sizeof(int32_t));
		for (s = 0; s < ref.n; ++s) sid.name[s] = ref.name[s], sid.len[s] = (int32_t)ref.len[s];
		for (s = 0; s < as.n; ++s) sid.name[ref.n + s] = as.name[s], sid.len[ref.n + s] = (int32_t)as.len[s];
		sn = (int64_t*)calloc(n_all, sizeof(int64_t));
		for (i = 0; i < T.n; ++i)
			if (T.a[i].rocc == 1 && T.a[i].cocc == 1) sn[T.a[i].csid]++, ++np;
		pt = (rb3_liftpt_t*)malloc((np > 0? np : 1) * sizeof(rb3_liftpt_t));
		sp = (const rb3_liftpt_t**)malloc(n_all * sizeof(rb3_liftpt_t*));
		{
			int64_t *fill = (int64_t*)calloc(n_all + 1, sizeof(int64_t));
			for (s = 0; s < n_all; ++s) fill[s + 1] = fill[s] + sn[s];
			for (s = 0; s < n_all; ++s) sp[s] = pt + fill[s];
			for (i = 0; i < T.n; ++i)
				if (T.a[i].rocc == 1 && T.a[i].cocc == 1) {
					rb3_liftpt_t *p = &pt[fill[T.a[i].csid]++];
					p->cpos = T.a[i].cpos, p->rpos = T.a[i].rpos, p->csid = T.a[i].csid, p->rsid = T.a[i].rsid;
				}
			free(fill);
		}
		fprintf(stderr, "[M::%s::%.1f] %s: %d sequences, %ld anchors unique in both\n", __func__, rb3_realtime() - t0, taxon, as.n, (long)np);
		{
			const uint8_t **seqs = (const uint8_t**)malloc(n_all * sizeof(uint8_t*));
			for (s = 0; s < ref.n; ++s) seqs[s] = ref.seq[s];
			for (s = 0; s < as.n; ++s) seqs[ref.n + s] = as.seq[s];
			ret |= qc_run(&opt, &sid, ref_pfx, sp, sn, any, seqs, n_threads, stdout, !header_done);
			free((void*)seqs);
		}
		fflush(stdout);
		header_done = 1;
		free(pt); free(sp); free(sn); free(sid.name); free(sid.len);
		qf_seqs_free(&as);
		free(taxon);
	}
	fprintf(stderr, "[M::%s::%.1f] done\n", __func__, rb3_realtime() - t0);
	if (qc_chain_fp) fclose(qc_chain_fp), qc_chain_fp = 0;
	if (qc_var_fp) fclose(qc_var_fp), qc_var_fp = 0;
	for (r = 0; r < ref.n; ++r) free(any[r]);
	free(any); free(T.a); free(T.slot); free(ref_taxon); free(ref_pfx);
	qf_seqs_free(&ref);
	return ret;
}
