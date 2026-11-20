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
#include <math.h>
#include "rb3priv.h"
#include "io.h"
#include "lift.h"
#include "ketopt.h"
#include "kthread.h"
#include <pthread.h>

static FILE *qc_chain_fp = 0; // --chains output (NULL = off)
static pthread_mutex_t qc_chain_lock = PTHREAD_MUTEX_INITIALIZER;
static const char *qc_cls_name[] = { "syn", "syn_inv", "other_chr", "scaffold" };

#define QC_SYN_FWD 0
#define QC_SYN_INV 1
#define QC_OTHER   2
#define QC_SCAF    3

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
} qc_ctx_t;

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
static void qc_chain_group(const qc_ctx_t *q, int32_t taxon, int32_t csid, int32_t rsid, int64_t n, const qc_pt_t *pts)
{
	const qc_opt_t *o = q->opt;
	int32_t r = q->seq_ref[rsid], s, i, j, k = o->k, n_ch = 0, m_ch = 0;
	qc_acc_t *a = &q->acc[(int64_t)taxon * q->n_ref + r];
	double *f[2];
	int32_t *p[2];
	uint8_t *used;
	qc_cand_t *cand;
	qc_chain_t *ch = 0, *blk = 0;
	int32_t n_blk = 0, m_blk = 0;
	int64_t bp_strand[2] = {0, 0}, nc = 0;
	int is_syn = (strcmp(qc_suffix(q->sid->name[csid]), qc_suffix(q->sid->name[rsid])) == 0);
	int is_scaf = (strncmp(qc_suffix(q->sid->name[csid]), "chr", 3) != 0);
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
		free(c->idx);
	}
	free(ch); free(used); free(cand);
	for (s = 0; s < 2; ++s) free(f[s]), free(p[s]);
}

static void qc_worker(void *data, long t, int tid) // one taxon per job
{
	qc_ctx_t *q = (qc_ctx_t*)data;
	int64_t sx, m_tmp = 0, m_pts = 0;
	rb3_liftpt_t *tmp = 0;
	qc_pt_t *pts = 0;
	(void)tid;
	for (sx = 0; sx < q->sid->n_seq; ++sx) {
		int64_t n, i, i0;
		const rb3_liftpt_t *lp;
		if (q->seq_taxon[sx] != t) continue;
		lp = q->sp[sx], n = q->sn[sx];
		if (n == 0) continue;
		if (n > m_tmp) m_tmp = n, tmp = (rb3_liftpt_t*)realloc(tmp, n * sizeof(rb3_liftpt_t));
		if (n > m_pts) m_pts = n, pts = (qc_pt_t*)realloc(pts, n * sizeof(qc_pt_t));
		memcpy(tmp, lp, n * sizeof(rb3_liftpt_t));
		qsort(tmp, n, sizeof(rb3_liftpt_t), qc_cmp_lpt);
		for (i0 = 0, i = 1; i <= n; ++i) {
			if (i == n || tmp[i].rsid != tmp[i0].rsid) {
				int64_t j;
				if (q->seq_ref[tmp[i0].rsid] >= 0) {
					for (j = i0; j < i; ++j) pts[j - i0].r = tmp[j].rpos, pts[j - i0].c = tmp[j].cpos;
					qc_chain_group(q, (int32_t)t, (int32_t)sx, tmp[i0].rsid, i - i0, pts);
				}
				i0 = i;
			}
		}
	}
	free(tmp); free(pts);
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
	fprintf(fp, "%s\t%s\t%ld\t%ld\t%ld\t%ld\t%.6f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.6f\t%ld\t%ld\t%.6f\t%.6f\n",
			taxon, chrom, (long)len, (long)n_anch, (long)a->n_pts, (long)mapped, len? (double)mapped / len : 0.0,
			(long)inv, (long)other, (long)scaf, (long)(len - mapped), (long)a->n_chain, (long)a->n_chain_inv,
			(long)a->ins_n, (long)a->ins_bp, (long)a->ins_max, (long)a->del_n, (long)a->del_bp, (long)a->del_max,
			(long)a->tight_bp, len? (double)a->tight_bp / len : 0.0, (long)a->tight_present, (long)a->tight_missing, fi, div);
}


/* Shared driver: sid = all sequences (reference ones start with ref_prefix), sp/sn = the
 * liftover points of each assembly sequence. any_in (optional, indexed by reference order)
 * marks the anchorable stride bins used as the divergence denominator; NULL derives it from
 * the points (panel mode: a bin is anchorable if any assembly carries its k-mer). */
static int qc_run(const qc_opt_t *opt, const rb3_sid_t *sid, const char *ref_prefix,
				  const rb3_liftpt_t **sp, const int64_t *sn, uint8_t **any_in, int n_threads, FILE *fp, int print_header)
{
	qc_ctx_t q;
	char **taxa = 0;
	int32_t m_taxa = 0, t, r;
	int64_t sx, plen, i, max_ivl = 0;
	qc_ivl_t *buf;
	uint8_t **any;

	memset(&q, 0, sizeof(q));
	q.opt = opt, q.sid = sid, q.sp = sp, q.sn = sn;
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
	kt_for(n_threads, qc_worker, &q, q.n_taxa);

	// report
	for (i = 0; i < (int64_t)q.n_taxa * q.n_ref; ++i)
		if (q.acc[i].n_ivl > max_ivl) max_ivl = q.acc[i].n_ivl;
	buf = (qc_ivl_t*)malloc((max_ivl + 1) * sizeof(qc_ivl_t));
	if (print_header) fprintf(fp, "taxon\tchrom\tref_len\tn_anchorable\tn_pts\tmapped_bp\tfrac_mapped\tinv_bp\tother_chr_bp\tscaf_bp\tgap_bp"
		   "\tn_chain\tn_chain_inv\tins_n\tins_bp_lb\tins_max\tdel_n\tdel_bp_lb\tdel_max"
		   "\ttight_bp\tfrac_tight\ttight_present\ttight_missing\tfrac_intact\tdiv_est\n");
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
			tot.n_pts += a->n_pts, tot.n_chain += a->n_chain, tot.n_chain_inv += a->n_chain_inv;
			tot.ins_n += a->ins_n, tot.ins_bp += a->ins_bp, tot.del_n += a->del_n, tot.del_bp += a->del_bp;
			if (a->ins_max > tot.ins_max) tot.ins_max = a->ins_max;
			if (a->del_max > tot.del_max) tot.del_max = a->del_max;
			tot.tight_bp += a->tight_bp, tot.tight_present += a->tight_present, tot.tight_missing += a->tight_missing;
		}
		qc_print_row(fp, taxa[t], "ALL", tl, tanch, &tot, tm, tinv, toth, tsc, opt->k);
	}
	free(buf);

	for (i = 0; i < (int64_t)q.n_taxa * q.n_ref; ++i) free(q.acc[i].ivl);
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
	ret = qc_run(&opt, sid, ref_prefix, sp, sn, 0, n_threads, stdout, 1);
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
		ret |= qc_run(&opt, &sid, ref_pfx, sp, sn, any, n_threads, stdout, !header_done);
		fflush(stdout);
		header_done = 1;
		free(pt); free(sp); free(sn); free(sid.name); free(sid.len);
		qf_seqs_free(&as);
		free(taxon);
	}
	fprintf(stderr, "[M::%s::%.1f] done\n", __func__, rb3_realtime() - t0);
	if (qc_chain_fp) fclose(qc_chain_fp), qc_chain_fp = 0;
	for (r = 0; r < ref.n; ++r) free(any[r]);
	free(any); free(T.a); free(T.slot); free(ref_taxon); free(ref_pfx);
	qf_seqs_free(&ref);
	return ret;
}
