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
 *     both orientations (minimap2-style DP with bounded look-back; a log-scaled
 *     penalty on |delta_asm - delta_ref| so multi-kb TE indels do not break a chain;
 *     syntenic chains bridge gaps up to --max-gap, 3 Mb by default, to match
 *     AnchorWave's multi-Mb deletions; off-chromosome chains use --max-gap-off).
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
	int64_t max_gap, max_gap_off, min_span, tight_gap, tight_diag;
} qc_opt_t;

typedef struct {
	const qc_opt_t *opt;
	const rb3_lift_t *lf;
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
	qc_chain_t *ch = 0;
	int64_t bp_strand[2] = {0, 0}, nc = 0;
	int is_syn = (strcmp(qc_suffix(q->sid->name[csid]), qc_suffix(q->sid->name[rsid])) == 0);
	int is_scaf = (strncmp(qc_suffix(q->sid->name[csid]), "chr", 3) != 0);
	// syntenic chains may bridge AnchorWave-scale deletions; chains to other chromosomes and
	// scaffolds use a stricter cap so sparse stray anchors are not bridged into long false spans
	int64_t max_gap = is_syn? o->max_gap : o->max_gap_off;

	a->n_pts += n;
	if (n < 2) return;
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
				sc = f[s][j] + (dr < k || dc < k? (dr < dc? dr : dc) : k) - qc_gap_pen(d);
				if (sc > best) best = sc, bj = j;
			}
			f[s][i] = best, p[s][i] = bj;
		}
	}
	// greedy extraction across both strands
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
			if (n_ch == m_ch) m_ch = m_ch? m_ch * 2 : 16, ch = (qc_chain_t*)realloc(ch, m_ch * sizeof(qc_chain_t));
			c = &ch[n_ch];
			c->idx = (int32_t*)malloc(cnt * sizeof(int32_t));
			c->n_idx = cnt, c->s = s, c->i_end = x;
			for (y = x; t > 0; y = p[s][y]) c->idx[--t] = y;
			c->rst = pts[c->idx[0]].r, c->ren = pts[x].r + k, c->n = cnt, c->bp = c->ren - c->rst;
			if (c->bp >= o->min_span) {
				bp_strand[s] += c->bp;
				++n_ch;
			} else free(c->idx);
		}
		for (y = x; y >= 0 && !used[y]; y = p[s][y]) used[y] = 1;
	}
	// classify and accumulate
	for (i = 0; i < n_ch; ++i) {
		qc_chain_t *c = &ch[i];
		int32_t cls, dom = bp_strand[1] > bp_strand[0]? 1 : 0;
		cls = !is_syn? (is_scaf? QC_SCAF : QC_OTHER) : (c->s == dom? QC_SYN_FWD : QC_SYN_INV);
		qc_add_ivl(a, c->rst, c->ren, cls);
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
		lp = rb3_lift_seq_points(q->lf, sx, &n);
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

int main_asmqc(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	static ko_longopt_t long_options[] = {
		{ "ref-prefix", ko_required_argument, 301 },
		{ "max-gap",    ko_required_argument, 302 },
		{ "min-span",   ko_required_argument, 303 },
		{ "tight-gap",  ko_required_argument, 304 },
		{ "tight-diag", ko_required_argument, 305 },
		{ "max-gap-off", ko_required_argument, 306 },
		{ 0, 0, 0 }
	};
	qc_opt_t opt;
	qc_ctx_t q;
	char *ref_prefix = 0, **taxa = 0;
	int32_t c, n_threads = 4, m_taxa = 0, t, r;
	int64_t sx, plen, i, max_ivl = 0;
	rb3_lift_t *lf;
	rb3_sid_t *sid;
	qc_ivl_t *buf;
	uint8_t **any;

	memset(&opt, 0, sizeof(opt));
	opt.k = 100, opt.stride = 2000, opt.h = 50, opt.min_anchor = 3;
	opt.max_gap = 3000000, opt.max_gap_off = 1000000, opt.min_span = 1000, opt.tight_gap = 1000, opt.tight_diag = 20;
	while ((c = ketopt(&o, argc, argv, 1, "k:s:t:h:n:", long_options)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 's') opt.stride = (int32_t)rb3_parse_num(o.arg);
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'h') opt.h = atoi(o.arg);
		else if (c == 'n') opt.min_anchor = atoi(o.arg);
		else if (c == 301) ref_prefix = o.arg;
		else if (c == 302) opt.max_gap = rb3_parse_num(o.arg);
		else if (c == 303) opt.min_span = rb3_parse_num(o.arg);
		else if (c == 304) opt.tight_gap = rb3_parse_num(o.arg);
		else if (c == 305) opt.tight_diag = rb3_parse_num(o.arg);
		else if (c == 306) opt.max_gap_off = rb3_parse_num(o.arg);
	}
	if (argc - o.ind < 2 || ref_prefix == 0) {
		fprintf(stderr, "Usage: ropebwt3 asmqc --ref-prefix=STR [options] <idx.lift> <idx.fmd.len.gz>\n");
		fprintf(stderr, "Options (k and stride are NOT stored in the .lift; pass the values used by `ropebwt3 lift`):\n");
		fprintf(stderr, "  -k INT          anchor k-mer length [%d]\n", opt.k);
		fprintf(stderr, "  -s NUM          anchor stride on the reference [%d]\n", opt.stride);
		fprintf(stderr, "  -t INT          threads (one assembly per job) [%d]\n", n_threads);
		fprintf(stderr, "  -h INT          chaining look-back in anchors [%d]\n", opt.h);
		fprintf(stderr, "  -n INT          min anchors per chain [%d]\n", opt.min_anchor);
		fprintf(stderr, "  --max-gap NUM   max reference/assembly gap inside a syntenic chain [%ld]\n", (long)opt.max_gap);
		fprintf(stderr, "  --max-gap-off NUM  same, for chains to other chromosomes / scaffolds [%ld]\n", (long)opt.max_gap_off);
		fprintf(stderr, "  --min-span NUM  min reference span of a kept chain [%ld]\n", (long)opt.min_span);
		fprintf(stderr, "  --tight-gap NUM / --tight-diag NUM  collinear-interval limits for the divergence estimate [%ld/%ld]\n",
				(long)opt.tight_gap, (long)opt.tight_diag);
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

	memset(&q, 0, sizeof(q));
	q.opt = &opt, q.lf = lf, q.sid = sid;
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

	// panel-anchorable reference stride bins: any assembly carries the anchor
	q.ref_nb = (int64_t*)calloc(q.n_ref, sizeof(int64_t));
	q.ref_cum = (int32_t**)calloc(q.n_ref, sizeof(int32_t*));
	any = (uint8_t**)calloc(q.n_ref, sizeof(uint8_t*));
	for (sx = 0; sx < sid->n_seq; ++sx)
		if (q.seq_ref[sx] >= 0) {
			r = q.seq_ref[sx];
			q.ref_nb[r] = sid->len[sx] / opt.stride + 1;
			any[r] = (uint8_t*)calloc(q.ref_nb[r], 1);
		}
	{
		int64_t n_off = 0, n_all = 0;
		for (sx = 0; sx < sid->n_seq; ++sx) {
			int64_t n;
			const rb3_liftpt_t *lp = rb3_lift_seq_points(lf, sx, &n);
			for (i = 0; i < n; ++i) {
				int32_t rr = lp[i].rsid >= 0 && lp[i].rsid < sid->n_seq? q.seq_ref[lp[i].rsid] : -1;
				if (rr < 0) continue;
				++n_all;
				if (lp[i].rpos % opt.stride != 0) ++n_off;
				if (lp[i].rpos / opt.stride < q.ref_nb[rr]) any[rr][lp[i].rpos / opt.stride] = 1;
			}
		}
		if (n_off > 0)
			fprintf(stderr, "[W::%s] %ld of %ld points are off the -s %d grid; is -s the stride the map was built with?\n",
					__func__, (long)n_off, (long)n_all, opt.stride);
		fprintf(stderr, "[M::%s] %ld liftover points, %d assemblies, %d reference sequences\n", __func__, (long)n_all, q.n_taxa, q.n_ref);
	}
	for (r = 0; r < q.n_ref; ++r) {
		q.ref_cum[r] = (int32_t*)malloc((q.ref_nb[r] + 1) * sizeof(int32_t));
		q.ref_cum[r][0] = 0;
		for (i = 0; i < q.ref_nb[r]; ++i) q.ref_cum[r][i + 1] = q.ref_cum[r][i] + any[r][i];
		free(any[r]);
	}
	free(any);

	q.acc = (qc_acc_t*)calloc((int64_t)q.n_taxa * q.n_ref, sizeof(qc_acc_t));
	kt_for(n_threads, qc_worker, &q, q.n_taxa);

	// report
	for (i = 0; i < (int64_t)q.n_taxa * q.n_ref; ++i)
		if (q.acc[i].n_ivl > max_ivl) max_ivl = q.acc[i].n_ivl;
	buf = (qc_ivl_t*)malloc((max_ivl + 1) * sizeof(qc_ivl_t));
	printf("taxon\tchrom\tref_len\tn_anchorable\tn_pts\tmapped_bp\tfrac_mapped\tinv_bp\tother_chr_bp\tscaf_bp\tgap_bp"
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
			qc_print_row(stdout, taxa[t], sid->name[sx], sid->len[sx], q.ref_cum[r][q.ref_nb[r]], a, m, inv, oth, sc, opt.k);
			tl += sid->len[sx], tanch += q.ref_cum[r][q.ref_nb[r]], tm += m, tinv += inv, toth += oth, tsc += sc;
			tot.n_pts += a->n_pts, tot.n_chain += a->n_chain, tot.n_chain_inv += a->n_chain_inv;
			tot.ins_n += a->ins_n, tot.ins_bp += a->ins_bp, tot.del_n += a->del_n, tot.del_bp += a->del_bp;
			if (a->ins_max > tot.ins_max) tot.ins_max = a->ins_max;
			if (a->del_max > tot.del_max) tot.del_max = a->del_max;
			tot.tight_bp += a->tight_bp, tot.tight_present += a->tight_present, tot.tight_missing += a->tight_missing;
		}
		qc_print_row(stdout, taxa[t], "ALL", tl, tanch, &tot, tm, tinv, toth, tsc, opt.k);
	}
	free(buf);

	for (i = 0; i < (int64_t)q.n_taxa * q.n_ref; ++i) free(q.acc[i].ivl);
	for (r = 0; r < q.n_ref; ++r) free(q.ref_cum[r]);
	for (t = 0; t < q.n_taxa; ++t) free(taxa[t]);
	free(q.acc); free(q.ref_cum); free(q.ref_nb); free(taxa);
	free(q.seq_taxon); free(q.seq_ref);
	rb3_sid_destroy(sid);
	rb3_lift_destroy(lf);
	return 0;
}
