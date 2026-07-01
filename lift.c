/* lift.c -- the "second SSA": a carrier->reference coordinate liftover for refmap.
 *
 * Motivation (see docs/results-maize.md, experiment E4): placing a query that is
 * absent from the reference by WALKING outward through carriers is slow and its
 * single-anchor fallback is inaccurate. Instead, precompute a coordinate map from
 * every carrier genome to the reference, built from unique shared anchors, and at
 * query time PROJECT a carrier hit to an approximate reference coordinate. Where
 * no confident projection exists we return "no position" (NULL) rather than guess,
 * which keeps precision.
 *
 * Build: slide a k-mer along the reference sequences; for every k-mer that occurs
 * few enough times (<= max_occ, i.e. single-copy-per-taxon) locate all its hits.
 * A hit in the reference paired with each hit in a carrier is a liftover point
 * (carrier_seq, carrier_pos) -> (ref_seq, ref_pos), stored in forward-strand
 * coordinates and sorted per carrier sequence.
 *
 * Project: given a carrier hit, gather liftover points within +-win bp, take the
 * majority reference sequence, pick orientation (+1/-1) by the tighter residual,
 * and return the median-intercept projection -- unless support is thin or the
 * residual MAD is large, in which case the locus is not collinear -> NULL.
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "fm-index.h"
#include "io.h"
#include "rb3priv.h"
#include "ketopt.h"
#include "kalloc.h"
#include "kthread.h"
#include "lift.h"

typedef struct {          // one liftover point (flat build array + per-seq arrays)
	int64_t cpos;         // carrier forward-strand start
	int64_t rpos;         // reference forward-strand start
	int32_t csid;         // carrier sequence index (in [0,n_seq))
	int32_t rsid;         // reference sequence index
} rb3_liftpt_t;

struct rb3_lift_s {
	int64_t n_seq;        // number of sequences (carrier keys)
	int64_t *off;         // off[sid]..off[sid+1] is the point range for carrier sid
	rb3_liftpt_t *pt;     // all points, sorted by (csid, cpos)
	int64_t n_pt;
};

/* ---- build ------------------------------------------------------------- */

static int lift_cmp(const void *a, const void *b)
{
	const rb3_liftpt_t *x = (const rb3_liftpt_t*)a, *y = (const rb3_liftpt_t*)b;
	if (x->csid != y->csid) return x->csid < y->csid? -1 : 1;
	if (x->cpos != y->cpos) return x->cpos < y->cpos? -1 : 1;
	return 0;
}

// forward-strand start of a length-rlen hit at (sid,pos)
static inline int64_t lift_fwd(const rb3_sid_t *sid, int64_t hsid, int64_t pos, int32_t rlen)
{
	int64_t clen = sid->len[hsid>>1];
	return (hsid & 1)? clen - (pos + rlen) : pos;
}

// Growable, multithreaded build context. Each of n_threads keeps its own output
// buffer + kalloc arena; anchors of one sequence are farmed out with kt_for and
// merged at finalize.
typedef struct {
	const rb3_fmi_t *f;
	const uint8_t *is_ref;
	int32_t k, stride, n_threads;
	int64_t max_occ;
	rb3_liftpt_t **buf;   // [tid] growable point buffer
	int64_t *bn, *bm;     // [tid] count/cap
	void **km;            // [tid] kalloc arena
	int64_t n_kmer;
	// per-sequence job passed to the workers:
	const uint8_t *seq;
	int64_t seq_len;
} lift_builder_t;

static lift_builder_t *lift_builder_init(const rb3_fmi_t *f, const uint8_t *is_ref, int32_t k, int32_t stride, int64_t max_occ, int32_t n_threads)
{
	int t;
	lift_builder_t *b = RB3_CALLOC(lift_builder_t, 1);
	b->f = f, b->is_ref = is_ref, b->k = k, b->stride = stride, b->max_occ = max_occ;
	b->n_threads = n_threads < 1? 1 : n_threads;
	b->buf = RB3_CALLOC(rb3_liftpt_t*, b->n_threads);
	b->bn = RB3_CALLOC(int64_t, b->n_threads);
	b->bm = RB3_CALLOC(int64_t, b->n_threads);
	b->km = RB3_CALLOC(void*, b->n_threads);
	for (t = 0; t < b->n_threads; ++t) b->km[t] = km_init();
	return b;
}

// Process one anchor (index i = anchor number within the current sequence).
static void lift_worker(void *data, long i, int tid)
{
	lift_builder_t *b = (lift_builder_t*)data;
	const rb3_fmi_t *f = b->f;
	void *km = b->km[tid];
	int32_t k = b->k, j, bad = 0;
	const uint8_t *q = b->seq + (int64_t)i * b->stride;
	rb3_sai_t ik, ok[RB3_ASIZE];
	rb3_pos_t pos[64];
	int64_t np, r, nref = 0, ref_sid = -1, ref_pos = -1, cap = b->max_occ < 64? b->max_occ : 64;
	if (q[k-1] < 1 || q[k-1] > 4) return;
	rb3_fmd_set_intv(f, q[k-1], &ik);
	if (ik.size == 0) return;
	for (j = k - 2; j >= 0; --j) {
		if (q[j] < 1 || q[j] > 4) { bad = 1; break; }
		rb3_fmd_extend(f, &ik, ok, 1);
		if (ok[q[j]].size == 0) { bad = 1; break; }
		ik = ok[q[j]];
	}
	if (bad || ik.size < 1 || ik.size > b->max_occ) return;
	np = rb3_ssa_multi(km, f, f->ssa, ik.x[0], ik.x[0] + ik.size, cap, pos);
	for (r = 0; r < np; ++r) {           // require exactly one reference occurrence
		int64_t sid = pos[r].sid, sx = sid >> 1;
		if (b->is_ref[sx]) ref_sid = sx, ref_pos = lift_fwd(f->sid, sid, pos[r].pos, k), ++nref;
	}
	if (nref != 1) return;
	for (r = 0; r < np; ++r) {           // pair each carrier hit with the reference hit
		int64_t sid = pos[r].sid, sx = sid >> 1;
		if (b->is_ref[sx]) continue;
		Kgrow(0, rb3_liftpt_t, b->buf[tid], b->bn[tid], b->bm[tid]);
		b->buf[tid][b->bn[tid]].csid = (int32_t)sx;
		b->buf[tid][b->bn[tid]].cpos = lift_fwd(f->sid, sid, pos[r].pos, k);
		b->buf[tid][b->bn[tid]].rsid = (int32_t)ref_sid;
		b->buf[tid][b->bn[tid]].rpos = ref_pos;
		b->bn[tid]++;
	}
}

// Add one sequence (already nt6-encoded): fan its anchors across the threads.
static void lift_builder_add(lift_builder_t *b, int64_t len, const uint8_t *s)
{
	int64_t n_anchor;
	if (len < b->k) return;
	n_anchor = (len - b->k) / b->stride + 1;
	b->n_kmer += n_anchor;
	b->seq = s, b->seq_len = len;
	kt_for(b->n_threads, lift_worker, b, n_anchor);
}

static rb3_lift_t *lift_builder_finalize(lift_builder_t *b)
{
	rb3_lift_t *lf = RB3_CALLOC(rb3_lift_t, 1);
	int64_t i, n_pt = 0, o = 0;
	int t;
	lf->n_seq = b->f->sid->n_seq;
	for (t = 0; t < b->n_threads; ++t) n_pt += b->bn[t];
	lf->n_pt = n_pt;
	lf->pt = RB3_MALLOC(rb3_liftpt_t, n_pt > 0? n_pt : 1);
	for (t = 0; t < b->n_threads; ++t) {          // merge per-thread buffers
		if (b->bn[t]) memcpy(lf->pt + o, b->buf[t], b->bn[t] * sizeof(rb3_liftpt_t));
		o += b->bn[t];
		free(b->buf[t]); km_destroy(b->km[t]);
	}
	qsort(lf->pt, lf->n_pt, sizeof(rb3_liftpt_t), lift_cmp);
	lf->off = RB3_CALLOC(int64_t, lf->n_seq + 1);
	for (i = 0; i < lf->n_pt; ++i) lf->off[lf->pt[i].csid + 1]++;
	for (i = 0; i < lf->n_seq; ++i) lf->off[i + 1] += lf->off[i];
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] scanned %ld k-mers, %ld liftover points\n",
				__func__, (long)b->n_kmer, (long)lf->n_pt);
	free(b->buf); free(b->bn); free(b->bm); free(b->km); free(b);
	return lf;
}

void rb3_lift_destroy(rb3_lift_t *lf)
{
	if (lf == 0) return;
	free(lf->pt); free(lf->off); free(lf);
}

/* ---- serialize --------------------------------------------------------- */

int rb3_lift_dump(const rb3_lift_t *lf, const char *fn)
{
	FILE *fp = fn && strcmp(fn, "-")? fopen(fn, "wb") : stdout;
	if (fp == 0) return -1;
	fwrite("LIFT\1", 1, 5, fp);
	fwrite(&lf->n_seq, 8, 1, fp);
	fwrite(&lf->n_pt, 8, 1, fp);
	fwrite(lf->off, 8, lf->n_seq + 1, fp);
	fwrite(lf->pt, sizeof(rb3_liftpt_t), lf->n_pt, fp);
	if (fp != stdout) fclose(fp);
	return 0;
}

rb3_lift_t *rb3_lift_restore(const char *fn)
{
	FILE *fp = fn && strcmp(fn, "-")? fopen(fn, "rb") : stdin;
	char magic[5];
	rb3_lift_t *lf;
	if (fp == 0) return 0;
	if (fread(magic, 1, 5, fp) != 5 || strncmp(magic, "LIFT\1", 5) != 0) { if (fp != stdin) fclose(fp); return 0; }
	lf = RB3_CALLOC(rb3_lift_t, 1);
	if (fread(&lf->n_seq, 8, 1, fp) != 1 || fread(&lf->n_pt, 8, 1, fp) != 1) { free(lf); if (fp != stdin) fclose(fp); return 0; }
	lf->off = RB3_MALLOC(int64_t, lf->n_seq + 1);
	lf->pt = RB3_MALLOC(rb3_liftpt_t, lf->n_pt);
	if (fread(lf->off, 8, lf->n_seq + 1, fp) != (size_t)(lf->n_seq + 1) ||
		fread(lf->pt, sizeof(rb3_liftpt_t), lf->n_pt, fp) != (size_t)lf->n_pt) {
		rb3_lift_destroy(lf); if (fp != stdin) fclose(fp); return 0;
	}
	if (fp != stdin) fclose(fp);
	return lf;
}

int64_t rb3_lift_n_seq(const rb3_lift_t *lf) { return lf ? lf->n_seq : 0; }

/* ---- project ----------------------------------------------------------- */

static int64_t i64_median(int64_t *a, int64_t n)
{
	// nth_element would be faster; n is small (window). Simple sort of a copy.
	int64_t i, j; // insertion sort
	for (i = 1; i < n; ++i) {
		int64_t v = a[i];
		for (j = i - 1; j >= 0 && a[j] > v; --j) a[j+1] = a[j];
		a[j+1] = v;
	}
	return a[n>>1];
}

// Project a carrier hit (csid, cpos) to a reference coordinate. Returns 1 and sets
// *out_rsid,*out_rpos on success; returns 0 (NULL slot) when not confidently
// collinear. win/max_mad in bp, min_support minimum anchors in window.
int rb3_lift_project(const rb3_lift_t *lf, void *km, int32_t csid, int64_t cpos,
					 int64_t win, int64_t max_mad, int32_t min_support,
					 int64_t *out_rsid, int64_t *out_rpos)
{
	int64_t lo = lf->off[csid], hi = lf->off[csid + 1];
	int64_t a, b, i, n, cnt, best_rsid, best_n, sign;
	int64_t *resP, *resM, madP, madM, med;
	rb3_liftpt_t *pt = lf->pt;
	if (hi - lo < min_support) return 0;
	// binary search the window [cpos-win, cpos+win] within [lo,hi)
	{
		int64_t l = lo, r = hi; // first index with cpos >= target-win
		int64_t t = cpos - win;
		while (l < r) { int64_t mid = (l + r) >> 1; if (pt[mid].cpos < t) l = mid + 1; else r = mid; }
		a = l;
		l = a, r = hi; t = cpos + win;
		while (l < r) { int64_t mid = (l + r) >> 1; if (pt[mid].cpos <= t) l = mid + 1; else r = mid; }
		b = l;
	}
	n = b - a;
	if (n < min_support) return 0;
	// majority reference sequence in the window
	best_rsid = -1, best_n = 0;
	for (i = a; i < b; ++i) {
		int64_t rs = pt[i].rsid; cnt = 0;
		int64_t j;
		for (j = a; j < b; ++j) if (pt[j].rsid == rs) ++cnt;
		if (cnt > best_n) best_n = cnt, best_rsid = rs;
	}
	if (best_n < min_support) return 0;
	// residuals for slope +1 (rpos-cpos) and -1 (rpos+cpos) over majority-chr points
	resP = Kmalloc(km, int64_t, best_n);
	resM = Kmalloc(km, int64_t, best_n);
	n = 0;
	for (i = a; i < b; ++i)
		if (pt[i].rsid == best_rsid) {
			resP[n] = pt[i].rpos - pt[i].cpos;
			resM[n] = pt[i].rpos + pt[i].cpos;
			++n;
		}
	{
		int64_t *cpP = Kmalloc(km, int64_t, n), *cpM = Kmalloc(km, int64_t, n);
		int64_t mP, mM;
		memcpy(cpP, resP, n * 8); memcpy(cpM, resM, n * 8);
		mP = i64_median(cpP, n); mM = i64_median(cpM, n);
		for (i = 0; i < n; ++i) { cpP[i] = llabs(resP[i] - mP); cpM[i] = llabs(resM[i] - mM); }
		madP = i64_median(cpP, n); madM = i64_median(cpM, n);
		kfree(km, cpP); kfree(km, cpM);
		if (madP <= madM) sign = 1, med = mP;
		else sign = -1, med = mM;
		if ((madP <= madM? madP : madM) > max_mad) { kfree(km, resP); kfree(km, resM); return 0; }
	}
	*out_rsid = best_rsid;
	*out_rpos = sign * cpos + med;
	kfree(km, resP); kfree(km, resM);
	return 1;
}

/* ---- the `lift` subcommand (build + dump) ------------------------------ */

int main_lift(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, k = 100, stride = 2000, n_threads = 4;
	int64_t max_occ = -1; // -1 = auto (#taxa)
	char *ref_prefix = 0, *out = "-";
	rb3_fmi_t f;
	rb3_lift_t *lf;
	lift_builder_t *bld;
	int64_t j, plen, n_ref = 0, n_taxa = 0;
	uint8_t *is_ref;
	rb3_seqio_t *fp;

	static ko_longopt_t lo[] = { { "ref-prefix", ko_required_argument, 301 }, { 0, 0, 0 } };
	while ((c = ketopt(&o, argc, argv, 1, "k:s:m:o:t:", lo)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 's') stride = rb3_parse_num(o.arg);
		else if (c == 'm') max_occ = atol(o.arg);
		else if (c == 'o') out = o.arg;
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 301) ref_prefix = o.arg;
	}
	if (argc - o.ind < 2 || ref_prefix == 0) {
		fprintf(stderr, "Usage: ropebwt3 lift --ref-prefix=STR [-k %d -s %d -t %d -m auto] <idx.fmd> <ref.fa> [...]\n", k, stride, n_threads);
		fprintf(stderr, "Builds the carrier->reference coordinate liftover (the \"second SSA\") used by\n");
		fprintf(stderr, "`refmap --lift`. Feed the reference genome FASTA(s); anchors are shared k-mers.\n");
		return 1;
	}
	if (rb3_fmi_load_all(&f, argv[o.ind], RB3_LOAD_ALL) < 0) return 1;
	if (f.ssa == 0 || f.sid == 0) { fprintf(stderr, "ERROR: lift needs .ssa and .len.gz\n"); return 1; }

	plen = strlen(ref_prefix);
	is_ref = RB3_CALLOC(uint8_t, f.sid->n_seq);
	for (j = 0; j < f.sid->n_seq; ++j)
		if (strncmp(f.sid->name[j], ref_prefix, plen) == 0) is_ref[j] = 1, ++n_ref;
	if (n_ref == 0) { fprintf(stderr, "ERROR: no sequence starts with '%s'\n", ref_prefix); return 1; }
	if (max_occ < 0) { // auto = number of distinct taxa (name prefixes before '_')
		for (j = 0; j < f.sid->n_seq; ++j) {
			const char *nm = f.sid->name[j], *us = strchr(nm, '_');
			int32_t pl = us? (int32_t)(us - nm) : (int32_t)strlen(nm), dup = 0, t;
			for (t = 0; t < j; ++t) {
				const char *n2 = f.sid->name[t], *u2 = strchr(n2, '_');
				int32_t l2 = u2? (int32_t)(u2 - n2) : (int32_t)strlen(n2);
				if (l2 == pl && strncmp(n2, nm, pl) == 0) { dup = 1; break; }
			}
			if (!dup) ++n_taxa;
		}
		max_occ = n_taxa;
	}
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] ref sequences: %ld; max_occ: %ld\n", __func__, (long)n_ref, (long)max_occ);

	bld = lift_builder_init(&f, is_ref, k, stride, max_occ, n_threads);
	for (j = o.ind + 1; j < argc; ++j) {
		const char *name; int64_t len; char *s;
		fp = rb3_seq_open(argv[j], 0);
		if (fp == 0) { fprintf(stderr, "ERROR: cannot open '%s'\n", argv[j]); return 1; }
		while ((s = rb3_seq_read1(fp, &len, &name)) != 0) {
			rb3_char2nt6(len, (uint8_t*)s);
			lift_builder_add(bld, len, (uint8_t*)s);
		}
		rb3_seq_close(fp);
	}
	lf = lift_builder_finalize(bld);
	rb3_lift_dump(lf, out);
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] wrote liftover: %ld points over %ld sequences to %s\n",
				__func__, (long)lf->n_pt, (long)lf->n_seq, out);
	rb3_lift_destroy(lf);
	free(is_ref);
	rb3_fmi_free(&f);
	return 0;
}
