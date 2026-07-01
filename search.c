#include <math.h>
#include "fm-index.h"
#include "align.h"
#include "rb3priv.h"
#include "io.h"
#include "ketopt.h"
#include "kthread.h"
#include "kalloc.h"
#include "lift.h"

typedef enum { RB3_SA_MEM_TG, RB3_SA_MEM_ORI, RB3_SA_SW, RB3_SA_HAPDIV, RB3_SA_REFMAP } rb3_search_algo_t;

#define RB3_WALK_CONSENSUS  0 // follow the base shared by the most carriers
#define RB3_WALK_STRICT     1 // stop walking at the first carrier disagreement
#define RB3_WALK_PERCARRIER 2 // one outward walk (and one result) per carrier genome

#define RB3_MF_NO_KALLOC   0x1
#define RB3_MF_WRITE_UNMAP 0x2
#define RB3_MF_WRITE_COV   0x4
#define RB3_MF_WRITE_ALL   0x8
#define RB3_MF_BOTH_DIR    0x10

typedef struct {
	uint32_t flag;
	int32_t n_threads, min_gap_len, hapdiv_k, hapdiv_w;
	int32_t max_pos;
	rb3_search_algo_t algo;
	int64_t min_occ, min_len, max_all_out;
	int64_t batch_size;
	char *ref_prefix;  // refmap: reference sequences are those whose name starts with this
	int32_t max_walk;  // refmap: max bases to walk outward along carriers, per flank
	int8_t walk_mode;  // refmap: RB3_WALK_*
	int64_t max_occ;   // refmap: reject reads/anchors occurring > max_occ times (0 = off; <0 = auto = #taxa)
	int64_t max_bracket; // refmap: reject a PLACED if |cR-cL| > max_bracket (0 = off)
	int8_t two_flank;  // refmap: require both flanks to anchor concordantly (1 = on)
	char *lift_fn;     // refmap: liftover file -> project carrier hits instead of walking
	int64_t lift_win, lift_mad; // refmap: liftover projection window / max residual MAD (bp)
	int32_t kmer_len, kmer_step, min_agree; // refmap: k-mer-agreement placement (0 = off)
	int64_t kmer_cluster;                   // refmap: cluster tolerance for agreeing k-mers (bp)
	rb3_swopt_t swo;
} rb3_mopt_t;

void rb3_mopt_init(rb3_mopt_t *opt)
{
	memset(opt, 0, sizeof(rb3_mopt_t));
	opt->n_threads = 4;
	opt->min_occ = 1;
	opt->min_len = 19;
	opt->hapdiv_k = 101;
	opt->hapdiv_w = 50;
	opt->batch_size = 100000000;
	opt->algo = RB3_SA_MEM_TG;
	opt->max_walk = 5000;
	opt->walk_mode = RB3_WALK_CONSENSUS;
	opt->max_occ = 0;     // off by default (E0 baseline)
	opt->max_bracket = 0; // off by default
	opt->two_flank = 0;   // off by default
	opt->lift_fn = 0;     // off by default (walk)
	opt->lift_win = 500000, opt->lift_mad = 200000;
	opt->kmer_len = 0;    // off by default (whole-read placement)
	opt->kmer_step = 15, opt->min_agree = 2, opt->kmer_cluster = 2000;
	rb3_swopt_init(&opt->swo);
}

typedef struct mp_tbuf_s {
	void *km;
	int32_t n_gap, m_gap;
	uint64_t *gap;
	rb3_sai_v mem; // this is allocated from km
} m_tbuf_t;

typedef struct {
	int64_t n_pos;
	rb3_sai_t mem;
	rb3_pos_t *pos;
} m_sai_pos_t;

typedef struct {
	char *name;
	uint8_t *seq;
	int64_t id, n_pos;
	int32_t len, n_mem, n_gap;
	uint64_t *gap;
	m_sai_pos_t *mem;
} m_seq_t;

typedef struct {
	const rb3_mopt_t *opt;
	int64_t id;
	rb3_fmi_t fmi;
	rb3_seqio_t *fp;
	uint8_t *is_ref; // refmap: is_ref[k]!=0 iff sequence k (in [0,n_seq)) belongs to the reference
	int64_t n_ref;   // refmap: number of reference sequences
	rb3_lift_t *lift; // refmap: carrier->reference liftover (NULL = walk)
} pipeline_t;

typedef struct {
	int32_t id, offset;
	rb3_hapdiv_t r;
} m_hapdiv_t;

// refmap: result of placing one query on the reference genome
#define RB3_RM_UNPLACED 0
#define RB3_RM_PLACED   1
#define RB3_RM_ONE_SIDE 2
#define RB3_RM_EXACT    3
#define RB3_RM_MULTI    4 // query occurs > max_occ times (repeat/retro): not placed

typedef struct refmap_rst_s {
	int8_t status;       // RB3_RM_*
	int8_t strand;       // reference strand (0 forward, 1 reverse); valid when placed/exact
	int32_t qlen;
	int32_t n_carrier;   // distinct carrier sequences seen (capped)
	int64_t ref_sid;     // reference sequence index (in [0,n_seq)); -1 if none
	int64_t cL, cR;      // reference coordinates bracketing the query; -1 if unknown
	int64_t ins_size;    // implied size inserted into the carrier relative to the reference
	int32_t n_car_list;  // number of carriers stored below
	rb3_pos_t *carriers; // a few carrier (sid,pos); allocated with RB3_MALLOC
	int32_t n_sub;       // per-carrier mode: number of sub-results (one per carrier)
	struct refmap_rst_s *sub; // per-carrier mode: one placement per carrier
	int32_t n_vote, agree, second, mapq; // --kmer: informative tiles, agreeing k-mers, runner-up, calibrated MAPQ
} refmap_rst_t;

typedef struct {
	const pipeline_t *p;
	int32_t n_seq, n_hapdiv;
	m_seq_t *seq;
	rb3_swrst_t *rst, *rst_rev;
	m_hapdiv_t *hapdiv;
	refmap_rst_t *refmap; // refmap: one entry per query
	m_tbuf_t *buf;
} step_t;

static void refmap_query(void *km, const pipeline_t *p, const m_seq_t *s, refmap_rst_t *r);
static int refmap_anchor_flank(void *km, const pipeline_t *p, const uint8_t *flank, int32_t flen, int side,
							   int64_t *out_sid, int *out_strand, int64_t *out_coord, int32_t *out_mlen);

static void worker_for_seq(void *data, long i, int tid)
{
	step_t *t = (step_t*)data;
	const pipeline_t *p = t->p;
	m_seq_t *s = &t->seq[i];
	m_tbuf_t *b = &t->buf[tid];
	if (rb3_dbg_flag & RB3_DBG_QNAME)
		fprintf(stderr, "Q\t%s\t%d\n", s->name, tid);
	rb3_char2nt6(s->len, s->seq);
	if (p->opt->algo == RB3_SA_SW) { // BWA-SW
		rb3_sw(b->km, &p->opt->swo, &p->fmi, s->len, s->seq, &t->rst[i]);
		if (t->rst_rev) {
			rb3_revcomp6(s->len, s->seq);
			rb3_sw(b->km, &p->opt->swo, &p->fmi, s->len, s->seq, &t->rst_rev[i]);
			rb3_revcomp6(s->len, s->seq);
		}
	} else if (p->opt->algo == RB3_SA_REFMAP) { // place the query on the reference genome
		refmap_query(b->km, p, s, &t->refmap[i]);
	} else { // MEM algorithms
		int32_t i;
		b->mem.n = 0;
		if (p->opt->algo == RB3_SA_MEM_TG)
			rb3_fmd_smem_TG(b->km, &p->fmi, s->len, s->seq, &b->mem, p->opt->min_occ, p->opt->min_len);
		else if (p->opt->algo == RB3_SA_MEM_ORI)
			rb3_fmd_smem(b->km, &p->fmi, s->len, s->seq, &b->mem, p->opt->min_occ, p->opt->min_len);
		s->n_mem = b->mem.n;
		s->mem = RB3_CALLOC(m_sai_pos_t, s->n_mem);
		for (i = 0; i < s->n_mem; ++i)
			s->mem[i].mem = b->mem.a[i];
		if (p->opt->min_gap_len > 0) { // find gaps not covered by MEMs
			int32_t last = 0;
			b->n_gap = 0;
			Kgrow(b->km, uint64_t, b->gap, b->mem.n + 1, b->m_gap);
			for (i = 0; i < b->mem.n; ++i) {
				int32_t st = b->mem.a[i].info>>32, en = (int32_t)b->mem.a[i].info;
				if (st > last) {
					if (st - last >= p->opt->min_gap_len)
						b->gap[b->n_gap++] = (uint64_t)last<<32 | st;
					last = en;
				} else last = last > en? last : en;
			}
			if (s->len - last >= p->opt->min_gap_len)
				b->gap[b->n_gap++] = (uint64_t)last<<32 | s->len;
			s->n_gap = b->n_gap;
			s->gap = RB3_MALLOC(uint64_t, s->n_gap);
			memcpy(s->gap, b->gap, s->n_gap * 8);
		} else if (p->opt->max_pos > 0) {
			#if 1 // faster algorithm
			rb3_pos_t *pos;
			pos = Kmalloc(b->km, rb3_pos_t, p->opt->max_pos);
			for (i = 0; i < s->n_mem; ++i) {
				m_sai_pos_t *q = &s->mem[i];
				q->n_pos = rb3_ssa_multi(b->km, &p->fmi, p->fmi.ssa, q->mem.x[0], q->mem.x[0] + q->mem.size, p->opt->max_pos, pos);
				q->pos = RB3_MALLOC(rb3_pos_t, q->n_pos);
				memcpy(q->pos, pos, sizeof(rb3_pos_t) * q->n_pos);
			}
			kfree(b->km, pos);
			#else // naive algorithm
			for (i = 0; i < s->n_mem; ++i) {
				m_sai_pos_t *q = &s->mem[i];
				int32_t j;
				q->n_pos = q->mem.size < p->opt->max_pos? q->mem.size : p->opt->max_pos;
				q->pos = RB3_MALLOC(rb3_pos_t, q->n_pos);
				for (j = 0; j < q->n_pos; ++j)
					q->pos[j].pos = rb3_ssa(&p->fmi, p->fmi.ssa, q->mem.x[0] + j, &q->pos[j].sid);
			}
			#endif
		}
	}
}

static void worker_for_hapdiv(void *data, long i, int tid)
{
	step_t *t = (step_t*)data;
	const pipeline_t *p = t->p;
	m_hapdiv_t *a = &t->hapdiv[i];
	rb3_hapdiv(t->buf[tid].km, &p->opt->swo, &p->fmi, p->opt->hapdiv_k, &t->seq[a->id].seq[a->offset], &a->r);
}

static inline void write_name(kstring_t *out, const m_seq_t *s)
{
	if (s->name) rb3_sprintf_lite(out, "%s", s->name);
	else rb3_sprintf_lite(out, "seq%ld", s->id + 1);
}

static void pos_stranded(const rb3_sid_t *sid, const rb3_pos_t *pos, int32_t rlen, int64_t *clen, int64_t *st, int64_t *en)
{
	*clen = sid->len[pos->sid>>1];
	if ((pos->sid & 1) == 0)
		*st = pos->pos, *en = pos->pos + rlen;
	else
		*st = *clen - (pos->pos + rlen), *en = *clen - pos->pos;
}

/********************************************************************
 * refmap: place a query on a designated reference genome by
 * walking outward through the carrier genomes to the breakpoints.
 ********************************************************************/

// Backward-search the whole query; on success *out holds its SA interval. Returns 1 on a full match.
static int refmap_query_interval(const rb3_fmi_t *f, int64_t len, const uint8_t *q, rb3_sai_t *out)
{
	rb3_sai_t ik, ok[RB3_ASIZE];
	int64_t i;
	if (len <= 0) return 0;
	rb3_fmd_set_intv(f, q[len-1], &ik);
	if (ik.size == 0) return 0;
	for (i = len - 2; i >= 0; --i) {
		if (q[i] == 0 || q[i] > 4) return 0; // sentinel or ambiguous base: give up on a full match
		rb3_fmd_extend(f, &ik, ok, 1);
		if (ok[q[i]].size == 0) return 0;
		ik = ok[q[i]];
	}
	*out = ik;
	return 1;
}

// Walk outward from interval Iq through the carriers, recording one base per step.
// side 0 = left (backward) flank, side 1 = right (forward) flank. On return, buf holds the flank
// in forward (5'->3') orientation; the breakpoint-proximal end (adjacent to the query) is the 3'
// end for the left flank and the 5' end (index 0) for the right flank. Returns the flank length.
// walk_mode picks the path: consensus follows the most-supported base; strict stops at the first
// carrier disagreement; per-carrier (mask != NULL) follows the single carrier marked in mask.
static int32_t refmap_extract_flank(void *km, const pipeline_t *p, const rb3_sai_t *Iq, int side, int32_t max_walk, int walk_mode, const uint8_t *mask, uint8_t *buf)
{
	const rb3_fmi_t *f = &p->fmi;
	rb3_sai_t I = *Iq, ok[RB3_ASIZE];
	int32_t n = 0, is_back = (side == 0);
	while (n < max_walk) {
		int c, best = -1;
		int64_t bestsz = 0;
		rb3_fmd_extend(f, &I, ok, is_back);
		if (mask) { // per-carrier: follow the base whose child still contains this carrier
			rb3_pos_t pos[1];
			for (c = 1; c <= 4; ++c) {
				rb3_sai_t *iv = is_back? &ok[c] : &ok[rb3_comp(c)];
				if (iv->size == 0) continue;
				if (rb3_ssa_multi_ref(km, f, f->ssa, iv->x[0], iv->x[0] + iv->size, mask, 1, 1<<16, pos) > 0) {
					best = c, bestsz = iv->size; break;
				}
			}
		} else {
			for (c = 1; c <= 4; ++c) { // consensus: most-supported neighboring base
				int64_t sz = is_back? ok[c].size : ok[rb3_comp(c)].size;
				if (sz > bestsz) bestsz = sz, best = c;
			}
		}
		if (best < 0 || bestsz == 0) break;                         // no carrier continues
		if (walk_mode == RB3_WALK_STRICT && bestsz != I.size) break; // carriers disagree
		buf[n++] = best;
		I = is_back? ok[best] : ok[rb3_comp(best)];
	}
	if (side == 0) { // reverse the left flank so buf is in forward orientation
		int32_t i;
		for (i = 0; i < n>>1; ++i) {
			uint8_t t = buf[i]; buf[i] = buf[n-1-i]; buf[n-1-i] = t;
		}
	}
	return n;
}

// Given the query interval and (optionally) a carrier mask, extract both flanks, re-anchor them in
// the reference and fill the placement result r (status/strand/ref_sid/cL/cR/ins_size).
static void refmap_place(void *km, const pipeline_t *p, const m_seq_t *s, const rb3_sai_t *Iq, const uint8_t *mask, refmap_rst_t *r)
{
	const rb3_mopt_t *o = p->opt;
	uint8_t *Lf, *Rf;
	int32_t Ll, Rl, mlenL = 0, mlenR = 0;
	int64_t sidL = -1, sidR = -1, coordL = -1, coordR = -1;
	int gotL, gotR, strandL = 0, strandR = 0;
	Lf = Kmalloc(km, uint8_t, o->max_walk);
	Rf = Kmalloc(km, uint8_t, o->max_walk);
	Ll = refmap_extract_flank(km, p, Iq, 0, o->max_walk, o->walk_mode, mask, Lf);
	Rl = refmap_extract_flank(km, p, Iq, 1, o->max_walk, o->walk_mode, mask, Rf);
	gotL = refmap_anchor_flank(km, p, Lf, Ll, 0, &sidL, &strandL, &coordL, &mlenL);
	gotR = refmap_anchor_flank(km, p, Rf, Rl, 1, &sidR, &strandR, &coordR, &mlenR);
	kfree(km, Lf); kfree(km, Rf);
	int concord = (gotL && gotR && sidL == sidR && strandL == strandR);
	if (concord && o->max_bracket > 0) { // E1 collinearity: the two anchors must bracket a small zone
		int64_t lo = coordL < coordR? coordL : coordR, hi = coordL < coordR? coordR : coordL;
		if (hi - lo > o->max_bracket) concord = 0; // anchors too far apart -> paralogous, not collinear
	}
	if (concord) {
		r->status = RB3_RM_PLACED, r->ref_sid = sidL, r->strand = strandL;
		r->cL = coordL < coordR? coordL : coordR;
		r->cR = coordL < coordR? coordR : coordL;
		r->ins_size = (int64_t)(Ll - mlenL) + s->len + (Rl - mlenR);
	} else if (!o->two_flank && (gotL || gotR)) { // E1: with two_flank, a lone anchor is not trusted
		r->status = RB3_RM_ONE_SIDE;
		if (gotL) r->ref_sid = sidL, r->strand = strandL, r->cL = coordL;
		else      r->ref_sid = sidR, r->strand = strandR, r->cR = coordR;
	}
}

// Re-anchor a flank in the reference. The reference-shared part of a flank lies at the end FAR
// from the query (it is contiguous with the query only in the carriers), so a plain SMEM would be
// swallowed by the longer carrier match that crosses the breakpoint. Instead we grow the match
// from the far end and keep the longest stretch still present in the reference:
//   side 0 = left flank  : longest reference-matching PREFIX  (forward extension from the 5' end)
//   side 1 = right flank : longest reference-matching SUFFIX  (backward extension from the 3' end)
// The breakpoint is the near edge of that reference stretch: the prefix end for the left flank
// (*out_coord = forward end), the suffix start for the right flank (*out_coord = forward start).
// *out_mlen returns the matched reference length (used to size the insertion).
static int refmap_anchor_flank(void *km, const pipeline_t *p, const uint8_t *flank, int32_t flen, int side,
							   int64_t *out_sid, int *out_strand, int64_t *out_coord, int32_t *out_mlen)
{
	const rb3_fmi_t *f = &p->fmi;
	rb3_sai_t ik, ok[RB3_ASIZE], best_iv;
	rb3_pos_t pos[1];
	int32_t i, best_len = 0, min_len = p->opt->min_len;
	int c;
	int64_t clen, rst, ren;
	if (flen <= 0) return 0;
	memset(&best_iv, 0, sizeof(best_iv));
	if (side == 1) { // right flank: extend the suffix leftward from the 3' end
		c = flank[flen-1];
		if (c < 1 || c > 4) return 0;
		rb3_fmd_set_intv(f, c, &ik);
		for (i = flen - 1; i >= 0; --i) {
			if (i < flen - 1) {
				c = flank[i];
				if (c < 1 || c > 4) break;
				rb3_fmd_extend(f, &ik, ok, 1);
				if (ok[c].size == 0) break;
				ik = ok[c];
			}
			if (rb3_ssa_multi_ref(km, f, f->ssa, ik.x[0], ik.x[0] + ik.size, p->is_ref, 1, 1<<16, pos) > 0)
				best_len = flen - i, best_iv = ik;
			else if (best_len > 0) break; // reference dropped out: the longest stretch ends here
		}
	} else { // left flank: extend the prefix rightward from the 5' end
		c = flank[0];
		if (c < 1 || c > 4) return 0;
		rb3_fmd_set_intv(f, c, &ik);
		for (i = 0; i < flen; ++i) {
			if (i > 0) {
				int cc = rb3_comp(flank[i]);
				if (flank[i] < 1 || flank[i] > 4) break;
				rb3_fmd_extend(f, &ik, ok, 0);
				if (ok[cc].size == 0) break;
				ik = ok[cc];
			}
			if (rb3_ssa_multi_ref(km, f, f->ssa, ik.x[0], ik.x[0] + ik.size, p->is_ref, 1, 1<<16, pos) > 0)
				best_len = i + 1, best_iv = ik;
			else if (best_len > 0) break;
		}
	}
	if (best_len < min_len) return 0; // too short to anchor confidently
	// E2: the anchor must be (near-)single-copy across taxa, not a repeat. best_iv is the longest
	// reference-matching stretch, so it has the smallest interval; if even that exceeds max_occ the
	// flank is stuck in a repeat and the located coordinate cannot be trusted.
	if (p->opt->max_occ > 0 && best_iv.size > p->opt->max_occ) return 0;
	if (rb3_ssa_multi_ref(km, f, f->ssa, best_iv.x[0], best_iv.x[0] + best_iv.size, p->is_ref, 1, 1<<16, pos) <= 0)
		return 0;
	pos_stranded(f->sid, &pos[0], best_len, &clen, &rst, &ren);
	*out_sid = pos[0].sid >> 1;
	*out_strand = pos[0].sid & 1;
	// The breakpoint is the inner edge of the reference-matching region: the 3' end of the prefix
	// for the left flank, the 5' end of the suffix for the right flank. On a reverse-strand hit the
	// flank's increasing coordinate runs against the reference's, so the two ends swap.
	if (side == 0) *out_coord = (pos[0].sid & 1)? rst : ren; // left flank: end of prefix
	else           *out_coord = (pos[0].sid & 1)? ren : rst; // right flank: start of suffix
	*out_mlen = best_len;
	return 1;
}

// Place a carrier-only query by PROJECTING its carrier hits to the reference via
// the liftover (the E4 "second SSA"), instead of walking outward. The majority
// reference sequence among the per-carrier projections wins; NULL projections
// (no confident collinear support) are dropped, and if none survive -> UNPLACED.
static void refmap_place_lift(void *km, const pipeline_t *p, const m_seq_t *s, refmap_rst_t *r)
{
	const rb3_fmi_t *f = &p->fmi;
	int64_t rsids[8], rposs[8], best_rsid = -1, med, v[8];
	int32_t n = 0, k, j, bestn, m;
	for (k = 0; k < r->n_car_list && n < 8; ++k) {
		int64_t clen, st, en, rsid, rpos;
		pos_stranded(f->sid, &r->carriers[k], s->len, &clen, &st, &en);
		if (rb3_lift_project(p->lift, km, (int32_t)(r->carriers[k].sid >> 1), st,
							 p->opt->lift_win, p->opt->lift_mad, 4, &rsid, &rpos))
			rsids[n] = rsid, rposs[n] = rpos, ++n;
	}
	if (n == 0) { r->status = RB3_RM_UNPLACED; return; }
	for (k = 0, bestn = 0; k < n; ++k) {           // majority reference sequence
		int32_t c = 0;
		for (j = 0; j < n; ++j) if (rsids[j] == rsids[k]) ++c;
		if (c > bestn) bestn = c, best_rsid = rsids[k];
	}
	for (k = 0, m = 0; k < n; ++k) if (rsids[k] == best_rsid) v[m++] = rposs[k];
	for (k = 1; k < m; ++k) { int64_t x = v[k]; for (j = k - 1; j >= 0 && v[j] > x; --j) v[j+1] = v[j]; v[j+1] = x; }
	med = v[m >> 1];
	r->status = RB3_RM_PLACED, r->ref_sid = best_rsid, r->strand = 0;
	r->cL = r->cR = med, r->ins_size = 0;
}

// One reference-coordinate vote from a k-mer.
typedef struct { int64_t rsid, rpos; int32_t kmer; } refmap_vote_t;

static int refmap_vote_cmp(const void *a, const void *b)
{
	const refmap_vote_t *x = (const refmap_vote_t*)a, *y = (const refmap_vote_t*)b;
	if (x->rsid != y->rsid) return x->rsid < y->rsid? -1 : 1;
	if (x->rpos != y->rpos) return x->rpos < y->rpos? -1 : 1;
	return 0;
}

// Map the k-mer-agreement signals to a Phred-scaled MAPQ. The number of agreeing
// k-mers is a well-calibrated, error-rate-robust precision predictor (measured on
// maize NAM, 0-1% substitution: 1->~0.77, 2->~0.87, 3->~0.92, 4->~0.95, 5->~0.96,
// >=6->~0.97). A runner-up cluster that ties the winner is an ambiguous locus and
// caps precision (measured ~0.85). MAPQ = -10*log10(1 - precision). The raw agree/
// second columns are emitted too, for callers who prefer their own thresholds.
static inline int refmap_kmer_mapq(int agree, int second)
{
	static const double prec[7] = { 0.0, 0.77, 0.87, 0.92, 0.95, 0.96, 0.97 };
	double p = prec[agree < 1? 1 : (agree > 6? 6 : agree)];
	if (second >= agree && p > 0.85) p = 0.85;   // tie with a runner-up locus -> ambiguous
	if (p > 0.99999) p = 0.99999;
	return (int)(-10.0 * log10(1.0 - p) + 0.499);
}

// Collect the reference-coordinate votes of one k-mer (only if it is a full-length
// exact, single-copy-per-taxon match): a reference hit votes directly; a carrier
// hit votes via the liftover projection. An error-containing k-mer has no
// full-length match and contributes nothing.
static void refmap_kmer_votes(void *km, const pipeline_t *p, const uint8_t *q, int32_t K, int32_t ki,
							  refmap_vote_t **votes, int64_t *nv, int64_t *mv, rb3_pos_t *pos, int64_t cap)
{
	const rb3_fmi_t *f = &p->fmi;
	const rb3_mopt_t *o = p->opt;
	rb3_sai_t Iq;
	int64_t np, j;
	if (!refmap_query_interval(f, K, q, &Iq) || Iq.size < 1 || Iq.size > cap) return;
	np = rb3_ssa_multi(km, f, f->ssa, Iq.x[0], Iq.x[0] + Iq.size, cap, pos);
	for (j = 0; j < np; ++j) {
		int64_t sx = pos[j].sid >> 1, clen, st, en, rsid, rpos;
		pos_stranded(f->sid, &pos[j], K, &clen, &st, &en); // st = forward start of the k-mer
		if (p->is_ref[sx]) rsid = sx, rpos = st;
		else if (!(p->lift && rb3_lift_project(p->lift, km, (int32_t)sx, st, o->lift_win, o->lift_mad, 4, &rsid, &rpos)))
			continue;
		Kgrow(km, refmap_vote_t, *votes, *nv, *mv);
		(*votes)[*nv].rsid = rsid, (*votes)[*nv].rpos = rpos, (*votes)[*nv].kmer = ki, (*nv)++;
	}
}

// Place a read from its k-mers: tile it, project each k-mer to the reference, and
// place only where >= min_agree DISTINCT k-mers agree on a locus (cluster within
// kmer_cluster on one reference sequence). Tolerates sequencing error (an
// error-free k-mer votes even when the whole read fails to match exactly) and
// raises precision (a lone paralogous k-mer is outvoted). See docs experiment E4.
static void refmap_query_kmer(void *km, const pipeline_t *p, const m_seq_t *s, refmap_rst_t *r)
{
	const rb3_mopt_t *o = p->opt;
	int32_t K = o->kmer_len, ki = 0;
	int64_t off, last_off, nv = 0, mv = 0, cap = o->max_occ > 0? o->max_occ : 8, i;
	int64_t best_support = 0, second_support = 0, best_rsid = -1, best_pos = -1;
	uint64_t all_mask = 0;
	refmap_vote_t *votes = 0;
	rb3_pos_t *pos = Kmalloc(km, rb3_pos_t, cap);

	memset(r, 0, sizeof(*r));
	r->status = RB3_RM_UNPLACED, r->qlen = s->len, r->ref_sid = -1, r->cL = r->cR = -1;
	if (s->len < K) { kfree(km, pos); return; }

	for (off = 0; off + K <= s->len; off += o->kmer_step)
		refmap_kmer_votes(km, p, s->seq + off, K, ki++, &votes, &nv, &mv, pos, cap);
	last_off = (s->len - K) / o->kmer_step * o->kmer_step;
	if (last_off != s->len - K)                     // ensure the 3' end is tiled
		refmap_kmer_votes(km, p, s->seq + (s->len - K), K, ki++, &votes, &nv, &mv, pos, cap);
	kfree(km, pos);

	// n_vote = distinct k-mers that cast any vote (informative tiles)
	for (i = 0; i < nv; ++i) if (votes[i].kmer < 64) all_mask |= 1ULL << votes[i].kmer;
	r->n_vote = __builtin_popcountll(all_mask);

	// cluster the votes; support = number of DISTINCT k-mers in the cluster. Track the
	// top-2 cluster supports (best drives the placement; second measures competition).
	qsort(votes, nv, sizeof(refmap_vote_t), refmap_vote_cmp);
	for (i = 0; i < nv; ) {
		int64_t j = i, rs = votes[i].rsid, start = votes[i].rpos, support = 0;
		uint64_t kmask = 0;
		for (; j < nv && votes[j].rsid == rs && votes[j].rpos - start <= o->kmer_cluster; ++j)
			if (votes[j].kmer < 64 && !(kmask >> votes[j].kmer & 1))
				kmask |= 1ULL << votes[j].kmer, ++support;
		if (support > best_support)
			second_support = best_support, best_support = support, best_rsid = rs, best_pos = votes[(i + j) >> 1].rpos;
		else if (support > second_support)
			second_support = support;
		i = j;
	}
	kfree(km, votes);
	r->agree = (int32_t)best_support, r->second = (int32_t)second_support;
	if (best_support >= o->min_agree) {
		r->status = RB3_RM_PLACED, r->ref_sid = best_rsid, r->strand = 0;
		r->cL = r->cR = best_pos, r->ins_size = 0, r->n_carrier = (int32_t)best_support;
		r->mapq = refmap_kmer_mapq((int32_t)best_support, (int32_t)second_support);
	}
}

static void refmap_query(void *km, const pipeline_t *p, const m_seq_t *s, refmap_rst_t *r)
{
	const rb3_fmi_t *f = &p->fmi;
	rb3_sai_t Iq;
	if (p->opt->kmer_len > 0) { refmap_query_kmer(km, p, s, r); return; }
	rb3_pos_t *pos;
	int64_t np, i;
	int car_seen[8], n_car = 0;

	memset(r, 0, sizeof(*r));
	r->status = RB3_RM_UNPLACED, r->strand = 0, r->qlen = s->len;
	r->ref_sid = -1, r->cL = r->cR = -1;

	if (!refmap_query_interval(f, s->len, s->seq, &Iq) || Iq.size == 0) {
		// no end-to-end match (e.g. a sequencing error): fall back to the longest exact core (SMEM)
		rb3_sai_v mem = {0,0,0};
		size_t bi;
		int32_t bestlen = 0;
		rb3_fmd_smem_TG(km, f, s->len, s->seq, &mem, 1, p->opt->min_len);
		for (bi = 0; bi < mem.n; ++bi) {
			int32_t st = mem.a[bi].info>>32, en = (int32_t)mem.a[bi].info;
			if (en - st > bestlen) bestlen = en - st, Iq = mem.a[bi];
		}
		kfree(km, mem.a);
		if (bestlen == 0 || Iq.size == 0) return; // nothing of the query occurs in any genome
	}

	// E2: an informative read maps at most once per taxon; a read occurring > max_occ times is a
	// repeat/retro and cannot be placed confidently (it would be reported at one arbitrary copy).
	if (p->opt->max_occ > 0 && Iq.size > p->opt->max_occ) {
		r->status = RB3_RM_MULTI;
		return;
	}

	// locate a sample of occurrences; separate reference hits (exact) from carriers
	pos = Kmalloc(km, rb3_pos_t, 64);
	np = rb3_ssa_multi(km, f, f->ssa, Iq.x[0], Iq.x[0] + Iq.size, 64, pos);
	for (i = 0; i < np; ++i) {
		if (p->is_ref[pos[i].sid>>1]) { // the query is present in the reference: report it directly
			int64_t clen, st, en;
			pos_stranded(f->sid, &pos[i], s->len, &clen, &st, &en);
			r->status = RB3_RM_EXACT, r->ref_sid = pos[i].sid>>1, r->strand = pos[i].sid&1;
			r->cL = st, r->cR = en, r->ins_size = 0;
			kfree(km, pos);
			return;
		}
	}
	r->carriers = RB3_CALLOC(rb3_pos_t, 8);
	for (i = 0; i < np && n_car < 8; ++i) { // keep a few distinct carrier sequences for reporting
		int32_t j, dup = 0;
		for (j = 0; j < n_car; ++j)
			if (car_seen[j] == (int)(pos[i].sid>>1)) { dup = 1; break; }
		if (dup) continue;
		car_seen[n_car] = pos[i].sid>>1;
		r->carriers[n_car++] = pos[i];
	}
	r->n_car_list = n_car, r->n_carrier = n_car;
	kfree(km, pos);

	if (p->lift) { // E4: project carrier hits to the reference instead of walking
		refmap_place_lift(km, p, s, r);
		return;
	}

	if (p->opt->walk_mode == RB3_WALK_PERCARRIER && n_car > 0) {
		// place each carrier separately, following that one carrier through divergences
		uint8_t *mask = RB3_CALLOC(uint8_t, f->sid->n_seq);
		int32_t k;
		r->sub = RB3_CALLOC(refmap_rst_t, n_car);
		r->n_sub = n_car;
		for (k = 0; k < n_car; ++k) {
			refmap_rst_t *sub = &r->sub[k];
			int64_t csid = r->carriers[k].sid >> 1;
			memset(sub, 0, sizeof(*sub));
			sub->status = RB3_RM_UNPLACED, sub->qlen = s->len, sub->ref_sid = -1, sub->cL = sub->cR = -1;
			sub->carriers = &r->carriers[k], sub->n_car_list = 1, sub->n_carrier = 1; // borrowed pointer
			mask[csid] = 1;
			refmap_place(km, p, s, &Iq, mask, sub);
			mask[csid] = 0;
		}
		free(mask);
		return;
	}
	// consensus / strict: a single placement from the shared carrier path
	refmap_place(km, p, s, &Iq, 0, r);
}

static void write_paf(kstring_t *out, const rb3_fmi_t *f, const rb3_swhit_t *h, const m_seq_t *s)
{
	int32_t k;
	write_name(out, s);
	rb3_sprintf_lite(out, "\t%d\t%d\t%d", s->len, h->qoff[0], h->qoff[0] + h->qlen);
	if (h->n_pos > 0) {
		int64_t sid = h->pos[0].sid, pos = h->pos[0].pos;
		if (f->sid) { // print with sequence names and lengths
			int64_t clen, st, en;
			pos_stranded(f->sid, &h->pos[0], h->rlen, &clen, &st, &en);
			rb3_sprintf_lite(out, "\t%c\t%s\t%ld\t%ld\t%ld", "+-"[sid&1], f->sid->name[sid>>1], (long)clen, st, en);
		} else {
			rb3_sprintf_lite(out, "\t+\t%ld\t*\t%ld\t%ld", sid, pos, pos + h->rlen); // always on the forward strand
		}
	} else {
		rb3_sprintf_lite(out, "\t*\t*\t%d\t*\t*", h->rlen);
	}
	rb3_sprintf_lite(out, "\t%d\t%d\t0", h->mlen, h->blen);
	rb3_sprintf_lite(out, "\tAS:i:%d\tqh:i:%d\trh:i:%ld\tcg:Z:", h->score, h->n_qoff, (long)(h->hi - h->lo));
	for (k = 0; k < h->n_cigar; ++k)
		rb3_sprintf_lite(out, "%d%c", h->cigar[k]>>4, "MIDNSHP=X"[h->cigar[k]&0xf]);
	rb3_sprintf_lite(out, "\tcs:Z:%s", h->cs);
	if (h->rseq) {
		rb3_sprintf_lite(out, "\trs:Z:");
		for (k = 0; k < h->rlen; ++k)
			rb3_sprintf_lite(out, "%c", "$ACGTN"[h->rseq[k]]);
	}
	if (h->n_pos > 1) {
		rb3_sprintf_lite(out, "\ta%c:Z:", f->sid? 'p' : 'q');
		for (k = 1; k < h->n_pos; ++k) {
			int64_t sid = h->pos[k].sid, pos = h->pos[k].pos;
			if (f->sid) {
				int64_t clen, st, en;
				pos_stranded(f->sid, &h->pos[k], h->rlen, &clen, &st, &en);
				rb3_sprintf_lite(out, "%s,%c,%ld;", f->sid->name[sid>>1], "+-"[sid&1], st);
			} else {
				rb3_sprintf_lite(out, "%ld,%ld;", sid, pos);
			}
		}
	}
	rb3_sprintf_lite(out, "\n");
}

static void write_all_hits(kstring_t *out, const m_seq_t *s, const rb3_swrst_t *r, char strand, int64_t max_all_out)
{
	int64_t n_out = 0, tot = 0;
	int32_t i;
	if (max_all_out <= 0) max_all_out = INT64_MAX;
	for (i = 0; i < r->n; ++i) tot += r->a[i].hi - r->a[i].lo;
	for (i = 0; i < r->n; ++i) {
		n_out += r->a[i].hi - r->a[i].lo;
		if (n_out >= max_all_out) break;
	}
	rb3_sprintf_lite(out, "QS\t");
	write_name(out, s);
	rb3_sprintf_lite(out, "\t%d\t%d\t%c\t%ld\t%ld\n", s->len, r->n, strand, n_out, tot);
	for (i = 0, n_out = 0; i < r->n; ++i) {
		const rb3_swhit_t *h = &r->a[i];
		rb3_sprintf_lite(out, "QH\t%ld\t%d\t%d\t%s\n", (long)(h->hi - h->lo), h->score, h->blen - h->mlen, h->cs);
		n_out += h->hi - h->lo;
		if (n_out >= max_all_out) break;
	}
	rb3_sprintf_lite(out, "//\n");
}

static void write_refmap1(kstring_t *out, const rb3_fmi_t *f, const m_seq_t *s, const refmap_rst_t *r, int kmer)
{
	static const char *status_str[5] = { "UNPLACED", "PLACED", "ONE_SIDE", "EXACT", "MULTI" };
	int32_t k;
	out->l = 0;
	write_name(out, s);
	rb3_sprintf_lite(out, "\t%d\t%s\t%d\t", r->qlen, status_str[(int)r->status], r->n_carrier);
	if (r->n_car_list > 0) { // carrier sequences (name:strand)
		for (k = 0; k < r->n_car_list; ++k) {
			int64_t sid = r->carriers[k].sid;
			rb3_sprintf_lite(out, "%s%s:%c", k? "," : "", f->sid->name[sid>>1], "+-"[sid&1]);
		}
	} else rb3_sprintf_lite(out, ".");
	if (r->ref_sid >= 0)
		rb3_sprintf_lite(out, "\t%s\t%c", f->sid->name[r->ref_sid], "+-"[r->strand & 1]);
	else
		rb3_sprintf_lite(out, "\t.\t.");
	if (r->cL >= 0) rb3_sprintf_lite(out, "\t%ld", (long)r->cL); else rb3_sprintf_lite(out, "\t.");
	if (r->cR >= 0) rb3_sprintf_lite(out, "\t%ld", (long)r->cR); else rb3_sprintf_lite(out, "\t.");
	if (r->status == RB3_RM_PLACED || r->status == RB3_RM_EXACT)
		rb3_sprintf_lite(out, "\t%ld\t%ld", (long)(r->cR - r->cL), (long)r->ins_size);
	else
		rb3_sprintf_lite(out, "\t.\t.");
	if (kmer) // --kmer confidence: informative tiles, agreeing k-mers, runner-up, calibrated MAPQ
		rb3_sprintf_lite(out, "\t%d\t%d\t%d\t%d", r->n_vote, r->agree, r->second, r->mapq);
	rb3_sprintf_lite(out, "\n");
}

static void write_refmap(step_t *t)
{
	const pipeline_t *p = t->p;
	const rb3_fmi_t *f = &p->fmi;
	int32_t j, k;
	kstring_t out = {0,0,0};
	for (j = 0; j < t->n_seq; ++j) {
		m_seq_t *s = &t->seq[j];
		refmap_rst_t *r = &t->refmap[j];
		int kmer = p->opt->kmer_len > 0;
		if (r->n_sub > 0) { // per-carrier mode: one line per carrier (sub->carriers borrows r->carriers)
			for (k = 0; k < r->n_sub; ++k) {
				write_refmap1(&out, f, s, &r->sub[k], kmer);
				fputs(out.s, stdout);
			}
			free(r->sub);
		} else {
			write_refmap1(&out, f, s, r, kmer);
			fputs(out.s, stdout);
		}
		free(r->carriers);
		free(s->seq);
		free(s->name);
	}
	free(out.s);
}

static void write_per_seq(step_t *t)
{
	const pipeline_t *p = t->p;
	int32_t i, j;
	kstring_t out = {0,0,0};
	for (j = 0; j < t->n_seq; ++j) {
		m_seq_t *s = &t->seq[j];
		free(s->seq);
		out.l = 0;
		if (p->opt->algo == RB3_SA_SW && (p->opt->flag & RB3_MF_WRITE_ALL)) { // write all hits in a compact format
			write_all_hits(&out, s, &t->rst[j], '+', p->opt->max_all_out);
			rb3_swrst_free(&t->rst[j]);
			if (t->rst_rev) {
				write_all_hits(&out, s, &t->rst_rev[j], '-', p->opt->max_all_out);
				rb3_swrst_free(&t->rst_rev[j]);
			}
			fputs(out.s, stdout);
		} else if (p->opt->algo == RB3_SA_SW) { // write PAF
			rb3_swrst_t *r = &t->rst[j];
			if (r->n > 0) { // mapped
				for (i = 0; i < r->n; ++i) {
					out.l = 0;
					write_paf(&out, &p->fmi, &r->a[i], s);
					fputs(out.s, stdout);
				}
			} else if (p->opt->flag & RB3_MF_WRITE_UNMAP) { // unmapped
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\n", s->len);
				fputs(out.s, stdout);
			}
			rb3_swrst_free(r);
		} else if (p->opt->min_gap_len > 0) { // output regions not covered by long MEMs
			for (i = 0; i < s->n_gap; ++i) {
				int32_t st = s->gap[i]>>32, en = (int32_t)s->gap[i];
				out.l = 0;
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t%d\t%d\n", st, en, s->len);
				fputs(out.s, stdout);
			}
		} else if (p->opt->flag & RB3_MF_WRITE_COV) { // output breadth of coverage
			int32_t st0 = 0, en0 = 0, cov = 0;
			for (i = 0; i < s->n_mem; ++i) {
				rb3_sai_t *q = &s->mem[i].mem;
				int32_t st = q->info>>32, en = (int32_t)q->info;
				if (st > en0) {
					cov += en0 - st0;
					st0 = st, en0 = en;
				} else en0 = en0 > en? en0 : en;
			}
			cov += en0 - st0;
			if (cov > 0) {
				out.l = 0;
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t%d\n", s->len, cov);
				fputs(out.s, stdout);
			}
		} else { // output long MEMs
			const rb3_fmi_t *f = &p->fmi;
			for (i = 0; i < s->n_mem; ++i) {
				m_sai_pos_t *r = &s->mem[i];
				rb3_sai_t *q = &r->mem;
				int32_t st = q->info>>32, en = (int32_t)q->info;
				out.l = 0;
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t%d\t%ld", st, en, (long)q->size);
				if (r->n_pos > 0) {
					int32_t j;
					rb3_sprintf_lite(&out, "\t%ld", r->n_pos);
					for (j = 0; j < r->n_pos; ++j) {
						rb3_pos_t *t = &r->pos[j];
						int64_t rlen = f->sid->len[t->sid>>1], pos;
						pos = t->sid&1? rlen - (t->pos + (en - st)) : t->pos;
						rb3_sprintf_lite(&out, "\t%s:%c:%ld", f->sid->name[t->sid>>1], "+-"[t->sid&1], pos);
					}
					free(r->pos);
				}
				rb3_sprintf_lite(&out, "\n");
				fputs(out.s, stdout);
			}
		}
		free(s->name); free(s->mem); free(s->gap);
	}
	free(out.s);
	free(t->rst);
	free(t->rst_rev);
}

static void write_hapdiv(step_t *t)
{
	int32_t j, ed;
	const m_hapdiv_t *p;
	kstring_t out = {0,0,0};
	for (j = 0; j < t->n_seq; ++j)
		free(t->seq[j].seq);
	if (t->n_hapdiv == 0) return;
	p = &t->hapdiv[0];
	for (j = 1; j <= t->n_hapdiv; ++j) {
		const m_hapdiv_t *q = t->hapdiv + j;
		if (j == t->n_hapdiv || p->id != q->id || memcmp(&p->r, &q->r, sizeof(p->r)) != 0) {
			m_seq_t *s = &t->seq[p->id];
			out.l = 0;
			write_name(&out, s);
			rb3_sprintf_lite(&out, "\t%d\t%d\t%d\t%d", p->offset, t->hapdiv[j-1].offset + t->p->opt->hapdiv_k, p->r.n_al, p->r.max_ed);
			for (ed = 0; ed <= RB2_SW_MAX_ED; ++ed)
				rb3_sprintf_lite(&out, "\t%d", p->r.n_hap[ed]);
			puts(out.s);
			p = q;
		}
	}
	for (j = 0; j < t->n_seq; ++j)
		free(t->seq[j].name);
	free(out.s);
	free(t->hapdiv);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	pipeline_t *p = (pipeline_t*)shared;
	step_t *t = (step_t*)in;
	int32_t i;
	if (step == 0) {
		const char *name;
		char *ss;
		int64_t len, tot = 0;
		int32_t n_seq = 0, m_seq = 0;
		m_seq_t *seq = 0;
		while ((ss = rb3_seq_read1(p->fp, &len, &name)) != 0) { // read sequences
			m_seq_t *s;
			RB3_GROW0(m_seq_t, seq, n_seq, m_seq);
			s = &seq[n_seq++];
			s->name = name? rb3_strdup(name) : 0;
			s->seq = (uint8_t*)rb3_strdup(ss);
			s->len = len;
			s->id = p->id++;
			s->mem = 0, s->n_mem = 0;
			tot += len;
			if (tot >= p->opt->batch_size)
				break;
		}
		if (n_seq > 0) { // construct a step_t object
			t = RB3_CALLOC(step_t, 1);
			t->p = p;
			t->seq = seq;
			t->n_seq = n_seq;
			if (p->opt->algo == RB3_SA_HAPDIV) { // the hapdiv mode
				int32_t j, n_hapdiv = 0;
				for (i = 0; i < n_seq; ++i)
					n_hapdiv += seq[i].len < p->opt->hapdiv_k? 0 : (seq[i].len - p->opt->hapdiv_k) / p->opt->hapdiv_w + 1;
				t->n_hapdiv = n_hapdiv;
				t->hapdiv = RB3_CALLOC(m_hapdiv_t, n_hapdiv);
				for (i = 0, n_hapdiv = 0; i < n_seq; ++i)
					for (j = 0; j + p->opt->hapdiv_k <= seq[i].len; j += p->opt->hapdiv_w)
						t->hapdiv[n_hapdiv].id = i, t->hapdiv[n_hapdiv++].offset = j;
				assert(n_hapdiv == t->n_hapdiv);
			} else if (p->opt->algo == RB3_SA_REFMAP) { // reference-placement mode
				t->refmap = RB3_CALLOC(refmap_rst_t, n_seq);
			} else { // per-sequence mode (sw, mem, gap and coverage)
				t->rst = RB3_CALLOC(rb3_swrst_t, n_seq);
				if (p->opt->flag & RB3_MF_BOTH_DIR)
					t->rst_rev = RB3_CALLOC(rb3_swrst_t, n_seq);
			}
			t->buf = RB3_CALLOC(m_tbuf_t, p->opt->n_threads);
			for (i = 0; i < p->opt->n_threads; ++i)
				t->buf[i].km = p->opt->flag & RB3_MF_NO_KALLOC? 0 : km_init();
			return t;
		}
	} else if (step == 1) {
		if (p->opt->algo == RB3_SA_HAPDIV)
			kt_for(p->opt->n_threads, worker_for_hapdiv, in, t->n_hapdiv);
		else
			kt_for(p->opt->n_threads, worker_for_seq, in, t->n_seq);
		return in;
	} else if (step == 2) {
		for (i = 0; i < p->opt->n_threads; ++i) {
			kfree(t->buf[i].km, t->buf[i].mem.a);
			km_destroy(t->buf[i].km);
		}
		free(t->buf);
		if (p->opt->algo == RB3_SA_HAPDIV)
			write_hapdiv(t);
		else if (p->opt->algo == RB3_SA_REFMAP)
			write_refmap(t), free(t->refmap);
		else
			write_per_seq(t);
		free(t->seq);
		if (rb3_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences\n", __func__, rb3_realtime(), rb3_percent_cpu(), t->n_seq);
		free(t);
	}
	return 0;
}

static ko_longopt_t long_options[] = {
	{ "no-ssa",          ko_no_argument,       301 },
	{ "seq",             ko_no_argument,       302 },
	{ "gap",             ko_required_argument, 303 },
	{ "cov",             ko_no_argument,       304 },
	{ "old-mem",         ko_no_argument,       305 },
	{ "all-e2e",         ko_no_argument,       306 },
	{ "ref-prefix",      ko_required_argument, 307 },
	{ "max-walk",        ko_required_argument, 308 },
	{ "walk-mode",       ko_required_argument, 309 },
	{ "max-occ",         ko_required_argument, 310 },
	{ "two-flank",       ko_no_argument,       311 },
	{ "max-bracket",     ko_required_argument, 312 },
	{ "lift",            ko_required_argument, 313 },
	{ "lift-win",        ko_required_argument, 314 },
	{ "lift-mad",        ko_required_argument, 315 },
	{ "kmer",            ko_required_argument, 316 },
	{ "kmer-step",       ko_required_argument, 317 },
	{ "min-agree",       ko_required_argument, 318 },
	{ "kmer-cluster",    ko_required_argument, 319 },
	{ "no-kalloc",       ko_no_argument,       501 },
	{ "dbg-dawg",        ko_no_argument,       502 },
	{ "dbg-sw",          ko_no_argument,       503 },
	{ "dbg-qname",       ko_no_argument,       504 },
	{ "dbg-bt",          ko_no_argument,       505 },
	{ 0, 0, 0 }
};

int main_search(int argc, char *argv[]) // "sw" and "mem" share the same CLI
{
	int32_t c, j, is_line = 0, ret, load_flag = 0, no_ssa = 0;
	rb3_mopt_t opt;
	pipeline_t p;
	ketopt_t o = KETOPT_INIT;

	rb3_mopt_init(&opt);
	p.opt = &opt, p.id = 0;
	while ((c = ketopt(&o, argc, argv, 1, "Ll:c:t:K:MdN:A:B:O:E:C:m:k:uj:ey:a:w:p:bg:", long_options)) >= 0) {
		if (c == 'L') is_line = 1;
		else if (c == 'a') opt.algo = RB3_SA_HAPDIV, opt.hapdiv_k = atoi(o.arg);
		else if (c == 'w') opt.algo = RB3_SA_HAPDIV, opt.hapdiv_w = atoi(o.arg);
		else if (c == 'd') opt.algo = RB3_SA_SW, load_flag |= RB3_LOAD_ALL;
		else if (c == 'l') opt.min_len = atol(o.arg);
		else if (c == 'c') opt.min_occ = atol(o.arg);
		else if (c == 'g') opt.max_all_out = atol(o.arg), opt.flag |= RB3_MF_WRITE_ALL, opt.swo.flag |= RB3_SWF_E2E, opt.swo.end_len = 1, no_ssa = 1;
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'K') opt.batch_size = rb3_parse_num(o.arg);
		else if (c == 'p') opt.max_pos = opt.swo.max_pos = atoi(o.arg);
		else if (c == 'N') opt.swo.n_best = atoi(o.arg);
		else if (c == 'M') load_flag |= RB3_LOAD_MMAP;
		else if (c == 'A') opt.swo.match = atoi(o.arg);
		else if (c == 'B') opt.swo.mis = atoi(o.arg);
		else if (c == 'O') opt.swo.gap_open = atoi(o.arg);
		else if (c == 'E') opt.swo.gap_ext = atoi(o.arg);
		else if (c == 'C') opt.swo.r2cache_size = rb3_parse_num(o.arg);
		else if (c == 'm') opt.swo.min_sc = atoi(o.arg);
		else if (c == 'k') opt.swo.end_len = atoi(o.arg);
		else if (c == 'j') opt.swo.min_mem_len = atoi(o.arg);
		else if (c == 'e') opt.swo.flag |= RB3_SWF_E2E, opt.swo.end_len = 1;
		else if (c == 'y') opt.swo.e2e_drop = atoi(o.arg);
		else if (c == 'u') opt.flag |= RB3_MF_WRITE_UNMAP;
		else if (c == 'b') opt.flag |= RB3_MF_BOTH_DIR;
		else if (c == 301) no_ssa = 1;
		else if (c == 302) opt.swo.flag |= RB3_SWF_KEEP_RS;
		else if (c == 303) opt.min_gap_len = rb3_parse_num(o.arg);
		else if (c == 304) opt.flag |= RB3_MF_WRITE_COV;
		else if (c == 305) opt.algo = RB3_SA_MEM_ORI;
		else if (c == 306) opt.flag |= RB3_MF_WRITE_ALL, opt.swo.flag |= RB3_SWF_E2E, opt.swo.end_len = 1, no_ssa = 1;
		else if (c == 307) opt.ref_prefix = o.arg;
		else if (c == 308) opt.max_walk = rb3_parse_num(o.arg);
		else if (c == 309) {
			if (strcmp(o.arg, "consensus") == 0) opt.walk_mode = RB3_WALK_CONSENSUS;
			else if (strcmp(o.arg, "strict") == 0) opt.walk_mode = RB3_WALK_STRICT;
			else if (strcmp(o.arg, "per-carrier") == 0) opt.walk_mode = RB3_WALK_PERCARRIER;
			else { fprintf(stderr, "ERROR: --walk-mode must be consensus, strict or per-carrier\n"); return 1; }
		}
		else if (c == 310) opt.max_occ = atol(o.arg);     // E2: occurrence cap; <0 = auto (#taxa)
		else if (c == 311) opt.two_flank = 1;             // E1: require both flanks to anchor
		else if (c == 312) opt.max_bracket = rb3_parse_num(o.arg); // E1: max |cR-cL| for a PLACED
		else if (c == 313) opt.lift_fn = o.arg;           // E4: liftover file (project, not walk)
		else if (c == 314) opt.lift_win = rb3_parse_num(o.arg);
		else if (c == 315) opt.lift_mad = rb3_parse_num(o.arg);
		else if (c == 316) opt.kmer_len = atoi(o.arg);    // k-mer-agreement placement (0 = off)
		else if (c == 317) opt.kmer_step = atoi(o.arg);
		else if (c == 318) opt.min_agree = atoi(o.arg);
		else if (c == 319) opt.kmer_cluster = rb3_parse_num(o.arg);
		else if (c == 501) opt.flag |= RB3_MF_NO_KALLOC;
		else if (c == 502) rb3_dbg_flag |= RB3_DBG_DAWG;
		else if (c == 503) rb3_dbg_flag |= RB3_DBG_SW;
		else if (c == 504) rb3_dbg_flag |= RB3_DBG_QNAME;
		else if (c == 505) rb3_dbg_flag |= RB3_DBG_BT;
		else {
			fprintf(stderr, "ERROR: unknown option\n");
			return 1;
		}
	}

	if (opt.min_gap_len > 0) opt.max_pos = 0;
	if (strcmp(argv[0], "sw") == 0) {
		opt.algo = RB3_SA_SW;
		if (!no_ssa) load_flag |= RB3_LOAD_ALL;
	} else if (strcmp(argv[0], "hapdiv") == 0) {
		opt.algo = RB3_SA_HAPDIV, opt.swo.end_len = 1;
	} else if (strcmp(argv[0], "mem") == 0) {
		if (opt.max_pos > 0)
			load_flag |= RB3_LOAD_ALL;
	} else if (strcmp(argv[0], "refmap") == 0) {
		opt.algo = RB3_SA_REFMAP;
		load_flag |= RB3_LOAD_ALL;
	}
	if (opt.algo == RB3_SA_HAPDIV)
		opt.swo.flag |= RB3_SWF_E2E | RB3_SWF_HAPDIV;

	if (argc - o.ind < 2) {
		fprintf(stdout, "Usage: ropebwt3 %s [options] <idx.fmr> <seq.fa> [...]\n", argv[0]);
		fprintf(stderr, "Options:\n");
		if (strcmp(argv[0], "mem") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -l INT      min MEM length [%ld]\n", (long)opt.min_len);
			fprintf(stderr, "  -c INT      min interval size [%ld]\n", (long)opt.min_occ);
			fprintf(stderr, "  --old-mem   use the original MEM algorithm (for testing)\n");
			fprintf(stderr, "  --gap=NUM   output regions >=NUM that are not covered by MEMs [%d]\n", opt.min_gap_len);
			fprintf(stderr, "  --cov       output breadth of coverage\n");
		}
		if (strcmp(argv[0], "refmap") == 0) {
			fprintf(stderr, "  --ref-prefix=STR  reference = sequences whose name starts with STR [required]\n");
			fprintf(stderr, "  --max-walk=NUM    max bases to walk outward along carriers per flank [%d]\n", opt.max_walk);
			fprintf(stderr, "  --walk-mode=STR   carrier path: consensus|strict|per-carrier [consensus]\n");
			fprintf(stderr, "  --max-occ=INT     drop reads/anchors occurring >INT times; <0 = auto (#taxa); 0 = off [%ld]\n", (long)opt.max_occ);
			fprintf(stderr, "  --two-flank       require both flanks to anchor concordantly (drop ONE_SIDE)\n");
			fprintf(stderr, "  --max-bracket=NUM with --two-flank, max |cR-cL| for a PLACED; 0 = off [%ld]\n", (long)opt.max_bracket);
			fprintf(stderr, "  --lift=FILE       project carrier hits via a `ropebwt3 lift` map instead of walking\n");
			fprintf(stderr, "  --lift-win=NUM    liftover projection window [%ld]\n", (long)opt.lift_win);
			fprintf(stderr, "  --lift-mad=NUM    liftover max residual MAD [%ld]\n", (long)opt.lift_mad);
			fprintf(stderr, "  --kmer=INT        place a read from INT-bp k-mers by agreement (0 = off, whole-read)\n");
			fprintf(stderr, "  --kmer-step=INT   k-mer tiling step [%d]\n", opt.kmer_step);
			fprintf(stderr, "  --min-agree=INT   min agreeing k-mers to place [%d]\n", opt.min_agree);
			fprintf(stderr, "  --kmer-cluster=NUM  agreeing k-mers must fall within NUM bp [%ld]\n", (long)opt.kmer_cluster);
			fprintf(stderr, "  -l INT      min anchor length when re-mapping a flank [%ld]\n", (long)opt.min_len);
		}
		if (strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -d          use BWA-SW for local alignment\n");
		}
		if (strcmp(argv[0], "hapdiv") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -a INT      annotate sliding INT-mers [%d]\n", opt.hapdiv_k);
			fprintf(stderr, "  -w INT      k-mer step size for annotation [%d]\n", opt.hapdiv_w);
		}
		if (strcmp(argv[0], "sw") == 0 || strcmp(argv[0], "hapdiv") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -N INT      keep up to INT hits per DAWG node [%d]\n", opt.swo.n_best);
			fprintf(stderr, "  -m INT      min alignment score [%d]\n", opt.swo.min_sc);
			fprintf(stderr, "  -A INT      match score [%d]\n", opt.swo.match);
			fprintf(stderr, "  -B INT      mismatch penalty [%d]\n", opt.swo.mis);
			fprintf(stderr, "  -O INT      gap open penalty [%d]\n", opt.swo.gap_open);
			fprintf(stderr, "  -E INT      gap extension penalty; a k-long gap costs O+k*E [%d]\n", opt.swo.gap_ext);
			fprintf(stderr, "  -C NUM      size of the ranking cache [%d]\n", opt.swo.r2cache_size);
			fprintf(stderr, "  -y INT      ignore secondary hits scored INT lower than the best [%d]\n", opt.swo.e2e_drop);
		}
		if (strcmp(argv[0], "sw") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -e          end-to-end mode (forcing -k to 1)\n");
			fprintf(stderr, "  -j INT      min MEM length to initiate alignment [%d]\n", opt.swo.min_mem_len);
			fprintf(stderr, "  -k INT      require INT-mer match at the end of alignment [%d]\n", opt.swo.end_len);
			fprintf(stderr, "  -b          align both strands (effective with --all-e2e)\n");
			fprintf(stderr, "  -u          write unmapped queries to PAF\n");
			fprintf(stderr, "  --seq       write reference sequence to the rs tag\n");
			fprintf(stderr, "  --all-e2e   write all end-to-end hits in a compact format (forcing -e)\n");
			fprintf(stderr, "  -g INT      cap the number of --all-e2e output to INT (forcing --all-e2e)\n");
			fprintf(stderr, "  --no-ssa    ignore the sampled suffix array\n");
		}
		fprintf(stderr, "  -t INT      number of threads [%d]\n", opt.n_threads);
		fprintf(stderr, "  -p INT      output up to INT positions [%d]\n", opt.max_pos);
		fprintf(stderr, "  -L          one sequence per line in the input\n");
		fprintf(stderr, "  -K NUM      query batch size [100m]\n");
		fprintf(stderr, "  -M          use mmap to load FMD\n");
		return 0;
	}

	ret = rb3_fmi_load_all(&p.fmi, argv[o.ind], load_flag);
	if (ret < 0) return 1;
	if (opt.max_pos > 0 && (p.fmi.ssa == 0 || p.fmi.sid == 0)) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: failed to load suffix array samples or sequence names/lengths\n");
		return 1;
	}
	if (!rb3_fmi_is_symmetric(&p.fmi)) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: BWT doesn't contain both strands\n");
		return 1;
	}
	p.is_ref = 0, p.n_ref = 0, p.lift = 0;
	if (opt.algo == RB3_SA_REFMAP) { // mark the reference sequences by name prefix
		int64_t k, plen;
		if (opt.ref_prefix == 0) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: refmap requires --ref-prefix\n");
			return 1;
		}
		if (p.fmi.ssa == 0 || p.fmi.sid == 0) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: refmap needs the sampled suffix array (.ssa) and sequence names (.len.gz)\n");
			return 1;
		}
		plen = strlen(opt.ref_prefix);
		p.is_ref = RB3_CALLOC(uint8_t, p.fmi.sid->n_seq);
		for (k = 0; k < p.fmi.sid->n_seq; ++k)
			if (strncmp(p.fmi.sid->name[k], opt.ref_prefix, plen) == 0)
				p.is_ref[k] = 1, p.n_ref++;
		if (p.n_ref == 0) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: no sequence name starts with '%s'\n", opt.ref_prefix);
			free(p.is_ref);
			return 1;
		}
		if (rb3_verbose >= 3)
			fprintf(stderr, "[M::%s] %ld of %ld sequences marked as reference\n", __func__, (long)p.n_ref, (long)p.fmi.sid->n_seq);
		// count distinct taxa (sequence-name prefixes before '_') for the auto --max-occ default
		if (opt.max_occ < 0) {
			int64_t k, m, n_taxa = 0;
			char **pre = RB3_CALLOC(char*, p.fmi.sid->n_seq);
			for (k = 0; k < p.fmi.sid->n_seq; ++k) {
				const char *nm = p.fmi.sid->name[k];
				const char *us = strchr(nm, '_');
				int32_t plen = us? (int32_t)(us - nm) : (int32_t)strlen(nm), dup = 0;
				for (m = 0; m < n_taxa; ++m)
					if ((int32_t)strlen(pre[m]) == plen && strncmp(pre[m], nm, plen) == 0) { dup = 1; break; }
				if (!dup) { pre[n_taxa] = RB3_MALLOC(char, plen + 1); memcpy(pre[n_taxa], nm, plen); pre[n_taxa][plen] = 0; ++n_taxa; }
			}
			for (m = 0; m < n_taxa; ++m) free(pre[m]);
			free(pre);
			opt.max_occ = n_taxa;
			if (rb3_verbose >= 3)
				fprintf(stderr, "[M::%s] auto --max-occ = %ld (distinct taxa)\n", __func__, (long)opt.max_occ);
		}
		if (opt.lift_fn) { // E4: load the carrier->reference liftover
			p.lift = rb3_lift_restore(opt.lift_fn);
			if (p.lift == 0) {
				if (rb3_verbose >= 1) fprintf(stderr, "ERROR: failed to load liftover '%s'\n", opt.lift_fn);
				free(p.is_ref);
				return 1;
			}
			if (rb3_verbose >= 3)
				fprintf(stderr, "[M::%s] loaded liftover over %ld sequences\n", __func__, (long)rb3_lift_n_seq(p.lift));
		}
	}
	if (opt.flag & RB3_MF_WRITE_ALL) {
		puts("CC\tQS  queryName  queryLen  numHap");
		puts("CC\tQH  refCount   score     editDist   cs   strand   nOut   totAln");
		puts("CC");
	}
	for (j = o.ind + 1; j < argc; ++j) {
		p.fp = rb3_seq_open(argv[j], is_line);
		if (p.fp == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to load the sequence file '%s'\n", argv[j]);
			break;
		}
		kt_pipeline(2, worker_pipeline, &p, 3);
		rb3_seq_close(p.fp);
	}
	rb3_fmi_free(&p.fmi);
	free(p.is_ref);
	rb3_lift_destroy(p.lift);
	return 0;
}
