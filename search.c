#include <math.h>
#include "fm-index.h"
#include "align.h"
#include "rb3priv.h"
#include "io.h"
#include "ketopt.h"
#include "kthread.h"
#include "kalloc.h"
#include "lift.h"
#include "ps4g.h"
#include "hitcount.h"
#include "khashl-km.h"

KHASHL_MAP_INIT(KH_LOCAL, rb3_name2sid_t, rb3_name2sid, kh_cstr_t, int32_t, kh_hash_str, kh_eq_str)

typedef enum { RB3_SA_MEM_TG, RB3_SA_MEM_ORI, RB3_SA_SW, RB3_SA_HAPDIV, RB3_SA_REFMAP, RB3_SA_CHAIN } rb3_search_algo_t;

#define RB3_WALK_CONSENSUS  0 // follow the base shared by the most assemblies
#define RB3_WALK_STRICT     1 // stop walking at the first assembly disagreement
#define RB3_WALK_PERASSEMBLY 2 // one outward walk (and one result) per assembly genome

#define RB3_MF_NO_KALLOC   0x1
#define RB3_MF_WRITE_UNMAP 0x2
#define RB3_MF_WRITE_COV   0x4
#define RB3_MF_WRITE_ALL   0x8
#define RB3_MF_BOTH_DIR    0x10

// Max distinct assembly (assembly) sequences kept per PLACED read, for both the
// refmap_place_lift() liftover projection and the PS4G/npy gamete set. Was
// hardcoded to 8 (an arbitrary suffix-array-traversal-order subset once a
// locus is shared by more assemblies than that -- see the Oh43 real-data eval,
// where it silently dropped the read's own true assembly from ~11% of sites).
// Sized to the rb3_ssa_multi() locate cap so it never truncates below the locate.
// This should scale with the #samples in the index (a locus can be shared by up to
// N assemblies); capped at 256 to keep counts byte-sized (Ed). For N>256 the extra
// assemblies are truncated (accepted). See docs/sample-relative-thresholds.md (#3).
#define RB3_RM_MAX_ASSEMBLY 256

typedef struct {
	uint32_t flag;
	int32_t n_threads, min_gap_len, hapdiv_k, hapdiv_w;
	int32_t max_pos;
	rb3_search_algo_t algo;
	int64_t min_occ, min_len, max_all_out;
	int64_t batch_size;
	char *ref_prefix;  // refmap: reference sequences are those whose name starts with this
	int32_t max_walk;  // refmap: max bases to walk outward along assemblies, per flank
	int8_t walk_mode;  // refmap: RB3_WALK_*
	int64_t max_occ;   // refmap: reject reads/anchors occurring > max_occ times (0 = off; <0 = auto = #taxa)
	int64_t max_bracket; // refmap: reject a PLACED if |cR-cL| > max_bracket (0 = off)
	int8_t two_flank;  // refmap: require both flanks to anchor concordantly (1 = on)
	int8_t allow_walk; // refmap: opt in to the DEPRECATED flank-walking path (no --lift)
	char *lift_fn;     // refmap: liftover file -> project assembly hits instead of walking
	int64_t lift_win, lift_mad; // refmap: liftover projection window / max residual MAD (bp)
	int32_t kmer_len, kmer_step, min_agree; // refmap: k-mer-agreement placement (0 = off)
	int64_t kmer_cluster;                   // refmap: cluster tolerance for agreeing k-mers (bp)
	char *ps4g_fn;      // refmap: write a PS4G file here (NULL = off)
	char *ps4g_per_read_fn;// refmap: write a per-read PS4G file here (NULL = off)
	char *npy_fn;       // refmap: write a numpy (.npy) training/inference array here (NULL = off)
	char *label_bed_fn; // refmap: diploid training labels (chrom start end sampleA [sampleB]), NULL = off
	int64_t bin_size;   // refmap: PS4G/npy position bin size in bp (default 256)
	int64_t max_intron;    // chain: reject a chain link whose unexplained ref jump (dRef-dQuery) exceeds this
	int32_t gap_intron;    // chain: reference gap (bp) that starts a new exon segment
	int32_t chain_max_occ; // chain: SMEMs with FM-interval size > this are uninformative (skipped)
	char *ref_fasta;       // chain: reference (or pangenome) FASTA -> GT-AG splice-site check (NULL = off)
	int32_t splice_min;    // chain: reference gap size (bp) above which the GT-AG check applies
	int32_t pav_min_len;   // chain --lift: min LONGEST assembly-only SMEM to emit a pav: row
	int32_t trim_polya;    // trim a terminal poly-A/poly-T run of >= this many bp (0 = off)
	int32_t pav_grid;      // chain --lift: snap the emitted pav: position to this grid (0 = off)
	int8_t diverged_rows;  // chain --lift: diverged rows: 0 = ordinary row, 1 = pav: row, 2 = drop
	int8_t npy_binary;  // refmap: --npy writes presence (1) instead of read counts (0 = off, counts)
	int64_t target_hits; // refmap: stop reading once this many PLACED+EXACT records have been written (0 = off)
	int8_t report_occ; // refmap: append the raw FM-index interval size (occurrence count) as an
	                    // extra output column (0 = off, opt-in -- see --report-occ)
	rb3_swopt_t swo;
} rb3_mopt_t;

void rb3_mopt_init(rb3_mopt_t *opt)
{
	memset(opt, 0, sizeof(rb3_mopt_t));
	opt->n_threads = 4;
	opt->min_occ = 1;
	opt->min_len = 31;    // plant pangenomes: 19 is too short to be specific (Ed, 2026-07-26)
	opt->max_intron = 500, opt->gap_intron = 30, opt->chain_max_occ = -1; // chain: chain_max_occ<0 = auto (2*#samples, cap 256)
	opt->ref_fasta = 0, opt->splice_min = 10; // chain: GT-AG splice check off unless --ref-fasta given
	opt->pav_min_len = 60;  // chain --lift: assembly-only specificity floor; see chain_emit_pav
	opt->trim_polya = 0;    // off by default; see m_trim_polya
	opt->pav_grid = 5000;   // chain --lift: pav is presence/absence; see chain_emit_pav (E3)
	opt->diverged_rows = 0;     // colinear rows have a real reference coordinate -> ordinary row
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
	opt->ps4g_fn = 0, opt->ps4g_per_read_fn = 0, opt->npy_fn = 0, opt->label_bed_fn = 0;
	opt->bin_size = 256;
	opt->npy_binary = 0; // off by default (write read counts, not presence/absence)
	opt->target_hits = 0; // off by default (read the whole input)
	opt->report_occ = 0;  // off by default (extra output column, opt-in via --report-occ)
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
	char **ref_seq;  // chain --ref-fasta: ref_seq[k] = uppercase ACGT sequence of reference seq k (NULL if not loaded)
	rb3_lift_t *lift; // refmap: assembly->reference liftover (NULL = walk)
	rb3_gtab_t *gtab;      // refmap --ps4g/--npy: sample (gamete) table, NULL unless requested
	rb3_ps4g_acc_t *ps4g_acc; // refmap --ps4g/--npy: accumulated per-read support events
	FILE *ps4g_per_read_fp;   // refmap --ps4g-per-read: per-read PS4G file, NULL unless requested
	rb3_bed_t *label_bed;  // refmap --label-bed: diploid training labels, NULL unless requested
	rb3_hitcount_t hitcount; // refmap --target-hits: PLACED+EXACT records written so far vs. the target
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
	int32_t n_assembly;   // distinct assembly sequences seen (capped)
	int64_t ref_sid;     // reference sequence index (in [0,n_seq)); -1 if none
	int64_t cL, cR;      // reference coordinates bracketing the query; -1 if unknown
	int64_t ins_size;    // implied size inserted into the assembly relative to the reference
	int32_t n_asm_list;  // number of assemblies stored below
	rb3_pos_t *assemblies; // a few assembly (sid,pos); allocated with RB3_MALLOC
	int32_t n_sub;       // per-assembly mode: number of sub-results (one per assembly)
	struct refmap_rst_s *sub; // per-assembly mode: one placement per assembly
	int32_t n_vote, agree, second, mapq; // --kmer: informative tiles, agreeing k-mers, runner-up, calibrated MAPQ
	int32_t *gametes, n_gametes; // PS4G/npy: sorted, deduped gamete (sample) indices supporting this call; allocated with RB3_MALLOC
	int64_t occ; // --report-occ: raw FM-index interval size for the matched query (whole read, or the
	             // matched SMEM core if there was no end-to-end match) -- 0 for UNPLACED and for
	             // --kmer-mode results (no single whole-query interval exists there). Pangenome-wide:
	             // counts exact matches to the read's given orientation anywhere across the reference
	             // + all assemblies (both strands are indexed, so a genuine inverted/palindromic repeat
	             // is correctly counted more than once -- this is real signal, not an artifact to
	             // divide away; empirically, a single-copy locus present in only one assembly reports
	             // occ=1, not 2 -- verified directly against known unique loci, see --report-occ help).
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

// Length of a terminal homopolymer run of `base`, scanning inward from one end. Scored scan:
// +1 per matching base, -RB3_TRIM_MM per mismatch, cut at the highest-scoring prefix. This
// absorbs sequencing error inside the tail without letting the scan run past the tail boundary
// and latch onto incidental A/T bases in genuine sequence -- a plain mismatch budget over-trims
// real sequence by several bp (verified: a 40 bp poly-T head cut 45 bases).
#define RB3_TRIM_MM   3   // mismatch penalty
#define RB3_TRIM_DROP 12  // stop once the score falls this far below the best (≈4 mismatches)
static int32_t trim_run(const uint8_t *seq, int32_t len, int32_t step, uint8_t base, int32_t min_run)
{
	int32_t i, k, sc = 0, best = 0, bestk = 0;
	for (k = 0; k < len; ++k) {
		i = step > 0? k : len - 1 - k;
		sc += (seq[i] == base)? 1 : -RB3_TRIM_MM;
		if (sc > best) best = sc, bestk = k + 1;
		else if (sc < best - RB3_TRIM_DROP) break;
	}
	return bestk >= min_run? bestk : 0;
}

// Trim a homopolymer tail: poly-A at the 3' end, or poly-T at the 5' end (the same tail on a
// reverse-complemented read). 3'-tag protocols prime on polyA and misprime on internal A-runs
// >=4 bp, so many reads carry a homopolymer tail with only a short informative remainder. The
// PAV specificity floor discards those reads (3'-tag: 2,363 rows at 20.5% source recall, mean
// longest A/T run 13.6 bp vs 4.0 bp for the high-recall population). A read that is entirely
// homopolymer trims to length 0 and yields no SMEMs, which is the correct outcome.
//
// MEASURED: this is a NO-OP for the chain/PAV pipeline, which is why it is off by default.
// Trimming cannot lengthen the informative core -- SMEM search already isolates it, since a
// polyA tail does not extend a genomic exact match unless the reference also carries A's there.
// On those same 2,363 reads: mean longest SMEM 42.6 bp untrimmed vs 42.7 bp trimmed, and 56
// reads lost outright. Whole-pipeline effect at --trim-polya=8 (3'-tag): reference rows
// 80,521 -> 80,287, pav source recall 76.2% -> 76.3%, ref source recall 96.5% -> 96.8%.
// The harmful polyA case is entirely covered by pav_low_complexity() instead. Kept because it
// is correct and cheap, and may matter for other read types or consumers.
// See pav_e1_findings_2026-07-26.md.
static void m_trim_polya(m_seq_t *s, int32_t min_run)
{
	int32_t cut;
	if (s->len <= 0) return;
	cut = trim_run(s->seq, s->len, -1, 1, min_run);          // poly-A at the 3' end
	if (cut > 0) s->len -= cut;
	if (s->len <= 0) { s->len = 0; return; }
	cut = trim_run(s->seq, s->len, 1, 4, min_run);           // poly-T at the 5' end
	if (cut > 0) { memmove(s->seq, s->seq + cut, s->len - cut); s->len -= cut; }
}

static void worker_for_seq(void *data, long i, int tid)
{
	step_t *t = (step_t*)data;
	const pipeline_t *p = t->p;
	m_seq_t *s = &t->seq[i];
	m_tbuf_t *b = &t->buf[tid];
	if (rb3_dbg_flag & RB3_DBG_QNAME)
		fprintf(stderr, "Q\t%s\t%d\n", s->name, tid);
	rb3_char2nt6(s->len, s->seq);
	if (p->opt->trim_polya > 0) m_trim_polya(s, p->opt->trim_polya);
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
		if (p->opt->algo == RB3_SA_MEM_TG || p->opt->algo == RB3_SA_CHAIN)
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
 * walking outward through the assembly genomes to the breakpoints.
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

// Walk outward from interval Iq through the assemblies, recording one base per step.
// side 0 = left (backward) flank, side 1 = right (forward) flank. On return, buf holds the flank
// in forward (5'->3') orientation; the breakpoint-proximal end (adjacent to the query) is the 3'
// end for the left flank and the 5' end (index 0) for the right flank. Returns the flank length.
// walk_mode picks the path: consensus follows the most-supported base; strict stops at the first
// assembly disagreement; per-assembly (mask != NULL) follows the single assembly marked in mask.
static int32_t refmap_extract_flank(void *km, const pipeline_t *p, const rb3_sai_t *Iq, int side, int32_t max_walk, int walk_mode, const uint8_t *mask, uint8_t *buf)
{
	const rb3_fmi_t *f = &p->fmi;
	rb3_sai_t I = *Iq, ok[RB3_ASIZE];
	int32_t n = 0, is_back = (side == 0);
	while (n < max_walk) {
		int c, best = -1;
		int64_t bestsz = 0;
		rb3_fmd_extend(f, &I, ok, is_back);
		if (mask) { // per-assembly: follow the base whose child still contains this assembly
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
		if (best < 0 || bestsz == 0) break;                         // no assembly continues
		if (walk_mode == RB3_WALK_STRICT && bestsz != I.size) break; // assemblies disagree
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

// Given the query interval and (optionally) a assembly mask, extract both flanks, re-anchor them in
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
// from the query (it is contiguous with the query only in the assemblies), so a plain SMEM would be
// swallowed by the longer assembly match that crosses the breakpoint. Instead we grow the match
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

// Place a assembly-only query by PROJECTING its assembly hits to the reference via
// the liftover (the E4 "second SSA"), instead of walking outward. The majority
// reference sequence among the per-assembly projections wins; NULL projections
// (no confident collinear support) are dropped, and if none survive -> UNPLACED.
static void refmap_place_lift(void *km, const pipeline_t *p, const m_seq_t *s, refmap_rst_t *r)
{
	const rb3_fmi_t *f = &p->fmi;
	int64_t rsids[RB3_RM_MAX_ASSEMBLY], rposs[RB3_RM_MAX_ASSEMBLY], best_rsid = -1, med, v[RB3_RM_MAX_ASSEMBLY];
	int32_t n = 0, k, j, bestn, m;
	for (k = 0; k < r->n_asm_list && n < RB3_RM_MAX_ASSEMBLY; ++k) {
		int64_t clen, st, en, rsid, rpos;
		pos_stranded(f->sid, &r->assemblies[k], s->len, &clen, &st, &en);
		if (rb3_lift_project(p->lift, km, (int32_t)(r->assemblies[k].sid >> 1), st,
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

// One reference-coordinate vote from a k-mer. sx = the physical source sequence
// (assembly/reference) this occurrence came from, before any liftover projection
// remaps the locus to rsid/rpos -- kept so the winning cluster's assemblies can be
// recovered for the PS4G/npy gameteSet (see refmap_query_kmer()).
typedef struct { int64_t rsid, rpos; int32_t kmer; int32_t sx; } refmap_vote_t;

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
// exact, single-copy-per-taxon match): a reference hit votes directly; a assembly
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
		(*votes)[*nv].rsid = rsid, (*votes)[*nv].rpos = rpos, (*votes)[*nv].kmer = ki;
		(*votes)[*nv].sx = (int32_t)sx, (*nv)++;
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
	int64_t best_i = -1, best_j = -1; // [best_i,best_j) into votes[]: the winning cluster,
	                                  // for PS4G/npy gameteSet collection below
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
			second_support = best_support, best_support = support, best_rsid = rs, best_pos = votes[(i + j) >> 1].rpos,
			best_i = i, best_j = j;
		else if (support > second_support)
			second_support = support;
		i = j;
	}
	r->agree = (int32_t)best_support, r->second = (int32_t)second_support;
	if (best_support >= o->min_agree) {
		r->status = RB3_RM_PLACED, r->ref_sid = best_rsid, r->strand = 0;
		r->cL = r->cR = best_pos, r->ins_size = 0, r->n_assembly = (int32_t)best_support;
		r->mapq = refmap_kmer_mapq((int32_t)best_support, (int32_t)second_support);
		// PS4G/npy: assemblies backing any vote in the winning cluster. Looser than
		// whole-read mode's gameteSet (one single longest-match unit) -- this is a
		// union across the read's several independent k-mer tiles -- but built the
		// same way whole-read mode's EXACT path does: no dedup here, rb3_ps4g_acc_add
		// sorts+dedupes on ingestion.
		if (p->gtab && best_i >= 0) {
			int32_t m = (int32_t)(best_j - best_i), k;
			r->gametes = RB3_MALLOC(int32_t, m);
			for (k = 0; k < m; ++k) r->gametes[k] = p->gtab->sid2g[votes[best_i + k].sx];
			r->n_gametes = m;
		}
	}
	kfree(km, votes);
}

static void refmap_query(void *km, const pipeline_t *p, const m_seq_t *s, refmap_rst_t *r)
{
	const rb3_fmi_t *f = &p->fmi;
	rb3_sai_t Iq;
	if (p->opt->kmer_len > 0) { refmap_query_kmer(km, p, s, r); return; }
	rb3_pos_t *pos;
	int64_t np, i;
	int asm_seen[RB3_RM_MAX_ASSEMBLY], n_asm = 0;

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
	r->occ = Iq.size; // --report-occ: raw interval size, before any --max-occ/assembly-cap truncation
	                  // below, so MULTI/EXACT/PLACED rows all carry the true occurrence count

	// E2: an informative read maps at most once per taxon; a read occurring > max_occ times is a
	// repeat/retro and cannot be placed confidently (it would be reported at one arbitrary copy).
	if (p->opt->max_occ > 0 && Iq.size > p->opt->max_occ) {
		r->status = RB3_RM_MULTI;
		return;
	}

	// locate a sample of occurrences; separate reference hits (exact) from assemblies
	pos = Kmalloc(km, rb3_pos_t, RB3_RM_MAX_ASSEMBLY);
	np = rb3_ssa_multi(km, f, f->ssa, Iq.x[0], Iq.x[0] + Iq.size, RB3_RM_MAX_ASSEMBLY, pos);
	{
		int ref_found = 0;
		for (i = 0; i < np; ++i) {
			if (!ref_found && p->is_ref[pos[i].sid>>1]) { // the query is present in the reference: report it directly
				int64_t clen, st, en;
				pos_stranded(f->sid, &pos[i], s->len, &clen, &st, &en);
				r->status = RB3_RM_EXACT, r->ref_sid = pos[i].sid>>1, r->strand = pos[i].sid&1;
				r->cL = st, r->cR = en, r->ins_size = 0;
				ref_found = 1; // keep scanning: PS4G/npy need every sample sharing this exact sequence
			}
		}
		if (ref_found) {
			if (p->gtab) { // equivalent hits: every occurrence of an EXACT read names a candidate parent
				int32_t tmp[RB3_RM_MAX_ASSEMBLY], m;
				for (i = 0, m = 0; i < np; ++i) tmp[m++] = p->gtab->sid2g[pos[i].sid>>1];
				r->gametes = RB3_MALLOC(int32_t, m);
				memcpy(r->gametes, tmp, m * sizeof(int32_t)); // rb3_ps4g_acc_add sorts+dedupes on ingestion
				r->n_gametes = m;
			}
			kfree(km, pos);
			return;
		}
	}
	r->assemblies = RB3_CALLOC(rb3_pos_t, RB3_RM_MAX_ASSEMBLY);
	for (i = 0; i < np && n_asm < RB3_RM_MAX_ASSEMBLY; ++i) { // keep distinct assembly sequences for reporting
		int32_t j, dup = 0;
		for (j = 0; j < n_asm; ++j)
			if (asm_seen[j] == (int)(pos[i].sid>>1)) { dup = 1; break; }
		if (dup) continue;
		asm_seen[n_asm] = pos[i].sid>>1;
		r->assemblies[n_asm++] = pos[i];
	}
	r->n_asm_list = n_asm, r->n_assembly = n_asm;
	kfree(km, pos);
	if (p->gtab && n_asm > 0) { // PS4G/npy: assembly samples backing a (possible) PLACED call below
		int32_t k;
		r->gametes = RB3_MALLOC(int32_t, n_asm);
		for (k = 0; k < n_asm; ++k) r->gametes[k] = p->gtab->sid2g[r->assemblies[k].sid>>1];
		r->n_gametes = n_asm;
	}

	if (p->lift) { // E4: project assembly hits to the reference instead of walking
		refmap_place_lift(km, p, s, r);
		return;
	}

	if (p->opt->walk_mode == RB3_WALK_PERASSEMBLY && n_asm > 0) {
		// place each assembly separately, following that one assembly through divergences
		uint8_t *mask = RB3_CALLOC(uint8_t, f->sid->n_seq);
		int32_t k;
		r->sub = RB3_CALLOC(refmap_rst_t, n_asm);
		r->n_sub = n_asm;
		for (k = 0; k < n_asm; ++k) {
			refmap_rst_t *sub = &r->sub[k];
			int64_t csid = r->assemblies[k].sid >> 1;
			memset(sub, 0, sizeof(*sub));
			sub->status = RB3_RM_UNPLACED, sub->qlen = s->len, sub->ref_sid = -1, sub->cL = sub->cR = -1;
			sub->occ = r->occ; // same overall query interval; each sub just reports a different assembly
			sub->assemblies = &r->assemblies[k], sub->n_asm_list = 1, sub->n_assembly = 1; // borrowed pointer
			if (p->gtab) {
				sub->gametes = RB3_MALLOC(int32_t, 1);
				sub->gametes[0] = p->gtab->sid2g[csid];
				sub->n_gametes = 1;
			}
			mask[csid] = 1;
			refmap_place(km, p, s, &Iq, mask, sub);
			mask[csid] = 0;
		}
		free(mask);
		return;
	}
	// consensus / strict: a single placement from the shared assembly path
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

static void write_refmap1(kstring_t *out, const rb3_fmi_t *f, const m_seq_t *s, const refmap_rst_t *r, int kmer, int report_occ)
{
	static const char *status_str[5] = { "UNPLACED", "PLACED", "ONE_SIDE", "EXACT", "MULTI" };
	int32_t k;
	out->l = 0;
	write_name(out, s);
	rb3_sprintf_lite(out, "\t%d\t%s\t%d\t", r->qlen, status_str[(int)r->status], r->n_assembly);
	if (r->n_asm_list > 0) { // assembly sequences (name:strand)
		for (k = 0; k < r->n_asm_list; ++k) {
			int64_t sid = r->assemblies[k].sid;
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
	if (report_occ) // --report-occ: raw FM-index interval size (pangenome-wide occurrence count; see refmap_rst_t.occ)
		rb3_sprintf_lite(out, "\t%ld", (long)r->occ);
	rb3_sprintf_lite(out, "\n");
}

// PS4G/npy: fold an EXACT or PLACED read's (ref position, supporting gametes) into the accumulator.
// Other statuses (UNPLACED/ONE_SIDE/MULTI) contribute no confident reference position and are skipped.
static void refmap_rst_accumulate(rb3_ps4g_acc_t *acc, const refmap_rst_t *r)
{
	if (acc && r->n_gametes > 0 && r->ref_sid >= 0 && (r->status == RB3_RM_EXACT || r->status == RB3_RM_PLACED))
		rb3_ps4g_acc_add(acc, r->ref_sid, r->cL, r->gametes, r->n_gametes);
}

// --ps4g-per-read: the per-read PS4G *before* aggregation -- one line per EXACT/PLACED
// read with the exact gameteSet it contributes, so a scorer can attribute each
// read's evidence strictly (no lookup into the collapsed PS4G). refContig is the
// bare contig part of the reference name (same first-'_' split as ps4g.c).
static void write_ps4g_read(FILE *fp, const rb3_fmi_t *f, const m_seq_t *s, const refmap_rst_t *r)
{
	int32_t k;
	const char *nm, *us;
	if (fp == 0 || r->n_gametes <= 0 || r->ref_sid < 0) return;
	if (r->status != RB3_RM_EXACT && r->status != RB3_RM_PLACED) return;
	nm = f->sid->name[r->ref_sid];
	us = strchr(nm, '_');
	fprintf(fp, "%s\t%s\t%ld\t", s->name? s->name : "?", us? us + 1 : nm, (long)r->cL);
	for (k = 0; k < r->n_gametes; ++k)
		fprintf(fp, "%s%d", k? "," : "", r->gametes[k]);
	fputc('\n', fp);
}

// --target-hits: count a written record toward the target if it's PLACED or EXACT.
// t->p is const (step_t's writers are meant to be read-only); the hit counter is
// the one deliberately mutable exception, so cast it away at this single call site
// rather than loosen const-ness everywhere step_t.p is used.
static void refmap_count_hit(const pipeline_t *p, int8_t status)
{
	rb3_hitcount_add(&((pipeline_t*)p)->hitcount, status == RB3_RM_PLACED || status == RB3_RM_EXACT);
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
		if (r->n_sub > 0) { // per-assembly mode: one line per assembly (sub->assemblies borrows r->assemblies)
			for (k = 0; k < r->n_sub; ++k) {
				write_refmap1(&out, f, s, &r->sub[k], kmer, p->opt->report_occ);
				fputs(out.s, stdout);
				refmap_rst_accumulate(p->ps4g_acc, &r->sub[k]);
				write_ps4g_read(p->ps4g_per_read_fp, f, s, &r->sub[k]);
				refmap_count_hit(p, r->sub[k].status);
				free(r->sub[k].gametes);
			}
			free(r->sub);
		} else {
			write_refmap1(&out, f, s, r, kmer, p->opt->report_occ);
			fputs(out.s, stdout);
			refmap_rst_accumulate(p->ps4g_acc, r);
			write_ps4g_read(p->ps4g_per_read_fp, f, s, r);
			refmap_count_hit(p, r->status);
		}
		free(r->assemblies);
		free(r->gametes);
		free(s->seq);
		free(s->name);
	}
	free(out.s);
}

/********************************************************************
 * chain: unite a read's SMEMs by colinear chaining + strict set
 * intersection, emitting a per-read PS4G row per exon segment.
 * Native port of rnaseq_ps4g/chain/chain_prototype.py, including the optional
 * GT-AG splice check when --ref-fasta loads the reference contigs. Reuses the MEM
 * SMEM+locate machinery: each s->mem[i] has its query span and located
 * occurrences; gametes come from gtab->sid2g, reference occurrences from
 * is_ref (pos_stranded gives the forward coordinate + strand).
 ********************************************************************/
#define RB3_CHAIN_GAP_NUM 2   // gap penalty numerator: 0.02 * 100 (score scaled x100 to stay integer)
#define RB3_CHAIN_SLACK   6   // ±bp window searched around a SMEM boundary for the canonical splice motif

/* Load reference contig sequences for the GT-AG splice check. Reads every record
 * of a reference (or whole-pangenome) FASTA and keeps, keyed by index sequence id,
 * an uppercase copy of the sequence for each contig marked in is_ref[]. Non-reference
 * records (other assemblies) and names absent from the index are ignored. Returns the
 * number of reference contigs filled; p->ref_seq[k] stays NULL for any not found. */
static int64_t chain_load_ref(pipeline_t *p, const char *fn)
{
	rb3_seqio_t *fp;
	rb3_name2sid_t *h;
	int64_t k, n_loaded = 0;
	const char *name;
	char *ss;
	int64_t len;
	int absent;
	p->ref_seq = RB3_CALLOC(char*, p->fmi.sid->n_seq);
	h = rb3_name2sid_init();                       // reference contig name -> index sid
	for (k = 0; k < p->fmi.sid->n_seq; ++k)
		if (p->is_ref[k]) {
			khint_t itr = rb3_name2sid_put(h, p->fmi.sid->name[k], &absent);
			kh_val(h, itr) = (int32_t)k;
		}
	fp = rb3_seq_open(fn, 0);
	if (fp == 0) { rb3_name2sid_destroy(h); free(p->ref_seq); p->ref_seq = 0; return -1; }
	while ((ss = rb3_seq_read1(fp, &len, &name)) != 0) {
		khint_t itr;
		int32_t sid;
		char *cp;
		int64_t i;
		if (name == 0) continue;
		itr = rb3_name2sid_get(h, name);
		if (itr == kh_end(h)) continue;            // not a reference contig -> skip
		sid = kh_val(h, itr);
		if (p->ref_seq[sid]) continue;             // first record for this name wins
		cp = RB3_MALLOC(char, len + 1);
		for (i = 0; i < len; ++i) { char c = ss[i]; cp[i] = (c >= 'a' && c <= 'z')? c - 32 : c; }
		cp[len] = 0;
		p->ref_seq[sid] = cp;
		++n_loaded;
	}
	rb3_seq_close(fp);
	rb3_name2sid_destroy(h);
	return n_loaded;
}

/* Does the reference intron spanning forward coordinates [lo,hi) have canonical
 * splice motifs? GT..AG in transcription orientation reads GT..AG on the forward
 * strand for a '+' gene and CT..AC (reverse complement) for a '-' gene. A SMEM's
 * end drifts a few bp from the true splice site (microhomology), so search a ±slack
 * window on each boundary, mirroring canonical_splice() in chain_prototype.py. */
static int chain_canonical_splice(const char *refseq, int64_t reflen, int8_t strand, int64_t lo, int64_t hi)
{
	const int slack = RB3_CHAIN_SLACK;
	int64_t d, a, dlo, dhi;
	const char *low_motif = strand? "CT" : "GT", *high_motif = strand? "AC" : "AG";
	if (hi - lo < 4) return 1;                     // too small to be a real intron -> allow
	dlo = lo - slack < 0? 0 : lo - slack;
	dhi = lo + slack < reflen - 1? lo + slack : reflen - 1;
	for (d = dlo; d <= dhi; ++d) {
		if (d + 1 >= reflen || refseq[d] != low_motif[0] || refseq[d+1] != low_motif[1]) continue;
		int64_t alo = d + 4 > hi - slack? d + 4 : hi - slack;
		int64_t ahi = hi + slack < reflen? hi + slack : reflen;
		for (a = alo; a <= ahi; ++a)
			if (a >= 2 && refseq[a-2] == high_motif[0] && refseq[a-1] == high_motif[1]) return 1;
	}
	return 0;
}

typedef struct {
	int32_t qs, qe;       // query span
	int64_t rpos;         // reference forward position of the chosen occurrence
	int32_t rsid;         // reference sequence index (for the contig name)
	int8_t  rstrand;      // 0 forward, 1 reverse
	int32_t n_ref;        // number of reference occurrences (>1 => ambiguous locus)
	int32_t *gam, n_gam;  // sorted, deduped gamete indices over ALL occurrences
	int64_t sc;           // DP: (bases - gap_penalty) * 100
	int32_t anch, prev;   // DP: #anchors and backtrack pointer
} chain_sm_t;

static int chain_cmp_i32(const void *a, const void *b)
{ int32_t x = *(const int32_t*)a, y = *(const int32_t*)b; return x < y? -1 : x > y? 1 : 0; }

static int chain_cmp_sm(const void *a, const void *b) // sort candidates by (qs, rpos)
{
	const chain_sm_t *x = (const chain_sm_t*)a, *y = (const chain_sm_t*)b;
	if (x->qs != y->qs) return x->qs < y->qs? -1 : 1;
	return x->rpos < y->rpos? -1 : x->rpos > y->rpos? 1 : 0;
}

static int32_t chain_isect(int32_t *a, int32_t na, const int32_t *b, int32_t nb) // a := sorted(a) ∩ sorted(b)
{
	int32_t i = 0, j = 0, k = 0;
	while (i < na && j < nb) {
		if (a[i] == b[j]) a[k++] = a[i], ++i, ++j;
		else if (a[i] < b[j]) ++i;
		else ++j;
	}
	return k;
}

// Tight colinearity tolerance (bp) for the PAV path: a assembly position whose flanking
// reference anchors disagree by more than this is inside an insertion -> use the nearest
// reference anchor (breakpoint) rather than a projected-through-the-insertion coordinate.
#define RB3_CHAIN_PAV_MAD 64
// Projection window (bp) for the PAV path: the breakpoint is the nearest anchor, so a
// modest window suffices and keeps the O(anchors^2) projection cheap (vs lift_win 500kb).
#define RB3_CHAIN_PAV_WIN 50000
#define RB3_CHAIN_PAV_NPROJ 16  // occurrences of the seed SMEM to project (majority + median consensus)
#define RB3_CHAIN_PAV_NSEED 8   // seed occurrences of the longest SMEM tried as a cluster anchor
// Carrier-space bound (bp) on a cluster's unexplained gap: two assembly-only SMEMs are
// "the same locus" only if their assembly offsets track their query offsets to within this.
#define RB3_CHAIN_PAV_CLUSTER 100000

typedef struct {
	int32_t qs, qe;        // query span
	int32_t *gam, n_gam;   // sorted, deduped gamete indices over all occurrences
	int32_t mi;            // index into s->mem
} pav_sm_t;

// Forward-strand assembly position of occurrence `t` of a SMEM of length `len`.
static inline int64_t pav_fwd_pos(const rb3_fmi_t *f, const rb3_pos_t *t, int32_t len)
{
	int32_t sidx = t->sid>>1;
	return (t->sid&1)? f->sid->len[sidx] - (t->pos + len) : t->pos;
}

// Max share (%) one base may occupy in the anchor before it is called low-complexity.
#define RB3_CHAIN_PAV_MAXBASE 80

// Is the query span [qs,qe) dominated by a single base? LENGTH IS NOT SPECIFICITY for a
// homopolymer: a 76 bp pure poly-A read has a 76 bp exact match, clears any length floor,
// and lands wherever some assembly happens to carry a long enough A-run. Measured: 772 such
// rows (741 full-length, 31 3'-tag) at **0.0%** source recall -- every one named the wrong
// assembly. 3'-tag libraries generate these by design (polyA priming), so gate on
// composition, not length. See pav_e1_findings_2026-07-26.md.
static int pav_low_complexity(const m_seq_t *s, int32_t qs, int32_t qe)
{
	int32_t cnt[7], i, c, best = 0, n = qe - qs;
	if (n <= 0) return 0;
	for (i = 0; i < 7; ++i) cnt[i] = 0;
	for (i = qs; i < qe; ++i) { c = s->seq[i]; if (c < 0 || c > 6) c = 6; ++cnt[c]; }
	for (i = 0; i < 7; ++i) if (cnt[i] > best) best = cnt[i];
	return best * 100 >= RB3_CHAIN_PAV_MAXBASE * n;
}

// Carrier-only (PAV) fallback for reads with NO reference-hitting SMEM: place the read
// at the nearest B73 breakpoint via the liftover, emitting one row `pav:<contig>` with
// the intersected assembly gamete set. Requires --lift. Presence/absence only (no
// internal resolution) -- see design/pav-assembly-coordinates-scope.md.
//
// The gamete sets are intersected only over SMEMs that are COLINEAR IN CARRIER SPACE.
// A read originates from one assembly, so its assembly-only SMEMs must be colinear in that
// assembly's coordinates; intersecting across unrelated loci (a PAV fragment plus an
// unrelated repeat fragment) silently drops the true source assembly. Measured before this
// clustering existed: 33% source dropout on 3'-tag RNAseq, concentrated in confidently-
// wrong singletons -- see RopeBWTGrits-eval/experiments/pav_e1_findings_2026-07-26.md.
// This mirrors what chain_emit does in reference space (colinear DP, then intersect over
// the winning chain only).
static void chain_emit_pav(const pipeline_t *p, const m_seq_t *s, void *km)
{
	const rb3_fmi_t *f = &p->fmi;
	const rb3_gtab_t *gt = p->gtab;
	int32_t i, j, k, n = 0, na = 0, *acc = 0, sd = -1, sd_len = 0;
	int32_t *clu = 0, n_clu = 0, *cand = 0, best_o = -1;
	int64_t best_sc = -1;
	pav_sm_t *cs;
	cs = RB3_CALLOC(pav_sm_t, s->n_mem);
	for (i = 0; i < s->n_mem; ++i) { // collect assembly-only informative SMEMs
		m_sai_pos_t *r = &s->mem[i];
		int32_t st = r->mem.info>>32, en = (int32_t)r->mem.info, ng = 0, nref = 0, *g, m;
		if (r->n_pos == 0 || (int64_t)r->mem.size > p->opt->chain_max_occ) continue;
		g = RB3_MALLOC(int32_t, r->n_pos);
		for (k = 0; k < r->n_pos; ++k) {
			int32_t sidx = r->pos[k].sid>>1;
			g[ng++] = gt->sid2g[sidx];
			if (p->is_ref[sidx]) nref++;
		}
		if (nref > 0) { free(g); continue; } // has a reference hit -> not a PAV SMEM
		qsort(g, ng, sizeof(int32_t), chain_cmp_i32);
		for (k = 0, m = 0; k < ng; ++k) if (m == 0 || g[k] != g[m-1]) g[m++] = g[k];
		cs[n].qs = st, cs[n].qe = en, cs[n].gam = g, cs[n].n_gam = m, cs[n].mi = i;
		if (en - st > sd_len) sd_len = en - st, sd = n;
		++n;
	}
	if (n == 0) { free(cs); return; } // no assembly-only SMEM
	// Cluster in assembly space: seed on the longest SMEM and, for each of its assembly
	// occurrences, keep the SMEMs whose assembly offset tracks their query offset. The
	// highest-scoring cluster (bases covered) wins; its seed occurrence is the anchor.
	clu = RB3_MALLOC(int32_t, n);
	cand = RB3_MALLOC(int32_t, n);
	{
		m_sai_pos_t *rs = &s->mem[cs[sd].mi];
		for (k = 0; k < rs->n_pos && k < RB3_CHAIN_PAV_NSEED; ++k) {
			int32_t sid_o = rs->pos[k].sid>>1, str_o = rs->pos[k].sid&1, nc = 0;
			int64_t cp_o = pav_fwd_pos(f, &rs->pos[k], sd_len), sc = sd_len;
			cand[nc++] = sd;
			for (j = 0; j < n; ++j) {
				m_sai_pos_t *rj;
				int32_t len_j, u, hit = 0;
				if (j == sd) continue;
				rj = &s->mem[cs[j].mi];
				len_j = cs[j].qe - cs[j].qs;
				for (u = 0; u < rj->n_pos && !hit; ++u) {
					int64_t dq, dc;
					if ((rj->pos[u].sid>>1) != sid_o || (rj->pos[u].sid&1) != str_o) continue;
					dq = cs[j].qs - cs[sd].qs;
					dc = pav_fwd_pos(f, &rj->pos[u], len_j) - cp_o;
					if (str_o) dc = -dc;                     // reverse: query fwd -> assembly fwd back
					if ((dq > 0 && dc < 0) || (dq < 0 && dc > 0)) continue;      // out of order
					if (llabs(dc - dq) > RB3_CHAIN_PAV_CLUSTER) continue;        // implausible jump
					hit = 1;
				}
				if (hit) cand[nc++] = j, sc += len_j;
			}
			if (sc > best_sc) {
				best_sc = sc, best_o = k, n_clu = nc;
				memcpy(clu, cand, nc * sizeof(int32_t));
			}
		}
	}
	// Two suppression rules, both calibrated on B97 vs the NAM PHG (pav_e1b):
	// (1) INCOHERENT read: the winning cluster must contain EVERY assembly-only SMEM. If any
	//     is left out, the read's SMEMs cannot be reconciled to one assembly locus (chimera,
	//     adapter, mispriming) and the emitted set is near-worthless -- source recall by
	//     cluster/total was 1/1 68%, 2/2 94%, but 1/2 17%, 1/3 16%, 1/4 8%, 1/5 1.5%.
	//     Before clustering existed, the strict all-SMEM intersection collapsed these to
	//     na==0 and dropped them by accident; that accident was doing real work.
	// (2) SPECIFICITY FLOOR: a short assembly-only match is not specific enough to name the
	//     assemblies it came from. Source recall vs cluster bases (3'-tag, single-SMEM):
	//     31-39bp 20%, 40-49 37%, 50-59 56%, 60-79 71%, 80-119 78%. `min_len` (31) is the
	//     floor for anchoring to the REFERENCE, where a colinear chain adds confirmation;
	//     the assembly-only path has no such confirmation and needs its own, longer floor.
	//     The floor is on the LONGEST single match (sd_len), not the cluster total: summing
	//     several short SMEMs clears any total-bases floor while no individual match is
	//     specific (a 7-SMEM class of 728 full-length rows did exactly that, at 0% recall).
	// (3) LOW COMPLEXITY: reject a homopolymer-dominated anchor -- long but not specific.
	if (best_o >= 0 && n_clu == n && sd_len >= p->opt->pav_min_len
		&& !pav_low_complexity(s, cs[sd].qs, cs[sd].qe)) {
		acc = RB3_MALLOC(int32_t, cs[clu[0]].n_gam);
		na = cs[clu[0]].n_gam;
		memcpy(acc, cs[clu[0]].gam, na * sizeof(int32_t));
		for (i = 1; i < n_clu; ++i) na = chain_isect(acc, na, cs[clu[i]].gam, cs[clu[i]].n_gam);
	}
	// CONSENSUS projection. Previously the breakpoint came from ONE occurrence of the seed
	// (best_o), with a 30 bp agreement check over the first 4. refmap_place_lift has always
	// done better: project EVERY assembly, take the majority reference sequence, then the MEDIAN
	// position among that majority. That asymmetry let reads at one locus anchor via different
	// assemblies and jump: E3 found 16% of loci whose breakpoints spread >5x the region's true
	// extent (median true span 4.9 kb, median spread 40.6 kb). Consensus + a dispersion test is
	// both the fix and the per-read form of "suppress the scattered ones".
	if (na > 0) {
		m_sai_pos_t *r = &s->mem[cs[sd].mi];
		int64_t prs[RB3_CHAIN_PAV_NPROJ], prp[RB3_CHAIN_PAV_NPROJ], v[RB3_CHAIN_PAV_NPROJ];
		int pmd[RB3_CHAIN_PAV_NPROJ];
		int32_t np = 0, nmaj = 0, bestn = 0, n_insertion = 0;
		int64_t best_rsid = -1, med = 0, disp = 0;
		for (k = 0; k < r->n_pos && np < RB3_CHAIN_PAV_NPROJ; ++k) {
			int64_t rsid, rpos; int mode;
			if (!rb3_lift_project_bp(p->lift, km, (int32_t)(r->pos[k].sid>>1),
									 pav_fwd_pos(f, &r->pos[k], sd_len),
									 RB3_CHAIN_PAV_WIN, RB3_CHAIN_PAV_MAD, 4, &rsid, &rpos, &mode)) continue;
			prs[np] = rsid, prp[np] = rpos, pmd[np] = mode, ++np;
		}
		for (k = 0; k < np; ++k) { // majority reference sequence
			int32_t c = 0;
			for (j = 0; j < np; ++j) if (prs[j] == prs[k]) ++c;
			if (c > bestn) bestn = c, best_rsid = prs[k];
		}
		for (k = 0; k < np; ++k) if (prs[k] == best_rsid) { v[nmaj++] = prp[k]; if (pmd[k] == 1) ++n_insertion; }
		for (k = 1; k < nmaj; ++k) { int64_t x = v[k]; for (j = k - 1; j >= 0 && v[j] > x; --j) v[j+1] = v[j]; v[j+1] = x; }
		if (nmaj > 0) {
			med = v[nmaj >> 1];
			for (k = 0; k < nmaj; ++k) { int64_t d = llabs(v[k] - med); if (d > disp) disp = d; }
		}
		// SUPPRESS when the occurrences do not agree on WHERE this is: require a majority of
		// projections on one reference sequence, and their spread within one grid cell. A read
		// whose own assemblies disagree by more than the emitted resolution cannot be placed.
		if (nmaj > 0 && bestn == np && disp <= (p->opt->pav_grid > 0? p->opt->pav_grid : p->opt->gap_intron)) {
			const char *nm = f->sid->name[best_rsid], *us = strchr(nm, '_');
			const char *contig = us? us + 1 : nm;
			int is_insertion = n_insertion * 2 >= nmaj;   // majority of projections were true breakpoints
			// A DIVERGED row is not a PAV: the sequence is in the reference, just too divergent
			// to share a SMEM, and the projection is a real colinear coordinate rather than a
			// flanking breakpoint. It therefore deserves the ordinary schema and its exact
			// position (no grid snapping, which exists only because a breakpoint is approximate).
			// Caveat measured on B97: colinear rows match mode1 on genomic (100.0%) and
			// full-length (97.7% vs 98.4%) source recall, but are 9.4 pts worse on 3'-tag
			// (75.0% vs 84.4%) -- those reads still have no reference SMEM, so the coordinate is
			// lifted from a assembly and is one step less direct. Hence the policy switch.
			if (!is_insertion) {
				if (p->opt->diverged_rows == 2) goto pav_done;          // drop
				if (p->opt->diverged_rows == 0) {                        // ordinary row, exact position
					printf("%s\t%s\t%ld\t", s->name? s->name : "?", contig, (long)med);
					for (i = 0; i < na; ++i) printf("%s%d", i? "," : "", acc[i]);
					putchar('\n');
					goto pav_done;
				}
			}
			// Snap to --pav-grid. The design calls PAV "presence/absence only (no internal
			// resolution)", but emission never enforced that, so single-base coordinates
			// promised a precision the method does not have and fragmented one insertion's
			// evidence across many bins (mean 6.3). Snapping collapses it -- safe because E3
			// showed source recall is 100% in every spread bucket, i.e. the assembly sets are
			// right and only the aggregation key moved. 5 kb -> 72% of insertions in one bin.
			int64_t bp = p->opt->pav_grid > 0? med / p->opt->pav_grid * p->opt->pav_grid : med;
			printf("%s\tpav:%s\t%ld\t", s->name? s->name : "?", contig, (long)bp);
			for (i = 0; i < na; ++i) printf("%s%d", i? "," : "", acc[i]);
			// 5th column = row class. 1 = INSERTION: absent from the reference, placed at
			// the nearest reference breakpoint, approximate (presence/absence only).
			// 0 = DIVERGED: present in the reference but too divergent to share a SMEM,
			// so the coordinate is colinear and exact. Only insertions reach here --
			// diverged rows are routed by --diverged-rows (default: ordinary row).
			printf("\t%d\n", is_insertion? 1 : 0);
		}
	}
pav_done:
	for (i = 0; i < n; ++i) free(cs[i].gam);
	free(acc); free(clu); free(cand); free(cs);
}

static void chain_emit(const pipeline_t *p, const m_seq_t *s, void *km)
{
	const rb3_fmi_t *f = &p->fmi;
	const rb3_gtab_t *gt = p->gtab;
	int32_t i, j, n = 0, nc;
	int8_t strand;
	int64_t spanp = 0, spanm = 0;
	chain_sm_t *cs;
	int32_t *chain, best;
	if (s->n_mem == 0) return;
	cs = RB3_CALLOC(chain_sm_t, s->n_mem);
	for (i = 0; i < s->n_mem; ++i) { // build informative candidates that hit the reference
		m_sai_pos_t *r = &s->mem[i];
		int32_t st = r->mem.info>>32, en = (int32_t)r->mem.info, k, ng = 0, nref = 0, rsid = -1;
		int64_t rpos = -1; int8_t rstr = 0;
		int32_t *g;
		if (r->n_pos == 0) continue;
		g = RB3_MALLOC(int32_t, r->n_pos);
		for (k = 0; k < r->n_pos; ++k) {
			rb3_pos_t *t = &r->pos[k];
			int32_t sidx = t->sid>>1;
			g[ng++] = gt->sid2g[sidx];
			if (p->is_ref[sidx]) {
				int64_t rlen = f->sid->len[sidx];
				int64_t fp = (t->sid&1)? rlen - (t->pos + (en - st)) : t->pos;
				if (nref == 0) rpos = fp, rstr = t->sid & 1, rsid = sidx;
				++nref;
			}
		}
		if (nref == 0 || (int64_t)r->mem.size > p->opt->chain_max_occ) { free(g); continue; } // uninformative or no ref hit
		qsort(g, ng, sizeof(int32_t), chain_cmp_i32);
		{ int32_t m = 0; for (k = 0; k < ng; ++k) if (m == 0 || g[k] != g[m-1]) g[m++] = g[k]; ng = m; }
		cs[n].qs = st, cs[n].qe = en, cs[n].rpos = rpos, cs[n].rsid = rsid, cs[n].rstrand = rstr;
		cs[n].n_ref = nref, cs[n].gam = g, cs[n].n_gam = ng;
		++n;
	}
	if (n == 0) { // no ref anchor -> try the PAV breakpoint path (which reads s->mem[].pos)
		free(cs);
		if (p->lift) chain_emit_pav(p, s, km);
		for (i = 0; i < s->n_mem; ++i) free(s->mem[i].pos);
		return;
	}
	for (i = 0; i < n; ++i) (cs[i].rstrand? &spanm : &spanp)[0] += cs[i].qe - cs[i].qs;
	strand = spanp >= spanm? 0 : 1;                       // dominant strand (ties -> '+')
	{ int32_t m = 0; for (i = 0; i < n; ++i) { if (cs[i].rstrand == strand) cs[m++] = cs[i]; else free(cs[i].gam); } n = m; }
	if (n == 0) { free(cs); for (i = 0; i < s->n_mem; ++i) free(s->mem[i].pos); return; }
	qsort(cs, n, sizeof(chain_sm_t), chain_cmp_sm);
	for (i = 0; i < n; ++i) cs[i].anch = 1, cs[i].sc = (int64_t)(cs[i].qe - cs[i].qs) * 100, cs[i].prev = -1;
	for (i = 0; i < n; ++i) { // colinear in-order DP: maximize (#anchors, bases - gap_penalty)
		int64_t ri = cs[i].rpos;
		for (j = 0; j < i; ++j) {
			int64_t rj = cs[j].rpos, dr = strand? rj - ri : ri - rj;
			int64_t dq = cs[i].qs - cs[j].qe, unexp;
			int32_t anch; int64_t sc;
			if (dr < 0 || cs[j].qs > cs[i].qs) continue;  // out of order -> not colinear
			if (dq < 0) dq = 0;
			unexp = dr - dq; if (unexp < 0) unexp = 0;    // intron length / paralog jump
			if (unexp > p->opt->max_intron) continue;     // implausible -> reject the link
			if (p->ref_seq && cs[i].rsid == cs[j].rsid && p->ref_seq[cs[j].rsid]) { // GT-AG splice-site check
				int32_t wj = cs[j].qe - cs[j].qs, wi = cs[i].qe - cs[i].qs;
				int64_t lo = strand? ri + wi : rj + wj, hi = strand? rj : ri; // forward intron bounds
				if (hi - lo >= p->opt->splice_min &&
					!chain_canonical_splice(p->ref_seq[cs[j].rsid], f->sid->len[cs[j].rsid], strand, lo, hi))
					continue;                             // non-canonical intron gap -> reject the link
			}
			anch = cs[j].anch + 1;
			sc = cs[j].sc + (int64_t)(cs[i].qe - cs[i].qs) * 100 - RB3_CHAIN_GAP_NUM * unexp;
			if (anch > cs[i].anch || (anch == cs[i].anch && sc > cs[i].sc))
				cs[i].anch = anch, cs[i].sc = sc, cs[i].prev = j;
		}
	}
	best = 0;
	for (i = 1; i < n; ++i)
		if (cs[i].anch > cs[best].anch || (cs[i].anch == cs[best].anch && cs[i].sc > cs[best].sc)) best = i;
	chain = RB3_MALLOC(int32_t, n);
	nc = 0;
	for (i = best; i >= 0; i = cs[i].prev) chain[nc++] = i;
	for (i = 0; i < nc/2; ++i) { int32_t tmp = chain[i]; chain[i] = chain[nc-1-i], chain[nc-1-i] = tmp; }
	{
		int ambiguous = 0;
		for (i = 0; i < nc; ++i) if (cs[chain[i]].n_ref > 1) ambiguous = 1; // maps to >1 locus -> suppress
		if (!ambiguous) {
			int32_t *acc = RB3_MALLOC(int32_t, cs[chain[0]].n_gam), na = cs[chain[0]].n_gam;
			int64_t *ivlo = RB3_MALLOC(int64_t, nc), *ivhi = RB3_MALLOC(int64_t, nc);
			memcpy(acc, cs[chain[0]].gam, na * sizeof(int32_t));
			for (i = 1; i < nc; ++i) na = chain_isect(acc, na, cs[chain[i]].gam, cs[chain[i]].n_gam);
			if (na > 0) { // segment forward intervals by intron-sized gaps, emit one row per segment
				const char *nm = f->sid->name[cs[chain[0]].rsid], *us = strchr(nm, '_');
				const char *contig = us? us + 1 : nm;
				int64_t lo, hi;
				for (i = 0; i < nc; ++i) ivlo[i] = cs[chain[i]].rpos, ivhi[i] = cs[chain[i]].rpos + (cs[chain[i]].qe - cs[chain[i]].qs);
				for (i = 0; i < nc - 1; ++i) // insertion sort ascending by start (nc is tiny)
					for (j = i + 1; j < nc; ++j)
						if (ivlo[j] < ivlo[i]) { int64_t t0 = ivlo[i], t1 = ivhi[i]; ivlo[i] = ivlo[j], ivhi[i] = ivhi[j], ivlo[j] = t0, ivhi[j] = t1; }
				lo = ivlo[0], hi = ivhi[0];
				for (i = 1; i < nc; ++i) {
					if (ivlo[i] - hi > p->opt->gap_intron) { // intron gap -> flush this exon segment
						printf("%s\t%s\t%ld\t", s->name? s->name : "?", contig, (long)lo);
						for (j = 0; j < na; ++j) printf("%s%d", j? "," : "", acc[j]);
						putchar('\n');
						lo = ivlo[i], hi = ivhi[i];
					} else if (ivhi[i] > hi) hi = ivhi[i];
				}
				printf("%s\t%s\t%ld\t", s->name? s->name : "?", contig, (long)lo);
				for (j = 0; j < na; ++j) printf("%s%d", j? "," : "", acc[j]);
				putchar('\n');
			}
			free(ivlo); free(ivhi);
			free(acc);
		}
	}
	for (i = 0; i < n; ++i) free(cs[i].gam);
	for (i = 0; i < s->n_mem; ++i) free(s->mem[i].pos);
	free(chain); free(cs);
}

static void write_per_seq(step_t *t)
{
	const pipeline_t *p = t->p;
	int32_t i, j;
	kstring_t out = {0,0,0};
	void *km = p->opt->algo == RB3_SA_CHAIN? km_init() : 0; // arena reused by the PAV lift projections
	for (j = 0; j < t->n_seq; ++j) {
		m_seq_t *s = &t->seq[j];
		free(s->seq);
		out.l = 0;
		if (p->opt->algo == RB3_SA_CHAIN) { // unite SMEMs -> per-read PS4G (per exon segment)
			chain_emit(p, s, km);
		} else if (p->opt->algo == RB3_SA_SW && (p->opt->flag & RB3_MF_WRITE_ALL)) { // write all hits in a compact format
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
	if (km) km_destroy(km);
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
		while (!rb3_hitcount_reached(&p->hitcount) && (ss = rb3_seq_read1(p->fp, &len, &name)) != 0) { // read sequences
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
	{ "ps4g",            ko_required_argument, 320 },
	{ "ps4g-per-read",      ko_required_argument, 326 },
	{ "npy",             ko_required_argument, 321 },
	{ "label-bed",       ko_required_argument, 322 },
	{ "bin-size",        ko_required_argument, 323 },
	{ "npy-binary",      ko_no_argument,       324 },
	{ "target-hits",     ko_required_argument, 325 },
	{ "report-occ",      ko_no_argument,       326 },
	{ "max-intron",      ko_required_argument, 330 },
	{ "gap-intron",      ko_required_argument, 331 },
	{ "chain-max-occ",   ko_required_argument, 332 },
	{ "ref-fasta",       ko_required_argument, 333 },
	{ "splice-min",      ko_required_argument, 334 },
	{ "pav-min-len",     ko_required_argument, 335 },
	{ "walk",            ko_no_argument,       336 },
	{ "trim-polya",      ko_required_argument, 337 },
	{ "pav-grid",        ko_required_argument, 338 },
	{ "diverged-rows",       ko_required_argument, 339 },
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
			else if (strcmp(o.arg, "per-assembly") == 0 || strcmp(o.arg, "per-carrier") == 0) opt.walk_mode = RB3_WALK_PERASSEMBLY;
			else { fprintf(stderr, "ERROR: --walk-mode must be consensus, strict or per-assembly\n"); return 1; }
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
		else if (c == 320) opt.ps4g_fn = o.arg;      // PS4G output (parents/gametes supporting each ref position)
		else if (c == 326) opt.ps4g_per_read_fn = o.arg; // per-read PS4G file (exact per-read gameteSet)
		else if (c == 321) opt.npy_fn = o.arg;       // numpy training/inference array
		else if (c == 322) opt.label_bed_fn = o.arg; // diploid training labels: chrom start end sampleA [sampleB]
		else if (c == 323) opt.bin_size = rb3_parse_num(o.arg); // PS4G/npy position bin size in bp
		else if (c == 324) opt.npy_binary = 1; // npy: write presence (1) instead of read counts
		else if (c == 325) opt.target_hits = rb3_parse_num(o.arg); // stop once this many PLACED+EXACT records are written
		else if (c == 326) opt.report_occ = 1; // append raw FM-index interval size (occurrence count) as an extra column
		else if (c == 330) opt.max_intron = rb3_parse_num(o.arg);    // chain: max unexplained ref jump per link
		else if (c == 331) opt.gap_intron = atoi(o.arg);             // chain: ref gap that starts a new exon segment
		else if (c == 332) opt.chain_max_occ = atoi(o.arg);          // chain: interval-size cap for an informative SMEM
		else if (c == 333) opt.ref_fasta = o.arg;                    // chain: reference FASTA -> GT-AG splice check
		else if (c == 334) opt.splice_min = atoi(o.arg);             // chain: min ref gap for the GT-AG check
		else if (c == 335) opt.pav_min_len = atoi(o.arg);            // chain --lift: assembly-only specificity floor
		else if (c == 336) opt.allow_walk = 1;                       // refmap: opt in to deprecated walking
		else if (c == 337) opt.trim_polya = atoi(o.arg);             // trim terminal poly-A/T runs >= INT bp
		else if (c == 338) opt.pav_grid = rb3_parse_num(o.arg);      // chain --lift: pav position grid
		else if (c == 339) {                                         // chain --lift: colinear-row policy
			if (strcmp(o.arg, "ordinary") == 0) opt.diverged_rows = 0;
			else if (strcmp(o.arg, "pav") == 0) opt.diverged_rows = 1;
			else if (strcmp(o.arg, "drop") == 0) opt.diverged_rows = 2;
			else { fprintf(stderr, "ERROR: --diverged-rows must be ordinary|pav|drop\n"); return 1; }
		}
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
	} else if (strcmp(argv[0], "chain") == 0) {
		opt.algo = RB3_SA_CHAIN;
		load_flag |= RB3_LOAD_ALL;
		if (opt.max_pos <= 0) opt.max_pos = opt.swo.max_pos = 64; // locate all occurrences of informative SMEMs
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
			fprintf(stderr, "  --lift=FILE       REQUIRED (standard): project assembly hits via a `ropebwt3 lift` map\n");
			fprintf(stderr, "  --walk            opt in to the DEPRECATED flank-walking path instead of --lift\n");
			fprintf(stderr, "  --max-walk=NUM    DEPRECATED (--walk only) max bases to walk per flank [%d]\n", opt.max_walk);
			fprintf(stderr, "  --walk-mode=STR   DEPRECATED (--walk only) consensus|strict|per-assembly [consensus]\n");
			fprintf(stderr, "  --max-occ=INT     drop reads/anchors occurring >INT times; <0 = auto (#taxa); 0 = off [%ld]\n", (long)opt.max_occ);
			fprintf(stderr, "  --two-flank       require both flanks to anchor concordantly (drop ONE_SIDE)\n");
			fprintf(stderr, "  --max-bracket=NUM with --two-flank, max |cR-cL| for a PLACED; 0 = off [%ld]\n", (long)opt.max_bracket);
			fprintf(stderr, "  --lift-win=NUM    liftover projection window [%ld]\n", (long)opt.lift_win);
			fprintf(stderr, "  --lift-mad=NUM    liftover max residual MAD [%ld]\n", (long)opt.lift_mad);
			fprintf(stderr, "  --kmer=INT        place a read from INT-bp k-mers by agreement (0 = off, whole-read)\n");
			fprintf(stderr, "  --kmer-step=INT   k-mer tiling step [%d]\n", opt.kmer_step);
			fprintf(stderr, "  --min-agree=INT   min agreeing k-mers to place [%d]\n", opt.min_agree);
			fprintf(stderr, "  --kmer-cluster=NUM  agreeing k-mers must fall within NUM bp [%ld]\n", (long)opt.kmer_cluster);
			fprintf(stderr, "  -l INT      min anchor length when re-mapping a flank [%ld]\n", (long)opt.min_len);
			fprintf(stderr, "  --ps4g=FILE       write PS4G v2.0 gamete-support counts (EXACT+PLACED reads)\n");
			fprintf(stderr, "  --ps4g-per-read=FILE write a per-read PS4G file (exact gameteSet per read)\n");
			fprintf(stderr, "  --npy=FILE        write a dense (bin x gamete+2) numpy training/inference array\n");
			fprintf(stderr, "  --label-bed=FILE  diploid training labels: chrom start end sampleA [sampleB]\n");
			fprintf(stderr, "  --bin-size=NUM    PS4G/npy reference position bin size in bp [%ld]\n", (long)opt.bin_size);
			fprintf(stderr, "  --npy-binary      npy: write presence (1) instead of read counts\n");
			fprintf(stderr, "  --target-hits=NUM stop once NUM PLACED/EXACT records are written (0 = off, read everything) [%ld]\n", (long)opt.target_hits);
			fprintf(stderr, "  --report-occ      append the raw FM-index interval size (occurrence count) as an extra\n");
			fprintf(stderr, "                    trailing column (0 = off, opt-in; pangenome-wide -- counts exact matches\n");
			fprintf(stderr, "                    to the reference + all assemblies together, not per-genome; 0 for --kmer mode)\n");
		}
		if (strcmp(argv[0], "chain") == 0) {
			fprintf(stderr, "  --ref-prefix=STR  reference = sequences whose name starts with STR [required]\n");
			fprintf(stderr, "  --max-intron=NUM  reject a chain link whose unexplained ref jump exceeds this [%ld]\n", (long)opt.max_intron);
			fprintf(stderr, "  --gap-intron=INT  reference gap starting a new exon segment [%d]\n", opt.gap_intron);
			fprintf(stderr, "  --chain-max-occ=INT  skip SMEMs whose FM interval exceeds this [auto: min(2*#samples,256)]\n");
			fprintf(stderr, "  --ref-fasta=FILE  reference/pangenome FASTA -> GT-AG splice check on intron-gap links\n");
			fprintf(stderr, "  --splice-min=INT  ref gap size above which the GT-AG check applies [%d]\n", opt.splice_min);
			fprintf(stderr, "  --lift=FILE       `ropebwt3 lift` map; enables PAV breakpoint anchoring: reads with\n");
			fprintf(stderr, "                    no reference SMEM emit one `pav:<contig>` row at the nearest\n");
			fprintf(stderr, "                    reference breakpoint, plus a 5th column (1 = true insertion,\n");
			fprintf(stderr, "                    0 = colinear = a divergent allele, not a PAV)\n");
			fprintf(stderr, "  --diverged-rows=STR  DIVERGED rows -- sequence that IS in the reference but too\n");
			fprintf(stderr, "                    divergent to share a SMEM, so the coordinate is colinear and exact\n");
			fprintf(stderr, "                    rather than a breakpoint: ordinary|pav|drop\n");
			fprintf(stderr, "                    [ordinary = emit as a normal row at its exact position]\n");
			fprintf(stderr, "  --pav-grid=NUM    snap the emitted pav: position to this grid, and require the\n");
			fprintf(stderr, "                    seed's assembly projections to agree within it (0 = off) [%d]\n", opt.pav_grid);
			fprintf(stderr, "  --pav-min-len=INT min LONGEST assembly-only SMEM to emit a pav: row; the\n");
			fprintf(stderr, "                    assembly-only path has no colinear reference confirmation, so it\n");
			fprintf(stderr, "                    needs a longer floor than -l [%d]\n", opt.pav_min_len);
			fprintf(stderr, "  --trim-polya=INT  trim a terminal poly-A (3') / poly-T (5') run of >=INT bp before\n");
			fprintf(stderr, "                    searching. MEASURED AS A NO-OP here (SMEMs already isolate the\n");
			fprintf(stderr, "                    genomic core); low-complexity anchors are gated instead [%d]\n", opt.trim_polya);
			fprintf(stderr, "  -l INT      min SMEM length [%ld]\n", (long)opt.min_len);
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
	p.is_ref = 0, p.n_ref = 0, p.lift = 0, p.ref_seq = 0;
	p.gtab = 0, p.ps4g_acc = 0, p.ps4g_per_read_fp = 0, p.label_bed = 0;
	rb3_hitcount_init(&p.hitcount, opt.target_hits);
	if (opt.algo == RB3_SA_REFMAP || opt.algo == RB3_SA_CHAIN) { // mark the reference sequences by name prefix
		int64_t k, plen;
		if (opt.ref_prefix == 0) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: %s requires --ref-prefix\n", argv[0]);
			return 1;
		}
		if (p.fmi.ssa == 0 || p.fmi.sid == 0) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: refmap needs the sampled suffix array (.ssa) and sequence names (.len.gz)\n");
			return 1;
		}
		// refmap: --lift is the standard resolution path. Flank walking is deprecated (it did
		// not work well) and used to be the SILENT default when --lift was omitted, so an
		// invocation that simply forgot --lift got the deprecated path and plausible-looking
		// output. Require an explicit choice instead of defaulting to the bad one.
		if (opt.algo == RB3_SA_REFMAP && opt.lift_fn == 0 && !opt.allow_walk) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: refmap needs --lift=FILE (the standard path; build it with `ropebwt3 lift`).\n"
								"       Flank walking is DEPRECATED and no longer the default; pass --walk to opt in.\n");
			return 1;
		}
		if (opt.algo == RB3_SA_REFMAP && opt.allow_walk && opt.lift_fn == 0 && rb3_verbose >= 1)
			fprintf(stderr, "WARNING: --walk selects the deprecated flank-walking path; prefer --lift=FILE.\n");
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
		{ // N = distinct taxa (name prefixes before '_'); drives sample-relative caps
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
			if (opt.max_occ < 0) { // refmap: auto occurrence cap = N
				opt.max_occ = n_taxa;
				if (rb3_verbose >= 3)
					fprintf(stderr, "[M::%s] auto --max-occ = %ld (distinct taxa)\n", __func__, (long)opt.max_occ);
			}
			if (opt.algo == RB3_SA_CHAIN) { // sample-relative caps: an FM interval counts BOTH strands, so
				// sequence present once per taxon has size ~2N; admit that (cap 256 to stay byte-sized).
				if (opt.chain_max_occ < 0) {
					int64_t c = 2 * n_taxa;
					opt.chain_max_occ = c < 256? (int32_t)c : 256;
				}
				if (opt.max_pos < opt.chain_max_occ) // locate all occurrences the filter admits
					opt.max_pos = opt.swo.max_pos = opt.chain_max_occ;
				if (rb3_verbose >= 3)
					fprintf(stderr, "[M::%s] chain: N=%ld -> --chain-max-occ=%d, max_pos=%d\n", __func__, (long)n_taxa, opt.chain_max_occ, opt.max_pos);
			}
		}
		if (opt.lift_fn) { // E4: load the assembly->reference liftover
			p.lift = rb3_lift_restore(opt.lift_fn);
			if (p.lift == 0) {
				if (rb3_verbose >= 1) fprintf(stderr, "ERROR: failed to load liftover '%s'\n", opt.lift_fn);
				free(p.is_ref);
				return 1;
			}
			if (rb3_verbose >= 3)
				fprintf(stderr, "[M::%s] loaded liftover over %ld sequences\n", __func__, (long)rb3_lift_n_seq(p.lift));
		}
		if (opt.algo == RB3_SA_CHAIN || opt.ps4g_fn || opt.npy_fn || opt.ps4g_per_read_fn) {
			p.gtab = rb3_gtab_build(p.fmi.sid); // gamete indices for PS4G, npy, the per-read PS4G file, and chain
			if (opt.algo == RB3_SA_CHAIN && opt.ref_fasta) { // load ref contigs for the GT-AG splice check
				int64_t n_load = chain_load_ref(&p, opt.ref_fasta);
				if (n_load < 0) {
					if (rb3_verbose >= 1) fprintf(stderr, "ERROR: failed to open --ref-fasta '%s'\n", opt.ref_fasta);
					free(p.is_ref);
					return 1;
				}
				if (rb3_verbose >= 3)
					fprintf(stderr, "[M::%s] loaded %ld of %ld reference contigs for the GT-AG splice check\n", __func__, (long)n_load, (long)p.n_ref);
			}
			if (opt.ps4g_fn || opt.npy_fn) {
				p.ps4g_acc = rb3_ps4g_acc_init(opt.bin_size);
				if (opt.label_bed_fn) {
					p.label_bed = rb3_bed_read(opt.label_bed_fn, p.gtab, p.fmi.sid, opt.ref_prefix);
					if (p.label_bed == 0) {
						if (rb3_verbose >= 1) fprintf(stderr, "ERROR: failed to read --label-bed '%s'\n", opt.label_bed_fn);
						free(p.is_ref);
						return 1;
					}
				}
			}
			if (opt.ps4g_per_read_fn) {
				p.ps4g_per_read_fp = fopen(opt.ps4g_per_read_fn, "w");
				if (p.ps4g_per_read_fp == 0) {
					if (rb3_verbose >= 1) fprintf(stderr, "ERROR: failed to write --ps4g-per-read '%s'\n", opt.ps4g_per_read_fn);
					free(p.is_ref);
					return 1;
				}
				fprintf(p.ps4g_per_read_fp, "readName\trefContig\trefPos\tgameteSet\n");
			}
		}
	} else if (opt.ps4g_fn || opt.npy_fn || opt.ps4g_per_read_fn || opt.label_bed_fn || opt.target_hits > 0) {
		if (rb3_verbose >= 1) fprintf(stderr, "ERROR: --ps4g/--npy/--label-bed/--target-hits only apply to refmap\n");
		return 1;
	}
	if (opt.algo == RB3_SA_CHAIN) puts("readName\trefContig\trefPos\tgameteSet");
	if (opt.flag & RB3_MF_WRITE_ALL) {
		puts("CC\tQS  queryName  queryLen  numHap");
		puts("CC\tQH  refCount   score     editDist   cs   strand   nOut   totAln");
		puts("CC");
	}
	for (j = o.ind + 1; j < argc; ++j) {
		if (rb3_hitcount_reached(&p.hitcount)) break; // --target-hits already satisfied; skip remaining input files
		p.fp = rb3_seq_open(argv[j], is_line);
		if (p.fp == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to load the sequence file '%s'\n", argv[j]);
			break;
		}
		kt_pipeline(2, worker_pipeline, &p, 3);
		rb3_seq_close(p.fp);
	}
	{
		int64_t tgt, n;
		if (rb3_hitcount_short(&p.hitcount, &tgt, &n) && rb3_verbose >= 1)
			fprintf(stderr, "WARNING: --target-hits=%ld requested but the input was exhausted after only %ld PLACED/EXACT records\n", (long)tgt, (long)n);
	}
	if (p.ps4g_acc) {
		kstring_t cmd = {0,0,0};
		int32_t k;
		rb3_sprintf_lite(&cmd, "ropebwt3");
		for (k = 0; k < argc; ++k) rb3_sprintf_lite(&cmd, " %s", argv[k]);
		rb3_ps4g_npy_finalize(p.ps4g_acc, p.gtab, p.fmi.sid, p.label_bed, opt.npy_binary, opt.ps4g_fn, opt.npy_fn, cmd.s);
		free(cmd.s);
		rb3_ps4g_acc_destroy(p.ps4g_acc);
	}
	if (p.ps4g_per_read_fp) fclose(p.ps4g_per_read_fp);
	rb3_bed_destroy(p.label_bed);
	rb3_gtab_destroy(p.gtab);
	if (p.ref_seq) {
		int64_t k;
		for (k = 0; k < p.fmi.sid->n_seq; ++k) free(p.ref_seq[k]);
		free(p.ref_seq);
	}
	rb3_fmi_free(&p.fmi);
	free(p.is_ref);
	rb3_lift_destroy(p.lift);
	return 0;
}
