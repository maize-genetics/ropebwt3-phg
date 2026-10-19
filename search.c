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
#include "pairing.h"
#include "errprof.h"
#include "khashl-km.h"

KHASHL_MAP_INIT(KH_LOCAL, rb3_name2sid_t, rb3_name2sid, kh_cstr_t, int32_t, kh_hash_str, kh_eq_str)

static int pos_cmp(const void *a, const void *b) // SMEM hits by (sequence/strand id, position)
{
	const rb3_pos_t *x = (const rb3_pos_t*)a, *y = (const rb3_pos_t*)b;
	if (x->sid != y->sid) return x->sid < y->sid? -1 : 1;
	return x->pos < y->pos? -1 : x->pos > y->pos? 1 : 0;
}

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
#define RB3_RM_MAX_ASSEMBLY 256 // refmap: default --max-locate (occurrences located per read)

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
	int32_t pav_agree;     // chain --lift: max spread among an assembly's projections to emit
	int8_t diverged_rows;  // chain --lift: diverged rows: 0 = ordinary row, 1 = pav: row, 2 = drop
	int8_t insertion_rows; // chain --lift: insertion rows: 0 = ordinary row (exact median, no snap),
	                       // 1 = pav: row (grid-snapped, DEFAULT -- unchanged existing behavior), 2 = drop
	int8_t npy_binary;  // refmap: --npy writes presence (1) instead of read counts (0 = off, counts)
	int64_t target_hits; // refmap: stop reading once this many PLACED+EXACT records have been written (0 = off)
	int8_t report_occ; // refmap: append the raw FM-index interval size (occurrence count) as an
	                    // extra output column (0 = off, opt-in -- see --report-occ)
	int8_t anchor_dist_npy; // refmap/chain --npy: also write a per-founder ternary read-sharing block
	                    // (match/diverged/deletion) and a per-founder distance-to-nearest-lift-anchor
	                    // block, widening --npy from (bin x gamete+2) to (bin x 3*gamete+2)
	                    // (0 = off, opt-in -- see --anchor-dist-npy; requires --lift)
	int64_t anchor_dist_thresh; // refmap/chain --anchor-dist-npy: distance (bp) at/under which a
	                    // non-matching founder reads as diverged rather than deletion
	int8_t support_best;   // chain: founder set = assemblies best explaining the read (1), or the
	                       // chain-wide intersection of reference-hitting SMEMs (0, default)
	int32_t support_tol;   // chain --support=best: keep assemblies within this many covered bases of the best
	int64_t support_win;   // chain --support=best --lift: an assembly-only SMEM counts only if it
	                       // projects within this many bp of the chain's reference span (0 = no check)
	int8_t dup_count;      // chain: a read named "<name>;x=N" stands for N identical reads (PS4G/npy counts)
	int8_t paired;         // --paired: reads are mate pairs (two files R1/R2, or one interleaved file); 0 = single-end
	double pe_mean, pe_sd; // --pe-isize=MEAN,SD: fixed insert-size distribution (pe_mean <= 0: estimate per batch)
	int64_t pe_slack;      // --pe-slack: bp added to the max insert size before two placed mates count as discordant
	int32_t max_locate;    // refmap: occurrences of a read located per read (the rest are not sampled) [256]
	int8_t pe_slack_set;   // --pe-slack given explicitly (else mem uses 0: SMEM positions are exact)
	int64_t pe_near_win;   // --pe-near-win: a rescued mate keeps founders whose own projection is within this
	                       // many bp of the anchor, same contig (0 = the concordance window) [100000]
	int8_t pe_trust_exact; // --pe-trust-exact: an unrescuable EXACT + PLACED pair moves the PLACED mate onto the EXACT one
	// refmap error-profile / span filters (all off by default; a rejected EXACT/PLACED read or pair is
	// reported as FILTERED and adds nothing to --ps4g/--npy/--ps4g-per-read/--target-hits)
	int32_t err_max_ev;    // --err-max-events: max mismatch events per read (per pair with --paired); -1 = off, -2 = auto
	int64_t err_calib;     // --err-calib: reads (pairs with --paired) in the calibration pre-pass [1000000]
	double err_alpha;      // --err-alpha: auto cutoff = smallest K with P(Poisson(mu_err) > K) <= alpha [0.001]
	char *err_hist_fn;     // --err-hist: write the calibration histogram + fit here (NULL = off)
	int8_t err_stats;      // --err-stats: append the events and reference-span columns to each row
	int64_t max_ref_span;  // --max-ref-span: max reference span of a read (pair with --paired), bp; 0 = off
	int8_t pe_drop_disc;   // --pe-drop-disc: discordant pairs are FILTERED instead of two single-end events
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
	opt->pav_agree = 5000;  // chain --lift: independent of the grid; see chain_emit_pav
	opt->diverged_rows = 0;     // colinear rows have a real reference coordinate -> ordinary row
	opt->insertion_rows = 1;    // pav: row, grid-snapped -- matches pre-existing chain_emit_pav behavior
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
	opt->anchor_dist_npy = 0;      // off by default (opt-in via --anchor-dist-npy)
	opt->anchor_dist_thresh = 2000; // matches `ropebwt3 lift`'s default anchor stride (-s 2000)
	opt->support_best = 0, opt->support_tol = 0, opt->support_win = 10000; // chain: default = intersection
	opt->max_locate = RB3_RM_MAX_ASSEMBLY;
	opt->pe_near_win = 100000; // lift projections are 2-20 kb noisy; the ~2.7 kb concordance window dropped true founders
	opt->paired = 0, opt->pe_mean = 0.0, opt->pe_sd = 0.0; // single-end by default; insert size estimated per batch
	opt->pe_slack = 2000;  // lift projections are typically within ~2kb (the --anchor-dist-thresh scale)
	opt->err_max_ev = -1, opt->err_calib = 1000000, opt->err_alpha = 1e-3; // error-profile filter off by default
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
	uint8_t *keep; // mem --paired: keep[j] != 0 iff pos[j] has a concordant mate hit (NULL if not computed)
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
	rb3_lift_ridx_t *lift_ridx; // refmap/chain --anchor-dist-npy: reference-position index over `lift`, NULL unless requested
	rb3_gtab_t *gtab;      // refmap/chain --ps4g/--npy: sample (gamete) table, NULL unless requested
	rb3_ps4g_acc_t *ps4g_acc; // refmap/chain --ps4g/--npy: accumulated per-read support events
	FILE *ps4g_per_read_fp;   // refmap --ps4g-per-read: per-read PS4G file, NULL unless requested
	rb3_bed_t *label_bed;  // refmap --label-bed: diploid training labels, NULL unless requested
	rb3_hitcount_t hitcount; // refmap --target-hits: PLACED+EXACT records written so far vs. the target
	rb3_seqio_t *fp2;        // --paired with two files: the R2 stream (NULL: interleaved, or single-end)
	rb3_isize_t isize;       // --paired: insert-size distribution in use (fixed, or the latest good batch estimate)
	int isize_fixed;         // --paired: isize came from --pe-isize; never re-estimated
	int64_t pe_cnt[5];       // --paired: pairs per RB3_PE_* class
	int pe_warned;           // --paired: the "too few pairs to estimate" warning was printed
	int8_t err_on;           // refmap: count mismatch events per read (any --err-* filter/stats flag)
	int32_t err_cut;         // refmap --err-max-events in force (fixed, or calibrated); -1 = off
	int64_t err_nfilt[3];    // refmap: reads (mates) FILTERED by events / span / --pe-drop-disc
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
#define RB3_RM_FILTERED 5 // EXACT/PLACED rejected by an opt-in read filter (--err-max-events/--max-ref-span/--pe-drop-disc)
#define RB3_SPAN_INF INT64_MAX // reference span of a read whose copies/projections lie on different sequences

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
	int32_t core_len; // length of the exact interval used (s->len, or the SMEM core)
	int8_t capped;    // the occurrence sample hit --max-locate: gametes may be incomplete
	int32_t n_proj;   // --paired: alternative reference positions for mate rescue -- per-assembly lift
	                  // projections (PLACED) or every reference copy (EXACT)
	int64_t *proj_rsid, *proj_pos; // allocated with RB3_MALLOC; NULL when n_proj == 0
	int32_t *proj_gid; // PLACED: gamete of each projection (NULL without --ps4g/--npy, and for EXACT)
	int8_t *proj_strand; // EXACT: strand of each reference copy (NULL for PLACED)
	int64_t occ; // --report-occ: raw FM-index interval size for the matched query (whole read, or the
	             // matched SMEM core if there was no end-to-end match) -- 0 for UNPLACED and for
	             // --kmer-mode results (no single whole-query interval exists there). Pangenome-wide:
	             // counts exact matches to the read's given orientation anywhere across the reference
	             // + all assemblies (both strands are indexed, so a genuine inverted/palindromic repeat
	             // is correctly counted more than once -- this is real signal, not an artifact to
	             // divide away; empirically, a single-copy locus present in only one assembly reports
	             // occ=1, not 2 -- verified directly against known unique loci, see --report-occ help).
	int16_t n_ev;  // --err-*: mismatch events against the best carrier of the exact core (-1 = not computed / no core)
	char filt;     // FILTERED reason: 'e' events, 'S' span, 'D' --pe-drop-disc; 0 = kept
	int64_t span;  // reference span of the call: EXACT = first-to-last reference copy, PLACED = spread of
	               // every assembly's lift projection; RB3_SPAN_INF across sequences, -1 = not set
} refmap_rst_t;

typedef struct {
	const pipeline_t *p;
	int32_t n_seq, n_hapdiv;
	m_seq_t *seq;
	rb3_swrst_t *rst, *rst_rev;
	m_hapdiv_t *hapdiv;
	refmap_rst_t *refmap; // refmap: one entry per query
	rb3_pe_res_t *pe;     // --paired: one entry per mate pair (seq[2i], seq[2i+1])
	int64_t pe_max_dist;  // --paired (mem/sw): concordance window used for this batch
	uint8_t **sw_keep;    // --paired (sw): per read, keep flags over its (hit, position) list; NULL = keep all
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
				// rb3_ssa_multi's hit order depends on the SA sampling rate (ssa -s); chain keeps the
				// first reference hit and breaks DP ties by order, so sort for an index-independent result
				qsort(pos, q->n_pos, sizeof(rb3_pos_t), pos_cmp);
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
static int refmap_place_debug = 0; // RB3_PLACE_DEBUG=1: per-assembly projection trace (PLD lines) on stderr

static void refmap_place_lift(void *km, const pipeline_t *p, const m_seq_t *s, refmap_rst_t *r)
{
	const rb3_fmi_t *f = &p->fmi;
	int32_t L = r->n_asm_list > 0? r->n_asm_list : 1; // at most --max-locate
	int64_t *rsids = Kmalloc(km, int64_t, L), *rposs = Kmalloc(km, int64_t, L), *v = Kmalloc(km, int64_t, L), best_rsid = -1, med;
	int32_t n = 0, k, j, bestn, m, *gids = Kmalloc(km, int32_t, L);
	for (k = 0; k < r->n_asm_list; ++k) {
		int64_t clen, st, en, rsid, rpos;
		pos_stranded(f->sid, &r->assemblies[k], s->len, &clen, &st, &en);
		if (rb3_lift_project(p->lift, km, (int32_t)(r->assemblies[k].sid >> 1), st,
							 p->opt->lift_win, p->opt->lift_mad, 4, &rsid, &rpos))
			gids[n] = p->gtab? p->gtab->sid2g[r->assemblies[k].sid >> 1] : -1, rsids[n] = rsid, rposs[n] = rpos, ++n;
		else rsid = -1, rpos = -1;
		if (refmap_place_debug)
			fprintf(stderr, "PLD\t%s\t%d\t%ld\t%s\t%ld\t%d\t%s\t%ld\n", s->name, r->core_len, (long)r->occ,
					f->sid->name[r->assemblies[k].sid >> 1], (long)st, rb3_lift_last_path,
					rsid >= 0? f->sid->name[rsid] : ".", (long)rpos);
	}
	if (n == 0) { r->status = RB3_RM_UNPLACED; goto lift_end; }
	if (p->opt->paired) { // keep every projection: a mate placed far from this one may find one near it
		r->n_proj = n;
		r->proj_rsid = RB3_MALLOC(int64_t, n), r->proj_pos = RB3_MALLOC(int64_t, n);
		memcpy(r->proj_rsid, rsids, n * sizeof(int64_t)), memcpy(r->proj_pos, rposs, n * sizeof(int64_t));
		if (p->gtab) r->proj_gid = RB3_MALLOC(int32_t, n), memcpy(r->proj_gid, gids, n * sizeof(int32_t));
	}
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
	r->span = bestn < n? RB3_SPAN_INF : v[m - 1] - v[0]; // spread of every projection (cf. --placed-max-spread: majority seq only)
lift_end:
	kfree(km, rsids), kfree(km, rposs), kfree(km, v), kfree(km, gids);
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

/*
 * --err-*: mismatch "events" of a read against the best-matching carrier of its exact core.
 * Starting from the core interval refmap places the read with (the whole read, or its longest
 * SMEM), extend it base by base to the right and then to the left in the bidirectional FMD
 * index. While the read's next base extends the interval, the read still matches at least one
 * indexed sequence end to end. When it does not, that is one event, resolved greedily as the
 * edit whose continuation matches the most following read bases (RB3_ERR_LOOK lookahead):
 *   substitution  -- carrier base b != read base (preferred on ties: the Illumina/wgsim error);
 *   read insertion -- skip the read base;
 *   read deletion  -- the carrier has an extra base b before the read base.
 * The final interval is non-empty, so the count is the number of edits that explain the WHOLE read
 * against ONE indexed sequence (a greedy upper bound on that minimum), never a mosaic of SMEMs
 * from different assemblies. If no base extends (the carrier ends: contig end), the overhang
 * is one event. Counting stops at RB3_ERR_CAP.
 */
#define RB3_ERR_LOOK 12

static inline int64_t err_step(const rb3_fmi_t *f, const rb3_sai_t *ik, int c, int back, rb3_sai_t *out)
{
	rb3_sai_t ok[RB3_ASIZE];
	if (c < 1 || c > 4) { out->size = 0; return 0; }
	rb3_fmd_extend(f, ik, ok, back);
	*out = ok[back? c : rb3_comp(c)];
	return out->size;
}

// read bases matched exactly from q[i] on in direction d (+1 right, -1 left), at most lim
static int32_t err_run(const rb3_fmi_t *f, rb3_sai_t ik, const uint8_t *q, int32_t len, int32_t i, int d, int32_t lim)
{
	int32_t n = 0;
	rb3_sai_t ok;
	for (; n < lim && i >= 0 && i < len; i += d, ++n) {
		if (err_step(f, &ik, q[i], d < 0, &ok) == 0) break;
		ik = ok;
	}
	return n;
}

static int32_t err_extend(const rb3_fmi_t *f, rb3_sai_t *ik, const uint8_t *q, int32_t len, int32_t i, int d, int32_t max_ev)
{
	int32_t ev = 0, back = d < 0;
	while (i >= 0 && i < len && ev < max_ev) {
		rb3_sai_t ok[RB3_ASIZE], best_ik, n2;
		int32_t b, rem, lim, best = -1, sc;
		rb3_fmd_extend(f, ik, ok, back);
		if (q[i] >= 1 && q[i] <= 4 && ok[back? q[i] : rb3_comp(q[i])].size > 0) { // still matching
			*ik = ok[back? q[i] : rb3_comp(q[i])], i += d;
			continue;
		}
		++ev;
		rem = d > 0? len - i - 1 : i; // read bases beyond q[i]
		lim = rem < RB3_ERR_LOOK? rem : RB3_ERR_LOOK;
		for (b = 1; b <= 4; ++b) { // substitution
			const rb3_sai_t *x = &ok[back? b : rb3_comp(b)];
			if (b == q[i] || x->size == 0) continue;
			sc = err_run(f, *x, q, len, i + d, d, lim);
			if (sc > best) best = sc, best_ik = *x;
		}
		if (best < 0) break; // the carrier ends here: the overhang is this one event
		if (rem > 0) { // read insertion: q[i] is extra
			sc = err_run(f, *ik, q, len, i + d, d, lim);
			if (sc > best) best = sc, best_ik = *ik;
		}
		if (rem > 0 && q[i] >= 1 && q[i] <= 4) { // read deletion: the carrier has an extra base b, then q[i]
			for (b = 1; b <= 4; ++b) {
				const rb3_sai_t *x = &ok[back? b : rb3_comp(b)];
				if (x->size == 0 || err_step(f, x, q[i], back, &n2) == 0) continue;
				sc = err_run(f, n2, q, len, i + d, d, lim);
				if (sc > best) best = sc, best_ik = n2;
			}
		}
		*ik = best_ik, i += d;
	}
	return ev;
}

// Events of read q[0..len) given the core interval ik over q[st..en) (st == 0 && en == len: 0)
static int32_t refmap_err_events(const rb3_fmi_t *f, int32_t len, const uint8_t *q, const rb3_sai_t *core, int32_t st, int32_t en)
{
	rb3_sai_t ik = *core;
	int32_t ev;
	if (st <= 0 && en >= len) return 0;
	ev = err_extend(f, &ik, q, len, en, +1, RB3_ERR_CAP);
	if (ev < RB3_ERR_CAP) ev += err_extend(f, &ik, q, len, st - 1, -1, RB3_ERR_CAP - ev);
	return ev < RB3_ERR_CAP? ev : RB3_ERR_CAP;
}

// The exact core refmap places a read with: the whole read, else its longest SMEM (first one on
// ties). Returns 0 if nothing of the read occurs in any genome; *st/*en = the core's read span.
static int refmap_core(void *km, const pipeline_t *p, const m_seq_t *s, rb3_sai_t *Iq, int32_t *st, int32_t *en)
{
	const rb3_fmi_t *f = &p->fmi;
	if (!refmap_query_interval(f, s->len, s->seq, Iq) || Iq->size == 0) {
		// no end-to-end match (e.g. a sequencing error): fall back to the longest exact core (SMEM)
		rb3_sai_v mem = {0,0,0};
		size_t bi;
		int32_t bestlen = 0;
		rb3_fmd_smem_TG(km, f, s->len, s->seq, &mem, 1, p->opt->min_len);
		for (bi = 0; bi < mem.n; ++bi) {
			int32_t a = mem.a[bi].info>>32, b = (int32_t)mem.a[bi].info;
			if (b - a > bestlen) bestlen = b - a, *Iq = mem.a[bi], *st = a, *en = b;
		}
		kfree(km, mem.a);
		if (bestlen == 0 || Iq->size == 0) return 0; // nothing of the query occurs in any genome
	} else *st = 0, *en = s->len;
	return 1;
}

static void refmap_query(void *km, const pipeline_t *p, const m_seq_t *s, refmap_rst_t *r)
{
	const rb3_fmi_t *f = &p->fmi;
	rb3_sai_t Iq;
	if (p->opt->kmer_len > 0) { refmap_query_kmer(km, p, s, r); return; }
	rb3_pos_t *pos;
	int64_t np, i;
	int32_t L = p->opt->max_locate, *asm_seen, n_asm = 0;

	memset(r, 0, sizeof(*r));
	r->status = RB3_RM_UNPLACED, r->strand = 0, r->qlen = s->len;
	r->ref_sid = -1, r->cL = r->cR = -1;

	r->n_ev = -1, r->span = -1;
	{
		int32_t st = 0, en = 0;
		if (!refmap_core(km, p, s, &Iq, &st, &en)) return; // nothing of the query occurs in any genome
		r->core_len = en - st;
		if (p->err_on) r->n_ev = refmap_err_events(f, s->len, s->seq, &Iq, st, en);
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
	pos = Kmalloc(km, rb3_pos_t, L);
	np = rb3_ssa_multi(km, f, f->ssa, Iq.x[0], Iq.x[0] + Iq.size, L, pos);
	r->capped = (np >= L);
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
		if (ref_found) { // span: first-to-last reference copy located (RB3_SPAN_INF across sequences)
			int64_t lo = INT64_MAX, hi = -1, rs = -1;
			for (i = 0; i < np && r->span != RB3_SPAN_INF; ++i) {
				int64_t clen, st, en;
				if (!p->is_ref[pos[i].sid>>1]) continue;
				if (rs >= 0 && rs != pos[i].sid>>1) { r->span = RB3_SPAN_INF; break; }
				pos_stranded(f->sid, &pos[i], s->len, &clen, &st, &en);
				rs = pos[i].sid>>1, lo = st < lo? st : lo, hi = en > hi? en : hi;
			}
			if (r->span != RB3_SPAN_INF) r->span = hi - lo;
		}
		if (ref_found && p->opt->paired) { // every reference copy: a mate can pick the right one
			int32_t m = 0;
			r->proj_rsid = RB3_MALLOC(int64_t, np), r->proj_pos = RB3_MALLOC(int64_t, np), r->proj_strand = RB3_MALLOC(int8_t, np);
			for (i = 0; i < np; ++i) {
				int64_t clen, st, en;
				if (!p->is_ref[pos[i].sid>>1]) continue;
				pos_stranded(f->sid, &pos[i], s->len, &clen, &st, &en);
				r->proj_rsid[m] = pos[i].sid>>1, r->proj_pos[m] = st, r->proj_strand[m++] = pos[i].sid&1;
			}
			r->n_proj = m;
		}
		if (ref_found) {
			if (p->gtab) { // equivalent hits: every occurrence of an EXACT read names a candidate parent
				r->gametes = RB3_MALLOC(int32_t, np > 0? np : 1); // rb3_ps4g_acc_add sorts+dedupes on ingestion
				for (i = 0; i < np; ++i) r->gametes[i] = p->gtab->sid2g[pos[i].sid>>1];
				r->n_gametes = (int32_t)np;
			}
			kfree(km, pos);
			return;
		}
	}
	r->assemblies = RB3_CALLOC(rb3_pos_t, L);
	asm_seen = Kmalloc(km, int32_t, L);
	for (i = 0; i < np && n_asm < L; ++i) { // keep distinct assembly sequences for reporting
		int32_t j, dup = 0;
		for (j = 0; j < n_asm; ++j)
			if (asm_seen[j] == (int)(pos[i].sid>>1)) { dup = 1; break; }
		if (dup) continue;
		asm_seen[n_asm] = pos[i].sid>>1;
		r->assemblies[n_asm++] = pos[i];
	}
	r->n_asm_list = n_asm, r->n_assembly = n_asm;
	kfree(km, asm_seen);
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

static void write_refmap1(kstring_t *out, const rb3_fmi_t *f, const m_seq_t *s, const refmap_rst_t *r, int kmer, int report_occ, int err_stats)
{
	static const char *status_str[6] = { "UNPLACED", "PLACED", "ONE_SIDE", "EXACT", "MULTI", "FILTERED" };
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
	if (err_stats) { // --err-stats: events, reference span (inf = across sequences, . = unset), filter reason
		rb3_sprintf_lite(out, "\t%d", r->n_ev);
		if (r->span == RB3_SPAN_INF) rb3_sprintf_lite(out, "\tinf");
		else if (r->span < 0) rb3_sprintf_lite(out, "\t.");
		else rb3_sprintf_lite(out, "\t%ld", (long)r->span);
		rb3_sprintf_lite(out, "\t%c", r->filt? r->filt : '.');
	}
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

// --paired: the single PS4G/npy event of a pair, at the anchor mate a's position, carrying the
// founders both mates share (rb3_pe_event_gametes). --ps4g-per-read gets one line, named after
// the fragment (the mate name without "/1" or "/2").
static void refmap_pair_event(const pipeline_t *p, const m_seq_t *s, const refmap_rst_t *a, const refmap_rst_t *b)
{
	int32_t *g, n, k;
	if (a->status == RB3_RM_FILTERED) return; // --err-max-events/--max-ref-span/--pe-drop-disc rejected the pair
	if (a->n_gametes <= 0 && b->n_gametes <= 0) return;
	g = RB3_MALLOC(int32_t, a->n_gametes + b->n_gametes);
	n = rb3_pe_event_gametes(a->gametes, a->n_gametes, a->capped, b->gametes, b->n_gametes, b->capped, g);
	if (n > 0 && p->ps4g_acc) rb3_ps4g_acc_add(p->ps4g_acc, a->ref_sid, a->cL, g, n);
	if (n > 0 && p->ps4g_per_read_fp) {
		const char *nm = p->fmi.sid->name[a->ref_sid], *us = strchr(nm, '_');
		fprintf(p->ps4g_per_read_fp, "%.*s\t%s\t%ld\t", s->name? rb3_pe_name_len(s->name) : 1, s->name? s->name : "?",
				us? us + 1 : nm, (long)a->cL);
		for (k = 0; k < n; ++k) fprintf(p->ps4g_per_read_fp, "%s%d", k? "," : "", g[k]);
		fputc('\n', p->ps4g_per_read_fp);
	}
	free(g);
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
				write_refmap1(&out, f, s, &r->sub[k], kmer, p->opt->report_occ, p->opt->err_stats);
				fputs(out.s, stdout);
				refmap_rst_accumulate(p->ps4g_acc, &r->sub[k]);
				write_ps4g_read(p->ps4g_per_read_fp, f, s, &r->sub[k]);
				refmap_count_hit(p, r->sub[k].status);
				free(r->sub[k].gametes);
			}
			free(r->sub);
		} else {
			write_refmap1(&out, f, s, r, kmer, p->opt->report_occ, p->opt->err_stats);
			if (t->pe) { // --paired: one extra column, the pair class (MOVED = this mate was rescued)
				const rb3_pe_res_t *e = &t->pe[j >> 1];
				out.l--; // drop the newline
				rb3_sprintf_lite(&out, "\t%s\n", e->cls == RB3_PE_RESC && e->rescued == (j & 1)? "MOVED" : rb3_pe_class_str[(int)e->cls]);
			}
			fputs(out.s, stdout);
			if (t->pe == 0 || t->pe[j >> 1].cls == RB3_PE_DISC) { // single-end, or a discordant mate: its own event
				refmap_rst_accumulate(p->ps4g_acc, r);
				write_ps4g_read(p->ps4g_per_read_fp, f, s, r);
			} else if (t->pe[j >> 1].anchor == (j & 1)) // one event per pair, at the anchor mate
				refmap_pair_event(p, s, r, &t->refmap[j ^ 1]);
			refmap_count_hit(p, r->status);
		}
	}
	for (j = 0; j < t->n_seq; ++j) { // after all rows: a pair event reads its mate's gametes
		refmap_rst_t *r = &t->refmap[j];
		free(r->assemblies);
		free(r->gametes);
		free(r->proj_rsid);
		free(r->proj_pos);
		free(r->proj_gid);
		free(r->proj_strand);
		free(t->seq[j].seq);
		free(t->seq[j].name);
	}
	free(out.s);
}

/********************************************************************
 * refmap --err-max-events / --max-ref-span / --pe-drop-disc
 ********************************************************************/
static inline int refmap_is_placed(const refmap_rst_t *r) { return r->status == RB3_RM_EXACT || r->status == RB3_RM_PLACED; }

static inline void refmap_filter1(pipeline_t *p, refmap_rst_t *r, char why, int which)
{
	if (!refmap_is_placed(r)) return;
	r->status = RB3_RM_FILTERED, r->filt = why, ++p->err_nfilt[which];
}

static inline int64_t refmap_se_span(const refmap_rst_t *r) // walk path / kmer: the bracket
{
	return r->span >= 0? r->span : r->cR >= r->cL && r->cL >= 0? r->cR - r->cL : 0;
}

static inline int32_t refmap_pair_ev(const refmap_rst_t *a, const refmap_rst_t *b) // no core counts as the cap
{
	int32_t e = (a->n_ev < 0? RB3_ERR_CAP : a->n_ev) + (b->n_ev < 0? RB3_ERR_CAP : b->n_ev);
	return e < RB3_ERR_CAP? e : RB3_ERR_CAP;
}

// Step 1, before pairing: single-end event + span filters, or the pair-sum event filter.
static void refmap_filter_pre(pipeline_t *p, step_t *t)
{
	int32_t j;
	if (!p->opt->paired) {
		for (j = 0; j < t->n_seq; ++j) {
			refmap_rst_t *r = &t->refmap[j];
			if (!refmap_is_placed(r)) continue;
			if (p->err_cut >= 0 && (r->n_ev < 0 || r->n_ev > p->err_cut)) refmap_filter1(p, r, 'e', 0);
			else if (p->opt->max_ref_span > 0 && refmap_se_span(r) > p->opt->max_ref_span) refmap_filter1(p, r, 'S', 1);
		}
	} else if (p->err_cut >= 0) {
		for (j = 0; j + 1 < t->n_seq; j += 2)
			if (refmap_pair_ev(&t->refmap[j], &t->refmap[j+1]) > p->err_cut)
				refmap_filter1(p, &t->refmap[j], 'e', 0), refmap_filter1(p, &t->refmap[j+1], 'e', 0);
	}
}

// Reference span of a resolved pair: first-to-last coordinate of both mates (CONC/RESC/DISC on one
// sequence), RB3_SPAN_INF for a DISC pair across sequences, the placed mate's own span for SOLO.
static int64_t refmap_pair_span(const rb3_pe_res_t *e, const refmap_rst_t *a, const refmap_rst_t *b)
{
	int64_t lo, hi;
	if (e->cls == RB3_PE_NONE) return -1;
	if (e->cls == RB3_PE_SOLO) return refmap_se_span(e->anchor? b : a);
	if (a->ref_sid != b->ref_sid) return RB3_SPAN_INF;
	lo = a->cL < b->cL? a->cL : b->cL, hi = a->cR > b->cR? a->cR : b->cR;
	return hi - lo;
}

// Step 1, after pairing: --pe-drop-disc and the pair span cap (both mates FILTERED).
static void refmap_filter_pairs(pipeline_t *p, step_t *t)
{
	int32_t i;
	for (i = 0; i < t->n_seq / 2; ++i) {
		const rb3_pe_res_t *e = &t->pe[i];
		refmap_rst_t *a = &t->refmap[2*i], *b = &t->refmap[2*i+1];
		if (p->opt->pe_drop_disc && e->cls == RB3_PE_DISC)
			refmap_filter1(p, a, 'D', 2), refmap_filter1(p, b, 'D', 2);
		else if (p->opt->max_ref_span > 0 && refmap_pair_span(e, a, b) > p->opt->max_ref_span)
			refmap_filter1(p, a, 'S', 1), refmap_filter1(p, b, 'S', 1);
	}
}

// Calibration pre-pass (--err-max-events=auto, --err-hist): events of the first opt->err_calib reads
// (pairs) of the input, exactly as the main pass counts them; fit the error model; set p->err_cut.
typedef struct {
	const pipeline_t *p;
	m_seq_t *seq;
	int16_t *ev;
	void **km;
} err_cal_t;

static void err_cal_worker(void *data, long i, int tid)
{
	err_cal_t *c = (err_cal_t*)data;
	m_seq_t *s = &c->seq[i];
	rb3_sai_t Iq;
	int32_t st = 0, en = 0;
	rb3_char2nt6(s->len, s->seq);
	if (c->p->opt->trim_polya > 0) m_trim_polya(s, c->p->opt->trim_polya);
	c->ev[i] = refmap_core(c->km[tid], c->p, s, &Iq, &st, &en)? refmap_err_events(&c->p->fmi, s->len, s->seq, &Iq, st, en) : -1;
}

static int refmap_err_calibrate(pipeline_t *p, const char *fn1, const char *fn2, int is_line)
{
	const rb3_mopt_t *o = p->opt;
	int32_t nm = o->paired? 2 : 1, chunk = 200000, k, cut, i;
	int64_t n_frag = 0, hist_r[RB3_ERR_CAP+1], hist_p[RB3_ERR_CAP+1];
	rb3_seqio_t *fp, *fp2 = 0;
	err_cal_t c;
	rb3_errfit_t fit_r, fit_p;
	memset(hist_r, 0, sizeof(hist_r)), memset(hist_p, 0, sizeof(hist_p));
	if ((fp = rb3_seq_open(fn1, is_line)) == 0 || (fn2 && (fp2 = rb3_seq_open(fn2, is_line)) == 0)) {
		if (rb3_verbose >= 1) fprintf(stderr, "ERROR: --err-*: failed to open '%s' for calibration\n", fp? fn2 : fn1);
		if (fp) rb3_seq_close(fp);
		return -1;
	}
	c.p = p, c.seq = RB3_CALLOC(m_seq_t, chunk), c.ev = RB3_CALLOC(int16_t, chunk), c.km = RB3_CALLOC(void*, o->n_threads);
	for (i = 0; i < o->n_threads; ++i) c.km[i] = km_init();
	while (n_frag < o->err_calib) {
		int32_t n = 0;
		while (n + nm <= chunk && n_frag < o->err_calib) {
			int64_t len;
			const char *name;
			char *ss;
			for (k = 0; k < nm; ++k) {
				if ((ss = rb3_seq_read1(k == 1 && fp2? fp2 : fp, &len, &name)) == 0) break;
				c.seq[n + k].seq = (uint8_t*)rb3_strdup(ss), c.seq[n + k].len = len;
			}
			if (k < nm) { for (i = 0; i < k; ++i) free(c.seq[n + i].seq); break; }
			n += nm, ++n_frag;
		}
		if (n == 0) break;
		kt_for(o->n_threads, err_cal_worker, &c, n);
		for (i = 0; i < n; ++i) {
			if (c.ev[i] >= 0) ++hist_r[c.ev[i] < RB3_ERR_CAP? c.ev[i] : RB3_ERR_CAP];
			free(c.seq[i].seq);
		}
		if (nm == 2)
			for (i = 0; i + 1 < n; i += 2) {
				refmap_rst_t a, b;
				a.n_ev = c.ev[i], b.n_ev = c.ev[i+1];
				++hist_p[refmap_pair_ev(&a, &b)];
			}
		if (n < chunk - nm + 1 && n_frag < o->err_calib) break; // input exhausted
	}
	for (i = 0; i < o->n_threads; ++i) km_destroy(c.km[i]);
	free(c.km), free(c.seq), free(c.ev);
	rb3_seq_close(fp);
	if (fp2) rb3_seq_close(fp2);
	rb3_err_fit(hist_r, RB3_ERR_CAP, &fit_r);
	if (nm == 2) rb3_err_fit(hist_p, RB3_ERR_CAP, &fit_p);
	{
		const rb3_errfit_t *fu = nm == 2? &fit_p : &fit_r;
		cut = fu->ok? rb3_err_cutoff(fu->mu0, o->err_alpha, RB3_ERR_CAP) : RB3_ERR_CAP - 1;
		if (o->err_max_ev == -2) p->err_cut = cut;
		if (rb3_verbose >= 3)
			fprintf(stderr, "[M::%s] %ld %s: error component pi0=%.4f mu0=%.4f, tail mu1=%.3f; P(Poisson(mu0) > K) <= %g at K=%d%s\n",
					__func__, (long)n_frag, nm == 2? "pairs" : "reads", fu->pi0, fu->mu0, fu->mu1, o->err_alpha, cut,
					o->err_max_ev == -2? " (in force)" : "");
	}
	if (o->err_hist_fn) {
		FILE *hf = fopen(o->err_hist_fn, "w");
		if (hf == 0) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: failed to open --err-hist file '%s'\n", o->err_hist_fn);
			return -1;
		}
		fprintf(hf, "#ropebwt3 refmap error profile: events = edits explaining the read against one carrier of its exact core\n");
		fprintf(hf, "#calibration %s=%ld\n", nm == 2? "pairs" : "reads", (long)n_frag);
		fprintf(hf, "unit\tk\tcount\tfrac\tcum_frac\terr_expected\tpost_err\n");
		rb3_err_write_hist(hf, "read", hist_r, RB3_ERR_CAP, &fit_r, o->err_alpha, rb3_err_cutoff(fit_r.mu0, o->err_alpha, RB3_ERR_CAP));
		if (nm == 2) rb3_err_write_hist(hf, "pair", hist_p, RB3_ERR_CAP, &fit_p, o->err_alpha, rb3_err_cutoff(fit_p.mu0, o->err_alpha, RB3_ERR_CAP));
		fclose(hf);
	}
	return 0;
}

/********************************************************************
 * --paired: mate pairing over one batch (refmap)
 ********************************************************************/
#define RB3_PE_MIN_PAIRS       25   // fewest EXACT-EXACT pairs for a batch's own insert-size estimate
#define RB3_PE_DEFAULT_MAX_INS 1000 // max insert used until some batch yields an estimate
#define RB3_PE_EST_WIN 10000 // mem/chain: largest span accepted into the insert-size sample

// Update p->isize from one batch's fragment spans (span[] is sorted in place); warn once while
// no batch has had enough pairs. Runs in pipeline step 1 (serialized across batches).
static void pe_update_isize(pipeline_t *p, int64_t *span, int32_t n_span)
{
	rb3_isize_t is;
	if (rb3_pe_isize_estimate(span, n_span, RB3_PE_MIN_PAIRS, &is)) {
		p->isize = is;
		if (rb3_verbose >= 3)
			fprintf(stderr, "[M::%s] insert size from %d pairs: mean %.1f, sd %.1f, max %ld\n", __func__,
					is.n, is.mean, is.sd, (long)is.max_ins);
	} else if (!p->isize.ok && !p->pe_warned && rb3_verbose >= 2) {
		fprintf(stderr, "WARNING: only %d confidently placed pairs in a batch (need %d) to estimate the insert size; "
						"using max insert %d until one does (or pass --pe-isize)\n", n_span, RB3_PE_MIN_PAIRS, RB3_PE_DEFAULT_MAX_INS);
		p->pe_warned = 1;
	}
}

static inline int8_t refmap_pe_rank(const refmap_rst_t *r)
{
	if (r->ref_sid < 0) return 0;
	return r->status == RB3_RM_EXACT? 2 : r->status == RB3_RM_PLACED? 1 : 0;
}

// Estimate the insert size from the batch's confidently placed pairs (both mates EXACT on one
// contig, opposite strands; span = fragment end - fragment start), then classify every pair
// and move rescued mates next to their anchor. Runs in pipeline step 1, which kt_pipeline
// serializes across batches, so updating p->isize / p->pe_cnt here is race-free.
static void refmap_pair_batch(pipeline_t *p, step_t *t)
{
	int32_t i, n_pair = t->n_seq / 2, n_span = 0;
	int64_t max_dist;
	t->pe = RB3_CALLOC(rb3_pe_res_t, n_pair > 0? n_pair : 1);
	if (!p->isize_fixed) {
		int64_t *span = RB3_MALLOC(int64_t, n_pair > 0? n_pair : 1);
		for (i = 0; i < n_pair; ++i) {
			const refmap_rst_t *a = &t->refmap[2*i], *b = &t->refmap[2*i+1];
			if (a->status != RB3_RM_EXACT || b->status != RB3_RM_EXACT) continue;
			if (a->ref_sid != b->ref_sid || a->strand == b->strand) continue;
			span[n_span++] = (a->cR > b->cR? a->cR : b->cR) - (a->cL < b->cL? a->cL : b->cL);
		}
		pe_update_isize(p, span, n_span);
		free(span);
	}
	max_dist = (p->isize.ok? p->isize.max_ins : RB3_PE_DEFAULT_MAX_INS) + p->opt->pe_slack;
	for (i = 0; i < n_pair; ++i) {
		refmap_rst_t *r[2];
		rb3_pe_mate_t m[2];
		int k;
		r[0] = &t->refmap[2*i], r[1] = &t->refmap[2*i+1];
		for (k = 0; k < 2; ++k) {
			m[k].rank = refmap_pe_rank(r[k]);
			m[k].rsid = r[k]->ref_sid, m[k].pos = r[k]->cL;
			m[k].n_proj = r[k]->n_proj, m[k].proj_rsid = r[k]->proj_rsid, m[k].proj_pos = r[k]->proj_pos;
			m[k].proj_gid = r[k]->proj_gid;
		}
		rb3_pe_resolve(m, max_dist, p->opt->pe_trust_exact, &t->pe[i]);
		if (t->pe[i].cls == RB3_PE_RESC) { // move the mate; keep only founders whose own copy projected there
			const rb3_pe_res_t *e = &t->pe[i];
			refmap_rst_t *q = r[(int)e->rescued];
			const rb3_pe_mate_t *a = &m[(int)e->anchor];
			if (q->proj_gid && q->gametes) {
				int64_t nw = p->opt->pe_near_win > 0? p->opt->pe_near_win : max_dist;
				int32_t ng = rb3_pe_near_gametes(&m[(int)e->rescued], a->rsid, a->pos, nw, q->gametes);
				if (ng > 0) q->n_gametes = ng; // q->gametes holds n_asm >= n_proj entries
			}
			if (q->status == RB3_RM_EXACT) { // moved to another reference copy: keep its length, take its strand
				int32_t c;
				for (c = 0; c < q->n_proj; ++c)
					if (q->proj_rsid[c] == e->resc_rsid && q->proj_pos[c] == e->resc_pos) { q->strand = q->proj_strand[c]; break; }
				q->cR = e->resc_pos + (q->cR - q->cL), q->cL = e->resc_pos, q->ref_sid = e->resc_rsid;
			} else q->ref_sid = e->resc_rsid, q->cL = q->cR = e->resc_pos;
		}
		++p->pe_cnt[(int)t->pe[i].cls];
	}
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
// chain --dup-count: a collapsed duplicate read is named "<name>;x=N" (N identical reads); it adds N
// to the PS4G/npy counts. Without the option, or without a valid suffix, every read counts once.
static int32_t chain_read_weight(const pipeline_t *p, const m_seq_t *s)
{
	const char *q;
	long n;
	if (!p->opt->dup_count || s->name == 0 || (q = strstr(s->name, ";x=")) == 0) return 1;
	n = strtol(q + 3, 0, 10);
	return n > 0 && n < INT32_MAX? (int32_t)n : 1;
}

// One chain output row. Single-end prints it at once; --paired collects a read's rows in a sink so
// both mates can be compared before anything is printed or counted.
typedef struct {
	int32_t rsid;     // reference contig
	int64_t pos, hi;  // row position (segment start / projected position) and segment end
	int32_t *acc, na; // founder (gamete) set, sorted
	int8_t pav;       // 1: "pav:"-prefixed, grid-snapped insertion row
	int8_t extra;     // 5th column of a pav: row (insertion/diverged flag), -1 if none
	int8_t lift;      // 1: from the lift-projection (PAV) path
} chain_row_t;

typedef struct {
	chain_row_t *a;
	int32_t n, m;
	int8_t strand;    // reference chain's strand
	int32_t bases;    // query bases in the reference chain
} chain_sink_t;

static void chain_print_row(const pipeline_t *p, const char *name, int name_len, const chain_row_t *r, int32_t w, const char *cls)
{
	const char *nm = p->fmi.sid->name[r->rsid], *us = strchr(nm, '_');
	int32_t i;
	if (p->ps4g_acc && r->na > 0) rb3_ps4g_acc_add_w(p->ps4g_acc, r->rsid, r->pos, r->acc, r->na, w);
	printf("%.*s\t%s%s\t%ld\t", name_len, name, r->pav? "pav:" : "", us? us + 1 : nm, (long)r->pos);
	for (i = 0; i < r->na; ++i) printf("%s%d", i? "," : "", r->acc[i]);
	if (r->extra >= 0) printf("\t%d", r->extra);
	if (cls) printf("\t%s", cls);
	putchar('\n');
}

// Emit (sk == NULL) or collect one row of read s.
static void chain_put_row(const pipeline_t *p, chain_sink_t *sk, const m_seq_t *s, int32_t rsid, int64_t pos, int64_t hi,
						  const int32_t *acc, int32_t na, int pav, int extra, int lift)
{
	chain_row_t r;
	r.rsid = rsid, r.pos = pos, r.hi = hi, r.na = na, r.pav = pav, r.extra = extra, r.lift = lift;
	if (sk == 0) {
		const char *nm = s->name? s->name : "?";
		r.acc = (int32_t*)acc;
		chain_print_row(p, nm, (int)strlen(nm), &r, chain_read_weight(p, s), 0);
		return;
	}
	r.acc = RB3_MALLOC(int32_t, na > 0? na : 1);
	memcpy(r.acc, acc, na * sizeof(int32_t));
	RB3_GROW(chain_row_t, sk->a, sk->n, sk->m);
	sk->a[sk->n++] = r;
}

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
static void chain_emit_pav(const pipeline_t *p, const m_seq_t *s, void *km, chain_sink_t *sk)
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
		// Agreement and resolution are separate concerns and now separate knobs. They used
		// to share --pav-grid, which made --pav-grid=0 read as "no constraint" while actually
		// being the STRICTEST setting: with no grid the tolerance fell back to --gap-intron
		// (30bp), so insertions carried by several assemblies were silently suppressed.
		if (nmaj > 0 && bestn == np && disp <= p->opt->pav_agree) {
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
					chain_put_row(p, sk, s, best_rsid, med, med, acc, na, 0, -1, 1);
					goto pav_done;
				}
			} else {
				// --insertion-rows: same routing as --diverged-rows, but for TRUE insertions. The
				// "ordinary" choice places the read inline on its real reference contig at the exact
				// (unsnapped) median -- this is MORE precise than the grid-snapped pav: row below,
				// which floors rather than centers (bp = med/grid*grid). It exists because the
				// downstream windower groups PS4G/npy rows by contig (ropebwt_npy_to_matrix.py), so a
				// `pav:chr1` row is a separate pseudo-chromosome disconnected from `chr1` -- "ordinary"
				// keeps insertion evidence positionally inline with everything else on the real contig.
				// DEFAULT stays `pav` (grid-snapped, prefixed): this branch changes nothing unless
				// explicitly requested.
				if (p->opt->insertion_rows == 2) goto pav_done;          // drop
				if (p->opt->insertion_rows == 0) {                       // ordinary row, exact median, no snap
					chain_put_row(p, sk, s, best_rsid, med, med, acc, na, 0, -1, 1);
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
			// PS4G/npy accumulate under the REAL ref_sid + snapped bp -- the accumulator has no
			// concept of the `pav:` contig prefix (it is stdout-only, added at
			// rb3_ps4g_npy_finalize/ps4g.c time from sid->name), so this shares bins with ordinary
			// rows on the same contig at the same coarse position. Known gap, not resolved here --
			// see the design note on --insertion-rows above; --insertion-rows=ordinary avoids it.
			// 5th column = row class. 1 = INSERTION: absent from the reference, placed at
			// the nearest reference breakpoint, approximate (presence/absence only).
			// 0 = DIVERGED: present in the reference but too divergent to share a SMEM,
			// so the coordinate is colinear and exact. Only insertions reach here --
			// diverged rows are routed by --diverged-rows (default: ordinary row).
			chain_put_row(p, sk, s, best_rsid, bp, bp, acc, na, 1, is_insertion? 1 : 0, 1);
		}
	}
pav_done:
	for (i = 0; i < n; ++i) free(cs[i].gam);
	free(acc); free(clu); free(cand); free(cs);
}

// chain --support=best: replace the chain-wide intersection by the assemblies that best explain
// the read. Placement is unchanged (the reference-hitting colinear chain); support is counted over
// the chain's SMEMs PLUS the read's assembly-only SMEMs (no reference occurrence -- the ones the
// intersection rule discards: a donor-allele read with a SNP vs the reference is matched longer by
// the assembly carrying that allele). An assembly-only SMEM occurrence counts for its gamete only
// if, with --lift and --support-win > 0, it projects onto the chain's reference contig within
// support_win bp of the chain's reference span (a paralog elsewhere does not count). Per gamete:
// the number of read bases covered by the counted SMEMs that occur in it; the set is every gamete
// within support_tol bases of the best. Non-chain SMEMs that do hit the reference are ignored, as
// in the default rule. Writes the sorted set to acc (capacity n_gamete) and returns its size.
static int32_t chain_support_best(const pipeline_t *p, const m_seq_t *s, void *km, const chain_sm_t *cs,
								  const int32_t *chain, int32_t nc, int64_t rlo, int64_t rhi, int32_t *acc)
{
	const rb3_fmi_t *f = &p->fmi;
	const rb3_gtab_t *gt = p->gtab;
	int32_t ng = gt->n_gamete, len = s->len, i, k, g, best = 0, na = 0;
	int32_t rsid0 = cs[chain[0]].rsid;
	uint8_t *mask;
	int32_t *cov;
	if (len <= 0) return 0;
	mask = RB3_CALLOC(uint8_t, (int64_t)ng * len);
	cov = RB3_CALLOC(int32_t, ng);
	for (i = 0; i < nc; ++i) { // the placement chain's SMEMs: every gamete they occur in
		const chain_sm_t *c = &cs[chain[i]];
		for (k = 0; k < c->n_gam; ++k)
			memset(mask + (int64_t)c->gam[k] * len + c->qs, 1, c->qe - c->qs);
	}
	for (i = 0; i < s->n_mem; ++i) { // assembly-only SMEMs near the chain locus
		const m_sai_pos_t *r = &s->mem[i];
		int32_t st = r->mem.info>>32, en = (int32_t)r->mem.info, nref = 0;
		if (r->n_pos == 0 || (int64_t)r->mem.size > p->opt->chain_max_occ) continue;
		for (k = 0; k < r->n_pos; ++k) if (p->is_ref[r->pos[k].sid>>1]) ++nref;
		if (nref > 0) continue;
		for (k = 0; k < r->n_pos; ++k) {
			const rb3_pos_t *t = &r->pos[k];
			g = gt->sid2g[t->sid>>1];
			if (g < 0 || g >= ng) continue;
			if (p->lift && p->opt->support_win > 0) {
				int64_t rsid, rpos;
				if (!rb3_lift_project(p->lift, km, (int32_t)(t->sid>>1), pav_fwd_pos(f, t, en - st),
									  p->opt->lift_win, p->opt->lift_mad, 4, &rsid, &rpos)) continue;
				if (rsid != rsid0 || rpos < rlo - p->opt->support_win || rpos > rhi + p->opt->support_win) continue;
			}
			memset(mask + (int64_t)g * len + st, 1, en - st);
		}
	}
	for (g = 0; g < ng; ++g) {
		const uint8_t *m = mask + (int64_t)g * len;
		for (k = 0; k < len; ++k) cov[g] += m[k];
		if (cov[g] > best) best = cov[g];
	}
	if (best > 0)
		for (g = 0; g < ng; ++g)
			if (cov[g] > 0 && cov[g] >= best - p->opt->support_tol) acc[na++] = g;
	free(mask); free(cov);
	return na;
}

static void chain_emit(const pipeline_t *p, const m_seq_t *s, void *km, chain_sink_t *sk)
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
		if (p->lift) chain_emit_pav(p, s, km, sk);
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
			if (cs[j].rsid != cs[i].rsid) continue;       // different contigs: positions are not comparable
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
			for (i = 0; i < nc; ++i) ivlo[i] = cs[chain[i]].rpos, ivhi[i] = cs[chain[i]].rpos + (cs[chain[i]].qe - cs[chain[i]].qs);
			if (p->opt->support_best) { // founder set = best-explaining assemblies (see chain_support_best)
				int64_t rlo = ivlo[0], rhi = ivhi[0];
				for (i = 1; i < nc; ++i) { if (ivlo[i] < rlo) rlo = ivlo[i]; if (ivhi[i] > rhi) rhi = ivhi[i]; }
				free(acc);
				acc = RB3_MALLOC(int32_t, gt->n_gamete);
				na = chain_support_best(p, s, km, cs, chain, nc, rlo, rhi, acc);
			}
			if (na > 0) { // segment forward intervals by intron-sized gaps, emit one row per segment
				int64_t lo, hi;
				if (sk) { // --paired: the chain's strand and size, for comparing mates
					sk->strand = strand;
					for (i = 0, sk->bases = 0; i < nc; ++i) sk->bases += cs[chain[i]].qe - cs[chain[i]].qs;
				}
				for (i = 0; i < nc - 1; ++i) // insertion sort ascending by start (nc is tiny)
					for (j = i + 1; j < nc; ++j)
						if (ivlo[j] < ivlo[i]) { int64_t t0 = ivlo[i], t1 = ivhi[i]; ivlo[i] = ivlo[j], ivhi[i] = ivhi[j], ivlo[j] = t0, ivhi[j] = t1; }
				lo = ivlo[0], hi = ivhi[0];
				for (i = 1; i < nc; ++i) {
					if (ivlo[i] - hi > p->opt->gap_intron) { // intron gap -> flush this exon segment
						chain_put_row(p, sk, s, cs[chain[0]].rsid, lo, hi, acc, na, 0, -1, 0);
						lo = ivlo[i], hi = ivhi[i];
					} else if (ivhi[i] > hi) hi = ivhi[i];
				}
				chain_put_row(p, sk, s, cs[chain[0]].rsid, lo, hi, acc, na, 0, -1, 0);
			}
			free(ivlo); free(ivhi);
			free(acc);
		}
	}
	for (i = 0; i < n; ++i) free(cs[i].gam);
	for (i = 0; i < s->n_mem; ++i) free(s->mem[i].pos);
	free(chain); free(cs);
}

/********************************************************************
 * --paired: chain
 ********************************************************************/
#define RB3_PE_LIFT_SLACK 2000 // extra gap allowed per lift-projected (PAV) mate: projections are ~kb-accurate

static void chain_sink_span(const chain_sink_t *sk, rb3_pe_chain_t *c)
{
	int32_t i;
	memset(c, 0, sizeof(*c));
	if (sk->n == 0) { c->kind = RB3_PE_CHAIN_NONE, c->rsid = -1; return; }
	c->kind = sk->a[0].lift? RB3_PE_CHAIN_LIFT : RB3_PE_CHAIN_REF;
	c->rsid = sk->a[0].rsid, c->lo = sk->a[0].pos, c->hi = sk->a[0].hi, c->strand = sk->strand, c->bases = sk->bases;
	for (i = 1; i < sk->n; ++i) {
		if (sk->a[i].pos < c->lo) c->lo = sk->a[i].pos;
		if (sk->a[i].hi > c->hi) c->hi = sk->a[i].hi;
	}
}

static void chain_write_paired(pipeline_t *p, step_t *t, void *km)
{
	int32_t i, j, n_pair = t->n_seq / 2, n_span = 0;
	chain_sink_t *sk = RB3_CALLOC(chain_sink_t, t->n_seq > 0? t->n_seq : 1);
	int64_t *span = RB3_MALLOC(int64_t, n_pair > 0? n_pair : 1), max_gap;
	t->pe = RB3_CALLOC(rb3_pe_res_t, n_pair > 0? n_pair : 1);
	for (j = 0; j < t->n_seq; ++j) chain_emit(p, &t->seq[j], km, &sk[j]); // collect, print nothing yet
	if (!p->isize_fixed) { // insert size: both mates single-segment reference chains, same contig, opposite strands
		for (i = 0; i < n_pair; ++i) {
			const chain_sink_t *a = &sk[2*i], *b = &sk[2*i+1];
			if (a->n != 1 || b->n != 1 || a->a[0].lift || b->a[0].lift) continue;
			if (a->a[0].rsid != b->a[0].rsid || a->strand == b->strand) continue;
			span[n_span] = (a->a[0].hi > b->a[0].hi? a->a[0].hi : b->a[0].hi) - (a->a[0].pos < b->a[0].pos? a->a[0].pos : b->a[0].pos);
			if (span[n_span] < RB3_PE_EST_WIN) ++n_span;
		}
		pe_update_isize(p, span, n_span);
	}
	max_gap = (p->isize.ok? p->isize.max_ins : RB3_PE_DEFAULT_MAX_INS) + p->opt->max_intron + p->opt->pe_slack;
	for (i = 0; i < n_pair; ++i) {
		const m_seq_t *s0 = &t->seq[2*i];
		chain_sink_t *k[2];
		rb3_pe_chain_t c[2];
		rb3_pe_res_t *e = &t->pe[i];
		const char *cls;
		k[0] = &sk[2*i], k[1] = &sk[2*i+1];
		chain_sink_span(k[0], &c[0]), chain_sink_span(k[1], &c[1]);
		rb3_pe_chain_resolve(c, max_gap, RB3_PE_LIFT_SLACK, e);
		++p->pe_cnt[(int)e->cls];
		cls = rb3_pe_class_str[(int)e->cls];
		if (e->cls == RB3_PE_DISC) { // each mate on its own, as single-end
			for (j = 0; j < 2; ++j) {
				const m_seq_t *s = &t->seq[2*i+j];
				const char *nm = s->name? s->name : "?";
				int32_t r;
				for (r = 0; r < k[j]->n; ++r) chain_print_row(p, nm, (int)strlen(nm), &k[j]->a[r], chain_read_weight(p, s), cls);
			}
		} else if (e->cls == RB3_PE_CONC || e->cls == RB3_PE_SOLO) { // one event: the anchor's first row
			const chain_sink_t *a = k[(int)e->anchor], *b = k[e->anchor ^ 1];
			const char *nm = s0->name? s0->name : "?";
			chain_row_t r = a->a[0];
			int32_t *g = RB3_MALLOC(int32_t, a->a[0].na + (b->n > 0? b->a[0].na : 0) + 1);
			r.na = rb3_pe_event_gametes(a->a[0].acc, a->a[0].na, 0, b->n > 0? b->a[0].acc : 0, b->n > 0? b->a[0].na : 0, 0, g);
			r.acc = g;
			chain_print_row(p, nm, rb3_pe_name_len(nm), &r, chain_read_weight(p, s0), cls);
			free(g);
		}
	}
	for (j = 0; j < t->n_seq; ++j) {
		int32_t r;
		for (r = 0; r < sk[j].n; ++r) free(sk[j].a[r].acc);
		free(sk[j].a);
	}
	free(sk); free(span);
}

static void write_per_seq(step_t *t)
{
	const pipeline_t *p = t->p;
	int32_t i, j;
	kstring_t out = {0,0,0};
	void *km = p->opt->algo == RB3_SA_CHAIN? km_init() : 0; // arena reused by the PAV lift projections
	if (p->opt->algo == RB3_SA_CHAIN && p->opt->paired) // both mates first, then one event per pair
		chain_write_paired((pipeline_t*)p, t, km);
	for (j = 0; j < t->n_seq; ++j) {
		m_seq_t *s = &t->seq[j];
		if (p->opt->algo != RB3_SA_CHAIN) free(s->seq); // chain's PAV path still reads s->seq (pav_low_complexity)
		out.l = 0;
		if (p->opt->algo == RB3_SA_CHAIN) { // unite SMEMs -> per-read PS4G (per exon segment)
			if (!p->opt->paired) chain_emit(p, s, km, 0);
			free(s->seq);
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
			const char *pcls = t->pe? rb3_pe_class_str[(int)t->pe[j >> 1].cls] : 0;
			uint8_t *keep = t->sw_keep? t->sw_keep[j] : 0; // --paired CONC: concordant positions only
			if (r->n > 0) { // mapped
				int32_t x = 0;
				for (i = 0; i < r->n; ++i) {
					rb3_swhit_t *h = &r->a[i];
					if (keep) { // move the concordant positions to the front; drop the hit if it has none
						int32_t k, m = 0;
						for (k = 0; k < h->n_pos; ++k, ++x)
							if (keep[x]) h->pos[m++] = h->pos[k];
						if (m == 0) continue;
						h->n_pos = m;
					}
					out.l = 0;
					write_paf(&out, &p->fmi, h, s);
					if (pcls) out.l--, rb3_sprintf_lite(&out, "\tpe:Z:%s\n", pcls);
					fputs(out.s, stdout);
				}
			} else if (p->opt->flag & RB3_MF_WRITE_UNMAP) { // unmapped
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0", s->len);
				if (pcls) rb3_sprintf_lite(&out, "\tpe:Z:%s", pcls);
				rb3_sprintf_lite(&out, "\n");
				fputs(out.s, stdout);
			}
			if (keep) free(keep);
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
			const rb3_pe_res_t *e = t->pe? &t->pe[j >> 1] : 0; // --paired: this read's pair
			int filt = e && e->cls == RB3_PE_CONC;            // concordant pair: print only concordant positions
			for (i = 0; i < s->n_mem; ++i) {
				m_sai_pos_t *r = &s->mem[i];
				rb3_sai_t *q = &r->mem;
				int32_t st = q->info>>32, en = (int32_t)q->info;
				out.l = 0;
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t%d\t%ld", st, en, (long)q->size);
				if (r->n_pos > 0) {
					int32_t j, nk = r->n_pos;
					if (filt) for (j = 0, nk = 0; j < r->n_pos; ++j) nk += r->keep && r->keep[j];
					rb3_sprintf_lite(&out, "\t%ld", (long)nk);
					for (j = 0; j < r->n_pos; ++j) {
						rb3_pos_t *t = &r->pos[j];
						int64_t rlen = f->sid->len[t->sid>>1], pos;
						if (filt && !(r->keep && r->keep[j])) continue;
						pos = t->sid&1? rlen - (t->pos + (en - st)) : t->pos;
						rb3_sprintf_lite(&out, "\t%s:%c:%ld", f->sid->name[t->sid>>1], "+-"[t->sid&1], pos);
					}
					free(r->pos);
				}
				free(r->keep);
				if (e) rb3_sprintf_lite(&out, "\t%s", rb3_pe_class_str[(int)e->cls]);
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
	free(t->sw_keep);
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

/********************************************************************
 * --paired: mate pairing over one batch (mem)
 ********************************************************************/

// A read's located SMEM hits, as implied read starts on each contig's forward strand, flattened
// in (SMEM, position) order; full_only: only SMEMs covering the whole read. Returns the count.
static int32_t mem_read_hits(const pipeline_t *p, const m_seq_t *s, int full_only, rb3_pe_hit_t *h)
{
	const rb3_fmi_t *f = &p->fmi;
	int32_t i, j, n = 0;
	for (i = 0; i < s->n_mem; ++i) {
		const m_sai_pos_t *r = &s->mem[i];
		int32_t st = r->mem.info>>32, en = (int32_t)r->mem.info;
		if (full_only && (st != 0 || en != s->len)) continue;
		for (j = 0; j < r->n_pos; ++j) {
			const rb3_pos_t *t = &r->pos[j];
			int64_t rlen = f->sid->len[t->sid>>1];
			int strand = t->sid & 1;
			int64_t pos = strand? rlen - (t->pos + (en - st)) : t->pos; // forward start of the SMEM
			h[n].sid = t->sid >> 1, h[n].strand = strand;
			h[n++].start = strand? pos - (s->len - en) : pos - st;
		}
	}
	return n;
}

// 1 if every occurrence of the read's longest SMEM was located (-p did not truncate it). A mate
// whose hits are only a sample cannot veto its partner's hits: the missing occurrence may be the
// partner's true locus.
static int mem_fully_located(const m_seq_t *s)
{
	int32_t i, best = -1, bl = -1;
	for (i = 0; i < s->n_mem; ++i) {
		int32_t l = (int32_t)s->mem[i].mem.info - (int32_t)(s->mem[i].mem.info>>32);
		if (l > bl) bl = l, best = i;
	}
	return best >= 0 && s->mem[best].n_pos >= s->mem[best].mem.size;
}

static int32_t mem_n_pos(const m_seq_t *s)
{
	int32_t i, n = 0;
	for (i = 0; i < s->n_mem; ++i) n += s->mem[i].n_pos;
	return n;
}

static void mem_pair_worker(void *data, long i, int tid)
{
	step_t *t = (step_t*)data;
	m_seq_t *s[2];
	rb3_pe_hit_t *h[2];
	uint8_t *k[2];
	int32_t n[2], c, nc;
	rb3_pe_res_t *e = &t->pe[i];
	s[0] = &t->seq[2*i], s[1] = &t->seq[2*i+1];
	for (c = 0; c < 2; ++c) {
		int32_t m = mem_n_pos(s[c]);
		h[c] = RB3_MALLOC(rb3_pe_hit_t, m > 0? m : 1);
		k[c] = RB3_CALLOC(uint8_t, m > 0? m : 1);
		n[c] = mem_read_hits(t->p, s[c], 0, h[c]);
	}
	nc = rb3_pe_concordant_hits(h[0], n[0], h[1], n[1], t->pe_max_dist, s[0]->len, s[1]->len, k[0], k[1], 0);
	for (c = 0; c < 2 && nc > 0; ++c)
		if (!mem_fully_located(s[c^1])) // the mate's hits are a -p sample: keep all of this read's positions
			memset(k[c], 1, n[c] > 0? n[c] : 1);
	e->anchor = e->rescued = -1;
	e->cls = nc > 0? RB3_PE_CONC : n[0] > 0 && n[1] > 0? RB3_PE_DISC : n[0] > 0 || n[1] > 0? RB3_PE_SOLO : RB3_PE_NONE;
	e->resc_support = nc; // mem: number of contigs carrying the fragment
	if (nc > 0) // scatter the flags back to each SMEM's positions (same flattened order as mem_read_hits)
		for (c = 0; c < 2; ++c) {
			int32_t a, b, x = 0;
			for (a = 0; a < s[c]->n_mem; ++a) {
				m_sai_pos_t *r = &s[c]->mem[a];
				if (r->n_pos <= 0) continue;
				r->keep = RB3_MALLOC(uint8_t, r->n_pos);
				for (b = 0; b < r->n_pos; ++b) r->keep[b] = k[c][x++];
			}
		}
	for (c = 0; c < 2; ++c) free(h[c]), free(k[c]);
}

/********************************************************************
 * --paired: mate pairing over one batch (sw)
 ********************************************************************/
// A read's located sw hits as implied read starts on each contig's forward strand, flattened in
// (hit, position) order; full_only: only end-to-end (whole-read) hits. Returns the count.
static int32_t sw_read_hits(const pipeline_t *p, const m_seq_t *s, const rb3_swrst_t *r, int full_only, rb3_pe_hit_t *h)
{
	int32_t i, k, n = 0;
	for (i = 0; i < r->n; ++i) {
		const rb3_swhit_t *x = &r->a[i];
		int32_t qs = x->n_qoff > 0? x->qoff[0] : 0, qe = qs + x->qlen;
		if (full_only && (qs != 0 || qe != s->len)) continue;
		for (k = 0; k < x->n_pos; ++k) {
			int64_t clen, st, en;
			int strand = x->pos[k].sid & 1;
			pos_stranded(p->fmi.sid, &x->pos[k], x->rlen, &clen, &st, &en);
			h[n].sid = x->pos[k].sid >> 1, h[n].strand = strand;
			h[n++].start = strand? st - (s->len - qe) : st - qs;
		}
	}
	return n;
}

static int32_t sw_n_pos(const rb3_swrst_t *r)
{
	int32_t i, n = 0;
	for (i = 0; i < r->n; ++i) n += r->a[i].n_pos;
	return n;
}

// 1 if every occurrence of the read's best-scoring hit was located (-p did not truncate it)
static int sw_fully_located(const rb3_swrst_t *r)
{
	int32_t i, b = -1;
	for (i = 0; i < r->n; ++i) if (b < 0 || r->a[i].score > r->a[b].score) b = i;
	return b >= 0 && r->a[b].n_pos >= r->a[b].hi - r->a[b].lo;
}

static void sw_pair_worker(void *data, long i, int tid)
{
	step_t *t = (step_t*)data;
	rb3_pe_hit_t *h[2];
	uint8_t *k[2];
	int32_t n[2], c, nc;
	rb3_pe_res_t *e = &t->pe[i];
	for (c = 0; c < 2; ++c) {
		int32_t m = sw_n_pos(&t->rst[2*i+c]);
		h[c] = RB3_MALLOC(rb3_pe_hit_t, m > 0? m : 1);
		k[c] = RB3_CALLOC(uint8_t, m > 0? m : 1);
		n[c] = sw_read_hits(t->p, &t->seq[2*i+c], &t->rst[2*i+c], 0, h[c]);
	}
	nc = rb3_pe_concordant_hits(h[0], n[0], h[1], n[1], t->pe_max_dist, t->seq[2*i].len, t->seq[2*i+1].len, k[0], k[1], 0);
	e->anchor = e->rescued = -1;
	e->cls = nc > 0? RB3_PE_CONC : n[0] > 0 && n[1] > 0? RB3_PE_DISC : n[0] > 0 || n[1] > 0? RB3_PE_SOLO : RB3_PE_NONE;
	e->resc_support = nc;
	for (c = 0; c < 2 && nc > 0; ++c)
		if (!sw_fully_located(&t->rst[2*i+(c^1)])) // the mate's positions are a -p sample: keep all of this read's
			memset(k[c], 1, n[c] > 0? n[c] : 1);
	for (c = 0; c < 2; ++c) {
		if (nc > 0) t->sw_keep[2*i+c] = k[c];
		else free(k[c]);
		free(h[c]);
	}
}

static void sw_pair_batch(pipeline_t *p, step_t *t)
{
	int32_t i, n_pair = t->n_seq / 2, n_span = 0;
	t->pe = RB3_CALLOC(rb3_pe_res_t, n_pair > 0? n_pair : 1);
	t->sw_keep = RB3_CALLOC(uint8_t*, t->n_seq > 0? t->n_seq : 1);
	if (!p->isize_fixed) { // sample: pairs whose mates both have end-to-end hits on a shared contig
		int64_t *span = RB3_MALLOC(int64_t, n_pair > 0? n_pair : 1);
		for (i = 0; i < n_pair; ++i) {
			int32_t ma = sw_n_pos(&t->rst[2*i]), mb = sw_n_pos(&t->rst[2*i+1]), na, nb;
			rb3_pe_hit_t *ha, *hb;
			int64_t sp;
			if (ma == 0 || mb == 0) continue;
			ha = RB3_MALLOC(rb3_pe_hit_t, ma), hb = RB3_MALLOC(rb3_pe_hit_t, mb);
			na = sw_read_hits(p, &t->seq[2*i], &t->rst[2*i], 1, ha), nb = sw_read_hits(p, &t->seq[2*i+1], &t->rst[2*i+1], 1, hb);
			if (rb3_pe_concordant_hits(ha, na, hb, nb, RB3_PE_EST_WIN, t->seq[2*i].len, t->seq[2*i+1].len, 0, 0, &sp) > 0 && sp > 0)
				span[n_span++] = sp;
			free(ha), free(hb);
		}
		pe_update_isize(p, span, n_span);
		free(span);
	}
	t->pe_max_dist = (p->isize.ok? p->isize.max_ins : RB3_PE_DEFAULT_MAX_INS) + p->opt->pe_slack;
	kt_for(p->opt->n_threads, sw_pair_worker, t, n_pair);
	for (i = 0; i < n_pair; ++i) ++p->pe_cnt[(int)t->pe[i].cls];
}

static void mem_pair_batch(pipeline_t *p, step_t *t)
{
	int32_t i, n_pair = t->n_seq / 2, n_span = 0;
	t->pe = RB3_CALLOC(rb3_pe_res_t, n_pair > 0? n_pair : 1);
	if (!p->isize_fixed) { // sample: pairs whose mates both have whole-read SMEM hits on a shared contig
		int64_t *span = RB3_MALLOC(int64_t, n_pair > 0? n_pair : 1);
		for (i = 0; i < n_pair; ++i) {
			const m_seq_t *a = &t->seq[2*i], *b = &t->seq[2*i+1];
			int32_t ma = mem_n_pos(a), mb = mem_n_pos(b), na, nb;
			rb3_pe_hit_t *ha, *hb;
			int64_t sp;
			if (ma == 0 || mb == 0) continue;
			ha = RB3_MALLOC(rb3_pe_hit_t, ma), hb = RB3_MALLOC(rb3_pe_hit_t, mb);
			na = mem_read_hits(p, a, 1, ha), nb = mem_read_hits(p, b, 1, hb);
			if (rb3_pe_concordant_hits(ha, na, hb, nb, RB3_PE_EST_WIN, a->len, b->len, 0, 0, &sp) > 0 && sp > 0)
				span[n_span++] = sp;
			free(ha), free(hb);
		}
		pe_update_isize(p, span, n_span);
		free(span);
	}
	t->pe_max_dist = (p->isize.ok? p->isize.max_ins : RB3_PE_DEFAULT_MAX_INS) + p->opt->pe_slack;
	kt_for(p->opt->n_threads, mem_pair_worker, t, n_pair);
	for (i = 0; i < n_pair; ++i) ++p->pe_cnt[(int)t->pe[i].cls];
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
		while (!rb3_hitcount_reached(&p->hitcount)) { // read sequences; --paired: whole mate pairs only
			int32_t k, nm = p->opt->paired? 2 : 1;
			for (k = 0; k < nm; ++k) {
				m_seq_t *s;
				if ((ss = rb3_seq_read1(k == 1 && p->fp2? p->fp2 : p->fp, &len, &name)) == 0) break;
				RB3_GROW0(m_seq_t, seq, n_seq, m_seq);
				s = &seq[n_seq++];
				s->name = name? rb3_strdup(name) : 0;
				s->seq = (uint8_t*)rb3_strdup(ss);
				s->len = len;
				s->id = p->id++;
				s->mem = 0, s->n_mem = 0;
				tot += len;
			}
			if (k == 0) { // end of input
				if (p->fp2 && rb3_seq_read1(p->fp2, &len, &name) != 0) {
					fprintf(stderr, "ERROR: --paired: the second file has more reads than the first\n");
					exit(1);
				}
				break;
			}
			if (k < nm) {
				fprintf(stderr, "ERROR: --paired: %s\n", p->fp2? "the first file has more reads than the second"
						: "odd number of reads in the interleaved file");
				exit(1);
			}
			if (nm == 2 && !rb3_pe_names_match(seq[n_seq-2].name, seq[n_seq-1].name)) {
				fprintf(stderr, "ERROR: --paired: reads '%s' and '%s' are not mates (names must match up to /1 /2)\n",
						seq[n_seq-2].name? seq[n_seq-2].name : "?", seq[n_seq-1].name? seq[n_seq-1].name : "?");
				exit(1);
			}
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
		if (p->opt->algo == RB3_SA_REFMAP && (p->err_cut >= 0 || p->opt->max_ref_span > 0))
			refmap_filter_pre(p, t);
		if (p->opt->paired && p->opt->algo == RB3_SA_REFMAP) {
			refmap_pair_batch(p, t);
			if (p->opt->pe_drop_disc || p->opt->max_ref_span > 0) refmap_filter_pairs(p, t);
		} else if (p->opt->paired && (p->opt->algo == RB3_SA_MEM_TG || p->opt->algo == RB3_SA_MEM_ORI))
			mem_pair_batch(p, t);
		else if (p->opt->paired && p->opt->algo == RB3_SA_SW)
			sw_pair_batch(p, t);
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
			write_refmap(t), free(t->refmap), free(t->pe);
		else
			write_per_seq(t), free(t->pe);
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
	{ "lift-local-win",  ko_required_argument, 360 },
	{ "lift-local-tol",  ko_required_argument, 361 },
	{ "lift-chain-nb",   ko_required_argument, 362 },
	{ "lift-chain-interp", ko_no_argument,     363 },
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
	{ "report-occ",      ko_no_argument,       341 },
	{ "max-intron",      ko_required_argument, 330 },
	{ "gap-intron",      ko_required_argument, 331 },
	{ "chain-max-occ",   ko_required_argument, 332 },
	{ "ref-fasta",       ko_required_argument, 333 },
	{ "splice-min",      ko_required_argument, 334 },
	{ "pav-min-len",     ko_required_argument, 335 },
	{ "walk",            ko_no_argument,       336 },
	{ "trim-polya",      ko_required_argument, 337 },
	{ "pav-grid",        ko_required_argument, 338 },
	{ "pav-agree",       ko_required_argument, 340 },
	{ "diverged-rows",       ko_required_argument, 339 },
	{ "insertion-rows",      ko_required_argument, 342 },
	{ "anchor-dist-npy",     ko_no_argument,       343 },
	{ "anchor-dist-thresh",  ko_required_argument, 344 },
	{ "support",             ko_required_argument, 345 },
	{ "support-tol",         ko_required_argument, 346 },
	{ "support-win",         ko_required_argument, 347 },
	{ "dup-count",           ko_no_argument,       348 },
	{ "paired",              ko_no_argument,       370 },
	{ "pe-isize",            ko_required_argument, 371 },
	{ "pe-slack",            ko_required_argument, 372 },
	{ "pe-trust-exact",      ko_no_argument,       373 },
	{ "max-locate",          ko_required_argument, 374 },
	{ "pe-near-win",         ko_required_argument, 375 },
	{ "err-max-events",      ko_required_argument, 380 },
	{ "err-calib",           ko_required_argument, 381 },
	{ "err-alpha",           ko_required_argument, 382 },
	{ "err-hist",            ko_required_argument, 383 },
	{ "err-stats",           ko_no_argument,       384 },
	{ "max-ref-span",        ko_required_argument, 385 },
	{ "pe-drop-disc",        ko_no_argument,       386 },
	{ "no-kalloc",       ko_no_argument,       501 },
	{ "dbg-dawg",        ko_no_argument,       502 },
	{ "dbg-sw",          ko_no_argument,       503 },
	{ "dbg-qname",       ko_no_argument,       504 },
	{ "dbg-bt",          ko_no_argument,       505 },
	{ 0, 0, 0 }
};

// --anchor-dist-npy/--anchor-dist-thresh usage, shared verbatim by refmap and chain: both
// commands feed the same PS4G/npy accumulator and rb3_ps4g_npy_finalize builds the block.
static void usage_anchor_dist(const rb3_mopt_t *opt)
{
	fprintf(stderr, "  --anchor-dist-npy with --npy, also write a per-founder ternary read-sharing block\n");
	fprintf(stderr, "                    (match=1/diverged=0/deletion=-1) and a per-founder distance-to-\n");
	fprintf(stderr, "                    nearest-lift-anchor block (bp, -1 = no anchor), widening --npy from\n");
	fprintf(stderr, "                    (bin x gamete+2) to (bin x 3*gamete+2); requires --lift; a\n");
	fprintf(stderr, "                    '<npy>.layout.tsv' sidecar records the column layout\n");
	fprintf(stderr, "  --anchor-dist-thresh=NUM  with --anchor-dist-npy, distance (bp) at/under which a\n");
	fprintf(stderr, "                    non-matching founder reads as diverged rather than deletion [%ld]\n", (long)opt->anchor_dist_thresh);
}

int main_search(int argc, char *argv[]) // "sw" and "mem" share the same CLI
{
	int32_t c, j, is_line = 0, ret, load_flag = 0, no_ssa = 0, ref_gamete = -1; // --anchor-dist-npy: gamete index of the reference genome, -1 if unknown
	rb3_mopt_t opt;
	pipeline_t p;
	ketopt_t o = KETOPT_INIT;
	int64_t lift_local_win_opt = 0, lift_local_tol_opt = 1000; // --lift-local-win/--lift-local-tol
	int32_t lift_chain_nb_opt = 0, lift_chain_interp_opt = 0; // --lift-chain-nb, --lift-chain-interp

	rb3_mopt_init(&opt);
	if (getenv("RB3_PLACE_DEBUG")) refmap_place_debug = atoi(getenv("RB3_PLACE_DEBUG"));
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
		else if (c == 360) lift_local_win_opt = rb3_parse_num(o.arg);
		else if (c == 361) lift_local_tol_opt = rb3_parse_num(o.arg);
		else if (c == 362) lift_chain_nb_opt = atoi(o.arg);
		else if (c == 363) lift_chain_interp_opt = 1;
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
		else if (c == 341) opt.report_occ = 1; // append raw FM-index interval size (occurrence count) as an extra column
		else if (c == 330) opt.max_intron = rb3_parse_num(o.arg);    // chain: max unexplained ref jump per link
		else if (c == 331) opt.gap_intron = atoi(o.arg);             // chain: ref gap that starts a new exon segment
		else if (c == 332) opt.chain_max_occ = atoi(o.arg);          // chain: interval-size cap for an informative SMEM
		else if (c == 333) opt.ref_fasta = o.arg;                    // chain: reference FASTA -> GT-AG splice check
		else if (c == 334) opt.splice_min = atoi(o.arg);             // chain: min ref gap for the GT-AG check
		else if (c == 335) opt.pav_min_len = atoi(o.arg);            // chain --lift: assembly-only specificity floor
		else if (c == 336) opt.allow_walk = 1;                       // refmap: opt in to deprecated walking
		else if (c == 337) opt.trim_polya = atoi(o.arg);             // trim terminal poly-A/T runs >= INT bp
		else if (c == 338) opt.pav_grid = rb3_parse_num(o.arg);      // chain --lift: pav position grid
		else if (c == 340) opt.pav_agree = rb3_parse_num(o.arg);     // chain --lift: projection agreement
		else if (c == 339) {                                         // chain --lift: colinear-row policy
			if (strcmp(o.arg, "ordinary") == 0) opt.diverged_rows = 0;
			else if (strcmp(o.arg, "pav") == 0) opt.diverged_rows = 1;
			else if (strcmp(o.arg, "drop") == 0) opt.diverged_rows = 2;
			else { fprintf(stderr, "ERROR: --diverged-rows must be ordinary|pav|drop\n"); return 1; }
		}
		else if (c == 342) {                                         // chain --lift: insertion-row policy
			if (strcmp(o.arg, "ordinary") == 0) opt.insertion_rows = 0;
			else if (strcmp(o.arg, "pav") == 0) opt.insertion_rows = 1;
			else if (strcmp(o.arg, "drop") == 0) opt.insertion_rows = 2;
			else { fprintf(stderr, "ERROR: --insertion-rows must be ordinary|pav|drop\n"); return 1; }
		}
		else if (c == 343) opt.anchor_dist_npy = 1;                  // npy: also write ternary sharing + lift-anchor distance blocks
		else if (c == 344) opt.anchor_dist_thresh = rb3_parse_num(o.arg); // diverged/deletion distance threshold (bp)
		else if (c == 345) {                                         // chain: founder-set rule
			if (strcmp(o.arg, "chain") == 0) opt.support_best = 0;
			else if (strcmp(o.arg, "best") == 0) opt.support_best = 1;
			else { fprintf(stderr, "ERROR: --support must be chain or best\n"); return 1; }
		}
		else if (c == 346) opt.support_tol = atoi(o.arg);            // chain --support=best: coverage tolerance (bases)
		else if (c == 347) opt.support_win = rb3_parse_num(o.arg);   // chain --support=best: lift proximity window
		else if (c == 348) opt.dup_count = 1;                        // chain: "<name>;x=N" reads count N times
		else if (c == 370) opt.paired = 1;                           // mate pairs: R1 R2, or one interleaved file
		else if (c == 371) {                                         // --paired: fixed insert size MEAN,SD
			if (sscanf(o.arg, "%lf,%lf", &opt.pe_mean, &opt.pe_sd) != 2 || opt.pe_mean <= 0 || opt.pe_sd < 0) {
				fprintf(stderr, "ERROR: --pe-isize must be MEAN,SD (MEAN > 0, SD >= 0)\n");
				return 1;
			}
		}
		else if (c == 372) opt.pe_slack = rb3_parse_num(o.arg), opt.pe_slack_set = 1; // --paired: discordance slack (bp)
		else if (c == 375) opt.pe_near_win = rb3_parse_num(o.arg);   // --paired: rescued-mate founder window
		else if (c == 380) opt.err_max_ev = strcmp(o.arg, "auto") == 0? -2 : atoi(o.arg); // refmap: events cutoff
		else if (c == 381) opt.err_calib = rb3_parse_num(o.arg);     // refmap: calibration reads (pairs)
		else if (c == 382) opt.err_alpha = atof(o.arg);              // refmap: auto-cutoff tail probability
		else if (c == 383) opt.err_hist_fn = o.arg;                  // refmap: calibration histogram sidecar
		else if (c == 384) opt.err_stats = 1;                        // refmap: events/span/reason columns
		else if (c == 385) opt.max_ref_span = rb3_parse_num(o.arg);  // refmap: reference span cap (bp)
		else if (c == 386) opt.pe_drop_disc = 1;                     // refmap --paired: drop DISC pairs
		else if (c == 374) opt.max_locate = atoi(o.arg);             // refmap: occurrences located per read
		else if (c == 373) opt.pe_trust_exact = 1;                   // --paired: move unrescuable PLACED mates onto EXACT ones
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
			fprintf(stderr, "  --paired    reads are mate pairs (<R1> <R2>, or one interleaved file); needs -p. A pair\n");
			fprintf(stderr, "              with hits on the same contig, opposite strands, implied starts within the\n");
			fprintf(stderr, "              insert is CONC and prints only those positions; a trailing column gives\n");
			fprintf(stderr, "              the pair class (CONC/DISC/SOLO/NONE). --pe-isize=MEAN,SD, --pe-slack=NUM [0]\n");
			fprintf(stderr, "              With --paired use a LARGE -p (start at -p 1024): a small cap locates only a\n");
			fprintf(stderr, "              sample of a repetitive SMEM's hits and misses the true locus (maize: -p 128\n");
			fprintf(stderr, "              91.4%% vs -p 1024 96.2%% true-founder recall; pairing keeps the sets small)\n");
		}
		if (strcmp(argv[0], "refmap") == 0) {
			fprintf(stderr, "  --ref-prefix=STR  reference = sequences whose name starts with STR [required]\n");
			fprintf(stderr, "  --lift=FILE       REQUIRED (standard): project assembly hits via a `ropebwt3 lift` map\n");
			fprintf(stderr, "  --walk            opt in to the DEPRECATED flank-walking path instead of --lift\n");
			fprintf(stderr, "  --max-walk=NUM    DEPRECATED (--walk only) max bases to walk per flank [%d]\n", opt.max_walk);
			fprintf(stderr, "  --walk-mode=STR   DEPRECATED (--walk only) consensus|strict|per-assembly [consensus]\n");
			fprintf(stderr, "  --max-occ=INT     drop reads/anchors occurring >INT times; <0 = auto (#taxa); 0 = off [%ld]\n", (long)opt.max_occ);
			fprintf(stderr, "  --max-locate=INT  locate at most INT occurrences of each read; the rest are not\n");
			fprintf(stderr, "                    sampled, so a repetitive read's gameteSet / reference copies may be\n");
			fprintf(stderr, "                    incomplete, and a read whose reference copy is not sampled is called\n");
			fprintf(stderr, "                    PLACED instead of EXACT (refmap's analogue of mem -p) [%d]\n", opt.max_locate);
			fprintf(stderr, "  --two-flank       require both flanks to anchor concordantly (drop ONE_SIDE)\n");
			fprintf(stderr, "  --max-bracket=NUM with --two-flank, max |cR-cL| for a PLACED; 0 = off [%ld]\n", (long)opt.max_bracket);
			fprintf(stderr, "  --lift-win=NUM    liftover projection window [%ld]\n", (long)opt.lift_win);
			fprintf(stderr, "  --lift-mad=NUM    liftover max residual MAD [%ld]\n", (long)opt.lift_mad);
			fprintf(stderr, "  --lift-local-win=NUM  project via nearest collinear flanking anchors within NUM bp first,\n"
			                "                    falling back to the --lift-win median (0 = off) [0]\n");
			fprintf(stderr, "  --lift-local-tol=NUM  collinearity tolerance for --lift-local-win [1000]\n");
			fprintf(stderr, "  --lift-chain-nb=INT   with --lift-local-win: drop anchors not collinear with one of their next\n"
			                "                    INT anchors, project via the nearest kept anchors (0 = off) [0]\n");
			fprintf(stderr, "  --lift-chain-interp   with --lift-chain-nb: interpolate across indels between the kept\n"
			                "                    anchors instead of extrapolating from the nearer one\n");
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
			usage_anchor_dist(&opt);
			fprintf(stderr, "  --paired          reads are mate pairs: <R1.fq> <R2.fq>, or one interleaved file. Each\n");
			fprintf(stderr, "                    pair adds ONE PS4G/npy event (the founders both mates share) at the\n");
			fprintf(stderr, "                    better-placed mate; a projected mate far from its mate is moved next\n");
			fprintf(stderr, "                    to it when one of its assembly projections lands there; stdout gains\n");
			fprintf(stderr, "                    a trailing pair-class column (NONE/SOLO/CONC/RESC/MOVED/DISC).\n");
			fprintf(stderr, "                    -p does not apply to refmap; raise --max-locate instead\n");
			fprintf(stderr, "  --pe-isize=MEAN,SD  --paired: fixed insert size [estimated per batch]\n");
			fprintf(stderr, "  --pe-slack=NUM    --paired: bp added to the max insert before mates are discordant [%ld]\n", (long)opt.pe_slack);
			fprintf(stderr, "  --pe-near-win=NUM --paired: a MOVED mate keeps the founders whose own projection lies\n");
			fprintf(stderr, "                    on its mate's contig within NUM bp (founder set only; the mate gap\n");
			fprintf(stderr, "                    is max insert + --pe-slack); 0 = that gap [%ld]\n", (long)opt.pe_near_win);
			fprintf(stderr, "  --pe-trust-exact  --paired: a discordant EXACT + PLACED pair with no rescue moves the\n");
			fprintf(stderr, "                    PLACED mate onto the EXACT one (default: keep them as two events)\n");
			fprintf(stderr, "  Error-profile / span filters (off by default; a rejected read, or both mates of a rejected\n");
			fprintf(stderr, "  pair, is reported as FILTERED and adds nothing to --ps4g/--npy/--ps4g-per-read):\n");
			fprintf(stderr, "  --err-max-events=INT|auto  max mismatch events (edits explaining the read against one\n");
			fprintf(stderr, "                    carrier of its exact core) per read, or per pair (sum) with --paired;\n");
			fprintf(stderr, "                    auto = calibrate on the first --err-calib reads/pairs (off)\n");
			fprintf(stderr, "  --err-calib=NUM   reads (pairs with --paired) in the calibration pre-pass [%ld]\n", (long)opt.err_calib);
			fprintf(stderr, "  --err-alpha=FLOAT auto: smallest K with P(Poisson(mu_err) > K) <= FLOAT [%g]\n", opt.err_alpha);
			fprintf(stderr, "  --err-hist=FILE   write the calibration histogram and error-model fit (TSV)\n");
			fprintf(stderr, "  --err-stats       append events, reference span and filter reason columns to each row\n");
			fprintf(stderr, "  --max-ref-span=NUM  max reference span (bp): EXACT first-to-last reference copy, PLACED\n");
			fprintf(stderr, "                    spread of all lift projections (another sequence = infinite); with\n");
			fprintf(stderr, "                    --paired, first-to-last coordinate of both mates (0 = off)\n");
			fprintf(stderr, "  --pe-drop-disc    --paired: FILTER discordant pairs instead of two single-end events\n");
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
			fprintf(stderr, "  --insertion-rows=STR  TRUE INSERTION rows -- same routing as --diverged-rows, but\n");
			fprintf(stderr, "                    for reads absent from the reference: ordinary|pav|drop\n");
			fprintf(stderr, "                    [pav = grid-snapped, prefixed `pav:<contig>` (default, unchanged)]\n");
			fprintf(stderr, "                    ordinary emits the read inline on its real contig at the exact\n");
			fprintf(stderr, "                    (unsnapped) median -- more precise than pav-grid snapping (which\n");
			fprintf(stderr, "                    floors, not centers) and keeps evidence positionally inline with\n");
			fprintf(stderr, "                    ordinary rows for windowing, at the cost of the `pav:` marker\n");
			fprintf(stderr, "  --pav-grid=NUM    snap the emitted pav: position to this grid; a breakpoint is\n");
			fprintf(stderr, "                    approximate, so this is the resolution actually claimed\n");
			fprintf(stderr, "                    (0 = emit the exact projected position) [%d]\n", opt.pav_grid);
			fprintf(stderr, "  --pav-agree=NUM   max spread among the assemblies' projections of the seed before\n");
			fprintf(stderr, "                    the read is suppressed as unplaceable [%d]\n", opt.pav_agree);
			fprintf(stderr, "  --pav-min-len=INT min LONGEST assembly-only SMEM to emit a pav: row; the\n");
			fprintf(stderr, "                    assembly-only path has no colinear reference confirmation, so it\n");
			fprintf(stderr, "                    needs a longer floor than -l [%d]\n", opt.pav_min_len);
			fprintf(stderr, "  --support=STR     founder set of a chained read: chain = intersection over the\n");
			fprintf(stderr, "                    reference-hitting chain's SMEMs (default, unchanged); best = the\n");
			fprintf(stderr, "                    assemblies covering the most read bases, counting the chain's SMEMs\n");
			fprintf(stderr, "                    plus assembly-only SMEMs (for diverged / out-of-panel alleles)\n");
			fprintf(stderr, "  --support-tol=INT --support=best: keep assemblies within INT covered bases of the best [%d]\n", opt.support_tol);
			fprintf(stderr, "  --support-win=NUM --support=best --lift: an assembly-only SMEM counts only if it\n");
			fprintf(stderr, "                    projects within NUM bp of the chain's reference span (0 = off) [%ld]\n", (long)opt.support_win);
			fprintf(stderr, "  --dup-count       a read named <name>;x=N stands for N identical reads: it adds N\n");
			fprintf(stderr, "                    to the PS4G/npy counts (collapse duplicates first, keep the counts)\n");
			fprintf(stderr, "  --trim-polya=INT  trim a terminal poly-A (3') / poly-T (5') run of >=INT bp before\n");
			fprintf(stderr, "                    searching. MEASURED AS A NO-OP here (SMEMs already isolate the\n");
			fprintf(stderr, "                    genomic core); low-complexity anchors are gated instead [%d]\n", opt.trim_polya);
			fprintf(stderr, "  -l INT      min SMEM length [%ld]\n", (long)opt.min_len);
			fprintf(stderr, "  --ps4g=FILE       write PS4G v2.0 gamete-support counts, aggregated over all rows\n");
			fprintf(stderr, "  --npy=FILE        write a dense (bin x gamete+2) numpy training/inference array\n");
			fprintf(stderr, "  --label-bed=FILE  diploid training labels: chrom start end sampleA [sampleB]\n");
			fprintf(stderr, "  --bin-size=NUM    PS4G/npy reference position bin size in bp [%ld]\n", (long)opt.bin_size);
			fprintf(stderr, "  --npy-binary      npy: write presence (1) instead of read counts\n");
			fprintf(stderr, "  --lift-local-win=NUM  --lift projection: nearest collinear flanking anchors within\n"
			                "                    NUM bp first (0 = off) [0]; --lift-local-tol=NUM tolerance [1000]\n");
			usage_anchor_dist(&opt);
			fprintf(stderr, "                    chain: a row's founder set is the read's chain-wide intersection;\n");
			fprintf(stderr, "                    a spliced read gives one row per exon segment, each judged at\n");
			fprintf(stderr, "                    its own start; use --insertion-rows=ordinary so insertion rows\n");
			fprintf(stderr, "                    are not --pav-grid snapped before the anchor-distance lookup\n");
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
			fprintf(stderr, "  --paired    mate pairs (<R1> <R2>, or one interleaved file); needs a LARGE -p (start at\n");
			fprintf(stderr, "              -p 1024: sw's -p is a per-read position budget, 0 locates one per hit). A pair\n");
			fprintf(stderr, "              with hits on the same contig, opposite strands, implied starts within the\n");
			fprintf(stderr, "              insert is CONC: its PAF lines keep only those positions. Every line gets a\n");
			fprintf(stderr, "              pe:Z: tag (CONC/DISC/SOLO/NONE). --pe-isize=MEAN,SD, --pe-slack=NUM [0]\n");
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

	if (opt.max_locate < 1) {
		if (rb3_verbose >= 1) fprintf(stderr, "ERROR: --max-locate must be >= 1\n");
		return 1;
	}
	if (opt.max_locate != RB3_RM_MAX_ASSEMBLY && opt.algo != RB3_SA_REFMAP) {
		if (rb3_verbose >= 1) fprintf(stderr, "ERROR: --max-locate applies to refmap only (mem/chain: use -p)\n");
		return 1;
	}
	if (opt.paired) { // check the paired-end setup before the (slow) index load
		int n_files = argc - o.ind - 1;
		const char *err = 0;
		int is_mem = (opt.algo == RB3_SA_MEM_TG || opt.algo == RB3_SA_MEM_ORI);
		int is_chain = (opt.algo == RB3_SA_CHAIN), is_sw = (opt.algo == RB3_SA_SW);
		if (opt.algo != RB3_SA_REFMAP && !is_mem && !is_chain && !is_sw) err = "--paired is implemented for refmap, mem, chain and sw";
		else if (is_sw && opt.max_pos <= 0) err = "sw --paired needs -p INT (positions to locate); use a large value, e.g. -p 1024";
		else if (is_sw && (opt.flag & RB3_MF_WRITE_ALL)) err = "sw --paired cannot be combined with --all-e2e/-g";
		else if (is_mem && opt.max_pos <= 0) err = "mem --paired needs -p INT (positions to locate per SMEM); use a large value, e.g. -p 1024";
		else if (is_mem && (opt.min_gap_len > 0 || (opt.flag & RB3_MF_WRITE_COV))) err = "mem --paired cannot be combined with --gap/--cov";
		else if ((is_mem || is_chain || is_sw) && opt.pe_trust_exact) err = "--pe-trust-exact applies to refmap only";
		else if (opt.walk_mode == RB3_WALK_PERASSEMBLY) err = "--paired does not support --walk-mode=per-assembly";
		else if (is_line) err = "--paired needs read names (FASTA/FASTQ), not -L";
		else if (n_files != 1 && n_files != 2) err = "--paired takes <R1.fq> <R2.fq> or one interleaved file";
		if (err) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: %s\n", err);
			return 1;
		}
		if ((is_mem || is_chain || is_sw) && !opt.pe_slack_set) opt.pe_slack = 0; // SMEM/chain/sw positions are exact
	} else if (opt.pe_mean > 0 || opt.pe_trust_exact || opt.pe_drop_disc) {
		if (rb3_verbose >= 1) fprintf(stderr, "ERROR: --pe-isize/--pe-trust-exact/--pe-drop-disc need --paired\n");
		return 1;
	}
	if (opt.err_max_ev != -1 || opt.err_hist_fn || opt.err_stats || opt.max_ref_span > 0 || opt.pe_drop_disc) {
		const char *err = 0;
		if (opt.algo != RB3_SA_REFMAP) err = "--err-*/--max-ref-span/--pe-drop-disc apply to refmap only";
		else if (opt.err_max_ev < -2) err = "--err-max-events takes an integer >= 0 or 'auto'";
		else if (opt.kmer_len > 0) err = "--err-*/--max-ref-span do not support --kmer";
		else if (opt.walk_mode == RB3_WALK_PERASSEMBLY) err = "--err-*/--max-ref-span do not support --walk-mode=per-assembly";
		else if (opt.err_calib < 1 || opt.err_alpha <= 0.0 || opt.err_alpha >= 1.0) err = "--err-calib must be >= 1 and 0 < --err-alpha < 1";
		else if ((opt.err_max_ev == -2 || opt.err_hist_fn) && (argc - o.ind < 2 || strcmp(argv[o.ind + 1], "-") == 0))
			err = "--err-max-events=auto/--err-hist read the input twice: give a file, not stdin";
		if (err) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: %s\n", err);
			return 1;
		}
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
	p.is_ref = 0, p.n_ref = 0, p.lift = 0, p.lift_ridx = 0, p.ref_seq = 0;
	p.gtab = 0, p.ps4g_acc = 0, p.ps4g_per_read_fp = 0, p.label_bed = 0;
	p.fp2 = 0, p.isize_fixed = 0, p.pe_warned = 0;
	memset(&p.isize, 0, sizeof(p.isize)), memset(p.pe_cnt, 0, sizeof(p.pe_cnt));
	if (opt.pe_mean > 0) { // --pe-isize: same max-insert rule as the estimator (mean + 4 sd)
		p.isize.ok = 1, p.isize_fixed = 1, p.isize.mean = opt.pe_mean, p.isize.sd = opt.pe_sd;
		p.isize.max_ins = (int64_t)(opt.pe_mean + 4.0 * opt.pe_sd + .5);
	}
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
			rb3_lift_set_local(lift_local_win_opt, lift_local_tol_opt);
			rb3_lift_set_chain(lift_chain_nb_opt, lift_chain_interp_opt);
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
			if (opt.npy_fn && opt.anchor_dist_npy) {
				int64_t k;
				if (p.lift == 0) {
					if (rb3_verbose >= 1) fprintf(stderr, "ERROR: --anchor-dist-npy requires --lift=FILE\n");
					free(p.is_ref);
					return 1;
				}
				for (k = 0; k < p.fmi.sid->n_seq; ++k) // any reference sequence's gamete is THE reference gamete
					if (p.is_ref[k]) { ref_gamete = p.gtab->sid2g[k]; break; }
				p.lift_ridx = rb3_lift_ridx_build(p.lift, p.gtab->sid2g, p.gtab->n_gamete);
				if (rb3_verbose >= 3)
					fprintf(stderr, "[M::%s] built lift reference-position index for --anchor-dist-npy (ref_gamete=%d)\n", __func__, ref_gamete);
				// chain shares this block unchanged: chain_emit/chain_emit_pav feed the same accumulator
				// (one event per emitted row: reference sequence, row start, the read's founder set) and
				// rb3_ps4g_npy_finalize derives the ternary/distance columns per (contig,bin,gameteSet)
				// row exactly as for refmap. So for a chained read: 1 = founder carries EVERY SMEM of the
				// winning chain (the chain-wide intersection, the same set stdout/PS4G report); each exon
				// segment is its own row, judged against that founder's lift anchors at the segment's own
				// start bin. Grid-snapped insertion rows (--insertion-rows=pav, the chain default) are
				// looked up at the snapped coordinate, up to --pav-grid bp from the real breakpoint.
				if (opt.algo == RB3_SA_CHAIN && opt.insertion_rows == 1 && opt.pav_grid > 0 && rb3_verbose >= 2)
					fprintf(stderr, "WARNING: chain --anchor-dist-npy with --insertion-rows=pav: insertion rows are snapped to "
									"--pav-grid=%d before the anchor-distance lookup; --insertion-rows=ordinary keeps them exact\n", opt.pav_grid);
			}
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
		if (rb3_verbose >= 1) fprintf(stderr, "ERROR: --ps4g/--npy/--label-bed/--target-hits only apply to refmap/chain\n");
		return 1;
	}
	if (opt.algo == RB3_SA_CHAIN) puts("readName\trefContig\trefPos\tgameteSet");
	if (opt.flag & RB3_MF_WRITE_ALL) {
		puts("CC\tQS  queryName  queryLen  numHap");
		puts("CC\tQH  refCount   score     editDist   cs   strand   nOut   totAln");
		puts("CC");
	}
	if (opt.algo == RB3_SA_REFMAP && opt.max_locate != RB3_RM_MAX_ASSEMBLY && opt.max_occ > 0
		&& opt.max_occ <= opt.max_locate && opt.kmer_len <= 0 && rb3_verbose >= 2)
		fprintf(stderr, "WARNING: --max-locate=%d has no effect: --max-occ=%ld already drops (MULTI) every read with more "
						"occurrences, so all occurrences of the reads kept are located\n", opt.max_locate, (long)opt.max_occ);
	p.err_on = (opt.err_max_ev != -1 || opt.err_hist_fn || opt.err_stats);
	p.err_cut = opt.err_max_ev >= 0? opt.err_max_ev : -1;
	memset(p.err_nfilt, 0, sizeof(p.err_nfilt));
	if (opt.algo == RB3_SA_REFMAP && (opt.err_max_ev == -2 || opt.err_hist_fn)) // calibration pre-pass
		if (refmap_err_calibrate(&p, argv[o.ind + 1], opt.paired && argc - o.ind - 1 == 2? argv[o.ind + 2] : 0, is_line) < 0)
			return 1;
	if (opt.paired) { // one pass over R1 (+ R2, or interleaved)
		p.fp = rb3_seq_open(argv[o.ind + 1], is_line);
		if (p.fp && argc - o.ind - 1 == 2) p.fp2 = rb3_seq_open(argv[o.ind + 2], is_line);
		if (p.fp == 0 || (argc - o.ind - 1 == 2 && p.fp2 == 0)) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to load the sequence file '%s'\n", p.fp == 0? argv[o.ind + 1] : argv[o.ind + 2]);
		} else kt_pipeline(2, worker_pipeline, &p, 3);
		if (p.fp) rb3_seq_close(p.fp);
		if (p.fp2) rb3_seq_close(p.fp2);
		if (rb3_verbose >= 3) {
			int64_t tot = 0;
			for (j = 0; j < 5; ++j) tot += p.pe_cnt[j];
			fprintf(stderr, "[M::%s] %ld pairs: %ld CONC, %ld RESC, %ld SOLO, %ld DISC, %ld NONE; insert size %s mean %.1f sd %.1f max %ld (+%ld slack)\n",
					__func__, (long)tot, (long)p.pe_cnt[RB3_PE_CONC], (long)p.pe_cnt[RB3_PE_RESC], (long)p.pe_cnt[RB3_PE_SOLO],
					(long)p.pe_cnt[RB3_PE_DISC], (long)p.pe_cnt[RB3_PE_NONE], p.isize_fixed? "(fixed)" : "(estimated)",
					p.isize.mean, p.isize.sd, (long)(p.isize.ok? p.isize.max_ins : RB3_PE_DEFAULT_MAX_INS), (long)opt.pe_slack);
		}
	}
	for (j = o.ind + 1; !opt.paired && j < argc; ++j) {
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
	if (opt.algo == RB3_SA_REFMAP && (p.err_cut >= 0 || opt.max_ref_span > 0 || opt.pe_drop_disc) && rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] FILTERED %s: %ld by events (> %d), %ld by span (> %ld bp), %ld by --pe-drop-disc\n", __func__,
				opt.paired? "mates" : "reads", (long)p.err_nfilt[0], p.err_cut, (long)p.err_nfilt[1], (long)opt.max_ref_span, (long)p.err_nfilt[2]);
	if (p.ps4g_acc) {
		kstring_t cmd = {0,0,0};
		int32_t k;
		rb3_sprintf_lite(&cmd, "ropebwt3");
		for (k = 0; k < argc; ++k) rb3_sprintf_lite(&cmd, " %s", argv[k]);
		rb3_ps4g_npy_finalize(p.ps4g_acc, p.gtab, p.fmi.sid, p.label_bed, opt.npy_binary,
							  p.lift_ridx, opt.anchor_dist_thresh, ref_gamete, opt.ps4g_fn, opt.npy_fn, cmd.s);
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
	rb3_lift_ridx_destroy(p.lift_ridx);
	rb3_lift_destroy(p.lift);
	return 0;
}
