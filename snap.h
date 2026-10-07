#ifndef RB3_SNAP_H
#define RB3_SNAP_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Base-exact refinement ("snap") of a projected interval end: align the reference
 * flank at that end to a window of the founder around the projected position.
 *
 * The aligner is a plain Gotoh affine-gap DP, global on the query (the reference
 * flank: every flank base is either aligned or deleted) and free at both ends of
 * the target (the founder window). Scores: +match, -mismatch, a gap of length L
 * costs gap_open + L*gap_ext (bwa-mem style). Sequences are nt6 codes (1..4 =
 * ACGT, 5 = N; N never matches). ropebwt3 ships no pairwise aligner (bwa-sw.c
 * aligns against the FM-index itself), hence this small one: the matrices are at
 * most F x (2R+2F) cells, e.g. 64 x 12k for the wide retry.
 *
 * Mapping to the truth convention (simulate_pangenome.py): the START of [s,e) is
 * the founder base aligned to the first flank base that is not a gap column (the
 * first reference base >= s present in the founder), the END is one past the
 * founder base aligned to the last such flank base. Founder-only bases next to the
 * end are never included (free target ends), and traceback prefers the diagonal,
 * so an equal-base insertion adjacent to the start maps s to its own (rightmost)
 * copy and one adjacent to the end is left out (the leftmost end is chosen).
 */

typedef struct {
	int32_t match, mismatch, gap_open, gap_ext; // positive costs, e.g. 1, 4, 6, 1
	int32_t min_aligned;  // minimum number of aligned (M) flank bases to accept
	double min_frac;      // ... and at least this fraction of the query aligned
	double min_id;        // minimum identity over the aligned core (first to last M column)
	int32_t margin;       // ambiguous if a placement >= excl bp away scores >= best - margin
	int32_t excl;         // exclusion zone around the best end position for the second best
	int32_t clip5;        // >= 0: the query may instead start anywhere, its skipped prefix
	                      // costing clip5 flat (a known-deleted stretch); < 0: global query
	int32_t no_seed;      // 1: always the full DP (no k-mer prefilter on wide windows; for testing)
	int32_t split;        // >= 0: also report the aligned target base of the last query base
	                      // before split and of the first from split on (t_left, t_right)
} rb3_snap_opt_t;

#define RB3_SNAP_OK      0
#define RB3_SNAP_AMBIG   1
#define RB3_SNAP_LOWID   2
#define RB3_SNAP_SHORT   3  // query (flank) or target (window) too short

typedef struct {
	int32_t score, second;     // best score; best score outside the exclusion zone (INT32_MIN if none)
	int32_t t_end;             // target end (exclusive) of the best alignment
	int32_t q_first, t_first;  // first aligned (M) column: query index, target index (q_first > 0
	                           // = that many leading query bases deleted, or clipped with clip5)
	int32_t q_last, t_last;    // last aligned (M) column
	int32_t t_left, t_right;   // with opt.split: target index of the last M column with query
	                           // index < split / the first with query index >= split (-1 if none)
	int32_t q_left, q_right;   // ... and their query indices
	int32_t n_aligned, n_match, n_core; // M columns, identical M columns, columns in [first M, last M]
	double identity;           // n_match / n_core
	int status;                // RB3_SNAP_*
} rb3_snap_res_t;

void rb3_snap_opt_init(rb3_snap_opt_t *o);

// Align q[0,n) (global) to t[0,m) (free ends). Fills r and returns r->status.
int rb3_snap_align(const uint8_t *q, int32_t n, const uint8_t *t, int32_t m, const rb3_snap_opt_t *o, rb3_snap_res_t *r);

#ifdef __cplusplus
}
#endif

#endif
