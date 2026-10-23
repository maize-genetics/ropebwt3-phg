#ifndef RB3_EDSTATS_H
#define RB3_EDSTATS_H

// Per-read edit-distance statistics against each assembly of the pangenome (refmap --ed-stats /
// --ed-hist). This file holds the pure part -- summaries of a per-assembly vector, the histogram
// accumulator and its writer, the npy writer -- so it can be unit tested without an index. The
// search itself (bounded-edit FM-index extension from exact read tiles) is refmap_ed in search.c.
//
// A per-assembly vector v[0..n_asm) holds, for one read (or one pair), the best edit distance
// found to each assembly: 0..cap, or RB3_ED_NA (none within the cap / not reachable from the
// read's seeds). Pairs: v = sum of the mates' values (NA if either mate is NA), cap = 2*cap.

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RB3_ED_NA      255
#define RB3_ED_MAX_CAP 60  // --ed-max upper limit; pair sums stay < RB3_ED_NA

// best = min(gbest, min over v) (NA if all NA; gbest = a best known without its assembly, or NA);
// n_best = assemblies at best; wrong = the closest DISTINCT wrong copy: min over assemblies with
// v > best + 1 (NA if none). Returns best.
int32_t rb3_ed_summary(const uint8_t *v, int32_t n, int32_t gbest, int32_t *n_best, int32_t *wrong);

// out[a] = a[a] + b[a] (NA if either is NA)
void rb3_ed_pair(const uint8_t *a, const uint8_t *b, int32_t n, uint8_t *out);

#define RB3_ED_UNIT_READ 0
#define RB3_ED_UNIT_PAIR 1

typedef struct {
	int32_t n_asm, cap[2];    // cap per unit: read cap, pair cap (2 * read cap)
	int64_t n[2];             // reads / pairs added
	int64_t *best[2], *asm_[2], *wrong[2]; // [cap+2]: 0..cap, cap+1 = NA
	int64_t *nbest[2];        // [n_asm+1]: assemblies at the best value (0 = no hit)
	int64_t n_trunc[2];       // reads whose search hit the budget (values above the reached cost are NA)
	int64_t n_rep[2];         // reads with at least one tile skipped as repetitive (--ed-max-occ)
} rb3_edhist_t;

rb3_edhist_t *rb3_edhist_init(int32_t n_asm, int32_t cap);
void rb3_edhist_destroy(rb3_edhist_t *h);
#define RB3_ED_F_TRUNC 0x1 // the search hit --ed-budget: values above the cost reached are NA
#define RB3_ED_F_REP   0x2 // a tile was skipped as repetitive (> --ed-max-occ occurrences)
#define RB3_ED_F_NOSEED 0x4 // no usable tile (read too short, or every tile skipped / contains N)
#define RB3_ED_F_LOC   0x8 // an aligned interval had more copies than --ed-max-loc: some were not located
#define RB3_ED_F_SKIP  0x10 // not in the --ed-frac subsample: not computed (all NA)
void rb3_edhist_add(rb3_edhist_t *h, int unit, const uint8_t *v, int32_t gbest, int flags);

// TSV: '#' header lines, then "unit measure k count frac cum_frac err_expected post_err" rows.
// measure = best | asm | wrong | nbest; the k = cap+1 row is "NA". The best rows carry the
// two-Poisson fit of errprof.c (component 0 = sequencing error; err_expected = expected count
// from sequencing error alone, post_err = posterior that a best hit at k edits is sequencing
// error); other measures print '.' there.
void rb3_edhist_write(FILE *fp, const rb3_edhist_t *h, const char *const *asm_name, double alpha);

// npy (uint8, C order, shape (rows, cols)): write a fixed-size header now, rewrite it with the
// final row count at the end (the file must be seekable).
#define RB3_ED_NPY_HDR 128
void rb3_ed_npy_header(FILE *fp, int64_t rows, int64_t cols);

#ifdef __cplusplus
}
#endif

#endif
