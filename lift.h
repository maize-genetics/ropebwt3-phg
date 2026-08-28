#ifndef RB3_LIFT_H
#define RB3_LIFT_H

#include <stdint.h>

// The assembly->reference coordinate liftover (the "second SSA"). Built by
// `ropebwt3 lift`, consulted by `refmap --lift`. See lift.c for the method.
typedef struct rb3_lift_s rb3_lift_t;

rb3_lift_t *rb3_lift_restore(const char *fn);
void rb3_lift_destroy(rb3_lift_t *lf);
int64_t rb3_lift_n_seq(const rb3_lift_t *lf);

// Project a assembly hit (csid, cpos) to a reference coordinate. Returns 1 and
// sets *out_rsid,*out_rpos on success; 0 (NULL) when not confidently collinear.
int rb3_lift_project(const rb3_lift_t *lf, void *km, int32_t csid, int64_t cpos,
					 int64_t win, int64_t max_mad, int32_t min_support,
					 int64_t *out_rsid, int64_t *out_rpos);

// Like rb3_lift_project, but when the locus is not colinear (a PAV/insertion between
// the flanking anchors) it returns the nearest reference anchor = the closest
// breakpoint, with *out_mode=1 (colinear projection sets *out_mode=0).
int rb3_lift_project_bp(const rb3_lift_t *lf, void *km, int32_t csid, int64_t cpos,
						int64_t win, int64_t max_mad, int32_t min_support,
						int64_t *out_rsid, int64_t *out_rpos, int *out_mode);

// A reference-position index over the same liftover points as rb3_lift_t, but
// grouped by (gamete, reference sequence) and sorted by reference position --
// rb3_lift_t's own pt[] is sorted by (assembly sequence, assembly position), which
// cannot answer "nearest anchor to reference position X for founder G" directly.
// Built once per run (see rb3_lift_ridx_build); the index only holds a permutation
// into the rb3_lift_t it was built from, so `lf` must outlive the returned index.
typedef struct rb3_lift_ridx_s rb3_lift_ridx_t;

// sid2g[seqIdx] = gamete (founder) index of sequence seqIdx, seqIdx in [0,n_seq)
// (rb3_gtab_t.sid2g from ps4g.h) -- needed because a liftover point's csid is an
// assembly SEQUENCE index (one founder owns many chromosomes/scaffolds), not a
// gamete index directly.
rb3_lift_ridx_t *rb3_lift_ridx_build(const rb3_lift_t *lf, const int32_t *sid2g, int32_t n_gamete);
void rb3_lift_ridx_destroy(rb3_lift_ridx_t *ridx);

// Distance (bp) from reference position (rsid,rpos) to the nearest liftover anchor
// belonging to `gamete` on that same reference sequence. Returns -1 if that gamete
// has no anchor at all on rsid (e.g. the reference gamete itself, which has zero
// liftover points by construction -- it never appears as an assembly/csid).
int64_t rb3_lift_nearest_ref(const rb3_lift_ridx_t *ridx, int32_t gamete, int64_t rsid, int64_t rpos);

// Ternary read-sharing state for one (founder, row): match (1) if the founder is in
// the row's matched-founder set, independent of distance; otherwise diverged (0) if
// dist is within thresh bp of an anchor, deletion (-1) if dist exceeds thresh or the
// founder has no anchor at all (dist < 0) on this reference sequence. A small, pure
// function so it can be unit-tested on its own, independent of the liftover index.
int8_t rb3_lift_ternary_state(int in_gameteset, int64_t dist, int64_t thresh);

int main_lift(int argc, char *argv[]);

#endif
