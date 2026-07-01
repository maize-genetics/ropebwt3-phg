#ifndef RB3_LIFT_H
#define RB3_LIFT_H

#include <stdint.h>

// The carrier->reference coordinate liftover (the "second SSA"). Built by
// `ropebwt3 lift`, consulted by `refmap --lift`. See lift.c for the method.
typedef struct rb3_lift_s rb3_lift_t;

rb3_lift_t *rb3_lift_restore(const char *fn);
void rb3_lift_destroy(rb3_lift_t *lf);
int64_t rb3_lift_n_seq(const rb3_lift_t *lf);

// Project a carrier hit (csid, cpos) to a reference coordinate. Returns 1 and
// sets *out_rsid,*out_rpos on success; 0 (NULL) when not confidently collinear.
int rb3_lift_project(const rb3_lift_t *lf, void *km, int32_t csid, int64_t cpos,
					 int64_t win, int64_t max_mad, int32_t min_support,
					 int64_t *out_rsid, int64_t *out_rpos);

int main_lift(int argc, char *argv[]);

#endif
