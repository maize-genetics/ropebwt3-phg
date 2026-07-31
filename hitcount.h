#ifndef RB3_HITCOUNT_H
#define RB3_HITCOUNT_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Generic, refmap-status-agnostic counter for "stop once N interesting
// records have been produced": the caller decides what counts as a "hit".
// target<=0 means the feature is off -- reached()/short() always return
// false regardless of n, so callers can unconditionally wire this in
// without a separate on/off branch.
typedef struct {
	int64_t target;
	int64_t n;
} rb3_hitcount_t;

void rb3_hitcount_init(rb3_hitcount_t *hc, int64_t target);

// Call once per output record; is_hit = does this record count toward the
// target. Safe to call concurrently with itself and with rb3_hitcount_reached
// from other threads (updates via an atomic add).
void rb3_hitcount_add(rb3_hitcount_t *hc, int is_hit);

// True once hc->n >= hc->target (target>0 required).
int rb3_hitcount_reached(const rb3_hitcount_t *hc);

// True if target>0 && n<target; fills *target_out/*n_out with the values at
// the time of the call. Does no I/O itself -- the caller formats/prints its
// own warning using the returned values.
int rb3_hitcount_short(const rb3_hitcount_t *hc, int64_t *target_out, int64_t *n_out);

#ifdef __cplusplus
}
#endif

#endif
