#include "hitcount.h"

void rb3_hitcount_init(rb3_hitcount_t *hc, int64_t target)
{
	hc->target = target;
	hc->n = 0;
}

void rb3_hitcount_add(rb3_hitcount_t *hc, int is_hit)
{
	if (is_hit) __atomic_fetch_add(&hc->n, 1, __ATOMIC_RELAXED);
}

int rb3_hitcount_reached(const rb3_hitcount_t *hc)
{
	if (hc->target <= 0) return 0; // feature off
	return __atomic_load_n(&hc->n, __ATOMIC_RELAXED) >= hc->target;
}

int rb3_hitcount_short(const rb3_hitcount_t *hc, int64_t *target_out, int64_t *n_out)
{
	int64_t n;
	if (hc->target <= 0) return 0; // feature off
	n = __atomic_load_n(&hc->n, __ATOMIC_RELAXED);
	if (n >= hc->target) return 0;
	*target_out = hc->target, *n_out = n;
	return 1;
}
