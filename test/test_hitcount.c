// Standalone unit tests for hitcount.c. Deterministic, no index building, no
// CLI. Run with `make test` from the repository root.
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "hitcount.h"

static int n_pass = 0, n_fail = 0;

#define CHECK(cond, msg) do { \
		if (cond) { ++n_pass; } \
		else { ++n_fail; fprintf(stderr, "FAIL: %s (%s:%d)\n", (msg), __FILE__, __LINE__); } \
	} while (0)

static void test_basic(void)
{
	rb3_hitcount_t hc;
	int64_t tgt = -1, n = -1;

	rb3_hitcount_init(&hc, 3);
	CHECK(hc.target == 3 && hc.n == 0, "init sets target and zeroes n");
	CHECK(!rb3_hitcount_reached(&hc), "not reached before any hits");

	rb3_hitcount_add(&hc, 0);
	CHECK(hc.n == 0, "add(is_hit=0) does not increment n");
	CHECK(!rb3_hitcount_reached(&hc), "still not reached after a non-hit");

	rb3_hitcount_add(&hc, 1);
	CHECK(hc.n == 1, "add(is_hit=1) increments n");
	CHECK(!rb3_hitcount_reached(&hc), "not reached at n=1 < target=3");

	rb3_hitcount_add(&hc, 1);
	rb3_hitcount_add(&hc, 1);
	CHECK(hc.n == 3, "n reaches 3 after 3 hits (non-hit adds were no-ops)");
	CHECK(rb3_hitcount_reached(&hc), "reached once n == target");

	rb3_hitcount_add(&hc, 1);
	CHECK(rb3_hitcount_reached(&hc), "still reached once n > target");

	CHECK(!rb3_hitcount_short(&hc, &tgt, &n), "short() is false once the target was met");
}

static void test_off(void)
{
	rb3_hitcount_t hc;
	int64_t tgt = -1, n = -1;
	int i;

	rb3_hitcount_init(&hc, 0); // target<=0 = feature off
	for (i = 0; i < 5; ++i) rb3_hitcount_add(&hc, 1);
	CHECK(hc.n == 5, "add() still counts even when the feature is off (harmless bookkeeping)");
	CHECK(!rb3_hitcount_reached(&hc), "reached() is always false when target<=0, regardless of n");
	CHECK(!rb3_hitcount_short(&hc, &tgt, &n), "short() is always false when target<=0 (no target was ever requested)");

	rb3_hitcount_init(&hc, -7); // negative also counts as "off"
	CHECK(!rb3_hitcount_reached(&hc), "reached() false for a negative target too");
	CHECK(!rb3_hitcount_short(&hc, &tgt, &n), "short() false for a negative target too");
}

static void test_short(void)
{
	rb3_hitcount_t hc;
	int64_t tgt = -1, n = -1;

	rb3_hitcount_init(&hc, 10);
	rb3_hitcount_add(&hc, 1);
	rb3_hitcount_add(&hc, 1);
	rb3_hitcount_add(&hc, 1);
	CHECK(rb3_hitcount_short(&hc, &tgt, &n), "short() true when n(3) < target(10)");
	CHECK(tgt == 10, "short() reports the requested target");
	CHECK(n == 3, "short() reports the actual count found");
}

#define N_THREADS 8
#define ADDS_PER_THREAD 20000

static void *adder_thread(void *arg)
{
	rb3_hitcount_t *hc = (rb3_hitcount_t*)arg;
	int i;
	for (i = 0; i < ADDS_PER_THREAD; ++i)
		rb3_hitcount_add(hc, 1);
	return 0;
}

// The reader (worker_pipeline step 0) and writer (write_refmap, step 2) run on
// different threads in production; this proves the atomic add never loses an
// update under real concurrent access from multiple threads at once.
static void test_concurrency(void)
{
	rb3_hitcount_t hc;
	pthread_t th[N_THREADS];
	int i;

	rb3_hitcount_init(&hc, N_THREADS * ADDS_PER_THREAD);
	for (i = 0; i < N_THREADS; ++i)
		pthread_create(&th[i], 0, adder_thread, &hc);
	for (i = 0; i < N_THREADS; ++i)
		pthread_join(th[i], 0);

	CHECK(hc.n == (int64_t)N_THREADS * ADDS_PER_THREAD,
		  "concurrent adds from multiple threads sum exactly, no lost updates");
	CHECK(rb3_hitcount_reached(&hc), "reached() true once the concurrent total hits the target");
}

int main(void)
{
	test_basic();
	test_off();
	test_short();
	test_concurrency();

	fprintf(stderr, "test_hitcount: %d passed, %d failed\n", n_pass, n_fail);
	return n_fail == 0? 0 : 1;
}
