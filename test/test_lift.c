// Standalone unit tests for the rb3_lift_ridx_* reference-position index and
// rb3_lift_ternary_state (lift.c). Deterministic, hand-computed expectations, no
// CLI, no real .lift build. Run with `make test` from the repository root.
//
// rb3_lift_t is opaque outside lift.c (lift.h), so unlike test_ps4g.c's rb3_sid_t
// fixture, this test can't hand-fill its fields -- instead it hand-writes a small
// binary .lift FILE matching the on-disk format (magic, n_seq, n_pt, off[], then
// n_pt 24-byte (cpos,rpos,csid,rsid) records) and calls the real rb3_lift_restore()
// on it, the same way test_ps4g.c writes a BED fixture and calls the real
// rb3_bed_read(). This exercises the actual file-format round trip.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "rb3priv.h"
#include "io.h"
#include "ps4g.h" // for rb3_gtab_build/rb3_gtab_t -- founder->gamete index, same as production
#include "lift.h"

// misc.o's (unused-by-us) rb3_init() references this dawg.c symbol; stub it out
// rather than linking in the rest of the DAWG matcher just to satisfy the linker
// (same trick as test_ps4g.c).
void rb3_bwtl_init(void) {}

static int n_pass = 0, n_fail = 0;

#define CHECK(cond, msg) do { \
		if (cond) { ++n_pass; } \
		else { ++n_fail; fprintf(stderr, "FAIL: %s (%s:%d)\n", (msg), __FILE__, __LINE__); } \
	} while (0)

// Mirrors rb3_liftpt_t's private layout in lift.c (int64 cpos, int64 rpos, int32
// csid, int32 rsid = 24 bytes, no padding). Kept local: the real struct isn't
// exposed via lift.h. sizeof() is asserted against 24 below as a format-drift guard.
typedef struct { int64_t cpos, rpos; int32_t csid, rsid; } fixture_pt_t;

/*******************************
 * fixture: sequences + gtab   *
 *******************************/

// 6 sequences: 2 reference chromosomes (B73_chr1/2, never carry lift points -- the
// lift builder skips reference hits) and 2 founders (A, B) each with a sequence on
// both chromosomes, so cross-chromosome and cross-founder isolation can be tested.
static rb3_sid_t *make_sid(void)
{
	rb3_sid_t *sid = RB3_CALLOC(rb3_sid_t, 1);
	const char *names[] = { "B73_chr1", "B73_chr2", "A_chr1", "A_chr2", "B_chr1", "B_chr2" };
	int64_t i;
	sid->n_seq = 6;
	sid->name = RB3_MALLOC(char*, sid->n_seq);
	sid->len = RB3_CALLOC(int32_t, sid->n_seq);
	for (i = 0; i < sid->n_seq; ++i) sid->name[i] = rb3_strdup(names[i]);
	return sid;
}

static void free_sid(rb3_sid_t *sid)
{
	int64_t i;
	for (i = 0; i < sid->n_seq; ++i) free(sid->name[i]);
	free(sid->name); free(sid->len); free(sid);
}

// Hand-write a small binary .lift file for the fixture above and return its path.
// Points are listed pre-sorted by (csid,cpos), matching what a real `ropebwt3 lift`
// build produces (rb3_lift_ridx_build itself doesn't rely on this order -- it scans
// every point directly -- but a realistic fixture should still look like a real file).
//
//   csid=2 (A_chr1): rpos 1000,2000,3000 on rsid=0 (chr1)         -- 3 points
//   csid=3 (A_chr2): none                                         -- A has ZERO anchors on chr2
//   csid=4 (B_chr1): rpos 1500,2500 on rsid=0 (chr1)               -- 2 points, offset from A's
//   csid=5 (B_chr2): rpos 1000 on rsid=1 (chr2)                    -- 1 point: a "trap" at the
//     exact position A is queried at on chr2, to prove a no-anchor answer for A there isn't
//     accidentally satisfied by B's unrelated anchor.
static const char *write_lift_fixture(const char *outdir)
{
	static char fn[512];
	fixture_pt_t pt[6] = {
		{100, 1000, 2, 0}, {200, 2000, 2, 0}, {300, 3000, 2, 0}, // csid=2 (A_chr1)
		{150, 1500, 4, 0}, {250, 2500, 4, 0},                     // csid=4 (B_chr1)
		{100, 1000, 5, 1},                                        // csid=5 (B_chr2)
	};
	int64_t off[7] = { 0, 0, 0, 3, 3, 5, 6 }; // off[csid]..off[csid+1], csid in [0,6)
	int64_t n_seq = 6, n_pt = 6;
	FILE *fp;
	CHECK(sizeof(fixture_pt_t) == 24, "fixture: local liftpt layout matches the documented 24-byte on-disk record (format-drift guard)");
	snprintf(fn, sizeof(fn), "%s/fixture.lift", outdir);
	fp = fopen(fn, "wb");
	fwrite("LIFT\1", 1, 5, fp);
	fwrite(&n_seq, 8, 1, fp);
	fwrite(&n_pt, 8, 1, fp);
	fwrite(off, 8, 7, fp);
	fwrite(pt, sizeof(fixture_pt_t), 6, fp);
	fclose(fp);
	return fn;
}

/*******************************
 * rb3_lift_ridx_build/destroy *
 *******************************/

static void test_ridx_build_null_safety(const rb3_lift_t *lf, const int32_t *sid2g, int32_t n_gamete)
{
	CHECK(rb3_lift_ridx_build(0, sid2g, n_gamete) == 0, "ridx_build: NULL lf -> NULL");
	CHECK(rb3_lift_ridx_build(lf, 0, n_gamete) == 0, "ridx_build: NULL sid2g -> NULL");
	CHECK(rb3_lift_ridx_build(lf, sid2g, 0) == 0, "ridx_build: n_gamete=0 -> NULL");
	CHECK(rb3_lift_ridx_build(lf, sid2g, -1) == 0, "ridx_build: n_gamete<0 -> NULL");
	rb3_lift_ridx_destroy(0); // must not crash
}

/*******************************
 * rb3_lift_nearest_ref        *
 *******************************/

// gamete indices, from rb3_gtab_build's alphabetical ordering over {A,B,B73}
#define G_A   0
#define G_B   1
#define G_B73 2
#define RSID_CHR1 0
#define RSID_CHR2 1

static void test_nearest_ref(const rb3_lift_ridx_t *ridx)
{
	CHECK(rb3_lift_nearest_ref(ridx, G_A, RSID_CHR1, 1000) == 0, "nearest_ref: A on an exact anchor (1000) -> dist 0");
	CHECK(rb3_lift_nearest_ref(ridx, G_A, RSID_CHR1, 1050) == 50, "nearest_ref: A, left anchor (1000) nearer than right (2000)");
	CHECK(rb3_lift_nearest_ref(ridx, G_A, RSID_CHR1, 1950) == 50, "nearest_ref: A, right anchor (2000) nearer than left (1000)");
	CHECK(rb3_lift_nearest_ref(ridx, G_A, RSID_CHR1, 500) == 500, "nearest_ref: A, before all anchors -- only a right neighbor (1000) exists");
	CHECK(rb3_lift_nearest_ref(ridx, G_A, RSID_CHR1, 3500) == 500, "nearest_ref: A, after all anchors -- only a left neighbor (3000) exists");
	CHECK(rb3_lift_nearest_ref(ridx, G_A, RSID_CHR2, 1000) == -1,
		  "nearest_ref: A has ZERO anchors on chr2 -- not satisfied by B's unrelated chr2 anchor at the same position (no cross-founder bleeding)");
	CHECK(rb3_lift_nearest_ref(ridx, G_B, RSID_CHR1, 1500) == 0, "nearest_ref: B on an exact anchor (1500) -> dist 0");
	CHECK(rb3_lift_nearest_ref(ridx, G_B, RSID_CHR1, 1000) == 500,
		  "nearest_ref: B's nearest own anchor from 1000 is 1500 (dist 500), NOT A's anchor which sits exactly at 1000 (would wrongly give dist 0) -- no cross-founder bleeding");
	CHECK(rb3_lift_nearest_ref(ridx, G_B73, RSID_CHR1, 12345) == -1, "nearest_ref: the reference gamete has zero anchors anywhere, by construction");

	CHECK(rb3_lift_nearest_ref(ridx, -1, RSID_CHR1, 1000) == -1, "nearest_ref: negative gamete -> -1");
	CHECK(rb3_lift_nearest_ref(ridx, 3, RSID_CHR1, 1000) == -1, "nearest_ref: gamete >= n_gamete(3) -> -1");
	CHECK(rb3_lift_nearest_ref(ridx, G_A, -1, 1000) == -1, "nearest_ref: negative rsid -> -1");
	CHECK(rb3_lift_nearest_ref(ridx, G_A, 6, 1000) == -1, "nearest_ref: rsid >= n_seq(6) -> -1");
	CHECK(rb3_lift_nearest_ref(0, G_A, RSID_CHR1, 1000) == -1, "nearest_ref: NULL ridx -> -1");
}

/*******************************
 * rb3_lift_ternary_state      *
 *******************************/

static void test_ternary_state(void)
{
	CHECK(rb3_lift_ternary_state(1, 12345, 0) == 1, "ternary_state: in_gameteset -> match(1), regardless of distance");
	CHECK(rb3_lift_ternary_state(1, -1, 0) == 1, "ternary_state: in_gameteset -> match(1), even with dist=-1 (no anchor)");
	CHECK(rb3_lift_ternary_state(0, 0, 0) == 0, "ternary_state: not in gameteset, dist==thresh(0) -> diverged(0) (boundary)");
	CHECK(rb3_lift_ternary_state(0, 2000, 2000) == 0, "ternary_state: not in gameteset, dist==thresh(2000) -> diverged(0) (boundary)");
	CHECK(rb3_lift_ternary_state(0, 2001, 2000) == -1, "ternary_state: not in gameteset, dist just over thresh -> deletion(-1)");
	CHECK(rb3_lift_ternary_state(0, -1, 2000) == -1, "ternary_state: not in gameteset, no anchor at all (dist=-1) -> deletion(-1)");
}

int main(int argc, char *argv[])
{
	const char *outdir = argc > 1? argv[1] : "test/output";
	const char *lift_fn;
	rb3_sid_t *sid;
	rb3_gtab_t *gtab;
	rb3_lift_t *lf;
	rb3_lift_ridx_t *ridx;

	sid = make_sid();
	gtab = rb3_gtab_build(sid);
	CHECK(gtab->n_gamete == 3, "gtab: 3 distinct samples (A,B,B73)");
	CHECK(gtab->sid2g[0] == G_B73 && gtab->sid2g[1] == G_B73, "gtab: B73_chr1/2 -> gamete B73");
	CHECK(gtab->sid2g[2] == G_A && gtab->sid2g[3] == G_A, "gtab: A_chr1/2 -> gamete A");
	CHECK(gtab->sid2g[4] == G_B && gtab->sid2g[5] == G_B, "gtab: B_chr1/2 -> gamete B");

	lift_fn = write_lift_fixture(outdir);
	lf = rb3_lift_restore(lift_fn);
	CHECK(lf != 0, "lift fixture: rb3_lift_restore parsed the hand-written file");
	CHECK(rb3_lift_n_seq(lf) == 6, "lift fixture: n_seq round-trips through the file");

	test_ridx_build_null_safety(lf, gtab->sid2g, gtab->n_gamete);

	ridx = rb3_lift_ridx_build(lf, gtab->sid2g, gtab->n_gamete);
	CHECK(ridx != 0, "ridx_build: built successfully from the fixture");
	test_nearest_ref(ridx);
	rb3_lift_ridx_destroy(ridx); // must not crash

	test_ternary_state();

	rb3_lift_destroy(lf);
	rb3_gtab_destroy(gtab);
	free_sid(sid);

	fprintf(stderr, "test_lift: %d passed, %d failed\n", n_pass, n_fail);
	return n_fail == 0? 0 : 1;
}
