// Standalone unit tests for pairing.c (paired-end helpers). Deterministic, no index
// building, no CLI. Run with `make test` from the repository root.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pairing.h"

static int n_pass = 0, n_fail = 0;

#define CHECK(cond, msg) do { \
		if (cond) { ++n_pass; } \
		else { ++n_fail; fprintf(stderr, "FAIL: %s (%s:%d)\n", (msg), __FILE__, __LINE__); } \
	} while (0)

static void test_names(void)
{
	CHECK(rb3_pe_name_len("r1/1") == 2, "name_len strips /1");
	CHECK(rb3_pe_name_len("r1/2") == 2, "name_len strips /2");
	CHECK(rb3_pe_name_len("r1/3") == 4, "name_len keeps /3");
	CHECK(rb3_pe_name_len("r1") == 2, "name_len: no suffix");
	CHECK(rb3_pe_name_len("/1") == 0, "name_len: suffix only");
	CHECK(rb3_pe_name_len("1") == 1, "name_len: single char");
	CHECK(rb3_pe_names_match("chr1_10_400_0:0:0_0:0:0_a/1", "chr1_10_400_0:0:0_0:0:0_a/2"), "wgsim mates match");
	CHECK(rb3_pe_names_match("frag7", "frag7"), "identical names (no suffix) match");
	CHECK(rb3_pe_names_match("frag7/1", "frag7"), "suffix on one side only still matches");
	CHECK(!rb3_pe_names_match("frag7/1", "frag8/2"), "different fragments do not match");
	CHECK(!rb3_pe_names_match("frag7/1", "frag77/2"), "prefix of another name does not match");
	CHECK(!rb3_pe_names_match(0, "x"), "NULL name never matches");
	CHECK(!rb3_pe_names_match("x", 0), "NULL name never matches (b)");
}

static void test_isize(void)
{
	rb3_isize_t is;
	int64_t span[100];
	int32_t i;

	// too few spans
	for (i = 0; i < 10; ++i) span[i] = 500;
	CHECK(rb3_pe_isize_estimate(span, 10, 25, &is) == 0 && is.ok == 0, "isize: fewer than min_n -> not ok");
	CHECK(rb3_pe_isize_estimate(span, 0, 0, &is) == 0, "isize: n=0 -> not ok");

	// constant spans: sd 0, IQR 0 -> max_ins == the span itself
	for (i = 0; i < 40; ++i) span[i] = 500;
	CHECK(rb3_pe_isize_estimate(span, 40, 25, &is) == 1 && is.ok, "isize: constant -> ok");
	CHECK(is.n == 40 && fabs(is.mean - 500.0) < 1e-9 && fabs(is.sd) < 1e-9 && is.max_ins == 500, "isize: constant mean/sd/max_ins");

	// 40 spans 480..519 (unsorted) + 2 gross outliers that the IQR filter must drop.
	// sorted: 3, 480..519, 50000 (sorted[k] = 479+k for k in 1..40). p25 = sorted[10] = 489,
	// p75 = sorted[31] = 510 -> IQR 21, keep [447, 552]; kept = 480..519: mean 499.5,
	// sd = sqrt((40^2-1)/12) = 11.5434. max_ins = max(round(499.5+4*11.5434) = 546, 510+63 = 573).
	for (i = 0; i < 40; ++i) span[i] = 480 + (i * 7) % 40;
	span[40] = 50000, span[41] = 3;
	CHECK(rb3_pe_isize_estimate(span, 42, 25, &is) == 1, "isize: noisy -> ok");
	CHECK(is.n == 40, "isize: both outliers dropped");
	CHECK(fabs(is.mean - 499.5) < 1e-9, "isize: mean of kept spans");
	CHECK(fabs(is.sd - sqrt((40.0 * 40.0 - 1.0) / 12.0)) < 1e-9, "isize: population sd of kept spans");
	CHECK(is.max_ins == 573, "isize: max_ins = max(mean+4sd, p75+3IQR)");
	CHECK(span[0] == 3 && span[41] == 50000, "isize: spans sorted in place");

	// wide distribution: mean + 4sd is the larger bound
	for (i = 0; i < 100; ++i) span[i] = (i < 50)? 300 : 700; // p25=300 p75=700 IQR=400; mean 500 sd 200
	CHECK(rb3_pe_isize_estimate(span, 100, 25, &is) == 1 && is.n == 100, "isize: bimodal kept");
	CHECK(is.max_ins == 1900, "isize: bimodal max_ins = max(500+800, 700+1200) = 1900");
}

static rb3_pe_mate_t mate(int rank, int64_t rsid, int64_t pos)
{
	rb3_pe_mate_t m;
	memset(&m, 0, sizeof(m));
	m.rank = rank, m.rsid = rsid, m.pos = pos;
	return m;
}

static void test_resolve(void)
{
	rb3_pe_mate_t m[2];
	rb3_pe_res_t r;
	int64_t prs[4], ppos[4], prs2[2], ppos2[2];

	m[0] = mate(0, -1, -1), m[1] = mate(0, -1, -1);
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_NONE && r.anchor == -1 && r.rescued == -1, "resolve: both unplaced -> NONE");

	m[0] = mate(0, -1, -1), m[1] = mate(1, 3, 5000);
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_SOLO && r.anchor == 1, "resolve: only mate 1 placed -> SOLO anchored on 1");
	m[0] = mate(2, 3, 5000), m[1] = mate(0, -1, -1);
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_SOLO && r.anchor == 0, "resolve: only mate 0 placed -> SOLO anchored on 0");

	m[0] = mate(2, 3, 5000), m[1] = mate(2, 3, 5900);
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_CONC && r.anchor == 0, "resolve: exact+exact within max_dist -> CONC anchored on mate 0");
	m[0] = mate(2, 3, 5000), m[1] = mate(2, 3, 6000);
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_CONC, "resolve: distance == max_dist is concordant");
	m[0] = mate(2, 3, 5000), m[1] = mate(2, 3, 6001);
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_DISC, "resolve: exact+exact one past max_dist -> DISC (no projections)");
	m[0] = mate(1, 3, 5000), m[1] = mate(2, 3, 4200);
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_CONC && r.anchor == 1, "resolve: CONC anchors on the exact mate");
	m[0] = mate(2, 3, 5000), m[1] = mate(2, 4, 5000);
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_DISC, "resolve: different contigs -> DISC");

	// projected mate 1 lands 80kb away, but 2 of its 4 assembly projections are near the exact mate
	m[0] = mate(2, 3, 100000), m[1] = mate(1, 3, 180000);
	prs[0] = 3, ppos[0] = 180000; // the (wrong) majority
	prs[1] = 3, ppos[1] = 100700;
	prs[2] = 3, ppos[2] = 100300;
	prs[3] = 5, ppos[3] = 100500; // right position, wrong contig: not counted
	m[1].n_proj = 4, m[1].proj_rsid = prs, m[1].proj_pos = ppos;
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_RESC && r.anchor == 0 && r.rescued == 1, "resolve: projected mate rescued next to exact mate");
	CHECK(r.resc_rsid == 3 && r.resc_support == 2, "resolve: rescue counts only same-contig projections in range");
	CHECK(r.resc_pos == 100700, "resolve: rescued position = upper median of the near projections");

	// projected mate 0 has nothing near exact mate 1 -> DISC; with trust_exact it moves onto mate 1
	m[0] = mate(1, 3, 100000), m[1] = mate(2, 3, 300000);
	m[0].n_proj = 4, m[0].proj_rsid = prs, m[0].proj_pos = ppos; // none near 300000
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_DISC, "resolve: projected mate without projections near the exact mate -> DISC");
	rb3_pe_resolve(m, 1000, 1, &r);
	CHECK(r.cls == RB3_PE_RESC && r.anchor == 1 && r.rescued == 0 && r.resc_pos == 300000 && r.resc_support == 0,
		  "resolve: trust_exact moves the projected mate onto the exact mate");
	{ // trust_exact never applies to two exact or two projected mates
		rb3_pe_mate_t e[2];
		e[0] = mate(2, 3, 5000), e[1] = mate(2, 3, 90000);
		rb3_pe_resolve(e, 1000, 1, &r);
		CHECK(r.cls == RB3_PE_DISC, "resolve: trust_exact leaves EXACT/EXACT discordant pairs alone");
		e[0] = mate(1, 3, 5000), e[1] = mate(1, 3, 90000);
		rb3_pe_resolve(e, 1000, 1, &r);
		CHECK(r.cls == RB3_PE_DISC, "resolve: trust_exact leaves PLACED/PLACED discordant pairs alone");
	}
	{ // tie-break among an exact read's reference copies: the copy next to the projected mate wins
		rb3_pe_mate_t e[2];
		int64_t crs[3] = { 3, 3, 4 }, cpos[3] = { 300000, 100400, 100000 };
		e[0] = mate(1, 3, 100000), e[1] = mate(2, 3, 300000); // exact mate reported at its first copy
		e[1].n_proj = 3, e[1].proj_rsid = crs, e[1].proj_pos = cpos;
		rb3_pe_resolve(e, 1000, 0, &r);
		CHECK(r.cls == RB3_PE_RESC && r.anchor == 0 && r.rescued == 1 && r.resc_pos == 100400 && r.resc_support == 1,
			  "resolve: exact mate moved to its reference copy next to the projected mate");
		// if the projected mate ALSO has a projection near the exact copy, moving the projected mate wins
		e[0].n_proj = 4, e[0].proj_rsid = prs, e[0].proj_pos = ppos; ppos[0] = 300200;
		rb3_pe_resolve(e, 1000, 0, &r);
		CHECK(r.cls == RB3_PE_RESC && r.anchor == 1 && r.rescued == 0 && r.resc_pos == 300200,
			  "resolve: moving the less confident mate is preferred over moving the exact one");
		ppos[0] = 180000;
		// two exact mates: one of them has a copy near the other -> moved (tie on rank, support decides)
		e[0] = mate(2, 3, 5000), e[1] = mate(2, 3, 300000);
		e[1].n_proj = 3, e[1].proj_rsid = crs, e[1].proj_pos = cpos; cpos[1] = 5300;
		rb3_pe_resolve(e, 1000, 0, &r);
		CHECK(r.cls == RB3_PE_RESC && r.anchor == 0 && r.rescued == 1 && r.resc_pos == 5300, "resolve: EXACT/EXACT tie-break between copies");
	}
	m[1].pos = 180500; // now projection 180000 of mate 0 is near the exact mate 1
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_RESC && r.anchor == 1 && r.rescued == 0 && r.resc_pos == 180000, "resolve: rescue toward the exact mate 1");

	// both projected and discordant: the mate with more agreeing projections is moved
	m[0] = mate(1, 3, 100000), m[1] = mate(1, 3, 180000);
	m[0].n_proj = 1, prs2[0] = 3, ppos2[0] = 179500, m[0].proj_rsid = prs2, m[0].proj_pos = ppos2; // 1 near mate 1
	m[1].n_proj = 4, m[1].proj_rsid = prs, m[1].proj_pos = ppos; // 2 near mate 0
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_RESC && r.anchor == 0 && r.rescued == 1 && r.resc_support == 2, "resolve: projected pair -> move the better-supported mate");
	m[1].n_proj = 2; // now mate 1 has only projection 180000 + 100700 -> 1 near mate 0: tie -> anchor 0
	rb3_pe_resolve(m, 1000, 0, &r);
	CHECK(r.cls == RB3_PE_RESC && r.anchor == 0 && r.resc_support == 1 && r.resc_pos == 100700, "resolve: support tie -> anchor on mate 0");
}

static int same(const int32_t *a, int32_t na, const int32_t *b, int32_t nb)
{
	return na == nb && memcmp(a, b, na * sizeof(int32_t)) == 0;
}

static void test_near_gametes(void)
{
	// mate b: 5 projections (founders 4,2,4,7,1). Near (contig 3, 100000 +- 1000): idx 1 (2), 2 (4), 3 (7, edge).
	int64_t prs[5] = { 3, 3, 3, 3, 9 }, ppos[5] = { 180000, 100700, 100300, 101000, 100000 };
	int32_t gid[5] = { 4, 2, 4, 7, 1 }, out[5], n;
	int32_t e_near[] = { 2, 4, 7 }, e_one[] = { 4 };
	rb3_pe_mate_t b = mate(1, 3, 180000);
	b.n_proj = 5, b.proj_rsid = prs, b.proj_pos = ppos;
	CHECK(rb3_pe_near_gametes(&b, 3, 100000, 1000, out) == 0, "near_gametes: no proj_gid -> 0");
	b.proj_gid = gid;
	n = rb3_pe_near_gametes(&b, 3, 100000, 1000, out);
	CHECK(same(out, n, e_near, 3), "near_gametes: same contig, within max_dist (inclusive), sorted+deduped");
	n = rb3_pe_near_gametes(&b, 3, 100000, 500, out);
	CHECK(same(out, n, e_one, 1), "near_gametes: narrower window keeps only founder 4 (100300)");
	CHECK(rb3_pe_near_gametes(&b, 3, 50000, 1000, out) == 0, "near_gametes: nothing near -> 0");
	CHECK(rb3_pe_near_gametes(&b, 5, 100000, 1000, out) == 0, "near_gametes: other contig -> 0");
	b.n_proj = 0;
	CHECK(rb3_pe_near_gametes(&b, 3, 100000, 1000, out) == 0, "near_gametes: no projections -> 0");
}

static void test_event_gametes(void)
{
	int32_t out[32], n;
	int32_t a[] = { 7, 2, 0, 2, 5 };   // unsorted, with a duplicate: {0,2,5,7}
	int32_t b[] = { 5, 9, 2 };         // {2,5,9}
	int32_t c[] = { 1, 3 };            // disjoint from a
	int32_t e_ab[] = { 2, 5 }, e_a[] = { 0, 2, 5, 7 }, e_b[] = { 2, 5, 9 };

	n = rb3_pe_event_gametes(a, 5, 0, b, 3, 0, out);
	CHECK(same(out, n, e_ab, 2), "gametes: intersection, sorted+deduped");
	n = rb3_pe_event_gametes(b, 3, 0, a, 5, 0, out);
	CHECK(same(out, n, e_ab, 2), "gametes: intersection is symmetric");
	n = rb3_pe_event_gametes(a, 5, 0, 0, 0, 0, out);
	CHECK(same(out, n, e_a, 4), "gametes: mate has no set -> anchor's set");
	n = rb3_pe_event_gametes(a, 5, 0, b, 3, 1, out);
	CHECK(same(out, n, e_a, 4), "gametes: mate's set capped -> anchor's set");
	n = rb3_pe_event_gametes(a, 5, 1, b, 3, 0, out);
	CHECK(same(out, n, e_b, 3), "gametes: anchor's set capped -> mate's set");
	n = rb3_pe_event_gametes(a, 5, 1, b, 3, 1, out);
	CHECK(same(out, n, e_a, 4), "gametes: both capped -> anchor's set");
	n = rb3_pe_event_gametes(a, 5, 0, c, 2, 0, out);
	CHECK(same(out, n, e_a, 4), "gametes: empty intersection -> anchor's set");
	n = rb3_pe_event_gametes(a, 5, 0, a, 5, 0, out);
	CHECK(same(out, n, e_a, 4), "gametes: self-intersection = deduped set");
	n = rb3_pe_event_gametes(0, 0, 0, b, 3, 0, out);
	CHECK(same(out, n, e_b, 3), "gametes: empty anchor set -> mate's set");
	CHECK(a[0] == 7 && b[0] == 5, "gametes: inputs are not modified");
}

int main(void)
{
	test_names();
	test_isize();
	test_resolve();
	test_event_gametes();
	test_near_gametes();
	printf("test_pairing: %d passed, %d failed\n", n_pass, n_fail);
	return n_fail > 0;
}
