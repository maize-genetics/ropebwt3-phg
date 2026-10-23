// Standalone unit tests for ps4g.c (no index building, no CLI). Deterministic,
// hand-computed expectations; run with `make test` from the repository root.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "rb3priv.h"
#include "ps4g.h"

// misc.o's (unused-by-us) rb3_init() references this dawg.c symbol; stub it out
// rather than linking in the rest of the DAWG matcher just to satisfy the linker.
void rb3_bwtl_init(void) {}

static int n_pass = 0, n_fail = 0;

#define CHECK(cond, msg) do { \
		if (cond) { ++n_pass; } \
		else { ++n_fail; fprintf(stderr, "FAIL: %s (%s:%d)\n", (msg), __FILE__, __LINE__); } \
	} while (0)

// evaluate a/b into locals first: a is often slurp_line(fp), which has the side
// effect of advancing the file; naively repeating (a) in the macro body would
// read multiple lines per check.
#define CHECK_STREQ(a, b, msg) do { \
		const char *_a = (a), *_b = (b); \
		if (_a && _b && strcmp(_a, _b) == 0) { ++n_pass; } \
		else { ++n_fail; fprintf(stderr, "FAIL: %s: got '%s', want '%s' (%s:%d)\n", (msg), _a?_a:"(null)", _b?_b:"(null)", __FILE__, __LINE__); } \
	} while (0)

static char *slurp_line(FILE *fp)
{
	static char buf[4096];
	if (fgets(buf, sizeof(buf), fp) == 0) return 0;
	buf[strcspn(buf, "\n")] = 0;
	return buf;
}

/*************************
 * rb3_gtab_build/find   *
 *************************/

static rb3_sid_t *make_sid(void)
{
	rb3_sid_t *sid = RB3_CALLOC(rb3_sid_t, 1);
	const char *names[] = { "B73_chr1", "B73_scaf_1", "A_chr1", "B_chr1", "C_chr1" };
	int64_t i;
	sid->n_seq = 5;
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

static void test_gtab(rb3_sid_t *sid, rb3_gtab_t *g)
{
	// sample names sorted ascending: "A" < "B" < "B73" < "C" (strcmp, "B" is a prefix of "B73")
	CHECK(g->n_gamete == 4, "gtab: 4 distinct samples (A,B,B73,C)");
	CHECK_STREQ(g->name[0], "A", "gtab: gamete 0 name");
	CHECK_STREQ(g->name[1], "B", "gtab: gamete 1 name");
	CHECK_STREQ(g->name[2], "B73", "gtab: gamete 2 name");
	CHECK_STREQ(g->name[3], "C", "gtab: gamete 3 name");

	CHECK(g->sid2g[0] == 2, "gtab: B73_chr1 -> gamete B73(2)");
	CHECK(g->sid2g[1] == 2, "gtab: B73_scaf_1 -> gamete B73(2) (scaffold shares the sample prefix)");
	CHECK(g->sid2g[2] == 0, "gtab: A_chr1 -> gamete A(0)");
	CHECK(g->sid2g[3] == 1, "gtab: B_chr1 -> gamete B(1)");
	CHECK(g->sid2g[4] == 3, "gtab: C_chr1 -> gamete C(3)");

	CHECK(rb3_gtab_find(g, "B73") == 2, "gtab_find: B73 found at index 2");
	CHECK(rb3_gtab_find(g, "A") == 0, "gtab_find: A found at index 0");
	CHECK(rb3_gtab_find(g, "nonexistent") == -1, "gtab_find: unknown sample returns -1");
	(void)sid;
}

/*************************
 * rb3_ps4g_acc_add      *
 *************************/

static void test_acc_dedup(void)
{
	rb3_ps4g_acc_t *acc = rb3_ps4g_acc_init(100);
	int32_t gametes[] = { 2, 0, 0, 3, 2 }; // unsorted with duplicates
	rb3_ps4g_acc_add(acc, 0, 250, gametes, 5); // bin = 250/100 = 2
	CHECK(acc->n_ev == 1, "acc_add: one event recorded");
	CHECK(acc->ev_ref_sid[0] == 0, "acc_add: ref_sid preserved");
	CHECK(acc->ev_bin[0] == 2, "acc_add: pos/bin_size computed correctly");
	CHECK(acc->ev_glen[0] == 3, "acc_add: duplicate gametes collapsed to 3 unique");
	CHECK(acc->garena[acc->ev_goff[0]+0] == 0 &&
		  acc->garena[acc->ev_goff[0]+1] == 2 &&
		  acc->garena[acc->ev_goff[0]+2] == 3, "acc_add: gamete list sorted ascending {0,2,3}");
	rb3_ps4g_acc_destroy(acc);
}

/*******************************************
 * full finalize round-trip: PS4G + npy    *
 *******************************************/

static void write_bed_fixture(const char *fn)
{
	FILE *fp = fopen(fn, "w");
	fprintf(fp, "chr1\t200\t300\tA\tB73\n");   // labels the (ref0,bin2) row: gA=A(0), gB=B73(2)
	fprintf(fp, "scaf_1\t0\t10\tC\n");         // labels the (ref1,bin0) row homozygous C(3)
	fclose(fp);
}

// hand-rolled .npy v1.0 reader: enough to validate shape/dtype/data for this test
static int32_t *read_npy_i4(const char *fn, int64_t *rows, int64_t *cols)
{
	FILE *fp = fopen(fn, "rb");
	uint8_t magic[6], ver[2], hlen_b[2];
	uint16_t hlen;
	char *hdr;
	int32_t *data;
	int64_t n;
	if (fp == 0) return 0;
	if (fread(magic, 1, 6, fp) != 6 || memcmp(magic, "\x93NUMPY", 6) != 0) { fclose(fp); return 0; }
	if (fread(ver, 1, 2, fp) != 2) { fclose(fp); return 0; }
	if (fread(hlen_b, 1, 2, fp) != 2) { fclose(fp); return 0; }
	hlen = (uint16_t)(hlen_b[0] | (hlen_b[1] << 8));
	hdr = RB3_MALLOC(char, hlen + 1);
	if (fread(hdr, 1, hlen, fp) != (size_t)hlen) { free(hdr); fclose(fp); return 0; }
	hdr[hlen] = 0;
	if (strstr(hdr, "'descr': '<i4'") == 0) { fprintf(stderr, "npy dtype is not <i4: %s\n", hdr); free(hdr); fclose(fp); return 0; }
	{
		char *p = strstr(hdr, "'shape': (");
		if (p == 0 || sscanf(p, "'shape': (%ld, %ld)", rows, cols) != 2) { free(hdr); fclose(fp); return 0; }
	}
	free(hdr);
	n = (*rows) * (*cols);
	data = RB3_MALLOC(int32_t, n);
	if (fread(data, sizeof(int32_t), n, fp) != (size_t)n) { free(data); fclose(fp); return 0; }
	fclose(fp);
	return data;
}

static void test_finalize(const char *outdir)
{
	rb3_sid_t *sid = make_sid();
	rb3_gtab_t *g = rb3_gtab_build(sid);
	rb3_ps4g_acc_t *acc = rb3_ps4g_acc_init(100);
	char ps4g_fn[512], npy_fn[512], bed_fn[512];
	rb3_bed_t *bed;
	FILE *fp;
	int32_t g1[] = { 2, 0, 0, 3, 2 }, g2[] = { 0, 2, 3 }, g3[] = { 1 }, g4[] = { 2 }, g5[] = { 1 };
	int32_t *mat;
	int64_t rows = 0, cols = 0;
	char npy_bin_fn[512];

	snprintf(ps4g_fn, sizeof(ps4g_fn), "%s/unit.ps4g", outdir);
	snprintf(npy_fn, sizeof(npy_fn), "%s/unit.npy", outdir);
	snprintf(npy_bin_fn, sizeof(npy_bin_fn), "%s/unit_binary.npy", outdir);
	snprintf(bed_fn, sizeof(bed_fn), "%s/unit_labels.bed", outdir);
	write_bed_fixture(bed_fn);

	// bin_size=100: pos 250,299 -> bin2 (same gameteSet {0,2,3}, aggregates to count=2);
	// pos 150 -> bin1 {1}; pos 210 -> bin2 but a DIFFERENT gameteSet {1} (tests that PS4G keeps
	// distinct gameteSets as separate rows while npy sums both into one (contig,bin) row);
	// ref_sid=1 pos 50 -> bin0 {2}, a different reference sequence entirely.
	rb3_ps4g_acc_add(acc, 0, 250, g1, 5);
	rb3_ps4g_acc_add(acc, 0, 299, g2, 3);
	rb3_ps4g_acc_add(acc, 0, 150, g3, 1);
	rb3_ps4g_acc_add(acc, 1, 50,  g4, 1);
	rb3_ps4g_acc_add(acc, 0, 210, g5, 1);
	CHECK(acc->n_ev == 5, "finalize: 5 events accumulated before aggregation");

	bed = rb3_bed_read(bed_fn, g, sid, "B73"); // no trailing '_': the real-world convention (see test_bed_contig_resolution)
	CHECK(bed != 0, "finalize: BED fixture parsed");
	CHECK(bed && bed->n_r == 2, "finalize: BED has 2 regions");

	rb3_ps4g_npy_finalize(acc, g, sid, bed, /*npy_binary=*/0, /*ridx=*/0, /*anchor_thresh=*/0, /*ref_gamete=*/-1, ps4g_fn, npy_fn, "unit-test-cmd");
	rb3_ps4g_npy_finalize(acc, g, sid, bed, /*npy_binary=*/1, /*ridx=*/0, /*anchor_thresh=*/0, /*ref_gamete=*/-1, 0, npy_bin_fn, "unit-test-cmd");

	// ---- PS4G file ----
	fp = fopen(ps4g_fn, "r");
	CHECK(fp != 0, "finalize: PS4G file was written");
	if (fp) {
		CHECK_STREQ(slurp_line(fp), "#PS4G", "ps4g: line 1");
		CHECK_STREQ(slurp_line(fp), "#version=2.0", "ps4g: line 2");
		CHECK_STREQ(slurp_line(fp), "#Command: unit-test-cmd", "ps4g: command line echoed verbatim");
		CHECK_STREQ(slurp_line(fp), "#TotalUniqueCounts: 5", "ps4g: total = sum of unique-row counts (1+1+2+1)");
		CHECK_STREQ(slurp_line(fp), "#gamete\tgameteIndex\tcount", "ps4g: gamete table header");
		CHECK_STREQ(slurp_line(fp), "#A\t0\t2", "ps4g: gamete A total count (in the {0,2,3} row, count 2)");
		CHECK_STREQ(slurp_line(fp), "#B\t1\t2", "ps4g: gamete B total count (bin1 row + bin2 {1} row)");
		CHECK_STREQ(slurp_line(fp), "#B73\t2\t3", "ps4g: gamete B73 total count (2 from {0,2,3} + 1 from ref1 {2})");
		CHECK_STREQ(slurp_line(fp), "#C\t3\t2", "ps4g: gamete C total count");
		CHECK_STREQ(slurp_line(fp), "gameteSet\trefContig\trefPosBinned\tcount", "ps4g: data section header");
		CHECK_STREQ(slurp_line(fp), "1\tchr1\t1\t1", "ps4g: row ref0/bin1 {B}");
		CHECK_STREQ(slurp_line(fp), "1\tchr1\t2\t1", "ps4g: row ref0/bin2 {B} (distinct gameteSet from the one below)");
		CHECK_STREQ(slurp_line(fp), "0,2,3\tchr1\t2\t2", "ps4g: row ref0/bin2 {A,B73,C} count 2");
		CHECK_STREQ(slurp_line(fp), "2\tscaf_1\t0\t1", "ps4g: row ref1/bin0 on a different contig, prefix stripped");
		CHECK(slurp_line(fp) == 0, "ps4g: no trailing rows");
		fclose(fp);
	}

	// ---- npy file (count mode): one row per PS4G row, i.e. per (contig,bin,gameteSet) --
	// NOT collapsed across the two distinct gameteSets ({B} and {A,B73,C}) that share bin2.
	// Losing that distinction would hide that {A,B73,C} co-occurred on the same reads.
	mat = read_npy_i4(npy_fn, &rows, &cols);
	CHECK(mat != 0, "finalize: npy file parsed");
	if (mat) {
		CHECK(rows == 4, "npy: 4 rows, matching the 4 PS4G rows exactly (gameteSets kept separate)");
		CHECK(cols == 4 + 2, "npy: n_gamete(4) + 2 label columns");
		// row0 = ref0/bin1 {B}count1 -> [0,1,0,0], unlabeled (pos=100 not in [200,300))
		CHECK(mat[0*cols+0]==0 && mat[0*cols+1]==1 && mat[0*cols+2]==0 && mat[0*cols+3]==0, "npy: row0 (bin1 {B}) counts");
		CHECK(mat[0*cols+4]==-1 && mat[0*cols+5]==-1, "npy: row0 unlabeled (-1,-1)");
		// row1 = ref0/bin2 {B}count1 -> [0,1,0,0]; labeled (pos=200 is in [200,300))
		CHECK(mat[1*cols+0]==0 && mat[1*cols+1]==1 && mat[1*cols+2]==0 && mat[1*cols+3]==0, "npy: row1 (bin2 {B}) counts, kept separate from row2's {A,B73,C}");
		CHECK(mat[1*cols+4]==0 && mat[1*cols+5]==2, "npy: row1 labeled (A=0, B73=2) from the chr1:200-300 BED region");
		// row2 = ref0/bin2 {A,B73,C}count2 -> [2,0,2,2]; same bin as row1, same label, different gameteSet
		CHECK(mat[2*cols+0]==2 && mat[2*cols+1]==0 && mat[2*cols+2]==2 && mat[2*cols+3]==2, "npy: row2 (bin2 {A,B73,C}) counts, not merged with row1");
		CHECK(mat[2*cols+4]==0 && mat[2*cols+5]==2, "npy: row2 labeled the same as row1 (same bin, independent of gameteSet)");
		// row3 = ref1/bin0 {B73}count1 -> [0,0,1,0]; homozygous C label
		CHECK(mat[3*cols+0]==0 && mat[3*cols+1]==0 && mat[3*cols+2]==1 && mat[3*cols+3]==0, "npy: row3 (ref1/bin0 {B73}) counts");
		CHECK(mat[3*cols+4]==3 && mat[3*cols+5]==3, "npy: row3 homozygous label (C,C) from a single-sample BED row");
		free(mat);
	}

	// ---- npy file (--npy-binary): identical nonzero pattern, but values clipped to 1 ----
	{
		int64_t rows2 = 0, cols2 = 0;
		int32_t *mat2 = read_npy_i4(npy_bin_fn, &rows2, &cols2);
		CHECK(mat2 != 0, "finalize: binary-mode npy file parsed");
		if (mat2) {
			CHECK(rows2 == 4 && cols2 == 6, "npy-binary: same shape as count mode");
			CHECK(mat2[0*cols2+1]==1, "npy-binary: row0 gamete B is 1 (was already count 1)");
			CHECK(mat2[2*cols2+0]==1 && mat2[2*cols2+1]==0 && mat2[2*cols2+2]==1 && mat2[2*cols2+3]==1,
				  "npy-binary: row2 {A,B73,C} clipped from count=2 to presence=1");
			CHECK(mat2[2*cols2+4]==0 && mat2[2*cols2+5]==2, "npy-binary: labels are unaffected by binary mode");
			free(mat2);
		}
	}

	// ---- companion TSVs ----
	{
		char aux_fn[600];
		snprintf(aux_fn, sizeof(aux_fn), "%s.bins.tsv", npy_fn);
		fp = fopen(aux_fn, "r");
		CHECK(fp != 0, "finalize: .bins.tsv written");
		if (fp) {
			CHECK_STREQ(slurp_line(fp), "row\tcontig\tbin", "bins.tsv: header");
			CHECK_STREQ(slurp_line(fp), "0\tchr1\t1", "bins.tsv: row 0");
			CHECK_STREQ(slurp_line(fp), "1\tchr1\t2", "bins.tsv: row 1");
			CHECK_STREQ(slurp_line(fp), "2\tchr1\t2", "bins.tsv: row 2 (bin repeats: a second gameteSet at the same bin)");
			CHECK_STREQ(slurp_line(fp), "3\tscaf_1\t0", "bins.tsv: row 3");
			fclose(fp);
		}
		snprintf(aux_fn, sizeof(aux_fn), "%s.gametes.tsv", npy_fn);
		fp = fopen(aux_fn, "r");
		CHECK(fp != 0, "finalize: .gametes.tsv written");
		if (fp) {
			CHECK_STREQ(slurp_line(fp), "gameteIndex\tsampleName", "gametes.tsv: header");
			CHECK_STREQ(slurp_line(fp), "0\tA", "gametes.tsv: row 0");
			CHECK_STREQ(slurp_line(fp), "1\tB", "gametes.tsv: row 1");
			CHECK_STREQ(slurp_line(fp), "2\tB73", "gametes.tsv: row 2");
			CHECK_STREQ(slurp_line(fp), "3\tC", "gametes.tsv: row 3");
			fclose(fp);
		}
	}

	rb3_bed_destroy(bed);
	rb3_ps4g_acc_destroy(acc);
	rb3_gtab_destroy(g);
	free_sid(sid);
}

// Regression test for a real bug: a raw name+strlen(ref_prefix) offset left a
// stray leading '_' ("_chr1") whenever --ref-prefix didn't itself include the
// separator -- exactly the convention this project's own scripts use
// (--ref-prefix=B73, not "B73_"). Two sequences share the bare contig name
// "chr1" here (a reference and a assembly) to also confirm the fix doesn't
// accidentally resolve a BED region to the wrong (non-reference) sequence.
static void test_bed_contig_resolution(const char *outdir)
{
	rb3_sid_t *sid = RB3_CALLOC(rb3_sid_t, 1);
	rb3_gtab_t *g;
	rb3_bed_t *bed;
	char bed_fn[512];
	FILE *fp;
	sid->n_seq = 2;
	sid->name = RB3_MALLOC(char*, 2);
	sid->len = RB3_CALLOC(int32_t, 2);
	sid->name[0] = rb3_strdup("B73_chr1");  // reference
	sid->name[1] = rb3_strdup("Oh43_chr1"); // assembly; same bare contig name "chr1"

	g = rb3_gtab_build(sid);
	snprintf(bed_fn, sizeof(bed_fn), "%s/bed_contig_resolution.bed", outdir);
	fp = fopen(bed_fn, "w");
	fprintf(fp, "chr1\t0\t100\tB73\n"); // bare contig name, no leading '_'
	fclose(fp);

	bed = rb3_bed_read(bed_fn, g, sid, "B73"); // no trailing '_'
	CHECK(bed != 0, "bed_resolve_contig: fixture parsed");
	CHECK(bed && bed->n_r == 1, "bed_resolve_contig: one region parsed");
	CHECK(bed && bed->n_r == 1 && bed->r[0].ref_sid == 0,
		  "bed_resolve_contig: 'chr1' with --ref-prefix=B73 (no '_') resolves to the reference sequence (0), not the assembly (1) sharing the same bare contig name, and not -1 (unresolved)");

	rb3_bed_destroy(bed);
	rb3_gtab_destroy(g);
	free(sid->name[0]); free(sid->name[1]); free(sid->name); free(sid->len); free(sid);
}

int main(int argc, char *argv[])
{
	const char *outdir = argc > 1? argv[1] : "test/output";
	rb3_sid_t *sid = make_sid();
	rb3_gtab_t *g = rb3_gtab_build(sid);

	test_gtab(sid, g);
	rb3_gtab_destroy(g);
	free_sid(sid);

	test_acc_dedup();
	test_finalize(outdir);
	test_bed_contig_resolution(outdir);

	fprintf(stderr, "test_ps4g: %d passed, %d failed\n", n_pass, n_fail);
	return n_fail == 0? 0 : 1;
}
