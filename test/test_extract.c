// Unit tests for coordinate-based extraction: the sampled inverse suffix array
// (isa.c: rb3_isa_from_ssa*, rb3_isa_extract, dump/restore, validation) and the
// reference->assembly inverse liftover projection (lift.c: rb3_lift_project_inv),
// plus a few hand-computed forward rb3_lift_project values on the same fixture.
//
// The FM-index is built IN MEMORY from known random sequences (the same
// libsais -> rld path `ropebwt3 build -d` uses), so every extracted base can be
// compared with the truth. Run with `make test` from the repository root.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "io.h"
#include "kalloc.h"
#include "ps4g.h"
#include "lift.h"
#include "isa.h"
#include "snap.h"

void rb3_bwtl_init(void) {} // misc.o's rb3_init() references this dawg.c symbol (see test_ps4g.c)

static int n_pass = 0, n_fail = 0;

#define CHECK(cond, msg) do { \
		if (cond) { ++n_pass; } \
		else { ++n_fail; fprintf(stderr, "FAIL: %s (%s:%d)\n", (msg), __FILE__, __LINE__); } \
	} while (0)

// Count-only variant for loops with thousands of cases: one PASS/FAIL per loop.
#define LOOP_CHECK(nbad, ntot, msg) do { \
		if ((nbad) == 0 && (ntot) > 0) ++n_pass; \
		else { ++n_fail; fprintf(stderr, "FAIL: %s: %ld of %ld cases wrong (%s:%d)\n", (msg), (long)(nbad), (long)(ntot), __FILE__, __LINE__); } \
	} while (0)

/*********************
 * in-memory fixture *
 *********************/

#define N_SEQ 6
// lengths cover: a 1-base sequence, sequences shorter than any sample spacing, and
// ones spanning many 2^shift buckets; seq 3 carries N runs (incl. at both ends).
static const int32_t seq_len[N_SEQ] = { 1, 7, 64, 1000, 5000, 3333 };
static char *seq_fwd[N_SEQ]; // ACGTN strings (truth), forward strand

static uint64_t rng_state = 11;
static uint32_t rng(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return (uint32_t)rng_state; }

static void make_seqs(void)
{
	int i, j;
	for (i = 0; i < N_SEQ; ++i) {
		seq_fwd[i] = RB3_MALLOC(char, seq_len[i] + 1);
		for (j = 0; j < seq_len[i]; ++j) seq_fwd[i][j] = "ACGT"[rng() & 3];
		seq_fwd[i][seq_len[i]] = 0;
	}
	// N runs in seq 3: at the start, the end, and two interior runs
	for (j = 0; j < 3; ++j) seq_fwd[3][j] = 'N';
	for (j = 400; j < 450; ++j) seq_fwd[3][j] = 'N';
	seq_fwd[3][700] = 'N';
	for (j = 995; j < 1000; ++j) seq_fwd[3][j] = 'N';
	// a repeat: copy 300 bases of seq 4 into seq 5 (multi-copy sequence)
	memcpy(seq_fwd[5] + 1000, seq_fwd[4] + 2000, 300);
}

static void revcomp_str(const char *s, int64_t n, char *out)
{
	int64_t i;
	for (i = 0; i < n; ++i) {
		char c = s[n - 1 - i];
		out[i] = c == 'A'? 'T' : c == 'C'? 'G' : c == 'G'? 'C' : c == 'T'? 'A' : 'N';
	}
	out[n] = 0;
}

// BWT string sid's truth: sid = 2j forward of seq j, 2j+1 its reverse complement
static void truth(int64_t sid, int64_t st, int64_t en, char *out)
{
	const char *s = seq_fwd[sid>>1];
	if ((sid & 1) == 0) { memcpy(out, s + st, en - st); out[en - st] = 0; }
	else {
		char *rc = RB3_MALLOC(char, seq_len[sid>>1] + 1);
		revcomp_str(s, seq_len[sid>>1], rc);
		memcpy(out, rc + st, en - st); out[en - st] = 0;
		free(rc);
	}
}

static void build_fmi(rb3_fmi_t *f)
{
	int64_t i, l = 0, tot = 0;
	char *s;
	rld_t *e;
	for (i = 0; i < N_SEQ; ++i) tot += 2 * (seq_len[i] + 1);
	s = RB3_MALLOC(char, tot);
	for (i = 0; i < N_SEQ; ++i) { // same layout as rb3_seq_read(..., is_for=1, is_rev=1)
		memcpy(s + l, seq_fwd[i], seq_len[i] + 1);
		rb3_char2nt6(seq_len[i], (uint8_t*)s + l);
		s[l + seq_len[i]] = 0;
		l += seq_len[i] + 1;
		memcpy(s + l, s + l - seq_len[i] - 1, seq_len[i] + 1);
		rb3_revcomp6(seq_len[i], (uint8_t*)s + l);
		l += seq_len[i] + 1;
	}
	rb3_build_sais(2 * N_SEQ, l, s, 1);
	e = rb3_enc_plain2rld(l, (uint8_t*)s, 3);
	free(s);
	memset(f, 0, sizeof(*f));
	rb3_fmi_init(f, e, 0);
	f->sid = 0;
}

static void nt6_to_str(const uint8_t *b, int64_t n, char *out)
{
	int64_t i;
	for (i = 0; i < n; ++i) out[i] = "$ACGTN"[b[i]];
	out[n] = 0;
}

// extract [st,en) of sid and compare with truth; returns 1 if equal
static int extract_ok(const rb3_fmi_t *f, const rb3_isa_t *isa, int64_t sid, int64_t st, int64_t en)
{
	uint8_t *b = RB3_MALLOC(uint8_t, en - st + 1);
	char *got = RB3_MALLOC(char, en - st + 1), *want = RB3_MALLOC(char, en - st + 1);
	int64_t r = rb3_isa_extract(f, isa, sid, st, en, b);
	int ok = 0;
	if (r >= en - st) {
		nt6_to_str(b, en - st, got);
		truth(sid, st, en, want);
		ok = strcmp(got, want) == 0;
	}
	free(b); free(got); free(want);
	return ok;
}

/*********
 * tests *
 *********/

static void test_fixture_numbering(const rb3_fmi_t *f)
{
	kstring_t s = {0,0,0};
	int64_t sid, nbad = 0;
	char *want = RB3_MALLOC(char, 6000);
	CHECK(f->acc[1] == 2 * N_SEQ, "fixture: 2 BWT strings (forward + reverse complement) per sequence");
	for (sid = 0; sid < 2 * N_SEQ; ++sid) { // `ropebwt3 get sid` == forward (even) / revcomp (odd) of seq sid>>1
		rb3_fmi_retrieve(f, sid, &s);
		truth(sid, 0, seq_len[sid>>1], want);
		if (strcmp(s.s, want) != 0) ++nbad;
	}
	LOOP_CHECK(nbad, 2 * N_SEQ, "numbering: get(sid) = seq sid>>1, forward for even sid, reverse complement for odd");
	free(s.s); free(want);
}

// every entry, LF-walked to the string start, must take exactly pos steps and land on
// the sentinel that r2i maps back to sid
static void test_entries_walk(const rb3_fmi_t *f, const rb3_ssa_t *sa, const rb3_isa_t *isa, const char *tag)
{
	int64_t sid, i, nbad = 0, nbad_sort = 0, nbad_sent = 0, n = 0;
	char msg[256];
	for (sid = 0; sid < isa->m; ++sid) {
		int64_t lo = isa->off[sid], hi = isa->off[sid+1];
		if (hi > lo) {
			if (isa->pos[hi-1] != (uint32_t)seq_len[sid>>1] || isa->row[hi-1] != sid) ++nbad_sent;
		}
		for (i = lo; i < hi; ++i) {
			int64_t k = isa->row[i], steps = 0, ok[RB3_ASIZE];
			int c;
			if (i > lo && isa->pos[i] <= isa->pos[i-1]) ++nbad_sort;
			while ((c = rb3_fmi_rank1a(f, k, ok)) > 0) k = f->acc[c] + ok[c], ++steps;
			// c == 0: k is the row of the sentinel preceding the string
			k = f->acc[0] + ok[0];
			if (steps != (int64_t)isa->pos[i] || (int64_t)sa->r2i[k] != sid) ++nbad;
			++n;
		}
	}
	snprintf(msg, sizeof msg, "%s: every entry LF-walks exactly pos steps to the start of its own string", tag); LOOP_CHECK(nbad, n, msg);
	snprintf(msg, sizeof msg, "%s: entries strictly increasing in pos within a string", tag); LOOP_CHECK(nbad_sort, n, msg);
	snprintf(msg, sizeof msg, "%s: last entry of each sampled string is the sentinel (pos=len, row=sid)", tag); CHECK(nbad_sent == 0, msg);
}

static void test_entry_density(const rb3_isa_t *isa, int shift)
{
	int64_t sid, nbad = 0;
	char msg[128];
	for (sid = 0; sid < isa->m; sid += 2) // at most one SSA entry per bucket + the sentinel
		if (isa->off[sid+1] - isa->off[sid] > (seq_len[sid>>1] >> shift) + 2) ++nbad;
	snprintf(msg, sizeof msg, "shift %d: at most one entry per 2^shift bucket (+ sentinel)", shift);
	CHECK(nbad == 0, msg);
}

static void test_extract_cases(const rb3_fmi_t *f, const rb3_isa_t *isa, int both, const char *tag)
{
	int64_t sid, st, en, i, nbad, n;
	char msg[256];
	int64_t step = both? 1 : 2;
	// exhaustive [st,en) for the short sequences
	nbad = n = 0;
	for (sid = 0; sid < 2 * N_SEQ; sid += step) {
		if (seq_len[sid>>1] > 64) continue;
		for (st = 0; st < seq_len[sid>>1]; ++st)
			for (en = st + 1; en <= seq_len[sid>>1]; ++en)
				nbad += !extract_ok(f, isa, sid, st, en), ++n;
	}
	snprintf(msg, sizeof msg, "%s: every [st,en) of the 1/7/64-bp sequences", tag); LOOP_CHECK(nbad, n, msg);
	// starts, ends, whole sequences, length 1 at both ends
	nbad = n = 0;
	for (sid = 0; sid < 2 * N_SEQ; sid += step) {
		int64_t L = seq_len[sid>>1];
		nbad += !extract_ok(f, isa, sid, 0, 1), ++n;
		nbad += !extract_ok(f, isa, sid, L - 1, L), ++n;
		nbad += !extract_ok(f, isa, sid, 0, L), ++n;
		if (L > 10) nbad += !extract_ok(f, isa, sid, 0, 10) + !extract_ok(f, isa, sid, L - 10, L), n += 2;
	}
	snprintf(msg, sizeof msg, "%s: sequence start, end, length 1 at both ends, whole sequence", tag); LOOP_CHECK(nbad, n, msg);
	// exactly at / across every entry position
	nbad = n = 0;
	for (sid = 0; sid < 2 * N_SEQ; sid += step) {
		int64_t L = seq_len[sid>>1];
		for (i = isa->off[sid]; i < isa->off[sid+1]; ++i) {
			int64_t p = isa->pos[i];
			if (p < L) nbad += !extract_ok(f, isa, sid, p, p + 1), ++n;          // the base AT the entry
			if (p > 0) nbad += !extract_ok(f, isa, sid, p - 1, p), ++n;          // ends exactly at the entry
			if (p > 0 && p < L) nbad += !extract_ok(f, isa, sid, p - 1, p + 1), ++n; // straddles it
			if (p >= 5 && p + 5 <= L) nbad += !extract_ok(f, isa, sid, p - 5, p + 5), ++n;
		}
	}
	snprintf(msg, sizeof msg, "%s: intervals at, ending at, and straddling every entry position", tag); LOOP_CHECK(nbad, n, msg);
	// random intervals, incl. long ones
	nbad = n = 0;
	for (i = 0; i < 3000; ++i) {
		int64_t L;
		sid = (rng() % (2 * N_SEQ)) & (both? ~0LL : ~1LL);
		L = seq_len[sid>>1];
		st = rng() % L;
		en = st + 1 + rng() % (i % 3 == 0? L - st : (L - st < 50? L - st : 50));
		nbad += !extract_ok(f, isa, sid, st, en), ++n;
	}
	snprintf(msg, sizeof msg, "%s: 3000 random intervals (1/3 up to the sequence end)", tag); LOOP_CHECK(nbad, n, msg);
	// N handling: the N runs of seq 3 come back as N
	{
		uint8_t b[60]; char s[61];
		int64_t r = rb3_isa_extract(f, isa, 6, 395, 455, b);
		nt6_to_str(b, 60, s);
		snprintf(msg, sizeof msg, "%s: an interior N run is extracted as N", tag);
		CHECK(r >= 60 && strncmp(s + 5, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", 50) == 0 && s[4] != 'N' && s[55] != 'N', msg);
		r = rb3_isa_extract(f, isa, 6, 0, 3, b);
		snprintf(msg, sizeof msg, "%s: leading N run", tag);
		CHECK(r >= 3 && b[0] == 5 && b[1] == 5 && b[2] == 5, msg);
	}
}

static void test_extract_errors(const rb3_fmi_t *f, const rb3_isa_t *isa_fwd)
{
	uint8_t b[16];
	CHECK(rb3_isa_extract(f, isa_fwd, 8, -1, 5, b) == RB3_ISA_ERR_RANGE, "error: negative start");
	CHECK(rb3_isa_extract(f, isa_fwd, 8, 0, 5001, b) == RB3_ISA_ERR_RANGE, "error: end past the sequence end");
	CHECK(rb3_isa_extract(f, isa_fwd, 8, 5, 5, b) == RB3_ISA_ERR_RANGE, "error: empty interval");
	CHECK(rb3_isa_extract(f, isa_fwd, 8, 6, 5, b) == RB3_ISA_ERR_RANGE, "error: start > end");
	CHECK(rb3_isa_extract(f, isa_fwd, 2 * N_SEQ, 0, 1, b) == RB3_ISA_ERR_SID, "error: sid out of range");
	CHECK(rb3_isa_extract(f, isa_fwd, -1, 0, 1, b) == RB3_ISA_ERR_SID, "error: negative sid");
	CHECK(rb3_isa_extract(f, isa_fwd, 9, 0, 1, b) == RB3_ISA_ERR_SID, "error: reverse string not sampled in a forward-only ISA");
	CHECK(rb3_isa_seq_len(isa_fwd, 8) == 5000 && rb3_isa_seq_len(isa_fwd, 9) == -1, "seq_len: from the sentinel entry; -1 if unsampled");
}

static int isa_equal(const rb3_isa_t *a, const rb3_isa_t *b)
{
	if (a == 0 || b == 0) return 0;
	if (a->shift != b->shift || a->ss != b->ss || a->strands != b->strands || a->m != b->m || a->n_ssa != b->n_ssa || a->n != b->n) return 0;
	return memcmp(a->off, b->off, (a->m + 1) * 8) == 0 && memcmp(a->pos, b->pos, a->n * 4) == 0 && memcmp(a->row, b->row, a->n * 8) == 0;
}

static void test_files(const rb3_fmi_t *f, const rb3_ssa_t *sa, const rb3_isa_t *isa, const char *outdir)
{
	char fn_ssa[512], fn_isa[512], fn_bad[512];
	int32_t badlen[N_SEQ];
	rb3_isa_t *x;
	FILE *fp;
	int verbose0 = rb3_verbose;
	snprintf(fn_ssa, sizeof fn_ssa, "%s/extract_fixture.ssa", outdir);
	snprintf(fn_isa, sizeof fn_isa, "%s/extract_fixture.isa", outdir);
	snprintf(fn_bad, sizeof fn_bad, "%s/extract_fixture_bad.ssa", outdir);
	CHECK(rb3_ssa_dump(sa, fn_ssa) == 0, "files: wrote the fixture .ssa");
	x = rb3_isa_from_ssa_file(fn_ssa, N_SEQ, seq_len, isa->shift, isa->strands);
	CHECK(isa_equal(x, isa), "files: streaming the .ssa file gives the same entries as the in-memory SSA");
	rb3_isa_destroy(x);
	CHECK(rb3_isa_dump(isa, fn_isa) == 0, "files: wrote the .isa");
	x = rb3_isa_restore(fn_isa);
	CHECK(isa_equal(x, isa), "files: .isa dump/restore round trip");
	CHECK(rb3_isa_check_fmi(x, f) == 0, "files: restored .isa validates against the FMD");
	rb3_isa_destroy(x);
	rb3_verbose = 0; // the rejections below print an ERROR line each; they are expected
	// wrong lengths (a .len.gz from another index): samples fall outside -> rejected
	memcpy(badlen, seq_len, sizeof badlen);
	badlen[4] = 100;
	CHECK(rb3_isa_from_ssa_file(fn_ssa, N_SEQ, badlen, isa->shift, isa->strands) == 0, "files: .len.gz shorter than the SSA samples is rejected");
	CHECK(rb3_isa_from_ssa_file(fn_ssa, N_SEQ - 1, seq_len, isa->shift, isa->strands) == 0, "files: #sequences != m/2 is rejected");
	CHECK(rb3_isa_from_ssa_file(fn_ssa, N_SEQ, seq_len, 31, isa->strands) == 0, "files: shift out of range is rejected");
	CHECK(rb3_isa_from_ssa_file(fn_ssa, N_SEQ, seq_len, 4, 0) == 0, "files: no strand selected is rejected");
	CHECK(rb3_isa_from_ssa_file("/nonexistent/x.ssa", N_SEQ, seq_len, 4, 1) == 0, "files: missing .ssa -> NULL");
	CHECK(rb3_isa_restore("/nonexistent/x.isa") == 0, "files: missing .isa -> NULL");
	// truncated .ssa
	{
		FILE *in = fopen(fn_ssa, "rb");
		char buf[1<<16];
		size_t n = fread(buf, 1, sizeof buf, in);
		fclose(in);
		fp = fopen(fn_bad, "wb"); fwrite(buf, 1, n - 16, fp); fclose(fp);
		CHECK(rb3_isa_from_ssa_file(fn_bad, N_SEQ, seq_len, 4, 1) == 0, "files: truncated .ssa is rejected");
		fp = fopen(fn_bad, "wb"); fwrite("XXXX", 1, 4, fp); fwrite(buf + 4, 1, n - 4, fp); fclose(fp);
		CHECK(rb3_isa_from_ssa_file(fn_bad, N_SEQ, seq_len, 4, 1) == 0, "files: wrong .ssa magic is rejected");
		CHECK(rb3_isa_restore(fn_bad) == 0, "files: a non-.isa file is rejected by rb3_isa_restore");
		in = fopen(fn_isa, "rb"); n = fread(buf, 1, sizeof buf, in); fclose(in);
		fp = fopen(fn_bad, "wb"); fwrite(buf, 1, n - 8, fp); fclose(fp);
		CHECK(rb3_isa_restore(fn_bad) == 0, "files: truncated .isa is rejected");
	}
	// an ISA from a different SSA sampling rate does not validate against this FMD
	{
		rb3_isa_t y = *isa;
		y.n_ssa += 1;
		CHECK(rb3_isa_check_fmi(&y, f) == -1, "validate: n_ssa mismatch with the FMD is detected");
		y = *isa; y.m += 2;
		CHECK(rb3_isa_check_fmi(&y, f) == -1, "validate: #strings mismatch with the FMD is detected");
	}
	rb3_verbose = verbose0;
}

/**********************************
 * inverse (reference->assembly) lift *
 **********************************/

typedef struct { int64_t cpos, rpos; int32_t csid, rsid; } fixture_pt_t;

static int fpt_cmp(const void *a, const void *b)
{
	const fixture_pt_t *x = (const fixture_pt_t*)a, *y = (const fixture_pt_t*)b;
	if (x->csid != y->csid) return x->csid < y->csid? -1 : 1;
	return x->cpos < y->cpos? -1 : x->cpos > y->cpos;
}

// sequences: 0 B73_chr1, 1 B73_chr2, 2 A_chr1, 3 B_chr1, 4 C_chr1, 5 C_ctg2
//   A_chr1: collinear, cpos = rpos + 50 for rpos 1000..10000, then a 500-bp
//           insertion in A: cpos = rpos + 550 for rpos 11000..20000
//   B_chr1: inverted, cpos = 30000 - rpos for rpos 1000..20000
//   C:      rpos 1000..10000 on C_chr1 (cpos = rpos), 11000..20000 on C_ctg2 (cpos = rpos - 11000)
//   B73_chr2 has no anchors at all
static const char *write_inv_fixture(const char *outdir)
{
	static char fn[512];
	fixture_pt_t pt[60];
	int64_t off[7] = {0}, n_seq = 6, n_pt = 0, r, i;
	FILE *fp;
	for (r = 1000; r <= 20000; r += 1000) {
		pt[n_pt].csid = 2, pt[n_pt].rsid = 0, pt[n_pt].rpos = r, pt[n_pt++].cpos = r + (r <= 10000? 50 : 550);
		pt[n_pt].csid = 3, pt[n_pt].rsid = 0, pt[n_pt].rpos = r, pt[n_pt++].cpos = 30000 - r;
		pt[n_pt].csid = r <= 10000? 4 : 5, pt[n_pt].rsid = 0, pt[n_pt].rpos = r, pt[n_pt++].cpos = r <= 10000? r : r - 11000;
	}
	qsort(pt, n_pt, sizeof(fixture_pt_t), fpt_cmp);
	for (i = 0; i < n_pt; ++i) off[pt[i].csid + 1]++;
	for (i = 0; i < n_seq; ++i) off[i + 1] += off[i];
	snprintf(fn, sizeof fn, "%s/extract_inv.lift", outdir);
	fp = fopen(fn, "wb");
	fwrite("LIFT\1", 1, 5, fp);
	fwrite(&n_seq, 8, 1, fp); fwrite(&n_pt, 8, 1, fp);
	fwrite(off, 8, n_seq + 1, fp);
	fwrite(pt, sizeof(fixture_pt_t), n_pt, fp);
	fclose(fp);
	return fn;
}

static void test_inverse_lift(const char *outdir)
{
	const char *names[] = { "B73_chr1", "B73_chr2", "A_chr1", "B_chr1", "C_chr1", "C_ctg2" };
	rb3_sid_t sid;
	rb3_gtab_t *g;
	rb3_lift_t *lf;
	rb3_lift_ridx_t *ridx;
	int32_t lens[6] = { 30000, 30000, 30000, 30000, 30000, 30000 }, gA, gB, gC;
	int64_t cs, cp, rs, rp;
	int sg, me;
	void *km = km_init();

	sid.n_seq = 6, sid.name = (char**)names, sid.len = lens, sid.tot_len = 0;
	g = rb3_gtab_build(&sid);
	gA = rb3_gtab_find(g, "A"), gB = rb3_gtab_find(g, "B"), gC = rb3_gtab_find(g, "C");
	lf = rb3_lift_restore(write_inv_fixture(outdir));
	CHECK(lf != 0, "inv lift: fixture restored");
	ridx = rb3_lift_ridx_build(lf, g->sid2g, g->n_gamete);

	// --- local projection on (refmap's --lift-local-win=20000 --lift-local-tol=1000) ---
	rb3_lift_set_local(20000, 1000);
	CHECK(rb3_lift_project_inv(ridx, km, gA, 0, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cs == 2 && cp == 1550 && sg == 1 && me == 1,
		  "inv local: bracketed collinear pair interpolates (A rpos 1500 -> A_chr1:1550, +)");
	CHECK(rb3_lift_project_inv(ridx, km, gA, 0, 1000, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cs == 2 && cp == 1050,
		  "inv local: a position exactly on an anchor maps to that anchor");
	CHECK(rb3_lift_project_inv(ridx, km, gA, 0, 10700, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cp == 10750 && me == 1,
		  "inv local: a 500-bp insertion within --lift-local-tol still brackets (interpolates on the left anchor's offset)");
	CHECK(rb3_lift_project_inv(ridx, km, gB, 0, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cs == 3 && cp == 28500 && sg == -1,
		  "inv local: an inverted founder projects with slope -1 (B rpos 1500 -> B_chr1:28500, -)");
	CHECK(rb3_lift_project_inv(ridx, km, gC, 0, 5500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cs == 4 && cp == 5500,
		  "inv local: founder C, left contig");
	CHECK(rb3_lift_project_inv(ridx, km, gC, 0, 15500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cs == 5 && cp == 4500,
		  "inv local: founder C, right contig");
	CHECK(rb3_lift_project_inv(ridx, km, gC, 0, 10500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cs == 4 && cp == 10500 && me == 2,
		  "inv local: flanks on different contigs are not collinear -> one-sided from the nearer (tie: left) anchor");
	rb3_lift_set_local(20000, 100);
	CHECK(rb3_lift_project_inv(ridx, km, gA, 0, 10700, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cp == 11250 && me == 2,
		  "inv local: with tol 100 the insertion breaks the bracket -> one-sided from the nearer right anchor (11550-300)");
	CHECK(rb3_lift_project_inv(ridx, km, gB, 0, 1500, 500000, 200000, 4, 31, &cs, &cp, &sg, &me) == 1 && cp == 28530 && sg == -1,
		  "inv local: on a slope -1 locus the anchor k-mer offset k-1 (=30) is added back");
	CHECK(rb3_lift_project_inv(ridx, km, gA, 0, 1500, 500000, 200000, 4, 31, &cs, &cp, &sg, &me) == 1 && cp == 1550 && sg == 1,
		  "inv local: no k-1 shift on a slope +1 locus");
	rb3_lift_set_local(20000, 3000); // tol > anchor gap: both slopes fit within tol
	CHECK(rb3_lift_project_inv(ridx, km, gB, 0, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cp == 28500 && sg == -1 && me == 1,
		  "inv local: with a tolerance admitting both slopes, the better-fitting one (-1) is taken");
	CHECK(rb3_lift_project(lf, km, 3, 28500, 500000, 200000, 4, &rs, &rp) == 1 && rs == 0 && rp == 2000,
		  "fwd local (unchanged refmap behaviour): slope +1 tested first, interpolation clamped to the bracket -> 2000");
	rb3_lift_set_local(20000, 100);
	// --- no anchors / bad input ---
	CHECK(rb3_lift_project_inv(ridx, km, gA, 1, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 0, "inv: no anchors on the reference contig -> 0");
	CHECK(rb3_lift_project_inv(ridx, km, -1, 0, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 0, "inv: bad gamete -> 0");
	CHECK(rb3_lift_project_inv(ridx, km, gA, 7, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 0, "inv: bad rsid -> 0");
	CHECK(rb3_lift_project_inv(0, km, gA, 0, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 0, "inv: NULL ridx -> 0");
	// --- local off: the wide-window median (refmap's default) ---
	rb3_lift_set_local(0, 1000);
	CHECK(rb3_lift_project_inv(ridx, km, gA, 0, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cs == 2 && cp == 2050 && sg == 1 && me == 3,
		  "inv median: offset = upper median of {50 x10, 550 x10} = 550 -> 2050 (the coarse wide-window fit)");
	CHECK(rb3_lift_project_inv(ridx, km, gA, 0, 1500, 600, 200000, 4, 1, &cs, &cp, &sg, &me) == 0,
		  "inv median: fewer than min_support anchors in the window -> 0");
	CHECK(rb3_lift_project_inv(ridx, km, gB, 0, 1500, 500000, 200000, 4, 1, &cs, &cp, &sg, &me) == 1 && cp == 28500 && sg == -1 && me == 3,
		  "inv median: inverted founder, slope -1 fit");
	CHECK(rb3_lift_project_inv(ridx, km, gA, 0, 1500, 500000, 100, 4, 1, &cs, &cp, &sg, &me) == 0,
		  "inv median: residual MAD above max_mad -> 0 (never a guess)");
	// --- the forward projector still gives the same answers on the same points ---
	CHECK(rb3_lift_project(lf, km, 2, 1550, 500000, 200000, 4, &rs, &rp) == 1 && rs == 0 && rp == 1500,
		  "fwd median: A_chr1:1550 -> offset = upper median of {-550 x10, -50 x10} = -50 -> B73_chr1:1500");
	rb3_lift_set_local(20000, 1000);
	CHECK(rb3_lift_project(lf, km, 2, 1550, 500000, 200000, 4, &rs, &rp) == 1 && rs == 0 && rp == 1500,
		  "fwd local: A_chr1:1550 -> B73_chr1:1500");
	CHECK(rb3_lift_project(lf, km, 3, 28500, 500000, 200000, 4, &rs, &rp) == 1 && rs == 0 && rp == 1500,
		  "fwd local: inverted B_chr1:28500 -> B73_chr1:1500");
	rb3_lift_set_local(0, 1000);

	rb3_lift_ridx_destroy(ridx);
	rb3_lift_destroy(lf);
	rb3_gtab_destroy(g);
	km_destroy(km);
}


/*****************************
 * snap: the flank aligner   *
 *****************************/

static void rand_bases(char *s, int n) { int i; for (i = 0; i < n; ++i) s[i] = "ACGT"[rng() & 3]; s[n] = 0; }

static int snap_str(const char *q, const char *t, const rb3_snap_opt_t *o, rb3_snap_res_t *r)
{
	int32_t n = (int32_t)strlen(q), m = (int32_t)strlen(t);
	uint8_t *qq = RB3_MALLOC(uint8_t, n + 1), *tt = RB3_MALLOC(uint8_t, m + 1);
	int st;
	memcpy(qq, q, n); memcpy(tt, t, m);
	rb3_char2nt6(n, qq); rb3_char2nt6(m, tt);
	st = rb3_snap_align(qq, n, tt, m, o, r);
	free(qq); free(tt);
	return st;
}

static char other_base(char c) { return c == 'A'? 'C' : 'A'; }
static char base_not_in(const char *s, int n) { const char *b = "ACGT"; int i, k; for (k = 0; k < 4; ++k) { for (i = 0; i < n && s[i] != b[k]; ++i) {} if (i == n) return b[k]; } return 'N'; }

static void test_snap_align(void)
{
	rb3_snap_opt_t o;
	rb3_snap_res_t r;
	char fl[65], L[101], R[101], t[1000], q[400], tmp[200];
	rb3_snap_opt_init(&o);
	rand_bases(fl, 64); rand_bases(L, 100); rand_bases(R, 100);
	// 1. exact
	sprintf(t, "%s%s%s", L, fl, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_first == 100 && r.t_last == 163 && r.q_first == 0 && r.q_last == 63 && r.identity == 1.0,
		  "snap: exact flank found at its offset");
	// 2. SNPs, incl. the very first and last flank base (present, so they still define the ends)
	strcpy(tmp, fl); tmp[0] = other_base(tmp[0]); tmp[20] = other_base(tmp[20]); tmp[63] = other_base(tmp[63]);
	sprintf(t, "%s%s%s", L, tmp, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_first == 100 && r.t_last == 163 && r.n_match == 61,
		  "snap: SNPs (incl. at both flank ends) do not move the ends");
	// 3. small deletion inside the flank (founder lacks flank[30,33))
	sprintf(t, "%s%.30s%s%s", L, fl, fl + 33, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_first == 100 && r.t_last == 100 + 60,
		  "snap: 3-bp deletion inside the flank (end shifts by -3)");
	// 4. small insertion inside the flank
	sprintf(t, "%s%.30sGTCA%s%s", L, fl, fl + 30, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_first == 100 && r.t_last == 167,
		  "snap: 4-bp insertion inside the flank (end shifts by +4)");
	// 5. the end base deleted: flank[0,3) absent; start = founder base of flank[3]
	strcpy(tmp, L); // no accidental slide of the deletion into the left flank (an equally valid placement)
	{ int i; char c = base_not_in(fl, 3); for (i = 95; i < 100; ++i) tmp[i] = c; }
	sprintf(t, "%s%s%s", tmp, fl + 3, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.q_first == 3 && r.t_first == 100,
		  "snap: deleted end bases -> the first present flank base (q_first=3) at its founder base");
	strcpy(tmp, R);
	{ int i; char c = base_not_in(fl + 61, 3); for (i = 0; i < 5; ++i) tmp[i] = c; }
	sprintf(t, "%s%.61s%s", L, fl, tmp);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.q_last == 60 && r.t_last == 160,
		  "snap: deleted last bases -> the last present flank base");
	// 6. insertion adjacent to the end is never included; an equal-base one maps to the own copy
	sprintf(t, "%sGGTTC%s%s", L, fl, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_first == 105, "snap: insertion right before the start is excluded");
	tmp[0] = fl[0], tmp[1] = 0;
	sprintf(t, "%s%s%s%s", L, tmp, fl, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_first == 101,
		  "snap: an equal-base insertion before the start: s maps to its own (rightmost) copy");
	tmp[0] = fl[63], tmp[1] = 0;
	sprintf(t, "%s%s%s%s", L, fl, tmp, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_last == 163 && r.t_end == 164,
		  "snap: an equal-base insertion after the end is left out (leftmost end)");
	// 7. repeat in the window -> ambiguous; a diverged second copy is not
	sprintf(t, "%s%s%.50s%s%s", L, fl, R, fl, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_AMBIG && r.second == r.score, "snap: two identical copies -> ambiguous");
	strcpy(tmp, fl); tmp[10] = other_base(tmp[10]); tmp[30] = other_base(tmp[30]); tmp[50] = other_base(tmp[50]);
	sprintf(t, "%s%s%.50s%s%s", L, fl, R, tmp, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_first == 100 && r.second <= r.score - 10,
		  "snap: a copy with 3 mismatches scores well below the best -> unique");
	tmp[50] = fl[50], tmp[30] = fl[30]; // one mismatch: best-5, within the margin
	sprintf(t, "%s%s%.50s%s%s", L, fl, R, tmp, R);
	CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_AMBIG, "snap: a copy with 1 mismatch is within --snap-margin -> ambiguous");
	// 8. too short
	CHECK(snap_str("ACGTACGTAC", t, &o, &r) == RB3_SNAP_SHORT, "snap: flank shorter than min_aligned -> short");
	CHECK(snap_str(fl, "ACGTACGTACGT", &o, &r) == RB3_SNAP_SHORT, "snap: window shorter than the required aligned length -> short");
	{ // a 20-bp flank (e.g. a short interval) at the very start of a contig window
		char f20[21]; memcpy(f20, fl, 20); f20[20] = 0;
		sprintf(t, "%s%s", fl, R);
		CHECK(snap_str(f20, t, &o, &r) == RB3_SNAP_OK && r.t_first == 0 && r.t_last == 19, "snap: short flank at the window (contig) start");
	}
	// 9. low identity: unrelated sequence
	rand_bases(tmp, 64);
	sprintf(t, "%s%s", L, R);
	CHECK(snap_str(tmp, t, &o, &r) == RB3_SNAP_LOWID, "snap: unrelated flank -> low identity");
	{ // 25% divergence -> below 0.85
		strcpy(tmp, fl);
		{ int i; for (i = 0; i < 64; i += 4) tmp[i] = other_base(tmp[i]); }
		sprintf(t, "%s%s%s", L, tmp, R);
		CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_LOWID, "snap: 25% mismatches -> low identity (fallback expected)");
		o.min_id = 0.7;
		CHECK(snap_str(fl, t, &o, &r) == RB3_SNAP_OK && r.t_first == 100, "snap: ... accepted with --snap-min-id 0.7");
		o.min_id = 0.85;
	}
	// 10. N in the flank never matches but does not break the placement
	strcpy(tmp, fl); tmp[40] = 'N';
	sprintf(t, "%s%s%s", L, fl, R);
	CHECK(snap_str(tmp, t, &o, &r) == RB3_SNAP_OK && r.t_first == 100 && r.n_match == 63, "snap: an N in the flank");
	// 11. clip5: a long deleted prefix (the walk-inward final alignment)
	rand_bases(q, 200); strcat(q, fl);
	sprintf(t, "%s%s%s", L, fl, R);
	{
		rb3_snap_opt_t c = o;
		c.clip5 = c.gap_open + c.gap_ext, c.min_frac = 0, c.min_aligned = 32;
		CHECK(snap_str(q, t, &c, &r) == RB3_SNAP_OK && r.q_first == 200 && r.t_first == 100,
			  "snap clip5: a 200-bp deleted prefix is clipped, start = founder base of the first present base");
		CHECK(snap_str(q, t, &o, &r) == RB3_SNAP_LOWID, "snap: without clip5 the same query needs half of it aligned -> low identity");
	}
	// 13. the k-mer prefilter on wide windows (m > 8n) gives the full DP's answer
	{
		int it, nbad = 0, n = 0;
		rb3_snap_opt_t full = o;
		rb3_snap_res_t r2;
		static char big[6000], v[200];
		full.no_seed = 1;
		for (it = 0; it < 300; ++it) {
			int pos = 200 + rng() % 5000, nd = rng() % 8, d, kind = rng() % 4;
			rand_bases(big, 5800);
			strcpy(v, fl);
			for (d = 0; d < nd; ++d) v[rng() % 64] = "ACGT"[rng() & 3];         // 0-7 random substitutions
			if (kind == 1) memmove(v + 30, v + 33, strlen(v + 33) + 1);          // 3-bp deletion
			if (kind == 2) { memmove(v + 33, v + 30, strlen(v + 30) + 1); memcpy(v + 30, "GTA", 3); } // insertion
			memcpy(big + pos, v, strlen(v));
			if (kind == 3) memcpy(big + (pos + 1500) % 5600, fl, 64);             // a second copy
			if (snap_str(fl, big, &o, &r) != snap_str(fl, big, &full, &r2)) ++nbad;
			else if (r.status == RB3_SNAP_OK && (r.t_first != r2.t_first || r.t_last != r2.t_last || r.score != r2.score)) ++nbad;
			++n;
		}
		LOOP_CHECK(nbad, n, "snap: seeded alignment on a 5.8 kb window == the full DP (300 random cases)");
	}
	// 12. identical flank on a reverse-complemented window (what extract does on '-' loci)
	{
		char rc[300]; int n;
		sprintf(t, "%s%s%s", L, fl, R);
		n = (int)strlen(t);
		revcomp_str(t, n, rc);   // the founder stores the locus inverted ...
		revcomp_str(rc, n, tmp); // ... extract revcomps the window back into reference orientation
		CHECK(snap_str(fl, tmp, &o, &r) == RB3_SNAP_OK && r.t_first == 100, "snap: reverse-complemented window, re-oriented");
	}
}

int main(int argc, char *argv[])
{
	const char *outdir = argc > 1? argv[1] : "test/output";
	rb3_fmi_t f;
	rb3_ssa_t *sa;
	int shifts[] = { 0, 2, 5, 12 }, si;

	make_seqs();
	build_fmi(&f);
	test_fixture_numbering(&f);
	sa = rb3_ssa_gen(&f, 3, 1); // one SA sample per 8 rows: many entries, irregular gaps
	CHECK(sa != 0 && sa->m == 2 * N_SEQ, "ssa: generated for the in-memory fixture");
	for (si = 0; si < 4; ++si) {
		char tag[64];
		rb3_isa_t *both = rb3_isa_from_ssa(sa, N_SEQ, seq_len, shifts[si], RB3_ISA_FWD | RB3_ISA_REV);
		rb3_isa_t *fwd = rb3_isa_from_ssa(sa, N_SEQ, seq_len, shifts[si], RB3_ISA_FWD);
		snprintf(tag, sizeof tag, "shift %d both strands", shifts[si]);
		CHECK(both != 0 && fwd != 0, "isa: built from the in-memory SSA");
		CHECK(rb3_isa_check_fmi(both, &f) == 0, "isa: validates against its FMD");
		test_entries_walk(&f, sa, both, tag);
		test_entry_density(both, shifts[si]);
		test_extract_cases(&f, both, 1, tag);
		snprintf(tag, sizeof tag, "shift %d forward only", shifts[si]);
		test_extract_cases(&f, fwd, 0, tag);
		test_extract_errors(&f, fwd);
		if (shifts[si] == 5) test_files(&f, sa, both, outdir);
		rb3_isa_destroy(both); rb3_isa_destroy(fwd);
	}
	{ // shift 0 with a dense SSA: every sampled SA position becomes an entry
		rb3_isa_t *x = rb3_isa_from_ssa(sa, N_SEQ, seq_len, 0, RB3_ISA_FWD | RB3_ISA_REV);
		CHECK(x->n == sa->n_ssa + 2 * N_SEQ, "shift 0: every SSA sample is an entry, plus one sentinel per string");
		rb3_isa_destroy(x);
	}
	test_inverse_lift(outdir);
	test_snap_align();

	rb3_ssa_destroy(sa);
	rb3_fmi_free(&f);
	for (si = 0; si < N_SEQ; ++si) free(seq_fwd[si]);
	fprintf(stderr, "test_extract: %d passed, %d failed\n", n_pass, n_fail);
	return n_fail == 0? 0 : 1;
}
