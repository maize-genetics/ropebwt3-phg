/* Unit tests for asmqc.c (asmqc / asmqc-fa). The module's helpers are static, so the source is
 * included directly; nothing here goes through the CLI or reads files.
 *
 * Covered:
 *  - qa_gapcost: AnchorWave's 2-piece affine gap cost min(4+2l, 80+l), including the crossover
 *  - qa_ksw: with the positive-match conversion, ksw2's alignment has exactly the optimal AnchorWave
 *    penalty -- checked against an independent unbanded 2-piece affine DP (match 0, mismatch 4) on
 *    300 random pairs; and replacement is preferred for unrelated sequence (check_rep)
 *  - qa_snp_classes: local SNP-density classes against a brute-force count
 *  - qa_align: SNP positions, a 3 bp deletion, a 5 kb deletion, a 4 kb insertion, a replaced block
 *  - qa_extend: extension crosses SNPs and stops where the sequences diverge, both directions
 *  - qf_build_table / qf_scan: rolling-hash anchor table -- uniqueness counted on both strands, windows
 *    with N skipped, assembly hits found on the reverse strand at the right forward position
 *  - qc_same_chrom / qc_is_chr, qc_union, qc_acc_merge, qf_taxon_from_path
 */
#include "../asmqc.c"

static int n_pass = 0, n_fail = 0;
#define CHECK(cond, msg) do { if (cond) ++n_pass; else { ++n_fail; fprintf(stderr, "FAIL: %s\n", msg); } } while (0)

static uint64_t rs = 88172645463325252ULL;
static uint32_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (uint32_t)rs; }
static void rand_seq(uint8_t *s, int n) { for (int i = 0; i < n; ++i) s[i] = 1 + rnd() % 4; } // nt6: 1..4 = ACGT

/* reference: minimum AnchorWave penalty (match 0, mismatch 4, gaps min over pieces of o+e*l), full DP */
static int64_t ref_penalty(const uint8_t *R, int lr, const uint8_t *A, int la)
{
	const int64_t INF = 1LL << 50, o1 = 4, e1 = 2, o2 = 80, e2 = 1;
	int64_t sz = (int64_t)(lr + 1) * (la + 1) * 8;
	int64_t *H = (int64_t*)malloc(sz), *E1 = (int64_t*)malloc(sz), *E2 = (int64_t*)malloc(sz);
	int64_t *F1 = (int64_t*)malloc(sz), *F2 = (int64_t*)malloc(sz), res;
#define I(i, j) ((int64_t)(i) * (la + 1) + (j))
	for (int i = 0; i <= lr; ++i) for (int j = 0; j <= la; ++j) {
		int64_t h = INF, e1v = INF, e2v = INF, f1v = INF, f2v = INF;
		if (i == 0 && j == 0) { H[I(0,0)] = 0; E1[I(0,0)] = E2[I(0,0)] = F1[I(0,0)] = F2[I(0,0)] = INF; continue; }
		if (j > 0) { // consume A (insertion)
			e1v = H[I(i,j-1)] + o1 + e1; if (E1[I(i,j-1)] + e1 < e1v) e1v = E1[I(i,j-1)] + e1;
			e2v = H[I(i,j-1)] + o2 + e2; if (E2[I(i,j-1)] + e2 < e2v) e2v = E2[I(i,j-1)] + e2;
		}
		if (i > 0) { // consume R (deletion)
			f1v = H[I(i-1,j)] + o1 + e1; if (F1[I(i-1,j)] + e1 < f1v) f1v = F1[I(i-1,j)] + e1;
			f2v = H[I(i-1,j)] + o2 + e2; if (F2[I(i-1,j)] + e2 < f2v) f2v = F2[I(i-1,j)] + e2;
		}
		if (i > 0 && j > 0) h = H[I(i-1,j-1)] + (R[i-1] == A[j-1]? 0 : 4);
		if (e1v < h) h = e1v;
		if (e2v < h) h = e2v;
		if (f1v < h) h = f1v;
		if (f2v < h) h = f2v;
		H[I(i,j)] = h, E1[I(i,j)] = e1v, E2[I(i,j)] = e2v, F1[I(i,j)] = f1v, F2[I(i,j)] = f2v;
	}
	res = H[I(lr,la)];
#undef I
	free(H); free(E1); free(E2); free(F1); free(F2);
	return res;
}

/* banded version for long pairs: |j - i| <= w (the band qa_ksw uses), rolling rows */
static int64_t ref_penalty_banded(const uint8_t *R, int lr, const uint8_t *A, int la, int w)
{
	const int64_t INF = 1LL << 50, o1 = 4, e1 = 2, o2 = 80, e2 = 1;
	int W = 2 * w + 1;
	int64_t *ph = (int64_t*)malloc(W * 8 * 5), *ch = (int64_t*)malloc(W * 8 * 5), *t, res;
	// row layout: [H | E1 | E2 | F1 | F2], index x = j - i + w
#define AT(row, k, x) (row)[(k) * W + (x)]
	for (int x = 0; x < W * 5; ++x) ph[x] = INF;
	for (int i = 0; i <= lr; ++i) {
		for (int x = 0; x < W * 5; ++x) ch[x] = INF;
		for (int x = 0; x < W; ++x) {
			int j = i + x - w;
			int64_t h = INF, a1 = INF, a2 = INF, b1 = INF, b2 = INF;
			if (j < 0 || j > la) continue;
			if (i == 0 && j == 0) { AT(ch, 0, x) = 0; continue; }
			if (j > 0 && x > 0) { // same row, j-1
				a1 = AT(ch, 0, x-1) + o1 + e1; if (AT(ch, 1, x-1) + e1 < a1) a1 = AT(ch, 1, x-1) + e1;
				a2 = AT(ch, 0, x-1) + o2 + e2; if (AT(ch, 2, x-1) + e2 < a2) a2 = AT(ch, 2, x-1) + e2;
			}
			if (i > 0 && x + 1 < W) { // previous row, same j -> x+1
				b1 = AT(ph, 0, x+1) + o1 + e1; if (AT(ph, 3, x+1) + e1 < b1) b1 = AT(ph, 3, x+1) + e1;
				b2 = AT(ph, 0, x+1) + o2 + e2; if (AT(ph, 4, x+1) + e2 < b2) b2 = AT(ph, 4, x+1) + e2;
			}
			if (i > 0 && j > 0) h = AT(ph, 0, x) + (R[i-1] == A[j-1]? 0 : 4);
			if (a1 < h) h = a1;
			if (a2 < h) h = a2;
			if (b1 < h) h = b1;
			if (b2 < h) h = b2;
			AT(ch, 0, x) = h, AT(ch, 1, x) = a1, AT(ch, 2, x) = a2, AT(ch, 3, x) = b1, AT(ch, 4, x) = b2;
		}
		t = ph, ph = ch, ch = t;
	}
	res = AT(ph, 0, la - lr + w);
#undef AT
	free(ph); free(ch);
	return res;
}

static int mutate(const uint8_t *s, int n, uint8_t *o, double snp, double indel, int big)
{
	int m = 0;
	for (int i = 0; i < n; ++i) {
		double r = (rnd() % 100000) / 100000.0;
		if (big && i == n / 2) { if (rnd() % 2) { for (int z = 0; z < 100; ++z) o[m++] = 1 + rnd() % 4; } else { i += 99; continue; } }
		if (r < indel / 2) { int L = 1 + rnd() % 5; for (int z = 0; z < L; ++z) o[m++] = 1 + rnd() % 4; o[m++] = s[i]; }
		else if (r < indel) { i += rnd() % 5; }
		else if (r < indel + snp) { o[m++] = 1 + (s[i] - 1 + 1 + rnd() % 3) % 4; }
		else o[m++] = s[i];
	}
	return m;
}

static int has_var(const qa_buf_t *b, int type, int64_t pos, int64_t tol, int64_t len, int64_t len_tol)
{
	for (int64_t i = 0; i < b->n_v; ++i) {
		const qa_var_t *v = &b->v[i];
		int64_t L = type == QA_VDEL || type == QA_VREP? v->rlen : v->alen;
		if (v->type == type && llabs(v->pos - pos) <= tol && (type == QA_VSNP || llabs(L - len) <= len_tol)) return 1;
	}
	return 0;
}

static void test_gapcost(void)
{
	CHECK(qa_gapcost(0) == 0, "gapcost(0) = 0");
	CHECK(qa_gapcost(1) == 6, "gapcost(1) = 4+2");
	CHECK(qa_gapcost(76) == 156 && 4 + 2 * 76 == 156 && 80 + 76 == 156, "gapcost crossover at 76 (both pieces 156)");
	CHECK(qa_gapcost(77) == 157, "gapcost(77) uses the long piece 80+77");
	CHECK(qa_gapcost(1000) == 1080, "gapcost(1000) = 80+1000");
}

static void test_ksw_optimal(void)
{
	uint8_t R[400], A[700];
	int bad = 0;
	qa_buf_t b; memset(&b, 0, sizeof(b));
	for (int t = 0; t < 300; ++t) {
		int lr = 40 + rnd() % 260, la;
		qa_stat_t st; memset(&st, 0, sizeof(st));
		rand_seq(R, lr);
		la = mutate(R, lr, A, 0.02 + (t % 5) * 0.05, (t % 3) * 0.02, t % 7 == 0);
		if (la == 0) { A[0] = 1; la = 1; }
		qa_ksw(R, lr, A, la, &b, &st, 0, 0);
		if (st.pen != ref_penalty(R, lr, A, la)) {
			if (!bad) fprintf(stderr, "  case %d: lr=%d la=%d asmqc penalty %ld, optimum %ld\n", t, lr, la, (long)st.pen, (long)ref_penalty(R, lr, A, la));
			++bad;
		}
		free(st.hist);
	}
	CHECK(bad == 0, "qa_ksw (converted scoring) reaches the optimal AnchorWave penalty on 300 random pairs");
	{ // long diverged pairs (20 kb, 25% SNPs, 2% indels): ksw2 with a zero match score mis-scores these
		uint8_t *X = (uint8_t*)malloc(20000), *Y = (uint8_t*)malloc(22000);
		int bad2 = 0;
		for (int t = 0; t < 4; ++t) {
			int lx = 19000 + (int)(rnd() % 1000), ly;
			qa_stat_t st; memset(&st, 0, sizeof(st));
			rand_seq(X, lx);
			ly = mutate(X, lx, Y, 0.25, 0.02, 0);
			qa_ksw(X, lx, Y, ly, &b, &st, 0, 0);
			int64_t opt = ref_penalty_banded(X, lx, Y, ly, abs(lx - ly) + QA_KSW_W);
			if (st.pen != opt) { if (!bad2) fprintf(stderr, "  long case %d: lr=%d la=%d asmqc %ld, banded optimum %ld\n", t, lx, ly, (long)st.pen, (long)opt); ++bad2; }
			free(st.hist);
		}
		CHECK(bad2 == 0, "qa_ksw reaches the banded optimum on 20 kb pairs at 25% divergence");
		free(X); free(Y);
	}
	{ // unrelated 2 kb sequences: deleting + inserting is cheaper than aligning
		uint8_t *X = (uint8_t*)malloc(2000), *Y = (uint8_t*)malloc(2000);
		qa_stat_t st; memset(&st, 0, sizeof(st));
		rand_seq(X, 2000), rand_seq(Y, 2000);
		CHECK(qa_ksw(X, 2000, Y, 2000, &b, &st, 0, 1) == 0 && st.snp == 0, "unrelated 2 kb: check_rep prefers replacement, adds nothing");
		mutate(X, 2000, Y, 0.05, 0.0, 0);
		CHECK(qa_ksw(X, 2000, Y, 2000, &b, &st, 0, 1) == 1 && st.snp > 50, "5%-diverged 2 kb: aligned under check_rep");
		free(st.hist); free(X); free(Y);
	}
	qa_buf_free(&b);
}

static void test_snp_classes(void)
{
	int64_t p[200]; int n = 0;
	for (int i = 0; i < 100; i += 2) p[n++] = 1000 + i;      // dense block: 50 SNPs in 100 bp
	for (int i = 0; i < 20; ++i) p[n++] = 5000 + 37 * i;     // sparse
	p[n++] = 20000;                                          // isolated
	qa_stat_t st; memset(&st, 0, sizeof(st));
	qa_snp_classes(p, n, &st);
	int64_t d10 = 0, d40 = 0;
	for (int i = 0; i < n; ++i) { int c = 0; for (int j = 0; j < n; ++j) if (j != i && llabs(p[j] - p[i]) <= 50) ++c; d10 += c >= 10; d40 += c >= 40; }
	CHECK(st.snp_d10 == d10 && st.snp_d40 == d40, "qa_snp_classes matches a brute-force neighbour count");
	CHECK(d40 > 0 && d10 > d40, "fixture has both density classes");
}

static void test_align(void)
{
	qa_buf_t b; memset(&b, 0, sizeof(b)); b.dump = 1;
	uint8_t *R = (uint8_t*)malloc(20000), *A = (uint8_t*)malloc(20000);
	qa_stat_t st;
	// SNPs + a 3 bp deletion inside a 2 kb gap
	rand_seq(R, 2000); memcpy(A, R, 2000);
	A[100] = 1 + (R[100] % 4); A[900] = 1 + (R[900] % 4);
	memmove(A + 1500, A + 1503, 500 - 3);
	memset(&st, 0, sizeof(st)); b.n_v = 0;
	qa_align(R, 2000, A, 1997, 0, &b, &st, 1000000);
	CHECK(st.snp == 2 && has_var(&b, QA_VSNP, 1000100, 0, 0, 0) && has_var(&b, QA_VSNP, 1000900, 0, 0, 0), "qa_align: both SNPs at their exact reference positions");
	CHECK(st.del_n == 1 && st.del_bp == 3 && has_var(&b, QA_VDEL, 1001500, 3, 3, 0), "qa_align: one 3 bp deletion at ~1,001,500");
	free(st.hist);
	// a 5 kb deletion between 1 kb flanks
	rand_seq(R, 7000); memcpy(A, R, 1000); memcpy(A + 1000, R + 6000, 1000);
	memset(&st, 0, sizeof(st)); b.n_v = 0;
	qa_align(R, 7000, A, 2000, 0, &b, &st, 0);
	CHECK(st.snp == 0 && st.del_bp == 5000 && has_var(&b, QA_VDEL, 1000, 20, 5000, 0), "qa_align: 5 kb deletion, no SNPs");
	free(st.hist);
	// a 4 kb insertion
	rand_seq(R, 2000); memcpy(A, R, 1000); rand_seq(A + 1000, 4000); memcpy(A + 5000, R + 1000, 1000);
	memset(&st, 0, sizeof(st)); b.n_v = 0;
	qa_align(R, 2000, A, 6000, 0, &b, &st, 0);
	CHECK(st.snp == 0 && st.ins_bp == 4000 && has_var(&b, QA_VINS, 1000, 20, 4000, 0), "qa_align: 4 kb insertion, no SNPs");
	free(st.hist);
	// a replaced 3 kb block between 1 kb flanks: reference 3 kb deleted, assembly 3 kb inserted
	rand_seq(R, 5000); memcpy(A, R, 1000); rand_seq(A + 1000, 3000); memcpy(A + 4000, R + 4000, 1000);
	memset(&st, 0, sizeof(st)); b.n_v = 0;
	qa_align(R, 5000, A, 5000, 0, &b, &st, 0);
	if (getenv("TEST_ASMQC_VERBOSE")) {
		fprintf(stderr, "  replaced-block case: del_bp=%ld ins_bp=%ld snp=%ld del_n=%ld ins_n=%ld pen=%ld; records:", (long)st.del_bp, (long)st.ins_bp, (long)st.snp, (long)st.del_n, (long)st.ins_n, (long)st.pen);
		for (int64_t z = 0; z < b.n_v && z < 12; ++z) fprintf(stderr, " %d@%ld(%ld/%ld)", b.v[z].type, (long)b.v[z].pos, (long)b.v[z].rlen, (long)b.v[z].alen);
		fprintf(stderr, "\n");
	}
	CHECK(st.del_bp >= 2950 && st.del_bp <= 3050 && st.ins_bp >= 2950 && st.ins_bp <= 3050 && st.snp < 30,
		  "qa_align: unrelated 3 kb block is deleted + inserted, not force-aligned");
	free(st.hist);
	free(R); free(A); qa_buf_free(&b);
}

static void test_extend(void)
{
	uint8_t *R = (uint8_t*)malloc(3000), *A = (uint8_t*)malloc(3000);
	int64_t er, ea;
	rand_seq(R, 3000); memcpy(A, R, 1000); rand_seq(A + 1000, 2000);
	for (int i = 50; i < 1000; i += 50) A[i] = 1 + (R[i] % 4); // SNP every 50 bp inside the shared part
	qa_extend(R, 3000, A, 3000, 0, &er, &ea);
	CHECK(er >= 950 && er <= 1010 && ea >= 950 && ea <= 1010, "qa_extend: crosses SNPs, stops at ~1000 where the sequences diverge");
	rand_seq(R, 3000); rand_seq(A, 2000); memcpy(A + 2000, R + 2000, 1000);
	qa_extend(R, 3000, A, 3000, 1, &er, &ea);
	CHECK(er >= 950 && er <= 1010 && ea >= 950 && ea <= 1010, "qa_extend (reverse): extends ~1000 from the end");
	free(R); free(A);
}

static void test_hash_table(void)
{
	const int k = 31, stride = 200, L = 4000;
	qf_seqs_t ref, as;
	qf_tab_t T;
	memset(&ref, 0, sizeof(ref)); memset(&as, 0, sizeof(as)); memset(&T, 0, sizeof(T));
	ref.n = 2; ref.m = 2;
	ref.name = (char**)calloc(2, sizeof(char*)); ref.len = (int64_t*)calloc(2, sizeof(int64_t)); ref.seq = (uint8_t**)calloc(2, sizeof(uint8_t*));
	ref.name[0] = strdup("REF_chr1"); ref.name[1] = strdup("REF_chr2");
	ref.len[0] = L; ref.len[1] = 1000;
	ref.seq[0] = (uint8_t*)malloc(L); ref.seq[1] = (uint8_t*)malloc(1000);
	rand_seq(ref.seq[0], L); rand_seq(ref.seq[1], 1000);
	memcpy(ref.seq[1] + 100, ref.seq[0] + 400, k);                                  // forward copy of the k-mer at 400
	for (int i = 0; i < k; ++i) ref.seq[1][600 + i] = 5 - ref.seq[0][800 + k - 1 - i]; // reverse-complement copy of 800
	ref.seq[0][1210] = 5;                                                           // N inside the window at 1200
	qf_build_table(&T, &ref, k, stride);
	qf_scan(&T, &ref, 0, 1, 2);
	int64_t got[4] = {-1, -1, -1, -1}, has1200 = 0;
	for (int64_t i = 0; i < T.n; ++i) {
		if (T.a[i].rsid != 0) continue;
		if (T.a[i].rpos == 0) got[0] = T.a[i].rocc;
		if (T.a[i].rpos == 400) got[1] = T.a[i].rocc;
		if (T.a[i].rpos == 800) got[2] = T.a[i].rocc;
		if (T.a[i].rpos == 1200) has1200 = 1;
	}
	CHECK(got[0] == 1, "anchor table: a unique k-mer occurs once");
	CHECK(got[1] == 2, "anchor table: a forward copy on another sequence makes it non-unique");
	CHECK(got[2] == 2, "anchor table: a reverse-complement copy counts too (both strands)");
	CHECK(!has1200, "anchor table: windows containing N are skipped");
	// assembly = reverse complement of REF_chr1: anchors found on the reverse strand at the right forward start
	as.n = as.m = 1; as.name = (char**)calloc(1, sizeof(char*)); as.len = (int64_t*)calloc(1, sizeof(int64_t)); as.seq = (uint8_t**)calloc(1, sizeof(uint8_t*));
	as.name[0] = strdup("CAND_chr1"); as.len[0] = L; as.seq[0] = (uint8_t*)malloc(L);
	for (int i = 0; i < L; ++i) { uint8_t c = ref.seq[0][L - 1 - i]; as.seq[0][i] = c >= 1 && c <= 4? 5 - c : c; }
	qf_scan(&T, &as, 2, 0, 2);
	int ok = 1, n_chk = 0;
	for (int64_t i = 0; i < T.n; ++i)
		if (T.a[i].rsid == 0 && T.a[i].rocc == 1) {
			++n_chk;
			if (!(T.a[i].cocc == 1 && T.a[i].csid == 2 && T.a[i].cpos == L - T.a[i].rpos - k)) ok = 0;
		}
	CHECK(ok && n_chk >= 15, "assembly scan: every unique anchor found once on the reverse strand, cpos = L - rpos - k");
	free(T.a); free(T.slot); qf_seqs_free(&ref); qf_seqs_free(&as);
}

static void test_names_and_intervals(void)
{
	CHECK(qc_same_chrom("chr01", "chr1") && qc_same_chrom("Chr1", "chr1") && qc_same_chrom("chr10", "CHR10"), "qc_same_chrom: case and leading zeros ignored");
	CHECK(!qc_same_chrom("chr1", "chr10") && !qc_same_chrom("chr1", "scaf_1"), "qc_same_chrom: chr1 != chr10, chr != scaffold");
	CHECK(qc_is_chr("chr5") && qc_is_chr("Chr5") && !qc_is_chr("scaf_31") && !qc_is_chr("Pt"), "qc_is_chr");
	CHECK(!strcmp(qc_suffix("Oh43_chr10"), "chr10") && !strcmp(qc_suffix("GEMS-0219_scaf_3"), "scaf_3"), "qc_suffix: label after the first '_'");
	{
		char *a = qf_taxon_from_path("/x/y/Oh7B_renamed.fa"), *b2 = qf_taxon_from_path("B73.fa.gz"), *c = qf_taxon_from_path("GEMS_0219.fasta");
		CHECK(!strcmp(a, "Oh7B") && !strcmp(b2, "B73") && !strcmp(c, "GEMS-0219"), "qf_taxon_from_path: extensions, _renamed and '_' handled");
		free(a); free(b2); free(c);
	}
	qc_acc_t x, y; qc_ivl_t buf[8];
	memset(&x, 0, sizeof(x)); memset(&y, 0, sizeof(y));
	qc_add_ivl(&x, 0, 10, QC_SYN_FWD); qc_add_ivl(&x, 5, 20, QC_SYN_FWD); qc_add_ivl(&x, 30, 40, QC_SYN_INV);
	CHECK(qc_union(&x, 1u << QC_SYN_FWD, buf) == 20, "qc_union: overlapping spans counted once");
	CHECK(qc_union(&x, 1u << QC_SYN_FWD | 1u << QC_SYN_INV, buf) == 30, "qc_union: class mask");
	y.n_pts = 7; y.ins_max = 50; y.del_max = 3; y.b.snp = 4; y.b.pen = 16; y.b.del_max = 9;
	qa_event(&y.b, 9, 1, 0, 0);
	qc_add_ivl(&y, 100, 200, QC_OTHER);
	x.ins_max = 10; x.b.snp = 1;
	qc_acc_merge(&x, &y);
	CHECK(x.n_pts == 7 && x.ins_max == 50 && x.b.snp == 5 && x.b.del_bp == 9 && x.b.hist && x.b.hist[QA_HIST + 9] == 1 && x.n_ivl == 4,
		  "qc_acc_merge: sums, maxima, histogram and intervals");
	free(x.ivl); free(y.ivl); free(x.b.hist); free(y.b.hist);
}

int main(void)
{
	test_gapcost();
	test_ksw_optimal();
	test_snp_classes();
	test_align();
	test_extend();
	test_hash_table();
	test_names_and_intervals();
	printf("test_asmqc: %d passed, %d failed\n", n_pass, n_fail);
	return n_fail? 1 : 0;
}
