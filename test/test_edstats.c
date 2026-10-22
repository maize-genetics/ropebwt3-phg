// Unit tests for edstats.c (edit-distance summaries, histograms, npy header), no index needed.
#include <stdio.h>
#include <string.h>
#include "edstats.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(cond, msg) do { if (cond) ++n_pass; else { ++n_fail; fprintf(stderr, "FAIL: %s\n", msg); } } while (0)
#define NA RB3_ED_NA

int main(void)
{
	int32_t nb, w, best;
	// summary: best, assemblies at best, closest distinct wrong copy (> best + 1)
	{
		uint8_t v[6] = { 2, 0, 1, 0, NA, 4 };
		best = rb3_ed_summary(v, 6, NA, &nb, &w);
		CHECK(best == 0 && nb == 2, "best 0 held by two assemblies");
		CHECK(w == 2, "wrong copy: the smallest value > best + 1, ignoring NA");
	}
	{
		uint8_t v[6] = { 3, 0, 1, 0, NA, 2 };
		rb3_ed_summary(v, 6, NA, &nb, &w);
		CHECK(w == 2, "wrong = 2 (the smallest value > best + 1)");
	}
	{
		uint8_t v[3] = { 1, 1, 2 };
		best = rb3_ed_summary(v, 3, NA, &nb, &w);
		CHECK(best == 1 && nb == 2 && w == NA, "no assembly > best+1: wrong is NA");
	}
	{
		uint8_t v[3] = { NA, NA, NA };
		best = rb3_ed_summary(v, 3, NA, &nb, &w);
		CHECK(best == NA && nb == 0 && w == NA, "all NA: best NA, nobody at best");
		best = rb3_ed_summary(v, 3, 0, &nb, &w);
		CHECK(best == 0 && nb == 0, "a known best whose assembly was not located: best 0, n_best 0");
	}
	{
		uint8_t v[3] = { 2, 3, NA };
		best = rb3_ed_summary(v, 3, 1, &nb, &w);
		CHECK(best == 1 && nb == 0 && w == 3, "gbest below every value: wrong measured from gbest");
	}
	// pairs: per-assembly sum, NA if either mate is NA
	{
		uint8_t a[4] = { 0, 1, NA, 3 }, b[4] = { 2, NA, 0, 3 }, o[4];
		rb3_ed_pair(a, b, 4, o);
		CHECK(o[0] == 2 && o[1] == NA && o[2] == NA && o[3] == 6, "pair sums");
	}
	// histogram: one read, then the writer
	{
		rb3_edhist_t *h = rb3_edhist_init(3, 4);
		uint8_t v1[3] = { 0, 2, NA }, v2[3] = { 1, 1, 9 }, p[3];
		FILE *fp = tmpfile();
		char line[512];
		int n_best_rows = 0, n_asm_rows = 0, n_pair = 0, n_hdr = 0;
		long cnt_best0 = -1, cnt_asm_na = -1, cnt_wrong2 = -1;
		rb3_edhist_add(h, RB3_ED_UNIT_READ, v1, NA, 0);
		rb3_edhist_add(h, RB3_ED_UNIT_READ, v2, NA, RB3_ED_F_TRUNC | 3 << 8);
		rb3_ed_pair(v1, v2, 3, p);
		rb3_edhist_add(h, RB3_ED_UNIT_PAIR, p, NA, 0);
		CHECK(h->n[0] == 2 && h->n[1] == 1 && h->n_trunc[0] == 1, "counts and truncated reads");
		CHECK(h->best[0][0] == 1 && h->best[0][1] == 1, "read best histogram");
		CHECK(h->asm_[0][5] == 2, "per-assembly NA and > cap (9) share the NA bin");
		CHECK(h->wrong[0][2] == 1 && h->wrong[0][5] == 1, "wrong copy: 2 for read 1, none for read 2");
		CHECK(h->nbest[0][1] == 1 && h->nbest[0][2] == 1, "assemblies at best: 1 and 2");
		CHECK(h->cap[1] == 8 && h->best[1][1] == 1, "pair cap is twice the read cap; pair best = 0+1");
		rb3_edhist_write(fp, h, 0, 1e-3);
		rewind(fp);
		while (fgets(line, sizeof(line), fp)) {
			char u[16], m[16], k[16];
			long c;
			if (line[0] == '#') { ++n_hdr; continue; }
			if (sscanf(line, "%15s %15s %15s %ld", u, m, k, &c) != 4 || strcmp(u, "unit") == 0) continue;
			if (strcmp(u, "pair") == 0) { ++n_pair; continue; }
			if (strcmp(m, "best") == 0) { ++n_best_rows; if (strcmp(k, "0") == 0) cnt_best0 = c; }
			if (strcmp(m, "asm") == 0) { ++n_asm_rows; if (strcmp(k, "NA") == 0) cnt_asm_na = c; }
			if (strcmp(m, "wrong") == 0 && strcmp(k, "2") == 0) cnt_wrong2 = c;
		}
		fclose(fp);
		CHECK(n_hdr >= 5, "writer: header lines incl. one #unit line per unit");
		CHECK(n_best_rows == 6 && n_asm_rows == 6, "writer: cap+2 rows (0..cap, NA) per measure");
		CHECK(cnt_best0 == 1 && cnt_asm_na == 2 && cnt_wrong2 == 1, "writer: counts round-trip");
		CHECK(n_pair == 3 * 10 + 4, "writer: pair block = 3 measures x (0..8, NA) + nbest 0..3");
		rb3_edhist_destroy(h);
	}
	// npy header: fixed 128 bytes, shape rewritten in place
	{
		FILE *fp = tmpfile();
		unsigned char buf[RB3_ED_NPY_HDR + 1];
		rb3_ed_npy_header(fp, 0, 27);
		rb3_ed_npy_header(fp, 1520864, 27);
		CHECK(ftell(fp) == 2 * RB3_ED_NPY_HDR, "header is exactly RB3_ED_NPY_HDR bytes");
		rewind(fp);
		CHECK(fread(buf, 1, RB3_ED_NPY_HDR, fp) == RB3_ED_NPY_HDR, "read header");
		buf[RB3_ED_NPY_HDR] = 0;
		CHECK(memcmp(buf, "\x93NUMPY\x01\x00", 8) == 0 && buf[8] + 256 * buf[9] == RB3_ED_NPY_HDR - 10, "magic, version, header length");
		CHECK(strstr((char*)buf + 10, "'descr': '|u1'") && strstr((char*)buf + 10, "'shape': (0, 27)"), "dtype uint8, shape");
		CHECK(buf[RB3_ED_NPY_HDR - 1] == '\n', "header ends with a newline");
		fclose(fp);
	}
	printf("test_edstats: %d passed, %d failed\n", n_pass, n_fail);
	return n_fail? 1 : 0;
}
