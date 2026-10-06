/* isa.c -- sampled inverse suffix array (entry points) derived from the .ssa, and
 * coordinate-based substring extraction from an FMD. See isa.h for the method. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "io.h"
#include "ketopt.h"
#include "isa.h"

/***********
 * builder *
 ***********/

typedef struct {
	int64_t m;
	int32_t shift, strands;
	const int32_t *len; // per sequence (BWT string sid has length len[sid>>1])
	int64_t *boff;      // m+1 bucket offsets
	uint32_t *bpos;     // per bucket: smallest sampled pos seen, UINT32_MAX = empty
	int64_t *brow;
	int64_t n_bad;      // samples with sid/pos out of range (a mismatched .len.gz)
} isa_bld_t;

static inline int isa_sid_used(int32_t strands, int64_t sid)
{
	return (strands >> (sid & 1)) & 1;
}

static int isa_bld_init(isa_bld_t *b, int64_t m, int64_t n_seq, const int32_t *len, int shift, int strands)
{
	int64_t sid, nb = 0;
	memset(b, 0, sizeof(*b));
	if (m != n_seq * 2 || shift < 0 || shift > 30 || (strands & 3) == 0) return -1;
	b->m = m, b->shift = shift, b->strands = strands & 3, b->len = len;
	b->boff = RB3_CALLOC(int64_t, m + 1);
	for (sid = 0; sid < m; ++sid) {
		b->boff[sid] = nb;
		if (isa_sid_used(b->strands, sid)) nb += ((int64_t)len[sid>>1] >> shift) + 1;
	}
	b->boff[m] = nb;
	b->bpos = RB3_MALLOC(uint32_t, nb > 0? nb : 1);
	b->brow = RB3_MALLOC(int64_t, nb > 0? nb : 1);
	if (b->bpos == 0 || b->brow == 0) { free(b->boff); free(b->bpos); free(b->brow); return -1; }
	memset(b->bpos, 0xff, nb * sizeof(uint32_t));
	return 0;
}

static inline void isa_bld_add(isa_bld_t *b, int64_t sid, int64_t pos, int64_t row)
{
	int64_t k;
	if (sid < 0 || sid >= b->m || pos < 0 || pos >= b->len[sid>>1]) { ++b->n_bad; return; }
	if (!isa_sid_used(b->strands, sid)) return;
	k = b->boff[sid] + (pos >> b->shift);
	if ((uint32_t)pos < b->bpos[k]) b->bpos[k] = (uint32_t)pos, b->brow[k] = row;
}

static rb3_isa_t *isa_bld_finish(isa_bld_t *b, int32_t ss, int64_t n_ssa)
{
	rb3_isa_t *isa;
	int64_t sid, k, n = 0;
	if (b->n_bad > 0) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: %ld SSA samples fall outside the sequence lengths; the .len.gz does not match the .ssa\n", (long)b->n_bad);
		free(b->boff); free(b->bpos); free(b->brow);
		return 0;
	}
	isa = RB3_CALLOC(rb3_isa_t, 1);
	isa->shift = b->shift, isa->ss = ss, isa->strands = b->strands, isa->m = b->m, isa->n_ssa = n_ssa;
	for (k = 0; k < b->boff[b->m]; ++k) n += (b->bpos[k] != UINT32_MAX);
	for (sid = 0; sid < b->m; ++sid) n += isa_sid_used(b->strands, sid);
	isa->n = n;
	isa->off = RB3_CALLOC(int64_t, b->m + 1);
	isa->pos = RB3_MALLOC(uint32_t, n > 0? n : 1);
	isa->row = RB3_MALLOC(int64_t, n > 0? n : 1);
	for (sid = 0, n = 0; sid < b->m; ++sid) {
		isa->off[sid] = n;
		if (!isa_sid_used(b->strands, sid)) continue;
		for (k = b->boff[sid]; k < b->boff[sid+1]; ++k)
			if (b->bpos[k] != UINT32_MAX)
				isa->pos[n] = b->bpos[k], isa->row[n++] = b->brow[k];
		isa->pos[n] = (uint32_t)b->len[sid>>1], isa->row[n++] = sid; // the sentinel: row sid ends string sid
	}
	isa->off[b->m] = n;
	free(b->boff); free(b->bpos); free(b->brow);
	return isa;
}

rb3_isa_t *rb3_isa_from_ssa(const rb3_ssa_t *sa, int64_t n_seq, const int32_t *len, int shift, int strands)
{
	isa_bld_t b;
	int64_t x, mask;
	if (sa == 0 || isa_bld_init(&b, sa->m, n_seq, len, shift, strands) < 0) return 0;
	mask = (1LL << sa->ms) - 1;
	for (x = 0; x < sa->n_ssa; ++x)
		isa_bld_add(&b, (int64_t)(sa->ssa[x] & mask), (int64_t)(sa->ssa[x] >> sa->ms), sa->m + (x << sa->ss));
	return isa_bld_finish(&b, sa->ss, sa->n_ssa);
}

rb3_isa_t *rb3_isa_from_ssa_file(const char *fn, int64_t n_seq, const int32_t *len, int shift, int strands)
{
	FILE *fp;
	char magic[4];
	uint32_t ss, ms;
	int64_t m, n_ssa, x = 0, mask;
	uint64_t *buf;
	size_t chunk = 1<<20, nr;
	isa_bld_t b;
	rb3_isa_t *isa;

	if ((fp = fopen(fn, "rb")) == 0) return 0;
	if (fread(magic, 1, 4, fp) != 4 || strncmp(magic, "SSA\1", 4) != 0
		|| fread(&ss, 4, 1, fp) != 1 || fread(&ms, 4, 1, fp) != 1
		|| fread(&m, 8, 1, fp) != 1 || fread(&n_ssa, 8, 1, fp) != 1
		|| fseeko(fp, (off_t)m * 8, SEEK_CUR) != 0) {
		fclose(fp);
		return 0;
	}
	if (isa_bld_init(&b, m, n_seq, len, shift, strands) < 0) { fclose(fp); return 0; }
	mask = (1LL << ms) - 1;
	buf = RB3_MALLOC(uint64_t, chunk);
	while (x < n_ssa && (nr = fread(buf, 8, (size_t)(n_ssa - x) < chunk? (size_t)(n_ssa - x) : chunk, fp)) > 0) {
		size_t i;
		for (i = 0; i < nr; ++i, ++x)
			isa_bld_add(&b, (int64_t)(buf[i] & mask), (int64_t)(buf[i] >> ms), m + (x << ss));
	}
	free(buf);
	fclose(fp);
	if (x != n_ssa) { // truncated file
		if (rb3_verbose >= 1) fprintf(stderr, "ERROR: '%s' is truncated (%ld of %ld samples)\n", fn, (long)x, (long)n_ssa);
		free(b.boff); free(b.bpos); free(b.brow);
		return 0;
	}
	isa = isa_bld_finish(&b, ss, n_ssa);
	return isa;
}

/*******
 * I/O *
 *******/

int rb3_isa_dump(const rb3_isa_t *isa, const char *fn)
{
	FILE *fp = fn && strcmp(fn, "-")? fopen(fn, "wb") : stdout;
	int32_t y[3];
	if (fp == 0) return -1;
	y[0] = isa->shift, y[1] = isa->ss, y[2] = isa->strands;
	fwrite("ISA\1", 1, 4, fp);
	fwrite(y, 4, 3, fp);
	fwrite(&isa->m, 8, 1, fp);
	fwrite(&isa->n_ssa, 8, 1, fp);
	fwrite(&isa->n, 8, 1, fp);
	fwrite(isa->off, 8, isa->m + 1, fp);
	fwrite(isa->pos, 4, isa->n, fp);
	fwrite(isa->row, 8, isa->n, fp);
	if (fp != stdout) { if (fclose(fp) != 0) return -1; }
	else fflush(fp);
	return 0;
}

rb3_isa_t *rb3_isa_restore(const char *fn)
{
	FILE *fp;
	char magic[4];
	int32_t y[3];
	rb3_isa_t *isa;
	if ((fp = fopen(fn, "rb")) == 0) return 0;
	if (fread(magic, 1, 4, fp) != 4 || strncmp(magic, "ISA\1", 4) != 0 || fread(y, 4, 3, fp) != 3) { fclose(fp); return 0; }
	isa = RB3_CALLOC(rb3_isa_t, 1);
	isa->shift = y[0], isa->ss = y[1], isa->strands = y[2];
	if (fread(&isa->m, 8, 1, fp) != 1 || fread(&isa->n_ssa, 8, 1, fp) != 1 || fread(&isa->n, 8, 1, fp) != 1
		|| isa->m < 0 || isa->n < 0) {
		free(isa); fclose(fp); return 0;
	}
	isa->off = RB3_MALLOC(int64_t, isa->m + 1);
	isa->pos = RB3_MALLOC(uint32_t, isa->n > 0? isa->n : 1);
	isa->row = RB3_MALLOC(int64_t, isa->n > 0? isa->n : 1);
	if (isa->off == 0 || isa->pos == 0 || isa->row == 0
		|| fread(isa->off, 8, isa->m + 1, fp) != (size_t)(isa->m + 1)
		|| fread(isa->pos, 4, isa->n, fp) != (size_t)isa->n
		|| fread(isa->row, 8, isa->n, fp) != (size_t)isa->n
		|| isa->off[isa->m] != isa->n) {
		rb3_isa_destroy(isa); fclose(fp); return 0;
	}
	fclose(fp);
	return isa;
}

void rb3_isa_destroy(rb3_isa_t *isa)
{
	if (isa == 0) return;
	free(isa->off); free(isa->pos); free(isa->row); free(isa);
}

int rb3_isa_check_fmi(const rb3_isa_t *isa, const rb3_fmi_t *f)
{
	int64_t n_ssa;
	if (isa == 0 || f == 0 || isa->m != f->acc[1] || isa->ss < 0 || isa->ss > 62) return -1;
	n_ssa = (f->acc[RB3_ASIZE] - f->acc[1] + (1LL<<isa->ss) - 1LL) >> isa->ss;
	return n_ssa == isa->n_ssa? 0 : -1;
}

int64_t rb3_isa_seq_len(const rb3_isa_t *isa, int64_t sid)
{
	if (isa == 0 || sid < 0 || sid >= isa->m || isa->off[sid+1] == isa->off[sid]) return -1;
	return isa->pos[isa->off[sid+1] - 1]; // the last entry is always the sentinel at pos = len
}

/**************
 * extraction *
 **************/

int64_t rb3_isa_extract(const rb3_fmi_t *f, const rb3_isa_t *isa, int64_t sid, int64_t st, int64_t en, uint8_t *buf)
{
	int64_t lo, hi, len, p, k, ok[RB3_ASIZE];
	len = rb3_isa_seq_len(isa, sid);
	if (len < 0) return RB3_ISA_ERR_SID;
	if (st < 0 || en > len || st >= en) return RB3_ISA_ERR_RANGE;
	lo = isa->off[sid], hi = isa->off[sid+1];
	while (lo < hi) { // first entry with pos >= en; exists since the last one has pos = len >= en
		int64_t mid = lo + ((hi - lo) >> 1);
		if ((int64_t)isa->pos[mid] < en) lo = mid + 1;
		else hi = mid;
	}
	k = isa->row[lo];
	for (p = isa->pos[lo] - 1; p >= st; --p) { // LF(k) is the row of the suffix at p; BWT[k] = base p
		int c = rb3_fmi_rank1a(f, k, ok);
		if (c <= 0) return RB3_ISA_ERR_WALK;
		if (p < en) buf[p - st] = (uint8_t)c;
		k = f->acc[c] + ok[c];
	}
	return isa->pos[lo] - st;
}

/***********************
 * the `isa` subcommand *
 ***********************/

int main_isa(int argc, char *argv[])
{
	int c, shift = 12, strands = RB3_ISA_FWD;
	ketopt_t o = KETOPT_INIT;
	char *fn_out = 0, *buf;
	rb3_sid_t *sid;
	rb3_isa_t *isa;
	double t0;

	while ((c = ketopt(&o, argc, argv, 1, "b:o:2", 0)) >= 0) {
		if (c == 'b') shift = atoi(o.arg);
		else if (c == 'o') fn_out = o.arg;
		else if (c == '2') strands = RB3_ISA_FWD | RB3_ISA_REV;
	}
	if (argc == o.ind) {
		fprintf(stderr, "Usage: ropebwt3 isa [options] <idx.fmd>\n");
		fprintf(stderr, "Derives the sampled inverse suffix array (extraction entry points) used by\n");
		fprintf(stderr, "`ropebwt3 extract` from <idx.fmd>.ssa and <idx.fmd>.len.gz (the BWT itself is not read).\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b INT     keep at most one entry per 2**INT bases of each sequence [%d]\n", shift);
		fprintf(stderr, "  -2         also sample reverse-complement strings (extract never needs them) [forward only]\n");
		fprintf(stderr, "  -o FILE    output file [<idx.fmd>.isa]\n");
		return 1;
	}
	if (shift < 0 || shift > 30) {
		fprintf(stderr, "ERROR: -b must be in [0,30]\n");
		return 1;
	}
	t0 = rb3_realtime();
	buf = RB3_CALLOC(char, strlen(argv[o.ind]) + 8);
	sid = rb3_sid_read(strcat(strcpy(buf, argv[o.ind]), ".len.gz"));
	if (sid == 0 || sid->n_seq == 0) {
		fprintf(stderr, "ERROR: failed to read sequence names/lengths from '%s'\n", buf);
		free(buf); rb3_sid_destroy(sid);
		return 1;
	}
	strcat(strcpy(buf, argv[o.ind]), ".ssa");
	isa = rb3_isa_from_ssa_file(buf, sid->n_seq, sid->len, shift, strands);
	if (isa == 0) {
		fprintf(stderr, "ERROR: failed to derive entry points from '%s' (missing, truncated, or #sequences != 2 x %ld in .len.gz)\n", buf, (long)sid->n_seq);
		free(buf); rb3_sid_destroy(sid);
		return 1;
	}
	if (fn_out == 0) fn_out = strcat(strcpy(buf, argv[o.ind]), ".isa");
	if (rb3_isa_dump(isa, fn_out) < 0) {
		fprintf(stderr, "ERROR: failed to write '%s'\n", fn_out);
		rb3_isa_destroy(isa); free(buf); rb3_sid_destroy(sid);
		return 1;
	}
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] %ld entry points (2^%d spacing, from 2^%d SSA) for %ld sequences in %.3f sec -> %s\n",
				__func__, (long)isa->n, isa->shift, isa->ss, (long)sid->n_seq, rb3_realtime() - t0, fn_out);
	rb3_isa_destroy(isa);
	rb3_sid_destroy(sid);
	free(buf);
	return 0;
}
