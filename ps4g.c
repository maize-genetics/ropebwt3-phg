#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "rb3priv.h"
#include "ps4g.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

/*******************
 * gamete table    *
 *******************/

static int32_t gtab_sample_len(const char *nm)
{
	const char *us = strchr(nm, '_');
	return us? (int32_t)(us - nm) : (int32_t)strlen(nm);
}

// Bare contig name for a sequence, e.g. "B73_chr1" -> "chr1": strip the sample
// prefix *and* its following '_' separator, using the same first-'_' split as
// gtab_sample_len/rb3_gtab_build. Deliberately independent of the exact
// --ref-prefix string a user passes (e.g. "B73" vs "B73_") -- a raw
// name+strlen(ref_prefix) offset left a stray leading '_' ("_chr1") whenever
// --ref-prefix didn't itself include the trailing underscore, which is the
// convention actually used in this project's own scripts (--ref-prefix=B73).
static const char *contig_suffix(const char *name)
{
	const char *us = strchr(name, '_');
	return us? us + 1 : name;
}

static int gtab_strp_cmp(const void *a, const void *b)
{
	return strcmp(*(char *const*)a, *(char *const*)b);
}

int32_t rb3_gtab_find(const rb3_gtab_t *g, const char *name)
{
	int32_t lo = 0, hi = g->n_gamete - 1;
	while (lo <= hi) {
		int32_t mid = (lo + hi) / 2;
		int c = strcmp(g->name[mid], name);
		if (c == 0) return mid;
		else if (c < 0) lo = mid + 1;
		else hi = mid - 1;
	}
	return -1;
}

rb3_gtab_t *rb3_gtab_build(const rb3_sid_t *sid)
{
	int64_t k;
	char **pre, **tmp;
	rb3_gtab_t *g;

	if (sid == 0 || sid->n_seq == 0) return 0;
	g = RB3_CALLOC(rb3_gtab_t, 1);
	pre = RB3_MALLOC(char*, sid->n_seq);
	for (k = 0; k < sid->n_seq; ++k) {
		int32_t plen = gtab_sample_len(sid->name[k]);
		pre[k] = RB3_MALLOC(char, plen + 1);
		memcpy(pre[k], sid->name[k], plen);
		pre[k][plen] = 0;
	}
	tmp = RB3_MALLOC(char*, sid->n_seq);
	memcpy(tmp, pre, sid->n_seq * sizeof(char*));
	qsort(tmp, sid->n_seq, sizeof(char*), gtab_strp_cmp);
	g->name = RB3_MALLOC(char*, sid->n_seq);
	g->n_gamete = 0;
	for (k = 0; k < sid->n_seq; ++k)
		if (k == 0 || strcmp(tmp[k], tmp[k-1]) != 0)
			g->name[g->n_gamete++] = rb3_strdup(tmp[k]);
	free(tmp);
	g->sid2g = RB3_MALLOC(int32_t, sid->n_seq);
	for (k = 0; k < sid->n_seq; ++k)
		g->sid2g[k] = rb3_gtab_find(g, pre[k]);
	for (k = 0; k < sid->n_seq; ++k) free(pre[k]);
	free(pre);
	return g;
}

void rb3_gtab_destroy(rb3_gtab_t *g)
{
	int32_t i;
	if (g == 0) return;
	for (i = 0; i < g->n_gamete; ++i) free(g->name[i]);
	free(g->name);
	free(g->sid2g);
	free(g);
}

/*******************
 * event accumulator *
 *******************/

rb3_ps4g_acc_t *rb3_ps4g_acc_init(int64_t bin_size)
{
	rb3_ps4g_acc_t *acc = RB3_CALLOC(rb3_ps4g_acc_t, 1);
	acc->bin_size = bin_size > 0? bin_size : 256;
	return acc;
}

static int rb3_i32_cmp(const void *a, const void *b)
{
	int32_t x = *(const int32_t*)a, y = *(const int32_t*)b;
	return x < y? -1 : x > y? 1 : 0;
}

void rb3_ps4g_acc_add(rb3_ps4g_acc_t *acc, int64_t ref_sid, int64_t pos, const int32_t *gametes, int32_t n_gametes)
{
	int32_t tmp[256], n, m, i;
	int64_t bin, ei, gi;
	if (ref_sid < 0 || n_gametes <= 0) return;
	n = n_gametes < 256? n_gametes : 256;
	memcpy(tmp, gametes, n * sizeof(int32_t));
	qsort(tmp, n, sizeof(int32_t), rb3_i32_cmp);
	m = 0;
	for (i = 0; i < n; ++i)
		if (i == 0 || tmp[i] != tmp[m-1]) tmp[m++] = tmp[i];
	bin = pos / acc->bin_size;

	ei = acc->n_ev;
	RB3_GROW(int64_t, acc->ev_ref_sid, ei, acc->m_ev);
	acc->ev_bin  = RB3_REALLOC(int64_t, acc->ev_bin,  acc->m_ev);
	acc->ev_goff = RB3_REALLOC(int32_t, acc->ev_goff, acc->m_ev);
	acc->ev_glen = RB3_REALLOC(int32_t, acc->ev_glen, acc->m_ev);
	acc->ev_ref_sid[ei] = ref_sid;
	acc->ev_bin[ei] = bin;
	acc->ev_goff[ei] = (int32_t)acc->n_garena;
	acc->ev_glen[ei] = m;
	acc->n_ev = ei + 1;

	gi = acc->n_garena;
	RB3_GROW(int32_t, acc->garena, gi + m - 1, acc->m_garena);
	memcpy(acc->garena + gi, tmp, m * sizeof(int32_t));
	acc->n_garena += m;
}

void rb3_ps4g_acc_destroy(rb3_ps4g_acc_t *acc)
{
	if (acc == 0) return;
	free(acc->ev_ref_sid); free(acc->ev_bin);
	free(acc->ev_goff); free(acc->ev_glen);
	free(acc->garena);
	free(acc);
}

/*******************
 * BED label reader *
 *******************/

static int64_t bed_resolve_contig(const rb3_sid_t *sid, const char *ref_prefix, const char *chrom)
{
	int64_t k;
	size_t plen = ref_prefix? strlen(ref_prefix) : 0;
	for (k = 0; k < sid->n_seq; ++k) { // try the stripped (PS4G-style) contig name, reference sequences only
		if (ref_prefix == 0 || strncmp(sid->name[k], ref_prefix, plen) != 0) continue; // must be a reference sequence, else e.g. "chr1" could match a assembly's own "Oh43_chr1"
		if (strcmp(contig_suffix(sid->name[k]), chrom) == 0) return k;
	}
	for (k = 0; k < sid->n_seq; ++k) // fall back to a literal sequence-name match
		if (strcmp(sid->name[k], chrom) == 0)
			return k;
	return -1;
}

static int bed_region_cmp(const void *a, const void *b)
{
	const rb3_bed_region_t *x = (const rb3_bed_region_t*)a, *y = (const rb3_bed_region_t*)b;
	if (x->ref_sid != y->ref_sid) return x->ref_sid < y->ref_sid? -1 : 1;
	if (x->start != y->start) return x->start < y->start? -1 : 1;
	return 0;
}

rb3_bed_t *rb3_bed_read(const char *fn, const rb3_gtab_t *gtab, const rb3_sid_t *sid, const char *ref_prefix)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t dret;
	rb3_bed_t *bed;
	int64_t m_r = 0;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	bed = RB3_CALLOC(rb3_bed_t, 1);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		char *fld[5]; // chrom, start, end, sampleA, sampleB
		int32_t nf = 0;
		char *p, *q;
		int64_t ref_sid;
		int32_t gA, gB;
		if (str.l == 0 || str.s[0] == '#') continue;
		if (strncmp(str.s, "track", 5) == 0 || strncmp(str.s, "browser", 7) == 0) continue;
		for (p = q = str.s; nf < 5; ++p) {
			if (*p == '\t' || *p == 0) {
				int done = (*p == 0);
				*p = 0;
				fld[nf++] = q;
				q = p + 1;
				if (done) break;
			}
		}
		if (nf < 4) continue; // need at least chrom, start, end, sampleA
		ref_sid = bed_resolve_contig(sid, ref_prefix, fld[0]);
		if (ref_sid < 0) {
			if (rb3_verbose >= 2) fprintf(stderr, "WARNING: --label-bed: contig '%s' not found in the index; skipped\n", fld[0]);
			continue;
		}
		gA = rb3_gtab_find(gtab, fld[3]);
		if (gA < 0) {
			if (rb3_verbose >= 2) fprintf(stderr, "WARNING: --label-bed: sample '%s' not found in the index; skipped\n", fld[3]);
			continue;
		}
		gB = nf >= 5 && fld[4][0]? rb3_gtab_find(gtab, fld[4]) : gA;
		if (gB < 0) gB = gA;
		RB3_GROW(rb3_bed_region_t, bed->r, bed->n_r, m_r);
		bed->r[bed->n_r].ref_sid = ref_sid;
		bed->r[bed->n_r].start = atol(fld[1]);
		bed->r[bed->n_r].end = atol(fld[2]);
		bed->r[bed->n_r].gA = gA;
		bed->r[bed->n_r].gB = gB;
		++bed->n_r;
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	qsort(bed->r, bed->n_r, sizeof(rb3_bed_region_t), bed_region_cmp);
	return bed;
}

void rb3_bed_destroy(rb3_bed_t *b)
{
	if (b == 0) return;
	free(b->r);
	free(b);
}

// binary search for a region on ref_sid whose [start,end) contains pos
static int bed_lookup(const rb3_bed_t *bed, int64_t ref_sid, int64_t pos, int32_t *gA, int32_t *gB)
{
	int64_t lo = 0, hi = bed->n_r - 1, best = -1;
	if (bed == 0 || bed->n_r == 0) return 0;
	while (lo <= hi) { // find the last region with (ref_sid,start) <= (ref_sid,pos)
		int64_t mid = (lo + hi) / 2;
		const rb3_bed_region_t *r = &bed->r[mid];
		int before = r->ref_sid < ref_sid || (r->ref_sid == ref_sid && r->start <= pos);
		if (before) { best = mid; lo = mid + 1; }
		else hi = mid - 1;
	}
	while (best >= 0 && bed->r[best].ref_sid == ref_sid) {
		if (bed->r[best].start <= pos && pos < bed->r[best].end) {
			*gA = bed->r[best].gA, *gB = bed->r[best].gB;
			return 1;
		}
		if (bed->r[best].end <= pos) break; // regions assumed non-overlapping; no earlier region can help
		--best;
	}
	return 0;
}

/*******************
 * finalize: sort + write PS4G / npy *
 *******************/

typedef struct { int64_t ref_sid, bin; int32_t goff, glen; } ev_t;

static const int32_t *g_sort_arena; // scratch for ev_cmp during the single-threaded finalize sort

static int ev_cmp(const void *a, const void *b)
{
	const ev_t *x = (const ev_t*)a, *y = (const ev_t*)b;
	int32_t i;
	if (x->ref_sid != y->ref_sid) return x->ref_sid < y->ref_sid? -1 : 1;
	if (x->bin != y->bin) return x->bin < y->bin? -1 : 1;
	if (x->glen != y->glen) return x->glen < y->glen? -1 : 1;
	for (i = 0; i < x->glen; ++i) {
		int32_t xa = g_sort_arena[x->goff + i], ya = g_sort_arena[y->goff + i];
		if (xa != ya) return xa < ya? -1 : 1;
	}
	return 0;
}

typedef struct { int64_t ref_sid, bin; int32_t goff, glen, count; } row_t; // one unique (ref_sid,bin,gameteSet)

static void write_npy_header(FILE *fp, int64_t rows, int64_t cols)
{
	char dict[256];
	int dict_len, pad, i;
	int64_t base_len, header_len;
	uint8_t ver[2] = {1, 0}, hlen[2];
	dict_len = snprintf(dict, sizeof(dict), "{'descr': '<i4', 'fortran_order': False, 'shape': (%ld, %ld), }", (long)rows, (long)cols);
	base_len = 10 + dict_len; // 6-byte magic + 2-byte version + 2-byte header-length field
	pad = (int)((64 - ((base_len + 1) % 64)) % 64);
	header_len = dict_len + pad + 1;
	hlen[0] = (uint8_t)(header_len & 0xff);
	hlen[1] = (uint8_t)((header_len >> 8) & 0xff);
	fwrite("\x93NUMPY", 1, 6, fp);
	fwrite(ver, 1, 2, fp);
	fwrite(hlen, 1, 2, fp);
	fwrite(dict, 1, dict_len, fp);
	for (i = 0; i < pad; ++i) fputc(' ', fp);
	fputc('\n', fp);
}

void rb3_ps4g_npy_finalize(rb3_ps4g_acc_t *acc, const rb3_gtab_t *gtab, const rb3_sid_t *sid,
							const rb3_bed_t *bed, int npy_binary,
							const char *ps4g_fn, const char *npy_fn, const char *cli_command)
{
	ev_t *ev;
	int64_t i, n_row = 0, m_row = 0;
	row_t *row = 0;
	int64_t *gamete_total;
	int64_t total_unique_counts = 0;

	if ((ps4g_fn == 0 && npy_fn == 0) || acc == 0) return;

	ev = RB3_MALLOC(ev_t, acc->n_ev);
	for (i = 0; i < acc->n_ev; ++i) {
		ev[i].ref_sid = acc->ev_ref_sid[i];
		ev[i].bin = acc->ev_bin[i];
		ev[i].goff = acc->ev_goff[i];
		ev[i].glen = acc->ev_glen[i];
	}
	g_sort_arena = acc->garena;
	qsort(ev, acc->n_ev, sizeof(ev_t), ev_cmp);

	// collapse runs of identical (ref_sid,bin,gameteSet) into unique rows with counts
	for (i = 0; i < acc->n_ev; ++i) {
		if (n_row > 0 && row[n_row-1].ref_sid == ev[i].ref_sid && row[n_row-1].bin == ev[i].bin &&
			row[n_row-1].glen == ev[i].glen &&
			memcmp(acc->garena + row[n_row-1].goff, acc->garena + ev[i].goff, ev[i].glen * sizeof(int32_t)) == 0) {
			row[n_row-1].count++;
		} else {
			RB3_GROW(row_t, row, n_row, m_row);
			row[n_row].ref_sid = ev[i].ref_sid, row[n_row].bin = ev[i].bin;
			row[n_row].goff = ev[i].goff, row[n_row].glen = ev[i].glen;
			row[n_row].count = 1;
			++n_row;
		}
	}
	free(ev);

	gamete_total = RB3_CALLOC(int64_t, gtab->n_gamete);
	for (i = 0; i < n_row; ++i) {
		int32_t k;
		total_unique_counts += row[i].count;
		for (k = 0; k < row[i].glen; ++k)
			gamete_total[acc->garena[row[i].goff + k]] += row[i].count;
	}

	if (ps4g_fn) {
		FILE *fp = fopen(ps4g_fn, "w");
		if (fp == 0) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: failed to write PS4G file '%s'\n", ps4g_fn);
		} else {
			int32_t g;
			fprintf(fp, "#PS4G\n#version=2.0\n");
			fprintf(fp, "#Command: %s\n", cli_command? cli_command : "");
			fprintf(fp, "#TotalUniqueCounts: %ld\n", (long)total_unique_counts);
			fprintf(fp, "#gamete\tgameteIndex\tcount\n");
			for (g = 0; g < gtab->n_gamete; ++g)
				fprintf(fp, "#%s\t%d\t%ld\n", gtab->name[g], g, (long)gamete_total[g]);
			fprintf(fp, "gameteSet\trefContig\trefPosBinned\tcount\n");
			for (i = 0; i < n_row; ++i) {
				int32_t k;
				const char *cname = contig_suffix(sid->name[row[i].ref_sid]);
				for (k = 0; k < row[i].glen; ++k)
					fprintf(fp, "%s%d", k? "," : "", acc->garena[row[i].goff + k]);
				fprintf(fp, "\t%s\t%ld\t%ld\n", cname, (long)row[i].bin, (long)row[i].count);
			}
			fclose(fp);
		}
	}

	if (npy_fn) {
		// One npy row per PS4G row (contig, bin, gameteSet) -- NOT collapsed across
		// gameteSets sharing a bin. Aggregating different gameteSets into one row
		// would discard exactly the co-occurrence information (which gametes were
		// jointly supported by the same reads) that the imputation model needs.
		int64_t cols = gtab->n_gamete + 2;
		int32_t *mat = RB3_CALLOC(int32_t, n_row * cols);
		FILE *fp;
		char *aux_fn;

		for (i = 0; i < n_row; ++i) {
			int32_t k, val = npy_binary? 1 : (int32_t)row[i].count;
			int32_t gA = -1, gB = -1;
			int64_t pos = row[i].bin * acc->bin_size;
			for (k = 0; k < row[i].glen; ++k)
				mat[i * cols + acc->garena[row[i].goff + k]] = val;
			if (bed) bed_lookup(bed, row[i].ref_sid, pos, &gA, &gB); // diploid training labels, -1 if unlabeled
			mat[i * cols + gtab->n_gamete] = gA;
			mat[i * cols + gtab->n_gamete + 1] = gB;
		}

		fp = fopen(npy_fn, "wb");
		if (fp == 0) {
			if (rb3_verbose >= 1) fprintf(stderr, "ERROR: failed to write numpy file '%s'\n", npy_fn);
		} else {
			write_npy_header(fp, n_row, cols);
			fwrite(mat, sizeof(int32_t), n_row * cols, fp);
			fclose(fp);
		}

		aux_fn = RB3_MALLOC(char, strlen(npy_fn) + 16);
		sprintf(aux_fn, "%s.bins.tsv", npy_fn);
		fp = fopen(aux_fn, "w");
		if (fp) { // a (contig,bin) can repeat across rows: each row is a distinct gameteSet at that bin
			fprintf(fp, "row\tcontig\tbin\n");
			for (i = 0; i < n_row; ++i)
				fprintf(fp, "%ld\t%s\t%ld\n", (long)i, contig_suffix(sid->name[row[i].ref_sid]), (long)row[i].bin);
			fclose(fp);
		}
		sprintf(aux_fn, "%s.gametes.tsv", npy_fn);
		fp = fopen(aux_fn, "w");
		if (fp) {
			int32_t g;
			fprintf(fp, "gameteIndex\tsampleName\n");
			for (g = 0; g < gtab->n_gamete; ++g)
				fprintf(fp, "%d\t%s\n", g, gtab->name[g]);
			fclose(fp);
		}
		free(aux_fn);
		free(mat);
	}

	free(gamete_total);
	free(row);
}
