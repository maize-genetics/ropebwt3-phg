#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pairing.h"

const char *rb3_pe_class_str[5] = { "NONE", "SOLO", "CONC", "RESC", "DISC" };

int rb3_pe_name_len(const char *name)
{
	int l = (int)strlen(name);
	if (l >= 2 && name[l-2] == '/' && (name[l-1] == '1' || name[l-1] == '2')) l -= 2;
	return l;
}

int rb3_pe_names_match(const char *a, const char *b)
{
	int la, lb;
	if (a == 0 || b == 0) return 0;
	la = rb3_pe_name_len(a), lb = rb3_pe_name_len(b);
	return la == lb && strncmp(a, b, la) == 0;
}

static int cmp_i64(const void *x, const void *y)
{
	int64_t a = *(const int64_t*)x, b = *(const int64_t*)y;
	return a < b? -1 : a > b;
}

static int cmp_i32(const void *x, const void *y)
{
	int32_t a = *(const int32_t*)x, b = *(const int32_t*)y;
	return a < b? -1 : a > b;
}

int rb3_pe_isize_estimate(int64_t *span, int32_t n, int32_t min_n, rb3_isize_t *out)
{
	int64_t p25, p75, iqr, lo, hi, cap;
	double s = 0.0, ss = 0.0;
	int32_t i, m = 0;
	memset(out, 0, sizeof(*out));
	if (n < min_n || n <= 0) return 0;
	qsort(span, n, sizeof(int64_t), cmp_i64);
	p25 = span[(int32_t)(.25 * n)], p75 = span[(int32_t)(.75 * n)];
	iqr = p75 - p25;
	lo = p25 - 2 * iqr, hi = p75 + 2 * iqr;
	for (i = 0; i < n; ++i)
		if (span[i] >= lo && span[i] <= hi) s += span[i], ++m;
	out->mean = s / m;
	for (i = 0; i < n; ++i)
		if (span[i] >= lo && span[i] <= hi) ss += (span[i] - out->mean) * (span[i] - out->mean);
	out->sd = sqrt(ss / m);
	out->n = m, out->ok = 1;
	out->max_ins = (int64_t)(out->mean + 4.0 * out->sd + .5);
	cap = p75 + 3 * iqr;
	if (out->max_ins < cap) out->max_ins = cap;
	return 1;
}

static int32_t sort_uniq(int32_t *a, int32_t n)
{
	int32_t i, m;
	if (n <= 1) return n;
	qsort(a, n, sizeof(int32_t), cmp_i32);
	for (i = 1, m = 1; i < n; ++i)
		if (a[i] != a[m-1]) a[m++] = a[i];
	return m;
}

static inline int64_t i64_abs(int64_t x) { return x < 0? -x : x; }

// Projections of mate b lying on the anchor's contig within max_dist of it; median of them.
static int32_t pe_near_proj(const rb3_pe_mate_t *a, const rb3_pe_mate_t *b, int64_t max_dist, int64_t *med)
{
	int64_t v[64];
	int32_t i, j, n = 0, tot = 0;
	for (i = 0; i < b->n_proj; ++i) {
		if (b->proj_rsid[i] != a->rsid || i64_abs(b->proj_pos[i] - a->pos) > max_dist) continue;
		++tot;
		if (n < 64) v[n++] = b->proj_pos[i];
	}
	if (n == 0) return 0;
	for (i = 1; i < n; ++i) { int64_t x = v[i]; for (j = i - 1; j >= 0 && v[j] > x; --j) v[j+1] = v[j]; v[j+1] = x; }
	*med = v[n >> 1];
	return tot;
}

void rb3_pe_resolve(const rb3_pe_mate_t m[2], int64_t max_dist, rb3_pe_res_t *res)
{
	int k, best_a = -1;
	int32_t best_sup = 0;
	int64_t best_pos = -1;
	memset(res, 0, sizeof(*res));
	res->anchor = res->rescued = -1, res->resc_rsid = res->resc_pos = -1;
	if (m[0].rank == 0 && m[1].rank == 0) { res->cls = RB3_PE_NONE; return; }
	if (m[0].rank == 0 || m[1].rank == 0) { res->cls = RB3_PE_SOLO, res->anchor = m[0].rank? 0 : 1; return; }
	if (m[0].rsid == m[1].rsid && i64_abs(m[0].pos - m[1].pos) <= max_dist) {
		res->cls = RB3_PE_CONC, res->anchor = m[1].rank > m[0].rank? 1 : 0;
		return;
	}
	// discordant: try moving a projected mate b next to an anchor a at least as confident.
	// An exact anchor beats a projected one; then more agreeing projections; then mate 0.
	for (k = 0; k < 2; ++k) {
		const rb3_pe_mate_t *a = &m[k], *b = &m[k^1];
		int64_t med;
		int32_t sup;
		if (b->rank != 1 || a->rank < b->rank) continue;
		if ((sup = pe_near_proj(a, b, max_dist, &med)) == 0) continue;
		if (best_a < 0 || a->rank > m[best_a].rank || (a->rank == m[best_a].rank && sup > best_sup))
			best_a = k, best_sup = sup, best_pos = med;
	}
	if (best_a < 0) { res->cls = RB3_PE_DISC; return; }
	res->cls = RB3_PE_RESC, res->anchor = best_a, res->rescued = best_a ^ 1;
	res->resc_rsid = m[best_a].rsid, res->resc_pos = best_pos, res->resc_support = best_sup;
}


int32_t rb3_pe_near_gametes(const rb3_pe_mate_t *b, int64_t rsid, int64_t pos, int64_t max_dist, int32_t *out)
{
	int32_t i, n = 0;
	if (b->proj_gid == 0) return 0;
	for (i = 0; i < b->n_proj; ++i)
		if (b->proj_rsid[i] == rsid && i64_abs(b->proj_pos[i] - pos) <= max_dist)
			out[n++] = b->proj_gid[i];
	return sort_uniq(out, n);
}

int32_t rb3_pe_event_gametes(const int32_t *a, int32_t na, int a_capped,
							 const int32_t *b, int32_t nb, int b_capped, int32_t *out)
{
	int32_t i, j, n, nu, ma, *tb;
	if (nb <= 0 || b_capped) { // nothing usable from the mate: the anchor's own set
		memcpy(out, a, na * sizeof(int32_t));
		return sort_uniq(out, na);
	}
	if (a_capped || na <= 0) { // the anchor's sample is incomplete (or empty): the mate's set
		memcpy(out, b, nb * sizeof(int32_t));
		return sort_uniq(out, nb);
	}
	memcpy(out, a, na * sizeof(int32_t));
	ma = sort_uniq(out, na);
	tb = (int32_t*)malloc(nb * sizeof(int32_t));
	memcpy(tb, b, nb * sizeof(int32_t));
	nu = sort_uniq(tb, nb);
	for (i = j = n = 0; i < ma && j < nu; ) { // in-place: out[n] is written at or behind out[i]
		if (out[i] < tb[j]) ++i;
		else if (out[i] > tb[j]) ++j;
		else out[n++] = out[i], ++i, ++j;
	}
	free(tb);
	if (n == 0) { // mates disagree on every founder: keep the anchor's evidence
		memcpy(out, a, na * sizeof(int32_t));
		return sort_uniq(out, na);
	}
	return n;
}
