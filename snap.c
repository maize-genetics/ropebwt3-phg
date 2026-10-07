/* snap.c -- affine-gap flank aligner used by `extract --lift` to snap projected
 * interval ends to the homologous base. See snap.h. */
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "rb3priv.h"
#include "snap.h"

#define SNAP_NEG_INF (INT32_MIN / 4)

void rb3_snap_opt_init(rb3_snap_opt_t *o)
{
	o->match = 1, o->mismatch = 4, o->gap_open = 6, o->gap_ext = 1;
	o->min_aligned = 16, o->min_frac = 0.5, o->min_id = 0.85, o->margin = 5, o->excl = 8, o->clip5 = -1, o->no_seed = 0, o->split = -1;
}

static int snap_align_dp(const uint8_t *q, int32_t n, const uint8_t *t, int32_t m, const rb3_snap_opt_t *o, rb3_snap_res_t *r)
{
	int32_t i, j, *H0, *H1, *F, best_j = -1, oe = o->gap_open + o->gap_ext;
	uint8_t *dir;
	char *col;
	int32_t n_col = 0, k, first = -1, last = -1, need;
	memset(r, 0, sizeof(*r));
	r->second = INT32_MIN, r->score = INT32_MIN, r->q_first = r->t_first = r->q_last = r->t_last = r->t_end = -1;
	r->t_left = r->t_right = r->q_left = r->q_right = -1;
	need = (int32_t)(o->min_frac * n);
	if (need < o->min_aligned) need = o->min_aligned;
	if (n <= 0 || m <= 0 || n < o->min_aligned || m < need) return (r->status = RB3_SNAP_SHORT);

	H0 = RB3_MALLOC(int32_t, m + 1), H1 = RB3_MALLOC(int32_t, m + 1), F = RB3_MALLOC(int32_t, m + 1);
	dir = RB3_MALLOC(uint8_t, (size_t)(n + 1) * (m + 1));
	for (j = 0; j <= m; ++j) H0[j] = 0, F[j] = SNAP_NEG_INF, dir[j] = 0; // free leading target
	for (i = 1; i <= n; ++i) {
		uint8_t *d = dir + (size_t)i * (m + 1);
		int32_t E = SNAP_NEG_INF, qc = q[i-1];
		// column 0: query bases [0,i) deleted (vertical gap from H[0][0] = 0)
		F[0] = i == 1? -oe : F[0] - o->gap_ext;
		H1[0] = F[0];
		d[0] = 2 | (i > 1? 8 : 0);
		for (j = 1; j <= m; ++j) {
			int32_t eo = H1[j-1] - oe, ee = E - o->gap_ext;
			int32_t fo = H0[j] - oe, fe = F[j] - o->gap_ext;
			int32_t h, sc = (qc == t[j-1] && qc >= 1 && qc <= 4)? o->match : -o->mismatch;
			uint8_t b = 0;
			if (ee > eo) E = ee, b |= 4; else E = eo;
			if (fe > fo) F[j] = fe, b |= 8; else F[j] = fo;
			h = H0[j-1] + sc; // prefer the diagonal on ties
			if (o->clip5 >= 0 && i > 1 && -o->clip5 + sc > h) h = -o->clip5 + sc, b |= 16; // start here, prefix clipped
			if (E > h) h = E, b = (b & ~19) | 1;
			if (F[j] > h) h = F[j], b = (b & ~19) | 2;
			H1[j] = h, d[j] = b;
		}
		{ int32_t *tmp = H0; H0 = H1; H1 = tmp; }
	}
	// H0 = last row: the score of the flank aligned ending at target position j
	for (j = 1; j <= m; ++j)
		if (best_j < 0 || H0[j] > H0[best_j]) best_j = j; // leftmost best end
	r->score = H0[best_j], r->t_end = best_j;
	for (j = 1; j <= m; ++j)
		if ((j - best_j >= o->excl || best_j - j >= o->excl) && H0[j] > r->second) r->second = H0[j];

	// traceback
	col = RB3_MALLOC(char, n + m + 2);
	{
		int state = 0; // 0 H, 1 E (horizontal: founder-only base), 2 F (vertical: deleted flank base)
		i = n, j = best_j;
		while (i > 0) {
			uint8_t b = dir[(size_t)i * (m + 1) + j];
			if (state == 0) {
				int src = b & 3;
				if (src == 0 && o->split >= 0) { // traceback runs backwards: the first M seen left of split is the last
					if (i - 1 >= o->split) r->t_right = j - 1, r->q_right = i - 1;
					else if (r->t_left < 0) r->t_left = j - 1, r->q_left = i - 1;
				}
				if (src == 0 && (b & 16)) { // the alignment starts here; query [0,i-1) clipped
					col[n_col++] = (q[i-1] == t[j-1] && q[i-1] >= 1 && q[i-1] <= 4)? '=' : 'X';
					if (last < 0) r->q_last = i - 1, r->t_last = j - 1, last = 1;
					r->q_first = i - 1, r->t_first = j - 1;
					break;
				} else if (src == 0) {
					col[n_col++] = (q[i-1] == t[j-1] && q[i-1] >= 1 && q[i-1] <= 4)? '=' : 'X';
					if (last < 0) r->q_last = i - 1, r->t_last = j - 1, last = 1;
					r->q_first = i - 1, r->t_first = j - 1;
					--i, --j;
				} else state = src;
			} else if (state == 1) {
				col[n_col++] = 'I';
				state = (b & 4)? 1 : 0;
				--j;
			} else {
				col[n_col++] = 'D';
				state = (b & 8)? 2 : 0;
				--i;
			}
		}
	}
	// stats over the core [first M, last M] (col[] is in reverse order)
	for (k = 0; k < n_col; ++k) if (col[k] == '=' || col[k] == 'X') { if (first < 0) first = k; last = k; }
	if (first >= 0) {
		int32_t n_open = 0;
		for (k = first; k <= last; ++k) {
			if (col[k] == '=' || col[k] == 'X') { ++r->n_aligned; r->n_match += col[k] == '='; }
			else if (k == first || col[k-1] != col[k]) ++n_open; // gap-compressed: one event per gap
			++r->n_core;
		}
		r->identity = (double)r->n_match / (r->n_aligned + n_open);
	}
	free(col); free(dir); free(H0); free(H1); free(F);
	if (r->n_aligned < need || r->identity < o->min_id) r->status = RB3_SNAP_LOWID;
	else if (r->second != INT32_MIN && r->second >= r->score - o->margin) r->status = RB3_SNAP_AMBIG;
	else r->status = RB3_SNAP_OK;
	return r->status;
}

/* Seed prefilter for wide windows (m > SNAP_SEED_RATIO * n): the full DP over a 12-25 kb
 * window costs ~n*m cells, and the wide retry / walk-inward steps run several of them
 * per end. Instead, collect the target diagonals sharing an exact k-mer (k = SNAP_SEED_K)
 * with the query, cut the target into the regions [diag - n, diag + 2n) around them,
 * and run the same DP on those regions joined by runs of N longer than any gap worth
 * opening (N never matches, a gap across costs > the best possible score), so the best
 * and second-best placements -- and thus the ambiguity call -- are those of the full
 * DP whenever every competitive placement shares a k-mer with the flank. A flank at
 * the default --snap-min-id with no exact SNAP_SEED_K-mer anywhere would need a
 * difference every <SNAP_SEED_K bases; such an alignment is reported as low identity. */
#define SNAP_SEED_K 10
#define SNAP_SEED_RATIO 8

typedef struct { int32_t c0, t0, len; } snap_seg_t; // concatenated start, target start, length

static int cmp_u64(const void *a, const void *b)
{
	uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
	return x < y? -1 : x > y;
}

int rb3_snap_align(const uint8_t *q, int32_t n, const uint8_t *t, int32_t m, const rb3_snap_opt_t *o, rb3_snap_res_t *r)
{
	const int k = SNAP_SEED_K;
	const uint32_t mask = (1U << 2*k) - 1;
	uint64_t *qk; // (kmer << 32) | query position
	int32_t nqk = 0, i, j, l, ns = 0, sep, cl = 0, *lo_hi, nlh = 0;
	uint32_t h;
	snap_seg_t *seg;
	uint8_t *ct;
	int st;
	if (o->no_seed || m <= SNAP_SEED_RATIO * n || n < k || n < o->min_aligned) return snap_align_dp(q, n, t, m, o, r);
	qk = RB3_MALLOC(uint64_t, n);
	for (i = 0, l = 0, h = 0; i < n; ++i) {
		if (q[i] < 1 || q[i] > 4) { l = 0, h = 0; continue; }
		h = (h << 2 | (q[i] - 1)) & mask;
		if (++l >= k) qk[nqk++] = (uint64_t)h << 32 | (uint32_t)(i - k + 1);
	}
	qsort(qk, nqk, 8, cmp_u64);
	lo_hi = RB3_MALLOC(int32_t, 2 * 64);
	{ int32_t m_lh = 64;
	for (j = 0, l = 0, h = 0; j < m; ++j) {
		int32_t a, b;
		if (t[j] < 1 || t[j] > 4) { l = 0, h = 0; continue; }
		h = (h << 2 | (t[j] - 1)) & mask;
		if (++l < k) continue;
		a = 0, b = nqk; // first query k-mer >= h
		while (a < b) { int32_t mid = (a + b) >> 1; if ((uint32_t)(qk[mid] >> 32) < h) a = mid + 1; else b = mid; }
		for (; a < nqk && (uint32_t)(qk[a] >> 32) == h; ++a) {
			int32_t diag = (j - k + 1) - (int32_t)(uint32_t)qk[a], s0 = diag - n, e0 = diag + 2 * n;
			if (s0 < 0) s0 = 0;
			if (e0 > m) e0 = m;
			if (nlh > 0 && s0 <= lo_hi[2*nlh-1]) { if (e0 > lo_hi[2*nlh-1]) lo_hi[2*nlh-1] = e0; } // diagonals arrive in near-sorted order
			else {
				if (nlh == m_lh) m_lh <<= 1, lo_hi = RB3_REALLOC(int32_t, lo_hi, 2 * m_lh);
				lo_hi[2*nlh] = s0, lo_hi[2*nlh+1] = e0, ++nlh;
			}
		}
	}
	}
	free(qk);
	if (nlh == 0) { // no shared k-mer anywhere: report as the DP would for an unrelated flank
		free(lo_hi);
		memset(r, 0, sizeof(*r));
		r->score = r->second = INT32_MIN, r->q_first = r->t_first = r->q_last = r->t_last = r->t_end = -1;
		return (r->status = RB3_SNAP_LOWID);
	}
	// merge (regions from different query k-mers may arrive out of order)
	{
		int32_t a, w = 0;
		for (a = 1; a < nlh; ++a) { // insertion sort by start; nlh is small
			int32_t s0 = lo_hi[2*a], e0 = lo_hi[2*a+1], b = a - 1;
			while (b >= 0 && lo_hi[2*b] > s0) lo_hi[2*b+2] = lo_hi[2*b], lo_hi[2*b+3] = lo_hi[2*b+1], --b;
			lo_hi[2*b+2] = s0, lo_hi[2*b+3] = e0;
		}
		for (a = 1; a < nlh; ++a) {
			if (lo_hi[2*a] <= lo_hi[2*w+1]) { if (lo_hi[2*a+1] > lo_hi[2*w+1]) lo_hi[2*w+1] = lo_hi[2*a+1]; }
			else ++w, lo_hi[2*w] = lo_hi[2*a], lo_hi[2*w+1] = lo_hi[2*a+1];
		}
		nlh = w + 1;
	}
	sep = 2 * n + o->gap_open + 8;
	for (i = 0; i < nlh; ++i) cl += lo_hi[2*i+1] - lo_hi[2*i] + (i? sep : 0);
	if (cl >= m) { free(lo_hi); return snap_align_dp(q, n, t, m, o, r); } // no saving
	ct = RB3_MALLOC(uint8_t, cl);
	seg = RB3_MALLOC(snap_seg_t, nlh);
	for (i = 0, cl = 0; i < nlh; ++i) {
		int32_t len = lo_hi[2*i+1] - lo_hi[2*i];
		if (i) { memset(ct + cl, 5, sep); cl += sep; }
		seg[ns].c0 = cl, seg[ns].t0 = lo_hi[2*i], seg[ns].len = len, ++ns;
		memcpy(ct + cl, t + lo_hi[2*i], len);
		cl += len;
	}
	free(lo_hi);
	st = snap_align_dp(q, n, ct, cl, o, r);
	{ // map concatenated target coordinates back
		int32_t *f[6];
		f[0] = &r->t_first, f[1] = &r->t_last, f[2] = &r->t_end, f[3] = &r->t_left, f[4] = &r->t_right, f[5] = 0;
		for (i = 0; f[i]; ++i) {
			int32_t x = *f[i], s2, off = f[i] == &r->t_end? 1 : 0; // t_end is exclusive
			if (x < 0) continue;
			for (s2 = ns - 1; s2 > 0 && seg[s2].c0 > x - off; --s2) {}
			*f[i] = x - seg[s2].c0 + seg[s2].t0;
		}
	}
	free(ct); free(seg);
	return st;
}
