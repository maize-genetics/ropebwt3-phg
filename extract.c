/* extract.c -- `ropebwt3 extract`: coordinate-based sequence extraction from an FMD
 * via the sampled inverse suffix array (isa.c), optionally after projecting a
 * reference interval into a target assembly with the refmap liftover (lift.c).
 *
 * Coordinates are 0-based, half-open [start,end) on the named sequence's FORWARD
 * strand (the sequence as it appeared in the FASTA given to `ropebwt3 build`).
 * Strand '-' returns the reverse complement of that same interval. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "io.h"
#include "ketopt.h"
#include "kalloc.h"
#include "kthread.h"
#include "lift.h"
#include "ps4g.h"
#include "isa.h"
#include "snap.h"

#define EXT_BATCH 4096

enum { ST_OK = 0, ST_INPUT_ERR = 1, ST_PROJ_FAIL = 2 };

typedef struct {
	char *qid;
	char *region;       // the original region string / TSV line (for messages)
	// direct mode, or the projected result
	int64_t tsid;       // target sequence index (.len.gz order), -1 if none
	int64_t st, en;
	int strand;         // +1 / -1
	// projected mode input
	int64_t rsid, rst, ren;
	char *asm_name;
	int clamped;        // the interval was clamped (input --clamp, or a projected end past the contig end)
	int meth[2];        // projection method of the start / end
	int snap[2];        // SNAP_* label of the start / end
	int raw_ok[2], raw_sign[2]; // raw (unsnapped) lift projection of each reference boundary
	int64_t raw_tsid[2], raw_pos[2]; // ... contig and forward fencepost (see print_proj_row)
	int64_t shift[2];   // bp the snap moved the founder base of the reference start / end (forward coords)
	int status;         // ST_*
	const char *msg;    // failure reason (static string or msgbuf), or 0
	char msgbuf[96];
	char *seq;          // extracted sequence (ACGTN), NUL-terminated
	int64_t steps;      // LF steps taken
} ext_query_t;

typedef struct {
	int32_t n_threads, line_width, clamp, project_only, use_mmap, min_support, lift_k;
	int32_t snap, snap_flank;
	int64_t snap_r1, snap_r2;
	rb3_snap_opt_t sopt;
	char *lift_fn, *ref_prefix, *to_asm, *isa_fn, *tsv_fn;
	int64_t lift_win, lift_mad, local_win, local_tol;
	int64_t max_span_excess;   // span check: fail if prj-ref > this AND prj > ratio*ref; 0 = off
	double max_span_ratio;
} ext_opt_t;

typedef struct {
	const ext_opt_t *opt;
	rb3_fmi_t f;
	rb3_sid_t *sid;
	rb3_isa_t *isa;
	int64_t *name_ord;  // sequence indices sorted by name
	rb3_lift_t *lift;
	rb3_lift_ridx_t *ridx;
	rb3_gtab_t *gtab;
	void **km;          // per-thread arena for the projection
	ext_query_t *q;     // current batch
} ext_ctx_t;

/***************
 * name lookup *
 ***************/

static const rb3_sid_t *g_sort_sid;
static int name_ord_cmp(const void *a, const void *b)
{
	return strcmp(g_sort_sid->name[*(const int64_t*)a], g_sort_sid->name[*(const int64_t*)b]);
}

static int64_t *name_index_build(const rb3_sid_t *sid)
{
	int64_t i, *ord = RB3_MALLOC(int64_t, sid->n_seq);
	for (i = 0; i < sid->n_seq; ++i) ord[i] = i;
	g_sort_sid = sid;
	qsort(ord, sid->n_seq, sizeof(int64_t), name_ord_cmp);
	return ord;
}

static int64_t name_lookup(const rb3_sid_t *sid, const int64_t *ord, const char *name)
{
	int64_t lo = 0, hi = sid->n_seq - 1;
	while (lo <= hi) {
		int64_t mid = lo + ((hi - lo) >> 1);
		int c = strcmp(sid->name[ord[mid]], name);
		if (c == 0) return ord[mid];
		else if (c < 0) lo = mid + 1;
		else hi = mid - 1;
	}
	return -1;
}

/***********
 * parsing *
 ***********/

// strict non-negative decimal integer, the whole string
static int parse_coord(const char *s, int64_t *x)
{
	char *end;
	long long v;
	if (s == 0 || *s < '0' || *s > '9') return -1;
	errno = 0;
	v = strtoll(s, &end, 10);
	if (errno || *end != 0) return -1;
	*x = v;
	return 0;
}

static void q_fail(ext_query_t *q, int status, const char *msg)
{
	if (q->status == ST_OK) q->status = status, q->msg = msg;
}

// CONTIG:START-END[:STRAND]; contig names may themselves contain ':'
static int parse_region(const char *reg, char **contig, int64_t *st, int64_t *en, int *strand)
{
	char *s = rb3_strdup(reg), *p, *dash;
	size_t l = strlen(s);
	*strand = 0;
	if (l >= 2 && s[l-2] == ':' && (s[l-1] == '+' || s[l-1] == '-')) {
		*strand = s[l-1] == '+'? 1 : -1;
		s[l-2] = 0;
	}
	if ((p = strrchr(s, ':')) == 0 || p == s) { free(s); return -1; }
	*p++ = 0;
	if ((dash = strchr(p, '-')) == 0) { free(s); return -1; }
	*dash++ = 0;
	if (parse_coord(p, st) < 0 || parse_coord(dash, en) < 0) { free(s); return -1; }
	*contig = s; // s now holds just the contig
	return 0;
}

static void split_tabs(char *line, char **col, int *n_col, int max_col)
{
	char *p = line;
	*n_col = 0;
	while (*n_col < max_col) {
		col[(*n_col)++] = p;
		if ((p = strchr(p, '\t')) == 0) break;
		*p++ = 0;
	}
}

// Fill the coordinate fields of a direct query; contig lookup + range checks.
static void q_set_direct(const ext_ctx_t *c, ext_query_t *q, const char *contig, int64_t st, int64_t en, int strand)
{
	int64_t len;
	q->strand = strand;
	q->tsid = name_lookup(c->sid, c->name_ord, contig);
	if (q->tsid < 0) { q_fail(q, ST_INPUT_ERR, "unknown contig"); return; }
	len = c->sid->len[q->tsid];
	if (c->opt->clamp) {
		int64_t s0 = st, e0 = en;
		if (st > len) st = len;
		if (en > len) en = len;
		q->clamped = (s0 != st || e0 != en);
	}
	q->st = st, q->en = en;
	if (st >= en) q_fail(q, ST_INPUT_ERR, en == st? "empty interval (start == end)" : "empty interval (start > end)");
	else if (en > len) q_fail(q, ST_INPUT_ERR, "interval extends past the contig end (use --clamp to clip)");
}

static void q_set_ref(const ext_ctx_t *c, ext_query_t *q, const char *ref_contig, int64_t st, int64_t en, const char *asm_name)
{
	int64_t len;
	q->asm_name = rb3_strdup(asm_name? asm_name : "");
	q->rsid = name_lookup(c->sid, c->name_ord, ref_contig);
	q->tsid = -1;
	if (q->rsid < 0) { q_fail(q, ST_INPUT_ERR, "unknown reference contig"); return; }
	if (c->opt->ref_prefix && strncmp(ref_contig, c->opt->ref_prefix, strlen(c->opt->ref_prefix)) != 0) {
		q_fail(q, ST_INPUT_ERR, "contig does not start with --ref-prefix"); return;
	}
	if (asm_name == 0 || *asm_name == 0) { q_fail(q, ST_INPUT_ERR, "no target assembly (5th TSV column or --to)"); return; }
	if (rb3_gtab_find(c->gtab, asm_name) < 0) { q_fail(q, ST_INPUT_ERR, "unknown target assembly"); return; }
	len = c->sid->len[q->rsid];
	if (c->opt->clamp) {
		int64_t s0 = st, e0 = en;
		if (st > len) st = len;
		if (en > len) en = len;
		q->clamped = (s0 != st || e0 != en);
	}
	q->rst = st, q->ren = en;
	if (st >= en) q_fail(q, ST_INPUT_ERR, en == st? "empty interval (start == end)" : "empty interval (start > end)");
	else if (en > len) q_fail(q, ST_INPUT_ERR, "interval extends past the contig end (use --clamp to clip)");
}

static void q_from_region(const ext_ctx_t *c, ext_query_t *q, const char *reg)
{
	char *contig = 0;
	int64_t st, en;
	int strand;
	memset(q, 0, sizeof(*q));
	q->qid = rb3_strdup(reg), q->region = rb3_strdup(reg), q->tsid = q->rsid = -1;
	if (parse_region(reg, &contig, &st, &en, &strand) < 0) {
		q_fail(q, ST_INPUT_ERR, "cannot parse region (want CONTIG:START-END[:+|-], 0-based half-open)");
		return;
	}
	if (c->lift) {
		if (strand != 0) q_fail(q, ST_INPUT_ERR, "a strand suffix is not allowed with --lift (the projection decides it)");
		else q_set_ref(c, q, contig, st, en, c->opt->to_asm);
	} else q_set_direct(c, q, contig, st, en, strand < 0? -1 : 1);
	free(contig);
}

// returns 1 if the line is a header (only allowed as the first data line)
static int q_from_tsv(const ext_ctx_t *c, ext_query_t *q, char *line, int first)
{
	char *col[8];
	int n_col;
	int64_t st, en;
	memset(q, 0, sizeof(*q));
	q->region = rb3_strdup(line), q->tsid = q->rsid = -1;
	split_tabs(line, col, &n_col, 8);
	q->qid = rb3_strdup(col[0]);
	if (n_col < 4) { q_fail(q, ST_INPUT_ERR, "TSV line has fewer than 4 columns"); return 0; }
	if (parse_coord(col[2], &st) < 0 || parse_coord(col[3], &en) < 0) {
		if (first) return 1;
		q_fail(q, ST_INPUT_ERR, "start/end is not a non-negative integer");
		return 0;
	}
	if (c->lift) {
		q_set_ref(c, q, col[1], st, en, n_col >= 5 && *col[4]? col[4] : c->opt->to_asm);
	} else {
		int strand = 1;
		if (n_col >= 5) {
			if (strcmp(col[4], "+") == 0) strand = 1;
			else if (strcmp(col[4], "-") == 0) strand = -1;
			else { q_fail(q, ST_INPUT_ERR, "strand column must be '+' or '-'"); return 0; }
		}
		q_set_direct(c, q, col[1], st, en, strand);
	}
	return 0;
}

/**************
 * projection *
 **************/

static const char *meth_name(int m)
{
	static const char *names[] = { "none", "local", "local1", "median", "breakpoint", "identity" };
	return m >= 0 && m <= 5? names[m] : "none";
}

enum { SNAP_OFF = 0, SNAP_SNAP, SNAP_WIDE, SNAP_AMBIG, SNAP_LOWID, SNAP_SHORT, SNAP_IDENTITY, SNAP_WALK,
	   SNAP_R_SNAP, SNAP_R_WIDE, SNAP_R_WALK };

static const char *snap_name(int m)
{
	static const char *names[] = { "proj", "snap", "snap_wide", "proj_ambiguous", "proj_lowid", "proj_short", "identity", "snap_walk",
								   "snap_x", "snap_x_wide", "snap_x_walk" };
	return m >= 0 && m <= 10? names[m] : "proj";
}

// An oriented founder window extracted once per end and reused by the wide retry and
// every walk step (each would otherwise pay its own LF walk of ~R2 bases).
typedef struct { int64_t lo, hi; uint8_t *buf; } snap_cache_t;

static void snap_cache_fill(const ext_ctx_t *c, ext_query_t *q, int64_t tsid, int sign, int64_t lo, int64_t hi, snap_cache_t *ca)
{
	int64_t L = c->sid->len[tsid], r;
	if (lo < 0) lo = 0;
	if (hi > L) hi = L;
	if (ca->buf || hi <= lo) return;
	ca->buf = RB3_MALLOC(uint8_t, hi - lo);
	if ((r = rb3_isa_extract(&c->f, c->isa, tsid * 2, sign > 0? lo : L - hi, sign > 0? hi : L - lo, ca->buf)) < 0) {
		free(ca->buf); ca->buf = 0; return;
	}
	q->steps += r;
	if (sign < 0) rb3_revcomp6(hi - lo, ca->buf);
	ca->lo = lo, ca->hi = hi;
}

// Extract the oriented founder window [wlo,whi) (clamped to [0,L)) of sequence tsid
// (from the cache when it covers it) and align q[0,nq) to it. Returns the status (RB3_SNAP_SHORT on an extraction
// failure); *wlo_out = the clamped window start.
static int snap_window(const ext_ctx_t *c, ext_query_t *q, int64_t tsid, int sign, int64_t wlo, int64_t whi,
					   const uint8_t *qs, int32_t nq, const rb3_snap_opt_t *so, rb3_snap_res_t *res, int64_t *wlo_out, int reversed,
					   const snap_cache_t *ca)
{
	int64_t L = c->sid->len[tsid], flo, fhi, r;
	uint8_t *win;
	if (wlo < 0) wlo = 0;
	if (whi > L) whi = L;
	*wlo_out = wlo;
	if (whi <= wlo) { memset(res, 0, sizeof(*res)); return (res->status = RB3_SNAP_SHORT); }
	win = RB3_MALLOC(uint8_t, whi - wlo);
	if (ca && ca->buf && ca->lo <= wlo && whi <= ca->hi) {
		memcpy(win, ca->buf + (wlo - ca->lo), whi - wlo);
	} else {
		flo = sign > 0? wlo : L - whi, fhi = sign > 0? whi : L - wlo;
		if ((r = rb3_isa_extract(&c->f, c->isa, tsid * 2, flo, fhi, win)) < 0) { free(win); memset(res, 0, sizeof(*res)); return (res->status = RB3_SNAP_SHORT); }
		q->steps += r;
		if (sign < 0) rb3_revcomp6(whi - wlo, win);
	}
	if (reversed) {
		int64_t j, m = whi - wlo;
		for (j = 0; j < m>>1; ++j) { uint8_t z = win[j]; win[j] = win[m-1-j]; win[m-1-j] = z; }
	}
	rb3_snap_align(qs, nq, win, (int32_t)(whi - wlo), so, res);
	if (reversed && res->t_first >= 0) { // map back: reversed first M = last M of the forward window
		int32_t m = (int32_t)(whi - wlo), tf = res->t_first, qf = res->q_first;
		res->t_first = m - 1 - res->t_last, res->t_last = m - 1 - tf;
		res->q_first = nq - 1 - res->q_last, res->q_last = nq - 1 - qf;
	}
	free(win);
	return res->status;
}

// The flank itself is deleted in the founder (a deletion longer than the flank): step
// the flank away from the boundary x (t = F, 2F, 4F, ... <= radius2, inside the
// reference contig) until it aligns, then align the whole reference stretch from x to
// that flank, so the last present base before x (left side, dir = 1) / the first
// present base from x on (right side, dir = 0) is found exactly. Oriented coordinates;
// p = the estimated homolog of the flank base next to x. Depends only on (x, dir, p),
// never on the query interval, so both pieces sharing a boundary get the same answer.
static int snap_walk(const ext_ctx_t *c, ext_query_t *q, int64_t rsid, int64_t x, int dir, int sign, int64_t tsid, int64_t p,
					 int64_t *u_out, const snap_cache_t *ca)
{
	const ext_opt_t *o = c->opt;
	int64_t t, tmax, F = o->snap_flank, Lr = c->sid->len[rsid], wlo, y = -1, r, ra, rb;
	rb3_snap_res_t res;
	uint8_t *rs, *qs;
	int status = RB3_SNAP_LOWID;
	for (tmax = 0, t = F; t <= o->snap_r2 && (dir? x - t - F >= 0 : x + t + F <= Lr); t <<= 1) tmax = t;
	if (tmax == 0) return status;
	// the reference stretch [x .. x +- tmax+F), extracted once
	ra = dir? x - tmax - F : x, rb = dir? x : x + tmax + F;
	rs = RB3_MALLOC(uint8_t, rb - ra);
	if ((r = rb3_isa_extract(&c->f, c->isa, rsid * 2, ra, rb, rs)) < 0) { free(rs); return RB3_SNAP_SHORT; }
	q->steps += r;
	for (t = F; t <= tmax; t <<= 1) {
		int64_t a = dir? x - t - F : x + t;
		// the stepped flank's homolog: near p if [x, x+-t) is deleted, up to t further otherwise
		if (!dir) status = snap_window(c, q, tsid, sign, p - o->snap_r2, p + t + F + o->snap_r2, rs + (a - ra), (int32_t)F, &o->sopt, &res, &wlo, 0, ca);
		else status = snap_window(c, q, tsid, sign, p - t - F - o->snap_r2, p + o->snap_r2 + 1, rs + (a - ra), (int32_t)F, &o->sopt, &res, &wlo, 0, ca);
		if (status == RB3_SNAP_OK) { y = wlo + (dir? res.t_last : res.t_first); break; }
		if (status == RB3_SNAP_AMBIG) break;
	}
	if (y < 0) { free(rs); return status; }
	// whole stretch [x .. stepped flank]: every base is aligned or deleted. Most of it
	// is deleted, so the deleted run next to x is clipped at a flat cost (as one gap
	// event) instead of letting random bases align across the founder window; the left
	// side is aligned reversed so that run is always the query prefix.
	{
		rb3_snap_opt_t so = o->sopt;
		int64_t a = dir? x - t - F : x, n = t + F, wlo2, whi2, j;
		so.min_frac = 0, so.min_aligned = (int32_t)(F / 2), so.clip5 = so.gap_open + so.gap_ext;
		qs = RB3_MALLOC(uint8_t, n);
		memcpy(qs, rs + (a - ra), n);
		if (!dir) wlo2 = y - t - o->snap_r1, whi2 = y + F + o->snap_r1;
		else wlo2 = y - F - o->snap_r1, whi2 = y + t + o->snap_r1 + 1;
		if (dir) for (j = 0; j < n>>1; ++j) { uint8_t z = qs[j]; qs[j] = qs[n-1-j]; qs[n-1-j] = z; }
		status = snap_window(c, q, tsid, sign, wlo2, whi2, qs, (int32_t)n, &so, &res, &wlo, dir, ca);
		free(qs);
		if (status == RB3_SNAP_OK) *u_out = !dir? wlo + res.t_first : wlo + res.t_last;
	}
	free(rs);
	return status;
}

// Snap one side of boundary x: dir = 1 aligns the left flank [x-F, x) and returns the
// homolog of the last present base before x; dir = 0 aligns the right flank [x, x+F)
// and returns the homolog of the first present base from x on (*u, oriented). p = the
// estimated homolog of the flank base next to x. Returns SNAP_SNAP/SNAP_WIDE/SNAP_WALK
// or the failure label. *edge = 1 if the flank base next to x (x-1 or x) is itself
// present (aligned), i.e. there is no deletion on that side of x.
// walk = 0: the +-R1 / +-R2 windows only; walk = 1: only the walk (both windows
// having given low identity).
static int snap_side(const ext_ctx_t *c, ext_query_t *q, int64_t rsid, int64_t x, int dir, int sign, int64_t tsid, int64_t p, int64_t *u, int *edge, int walk,
					 snap_cache_t *ca)
{
	const ext_opt_t *o = c->opt;
	int64_t L = c->sid->len[tsid], Lr = c->sid->len[rsid], a, b, radii[2], r, wlo;
	int32_t F = o->snap_flank, nq, k;
	int lab = SNAP_SHORT;
	uint8_t *fl;
	rb3_snap_res_t res;
	*u = -1, *edge = 0;
	a = dir? (x - F > 0? x - F : 0) : x;
	b = dir? x : (x + F < Lr? x + F : Lr);
	nq = (int32_t)(b - a);
	if (nq < o->sopt.min_aligned) return SNAP_SHORT;
	fl = RB3_MALLOC(uint8_t, nq);
	if ((r = rb3_isa_extract(&c->f, c->isa, rsid * 2, a, b, fl)) < 0) { free(fl); return SNAP_SHORT; }
	q->steps += r;
	if (p < 0) p = 0;
	if (p > L - 1) p = L - 1;
	radii[0] = o->snap_r1, radii[1] = o->snap_r2;
	if (walk) { // the flank itself may be deleted
		lab = SNAP_LOWID;
		if (o->snap_r2 > 0) {
			int st;
			if (!dir) snap_cache_fill(c, q, tsid, sign, p - radii[1] - F, p + 2 * radii[1] + 2 * F + 1, ca);
			else snap_cache_fill(c, q, tsid, sign, p - 2 * radii[1] - 2 * F, p + radii[1] + F + 1, ca);
			st = snap_walk(c, q, rsid, x, dir, sign, tsid, p, u, ca); // *edge stays 0: the next base is deleted
			if (st == RB3_SNAP_OK) lab = SNAP_WALK;
			else if (st == RB3_SNAP_AMBIG) lab = SNAP_AMBIG;
		}
		free(fl);
		return lab;
	}
	for (k = 0; k < 2; ++k) {
		int st;
		if (k == 1 && radii[1] <= radii[0]) break;
		if (k == 1) { // one extraction covering the wide window and every walk step (t <= R2)
			if (!dir) snap_cache_fill(c, q, tsid, sign, p - radii[1] - F, p + 2 * radii[1] + 2 * F + 1, ca);
			else snap_cache_fill(c, q, tsid, sign, p - 2 * radii[1] - 2 * F, p + radii[1] + F + 1, ca);
		}
		st = snap_window(c, q, tsid, sign, p - radii[k] - F, p + radii[k] + F + 1, fl, nq, &o->sopt, &res, &wlo, 0, ca);
		if (st == RB3_SNAP_OK) {
			*u = wlo + (dir? res.t_last : res.t_first);
			*edge = dir? res.q_last == nq - 1 : res.q_first == 0;
			lab = k == 0? SNAP_SNAP : SNAP_WIDE;
			break;
		}
		lab = st == RB3_SNAP_AMBIG? SNAP_AMBIG : st == RB3_SNAP_LOWID? SNAP_LOWID : SNAP_SHORT;
		if (st != RB3_SNAP_LOWID) break; // a wider window cannot resolve an ambiguity
	}
	free(fl);
	return lab;
}

// label of a side placed by the flank on the OTHER side of the boundary than its own
// (a start placed by its left flank, an end by its right flank)
static int snap_cross(int lab)
{
	return lab == SNAP_SNAP? SNAP_R_SNAP : lab == SNAP_WIDE? SNAP_R_WIDE : SNAP_R_WALK;
}

// The founder fencepost (oriented) of reference boundary x, between bases x-1 and x.
// An END at x is one past the homolog of the last reference base before x present in
// the founder (left flank [x-F, x)); a START at x is the homolog of the first present
// base from x on (right flank [x, x+F)). Both are derived from the same placements of
// the two flanks of x, in this order, so the pieces [a,x) and [x,b) always agree:
//  1. both flanks placed within +-R1/R2, x-1 and x both present (no deletion at x) and
//     the right flank not overlapping the left: end = after the left, start = the
//     right; any gap is founder sequence inserted exactly at x (truth convention:
//     in neither piece);
//  2. otherwise one flank placed (left preferred): end = start = that junction, so a
//     deletion containing x maps both pieces to the same founder base, and founder-
//     only (replacement) bases at x go to the right piece;
//  3. neither placed: walk the left flank, then the right flank, away from x (the
//     flank itself deleted); end = start = that junction;
//  4. otherwise the lift projection (u_raw).
// Depends only on x (never on the interval). Labels: snap/snap_wide/snap_walk = the
// side's own flank (start: right, end: left); snap_x* = the flank across the boundary.
// Joint placement of boundary x: align the reference [x-F, x+F) -- both flanks as one
// query -- around the projection, so a deletion or insertion at x is an INTERNAL gap
// of an affine alignment and is placed by the whole context, rather than a terminal
// gap of a one-sided flank (where a flank's deleted tail is cheaper to align as a
// mismatch plus chance matches than to leave as a gap). Returns 1 on an accepted
// alignment with *tl = homolog of the last present base before x (-1 if none in the
// query), *tr = homolog of the first present base from x on (-1 if none), *both = x-1
// and x are both present (no deletion at x); oriented coordinates. *lab = SNAP_SNAP
// (+-R1) / SNAP_WIDE (+-R2) on success, else the failure label.
static int snap_joint(const ext_ctx_t *c, ext_query_t *q, int64_t x, int sign, int64_t tsid, int64_t p,
					  int64_t *tl, int64_t *tr, int *both, int *lab, snap_cache_t *ca)
{
	const ext_opt_t *o = c->opt;
	int64_t L = c->sid->len[tsid], Lr = c->sid->len[q->rsid], a, b, r, wlo, radii[2];
	int32_t F = o->snap_flank, nq, k;
	rb3_snap_opt_t so = o->sopt;
	rb3_snap_res_t res;
	uint8_t *fl;
	int ok = 0;
	*tl = *tr = -1, *both = 0, *lab = SNAP_SHORT;
	a = x - F > 0? x - F : 0, b = x + F < Lr? x + F : Lr;
	nq = (int32_t)(b - a);
	if (nq < so.min_aligned) return 0;
	fl = RB3_MALLOC(uint8_t, nq);
	if ((r = rb3_isa_extract(&c->f, c->isa, q->rsid * 2, a, b, fl)) < 0) { free(fl); return 0; }
	q->steps += r;
	so.split = (int32_t)(x - a);
	if (p < 0) p = 0;
	if (p > L - 1) p = L - 1;
	radii[0] = o->snap_r1, radii[1] = o->snap_r2;
	for (k = 0; k < 2; ++k) {
		int st;
		if (k == 1 && radii[1] <= radii[0]) break;
		if (k == 1) snap_cache_fill(c, q, tsid, sign, p - 2 * radii[1] - 2 * F, p + radii[1] + F + 1, ca);
		st = snap_window(c, q, tsid, sign, p - radii[k] - F, p + radii[k] + F + 1, fl, nq, &so, &res, &wlo, 0, ca);
		if (st == RB3_SNAP_OK && (res.t_left >= 0 || res.t_right >= 0)) {
			if (res.t_left >= 0) *tl = wlo + res.t_left;
			if (res.t_right >= 0) *tr = wlo + res.t_right;
			*both = res.q_left == so.split - 1 && res.q_right == so.split;
			*lab = k == 0? SNAP_SNAP : SNAP_WIDE;
			ok = 1;
			break;
		}
		*lab = st == RB3_SNAP_AMBIG? SNAP_AMBIG : st == RB3_SNAP_LOWID? SNAP_LOWID : SNAP_SHORT;
		if (st != RB3_SNAP_LOWID) break;
	}
	free(fl);
	return ok;
}

static int snap_boundary(const ext_ctx_t *c, ext_query_t *q, int64_t x, int is_start, int sign, int64_t tsid, int64_t u_raw, int64_t *U)
{
	int64_t uL = -1, uR = -1;
	int labL = SNAP_SHORT, labR = SNAP_SHORT, eL = 0, eR = 0, lab;
	snap_cache_t caL = { 0, 0, 0 }, caR = { 0, 0, 0 }; // kept for the walks
	*U = u_raw;
	{ // 0. both flanks as one query (see snap_joint)
		int64_t tl, tr;
		int both, jl;
		// A deletion at x (x-1 or x not present): the junction, for both pieces. With x-1 and
		// x both present, any founder bases between them are an insertion at x, which the
		// two one-sided flanks below place exactly (their free target ends leave it out),
		// whereas in the joint alignment an insertion's I gap can slide by equal bases.
		if (snap_joint(c, q, x, sign, tsid, u_raw, &tl, &tr, &both, &jl, &caL) && !both) {
			*U = tl >= 0? tl + 1 : tr;
			free(caL.buf);
			return is_start && tl >= 0? snap_cross(jl) : !is_start && tl < 0? snap_cross(jl) : jl;
		}
	}
	if (x > 0) labL = snap_side(c, q, q->rsid, x, 1, sign, tsid, u_raw - 1, &uL, &eL, 0, &caL);
	if (is_start || uL < 0) labR = snap_side(c, q, q->rsid, x, 0, sign, tsid, u_raw, &uR, &eR, 0, &caR);
	if (uL >= 0) {
		if (is_start && uR >= 0 && eL && eR && uR >= uL + 1) *U = uR, lab = labR;
		else *U = uL + 1, lab = is_start? snap_cross(labL) : labL;
	} else if (uR >= 0) {
		*U = uR, lab = is_start? labR : snap_cross(labR);
	} else {
		lab = labL == SNAP_SHORT? labR : labL;
		if (labL == SNAP_LOWID) {
			labL = snap_side(c, q, q->rsid, x, 1, sign, tsid, u_raw - 1, &uL, &eL, 1, &caL);
			if (uL >= 0) *U = uL + 1, lab = is_start? snap_cross(labL) : labL;
		}
		if (uL < 0 && labR == SNAP_LOWID) {
			labR = snap_side(c, q, q->rsid, x, 0, sign, tsid, u_raw, &uR, &eR, 1, &caR);
			if (uR >= 0) *U = uR, lab = is_start? labR : snap_cross(labR);
		}
		if (uL < 0 && uR < 0) lab = labL == SNAP_SHORT? labR : labL;
	}
	free(caL.buf); free(caR.buf);
	return lab;
}

static void q_project(const ext_ctx_t *c, ext_query_t *q, void *km)
{
	const ext_opt_t *o = c->opt;
	int32_t g = rb3_gtab_find(c->gtab, q->asm_name);
	int64_t cs[2], cp[2], len, U[2], x[2];
	int sg[2], i;
	if (g == c->gtab->sid2g[q->rsid]) { // the reference itself: identity
		q->tsid = q->rsid, q->st = q->rst, q->en = q->ren, q->strand = 1, q->meth[0] = q->meth[1] = 5;
		q->raw_ok[0] = q->raw_ok[1] = 1, q->raw_tsid[0] = q->raw_tsid[1] = q->rsid, q->raw_sign[0] = q->raw_sign[1] = 1;
		q->raw_pos[0] = q->rst, q->raw_pos[1] = q->ren;
		if (o->snap) q->snap[0] = q->snap[1] = SNAP_IDENTITY;
		return;
	}
	// each boundary x is projected through the base before it (x-1; base 0 for x = 0),
	// so the end of [a,m) and the start of [m,b) are the same computation
	x[0] = q->rst, x[1] = q->ren;
	for (i = 0; i < 2; ++i) { // both ends, so that a failed row still reports the raw end that projected
		int64_t rp = x[i] > 0? x[i] - 1 : 0;
		if (rb3_lift_project_inv(c->ridx, km, g, q->rsid, rp, o->lift_win, o->lift_mad, o->min_support, o->lift_k, &cs[i], &cp[i], &sg[i], &q->meth[i])) {
			int64_t Lc = c->sid->len[cs[i]], uc = sg[i] > 0? cp[i] : Lc - 1 - cp[i], u = x[i] > 0? uc + 1 : uc;
			q->raw_ok[i] = 1, q->raw_tsid[i] = cs[i], q->raw_sign[i] = sg[i];
			q->raw_pos[i] = sg[i] > 0? u : Lc - u; // forward fencepost
		}
	}
	for (i = 0; i < 2; ++i)
		if (!q->raw_ok[i]) {
			q_fail(q, ST_PROJ_FAIL, i == 0? "no confident projection for the start" : "no confident projection for the end");
			return;
		}
	if (cs[0] != cs[1]) { q_fail(q, ST_PROJ_FAIL, "start and end project to different contigs"); q->tsid = -1; return; }
	if (sg[0] != sg[1]) { q_fail(q, ST_PROJ_FAIL, "start and end project with different orientations"); return; }
	q->tsid = cs[0], q->strand = sg[0];
	len = c->sid->len[q->tsid];
	for (i = 0; i < 2; ++i) { // oriented fenceposts (u = L-1-cp on '-')
		int64_t uc = q->strand > 0? cp[i] : len - 1 - cp[i], uraw = x[i] > 0? uc + 1 : uc;
		U[i] = uraw;
		if (o->snap) {
			q->snap[i] = snap_boundary(c, q, x[i], i == 0, q->strand, q->tsid, uraw, &U[i]);
			q->shift[i] = q->strand > 0? U[i] - uraw : uraw - U[i]; // in forward coordinates
		}
	}
	if (q->strand > 0) q->st = U[0], q->en = U[1];
	else q->st = len - U[1], q->en = len - U[0];
	if (q->en == q->st) { q_fail(q, ST_PROJ_FAIL, "interval deleted in the founder (empty projection)"); return; }
	if (q->en < q->st) { q_fail(q, ST_PROJ_FAIL, "projected end precedes projected start"); return; }
	if (q->st < 0) q->st = 0, q->clamped = 1;
	if (q->st >= len) { q_fail(q, ST_PROJ_FAIL, "projected start is past the contig end"); return; }
	if (q->en > len) q->en = len, q->clamped = 1;
	// sanity: a span far longer than the reference interval is a coarse (wide-median)
	// projection gone wrong, not a real insertion (same rule as the Python builder)
	if (o->max_span_excess > 0 && o->max_span_ratio > 0) {
		int64_t ref_len = q->ren - q->rst, prj_len = q->en - q->st;
		if (prj_len - ref_len > o->max_span_excess && (double)prj_len > o->max_span_ratio * (double)ref_len) {
			snprintf(q->msgbuf, sizeof(q->msgbuf), "projected span implausibly long (%ld bp for %ld bp)", (long)prj_len, (long)ref_len);
			q_fail(q, ST_PROJ_FAIL, q->msgbuf);
			return;
		}
	}
}

/**************
 * extraction *
 **************/

static void q_extract(const ext_ctx_t *c, ext_query_t *q)
{
	int64_t n = q->en - q->st, i, r;
	uint8_t *buf = RB3_MALLOC(uint8_t, n + 1);
	r = rb3_isa_extract(&c->f, c->isa, q->tsid * 2, q->st, q->en, buf);
	if (r < 0) {
		free(buf);
		q_fail(q, ST_INPUT_ERR, r == RB3_ISA_ERR_WALK? "LF walk hit a sentinel: the .isa does not match this index" : "extraction failed (range/strand not covered by the .isa)");
		return;
	}
	q->steps += r;
	if (q->strand < 0) rb3_revcomp6(n, buf);
	for (i = 0; i < n; ++i) buf[i] = "$ACGTN"[buf[i] < 6? buf[i] : 5];
	buf[n] = 0;
	q->seq = (char*)buf;
}

static void ext_worker(void *data, long i, int tid)
{
	ext_ctx_t *c = (ext_ctx_t*)data;
	ext_query_t *q = &c->q[i];
	if (q->status != ST_OK) return;
	if (c->lift) q_project(c, q, c->km[tid]);
	if (q->status != ST_OK || c->opt->project_only) return;
	q_extract(c, q);
}

/**********
 * output *
 **********/

static void print_proj_header(void)
{
	printf("qid\tref_contig\tref_start\tref_end\tassembly\tcontig\tstart\tend\tstrand\tstatus\tmethod\tsnap\tsnap_shift\traw_start\traw_end\n");
}

// raw_start / raw_end: the unsnapped lift projection of the reference boundary
// ref_start / ref_end, as CONTIG:POS:STRAND, '.' if that end did not project. POS is a
// forward-strand fencepost on CONTIG: boundary x is projected through reference base
// x-1 (base 0 for x = 0); on '+' POS = homolog(x-1)+1 (homolog(0) for x = 0), on '-'
// POS = homolog(x-1) (homolog(0)+1 for x = 0). So with both ends on one contig and
// strand the raw interval is [raw_start, raw_end) on '+' and [raw_end, raw_start) on
// '-' -- exactly what --no-snap reports when the row succeeds.
static void print_raw(const ext_ctx_t *c, const ext_query_t *q, int i)
{
	if (q->raw_ok[i]) printf("%s:%ld:%c", c->sid->name[q->raw_tsid[i]], (long)q->raw_pos[i], q->raw_sign[i] > 0? '+' : '-');
	else putchar('.');
}

static void print_proj_row(const ext_ctx_t *c, const ext_query_t *q)
{
	const char *ref = q->rsid >= 0? c->sid->name[q->rsid] : ".";
	printf("%s\t%s\t", q->qid, ref);
	if (q->rsid >= 0 && q->status != ST_INPUT_ERR) printf("%ld\t%ld\t", (long)q->rst, (long)q->ren);
	else printf(".\t.\t");
	printf("%s\t", q->asm_name && *q->asm_name? q->asm_name : ".");
	if (q->status == ST_OK)
		printf("%s\t%ld\t%ld\t%c\t%s\t%s,%s\t%s,%s\t%ld,%ld\t", c->sid->name[q->tsid], (long)q->st, (long)q->en, q->strand > 0? '+' : '-',
			   q->clamped? "ok_clamped" : "ok", meth_name(q->meth[0]), meth_name(q->meth[1]),
			   snap_name(q->snap[0]), snap_name(q->snap[1]), (long)q->shift[0], (long)q->shift[1]);
	else
		printf(".\t.\t.\t.\t%s: %s\t%s,%s\t%s,%s\t%ld,%ld\t", q->status == ST_INPUT_ERR? "input_error" : "projection_failed", q->msg,
			   meth_name(q->meth[0]), meth_name(q->meth[1]), snap_name(q->snap[0]), snap_name(q->snap[1]), (long)q->shift[0], (long)q->shift[1]);
	print_raw(c, q, 0); putchar('\t'); print_raw(c, q, 1); putchar('\n');
}

static void print_fasta(const ext_ctx_t *c, const ext_query_t *q)
{
	int64_t n = q->en - q->st, i;
	printf(">%s", q->qid);
	if (c->lift)
		printf(" %s:%ld-%ld:%c ref=%s:%ld-%ld", c->sid->name[q->tsid], (long)q->st, (long)q->en, q->strand > 0? '+' : '-',
			   c->sid->name[q->rsid], (long)q->rst, (long)q->ren);
	if (c->lift) printf(" snap=%s,%s", snap_name(q->snap[0]), snap_name(q->snap[1]));
	if (q->clamped) printf(" clamped");
	putchar('\n');
	if (c->opt->line_width <= 0) { fwrite(q->seq, 1, n, stdout); putchar('\n'); }
	else for (i = 0; i < n; i += c->opt->line_width) {
		int64_t l = n - i < c->opt->line_width? n - i : c->opt->line_width;
		fwrite(q->seq + i, 1, l, stdout); putchar('\n');
	}
}

static void q_free(ext_query_t *q)
{
	free(q->qid); free(q->region); free(q->asm_name); free(q->seq);
}

/********
 * main *
 ********/

static int usage_extract(FILE *fp, const ext_opt_t *o)
{
	fprintf(fp, "Usage: ropebwt3 extract [options] <idx.fmd> [CONTIG:START-END[:+|-] ...]\n");
	fprintf(fp, "Extracts sequence by coordinate. Coordinates are 0-BASED, HALF-OPEN [START,END) on the\n");
	fprintf(fp, "forward strand of the named sequence (names as in <idx.fmd>.len.gz, e.g. Oh43_chr1);\n");
	fprintf(fp, "':-' returns the reverse complement of the same interval. Output is FASTA, '>qid' (the\n");
	fprintf(fp, "region string on the command line), bases in upper case ACGTN.\n");
	fprintf(fp, "Queries:\n");
	fprintf(fp, "  -f FILE        TSV queries ('-' = stdin): qid contig start end [strand]; with --lift:\n");
	fprintf(fp, "                 qid ref_contig ref_start ref_end [assembly] (extra columns ignored);\n");
	fprintf(fp, "                 one header line allowed, '#' lines skipped\n");
	fprintf(fp, "  --clamp        clip intervals to [0,contig length) instead of rejecting them\n");
	fprintf(fp, "Index:\n");
	fprintf(fp, "  -I FILE        entry points from `ropebwt3 isa` [<idx.fmd>.isa; if absent, derived\n");
	fprintf(fp, "                 in memory from <idx.fmd>.ssa]\n");
	fprintf(fp, "  -M             load the FMD via mmap\n");
	fprintf(fp, "  -t INT         number of threads [%d]\n", o->n_threads);
	fprintf(fp, "Projection (reference interval -> founder assembly, the refmap liftover):\n");
	fprintf(fp, "  --lift=FILE    `ropebwt3 lift` map; queries are then REFERENCE intervals\n");
	fprintf(fp, "  --ref-prefix=STR  reference contig name prefix, e.g. B73 (validates ref contigs)\n");
	fprintf(fp, "  --to=ASM       target assembly (name prefix before '_') for region queries / TSV rows\n");
	fprintf(fp, "                 without an assembly column\n");
	fprintf(fp, "  --project-only print the projected interval as TSV, do not extract\n");
	fprintf(fp, "  --lift-win=NUM, --lift-mad=NUM  projection window / max residual MAD [%ld, %ld]\n", (long)o->lift_win, (long)o->lift_mad);
	fprintf(fp, "  --lift-local-win=NUM  nearest collinear flanking anchors within NUM bp first; 0 = off,\n");
	fprintf(fp, "                 the wide median only (refmap's built-in default) [%ld]\n", (long)o->local_win);
	fprintf(fp, "  --lift-local-tol=NUM  collinearity tolerance for --lift-local-win [%ld]\n", (long)o->local_tol);
	fprintf(fp, "  --lift-min-support=INT  min anchors for the median fit (refmap uses 4) [%d]\n", o->min_support);
	fprintf(fp, "  --lift-k=INT   the -k the .lift was built with (not stored in the file); only\n");
	fprintf(fp, "                 shifts projections onto inverted ('-') loci by k-1 [%d]\n", o->lift_k);
	fprintf(fp, "  --max-span-excess=NUM, --max-span-ratio=FLOAT  fail a projection (after snapping) whose\n");
	fprintf(fp, "                 span exceeds the reference interval by > NUM bp AND is > FLOAT x as long;\n");
	fprintf(fp, "                 0 for either turns the check off [%ld, %.1f]\n", (long)o->max_span_excess, o->max_span_ratio);
	fprintf(fp, "Snapping projected ends to the homologous base (with --lift):\n");
	fprintf(fp, "  --no-snap      keep the lift projection as is [snap on]\n");
	fprintf(fp, "  --snap-flank=INT  reference flank aligned at each end [%d]\n", o->snap_flank);
	fprintf(fp, "  --snap-radius=NUM, --snap-radius2=NUM  founder window +- around the projection; the\n");
	fprintf(fp, "                 second (0 = off) is tried when the first gives low identity [%ld, %ld]\n", (long)o->snap_r1, (long)o->snap_r2);
	fprintf(fp, "  --snap-min-id=FLOAT  min gap-compressed identity of the flank alignment [%.2f]\n", o->sopt.min_id);
	fprintf(fp, "  --snap-margin=INT  ambiguous if another placement scores within INT of the best [%d]\n", o->sopt.margin);
	fprintf(fp, "Output:\n");
	fprintf(fp, "  -w INT         FASTA line width, 0 = one line [%d]\n", o->line_width);
	fprintf(fp, "Exit status 1 if any query had an input error (unparsable, unknown contig/assembly,\n");
	fprintf(fp, "out of range, empty); projection failures are reported but do not set it.\n");
	return fp == stdout? 0 : 1;
}

int main_extract(int argc, char *argv[])
{
	ext_opt_t opt;
	ext_ctx_t ctx;
	ketopt_t o = KETOPT_INIT;
	int c, ret = 0, i;
	char *buf;
	ext_query_t *qs = 0;
	int64_t n_q = 0, m_q = 0, k, n_input_err = 0, n_proj_fail = 0, tot_steps = 0, tot_bases = 0;
	double t_load, t0;
	static ko_longopt_t lo[] = {
		{ "lift",             ko_required_argument, 401 },
		{ "ref-prefix",       ko_required_argument, 402 },
		{ "to",               ko_required_argument, 403 },
		{ "project-only",     ko_no_argument,       404 },
		{ "lift-win",         ko_required_argument, 405 },
		{ "lift-mad",         ko_required_argument, 406 },
		{ "lift-local-win",   ko_required_argument, 407 },
		{ "lift-local-tol",   ko_required_argument, 408 },
		{ "lift-min-support", ko_required_argument, 409 },
		{ "clamp",            ko_no_argument,       410 },
		{ "lift-k",           ko_required_argument, 411 },
		{ "snap",             ko_no_argument,       412 },
		{ "no-snap",          ko_no_argument,       413 },
		{ "snap-flank",       ko_required_argument, 414 },
		{ "snap-radius",      ko_required_argument, 415 },
		{ "snap-radius2",     ko_required_argument, 416 },
		{ "snap-min-id",      ko_required_argument, 417 },
		{ "snap-margin",      ko_required_argument, 418 },
		{ "max-span-excess",  ko_required_argument, 419 },
		{ "max-span-ratio",   ko_required_argument, 420 },
		{ 0, 0, 0 }
	};

	memset(&opt, 0, sizeof(opt));
	opt.n_threads = 1, opt.line_width = 0, opt.min_support = 4, opt.lift_k = 100;
	opt.max_span_excess = 50000, opt.max_span_ratio = 3.0;
	opt.snap = 1, opt.snap_flank = 64, opt.snap_r1 = 200, opt.snap_r2 = 6000;
	rb3_snap_opt_init(&opt.sopt);
	// refmap's projection defaults, except the local bracket is ON (20 kb, the adopted
	// refmap setting): without it the +-500 kb median offset is kb-scale off for an
	// interval (see README); --lift-local-win=0 gives refmap's built-in default exactly
	opt.lift_win = 500000, opt.lift_mad = 200000, opt.local_win = 20000, opt.local_tol = 1000;
	while ((c = ketopt(&o, argc, argv, 1, "f:I:Mt:w:", lo)) >= 0) {
		if (c == 'f') opt.tsv_fn = o.arg;
		else if (c == 'I') opt.isa_fn = o.arg;
		else if (c == 'M') opt.use_mmap = 1;
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'w') opt.line_width = atoi(o.arg);
		else if (c == 401) opt.lift_fn = o.arg;
		else if (c == 402) opt.ref_prefix = o.arg;
		else if (c == 403) opt.to_asm = o.arg;
		else if (c == 404) opt.project_only = 1;
		else if (c == 405) opt.lift_win = rb3_parse_num(o.arg);
		else if (c == 406) opt.lift_mad = rb3_parse_num(o.arg);
		else if (c == 407) opt.local_win = rb3_parse_num(o.arg);
		else if (c == 408) opt.local_tol = rb3_parse_num(o.arg);
		else if (c == 409) opt.min_support = atoi(o.arg);
		else if (c == 410) opt.clamp = 1;
		else if (c == 411) opt.lift_k = atoi(o.arg);
		else if (c == 412) opt.snap = 1;
		else if (c == 413) opt.snap = 0;
		else if (c == 414) opt.snap_flank = atoi(o.arg);
		else if (c == 415) opt.snap_r1 = rb3_parse_num(o.arg);
		else if (c == 416) opt.snap_r2 = rb3_parse_num(o.arg);
		else if (c == 417) opt.sopt.min_id = atof(o.arg);
		else if (c == 418) opt.sopt.margin = atoi(o.arg);
		else if (c == 419) opt.max_span_excess = rb3_parse_num(o.arg);
		else if (c == 420) opt.max_span_ratio = atof(o.arg);
		else if (c == '?' || c == ':') {
			fprintf(stderr, "ERROR: unknown or incomplete option '%s'\n", argv[o.i - 1]);
			return 1;
		}
	}
	if (argc == o.ind) return usage_extract(stdout, &opt);
	if (opt.tsv_fn == 0 && argc - o.ind < 2) {
		fprintf(stderr, "ERROR: no queries (give regions or -f FILE)\n");
		return 1;
	}
	if (opt.project_only && opt.lift_fn == 0) { fprintf(stderr, "ERROR: --project-only needs --lift\n"); return 1; }
	if (opt.lift_fn == 0 && (opt.to_asm || opt.ref_prefix)) { fprintf(stderr, "ERROR: --to/--ref-prefix need --lift\n"); return 1; }
	if (opt.lift_fn && opt.ref_prefix == 0) { fprintf(stderr, "ERROR: --lift needs --ref-prefix\n"); return 1; }
	if (opt.n_threads < 1) opt.n_threads = 1;
	if (opt.lift_fn == 0) opt.snap = 0;
	if (opt.max_span_excess < 0 || opt.max_span_ratio < 0) {
		fprintf(stderr, "ERROR: --max-span-excess/--max-span-ratio must be >= 0 (0 = check off)\n");
		return 1;
	}
	if (opt.snap && (opt.snap_flank < opt.sopt.min_aligned || opt.snap_r1 < 0 || opt.snap_r2 < 0)) {
		fprintf(stderr, "ERROR: --snap-flank must be >= %d and the radii >= 0\n", opt.sopt.min_aligned);
		return 1;
	}

	memset(&ctx, 0, sizeof(ctx));
	ctx.opt = &opt;
	buf = RB3_CALLOC(char, strlen(argv[o.ind]) + 8);
	t0 = rb3_realtime();
	// sequence names/lengths
	ctx.sid = rb3_sid_read(strcat(strcpy(buf, argv[o.ind]), ".len.gz"));
	if (ctx.sid == 0 || ctx.sid->n_seq == 0) {
		fprintf(stderr, "ERROR: failed to read sequence names/lengths from '%s'\n", buf);
		ret = 1; goto end_extract;
	}
	ctx.name_ord = name_index_build(ctx.sid);
	// the BWT (not needed for --project-only without snapping)
	if (!opt.project_only || opt.snap) {
		rb3_fmi_restore(&ctx.f, argv[o.ind], opt.use_mmap);
		if (ctx.f.e == 0 && ctx.f.r == 0) {
			fprintf(stderr, "ERROR: failed to load the FM-index '%s'\n", argv[o.ind]);
			ret = 1; goto end_extract;
		}
		if (ctx.f.acc[1] != ctx.sid->n_seq * 2) {
			fprintf(stderr, "ERROR: the index has %ld strings but .len.gz lists %ld sequences (need 2 strings each: build without -F/-R)\n",
					(long)ctx.f.acc[1], (long)ctx.sid->n_seq);
			ret = 1; goto end_extract;
		}
		if (rb3_verbose >= 3) fprintf(stderr, "[M::%s] loaded the FM-index in %.3f sec\n", __func__, rb3_realtime() - t0);
		// entry points
		t_load = rb3_realtime();
		if (opt.isa_fn) ctx.isa = rb3_isa_restore(opt.isa_fn);
		else {
			FILE *fp;
			strcat(strcpy(buf, argv[o.ind]), ".isa");
			if ((fp = fopen(buf, "rb")) != 0) {
				fclose(fp);
				ctx.isa = rb3_isa_restore(buf);
			} else {
				strcat(strcpy(buf, argv[o.ind]), ".ssa");
				if (rb3_verbose >= 3) fprintf(stderr, "[M::%s] no .isa; deriving entry points from '%s' in memory\n", __func__, buf);
				ctx.isa = rb3_isa_from_ssa_file(buf, ctx.sid->n_seq, ctx.sid->len, 12, RB3_ISA_FWD);
			}
		}
		if (ctx.isa == 0) {
			fprintf(stderr, "ERROR: failed to load/derive the entry points (.isa, or .ssa to derive it from)\n");
			ret = 1; goto end_extract;
		}
		if (rb3_isa_check_fmi(ctx.isa, &ctx.f) < 0) {
			fprintf(stderr, "ERROR: the entry points were not built from this index's .ssa\n");
			ret = 1; goto end_extract;
		}
		for (k = 0; k < ctx.sid->n_seq; ++k)
			if (rb3_isa_seq_len(ctx.isa, k * 2) != ctx.sid->len[k]) break;
		if (k < ctx.sid->n_seq) {
			fprintf(stderr, "ERROR: the entry points disagree with .len.gz on the length of '%s'\n", ctx.sid->name[k]);
			ret = 1; goto end_extract;
		}
		if (rb3_verbose >= 3) fprintf(stderr, "[M::%s] %ld entry points ready in %.3f sec\n", __func__, (long)ctx.isa->n, rb3_realtime() - t_load);
	}
	// the liftover
	if (opt.lift_fn) {
		t_load = rb3_realtime();
		rb3_lift_set_local(opt.local_win, opt.local_tol);
		ctx.lift = rb3_lift_restore(opt.lift_fn);
		if (ctx.lift == 0) { fprintf(stderr, "ERROR: failed to load liftover '%s'\n", opt.lift_fn); ret = 1; goto end_extract; }
		if (rb3_lift_n_seq(ctx.lift) != ctx.sid->n_seq) {
			fprintf(stderr, "ERROR: liftover covers %ld sequences, .len.gz lists %ld\n", (long)rb3_lift_n_seq(ctx.lift), (long)ctx.sid->n_seq);
			ret = 1; goto end_extract;
		}
		ctx.gtab = rb3_gtab_build(ctx.sid);
		ctx.ridx = rb3_lift_ridx_build(ctx.lift, ctx.gtab->sid2g, ctx.gtab->n_gamete);
		if (opt.to_asm && rb3_gtab_find(ctx.gtab, opt.to_asm) < 0) {
			fprintf(stderr, "ERROR: --to '%s' is not an assembly in the index\n", opt.to_asm);
			ret = 1; goto end_extract;
		}
		if (rb3_verbose >= 3) fprintf(stderr, "[M::%s] loaded + indexed the liftover in %.3f sec\n", __func__, rb3_realtime() - t_load);
	}
	ctx.km = RB3_CALLOC(void*, opt.n_threads);
	for (i = 0; i < opt.n_threads; ++i) ctx.km[i] = km_init();

	// read the queries (parsing is cheap; extraction below runs in bounded batches)
	for (i = o.ind + 1; i < argc; ++i) {
		RB3_GROW(ext_query_t, qs, n_q, m_q);
		q_from_region(&ctx, &qs[n_q++], argv[i]);
	}
	if (opt.tsv_fn) {
		FILE *fp = strcmp(opt.tsv_fn, "-")? fopen(opt.tsv_fn, "r") : stdin;
		char *line = 0;
		size_t cap = 0;
		ssize_t l;
		int first = 1;
		if (fp == 0) { fprintf(stderr, "ERROR: cannot open '%s'\n", opt.tsv_fn); ret = 1; goto end_extract; }
		while ((l = getline(&line, &cap, fp)) >= 0) {
			while (l > 0 && (line[l-1] == '\n' || line[l-1] == '\r')) line[--l] = 0;
			if (l == 0 || line[0] == '#') continue;
			RB3_GROW(ext_query_t, qs, n_q, m_q);
			if (q_from_tsv(&ctx, &qs[n_q], line, first)) q_free(&qs[n_q]);
			else ++n_q;
			first = 0;
		}
		free(line);
		if (fp != stdin) fclose(fp);
	}

	t0 = rb3_realtime();
	if (opt.project_only) print_proj_header();
	for (k = 0; k < n_q; k += EXT_BATCH) {
		int64_t j, nb = n_q - k < EXT_BATCH? n_q - k : EXT_BATCH;
		ctx.q = qs + k;
		kt_for(opt.n_threads, ext_worker, &ctx, nb);
		for (j = 0; j < nb; ++j) {
			ext_query_t *q = &ctx.q[j];
			if (q->status == ST_INPUT_ERR) {
				++n_input_err;
				fprintf(stderr, "[E::%s] query '%s': %s\n", __func__, q->region, q->msg);
			} else if (q->status == ST_PROJ_FAIL) {
				++n_proj_fail;
				fprintf(stderr, "[W::%s] query '%s': projection failed: %s\n", __func__, q->qid, q->msg);
			}
			if (opt.project_only) print_proj_row(&ctx, q);
			else if (q->status == ST_OK) {
				print_fasta(&ctx, q);
				tot_steps += q->steps, tot_bases += q->en - q->st;
			}
			q_free(q);
		}
	}
	fflush(stdout);
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] %ld queries (%ld input errors, %ld projection failures); %ld bases via %ld LF steps in %.3f sec\n",
				__func__, (long)n_q, (long)n_input_err, (long)n_proj_fail, (long)tot_bases, (long)tot_steps, rb3_realtime() - t0);
	if (n_input_err > 0) ret = 1;

end_extract:
	if (ctx.km) { for (i = 0; i < opt.n_threads; ++i) km_destroy(ctx.km[i]); free(ctx.km); }
	free(qs);
	free(buf);
	free(ctx.name_ord);
	rb3_isa_destroy(ctx.isa);
	rb3_lift_ridx_destroy(ctx.ridx);
	rb3_lift_destroy(ctx.lift);
	rb3_gtab_destroy(ctx.gtab);
	if (ctx.f.e || ctx.f.r) rb3_fmi_free(&ctx.f);
	rb3_sid_destroy(ctx.sid);
	return ret;
}
