#ifndef RB3_PAIRING_H
#define RB3_PAIRING_H

// Paired-end helpers shared by the alignment modes (refmap first; chain and mem/sw later).
// Everything here is pure (no index, no I/O) so it can be unit tested in isolation; the
// mode-specific glue in search.c converts its per-read results into rb3_pe_mate_t.

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Pair classes (one per pair; both mates' stdout rows carry it).
#define RB3_PE_NONE 0 // neither mate has a confident reference position: no event
#define RB3_PE_SOLO 1 // exactly one mate placed: one event at that mate
#define RB3_PE_CONC 2 // both placed, same contig, within max_dist: one event at the anchor mate
#define RB3_PE_RESC 3 // discordant, but one mate has an alternative position near the other: moved there, one event
#define RB3_PE_DISC 4 // discordant and not rescuable: each mate is its own single-end event

extern const char *rb3_pe_class_str[5];

// Insert-size distribution (fragment span, bp).
typedef struct {
	int32_t n;       // pairs the estimate is based on
	int ok;          // 1 if the estimate is usable
	double mean, sd;
	int64_t max_ins; // largest fragment span still considered concordant
} rb3_isize_t;

// One mate, as the pair resolver sees it.
typedef struct {
	int8_t rank;             // 2 = exact reference hit, 1 = projected (PLACED), 0 = no confident position
	int64_t rsid, pos;       // reference sequence id and position; ignored when rank == 0
	int32_t n_proj;          // alternative positions: per-assembly lift projections (PLACED) or all
	                         // reference copies (EXACT); may be 0
	const int64_t *proj_rsid, *proj_pos;
	const int32_t *proj_gid; // gamete (founder) of each projection; may be NULL
} rb3_pe_mate_t;

typedef struct {
	int8_t cls;              // RB3_PE_*
	int8_t anchor;           // mate (0/1) whose position the pair event uses; -1 for NONE/DISC
	int8_t rescued;          // mate (0/1) moved next to the anchor (RESC only); -1 otherwise
	int64_t resc_rsid, resc_pos; // new position of the rescued mate
	int32_t resc_support;    // how many of its projections agreed with the anchor
} rb3_pe_res_t;

// Length of a read name without a trailing "/1" or "/2".
int rb3_pe_name_len(const char *name);

// 1 if two names are mates (equal after dropping "/1" / "/2"); a NULL name never matches.
int rb3_pe_names_match(const char *a, const char *b);

// bwa-style estimate from n fragment spans (span[] is sorted in place): drop spans outside
// [p25 - 2*IQR, p75 + 2*IQR], then mean/sd of the rest; max_ins = mean + 4*sd, and at least
// p75 + 3*IQR. Returns out->ok = 1 iff at least min_n spans were given.
int rb3_pe_isize_estimate(int64_t *span, int32_t n, int32_t min_n, rb3_isize_t *out);

// Classify a pair. max_dist: largest |pos0 - pos1| on one contig that is still concordant.
// Discordant pairs are rescued by moving one mate to the median of its alternative positions
// within max_dist of the other mate, preferring to move the less confident mate. trust_exact:
// an unrescuable EXACT + PLACED pair moves the PLACED mate onto the EXACT one (support 0).
void rb3_pe_resolve(const rb3_pe_mate_t m[2], int64_t max_dist, int trust_exact, rb3_pe_res_t *res);

// Founders of mate b whose own projection lies on contig rsid within max_dist of pos: sorted,
// deduplicated, into out (size >= b->n_proj). A rescued mate keeps only these, so an assembly
// whose copy of the read sits elsewhere (a TE / paralog copy) does not count as sharing the
// fragment. Returns 0 when b has no projections or no proj_gid.
int32_t rb3_pe_near_gametes(const rb3_pe_mate_t *b, int64_t rsid, int64_t pos, int64_t max_dist, int32_t *out);

// Founder (gamete) set of a single pair event: sorted, deduplicated intersection of the
// anchor's set a[] and the mate's set b[]. Falls back to a[] (sorted, deduped) when b[] is
// empty or capped (an incomplete occurrence sample cannot veto founders), and to b[] when
// only a[] is capped; an empty intersection also falls back to a[]. out must hold na + nb
// entries. Returns the output size.
int32_t rb3_pe_event_gametes(const int32_t *a, int32_t na, int a_capped,
							 const int32_t *b, int32_t nb, int b_capped, int32_t *out);

#ifdef __cplusplus
}
#endif

#endif
