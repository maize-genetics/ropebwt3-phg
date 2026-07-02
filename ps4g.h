#ifndef RB3_PS4G_H
#define RB3_PS4G_H

#include <stdint.h>
#include "io.h" // for rb3_sid_t

#ifdef __cplusplus
extern "C" {
#endif

// gamete table: one gamete per sample (sequence-name prefix before the first '_'),
// sorted by name so indices are stable and comparable across runs.
typedef struct {
	int32_t n_gamete;
	char **name;    // name[g] = sample name, sorted ascending
	int32_t *sid2g; // sid2g[seqIdx] = gamete index of sequence seqIdx (seqIdx in [0,sid->n_seq))
} rb3_gtab_t;

rb3_gtab_t *rb3_gtab_build(const rb3_sid_t *sid);
void rb3_gtab_destroy(rb3_gtab_t *g);
int32_t rb3_gtab_find(const rb3_gtab_t *g, const char *name); // -1 if absent

// accumulator: one event per contributing read (EXACT or PLACED), recording the
// reference sequence + position it was binned at and the set of gametes it supports.
typedef struct {
	int64_t bin_size;
	int64_t *ev_ref_sid, *ev_bin; // parallel arrays, one entry per event
	int32_t *ev_goff, *ev_glen;   // offset/length into garena for this event's gamete list
	int64_t n_ev, m_ev;
	int32_t *garena;
	int64_t n_garena, m_garena;
} rb3_ps4g_acc_t;

rb3_ps4g_acc_t *rb3_ps4g_acc_init(int64_t bin_size);
void rb3_ps4g_acc_add(rb3_ps4g_acc_t *acc, int64_t ref_sid, int64_t pos, const int32_t *gametes, int32_t n_gametes);
void rb3_ps4g_acc_destroy(rb3_ps4g_acc_t *acc);

// BED-based diploid training labels: chrom, start, end, sampleA[, sampleB]
typedef struct {
	int64_t ref_sid; // resolved reference sequence index
	int64_t start, end;
	int32_t gA, gB; // gamete indices; gB == gA when the BED gave one label (homozygous)
} rb3_bed_region_t;

typedef struct {
	rb3_bed_region_t *r;
	int64_t n_r;
} rb3_bed_t;

rb3_bed_t *rb3_bed_read(const char *fn, const rb3_gtab_t *gtab, const rb3_sid_t *sid, const char *ref_prefix);
void rb3_bed_destroy(rb3_bed_t *b);

// write both outputs from one sorted pass over the accumulated events; either
// ps4g_fn or npy_fn may be NULL to skip that output. The npy matrix has one row
// per (contig,bin,gameteSet) -- the same granularity as the PS4G data rows, not
// collapsed across gameteSets sharing a bin -- so gamete co-occurrence within a
// read is never lost. npy_binary writes presence (1) instead of the read count.
void rb3_ps4g_npy_finalize(rb3_ps4g_acc_t *acc, const rb3_gtab_t *gtab, const rb3_sid_t *sid,
							const char *ref_prefix, const rb3_bed_t *bed, int npy_binary,
							const char *ps4g_fn, const char *npy_fn, const char *cli_command);

#ifdef __cplusplus
}
#endif

#endif
