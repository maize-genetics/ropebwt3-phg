#ifndef RB3_ISA_H
#define RB3_ISA_H

#include <stdint.h>
#include "fm-index.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Sampled inverse suffix array ("entry points") for coordinate-based extraction.
 *
 * An entry (sid, pos, row) says: BWT row `row` is the suffix of BWT string `sid`
 * that starts at 0-based offset `pos`. LF-walking from that row yields the bases
 * at pos-1, pos-2, ..., 0 of that string, so the substring [st,en) is recovered by
 * starting at the first entry with pos >= en and walking pos-st steps.
 *
 * Entries are not computed by a new pass over the BWT: `ropebwt3 ssa` already
 * stores (row -> sid,pos) for every row with (row - acc[1]) % 2^ss == 0, so the
 * entries are exactly the inverted .ssa samples, thinned to at most one per
 * 2^shift-base bucket of each string (the sample with the smallest pos in the
 * bucket), plus one free exact entry per string: (pos = len, row = sid), the
 * string's own sentinel row (ropebwt3 numbers sentinels so that row sid ends
 * string sid; this is what `ropebwt3 get sid` walks from).
 *
 * BWT string numbering (`ropebwt3 build`, any of its input orders): sequence j of
 * the .len.gz file (0-based, the order the FASTA records were read) is BWT string
 * 2j (forward strand) and 2j+1 (its reverse complement). Both have length len[j].
 */
typedef struct {
	int32_t shift;   // entry spacing: at most one SSA-derived entry per 2^shift bases
	int32_t ss;      // the source .ssa's sampling shift (validation)
	int32_t strands; // bit 0: forward strings (even sid) sampled; bit 1: reverse (odd sid)
	int64_t m;       // number of BWT strings (== acc[1] == 2 * #sequences)
	int64_t n_ssa;   // the source .ssa's n_ssa (validation against the FMD)
	int64_t n;       // total number of entries
	int64_t *off;    // m+1: entries of BWT string sid are [off[sid], off[sid+1]), sorted by pos
	uint32_t *pos;   // n: 0-based offset on the string (== string length for the sentinel entry)
	int64_t *row;    // n: BWT row of the suffix starting at pos
} rb3_isa_t;

#define RB3_ISA_FWD 1
#define RB3_ISA_REV 2

// error codes of rb3_isa_extract()
#define RB3_ISA_ERR_SID    (-1) // sid out of range or its strand not sampled
#define RB3_ISA_ERR_RANGE  (-2) // st < 0, en > len, or st >= en
#define RB3_ISA_ERR_WALK   (-3) // LF walk hit a sentinel early: ISA does not match this FMD

// Build from an in-memory SSA. len[j] = length of sequence j (j < n_seq, from .len.gz).
rb3_isa_t *rb3_isa_from_ssa(const rb3_ssa_t *sa, int64_t n_seq, const int32_t *len, int shift, int strands);
// Same, but streams a .ssa file (only the header and the ssa[] array are read; the
// r2i[] table is skipped), so the 8-byte-per-sample SSA is never held in memory.
rb3_isa_t *rb3_isa_from_ssa_file(const char *fn, int64_t n_seq, const int32_t *len, int shift, int strands);

int rb3_isa_dump(const rb3_isa_t *isa, const char *fn);
rb3_isa_t *rb3_isa_restore(const char *fn);
void rb3_isa_destroy(rb3_isa_t *isa);

// 0 if the ISA was built from this FMD's SSA (same #strings and same #SSA samples
// for the FMD's length), else -1.
int rb3_isa_check_fmi(const rb3_isa_t *isa, const rb3_fmi_t *f);

// Length of BWT string sid as recorded by its sentinel entry, or -1.
int64_t rb3_isa_seq_len(const rb3_isa_t *isa, int64_t sid);

// Extract [st,en) (0-based, half-open) of BWT string sid into buf[0..en-st) as nt6
// codes (1..5 = A,C,G,T,N), in the string's own 5'->3' order. Returns the number
// of LF steps taken (>= en-st) or a negative RB3_ISA_ERR_* code.
int64_t rb3_isa_extract(const rb3_fmi_t *f, const rb3_isa_t *isa, int64_t sid, int64_t st, int64_t en, uint8_t *buf);

int main_isa(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
