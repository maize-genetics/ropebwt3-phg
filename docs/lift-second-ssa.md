# The lift "second SSA" — how the carrier→reference liftover is built

`ropebwt3 lift` builds a **second SSA**: a per-carrier→reference coordinate map
that lets `refmap --lift` place a read by *projecting* a carrier hit to the
reference instead of walking outward through the carriers. This document explains
the construction. For why it exists and the measured payoff (precision 57→95.5%,
recall 55→78%, ~15× faster) see [`results-maize.md`](results-maize.md) experiment
E4; the code is `lift.c` / `lift.h`.

Key idea in one line: **slide k-mers along the reference genome, and for every
k-mer that is single-copy across the pangenome, ask the *primary* SSA "where does
this land in each genome?" — each answer is one carrier↔reference coordinate
pair.** The second SSA is bootstrapped from the first.

## Inputs (`main_lift`)

- The FM-index with its **`.ssa` and `.len.gz`** loaded (`rb3_fmi_load_all`): we
  reuse the existing sampled suffix array to *locate*, and the sequence
  names/lengths to interpret coordinates.
- `--ref-prefix=B73` → an `is_ref[]` bitmap marking which sequences are the
  reference.
- `-m auto` → counts distinct taxa from name prefixes (`B73`, `B97`, …) = the
  single-copy-per-taxon cap (`max_occ`, 4 here).
- The reference FASTA (`B73.chr.fa.gz`) as the source of anchor k-mers.

## Per-anchor pipeline (`lift_worker`, parallel via `kt_for`)

Each thread takes one anchor — a length-`k` (=100) window at offset `i·stride`
(stride=2000) along a reference chromosome — and does four things:

1. **Backward-search the k-mer** through the BWT (`rb3_fmd_set_intv` + a loop of
   `rb3_fmd_extend`, is_back=1). This yields the k-mer's SA interval; `ik.size` is
   its total number of occurrences across all genomes (both strands).
2. **Uniqueness filter.** If `ik.size > max_occ`, skip — that's a retro/repeat
   (the same reads E2 discards). Keep only k-mers occurring at most once per taxon.
3. **Locate every occurrence** with the *primary* SSA: `rb3_ssa_multi(...)` walks
   each of the ≤`max_occ` interval rows to a sampled position, returning
   `(sid, pos)` per hit. `sid>>1` is the sequence index, `sid&1` the strand.
4. **Pair carrier↔reference.** Split hits by `is_ref[sid>>1]`. Require **exactly
   one reference occurrence** (if the anchor is duplicated in the reference it is
   ambiguous → skip). Then for **each carrier occurrence**, emit one point
   `(carrier_seq, carrier_pos) → (reference_seq, reference_pos)`. Both coordinates
   are normalized to **forward-strand starts** (`lift_fwd`: a reverse hit becomes
   `chr_len − (pos + k)`), so orientation is uniform downstream.

So one shared reference window present in B97 and CML247 yields two points:
`B97:x → B73:p` and `CML247:y → B73:p`. On the maize index: **1,065,928 k-mers
scanned → 897,784 liftover points** (the rest are repeats, or reference-unique
with no carrier).

## Assembly into a queryable map (`lift_builder_finalize`)

- Concatenate the per-thread point buffers into one array.
- **Sort by `(carrier_seq, carrier_pos)`** so each carrier chromosome's points are
  in coordinate order.
- Build a CSR-style **offset index** `off[]`: `off[sid]..off[sid+1]` is the point
  range for carrier sequence `sid`. Projection is then a fast binary search —
  jump to the carrier's slice and bisect the window around the query position.

This sorted-array-plus-offset-index *is* the second SSA: a per-carrier coordinate
table, parallel in spirit to how `.ssa` is a per-row sampled table.

## On-disk format (`rb3_lift_dump`)

```
"LIFT\1"                      magic
int64  n_seq                  number of sequences
int64  n_pt                   number of liftover points
int64  off[n_seq+1]           CSR offsets per carrier sequence
{ int64 cpos; int64 rpos; int32 csid; int32 rsid } pt[n_pt]   sorted points
```

21 MB for the maize index; `rb3_lift_restore` reads it straight back at
`refmap --lift` startup.

## Why it is built this way

- **Slide over the reference, not the carriers.** Every pair needs a reference
  coordinate, so anchoring on the reference guarantees it and touches one genome
  instead of all; by symmetry it captures the same shared anchors.
- **Single-copy-per-taxon (`max_occ`)** keeps only trustworthy anchors; retros
  never enter the map, so they cannot manufacture false collinearity at query time.
- **Only shared sequence gets points.** Carrier-specific insertions have no nearby
  reference anchor — correct: the query-time projection interpolates across the gap
  from flanking shared anchors, and returns NULL (→ `UNPLACED`) when the flanks
  disagree.
- **Reuses the primary SSA + FM-index** for both search and locate, so there is no
  new heavy machinery — the second SSA is a thin, sorted coordinate table derived
  from the first.

## The query side (mirror image, `rb3_lift_project`)

Given a carrier hit `(sid, pos)`: bisect `off[sid]`'s slice to gather points within
`±lift_win` bp; take the majority reference sequence; decide orientation (+1/−1)
by whichever residual `rpos∓cpos` is tighter; project with the median intercept.
Return NULL unless there is enough collinear support and the residual **MAD ≤
lift_mad** — the MAD gate is what protects precision in rearranged/insertion
regions.

## Build & use

```sh
ropebwt3 lift --ref-prefix=B73 -k 100 -s 2000 -t 20 -o nam4.lift nam4.fmd B73.chr.fa.gz
ropebwt3 refmap --ref-prefix=B73 --max-occ=-1 --lift nam4.lift -t 20 nam4.fmd reads.fa
```
