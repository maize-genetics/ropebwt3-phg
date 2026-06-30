# Handoff — `refmap`

Status: **implemented, tested, builds clean, ASAN-clean.** Merged to `master`
and published at `maize-genetics/ropebwt3-phg` (a fork of `lh3/ropebwt3`).

## What it does

`ropebwt3 refmap` places a query on a designated reference genome within a
multi-genome (pangenome) index, including when the query is absent from the
reference because it falls inside an insertion / reference-side deletion. See
`docs/plan.md` for the design and `docs/experiments.md` + `docs/results.md` for
testing.

The algorithm in one picture — [`docs/refmap-walk.svg`](refmap-walk.svg): the
query never touches the reference, so refmap locates it in the carriers, walks
outward through a carrier to the insertion breakpoints, and re-anchors the shared
flanks (longest reference-matching prefix / suffix) back onto the reference.

## Build

Plain `make` fails on macOS: Apple `clang` does not accept `-fopenmp` (used only
by `libsais`). Build with:

```sh
make omp=0            # normal optimized build
make omp=0 asan=1     # AddressSanitizer build for debugging
```

On Linux with real gcc, plain `make` works.

## Run

```sh
# Index must have the sampled suffix array and sequence names:
ropebwt3 ssa -o pan.fmd.ssa -s8 pan.fmd
cat *.fa | seqtk comp | cut -f1,2 | gzip > pan.fmd.len.gz
ropebwt3 refmap --ref-prefix=B73 pan.fmd query.fa
```

Options: `--ref-prefix=STR` (required), `--max-walk=NUM` (default 5000),
`--walk-mode=consensus|strict|per-carrier`, `-l INT` (min anchor length),
`-t INT` (threads).

## Code map

All new code is gated behind the `refmap` subcommand; `mem` / `sw` / `hapdiv`
are untouched.

| File | Change |
|---|---|
| `ssa.c` | `rb3_ssa_multi_ref()` + `ssa_add_intv_ref()` — reference-filtered locate with a scan budget. Added `is_ref`/`scan_left` to `ssa_aux_t`. |
| `fm-index.h` | declaration of `rb3_ssa_multi_ref()`. |
| `search.c` | the feature (see below). |
| `main.c` | `refmap` in the command list and dispatch. |
| `README.md` | "Placing a query on a reference genome" section. |

### `search.c` walkthrough

* `rb3_mopt_t` gains `ref_prefix`, `max_walk`, `walk_mode`; defaults in
  `rb3_mopt_init`. New algo `RB3_SA_REFMAP`. Walk modes `RB3_WALK_*`.
* `pipeline_t` gains `is_ref` (per-sequence reference bitmap, indexed by
  `sid>>1`) and `n_ref`; built in `main_search` after `rb3_fmi_load_all`.
* `refmap_rst_t` — per-query result (status, strand, ref_sid, cL, cR, ins_size,
  carriers, and `sub`/`n_sub` for per-carrier mode). One per query in
  `step_t.refmap`.
* `worker_for_seq` → `refmap_query` for `RB3_SA_REFMAP`.
* `refmap_query_interval` — backward-search the full query to its SA interval
  (`rb3_fmd_set_intv` + `rb3_fmd_extend` is_back=1).
* `refmap_query` — locate carriers (`rb3_ssa_multi`, then `is_ref` split). If a
  reference hit is in the interval → `EXACT`. Else extract + anchor flanks. On
  no end-to-end match, falls back to the longest SMEM core.
* `refmap_extract_flank` — outward walk. `mask != NULL` ⇒ per-carrier (follow the
  child interval still containing that carrier, via `rb3_ssa_multi_ref` with a
  one-bit mask). Otherwise consensus (largest child) or strict (stop when the
  chosen child is smaller than the parent = a carrier disagreed). Returns the
  flank in forward orientation.
* `refmap_anchor_flank` — grows a match from the flank end *far* from the query
  and keeps the longest stretch still present in the reference (prefix for the
  left flank via forward extension, suffix for the right flank via backward
  extension), then locates it. Returns the inner-edge reference coordinate
  (strand-aware) and matched length.
* `refmap_place` — extract both flanks, anchor, fill the result
  (PLACED / ONE_SIDE).
* `write_refmap` / `write_refmap1` — emit the line(s).

## Gotchas for the next person

* **Coordinates are strand-aware.** `pos_stranded()` converts a located
  `rb3_pos_t` to a forward-strand reference interval. On a reverse-strand hit the
  breakpoint edge swaps (handled in `refmap_anchor_flank`).
* **Microhomology** at a breakpoint makes `cL`/`cR` differ by a few bp — this is
  correct, not a bug (the inserted sequence shares terminal bases with the
  flanks). `refSpan = cR - cL` reports the ambiguity zone.
* **Per-carrier `sub[].carriers` is a borrowed pointer** into the parent
  `carriers` array; `write_refmap` frees the parent once and never the borrowed
  copies. Memory ownership confirmed under ASAN.
* Anchoring does base-by-base reference-membership checks (O(flank)). Fine for
  the kb-scale insertions targeted here; for very long flanks at full
  pangenome scale, switch to exponential+binary search on the anchor length.

## Open questions / next investigations

These are research directions for the next person, not yet started. The first
three are the priorities flagged after the initial implementation.

### 1. Speed and the locate tradeoff on real genomes

The synthetic tests are tiny; the real cost is unknown on maize-scale data
(chromosome-length sequences × ~24 genomes). The likely hot spot is
`refmap_anchor_flank`, which today does a base-by-base reference-membership check
(`rb3_ssa_multi_ref` per extension step) — `O(flank length)` locate calls, each a
BWT walk bounded by `max_scan`. The knobs that trade speed against memory/quality:

* **SSA sample rate** (`ssa -s`): denser sampling → faster `locate`, larger
  `.ssa`. ropebwt3 itself notes locate is its slow operation.
* **`--max-walk`**: caps walk length; larger handles bigger SVs but costs more.
* **`max_scan`** budget inside `rb3_ssa_multi_ref` (currently `1<<16`).
* **Anchoring search**: replace the linear per-step membership scan with
  exponential + binary search on the anchor length → `O(log flank)` locates.

Needs a real benchmark: build a true multi-genome index, time `refmap` across
query sets, and profile where time goes (walk vs. membership-locate vs. the final
locate). Record index/SSA/`.len.gz` sizes and peak RSS. **These tradeoffs need
careful, measured consideration before scaling** — don't tune blind.

### 2. Is the BWT global, or can it be biased toward synteny?

The BWT is **global**: ropebwt3 concatenates all sequences (both strands) and
sorts suffixes lexicographically. It has no notion of synteny or position — a
match is found wherever the exact string occurs, in any genome. Build ordering
options (`build -s` RLO, `-r` RCLO) change run-length compression but **not** which
matches exist. So you cannot make the index itself "syntenic." Emphasis on
syntenic regions would have to live in a layer *around* the BWT, e.g.:

* constrain carriers/anchors to be **collinear** with the reference anchor
  (chain the located flank positions; reject off-diagonal carriers) — a
  query-time synteny constraint;
* partition genomes into synteny blocks and index per block;
* filter carriers by chromosome via name prefix before walking.

Open question worth scoping: does adding a collinearity check in `refmap_place`
(compare each carrier's flank coordinates against the reference anchor) improve
placement in repeat-rich / paralogous regions?

### 3. Quantify the "messiness" of the sequence among the query's carriers

The outward walk already sees, at every backward-extension step, the distribution
of child-interval sizes — i.e. how many carriers agree on the next base
(`refmap_extract_flank`: consensus picks the largest child, strict stops at the
first disagreement). That signal can be turned into a **divergence score** for the
locus, for example:

* fraction of walk steps where carriers disagreed (consensus ≠ unanimous);
* per-step base-distribution entropy, summarised over the flank;
* number of distinct carrier alleles (from `--walk-mode=per-carrier`);
* **spread of `cL`/`cR` across carriers** in per-carrier mode — consistent
  breakpoints = clean locus, scattered = messy/ambiguous.

Concretely: add an output column (e.g. divergent-step count or allele count) fed
from the size distribution in `refmap_extract_flank`, and/or compare the
consensus result against the per-carrier results to flag noisy loci.

### Also

* Emit BED/PAF-style output; multi-locus handling for large rearrangements.
* Offer `refmap` upstream to `lh3/ropebwt3` as a PR if of general interest.
