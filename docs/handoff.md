# Handoff — `refmap`

Status: **implemented, tested, builds clean, ASAN-clean.** On branch
`feature/refmap`.

## What it does

`ropebwt3 refmap` places a query on a designated reference genome within a
multi-genome (pangenome) index, including when the query is absent from the
reference because it falls inside an insertion / reference-side deletion. See
`docs/plan.md` for the design and `docs/experiments.md` + `docs/results.md` for
testing.

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

## Suggested next steps

1. Validate on a real maize pangenome index (24 genomes); profile `refmap` time.
2. If slow, optimize `refmap_anchor_flank` membership search (see above).
3. Optional: emit BED/PAF-style output; multi-locus handling for rearrangements.
