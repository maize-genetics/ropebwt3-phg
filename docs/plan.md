# Plan — `refmap`: place a query on a reference genome within a pangenome

## Context / motivation

ropebwt3 indexes many genomes of a species together in one BWT/FM-index. For a
plant / maize pangenome (~24 genomes) we designate **one** genome as the
reference and want to answer:

> Given a query sequence (e.g. ~150 bp), where does it sit on the reference
> genome — even when the query is **absent** from the reference because it lies
> inside an insertion (equivalently, a deletion on the reference side)?

The existing `mem` and `sw` commands cannot do this. Every exact match of such a
query occurs only in the *other* genomes (the carriers), never in the reference,
so no amount of locating that query yields a reference coordinate.

## Key idea

1. Locate the query in the genomes that carry it (the carriers).
2. **Walk outward** along the carriers (BWT backward-extension) past the query
   until each flank crosses the insertion breakpoint and rejoins
   reference-shared sequence.
3. Re-anchor each flank in the reference. The reference-shared part of a flank
   sits at the end *far* from the query, so a plain SMEM is swallowed by the
   longer carrier match that crosses the breakpoint. Instead we take the
   **longest reference-matching prefix** (left flank) / **suffix** (right flank).
4. The two reference anchors bracket the query's approximate location; their
   inner edges are the insertion breakpoints.

## Decisions (agreed with user)

| Question | Decision |
|---|---|
| How to designate the reference | by sequence-name **prefix** (`--ref-prefix=STR`) |
| "Nearest" when the query is absent from the reference | walk outward through carriers to the breakpoints (not chain query-internal SMEMs — the deletion case has none) |
| How far to walk | until anchored, capped at `--max-walk` (default 5000), configurable |
| Carrier divergence in flanks | `--walk-mode` = `consensus` (default) / `strict` / `per-carrier` (all three) |
| Output | one per-query reference interval line |

## Implementation summary

New subcommand `refmap`, reusing the existing `main_search` pipeline.

* **`ssa.c` / `fm-index.h`** — `rb3_ssa_multi_ref()`: a reference-filtered
  variant of `rb3_ssa_multi()` (BWT walk that only returns hits in reference
  sequences, with a `max_scan` budget so confirming a rare/absent reference hit
  stays cheap).
* **`search.c`** — the bulk:
  * options `--ref-prefix`, `--max-walk`, `--walk-mode`; `RB3_SA_REFMAP` algo.
  * `is_ref[]` bitmap built from sequence names after index load.
  * `refmap_query_interval()` — backward-search the whole query (longest-SMEM
    fallback if it doesn't match a carrier end to end).
  * `refmap_extract_flank()` — outward consensus / strict / per-carrier walk.
  * `refmap_anchor_flank()` — longest reference-matching prefix/suffix + locate.
  * `refmap_place()` / `refmap_query()` — orchestrate; `write_refmap()` output.
* **`main.c`** — register/dispatch `refmap`.
* **`README.md`** — user-facing section.

## Output format

One tab-separated line per query (one per carrier in `per-carrier` mode):

```
qname  qlen  status  nCarrier  carriers  refName  refStrand  cL  cR  refSpan  insSize
```

`status` ∈ `PLACED` | `ONE_SIDE` | `EXACT` | `UNPLACED`. `cL`/`cR` bracket the
query on the reference (equal for a clean breakpoint; a few bp apart with
microhomology). `EXACT` means the query is in the reference itself.

## Out of scope (future)

* Base-by-base reference-membership checks during anchoring are O(flank);
  replace with exponential + binary search for very long flanks at scale.
* Multiple chains per query (large structural rearrangements).
* Regex / multi-prefix reference selection.
* `insSize` is approximate when a query is placed via a partial (SMEM) core.
