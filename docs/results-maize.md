# Results — `refmap` maize-pangenome benchmark

Measured runs of the benchmark in [`benchmark-maize.md`](benchmark-maize.md):
4 NAM founders (B73 ref + B97/Ki3/CML247), one FMD index, **100,000** 150 bp
reads (25k/line, seed 7), `-t 20`. The scoreboards accumulate one row per variant
(E0 baseline measured; E1–E3 to come).

## Scoreboard — sensitivity (tol ±5 Mb)

precision = correct / placed; recall = correct / 100,000; placed = reads asserting
a coordinate (not UNPLACED/MULTI). Canonical numbers from `analyze.py` ALL row.

| ID | Variant | placed | correct | wrong-chr | **precision** | recall |
|---|---|---:|---:|---:|---:|---:|
| **E0** | baseline | 95,218 | 54,545 | 35,768 | **57.3%** | 54.6% |
| **E1** | `--two-flank` | 66,555 | 49,736 | 14,263 | **74.7%** | 49.7% |
| **E2** | `--max-occ=4` | 56,858 | 50,083 | 5,219 | **88.1%** | 50.1% |
| **E3** | `--two-flank --max-occ=4` | 47,337 | 45,459 | **920** | **96.0%** | 45.5% |
| **E4** | `--lift` synteny prior (C, integrated) | 81,536 | 77,893 | 1,417 | **95.5%** | **77.9%** |

**Takeaways.** Both fixes do exactly what the diagnosis predicted, and stack:
- **E1** removes `ONE_SIDE`: wrong-chr 35.8k → 14.3k (−60%), precision 57→75%, at
  ~10% recall cost.
- **E2** removes retro/repeat reads (18,062 flagged `MULTI`): wrong-chr → 5.2k
  (−85%) at almost **no** recall cost (54.5→50.1%) — repeats were nearly all
  wrong anyway, so filtering them is close to free. Best precision-per-recall.
- **E3** (both): wrong-chr **35,768 → 920 (−97%)**, **precision 96.0%**, recall
  45.5%. The right operating point for genotyping (precision ≫ recall).
- `--max-occ=4` (= #taxa) also trims B73 control 86.6→81.3% — a few correctly
  placed shared loci (occ 5–8) get flagged `MULTI`; a slightly higher threshold
  would recover them. **Sweep `--max-occ ∈ {4,6,8,…}` is the obvious next tuning.**
- **E4 (synteny prior) is the breakthrough, now implemented in C.** It replaces
  the slow walk with a carrier→B73 coordinate projection (the "second SSA", built
  once by `ropebwt3 lift`). Integrated `refmap --lift`: **precision 95.5%, recall
  77.9% (+32 pts over E3), and 2.76 s for 100k reads (~15× faster than the walk,
  essentially `mem` speed).** ASAN-clean. The liftover build is a one-time
  index-side cost (21.6 s over B73, 21 MB file), analogous to `.ssa`. How the
  liftover is built: [`lift-second-ssa.md`](lift-second-ssa.md). Method +
  prototype sweeps in
  [`../experiments/ref-sensitivity/e4/RESULTS.md`](../experiments/ref-sensitivity/e4/RESULTS.md).
  Recall is now capped by the ~18% genuinely-repetitive (`MULTI`) reads.

## Scoreboard — speed (100k reads, `-t 20`)

| ID | Tool / variant | Wall | Reads/sec | CPU% | CPU-s | Peak RSS | vs `mem` |
|---|---|---:|---:|---:|---:|---:|---:|
| — | `mem` (SMEM+locate) | 1.16 s | 86,207 | 973% | 11.3 | 2.1 GB | 1× |
| — | `sw` (local align) | 9.44 s | 10,593 | 1816% | 171.4 | 2.8 GB | 8.1× |
| **E0** | `refmap` baseline | 40.44 s | 2,473 | 1948% | 787.8 | 2.5 GB | **35× wall / 70× CPU** |
| **E1** | `refmap --two-flank` | 39.50 s | 2,532 | — | — | — | 34× |
| **E2** | `refmap --max-occ=4` | 37.44 s | 2,671 | — | — | — | 32× |
| **E3** | `refmap` E1+E2 | 37.47 s | 2,669 | — | — | — | 32× |
| **E4** | `refmap --lift` (C) | **2.76 s** | **36,232** | — | — | 2.6 GB | **1.03×** |

Speed read: `refmap` is ~35× slower than `mem` by wall (≈70× by CPU-seconds, since
`mem` finishes too fast at 1.16 s to fully use 20 threads) and ~4.3× slower than
`sw`. 60% of reads are `EXACT` (fast `mem`-like path); the wall time is paid by
the ~40% that walk + anchor, i.e. ≈1 ms/read on the placement path. **The E1–E3
precision fixes are essentially speed-neutral** (E2/E3 slightly faster: `MULTI`
reads skip the walk). The walk itself is untouched — the real speed lever is the
planned **E4 synteny prior** (use the carrier coordinate the SSA already returns
to skip/shorten the walk and reject off-diagonal carriers).

## `--max-occ` sweep (knee finding)

precision = correct/placed; recall = correct/100k; B73ctrl = B73-control %correct.

| --max-occ | solo prec | solo recall | +two-flank prec | +two-flank recall | B73 ctrl |
|---:|---:|---:|---:|---:|---:|
| 3 | 85.9% | 37.5% | 96.9% | 31.6% | 68.7% |
| **4** (=#taxa) | 84.2% | **54.0%** | **96.0%** | 45.5% | 81.3% |
| 6 | 76.6% | 55.3% | 94.2% | 47.1% | 83.3% |
| 8 | 72.5% | 55.8% | 92.6% | 48.0% | 84.4% |
| 12 | 68.5% | 55.9% | 90.4% | 48.6% | 85.2% |
| 16 | 66.5% | 56.0% | 89.0% | 48.9% | 85.6% |

**Verdict: `--max-occ = n_taxa` (auto default) is the knee.** Recall plateaus at
occ=4 (solo 54.0% → only 56.0% at occ=16) while precision falls steadily
(84→66%) — reads occurring 5–16× are overwhelmingly repeats that place wrong, so
raising the cap buys wrong placements, not recall. The earlier guess that a
higher cap would recover recall was wrong. The one real lever: with `--two-flank`,
occ 6–8 trades ~3 pt precision for ~2.5 pt recall and recovers a little B73
control (81→84%) — a recall-leaning operating point if wanted. Default stays
`auto` (=4 here).

---

# E0 baseline — full result

## Sensitivity headline (`data/report.txt`)

| line | n | correct | %correct | wrong_chr | far(>tol) | unplaced |
|---|---:|---:|---:|---:|---:|---:|
| B73 (control) | 25,000 | 21,660 | **86.64%** | 2,982 | 358 | 0 |
| B97 | 25,000 | 11,715 | 46.86% | 10,325 | 1,504 | 1,456 |
| Ki3 | 25,000 | 10,626 | 42.50% | 11,060 | 1,670 | 1,644 |
| CML247 | 25,000 | 10,544 | 42.18% | 11,401 | 1,373 | 1,682 |
| **ALL** | 100,000 | 54,545 | **54.55%** | 35,768 | 4,905 | 4,782 |

Status mix (all reads): **60,089 EXACT · 32,719 ONE_SIDE · 6,466 PLACED · 726 UNPLACED.**

## Diagnosis — where the error comes from

The error is not diffuse; it concentrates in two explainable places.

### 1. Read multiplicity (repeats) — the control deficit is entirely repeats

B73 reads by occurrence count (exact full-length matches in the index):

| occ bucket | n | %correct | %wrong-chr |
|---|---:|---:|---:|
| 1 (unique) | 7,499 | **100.0%** | 0.0% |
| 2–5 | 13,326 | 98.2% | 1.7% |
| 6–20 | 1,696 | 51.4% | 43.3% |
| 21–100 | 1,137 | 13.0% | 77.7% |
| >100 | 1,342 | 4.4% | 85.0% |

Reads drawn from B73 place **perfectly when unique** and ~98% when low-copy. The
entire 13% control deficit is reads with occ ≥ 6 (retros), where `refmap` reports
one arbitrary copy. Only **31% of 150 bp reads are unique**; 52% are occ 2–5.
`refmap` emits no multiplicity/MAPQ to flag the repeats.

### 2. The `ONE_SIDE` single-anchor path — 60% of all errors

Carrier reads (75k) by status:

| status | % of carriers | %correct | %wrong-chr |
|---|---:|---:|---:|
| EXACT (shared with B73) | 46.8% | 65.0% | 29.6% |
| **PLACED (both flanks anchor)** | 8.6% | **81.3%** | 13.9% |
| **ONE_SIDE (one flank anchors)** | 43.6% | **14.7%** | **65.7%** |
| UNPLACED | 1.0% | — | — |

When **both** flanks anchor (`PLACED`) refmap is 81% correct — the coordinate math
is sound. The failure is `ONE_SIDE`: a single-flank anchor with no concordance
check lands on the wrong chromosome 2/3 of the time, and it is 44% of carrier
reads. **`ONE_SIDE` alone contributes 21,505 of the 35,768 wrong-chr errors (60%).**

## Conclusions → motivated experiments

1. **Lone single-flank anchors are untrustworthy.** → **E1**: require two-flank
   concordance (both anchored, same chr, collinear); demote `ONE_SIDE`.
   Code: `refmap_place` (`search.c:282`).
2. **Anchors and reads should be single-copy-per-taxon, not just matching.** →
   **E2**: `--max-occ` (default = #taxa); extend anchors until reference interval
   `< max-occ`; flag reads with whole-query occ `≥ max-occ` as `MULTI`.
   Code: `refmap_anchor_flank` (`search.c:317`), `refmap_query_interval`.
3. Both are **precision-over-recall** trades — appropriate for genotyping. E1–E3
   will report the precision gain, the yield cost, and any speed change (the
   uniqueness early-exit may also shorten anchoring).

Caveat: `wrong_chr` is genuine error; the smaller `far(>tol)` bucket and the
~1–2 Mb median carrier distances partly reflect *real* inter-genome coordinate
drift, so the distance-based numbers modestly under-state true accuracy.
