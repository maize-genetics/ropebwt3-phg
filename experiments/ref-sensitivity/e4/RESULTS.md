# E4 prototype — synteny-prior placement via carrier→B73 liftover

Goal: replace refmap's slow outward **walk** with a coordinate **projection**.
A carrier hit's exact coordinate comes free from the SSA at locate time; project
it to an approximate B73 coordinate through a precomputed carrier→B73 liftover,
and mark **NULL** (→ UNPLACED) wherever no confident projection exists, to keep
precision (per the "two SSA + null slots" design).

This is a **Python prototype** over `ropebwt3 mem -p` output (not yet the
integrated C "second SSA"), built only from the existing index — no external
aligner, no synteny data.

## Pipeline
1. `gen_anchors.py` — sample anchor windows along each carrier (stride S, len L).
2. `mem -p` the anchors → for every single-copy-per-taxon anchor with one B73
   hit, a pair (carrier_pos ↔ B73_pos).
3. `liftover.py` — sort pairs per (carrier,chr); `project()` interpolates between
   flanking anchors when they agree on chr and are locally collinear (slope≈1),
   else returns NULL.
4. `place_reads.py` — EXACT if a B73 hit; else project carrier hits, place at the
   agreed chr/median, else UNPLACED; MULTI if occ > max_occ.

## Liftover self-test (leave-one-out, stride 10k, L=100, max_occ=4)
- 179,517 usable liftover points from 639,501 anchors (28.1%).
- Projection: 85.7% committed, **14.3% NULL**; of committed, **99.76% within 5 Mb**.
- |Δ B73| median **7 bp**, p90 5.3 kb, p99 44 kb. Interpolation is very accurate.

## Read placement vs refmap (100k reads, tol ±5 Mb)

| variant | precision | recall | wrong-chr | placement speed |
|---|---:|---:|---:|---:|
| E0 baseline refmap | 57.3% | 54.6% | 35,768 | 40 s (walk) |
| E3 two-flank+max-occ=4 | 96.0% | 45.5% | 920 | 37 s (walk) |
| **E4 liftover (stride 10k, occ4)** | **95.9%** | **63.7%** | 1,068 | ~3 s locate + project |

By status: EXACT 42,242 (96.7% correct) · PLACED 24,218 (**94.5% correct**) ·
MULTI 18,062 (filtered) · UNPLACED 15,045 (null/disagree) · AMBIG 433.

**Takeaway:** the synteny prior beats the walk on precision, recall, AND speed at
once. The carrier-projected PLACED bucket is 94.5% correct vs refmap's 14.7%
ONE_SIDE. Recall +18 pts over E3 at the same ~96% precision, ~15× faster.

## Parameter experiments

### Two-anchor projection: stride × max_occ
- `max_occ=4` (=#taxa) is the knee for E4 too: precision 95.9→94.4% as occ 4→8,
  and recall *falls* (63.7→59.7%). Same verdict as refmap.
- Counter-intuitively, **sparser anchors gave higher recall** with the naive
  two-anchor projection (stride 10k 63.7% > 5k 60.2% > 2k 55.9%): a short baseline
  makes one local indel blow up the apparent slope → spurious NULL. The fix is a
  windowed fit, not sparse anchors.

### Robust windowed projection (the win)
`project_robust`: gather anchors within ±`win` of the query, take the majority
B73 chr, choose orientation (+1/−1) by the tighter residual `b∓c`, project with
the median intercept; NULL unless ≥`min_support` anchors and residual MAD ≤
`max_mad`. Exploits the within-species **slope ≈ ±1** prior and the MAD gate
protects precision.

| projection (occ4) | recall | precision |
|---|---:|---:|
| two-anchor, stride 10k | 63.7% | 95.9% |
| two-anchor, stride 2k | 55.9% | 96.0% |
| **robust, stride 2k, win 500k, mad 200k** | **77.6%** | **95.9%** |

Robust is **stride-insensitive** (now denser helps: 10k 76.2% → 2k 77.6%) and
flat across win∈{0.5–1 M}, mad∈{0.2–0.5 M} (recall 77.6–77.9%, precision ~95.9%).
With robust projection only ~1,000 reads NULL out — recall is now capped by the
**18,062 `MULTI`** (occ>4 repeat) reads, which are genuinely unplaceable.

## Final comparison (100k reads, tol ±5 Mb)

| variant | precision | recall | wrong-chr | placement |
|---|---:|---:|---:|---:|
| E0 baseline refmap | 57.3% | 54.6% | 35,768 | 40 s walk |
| E3 two-flank+max-occ=4 | 96.0% | 45.5% | 920 | 37 s walk |
| E4 liftover (two-anchor) | 95.9% | 63.7% | 1,068 | ~3 s locate |
| **E4 liftover (robust)** | **95.9%** | **77.6%** | 1,177 | ~3 s locate |

**E4 robust dominates:** same ~96% precision as the best walk-based variant, but
**+32 pts recall** over E3 and ~15× faster. The carrier-projected `PLACED` bucket
is ~94% correct vs refmap's 14.7% `ONE_SIDE`.

## Robustness checks
- **Speed (measured):** projection step (load 3.2 M-point liftover + place 100k)
  = 5.45 s single-thread Python, 237 MB RSS; with `mem -p` locate (2.67 s, 20 thr)
  end-to-end ≈ 8 s vs refmap 40 s (~5×). The liftover load is one-time (like
  `.ssa`); marginal projection is cheap, so a C "second SSA" should near `mem`.
- **Generalization (fresh seed):** seed-13 reads (unseen by the anchors) →
  precision 96.0%, recall 77.8%, B73ctrl 80.6% — matches seed-7
  (95.9% / 77.6% / 80.7%). Not overfit.

## Read-length sensitivity (E4 robust, occ4, stride 2k)

| read len | MULTI (repeat) | recall | precision | B73 ctrl |
|---:|---:|---:|---:|---:|
| 75 bp | 39.3% | 56.7% | 95.0% | — |
| 100 bp | 29.5% | 66.3% | 95.4% | — |
| 150 bp | 18.1% | 77.6% | 95.9% | 80.7% |
| 250 bp | 8.6% | 87.4% | 96.5% | 90.7% |

Recall is governed by the repeat (`MULTI`) fraction, which shrinks with read
length; **precision is rock-stable (~95–96%, even improving) at every length.**
Practical implication: longer reads buy large recall gains (56.7→87.4% from
75→250 bp) at no precision cost — use the longest reads available.

## Locked prototype config
`--max-occ 4` (auto=#taxa), anchors stride 2000 len 100, `--robust --win 500000
--max-mad 200000`.

## Precision under sequencing error: 75 bp k-mer agreement

Real reads carry errors, so a 150 bp read rarely matches the index exactly; the
whole-read placement then falls back to a short, less-unique core and loses
precision. Splitting the read into overlapping 75 bp k-mers and requiring
**≥ N of them to agree** on a locus (only full-length exact k-mer hits count, so
an error-containing k-mer simply contributes nothing) tolerates error and raises
precision. Scripts: `add_errors.py`, `gen_kmers.py`, `place_kmers.py`,
`kmer_exp.sh`.

Precision by substitution rate (whole-read `--lift` vs 75 bp, step 15, agree ≥ 2):

| error | whole-read prec / recall | k-mer≥2 prec / recall |
|---:|---:|---:|
| 0% | 95.5% / 77.9% | 96.1% / 68.6% |
| 1% | 85.1% / 65.3% | 91.5% / 39.5% |
| 2% | 75.7% / 55.6% | 88.1% / 19.9% |
| 3% | 68.4% / 48.5% | 85.2% / 12.0% |

**The k-mer precision advantage grows with error rate** (+0.6 → +16.8 pts from
0→3%): whole-read degrades fast as exact 150 bp matches vanish, while k-mer
agreement holds precision high. The cost is recall (a conservative, high-precision
operating point). More agreement (≥3) or denser tiling trades recall for still
higher precision. Good fit for genotyping, where a wrong call is costlier than a
missing one.

**Integrated in C** as `refmap --kmer 75 [--kmer-step --min-agree --kmer-cluster]`
(with `--lift`): tiles each read, projects each k-mer, and places only where
`--min-agree` distinct k-mers agree. Matches the prototype (error-free
96.0%/68.5%, 1% error 91.5%/39.5%), ~4-5 s / 100k reads, ASAN-clean.

## Integrated as the C "second SSA" (DONE)
The prototype was integrated into ropebwt3: `ropebwt3 lift` builds the
carrier→reference map and `refmap --lift` projects at the locate step (no walk),
with NULL slots → UNPLACED. Measured `refmap --lift`: precision 95.5%, recall
77.9%, 2.76 s / 100k reads (~15×), ASAN-clean. Construction walkthrough:
[`../../../docs/lift-second-ssa.md`](../../../docs/lift-second-ssa.md).

