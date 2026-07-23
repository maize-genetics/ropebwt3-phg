# Stage-0 exploration results

Stock `ropebwt3 refmap` baseline on the Module-0 simulator (5-genotype synthetic
pangenome, B73 source). Scoring per `eval/score.py` (strict per-read attribution).
These measurements motivate `NEXT_STEPS.md` (uniting a read's SMEMs).

## 1. Minimum-SMEM-length (`-l`) is the wrong knob (HANDOFF decision 1)

8k within-exon reads, seed 11. `-l` gates the exact core used to place an
error-containing read; error-free reads match end-to-end and are unaffected.

**1% substitution error**

| `-l` | recall | misplaced | unmapped | spurious founders |
|--:|--:|--:|--:|--:|
| 19 (default) | 99.4% | 0.6% | 0.0% | 31.1% |
| 30 | 99.4% | 0.6% | 0.0% | 31.1% |
| 40 | 99.3% | 0.6% | 0.2% | 31.0% |
| 50 | 97.8% | 0.6% | 1.6% | 30.3% |
| 60 | 93.3% | 0.5% | 6.2% | 28.2% |

**2% substitution error**

| `-l` | recall | misplaced | unmapped | spurious founders |
|--:|--:|--:|--:|--:|
| 19 (default) | 98.8% | 1.2% | 0.0% | 49.4% |
| 40 | 96.2% | 1.2% | 2.6% | 48.5% |
| 50 | 88.8% | 1.1% | 10.1% | 46.1% |
| 60 | 76.5% | 0.9% | 22.5% | 41.7% |

**Finding.** Raising `-l` barely improves specificity — at 1% error `spurious
founders` fall only 31→28% and `misplaced` is already ~0.6%; at 2% error spurious
falls only 49→42% — but it sharply raises `unmapped` (to 6% at 1% / **22%** at 2%
error) and craters `recall` (to 93% / **77%**) as the floor climbs past ~40. So the
single-longest-core design has a **specificity floor `-l` cannot break**: spurious
founders come from the core matching *other founders over its span* (conserved
exons), not from short off-target hits, so a longer floor just discards reads.
Reducing spurious without losing coverage needs the read's founder set constrained
by **all** its exact segments (chaining across error gaps), not one core.
Practical default: keep `-l` low (~20–30); do not raise it to chase specificity.

## 2. Exon-spanning (junction) reads are badly served (motivates splice segmentation)

8k reads, `--max-junctions 2`, 1% error → 48% within-exon, 49% one-junction, 3%
two-junction. `eval/analyze_junctions.py`:

| class | reads | segs/read | recall (near-exon) | segment coverage cap | valid emitted base |
|--|--:|--:|--:|--:|--:|
| within-exon | 3857 | 1.00 | 99.5% | ≤100% | 58.3% |
| one-junction | 3932 | 2.00 | 75.7% | ≤50% | 8.5% |
| two-junction | 211 | 3.00 | 95.7% | ≤33% | 5.4% |

Overall score on this 51%-junction set vs a within-exon set: recall 99.4→**87.7%**,
misplaced 0.6→**12.3%**.

**Two distinct failures on spliced reads:**
- **Structural** — stock refmap emits *one* PS4G position per read, so a read
  spanning N exons covers at most 1 of them: per-segment coverage is capped at
  **1/N** (≤50% one-junction, ≤33% two-junction). The read's other exons get no
  evidence (this also shows up as inflated `missed-IBS` in the overall score,
  since only one exon's founders are emitted).
- **Coordinate** — the whole-read `[cL,cR)` is *extrapolated* from the core as if
  the read were contiguous, which is invalid across a junction. So the one emitted
  base usually lands **off** the true exons: valid-base rate collapses from 58%
  (within-exon) to 8.5% (one-junction). Junction reads deposit PS4G rows at
  invalid (often intronic) positions, and 24% miss every true exon entirely.

Note the 58% valid-base for *within-exon* reads: the same extrapolation shifts
`-`-strand error reads a few–tens of bp below the true start, so even clean reads
often don't land on the exact start base (interval-tolerant recall is still 99.5%,
but single-base PS4G positions are imprecise). Emitting each SMEM at its **own**
start base (segmentation) fixes this too.

## Conclusion

Both experiments point to the same fix: **stop relying on one core per read.**
Collect all of a read's SMEMs, chain the colinear ones, refine the founder set
across the chain, and emit **per exon segment** at each segment's own start base.
This lowers `spurious` (set constrained by the whole read), keeps `-l` low
(chaining supplies specificity), covers **all** exon segments of spliced reads,
and puts every row at a valid base. Design in `NEXT_STEPS.md`.
