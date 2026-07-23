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

## 3. Chaining vs stock — full sweep (set-based colinear chaining)

`chain/chain_prototype.py` (unite each read's SMEMs over `ropebwt3 mem -p`, strict
colinear intersection) vs stock refmap, scored by `eval/score.py`. 8k within-exon
reads/rate, seed 11. **S** = stock, **C** = chaining:

| error | recall S→C | misplaced S→C | spurious/read S→C | resolution S→C | exact-set S→C |
|--:|--:|--:|--:|--:|--:|
| 0.000 | 100→100% | 0.0→0.0% | 0.00→0.00 | 1.000→1.000 | 100→100% |
| 0.005 | 99.7→99.8% | 0.3→**0.0%** | 0.22→**0.02** | 0.927→0.992 | 82→98% |
| 0.010 | 99.4→99.9% | 0.6→**0.0%** | 0.40→**0.06** | 0.869→0.979 | 69→95% |
| 0.020 | 98.8→99.9% | 1.2→**0.0%** | 0.68→**0.13** | 0.791→0.953 | 51→89% |
| 0.030 | 98.7→100%  | 1.3→**0.0%** | 0.90→**0.23** | 0.735→0.920 | 39→81% |

Chaining wins on **every** metric at every error rate: recall equal-or-higher,
**misplaced → 0%** (colinearity + the reference-span bound kill mis-mapping),
**spurious founders 3–9× lower**, resolution/exact-set far higher; founder-dropout
and missed-IBS stay 0. The true founder survives every intersection by
construction, so specificity rises with no recall cost. Junction reads: the same
chainer covers **both** exons at valid bases (unit-tested; r000001 → `{0,1,2,4}` at
both exon starts) — a proper junction sweep needs the whole-read (intersected)
oracle in `score.py`, a small TODO.

Run it: `sh rnaseq_ps4g/run_stage0.sh` now emits both `score.txt` (stock) and
`score.chain.txt` (chaining). Speed: the prototype's `mem -p` is locate-bound, so
`run_stage0.sh` builds the SSA dense (`-s4`); the chaining itself is instant, and a
C `--smem-out` would drop the separate `mem` pass entirely.

## 4. Chaining on junction (spliced) reads

8k reads, `--max-junctions 2` (≈48% within-exon, 49% one-junction, 3% two-junction),
scored against the **whole-read oracle** (intersection of per-segment oracle sets —
a spliced read's founder must match every exon). **S** = stock, **C** = chaining.

Set quality (all reads):

| error | recall S→C | misplaced S→C | spurious/read S→C | resolution S→C |
|--:|--:|--:|--:|--:|
| 0.00 | 84.2→91.4% | 15.8→**0.0%** | 0.14→**0.06** | 0.933→0.952 |
| 0.01 | 87.7→97.0% | 12.3→**0.0%** | 0.52→**0.17** | 0.825→0.923 |
| 0.02 | 90.1→98.8% |  9.9→**0.0%** | 0.80→**0.26** | 0.754→0.900 |

Per-segment coverage (fraction of true exon segments getting a valid emitted base),
at 1% error:

| read class | stock | chaining |
|--|--:|--:|
| within-exon | 58.3% | **100%** |
| one-junction | 8.5% | **67.7%** |
| two-junction | 5.4% | **50.2%** |

Chaining emits one row per exon, so a spliced read covers **all** its exons at valid
bases (stock covers ≤1/N, and its extrapolated base usually lands off the exons —
even within-exon it hits the true base only 58% of the time). Chaining also
converts stock's false **misplaced** (12% at 1% error, from cross-junction
coordinate extrapolation) into honest recall/unmapped, and cuts spurious founders.
Remaining chaining gaps: exon fragments below the 19 bp SMEM floor can't be covered
(so one/two-junction coverage caps below 100%), and a few hard reads go `unmapped`
rather than mis-mapped — the safe direction (lost evidence, not false evidence).

## 5. Adversarial fixtures + intron-aware gap penalty

Fixtures added to the simulator (default off): `sim_genomes --tandem-dup`
(identical adjacent gene copy) and `--large-intron-kb` (one gene with a >2 kb
intron); `sim_rnaseq --chimera` and `--intron-retention` (labeled in a `read_class`
truth column). The chainer's hard `max-ref-span` cap was replaced by a penalty on
the **unexplained reference jump** (`Δref − Δquery`, i.e. the intron length),
folded into the `(#anchors, bases − penalty)` DP, plus a `max-intron` reject and an
**ambiguity flag** (a SMEM with >1 reference occurrence → the read maps to multiple
loci → suppress, don't guess).

Measured on a fixtured genome (6k reads, 1% error):

| fixture | old behaviour | now |
|--|--|--|
| **large intron** (gene5, 5 kb) | rejected by the 2 kb cap → under-served | 155/264 junction reads chain **both exons** across the 5 kb intron |
| **tandem duplication** (gene0≡gene0dup) | placed at an arbitrary copy (false position) | **0/3181** placed — all flagged **ambiguous** (no false evidence) |

Within-exon results are unchanged (penalty ≈ 0 for contiguous reads: recall 99.9%,
misplaced 0.0%, spurious 5.0%).

### Chimera false-evidence (the honest cost of chaining)

Chimeras (two-transcript fusions, `read_class=chimera`) are artifacts from no real
founder, so a chain that **bridges the two loci** fabricates a fusion/linkage that
does not exist. Chaining's uniting behaviour makes this *worse* than stock, and it
is tunable by `--max-intron`. 8k reads, 50% chimera, 1% error:

| | placed (any) | **fabricated fusion** (chain bridges >10 kb) |
|--|--:|--:|
| stock refmap | 93% | **0%** (one position/read — can't span two loci) |
| chaining, `--max-intron 200000` | 99.9% | **18.3%** |
| chaining, `--max-intron 20000` (new default) | 99.9% | **0.0%** |

Fabricated-fusion rate vs `--max-intron`: 0% up to 20 kb, then 15% @100 kb, 18% @200 kb.
Genes here are ≥30 kb apart and the large-intron fixture is 5 kb, so **max-intron in
~[6 kb, 20 kb] rejects every cross-gene fusion while still placing the 5 kb intron** —
the old 200 kb default was simply too loose. The default is now **20 kb**; raise it
for large-intron organisms (at the cost of more chimeric fusions).

`--max-intron` cannot catch a chimera fusing two genes *closer* than the cap; the
real discriminator there is **canonical splice sites (GT–AG)**, which needs real
sequence (the sim uses random bases) — a future filter. The `read_class` label is
what lets us keep this rate honest.

## Conclusion

Both experiments point to the same fix: **stop relying on one core per read.**
Collect all of a read's SMEMs, chain the colinear ones, refine the founder set
across the chain, and emit **per exon segment** at each segment's own start base.
This lowers `spurious` (set constrained by the whole read), keeps `-l` low
(chaining supplies specificity), covers **all** exon segments of spliced reads,
and puts every row at a valid base. Design in `NEXT_STEPS.md`.
