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
| **large intron** (gene5, 5 kb) | rejected by a hard 2 kb cap | placed across the 5 kb intron **when `--max-intron ≥ 5 kb`** (155/264 both exons); at the conservative 500 bp default it is intentionally not linked (GT–AG makes raising the cap chimera-safe) |
| **tandem duplication** (gene0≡gene0dup) | placed at an arbitrary copy (false position) | **0/3181** placed — all flagged **ambiguous** (no false evidence) |

Within-exon results are unchanged (penalty ≈ 0 for contiguous reads: recall 99.9%,
misplaced 0.0%, spurious 5.0%).

### Chimera false-evidence (the honest cost of chaining)

Chimeras (two-transcript fusions, `read_class=chimera`) are artifacts from no real
founder, so a chain that **bridges the two loci** fabricates a fusion/linkage that
does not exist. Chaining's uniting behaviour makes this *worse* than stock. Two
defenses gate it: a `--max-intron` cap on the reference jump, and a **GT–AG
splice-site check**. 8k reads, 50% chimera, 1% error:

| | placed (any) | **fabricated fusion** (chain bridges >10 kb) |
|--|--:|--:|
| stock refmap | 93% | **0%** (one position/read — can't span two loci) |
| chaining, `--max-intron 200000`, no GT–AG | 99.9% | **18.3%** |
| chaining, `--max-intron 500` (default) | 99.9% | **0.0%** |

Genes here are ≥30 kb apart, so a tight cap alone kills every cross-gene fusion; the
default is **500 bp** (`--max-intron`) — conservative, trading large-intron *power*
for chimera *specificity* per the operating choice. Fabricated fusion vs cap (no
GT–AG): 0% ≤20 kb, 15% @100 kb, 18% @200 kb.

### GT–AG splice-site validation

`sim_genomes` now writes **canonical splice motifs** into every reference intron
(GT‥AG in transcription orientation — genomic `GT‥AG` for `+` genes, `CT‥AC` for
`−`; 17/17 introns verified). The chainer (`--ref-fasta`) then requires an
intron-sized chain link to have canonical boundaries, searching a small ±`slack`
window because a SMEM end drifts a few bp from the true splice site (microhomology /
error near the junction).

The payoff is that GT–AG lets you **raise `--max-intron`** to place large real
introns while still rejecting far chimeras — a lever the cap alone can't give (the
cap trades away *all* large introns). At `--max-intron 200000`, GT–AG cuts fabricated
fusion from 18.3% toward zero; `slack` trades junction recovery against chimera
leakage (junction "both-exon" baseline without GT–AG is 45.8%):

| `slack` | fabricated fusion | normal-junction both-exon |
|--:|--:|--:|
| 2 | 2.8% | 31.7% |
| 6 (default) | 8.8% | 38.5% |
| 12 | 13.8% | 43.8% |

At the **tight default cap (500 bp)** GT–AG is a *redundant* second gate on
chimeras (the cap already gives 0% fusion) — its value is unlocked when the cap is
raised. Crucially, GT–AG only ever rejects the intron **link** (the read still
places at one exon; the second exon is dropped → benign `missed-IBS`/coverage),
so it **never** causes founder-dropout or misplacement — the safe direction.

Scored at the final defaults (`--max-intron 500`, GT–AG on, canonical genome, 1%
error), the imputation-critical metrics stay perfect:

| set | recall | founder-dropout | misplaced | spurious | missed-IBS |
|--|--:|--:|--:|--:|--:|
| within-exon (8k) | 99.9% | **0.0%** | **0.0%** | 6.5% | 0.0% |
| junction, `--max-junctions 2` (8k) | 97.0% | **0.0%** | **0.0%** | 16.3% | 7.4% |

## 6. Native C chaining (`ropebwt3 chain`)

The prototype chainer is now implemented natively in `search.c` as a new subcommand
`ropebwt3 chain`, reusing the SMEM+locate machinery (`rb3_fmd_smem_TG` +
`rb3_ssa_multi`) and the `gtab` gamete table — no separate `mem` process, no Python.
It runs the same colinear in-order DP + strict intersection + per-exon-segment
emission + `max-intron` cap + ambiguity flag as `chain_prototype.py`, emitting the
identical per-read PS4G (`readName refContig refPos gameteSet`).

**Accuracy** — byte-identical to the Python prototype on 5k mixed junction reads
(**6453/6453 rows match, 0 diffs**) at matched parameters (`-l 19`, `--max-occ 5`,
`--max-intron 500`, `--gap-intron 30`). Scored on 100k reads: recall 97.6%,
**founder-dropout 0.0%, misplaced 0.0%**, unmapped 2.4% — matching the Python metrics.

**Speed** — 100k RNAseq reads, 16 threads, single-chr pangenome; compute time
(`worker_pipeline`, excludes index load), best of 3:

| tool | time | vs refmap |
|--|--:|--:|
| `ropebwt3 chain` (native) | **0.49 s** | 0.98× |
| `ropebwt3 refmap --lift --ps4g` (stock) | 0.50 s | 1.00× |
| `mem -p64 \| chain_prototype.py` (Python) | 1.33 s | 2.7× |

Native chaining is **as fast as stock refmap** — uniting a read's SMEMs adds no
measurable cost over placing one — and ~2.7× faster than the Python pipeline, while
emitting the richer united, per-segment PS4G. Run it:
`ropebwt3 chain --ref-prefix B73 idx.fmd reads.fq > reads.ps4g`.

**Not yet in C**: the GT–AG splice check (needs reference-sequence access at the
intron boundaries); `ropebwt3 chain` equals `chain_prototype.py` *without*
`--ref-fasta`. Adding it (load the reference contigs, look up the boundary motifs
with the same ±slack search) is the next step.

## Conclusion

Both experiments point to the same fix: **stop relying on one core per read.**
Collect all of a read's SMEMs, chain the colinear ones, refine the founder set
across the chain, and emit **per exon segment** at each segment's own start base.
This lowers `spurious` (set constrained by the whole read), keeps `-l` low
(chaining supplies specificity), covers **all** exon segments of spliced reads,
and puts every row at a valid base. Design in `NEXT_STEPS.md`.
