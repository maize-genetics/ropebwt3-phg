# rnaseq_ps4g — RNAseq → PS4G for GRITS imputation

Adapts the WGS RopeBWT3-RefMap pipeline to map RNAseq reads and emit **PS4G**
(the per-reference-position set of pangenome haplotypes/"gametes" a read is
exactly consistent with) for GRITS imputation. See `../HANDOFF.md` for the full
spec.

**Module 0 (this directory, the first deliverable): a synthetic RNAseq read
simulator with exact ground truth, plus a Tier-A scoring harness.** Mapping
sensitivity is unmeasurable without labeled truth; the simulator makes truth
exact by construction so every tuning knob can be scored against known answers.

## Layout

```
sim/
  simlib.py       FASTA/GFF I/O, EditMap (exact genotype<->reference coords),
                  gamete indexing (matches ropebwt3 ps4g.c), oracle sets
  sim_genomes.py  synthetic pangenome: conserved exons + churned introns/
                  intergenic + PAV  ->  FASTAs, GFF, exact coord maps, gamete map
  sim_rnaseq.py   splice mRNA -> expression-weighted fragments -> 150 bp reads
                  + substitution error  ->  FASTQ, per-read truth TSV, answer-key BED
eval/
  score.py        Tier-A: score a refmap/PS4G run vs the oracle (recall, resolution)
tests/
  test_sim.py     divergence targets, coord-map integrity, oracle self-consistency
run_stage0.sh     end-to-end: sim -> index -> stock refmap --ps4g --bin-size 1 -> score
```

## Quick start

```sh
# build the aligner once (repo root)
make omp=0

# whole Stage-0 pipeline into an output dir
sh rnaseq_ps4g/run_stage0.sh /workdir/esb33/stage0 20000 0.01

# tests
python3 -m pytest rnaseq_ps4g/tests/     # or: python3 rnaseq_ps4g/tests/test_sim.py
```

`sim_genomes.py` / `sim_rnaseq.py` / `score.py` are standalone (`--help` each);
write all outputs to an explicit `--outdir`, never into the repo.

## Key design points

- **No binning.** Positions are the exact reference read-start base
  (`refmap --bin-size 1`); GRITS `build_answer_key(bin_size=1)` matches.
- **Gamete index by construction.** A gamete is the sequence name before the
  first `_` (`B73_chr1` -> `B73`); indices are assigned in sorted unique-name
  order — the same rule as `ps4g.c`, so simulator, aligner and GRITS agree.
- **Exact truth.** Each genotype carries an exact block map to the reference
  (`simlib.EditMap`), so read reference coordinates and oracle sets are exact.
- **Mutations in lowercase.** Non-ancestral bases (substitutions, insertions,
  duplications) are lowercased in the genotype FASTAs so they are visible on
  inspection; the reference (ancestral) stays uppercase. `ropebwt3` folds case
  and `read_fasta` uppercases, so this is cosmetic — index, reads and oracle are
  all case-normalized.
- **Oracle set.** The maximal gamete set for a read segment = genotypes whose
  reference projection is identical to the source over the segment's orthologous
  interval (equals the index's exact-match set for colinear exon reads).
- **Scoring from the individual's perspective (HANDOFF 4.4).** Reads come from one
  gamete now (homozygous), a het individual (two gametes) later. `score.py` puts
  every read in exactly one **imputation outcome** (these sum to 100%):
  - `recall` — placed at the correct locus and the true founder is in the set;
  - `founder-dropout` — placed right but the true founder is **missing**;
  - `misplaced` — placed at a locus overlapping none of the true segments
    (evidence at the wrong place);
  - `unmapped` — no confident placement (UNPLACED/MULTI).

  **`founder-dropout` and `misplaced` are the errors that hurt imputation** — the
  true founder loses support, or a wrong locus gains false support. Everything
  else is *secondary specificity*, not a founder-path error.

  For `recall` reads (the true founder is already present), the *other* founders
  in the set are split: `IBS founders` = non-source founders whose sequence is
  **identical over the read span** (identity-by-state — inherent ambiguity, not an
  error); `spurious founders` = founders added only because sequencing error
  shortened the exact match so it also matched them (not truly consistent with the
  full read — a specificity cost the CRF's multihot absorbs; it does **not** drop
  the true founder); `missed-IBS` = IBS founders we failed to emit; plus
  `resolution` and `exact`%. Placement uses interval overlap because refmap
  extrapolates a full-read `[cL,cR)` from the exact core.
- **Strict per-read attribution.** The emitted set for a read is its *exact*
  contribution, read from a **per-read PS4G file** (`refmap --ps4g-per-read`,
  one row per EXACT/PLACED read with its exact gameteSet) — not a lookup into the
  aggregated PS4G, so co-located reads never mix. This makes `missed-IBS` a true
  **0** (a read's core-based set is always a superset of its full-read IBS set)
  and counts every spurious founder (the aggregated lookup under-counted them by
  picking a smaller neighbouring set).

## Stage-0 defaults & result

Milestone 0: B73-only source, substitution-only error, single-end, within-exon
reads (`--max-junctions 0`). Full generative structure (multi-exon splicing,
indels/dups, PAV, other sources, paired-end, adversarial artifacts) is present
behind flags and off by default; later stages flip them on.

Stock `refmap` baseline (8k reads, seed 11), all-reads, vs substitution error
rate. The **imputation-critical** errors — `founder-dropout` and `misplaced` —
stay near zero; `founder-dropout` and `unmapped` are **0.0% at every rate**. What
grows with error is the *secondary* `spurious founders` breadth (strict per-read
attribution, so `missed-IBS` is a true 0.0%):

| error | recall | founder-dropout | misplaced | spurious founders (secondary) | IBS founders |
|------:|-------:|:---------------:|:---------:|------------------------------:|-------------:|
| 0.000 | 100.0% | 0.0% | 0.0% | 0.0%  (0.00/rd) | 1.12 /rd |
| 0.005 |  99.7% | 0.0% | 0.3% | 17.7% (0.22/rd) | 1.12 /rd |
| 0.010 |  99.4% | 0.0% | 0.6% | 31.1% (0.40/rd) | 1.12 /rd |
| 0.020 |  98.8% | 0.0% | 1.2% | 49.4% (0.68/rd) | 1.12 /rd |
| 0.030 |  98.7% | 0.0% | 1.3% | 61.0% (0.90/rd) | 1.13 /rd |

Read this the right way: stock refmap **almost never commits the errors that hurt
imputation** — the true founder is essentially always recovered (`founder-dropout`
≈ 0) and evidence rarely lands at the wrong locus (`misplaced` ≤ 1.3%). The
big-looking `spurious founders` number (up to 61%) is the *benign* one: the true
founder is still in the set, there are just extra founders that error made
consistent over a shortened match — honest ambiguity the CRF's multihot absorbs,
alongside the ~1.1/read of genuine `IBS founders`. A visual report is at
`<outdir>/report.html`.

## Deferred / possible options

- **SMEM length in the PS4G.** Record the exact-match (SMEM/core) length per row
  as PS4G confidence (longer core = more specific, higher-confidence evidence).
  The length is already computed in `search.c` (`s->len` for a whole-read exact
  match; `bestlen`, the longest core, in the error fallback around line 571).
  Two implementation paths, deferred until the Section-5.0 gap report justifies a
  C change (HANDOFF: prefer Python post-processing over touching RopeBWT):
  1. *C column* — add `ev_len` to `rb3_ps4g_acc_t` (`ps4g.h`), thread the length
     through `rb3_ps4g_acc_add`/`refmap_rst_accumulate` (`search.c`), aggregate
     on collapse (max core length per `(pos,gameteSet)` row), and append a
     `smemLen` column in `rb3_ps4g_npy_finalize` (GRITS parses by column name, so
     an extra column is compatible).
  2. *Python post-process* — a thin joiner over `ropebwt3 mem` SMEM output (which
     already reports match lengths) that annotates the PS4G rows, no C change.
