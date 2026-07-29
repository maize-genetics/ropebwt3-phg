# Benchmark — `refmap` at maize-pangenome scale

The synthetic locus in [`experiments.md`](experiments.md) proves correctness on a
known breakpoint. This benchmark measures **sensitivity** and **speed** on real
data: four NAM maize founders indexed together, with 100k synthetic-but-grounded
reads whose true origin is known. Results (and a running scoreboard across
variants) live in [`results-maize.md`](results-maize.md).

Driver script: [`../experiments/ref-sensitivity/server-run.sh`](../experiments/ref-sensitivity/server-run.sh).

## Data

Four NAM founder assemblies from MaizeGDB, chr1–chr10 only, renamed `<LINE>_chrN`:

| Line | role | chr1–10 size (fwd) |
|---|---|---|
| B73 | **reference** (`--ref-prefix=B73`) | 2.13 Gbp |
| B97 | carrier | ~2.14 Gbp |
| Ki3 | carrier | ~2.14 Gbp |
| CML247 | carrier | ~2.14 Gbp |

Total input: **~8.5 Gbp forward → 17.10 G both-strand symbols** (per-genome
symbol counts from `build.log`: 4.264 + 4.270 + 4.279 + 4.286 G).

## Index build (one big-RAM box, 20 threads)

```sh
NCORES=20 MEM_BATCH=60G sh experiments/ref-sensitivity/server-run.sh
```

| Step | Command (key flags) | Wall | Peak RAM | Output |
|---|---|---:|---:|---|
| build FMD | `build -d -t20 -p20 -m60G` | 32:23 | 38.8 GB | `nam4.fmd` 2.13 GB |
| sampled SA | `ssa -t20 -s8` | 9:24 | 3.0 GB | `nam4.fmd.ssa` 510 MiB |
| len table | (awk over the FASTAs) | <1 min | — | `nam4.fmd.len.gz` (40 seqs) |

Notes:
- `-m 60G` exceeds the 17.1 G both-strand input, so libsais runs per-file and
  merges across the 4 files (4 batches); RAM peak 38.8 GB ≪ box.
- `ssa` is **kthread-parallel** (`-t`), not OpenMP. At the default `-t4` it was
  the pipeline long pole (>40 min); `-t20` cuts it to 9:24. Size is a flat
  8 bytes × (symbols / 2⁸) ≈ 510 MiB, deterministic from `-s8`.

## Queries (deterministic)

```sh
gen_queries.py --n 25000 --len 150 --seed 7 \
  B73:… B97:… Ki3:… CML247:…        # --n is PER GENOME
```

→ **100,000 reads** (25k each × 4 lines), 150 bp, no `N`. The read name encodes
ground truth `LINE|chr|pos0|len`; `truth.tsv` is the join key.

## Metrics

### Sensitivity (`analyze.py`, tol = ±5 Mb)

A placement is **correct** when refmap's reference chromosome equals the read's
source chromosome **and** `min(|cL−src|,|cR−src|) ≤ tol`. NAM founders are
largely collinear with B73, so the source coordinate is a *proxy* for the
expected B73 coordinate and the tolerance absorbs indel drift.

- `%correct` — correct / reads (per line and overall).
- `wrong_chr` — placed on the wrong chromosome (**genuine error**; drift cannot
  move a read to another chromosome).
- `far(>tol)` — right chromosome, > tol away (partly real inter-genome drift, so
  this bucket slightly *under*-counts true accuracy).
- **precision** = correct / (reads − UNPLACED); **recall/yield** = placed / reads.
  For genotyping, precision is the figure of merit.

Report also stratifies by **status** (EXACT / PLACED / ONE_SIDE / UNPLACED) and
by **read multiplicity** (occurrence count from `ropebwt3 mem` col 4 = SA-interval
size), which is the single best predictor of correctness (see results).

### Speed (`/usr/bin/time -v`, `-t 20`, same 100k reads + index)

Compare `refmap` against the regular mappers `mem` (SMEM+locate) and `sw` (local
align). Record wall, reads/sec, CPU% (scaling), CPU-seconds, peak RSS. `refmap`'s
slowdown over `mem`/`sw` quantifies the cost of the walk + flank-anchoring.

## Experiment matrix

`E0` is the current baseline (measured). `E1–E3` test the precision fixes
motivated by the baseline diagnosis (see [`results-maize.md`](results-maize.md)),
each behind a flag, re-run on the **same** index + 100k reads.

| ID | Variant | Flag | Hypothesis | Result |
|---|---|---|---|---|
| **E0** | baseline | — | precision low, `ONE_SIDE` dominates error | prec 57.3% |
| **E1** | two-flank concordance | `--two-flank` (+ `--max-bracket`) | kills the bucket that is ~60% of errors; precision ↑, yield ↓ | **prec 74.7%**, wrong-chr −60% ✓ |
| **E2** | uniqueness filter | `--max-occ=N` (`-1`=auto=#taxa) | removes retro/repeat mis-placements; precision ↑ at low recall cost | **prec 88.1%**, wrong-chr −85%, ~free recall ✓ |
| **E3** | E1 + E2 | both | best precision; measure recall + speed | **prec 96.0%**, wrong-chr −97% ✓ |
| **E4** | synteny prior (planned) | _tbd_ | use the carrier coordinate the SSA already returns + a coarse per-carrier→B73 offset map to skip/shorten the walk and reject off-diagonal carriers | speed ↑↑ (attacks the walk) and precision ↑ |

Implemented flags (E1–E3), all default-off so E0 is reproducible:
- `--two-flank` — require both flanks to anchor concordantly (same chr + strand);
  a lone anchor yields `UNPLACED` instead of `ONE_SIDE`.
- `--max-bracket=NUM` — with `--two-flank`, reject a `PLACED` whose two anchors
  bracket > NUM bp (collinearity guard; 0 = off).
- `--max-occ=N` — drop reads (status `MULTI`) and flank anchors occurring > N
  times; `N<0` = auto = number of taxa (distinct name prefixes); 0 = off.

E4 is motivated by the speed analysis: E1–E3 are precision fixes and leave the
walk (the 35× cost) intact. The SSA returns the exact *carrier* coordinate for
free at the locate step; projecting it through a small precomputed
carrier→reference offset map gives an approximate B73 band without walking, and
simultaneously enables an off-diagonal (paralog) rejection — see
[`handoff.md`](handoff.md) open-question #2.

**E4 is now implemented** as the `lift` subcommand (the "second SSA") + `refmap
--lift` (`lift.c`, `search.c`). Workflow:

```sh
# build the carrier->reference liftover once (like ssa), from the reference FASTA
ropebwt3 lift --ref-prefix=B73 -k 100 -s 2000 -t 20 -o nam4.lift nam4.fmd B73.chr.fa.gz
# place reads by projection instead of walking
ropebwt3 refmap --ref-prefix=B73 --max-occ=-1 --lift nam4.lift -t 20 nam4.fmd queries.fa
```

`lift` slides k-mers over the reference, and for each single-copy-per-taxon k-mer
(≤ `-m`, auto = #taxa) records (carrier_pos ↔ reference_pos) pairs. At query time
`refmap --lift` projects a carrier hit with a windowed slope=±1 robust fit
(`--lift-win`, `--lift-mad`), returning `PLACED` or, where no confident collinear
support exists, leaving the read `UNPLACED` (the NULL slot). Measured: **95.5%
precision, 77.9% recall, 2.76 s / 100k reads** (see results-maize.md E4 row).
How the liftover is constructed, step by step: [`lift-second-ssa.md`](lift-second-ssa.md).

**Biological basis for E2:** ~40% of sequence is not shared between two maize
lines, but retrotransposons are shared and high-copy. A single-copy (informative)
locus contributes at most one hit per taxon, so an informative read maps
**< (number of taxa)** times; a retro maps ≫ taxa. Occurrence count therefore
separates informative loci from repeats, and the threshold scales with the panel
(≈4 here, ≈24 for the full NAM pangenome).

Success criteria per variant: report **precision, wrong-chr rate, yield**, and
**reads/sec**; a variant "wins" if it raises precision materially without
collapsing yield, at acceptable speed.

## Reproduce

```sh
# 0. build tool + index + ssa + len + queries  (≈45 min on 20 cores)
NCORES=20 MEM_BATCH=60G sh experiments/ref-sensitivity/server-run.sh

# 1. baseline sensitivity already emitted as data/report.txt; speed:
D=experiments/ref-sensitivity/data
for cmd in "refmap --ref-prefix=B73" mem sw; do
  /usr/bin/time -v ./ropebwt3 $cmd -t 20 $D/nam4.fmd $D/queries.fa \
    > $D/${cmd%% *}.out 2> $D/${cmd%% *}.time
done
python3 experiments/ref-sensitivity/analyze.py \
  --refmap $D/refmap.out --truth $D/truth.tsv --tol 5000000
```
