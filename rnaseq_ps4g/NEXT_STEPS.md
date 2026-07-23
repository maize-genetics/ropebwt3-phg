# Next steps: uniting information across a read's SMEMs (chaining / splice segmentation)

Status: plan for review (Ed asked to start this while the min_len sweep ran).
Motivating measurements live in `RESULTS.md`.

## Why

An RNAseq read fragments into several exact matches (SMEMs) for two reasons that
look identical to the index: **sequencing errors** and **splice junctions**. Both
break the exact match. Stock `refmap` keeps only the *single longest* SMEM (its
"core") and emits one placement per read. That throws away everything the read's
*other* SMEMs know, and it creates the two failure modes the Stage-0 metrics
already isolate:

1. **misplaced + spurious founders** — a lone short core is not specific: it
   matches paralogous loci (misplacement) and matches extra founders over its
   shortened span (spurious founders). Raising the min-SMEM-length floor (`-l`)
   buys specificity but only by dropping reads to `unmapped` (see the RESULTS
   sweep) — you cannot win both with a single core.
2. **lost segment coverage on junction reads** — a spliced read covers N exon
   segments; the longest core is in one exon, so only ~1/N of the read's exonic
   evidence is emitted (RESULTS: junction reads get ~1/N segment coverage). The
   other exons — and every read that spans a junction — are under-served.

**Uniting the SMEMs fixes all three at once**, and decouples the min_len tradeoff:
you can keep short SMEMs (high coverage) because chaining supplies the specificity
a single short core lacks.

## The idea (the whiteboard "intersection of hits")

For each read, instead of one core:

1. **Collect all SMEMs** (both strands). Per SMEM: query span `[qs,qe)`, its SA
   interval → the set of gametes carrying it, and each occurrence's reference
   coordinate (RefMap/lift projection).
2. **Chain** SMEMs whose reference projections are **colinear**: same contig, same
   strand, monotonic in both query and reference, gaps consistent with introns.
   Each colinear cluster is one exon block; a chain across an intron-sized gap is a
   spliced read; non-colinear/chimeric sets are flagged, not emitted as one locus.
3. **Refine the gamete set along the chain.** The true founder is exact over the
   *union* of the read's cores, so combine the per-SMEM gamete sets down the chain
   (intersection at colinear positions). A founder that matches one short SMEM but
   not the adjacent one is dropped → **spurious founders fall**, while the true
   founder (consistent everywhere) and genuine IBS founders survive.
4. **Emit per segment.** Each exon cluster emits its own PS4G rows at its own exact
   reference start base; intronic positions get nothing → **junction reads cover
   all N segments**, zero intronic rows (HANDOFF 2.2/2.3).

### Which metric each win moves
| win | mechanism | metric |
|--|--|--|
| disambiguation | colinear chain confirms the locus, rejects off-target short-SMEM hits | `misplaced` ↓; lets `-l` stay low without `unmapped` ↑ |
| set refinement | founder must be exact across the chain, not one short core | `spurious founders` ↓ (toward the IBS-only floor) |
| per-segment emission | every colinear cluster emits its own rows | junction `segment coverage` → ~100% |
| (preserved) | the true founder is exact over every core it spans | `recall` maintained, `founder-dropout` stays ~0 |

## Implementation (prefer Python post-processing; minimal C)

We already emit a per-read PS4G file (`--ps4g-per-read`) with the *one* chosen
placement. Chaining needs **all SMEMs per read** with their gamete sets and
reference projections. Two build options:

- **A — Python prototype over `ropebwt3 mem` (fastest to stand up).** `mem` already
  emits every SMEM (query span + SA interval). Reuse the existing E4 projection
  (`experiments/ref-sensitivity/e4/liftover.py: project_robust`) to place each
  SMEM occurrence on the reference, map SA interval → gametes via `sid2g`, then
  chain + refine + emit in Python. No C change; leans on code already in the repo.
- **B — extend the C per-read PS4G file to per-SMEM (cleaner, later).** Add
  `refmap --smem-out=FILE` emitting one row per (read, SMEM): `read, qStart, qEnd,
  matchLen, refContig, refPos, gameteSet`. Then the Python chainer consumes that
  directly (no separate `mem` run / re-projection). This is the natural successor
  to `--ps4g-per-read` and reuses the same `write_ps4g_read` plumbing in `search.c`.

Recommended: build **A** first to validate the algorithm and the metric gains on
the simulator, then promote the hot path to **B**.

### New module layout
```
rnaseq_ps4g/chain/
  collect.py    # read SMEMs (mem output or --smem-out) -> per-read SMEM list
  project.py    # SMEM occurrences -> reference coords (reuse e4 liftover)
  chain.py      # colinear clustering + intron-gap segmentation + chimera flag
  refine.py     # combine per-SMEM gamete sets along a chain
  emit.py       # per-segment, exact-base PS4G rows (+ per-read PS4G file)
```
`eval/score.py` and `eval/analyze_junctions.py` are already the scorers; the new
emitter just replaces stock refmap's PS4G/per-read file in `run_stage0.sh`.

## Validation (same simulator, honest deltas)
Re-run the RESULTS sweeps against the chaining emitter and expect, vs the stock
baseline: `misplaced` ↓, `spurious founders` ↓ (toward IBS-only), junction
`segment coverage` → ~100%, `recall` unchanged / `founder-dropout` ~0, and the
recall-vs-`-l` curve flattened (short SMEMs safe). Unit fixtures (HANDOFF 6):
within-exon, one-junction, two-junction (small middle exon → 3 SMEMs), PAV-gene,
tandem-duplicate, highly-divergent, chimera, intron-retention.

## Open questions for Ed
- Intron-gap bounds for "same-strand colinear = spliced": fixed cap vs annotation
  (junction-augmented index, HANDOFF decision 5)?
- Set-refinement operator: strict intersection across chain SMEMs, or a
  vote/consistency threshold (robust to one bad SMEM)?
- Do we want per-SMEM emission (option B) now, or is the `mem`-based prototype
  enough to set the min_len floor and quantify the chaining gain first?

## Backlog

- **Remove `--kmer`** (Ed, longer-term). The `--kmer` mode (`refmap_query_kmer` /
  `refmap_kmer_votes` in `search.c`, options `--kmer/--kmer-step/--min-agree/
  --kmer-cluster`, and the MAPQ calibration) is a separate positional k-mer-
  agreement *placement* path that emits confidence columns but **no PS4G/npy**.
  Set-based colinear chaining supersedes its purpose (multi-seed agreement for
  specificity) *and* produces refined PS4G founder sets. Once chaining lands and
  is validated, delete `--kmer` and its options to cut surface area.

## Progress (prototype, this session)

- `chain/chain_prototype.py`: SMEM collect from `mem -p` → in-order colinear DP
  chain scored by (anchor-count, bases) → spatial-compactness bound
  (`--max-ref-span`, rejects distant chimeras) → strict whole-read set
  intersection → per-exon-segment emission (per-read PS4G file schema).
  Handles negative strand and junctions; excludes carrier-only/chimeric SMEMs.
- `sim/sim_rnaseq.py`: added `--library sense|antisense|unstranded`.
- `tests/test_chain.py`: 5 unit tests (within-exon ±, junction linkage, chimera
  exclusion, out-of-order rejection). Verified on real reads r000000/1/2:
  r000002 emits `{0}` (spurious W22 removed vs stock `{0,4}`); junction r000001
  covers both exons with the correct intersection `{0,1,2,4}`.
- TODO: score.py per-read PS4G file-only + multi-row-per-read (junction per-segment
  scoring); carrier-only SMEM lift projection (non-B73 sources / PAV);
  extrapolate segment start to the read-start base; then the RESULTS.md sweeps.
