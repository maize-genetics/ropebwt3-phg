# asmqc: fast assembly-vs-reference QC

`ropebwt3 asmqc` and `ropebwt3 asmqc-fa` compare genome assemblies with a reference (e.g. B73) and report
the same kind of numbers PHGv2 computes from an AnchorWave alignment — mapped / identical / deleted /
inserted bp, SNPs, indels — plus structural flags (low coverage, inversions, sequence on other chromosomes
or unplaced scaffolds). The intended use is to **screen candidate assemblies before building a pangenome
index**: run every candidate against the reference, have someone review the flags, and build the index from
the lines that pass.

| | `asmqc-fa` | AnchorWave + gVCF (PHGv2 `align-assemblies`) |
|---|---|---|
| Input | reference FASTA + assembly FASTA(s) | same |
| Time per maize line, 24 threads | ~2.3 min | 90–108 min |
| Peak memory | ~7.5 GB | ~278 GB |

## Two modes

**`asmqc-fa` — from FASTA (recommended for screening).** No index is built. The reference's unique
k-mers at a fixed stride are tabulated once per run; each assembly is scanned for them, anchors are chained,
and the sequence between anchors is aligned.

```
ropebwt3 asmqc-fa -t 24 [--chains chains.tsv] [--variants variants.tsv] B73.fa line1.fa line2.fa ... > qc.tsv
```

Line names come from the FASTA file names (text before the first `.`, `_renamed` dropped). Sequence
names may be `chr1`, `chr01`, `Chr1` or `<line>_chr1`; chromosomes are matched by number.

**`asmqc` — from an existing panel index's liftover map (anchors only).** Uses the `.lift` built by
`ropebwt3 lift`; gives the structural metrics in seconds for every assembly already in an index, but no
base-level numbers (no SNPs).

```
ropebwt3 asmqc --ref-prefix=B73 -k 100 -s 200 -t 16 index.lift index.fmd.len.gz > qc.tsv
```

`-k` and `-s` must match the values the `.lift` was built with (they are not stored in the file).

## Method

1. **Anchors.** Reference k-mers (`-k`, default 100) sampled every `-s` bp (200 in `asmqc-fa`) that occur
   exactly once in the reference and once in the assembly (both strands, verified base by base).
2. **Chaining.** Per assembly sequence × reference sequence, in two stages: blocks along a near-constant
   diagonal (`--block-diag`, 10 kb), then blocks joined across gaps up to `--max-gap` (3 Mb, for
   AnchorWave-scale deletions) or `--max-gap-off` (1 Mb) for chains to other chromosomes / scaffolds.
   Chains to the matching chromosome are syntenic; the minority orientation is reported as inverted.
3. **Base level (`asmqc-fa` only).** Between consecutive anchors, reference and assembly sequence are
   aligned with **AnchorWave's gap model** — match 0, mismatch 4, 2-piece affine gaps min(4+2l, 80+l) —
   using ksw2 (vendored from minimap2; see `KSW2_NOTICE.txt`). ksw2 needs a positive match score, so the
   penalties are converted to an exactly equivalent scoring (+2 / −2 / open 4 extend 1 / open 80 extend 0).
   Each gap is first extended from both flanking anchors until similarity collapses (these breakpoints act
   as extra anchors); the pieces are re-aligned, and the middle is a clean indel, is split on shared
   k-mers, or is whichever of "align" vs "delete reference + insert assembly" scores better under the same
   model. Like AnchorWave, every reference base is aligned at most once: anchors and balanced gaps first
   (inverted chains, then forward chains by size), lopsided gaps second.

## Output (`qc.tsv`)

One row per assembly × reference sequence plus an `ALL` row per assembly.

| Column | Meaning |
|---|---|
| `frac_mapped`, `mapped_bp`, `gap_bp` | reference inside syntenic chains (≈ PHG `percentMappedToRef`) / outside |
| `inv_bp`, `other_chr_bp`, `scaf_bp` | reference covered by inverted chains / chains to another chromosome / to unplaced scaffolds |
| `numSNPs`, `numIndels`, `numIns`, `numDel` | AnchorWave-named base-level counts (`asmqc-fa`) |
| `numBasesInserted`, `numBasesDeleted` | inserted / deleted bp |
| `percentRefAligned`, `percentRefDeleted`, `percentRefInserted` | fractions of the reference length; identical bp = `percentRefAligned`·length − `numSNPs` (≈ PHG `percentIdentityWithRef`) |
| `mean/lowerQuartile/median/upperQuartileIndelSize`, `largestIndel` | indel size distribution |
| `numReplacements` | segments with no alignment, counted as one deletion + one insertion |
| `numSNPs_localMismatch10/40` | SNPs with ≥10 / ≥40 other SNPs within ±50 bp |
| `ins_bp_lb`, `del_bp_lb`, `tight_bp`, `div_est`, … | anchor-only estimates (also in lift mode) |

`--chains FILE`: every kept chain — reference span, assembly span, strand, anchor count, class
(`syn`, `syn_inv`, `other_chr`, `scaffold`). `--variants FILE` (`asmqc-fa`): every SNP, INS, DEL and REP
(replacement) with its 1-based reference position.

## Review workflow

```
scripts/asmqc_summary.py --qc qc.tsv --chains chains.tsv [--anchorwave metrics.tsv] --out-tsv summary.tsv --out-json summary.json
scripts/asmqc_vs_anchorwave.py --qc qc.tsv --anchorwave metrics.tsv
```

`asmqc_summary.py` gives one row per line with **check** flags (a chromosome < 90% mapped; ≥ 10 Mb
inverted, on other chromosomes or on scaffolds; genome < 97% mapped) and **note** flags (90–95% mapped,
5–10 Mb of the above, unusually high divergence). An inversion whose breakpoints recur (±1 Mb) in another
line of the batch is reported as a shared inversion — a note, not a check — since independent assemblies do
not share errors at identical breakpoints. Screen lines together so these can be recognised.
`--anchorwave` accepts `phg calc-vcf-metrics` output (or the older PHG metrics CSV) for side-by-side
review.

## Agreement with AnchorWave

24 NAM founders, chr1–10, against AnchorWave gVCFs of the **same assembly files** scored with
`phg calc-vcf-metrics` (ratios are asmqc / AnchorWave):

| Metric | ratio (range over lines) | r across lines | r across 240 chromosome rows |
|---|---|---|---|
| Mapped bp | 1.00 (0.99–1.01) | — (AnchorWave ≈ 99% for every line) | 0.998 |
| Identical bp | 1.02 (1.00–1.03) | 0.984 | 0.997 |
| Inserted bp | 1.02 (1.00–1.03) | 0.991 | 0.994 |
| Deleted bp | 1.01 (0.99–1.02) | 0.991 | 0.995 |
| Insertion + deletion events | 0.84 (0.82–0.86) | 0.984 | 0.990 |
| SNPs | 0.75 (0.72–0.76) | 0.949 | 0.969 |

**Why SNPs are ~0.75×.** On Oh43 chr10, ~90% of the AnchorWave SNPs asmqc does not call lie in reference
sequence asmqc reports as deleted or replaced. Mapping those reference segments to the whole Oh43 genome
(minimap2, evaluation only): for 87% of those SNPs the sequence's ≥90%-identity copy is elsewhere in
Oh43 (> 200 kb away or another chromosome) — it is a TE-type sequence absent at that locus, which
AnchorWave aligns to different nearby sequence anyway (often a diverged TE copy, sometimes unrelated
sequence at ~70% mismatch). Only ~2% are present at the same locus and missed by asmqc. asmqc's SNP count
is therefore a "same-locus" count; AnchorWave's also includes differences between different copies placed
side by side. The ratio is stable across lines, so lines rank the same.

Do not compare against PHG metrics files produced from other AnchorWave runs: a fresh AnchorWave run on
B97 gives 22.4 M SNPs where an older metrics CSV for B97 gives 18.9 M.

## Experimental options (off by default)

Each was tested on Oh43 chr10 against AnchorWave and left off because it did not help; enabled by
environment variable.

| Variable | What it does | Result |
|---|---|---|
| `QA_CHAIN=1` (`QA_CK`, `QA_COCC`, `QA_CMIN`) | seed gaps with all shared k-mers incl. repeats, best minimap2-style chain | no more shared SNPs at k = 11–19 |
| `QA_RESCUE=1` | before calling a replacement, align the best diagonal found by 11-mer voting | +0.1% shared SNPs |
| `QA_SOFT=1` (`=2` prints stats) | drop an anchor when aligning across it has lower AnchorWave penalty | 47 of 291k anchors dropped; +0.04% |

## Tests

`test/run_asmqc.sh` (part of `make test`) builds a synthetic reference/assembly pair
(`test/make_asmqc_fixture.py`) with planted SNPs, small indels, a replacement, a 5 kb insertion, a
20 kb deletion, an inversion, a translocation and a solo LTR, and checks both modes and the summary script.
