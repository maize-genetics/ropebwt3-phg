# Ref-sensitivity benchmark

**Question.** When we take random short sequences from maize NAM founder genomes
and ask `ropebwt3 refmap` to place them on the B73 reference, what fraction land
at approximately the right B73 coordinate (same chromosome, within ±5 Mb of the
source locus)?

This is the first real maize-scale test of `refmap` — the synthetic locus in
`docs/experiments.md` is tiny, and `docs/handoff.md` ("Open question #1") flags
that speed and placement accuracy on real genomes are unknown.

## Design

- **Index:** 4 NAM genomes (chr1–chr10 only), B73 = reference.
  - B73 `Zm-B73-REFERENCE-NAM-5.0`, B97/Ki3/CML247 `...-NAM-1.0` (MaizeGDB).
  - Sequences renamed `<LINE>_chrN` so `--ref-prefix=B73` selects the reference
    and every genome's chromosomes are uniquely named.
- **Queries:** 100,000 random 150 bp windows, 25,000 per line, no `N`, uniform
  over chromosome length. B73's 25k are a **control** (should be `EXACT`, Δ≈0);
  the 75k from B97/Ki3/CML247 are the cross-genome test. Each query name encodes
  ground truth: `LINE|chr|pos0|len`.
- **Placement:** `refmap --ref-prefix=B73 -t 12` (default `--max-walk=5000`,
  consensus mode).
- **Scoring:** a query is **correct** when refmap's `refName` is on the same
  chromosome as the source and the nearest breakpoint edge (`cL`/`cR`) is within
  ±5 Mb of the source position. Founders are largely colinear with B73, so the
  source coordinate proxies the expected B73 coordinate; ±5 Mb absorbs indel drift.

## Caveats baked into the analysis

- **Mostly `EXACT`.** A random 150 bp founder window is usually also in B73, so
  refmap returns `EXACT` (B73's direct coordinate); only insertion-borne windows
  exercise the `PLACED` outward-walk. The report splits accuracy by status.
- **Maize is ~85% TEs.** A 150 bp query can hit many B73 loci; refmap reports
  one. Some "wrong" placements are repeat multi-mapping, not algorithm error —
  the B73 control quantifies that floor.
- **Colinearity drift** can push a genuinely-correct placement past ±5 Mb near
  chromosome ends, so the full |Δpos| distribution is reported, not just pass/fail.

## Reproduce

From the repository root, after `make omp=0`:

```sh
sh experiments/ref-sensitivity/run.sh
```

This downloads (~2.6 GB), builds the index, generates queries, runs refmap, and
prints the report (also saved to `data/report.txt`). Large artifacts land in
`experiments/ref-sensitivity/data/` (gitignored). Needs ~12 GB free disk and a
few GB RAM; the index build and the 100k-query run are the slow steps.

## Files

| File | Purpose |
|---|---|
| `run.sh` | end-to-end pipeline (download → index → query → score) |
| `gen_queries.py` | draw N random fixed-length queries per genome + truth table |
| `analyze.py` | join refmap output to truth; ±5 Mb buckets + |Δ| distribution |
| `results.md` | measured results, timings, and interpretation |
| `data/` | genomes, index, queries, refmap output (gitignored) |
