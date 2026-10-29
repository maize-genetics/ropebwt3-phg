# ropebwt3 for GRITS: what was added and how it works

A plain-language summary of everything in this branch. Each section says what a feature is for, how it works in a
few sentences, and which flags turn it on. Every new option is **off by default**: with none of them set, output
is byte-identical to what it was before. The detailed design notes, test results and measurements are in the
`docs/*/PLAN.md` files linked at the end; you do not need them to use the tool.

## The big picture

GRITS imputes a sample's genome as a mosaic of known founder genomes. ropebwt3 does the first step: for each
sequencing read it works out (1) where the read sits on the B73 reference and (2) which founders the read matches.
`refmap` writes this as a table of reference positions × founders (`--ps4g`, `--npy`), which the GRITS model reads.

All founders plus B73 live in one FM-index. A **lift** file (`ropebwt3 lift`) records "anchors": short sequences
that occur once in B73 and at a matching place in a founder. Anchors are what let a position in a founder be
translated into a B73 position.

Every read gets a status:

- **EXACT**: the read matches B73 itself.
- **PLACED**: the read matches only founders, so its B73 position comes from translating the founder position
  through the lift.
- **MULTI**: too many copies (a repeat); not used.
- **FILTERED**: removed by one of the optional filters below; not used.

Only EXACT and PLACED reads go into the GRITS input.

## 1. Better placement of founder-only reads (the "PLACED fix")

**Problem.** A PLACED read was translated to B73 using only the two anchors right next to it. One bad anchor (for
example a repeat that happens to look unique in B73) broke the translation, and the read fell back to a coarse
estimate that was often tens of kilobases off. About half of PLACED reads ended up more than 20 kb away from where
they belong.

**What it does now.** It uses only anchors that agree with their neighbours (a consistent "chain"), and when two
good anchors disagree because of an insertion or deletion between them, it interpolates between them instead of
guessing from one side. Misplaced PLACED reads drop from about 45% to about 27%.

**Flags:** `--lift-local-win=20000 --lift-chain-nb=3 --lift-chain-interp`, with a lift built at `-k 50 -s 200`.

## 2. Paired-end reads (`--paired`)

**What it does.** It reads R1/R2 files (or one interleaved file) and treats the two reads of a pair as one
fragment. If one mate is placed confidently and the other lands far away, the far mate is moved next to its
partner when one of its founder copies is there. Each pair adds **one** entry to the output, credited to the
founders that both mates support. Works in `refmap`, `chain`, `mem` and `sw`.

**Note for GRITS.** In our tests paired input helped in-panel hybrids but did not improve the model on genomes
outside the panel, so the current GRITS pipeline uses single-end mode. Paired mode is there for those who want it.

## 3. Reads with too many differences are removed (the edit filter)

**Problem.** A read that only partly matches every genome gets placed by its longest exact piece, and those reads
are often placed wrongly. Sequencing errors alone almost never put more than one difference in a 150 bp read.

**What it does.** For each read it walks base by base from the exact piece out to both ends, against one genome,
and counts every difference: a wrong base costs 1, a k-base insertion or deletion costs k, and every base hanging
past the end of the genome costs 1 (nothing is clipped for free). Reads with more differences than the limit are
FILTERED.

**Flags:** `--err-noclip --err-max-events=1`, which keeps reads with at most one difference.
`--err-max-events` sets the limit. `--err-noclip` switches on the counting described above. Without it,
an older, looser count is used: everything hanging past a contig end counts as a single difference
however long it is, and an insertion or deletion longer than 1 bp is not recognised as one event.
Optional extras:
`--err-stats` adds the counts to each output row, `--err-cigar` writes the read-vs-genome alignment (for example
`75=1X74=`), `--min-match-frac=0.99` filters on the fraction of exactly matching bases, `--err-hist=FILE` writes a
histogram of differences for the whole run.

## 4. Reads whose matches disagree are removed (preset E)

**Problem.** Some reads match in several places, or their founder copies point to different places on B73.

**What it does.** It drops a read if any of these is true:

- **its exact match is too short** to trust (`--min-core-len=100`);
- **EXACT read with more than one copy in B73** (`--exact-max-ref-hits=1`). Before, one copy was picked;
- **EXACT read whose founder copies point somewhere else on B73,** more than 100 kb away (`--exact-min-asm-agree=1`);
- **PLACED read whose founder copies land on different chromosomes or more than 1 Mb apart**
  (`--placed-min-chr-agree=1 --placed-max-spread=1000000`).

`--filter-stats` adds the numbers behind each decision to the output.

## 5. Edit-distance statistics (`--ed-stats`, `--ed-hist`)

**What it does.** For each read it reports how many differences separate it from **each** founder genome, plus the
best over all founders. `--ed-hist` writes summary distributions (best founder, every founder, closest wrong
founder) that help choose a filter cutoff. It is slow on every read, so use `--ed-frac=0.1` to measure a 10%
sample; the distributions come out the same.

## 6. Pulling sequence out by coordinate (`extract`, `isa`)

**What it does.** It reads a stretch of any genome straight from the index, without the FASTA. Ask for a B73
interval and a founder (`--lift --to=FOUNDER`), and it translates both ends into that founder, nudges each end onto
the exact matching base, and returns the founder's sequence. Pieces cut at the same point join exactly. `isa`
builds the small start-point table it needs, once per index. GRITS uses this to rebuild imputed sequence for its
gene-level accuracy check.

## 7. `chain` improvements (RNA-seq)

`--support=best` credits a read to the founders that cover the most of it, instead of always including B73.
Without it, reads from a donor that differs from every assembly were wrongly credited to B73. `--dup-count` lets
you collapse identical reads first and keep their counts.

## The recommended GRITS command

```sh
ropebwt3 lift --ref-prefix=B73 -k 50 -s 200 -t 16 -o idx.k50s200.lift idx.fmd B73.fa  # once per index
ropebwt3 isa idx.fmd                                                                   # once per index (for extract)

ropebwt3 refmap --ref-prefix=B73 --max-occ=-1 -l 19 --lift=idx.k50s200.lift -t 16 \
  --anchor-dist-npy --anchor-dist-thresh=2000 --ps4g=out.ps4g --npy=out.npy \
  --lift-local-win=20000 --lift-chain-nb=3 --lift-chain-interp \
  --placed-min-chr-agree=1 --exact-min-asm-agree=1 --exact-max-ref-hits=1 --placed-max-spread=1000000 --min-core-len=100 \
  --err-noclip --err-max-events=1 \
  idx.fmd reads.fq > out.tsv
```

The four lines of options are: the base settings and outputs; the placement fix (section 1); the preset E read
filters (section 4); the edit filter (section 3). On simulated maize lines outside the reference panel, this command cut
the GRITS model's SNP error at 1x coverage by about 20-30% compared with before these changes (for example
held-out hybrids 11.2% → 7.9%).

## Known limitations and ideas for later

- **Whole reads are dropped, not individual founders.** When a read's founder copies disagree, the filters remove
  the whole read. Keeping the read but removing only the founders whose copy is far away (or placing the read
  separately for each founder) may keep more information. Not yet tested.
- **No read-length-relative match length.** `-l` and `--min-core-len` are fixed numbers of bases. With variable
  read lengths it may be better to say "the whole read", or "the read minus one base". `--err-max-events=0` gives
  "whole read matches one genome exactly", which is close. Not yet tested as an `-l` option.
- **Multi-copy EXACT reads are labelled FILTERED, not MULTI,** and only when the filter is on.
- **One flaky check.** A `--target-hits` integration check depends on thread timing and fails occasionally. This
  was already true before this branch.

## Where the details are

| Feature | Design notes and results |
|---|---|
| Paired-end | `docs/paired-end/PLAN.md` |
| Edit filter (counts) | `docs/err-profile/PLAN.md` |
| Edit filter (no-clip, match fraction, CIGAR) | `docs/edit-match-frac/PLAN.md` |
| Edit-distance statistics | `docs/ed-stats/PLAN.md` |
| Extract and isa | README section "Extracting sequence by coordinate" |
| Everything, as flags | README section "GRITS pipeline" |
