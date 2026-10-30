# `extract` and `isa`: technical details

The README section [Extracting sequence by coordinate](../../README.md#extract) explains how to use
`extract`. This page keeps the full technical description: exact coordinate rules, how entry points
are chosen, how the B73-to-founder projection and end snapping work, the extra `--project-only`
columns, and the measurements behind the defaults.

**Coordinates.** 0-based, half-open `[start,end)` on the forward strand of the named
sequence, i.e. the sequence as it appeared in the FASTA given to `ropebwt3 build`
(names as in `.len.gz`, e.g. `Oh43_chr1`). This is NOT samtools' convention:
`CTG:S-E` here is `samtools faidx CTG:S+1-E`. Strand `-` returns the reverse
complement of the same interval (`samtools faidx -i`). Bases come back upper
case `ACGTN`: the BWT stores no case, and any non-ACGT base was stored as `N` at
build time. Region strings are `CONTIG:START-END[:+|-]` (the contig may itself
contain `:`).

**Errors and clamping.** Unknown contig, unparsable region, `start >= end`, or `end`
past the contig end are errors: the query is reported on stderr with its reason,
no record is written for it, the other queries still are, and the exit status is
1. Nothing is clamped silently: `--clamp` clips to the contig and marks the FASTA
header with ` clamped` (an interval that is empty after clipping is still an error).

**Output.** FASTA, header `>qid` (TSV) or the region string as given, one line per
sequence unless `-w INT`. With `--lift` the header is `>qid CONTIG:START-END:STRAND
ref=REF:START-END`.

**Entry points (`ropebwt3 isa`).** `ropebwt3 ssa -s S` keeps the suffix-array value
of every BWT row with `(row - #sequences) % 2^S == 0`, i.e. it samples by *row*,
which lands on effectively random text positions (mean spacing 2^S, geometric
gaps). Inverting these (row -> sequence, offset) pairs gives exactly the entry
points needed, so no new pass over the BWT is made: `isa` streams the `.ssa`
array, keeps per sequence the sample with the smallest offset in each 2^b-base
bucket (`-b`, default 12 = 4 kb), and adds one exact entry per sequence for free,
its end: row `sid` is the sentinel row that ends BWT string `sid` (what `get sid`
walks from). Extraction binary-searches the first entry at offset >= end and
LF-walks `entry - start` steps. On the 25-founder maize index (55.2 Gbp forward,
3.4 GB `.ssa`, S=8), `isa -b 12` takes 6 s and 0.32 GB RAM and writes 13.5 M
entries (162 MB); `-b 10` 7 s, 1.25 GB, 52.9 M entries (636 MB). `extract` loads
the `.isa` in 0.05 s; without one it derives the same entries from the `.ssa` in
memory (6 s per run), so `isa` is worth running once per index.

**Strands.** `build` stores sequence j of `.len.gz` (FASTA order) as BWT string 2j
(forward) and 2j+1 (reverse complement); `ropebwt3 get` numbers by BWT string.
`isa` samples forward strings only (`-2` adds the reverse ones, which `extract`
never needs): a `-` query is the forward interval, reverse-complemented.

**Projection (`--lift`).** Queries are then *reference* intervals (`--ref-prefix`
contigs), TSV `qid ref_contig ref_start ref_end [assembly]` (extra columns ignored)
or region strings plus `--to=ASM` (the name prefix before `_`). The first and last
reference base are each projected into the assembly with the same projector refmap
uses (`lift.c`, run with the reference as the key axis over the founder's anchors
on that reference contig: the nearest collinear flanking anchors within
`--lift-local-win`, else the `--lift-win` median/MAD fit). The interval is
`[proj(first), proj(last)+1)` on `+`, or `[proj(last), proj(first)+1)` on `-` for an
inverted locus (the sequence is then reverse-complemented, so it reads in
reference orientation). If either end does not project, or the ends land on
different contigs or orientations, or the end precedes the start, the query fails
with that reason (stderr, and the `status` column of `--project-only`), never a
guess; projection failures do not set the exit status. `--project-only` prints
`qid ref_contig ref_start ref_end assembly contig start end strand status method snap snap_shift`,
where `method` is the start/end projector (`local` bracket, `local1` one-sided,
`median`, `identity` for the reference itself). Differences from refmap: the local
bracket is ON by default (`--lift-local-win=20000`; `0` gives refmap's built-in
default, the wide median only, which for an interval is typically kb off); the
bracket picks the better-fitting slope rather than +1 first (refmap's forward
projector is unchanged); and on `-` loci the result is shifted by k-1 so a base
lands on its homologous base (`--lift-k`, default 100 = `lift`'s default; set it
if the `.lift` was built with another `-k`, it is not stored in the file).
`--lift-win=0` disables the median fallback entirely (only bracketed projections).

**Span sanity check.** After projection and snapping, a query FAILS when its projected
span is both more than `--max-span-excess` (50000) bp longer than the reference
interval AND more than `--max-span-ratio` (3.0) times as long; 0 for either option
turns the check off. It applies to clamped (`ok_clamped`) results too, and is reported
like any projection failure (`projection_failed: projected span implausibly long (PRJ bp
for REF bp)` in the status column and on stderr, no FASTA record, exit status unchanged).
The rule is shared with the GRITS Python builder. Such spans come from the coarse
wide-median fallback: on 1,000 random B73 8 kb windows into the maize founders it
removes 21 of 890 OK projections (19 involve a `median` end), and no OK projection
more than 50 kb longer than its reference interval remains.

**Snapping projected ends to the homologous base (default; `--no-snap` turns it off).**
The lift only knows anchors every few hundred bp, so a projected end is interpolated
and is off by any indel between it and its anchor. Each interval end is therefore a
*boundary* x of the reference (the fencepost between bases x-1 and x; the start s or
the end e), projected through base x-1 and refined by aligning the reference flanks of
x (left `[x-F, x)`, right `[x, x+F)`, read from the index itself; `--snap-flank`,
F = 64) to the founder within +-`--snap-radius` (200) of the projection (in reference
orientation; on `-` loci the window is reverse-complemented). The aligner (`snap.c`;
ropebwt3 had no pairwise aligner) is an affine-gap DP, global on the flank and free at
both ends of the window (match 1, mismatch 4, gap 6+1/bp). A flank placement is
accepted if >= half of it aligns at gap-compressed identity >= `--snap-min-id` (0.85)
and no placement >= 8 bp away scores within `--snap-margin` (5, about one mismatch) of
the best; at low identity it is retried within +-`--snap-radius2` (6000; 0 = off),
which catches large indels between x and its anchors.

As `simulate_pangenome.py` defines the truth, an END at x is one past the founder base
aligned to the last reference base before x present in the founder (left flank), a
START at x is the founder base aligned to the first present base from x on (right
flank). Both come from the same placements around x, so two pieces `[a,x)` and `[x,b)`
projected into the same founder always meet at x:
0. first, the two flanks are aligned as ONE query `[x-F, x+F)` (+-R1, then +-R2), so a
   deletion at x is an internal gap of the affine alignment, placed by the whole
   context. (A one-sided flank must align its deleted tail as a terminal gap, 6+1/bp,
   which loses to a mismatch plus chance matches: boundaries inside small deletions
   came out 1-6 bp past the junction.) If x-1 or x is not present there (a deletion
   at x), end = start = one past the last present base before x (or, with none in the
   query, the first present base from x on). If both are present, steps 1-4 decide
   (one-sided flanks leave an insertion at x out exactly; in a joint alignment an
   insertion's gap can slide over equal bases);
1. both flanks placed, bases x-1 and x both present (no deletion at x), right not
   overlapping left: end = after the left flank, start = the right flank; anything in
   between is founder sequence inserted exactly at x, which belongs to neither piece;
2. otherwise the flank that was placed (left preferred) defines BOTH: a deletion
   containing x maps both pieces to the junction, and founder-only bases replacing
   the reference there go to the right piece;
3. neither placed: the flank itself may be deleted (an indel longer than F): step the
   left flank, then the right flank, away from x (F, 2F, 4F, ... <= radius2, anywhere
   on the reference contig) until it aligns, then align the whole stretch from x to it
   with the deleted run next to x clipped; again one junction for both pieces;
4. otherwise the lift projection is kept (the same for both pieces, since both
   project base x-1; also with `--no-snap`).
Nothing depends on the interval itself, only on x and the founder. An interval whose
two boundaries meet (it lies entirely in a deletion) fails as `interval deleted in the
founder (empty projection)`.

`--project-only` adds four columns AFTER the existing ones, `snap` (start,end labels:
`snap`, `snap_wide`, `snap_walk` = placed by the end's own flank (start: right, end:
left); `snap_x`, `snap_x_wide`, `snap_x_walk` = placed by the flank across the boundary
(a junction, or a fallback); `proj_ambiguous`, `proj_lowid`, `proj_short` (flank < 16
bp) = projection kept; `identity`; `proj` with `--no-snap`), `snap_shift` (bp each
end's fencepost moved, forward coordinates), and `raw_start`, `raw_end`; the FASTA
header gets ` snap=START,END`.

`raw_start` / `raw_end` are the unsnapped lift projections of the reference boundaries
ref_start / ref_end, on EVERY row (also failed ones, so another tool can snap or rescue
them): `CONTIG:POS:STRAND`, or `.` if that end did not project (no confident lift
projection; also `.` on input errors). Each end carries its own contig and strand,
since a row can fail exactly because they differ. POS is a fencepost in FORWARD
coordinates of CONTIG, from projecting reference base x-1 for boundary x (base 0 for
x = 0) to its homolog h: on `+`, POS = h+1 (h for x = 0); on `-`, POS = h (h+1 for
x = 0). So when both ends share contig and strand, the raw interval is `[raw_start,
raw_end)` on `+` and `[raw_end, raw_start)` on `-` -- exactly what `--no-snap` reports
when the row succeeds (POS is not clamped; it can lie outside the contig, which
`--no-snap` would report as clamped or past the end). The reference itself gives
`REF_CONTIG:ref_start:+` / `REF_CONTIG:ref_end:+`. The first 10 columns are unchanged
(the GRITS TableProvider reads them).
F = 64 makes a chance hit in a 12 kb window negligible while staying short enough that
most flanks carry no indel; windows wider than 8F are first cut to the regions sharing
an exact 10-mer with the flank (the same answer as the full DP; tested), which keeps the
wide retry cheap. Snapping needs the FMD and `.isa` even with `--project-only`.

Measured on the simulated pangenome (`test/simulate_pangenome.py`: 3 x 150-300 kb
chromosomes, 5 founders with 1% SNPs, small and 0.5-5 kb indels, a 6-copy 2 kb
element; 500 random 8 kb B73 windows) with an `s200` lift:

| | ends exact | intervals exact |
|---|---|---|
| `--lift-local-win=0` (refmap's wide median), no snap | 2.6% | 0.2% |
| local bracket, `--no-snap` | 65.6% | 44.6% |
| local bracket + snap (default) | 99.4% | 98.8% |

The 6 remaining ends are all deletion-boundary placements: 6 are equivalent alignments
(the bases separating the two placements are identical, so the sequence cannot tell
them apart; 1-3 bp). Splitting each window at its midpoint and projecting both
halves into the same founder (`test/check_extract_junctions.py`), 500/500 junctions
equal the truth gap (499 exactly 0, 1 an insertion at the midpoint). Boundaries placed
strictly inside every >= 2 bp deletion of the fixture and exactly at every insertion
(`test/check_extract_boundaries.py`, pieces of 3 kb either side): deletions 789/1,060 at
the true junction plus 235 at an equivalent placement of the same deletion (96.6%; the
one-sided flanks before the joint step: 142 + 178), insertions 1,549/1,554 exact. On the 25-founder maize index (200 random B73 8 kb windows,
`_v2_s200.lift`) 66% of ends snap (`snap` 198, `snap_wide` 16, `snap_walk` 26 of 364
projected ends); of the 145 `proj_lowid`, 117 had only the coarse wide-median
projection (no collinear anchors within 20 kb, i.e. non-syntenic sequence), where the
flank is genuinely not near the projection.

**Speed.** 1,000 random 8 kb windows across all 25 maize assemblies: 6.5 s on one
thread (6.5 ms per window, ~10.0 M LF steps), 1.5 s with `-t4`; plus 1.6 s to
load the FMD (5.0 GB RSS; `-M` mmaps it instead) and 0.05 s for the `.isa`. With
`--lift` (1,000 random B73 windows): 6.4 s without snapping, 27.6 s with (`-t4`: 1.7 s
vs 7.3 s); the extra time is mostly the LF walks for the flanks and windows (39.6 M vs
12.2 M steps), concentrated in the ends that fall back to the wide/walk search. With
the shared boundary computation and the joint two-flank step the same run takes 47.0 s
(`-t4`: 12.5 s; 71.6 M LF steps).
