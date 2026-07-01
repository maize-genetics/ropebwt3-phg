#!/bin/sh
# Build a small pangenome index and place queries on the B73 reference.
# Run from the repository root after building the tool (`make omp=0` on macOS).
set -e

RB=./ropebwt3
DIR=docs/examples
IDX=$DIR/pangenome.fmd

# 1. Build the BWT (FMD format). `build -d` dumps the static FMD index.
$RB build -d -o "$IDX" "$DIR/pangenome.fa"

# 2. Sampled suffix array (needed to turn matches into coordinates).
$RB ssa -s8 -o "$IDX.ssa" "$IDX"

# 3. Sequence names + lengths sidecar (needed to print names/coordinates).
#    seqtk works too: `seqtk comp pangenome.fa | cut -f1,2 | gzip > pangenome.fmd.len.gz`
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}
     END{if(n)print n"\t"l}' "$DIR/pangenome.fa" | gzip > "$IDX.len.gz"

# 4. RECOMMENDED: build the lift "second SSA" (carrier->reference coordinate map)
#    once, then place queries by projecting carrier hits. Tiny data needs a small
#    k/stride; on real genomes the defaults (-k 100 -s 2000) and the reference
#    FASTA (e.g. B73.fa) are appropriate. See docs/lift-second-ssa.md.
$RB lift --ref-prefix=B73 -k 31 -s 10 -o "$DIR/pangenome.lift" "$IDX" "$DIR/pangenome.fa"
echo "# refmap --lift (project; fast, recommended):"
$RB refmap --ref-prefix=B73 --max-occ=-1 --lift "$DIR/pangenome.lift" "$IDX" "$DIR/queries.fa"

# 5. Fallback (no lift build): the outward walk. It brackets the insertion
#    breakpoint (cL/cR at the A|B junction) and reports the inserted size, where
#    the projection instead returns one approximate coordinate.
echo "# refmap walk (bracket the breakpoint; no setup):"
$RB refmap --ref-prefix=B73 "$IDX" "$DIR/queries.fa"
