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

# 4. Place the queries on the reference (sequences whose name starts with "B73").
$RB refmap --ref-prefix=B73 "$IDX" "$DIR/queries.fa"
