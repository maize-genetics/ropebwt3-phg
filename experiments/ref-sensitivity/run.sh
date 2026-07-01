#!/bin/sh
# Ref-sensitivity experiment: how well does `ropebwt3 refmap` anchor random
# pangenome queries back to the correct B73 coordinate?
#
# Run from the repository root after `make omp=0`. Heavy: downloads ~2.6 GB,
# builds an ~8.8 Gbp index, runs 100k queries. See README.md.
set -e

RB=./ropebwt3
DIR=experiments/ref-sensitivity
DATA=$DIR/data
mkdir -p "$DATA"

# Lines: B73 is the reference; the other three are carriers.
B73=Zm-B73-REFERENCE-NAM-5.0
B97=Zm-B97-REFERENCE-NAM-1.0
KI3=Zm-Ki3-REFERENCE-NAM-1.0
CML=Zm-CML247-REFERENCE-NAM-1.0
BASE=https://download.maizegdb.org

# 1. Download (chr-level genome assemblies, gzipped)
for L in $B73 $B97 $KI3 $CML; do
  [ -f "$DATA/$L.fa.gz" ] || curl -sL --retry 3 -o "$DATA/$L.fa.gz" "$BASE/$L/$L.fa.gz"
done

# 2. Rename to <LINE>_chrN and keep only chr1..chr10 (drop scaffolds), gzipped.
rename() {  # <src.gz> <prefix>
  gzip -cd "$1" | awk -v p="$2" '
    /^>/ { n=substr($1,2); if (n ~ /^chr([1-9]|10)$/){k=1; print ">"p"_"n} else k=0; next }
    k { print }' | gzip > "$DATA/$2.chr.fa.gz"
}
rename "$DATA/$B73.fa.gz" B73
rename "$DATA/$B97.fa.gz" B97
rename "$DATA/$KI3.fa.gz" Ki3
rename "$DATA/$CML.fa.gz" CML247

# 3. Build the index triad (.fmd + .ssa + .len.gz).
IDX=$DATA/nam4.fmd
$RB build -d -o "$IDX" "$DATA/B73.chr.fa.gz" "$DATA/B97.chr.fa.gz" \
                       "$DATA/Ki3.chr.fa.gz" "$DATA/CML247.chr.fa.gz"
$RB ssa -s8 -o "$IDX.ssa" "$IDX"
for L in B73 B97 Ki3 CML247; do gzip -cd "$DATA/$L.chr.fa.gz"; done | \
  awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}
       END{if(n)print n"\t"l}' | gzip > "$IDX.len.gz"

# 4. Draw 100k random 150 bp queries (25k per line; B73 = control).
python3 "$DIR/gen_queries.py" --n 25000 --len 150 --seed 7 \
  --out "$DATA/queries.fa" --truth "$DATA/truth.tsv" \
  B73:"$DATA/B73.chr.fa.gz" B97:"$DATA/B97.chr.fa.gz" \
  Ki3:"$DATA/Ki3.chr.fa.gz" CML247:"$DATA/CML247.chr.fa.gz"

# 5. Place them on B73.
$RB refmap --ref-prefix=B73 -t 12 "$IDX" "$DATA/queries.fa" > "$DATA/refmap.out"

# 6. Score.
python3 "$DIR/analyze.py" --refmap "$DATA/refmap.out" --truth "$DATA/truth.tsv" \
  --tol 5000000 | tee "$DATA/report.txt"
