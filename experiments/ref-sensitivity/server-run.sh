#!/bin/sh
# Ref-sensitivity, full 4-genome build on a big-RAM Linux server.
#
# ropebwt3 is CPU + RAM bound (no GPU use). The win on a server is enough RAM to
# build the whole index in ONE libsais pass (no slow incremental merge). To fit
# 4 maize genomes (~17 G both-strand symbols) in a single batch you need roughly
# 78 GB free RAM; set MEM_BATCH below to comfortably exceed the total so the
# build never spills into a merge. On Linux, plain `make` enables OpenMP, so
# `-p $NCORES` multithreads the libsais construction.
#
# Usage (from the repo root, after `git clone -b ref-sensitivity ...`):
#     sh experiments/ref-sensitivity/server-run.sh
#
# Tunables (env vars):
#     NCORES     threads               [all cores]
#     MEM_BATCH  build batch size      [40G]  (must exceed total input symbols)
set -e

NCORES="${NCORES:-$(nproc 2>/dev/null || echo 16)}"
MEM_BATCH="${MEM_BATCH:-40G}"
DIR=experiments/ref-sensitivity
DATA=$DIR/data
RB=./ropebwt3
mkdir -p "$DATA"
echo "cores=$NCORES  batch=$MEM_BATCH"

# 0. Build ropebwt3 (Linux gcc => OpenMP on by default).
[ -x "$RB" ] || make -j"$NCORES"

# 1. Download the 4 NAM genomes (public MaizeGDB; ~2.6 GB).
B73=Zm-B73-REFERENCE-NAM-5.0; B97=Zm-B97-REFERENCE-NAM-1.0
KI3=Zm-Ki3-REFERENCE-NAM-1.0; CML=Zm-CML247-REFERENCE-NAM-1.0
BASE=https://download.maizegdb.org
for L in $B73 $B97 $KI3 $CML; do
  [ -f "$DATA/$L.fa.gz" ] || curl -sL --retry 3 -o "$DATA/$L.fa.gz" "$BASE/$L/$L.fa.gz"
done

# 2. Rename to <LINE>_chrN, keep chr1..chr10 only.
rename() { gzip -cd "$1" | awk -v p="$2" '
  /^>/ { n=substr($1,2); if (n ~ /^chr([1-9]|10)$/){k=1; print ">"p"_"n} else k=0; next }
  k { print }' | gzip > "$DATA/$2.chr.fa.gz"; }
[ -f "$DATA/B73.chr.fa.gz" ]    || rename "$DATA/$B73.fa.gz" B73
[ -f "$DATA/B97.chr.fa.gz" ]    || rename "$DATA/$B97.fa.gz" B97
[ -f "$DATA/Ki3.chr.fa.gz" ]    || rename "$DATA/$KI3.fa.gz" Ki3
[ -f "$DATA/CML247.chr.fa.gz" ] || rename "$DATA/$CML.fa.gz" CML247

# 3. SINGLE-BATCH build (no merge) + ssa + len. This is the whole point of the
#    big-RAM box: -m exceeds total input, so libsais runs once on everything.
IDX=$DATA/nam4.fmd
/usr/bin/time -v "$RB" build -d -t "$NCORES" -p "$NCORES" -m "$MEM_BATCH" -o "$IDX" \
  "$DATA/B73.chr.fa.gz" "$DATA/B97.chr.fa.gz" "$DATA/Ki3.chr.fa.gz" "$DATA/CML247.chr.fa.gz" \
  2> "$DATA/build.log" || { echo "build failed; see $DATA/build.log"; exit 1; }
"$RB" ssa -t "$NCORES" -s8 -o "$IDX.ssa" "$IDX"
for L in B73 B97 Ki3 CML247; do gzip -cd "$DATA/$L.chr.fa.gz"; done | \
  awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}
       END{if(n)print n"\t"l}' | gzip > "$IDX.len.gz"
echo "index built:"; ls -la "$IDX"*

# 4. Queries (deterministic: seed 7 reproduces the same 100k on any machine).
python3 "$DIR/gen_queries.py" --n 25000 --len 150 --seed 7 \
  --out "$DATA/queries.fa" --truth "$DATA/truth.tsv" \
  B73:"$DATA/B73.chr.fa.gz" B97:"$DATA/B97.chr.fa.gz" \
  Ki3:"$DATA/Ki3.chr.fa.gz" CML247:"$DATA/CML247.chr.fa.gz"

# 5. Place on B73, then score.
/usr/bin/time -v "$RB" refmap --ref-prefix=B73 -t "$NCORES" "$IDX" "$DATA/queries.fa" \
  > "$DATA/refmap.out" 2> "$DATA/refmap.time"
python3 "$DIR/analyze.py" --refmap "$DATA/refmap.out" --truth "$DATA/truth.tsv" \
  --tol 5000000 | tee "$DATA/report.txt"
