#!/bin/sh
# Module 0 end-to-end (Stage 0 / baseline 5.0): build a synthetic pangenome,
# simulate RNAseq reads with exact truth, index them, run stock ropebwt3 refmap
# (no aligner code) with single-base PS4G, and score against the oracle.
#
# Usage:  sh rnaseq_ps4g/run_stage0.sh [outdir] [nreads] [error]
#         outdir defaults to /workdir/esb33/rnaseq-stage0
# Requires:  ropebwt3 built at the repo root (make omp=0).
set -e

REPO=$(cd "$(dirname "$0")/.." && pwd)
RB="$REPO/ropebwt3"
SIM="$REPO/rnaseq_ps4g/sim"
EVAL="$REPO/rnaseq_ps4g/eval"

OUT=${1:-/workdir/esb33/rnaseq-stage0}
NREADS=${2:-20000}
ERROR=${3:-0.01}
G="$OUT/genomes"
R="$OUT/reads"
mkdir -p "$G" "$R"

[ -x "$RB" ] || { echo "build ropebwt3 first: (cd $REPO && make omp=0)"; exit 1; }

echo "# 1. synthetic pangenome + GFF + truth maps"
python3 "$SIM/sim_genomes.py" --outdir "$G" --seed 7

echo "# 2. RNAseq reads + per-read truth + answer-key BED (Stage-0 defaults)"
python3 "$SIM/sim_rnaseq.py" --simdir "$G" --outdir "$R" \
    --nreads "$NREADS" --seed 11 --error "$ERROR"

echo "# 3. build the RopeBWT3-RefMap index over the pangenome"
IDX="$G/idx.fmd"
"$RB" build -d -o "$IDX" "$G/pangenome.fa"
"$RB" ssa -s16 -o "$IDX.ssa" "$IDX"
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}
     END{if(n)print n"\t"l}' "$G/pangenome.fa" | gzip > "$IDX.len.gz"
"$RB" lift --ref-prefix=B73 -k 61 -s 500 -o "$G/idx.lift" "$IDX" "$G/pangenome.fa"

echo "# 4. stock refmap -> single-base PS4G + per-read PS4G file; table -> reads.refmap"
"$RB" refmap --ref-prefix=B73 --max-occ=-1 --lift "$G/idx.lift" \
    --ps4g "$R/out.ps4g" --ps4g-per-read "$R/reads.ps4g" --bin-size 1 \
    "$IDX" "$R/reads.fq" > "$R/reads.refmap"

echo "# 5. Tier-A scoring vs oracle truth  (also written to $OUT/score.txt)"
python3 "$EVAL/score.py" --truth "$R/truth.tsv" --refmap "$R/reads.refmap" \
    --ps4g "$R/out.ps4g" --ps4g-per-read "$R/reads.ps4g" \
    --gametes "$G/gametes.tsv" | tee "$OUT/score.txt"
