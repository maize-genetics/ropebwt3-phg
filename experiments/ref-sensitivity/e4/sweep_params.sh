#!/bin/sh
# E4 parameter exploration: anchor stride x max_occ. Appends to param_log.txt.
set -e
cd /home/esb33/ropebwt3-phg
DATA=experiments/ref-sensitivity/data; E4=experiments/ref-sensitivity/e4
RB=./ropebwt3
LOG=$E4/param_log.txt
: > "$LOG"
echo "E4 parameter sweep  $(date)" >> "$LOG"

for stride in 10000 5000 2000; do
  AF=$E4/anchors.s$stride.fa; AM=$E4/anchors.s$stride.mem
  if [ "$stride" = 10000 ] && [ -f "$E4/anchors.mem" ]; then
    AM=$E4/anchors.mem   # reuse the stride-10k run already done
  else
    python3 "$E4/gen_anchors.py" --len 100 --stride $stride --out "$AF" \
      B97:"$DATA/B97.chr.fa.gz" Ki3:"$DATA/Ki3.chr.fa.gz" CML247:"$DATA/CML247.chr.fa.gz" 2>>"$LOG"
    "$RB" mem -p 8 -t 20 "$DATA/nam4.fmd" "$AF" > "$AM" 2>/dev/null
  fi
  nanch=$(wc -l < "$AM")
  for occ in 4 6 8; do
    echo "=== stride=$stride anchors=$nanch max_occ=$occ ===" >> "$LOG"
    python3 "$E4/place_reads.py" --anchors-mem "$AM" --reads-mem "$E4/queries.mem" \
      --truth "$DATA/truth.tsv" --max-occ $occ 2>/dev/null >> "$LOG"
    echo "" >> "$LOG"
  done
done
echo "PARAM SWEEP DONE" >> "$LOG"
