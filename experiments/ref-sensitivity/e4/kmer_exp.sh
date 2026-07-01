#!/bin/sh
# Precision under sequencing error: whole-read refmap --lift vs 75bp k-mer agreement.
set -e
cd /home/esb33/ropebwt3-phg
E4=experiments/ref-sensitivity/e4; DATA=experiments/ref-sensitivity/data
RB=./ropebwt3; LOG=$E4/kmer_exp.log
: > "$LOG"

score_whole() { # $1 = refmap.out
python3 - "$1" "$DATA/truth.tsv" <<'PY'
import sys
truth={}
for ln in open(sys.argv[2]):
    if ln.startswith("qname"): continue
    q,line,c,p,L=ln.rstrip().split("\t"); truth[q]=(line,c,int(p))
n=cor=wr=pl=0
for ln in open(sys.argv[1]):
    F=ln.rstrip("\n").split("\t")
    if len(F)<11 or F[0] not in truth: continue
    q,_,st,_,_,refn,_,cL,cR,_,_=F[:11]; line,sc,sp=truth[q]; n+=1
    rc=None if refn=="." else refn.split("_",1)[1]
    if st in ("UNPLACED","MULTI") or rc is None or cL==".": continue
    pl+=1
    if rc!=sc: wr+=1
    elif abs(int(cL)-sp)<=5e6: cor+=1
print(f"  WHOLE-READ   placed={pl:6d} correct={cor:6d} wrong={wr:5d}  PREC={100*cor/pl:.1f}% RECALL={100*cor/n:.1f}%")
PY
}

for err in 0 0.01 0.02 0.03; do
  echo "===== substitution rate = $err =====" >> "$LOG"
  R=$E4/q.err$err.fa
  if [ "$err" = 0 ]; then cp "$DATA/queries.fa" "$R";
  else python3 "$E4/add_errors.py" --in "$DATA/queries.fa" --out "$R" --rate $err --seed 1 2>>"$LOG"; fi
  # whole read
  "$RB" refmap --ref-prefix=B73 --max-occ=-1 --lift "$DATA/nam4.lift" -t 20 "$DATA/nam4.fmd" "$R" > "$E4/rm.err$err.out" 2>/dev/null
  score_whole "$E4/rm.err$err.out" >> "$LOG"
  # k-mers, step 15 (~6 tiles/read)
  K=$E4/km.err$err.fa
  python3 "$E4/gen_kmers.py" --in "$R" --out "$K" --k 75 --step 15 2>/dev/null
  "$RB" mem -p 8 -t 20 "$DATA/nam4.fmd" "$K" > "$E4/km.err$err.mem" 2>/dev/null
  for ma in 1 2 3; do
    r=$(python3 "$E4/place_kmers.py" --anchors-mem "$E4/anchors.s2000.mem" --kmers-mem "$E4/km.err$err.mem" \
      --truth "$DATA/truth.tsv" --k 75 --max-occ 4 --min-agree $ma --cluster-tol 2000 2>/dev/null | sed -n '2p')
    echo "  KMER(75) agree>=$ma  $r" >> "$LOG"
  done
  rm -f "$K" "$E4/km.err$err.mem" "$R"
done
echo "DONE" >> "$LOG"
