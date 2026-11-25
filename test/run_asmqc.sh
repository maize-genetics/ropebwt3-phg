#!/bin/bash
# End-to-end test for `ropebwt3 asmqc` (lift mode) and `ropebwt3 asmqc-fa` (FASTA mode) on a synthetic
# reference/assembly pair with planted events at known positions (test/make_asmqc_fixture.py):
# SNPs, small indels, a replacement, a large insertion and deletion, an inversion, a translocation,
# and a solo LTR against a full LTR element. Also smoke-tests scripts/asmqc_summary.py.
#
# Run via `make test` from the repository root, or `cd test && ./run_asmqc.sh` after `make`.
set -u
cd "$(dirname "$0")"
RB=../ropebwt3
TMP=tmp_asmqc
PY=python3
rm -rf "$TMP"; mkdir -p "$TMP"
n_pass=0; n_fail=0
check() { if [ "$1" = 1 ]; then echo "ok   - $2"; n_pass=$((n_pass+1)); else echo "FAIL - $2"; n_fail=$((n_fail+1)); fi; }

$PY make_asmqc_fixture.py "$TMP" || { echo "FAIL: fixture"; exit 1; }

# ---- FASTA mode ----
"$RB" asmqc-fa -t 2 --variants "$TMP/var.tsv" --chains "$TMP/chains.tsv" "$TMP/REF.fa" "$TMP/CAND.fa" > "$TMP/fa.tsv" 2> "$TMP/fa.log" \
	|| { echo "FAIL: asmqc-fa exited non-zero"; cat "$TMP/fa.log" >&2; exit 1; }

eval "$($PY - "$TMP" <<'EOF'
import csv, sys, bisect
d = sys.argv[1]
rows = {r['chrom']: r for r in csv.DictReader(open(f'{d}/fa.tsv'), delimiter='\t')}
truth = [l.rstrip('\n').split('\t') for l in open(f'{d}/truth.tsv')]
var = [l.rstrip('\n').split('\t') for l in open(f'{d}/var.tsv') if not l.startswith('#')]
snp1 = {int(v[2]) for v in var if v[1] == 'REF_chr1' and v[3] == 'SNP'}
t_snp = {int(t[2]) for t in truth if t[0] == 'SNP'}
called = {p for p in snp1 if p <= 400000}
ok = lambda c: 1 if c else 0
def near(typ, pos, length, tol_pos=50, tol_len=100):
    for v in var:
        if v[1] == 'REF_chr1' or v[1] == 'REF_chr3':
            if v[3] == typ and abs(int(v[2]) - pos) <= tol_pos:
                L = int(v[4]) if typ in ('DEL', 'REP') else int(v[5])
                if abs(L - length) <= tol_len: return 1
    return 0
small = [t for t in truth if t[0] in ('INS3', 'DEL3')]
found = sum(1 for t in small if any(v[3] == ('INS' if t[0] == 'INS3' else 'DEL') and v[1] == 'REF_chr1'
                                    and abs(int(v[2]) - int(t[2])) <= 5 for v in var))
c1, c2, c3 = rows['REF_chr1'], rows['REF_chr2'], rows['REF_chr3']
snp3 = sum(1 for v in var if v[1] == 'REF_chr3' and v[3] == 'SNP' and 60000 < int(v[2]) <= 61600)
print(f"SNP_RECALL={ok(len(called & t_snp) >= 0.99 * len(t_snp))}")
print(f"SNP_FP={ok(len(called - t_snp) <= 2)}")
print(f"SMALL_INDELS={ok(found >= 17)}")
print(f"REP={ok(near('REP', 700001, 3000))}")
print(f"INS5K={ok(near('INS', 800001, 5000))}")
print(f"DEL20K={ok(near('DEL', 900001, 20000))}")
print(f"INV={ok(abs(int(c1['inv_bp']) - 100000) <= 5000)}")
print(f"TRANSLOC={ok(abs(int(c2['other_chr_bp']) - 100000) <= 5000)}")
print(f"SOLO_LTR={ok((near('DEL', 60001 + 1500, 28500, tol_pos=200, tol_len=300) or near('REP', 60001 + 1500, 28500, tol_pos=200, tol_len=300)) and 100 <= snp3 <= 200)}")
print(f"N_SNP_COL={ok(int(c1['numSNPs']) == len(snp1))}")
EOF
)"
check "$SNP_RECALL" "asmqc-fa recovers >= 99% of planted SNP positions"
check "$SNP_FP" "asmqc-fa: at most 2 false SNPs in the SNP block"
check "$SMALL_INDELS" "asmqc-fa: >= 17 of 19 planted 3 bp indels at the right place"
check "$REP" "asmqc-fa: 3 kb replacement reported as REP 3000/3000 at 700,001"
check "$INS5K" "asmqc-fa: 5 kb insertion at 800,001"
check "$DEL20K" "asmqc-fa: 20 kb deletion at 900,001"
check "$INV" "asmqc-fa: ~100 kb inverted on chr1"
check "$TRANSLOC" "asmqc-fa: ~100 kb of chr2 on another chromosome"
check "$SOLO_LTR" "asmqc-fa: solo LTR aligned (100-200 SNPs), the other ~28.5 kb of the element deleted/replaced"
check "$N_SNP_COL" "asmqc-fa: numSNPs column equals the --variants SNP count"

# ---- lift mode: build a two-genome index, lift, asmqc ----
sed 's/^>/>REF_/' "$TMP/REF.fa" > "$TMP/REFp.fa"
cat "$TMP/REFp.fa" "$TMP/CAND.fa" > "$TMP/pair.fa"
"$RB" build -d -o "$TMP/pair.fmd" "$TMP/pair.fa" >/dev/null 2>"$TMP/build.log" || { echo "FAIL: build"; cat "$TMP/build.log" >&2; exit 1; }
"$RB" ssa -s8 -o "$TMP/pair.fmd.ssa" "$TMP/pair.fmd" >/dev/null 2>"$TMP/ssa.log" || { echo "FAIL: ssa"; exit 1; }
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)} END{if(n)print n"\t"l}' "$TMP/pair.fa" | gzip > "$TMP/pair.fmd.len.gz"
"$RB" lift --ref-prefix=REF -k 100 -s 200 -t 2 -o "$TMP/pair.lift" "$TMP/pair.fmd" "$TMP/REFp.fa" >/dev/null 2>"$TMP/lift.log" || { echo "FAIL: lift"; cat "$TMP/lift.log" >&2; exit 1; }
"$RB" asmqc --ref-prefix=REF -k 100 -s 200 -t 2 "$TMP/pair.lift" "$TMP/pair.fmd.len.gz" > "$TMP/lift.tsv" 2>"$TMP/asmqc.log" || { echo "FAIL: asmqc"; cat "$TMP/asmqc.log" >&2; exit 1; }
eval "$($PY - "$TMP" <<'EOF'
import csv, sys
r = {x['chrom']: x for x in csv.DictReader(open(f'{sys.argv[1]}/lift.tsv'), delimiter='\t')}
f = {x['chrom']: x for x in csv.DictReader(open(f'{sys.argv[1]}/fa.tsv'), delimiter='\t')}
print(f"LIFT_MAPPED={1 if abs(float(r['REF_chr2']['frac_mapped']) - 500000/600000) < 0.02 else 0}")
print(f"LIFT_AGREES={1 if all(abs(float(r[c]['frac_mapped']) - float(f[c]['frac_mapped'])) < 0.01 for c in ('REF_chr1','REF_chr2','REF_chr3')) else 0}")
EOF
)"
check "$LIFT_MAPPED" "asmqc (lift mode): chr2 mapped fraction ~5/6 after the 100 kb translocation"
check "$LIFT_AGREES" "asmqc (lift mode) and asmqc-fa agree on mapped fractions (within 1%)"

# ---- summary script ----
$PY ../scripts/asmqc_summary.py --qc "$TMP/fa.tsv" --chains "$TMP/chains.tsv" --out-tsv "$TMP/summary.tsv" --out-json "$TMP/summary.json" 2>"$TMP/summary.log"
check "$( [ $? = 0 ] && grep -q 'chr2: only 83.3% of B73 in syntenic chains' "$TMP/summary.tsv" && echo 1 || echo 0 )" "asmqc_summary.py runs and flags chr2 (83.3% mapped after the translocation)"

echo "asmqc: $n_pass passed, $n_fail failed"
[ "$n_fail" = 0 ]
