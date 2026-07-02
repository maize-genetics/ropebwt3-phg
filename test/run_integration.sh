#!/bin/bash
# End-to-end integration test for `ropebwt3 refmap`, including the new
# --ps4g/--npy/--label-bed outputs. Builds a tiny synthetic pangenome index
# (docs/examples/pangenome.fa) and checks refmap's output against values
# captured from a known-good run; see docs/examples/refmap.out for the
# original (pre-PS4G) golden output this reuses.
#
# Run via `make test` from the repository root, or directly as
# `cd test && ./run_integration.sh` after `make` has built ../ropebwt3.
set -u
cd "$(dirname "$0")"

RB=../ropebwt3
DIR=../docs/examples
TMP=tmp
IDX="$TMP/pangenome.fmd"

n_pass=0
n_fail=0

fail() { n_fail=$((n_fail+1)); echo "FAIL: $1" >&2; }
pass() { n_pass=$((n_pass+1)); }

check_eq() { # check_eq <desc> <actual> <expected>
	if [ "$2" = "$3" ]; then pass; else fail "$1: got '$2' want '$3'"; fi
}

check_line_in() { # check_line_in <desc> <file> <exact line>
	if grep -qxF -- "$3" "$2" 2>/dev/null; then pass; else fail "$1: no line '$3' in $2"; fi
}

if [ ! -x "$RB" ]; then
	echo "FAIL: $RB not built; run 'make' first" >&2
	exit 1
fi

rm -rf "$TMP"
mkdir -p "$TMP"

# --- build the tiny pangenome index (same steps as docs/examples/run.sh) ---
"$RB" build -d -o "$IDX" "$DIR/pangenome.fa" >/dev/null 2>"$TMP/build.log" || { echo "FAIL: build failed"; cat "$TMP/build.log" >&2; exit 1; }
"$RB" ssa -s8 -o "$IDX.ssa" "$IDX" >/dev/null 2>"$TMP/ssa.log" || { echo "FAIL: ssa failed"; cat "$TMP/ssa.log" >&2; exit 1; }
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}
     END{if(n)print n"\t"l}' "$DIR/pangenome.fa" | gzip > "$IDX.len.gz"

# --- 1. baseline TSV output is unchanged by the new code ---
"$RB" refmap --ref-prefix=B73 -t1 "$IDX" "$DIR/queries.fa" > "$TMP/baseline.tsv" 2>"$TMP/baseline.log"
check_eq "refmap exits 0 (baseline)" "$?" "0"
check_line_in "baseline: ins_query PLACED, 6 carriers" "$TMP/baseline.tsv" \
	"$(printf 'ins_query\t150\tPLACED\t6\tOh43_chr1:+,B97_chr1:+,CML247_chr1:+,Mo17_chr1:+,Ki3_chr1:+,Tx303_chr1:+\tB73_chr1\t+\t300\t300\t0\t500')"
check_line_in "baseline: ref_query EXACT" "$TMP/baseline.tsv" \
	"$(printf 'ref_query\t150\tEXACT\t0\t.\tB73_chr1\t+\t250\t400\t150\t0')"
check_line_in "baseline: unrelated UNPLACED" "$TMP/baseline.tsv" \
	"$(printf 'unrelated\t120\tUNPLACED\t0\t.\t.\t.\t.\t.\t.\t.')"

# --- 2. --ps4g / --npy, no labels ---
"$RB" refmap --ref-prefix=B73 --ps4g "$TMP/out.ps4g" --npy "$TMP/out.npy" -t1 "$IDX" "$DIR/queries.fa" \
	> "$TMP/out.tsv" 2>"$TMP/out.log"
rc=$?
check_eq "refmap exits 0 (--ps4g --npy)" "$rc" "0"

check_line_in "ps4g: header line 1" "$TMP/out.ps4g" "#PS4G"
check_line_in "ps4g: header line 2" "$TMP/out.ps4g" "#version=2.0"
check_line_in "ps4g: gamete table header" "$TMP/out.ps4g" "$(printf '#gamete\tgameteIndex\tcount')"
check_line_in "ps4g: gametes sorted, B73 first (0)" "$TMP/out.ps4g" "$(printf '#B73\t0\t1')"
check_line_in "ps4g: 6 non-reference carriers, alphabetical from index 1" "$TMP/out.ps4g" "$(printf '#B97\t1\t3')"
check_line_in "ps4g: data section header" "$TMP/out.ps4g" "$(printf 'gameteSet\trefContig\trefPosBinned\tcount')"
check_line_in "ps4g: ref_query's EXACT hit contributes gamete B73 only (not equivalent to any carrier)" "$TMP/out.ps4g" \
	"$(printf '0\t_chr1\t0\t1')"
check_line_in "ps4g: the 3 PLACED reads at the insertion share one 6-carrier gameteSet, count aggregated to 3" "$TMP/out.ps4g" \
	"$(printf '1,2,3,4,5,6\t_chr1\t1\t3')"

# invariant: #TotalUniqueCounts must equal the sum of the data rows' count column
declare_total=$(grep '^#TotalUniqueCounts:' "$TMP/out.ps4g" | awk '{print $2}')
row_total=$(awk 'BEGIN{FS="\t"} /^gameteSet/{p=1;next} p{s+=$4} END{print s+0}' "$TMP/out.ps4g")
check_eq "ps4g: #TotalUniqueCounts equals sum of data-row counts" "$declare_total" "$row_total"

# invariant: no PS4G row lists a gamete index that isn't in the gamete table.
# (Deliberately pure awk, not grep -c with a \t pattern: plain POSIX/GNU grep
# does not expand \t to a tab in its default BRE mode, unlike awk's ERE.)
bad_idx=$(awk 'BEGIN{FS="\t"; n=0}
	/^#PS4G$/{next} /^#version=/{next} /^#Command:/{next} /^#TotalUniqueCounts:/{next}
	/^#gamete\tgameteIndex\tcount$/{next}
	/^#/{n++; next}
	/^gameteSet/{p=1; next}
	p{ split($1,a,","); for (i in a) if (a[i]+0 >= n) print "bad:"a[i] }
' "$TMP/out.ps4g")
check_eq "ps4g: every gameteSet index is within the gamete table" "$bad_idx" ""

# --- npy: shape/dtype readable without a numpy dependency (ASCII header is self-describing) ---
npy_header=$(head -c 200 "$TMP/out.npy" | tr -d '\0')
case "$npy_header" in
	*'<i4'*) pass ;;
	*) fail "npy: header does not declare dtype <i4: $npy_header" ;;
esac
shape=$(echo "$npy_header" | grep -oE "'shape': \([0-9]+, [0-9]+\)" | grep -oE '[0-9]+' | tr '\n' ' ')
check_eq "npy: shape is 2 bins x 9 columns (7 gametes + 2 labels)" "$shape" "2 9 "

check_line_in "npy: gametes.tsv lists all 7 samples, B73 first" "$TMP/out.npy.gametes.tsv" "$(printf '0\tB73')"
check_line_in "npy: bins.tsv row for the insertion locus (bin 1, contig stripped)" "$TMP/out.npy.bins.tsv" "$(printf '1\t_chr1\t1')"

# --- 2b. --npy-binary: same locations, but presence (1) instead of the read count ---
"$RB" refmap --ref-prefix=B73 --npy "$TMP/binary.npy" --npy-binary -t1 "$IDX" "$DIR/queries.fa" \
	> /dev/null 2>"$TMP/binary.log"
check_eq "refmap exits 0 (--npy-binary)" "$?" "0"
hdr_len_bin=$(python3 - "$TMP/binary.npy" <<'EOF' 2>/dev/null
import sys
with open(sys.argv[1], 'rb') as f:
    f.read(8)
    hlen = int.from_bytes(f.read(2), 'little')
    print(10 + hlen)
EOF
)
if [ -n "$hdr_len_bin" ]; then
	# row1 = the insertion-locus row (6-carrier gameteSet); B97 is gamete column 1
	cnt_mode=$(od -An -tu4 -j $(( hdr_len_bin + (1*9+1)*4 )) -N4 "$TMP/out.npy" | tr -d ' ')
	bin_mode=$(od -An -tu4 -j $(( hdr_len_bin + (1*9+1)*4 )) -N4 "$TMP/binary.npy" | tr -d ' ')
	check_eq "npy: default mode keeps the read count (3 reads placed there)" "$cnt_mode" "3"
	check_eq "npy: --npy-binary clips the same cell to presence (1)" "$bin_mode" "1"
else
	echo "SKIP: --npy-binary cell check needs python3 to locate the .npy data offset" >&2
fi

# --- 3. --label-bed: diploid training labels ---
cat > "$TMP/labels.bed" <<'EOF'
_chr1	0	100	B73
_chr1	200	400	Ki3	Mo17
EOF
"$RB" refmap --ref-prefix=B73 --npy "$TMP/labeled.npy" --label-bed "$TMP/labels.bed" -t1 "$IDX" "$DIR/queries.fa" \
	> /dev/null 2>"$TMP/labeled.log"
check_eq "refmap exits 0 (--label-bed)" "$?" "0"

read_col() { # read_col <npy> <row> <col> <ncols>  (best-effort raw int32 read via od)
	od -An -tu4 -j $(( 128 + ($2 * $4 + $3) * 4 )) -N4 "$1" | tr -d ' '
}
# both rows have 9 columns (7 gametes + 2 labels); label cols are index 7,8.
# NOTE: npy data starts after a header whose length varies with shape-string
# width, so we locate it by searching for the newline that ends the ASCII
# header instead of assuming a fixed offset.
hdr_len=$(python3 - "$TMP/labeled.npy" <<'EOF' 2>/dev/null
import sys
with open(sys.argv[1], 'rb') as f:
    f.read(8)
    hlen = int.from_bytes(f.read(2), 'little')
    print(10 + hlen)
EOF
)
if [ -n "$hdr_len" ]; then
	g0=$(od -An -tu4 -j $(( hdr_len + (0*9+7)*4 )) -N4 "$TMP/labeled.npy" | tr -d ' ')
	g1=$(od -An -tu4 -j $(( hdr_len + (0*9+8)*4 )) -N4 "$TMP/labeled.npy" | tr -d ' ')
	k0=$(od -An -tu4 -j $(( hdr_len + (1*9+7)*4 )) -N4 "$TMP/labeled.npy" | tr -d ' ')
	k1=$(od -An -tu4 -j $(( hdr_len + (1*9+8)*4 )) -N4 "$TMP/labeled.npy" | tr -d ' ')
	check_eq "npy: bin0 label columns are homozygous B73 (0,0) from the single-sample BED row" "$g0 $g1" "0 0"
	check_eq "npy: bin1 label columns are diploid Ki3/Mo17 (3,4) from the two-sample BED row" "$k0 $k1" "3 4"
else
	echo "SKIP: label column check needs python3 to locate the .npy data offset" >&2
fi

# --- optional deeper check: load with numpy if any interpreter on this box has it ---
NPY_PY=""
for cand in python3 /home/zrm22/mambaforge/envs/ml-impute-env/bin/python /home/zrm22/mambaforge/envs/phg-ml/bin/python; do
	if command -v "$cand" >/dev/null 2>&1 && "$cand" -c "import numpy" >/dev/null 2>&1; then NPY_PY="$cand"; break; fi
done
if [ -n "$NPY_PY" ]; then
	"$NPY_PY" - "$TMP/out.npy" "$TMP/labeled.npy" <<'EOF'
import sys, numpy as np
a = np.load(sys.argv[1])
assert a.shape == (2, 9) and a.dtype == np.int32, f"unexpected out.npy shape/dtype: {a.shape} {a.dtype}"
assert (a[:, -2:] == -1).all(), "out.npy (no --label-bed) should have all-(-1,-1) label columns"
assert a[0].tolist() == [1,0,0,0,0,0,0,-1,-1], f"unexpected out.npy row0: {a[0].tolist()}"
assert a[1].tolist() == [0,3,3,3,3,3,3,-1,-1], f"unexpected out.npy row1: {a[1].tolist()}"
b = np.load(sys.argv[2])
assert b[0,-2:].tolist() == [0,0] and b[1,-2:].tolist() == [3,4], f"unexpected labeled.npy labels: {b[:,-2:].tolist()}"
print("PYNPY_OK")
EOF
	if [ $? -eq 0 ]; then pass; else fail "numpy deep-check ($NPY_PY) raised an assertion"; fi
else
	echo "SKIP: no python3 with numpy found on PATH; skipping the numpy load/shape/value deep-check" >&2
fi

# --- memory safety, if valgrind is available ---
if command -v valgrind >/dev/null 2>&1; then
	valgrind --error-exitcode=99 --leak-check=full -q \
		"$RB" refmap --ref-prefix=B73 --ps4g "$TMP/vg.ps4g" --npy "$TMP/vg.npy" --label-bed "$TMP/labels.bed" -t1 "$IDX" "$DIR/queries.fa" \
		> /dev/null 2>"$TMP/valgrind.log"
	if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks; see $TMP/valgrind.log"; cat "$TMP/valgrind.log" >&2; fi
else
	echo "SKIP: valgrind not installed; skipping the memory-safety pass" >&2
fi

echo "run_integration: $n_pass passed, $n_fail failed"
[ "$n_fail" -eq 0 ]
