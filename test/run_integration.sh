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
# These sections exercise the DEPRECATED flank-walking path and are pinned to golden output
# captured from it, so they pass --walk explicitly (refmap now requires --lift or --walk rather
# than silently defaulting to walking). Section 6 covers the supported --lift path; the two do
# not agree on insertion placement, which is documented there.
"$RB" refmap --ref-prefix=B73 --walk -t1 "$IDX" "$DIR/queries.fa" > "$TMP/baseline.tsv" 2>"$TMP/baseline.log"
check_eq "refmap exits 0 (baseline)" "$?" "0"
check_line_in "baseline: ins_query PLACED, 6 carriers" "$TMP/baseline.tsv" \
	"$(printf 'ins_query\t150\tPLACED\t6\tOh43_chr1:+,B97_chr1:+,CML247_chr1:+,Mo17_chr1:+,Ki3_chr1:+,Tx303_chr1:+\tB73_chr1\t+\t300\t300\t0\t500')"
check_line_in "baseline: ref_query EXACT" "$TMP/baseline.tsv" \
	"$(printf 'ref_query\t150\tEXACT\t0\t.\tB73_chr1\t+\t250\t400\t150\t0')"
check_line_in "baseline: unrelated UNPLACED" "$TMP/baseline.tsv" \
	"$(printf 'unrelated\t120\tUNPLACED\t0\t.\t.\t.\t.\t.\t.\t.')"

# --- 2. --ps4g / --npy, no labels ---
"$RB" refmap --ref-prefix=B73 --walk --ps4g "$TMP/out.ps4g" --npy "$TMP/out.npy" -t1 "$IDX" "$DIR/queries.fa" \
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
	"$(printf '0\tchr1\t0\t1')"
check_line_in "ps4g: the 3 PLACED reads at the insertion share one 6-carrier gameteSet, count aggregated to 3" "$TMP/out.ps4g" \
	"$(printf '1,2,3,4,5,6\tchr1\t1\t3')"
# refContig must be the bare contig name ("chr1"), never a stray leading '_' left
# over from a --ref-prefix that (like the real-world --ref-prefix=B73 used here
# and in the project's own scripts) doesn't itself include the separator.
check_eq "ps4g: no data row's refContig has a stray leading underscore" \
	"$(awk 'BEGIN{FS="\t"} /^gameteSet/{p=1;next} p && $2 ~ /^_/{print $2}' "$TMP/out.ps4g")" ""

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
check_line_in "npy: bins.tsv row for the insertion locus (bin 1, contig stripped)" "$TMP/out.npy.bins.tsv" "$(printf '1\tchr1\t1')"

# --- 2b. --npy-binary: same locations, but presence (1) instead of the read count ---
"$RB" refmap --ref-prefix=B73 --walk --npy "$TMP/binary.npy" --npy-binary -t1 "$IDX" "$DIR/queries.fa" \
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
# chrom is the bare contig name ("chr1", not "_chr1"): this exercises
# bed_resolve_contig's contig_suffix-based match against --ref-prefix=B73
# (no trailing underscore), the same convention the project's own scripts use.
cat > "$TMP/labels.bed" <<'EOF'
chr1	0	100	B73
chr1	200	400	Ki3	Mo17
EOF
"$RB" refmap --ref-prefix=B73 --walk --npy "$TMP/labeled.npy" --label-bed "$TMP/labels.bed" -t1 "$IDX" "$DIR/queries.fa" \
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

# --- generic PS4G<->npy correspondence check (verify_ps4g_npy.py), same numpy interpreter ---
if [ -n "$NPY_PY" ]; then
	if "$NPY_PY" verify_ps4g_npy.py "$TMP/out.ps4g" "$TMP/out.npy" > "$TMP/verify.log" 2>&1; then
		pass
	else
		fail "verify_ps4g_npy.py found a mismatch between out.ps4g and out.npy"
		cat "$TMP/verify.log" >&2
	fi
else
	echo "SKIP: verify_ps4g_npy.py needs a python3 with numpy" >&2
fi

# --- 4. --target-hits: stop once N PLACED/EXACT records are written ---
# NOTE on exact counts below: the outer read/process/write pipeline
# (kt_pipeline) always runs with 2 scheduling workers regardless of -t (-t only
# controls the per-batch compute parallelism inside each step), and its
# worker-ordering rule guarantees round N+1's read can start as soon as round
# N leaves its read step -- i.e. before round N's write has updated the hit
# counter. With -K set small enough to force ~1 record per round, this makes
# the overshoot an exact, deterministic "+1 record" every time (verified by
# running the same case 5x and getting an identical count each time), not a
# vague "a few more". `-t1` doesn't remove this overshoot (it's not about
# compute threading) but is used below anyway for a clean, singlethreaded run.

# 4a. PLACED counts toward the target too, not just EXACT (existing fixture:
# ins_query=PLACED, ins_query_rc=PLACED, ref_query=EXACT, in that file order).
"$RB" refmap --ref-prefix=B73 --walk --target-hits=2 -K 1 -t1 "$IDX" "$DIR/queries.fa" \
	> "$TMP/target2.tsv" 2>"$TMP/target2.log"
check_eq "refmap exits 0 (--target-hits=2)" "$?" "0"
check_eq "target-hits=2: exactly 3 rows (2 PLACED to hit the target + 1 deterministic overshoot record)" \
	"$(wc -l < "$TMP/target2.tsv" | tr -d ' ')" "3"
check_eq "target-hits=2: stops right after ref_query, never reaches ins_query_1snp/unrelated" \
	"$(cut -f1 "$TMP/target2.tsv" | tr '\n' ',')" "ins_query,ins_query_rc,ref_query,"

# 4b. Tight, hand-computed fixture: 30 records alternating a known-EXACT 60bp
# chunk (e1..e15, lifted from pangenome.fa's B73_chr1) and the known-UNPLACED
# junk sequence already used as queries.fa's "unrelated" read (u1..u15).
EXACT_SEQ=$(awk 'NR==2' "$DIR/pangenome.fa")
JUNK_SEQ=$(awk '/^>unrelated/{getline; print}' "$DIR/queries.fa")
: > "$TMP/target_hits.fa"
for i in $(seq 1 15); do
	printf '>e%d\n%s\n' "$i" "$EXACT_SEQ" >> "$TMP/target_hits.fa"
	printf '>u%d\n%s\n' "$i" "$JUNK_SEQ" >> "$TMP/target_hits.fa"
done

"$RB" refmap --ref-prefix=B73 --walk --target-hits=5 -K 10 -t1 "$IDX" "$TMP/target_hits.fa" \
	> "$TMP/target5.tsv" 2>"$TMP/target5.log"
check_eq "refmap exits 0 (--target-hits=5, new fixture)" "$?" "0"
check_eq "target-hits=5: exactly 10 rows (5th EXACT is record 9; +1 deterministic overshoot = record 10)" \
	"$(wc -l < "$TMP/target5.tsv" | tr -d ' ')" "10"
check_eq "target-hits=5: rows are e1,u1,...,e5,u5 in order, nothing from e6/u6 onward" \
	"$(cut -f1 "$TMP/target5.tsv" | tr '\n' ',')" "e1,u1,e2,u2,e3,u3,e4,u4,e5,u5,"
check_eq "target-hits=5: all 5 EXACT rows are actually status EXACT" \
	"$(awk -F'\t' '$1 ~ /^e/{print $3}' "$TMP/target5.tsv" | sort -u | tr '\n' ',')" "EXACT,"
check_eq "target-hits=5: all 5 junk rows are actually status UNPLACED (not filtered out)" \
	"$(awk -F'\t' '$1 ~ /^u/{print $3}' "$TMP/target5.tsv" | sort -u | tr '\n' ',')" "UNPLACED,"

# 4c. Shortfall: requesting more PLACED/EXACT than exist must not error, must
# still process the whole input, and must warn with the actual count found.
"$RB" refmap --ref-prefix=B73 --walk --target-hits=100 -K 10 -t1 "$IDX" "$TMP/target_hits.fa" \
	> "$TMP/target100.tsv" 2>"$TMP/target100.log"
check_eq "refmap exits 0 even when the target is never reached (not an error)" "$?" "0"
check_eq "target-hits=100 (unreachable, only 15 EXACT exist): all 30 input records still processed" \
	"$(wc -l < "$TMP/target100.tsv" | tr -d ' ')" "30"
check_line_in "target-hits=100: WARNING reports the requested target and the actual count found" "$TMP/target100.log" \
	"WARNING: --target-hits=100 requested but the input was exhausted after only 15 PLACED/EXACT records"

# --- 5. the lift path, chaining, and PAV breakpoint anchoring ---
# --lift is the supported resolution path (walking is deprecated); `chain` unites a read's SMEMs
# for RNAseq; and `chain --lift` adds the PAV fallback for reads with no reference-hitting SMEM.
# The fixture is ideal for this: ins_query is 150 bp from the middle of a 500 bp insertion that
# every non-B73 assembly carries and the reference lacks.
LIFT="$TMP/pangenome.lift"
# -k/-s must be far below the defaults (100/2000): the whole fixture reference is 600 bp.
"$RB" lift --ref-prefix=B73 -k 31 -s 50 -o "$LIFT" "$IDX" "$DIR/pangenome.fa" \
	>/dev/null 2>"$TMP/lift.log"
check_eq "lift builds on the fixture" "$?" "0"
check_eq "lift file is non-empty" "$([ -s "$LIFT" ] && echo yes || echo no)" "yes"

# refmap must now be told which resolution path to use rather than silently walking.
"$RB" refmap --ref-prefix=B73 -t1 "$IDX" "$DIR/queries.fa" >/dev/null 2>"$TMP/nopath.log"
check_eq "refmap without --lift/--walk exits non-zero" "$([ $? -ne 0 ] && echo yes || echo no)" "yes"
if grep -q -- "--lift=FILE" "$TMP/nopath.log"; then pass; else fail "refmap error names --lift as the standard path"; fi

# refmap --lift agrees with the walk golden on STATUS and CARRIERS.
"$RB" refmap --ref-prefix=B73 --lift="$LIFT" -t1 "$IDX" "$DIR/queries.fa" \
	> "$TMP/lift.tsv" 2>"$TMP/liftrun.log"
check_eq "refmap --lift exits 0" "$?" "0"
check_eq "refmap --lift: ins_query PLACED with 6 carriers" \
	"$(awk -F'\t' '$1=="ins_query"{print $3"/"$4}' "$TMP/lift.tsv")" "PLACED/6"
check_line_in "refmap --lift: ref_query EXACT, unchanged from walk" "$TMP/lift.tsv" \
	"$(printf 'ref_query\t150\tEXACT\t0\t.\tB73_chr1\t+\t250\t400\t150\t0')"
check_eq "refmap --lift: unrelated still UNPLACED" \
	"$(awk -F'\t' '$1=="unrelated"{print $3}' "$TMP/lift.tsv")" "UNPLACED"
# NOTE: walk and lift do NOT agree on where insertion content lands. Walk reports the insertion
# point and size (300 / ins_size 500); lift reports a projected colinear coordinate with no size,
# because --lift-mad defaults to 200 kb and nothing in a 1.1 kb fixture exceeds it. Recorded, not
# asserted as equal -- `chain --lift` below is the path that resolves insertions properly.

# `chain` alone is reference-anchored: a read living entirely in the insertion has no reference
# SMEM and is dropped. Only the read that spans reference sequence survives.
"$RB" chain --ref-prefix=B73 -t1 "$IDX" "$DIR/queries.fa" > "$TMP/chain.tsv" 2>"$TMP/chain.log"
check_eq "chain exits 0" "$?" "0"
check_line_in "chain: ref_query placed at the exact reference base, B73 only" "$TMP/chain.tsv" \
	"$(printf 'ref_query\tchr1\t250\t0')"
check_eq "chain without --lift emits no pav: rows" \
	"$(grep -c 'pav:' "$TMP/chain.tsv")" "0"
check_eq "chain without --lift drops the insertion reads" \
	"$(grep -cE '^ins_query' "$TMP/chain.tsv")" "0"

# `chain --lift` rescues them via PAV breakpoint anchoring. --pav-grid=0 disables the 5 kb
# presence/absence grid, which is wider than this entire 600 bp fixture contig.
"$RB" chain --ref-prefix=B73 --lift="$LIFT" --pav-grid=0 -t1 "$IDX" "$DIR/queries.fa" \
	> "$TMP/pav.tsv" 2>"$TMP/pav.log"
check_eq "chain --lift exits 0" "$?" "0"
# class column: 1 = INSERTION (absent from the reference, placed at the nearest breakpoint),
# 0 = DIVERGED (present in the reference but too divergent to share a SMEM).
check_line_in "chain --lift: ins_query -> pav row, all 6 non-reference assemblies, class insertion" \
	"$TMP/pav.tsv" "$(printf 'ins_query\tpav:chr1\t250\t1,2,3,4,5,6\t1')"
check_line_in "chain --lift: reverse complement placed identically" \
	"$TMP/pav.tsv" "$(printf 'ins_query_rc\tpav:chr1\t250\t1,2,3,4,5,6\t1')"
check_line_in "chain --lift: one SNP inside the insertion still places" \
	"$TMP/pav.tsv" "$(printf 'ins_query_1snp\tpav:chr1\t250\t1,2,3,4,5,6\t1')"
check_line_in "chain --lift: ref_query is an ordinary row, not pav:" "$TMP/pav.tsv" \
	"$(printf 'ref_query\tchr1\t250\t0')"
check_eq "chain --lift: pav rows never contain the reference assembly (index 0)" \
	"$(awk -F'\t' '$2 ~ /^pav:/ {n=split($4,g,","); for(i=1;i<=n;i++) if(g[i]=="0") bad++} END{print bad+0}' "$TMP/pav.tsv")" "0"
check_eq "chain --lift: unrelated read produces no row" \
	"$(grep -c '^unrelated' "$TMP/pav.tsv")" "0"

# the specificity floor: no assembly-only match reaches 200 bp in 150 bp reads, so all pav rows go
"$RB" chain --ref-prefix=B73 --lift="$LIFT" --pav-grid=0 --pav-min-len=200 -t1 "$IDX" "$DIR/queries.fa" \
	> "$TMP/pavmin.tsv" 2>/dev/null
check_eq "--pav-min-len above the read length suppresses every pav row" \
	"$(grep -c 'pav:' "$TMP/pavmin.tsv")" "0"
check_line_in "--pav-min-len does not affect reference-anchored rows" "$TMP/pavmin.tsv" \
	"$(printf 'ref_query\tchr1\t250\t0')"

# the 5 kb default grid snaps the breakpoint; on this fixture that collapses it to 0
"$RB" chain --ref-prefix=B73 --lift="$LIFT" -t1 "$IDX" "$DIR/queries.fa" > "$TMP/pavgrid.tsv" 2>/dev/null
check_eq "default --pav-grid snaps the breakpoint (250 -> 0 at 5 kb)" \
	"$(awk -F'\t' '$1=="ins_query"{print $3}' "$TMP/pavgrid.tsv")" "0"

# --diverged-rows policy switch is accepted; this fixture yields only insertions, so dropping
# diverged rows must not change the output
"$RB" chain --ref-prefix=B73 --lift="$LIFT" --pav-grid=0 --diverged-rows=drop -t1 "$IDX" \
	"$DIR/queries.fa" > "$TMP/pavdrop.tsv" 2>/dev/null
check_eq "--diverged-rows=drop accepted" "$?" "0"
check_eq "--diverged-rows=drop leaves the all-insertion output unchanged" \
	"$(cmp -s "$TMP/pav.tsv" "$TMP/pavdrop.tsv" && echo same || echo differ)" "same"
"$RB" chain --ref-prefix=B73 --lift="$LIFT" --diverged-rows=bogus -t1 "$IDX" "$DIR/queries.fa" \
	>/dev/null 2>"$TMP/badarg.log"
check_eq "--diverged-rows rejects an invalid value" "$([ $? -ne 0 ] && echo yes || echo no)" "yes"

# --- 6. chain --anchor-dist-npy: the same ternary + distance block refmap writes ---
# chain feeds the same PS4G/npy accumulator as refmap, so --anchor-dist-npy must give the same
# (3*n_gamete+2)-column layout and sidecars. A synthetic spliced read (B73[0,60) + B73[400,460),
# a 340 bp "intron") checks that each exon segment becomes its own row at its own start bin,
# carrying the read's chain-wide founder set.
B73_SEQ=$(awk '/^>B73_chr1/{f=1;next}/^>/{f=0}f' "$DIR/pangenome.fa" | tr -d '\n')
{ cat "$DIR/queries.fa"; printf '>spliced_query\n%s%s\n' "${B73_SEQ:0:60}" "${B73_SEQ:400:60}"; } > "$TMP/chain_q.fa"

"$RB" chain --ref-prefix=B73 --npy "$TMP/cad_nolift.npy" --anchor-dist-npy -t1 "$IDX" "$TMP/chain_q.fa" \
	>/dev/null 2>"$TMP/cad_nolift.log"
check_eq "chain --anchor-dist-npy without --lift exits non-zero" "$([ $? -ne 0 ] && echo yes || echo no)" "yes"
check_line_in "chain --anchor-dist-npy without --lift names the missing option (same error as refmap)" \
	"$TMP/cad_nolift.log" "ERROR: --anchor-dist-npy requires --lift=FILE"

"$RB" chain --ref-prefix=B73 --lift="$LIFT" --pav-grid=0 --ps4g "$TMP/cad_plain.ps4g" --npy "$TMP/cad_plain.npy" \
	-t1 "$IDX" "$TMP/chain_q.fa" > "$TMP/cad_plain.tsv" 2>/dev/null
check_eq "chain --npy (plain) exits 0" "$?" "0"
"$RB" chain --ref-prefix=B73 --lift="$LIFT" --pav-grid=0 --ps4g "$TMP/cad.ps4g" --npy "$TMP/cad.npy" \
	--anchor-dist-npy --anchor-dist-thresh=2000 -t1 "$IDX" "$TMP/chain_q.fa" > "$TMP/cad.tsv" 2>/dev/null
check_eq "chain --anchor-dist-npy exits 0" "$?" "0"
check_eq "chain --anchor-dist-npy leaves stdout unchanged" \
	"$(cmp -s "$TMP/cad_plain.tsv" "$TMP/cad.tsv" && echo same || echo differ)" "same"
check_eq "chain --anchor-dist-npy leaves the PS4G file unchanged (modulo #Command)" \
	"$(cmp -s <(grep -v '^#Command' "$TMP/cad_plain.ps4g") <(grep -v '^#Command' "$TMP/cad.ps4g") && echo same || echo differ)" "same"
check_eq "chain --anchor-dist-npy leaves bins.tsv unchanged" \
	"$(cmp -s "$TMP/cad_plain.npy.bins.tsv" "$TMP/cad.npy.bins.tsv" && echo same || echo differ)" "same"
check_line_in "chain spliced read: exon 1 row at its own start" "$TMP/cad.tsv" "$(printf 'spliced_query\tchr1\t0\t0,1,2,3,4,5,6')"
check_line_in "chain spliced read: exon 2 row at its own start" "$TMP/cad.tsv" "$(printf 'spliced_query\tchr1\t400\t0,1,2,3,4,5,6')"
cad_header=$(head -c 200 "$TMP/cad.npy" | tr -d '\0')
cad_shape=$(echo "$cad_header" | grep -oE "'shape': \([0-9]+, [0-9]+\)" | grep -oE '[0-9]+' | tr '\n' ' ')
check_eq "chain --anchor-dist-npy: 4 rows x 23 columns (3*7 gametes + 2 labels)" "$cad_shape" "4 23 "
check_line_in "chain --anchor-dist-npy writes the layout sidecar (ternary block)" "$TMP/cad.npy.layout.tsv" \
	"$(printf 'ternary\t7\t14\tmatch(1)/diverged(0)/deletion(-1) per founder')"
check_eq "plain chain --npy writes no layout sidecar" "$([ -e "$TMP/cad_plain.npy.layout.tsv" ] && echo yes || echo no)" "no"

"$RB" refmap --ref-prefix=B73 --lift="$LIFT" --npy "$TMP/rad.npy" --anchor-dist-npy -t1 "$IDX" "$TMP/chain_q.fa" \
	>/dev/null 2>/dev/null
check_eq "refmap --anchor-dist-npy exits 0 on the same input" "$?" "0"
check_eq "chain and refmap write the same layout sidecar" \
	"$(cmp -s "$TMP/rad.npy.layout.tsv" "$TMP/cad.npy.layout.tsv" && echo same || echo differ)" "same"

if [ -n "$NPY_PY" ]; then
	"$NPY_PY" - "$TMP/cad_plain.npy" "$TMP/cad.npy" "$TMP/rad.npy" <<'EOF'
import sys, numpy as np
plain, wide, rwide = (np.load(f) for f in sys.argv[1:4])
K = 7
assert wide.dtype == np.int32 and wide.shape == (4, 3 * K + 2), wide.shape
assert rwide.dtype == wide.dtype and rwide.shape[1] == wide.shape[1], (rwide.shape, wide.shape)
# the count block and the label columns are exactly the plain --npy output
assert (wide[:, :K] == plain[:, :K]).all() and (wide[:, -2:] == plain[:, -2:]).all()
for a in (wide, rwide):  # refmap's definition, applied cell by cell, holds for both commands
    cnt, tern, dist = a[:, :K], a[:, K:2*K], a[:, 2*K:3*K]
    assert set(np.unique(tern)) <= {-1, 0, 1}
    assert ((tern == 1) == (cnt > 0)).all(), "ternary match must be exactly the row's founder set"
    off = tern != 1
    assert ((tern == 0) == (off & (dist >= 0) & (dist <= 2000))).all()
    assert ((tern == -1) == (off & ((dist < 0) | (dist > 2000)))).all()
    assert (dist[:, 0] == 0).all(), "reference gamete distance is 0 by construction"
# rows: bin0 {B73} (ref_query), bin0 {1..6} (insertion reads), bin0 {all} and bin1 {all} (spliced read)
assert wide[:, K:2*K].tolist() == [[1,0,0,0,0,0,0], [0,1,1,1,1,1,1], [1]*7, [1]*7], wide[:, K:2*K].tolist()
assert wide[2, :K].tolist() == [1]*7 and wide[3, :K].tolist() == [1]*7, "one count per exon segment"
print("CHAIN_ANCHOR_DIST_OK")
EOF
	if [ $? -eq 0 ]; then pass; else fail "chain --anchor-dist-npy numpy check raised an assertion"; fi
else
	echo "SKIP: chain --anchor-dist-npy numpy check needs a python3 with numpy" >&2
fi

# --- 7. chain never links SMEMs on different reference contigs ---
# Two unrelated reference contigs: B73_chr1 (copied into Oh43) and B73_chr2 (copied into Mo17). The
# read is 80 bp of chr1 at 1000 followed by 80 bp of chr2 at 1020, so the two SMEMs sit at
# similar coordinates on different contigs. Linking them would emit one chr1 row whose founder set is
# {B73,Oh43} ∩ {B73,Mo17} = {B73}; the colinear DP must instead keep them apart, leaving the best
# single-contig chain. With this seed the chr2 SMEM extends 2 bp left into the chr1 half (the random
# sequences share those bases), so it is the longer one and wins: chr2 at 1018, B73 + Mo17 = 0,1.
awk 'BEGIN{srand(7); split("ACGT", b, ""); for (c = 1; c <= 2; ++c) { s = ""; for (i = 0; i < 3000; ++i) s = s b[int(rand()*4)+1]; seq[c] = s }
     print ">B73_chr1"; print seq[1]; print ">B73_chr2"; print seq[2]; print ">Oh43_chr1"; print seq[1]; print ">Mo17_chr2"; print seq[2]
     print substr(seq[1], 1001, 80) substr(seq[2], 1021, 80) > "'"$TMP"'/xcontig_query.seq"}' > "$TMP/xcontig.fa"
printf '>xcontig\n%s\n' "$(cat "$TMP/xcontig_query.seq")" > "$TMP/xcontig_query.fa"
"$RB" build -d -o "$TMP/xcontig.fmd" "$TMP/xcontig.fa" >/dev/null 2>"$TMP/xcontig_build.log" &&
	"$RB" ssa -s8 -o "$TMP/xcontig.fmd.ssa" "$TMP/xcontig.fmd" >/dev/null 2>>"$TMP/xcontig_build.log"
check_eq "cross-contig fixture builds" "$?" "0"
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}
     END{if(n)print n"\t"l}' "$TMP/xcontig.fa" | gzip > "$TMP/xcontig.fmd.len.gz"
"$RB" chain --ref-prefix=B73 -t1 "$TMP/xcontig.fmd" "$TMP/xcontig_query.fa" > "$TMP/xcontig.tsv" 2>"$TMP/xcontig.log"
check_eq "chain exits 0 (cross-contig read)" "$?" "0"
check_eq "chain: one row for the cross-contig read" "$(grep -c '^xcontig' "$TMP/xcontig.tsv")" "1"
check_line_in "chain: SMEMs on different contigs are not chained (best single-contig chain, B73 + Mo17)" \
	"$TMP/xcontig.tsv" "$(printf 'xcontig\tchr2\t1018\t0,1')"
check_eq "chain: no row with the cross-contig intersection {B73} alone" \
	"$(awk -F'\t' '$1=="xcontig" && $4=="0"' "$TMP/xcontig.tsv" | wc -l | tr -d ' ')" "0"

# --- 8. chain --support=best: founder set = the assemblies that best explain the read ---
# B73_chr1 and Mo17_chr1 are identical random sequence; Zd_chr1 is the same with SNPs at 1130 and
# 1160 -- a diverged allele. The donor read is Zd's [1100,1250) with one extra base change at read
# offset 110 (the donor is a different accession than Zd), so it matches no assembly end to end:
# SMEMs [0,110) in Zd only and [111,150) in all three (B73's own match [61,110) is contained in the
# Zd SMEM, so it is not super-maximal). The default rule keeps only reference-hitting SMEMs, so the
# read is called {B73,Mo17,Zd}; --support=best counts the Zd-only SMEM (150 bases for Zd vs 39) and
# calls Zd alone, at the same placement. A B73-allele read
# from the same region stays {B73,Mo17}; a read outside the diverged region stays {all three}.
awk 'BEGIN{srand(11); split("ACGT", b, ""); s = ""; for (i = 0; i < 3000; ++i) s = s b[int(rand()*4)+1]
     z = s; for (p = 1130; p <= 1160; p += 30) { c = substr(z, p+1, 1); n = (c == "A")? "C" : "A"; z = substr(z, 1, p) n substr(z, p+2) }
     print ">B73_chr1"; print s; print ">Mo17_chr1"; print s; print ">Zd_chr1"; print z
     d = substr(z, 1101, 150); c = substr(d, 111, 1); n = (c == "G")? "T" : "G"; d = substr(d, 1, 110) n substr(d, 112)
     print ">donor_read\n" d > "'"$TMP"'/sup_q.fa"
     print ">b73_read\n" substr(s, 1101, 150) >> "'"$TMP"'/sup_q.fa"
     print ">shared_read\n" substr(s, 2001, 150) >> "'"$TMP"'/sup_q.fa"}' > "$TMP/sup.fa"
"$RB" build -d -o "$TMP/sup.fmd" "$TMP/sup.fa" >/dev/null 2>"$TMP/sup_build.log" &&
	"$RB" ssa -s8 -o "$TMP/sup.fmd.ssa" "$TMP/sup.fmd" >/dev/null 2>>"$TMP/sup_build.log"
check_eq "support fixture builds" "$?" "0"
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}
     END{if(n)print n"\t"l}' "$TMP/sup.fa" | gzip > "$TMP/sup.fmd.len.gz"
"$RB" lift --ref-prefix=B73 -k 31 -s 20 -o "$TMP/sup.lift" "$TMP/sup.fmd" "$TMP/sup.fa" >/dev/null 2>"$TMP/sup_lift.log"
check_eq "support fixture lift builds" "$?" "0"
SUPC="chain --ref-prefix=B73 -l 19 --lift=$TMP/sup.lift -t1"
"$RB" $SUPC "$TMP/sup.fmd" "$TMP/sup_q.fa" > "$TMP/sup_default.tsv" 2>"$TMP/sup_default.log"
check_eq "chain (default support) exits 0" "$?" "0"
"$RB" $SUPC --support=chain "$TMP/sup.fmd" "$TMP/sup_q.fa" > "$TMP/sup_chain.tsv" 2>/dev/null
check_eq "--support=chain is byte-identical to the default" "$(cmp -s "$TMP/sup_default.tsv" "$TMP/sup_chain.tsv" && echo same)" "same"
check_eq "default rule: the donor read's set includes B73 (the bias --support=best fixes)" \
	"$(awk -F'\t' '$1=="donor_read"{print $4}' "$TMP/sup_default.tsv")" "0,1,2"
"$RB" $SUPC --support=best "$TMP/sup.fmd" "$TMP/sup_q.fa" > "$TMP/sup_best.tsv" 2>"$TMP/sup_best.log"
check_eq "chain --support=best exits 0" "$?" "0"
check_eq "--support=best: donor read -> Zd alone (gamete 2), same placement as default" \
	"$(awk -F'\t' '$1=="donor_read"{print $2":"$3":"$4}' "$TMP/sup_best.tsv")" \
	"$(awk -F'\t' '$1=="donor_read"{print $2":"$3}' "$TMP/sup_default.tsv"):2"
check_eq "--support=best: B73-allele read -> B73 + Mo17" \
	"$(awk -F'\t' '$1=="b73_read"{print $4}' "$TMP/sup_best.tsv")" "0,1"
check_eq "--support=best: read in shared sequence -> all three" \
	"$(awk -F'\t' '$1=="shared_read"{print $4}' "$TMP/sup_best.tsv")" "0,1,2"
"$RB" $SUPC --support=best --support-tol=200 "$TMP/sup.fmd" "$TMP/sup_q.fa" > "$TMP/sup_tol.tsv" 2>/dev/null
check_eq "--support-tol larger than the read keeps every assembly with any support" \
	"$(awk -F'\t' '$1=="donor_read"{print $4}' "$TMP/sup_tol.tsv")" "0,1,2"
"$RB" $SUPC --support=bogus "$TMP/sup.fmd" "$TMP/sup_q.fa" >/dev/null 2>&1
check_eq "--support=bogus exits non-zero" "$([ $? -ne 0 ] && echo yes || echo no)" "yes"

# --- 9. chain output does not depend on the suffix-array sampling rate (hit order) ---
# rb3_ssa_multi returns an interval's occurrences in an order that depends on `ssa -s`; chain keeps
# the first reference hit and breaks ties by order, so the hits are sorted before use. Same index,
# SA sampled at 1/256 and 1/2: chain output (stdout, PS4G and npy) must be byte-identical.
cp "$TMP/pangenome.fmd" "$TMP/ss2.fmd"; cp "$IDX.len.gz" "$TMP/ss2.fmd.len.gz"
"$RB" ssa -s1 -o "$TMP/ss2.fmd.ssa" "$TMP/ss2.fmd" >/dev/null 2>"$TMP/ss2.log"
for ix in "$IDX" "$TMP/ss2.fmd"; do
	tag=$(basename "$ix" .fmd)
	"$RB" chain --ref-prefix=B73 -l 19 --lift="$LIFT" --support=best --ps4g "$TMP/ho_$tag.ps4g" --npy "$TMP/ho_$tag.npy" -t1 \
		"$ix" "$DIR/queries.fa" > "$TMP/ho_$tag.tsv" 2>/dev/null
done
check_eq "chain stdout identical across SA sampling rates" "$(cmp -s "$TMP/ho_pangenome.tsv" "$TMP/ho_ss2.tsv" && echo same)" "same"
check_eq "chain npy identical across SA sampling rates" "$(cmp -s "$TMP/ho_pangenome.npy" "$TMP/ho_ss2.npy" && echo same)" "same"
check_eq "chain PS4G identical across SA sampling rates (minus the command line)" \
	"$(diff <(grep -v '^#Command' "$TMP/ho_pangenome.ps4g") <(grep -v '^#Command' "$TMP/ho_ss2.ps4g") >/dev/null && echo same)" "same"

# --- 10. chain --dup-count: a collapsed duplicate read named <name>;x=N counts N times ---
# Every query three times (all reads) versus once each named ";x=3": identical PS4G and npy.
awk '/^>/{n=substr($1,2); next}{s[n]=s[n] $0; o[++k]=n} END{for (i=1;i<=k;++i) if (!(o[i] in d)) {d[o[i]]=1; for (r=0;r<3;++r) print ">" o[i] "_" r "\n" s[o[i]]; print ">" o[i] ";x=3\n" s[o[i]] > "'"$TMP"'/dup_collapsed.fa"}}' \
	"$DIR/queries.fa" > "$TMP/dup_full.fa"
"$RB" chain --ref-prefix=B73 -l 19 --lift="$LIFT" --insertion-rows=ordinary --ps4g "$TMP/dup_full.ps4g" --npy "$TMP/dup_full.npy" --anchor-dist-npy -t1 \
	"$IDX" "$TMP/dup_full.fa" > /dev/null 2>&1
"$RB" chain --ref-prefix=B73 -l 19 --lift="$LIFT" --insertion-rows=ordinary --dup-count --ps4g "$TMP/dup_coll.ps4g" --npy "$TMP/dup_coll.npy" --anchor-dist-npy -t1 \
	"$IDX" "$TMP/dup_collapsed.fa" > /dev/null 2>&1
check_eq "--dup-count: collapsed reads give the same npy as the duplicated input" \
	"$(cmp -s "$TMP/dup_full.npy" "$TMP/dup_coll.npy" && echo same)" "same"
check_eq "--dup-count: same PS4G counts (minus the command line)" \
	"$(diff <(grep -v '^#Command' "$TMP/dup_full.ps4g") <(grep -v '^#Command' "$TMP/dup_coll.ps4g") >/dev/null && echo same)" "same"
"$RB" chain --ref-prefix=B73 -l 19 --lift="$LIFT" --insertion-rows=ordinary --ps4g "$TMP/dup_nocount.ps4g" -t1 \
	"$IDX" "$TMP/dup_collapsed.fa" > /dev/null 2>&1
check_eq "without --dup-count the ;x=N suffix is ignored (counts differ from the duplicated input)" \
	"$(diff <(grep -v '^#Command' "$TMP/dup_full.ps4g") <(grep -v '^#Command' "$TMP/dup_nocount.ps4g") >/dev/null && echo same || echo differ)" "differ"

# --- 11. refmap --paired: mate pairs -> one PS4G/npy event per pair ---
# Fixture docs/examples/pe/ (generated by make_pe_fixture.py; see its docstring for the design):
# 4 genomes (gametes B73=0 Ki3=1 Mo17=2 Oh43=3), 40 clean pairs for the insert-size estimate and
# seven hand-built pairs P1..P7 covering CONC (with a founder-set intersection), SOLO, DISC, RESC,
# a reference-copy tie-break (P6) and --pe-trust-exact (P7).
PE="$DIR/pe"
PIDX="$TMP/pe.fmd"; PLIFT="$TMP/pe.lift"
"$RB" build -d -o "$PIDX" "$PE/pangenome.fa" >/dev/null 2>"$TMP/pe_build.log"
"$RB" ssa -s8 -o "$PIDX.ssa" "$PIDX" >/dev/null 2>>"$TMP/pe_build.log"
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)} END{if(n)print n"\t"l}' "$PE/pangenome.fa" | gzip > "$PIDX.len.gz"
"$RB" lift --ref-prefix=B73 -k 31 -s 50 -o "$PLIFT" "$PIDX" "$PE/pangenome.fa" >/dev/null 2>>"$TMP/pe_build.log"
check_eq "paired fixture: index + lift build" "$?" "0"
PRM="refmap --ref-prefix=B73 --lift=$PLIFT -t2"

"$RB" $PRM --paired --ps4g="$TMP/pe.ps4g" --npy="$TMP/pe.npy" --ps4g-per-read="$TMP/pe.pr" "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > "$TMP/pe.tsv" 2>"$TMP/pe.log"
check_eq "--paired (R1 R2) exits 0" "$?" "0"
check_eq "--paired: one row per mate (47 pairs)" "$(wc -l < "$TMP/pe.tsv" | tr -d ' ')" "94"
check_eq "--paired: every row gains exactly one trailing column (12)" \
	"$(awk -F'\t' '{print NF}' "$TMP/pe.tsv" | sort -u | tr '\n' ' ')" "12 "
pe_class() { awk -F'\t' -v n="$1" '$1==n{print $3"/"$NF}' "$TMP/pe.tsv"; }
check_eq "P1: projected mate + exact mate within the insert -> CONC" "$(pe_class P1/1) $(pe_class P1/2)" "PLACED/CONC EXACT/CONC"
check_eq "P2: two exact mates -> CONC" "$(pe_class P2/1) $(pe_class P2/2)" "EXACT/CONC EXACT/CONC"
check_eq "P3: one mate unplaceable -> SOLO" "$(pe_class P3/1) $(pe_class P3/2)" "EXACT/SOLO UNPLACED/SOLO"
check_eq "P4: exact mates 17 kb apart -> DISC" "$(pe_class P4/1) $(pe_class P4/2)" "EXACT/DISC EXACT/DISC"
check_eq "P5: projected mate rescued next to its exact mate -> MOVED/RESC" "$(pe_class P5/1) $(pe_class P5/2)" "PLACED/MOVED EXACT/RESC"
check_eq "P5: the moved mate's position is Mo17's projection (9200), not Ki3's TE copy (16100)" \
	"$(awk -F'\t' '$1=="P5/1"{print $8}' "$TMP/pe.tsv")" "9200"
check_eq "--paired: all 40 clean pairs CONC" "$(grep -c '^iso.*	CONC$' "$TMP/pe.tsv")" "80"
check_eq "P6: exact mate with two reference copies ends at the copy next to its mate (11000, not 19000)" \
	"$(awk -F'\t' '$1=="P6/1"{print $8}' "$TMP/pe.tsv")" "11000"
check_eq "P7: projected mate with nothing near its exact mate -> DISC by default" "$(pe_class P7/1) $(pe_class P7/2)" "PLACED/DISC EXACT/DISC"
if grep -q 'insert size from 4[12] pairs: mean 399.5' "$TMP/pe.log"; then pass; else fail "--paired: insert size estimate (mean 399.5) not logged"; fi
check_eq "--paired: summary = 44 CONC+RESC (P6 lands in either), 1 SOLO, 2 DISC, 0 NONE" \
	"$(sed -n 's/.*] 47 pairs: \([0-9]*\) CONC, \([0-9]*\) RESC, \([0-9]*\) SOLO, \([0-9]*\) DISC, \([0-9]*\) NONE.*/\1+\2 \3 \4 \5/p' "$TMP/pe.log" | awk '{split($1,a,"+"); print a[1]+a[2], $2, $3, $4}')" "44 1 2 0"

# PS4G: one event per non-discordant pair, two per discordant one -> 40 + 5 + 2*2 = 49 counts
check_eq "--paired ps4g: #TotalUniqueCounts = 40 clean + P1,P2,P3,P5,P6 + 2 each for discordant P4, P7 = 49" \
	"$(awk '/^#TotalUniqueCounts:/{print $2}' "$TMP/pe.ps4g")" "49"
check_eq "--paired ps4g: data-row counts sum to the total" \
	"$(awk -F'\t' 'f{s+=$4} /^gameteSet/{f=1} END{print s}' "$TMP/pe.ps4g")" "49"
check_line_in "P1 event: founders both mates share = {Mo17} only, at the exact mate (10300 -> bin 40)" "$TMP/pe.ps4g" "$(printf '2\tchr1\t40\t1')"
check_line_in "P5 event: rescued mate keeps only Mo17 (Ki3's copy is elsewhere), bin 37" "$TMP/pe.ps4g" "$(printf '2\tchr1\t37\t1')"
check_line_in "per-read: P1 is one line named after the fragment" "$TMP/pe.pr" "$(printf 'P1\tchr1\t10300\t2')"
check_line_in "per-read: P3 (SOLO) uses the placed mate's set" "$TMP/pe.pr" "$(printf 'P3\tchr1\t14000\t0,1,2,3')"
check_eq "per-read: discordant P4 gives one line per mate" "$(grep -c '^P4/[12]	' "$TMP/pe.pr")" "2"
check_eq "per-read: 49 events (+ header)" "$(wc -l < "$TMP/pe.pr" | tr -d ' ')" "50"

# single-end run of the same reads, for contrast: P1 then gives two separate events
"$RB" $PRM --ps4g="$TMP/pe_se.ps4g" "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > "$TMP/pe_se.tsv" 2>/dev/null
check_eq "without --paired: 11 columns, files read one after another" \
	"$(awk -F'\t' '{print NF}' "$TMP/pe_se.tsv" | sort -u | tr '\n' ' ')/$(wc -l < "$TMP/pe_se.tsv" | tr -d ' ')" "11 /94"
check_eq "without --paired: no P1 event has the Mo17-only set at bin 40" \
	"$(grep -c "^2	chr1	40	" "$TMP/pe_se.ps4g")" "0"
check_eq "without --paired: one event per placed read (93 = 94 reads minus P3's junk mate)" \
	"$(awk '/^#TotalUniqueCounts:/{print $2}' "$TMP/pe_se.ps4g")/$(awk -F'\t' '$3=="PLACED"||$3=="EXACT"' "$TMP/pe_se.tsv" | wc -l | tr -d ' ')" "93/93"

# interleaved input == two files
"$RB" $PRM --paired --ps4g="$TMP/pe_il.ps4g" --npy="$TMP/pe_il.npy" "$PIDX" "$PE/interleaved.fa" > "$TMP/pe_il.tsv" 2>/dev/null
check_eq "--paired interleaved: same stdout as R1 R2" "$(cmp -s "$TMP/pe.tsv" "$TMP/pe_il.tsv" && echo same)" "same"
check_eq "--paired interleaved: same PS4G (minus the command line)" \
	"$(diff <(grep -v '^#Command' "$TMP/pe.ps4g") <(grep -v '^#Command' "$TMP/pe_il.ps4g") >/dev/null && echo same)" "same"
check_eq "--paired interleaved: same npy" "$(cmp -s "$TMP/pe.npy" "$TMP/pe_il.npy" && echo same)" "same"
# small batches (-K) must not split a pair; output identical
"$RB" $PRM --paired -K 250 "$PIDX" "$PE/interleaved.fa" > "$TMP/pe_k.tsv" 2>"$TMP/pe_k.log"
check_eq "--paired -K 250 (pair-sized batches): same stdout classes" \
	"$(cut -f1-11 "$TMP/pe_k.tsv" | cmp -s - <(cut -f1-11 "$TMP/pe.tsv") && echo same)" "same"
if grep -q 'WARNING: only' "$TMP/pe_k.log"; then pass; else fail "--paired: tiny batches warn that the insert size cannot be estimated"; fi

# --pe-trust-exact: P7's projected mate is moved onto its exact mate; nothing else changes
"$RB" $PRM --paired --pe-trust-exact --ps4g="$TMP/pe_te.ps4g" "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > "$TMP/pe_te.tsv" 2>/dev/null
check_eq "--pe-trust-exact: P7 -> PLACED/MOVED at 6000 + EXACT/RESC" \
	"$(awk -F'\t' '$1=="P7/1"{print $3"/"$NF"@"$8} $1=="P7/2"{print $3"/"$NF}' "$TMP/pe_te.tsv" | tr '\n' ' ')" "PLACED/MOVED@6000 EXACT/RESC "
check_eq "--pe-trust-exact: every other row unchanged" \
	"$(diff <(grep -v '^P7/' "$TMP/pe.tsv") <(grep -v '^P7/' "$TMP/pe_te.tsv") >/dev/null && echo same)" "same"
check_eq "--pe-trust-exact: P7 is now one event (48)" "$(awk '/^#TotalUniqueCounts:/{print $2}' "$TMP/pe_te.ps4g")" "48"
pe_err_te() { "$RB" $PRM --pe-trust-exact "$PIDX" "$PE/r1.fa" >/dev/null 2>"$TMP/pe_te_err.log"; }
pe_err_te; if [ $? -ne 0 ] && grep -q 'need --paired' "$TMP/pe_te_err.log"; then pass; else fail "--pe-trust-exact without --paired must fail"; fi

# fixed insert size
"$RB" $PRM --paired --pe-isize=400,10 "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > "$TMP/pe_fix.tsv" 2>"$TMP/pe_fix.log"
check_eq "--pe-isize: exits 0, same classes" "$(cmp -s "$TMP/pe.tsv" "$TMP/pe_fix.tsv" && echo same)" "same"
if grep -q 'insert size (fixed) mean 400.0 sd 10.0 max 440' "$TMP/pe_fix.log"; then pass; else fail "--pe-isize: fixed insert size not reported"; fi
"$RB" $PRM --paired --pe-slack=0 --pe-isize=400,10 "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > "$TMP/pe_s0.tsv" 2>/dev/null
check_eq "--pe-slack=0: P1's projected mate (100 bp off) is still concordant at max 440" "$(awk -F'\t' '$1=="P1/1"{print $NF}' "$TMP/pe_s0.tsv")" "CONC"

# input errors (each must exit non-zero with a message, before any output)
pe_err() { # pe_err <desc> <expected message fragment> args...
	local d="$1" m="$2"; shift 2
	"$RB" "$@" > "$TMP/pe_err.out" 2> "$TMP/pe_err.log"
	if [ $? -ne 0 ] && grep -q -- "$m" "$TMP/pe_err.log"; then pass; else fail "$d (want non-zero exit + '$m')"; cat "$TMP/pe_err.log" >&2; fi
}
printf '>x/2\nACGT\n>y/2\nACGT\n' > "$TMP/r2_badname.fa"; head -4 "$PE/r1.fa" > "$TMP/r1_two.fa"
pe_err "--paired: mate names differ" "are not mates" $PRM --paired "$PIDX" "$TMP/r1_two.fa" "$TMP/r2_badname.fa"
head -6 "$PE/r2.fa" > "$TMP/r2_three.fa"
pe_err "--paired: R2 longer than R1" "second file has more reads" $PRM --paired "$PIDX" "$TMP/r1_two.fa" "$TMP/r2_three.fa"
pe_err "--paired: R1 longer than R2" "first file has more reads" $PRM --paired "$PIDX" "$TMP/r2_three.fa" "$TMP/r1_two.fa"
head -6 "$PE/interleaved.fa" > "$TMP/il_odd.fa"
pe_err "--paired: odd interleaved file" "odd number of reads" $PRM --paired "$PIDX" "$TMP/il_odd.fa"
pe_err "--paired: three input files" "R1.fq> <R2.fq> or one interleaved" $PRM --paired "$PIDX" "$PE/r1.fa" "$PE/r2.fa" "$PE/r1.fa"
pe_err "--paired: chain not implemented yet" "refmap and mem only" chain --ref-prefix=B73 --paired "$PIDX" "$PE/r1.fa" "$PE/r2.fa"
pe_err "--paired: -L has no read names" "not -L" $PRM --paired -L "$PIDX" "$PE/interleaved.fa"
pe_err "--pe-isize without --paired" "need --paired" $PRM --pe-isize=400,10 "$PIDX" "$PE/r1.fa"
pe_err "--pe-isize malformed" "MEAN,SD" $PRM --paired --pe-isize=400 "$PIDX" "$PE/r1.fa" "$PE/r2.fa"

if [ -n "$NPY_PY" ]; then
	if "$NPY_PY" verify_ps4g_npy.py "$TMP/pe.ps4g" "$TMP/pe.npy" > "$TMP/pe_verify.log" 2>&1; then pass
	else fail "verify_ps4g_npy.py: --paired PS4G and npy disagree"; cat "$TMP/pe_verify.log" >&2; fi
fi

# --- 12. mem --paired: per-assembly concordance of the mates' SMEM hits ---
# Same fixture. A pair is CONC when both mates hit one contig on opposite strands with implied
# read starts within the insert; only those positions are printed. Fixture spans differ per
# genome, so expectations below are per contig (see make_pe_fixture.py for coordinates).
"$RB" mem --paired -p 64 -t2 "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > "$TMP/mem_pe.tsv" 2>"$TMP/mem_pe.log"
check_eq "mem --paired exits 0" "$?" "0"
mem_cls() { awk -F'\t' -v n="$1" '$1==n{print $NF}' "$TMP/mem_pe.tsv" | sort -u | tr '\n' ' '; }
mem_pos() { awk -F'\t' -v n="$1" '$1==n{for(i=6;i<NF;++i) print $i}' "$TMP/mem_pe.tsv" | sort | tr '\n' ' '; }
check_eq "mem P2: CONC on all four genomes" "$(mem_cls P2/1)$(mem_cls P2/2)" "CONC CONC "
check_eq "mem P5: CONC, only Mo17 carries both mates (Ki3's TE copy is elsewhere)" "$(mem_pos P5/1)/$(mem_pos P5/2)" "Mo17_chr1:+:9600 /Mo17_chr1:-:10000 "
check_eq "mem P6: the duplicated mate keeps only the copy next to its mate (4 of 8 positions)" \
	"$(awk -F'\t' '$1=="P6/1"{print $4"/"$5}' "$TMP/mem_pe.tsv")" "8/4"
check_eq "mem P6: kept positions are the 11000 copy" "$(mem_pos P6/1)" "B73_chr1:+:11000 Ki3_chr1:+:12000 Mo17_chr1:+:12400 Oh43_chr1:+:11000 "
check_eq "mem P1: mates 900 bp apart on Mo17 (> max insert) -> DISC, positions unfiltered" "$(mem_cls P1/1)$(awk -F'\t' '$1=="P1/2"{print $5}' "$TMP/mem_pe.tsv")" "DISC 3"
check_eq "mem P3: junk mate has no SMEM -> SOLO" "$(mem_cls P3/1)/$(grep -c '^P3/2' "$TMP/mem_pe.tsv")" "SOLO /0"
check_eq "mem P4, P7: DISC" "$(mem_cls P4/1)$(mem_cls P7/2)" "DISC DISC "
check_eq "mem: 40 clean pairs all CONC" "$(awk -F'\t' '$1~/^iso/ && $NF!="CONC"' "$TMP/mem_pe.tsv" | wc -l | tr -d ' ')" "0"
if grep -q '47 pairs: 43 CONC, 0 RESC, 1 SOLO, 3 DISC, 0 NONE' "$TMP/mem_pe.log"; then pass; else fail "mem --paired: pair-class summary"; fi
if grep -q 'insert size from 42 pairs: mean 400.7' "$TMP/mem_pe.log"; then pass; else fail "mem --paired: insert size from whole-read SMEM pairs"; fi
"$RB" mem --paired -p 64 -t2 "$PIDX" "$PE/interleaved.fa" > "$TMP/mem_il.tsv" 2>/dev/null
check_eq "mem --paired: interleaved == R1 R2" "$(cmp -s "$TMP/mem_pe.tsv" "$TMP/mem_il.tsv" && echo same)" "same"
"$RB" mem -p 64 -t2 "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > "$TMP/mem_se.tsv" 2>/dev/null
check_eq "mem without --paired: the same SMEM rows (files read one after another, so sort)" \
	"$(cut -f1-4 "$TMP/mem_se.tsv" | sort | cmp -s - <(cut -f1-4 "$TMP/mem_pe.tsv" | sort) && echo same)" "same"
check_eq "mem without --paired: P6/1 lists both copies (8 positions)" "$(awk -F'\t' '$1=="P6/1"{print $5}' "$TMP/mem_se.tsv")" "8"
"$RB" mem --paired --pe-isize=400,10 -p 64 -t2 "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > /dev/null 2>"$TMP/mem_fix.log"
if grep -q 'insert size (fixed) mean 400.0 sd 10.0 max 440 (+0 slack)' "$TMP/mem_fix.log"; then pass; else fail "mem --pe-isize: fixed insert, no slack by default"; fi
# -p 2 locates only 2 of each SMEM's 4-8 occurrences: no mate is fully located, so nothing is filtered
"$RB" mem --paired -p 2 -t2 "$PIDX" "$PE/r1.fa" "$PE/r2.fa" > "$TMP/mem_p2.tsv" 2>/dev/null
check_eq "mem --paired -p 2: sampled mates cannot veto -- every row keeps all its located positions" \
	"$(awk -F'\t' 'NF>5 && $5!=NF-6' "$TMP/mem_p2.tsv" | wc -l | tr -d ' ')/$([ "$(grep -c 'CONC$' "$TMP/mem_p2.tsv")" -gt 0 ] && echo conc)" "0/conc"
pe_err "mem --paired without -p" "needs -p INT" mem --paired "$PIDX" "$PE/r1.fa" "$PE/r2.fa"
pe_err "mem --paired with --cov" "cannot be combined" mem --paired --cov -p 64 "$PIDX" "$PE/r1.fa" "$PE/r2.fa"
pe_err "mem --pe-trust-exact" "refmap only" mem --paired --pe-trust-exact -p 64 "$PIDX" "$PE/r1.fa" "$PE/r2.fa"
pe_err "sw --paired not implemented yet" "refmap and mem only" sw --paired "$PIDX" "$PE/r1.fa" "$PE/r2.fa"

# --- memory safety, if valgrind is available ---
if command -v valgrind >/dev/null 2>&1; then
	valgrind --error-exitcode=99 --leak-check=full -q \
		"$RB" refmap --ref-prefix=B73 --walk --ps4g "$TMP/vg.ps4g" --npy "$TMP/vg.npy" --label-bed "$TMP/labels.bed" -t1 "$IDX" "$DIR/queries.fa" \
		> /dev/null 2>"$TMP/valgrind.log"
	if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks; see $TMP/valgrind.log"; cat "$TMP/valgrind.log" >&2; fi

	# --target-hits with real multi-threaded compute (-t4), exercising the
	# atomic hit counter under actual concurrent access.
	valgrind --error-exitcode=99 --leak-check=full -q \
		"$RB" refmap --ref-prefix=B73 --walk --target-hits=5 -K 10 -t4 "$IDX" "$TMP/target_hits.fa" \
		> /dev/null 2>"$TMP/valgrind_target.log"
	if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks on --target-hits; see $TMP/valgrind_target.log"; cat "$TMP/valgrind_target.log" >&2; fi

	# chain --support=best with --lift and the npy accumulator (best-support masks, lift projections)
	valgrind --error-exitcode=99 --leak-check=full -q \
		"$RB" chain --ref-prefix=B73 -l 19 --lift="$TMP/sup.lift" --support=best --npy "$TMP/vg_sup.npy" --anchor-dist-npy -t2 \
		"$TMP/sup.fmd" "$TMP/sup_q.fa" > /dev/null 2>"$TMP/valgrind_sup.log"
	if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks on chain --support=best; see $TMP/valgrind_sup.log"; cat "$TMP/valgrind_sup.log" >&2; fi

	# --paired: projections, rescue, pair events and the per-read file
	valgrind --error-exitcode=99 --leak-check=full -q \
		"$RB" refmap --ref-prefix=B73 --lift="$TMP/pe.lift" --paired --ps4g "$TMP/vg_pe.ps4g" --npy "$TMP/vg_pe.npy" --ps4g-per-read "$TMP/vg_pe.pr" -K 2000 -t2 \
		"$TMP/pe.fmd" "$DIR/pe/r1.fa" "$DIR/pe/r2.fa" > /dev/null 2>"$TMP/valgrind_pe.log"
	if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks on refmap --paired; see $TMP/valgrind_pe.log"; cat "$TMP/valgrind_pe.log" >&2; fi
	valgrind --error-exitcode=99 --leak-check=full -q \
		"$RB" mem --paired -p 64 -K 2000 -t2 "$TMP/pe.fmd" "$DIR/pe/r1.fa" "$DIR/pe/r2.fa" > /dev/null 2>"$TMP/valgrind_mem_pe.log"
	if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks on mem --paired; see $TMP/valgrind_mem_pe.log"; cat "$TMP/valgrind_mem_pe.log" >&2; fi
else
	echo "SKIP: valgrind not installed; skipping the memory-safety pass" >&2
fi

echo "run_integration: $n_pass passed, $n_fail failed"
[ "$n_fail" -eq 0 ]
