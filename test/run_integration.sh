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

# --- 7. isa + extract: coordinate-based sequence extraction ---
# Truth for the tiny fixture comes straight from the FASTA (awk substr), so these checks
# need no python. Coordinates are 0-based half-open; ':-' is the reverse complement.
fa_seq() { # fa_seq <fasta> <name> -> the whole sequence on one line
	awk -v n="$2" '/^>/{p=(substr($1,2)==n); next} p{printf "%s", $0} END{print ""}' "$1"
}
fa_sub() { # fa_sub <fasta> <name> <start> <end> [-]
	local s; s=$(fa_seq "$1" "$2" | awk -v a="$3" -v b="$4" '{print substr($0, a+1, b-a)}')
	if [ "${5:-+}" = "-" ]; then echo "$s" | rev | tr ACGTN TGCAN; else echo "$s"; fi
}
fa_rec() { # fa_rec <fasta-out> <qid> -> the record's sequence joined onto one line
	awk -v q="$2" '/^>/{p=(substr($1,2)==q); next} p{printf "%s", $0} END{print ""}' "$1"
}
FA="$DIR/pangenome.fa"

"$RB" isa -b 3 "$IDX" 2>"$TMP/isa.log"
check_eq "isa exits 0" "$?" "0"
check_eq "isa writes <idx>.isa by default" "$([ -s "$IDX.isa" ] && echo yes || echo no)" "yes"
"$RB" isa 2>/dev/null; check_eq "isa without an index prints usage and exits 1" "$?" "1"

REGS="B73_chr1:0-1 B73_chr1:599-600 B73_chr1:0-600 Oh43_chr1:100-1100 Oh43_chr1:1099-1100 Mo17_chr1:250-260:- Tx303_chr1:1099-1100:- Ki3_chr1:0-1100:- CML247_chr1:290-310:+ B97_chr1:512-777"
"$RB" extract "$IDX" $REGS > "$TMP/ext.fa" 2>"$TMP/ext.log"
check_eq "extract exits 0 (direct regions)" "$?" "0"
check_eq "extract: one record per region" "$(grep -c '^>' "$TMP/ext.fa")" "10"
nbad=0
for r in $REGS; do
	strand=+; reg=$r
	case "$r" in *:-) strand=-; reg=${r%:-} ;; *:+) reg=${r%:+} ;; esac
	ctg=${reg%:*}; rng=${reg##*:}
	want=$(fa_sub "$FA" "$ctg" "${rng%-*}" "${rng#*-}" "$strand")
	got=$(fa_rec "$TMP/ext.fa" "$r")
	[ "$got" = "$want" ] || { nbad=$((nbad+1)); echo "  mismatch $r: got $got want $want" >&2; }
done
check_eq "extract: every region equals the FASTA substring (revcomp for '-')" "$nbad" "0"
check_eq "extract: '-' of a full sequence is its reverse complement" "$(fa_rec "$TMP/ext.fa" Ki3_chr1:0-1100:-)" "$(fa_seq "$FA" Ki3_chr1 | rev | tr ACGT TGCA)"
check_eq "extract: '+' suffix is accepted and equals no suffix" "$(fa_rec "$TMP/ext.fa" CML247_chr1:290-310:+)" "$(fa_sub "$FA" CML247_chr1 290 310)"

# errors: reported per query on stderr, other queries still printed, exit status 1
"$RB" extract "$IDX" B73_chr1:0-5 Nope_chr1:0-5 B73_chr1:0-601 B73_chr1:10-10 B73_chr1:10-5 B73_chr1:abc B73_chr1:5-10:x B73_chr1:1-3 \
	> "$TMP/ext_err.fa" 2>"$TMP/ext_err.log"
check_eq "extract: input errors -> exit status 1" "$?" "1"
check_eq "extract: valid queries still written alongside errors" "$(grep -c '^>' "$TMP/ext_err.fa")" "2"
check_eq "extract: unknown contig reported" "$(grep -c "query 'Nope_chr1:0-5': unknown contig" "$TMP/ext_err.log")" "1"
check_eq "extract: past-the-end reported (no silent clamp)" "$(grep -c "query 'B73_chr1:0-601': interval extends past the contig end" "$TMP/ext_err.log")" "1"
check_eq "extract: start == end reported" "$(grep -c "'B73_chr1:10-10': empty interval (start == end)" "$TMP/ext_err.log")" "1"
check_eq "extract: start > end reported" "$(grep -c "'B73_chr1:10-5': empty interval (start > end)" "$TMP/ext_err.log")" "1"
check_eq "extract: unparsable regions reported" "$(grep -c 'cannot parse region' "$TMP/ext_err.log")" "2"
"$RB" extract --clamp "$IDX" B73_chr1:590-700 B73_chr1:700-800 > "$TMP/ext_clamp.fa" 2>"$TMP/ext_clamp.log"
check_eq "extract --clamp: an interval that is empty after clamping is still an error" "$?" "1"
check_line_in "extract --clamp: header marks a clamped interval" "$TMP/ext_clamp.fa" ">B73_chr1:590-700 clamped"
check_eq "extract --clamp: clipped to the contig end" "$(fa_rec "$TMP/ext_clamp.fa" B73_chr1:590-700)" "$(fa_sub "$FA" B73_chr1 590 600)"
"$RB" extract "$IDX" 2>/dev/null >/dev/null; check_eq "extract with no queries -> exit 1" "$?" "1"

# TSV input: header, comment, strand column, bad strand
printf 'qid\tcontig\tstart\tend\tstrand\n# a comment\nq1\tOh43_chr1\t295\t810\t-\nq2\tB73_chr1\t5\t9\t+\nq3\tB73_chr1\t5\t9\n' > "$TMP/q.tsv"
"$RB" extract -f "$TMP/q.tsv" "$IDX" > "$TMP/ext_tsv.fa" 2>/dev/null
check_eq "extract -f: exits 0" "$?" "0"
check_eq "extract -f: header line and comment skipped, 3 records" "$(grep -c '^>' "$TMP/ext_tsv.fa")" "3"
check_eq "extract -f: '-' strand column" "$(fa_rec "$TMP/ext_tsv.fa" q1)" "$(fa_sub "$FA" Oh43_chr1 295 810 -)"
check_eq "extract -f: strand column optional (+)" "$(fa_rec "$TMP/ext_tsv.fa" q3)" "$(fa_sub "$FA" B73_chr1 5 9)"
printf 'q1\tB73_chr1\t5\t9\t*\nq2\tB73_chr1\tx\t9\t+\n' | "$RB" extract -f - "$IDX" > /dev/null 2>"$TMP/ext_tsv_bad.log"
check_eq "extract -f -: bad strand / non-numeric start after line 1 -> exit 1" "$?" "1"
check_eq "extract -f: bad strand reported" "$(grep -c "strand column must be" "$TMP/ext_tsv_bad.log")" "1"

# many random queries: -t4 == -t1, line wrapping, in-memory ISA == .isa file, -I
awk 'BEGIN{srand(7); split("B73_chr1 Oh43_chr1 Mo17_chr1 Ki3_chr1",c," "); L["B73_chr1"]=600;
	for(i=0;i<300;i++){ n=c[int(rand()*4)+1]; len=(n=="B73_chr1")?600:1100; a=int(rand()*len); b=a+1+int(rand()*(len-a));
	print "r"i"\t"n"\t"a"\t"b"\t"(rand()<.5?"+":"-") }}' > "$TMP/rand.tsv"
"$RB" extract -t1 -f "$TMP/rand.tsv" "$IDX" > "$TMP/rand1.fa" 2>/dev/null
"$RB" extract -t4 -f "$TMP/rand.tsv" "$IDX" > "$TMP/rand4.fa" 2>/dev/null
check_eq "extract: -t4 output identical to -t1 (300 random queries)" "$(cmp -s "$TMP/rand1.fa" "$TMP/rand4.fa" && echo same)" "same"
nbad=$(awk -F'\t' '{print $1, $2, $3, $4, $5}' "$TMP/rand.tsv" | head -60 | while read q c a b s; do
	[ "$(fa_rec "$TMP/rand1.fa" "$q")" = "$(fa_sub "$FA" "$c" "$a" "$b" "$s")" ] || echo bad; done | wc -l | tr -d ' ')
check_eq "extract: random queries equal the FASTA (first 60 checked)" "$nbad" "0"
"$RB" extract -w 7 -f "$TMP/rand.tsv" "$IDX" > "$TMP/randw.fa" 2>/dev/null
check_eq "extract -w 7: no sequence line longer than 7" "$(awk '!/^>/ && length($0)>7' "$TMP/randw.fa" | wc -l | tr -d ' ')" "0"
check_eq "extract -w 7: same records once unwrapped" "$(awk '/^>/{if(s)print s; print; s=""; next}{s=s $0} END{print s}' "$TMP/randw.fa" | md5sum)" "$(md5sum < "$TMP/rand1.fa")"
mv "$IDX.isa" "$TMP/moved.isa"
"$RB" extract -f "$TMP/rand.tsv" "$IDX" > "$TMP/rand_lazy.fa" 2>/dev/null
check_eq "extract: no .isa -> entry points derived from the .ssa in memory, same output" "$(cmp -s "$TMP/rand1.fa" "$TMP/rand_lazy.fa" && echo same)" "same"
"$RB" extract -I "$TMP/moved.isa" -f "$TMP/rand.tsv" "$IDX" > "$TMP/rand_I.fa" 2>/dev/null
check_eq "extract -I FILE: same output" "$(cmp -s "$TMP/rand1.fa" "$TMP/rand_I.fa" && echo same)" "same"
"$RB" extract -I /nonexistent.isa "$IDX" B73_chr1:0-5 > /dev/null 2>&1
check_eq "extract -I missing file -> exit 1" "$?" "1"
mv "$TMP/moved.isa" "$IDX.isa"
check_eq "get is unchanged: sequence 0 is B73_chr1" "$("$RB" get "$IDX" 0 2>/dev/null | tail -1)" "$(fa_seq "$FA" B73_chr1)"
check_eq "get is unchanged: sequence 1 is its reverse complement" "$("$RB" get "$IDX" 1 2>/dev/null | tail -1)" "$(fa_seq "$FA" B73_chr1 | rev | tr ACGT TGCA)"

# --lift projection on the tiny fixture: every founder = B73 + a 500 bp insertion at 300
"$RB" extract --lift="$LIFT" --ref-prefix=B73 --to Oh43 --project-only "$IDX" B73_chr1:100-250 B73_chr1:250-400 B73_chr1:350-500 \
	> "$TMP/proj.tsv" 2>"$TMP/proj.log"
check_eq "extract --project-only exits 0" "$?" "0"
check_line_in "project-only: TSV header" "$TMP/proj.tsv" "$(printf 'qid\tref_contig\tref_start\tref_end\tassembly\tcontig\tstart\tend\tstrand\tstatus\tmethod\tsnap\tsnap_shift')"
check_line_in "project-only: collinear stretch before the insertion maps 1:1" "$TMP/proj.tsv" \
	"$(printf 'B73_chr1:100-250\tB73_chr1\t100\t250\tOh43\tOh43_chr1\t100\t250\t+\tok\tlocal,local\tsnap,snap\t0,0')"
check_line_in "project-only: an interval spanning the 500 bp insertion grows by 500" "$TMP/proj.tsv" \
	"$(printf 'B73_chr1:250-400\tB73_chr1\t250\t400\tOh43\tOh43_chr1\t250\t900\t+\tok\tlocal,local\tsnap,snap\t0,0')"
check_line_in "project-only: after the insertion, offset +500" "$TMP/proj.tsv" \
	"$(printf 'B73_chr1:350-500\tB73_chr1\t350\t500\tOh43\tOh43_chr1\t850\t1000\t+\tok\tlocal,local\tsnap,snap\t0,0')"
"$RB" extract --lift="$LIFT" --ref-prefix=B73 --to B73 --project-only "$IDX" B73_chr1:10-20 > "$TMP/proj_id.tsv" 2>/dev/null
check_line_in "project-only: the reference itself is the identity" "$TMP/proj_id.tsv" \
	"$(printf 'B73_chr1:10-20\tB73_chr1\t10\t20\tB73\tB73_chr1\t10\t20\t+\tok\tidentity,identity\tidentity,identity\t0,0')"
"$RB" extract --lift="$LIFT" --ref-prefix=B73 --to Oh43 "$IDX" B73_chr1:250-400 > "$TMP/proj.fa" 2>/dev/null
check_line_in "extract --lift: FASTA header carries the projected and the reference interval" "$TMP/proj.fa" ">B73_chr1:250-400 Oh43_chr1:250-900:+ ref=B73_chr1:250-400 snap=snap,snap"
check_eq "extract --lift: sequence = the founder at the projected interval" "$(fa_rec "$TMP/proj.fa" B73_chr1:250-400)" "$(fa_sub "$FA" Oh43_chr1 250 900)"
"$RB" extract --lift="$LIFT" --ref-prefix=B73 --to Oh43 --lift-local-win=0 --project-only "$IDX" B73_chr1:100-250 > "$TMP/proj_med.tsv" 2>/dev/null
check_eq "extract --lift-local-win=0: falls back to refmap's wide median" "$(awk -F'\t' 'NR==2{print $11}' "$TMP/proj_med.tsv")" "median,median"
"$RB" extract --lift="$LIFT" --ref-prefix=B73 --to Nope "$IDX" B73_chr1:0-5 > /dev/null 2>&1
check_eq "extract --to unknown assembly -> exit 1" "$?" "1"
"$RB" extract --project-only "$IDX" B73_chr1:0-5 > /dev/null 2>&1
check_eq "extract --project-only without --lift -> exit 1" "$?" "1"
"$RB" extract --lift="$LIFT" --to Oh43 "$IDX" B73_chr1:0-5 > /dev/null 2>&1
check_eq "extract --lift without --ref-prefix -> exit 1" "$?" "1"
"$RB" extract --lift="$LIFT" --ref-prefix=B73 --to Oh43 --project-only "$IDX" B73_chr1:0-5:- Oh43_chr1:0-5 > "$TMP/proj_bad.tsv" 2>"$TMP/proj_bad.log"
check_eq "extract --lift: strand suffix / non-reference contig -> exit 1" "$?" "1"
check_eq "extract --lift: input errors still get a TSV row" "$(awk -F'\t' 'NR>1 && $10 ~ /^input_error/' "$TMP/proj_bad.tsv" | wc -l | tr -d ' ')" "2"

# inverted and split founders: projection decides the strand, and refuses to guess
awk 'BEGIN{srand(3); s=""; for(i=0;i<3000;i++) s=s substr("ACGT",int(rand()*4)+1,1);
	print ">B73_chr1"; print s; print ">Fwd_chr1"; print s; print ">Spl_a"; print substr(s,1,1500); print ">Spl_b"; print substr(s,1501);
	t=""; for(i=0;i<3000;i++) t=t substr("ACGT",int(rand()*4)+1,1); print ">Far_chr1"; print t}' > "$TMP/inv.fa"
fa_seq "$TMP/inv.fa" B73_chr1 | rev | tr ACGT TGCA | awk '{print ">Inv_chr1"; print}' >> "$TMP/inv.fa"
"$RB" build -d -o "$TMP/inv.fmd" "$TMP/inv.fa" 2>/dev/null && "$RB" ssa -s4 -o "$TMP/inv.fmd.ssa" "$TMP/inv.fmd" 2>/dev/null
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)} END{if(n)print n"\t"l}' "$TMP/inv.fa" | gzip > "$TMP/inv.fmd.len.gz"
"$RB" lift --ref-prefix=B73 -k 31 -s 50 -o "$TMP/inv.lift" "$TMP/inv.fmd" "$TMP/inv.fa" 2>/dev/null
"$RB" isa -b 4 "$TMP/inv.fmd" 2>/dev/null
printf 'qid\tref_contig\tref_start\tref_end\tassembly\nfwd\tB73_chr1\t1000\t1500\tFwd\ninv\tB73_chr1\t1000\t1500\tInv\nsplit\tB73_chr1\t1000\t2000\tSpl\nfar\tB73_chr1\t1000\t1500\tFar\n' > "$TMP/inv_q.tsv"
"$RB" extract --lift="$TMP/inv.lift" --lift-k 31 --ref-prefix=B73 --project-only -f "$TMP/inv_q.tsv" "$TMP/inv.fmd" > "$TMP/inv_proj.tsv" 2>"$TMP/inv_proj.log"
check_eq "extract --project-only: projection failures do not set the exit status" "$?" "0"
check_line_in "projection: an identical founder maps 1:1" "$TMP/inv_proj.tsv" "$(printf 'fwd\tB73_chr1\t1000\t1500\tFwd\tFwd_chr1\t1000\t1500\t+\tok\tlocal,local\tsnap,snap\t0,0')"
check_line_in "projection: a reverse-complemented founder maps to [L-end, L-start) on '-'" "$TMP/inv_proj.tsv" "$(printf 'inv\tB73_chr1\t1000\t1500\tInv\tInv_chr1\t1500\t2000\t-\tok\tlocal,local\tsnap,snap\t0,0')"
check_eq "projection: ends on different contigs -> explicit failure, no guess" "$(awk -F'\t' '$1=="split"{print $6"|"$10}' "$TMP/inv_proj.tsv")" ".|projection_failed: start and end project to different contigs"
check_eq "projection: no anchors -> explicit failure" "$(awk -F'\t' '$1=="far"{print $6"|"$10}' "$TMP/inv_proj.tsv")" ".|projection_failed: no confident projection for the start"
check_eq "projection failures are also reported on stderr" "$(grep -c 'projection failed' "$TMP/inv_proj.log")" "2"
"$RB" extract --lift="$TMP/inv.lift" --lift-k 31 --ref-prefix=B73 -f "$TMP/inv_q.tsv" "$TMP/inv.fmd" > "$TMP/inv_proj.fa" 2>/dev/null
check_eq "extract --lift: only the projected queries get a FASTA record" "$(grep -c '^>' "$TMP/inv_proj.fa")" "2"
check_eq "extract --lift: an inverted founder's '-' extraction reads in reference orientation" \
	"$(fa_rec "$TMP/inv_proj.fa" inv)" "$(fa_sub "$TMP/inv.fa" B73_chr1 1000 1500)"
# snapping repairs a wrong --lift-k on the inverted founder; --no-snap shows the raw projection
"$RB" extract --lift="$TMP/inv.lift" --lift-k 1 --ref-prefix=B73 --project-only -f "$TMP/inv_q.tsv" "$TMP/inv.fmd" > "$TMP/inv_k1.tsv" 2>/dev/null
check_eq "snap: an inverted locus projected with the wrong k is snapped back to the homologous bases" \
	"$(awk -F'\t' '$1=="inv"{print $6":"$7"-"$8":"$9" "$12" "$13}' "$TMP/inv_k1.tsv")" "Inv_chr1:1500-2000:- snap,snap 30,30"
"$RB" extract --lift="$TMP/inv.lift" --lift-k 1 --no-snap --ref-prefix=B73 --project-only -f "$TMP/inv_q.tsv" "$TMP/inv.fmd" > "$TMP/inv_k1_nosnap.tsv" 2>/dev/null
check_eq "--no-snap: the raw projection (off by k-1 = 30 here), labelled proj" \
	"$(awk -F'\t' '$1=="inv"{print $6":"$7"-"$8":"$9" "$12" "$13}' "$TMP/inv_k1_nosnap.tsv")" "Inv_chr1:1470-1970:- proj,proj 0,0"
check_eq "--no-snap: still 13 columns" "$(awk -F'\t' 'NR>1{print NF}' "$TMP/inv_k1_nosnap.tsv" | sort -u)" "13"
"$RB" extract --lift="$TMP/inv.lift" --snap-flank=8 --ref-prefix=B73 --project-only -f "$TMP/inv_q.tsv" "$TMP/inv.fmd" > /dev/null 2>&1
check_eq "--snap-flank below the minimum aligned length -> exit 1" "$?" "1"
# projected-span sanity check (--max-span-excess / --max-span-ratio, after snapping):
# Ins = B73[0,50k) + 60 kb insertion + B73[50k,100k); Clp = the first 112 kb of Ins
awk 'BEGIN{srand(5); s=""; for(i=0;i<100000;i++) s=s substr("ACGT",int(rand()*4)+1,1); x=""; for(i=0;i<60000;i++) x=x substr("ACGT",int(rand()*4)+1,1);
	ins=substr(s,1,50000) x substr(s,50001); print ">B73_chr1"; print s; print ">Ins_chr1"; print ins; print ">Clp_chr1"; print substr(ins,1,112000)}' > "$TMP/span.fa"
"$RB" build -d -o "$TMP/span.fmd" "$TMP/span.fa" 2>/dev/null && "$RB" ssa -s4 -o "$TMP/span.fmd.ssa" "$TMP/span.fmd" 2>/dev/null
awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)} END{if(n)print n"\t"l}' "$TMP/span.fa" | gzip > "$TMP/span.fmd.len.gz"
"$RB" lift --ref-prefix=B73 -k 31 -s 50 -o "$TMP/span.lift" "$TMP/span.fmd" "$TMP/span.fa" 2>/dev/null
"$RB" isa -b 6 "$TMP/span.fmd" 2>/dev/null
span() { # span <extra options...> -> "status-or-interval" for the two Ins queries and the Clp one
	"$RB" extract --lift="$TMP/span.lift" --lift-k 31 --ref-prefix=B73 --project-only "$@" \
		-f "$TMP/span_q.tsv" "$TMP/span.fmd" 2>/dev/null | awk -F'\t' 'NR>1{printf "%s=%s ", $1, ($10 ~ /^ok/? $6":"$7"-"$8":"$10 : $10)}'
}
printf 'trip\tB73_chr1\t49000\t51000\tIns\nplausible\tB73_chr1\t30000\t70000\tIns\nclamped\tB73_chr1\t49000\t60000\tClp\n' > "$TMP/span_q.tsv"
check_eq "span check: 2 kb -> 62 kb fails, 40 kb -> 100 kb (big but plausible insertion) passes, a clamped 63 kb span fails" "$(span)" \
	"trip=projection_failed: projected span implausibly long (62000 bp for 2000 bp) plausible=Ins_chr1:30000-130000:ok clamped=projection_failed: projected span implausibly long (63000 bp for 11000 bp) "
check_eq "span check: --max-span-excess=0 turns it off" "$(span --max-span-excess=0)" \
	"trip=Ins_chr1:49000-111000:ok plausible=Ins_chr1:30000-130000:ok clamped=Clp_chr1:49000-112000:ok_clamped "
check_eq "span check: --max-span-ratio=0 turns it off" "$(span --max-span-ratio=0)" \
	"trip=Ins_chr1:49000-111000:ok plausible=Ins_chr1:30000-130000:ok clamped=Clp_chr1:49000-112000:ok_clamped "
check_eq "span check: BOTH conditions needed (ratios 31/5.7 < 40; excesses 60/52 kb < 70 kb)" "$(span --max-span-ratio=40) $(span --max-span-excess=70000)" \
	"trip=Ins_chr1:49000-111000:ok plausible=Ins_chr1:30000-130000:ok clamped=Clp_chr1:49000-112000:ok_clamped  trip=Ins_chr1:49000-111000:ok plausible=Ins_chr1:30000-130000:ok clamped=Clp_chr1:49000-112000:ok_clamped "
check_eq "span check: also applies with --no-snap" "$(span --no-snap | cut -d' ' -f1-2)" \
	"trip=projection_failed: projected"
"$RB" extract --lift="$TMP/span.lift" --lift-k 31 --ref-prefix=B73 -f "$TMP/span_q.tsv" "$TMP/span.fmd" > "$TMP/span.fa.out" 2>"$TMP/span.log"
check_eq "span check: exit status unchanged (0), only the plausible query gets a record" "$?:$(grep '^>' "$TMP/span.fa.out" | cut -d' ' -f1 | tr '\n' ' ')" "0:>plausible "
check_eq "span check: failures reported on stderr" "$(grep -c 'projection failed: projected span implausibly long' "$TMP/span.log")" "2"
"$RB" extract --lift="$TMP/span.lift" --ref-prefix=B73 --max-span-ratio=-1 -f "$TMP/span_q.tsv" "$TMP/span.fmd" > /dev/null 2>&1
check_eq "span check: a negative value is an error" "$?" "1"

"$RB" extract -I "$TMP/inv.fmd.isa" "$IDX" B73_chr1:0-5 > /dev/null 2>"$TMP/ext_wrongisa.log"
check_eq "extract: an .isa from another index is rejected" "$(grep -c 'not built from this index' "$TMP/ext_wrongisa.log")" "1"

# the shared simulated fixture (simulate_pangenome.py), scaled down: SNPs, small and
# large indels, a multi-copy element; needs numpy for the simulator and the checker
if [ -n "$NPY_PY" ]; then
	SIM="$TMP/sim"
	"$NPY_PY" simulate_pangenome.py --out "$SIM" --seed 7 --chrom-lens 40000 25000 --founders Oh43 Mo17 Ki3 \
		--n-direct 150 --n-ref 150 --window 2000 --te-copies 3 --large-rate 5e-5 > "$TMP/sim.log" 2>&1
	check_eq "simulated fixture generated" "$?" "0"
	"$RB" build -d -o "$SIM/p.fmd" "$SIM/pangenome.fa" 2>/dev/null && "$RB" ssa -s4 -o "$SIM/p.fmd.ssa" "$SIM/p.fmd" 2>/dev/null
	awk '/^>/{if(n)print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)} END{if(n)print n"\t"l}' "$SIM/pangenome.fa" | gzip > "$SIM/p.fmd.len.gz"
	"$RB" isa -b 6 "$SIM/p.fmd" 2>/dev/null
	"$RB" lift --ref-prefix=B73 -s 200 -o "$SIM/p.lift" "$SIM/p.fmd" "$SIM/B73.fa" 2>/dev/null
	"$RB" extract -f "$SIM/direct_queries.tsv" "$SIM/p.fmd" > "$SIM/direct.fa" 2>/dev/null
	check_eq "sim: extract exits 0 on the direct queries (incl. edge cases)" "$?" "0"
	"$RB" extract --lift="$SIM/p.lift" --ref-prefix=B73 --project-only -f "$SIM/ref_queries.tsv" "$SIM/p.fmd" > "$SIM/proj.tsv" 2>/dev/null
	"$RB" extract --lift="$SIM/p.lift" --ref-prefix=B73 --no-snap --project-only -f "$SIM/ref_queries.tsv" "$SIM/p.fmd" > "$SIM/proj_nosnap.tsv" 2>/dev/null
	"$RB" extract --lift="$SIM/p.lift" --ref-prefix=B73 --no-snap -f "$SIM/ref_queries.tsv" "$SIM/p.fmd" > "$SIM/proj_nosnap.fa" 2>/dev/null
	"$RB" extract --lift="$SIM/p.lift" --ref-prefix=B73 -f "$SIM/ref_queries.tsv" "$SIM/p.fmd" > "$SIM/proj.fa" 2>/dev/null
	# hard: every direct query exact; every projected record exact at its projected interval;
	# without snapping (lift anchor density limits it): >= 85% within 100 bp at both ends;
	# with snapping (default): >= 90% of intervals base-exact and >= 97% within 5 bp
	if "$NPY_PY" check_extract_fixture.py "$SIM" "$SIM/direct.fa" "$SIM/proj_nosnap.tsv" "$SIM/proj_nosnap.fa" --min-proj-within 0.85 100 > "$TMP/sim_check_nosnap.log" 2>&1; then
		pass
	else
		fail "simulated fixture check (--no-snap) failed; see $TMP/sim_check_nosnap.log"; cat "$TMP/sim_check_nosnap.log" >&2
	fi
	if "$NPY_PY" check_extract_fixture.py "$SIM" "$SIM/direct.fa" "$SIM/proj.tsv" "$SIM/proj.fa" --min-proj-exact 0.90 --min-proj-within 0.97 5 > "$TMP/sim_check.log" 2>&1; then
		pass
	else
		fail "simulated fixture check (snap) failed; see $TMP/sim_check.log"; cat "$TMP/sim_check.log" >&2
	fi
else
	echo "SKIP: simulated-pangenome extract check needs a python3 with numpy" >&2
fi

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

	# extract: multi-threaded direct extraction, and --lift projection incl. failures
	valgrind --error-exitcode=99 --leak-check=full -q \
		"$RB" extract -t2 -f "$TMP/rand.tsv" "$IDX" > /dev/null 2>"$TMP/valgrind_extract.log"
	if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks on extract; see $TMP/valgrind_extract.log"; fi
	valgrind --error-exitcode=99 --leak-check=full -q \
		"$RB" extract --lift="$TMP/inv.lift" --lift-k 31 --ref-prefix=B73 -f "$TMP/inv_q.tsv" "$TMP/inv.fmd" > /dev/null 2>"$TMP/valgrind_extract_lift.log"
	if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks on extract --lift; see $TMP/valgrind_extract_lift.log"; fi
	if [ -s "$TMP/sim/p.lift" ]; then # snapping incl. the wide retry and walk-inward paths
		valgrind --error-exitcode=99 --leak-check=full -q \
			"$RB" extract --lift="$TMP/sim/p.lift" --ref-prefix=B73 -t2 -f "$TMP/sim/ref_queries.tsv" "$TMP/sim/p.fmd" > /dev/null 2>"$TMP/valgrind_extract_snap.log"
		if [ $? -eq 0 ]; then pass; else fail "valgrind reported errors or leaks on extract snapping; see $TMP/valgrind_extract_snap.log"; fi
	fi
else
	echo "SKIP: valgrind not installed; skipping the memory-safety pass" >&2
fi

echo "run_integration: $n_pass passed, $n_fail failed"
[ "$n_fail" -eq 0 ]
