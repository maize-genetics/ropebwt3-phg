#!/bin/sh
# Build an OpenMP-enabled `ropebwt3-omp` on macOS, where Apple clang has no
# bundled OpenMP and (on recent macOS) Homebrew gcc / brew may be unusable.
#
# Strategy: compile with Apple clang (which targets the current macOS SDK
# natively) against an LLVM OpenMP runtime (omp.h + libomp.dylib). We locate
# that runtime from an active conda/miniconda install (conda-forge ships
# `llvm-openmp`); override with OMP_PREFIX=/path if yours lives elsewhere.
#
# This parallelizes the slow libsais BWT-construction step (`build -p N`);
# the plain `make omp=0` binary leaves it single-threaded. Linux users with
# real gcc should just use `make` instead of this script.
#
# Usage:  sh build-omp.sh           # autodetect conda libomp
#         OMP_PREFIX=/x sh build-omp.sh
set -e

PREFIX="${OMP_PREFIX:-${CONDA_PREFIX:-$HOME/Code/miniconda3}}"
if [ ! -f "$PREFIX/include/omp.h" ] || [ ! -f "$PREFIX/lib/libomp.dylib" ]; then
  echo "error: omp.h / libomp.dylib not found under $PREFIX" >&2
  echo "       set OMP_PREFIX=/path/to/llvm-openmp (e.g. a conda env), or" >&2
  echo "       conda install -c conda-forge llvm-openmp" >&2
  exit 1
fi
echo "using OpenMP runtime: $PREFIX"

SELF=$(cd "$(dirname "$0")" && pwd)
B="${TMPDIR:-/tmp}/ropebwt3-omp-build"
rm -rf "$B" && mkdir -p "$B"
cp "$SELF"/*.c "$SELF"/*.h "$SELF/Makefile" "$B"/

# omp=0 stops the Makefile adding bare `-fopenmp` (Apple clang rejects it); we
# inject the clang-compatible OpenMP flags via CPPFLAGS/LIBS instead. An rpath
# to $PREFIX/lib lets the binary find libomp.dylib at runtime.
( cd "$B" && make omp=0 CC=clang -j8 \
    CPPFLAGS="-DLIBSAIS_OPENMP -Xpreprocessor -fopenmp -I$PREFIX/include" \
    LIBS="-lpthread -lz -lm -L$PREFIX/lib -lomp -Wl,-rpath,$PREFIX/lib" )

cp "$B/ropebwt3" "$SELF/ropebwt3-omp"
echo "built: $SELF/ropebwt3-omp"
"$SELF/ropebwt3-omp" version 2>&1 | head -1
echo "run with e.g.:  ./ropebwt3-omp build -d -t 12 -p 12 -o idx.fmd in1.fa.gz in2.fa.gz ..."
