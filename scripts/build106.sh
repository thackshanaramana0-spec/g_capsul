#!/bin/bash
# Build the in-process encoder (stage 106).
#
# Recorded as a script because it had only ever existed as a shell-history
# one-liner -- which is exactly how -fopenmp went missing. Without that flag the
# single `#pragma omp parallel for` in the sweep is silently discarded and the
# largest stage runs serial; linking it correctly was worth -45% wall time.
set -e
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-/tmp/best106}"
g++ -O3 -march=native -std=c++17 -pthread -fopenmp -o "$OUT" \
    "$HERE/stages/106_inprocess.cpp" \
    "$HERE/thirdparty/ppmd/Ppmd7.c" "$HERE/thirdparty/ppmd/Ppmd7Enc.c" \
    "$HERE/thirdparty/ppmd/Ppmd7Dec.c" "$HERE/thirdparty/ppmd/Alloc.c" \
    "$HERE/thirdparty/ppmd/CpuArch.c" \
    "$HERE/thirdparty/fse/fse_compress.c" "$HERE/thirdparty/fse/fse_decompress.c" \
    "$HERE/thirdparty/fse/huf_compress.c" "$HERE/thirdparty/fse/huf_decompress.c" \
    "$HERE/thirdparty/fse/entropy_common.c" "$HERE/thirdparty/fse/hist.c" \
    "$HERE/thirdparty/fse/debug.c" \
    -I"$HERE" -I"$HERE/include" -I"$HERE/thirdparty/ppmd" -I"$HERE/thirdparty/fse" -llzma
echo "built $OUT"
