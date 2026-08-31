#!/bin/bash
# Build the in-process encoder. Recorded as a script because it had only ever
# existed as a shell-history one-liner, which is not reproducible.
set -e
OUT="${1:-/tmp/best106}"
g++ -O3 -march=native -std=c++17 -pthread -fopenmp -o "$OUT" \
    106_inprocess.cpp \
    thirdparty/ppmd/Ppmd7.c thirdparty/ppmd/Ppmd7Enc.c thirdparty/ppmd/Ppmd7Dec.c \
    thirdparty/ppmd/Alloc.c thirdparty/ppmd/CpuArch.c \
    thirdparty/fse/fse_compress.c thirdparty/fse/fse_decompress.c \
    thirdparty/fse/huf_compress.c thirdparty/fse/huf_decompress.c \
    thirdparty/fse/entropy_common.c thirdparty/fse/hist.c thirdparty/fse/debug.c \
    -Ithirdparty/ppmd -Ithirdparty/fse -llzma
echo "built $OUT"
