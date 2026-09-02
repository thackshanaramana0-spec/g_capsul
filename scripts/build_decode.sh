#!/bin/bash
# Build the archive decoder. Mirrors build106.sh, including the gcc/g++ split
# for the vendored htscodecs C sources.
set -e
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-/tmp/capsule_decode}"
HTSOBJ="$(mktemp -d)"
gcc -O3 -I"$HERE/thirdparty/htscodecs" -c "$HERE/thirdparty/htscodecs/fqzcomp_qual.c" -o "$HTSOBJ/fqzcomp_qual.o"
gcc -O3 -I"$HERE/thirdparty/htscodecs" -c "$HERE/thirdparty/htscodecs/utils.c"        -o "$HTSOBJ/utils.o"
g++ -O3 -march=native -std=c++17 -pthread -o "$OUT" \
    "$HERE/stages/capsule_decode.cpp" \
    "$HERE/thirdparty/ppmd/Ppmd7.c" "$HERE/thirdparty/ppmd/Ppmd7Enc.c" \
    "$HERE/thirdparty/ppmd/Ppmd7Dec.c" "$HERE/thirdparty/ppmd/Alloc.c" \
    "$HERE/thirdparty/ppmd/CpuArch.c" \
    "$HERE/thirdparty/fse/fse_compress.c" "$HERE/thirdparty/fse/fse_decompress.c" \
    "$HERE/thirdparty/fse/huf_compress.c" "$HERE/thirdparty/fse/huf_decompress.c" \
    "$HERE/thirdparty/fse/entropy_common.c" "$HERE/thirdparty/fse/hist.c" \
    "$HERE/thirdparty/fse/debug.c" \
    "$HTSOBJ/fqzcomp_qual.o" "$HTSOBJ/utils.o" \
    -I"$HERE" -I"$HERE/include" -I"$HERE/thirdparty/ppmd" -I"$HERE/thirdparty/fse" \
    -I"$HERE/thirdparty/htscodecs" -llzma -lm
rm -rf "$HTSOBJ"
echo "built $OUT"
