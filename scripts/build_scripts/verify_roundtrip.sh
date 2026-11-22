#!/usr/bin/env bash
# Ordered, byte-for-byte FASTQ round trip: compress, decompress, rebuild the
# four-line records, and compare the rebuilt file with the input using cmp.
# Unlike a comparison of sorted (name, sequence, quality) tuples, this checks read
# order, line 3 and every byte of the file.
#
#   usage: verify_roundtrip.sh <reads.fq> [workdir]
#   env:   BEST=<encoder binary>  DEC=<capsule_decode binary>   (default: build/ then bin/)
#
# Line 3 is rebuilt from the archive's layout: '+' alone, or '+' followed by the
# header. The mode is read from the first input record, and the encoder refuses
# inputs where the two are mixed, so one mode describes the whole file.
# Exit status: 0 byte-identical, 1 differs, 2 the encoder refused the input.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IN="${1:?usage: verify_roundtrip.sh reads.fq [workdir]}"
W="${2:-$(mktemp -d)}"; mkdir -p "$W"
BEST="${BEST:-$(ls "$HERE/build/best106" "$HERE/bin/best106" 2>/dev/null | head -1)}"
DEC="${DEC:-$(ls "$HERE/build/capsule_decode" "$HERE/bin/capsule_decode" 2>/dev/null | head -1)}"
[ -x "${BEST:-}" ] && [ -x "${DEC:-}" ] || { echo "build the encoder and decoder first (see README)"; exit 2; }

CAPS_NAMES=1 CAPS_QUAL=1 INPUT="$IN" ARCHIVE="$W/rt.capsule" BEST="$BEST" \
    bash "$HERE/scripts/build_scripts/encode_adaptive.sh" >"$W/enc.log" 2>&1
if [ ! -s "$W/rt.capsule" ]; then
    echo "REFUSED   the encoder did not write an archive"
    grep -a -h -E "FATAL|OUT OF SCOPE" "$W"/rt.capsule.log* "$W/enc.log" 2>/dev/null | head -3
    exit 2
fi
mkdir -p "$W/dec"
"$DEC" "$W/rt.capsule" "$W/dec" "$W/dec/reads.fq" >"$W/dec.log" 2>&1
l3=$(sed -n 3p "$IN")
if [ "$l3" = "+" ]; then
    paste -d'\n' <(sed 's/^/@/' "$W/dec/reads.fq.names") "$W/dec/reads.fq" \
                 <(sed 's/.*/+/' "$W/dec/reads.fq") "$W/dec/reads.fq.qual" > "$W/rebuilt.fq"
else
    paste -d'\n' <(sed 's/^/@/' "$W/dec/reads.fq.names") "$W/dec/reads.fq" \
                 <(sed 's/^/+/' "$W/dec/reads.fq.names") "$W/dec/reads.fq.qual" > "$W/rebuilt.fq"
fi
if cmp -s "$IN" "$W/rebuilt.fq"; then
    echo "BYTE-IDENTICAL  archive=$(stat -c%s "$W/rt.capsule") bytes  records=$(( $(wc -l < "$IN") / 4 ))  input=$IN"
else
    echo "DIFFERS  first difference: $(cmp "$IN" "$W/rebuilt.fq" 2>&1 | head -1)"; exit 1
fi
