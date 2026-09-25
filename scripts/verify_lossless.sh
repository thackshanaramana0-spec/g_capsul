#!/bin/bash
# Encode a FASTQ, then decode the streams back and compare against the ORIGINAL
# file's sequence column. This is a real round trip, not a coder-level check:
# the reference is the input file itself.
#
#   scripts/verify_lossless.sh reads.fq [workdir]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN="${1:?usage: verify_lossless.sh reads.fq [workdir]}"
W="${2:-$(mktemp -d)}"; mkdir -p "$W"
# Always build from THIS tree unless the caller pins BEST explicitly. Reusing
# whatever /tmp/best106 happens to be silently verified a stale binary once.
BEST="${BEST:-$W/best106}"
[ -x "$BEST" ] || bash "$HERE/scripts/build106.sh" "$BEST" >/dev/null

LMAX=$(awk 'NR%4==2{if(length($0)>m)m=length($0)} NR>400000{exit} END{print m}' "$IN")
C1=$(( LMAX / 13 )); [ "$C1" -lt 6 ] && C1=6
C2=$(( LMAX / 5  )); [ "$C2" -lt 6 ] && C2=6
M1=16; M2=$(( LMAX * 35 / 100 )); [ "$M2" -lt 16 ] && M2=16

cd "$W"
CANDIDATES="${C1}:${M1},${C1}:${M2},${C2}:${M1},${C2}:${M2}" ARCHIVE="$W/out.arc" \
  DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 VERIFY_DUMP=1 \
  "$BEST" "$IN" 3 16 16 22 16 16 1 24 64 1 >/dev/null 2>"$W/enc.log"
sed -n 's/.*\[a3\] chose //p' "$W/enc.log" || true

read -r PG MAIN < "$W/pg_params.txt"
python3 "$HERE/scripts/decode_105.py" "$W" "$PG" "$MAIN" "$W/dec.seq" >"$W/dec.log" 2>&1

awk 'NR%4==2' "$IN" > "$W/orig.seq"
if cmp -s "$W/orig.seq" "$W/dec.seq"; then
    echo "LOSSLESS  archive=$(stat -c%s "$W/out.arc") bytes  input=$IN"
else
    echo "LOSSY     orig=$(wc -l <"$W/orig.seq") lines  decoded=$(wc -l <"$W/dec.seq") lines"; exit 1
fi
