#!/bin/bash
# THE test: compress to ONE file, then decompress using ONLY that file.
#
# The old check decoded the encoder's dumped intermediate streams, so a stream
# missing from the ARCHIVE could not be detected -- which is exactly how the
# reference destinations, lengths and RC flags went unnoticed. This reads the
# archive and nothing else.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN="${1:?usage: capsule_roundtrip.sh reads.fq [workdir]}"
W="${2:-$(mktemp -d)}"; mkdir -p "$W/out"
ENC="${ENC:-/tmp/best106}"; DEC="${DEC:-/tmp/capsule_d}"

LMAX=$(awk 'NR%4==2{if(length($0)>m)m=length($0)} NR>400000{exit} END{print m}' "$IN")
MM="${MAXMAP:-$((LMAX/13))}"; [ "$MM" -lt 6 ] && MM=6
MO="${MINOV:-16}"

cd "$W"
DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 MAXMAP=$MM MINOV=$MO ARCHIVE="$W/reads.capsule" \
    "$ENC" "$IN" 3 "$MO" 16 22 16 16 1 24 64 1 >/dev/null 2>enc.log
SZ=$(stat -c%s reads.capsule)

# decompress from the archive alone
"$DEC" reads.capsule out >dec.log 2>&1
read -r PG MAIN < out/pg_params.txt
python3 "$HERE/scripts/decode_105.py" out "$PG" "$MAIN" out/reads.seq >>dec.log 2>&1

awk 'NR%4==2' "$IN" > orig.seq
if cmp -s orig.seq out/reads.seq; then
    echo "LOSSLESS  archive=$SZ B  reads=$(wc -l < orig.seq)  input=$IN"
else
    echo "FAILED    archive=$SZ B  input=$IN"; exit 1
fi
