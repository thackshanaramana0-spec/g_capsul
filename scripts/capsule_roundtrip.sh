#!/bin/bash
# compress to ONE file, then decompress using ONLY that file. Single command
# each way -- the read reconstruction is in the decoder now, not a Python step.
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
"$DEC" reads.capsule out out/reads.seq >dec.log 2>&1
awk 'NR%4==2' "$IN" > orig.seq
if cmp -s orig.seq out/reads.seq; then
    echo "LOSSLESS  archive=$SZ B  reads=$(wc -l < orig.seq)  input=$IN"
else
    echo "FAILED    archive=$SZ B  input=$IN"; exit 1
fi
