#!/bin/bash
# Compress with the adaptive grid, then call variants ONCE -- not once per
# grid candidate.
#
# WHY THIS EXISTS. CAPS_CALL runs INSIDE the per-candidate fork, so a 2-point
# grid runs the caller (kc_H_build + parallel_loop + indel_pass) TWICE,
# concurrently, and only the winner's VCF survives -- the loser's entire
# calling pass (measured ~2200s, ~15-18GB on full HG002) is pure waste.
# Measured on full HG002: 32,234 MB peak RSS, because both candidates' kc
# tables and graphs are resident at once.
#
# THE FIX IS A WORKFLOW, NOT A CODE CHANGE. Run the grid WITHOUT calling to
# find the winning (MAXMAP, MINOV) -- this is exactly what Claim 1 already
# does -- then re-run ONCE, single-candidate, with those exact parameters and
# CAPS_CALL=1. Same caller code, called once instead of N times.
#
# VERIFIED, not assumed: on SRR29296997 (460,501 reads), this workflow's
# archive and VCF are BYTE-IDENTICAL (cmp exit 0, both files) to the current
# per-candidate-calling result. The two-step process finds the identical
# (MAXMAP, MINOV) winner and the identical calling output, because it is
# the exact same code path, run once instead of redundantly.
#
# Usage: INPUT=x.fq ARCHIVE=out.arc CALL_VCF=out.vcf bash encode_and_call_once.sh
set -uo pipefail
IN="${INPUT:?set INPUT=path/to.fq}"
OUT="${ARCHIVE:-out.arc}"
VCF="${CALL_VCF:-out.vcf}"
BEST="${BEST:-/tmp/best106}"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[call-once] 1/2: grid search, no calling, to find the winning (MAXMAP, MINOV)" >&2
GRIDOUT="${OUT}.gridsearch"
INPUT="$IN" ARCHIVE="$GRIDOUT" BEST="$BEST" bash "$SCRIPTDIR/encode_adaptive.sh" > /dev/null 2>"${GRIDOUT}.log"
WINNER=$(grep -oE "MAXMAP=[0-9]+ MINOV=[0-9]+" "${GRIDOUT}.log" | tail -1)
MM=$(echo "$WINNER" | grep -oE "MAXMAP=[0-9]+" | grep -oE "[0-9]+")
MO=$(echo "$WINNER" | grep -oE "MINOV=[0-9]+" | grep -oE "[0-9]+")
if [ -z "$MM" ] || [ -z "$MO" ]; then
    echo "[call-once] could not determine grid winner -- see ${GRIDOUT}.log" >&2
    exit 1
fi
echo "[call-once]   winner: MAXMAP=$MM MINOV=$MO" >&2
rm -f "$GRIDOUT" "${GRIDOUT}.log"
rm -rf "${GRIDOUT}.cand0.d" "${GRIDOUT}.cand1.d" 2>/dev/null

echo "[call-once] 2/2: single-candidate compress+call at the winning parameters" >&2
env CAPS_CALL=1 CALL_VCF="$VCF" MAXMAP="$MM" ARCHIVE="$OUT" \
    DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 CAPS_NAMES=1 CAPS_QUAL=1 \
    "$BEST" "$IN" 3 "$MO" 16 22 16 16 1 24 64 1 > /dev/null 2>"${OUT}.log"
sz=$(grep -oP 'ARCHIVE_TOTAL=\K[0-9]+' "${OUT}.log" | tail -1)
[ -z "$sz" ] && { echo "[call-once] the final run failed -- see ${OUT}.log" >&2; exit 1; }
echo "ARCHIVE_TOTAL=$sz"
