#!/bin/bash
# Adaptive mapping-ceiling selection.
#
# VERIFIED FINDING (PLAN_PG_VOLUME.md 1.4): the optimal mapping ceiling varies
# 5.6x across datasets (halo ~8 = L/19, sulfo ~45 = L/5.6) while read length
# varies only 1.66x, so NO constant and no fixed ratio of read length can
# express it -- including the MAXMAP=Lmax/13 currently compiled in. The driver
# is second-region compressibility (MEM-second removal 93.1% on halo vs 50.3%
# on sulfo), which is not knowable until after the decision it governs.
#
# That circularity is why this is resolved by MEASUREMENT rather than a
# formula: encode at each candidate ceiling and keep the smallest archive.
# Non-regressing by construction -- the result is never worse than the better
# candidate, on any input.
#
# COST: one extra encode pass. The expensive prefix stages (load, round 1,
# round 2, chain emission) are repeated because MAXMAP is consumed inside the
# mapping loop and everything after depends on it. Sharing that prefix is a
# pure-runtime refactor that cannot change the output, and is the obvious
# follow-up.
#
# Usage: INPUT=x.fq ARCHIVE=out.arc bash encode_adaptive.sh [extra args...]
set -uo pipefail
IN="${INPUT:?set INPUT=path/to.fq}"
OUT="${ARCHIVE:-out.arc}"
BEST="${BEST:-/tmp/best106}"
ARGS="${ENC_ARGS:-3 16 16 22 16 16 1 24 64 1}"

# Candidates: the compiled-in ratio, and a permissive one. Two points bracket
# the observed optima (L/19 .. L/5.6) closely enough that min() lands on the
# better side for every dataset measured; add more only with evidence.
LMAX=$(awk 'NR%4==2{if(length($0)>m)m=length($0)} NR>400000{exit} END{print m}' "$IN")
C1=$(( LMAX / 13 )); [ "$C1" -lt 6 ] && C1=6
C2=$(( LMAX / 5  )); [ "$C2" -lt 6 ] && C2=6

# MINOV candidates. The optimal overlap floor is ALSO dataset-dependent and
# ALSO not a fixed fraction -- measured optima 0.10L (P. aeruginosa), 0.30L
# (H. salinarum), 0.45L (S. acidocaldarius) -- and it INTERACTS with the
# ceiling, so neither can be chosen independently. At MAXMAP=50 on sulfo,
# raising MINOV 16 -> 113 cuts mem_triples 35% (141,115 -> 91,690) because a
# tighter chain floor yields a smaller pg, which MEM then has less work to
# collapse. An earlier test found MINOV worth only 1.6% -- that was measured
# at the wrong ceiling (MAXMAP=19) and understated it.
M1=16
M2=$(( LMAX * 35 / 100 )); [ "$M2" -lt 16 ] && M2=16

best_sz=""; best_c=""; best_m=""
for C in "$C1" "$C2"; do
  for M in "$M1" "$M2"; do
    tag="${C}_${M}"
    MAXMAP="$C" ARCHIVE="${OUT}.cand${tag}" \
        DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 "$BEST" "$IN" 3 "$M" 16 22 16 16 1 24 64 1 \
        > /dev/null 2>"${OUT}.log${tag}"
    sz=$(grep -oP 'ARCHIVE_TOTAL=\K[0-9]+' "${OUT}.log${tag}" 2>/dev/null)
    [ -z "$sz" ] && continue
    echo "  [adaptive] MAXMAP=$C MINOV=$M -> $sz B" >&2
    if [ -z "$best_sz" ] || [ "$sz" -lt "$best_sz" ]; then
        best_sz="$sz"; best_c="$C"; best_m="$M"; fi
  done
done

[ -z "$best_sz" ] && { echo "[adaptive] all candidates failed" >&2; exit 1; }
mv -f "${OUT}.cand${best_c}_${best_m}" "$OUT"
rm -f "${OUT}.cand"*
cp -f "${OUT}.log${best_c}_${best_m}" "${OUT}.log"; rm -f "${OUT}.log"[0-9]*
echo "  [adaptive] chose MAXMAP=$best_c MINOV=$best_m -> $best_sz B" >&2
echo "ARCHIVE_TOTAL=$best_sz"
