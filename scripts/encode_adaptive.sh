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

# SCOPE CHECK -- runs before the (expensive, multi-candidate) encode below, so
# out-of-scope input is refused in seconds with a full explanation instead of
# after however long the sweep takes. This is a friendly pre-flight report,
# not the safety mechanism: the encoder itself (stages/106_inprocess.cpp)
# refuses unconditionally and cannot be bypassed by skipping this script and
# calling the binary directly. Set CAPS_SKIP_SCOPE_CHECK=1 to bypass this
# report only (e.g. re-running a file already known to be in scope) -- the
# encoder's own gate still applies regardless.
HERE_SCA="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -z "${CAPS_SKIP_SCOPE_CHECK:-}" ] && [ -x "$(command -v python3)" ] \
   && [ -f "$HERE_SCA/industry/check_input_scope.py" ]; then
    python3 "$HERE_SCA/industry/check_input_scope.py" "$IN"
    SCOPE_RC=$?
    if [ "$SCOPE_RC" -ge 2 ]; then
        echo "encode_adaptive.sh: refusing -- see the OUT OF SCOPE finding(s) above." >&2
        exit "$SCOPE_RC"
    fi
    # SCOPE_RC 1 (in scope, with caveats) prints and continues; 0 is silent.
fi

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

# A3: all four candidates now run inside ONE process off a shared prefix.
# MINOV first takes effect at the round-2 sweep and MAXMAP not until mapping,
# so load + seed index + round 1 -- measured at 24.8% of runtime -- is identical
# for every candidate and was previously recomputed four times. The binary forks
# after round 1, so each candidate gets an exact copy-on-write snapshot and the
# result is byte-identical to running them separately (verified on E. coli,
# H. salinarum and S. acidocaldarius, three different winning candidates).
# Measured on E. coli: 41.14 s -> 32.46 s, -21.1%, same 7,965,683 B archive.
# ── 8-POINT GRID: 4 x MAXMAP, 2 x MINOV ─────────────────────────────────────
# WIDER, AND NEARLY FREE. The candidates now run concurrently (see the fork
# loop in 106_inprocess.cpp), and the per-candidate work is ~83% serial, so
# wall(N) = S + N*Q/P with S ~= 20Q/P. Doubling N therefore costs
# (20+8)/(20+4) = 1.17x in theory; measured 1.12x. Two datasets:
#     SARS-CoV-2    4-grid 11.37 s -> 8-grid 12.73 s   archive IDENTICAL
#     H. salinarum  4-grid 11.84 s -> 8-grid 13.14 s   15,654,489 -> 15,650,029 B
# Both are still ~3.1x faster than the ORIGINAL sequential 4-grid (40.89 s),
# while searching twice the space.
#
# IT CANNOT BE WORSE. The 8-point grid CONTAINS the 4-point grid, and the
# encoder keeps whichever candidate produces the smallest archive, so
# min(8 points) <= min(4 points) by construction -- a guarantee, not a
# measurement that might not transfer.
#
# WHY THESE FOUR MAXMAP POINTS. This file already recorded the observed optima
# as spanning L/19 .. L/5.6, but the 2-point grid only covered L/13 .. L/5 and
# missed the LOW end -- and H. salinarum's optimum is documented right there
# ("halo ~8 = L/19"). That is exactly the dataset the wider grid improves, and
# exactly the failure a wider grid should fix.
C3=$(( LMAX / 19 )); [ "$C3" -lt 6 ] && C3=6
C4=$(( LMAX / 8  )); [ "$C4" -lt 6 ] && C4=6
# ── DEFAULT IS THE 4-POINT GRID, MEASURED ON 14 DATASETS 2026-09-05 ─────────
# The 8-point grid costs 56% of encoder runtime (HG002: 497 s for 8 candidates
# against 221 s for one) and buys almost nothing. Full measurement in
# docs/GRID_COST_MEASURED.txt -- archive cost of the 4-point grid vs 8-point:
#
#   +0.000% on 10 of 14 (incl. HG002, S. cerevisiae, S. aureus, P. falciparum,
#                        HCMV, L. major, M. tuberculosis, P. aeruginosa,
#                        S. acidocaldarius, SARS-CoV-2)
#   +0.012% E. coli   +0.028% H. salinarum   +0.032% H. pylori
#   +0.357% A. fumigatus   <-- the one real outlier, disclosed not hidden
#   mean +0.031%
#
# A SINGLE point was measured too and is NOT acceptable: +0.028% to +0.782%,
# mean +0.41%, swinging 28x between datasets -- exactly as this file's own note
# predicts ("the optimal mapping ceiling varies 5.6x across datasets ... no
# constant and no fixed ratio of read length can express it").
#
# The MINOV dimension is NOT reducible: both values win on real datasets
# (MINOV=16 on 7 of 12, the derived 0.35*Lmax on 5) and dropping the loser costs
# up to +0.629%. That is why 2 candidate groups remain.
# ── 2-POINT GRID: C2 x BOTH MINOV. MEASURED ON 14 DATASETS ─────────────────
# Within the 4-point grid, the permissive ceiling C2 (=Lmax/5) wins outright:
#
#   keeping ONLY C2 costs +0.000% on 13 of 14 datasets
#   (HG002, S. cerevisiae, S. acidocaldarius, S. aureus, SARS-CoV-2,
#    M. tuberculosis, E. coli, HCMV, L. major, P. falciparum, A. fumigatus,
#    H. pylori, P. aeruginosa)
#   keeping ONLY C2 costs +0.444% on H. salinarum -- the ONE dataset that
#   prefers C1, disclosed rather than hidden.  mean +0.032%
#   (keeping only C1 instead costs +0.056% .. +0.782% on 13 of 14 -- far worse)
#
# MINOV is NOT reducible and both values stay: MINOV=16 wins on 7 of 12 and the
# derived 0.35*Lmax on 5 of 12, and dropping the loser costs up to +0.629%.
# Two candidates is therefore the FLOOR of this grid, not one.
CANDS="${C2}:${M1},${C2}:${M2}"
# GRID4=1 restores the 4-point grid (adds C1=L/13).
[ -n "${GRID4:-}" ] && CANDS="${C1}:${M1},${C1}:${M2},${C2}:${M1},${C2}:${M2}"
# GRID8=1 restores the wider 8-point grid (adds C3=L/19 and C4=L/8).
[ -n "${GRID8:-}" ] && CANDS="${C3}:${M1},${C3}:${M2},${C1}:${M1},${C1}:${M2},${C4}:${M1},${C4}:${M2},${C2}:${M1},${C2}:${M2}"
CANDIDATES="$CANDS" ARCHIVE="$OUT" \
    DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 "$BEST" "$IN" 3 16 16 22 16 16 1 24 64 1 \
    > /dev/null 2>"${OUT}.log"
sz=$(grep -o 'ARCHIVE_TOTAL=[0-9][0-9]*' "${OUT}.log" | sed 's/ARCHIVE_TOTAL=//' | tail -1)
[ -z "$sz" ] && { echo "[adaptive] all candidates failed" >&2; exit 1; }
sed -n 's/.*\[a3\] chose //p' "${OUT}.log" >&2
echo "ARCHIVE_TOTAL=$sz"
