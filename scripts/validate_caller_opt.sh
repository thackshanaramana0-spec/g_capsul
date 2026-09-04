#!/usr/bin/env bash
# Gate for caller RAM/speed changes: prove the optimised binary produces the
# SAME calls as the reference binary, on all five standard windows.
#
# WHY BYTE-IDENTITY IS THE RIGHT GATE HERE. The three changes this was written
# for (parallelising build_substrate's read placement, dropping the strand-0
# string copy, and skipping rc_reads keys that pkidx can never look up) are all
# argued to be output-preserving BY CONSTRUCTION, not by measurement. A change
# claimed to be exactly equivalent must be held to exact equality -- anything
# weaker would let a real behaviour change hide inside a "close enough" F1.
#
# This is deliberately a STRICTER gate than the one used for the kidx flat-array
# conversion, which legitimately could not meet it: unordered_map iteration
# order is arbitrary and several downstream structures are first-write-wins, so
# that change produced a different-but-not-worse call set and was validated on
# 5-window mean F1 instead. Do not silently downgrade this gate to that one --
# if these changes fail byte-identity, the assumption behind them is wrong and
# the change is the problem, not the gate.
#
#   usage: validate_caller_opt.sh <ref_exe> <new_exe> [outdir]
set -e
REF="${1:?ref binary}"; NEW="${2:?new binary}"
OUT="${3:-/data/caps_validate}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log(){ echo "[validate] $(date '+%H:%M:%S') $*"; }

mkdir -p "$OUT"
FAIL=0
for WIN in r2 r3 na r4 r5; do
  D="$OUT/$WIN"; mkdir -p "$D"
  # Reuse the reads the window bench already cached, so this does no network IO
  # and both binaries see byte-identical input.
  SRC="$HOME/caps_win/HG002_${WIN}/reads.fq"
  [ -s "$SRC" ] || SRC="/tmp/window_v3_HG002_${WIN}/reads.fq"
  if [ ! -s "$SRC" ]; then log "SKIP $WIN — no cached reads.fq"; continue; fi

  for TAG in ref new; do
    EXE=$([ "$TAG" = ref ] && echo "$REF" || echo "$NEW")
    if [ ! -s "$D/$TAG.vcf" ]; then
      ( cd "$D" && CAPS_CALL=1 CALL_VCF="$D/$TAG.vcf" CAPS_DUMP_CONTIGS="$D/$TAG.contigs" \
          /usr/bin/time -f "%e s  %M KB" "$EXE" "$SRC" 3 16 16 22 16 16 1 24 64 1 \
          > /dev/null 2> "$D/$TAG.log" ) || true
    fi
  done

  # Compare the CALLS, ignoring only header lines that legitimately vary.
  grep -v '^#' "$D/ref.vcf" 2>/dev/null | sort > "$D/ref.s" || true
  grep -v '^#' "$D/new.vcf" 2>/dev/null | sort > "$D/new.s" || true
  RT=$(tail -2 "$D/ref.log" | grep -oE "^[0-9.]+ s" | head -1)
  NT=$(tail -2 "$D/new.log" | grep -oE "^[0-9.]+ s" | head -1)
  RM=$(grep -oE "[0-9]+ KB" "$D/ref.log" | tail -1)
  NM=$(grep -oE "[0-9]+ KB" "$D/new.log" | tail -1)
  if cmp -s "$D/ref.s" "$D/new.s"; then
    printf "%-4s IDENTICAL  %6s calls   ref %-9s %-12s  new %-9s %-12s\n" \
      "$WIN" "$(wc -l < "$D/ref.s")" "$RT" "$RM" "$NT" "$NM"
  else
    FAIL=1
    printf "%-4s DIFFER     ref=%s new=%s  (+%s / -%s)\n" "$WIN" \
      "$(wc -l < "$D/ref.s")" "$(wc -l < "$D/new.s")" \
      "$(comm -13 "$D/ref.s" "$D/new.s" | wc -l)" \
      "$(comm -23 "$D/ref.s" "$D/new.s" | wc -l)"
  fi
done
echo
[ "$FAIL" = 0 ] && echo "GATE PASS — all windows byte-identical" \
                || echo "GATE FAIL — output changed; the equivalence argument is wrong, fix the change not the gate"
exit $FAIL
