#!/usr/bin/env bash
# Seed-floor sweep against the ORIENTATION-CORRECTED baseline.
#
# A first pass reported 0.9448 recall at every floor, but most of that shortfall
# was a bug in the measurement: reads are returned in pseudogenome orientation
# and were being compared as raw strings. Canonicalising both sides moves the
# baseline to 0.9789. The question is whether the remaining 2.1% is tolerance.
set -u
cd ~/t34real
export CAPS_SPANS=1
R=/mnt/c/Users/Lenovo/OneDrive/Desktop/g_capsul/g_capsul
head -50 bl_het.tsv > r50.tsv

printf '%9s %7s %10s %9s %9s\n' MINSEED k_max recall missed wall_s
for cfg in "12 2" "8 4" "6 5" "5 7"; do
  set -- $cfg; MS=$1; K=$2
  rm -rf qrec; mkdir -p qrec
  t0=$(date +%s.%N)
  while IFS=$'\t' read -r POS REF ALT U D; do
    CAPS_QUERY_MINSEED=$MS CAPS_QUERY_MM=$K ~/bin/ms_decode query realreads.capsule \
        qrec/up_$POS.fa   "$U" >/dev/null 2>&1
    CAPS_QUERY_MINSEED=$MS CAPS_QUERY_MM=$K ~/bin/ms_decode query realreads.capsule \
        qrec/down_$POS.fa "$D" >/dev/null 2>&1
  done < r50.tsv
  t1=$(date +%s.%N)
  out=$(QDIR=qrec PROBES=r50.tsv python3 "$R/scripts/t34_realdata/exact_match_recall.py" 50 2>/dev/null)
  rec=$(echo "$out" | grep 'recall' | head -1 | awk '{print $NF}')
  mis=$(echo "$out" | grep 'MISSED' | awk '{print $NF}')
  printf '%9s %7s %10s %9s %9s\n' "$MS" "$K" "$rec" "$mis" \
     "$(echo "$t1 $t0" | awk '{printf "%.1f",$1-$2}')"
done
