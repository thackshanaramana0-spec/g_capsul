#!/usr/bin/env bash
# Mismatch-tolerance sweep on REAL reads.
#
# The nine sites the mm=2 arm still misses have 3.5x the local GIAB variant
# density of the sites that succeed, and the search returns 3.2 occurrences
# there against 14.5 at successful sites. That says the alternate-haplotype
# contig differs from the probe in MORE than two places, so k=2 cannot reach
# it. If that reading is right, raising k recovers those sites.
#
# It cannot be a free win: a larger k also admits more spurious matches, so the
# homozygous control is swept alongside at every k. The useful quantity is the
# SEPARATION between the two, not the het number on its own.
set -u
cd ~/t34real
export CAPS_SPANS=1
DEC=~/bin/capsule_decode
ARC=realreads.capsule
SC=~/gc/scripts/score_locus_fidelity_v2.py

run_set(){   # $1 probes  $2 outdir  $3 k
  rm -rf "$2"; mkdir -p "$2"
  cat "$1" | xargs -P 6 -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r POS R A P V <<< "{}"
    CAPS_QUERY_MM=0 '"$DEC"' query '"$ARC"' '"$2"'/s0_$POS.fa "$P" >/dev/null 2>'"$2"'/s0_$POS.occ
    CAPS_QUERY_MM='"$3"' '"$DEC"' query '"$ARC"' '"$2"'/s1_$POS.fa "$P" >/dev/null 2>'"$2"'/s1_$POS.occ
    '"$DEC"' query '"$ARC"' '"$2"'/c_$POS.fa $((V-100))-$((V+100)) >/dev/null 2>&1'
}

cell(){      # $1 probes  $2 qdir  $3 label  -> prints "both"
  cp "$1" probes.tsv
  rm -rf q; cp -r "$2" q
  python3 "$SC" "$3" 2 | grep '^T34V2' | cut -d, -f15
}

printf '%6s %14s %22s %12s\n' k "het both/400" "homozygous false/400" separation
for K in 1 2 3 4 6 8; do
  run_set probes_real.tsv     "qk$K"  "$K"
  run_set probes_real_neg.tsv "qnk$K" "$K"
  H=$(cell probes_real.tsv     "qk$K"  HET)
  G=$(cell probes_real_neg.tsv "qnk$K" NEG)
  printf '%6s %14s %22s %12s\n' "$K" "$H" "$G" "$((H-G))"
done
