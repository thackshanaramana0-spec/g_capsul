#!/usr/bin/env bash
# The multi-probe change must be an EFFICIENCY change only. If it alters a
# single occurrence it is a behaviour change and has to be treated as one.
#
# Check: for 100 real sites, the occurrences found by ONE invocation with two
# comma-separated probes must equal the union of the occurrences found by TWO
# separate invocations. Compared as sorted sets of (start,end,strand).
set -u
cd ~/t34real
export CAPS_SPANS=1
OLD=~/bin/capsule_decode      # one probe per invocation
NEW=~/bin/mp_decode           # comma-separated probes
MM=2

head -100 bl_het.tsv > mp100.tsv
rm -rf mpA mpB && mkdir -p mpA mpB

while IFS=$'\t' read -r POS R A U D; do
  CAPS_QUERY_MM=$MM $OLD query realreads.capsule mpA/u_$POS.fa "$U" >/dev/null 2>mpA/u_$POS.occ
  CAPS_QUERY_MM=$MM $OLD query realreads.capsule mpA/d_$POS.fa "$D" >/dev/null 2>mpA/d_$POS.occ
done < mp100.tsv

while IFS=$'\t' read -r POS R A U D; do
  CAPS_QUERY_MM=$MM $NEW query realreads.capsule mpB/$POS.fa "$U,$D" >/dev/null 2>mpB/$POS.occ
done < mp100.tsv

same=0; diff=0
while IFS=$'\t' read -r POS R A U D; do
  a=$( { grep -oE 'occ [0-9]+ [0-9]+ [+-]' mpA/u_$POS.occ;
         grep -oE 'occ [0-9]+ [0-9]+ [+-]' mpA/d_$POS.occ; } | sort -u )
  b=$( grep -oE 'occ [0-9]+ [0-9]+ [+-]' mpB/$POS.occ | sort -u )
  if [ "$a" == "$b" ]; then same=$((same+1)); else diff=$((diff+1)); fi
done < mp100.tsv
echo "occurrence sets: identical=$same differing=$diff"

echo
echo "=== speed: 100 sites, two probes each ==="
/usr/bin/time -f "  two invocations : %e s  peak %M kB" bash -c '
  while IFS=$'"'"'\t'"'"' read -r POS R A U D; do
    CAPS_QUERY_MM=2 '"$OLD"' query realreads.capsule /tmp/x.fa "$U" >/dev/null 2>&1
    CAPS_QUERY_MM=2 '"$OLD"' query realreads.capsule /tmp/x.fa "$D" >/dev/null 2>&1
  done < mp100.tsv'
/usr/bin/time -f "  one invocation  : %e s  peak %M kB" bash -c '
  while IFS=$'"'"'\t'"'"' read -r POS R A U D; do
    CAPS_QUERY_MM=2 '"$NEW"' query realreads.capsule /tmp/x.fa "$U,$D" >/dev/null 2>&1
  done < mp100.tsv'
