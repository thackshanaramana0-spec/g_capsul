#!/usr/bin/env bash
# FINAL REGRESSION GATES.
#
# The multi-probe change altered the stderr format and added a dedup on the
# range list, so byte-identity of the ANSWER has to be re-proved rather than
# reasoned about.
set -u
cd ~/t34real
export CAPS_SPANS=1
SP=/mnt/c/Users/Lenovo/AppData/Local/Temp/claude/C--Users-Lenovo-OneDrive-Desktop-g-capsul-g-capsul/0f9f1ca8-c5b5-4ca9-87fd-b59b9d046590/scratchpad

echo "=== GATE 1: default path .fa byte-identical to git HEAD, 200 sites ==="
rm -rf gA gB; mkdir -p gA gB
head -200 probes_real.tsv > g200.tsv
while IFS=$'\t' read -r POS R A P V; do
  ~/bin/capsule_decode_base query realreads.capsule gA/$POS.fa "$P" >/dev/null 2>&1
  ~/bin/mp_decode            query realreads.capsule gB/$POS.fa "$P" >/dev/null 2>&1
done < g200.tsv
s=0; d=0
for f in gA/*.fa; do b=$(basename "$f"); if cmp -s "$f" "gB/$b"; then s=$((s+1)); else d=$((d+1)); fi; done
echo "  identical=$s differing=$d"

echo "=== GATE 2: k=0 identical to exact std::string::find ==="
g++ -O2 -std=c++17 -o /tmp/tm "$SP/test_matcher.cpp" && /tmp/tm

echo "=== GATE 3: archive LOSSLESS through the new binary ==="
rm -rf decg; mkdir -p decg
~/bin/mp_decode realreads.capsule decg decg/reads.out >/dev/null 2>&1
cmp -s real_orig_seqs.txt decg/reads.out && echo "  LOSSLESS byte-identical" || echo "  FAIL"

echo "=== GATE 4: headline bilateral reproduces ==="
python3 ~/gc/scripts/score_bilateral.py 2>&1 | grep -E "BILATERAL \(union\)|^sites:"

echo "=== GATE 5: multi-probe == two invocations, occurrence sets ==="
same=0; diff=0
while IFS=$'\t' read -r POS R A U D; do
  a=$( { CAPS_QUERY_MM=2 ~/bin/mp_decode query realreads.capsule /tmp/a1.fa "$U" 2>&1 >/dev/null
         CAPS_QUERY_MM=2 ~/bin/mp_decode query realreads.capsule /tmp/a2.fa "$D" 2>&1 >/dev/null; } \
       | grep -oE 'occ [0-9]+ [0-9]+ [+-]' | sort -u )
  b=$( CAPS_QUERY_MM=2 ~/bin/mp_decode query realreads.capsule /tmp/b1.fa "$U,$D" 2>&1 >/dev/null \
       | grep -oE 'occ [0-9]+ [0-9]+ [+-]' | sort -u )
  if [ "$a" == "$b" ]; then same=$((same+1)); else diff=$((diff+1)); fi
done < <(head -50 bl_het.tsv)
echo "  identical=$same differing=$diff"
