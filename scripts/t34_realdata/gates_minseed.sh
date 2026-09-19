#!/usr/bin/env bash
# CAPS_QUERY_MINSEED must be inert unless set.
set -u
cd ~/t34real
export CAPS_SPANS=1

echo "=== GATE A: default path byte-identical to git HEAD, 200 sites ==="
rm -rf mA mB; mkdir -p mA mB
head -200 probes_real.tsv > m200.tsv
while IFS=$'\t' read -r POS R A P V; do
  ~/bin/capsule_decode_base query realreads.capsule mA/$POS.fa "$P" >/dev/null 2>&1
  ~/bin/ms_decode            query realreads.capsule mB/$POS.fa "$P" >/dev/null 2>&1
done < m200.tsv
s=0; d=0
for f in mA/*.fa; do b=$(basename "$f"); cmp -s "$f" "mB/$b" && s=$((s+1)) || d=$((d+1)); done
echo "  identical=$s differing=$d"

echo "=== GATE B: unset MINSEED == explicit 12 ==="
rm -rf mC; mkdir -p mC
while IFS=$'\t' read -r POS R A P V; do
  CAPS_QUERY_MINSEED=12 CAPS_QUERY_MM=2 ~/bin/ms_decode query realreads.capsule mC/$POS.fa "$P" >/dev/null 2>&1
done < m200.tsv
rm -rf mD; mkdir -p mD
while IFS=$'\t' read -r POS R A P V; do
  CAPS_QUERY_MM=2 ~/bin/ms_decode query realreads.capsule mD/$POS.fa "$P" >/dev/null 2>&1
done < m200.tsv
s=0; d=0
for f in mC/*.fa; do b=$(basename "$f"); cmp -s "$f" "mD/$b" && s=$((s+1)) || d=$((d+1)); done
echo "  identical=$s differing=$d"

echo "=== GATE C: bilateral headline unchanged with the new binary ==="
rm -rf qbl2; mkdir -p qbl2
cat bl_het.tsv | xargs -P 3 -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r POS R A U D <<< "{}"
  CAPS_QUERY_MM=2 ~/bin/ms_decode query realreads.capsule qbl2/up_$POS.fa   "$U" >/dev/null 2>qbl2/up_$POS.occ
  CAPS_QUERY_MM=2 ~/bin/ms_decode query realreads.capsule qbl2/down_$POS.fa "$D" >/dev/null 2>qbl2/down_$POS.occ'
s=0; d=0
for f in qbl_het/*.fa; do b=$(basename "$f"); cmp -s "$f" "qbl2/$b" && s=$((s+1)) || d=$((d+1)); done
echo "  query outputs identical=$s differing=$d"

echo "=== GATE D: archive still LOSSLESS through the new binary ==="
rm -rf decm; mkdir -p decm
~/bin/ms_decode realreads.capsule decm decm/reads.out >/dev/null 2>&1
cmp -s real_orig_seqs.txt decm/reads.out && echo "  LOSSLESS byte-identical" || echo "  FAIL"
rm -rf decm mA mB mC mD qbl2
