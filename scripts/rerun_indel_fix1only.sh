#!/bin/bash
# FIX-1-ONLY: correct multi-allelic CALL classification, but keep the ORIGINAL
# truth definition (het() on raw truth, so GT=1/2 sites stay excluded exactly as
# they were in every previously-published number). This isolates the pure
# bug-fix from the benchmark-definition change, so the comparison cannot be
# accused of redefining the task in CAPSULE's favour.
set -u
REF=~/refs/chr20.fa; SDF=~/refs/chr20.sdf
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
dedup() { awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++'; }
f1() { rm -rf /tmp/_ev3; rtg vcfeval -b "$1" -c "$2" -t "$SDF" --squash-ploidy \
        --bed-regions "$3" -o /tmp/_ev3 >/dev/null 2>&1
  [ -f /tmp/_ev3/summary.txt ] && awk -v t="$4" 'NR>2{printf "%s TP=%s FP=%s FN=%s P=%.3f R=%.3f F1=%.3f\n",t,$3,$4,$5,$6,$7,$8}' /tmp/_ev3/summary.txt || echo "$4 FAILED"
}
for d in ~/bench2/*/; do
  W=$(basename "$d"); cd "$d" || continue
  [ -f cur_lift.vcf ] && [ -f o_t.vcf.gz ] || { echo "### $W SKIP"; continue; }
  echo "### $W (ORIGINAL truth = $(zcat o_t.vcf.gz|grep -vc '^#') het-indels)"
  # current caller calls, normalised+classified correctly, vs ORIGINAL truth
  (grep '^#' cur_lift.vcf; grep -v '^#' cur_lift.vcf|sort -k2,2n) > f1.pre.vcf
  bcftools norm -f "$REF" -m -any f1.pre.vcf 2>/dev/null | bcftools sort 2>/dev/null > f1.n.vcf || cp f1.pre.vcf f1.n.vcf
  indel < f1.n.vcf | dedup > f1.s.vcf
  bgzip -cf f1.s.vcf > f1.vcf.gz; tabix -f -p vcf f1.vcf.gz 2>/dev/null
  f1 o_t.vcf.gz f1.vcf.gz regions.bed "  CAPSULE(fix1,orig-truth) "
  f1 o_t.vcf.gz d_ind.vcf.gz regions.bed "  DISCO(orig-truth)        "
done
