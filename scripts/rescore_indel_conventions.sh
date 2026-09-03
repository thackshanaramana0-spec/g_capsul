#!/bin/bash
# Re-score all 8 het-indel evaluations under BOTH the old and the corrected
# truth/classification convention, for BOTH tools, from already-existing call
# sets (no tool is re-run, so nothing about either caller changes).
REF=~/refs/chr20.fa; SDF=~/refs/chr20.sdf
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
het()   { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1];
          if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
dedup() { awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++'; }

f1() { # $1=truth.gz $2=calls.gz $3=regions $4=tag
  rm -rf /tmp/_ev; rtg vcfeval -b "$1" -c "$2" -t "$SDF" --squash-ploidy \
      --bed-regions "$3" -o /tmp/_ev >/dev/null 2>&1
  if [ -f /tmp/_ev/summary.txt ]; then
    awk -v t="$4" 'NR>2{printf "%s TP=%s FP=%s FN=%s P=%.3f R=%.3f F1=%.3f\n",t,$3,$4,$5,$6,$7,$8}' /tmp/_ev/summary.txt
  else echo "$4 FAILED"; fi
}

printf "%-12s %-38s %-38s\n" "WINDOW" "OLD convention" "NEW convention"
for d in ~/bench2/*/; do
  W=$(basename "$d"); cd "$d" || continue
  [ -f truth.vcf ] && [ -f regions.bed ] && [ -f caps_lift.vcf ] && [ -f d_ind.vcf.gz ] || continue

  # ---- OLD truth: filter first, then norm (with the original silent fallback)
  het < truth.vcf | indel > o_t.v.vcf
  (grep '^#' o_t.v.vcf; grep -v '^#' o_t.v.vcf|sort -k2,2n) > o_t.s.vcf
  bcftools norm -f "$REF" -m -any o_t.s.vcf 2>/dev/null | bcftools sort 2>/dev/null | bgzip > o_t.vcf.gz \
    || bgzip -c o_t.s.vcf > o_t.vcf.gz
  tabix -f -p vcf o_t.vcf.gz 2>/dev/null

  # ---- NEW truth: norm first (splits multi-allelic), then filter, then dedupe
  (grep '^#' truth.vcf; grep -v '^#' truth.vcf|sort -k2,2n) > n_t.pre.vcf
  bcftools norm -f "$REF" -m -any n_t.pre.vcf 2>/dev/null | bcftools sort 2>/dev/null > n_t.n.vcf \
    || cp n_t.pre.vcf n_t.n.vcf
  het < n_t.n.vcf | indel | dedup > n_t.s.vcf
  bgzip -cf n_t.s.vcf > n_t.vcf.gz; tabix -f -p vcf n_t.vcf.gz 2>/dev/null

  # ---- CAPSULE calls, OLD convention (filter then norm-with-fallback)
  indel < caps_lift.vcf > o_c.v.vcf
  (grep '^#' o_c.v.vcf; grep -v '^#' o_c.v.vcf|sort -k2,2n) > o_c.s.vcf
  bcftools norm -f "$REF" -m -any o_c.s.vcf 2>/dev/null | bcftools sort 2>/dev/null | bgzip > o_c.vcf.gz \
    || bgzip -c o_c.s.vcf > o_c.vcf.gz
  tabix -f -p vcf o_c.vcf.gz 2>/dev/null

  # ---- CAPSULE calls, NEW convention (norm first, then filter, then dedupe)
  (grep '^#' caps_lift.vcf; grep -v '^#' caps_lift.vcf|sort -k2,2n) > n_c.pre.vcf
  bcftools norm -f "$REF" -m -any n_c.pre.vcf 2>/dev/null | bcftools sort 2>/dev/null > n_c.n.vcf \
    || cp n_c.pre.vcf n_c.n.vcf
  indel < n_c.n.vcf | dedup > n_c.s.vcf
  bgzip -cf n_c.s.vcf > n_c.vcf.gz; tabix -f -p vcf n_c.vcf.gz 2>/dev/null

  # disco: emits no multi-allelic records, and bcftools norm fails on its VCF
  # (strand-dependent anchoring), so its prepared indel set is the same either way.
  echo "### $W  (truth old=$(zcat o_t.vcf.gz|grep -vc '^#')  new=$(zcat n_t.vcf.gz|grep -vc '^#'))"
  f1 o_t.vcf.gz o_c.vcf.gz regions.bed "  CAPS-old  "
  f1 n_t.vcf.gz n_c.vcf.gz regions.bed "  CAPS-new  "
  f1 o_t.vcf.gz d_ind.vcf.gz regions.bed "  DISCO-old "
  f1 n_t.vcf.gz d_ind.vcf.gz regions.bed "  DISCO-new "
done
