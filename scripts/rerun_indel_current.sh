#!/bin/bash
# Run the CURRENT caller on all 8 cached windows and score under the corrected
# convention, against the same truth already built by rescore_all.sh (n_t.vcf.gz).
# Disco is scored against that identical truth for comparison. Nothing is re-run
# for disco (its call set is unchanged); only CAPSULE is regenerated, because the
# caps_lift.vcf in these dirs is stale relative to the current caller.
set -u
REF=~/refs/chr20.fa; SDF=~/refs/chr20.sdf
SC=/root/arcs-clean/c_star_pg_advance/scripts
BEST=/tmp/capsule_bin/best106
[ -x "$BEST" ] || bash "$SC/build106.sh" "$BEST" >/dev/null 2>&1
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
dedup() { awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++'; }
f1() { rm -rf /tmp/_ev2; rtg vcfeval -b "$1" -c "$2" -t "$SDF" --squash-ploidy \
        --bed-regions "$3" -o /tmp/_ev2 >/dev/null 2>&1
  [ -f /tmp/_ev2/summary.txt ] && awk -v t="$4" 'NR>2{printf "%s TP=%s FP=%s FN=%s P=%.3f R=%.3f F1=%.3f\n",t,$3,$4,$5,$6,$7,$8}' /tmp/_ev2/summary.txt || echo "$4 FAILED"
}
for d in ~/bench2/*/; do
  W=$(basename "$d"); cd "$d" || continue
  [ -f reads.fq ] && [ -f n_t.vcf.gz ] && [ -f regions.bed ] || { echo "### $W SKIP (missing inputs)"; continue; }
  echo "### $W"
  CAPS_CALL=1 CALL_VCF="$d/cur_calls.vcf" CAPS_DUMP_CONTIGS="$d/cur_contigs.tsv" \
    "$BEST" reads.fq 3 16 16 22 16 16 1 24 64 1 >/dev/null 2>cur_caps.log
  cp cur_contigs.tsv cur_contigs.fa 2>/dev/null
  [ -s "$REF.bwt" ] || bwa index "$REF" 2>/dev/null
  bwa mem -t4 "$REF" cur_contigs.fa 2>/dev/null > cur_c2r.sam
  python3 "$SC/lift_vcf.py" cur_calls.vcf cur_c2r.sam "$REF" 20 cur_lift.vcf cur_contigs.fa >/dev/null 2>&1
  (grep '^#' cur_lift.vcf; grep -v '^#' cur_lift.vcf|sort -k2,2n) > cur.pre.vcf
  bcftools norm -f "$REF" -m -any cur.pre.vcf 2>/dev/null | bcftools sort 2>/dev/null > cur.n.vcf || cp cur.pre.vcf cur.n.vcf
  indel < cur.n.vcf | dedup > cur.s.vcf
  bgzip -cf cur.s.vcf > cur.vcf.gz; tabix -f -p vcf cur.vcf.gz 2>/dev/null
  f1 n_t.vcf.gz cur.vcf.gz regions.bed "  CAPSULE(current,new) "
  f1 n_t.vcf.gz d_ind.vcf.gz regions.bed "  DISCO(new)           "
done
