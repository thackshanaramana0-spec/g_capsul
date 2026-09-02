#!/usr/bin/env bash
# CAPSULE version of the outer ARCS project's run_polyploid_bench.sh.
#   usage: run_polyploid_bench_capsule.sh <workdir> <capsule_exe> <scripts_dir> <ploidy>
set -e
export PATH=~/miniconda3/bin:$PATH
WD="$1"; CAPS="$2"; SC="$3"; K="${4:-3}"; CHROM=chrPLD
cd "$WD"
export CAPS_CALL=1 CALL_VCF="$WD/calls.vcf" CAPS_DUMP_CONTIGS="$WD/contigs.tsv" CAPS_PLOIDY=$K
"$CAPS" reads.fq 3 16 16 22 16 16 1 24 64 1 > /dev/null 2>&1 || true
grep -E "CAPS-CALL" "$WD"/*.log 2>/dev/null || true
cp contigs.tsv contigs.fa   # CAPS_DUMP_CONTIGS already emits FASTA (see run_window_bench_capsule.sh)
[ -s ref.fa.bwt ] || bwa index ref.fa 2>/dev/null
bwa mem -t4 ref.fa contigs.fa 2>/dev/null > c2r.sam
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam ref.fa $CHROM lifted.vcf contigs.fa
rm -rf sdf; rtg format -o sdf ref.fa >/dev/null 2>&1
prep(){ (grep '^#' "$1"; grep -v '^#' "$1"|sort -k2,2n) | bgzip > "$2.vcf.gz"; tabix -f -p vcf "$2.vcf.gz"; }
prep truth.vcf  truth
prep lifted.vcf calls
rm -rf e_pld
rtg vcfeval -b truth.vcf.gz -c calls.vcf.gz -t sdf --squash-ploidy --bed-regions regions.bed -o e_pld >/dev/null 2>&1 || true
echo "======== CAPSULE POLYPLOID (k=$K) reference-free SNV — rtg vcfeval (--squash-ploidy) ========"
awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{printf "k=%s  TP=%s FP=%s FN=%s  P=%.3f R=%.3f F1=%.3f\n","'$K'",tp,fp,fn,p,r,f}' e_pld/summary.txt
echo "truth het SNV: $(zcat truth.vcf.gz|grep -vc '^#')   multiallelic truth: $(zcat truth.vcf.gz|grep -v '^#'|awk -F'\t' '$5~/,/'|wc -l)   calls: $(zcat calls.vcf.gz|grep -vc '^#')"
