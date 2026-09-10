#!/usr/bin/env bash
# CAPSULE version of the outer ARCS project's run_indel_bench.sh: synthetic
# diploid indel benchmark, scored by rtg vcfeval. Only the caller invocation
# differs from the ARCS original (env-var interface instead of `arcs call`);
# lift/scoring below is copied unchanged.
#   usage: run_indel_bench_capsule.sh <workdir> <capsule_exe> <scripts_dir>
set -e
WD="$1"; CAPS="$2"; SC="$3"
CHROM=chrSIM
cd "$WD"

# 1. CAPSULE reference-free call (fixed MAXMAP/MINOV — NOT encode_adaptive.sh's
#    4-candidate sweep, which would race 4 writers on the same CALL_VCF path;
#    see docs/_removable/CLAIM2_BUILDER.md). Args after reads.fq are the standard
#    ENC_ARGS positional set (3 16 16 22 16 16 1 24 64 1).
export CAPS_CALL=1 CALL_VCF="$WD/calls.vcf" CAPS_DUMP_CONTIGS="$WD/contigs.tsv"
"$CAPS" reads.fq 3 16 16 22 16 16 1 24 64 1 > /dev/null 2> capsule_call.log
grep -E "CAPS-CALL" capsule_call.log || true
cp contigs.tsv contigs.fa   # CAPS_DUMP_CONTIGS already emits FASTA (see run_window_bench_capsule.sh)

# 2. Place contigs on the reference with BWA (evaluation-only coordinate lift)
bwa index ref.fa 2> bwa_index.log
bwa mem -t 4 ref.fa contigs.fa 2> bwa_mem.log > c2r.sam

# 3. Lift contig-coordinate calls to genome coordinates
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam ref.fa $CHROM lifted.vcf contigs.fa

# 4. Split truth + calls into SNV-only and INDEL-only, score each with rtg vcfeval
#
# TWO BUGS FIXED 2026-09-03, found in a full script-suite audit
# (docs/HET_INDEL_FRESH_SCAN.md Finding 4; same audit that found the
# identical classify-before-normalize defect in run_giab_indel_capsule.sh
# and run_window_bench_capsule.sh):
#
# BUG 1 (this script only): `score()` read TP from rtg's summary.txt column
# $2 ("True-pos-baseline") while every other script in this repo correctly
# reads column $3 ("True-pos-call") -- confirmed against a real
# summary.txt header: "Threshold True-pos-baseline True-pos-call
# False-pos False-neg Precision Sensitivity F-measure". These two TP counts
# are USUALLY equal but are logically distinct columns; using the wrong one
# is a real defect even where it happens not to have changed a number yet.
#
# BUG 2 (shared with the other two scripts): classification ran BEFORE
# normalisation, on the raw ALT field. A multi-allelic call (ALT="C,A",
# string length 3) is misfiled as an INDEL by `length($5)!=1`; bcftools
# norm then splits it into SNV-shaped rows sitting in the indel call set as
# false positives. Fixed the same way: normalise first, classify, dedupe.
snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
dedup() { awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++'; }

rtg format -o sdf ref.fa > /dev/null 2>&1 || true
NORM=${NORM:-1}
prep() {
  if [ "$NORM" = 1 ] && command -v bcftools >/dev/null 2>&1; then
    (grep '^#' "$1"; grep -v '^#' "$1" | sort -k2,2n) > "$2.pre.vcf"
    bcftools norm -f ref.fa -c s "$2.pre.vcf" 2>/dev/null | bcftools sort 2>/dev/null > "$2.n.vcf" || cp "$2.pre.vcf" "$2.n.vcf"
  else
    cp "$1" "$2.n.vcf"
  fi
  $3 < "$2.n.vcf" | dedup > "$2.vcf"
  bgzip -f "$2.vcf"; tabix -f -p vcf "$2.vcf.gz"; }

score() {
  rm -rf "eval_$1"
  rtg vcfeval -b "$2" -c "$3" -t sdf --squash-ploidy --bed-regions regions.bed \
      -o "eval_$1" > /dev/null 2>&1 || true
  if [ -f "eval_$1/summary.txt" ]; then
    awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{printf "%-6s TP=%s FP=%s FN=%s  P=%.3f R=%.3f F1=%.3f\n","'$1'",tp,fp,fn,p,r,f}' "eval_$1/summary.txt"
  else echo "$1: vcfeval produced no summary (see eval_$1)"; fi
}

prep truth.vcf  truth_snv   snv
prep truth.vcf  truth_indel indel
prep lifted.vcf call_snv    snv
prep lifted.vcf call_indel  indel

echo "=================== CAPSULE synthetic indel — rtg vcfeval (gold-standard) ==================="
score SNV   truth_snv.vcf.gz   call_snv.vcf.gz
score INDEL truth_indel.vcf.gz call_indel.vcf.gz
echo "==============================================================================================="
