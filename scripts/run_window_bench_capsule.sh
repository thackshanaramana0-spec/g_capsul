#!/usr/bin/env bash
# CAPSULE reference-free het-SNV/indel benchmark on ONE chr20 window of ONE
# GIAB individual, scored by rtg vcfeval (the GA4GH engine hap.py wraps) —
# the same methodology and the same five windows the outer ARCS project used,
# so numbers are directly comparable to its published 0.936 het-SNV F1.
#
# Reports precision / recall / F1 for SNV and INDEL separately, like the
# DiscoSNP++ and Kmer2SNP papers do.
#
# WHY A WINDOW: full chr20 at 30x is ~4 GB and ~40 min per call. A 400 kb
# window is ~70k reads and runs in about a minute, which is what makes an
# iterate-debug-rerun loop possible at all.
#
# WHY THIS DOES NOT OVERFIT: r2 is the DESIGNATED TUNING WINDOW (the outer
# project tuned its frozen parameters on it once, years of runs ago) and
# r3/na/r4/r5 are held out. The caller's parameters are frozen constants
# ported verbatim (HDMAX/MAF/DHI/KHI/MC/TRI/HALF in include/caps_caller.h)
# and MUST NOT be edited in response to what any window scores. Iterate on
# r2 for CORRECTNESS only; report the held-out windows.
#
#   usage: run_window_bench_capsule.sh <capsule_exe> <scripts_dir> <ref.fa> \
#                                      <individual> <window> [outdir]
#     individual: HG002 | HG003 | HG004 | HG005
#     window:     r2 | r3 | na | r4 | r5   (or explicit LO:HI)
set -e
T_START=$(date +%s)
log() { echo "[window] $(date '+%H:%M:%S') $*"; }

CAPS="$1"; SC="$2"; REF="$3"; IND="${4:-HG002}"; WIN="${5:-r2}"
OUT="${6:-$HOME/caps_win/${IND}_${WIN}}"
CHROM=20
TARGET_COV="${TARGET_COV:-30}"     # standardized depth (source BAMs are 300x)

# ── the five windows: r2 is the tuning window, the rest are held out ────────
case "$WIN" in
  r2) LO=3000000; HI=3400000;;
  r3) LO=4000000; HI=4400000;;
  na) LO=2000000; HI=2400000;;
  r4) LO=5000000; HI=5400000;;
  r5) LO=6000000; HI=6400000;;
  *)  LO="${WIN%%:*}"; HI="${WIN##*:}";;
esac
REGION="$CHROM:$LO-$HI"

# ── public GIAB 300x GRCh37 BAMs (streamed by region; nothing pre-downloaded) ─
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data
case "$IND" in
  HG002) BAM="$G/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x_chr20.bam";;
  HG003) BAM="$G/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam";;
  HG004) BAM="$G/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam";;
  HG005) BAM="$G/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.hs37d5.300x.bam";;
  *) echo "unknown individual $IND" >&2; exit 1;;
esac
TRUTH="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
[ -s "$BED" ] || BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.bed"

mkdir -p "$OUT"; cd "$OUT"
echo "=== $IND $WIN ($REGION), target ${TARGET_COV}x ==="

# ── 1. Stream the window and downsample 300x -> TARGET_COV ─────────────────
FRAC=$(awk -v t="$TARGET_COV" 'BEGIN{printf "%.4f", t/300}')
if [ ! -s reads.fq ]; then
  samtools view -h -s "$FRAC" "$BAM" "$REGION" 2>/dev/null \
    | samtools fastq -n - 2>/dev/null > reads.fq
fi
NR_=$(( $(wc -l < reads.fq) / 4 ))
echo "reads: $NR_"
[ "$NR_" -gt 1000 ] || { echo "too few reads streamed — network or region problem" >&2; exit 1; }

# ── 2. CAPSULE reference-free call (frozen params; single fixed candidate) ──
export CAPS_CALL=1 CALL_VCF="$OUT/calls.vcf" CAPS_DUMP_CONTIGS="$OUT/contigs.tsv"
"$CAPS" reads.fq 3 16 16 22 16 16 1 24 64 1 > /dev/null 2> capsule.log || true
grep -E "CAPS-CALL" capsule.log || true
[ -s contigs.tsv ] || { echo "no contigs dumped" >&2; exit 1; }
# CAPS_DUMP_CONTIGS already writes FASTA (">contig_N\n<seq>"), unlike ARCS's
# caller which writes TSV and needs an awk conversion. Copy, do not convert —
# running the TSV awk over FASTA silently produces a garbage "contigs.fa"
# whose records all fail to lift (symptom: "lifted 0 calls").
cp contigs.tsv contigs.fa

# ── 3. Place contigs on the reference (EVALUATION ONLY — not part of calling) ─
bwa mem -t "$(nproc)" "$REF" contigs.fa 2>/dev/null > c2r.sam

# ── 4. Lift contig-coordinate calls to genome coordinates ──────────────────
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam "$REF" $CHROM lifted.vcf contigs.fa

# ── 5. Truth for this window: het-only (reference-free sees only het), inside
#       the GIAB confident regions, split SNV / INDEL ───────────────────────
tabix -h "$TRUTH" "$REGION" 2>/dev/null | awk -v OFS='\t' '
  /^##contig/{next} /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
  /^#/{print; next} {print}' > truth.vcf
awk -v c=$CHROM -v lo=$LO -v hi=$HI '$1==c && $3>lo && $2<hi{
  s=($2>lo?$2:lo); e=($3<hi?$3:hi); if(e>s) print c"\t"s"\t"e}' "$BED" > regions.bed
echo "confident bp in window: $(awk '{n+=$3-$2}END{print n+0}' regions.bed)"

# SNV/INDEL classification. NOTE the ordering below: these run AFTER
# `bcftools norm -m -any` has split multi-allelic records, so $5 is always a
# SINGLE allele here and a plain length test is correct.
#
# BUG FIXED 2026-09-03 (docs/HET_INDEL_FRESH_SCAN.md): these filters used to
# run BEFORE normalisation, on the raw ALT field. A genuine multi-allelic SNV
# (`T -> C,A`, which CAPSULE emits natively and is its own T5.2 capability)
# has an ALT *string* of length 3, so `length($5)!=1` classified it as an
# INDEL; `bcftools norm -m -any` then split it into two SNV-shaped rows that
# sat in the INDEL call set and scored as indel false positives. Measured on
# HG002 r2: 4 of 16 indel FPs were multi-allelic SNVs misfiled this way.
# The bug could only ever penalise CAPSULE, because CAPSULE is the only tool
# in this comparison that emits true multi-allelic records at all --
# DiscoSNP++ emits separate biallelic rows (docs/CLAIM2_TABLES_AND_INDEL_SCAN.md
# T5.2), which the old filter classified correctly by accident.
# Normalising first also collapses records that only become identical after
# left-alignment -- a second observed FP source (two rows differing solely in
# reference-inherited case, both normalising onto 3041361).
snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
het()   { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1];
          if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
# normalise (split multi-allelics + left-align) FIRST, then classify, then
# de-duplicate rows that normalisation made identical. Applied identically to
# truth and to calls, so the comparison stays symmetric.
prep(){ (grep '^#' "$1"; grep -v '^#' "$1"|sort -k2,2n) > "$2.pre.vcf"
  bcftools norm -f "$REF" -m -any "$2.pre.vcf" 2>/dev/null | bcftools sort 2>/dev/null > "$2.n.vcf" \
    || cp "$2.pre.vcf" "$2.n.vcf"
  eval "$3" < "$2.n.vcf" > "$2.v.vcf"
  awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++' "$2.v.vcf" > "$2.s.vcf"
  bgzip -c "$2.s.vcf" > "$2.vcf.gz"
  tabix -f -p vcf "$2.vcf.gz"; }
prep truth.vcf  t_snv 'het | snv'
prep truth.vcf  t_ind 'het | indel'
prep lifted.vcf c_snv snv
prep lifted.vcf c_ind indel

# ── 6. Score with rtg vcfeval ──────────────────────────────────────────────
[ -d "$HOME/refs/chr20.sdf" ] && SDF="$HOME/refs/chr20.sdf" || { rm -rf sdf; rtg format -o sdf "$REF" >/dev/null 2>&1; SDF=sdf; }
score(){ rm -rf "e_$1"
  rtg vcfeval -b "$2" -c "$3" -t "$SDF" --squash-ploidy --bed-regions regions.bed \
      -o "e_$1" >/dev/null 2>&1 || true
  if [ -f "e_$1/summary.txt" ]; then
    awk -v n="$1" 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{
      printf "%-6s TP=%-5s FP=%-5s FN=%-5s  P=%.3f  R=%.3f  F1=%.3f\n",n,tp,fp,fn,p,r,f}' "e_$1/summary.txt"
  else echo "$1: no summary (see e_$1)"; fi; }

echo "======== $IND $WIN — CAPSULE reference-free — rtg vcfeval (GA4GH) ========"
score SNV   t_snv.vcf.gz c_snv.vcf.gz
score INDEL t_ind.vcf.gz c_ind.vcf.gz
echo "truth het-SNV: $(zcat t_snv.vcf.gz|grep -vc '^#')  truth het-indel: $(zcat t_ind.vcf.gz|grep -vc '^#')"
echo "capsule SNV calls: $(zcat c_snv.vcf.gz|grep -vc '^#')  indel calls: $(zcat c_ind.vcf.gz|grep -vc '^#')"
echo "==========================================================================="
log "Total elapsed: $(( $(date +%s) - T_START ))s"
