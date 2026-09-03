#!/usr/bin/env bash
# T5.3 — real tetraploid benchmark, one command.
#
# Method: Cooke, Wedge & Lunter, "Benchmarking small-variant genotyping in
# polyploids", Genome Research 2022 (PMC8805713). Two real diploid GIAB
# individuals' real reads are concatenated into one 4-copy sample; the truth
# is the union of their real GIAB v4.2.1 calls. Nothing simulated -- every
# read and every truth allele is real data from a real person. See
# docs/CLAIM2_TABLES_AND_INDEL_SCAN.md sec T5.3 for the full writeup.
#
# This replays, as a single script, the commands that were run by hand in
# /root/tetraploid_bench/ earlier this project (docs/HET_INDEL_FRESH_SCAN.md,
# results/claim2/t5_3_tetraploid.csv) -- that CSV's numbers should reproduce
# from this script within measurement noise.
#
#   usage: run_tetraploid_bench_capsule.sh <capsule_exe> <scripts_dir> <ref.fa> \
#                                          [IND_A] [IND_B] [PLOIDY] [CHROM:LO-HI] [outdir]
#   IND_A, IND_B: two of HG002/HG003/HG004/HG005 (default HG003 HG004)
#   PLOIDY: ploidy of the mixed sample = 2 * (number of individuals mixed), default 4
#   region: default chr20:3000000-3400000, the standard evaluation window
set -euo pipefail
export PATH=~/miniconda3/bin:$PATH
CAPS="$1"; SC="$2"; REF="$3"
INDA="${4:-HG003}"; INDB="${5:-HG004}"; PLOIDY="${6:-4}"
REGIONSPEC="${7:-20:3000000-3400000}"
OUT="${8:-$HOME/tetra_win/${INDA}_${INDB}}"
CHROM="${REGIONSPEC%%:*}"; RANGE="${REGIONSPEC#*:}"; LO="${RANGE%-*}"; HI="${RANGE#*-}"
TARGET_COV=30

log()  { echo "[tetraploid] $*"; }
mkdir -p "$OUT"; cd "$OUT"
log "mixing $INDA + $INDB reads at ${TARGET_COV}x each -> ${PLOIDY}-copy sample, region $CHROM:$LO-$HI"

# ── 1. Stream both individuals' windows from GIAB public S3/HTTP BAMs, mix ──
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data
declare -A BAM
BAM[HG002]="$G/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x_chr20.bam"
BAM[HG003]="$G/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam"
BAM[HG004]="$G/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam"
BAM[HG005]="$G/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.hs37d5.300x.bam"
FRAC=$(awk -v t="$TARGET_COV" 'BEGIN{printf "%.4f", t/300}')
REGION="$CHROM:$LO-$HI"
for IND in "$INDA" "$INDB"; do
    [ -s "r_${IND}.fq" ] || samtools view -h -s "$FRAC" "${BAM[$IND]}" "$REGION" 2>/dev/null \
        | samtools fastq -n - 2>/dev/null > "r_${IND}.fq"
    log "$IND window reads: $(( $(wc -l < "r_${IND}.fq") / 4 ))"
done
cat "r_${INDA}.fq" "r_${INDB}.fq" > reads.fq
log "mixed reads: $(( $(wc -l < reads.fq) / 4 ))"

# ── 2. Real per-individual GIAB truth, chr20 only ────────────────────────────
for IND in "$INDA" "$INDB"; do
    TRUTH="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
    [ -s "$TRUTH" ] || { echo "FAIL: $TRUTH not found -- see docs/SERVER_SETUP_AND_DOWNLOADS.md sec 3" >&2; exit 1; }
    bcftools view -r "$CHROM" "$TRUTH" -Oz -o "truth_${IND}_chr.vcf.gz" 2>/dev/null
    tabix -f -p vcf "truth_${IND}_chr.vcf.gz"
    bcftools norm -f "$REF" -m -any "truth_${IND}_chr.vcf.gz" -Oz -o "truth_${IND}_norm.vcf.gz" 2>/dev/null
    tabix -f -p vcf "truth_${IND}_norm.vcf.gz"
done

# ── 3. Confident regions = intersection of both individuals' BEDs ───────────
BEDA="$HOME/giab_truth/${INDA}_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
BEDB="$HOME/giab_truth/${INDB}_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
[ -s "$BEDA" ] || BEDA="$HOME/giab_truth/${INDA}_GRCh37_1_22_v4.2.1_benchmark.bed"
[ -s "$BEDB" ] || BEDB="$HOME/giab_truth/${INDB}_GRCh37_1_22_v4.2.1_benchmark.bed"
bedtools intersect -a "$BEDA" -b "$BEDB" | awk -v c="$CHROM" '$1==c' > tetra_confident_full.bed
awk -v c="$CHROM" -v lo="$LO" -v hi="$HI" '$1==c && $3>lo && $2<hi{
  s=($2>lo?$2:lo); e=($3<hi?$3:hi); if(e>s) print c"\t"s"\t"e}' tetra_confident_full.bed > regions.bed
log "confident bp in window: $(awk '{n+=$3-$2}END{print n+0}' regions.bed)"

# ── 4. Build the tetraploid union truth (v2: normalise first, drop nothing) ──
python3 "$SC/build_tetraploid_truth_v2.py" "truth_${INDA}_norm.vcf.gz" "truth_${INDB}_norm.vcf.gz" tetra_truth_v2.vcf
awk -v lo="$LO" -v hi="$HI" '/^#/{print;next} $2>lo && $2<hi' tetra_truth_v2.vcf > tetra_truth_v2_win.vcf
log "truth sites in window: $(grep -vc '^#' tetra_truth_v2_win.vcf)"

# ── 5. CAPSULE reference-free call at the mixed sample's ploidy ─────────────
export CAPS_CALL=1 CAPS_PLOIDY="$PLOIDY" CALL_VCF="$OUT/calls.vcf" CAPS_DUMP_CONTIGS="$OUT/contigs.tsv"
"$CAPS" reads.fq 3 16 16 22 16 16 1 24 64 1 >/dev/null 2>capsule_call.log
grep -E "CAPS-CALL" capsule_call.log || true
cp contigs.tsv contigs.fa
[ -s "$REF.bwt" ] || bwa index "$REF" 2>/dev/null
bwa mem -t"$(nproc)" "$REF" contigs.fa 2>/dev/null > c2r.sam
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam "$REF" "$CHROM" lifted.vcf contigs.fa

# ── 6. DiscoSNP++ on the identical mixed reads ───────────────────────────────
if command -v run_discoSnp++.sh >/dev/null; then
    echo "reads.fq" > fof.txt
    run_discoSnp++.sh -r fof.txt -k 31 -c 3 -D 100 -P 3 -b 0 -G "$REF" -T >disco.log 2>&1 || true
    D=$(ls discoRes_*coherent.vcf 2>/dev/null | head -1)
else
    D=""
fi
if [ -n "$D" ]; then
    (echo '##fileformat=VCFv4.2'; echo "##contig=<ID=$CHROM,length=63025520>";
     echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';
     printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
     grep -v "^#" "$D" | awk -F'\t' -v OFS='\t' -v c="$CHROM" '$1==c{
        split($10,g,":"); gt=g[1]; gsub(/\|/,"/",gt);
        print $1,$2-1,".",$4,$5,30,"PASS",".","GT",gt}' | sort -k2,2n) > d_raw.vcf
else
    log "SKIP: DiscoSNP++ not found on PATH -- CAPSULE-only run"
    printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n" > d_raw.vcf
fi

# ── 7. Score both tools, same corrected convention as run_window_bench_capsule.sh
SDF="$OUT/sdf"; [ -d "$SDF" ] || rtg format -o "$SDF" "$REF" >/dev/null 2>&1
snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
dedup() { awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++'; }
prep(){ (grep '^#' "$1"; grep -v '^#' "$1"|sort -k2,2n) > "$2.pre.vcf"
  bcftools norm -f "$REF" -m -any "$2.pre.vcf" 2>/dev/null | bcftools sort 2>/dev/null > "$2.n.vcf" || cp "$2.pre.vcf" "$2.n.vcf"
  eval "$3" < "$2.n.vcf" | dedup > "$2.s.vcf"; bgzip -cf "$2.s.vcf" > "$2.vcf.gz"; tabix -f -p vcf "$2.vcf.gz"; }
prep tetra_truth_v2_win.vcf t_snv 'snv'
prep tetra_truth_v2_win.vcf t_ind 'indel'
prep lifted.vcf              c_snv 'snv'
prep lifted.vcf              c_ind 'indel'
prep d_raw.vcf                d_snv 'snv'
prep d_raw.vcf                d_ind 'indel'

score(){ rm -rf "e_$1"; rtg vcfeval -b "$2" -c "$3" -t "$SDF" --squash-ploidy --bed-regions regions.bed -o "e_$1" >/dev/null 2>&1 || true
  [ -f "e_$1/summary.txt" ] && awk -v n="$1" 'NR>2{printf "%-14s TP=%-5s FP=%-5s FN=%-5s P=%.3f R=%.3f F1=%.3f\n",n,$3,$4,$5,$6,$7,$8}' "e_$1/summary.txt" || echo "$1: no summary"; }

echo "===== T5.3 TETRAPLOID ($INDA+$INDB, ploidy=$PLOIDY) ====="
score "CAPSULE_SNV"   t_snv.vcf.gz c_snv.vcf.gz
score "DiscoSNP++_SNV" t_snv.vcf.gz d_snv.vcf.gz
score "CAPSULE_INDEL" t_ind.vcf.gz c_ind.vcf.gz
score "DiscoSNP++_INDEL" t_ind.vcf.gz d_ind.vcf.gz
echo "=========================================================="
log "Results in: $OUT"
