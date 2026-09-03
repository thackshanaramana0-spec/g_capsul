#!/usr/bin/env bash
# T5.2 — real multi-allelic benchmark, one command.
#
# KEY FINDING THIS REPLICATES (docs/POLYPLOID_BENCHMARK.md sec 1/1b): pooling
# multiple GIAB individuals does NOT create multi-allelic sites -- human SNPs
# are essentially always biallelic in the population (measured: 1 site in 763
# pooling all four samples). The real source of multi-allelic sites is a
# SINGLE diploid individual's own GT=1/2 truth records, where the two
# haplotypes carry two DIFFERENT non-reference alleles (REF=C, ALT=A,G) --
# HG002 alone has 952 such sites on chr20. This script uses ONE individual,
# not a pool.
#
#   usage: run_multiallelic_bench_capsule.sh <capsule_exe> <scripts_dir> <ref.fa> \
#                                            [IND] [CHROM:LO-HI] [outdir]
#   default: HG002, chr20:1000000-6000000 (the primary region;
#            re-run with chr20:6000000-26000000 for the independent check)
set -euo pipefail
export PATH=~/miniconda3/bin:$PATH
CAPS="$1"; SC="$2"; REF="$3"
IND="${4:-HG002}"; REGIONSPEC="${5:-20:1000000-6000000}"
OUT="${6:-$HOME/multiallelic_bench/${IND}_$(echo "$REGIONSPEC" | tr ':' '_')}"
CHROM="${REGIONSPEC%%:*}"; RANGE="${REGIONSPEC#*:}"; LO="${RANGE%-*}"; HI="${RANGE#*-}"
TARGET_COV=30

log() { echo "[multiallelic] $(date '+%H:%M:%S') $*"; }
T_START=$(date +%s)
mkdir -p "$OUT"; cd "$OUT"
log "$IND, region $CHROM:$LO-$HI, target ${TARGET_COV}x, single individual (no pooling)"

# ── 1. Stream this individual's real reads for the region ──────────────────
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data
declare -A BAM
BAM[HG002]="$G/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x_chr20.bam"
BAM[HG003]="$G/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam"
BAM[HG004]="$G/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam"
BAM[HG005]="$G/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.hs37d5.300x.bam"
FRAC=$(awk -v t="$TARGET_COV" 'BEGIN{printf "%.4f", t/300}')
REGION="$CHROM:$LO-$HI"
[ -s reads.fq ] || samtools view -h -s "$FRAC" "${BAM[$IND]}" "$REGION" 2>/dev/null \
    | samtools fastq -n - 2>/dev/null > reads.fq
log "reads: $(( $(wc -l < reads.fq) / 4 ))"

# ── 2. Real truth for this individual, het sites (includes GT=1/2 naturally)
TRUTH="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
[ -s "$BED" ] || BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.bed"
[ -s "$TRUTH" ] || { echo "FAIL: $TRUTH not found -- see docs/SERVER_SETUP_AND_DOWNLOADS.md sec 3" >&2; exit 1; }
bcftools view -r "$REGION" "$TRUTH" -Oz -o truth_raw.vcf.gz 2>/dev/null
tabix -f -p vcf truth_raw.vcf.gz
awk -v c="$CHROM" -v lo="$LO" -v hi="$HI" '$1==c && $3>lo && $2<hi{
  s=($2>lo?$2:lo); e=($3<hi?$3:hi); if(e>s) print c"\t"s"\t"e}' "$BED" > regions.bed
log "confident bp in region: $(awk '{n+=$3-$2}END{print n+0}' regions.bed)"
N_MULTI=$(zcat truth_raw.vcf.gz | grep -v '^#' | awk -F'\t' '$5~/,/' | wc -l)
log "truth multi-allelic sites in region: $N_MULTI"

# ── 3. CAPSULE reference-free call (standard diploid, no ploidy override) ──
export CAPS_CALL=1 CALL_VCF="$OUT/calls.vcf" CAPS_DUMP_CONTIGS="$OUT/contigs.tsv"
"$CAPS" reads.fq 3 16 16 22 16 16 1 24 64 1 >/dev/null 2>capsule_call.log
grep -E "CAPS-CALL" capsule_call.log || true
cp contigs.tsv contigs.fa
[ -s "$REF.bwt" ] || bwa index "$REF" 2>/dev/null
bwa mem -t"$(nproc)" "$REF" contigs.fa 2>/dev/null > c2r.sam
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam "$REF" "$CHROM" lifted.vcf contigs.fa

# ── 4. DiscoSNP++ on the identical reads ─────────────────────────────────────
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

# ── 5. Extract truth's multi-allelic sites specifically, then check recovery
zcat truth_raw.vcf.gz | awk -F'\t' '/^#/{print;next} $5~/,/' > truth_multiallelic.vcf
N_TRUTH_MULTI=$(grep -vc '^#' truth_multiallelic.vcf)
log "scoring recovery of the $N_TRUTH_MULTI truth multi-allelic sites specifically"

recovered() {
    # $1 = calls file (lifted.vcf or d_raw.vcf).
    #
    # CAVEAT, stated plainly: this counts truth multi-allelic POSITIONS where
    # the caller emitted ANY call at that exact position -- it is a proxy for
    # "attempted this site," not the stricter check the original 11/18 vs
    # 0/18 analysis used (whether BOTH non-reference alleles were correctly
    # captured, either as one native multi-allelic record or two matching
    # biallelic ones). That original check was done as ad-hoc analysis, not a
    # saved script, so this rewrite cannot be verified byte-for-byte against
    # it. Treat this script's output as indicative, re-derive the strict
    # allele-level check by hand from lifted.vcf/d_raw.vcf against
    # truth_multiallelic.vcf before quoting a number in the paper.
    local CALLS="$1"
    awk -F'\t' 'NR==FNR{if($1!~/^#/) pos[$2]=1; next}
        $1!~/^#/ && ($2 in pos){print $2}' truth_multiallelic.vcf "$CALLS" \
        | sort -u | wc -l
}
CAPS_HIT=$(recovered lifted.vcf)
DISCO_HIT=$(recovered d_raw.vcf)

echo "===== T5.2 MULTI-ALLELIC ($IND, $CHROM:$LO-$HI) ====="
echo "truth multi-allelic sites: $N_TRUTH_MULTI"
echo "CAPSULE    sites with a call at that position: $CAPS_HIT / $N_TRUTH_MULTI"
echo "DiscoSNP++ sites with a call at that position: $DISCO_HIT / $N_TRUTH_MULTI"
echo "======================================================="
log "Results in: $OUT"
log "Total elapsed: $(( $(date +%s) - T_START ))s"
