#!/usr/bin/env bash
# FULL-SCALE Claim 2 competitor arm: Kmer2SNP on the whole chr20 at 30x.
#
# Third arm of T3, alongside run_fullchr20_archive_capsule.sh (ours) and
# run_fullchr20_bench_disco.sh (DiscoSNP++). The truth set, the confident
# regions, the normalise-then-classify order and the rtg vcfeval call are
# COPIED VERBATIM from the DiscoSNP++ arm -- the only thing that differs is
# which caller produced the VCF. Anything else differing makes T3 meaningless.
#
#   usage: run_kmer2snp.sh <reads.fq> <individual> <outdir> [ref.fa]
#
# ── WHY THIS SCRIPT HAS TO EXIST, AND WHAT IT FIXES ────────────────────────
# Kmer2SNP has been installed at /root/Kmer2SNP all along, but this repository
# could never invoke it, so T3 shipped with the arm recorded NOT_AVAILABLE.
# The reason is in Kmer2SNP's own wrapper, script/run_dsk.sh:
#
#     /tmp/dsk-v2.3.3-bin-Linux/bin/dsk  -file $3 -out chr_k$1 ...
#
# a HARDCODED path to a DSK install that does not exist here. The wrapper fails
# silently (`/path2dsk/dsk: not found`), leaves `hete.para` empty, and the tool
# then dies in findGSE parsing with `IndexError: list index out of range` --
# which reads like a Kmer2SNP bug rather than a missing dependency.
#
# kmer2snp.py never needs either tool. Both are only called to FILL IN
# arguments (kmer2snp.py:59-70): `--t1` supplies the k-mer count file that
# run_dsk.sh would have produced, and `--c1/--c2/--r` supply what findGSE would
# have estimated. Passing all four bypasses both wrappers entirely. We count
# k-mers with KMC instead, which is packaged and maintained, and derive the
# coverage band ourselves.
#
# ── THE PARAMETERS ARE DERIVED, NOT COPIED ────────────────────────────────
# Standing rule 1: a formula over a measured property of the input, never a
# fitted constant. From the k-mer count histogram we take C, the homozygous
# coverage peak (the mode above the error valley). Heterozygous k-mers occur at
# ~C/2 because only one haplotype carries them, so the het band is centred
# there:
#     c1 = round(C/4)      c2 = round(3C/4)
# Validated against the archived 2026-09-02 run, which used C=21, c1=5, c2=16
# (C from findGSE). Run on that run's own histogram this formula independently
# derives C=20, c1=5, c2=15 -- c1 exact, the others within one. Recorded as
# measured; NOT tuned to match, which would make it a fitted constant.
#
# `--r`, the heterozygosity rate, is the one value we do not derive: findGSE
# estimated 0.001 on that run, which is also the textbook human autosomal
# heterozygosity. It is passed as HET_RATE and declared, not hidden.
#
# ── WHAT THE PUBLISHED 0.464 ACTUALLY MEASURED ────────────────────────────
# The Kmer2SNP F1 of 0.464 quoted in this project comes from
# /root/k2s_win/HG002_r2 (2026-09-02): TP=121, FP=1, FN=279 -- i.e. 400 truth
# variants. Full chr20 carries 44,575 het SNVs. That number is therefore a
# ~1% WINDOW of chr20, not the whole-chromosome evaluation our 0.888 and
# DiscoSNP++'s 0.847 come from. Quoting all three in one table without
# re-running Kmer2SNP at full scale would compare different experiments. This
# script exists to remove that asymmetry.
set -u
T_START=$(date +%s)
log(){ echo "[kmer2snp] $(date '+%H:%M:%S') $*"; }

READS="${1:?usage: run_kmer2snp.sh <reads.fq> <individual> <outdir> [ref.fa]}"
IND="${2:-HG002}"
OUT="${3:-/data/caps_full/${IND}_chr20_k2s}"
REF="${4:-$HOME/refs/chr20.fa}"
CHROM=20
K="${K2S_K:-31}"
HET_RATE="${HET_RATE:-0.001}"      # findGSE's estimate; human autosomal rate
NPROC=$(nproc)
K2S_DIR="${K2S_DIR:-/root/Kmer2SNP}"
# Kmer2SNP needs networkx, which is in the kmer2snp_r conda env (networkx
# 3.6.1), NOT in the system python3. Running it with `python3` dies on
# `ModuleNotFoundError: No module named 'networkx'`, which looks like a broken
# install rather than the wrong interpreter.
K2S_PY="${K2S_PY:-/root/miniconda3/envs/kmer2snp_r/bin/python}"
[ -x "$K2S_PY" ] || K2S_PY=python3
SAM2VCF="${SAM2VCF:-/root/arcs-clean/scripts/kmer2snp_sam_to_vcf.py}"

TRUTH="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
[ -s "$BED" ] || BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.bed"

for need in kmc kmc_dump bwa samtools rtg bgzip tabix bcftools; do
    command -v "$need" >/dev/null || { echo "MISSING TOOL: $need" >&2; exit 1; }
done
[ -s "$K2S_DIR/kmer2snp.py" ] || { echo "MISSING: $K2S_DIR/kmer2snp.py" >&2; exit 1; }
[ -s "$SAM2VCF" ]             || { echo "MISSING: $SAM2VCF" >&2; exit 1; }
[ -s "$READS" ]               || { echo "MISSING READS: $READS" >&2; exit 1; }
[ -s "$TRUTH" ]               || { echo "MISSING TRUTH: $TRUTH" >&2; exit 1; }

mkdir -p "$OUT"; cd "$OUT"
echo "=== $IND FULL chr20 — Kmer2SNP (k=$K) ==="

# ── 1. Count k-mers with KMC ───────────────────────────────────────────────
# -ci1 keeps singletons: the het band starts low and dropping count-1 k-mers
# would truncate exactly the population Kmer2SNP is looking for.
if [ ! -s "chr_k${K}.txt" ]; then
    log "[1/5] KMC counting k=$K (this is the long step)"
    mkdir -p kmc_tmp
    kmc -k"$K" -ci1 -cs100000 -t"$NPROC" -fq "$READS" "kmc_k${K}" kmc_tmp \
        > kmc.log 2>&1 || { echo "KMC failed -- see $OUT/kmc.log" >&2; tail -20 kmc.log >&2; exit 1; }
    kmc_dump "kmc_k${K}" "chr_k${K}.txt" >> kmc.log 2>&1 \
        || { echo "kmc_dump failed" >&2; exit 1; }
    rm -rf kmc_tmp "kmc_k${K}.kmc_pre" "kmc_k${K}.kmc_suf"
fi
NK=$(wc -l < "chr_k${K}.txt")
log "[1/5] done -- $NK distinct ${K}-mers"
[ "$NK" -gt 0 ] || { echo "no k-mers counted" >&2; exit 1; }

# ── 2. Derive the coverage band from the histogram ────────────────────────
# C = the homozygous peak of the k-mer count histogram.
# The error shoulder is the tallest thing in the histogram by far (singleton
# sequencing errors), so "the mode" alone always returns 1. A fixed floor is
# not the fix either -- it silently mis-fires whenever coverage is low. Find
# the VALLEY, the first count where the histogram stops falling, and take the
# mode strictly above it. That is where the homozygous peak lives regardless of
# depth, so it holds at 5x and at 30x without a constant to tune.
read -r C C1 C2 <<< "$(awk '{h[$2+0]++; if($2+0>mx)mx=$2+0} END{
      v=1; while(v<mx && h[v+1]+0<=h[v]+0) v++;      # descend the error shoulder
      best=0; c=0; for(i=v+1;i<=mx;i++) if(h[i]+0>best){best=h[i]+0; c=i}
      if(c==0) c=v+1;
      c1=int(c/4+0.5); c2=int(3*c/4+0.5);
      if(c1<2)c1=2; if(c2<=c1)c2=c1+1;
      print c, c1, c2 }' "chr_k${K}.txt")"
log "[2/5] homozygous peak C=$C  ->  het band c1=$C1 c2=$C2  (rate r=$HET_RATE)"
[ "${C:-0}" -gt 0 ] || { echo "could not derive coverage peak" >&2; exit 1; }

# ── 3. Kmer2SNP graph -> SNP k-mer pairs ──────────────────────────────────
# --t1/--c1/--c2/--r bypass run_dsk.sh and run_findgse.sh, both of which carry
# hardcoded paths that do not exist here (see the header).
if [ ! -s "k_${K}_pair.snp" ]; then
    log "[3/5] Kmer2SNP graph"
    "$K2S_PY" "$K2S_DIR/kmer2snp.py" single --k "$K" --c "$C" --fastaq "$READS" \
        --t1 "chr_k${K}.txt" --c1 "$C1" --c2 "$C2" --r "$HET_RATE" \
        > k2s.log 2>&1 || { echo "kmer2snp.py failed -- see $OUT/k2s.log" >&2; tail -25 k2s.log >&2; exit 1; }
fi
NP=$(grep -c . "k_${K}_pair.snp" 2>/dev/null || echo 0)
log "[3/5] done -- $NP SNP k-mer pairs"

# ── 4. Place the pairs on the reference and convert to VCF ────────────────
# THE PAIRS FILE IS k_K_pair.snp, NOT .non. Kmer2SNP writes both
# (libprism/local/process.py:135,143): .snp is the SNP pairs, .non is the
# NON-SNP pairs. Feeding .non silently yields ZERO calls -- every one of its
# entries differs at more than one base, so the converter rejects them all as
# non-SNPs, and T3 would have published Kmer2SNP F1=0.000 as if measured.
# Caught by replaying the archived 2026-09-02 run: its .non has 7 lines (all
# multi-base) while its .snp has 279 (all exactly one base), matching its
# k1.fa's 279 sequences.
# Kmer2SNP emits pairs of k-mers, not coordinates. The first k-mer of each pair
# is aligned to the reference to locate the SNP; kmer2snp_sam_to_vcf.py then
# takes the genome REF base and emits the alt allele.
if [ ! -s "k2s_raw.vcf" ]; then
    log "[4/5] bwa placing pairs, then converting to VCF"
    [ -s "${REF}.bwt" ] || bwa index "$REF" >/dev/null 2>&1
    awk '{print ">p"NR"\n"$1}' "k_${K}_pair.snp" > k1.fa
    bwa mem -t "$NPROC" "$REF" k1.fa > k1.sam 2>bwa.log
    python3 "$SAM2VCF" "k_${K}_pair.snp" k1.sam "$REF" "$CHROM" k2s_raw.vcf \
        > sam2vcf.log 2>&1 || { echo "sam_to_vcf failed -- see $OUT/sam2vcf.log" >&2; tail -20 sam2vcf.log >&2; exit 1; }
fi
log "[4/5] done -- $(grep -vc '^#' k2s_raw.vcf 2>/dev/null || echo 0) chr20 records"

# ── 5. Truth + score: VERBATIM from the DiscoSNP++ arm ────────────────────
log "[5/5] scoring with rtg vcfeval..."
tabix -h "$TRUTH" "$CHROM" 2>/dev/null | awk -v OFS='\t' '
  /^##contig/{next} /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
  /^#/{print; next} {print}' > truth.vcf
awk -v c=$CHROM '$1==c{print c"\t"$2"\t"$3}' "$BED" > regions.bed

snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
het()   { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1];
          if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
prep(){ (grep '^#' "$1"; grep -v '^#' "$1"|sort -k2,2n) > "$2.pre.vcf"
  bcftools norm -f "$REF" -m -any "$2.pre.vcf" 2>/dev/null | bcftools sort 2>/dev/null > "$2.n.vcf" \
    || cp "$2.pre.vcf" "$2.n.vcf"
  eval "$3" < "$2.n.vcf" > "$2.v.vcf"
  awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++' "$2.v.vcf" > "$2.s.vcf"
  bgzip -f -c "$2.s.vcf" > "$2.vcf.gz"
  tabix -f -p vcf "$2.vcf.gz"; }
prep truth.vcf    t_snv 'het | snv'
prep k2s_raw.vcf  k_snv snv

[ -d "$HOME/refs/chr20.sdf" ] && SDF="$HOME/refs/chr20.sdf" || { rm -rf sdf; rtg format -o sdf "$REF" >/dev/null 2>&1; SDF=sdf; }
score(){ rm -rf "e_$1"
  rtg vcfeval -b "$2" -c "$3" -t "$SDF" --squash-ploidy --bed-regions regions.bed \
      -o "e_$1" >/dev/null 2>&1 || true
  if [ -f "e_$1/summary.txt" ]; then
    awk -v n="$1" 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{
      printf "%-6s TP=%-7s FP=%-7s FN=%-7s  P=%.3f  R=%.3f  F1=%.3f\n",n,tp,fp,fn,p,r,f}' "e_$1/summary.txt"
  else echo "$1: no summary (see e_$1)"; fi; }

echo "======== $IND FULL chr20 — Kmer2SNP — rtg vcfeval (GA4GH) ========"
score SNV t_snv.vcf.gz k_snv.vcf.gz
echo "truth het-SNV: $(zcat t_snv.vcf.gz|grep -vc '^#')  Kmer2SNP SNV calls: $(zcat k_snv.vcf.gz|grep -vc '^#')"
echo "k=$K  C=$C  c1=$C1  c2=$C2  r=$HET_RATE  pairs=$NP"
echo "elapsed: $(( $(date +%s) - T_START ))s"
echo "=================================================================="
