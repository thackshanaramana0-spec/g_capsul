#!/usr/bin/env bash
# FULL-SCALE Claim 2 competitor arm: DiscoSNP++ on the whole chr20 at 30x.
#
# Companion to run_fullchr20_bench_capsule.sh. Same reads, same truth, same
# confident regions, same normalise-then-classify order, same rtg vcfeval
# call — the ONLY difference is which caller produced the VCF. Anything else
# differing would make the head-to-head meaningless.
#
# The DiscoSNP++ invocation is copied verbatim from the window benchmarks
# (run_multiallelic_bench_capsule.sh, run_tetraploid_bench_capsule.sh) so the
# competitor is not accidentally tuned differently at full scale than it was
# in the window results this is meant to be compared against:
#     run_discoSnp++.sh -r fof.txt -k 31 -c 3 -D 100 -P 3 -b 0 -G <ref> -T
# `-G` is mandatory per DEVNOTES.md; `-T` maps the bubbles back to the reference
# so the output carries genome coordinates. The `$2-1` in the POS conversion
# is DiscoSNP++'s documented off-by-one, corrected everywhere in this project.
#
#   usage: run_fullchr20_bench_disco.sh <ref.fa> <reads.fq> [individual] [outdir]
set -e
T_START=$(date +%s)
log() { echo "[disco-full] $(date '+%H:%M:%S') $*"; }

REF="$1"; READS="$2"; IND="${3:-HG002}"
OUT="${4:-/data/caps_full/${IND}_chr20_disco}"
CHROM=20

TRUTH="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
[ -s "$BED" ] || BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.bed"
command -v run_discoSnp++.sh >/dev/null || { echo "run_discoSnp++.sh not on PATH" >&2; exit 1; }

mkdir -p "$OUT"; cd "$OUT"
echo "=== $IND FULL chr20 — DiscoSNP++ ==="

# ── 1. DiscoSNP++ ──────────────────────────────────────────────────────────
log "[1/4] running DiscoSNP++ on $READS (this is the long step)"
echo "$READS" > fof.txt
D=$(ls discoRes_*coherent.vcf 2>/dev/null | head -1)
if [ -z "$D" ]; then
  /usr/bin/time -v run_discoSnp++.sh -r fof.txt -k 31 -c 3 -D 100 -P 3 -b 0 \
      -G "$REF" -T > disco.log 2>&1 || true
  D=$(ls discoRes_*coherent.vcf 2>/dev/null | head -1)
fi
[ -n "$D" ] || { echo "DiscoSNP++ produced no coherent VCF — see $OUT/disco.log" >&2; exit 1; }
grep -E "Maximum resident|Elapsed" disco.log || true
log "[1/4] done -- $D ($(grep -vc '^#' "$D") raw records)"

# ── 2. Convert to a comparable VCF (POS off-by-one corrected) ──────────────
log "[2/4] converting DiscoSNP++ output to comparable VCF..."
(echo '##fileformat=VCFv4.2'; echo "##contig=<ID=$CHROM,length=63025520>";
 echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';
 printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
 grep -v "^#" "$D" | awk -F'\t' -v OFS='\t' -v c="$CHROM" '$1==c{
    split($10,g,":"); gt=g[1]; gsub(/\|/,"/",gt);
    print $1,$2-1,".",$4,$5,30,"PASS",".","GT",gt}' | sort -k2,2n) > d_raw.vcf
log "[2/4] done -- $(grep -vc '^#' d_raw.vcf) chr20 records"

# ── 3. Truth: whole chr20, het-only, inside GIAB confident regions ─────────
log "[3/4] fetching GIAB truth for all of chr20..."
tabix -h "$TRUTH" "$CHROM" 2>/dev/null | awk -v OFS='\t' '
  /^##contig/{next} /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
  /^#/{print; next} {print}' > truth.vcf
awk -v c=$CHROM '$1==c{print c"\t"$2"\t"$3}' "$BED" > regions.bed
log "[3/4] done -- $(awk '{n+=$3-$2}END{print n+0}' regions.bed) confident bp on chr20"

snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
het()   { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1];
          if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
prep(){ (grep '^#' "$1"; grep -v '^#' "$1"|sort -k2,2n) > "$2.pre.vcf"
  bcftools norm -f "$REF" -m -any "$2.pre.vcf" 2>/dev/null | bcftools sort 2>/dev/null > "$2.n.vcf" \
    || cp "$2.pre.vcf" "$2.n.vcf"
  eval "$3" < "$2.n.vcf" > "$2.v.vcf"
  awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++' "$2.v.vcf" > "$2.s.vcf"
  bgzip -f -c "$2.s.vcf" > "$2.vcf.gz"
  tabix -f -p vcf "$2.vcf.gz"; }
prep truth.vcf t_snv 'het | snv'
prep truth.vcf t_ind 'het | indel'
prep d_raw.vcf d_snv snv
prep d_raw.vcf d_ind indel

# ── 4. Score with rtg vcfeval ──────────────────────────────────────────────
log "[4/4] scoring SNV+INDEL with rtg vcfeval..."
[ -d "$HOME/refs/chr20.sdf" ] && SDF="$HOME/refs/chr20.sdf" || { rm -rf sdf; rtg format -o sdf "$REF" >/dev/null 2>&1; SDF=sdf; }
score(){ rm -rf "e_$1"
  rtg vcfeval -b "$2" -c "$3" -t "$SDF" --squash-ploidy --bed-regions regions.bed \
      -o "e_$1" >/dev/null 2>&1 || true
  if [ -f "e_$1/summary.txt" ]; then
    awk -v n="$1" 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{
      printf "%-6s TP=%-7s FP=%-7s FN=%-7s  P=%.3f  R=%.3f  F1=%.3f\n",n,tp,fp,fn,p,r,f}' "e_$1/summary.txt"
  else echo "$1: no summary (see e_$1)"; fi; }

echo "======== $IND FULL chr20 — DiscoSNP++ — rtg vcfeval (GA4GH) ========"
score SNV   t_snv.vcf.gz d_snv.vcf.gz
score INDEL t_ind.vcf.gz d_ind.vcf.gz
echo "truth het-SNV: $(zcat t_snv.vcf.gz|grep -vc '^#')  truth het-indel: $(zcat t_ind.vcf.gz|grep -vc '^#')"
echo "disco SNV calls: $(zcat d_snv.vcf.gz|grep -vc '^#')  indel calls: $(zcat d_ind.vcf.gz|grep -vc '^#')"
echo "===================================================================="
log "Total elapsed: $(( $(date +%s) - T_START ))s"
