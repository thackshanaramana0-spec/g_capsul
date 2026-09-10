#!/usr/bin/env bash
# CAPSULE version of the outer ARCS project's run_giab_indel.sh: REAL-GIAB
# reference-free SNV+indel benchmark, HG002 chr20:2.0-2.4M, ~30x, scored by
# rtg vcfeval against GIAB v4.2.1 truth restricted to het + the confident BED.
# Only the caller invocation differs from the ARCS original.
#   usage: run_giab_indel_capsule.sh <capsule_exe> <scripts_dir> <reads.fq> <ref.fa>
set -e
T_START=$(date +%s)
log() { echo "[giab] $(date '+%H:%M:%S') $*"; }
export PATH=~/miniconda3/bin:$PATH
CAPS="$1"; SC="$2"; READS="$3"; REF="$4"
CHROM=20; RLO=2000000; RHI=2400000; REGION="$CHROM:$RLO-$RHI"
WD=~/giab_indel_capsule; mkdir -p "$WD"; cd "$WD"
FTP=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37
VCFURL="$FTP/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BEDURL="$FTP/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"

# 1. CAPSULE reference-free call + dump contigs (fixed MAXMAP/MINOV, not the
#    encode_adaptive.sh sweep -- see docs/_removable/CLAIM2_BUILDER.md for why).
log "[1/6] running CAPSULE caller (CAPS_CALL=1) on $READS..."
export CAPS_CALL=1 CALL_VCF="$WD/calls.vcf" CAPS_DUMP_CONTIGS="$WD/contigs.tsv"
if [ -s calls.vcf ]; then
  log "  calls.vcf already present, reusing (cached from a prior run)"
else
  "$CAPS" "$READS" 3 16 16 22 16 16 1 24 64 1 > /dev/null 2> capsule_call.log
fi
grep -E "CAPS-CALL" capsule_call.log || true
cp contigs.tsv contigs.fa   # CAPS_DUMP_CONTIGS already emits FASTA (see run_window_bench_capsule.sh)
log "[1/6] done -- $(grep -c '^>' contigs.fa) contigs"

# 2. Place contigs (eval-only coordinate lift)
log "[2/6] aligning contigs to $REF with bwa mem (eval-only coordinate lift)..."
[ -s "$REF.bwt" ] || bwa index "$REF" 2>/dev/null
bwa mem -t4 "$REF" contigs.fa 2>/dev/null > c2r.sam
log "[2/6] done -- c2r.sam written"

# 3. Lift contig-coordinate calls to genome coords
log "[3/6] lifting contig-coordinate calls to genome coordinates..."
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam "$REF" $CHROM lifted.vcf contigs.fa
log "[3/6] done"

# 4. Truth for the region (remote tabix), header patched with contig length
log "[4/6] fetching GIAB truth for $REGION..."
tabix -h "$VCFURL" "$REGION" 2>/dev/null | awk -v OFS='\t' '
  /^##contig/{next} /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
  /^#/{print; next} {print}' > truth.vcf
log "[4/6] done"

# 5. Confident BED subset to the eval region (download once)
log "[5/6] fetching confident BED, restricting to region..."
[ -s giab_conf.bed ] || curl -s "$BEDURL" -o giab_conf.bed
awk -v c=$CHROM -v lo=$RLO -v hi=$RHI '$1==c && $3>lo && $2<hi{
  s=($2>lo?$2:lo); e=($3<hi?$3:hi); if(e>s) print c"\t"s"\t"e}' giab_conf.bed > regions.bed
log "[5/6] done -- $(awk '{n+=$3-$2}END{print n+0}' regions.bed) confident bp in region"

# 6. rtg SDF + split SNV/indel + score each restricted to the confident region
#
# BUG FIXED 2026-09-03 (same defect as run_window_bench_capsule.sh, see
# docs/HET_INDEL_FRESH_SCAN.md Finding 4): this used to classify BEFORE
# normalising, on the raw ALT field. A multi-allelic SNV (ALT="C,A", string
# length 3) was misfiled as an INDEL by `length($5)!=1`, then bcftools norm
# split it into SNV-shaped rows sitting in the indel call set as false
# positives -- a bug that could only ever penalise CAPSULE, since it is the
# only tool here that emits true multi-allelic records. Fixed by normalising
# FIRST (splitting multi-allelics + left-aligning), THEN classifying, THEN
# de-duplicating rows normalisation made identical -- applied identically to
# truth and calls.
log "[6/6] normalising, classifying, and scoring SNV+INDEL with rtg vcfeval..."
rm -rf sdf; rtg format -o sdf "$REF" > /dev/null 2>&1
snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
het()   { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1]; if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
dedup() { awk -F'\t' '/^#/{print;next} !seen[$1"\t"$2"\t"toupper($4)"\t"toupper($5)]++'; }
prep(){ (grep '^#' "$1"; grep -v '^#' "$1"|sort -k2,2n) > "$2.pre.vcf"
  bcftools norm -f "$REF" -m -any "$2.pre.vcf" 2>/dev/null | bcftools sort 2>/dev/null > "$2.n.vcf" \
    || cp "$2.pre.vcf" "$2.n.vcf"
  eval "$3" < "$2.n.vcf" | dedup > "$2.s.vcf"
  bgzip -cf "$2.s.vcf" > "$2.vcf.gz"; tabix -f -p vcf "$2.vcf.gz"; }
prep truth.vcf  t_snv   'het | snv'
prep truth.vcf  t_ind   'het | indel'
prep lifted.vcf c_snv   snv
prep lifted.vcf c_ind   indel
score(){ rm -rf "e_$1"; rtg vcfeval -b "$2" -c "$3" -t sdf --squash-ploidy \
    --bed-regions regions.bed -o "e_$1" >/dev/null 2>&1 || true
  if [ -f "e_$1/summary.txt" ]; then
    awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{printf "%-6s TP=%s FP=%s FN=%s  P=%.3f R=%.3f F1=%.3f\n","'$1'",tp,fp,fn,p,r,f}' "e_$1/summary.txt"
  else echo "$1: no summary"; fi; }
log "[6/6] done -- scoring"
echo "======== REAL-GIAB HG002 chr20:2.0-2.4M — CAPSULE — rtg vcfeval (gold-standard) ========"
score SNV   t_snv.vcf.gz c_snv.vcf.gz
score INDEL t_ind.vcf.gz c_ind.vcf.gz
echo "=========================================================================================="
echo "truth het-indel in region: $(zcat t_ind.vcf.gz|grep -vc '^#')   capsule indel calls: $(zcat c_ind.vcf.gz|grep -vc '^#')"
log "Total elapsed: $(( $(date +%s) - T_START ))s"
