#!/usr/bin/env bash
# FULL-SCALE Claim 2 benchmark: whole chr20 at 30x, not a 400 kb window.
#
# WHY THIS EXISTS. Every Claim 2 number in the project so far
# (docs/CLAIM2_FINAL_VERDICT.md, CLAIM2_TABLES_AND_INDEL_SCAN.md) was measured
# on a ~75k-read chr20 window, while every Claim 1 number is full-scale. That
# asymmetry is the largest single exposure in the paper: a reviewer can ask
# "do the window numbers hold at 12.6M reads?" and until this script runs, the
# honest answer is "unknown". A prior attempt to run full scale was killed at
# 78 GB RSS; the caller RAM work (ridx off by default, chunked k-mer counting,
# contig-batched pileup, flat sorted kidx) is what makes it feasible now.
#
# METHODOLOGY IS DELIBERATELY IDENTICAL to run_window_bench_capsule.sh --
# same caller invocation, same frozen parameters, same lift, the same
# normalise-then-classify order (the 2026-09-03 multi-allelic scoring fix),
# the same rtg vcfeval call. ONLY the input scale and the region change. If
# any step here differs from the window script, the two are not comparable and
# the comparison is the entire point.
#
#   usage: run_fullchr20_archive_capsule.sh <decoder_exe> <scripts_dir> <ref.fa> \
#                                           <in.capsule> [individual] [outdir]
#
# ARCHIVE PATH VARIANT. Identical to run_fullchr20_bench_capsule.sh in every
# respect -- same lift, same normalise-then-classify order, same rtg vcfeval --
# EXCEPT that the caller reads the ARCHIVE instead of the FASTQ:
#
#     capsule_decode call <in.capsule> <out.vcf> <workdir>
#
# That is the architecture Claim 2 actually asserts: variants served from the
# compressed archive, with no FASTQ present. The FASTQ runner measures the
# encoder doing assembly and calling in one pass, which is a different (and
# much heavier) operation -- 381 s / 33.3 GB against 122 s / 15.2 GB measured
# on the same chr20 data.
set -e
T_START=$(date +%s)
log() { echo "[full] $(date '+%H:%M:%S') $*"; }

CAPS="$1"; SC="$2"; REF="$3"; READS="$4"; IND="${5:-HG002}"
OUT="${6:-/data/caps_full/${IND}_chr20}"
CHROM=20

TRUTH="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
[ -s "$BED" ] || BED="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.bed"
[ -s "$READS" ] || { echo "no archive: $READS" >&2; exit 1; }

mkdir -p "$OUT"; cd "$OUT"
NR_=$(( $(wc -l < "$READS") / 4 ))
echo "=== $IND FULL chr20 — $NR_ reads ==="

# ── 1. CAPSULE reference-free call, full scale ─────────────────────────────
# RSS is sampled every 10 s into rss.log rather than relying on /usr/bin/time
# alone, because the shape of the curve (which phase peaks) is what the RAM
# work needs to be judged on, not just the maximum.
log "[1/5] running CAPSULE caller FROM THE ARCHIVE — no FASTQ is read"
# CAPS_DUMP_CONTIGS is what the lift needs to map contig-space calls to the
# genome; it is set by the harness, not by production use, and it is the reason
# the caller emits contigs.tsv at all here.
export CAPS_DUMP_CONTIGS="$OUT/contigs.tsv"
if [ ! -s contigs.tsv ]; then
  # CAPS_CALL_INDELS=1 selects the full SNV+indel caller. Without it the
  # decoder takes the graph-only path, which is a different configuration and
  # must not be mixed into one measurement.
  /usr/bin/time -v env CAPS_CALL_INDELS=1 "$CAPS" call "$READS" "$OUT/calls.vcf" "$OUT/callwk" \
      > /dev/null 2> capsule.log || true
  rm -rf "$OUT/callwk"
else
  log "  contigs.tsv already present, reusing"
fi
if grep -q "ARCHIVE LACKS" capsule.log 2>/dev/null; then
  echo "FATAL: this archive cannot serve the call operation. Re-compress with CAPS_CALL=1" >&2
  grep "ARCHIVE LACKS" capsule.log >&2
  exit 1
fi
grep -E "CAPS-CALL|CAPS-CALL-TIMING|HETSCAN|Maximum resident|Elapsed" capsule.log || true
[ -s contigs.tsv ] || { echo "no contigs dumped — see $OUT/capsule.log" >&2; exit 1; }
cp contigs.tsv contigs.fa
log "[1/5] done -- $(grep -c '^>' contigs.fa) contigs"

# ── 2. Place contigs on the reference (EVALUATION ONLY — not part of calling) ─
log "[2/5] aligning contigs to $REF with bwa mem (eval-only coordinate lift)..."
[ -s c2r.sam ] || bwa mem -t "$(nproc)" "$REF" contigs.fa 2>/dev/null > c2r.sam
log "[2/5] done -- $(du -h c2r.sam | cut -f1) c2r.sam"

# ── 3. Lift contig-coordinate calls to genome coordinates ──────────────────
log "[3/5] lifting contig-coordinate calls to genome coordinates..."
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam "$REF" $CHROM lifted.vcf contigs.fa
log "[3/5] done"

# ── 4. Truth: whole chr20, het-only, inside GIAB confident regions ─────────
log "[4/5] fetching GIAB truth for all of chr20, restricting to confident regions..."
tabix -h "$TRUTH" "$CHROM" 2>/dev/null | awk -v OFS='\t' '
  /^##contig/{next} /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
  /^#/{print; next} {print}' > truth.vcf
awk -v c=$CHROM '$1==c{print c"\t"$2"\t"$3}' "$BED" > regions.bed
log "[4/5] done -- $(awk '{n+=$3-$2}END{print n+0}' regions.bed) confident bp on chr20"

# Classification runs AFTER `bcftools norm -m -any`, so $5 is a single allele
# and a plain length test is correct. See run_window_bench_capsule.sh for the
# full account of the multi-allelic misfiling bug this ordering fixes.
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
prep truth.vcf  t_snv 'het | snv'
prep truth.vcf  t_ind 'het | indel'
prep lifted.vcf c_snv snv
prep lifted.vcf c_ind indel

# ── 5. Score with rtg vcfeval ──────────────────────────────────────────────
log "[5/5] scoring SNV+INDEL with rtg vcfeval..."
[ -d "$HOME/refs/chr20.sdf" ] && SDF="$HOME/refs/chr20.sdf" || { rm -rf sdf; rtg format -o sdf "$REF" >/dev/null 2>&1; SDF=sdf; }
score(){ rm -rf "e_$1"
  rtg vcfeval -b "$2" -c "$3" -t "$SDF" --squash-ploidy --bed-regions regions.bed \
      -o "e_$1" >/dev/null 2>&1 || true
  if [ -f "e_$1/summary.txt" ]; then
    awk -v n="$1" 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{
      printf "%-6s TP=%-7s FP=%-7s FN=%-7s  P=%.3f  R=%.3f  F1=%.3f\n",n,tp,fp,fn,p,r,f}' "e_$1/summary.txt"
  else echo "$1: no summary (see e_$1)"; fi; }

echo "======== $IND FULL chr20 — CAPSULE reference-free — rtg vcfeval (GA4GH) ========"
score SNV   t_snv.vcf.gz c_snv.vcf.gz
score INDEL t_ind.vcf.gz c_ind.vcf.gz
echo "truth het-SNV: $(zcat t_snv.vcf.gz|grep -vc '^#')  truth het-indel: $(zcat t_ind.vcf.gz|grep -vc '^#')"
echo "capsule SNV calls: $(zcat c_snv.vcf.gz|grep -vc '^#')  indel calls: $(zcat c_ind.vcf.gz|grep -vc '^#')"
echo "================================================================================"
log "Total elapsed: $(( $(date +%s) - T_START ))s"
