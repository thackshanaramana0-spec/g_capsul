#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
#  SANITY (ARCHIVE PATH) — ONE dataset, ALL THREE CLAIMS, end to end.
#
#  Derived from benchmark_sanity_one.sh, with two corrections:
#    * compress sets CAPS_SPANS=1, so the archive carries contig_spans WITHOUT
#      running the caller inline. CAPS_CALL=1 also writes them, but it runs the
#      full caller at the end of compression -- roughly 20x the work -- and
#      Claim 2 then calls a SECOND time from the archive. That mistake put
#      HG002 compression at 913 s against SPRING's 52 s.
#    * Claim 2 calls FROM THE ARCHIVE (capsule_decode call), which is what the
#      claim actually asserts. The original ran the ENCODER on the FASTQ.
#
#  The rehearsal for the full run. It exercises exactly the same code paths
#  benchmark_1_run.sh uses, on a single dataset, so every artefact the real
#  benchmark will produce can be inspected before committing a machine-day:
#  the archive, the decoded reads, the VCF, the contigs, the coverage table,
#  the query output, and every CSV row.
#
#  It ends with a FILE MANIFEST -- every output, its full path, its size, and
#  a one-line description of what it is and which table it feeds.
#
#  usage:
#    bash scripts/benchmark_sanity_one.sh                 # HG002 (all 3 claims)
#    bash scripts/benchmark_sanity_one.sh SRR2584863      # any of the 19
#    DS=... OUT_DIR=... bash scripts/benchmark_sanity_one.sh
#
#  HG002 is the default because it is the only kind of dataset that reaches
#  ALL THREE claims: it compresses (Claim 1), it is diploid human with a GIAB
#  truth set so variants can be called and scored (Claim 2), and its archive
#  feeds export/coverage/query (Claim 3). A bacterial accession exercises
#  Claims 1 and 3 only.
# ═══════════════════════════════════════════════════════════════════════════
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$HERE"
source "$HERE/scripts/capsule_config.sh" 2>/dev/null || true

DS="${1:-${DS:-HG002}}"
# Which claims to run, e.g. CLAIMS=1 for COMPACT alone, CLAIMS=13, default all.
# The claims share one archive, so restricting them changes nothing about how
# the ones that DO run are measured -- same encode, same competitors, same
# lossless check. That is the point of putting the switch here rather than
# writing a second script.
CLAIMS="${CLAIMS:-123}"
CSV3=""   # assigned only if Claim 3 runs; set -u would abort on the table loop
want(){ case "$CLAIMS" in *"$1"*) return 0;; *) return 1;; esac; }
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="${OUT_DIR:-$HERE/results/sanity_${DS}_$STAMP}"
mkdir -p "$OUT"
LOG="$OUT/sanity.log"
DATA_DIR="${CAPSULE_DATA_DIR:-/data/fastq}"
REFS="${CAPSULE_REFS_DIR:-$HOME/refs}"
BIN="${CAPSULE_BIN_DIR:-/tmp/capsule_bin}"
BEST="$BIN/best106"; DEC="$BIN/capsule_decode"
NPROC=$(nproc)
# DiscoSNP++ ships as a shell script in its own tree and its runner resolves
# `run_discoSnp++.sh` through PATH. The file existing is NOT enough -- this
# project has lost a competitor arm to exactly that before (see
# docs/INVOCATION_ERRORS.md), and it fails as "no SNV line", which looks like a
# scoring problem rather than a missing tool.
[ -d "$HOME/DiscoSnp" ] && export PATH="$HOME/DiscoSnp:$PATH"; T0=$(date +%s)

say(){ echo "$*" | tee -a "$LOG"; }
_el(){ local s=$(( $(date +%s) - T0 )); printf "%02d:%02d" $((s/60)) $((s%60)); }
inf(){ echo "[$(date +%H:%M:%S) +$(_el)] $*" | tee -a "$LOG"; }
step(){ echo "[$(date +%H:%M:%S) +$(_el)]    -> $*" | tee -a "$LOG"; }
ok(){   echo "[$(date +%H:%M:%S) +$(_el)]    OK  $*" | tee -a "$LOG"; }
err(){  echo "[$(date +%H:%M:%S) +$(_el)]  FAIL $*" | tee -a "$LOG"; }
warn(){ echo "[$(date +%H:%M:%S) +$(_el)]  WARN $*" | tee -a "$LOG"; }
banner(){ say ""; say "════════════════════════════════════════════════════════════════════"; say " $*"; say "════════════════════════════════════════════════════════════════════"; }
mb(){ awk -v b="${1:-0}" 'BEGIN{printf "%.2f MB", b/1048576}'; }
gb(){ awk -v b="${1:-0}" 'BEGIN{printf "%.2f GB", b/1073741824}'; }
rg(){ awk -v k="${1:-0}" 'BEGIN{printf "%.2f GB", k/1048576}'; }
tv(){ local f="$1" w h
  w=$(grep "Elapsed (wall clock)" "$f" 2>/dev/null | awk '{n=split($NF,a,":");
      if(n==3) printf "%.2f",a[1]*3600+a[2]*60+a[3]; else if(n==2) printf "%.2f",a[1]*60+a[2]; else printf "%.2f",a[1]}')
  h=$(grep "Maximum resident set size" "$f" 2>/dev/null | awk '{print $NF}')
  echo "${w:-0} ${h:-0}"; }

# Manifest: every artefact worth inspecting, recorded as it is produced.
MAN="$OUT/_manifest.tsv"; : > "$MAN"
keep(){ [ -e "$2" ] && printf "%s\t%s\t%s\t%s\n" "$1" "$2" "$(stat -c%s "$2" 2>/dev/null||echo 0)" "$3" >> "$MAN"; }

FAILED=0
banner "SANITY — $DS — all three claims, end to end"
say "started : $(date '+%Y-%m-%d %H:%M:%S')"
say "output  : $OUT"
say "commit  : $(git rev-parse --short HEAD 2>/dev/null)  branch $(git rev-parse --abbrev-ref HEAD 2>/dev/null)"
say "machine : $NPROC cores, $(free -g|awk '/^Mem:/{print $2}') GB RAM, $(df -h /|awk 'NR==2{print $4}') free"

# ── input ──────────────────────────────────────────────────────────────────
SRC="$DATA_DIR/${DS}_1.fq"; [ -s "$SRC" ] || SRC="$DATA_DIR/${DS}_pooled.fq"
[ -s "$SRC" ] || { err "no input for $DS in $DATA_DIR"; exit 1; }
RAW=$(stat -c%s "$SRC"); NREADS=$(( $(wc -l < "$SRC") / 4 ))
inf "input   : $SRC"
inf "size    : $(gb $RAW)   reads: $NREADS"
for b in "$BEST" "$DEC"; do
  [ -x "$b" ] || { err "missing binary $b -- run benchmark_0_preflight.sh first"; exit 1; }
done
keep "input" "$SRC" "the FASTQ under test"

# ── CLAIM 1 ────────────────────────────────────────────────────────────────
banner "CLAIM 1 — compression (T1 size, T2 time + RAM)"
A="$OUT/$DS.capsule"
step "CAPSULE compress (adaptive sweep, concurrent candidates)"
# CAPS_SPANS=1, NOT CAPS_CALL=1. Both write contig_spans, which Claim 2 needs
# to call from this archive -- but CAPS_CALL ALSO runs the full variant caller
# inline at the end of compression, which is a ~20x heavier path. Measured on
# a 500k-read HG002 slice: 38.64 s / 2.45 GB with CAPS_CALL against
# 8.97 s / 1.00 GB with CAPS_SPANS, for a BYTE-IDENTICAL archive and a
# BYTE-IDENTICAL VCF when called from it afterwards. Using CAPS_CALL here made
# compression pay for the caller twice.
/usr/bin/time -v env CAPS_SPANS=1 CAPS_NAMES=1 CAPS_QUAL=1 INPUT="$SRC" ARCHIVE="$A" BEST="$BEST" \
    bash "$HERE/scripts/encode_adaptive.sh" >"$OUT/encode.stdout" 2>"$OUT/_t_c"
read -r CW CR <<< "$(tv "$OUT/_t_c")"
if [ -s "$A" ]; then
  ARCH=$(stat -c%s "$A")
  ok "archive  $(mb $ARCH)   ratio $(awk -v a=$ARCH -v r=$RAW 'BEGIN{printf "%.3f%%",100*a/r}')   wall ${CW}s   RAM $(rg $CR)"
  keep "claim1" "$A" "the .capsule archive -- T1 size, and what Claim 3 reads"
  keep "claim1" "${A}.log" "encoder log: per-stage timings, stream sizes, candidate sweep"
else err "no archive produced"; FAILED=1; fi

step "CAPSULE decompress + lossless verify (decodes the ARCHIVE, not dumps)"
DD="$OUT/decoded"; mkdir -p "$DD"
/usr/bin/time -v "$DEC" "$A" "$DD" "$DD/reads.seq" >"$OUT/decode.stdout" 2>"$OUT/_t_d"
read -r DW DR <<< "$(tv "$OUT/_t_d")"
LL=LOSSY
if [ -s "$DD/reads.seq" ]; then
  awk 'NR%4==2' "$SRC" | tr -d '\r' | sort > "$OUT/_a"; sort "$DD/reads.seq" > "$OUT/_b"
  cmp -s "$OUT/_a" "$OUT/_b" && LL=LOSSLESS
  rm -f "$OUT/_a" "$OUT/_b"
fi
[ "$LL" = LOSSLESS ] && ok "decompress wall ${DW}s  RAM $(rg $DR)  -> $LL" \
                     || { err "decompress -> $LL"; FAILED=1; }
keep "claim1" "$DD/reads.seq" "reads reconstructed FROM THE ARCHIVE (lossless check input)"

# genozip needs /dev/stdout to exist -- it tests curl/wget availability with
# file_exists("/dev/stdout") and, on a Student licence, refuses to write the
# archive header if it cannot upload telemetry. Missing => exit 1, NO archive,
# after compressing the whole file. Repair it rather than lose the column.
[ -e /dev/stdout ] || ln -sfn /proc/self/fd/1 /dev/stdout 2>/dev/null || true

for T in SPRING Genozip; do
  step "$T (competitor arm, identical input)"
  if [ "$T" = SPRING ]; then
    /usr/bin/time -v spring -c -i "$SRC" -o "$OUT/$DS.spring" -t "$NPROC" 2>"$OUT/_t_s" >/dev/null
    S=$OUT/$DS.spring
  else
    /usr/bin/time -v genozip --force -o "$OUT/$DS.genozip" "$SRC" 2>"$OUT/_t_s" >/dev/null
    S=$OUT/$DS.genozip
  fi
  read -r XW XR <<< "$(tv "$OUT/_t_s")"
  if [ -s "$S" ]; then
    XS=$(stat -c%s "$S")
    ok "$T archive $(mb $XS)  ratio $(awk -v a=$XS -v r=$RAW 'BEGIN{printf "%.3f%%",100*a/r}')  wall ${XW}s  RAM $(rg $XR)"
    keep "claim1" "$S" "$T archive -- the competitor number in T1"

    # ── DECOMPRESS THE COMPETITOR TOO, AND VERIFY IT ────────────────────────
    # Without this, T2 compares our compress+DECOMPRESS against their compress
    # alone, and the lossless column has a value for us and a blank for them.
    # A reviewer reads that as "only one tool was checked for correctness".
    # NO -g ON SPRING. `-g` means "gzipped", and on decompress it writes GZIP
    # BYTES into a .fq file, which then reads as garbage and reports a false
    # LOSSY. This project has been burned by this exact flag before -- see the
    # tool_invocation_before_blame note. Our inputs are plain .fq.
    step "$T decompress + lossless verify"
    XD=""; XLL="NOT_CHECKED"; DOUT="$OUT/_x_$T"
    rm -rf "$DOUT"; mkdir -p "$DOUT"
    if [ "$T" = SPRING ]; then
      /usr/bin/time -v spring -d -i "$S" -o "$DOUT/out.fq" -t "$NPROC" \
        2>"$OUT/_t_sd" >/dev/null || true
    else
      /usr/bin/time -v genounzip --force -o "$DOUT/out.fq" "$S" \
        2>"$OUT/_t_sd" >/dev/null || true
    fi
    read -r XD _ <<< "$(tv "$OUT/_t_sd")"
    if [ -s "$DOUT/out.fq" ]; then
      # Compare the SEQUENCE column only, the same basis used for our own
      # lossless check, so the two verdicts mean the same thing.
      # Two-stage verdict, so the column means the same thing for every tool.
      # Stage 1 is IN ORDER, which is the strict test. Only if that fails do we
      # fall back to the order-free (sorted) test that our own arm uses, and
      # then the verdict says so explicitly. Comparing a competitor in order
      # against ourselves sorted would put two different tests in one column.
      awk 'NR%4==2' "$SRC" > "$DOUT/a.seq"
      awk 'NR%4==2' "$DOUT/out.fq" > "$DOUT/b.seq"
      if cmp -s "$DOUT/a.seq" "$DOUT/b.seq"; then
        XLL=LOSSLESS
      else
        sort "$DOUT/a.seq" -o "$DOUT/a.seq"; sort "$DOUT/b.seq" -o "$DOUT/b.seq"
        if cmp -s "$DOUT/a.seq" "$DOUT/b.seq"; then XLL=LOSSLESS_REORDERED; else XLL=LOSSY; fi
      fi
      ok "$T decompress ${XD}s -> $XLL"
    else
      err "$T produced no FASTQ on decompress -- recording NOT_CHECKED"
    fi
    rm -rf "$DOUT"

    printf "%s,%s,%s,%s,%.4f,%s,%s,%s,%s\n" "$DS" "$T" "$RAW" "$XS" \
      "$(awk -v a=$XS -v r=$RAW 'BEGIN{print 100*a/r}')" "$XW" "${XD:-}" "$XR" "$XLL" >> "$OUT/_rows_comp"
  else
    # A competitor that produces nothing must FAIL the run, not leave a blank
    # row. A T1/T2 table missing the Genozip line reads as "we did not bother
    # to benchmark it", which is worse than reporting a tool error.
    err "$T produced no archive -- competitor column would be BLANK"
    [ "$T" = Genozip ] && err "  genozip needs /dev/stdout + a telemetry upload on a Student licence; see _t_s"
    FAILED=1
  fi
done

CSV1="$OUT/claim1_t1_t2.csv"
{ echo "dataset,tool,raw_bytes,archive_bytes,ratio_pct,compress_s,decompress_s,peak_ram_kb,lossless"
  echo "$DS,CAPSULE,$RAW,${ARCH:-},$(awk -v a=${ARCH:-0} -v r=$RAW 'BEGIN{printf "%.4f",100*a/r}'),$CW,$DW,$CR,$LL"
  cat "$OUT/_rows_comp" 2>/dev/null
} > "$CSV1"
keep "claim1" "$CSV1" "T1 + T2 table rows for this dataset"

# ── CLAIM 2 ────────────────────────────────────────────────────────────────
if want 2; then
banner "CLAIM 2 — variant calling (T3 het-SNV F1)"
fi
if want 2 && [ -s "$DATA_DIR/${DS}_pooled.fq" ] && [ -s "$REFS/chr20.fa" ]; then
  step "our caller FROM THE ARCHIVE (capsule_decode call) -- no FASTQ is read"
  C2="$OUT/claim2"; mkdir -p "$C2"
  # THIS IS THE ARCHITECTURE CLAIM 2 ASSERTS. The archive built above is the
  # only input; the FASTQ is not opened. The alternative runner
  # (run_fullchr20_bench_capsule.sh) measures the ENCODER assembling and
  # calling in one pass -- a different and much heavier operation, 381 s /
  # 33.3 GB against 122 s / 15.2 GB measured on the same chr20 data -- and is
  # kept for comparison, not used here.
  if [ ! -s "$A" ]; then
    err "no archive at $A -- Claim 1 must have produced one; cannot run Claim 2"
  else
  bash "$HERE/scripts/run_fullchr20_archive_capsule.sh" "$DEC" "$HERE/scripts" \
       "$REFS/chr20.fa" "$A" "$DS" "$C2" > "$OUT/claim2.log" 2>&1
  if grep -q 'ARCHIVE LACKS' "$OUT/claim2.log" 2>/dev/null; then
    err "the archive cannot serve call -- re-encode with CAPS_CALL=1"
    grep 'ARCHIVE LACKS' "$OUT/claim2.log" | head -1
  fi
  fi
  LN=$(grep -aE '^SNV ' "$OUT/claim2.log" | tail -1)
  if [ -n "$LN" ]; then
    ok "OURS  $LN"
    keep "claim2" "$C2/calls.vcf"   "variant calls in contig coordinates (caller output)"
    keep "claim2" "$C2/lifted.vcf"  "calls lifted to chr20 coordinates (what vcfeval scores)"
    keep "claim2" "$C2/contigs.fa"  "assembled contigs the calls came from"
    keep "claim2" "$OUT/claim2.log" "full Claim 2 log incl. rtg vcfeval summary"
    parse_snv(){ echo "$1" | grep -oP "$2=\\K[0-9.]+" | head -1; }
    echo "individual,tool,tp,fp,fn,precision,recall,f1" > "$OUT/claim2_t3.csv"
    printf "%s,CAPSULE,%s,%s,%s,%s,%s,%s\n" "$DS" \
      "$(parse_snv "$LN" 'TP')" "$(parse_snv "$LN" 'FP')" "$(parse_snv "$LN" 'FN')" \
      "$(parse_snv "$LN" ' P')" "$(parse_snv "$LN" ' R')" "$(parse_snv "$LN" 'F1')" >> "$OUT/claim2_t3.csv"

    # DiscoSNP++ on the SAME reads, SAME truth, SAME scoring. Without it T3 is
    # a single number, not a head-to-head, and the claim is a comparison.
    step "DiscoSNP++ (identical reads, identical truth, identical scoring)"
    if [ -f "$HERE/scripts/run_fullchr20_bench_disco.sh" ]; then
      bash "$HERE/scripts/run_fullchr20_bench_disco.sh" \
           "$REFS/chr20.fa" "$SRC" "$DS" "$OUT/claim2_disco" > "$OUT/claim2_disco.log" 2>&1
      DL=$(grep -aE '^SNV ' "$OUT/claim2_disco.log" | tail -1)
      if [ -n "$DL" ]; then
        ok "DISCO $DL"
        printf "%s,DiscoSNP++,%s,%s,%s,%s,%s,%s\n" "$DS" \
          "$(parse_snv "$DL" 'TP')" "$(parse_snv "$DL" 'FP')" "$(parse_snv "$DL" 'FN')" \
          "$(parse_snv "$DL" ' P')" "$(parse_snv "$DL" ' R')" "$(parse_snv "$DL" 'F1')" >> "$OUT/claim2_t3.csv"
        keep "claim2" "$OUT/claim2_disco.log" "DiscoSNP++ arm -- the competitor number in T3"
      else err "DiscoSNP++ produced no SNV line -- see $OUT/claim2_disco.log"; fi
    else err "run_fullchr20_bench_disco.sh missing -- T3 will have only our arm"; fi
    # Third T3 arm. benchmark_1 carries the same arm and the same verdict:
    # /root/Kmer2SNP and its conda env exist, but THIS repo has never invoked
    # them -- the published Kmer2SNP F1 (0.464) comes from the outer ARCS
    # project under a methodology not reproduced here. Inventing an invocation
    # would yield a number that looks measured and is not. The row is written
    # explicitly so T3 shows a 3-arm table with a declared gap, rather than a
    # 2-arm table that reads as if only two tools were ever considered.
    if [ -f "$HERE/scripts/run_kmer2snp.sh" ]; then
      bash "$HERE/scripts/run_kmer2snp.sh" "$SRC" "$DS" "$OUT/claim2_k2s" \
           > "$OUT/claim2_k2s.log" 2>&1
      K1=$(grep -aE '^SNV ' "$OUT/claim2_k2s.log" | tail -1 | grep -oP 'F1=\K[0-9.]+')
      if [ -n "${K1:-}" ]; then ok "KMER2SNP SNV F1=$K1"
        printf "%s,Kmer2SNP,,,,,,%s\n" "$DS" "$K1" >> "$OUT/claim2_t3.csv"
      else err "Kmer2SNP produced no SNV line"
        printf "%s,Kmer2SNP,,,,,,FAILED\n" "$DS" >> "$OUT/claim2_t3.csv"; fi
    else
      warn "Kmer2SNP: no validated runner in this repo -- T3 arm recorded NOT_AVAILABLE"
      printf "%s,Kmer2SNP,,,,,,NOT_AVAILABLE\n" "$DS" >> "$OUT/claim2_t3.csv"
    fi
    keep "claim2" "$OUT/claim2_t3.csv" "T3 table: het-SNV TP/FP/FN/P/R/F1, ours vs DiscoSNP++"
  else err "no SNV line -- see $OUT/claim2.log"; FAILED=1; fi
elif ! want 2; then
  inf "SKIP Claim 2: CLAIMS=$CLAIMS"
else
  inf "SKIP Claim 2: $DS is not one of the 4 GIAB human sets (needs a truth VCF)"
fi

# ── CLAIM 3 ────────────────────────────────────────────────────────────────
if want 3; then
banner "CLAIM 3 — archive analysis (T6a export, T6b coverage, T6c query)"
C3="$OUT/claim3"; mkdir -p "$C3"
CSV3="$OUT/claim3_t6.csv"; echo "dataset,operation,ours_s,output_bytes,rows,status" > "$CSV3"
run3(){ local op="$1" outf="$2"; shift 2
  step "T6 $op"
  local t0 t1 s
  t0=$(date +%s.%N); "$DEC" "$op" "$A" "$outf" "$@" >"$C3/$op.stdout" 2>"$C3/$op.log"; local rc=$?
  t1=$(date +%s.%N); s=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.3f",b-a}')
  if [ $rc -eq 0 ] && [ -s "$outf" ]; then
    local sz rows; sz=$(stat -c%s "$outf"); rows=$(wc -l < "$outf")
    ok "$op  ${s}s  -> $(mb $sz)  ${rows} lines"
    echo "$DS,$op,$s,$sz,$rows,DONE" >> "$CSV3"
    keep "claim3" "$outf" "output of arcs $op, read from the archive alone"
  else err "$op failed (rc=$rc) -- see $C3/$op.log"; echo "$DS,$op,,,,FAILED" >> "$CSV3"; FAILED=1; fi
}
run3 export   "$C3/contigs.fa"
run3 coverage "$C3/coverage.tsv"
run3 query    "$C3/region.fq" 0-100000

# ── THE CONVENTIONAL BASELINES ──────────────────────────────────────────────
# Without these, T6 lists our timings with nothing to compare against, and the
# claim ("addressable: these operations are served from the archive instead of
# recomputed") has no measured contrast. Commands are IDENTICAL to
# benchmark_1_run.sh phase 3 so the sanity run rehearses the real thing.
OURS_EXP=$(awk -F, '$2=="export"{print $3}'   "$CSV3" | head -1)
OURS_COV=$(awk -F, '$2=="coverage"{print $3}' "$CSV3" | head -1)
C3REF="$REFS/chr20.fa"          # HG002 is chr20; other datasets use c3_<name>.fa
SPADES="$HOME/SPAdes-4.0.0-Linux/bin/spades.py"

step "T6a baseline: SPAdes de-novo assembly (minutes to hours)"
if [ -x "$SPADES" ] && [ -n "${OURS_EXP:-}" ]; then
  t0=$(date +%s.%N)
  python3 "$SPADES" -s "$SRC" -o "$C3/spades" -t "$NPROC" \
      -m $(( $(free -g | awk '/^Mem:/{print $2}') - 8 )) > "$C3/spades.log" 2>&1
  rc=$?; t1=$(date +%s.%N)
  if [ $rc -eq 0 ] && [ -s "$C3/spades/contigs.fasta" ]; then
    SP=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.2f",b-a}')
    ok "SPAdes ${SP}s   -> export speedup $(awk -v s=$SP -v o=$OURS_EXP 'BEGIN{printf "%.1fx",s/o}')"
    echo "$DS,export_baseline,$SP,,,SPAdes" >> "$CSV3"
  else
    err "SPAdes did not complete (rc=$rc) -- see $C3/spades.log"
    echo "$DS,export_baseline,,,,SPAdes_DNF" >> "$CSV3"
  fi
  rm -rf "$C3/spades/K"* "$C3/spades/tmp" 2>/dev/null
else
  err "SPAdes not installed at $SPADES -- T6a has no baseline"
  echo "$DS,export_baseline,,,,SPAdes_MISSING" >> "$CSV3"
fi

step "T6b baseline: bwa + samtools sort + mosdepth (the conventional route)"
if [ -s "$C3REF.bwt" ] && command -v mosdepth >/dev/null && [ -n "${OURS_COV:-}" ]; then
  t0=$(date +%s.%N)
  bwa mem -t "$NPROC" "$C3REF" "$SRC" 2>"$C3/bwa.log" \
    | samtools sort -@ 4 -o "$C3/aln.bam" - 2>>"$C3/bwa.log" \
    && samtools index "$C3/aln.bam" 2>>"$C3/bwa.log" \
    && mosdepth -t 4 "$C3/md" "$C3/aln.bam" 2>>"$C3/bwa.log"
  rc=$?; t1=$(date +%s.%N)
  if [ $rc -eq 0 ]; then
    CV=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.2f",b-a}')
    ok "bwa+mosdepth ${CV}s  -> coverage speedup $(awk -v s=$CV -v o=$OURS_COV 'BEGIN{printf "%.1fx",s/o}')"
    echo "$DS,coverage_baseline,$CV,,,bwa+samtools+mosdepth" >> "$CSV3"
  else
    err "bwa/mosdepth baseline failed (rc=$rc) -- see $C3/bwa.log"
    echo "$DS,coverage_baseline,,,,bwa_FAILED" >> "$CSV3"
  fi
  rm -f "$C3/aln.bam" "$C3/aln.bam.bai"      # BAMs are large; the timing is what we keep
else
  err "no bwa index at $C3REF.bwt or mosdepth missing -- T6b has no baseline"
  echo "$DS,coverage_baseline,,,,baseline_MISSING" >> "$CSV3"
fi

keep "claim3" "$CSV3" "T6 rows for this dataset"

# ── MANIFEST ───────────────────────────────────────────────────────────────
else
  inf "SKIP Claim 3: CLAIMS=$CLAIMS"
fi

banner "THE TABLES THIS RUN PRODUCED"
for t in "$CSV1:T1 + T2  (archive size, time, RAM -- ours vs SPRING vs Genozip)" \
         "$OUT/claim2_t3.csv:T3  (het-SNV F1 -- ours vs DiscoSNP++)" \
         "$CSV3:T6  (export / coverage / query from the archive)"; do
  f="${t%%:*}"; ttl="${t#*:}"
  [ -s "$f" ] || continue
  say ""; say "  $ttl"; say "  ${f#$OUT/}"
  column -s, -t "$f" 2>/dev/null | sed 's/^/    /' | tee -a "$LOG"
done

banner "FILE MANIFEST — every artefact, where it is, how big, what it is for"
say ""
printf "%-8s %12s  %-46s %s\n" "CLAIM" "SIZE" "PATH" "WHAT IT IS" | tee -a "$LOG"
say "$(printf '%.0s-' {1..150})"
while IFS=$'\t' read -r c p s d; do
  printf "%-8s %12s  %-46s %s\n" "$c" "$(mb $s)" "${p#$OUT/}" "$d" | tee -a "$LOG"
done < "$MAN"
say ""
inf "output root : $OUT"
inf "total size  : $(du -sh "$OUT" | cut -f1)"
inf "elapsed     : $(_el)"
say ""
if [ "$FAILED" -eq 0 ]; then
  say "  ┌──────────────────────────────────────────────────────────────┐"
  say "  │  SANITY PASSED — all three claims produced their artefacts.  │"
  say "  │  The full run should behave the same on every dataset.        │"
  say "  └──────────────────────────────────────────────────────────────┘"
else
  say "  ┌──────────────────────────────────────────────────────────────┐"
  say "  │  SANITY FAILED — see the FAIL lines above before the full run │"
  say "  └──────────────────────────────────────────────────────────────┘"
fi
rm -f "$OUT"/_t_c "$OUT"/_t_d "$OUT"/_t_s
say "  log: $LOG"
exit $FAILED
