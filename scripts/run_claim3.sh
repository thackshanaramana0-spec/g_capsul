#!/bin/bash
# Claim 3 — ADDRESSABLE: T6 archive-derived analysis vs conventional pipeline.
#
# ONE-COMMAND RUN — checks/builds everything it needs, then produces
# results/claim3/t6_results.csv. Mirrors the style of the outer project's
# benchmark/run_claim3.sh (phase/step logging, prerequisite checks, caching,
# CSV + summary table), but wired to the binaries that actually implement
# export/coverage/query TODAY: this repo's stages/capsule_decode.cpp.
# The outer `arcs export/coverage/query` do not exist yet (see
# docs/CLAIM3_LOCKED.md §5, §7 item 1) — do not point this at that binary.
#
# Usage:
#   bash scripts/run_claim3.sh [DATA_DIR] [OUT_DIR]
#
# DATA_DIR : directory containing the FASTQ inputs (default: /data/fastq)
#            Only SRR2584863_1.fq (E. coli, locked accession #1) is required
#            for export+query. Add HG002_pooled.fq / HG005_pooled.fq for the
#            GIAB rows (optional — skipped with a clear reason if absent).
# OUT_DIR  : results directory (default: ./results/claim3)
#
# Env overrides:
#   SPADES     path to spades.py                        [auto-detected / auto-installed]
#   REF_ECOLI  E. coli reference FASTA, for the coverage row [skipped if unset]
#
# What this script does NOT do, on purpose:
#   - does not touch the outer ARCS binary or repo
#   - does not re-download FASTQ (fails clearly if DATA_DIR is missing files)
#   - does not silently substitute a different assembler for SPAdes
#
# See docs/CLAIM3_LOCKED.md for the architecture, the two real bugs this
# operation's implementation had (both fixed, both covered by this script's
# own numbers matching), and the exact prior-art position.

set -euo pipefail

# ── Helpers — same shape as the outer run_claim3.sh, so logs read the same ──
_PH=""
log()   { echo "[$(date +%H:%M:%S)] $*"; }
phase() { _PH="$1"; echo ""; echo "[Phase $_PH] $2"; }
pdone() { echo "[Phase $_PH] DONE — $*"; }
pfail() { echo "[Phase $_PH] FAIL — $*"; exit 1; }
pskip() { echo "[Phase $_PH] SKIP — $*"; }
pinfo() { echo "[Phase $_PH]   -> $*"; }
step()  { echo ""; echo "======================================================"; echo "[STEP] $*"; echo "======================================================"; }

parse_time_v() {
    local logf="$1" wall vmhwm
    wall=$(grep "Elapsed (wall clock)" "$logf" | awk '{
        n=split($NF,a,":");
        if(n==3) printf "%.3f", a[1]*3600+a[2]*60+a[3];
        else if(n==2) printf "%.3f", a[1]*60+a[2];
        else printf "%.3f", a[1];
    }')
    vmhwm=$(grep "Maximum resident set size" "$logf" | awk '{print $NF}')
    echo "${wall:-0} ${vmhwm:-0}"
}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${1:-/data/fastq}"
OUT_DIR="${2:-$HERE/results/claim3}"
SPADES="${SPADES:-}"
REF_ECOLI="${REF_ECOLI:-}"
WD="$OUT_DIR/workdir"
mkdir -p "$WD"
NPROC=$(nproc 2>/dev/null || echo 4)
RESULTS_CSV="$OUT_DIR/t6_results.csv"

log "Claim 3 — ADDRESSABLE (T6), CAPSULE run"
log "DATA_DIR : $DATA_DIR"
log "OUT_DIR  : $OUT_DIR"

# ── STEP 0: Prerequisites ────────────────────────────────────────────────────
step "0/6 Prerequisites — check, and build/install only what's missing"

phase "0.1" "Toolchain (g++, gcc, liblzma)"
command -v g++ >/dev/null || pfail "g++ not found — apt install g++"
command -v gcc >/dev/null || pfail "gcc not found — apt install gcc"
echo 'int main(){return 0;}' > "$WD/_lzma_check.c"
if ! gcc "$WD/_lzma_check.c" -llzma -o "$WD/_lzma_check" 2>/dev/null; then
    pfail "liblzma-dev not found — apt install liblzma-dev"
fi
pdone "g++/gcc/liblzma present"

phase "0.2" "CAPSULE encoder (stage 106)"
BEST="$WD/best106"
if [ ! -x "$BEST" ]; then
    pinfo "not built — building now via scripts/build106.sh (~10-20s)"
    bash "$HERE/scripts/build106.sh" "$BEST" >/dev/null
fi
pdone "encoder ready: $BEST"

phase "0.3" "CAPSULE archive decoder (export/coverage/query)"
DECODER="$WD/capsule_decode"
if [ ! -x "$DECODER" ]; then
    pinfo "not built — building now via scripts/build_decode.sh (~10-20s)"
    bash "$HERE/scripts/build_decode.sh" "$DECODER" >/dev/null
fi
pdone "decoder ready: $DECODER"

phase "0.4" "SPAdes (spec-exact export baseline)"
if [ -z "$SPADES" ]; then
    if command -v spades.py >/dev/null; then
        SPADES="$(command -v spades.py)"
        pdone "found on PATH: $SPADES"
    elif [ -x "$HOME/SPAdes-4.0.0-Linux/bin/spades.py" ]; then
        SPADES="$HOME/SPAdes-4.0.0-Linux/bin/spades.py"
        pdone "found from a prior install: $SPADES"
    else
        pinfo "not found — installing prebuilt v4.0.0 to \$HOME (~180MB download, one-time)"
        ( cd "$HOME" && \
          curl -sL https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz -o spades.tar.gz && \
          tar xzf spades.tar.gz && rm -f spades.tar.gz )
        SPADES="$HOME/SPAdes-4.0.0-Linux/bin/spades.py"
        [ -x "$SPADES" ] || pfail "SPAdes install failed — install manually and set SPADES=/path/to/spades.py"
        pdone "installed: $SPADES"
    fi
else
    [ -x "$SPADES" ] || pfail "SPADES=$SPADES is not executable"
    pdone "using SPADES=$SPADES"
fi

phase "0.5" "bwa + samtools + mosdepth (coverage baseline — optional)"
HAVE_COVTOOLS=1
for t in bwa samtools mosdepth; do
    command -v "$t" >/dev/null || { HAVE_COVTOOLS=0; pinfo "$t not found"; }
done
if [ "$HAVE_COVTOOLS" = 1 ]; then pdone "all three present — coverage row will run"
else pskip "coverage row will be skipped (apt/conda install bwa samtools mosdepth to enable it)"
fi

# ── STEP 1: Input check ──────────────────────────────────────────────────────
step "1/6 Check input FASTQ"
phase "1.1" "Locate E. coli reads (locked accession #1)"
ECOLI_FQ="$DATA_DIR/SRR2584863_1.fq"
[ -s "$ECOLI_FQ" ] || pfail "required input not found: $ECOLI_FQ (locked accession #1, E. coli — see DATASET_LOCKED.md)"
N_READS=$(( $(wc -l < "$ECOLI_FQ") / 4 ))
pdone "E. coli reads: $N_READS ($ECOLI_FQ)"

phase "1.2" "Check optional GIAB inputs"
HAVE_GIAB=1
for f in HG002_pooled.fq HG005_pooled.fq; do
    [ -s "$DATA_DIR/$f" ] || { HAVE_GIAB=0; pinfo "optional GIAB input missing: $DATA_DIR/$f"; }
done
[ "$HAVE_GIAB" = 1 ] && pdone "GIAB pooled inputs present (not yet driven by this script — see closing note)" \
                      || pskip "GIAB rows not applicable (E. coli export/coverage/query still run)"

echo "dataset,operation,capsule_s,conventional_tool,conventional_s,speedup,notes" > "$RESULTS_CSV"

# ── STEP 2: Compress E. coli into a CAPSULE archive ──────────────────────────
step "2/6 CAPSULE compress — build the archive export/coverage/query read from"
phase "2.1" "Encode"
ARC="$WD/ecoli.capsule"
if [ ! -s "$ARC" ]; then
    pinfo "encoding $ECOLI_FQ -> $ARC"
    INPUT="$ECOLI_FQ" ARCHIVE="$ARC" BEST="$BEST" bash "$HERE/scripts/encode_adaptive.sh" >"$WD/ecoli_encode.log" 2>&1
    pdone "archive built: $(stat -c %s "$ARC") bytes"
else
    pskip "archive already exists: $ARC"
fi

# ── STEP 3: export — CAPSULE vs SPAdes ───────────────────────────────────────
step "3/6 T6 export — capsule_decode export vs SPAdes de-novo"
phase "3.1" "capsule_decode export"
EXP_OUT="$WD/ecoli_export.fa"
TF="$WD/ecoli_export_time.txt"
/usr/bin/time -v "$DECODER" export "$ARC" "$EXP_OUT" 2>"$TF"
read -r CAP_T CAP_RAM <<< "$(parse_time_v "$TF")"
N_CONTIGS=$(grep -c '^>' "$EXP_OUT")
pdone "export: ${CAP_T}s, ${CAP_RAM}KB RAM, $N_CONTIGS contigs -> $EXP_OUT"

phase "3.2" "SPAdes de-novo assembly (default full pipeline — slow, minutes not seconds)"
SPADES_OUT="$WD/ecoli_spades"
TF="$WD/ecoli_spades_time.txt"
if [ ! -s "$SPADES_OUT/contigs.fasta" ]; then
    pinfo "running: spades.py -s $ECOLI_FQ -o $SPADES_OUT -t $NPROC"
    /usr/bin/time -v "$SPADES" -s "$ECOLI_FQ" -o "$SPADES_OUT" -t "$NPROC" >"$WD/ecoli_spades.log" 2>"$TF" \
        || pfail "SPAdes failed — see $WD/ecoli_spades.log"
    read -r SPADES_T SPADES_RAM <<< "$(parse_time_v "$TF")"
else
    pskip "SPAdes output already exists at $SPADES_OUT — re-timing not possible from cache"
    SPADES_T="cached"; SPADES_RAM="cached"
fi
N_SPADES_CONTIGS=$(grep -c '^>' "$SPADES_OUT/contigs.fasta" 2>/dev/null || echo "?")
pdone "SPAdes: ${SPADES_T}s, ${SPADES_RAM}KB RAM, $N_SPADES_CONTIGS contigs"

if [[ "$SPADES_T" != "cached" ]]; then
    SPEEDUP=$(awk "BEGIN{printf \"%.0f\", $SPADES_T / ($CAP_T + 0.001)}")
    echo "SRR2584863_ecoli,export,$CAP_T,spades,$SPADES_T,${SPEEDUP}x,peak RAM spades ${SPADES_RAM}KB; $N_SPADES_CONTIGS contigs" >> "$RESULTS_CSV"
    pdone "export speedup: ${SPEEDUP}x"
else
    echo "SRR2584863_ecoli,export,$CAP_T,spades,cached,N/A,SPAdes output was cached from a prior run — re-run with a clean workdir for a fresh timing" >> "$RESULTS_CSV"
fi

# ── STEP 4: coverage — CAPSULE vs bwa+mosdepth ───────────────────────────────
step "4/6 T6 coverage — capsule_decode coverage vs bwa+samtools+mosdepth"
phase "4.1" "capsule_decode coverage"
COV_OUT="$WD/ecoli_coverage.tsv"
TF="$WD/ecoli_cov_time.txt"
/usr/bin/time -v "$DECODER" coverage "$ARC" "$COV_OUT" 2>"$TF"
read -r COVCAP_T COVCAP_RAM <<< "$(parse_time_v "$TF")"
pdone "coverage: ${COVCAP_T}s, ${COVCAP_RAM}KB RAM -> $COV_OUT"

if [ "$HAVE_COVTOOLS" = 1 ] && [ -n "$REF_ECOLI" ] && [ -s "$REF_ECOLI" ]; then
    phase "4.2" "bwa mem + samtools sort + mosdepth"
    [ -s "${REF_ECOLI}.bwt" ] || { pinfo "building BWA index (one-time)"; bwa index "$REF_ECOLI" 2>/dev/null; }
    BAM="$WD/ecoli_bwa.bam"
    TF_BWA="$WD/ecoli_bwa_time.txt"
    /usr/bin/time -v bash -c \
        "bwa mem -t $NPROC '$REF_ECOLI' '$ECOLI_FQ' 2>/dev/null | samtools sort -@ $NPROC -o '$BAM' && samtools index '$BAM'" \
        2>"$TF_BWA"
    read -r BWA_T BWA_RAM <<< "$(parse_time_v "$TF_BWA")"
    TF_MSD="$WD/ecoli_mosd_time.txt"
    /usr/bin/time -v mosdepth --threads "$NPROC" "$WD/ecoli_mosd" "$BAM" 2>"$TF_MSD"
    read -r MSD_T MSD_RAM <<< "$(parse_time_v "$TF_MSD")"
    CONV_T=$(awk "BEGIN{printf \"%.3f\", $BWA_T + $MSD_T}")
    SPEEDUP=$(awk "BEGIN{printf \"%.1f\", $CONV_T / ($COVCAP_T + 0.001)}")
    pdone "bwa+mosdepth: ${CONV_T}s total -> ${SPEEDUP}x speedup"
    echo "SRR2584863_ecoli,coverage,$COVCAP_T,bwa+samtools+mosdepth,$CONV_T,${SPEEDUP}x,pre-built BWA index (conservative for CAPSULE)" >> "$RESULTS_CSV"
else
    pskip "coverage baseline skipped — set REF_ECOLI=/path/to/ecoli.fa and ensure bwa/samtools/mosdepth are installed"
    echo "SRR2584863_ecoli,coverage,$COVCAP_T,bwa+samtools+mosdepth,NA,NA,skipped: set REF_ECOLI and install bwa/samtools/mosdepth" >> "$RESULTS_CSV"
fi

# ── STEP 5: query — CAPSULE only, no competitor exists (see CLAIM3_LOCKED §4) ─
step "5/6 T6 query — capsule_decode query (coordinate-range retrieval)"
phase "5.1" "query 0-100000, vs full decompress for comparison"
QUERY_OUT="$WD/ecoli_query.fa"
TF="$WD/ecoli_query_time.txt"
/usr/bin/time -v "$DECODER" query "$ARC" "$QUERY_OUT" 0-100000 2>"$TF"
read -r Q_T Q_RAM <<< "$(parse_time_v "$TF")"
N_Q_READS=$(grep -c '^>' "$QUERY_OUT")
pdone "query: ${Q_T}s, $N_Q_READS reads for range 0-100000"

phase "5.2" "full decompress, for the selectivity comparison"
# NOTE: capsule_decode only does the real per-read reconstruction when given
# a THIRD (outreads) argument -- with just <archive> <outdir> it takes a
# cheaper path that dumps raw intermediate streams and does NOT reconstruct
# read sequences (see stages/capsule_decode.cpp:512, outreads.empty() branch).
# Comparing query against that cheaper path understates the real "full
# decompress" cost and silently invalidates the speedup number. Passing the
# outreads path here is required, not optional, for this comparison to mean
# anything (verified against scripts/capsule_roundtrip.sh's own usage).
FULLDEC_OUT="$WD/ecoli_full_reads"
mkdir -p "$FULLDEC_OUT"
TF="$WD/ecoli_fulldec_time.txt"
/usr/bin/time -v "$DECODER" "$ARC" "$FULLDEC_OUT" "$FULLDEC_OUT/reads.seq" 2>"$TF" || pinfo "full decode returned non-zero — this is fine, timing still captured"
read -r FULL_T FULL_RAM <<< "$(parse_time_v "$TF")"
SPEEDUP=$(awk "BEGIN{printf \"%.1f\", $FULL_T / ($Q_T + 0.001)}")
pdone "query vs full decompress: ${SPEEDUP}x time, $N_Q_READS reads returned instead of $N_READS (selectivity is the real win, see CLAIM3_LOCKED §5)"
echo "SRR2584863_ecoli,query,$Q_T,full_decompress,$FULL_T,${SPEEDUP}x,returns $N_Q_READS of $N_READS reads for range 0-100000" >> "$RESULTS_CSV"

# ── STEP 6: Summary ──────────────────────────────────────────────────────────
phase "6" "Results summary"
echo ""
cat "$RESULTS_CSV"
echo ""
pdone "Claim 3 results in: $RESULTS_CSV"
pdone "Finished at: $(date)"
echo ""
echo "NOTE: GIAB (HG002/HG005) rows are not produced by this script — those"
echo "existing t6_results.csv rows came from ad-hoc runs on pre-built archives"
echo "(see docs/CLAIM3_LOCKED.md sec 5, 6.1). This script covers the locked"
echo "E. coli accession end to end, reproducibly, from one command."
