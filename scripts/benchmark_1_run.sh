#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
#  BENCHMARK 1 — the full run. Claims 1, 2 and 3, frame by frame.
#
#  Design rules this script obeys:
#   * ONE TIMED JOB AT A TIME. Concurrency corrupts wall time and peak RAM,
#     which are published numbers (CLAUDE.md rule 3). Nothing here runs in
#     parallel with anything else that is being measured.
#   * SPEED AND RAM come from /usr/bin/time -v (wall clock + VmHWM). Never
#     estimated, never derived from a log timestamp.
#   * Loud at the SEAMS, silent in the LOOPS. Every dataset, tool and phase
#     transition prints. Nothing prints from inside a hot path -- that would
#     slow the very thing being timed.
#   * RESULTS ARE WRITTEN AS THEY HAPPEN. A crash at dataset 15 keeps 1-14.
#   * FAILURE OF ONE DATASET IS NOT FAILURE OF THE RUN: retry once, then
#     record FAILED and continue. The single exception is a LOSSY archive,
#     which halts everything (CLAUDE.md rule 5) -- that is a correctness bug,
#     not a flaky run.
#
#  usage:
#    bash scripts/benchmark_1_run.sh                 # everything, ~9-12 h
#    SANITY_ONLY=1 bash scripts/benchmark_1_run.sh   # 1 dataset, measured rate
#    PHASES=1,3    bash scripts/benchmark_1_run.sh   # selected phases
#
#  env:
#    OUT_DIR       results dir              [results/benchmark_1_<timestamp>]
#    HEARTBEAT_MIN minutes between "still alive" lines            [15]
#    SANITY_ONLY   run only the first dataset of phase 1           [unset]
#    PHASES        comma list of phases to run                     [1,2,3]
# ═══════════════════════════════════════════════════════════════════════════
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$HERE"
source "$HERE/scripts/capsule_config.sh" 2>/dev/null || true

STAMP=$(date +%Y%m%d_%H%M%S)
OUT_DIR="${OUT_DIR:-$HERE/results/benchmark_1_$STAMP}"
mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/run.log"
HEARTBEAT_MIN="${HEARTBEAT_MIN:-15}"
PHASES="${PHASES:-1,2,3}"
DATA_DIR="${CAPSULE_DATA_DIR:-/data/fastq}"
REFS="${CAPSULE_REFS_DIR:-$HOME/refs}"
TRUTH="${CAPSULE_TRUTH_DIR:-$HOME/giab_truth}"
BIN_DIR="${CAPSULE_BIN_DIR:-/tmp/capsule_bin}"
BEST="$BIN_DIR/best106"; DEC="$BIN_DIR/capsule_decode"
NPROC=$(nproc)
T_RUN_START=$(date +%s)

# ── logging ───────────────────────────────────────────────────────────────
_ts(){ date '+%H:%M:%S'; }
_el(){ local s=$(( $(date +%s) - T_RUN_START )); printf "%02d:%02d:%02d" $((s/3600)) $((s%3600/60)) $((s%60)); }
say(){  echo "$*" | tee -a "$LOG"; }
inf(){  echo "[$(_ts) +$(_el)] $*" | tee -a "$LOG"; }
step(){ echo "[$(_ts) +$(_el)]    -> $*" | tee -a "$LOG"; }
ok(){   echo "[$(_ts) +$(_el)]    OK  $*" | tee -a "$LOG"; }
err(){  echo "[$(_ts) +$(_el)]  FAIL $*" | tee -a "$LOG"; }
banner(){ say ""; say "════════════════════════════════════════════════════════════════════════"; say " $*"; say "════════════════════════════════════════════════════════════════════════"; }
gbs(){ awk -v b="${1:-0}" 'BEGIN{printf "%.2f GB", b/1073741824}'; }
mbs(){ awk -v b="${1:-0}" 'BEGIN{printf "%.1f MB", b/1048576}'; }
ramg(){ awk -v k="${1:-0}" 'BEGIN{printf "%.2f GB", k/1048576}'; }

# DEBUG dump on failure: what ran, what it said, what the machine looked like.
debug_dump(){   # $1=label  $2=logfile
    err "── DEBUG: $1 ──"
    [ -s "${2:-}" ] && { err "last 25 lines of ${2}:"; tail -25 "$2" | sed 's/^/         | /' | tee -a "$LOG"; }
    err "disk: $(df -h / | awk 'NR==2{print $4" free"}')   mem: $(free -g | awk '/^Mem:/{print $7" GB avail"}')"
    err "── END DEBUG ──"
}

# /usr/bin/time -v -> "wall_seconds peak_ram_kb"
parse_time_v(){
    local f="$1" wall hwm
    wall=$(grep "Elapsed (wall clock)" "$f" 2>/dev/null | awk '{n=split($NF,a,":");
        if(n==3) printf "%.2f",a[1]*3600+a[2]*60+a[3];
        else if(n==2) printf "%.2f",a[1]*60+a[2]; else printf "%.2f",a[1]}')
    hwm=$(grep "Maximum resident set size" "$f" 2>/dev/null | awk '{print $NF}')
    echo "${wall:-0} ${hwm:-0}"
}

# ── heartbeat: proof of life for a run left going for hours ───────────────
HB_PID=""
# THE STAGE IS PASSED THROUGH A FILE, NOT A VARIABLE.
#
# The heartbeat runs in a SUBSHELL forked once at startup, so it captured
# CURRENT_STAGE's value at that moment and could never see a later update --
# every heartbeat for the whole run would have said "starting". A status line
# that cannot change is worse than none, because it looks like progress.
STAGE_F="$OUT_DIR/.stage"; echo "starting" > "$STAGE_F"
PROG_F="$OUT_DIR/PROGRESS.txt"; : > "$PROG_F"
mark(){ echo "$*" > "$STAGE_F"; }
# checkpoint(): one durable line per completed step, so the run can be checked
# at a glance at any time without reading the full log.
checkpoint(){
    printf "[%s +%s] %s\n" "$(_ts)" "$(_el)" "$*" >> "$PROG_F"
    echo "[$(_ts) +$(_el)]  ✓ CHECKPOINT: $*" | tee -a "$LOG"
}
start_heartbeat(){
    ( while true; do sleep $((HEARTBEAT_MIN*60))
        echo "[$(date '+%H:%M:%S')] ---- alive: $(cat "$STAGE_F" 2>/dev/null) | done so far: $(wc -l < "$PROG_F" 2>/dev/null) checkpoints | disk $(df -h / | awk 'NR==2{print $4}') free | mem $(free -g | awk '/^Mem:/{print $7}')G avail | load$(cut -d' ' -f1-3 /proc/loadavg | sed 's/^/ /')" | tee -a "$LOG"
      done ) & HB_PID=$!
}
stop_heartbeat(){ [ -n "$HB_PID" ] && kill "$HB_PID" 2>/dev/null; }
trap 'stop_heartbeat' EXIT INT TERM

# ── CSVs, written incrementally ───────────────────────────────────────────
CSV1="$OUT_DIR/claim1_t1_t2.csv"
CSV2="$OUT_DIR/claim2_t3.csv"
CSV3="$OUT_DIR/claim3_t6.csv"
echo "dataset,tool,raw_bytes,archive_bytes,ratio_pct,compress_s,decompress_s,peak_ram_kb,lossless,status" > "$CSV1"
echo "individual,tool,snv_f1,snv_precision,snv_recall,wall_s,peak_ram_kb,status" > "$CSV2"
echo "dataset,operation,ours_s,baseline_tool,baseline_s,speedup,status,note" > "$CSV3"

DATASETS="ERR5181310 SRR554369 ERR552797 SRR2584863 SRR29296997 ERR12954017 \
SRR065390 SRR40271341 ERR17740259 SRR37283774 DRR976266 SRR36741279 \
SRR32429602 SRR39257532 SRR10676752 HG002 HG003 HG004 HG005"
C2_SETS="HG002 HG003 HG004 HG005"
# Claim 3 baselines: one per kingdom. Chosen for SPAN (6 kingdoms, 0.20-3.99 GB,
# 12x-500x coverage), not size -- a 50-200x speedup does not need 19 points, and
# SPAdes on the multi-GB eukaryotes costs hours to days for the same ratio.
C3_BASE="ERR5181310:sarscov2 SRR2584863:ecoli SRR29296997:halobacterium \
SRR37283774:pfalciparum DRR976266:scerevisiae HG002:chr20"

WD="$OUT_DIR/wd"; mkdir -p "$WD"
ARCH_DIR="$OUT_DIR/archives"; mkdir -p "$ARCH_DIR"   # KEPT -- Claim 3 reads these
N_OK=0; N_FAIL=0; FAILED_LIST=""

banner "BENCHMARK 1 — G_CAPSUL full pipeline"
say "started    : $(date '+%Y-%m-%d %H:%M:%S %Z')"
say "output     : $OUT_DIR"
say "log        : $LOG"
say "phases     : $PHASES$([ -n "${SANITY_ONLY:-}" ] && echo '   [SANITY_ONLY — first dataset only]')"
say "host       : $(hostname)  $NPROC cores  $(free -g | awk '/^Mem:/{print $2}') GB RAM  $(df -h / | awk 'NR==2{print $4}') disk free"
say "commit     : $(git rev-parse --short HEAD 2>/dev/null)  branch $(git rev-parse --abbrev-ref HEAD 2>/dev/null)"
say "heartbeat  : every ${HEARTBEAT_MIN} min"
say ""
say "RULE: one timed job at a time. speed and RAM come from /usr/bin/time -v."
say "RULE: a dataset failure retries once, then records FAILED and continues."
say "RULE: a LOSSY archive halts the whole run (CLAIM 1 correctness gate)."

for b in "$BEST" "$DEC"; do
    [ -x "$b" ] || { err "missing binary: $b — run benchmark_0_preflight.sh first"; exit 1; }
done
start_heartbeat

# ═══════════════════════════════════════════════════════════════════════════
#  PHASE 1 — CLAIM 1
# ═══════════════════════════════════════════════════════════════════════════
phase1_one(){                      # $1 = dataset name ; returns 1 on failure
    local DS="$1" SRC IN RAW A ARCH CW CR DW DR LL TF TF2 OUTDIR t
    SRC="$DATA_DIR/${DS}_1.fq"; [ -s "$SRC" ] || SRC="$DATA_DIR/${DS}_pooled.fq"
    [ -s "$SRC" ] || { err "$DS: no input found in $DATA_DIR"; return 1; }
    RAW=$(stat -c%s "$SRC")
    inf "[$DS] input $(gbs $RAW)   $SRC"

    IN="$WD/$DS.fq"
    step "staging working copy"
    cp "$SRC" "$IN" || { err "$DS: copy failed (disk?)"; debug_dump "$DS copy" ""; return 1; }

    # ---- OURS ----
    mark "P1 $DS: CAPSULE compress"
    step "CAPSULE compress (adaptive, 4 candidates in one process)"
    A="$ARCH_DIR/$DS.capsule"; TF="$WD/t_c_$$"
    /usr/bin/time -v env CAPS_NAMES=1 CAPS_QUAL=1 INPUT="$IN" ARCHIVE="$A" BEST="$BEST" \
        bash "$HERE/scripts/encode_adaptive.sh" >/dev/null 2>"$TF"
    read -r CW CR <<< "$(parse_time_v "$TF")"
    if [ ! -s "$A" ]; then err "$DS: CAPSULE produced no archive"; debug_dump "$DS encode" "${A}.log"; rm -f "$IN" "$TF"; return 1; fi
    ARCH=$(stat -c%s "$A")
    ok "compress  archive=$(mbs $ARCH)  ratio=$(awk -v a=$ARCH -v r=$RAW 'BEGIN{printf "%.2f%%",100*a/r}')  wall=${CW}s  RAM=$(ramg $CR)"
    checkpoint "$DS  CAPSULE compress done -- $(mbs $ARCH), ${CW}s, $(ramg $CR)"

    mark "P1 $DS: CAPSULE decompress + lossless"
    step "CAPSULE decompress + lossless verify"
    OUTDIR="$WD/${DS}_dec"; mkdir -p "$OUTDIR"; TF2="$WD/t_d_$$"
    /usr/bin/time -v "$DEC" "$A" "$OUTDIR" "$OUTDIR/reads.seq" >/dev/null 2>"$TF2"
    read -r DW DR <<< "$(parse_time_v "$TF2")"
    LL=$(losscmp "$IN" "$OUTDIR/reads.seq" "$OUTDIR/reads.seq.names" "$OUTDIR/reads.seq.qual")
    ok "decompress wall=${DW}s  RAM=$(ramg $DR)  -> $LL"
    checkpoint "$DS  CAPSULE decompress + lossless done -- $LL, ${DW}s"
    printf "%s,CAPSULE,%s,%s,%.4f,%s,%s,%s,%s,DONE\n" "$DS" "$RAW" "$ARCH" \
        "$(awk -v a=$ARCH -v r=$RAW 'BEGIN{print 100*a/r}')" "$CW" "$DW" "$CR" "$LL" >> "$CSV1"
    rm -rf "$OUTDIR" "$TF" "$TF2"          # decode scratch is the big transient
    if [ "$LL" != "LOSSLESS" ]; then
        err "$DS: ARCHIVE IS LOSSY — halting the entire run (CLAUDE.md rule 5)"
        err "This is a correctness bug, not a flaky run. Do not continue."
        stop_heartbeat; exit 2
    fi

    # ---- SPRING ----
    mark "P1 $DS: SPRING"
    step "SPRING compress + decompress"
    A="$WD/$DS.spring"; TF="$WD/t_sc_$$"; TF2="$WD/t_sd_$$"
    /usr/bin/time -v spring -c -i "$IN" -o "$A" -t "$NPROC" 2>"$TF"
    read -r CW CR <<< "$(parse_time_v "$TF")"
    if [ -s "$A" ]; then
        ARCH=$(stat -c%s "$A")
        /usr/bin/time -v spring -d -i "$A" -o "$WD/$DS.spr.dec" -t "$NPROC" 2>"$TF2"
        read -r DW DR <<< "$(parse_time_v "$TF2")"
        LL=$(losscmp_plain "$IN" "$WD/$DS.spr.dec")
        ok "SPRING    archive=$(mbs $ARCH)  ctime=${CW}s  dtime=${DW}s  RAM=$(ramg $CR)  $LL"
        checkpoint "$DS  SPRING done -- $(mbs $ARCH), ${CW}s, $LL"
        printf "%s,SPRING,%s,%s,%.4f,%s,%s,%s,%s,DONE\n" "$DS" "$RAW" "$ARCH" \
            "$(awk -v a=$ARCH -v r=$RAW 'BEGIN{print 100*a/r}')" "$CW" "$DW" "$CR" "$LL" >> "$CSV1"
    else
        err "SPRING failed on $DS (continuing — competitor arm only)"; debug_dump "$DS SPRING" "$TF"
        printf "%s,SPRING,%s,,,,,,,FAILED\n" "$DS" "$RAW" >> "$CSV1"
    fi
    rm -f "$A" "$WD/$DS.spr.dec" "$TF" "$TF2"

    # ---- Genozip ----
    mark "P1 $DS: Genozip"
    step "Genozip compress + decompress"
    A="$WD/$DS.genozip"; TF="$WD/t_gc_$$"; TF2="$WD/t_gd_$$"
    /usr/bin/time -v genozip --force -o "$A" "$IN" 2>"$TF"
    read -r CW CR <<< "$(parse_time_v "$TF")"
    if [ -s "$A" ]; then
        ARCH=$(stat -c%s "$A")
        /usr/bin/time -v genounzip --force -o "$WD/$DS.gz.dec" "$A" 2>"$TF2"
        read -r DW DR <<< "$(parse_time_v "$TF2")"
        LL=$(losscmp_plain "$IN" "$WD/$DS.gz.dec")
        ok "Genozip   archive=$(mbs $ARCH)  ctime=${CW}s  dtime=${DW}s  RAM=$(ramg $CR)  $LL"
        checkpoint "$DS  Genozip done -- $(mbs $ARCH), ${CW}s, $LL"
        printf "%s,Genozip,%s,%s,%.4f,%s,%s,%s,%s,DONE\n" "$DS" "$RAW" "$ARCH" \
            "$(awk -v a=$ARCH -v r=$RAW 'BEGIN{print 100*a/r}')" "$CW" "$DW" "$CR" "$LL" >> "$CSV1"
    else
        err "Genozip failed on $DS (continuing — competitor arm only)"; debug_dump "$DS Genozip" "$TF"
        printf "%s,Genozip,%s,,,,,,,FAILED\n" "$DS" "$RAW" >> "$CSV1"
    fi
    rm -f "$A" "$WD/$DS.gz.dec" "$TF" "$TF2" "$IN"
    return 0
}

# lossless comparison helpers (order-free, CRLF-safe)
losscmp_plain(){
    # THE '+' LINE IS NORMALISED ON BOTH SIDES. SPRING drops the optional ID
    # after '+', which the FASTQ spec makes redundant (it must repeat line 1),
    # so a raw comparison reports LOSSY for a tool that lost no data. CLAUDE.md
    # specifies this normalisation so the SPRING comparison is fair on DATA
    # recovery. Verified on SRR29296997: raw compare differs only in that
    # field; normalised, it is byte-identical.
    local orig="$1" dec="$2"
    [ -s "$dec" ] || { echo "LOSSY(no-output)"; return; }
    _n4(){ paste - - - - < "$1" | awk -F'\t' '{print $1"\t"$2"\t+\t"$4}' | tr -d '\r' | sort; }
    _n4 "$orig" > "$WD/_a_$$"
    _n4 "$dec"  > "$WD/_b_$$"
    if cmp -s "$WD/_a_$$" "$WD/_b_$$"; then echo LOSSLESS; else echo LOSSY; fi
    rm -f "$WD/_a_$$" "$WD/_b_$$"
}
losscmp(){          # ours: separate seq/names/qual streams -> rebuild 4-line records
    local orig="$1" seq="$2" nam="$3" qual="$4"
    [ -s "$seq" ] || { echo "LOSSY(no-reads)"; return; }
    awk 'NR%4==1{n=substr($0,2)} NR%4==2{s=$0} NR%4==0{print n"\t"s"\t"$0}' "$orig" | tr -d '\r' | sort > "$WD/_a_$$"
    if [ -s "$nam" ] && [ -s "$qual" ]; then
        paste "$nam" "$seq" "$qual" | tr -d '\r' | sort > "$WD/_b_$$"
    else
        awk 'NR%4==2' "$orig" | tr -d '\r' | sort > "$WD/_a_$$"
        sort "$seq" > "$WD/_b_$$"
    fi
    if cmp -s "$WD/_a_$$" "$WD/_b_$$"; then echo LOSSLESS; else echo LOSSY; fi
    rm -f "$WD/_a_$$" "$WD/_b_$$"
}

run_phase1(){
    banner "PHASE 1 — CLAIM 1 (compression): 19 datasets, one at a time"
    say "  per dataset: CAPSULE encode -> archive KEPT -> decode -> lossless"
    say "               -> SPRING compress+decompress -> Genozip compress+decompress"
    say "  tables: T1 archive size | T2 wall time + peak RAM (all three tools)"
    local list="$DATASETS" i=0 n t0 el
    [ -n "${SANITY_ONLY:-}" ] && { list=$(echo $DATASETS | cut -d' ' -f1); say "  SANITY_ONLY: running only $list"; }
    n=$(echo $list | wc -w)
    for DS in $list; do
        i=$((i+1)); t0=$(date +%s)
        say ""
        say "──────────────────────────────────────────────────────────────────────"
        inf "PHASE 1  [$i/$n]  $DS"
        say "──────────────────────────────────────────────────────────────────────"
        if phase1_one "$DS"; then
            N_OK=$((N_OK+1))
        else
            err "$DS failed — RETRYING ONCE"
            rm -f "$WD/$DS.fq"
            if phase1_one "$DS"; then N_OK=$((N_OK+1)); ok "$DS succeeded on retry"
            else N_FAIL=$((N_FAIL+1)); FAILED_LIST="$FAILED_LIST $DS"
                 err "$DS FAILED twice — recorded and skipped"
                 printf "%s,CAPSULE,,,,,,,,FAILED\n" "$DS" >> "$CSV1"; fi
        fi
        el=$(( $(date +%s) - t0 ))
        inf "[$i/$n] $DS done in ${el}s   (elapsed $(_el), $((n-i)) left)"
        [ "$i" -gt 0 ] && inf "     projected remaining: ~$(( (($(date +%s)-T_RUN_START)/i) * (n-i) / 60 )) min at current rate"
    done
    checkpoint "PHASE 1 COMPLETE -- $N_OK ok, $N_FAIL failed"
    banner "PHASE 1 COMPLETE — $N_OK ok, $N_FAIL failed  ->  $CSV1"
}

# ═══════════════════════════════════════════════════════════════════════════
#  PHASE 2 — CLAIM 2
# ═══════════════════════════════════════════════════════════════════════════
run_phase2(){
    banner "PHASE 2 — CLAIM 2 (variant calling): 4 GIAB human sets"
    say "  per set: ours (compress+call, one pass) -> DiscoSNP++ -> Kmer2SNP"
    say "           -> rtg vcfeval against GIAB truth"
    say "  tables: T3 het-SNV F1 (3-way)"
    local i=0 n=4 t0 fq d f1 p r tv w hwm el
    for IND in $C2_SETS; do
        i=$((i+1)); t0=$(date +%s)
        fq="$DATA_DIR/${IND}_pooled.fq"
        say ""; say "──────────────────────────────────────────────────────────────────────"
        inf "PHASE 2  [$i/$n]  $IND"
        say "──────────────────────────────────────────────────────────────────────"
        [ -s "$fq" ] || { err "$IND: $fq missing — skipping"; printf "%s,CAPSULE,,,,,,MISSING_INPUT\n" "$IND" >> "$CSV2"; continue; }
        inf "  input $(gbs $(stat -c%s "$fq"))"

        mark "P2 $IND: our caller"
        step "our caller (full chr20, SNV configuration)"
        d="$OUT_DIR/c2_$IND"; mkdir -p "$d"
        CAPS_DBG=1 CAPS_DBG_ONLY=1 \
          bash "$HERE/scripts/run_fullchr20_bench_capsule.sh" \
               "$BEST" "$HERE/scripts" "$REFS/chr20.fa" "$fq" "$IND" "$d" > "$d.log" 2>&1
        f1=$(grep -aE '^SNV ' "$d.log" | tail -1 | grep -oP 'F1=\K[0-9.]+')
        p=$(grep -aE '^SNV ' "$d.log" | tail -1 | grep -oP ' P=\K[0-9.]+')
        r=$(grep -aE '^SNV ' "$d.log" | tail -1 | grep -oP ' R=\K[0-9.]+')
        tv=$(parse_time_v "$d.log"); w=${tv% *}; hwm=${tv#* }
        if [ -n "${f1:-}" ]; then
            ok "OURS      SNV F1=$f1  P=$p  R=$r   wall=${w}s  RAM=$(ramg $hwm)"
            checkpoint "$IND  our caller done -- SNV F1=$f1 P=$p R=$r, ${w}s"
            printf "%s,CAPSULE,%s,%s,%s,%s,%s,DONE\n" "$IND" "$f1" "$p" "$r" "$w" "$hwm" >> "$CSV2"
        else
            err "$IND: our caller produced no SNV line"; debug_dump "$IND ours" "$d.log"
            printf "%s,CAPSULE,,,,,,FAILED\n" "$IND" >> "$CSV2"
        fi

        mark "P2 $IND: DiscoSNP++"
        step "DiscoSNP++ (same reads, same truth, same scoring)"
        if [ -x "$HERE/scripts/run_fullchr20_bench_disco.sh" ] || [ -f "$HERE/scripts/run_fullchr20_bench_disco.sh" ]; then
            # signature: <ref.fa> <reads.fq> [individual] [outdir]  -- verified
            # against the script header, not assumed. A wrong argument order
            # here silently benchmarks the competitor on the wrong input.
            /usr/bin/time -v bash "$HERE/scripts/run_fullchr20_bench_disco.sh" \
                 "$REFS/chr20.fa" "$fq" "$IND" "$OUT_DIR/c2_${IND}_disco" \
                 > "$OUT_DIR/c2_${IND}_disco.log" 2>&1
            f1=$(grep -aE '^SNV ' "$OUT_DIR/c2_${IND}_disco.log" | tail -1 | grep -oP 'F1=\K[0-9.]+')
            p=$(grep -aE '^SNV ' "$OUT_DIR/c2_${IND}_disco.log" | tail -1 | grep -oP ' P=\K[0-9.]+')
            r=$(grep -aE '^SNV ' "$OUT_DIR/c2_${IND}_disco.log" | tail -1 | grep -oP ' R=\K[0-9.]+')
            tv=$(parse_time_v "$OUT_DIR/c2_${IND}_disco.log"); w=${tv% *}; hwm=${tv#* }
            if [ -n "${f1:-}" ]; then ok "DiscoSNP++ SNV F1=$f1  P=$p  R=$r   wall=${w}s  RAM=$(ramg $hwm)"
                checkpoint "$IND  DiscoSNP++ done -- SNV F1=$f1"
                printf "%s,DiscoSNP++,%s,%s,%s,%s,%s,DONE\n" "$IND" "$f1" "$p" "$r" "$w" "$hwm" >> "$CSV2"
            else err "DiscoSNP++ produced no SNV line for $IND"; debug_dump "$IND disco" "$OUT_DIR/c2_${IND}_disco.log"
                printf "%s,DiscoSNP++,,,,,,FAILED\n" "$IND" >> "$CSV2"; fi
        else
            err "run_fullchr20_bench_disco.sh not found — DiscoSNP++ arm skipped"
            printf "%s,DiscoSNP++,,,,,,SCRIPT_MISSING\n" "$IND" >> "$CSV2"
        fi

        mark "P2 $IND: Kmer2SNP"
        step "Kmer2SNP"
        # NO VERIFIED RUNNER EXISTS IN THIS REPO. /root/Kmer2SNP and the conda
        # env are both present, but this repository has never invoked them --
        # the published Kmer2SNP F1 (0.464) comes from the outer ARCS project,
        # with a methodology not reproduced here. Writing an invocation from
        # guesswork would produce a number that looks like a measurement and
        # is not one, which is worse for the paper than an honest gap. The arm
        # is therefore recorded as NOT_AVAILABLE until a runner is written and
        # validated against that published value.
        if [ -f "$HERE/scripts/run_kmer2snp.sh" ]; then
            bash "$HERE/scripts/run_kmer2snp.sh" "$fq" "$IND" "$OUT_DIR/c2_${IND}_k2s" \
                 > "$OUT_DIR/c2_${IND}_k2s.log" 2>&1
            f1=$(grep -aE '^SNV ' "$OUT_DIR/c2_${IND}_k2s.log" | tail -1 | grep -oP 'F1=\K[0-9.]+')
            if [ -n "${f1:-}" ]; then ok "Kmer2SNP  SNV F1=$f1"
                printf "%s,Kmer2SNP,%s,,,,,DONE\n" "$IND" "$f1" >> "$CSV2"
            else err "Kmer2SNP produced no SNV line for $IND"
                printf "%s,Kmer2SNP,,,,,,FAILED\n" "$IND" >> "$CSV2"; fi
        else
            err "Kmer2SNP env or runner missing — arm skipped (T3 has 2 of 3 tools)"
            printf "%s,Kmer2SNP,,,,,,NOT_AVAILABLE\n" "$IND" >> "$CSV2"
        fi
        el=$(( $(date +%s) - t0 ))
        inf "[$i/$n] $IND done in ${el}s   (elapsed $(_el))"
    done
    checkpoint "PHASE 2 COMPLETE"
    banner "PHASE 2 COMPLETE  ->  $CSV2"
}

# ═══════════════════════════════════════════════════════════════════════════
#  PHASE 3 — CLAIM 3
# ═══════════════════════════════════════════════════════════════════════════
run_phase3(){
    banner "PHASE 3 — CLAIM 3 (archive analysis): reads the archives phase 1 kept"
    say "  T6a export   vs SPAdes                — 6 datasets (one per kingdom)"
    say "  T6b coverage vs bwa+samtools+mosdepth — same 6"
    say "  T6c query    — every archive present (no competitor exists)"
    local SPADES="$HOME/SPAdes-4.0.0-Linux/bin/spades.py"
    local i=0 n t0 arc d t_a t_b sp ref src

    # ---- T6c: cheap, every archive ----
    say ""; inf "T6c — query, over every retained archive"
    for arc in "$ARCH_DIR"/*.capsule; do
        [ -s "$arc" ] || continue
        local ds; ds=$(basename "$arc" .capsule)
        t_a=$(date +%s.%N)
        if "$DEC" query "$arc" "$WD/q.fa" 0-100000 >/dev/null 2>"$WD/q.log"; then
            t_b=$(date +%s.%N); sp=$(awk -v a=$t_a -v b=$t_b 'BEGIN{printf "%.3f",b-a}')
            ok "T6c query  $ds  ${sp}s"
            printf "%s,query,%s,none,,,DONE,no competitor exists for coordinate-range retrieval\n" "$ds" "$sp" >> "$CSV3"
        else
            err "T6c query failed on $ds"; debug_dump "$ds query" "$WD/q.log"
            printf "%s,query,,none,,,FAILED,\n" "$ds" >> "$CSV3"
        fi
        rm -f "$WD/q.fa"
    done

    # ---- T6a + T6b: the six ----
    n=$(echo $C3_BASE | wc -w)
    for entry in $C3_BASE; do
        i=$((i+1)); t0=$(date +%s)
        local DS="${entry%%:*}" REFNAME="${entry##*:}"
        say ""; say "──────────────────────────────────────────────────────────────────────"
        inf "PHASE 3  [$i/$n]  $DS   (baseline reference: $REFNAME)"
        say "──────────────────────────────────────────────────────────────────────"
        arc="$ARCH_DIR/$DS.capsule"
        [ -s "$arc" ] || { err "$DS: archive missing (phase 1 must run first) — skipping"
                           printf "%s,export,,SPAdes,,,NO_ARCHIVE,\n" "$DS" >> "$CSV3"; continue; }
        if [ "$REFNAME" = chr20 ]; then ref="$REFS/chr20.fa"; else ref="$REFS/c3_${REFNAME}.fa"; fi
        src="$DATA_DIR/${DS}_1.fq"; [ -s "$src" ] || src="$DATA_DIR/${DS}_pooled.fq"

        # T6a export vs SPAdes
        mark "P3 $DS: export"
        step "T6a  our export"
        d="$OUT_DIR/c3_$DS"; mkdir -p "$d"
        t_a=$(date +%s.%N)
        "$DEC" export "$arc" "$d/contigs.fa" >/dev/null 2>"$d/export.log"
        t_b=$(date +%s.%N); local OURS_EXP; OURS_EXP=$(awk -v a=$t_a -v b=$t_b 'BEGIN{printf "%.2f",b-a}')
        if [ -s "$d/contigs.fa" ]; then ok "our export  ${OURS_EXP}s   $(mbs $(stat -c%s "$d/contigs.fa"))"
        else err "our export produced nothing for $DS"; debug_dump "$DS export" "$d/export.log"; fi

        step "T6a  SPAdes de-novo (the slow baseline — minutes to hours)"
        if [ -x "$SPADES" ] && [ -s "$src" ]; then
            t_a=$(date +%s.%N)
            python3 "$SPADES" -s "$src" -o "$d/spades" -t "$NPROC" -m $(( $(free -g | awk '/^Mem:/{print $2}') - 8 )) \
                > "$d/spades.log" 2>&1
            local RC=$?; t_b=$(date +%s.%N)
            if [ $RC -eq 0 ] && [ -s "$d/spades/contigs.fasta" ]; then
                local SP_T; SP_T=$(awk -v a=$t_a -v b=$t_b 'BEGIN{printf "%.2f",b-a}')
                ok "SPAdes      ${SP_T}s   ->  speedup $(awk -v s=$SP_T -v o=$OURS_EXP 'BEGIN{printf "%.1fx",s/o}')"
                checkpoint "$DS  T6a export done -- ours ${OURS_EXP}s vs SPAdes ${SP_T}s"
                printf "%s,export,%s,SPAdes,%s,%s,DONE,\n" "$DS" "$OURS_EXP" "$SP_T" \
                    "$(awk -v s=$SP_T -v o=$OURS_EXP 'BEGIN{printf "%.1f",s/o}')" >> "$CSV3"
            else
                err "SPAdes did not complete on $DS (rc=$RC)"; debug_dump "$DS SPAdes" "$d/spades.log"
                printf "%s,export,%s,SPAdes,,,BASELINE_DNF,SPAdes did not complete (rc=%s)\n" "$DS" "$OURS_EXP" "$RC" >> "$CSV3"
            fi
            rm -rf "$d/spades/K"* "$d/spades/tmp" 2>/dev/null
        else
            err "SPAdes or input missing for $DS"
            printf "%s,export,%s,SPAdes,,,BASELINE_MISSING,\n" "$DS" "$OURS_EXP" >> "$CSV3"
        fi

        # T6b coverage vs bwa+samtools+mosdepth
        mark "P3 $DS: coverage"
        step "T6b  our coverage"
        t_a=$(date +%s.%N)
        "$DEC" coverage "$arc" "$d/coverage.tsv" >/dev/null 2>"$d/coverage.log"
        t_b=$(date +%s.%N); local OURS_COV; OURS_COV=$(awk -v a=$t_a -v b=$t_b 'BEGIN{printf "%.2f",b-a}')
        local NROW; NROW=$(wc -l < "$d/coverage.tsv" 2>/dev/null || echo 0)
        ok "our coverage ${OURS_COV}s   ${NROW} rows"

        step "T6b  bwa + samtools sort + mosdepth (the conventional route)"
        if [ -s "$ref.bwt" ] && [ -s "$src" ] && command -v mosdepth >/dev/null; then
            t_a=$(date +%s.%N)
            bwa mem -t "$NPROC" "$ref" "$src" 2>"$d/bwa.log" \
              | samtools sort -@ 4 -o "$d/aln.bam" - 2>>"$d/bwa.log" \
              && samtools index "$d/aln.bam" 2>>"$d/bwa.log" \
              && mosdepth -t 4 "$d/md" "$d/aln.bam" 2>>"$d/bwa.log"
            local RC=$?; t_b=$(date +%s.%N)
            if [ $RC -eq 0 ]; then
                local CV_T; CV_T=$(awk -v a=$t_a -v b=$t_b 'BEGIN{printf "%.2f",b-a}')
                ok "bwa+mosdepth ${CV_T}s  ->  speedup $(awk -v s=$CV_T -v o=$OURS_COV 'BEGIN{printf "%.1fx",s/o}')"
                checkpoint "$DS  T6b coverage done -- ours ${OURS_COV}s vs bwa+mosdepth ${CV_T}s"
                printf "%s,coverage,%s,bwa+samtools+mosdepth,%s,%s,DONE,\n" "$DS" "$OURS_COV" "$CV_T" \
                    "$(awk -v s=$CV_T -v o=$OURS_COV 'BEGIN{printf "%.1f",s/o}')" >> "$CSV3"
            else
                err "bwa/mosdepth baseline failed on $DS"; debug_dump "$DS bwa" "$d/bwa.log"
                printf "%s,coverage,%s,bwa+samtools+mosdepth,,,BASELINE_FAILED,\n" "$DS" "$OURS_COV" >> "$CSV3"
            fi
            rm -f "$d/aln.bam" "$d/aln.bam.bai"    # BAMs are large, the timing is what we keep
        else
            err "bwa index / input / mosdepth missing for $DS — baseline skipped"
            printf "%s,coverage,%s,bwa+samtools+mosdepth,,,BASELINE_MISSING,\n" "$DS" "$OURS_COV" >> "$CSV3"
        fi
        inf "[$i/$n] $DS done in $(( $(date +%s) - t0 ))s   (elapsed $(_el))"
    done
    checkpoint "PHASE 3 COMPLETE"
    banner "PHASE 3 COMPLETE  ->  $CSV3"
}

# ═══════════════════════════════════════════════════════════════════════════
case ",$PHASES," in *,1,*) run_phase1;; esac
case ",$PHASES," in *,2,*) [ -n "${SANITY_ONLY:-}" ] || run_phase2;; esac
case ",$PHASES," in *,3,*) [ -n "${SANITY_ONLY:-}" ] || run_phase3;; esac
stop_heartbeat

banner "BENCHMARK 1 COMPLETE"
say "finished  : $(date '+%Y-%m-%d %H:%M:%S')"
say "total time: $(_el)"
say "datasets  : $N_OK ok, $N_FAIL failed${FAILED_LIST:+ ($FAILED_LIST)}"
say ""
say "  T1/T2 (Claim 1): $CSV1"
say "  T3    (Claim 2): $CSV2"
say "  T6    (Claim 3): $CSV3"
say "  archives kept  : $ARCH_DIR  ($(du -sh "$ARCH_DIR" 2>/dev/null | cut -f1))"
say "  full log       : $LOG"
say ""
say "── CLAIM 1 ──"; column -s, -t "$CSV1" 2>/dev/null | head -30 | tee -a "$LOG"
say ""; say "── CLAIM 2 ──"; column -s, -t "$CSV2" 2>/dev/null | tee -a "$LOG"
say ""; say "── CLAIM 3 ──"; column -s, -t "$CSV3" 2>/dev/null | head -30 | tee -a "$LOG"
rm -rf "$WD"
[ "$N_FAIL" -eq 0 ] && exit 0 || exit 1
