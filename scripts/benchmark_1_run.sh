#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
#  BENCHMARK 1 — the full run. Claims 1, 2 and 3, frame by frame.
#
#  Design rules this script obeys:
#   * ONE TIMED JOB AT A TIME. Concurrency corrupts wall time and peak RAM,
#     which are published numbers (BTR_NOTES.md rule 3). Nothing here runs in
#     parallel with anything else that is being measured.
#   * SPEED AND RAM come from /usr/bin/time -v (wall clock + VmHWM). Never
#     estimated, never derived from a log timestamp.
#   * Loud at the SEAMS, silent in the LOOPS. Every dataset, tool and phase
#     transition prints. Nothing prints from inside a hot path -- that would
#     slow the very thing being timed.
#   * RESULTS ARE WRITTEN AS THEY HAPPEN. A crash at dataset 15 keeps 1-14.
#   * FAILURE OF ONE DATASET IS NOT FAILURE OF THE RUN: retry once, then
#     record FAILED and continue. The single exception is a LOSSY archive,
#     which halts everything (BTR_NOTES.md rule 5) -- that is a correctness bug,
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
# DiscoSNP++ ships as a shell script in its own tree and its runner resolves
# `run_discoSnp++.sh` through PATH. The file existing is NOT enough -- this
# project has lost a competitor arm to exactly that before (see
# docs/INVOCATION_ERRORS.md), and it fails as "no SNV line", which looks like a
# scoring problem rather than a missing tool.
[ -d "$HOME/DiscoSnp" ] && export PATH="$HOME/DiscoSnp:$PATH"
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
    # TAKE THE LAST BLOCK ONLY. Some arms redirect /usr/bin/time -v into the
    # SAME log as the tool's own output, and those inner scripts run time -v
    # too -- so the file holds several "Elapsed" lines. Without tail -1, awk
    # printed one number per match and they concatenated: DiscoSNP++ reported
    # "wall=78.3384.97s", i.e. 78.33 and 84.97 glued together, which would have
    # gone into the published T3 table as a single bogus figure. The outer
    # timer finishes last, so the last block is the one that timed the whole
    # invocation. Peak RAM takes the max, not the last, since the heaviest
    # child is the honest answer for a pipeline.
    local f="$1" wall hwm
    wall=$(grep "Elapsed (wall clock)" "$f" 2>/dev/null | tail -1 | awk '{n=split($NF,a,":");
        if(n==3) printf "%.2f",a[1]*3600+a[2]*60+a[3];
        else if(n==2) printf "%.2f",a[1]*60+a[2]; else printf "%.2f",a[1]}')
    hwm=$(grep "Maximum resident set size" "$f" 2>/dev/null | awk '{if($NF+0>m)m=$NF+0}END{printf "%d",m}')
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
# ── projection from benchmark_0, so a wrong number is visible at minute 19
# and not at hour 9. Absent file = projections simply not shown; never fatal,
# and never used to judge a result (BTR_NOTES.md rule 6: measured wins).
PROJ_F="${PROJ_FILE:-$HERE/results/BENCHMARK_0_PROJECTION.tsv}"
proj_ds(){   # $1=accession -> "ours_c ours_d spring_c geno_c", empty if unknown
    [ -s "$PROJ_F" ] || return 0
    awk -v d="$1" '$1=="DS" && $2==d {print $4,$5,$6,$7; exit}' "$PROJ_F"; }
proj_key(){ [ -s "$PROJ_F" ] || return 0; awk -v k="$1" '$1==k{print $2; exit}' "$PROJ_F"; }
vs_proj(){   # $1=actual seconds  $2=projected seconds
    [ -n "${2:-}" ] || return 0
    awk -v a="$1" -v p="$2" 'BEGIN{ if(p<=0){print "";exit}
        r=a/p; printf "  [proj %.1fs, %.2fx %s]", p, r, (r>1.5)?"SLOWER THAN PROJECTED":((r<0.67)?"faster":"on track") }'; }
CSV1="$OUT_DIR/claim1_T1.1_T1.2.csv"
CSV2="$OUT_DIR/claim2_T2.1_snv.csv"
CSV3="$OUT_DIR/claim3_T3.1_T3.2_T3.3.csv"
echo "dataset,tool,raw_bytes,archive_bytes,ratio_pct,compress_s,decompress_s,peak_ram_kb,lossless,status" > "$CSV1"
echo "individual,tool,tp,fp,fn,precision,recall,f1,wall_s,peak_ram_kb,status" > "$CSV2"
echo "table,dataset,operation,ours_s,ours_peak_ram_kb,output_bytes,rows,baseline_tool,baseline_s,baseline_peak_ram_kb,speedup,status,note" > "$CSV3"
CSV4="$OUT_DIR/claim2_T2.2_coverage_sweep.csv"
CSV5="$OUT_DIR/claim2_T2.3_indel.csv"
CSV6="$OUT_DIR/claim2_T2.4_multiallelic.csv"
CSV7="$OUT_DIR/claim2_T2.5_tetraploid.csv"
CSV8="$OUT_DIR/claim3_T3.4_locus_fidelity.csv"
echo "individual,depth_x,reads,archive_bytes,tp,fp,fn,precision,recall,f1,wall_s,peak_ram_kb,status" > "$CSV4"
echo "individual,tool,tp,fp,fn,precision,recall,f1,wall_s,peak_ram_kb,status" > "$CSV5"
echo "individual,region,truth_multiallelic_sites,capsule_anycall,disco_anycall,capsule_both_alleles,disco_both_alleles,capsule_rate,disco_rate,wall_s,peak_ram_kb,status" > "$CSV6"
echo "pair,ploidy,region,tool,class,tp,fp,fn,precision,recall,f1,wall_s,peak_ram_kb,status" > "$CSV7"
echo "individual,sites,coord_both,coord_one,coord_none,content_both,content_one,content_none,ref_reads,alt_reads,status" > "$CSV8"

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
    # CAPS_SPANS=1 writes `contig_spans`, which Phase 2 needs to call from the
    # archive. NOT CAPS_CALL=1: that also runs the full caller inline during
    # compression (~20x heavier). Measured: 38.64 s / 2.45 GB vs 8.97 s /
    # 1.00 GB on the same input, for a byte-identical archive.
    /usr/bin/time -v env CAPS_SPANS=1 CAPS_NAMES=1 CAPS_QUAL=1 INPUT="$IN" ARCHIVE="$A" BEST="$BEST" \
        bash "$HERE/scripts/encode_adaptive.sh" >/dev/null 2>"$TF"
    read -r CW CR <<< "$(parse_time_v "$TF")"
    if [ ! -s "$A" ]; then err "$DS: CAPSULE produced no archive"; debug_dump "$DS encode" "${A}.log"; rm -f "$IN" "$TF"; return 1; fi
    ARCH=$(stat -c%s "$A")
    read -r PJC PJD PJS PJG <<< "$(proj_ds "$DS")"
    ok "compress  archive=$(mbs $ARCH)  ratio=$(awk -v a=$ARCH -v r=$RAW 'BEGIN{printf "%.2f%%",100*a/r}')  wall=${CW}s  RAM=$(ramg $CR)$(vs_proj "$CW" "${PJC:-}")"
    checkpoint "$DS  CAPSULE compress done -- $(mbs $ARCH), ${CW}s, $(ramg $CR)"

    mark "P1 $DS: CAPSULE decompress + lossless"
    step "CAPSULE decompress + lossless verify"
    OUTDIR="$WD/${DS}_dec"; mkdir -p "$OUTDIR"; TF2="$WD/t_d_$$"
    /usr/bin/time -v "$DEC" "$A" "$OUTDIR" "$OUTDIR/reads.seq" >/dev/null 2>"$TF2"
    read -r DW DR <<< "$(parse_time_v "$TF2")"
    LL=$(losscmp "$IN" "$OUTDIR/reads.seq" "$OUTDIR/reads.seq.names" "$OUTDIR/reads.seq.qual")
    ok "decompress wall=${DW}s  RAM=$(ramg $DR)  -> $LL$(vs_proj "$DW" "${PJD:-}")"
    checkpoint "$DS  CAPSULE decompress + lossless done -- $LL, ${DW}s"
    printf "%s,CAPSULE,%s,%s,%.4f,%s,%s,%s,%s,DONE\n" "$DS" "$RAW" "$ARCH" \
        "$(awk -v a=$ARCH -v r=$RAW 'BEGIN{print 100*a/r}')" "$CW" "$DW" "$CR" "$LL" >> "$CSV1"
    rm -rf "$OUTDIR" "$TF" "$TF2"          # decode scratch is the big transient
    if [ "$LL" != "LOSSLESS" ]; then
        err "$DS: ARCHIVE IS LOSSY — halting the entire run (BTR_NOTES.md rule 5)"
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
        ok "SPRING    archive=$(mbs $ARCH)  ctime=${CW}s  dtime=${DW}s  RAM=$(ramg $CR)  $LL$(vs_proj "$CW" "${PJS:-}")"
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
    # genozip tests curl/wget availability with file_exists("/dev/stdout") and,
    # on a Student licence, will not write the archive header unless it can
    # upload telemetry. If /dev/stdout is missing it exits 1 with NO archive
    # after compressing the whole file -- a silently blank competitor column.
    [ -e /dev/stdout ] || ln -sfn /proc/self/fd/1 /dev/stdout 2>/dev/null || true
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
    # so a raw comparison reports LOSSY for a tool that lost no data. BTR_NOTES.md
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
    say "  tables: T1.1 archive size | T1.2 wall time + peak RAM (all three tools)"
    local list="$DATASETS" i=0 n t0 el
    [ -n "${SANITY_ONLY:-}" ] && { list="${SANITY_DS:-$(echo $DATASETS | cut -d' ' -f1)}"; say "  SANITY_ONLY: running only $list"; }
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
    P1P=$(proj_key PHASE1_ALL_S)
    [ -n "${P1P:-}" ] && say "  phase 1 wall: $(_el)   projected ${P1P}s"
    banner "PHASE 1 COMPLETE — $N_OK ok, $N_FAIL failed  ->  $CSV1"
}

# ═══════════════════════════════════════════════════════════════════════════
#  PHASE 2 — CLAIM 2
# ═══════════════════════════════════════════════════════════════════════════
run_phase2(){
    banner "PHASE 2 — CLAIM 2 (variant calling): 4 GIAB human sets"
    # A rehearsal that skips a phase does not rehearse it. Under SANITY_ONLY we
    # still ENTER phase 2, restricted to the GIAB sets whose archive phase 1
    # actually produced -- so the calling path is exercised when it can be, and
    # says so loudly when it cannot, instead of being silently jumped over.
    if [ -n "${SANITY_ONLY:-}" ]; then
        local have=""
        for _i in $C2_SETS; do [ -s "$ARCH_DIR/$_i.capsule" ] && have="$have $_i"; done
        if [ -z "$have" ]; then
            err "SANITY_ONLY: no GIAB archive from phase 1 -- the Claim 2 path is NOT rehearsed."
            err "  To rehearse it end to end:  SANITY_ONLY=1 SANITY_DS=HG002 bash scripts/benchmark_1_run.sh"
            return 0
        fi
        C2_SETS="$(echo $have | cut -d' ' -f1)"
        say "  SANITY_ONLY: only $C2_SETS"
    fi
    say "  per set: ours (compress+call, one pass) -> DiscoSNP++ -> Kmer2SNP"
    say "           -> rtg vcfeval against GIAB truth"
    say "  tables: T2.1 het-SNV F1 (3-way)"
    local i=0 n=4 t0 fq d f1 p r tv w hwm el
    for IND in $C2_SETS; do
        i=$((i+1)); t0=$(date +%s)
        fq="$DATA_DIR/${IND}_pooled.fq"
        say ""; say "──────────────────────────────────────────────────────────────────────"
        inf "PHASE 2  [$i/$n]  $IND"
        say "──────────────────────────────────────────────────────────────────────"
        [ -s "$fq" ] || { err "$IND: $fq missing — skipping"; printf "%s,CAPSULE,,,,,,,,,MISSING_INPUT\n" "$IND" >> "$CSV2"; continue; }
        inf "  input $(gbs $(stat -c%s "$fq"))"

        mark "P2 $IND: our caller (FROM THE ARCHIVE)"
        step "our caller: capsule_decode call <archive> -- no FASTQ is read"
        d="$OUT_DIR/c2_$IND"; mkdir -p "$d"
        # THE ARCHIVE PATH IS WHAT CLAIM 2 ASSERTS. Phase 1 already built and
        # KEPT this archive with CAPS_CALL=1, so calling from it here is both
        # the honest measurement and free of a second compression pass.
        # The FASTQ runner (run_fullchr20_bench_capsule.sh) measures the encoder
        # assembling and calling in one pass -- a different, much heavier
        # operation -- and is retained for comparison, not used here.
        c2arc="$ARCH_DIR/$IND.capsule"
        if [ ! -s "$c2arc" ]; then
            err "$IND: no archive at $c2arc -- Phase 1 must run first (it is KEPT for this)"
            printf "%s,CAPSULE,,,,,,,,,NO_ARCHIVE\n" "$IND" >> "$CSV2"
            continue
        fi
        bash "$HERE/scripts/run_fullchr20_archive_capsule.sh" \
               "$DEC" "$HERE/scripts" "$REFS/chr20.fa" "$c2arc" "$IND" "$d" > "$d.log" 2>&1
        SL=$(grep -aE '^SNV ' "$d.log" | tail -1)
        f1=$(echo "$SL" | grep -oP 'F1=\K[0-9.]+'); p=$(echo "$SL" | grep -oP ' P=\K[0-9.]+')
        r=$(echo "$SL" | grep -oP ' R=\K[0-9.]+');  tp=$(echo "$SL" | grep -oP 'TP=\K[0-9]+')
        fp=$(echo "$SL" | grep -oP 'FP=\K[0-9]+');  fn=$(echo "$SL" | grep -oP 'FN=\K[0-9]+')
        tv=$(parse_time_v "$d.log"); w=${tv% *}; hwm=${tv#* }
        if [ -n "${f1:-}" ]; then
            ok "OURS      SNV F1=$f1  P=$p  R=$r   wall=${w}s  RAM=$(ramg $hwm)"
            checkpoint "$IND  our caller done -- SNV F1=$f1 P=$p R=$r, ${w}s"
            printf "%s,CAPSULE,%s,%s,%s,%s,%s,%s,%s,%s,DONE\n" "$IND" \
                "${tp:-}" "${fp:-}" "${fn:-}" "$p" "$r" "$f1" "$w" "$hwm" >> "$CSV2"
            # T2.3: both runners already score INDEL in the same pass -- the line
            # was being printed and thrown away. No extra compute.
            IL=$(grep -aE '^INDEL ' "$d.log" | tail -1)
            if [ -n "${IL:-}" ]; then
                printf "%s,CAPSULE,%s,%s,%s,%s,%s,%s,%s,%s,DONE\n" "$IND" \
                    "$(echo "$IL"|grep -oP 'TP=\K[0-9]+')" "$(echo "$IL"|grep -oP 'FP=\K[0-9]+')" \
                    "$(echo "$IL"|grep -oP 'FN=\K[0-9]+')" "$(echo "$IL"|grep -oP ' P=\K[0-9.]+')" \
                    "$(echo "$IL"|grep -oP ' R=\K[0-9.]+')" "$(echo "$IL"|grep -oP 'F1=\K[0-9.]+')" \
                    "$w" "$hwm" >> "$CSV5"
                ok "OURS      INDEL $(echo "$IL"|grep -oP 'F1=\K[0-9.]+')"
            else printf "%s,CAPSULE,,,,,,,,,NO_INDEL_LINE\n" "$IND" >> "$CSV5"; fi
        else
            err "$IND: our caller produced no SNV line"; debug_dump "$IND ours" "$d.log"
            printf "%s,CAPSULE,,,,,,,,,FAILED\n" "$IND" >> "$CSV2"
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
            SL=$(grep -aE '^SNV ' "$OUT_DIR/c2_${IND}_disco.log" | tail -1)
            f1=$(echo "$SL" | grep -oP 'F1=\K[0-9.]+'); p=$(echo "$SL" | grep -oP ' P=\K[0-9.]+')
            r=$(echo "$SL" | grep -oP ' R=\K[0-9.]+');  tp=$(echo "$SL" | grep -oP 'TP=\K[0-9]+')
            fp=$(echo "$SL" | grep -oP 'FP=\K[0-9]+');  fn=$(echo "$SL" | grep -oP 'FN=\K[0-9]+')
            tv=$(parse_time_v "$OUT_DIR/c2_${IND}_disco.log"); w=${tv% *}; hwm=${tv#* }
            if [ -n "${f1:-}" ]; then ok "DiscoSNP++ SNV F1=$f1  P=$p  R=$r   wall=${w}s  RAM=$(ramg $hwm)"
                checkpoint "$IND  DiscoSNP++ done -- SNV F1=$f1"
                printf "%s,DiscoSNP++,%s,%s,%s,%s,%s,%s,%s,%s,DONE\n" "$IND" \
                    "${tp:-}" "${fp:-}" "${fn:-}" "$p" "$r" "$f1" "$w" "$hwm" >> "$CSV2"
                IL=$(grep -aE '^INDEL ' "$OUT_DIR/c2_${IND}_disco.log" | tail -1)
                if [ -n "${IL:-}" ]; then
                    printf "%s,DiscoSNP++,%s,%s,%s,%s,%s,%s,%s,%s,DONE\n" "$IND" \
                        "$(echo "$IL"|grep -oP 'TP=\K[0-9]+')" "$(echo "$IL"|grep -oP 'FP=\K[0-9]+')" \
                        "$(echo "$IL"|grep -oP 'FN=\K[0-9]+')" "$(echo "$IL"|grep -oP ' P=\K[0-9.]+')" \
                        "$(echo "$IL"|grep -oP ' R=\K[0-9.]+')" "$(echo "$IL"|grep -oP 'F1=\K[0-9.]+')" \
                        "$w" "$hwm" >> "$CSV5"
                    ok "DiscoSNP++ INDEL $(echo "$IL"|grep -oP 'F1=\K[0-9.]+')"
                else printf "%s,DiscoSNP++,,,,,,,,,NO_INDEL_LINE\n" "$IND" >> "$CSV5"; fi
            else err "DiscoSNP++ produced no SNV line for $IND"; debug_dump "$IND disco" "$OUT_DIR/c2_${IND}_disco.log"
                printf "%s,DiscoSNP++,,,,,,,,,FAILED\n" "$IND" >> "$CSV2"; fi
        else
            err "run_fullchr20_bench_disco.sh not found — DiscoSNP++ arm skipped"
            printf "%s,DiscoSNP++,,,,,,,,,SCRIPT_MISSING\n" "$IND" >> "$CSV2"
        fi

        mark "P2 $IND: Kmer2SNP"
        step "Kmer2SNP"
        # RUNNER EXISTS AND IS VALIDATED as of 2026-09-08 (scripts/run_kmer2snp.sh).
        # Kmer2SNP's own DSK/findGSE wrappers hardcode paths that do not exist
        # here; the runner bypasses them by passing --t1/--c1/--c2/--r directly,
        # counting k-mers with KMC and deriving the coverage band from the
        # histogram. Verified two ways: the converter reproduces the archived
        # 2026-09-02 VCF byte-identically (131 records), and a full chr20 run
        # gives F1=0.464 (TP=13565 FP=290 FN=31010), matching the window figure
        # this project has been quoting (0.4636) at 111x the scale.
        # The branch below stays conditional: if the runner is ever absent the
        # arm records NOT_AVAILABLE rather than silently leaving T3 two-armed.
        if [ -f "$HERE/scripts/run_kmer2snp.sh" ]; then
            /usr/bin/time -v bash "$HERE/scripts/run_kmer2snp.sh" "$fq" "$IND" "$OUT_DIR/c2_${IND}_k2s" \
                 > "$OUT_DIR/c2_${IND}_k2s.log" 2>&1
            SL=$(grep -aE '^SNV ' "$OUT_DIR/c2_${IND}_k2s.log" | tail -1)
            f1=$(echo "$SL" | grep -oP 'F1=\K[0-9.]+'); p=$(echo "$SL" | grep -oP ' P=\K[0-9.]+')
            r=$(echo "$SL" | grep -oP ' R=\K[0-9.]+');  tp=$(echo "$SL" | grep -oP 'TP=\K[0-9]+')
            fp=$(echo "$SL" | grep -oP 'FP=\K[0-9]+');  fn=$(echo "$SL" | grep -oP 'FN=\K[0-9]+')
            tv=$(parse_time_v "$OUT_DIR/c2_${IND}_k2s.log"); w=${tv% *}; hwm=${tv#* }
            if [ -n "${f1:-}" ]; then ok "Kmer2SNP  SNV F1=$f1  P=$p  R=$r   wall=${w}s  RAM=$(ramg $hwm)"
                checkpoint "$IND  Kmer2SNP done -- SNV F1=$f1"
                printf "%s,Kmer2SNP,%s,%s,%s,%s,%s,%s,%s,%s,DONE\n" "$IND" \
                    "${tp:-}" "${fp:-}" "${fn:-}" "$p" "$r" "$f1" "$w" "$hwm" >> "$CSV2"
                # Kmer2SNP emits SNP k-mer PAIRS only; it has no indel model at
                # all. Recorded explicitly so T5's blank is a property of the
                # tool, not a gap in this benchmark.
                printf "%s,Kmer2SNP,,,,,,,,,NOT_APPLICABLE_SNP_ONLY\n" "$IND" >> "$CSV5"
            else err "Kmer2SNP produced no SNV line for $IND"
                printf "%s,Kmer2SNP,,,,,,,,,FAILED\n" "$IND" >> "$CSV2"; fi
        else
            err "Kmer2SNP env or runner missing — arm skipped (T3 has 2 of 3 tools)"
            printf "%s,Kmer2SNP,,,,,,,,,NOT_AVAILABLE\n" "$IND" >> "$CSV2"
        fi
        el=$(( $(date +%s) - t0 ))
        inf "[$i/$n] $IND done in ${el}s   (elapsed $(_el))"
    done
    # ═══ T2.2 — COVERAGE SWEEP (HG002 only) ══════════════════════════════════
    # Claim 2 is T2.1 + T2.2 + T2.3 (+T2.4/T2.5). T2.2 asks how F1 holds up as depth
    # falls, which is the question a reviewer asks of any k-mer/graph caller:
    # does it only work at luxurious coverage? It needs REAL runs -- the reads
    # are subsampled, re-compressed, and called from the resulting archive by
    # exactly the path T3 uses, so the only variable is depth.
    #
    # HG002_pooled.fq is standardised to 30x (DATASET_LOCKED.md), so the
    # fraction for a target depth is target/30 -- a formula over a declared
    # property of the input, not a fitted constant. 30x itself is not re-run:
    # it IS the T3 row, and re-running it would spend 6 minutes to reproduce a
    # number we already have.
    local T4_IND=HG002 BASE_DEPTH=30
    local t4fq="$DATA_DIR/${T4_IND}_pooled.fq"
    if [ -s "$t4fq" ] && case ",$PHASES," in *,2,*) true;; *) false;; esac; then
        banner "T2.2 — coverage sweep ($T4_IND, ${BASE_DEPTH}x source)"
        local BASE_READS; BASE_READS=$(( $(wc -l < "$t4fq") / 4 ))
        # carry the 30x row over from T3 so the sweep is complete in one table
        local r30; r30=$(awk -F, -v i="$T4_IND" '$1==i && $2=="CAPSULE" && $11=="DONE"{print $3","$4","$5","$6","$7","$8","$9","$10; exit}' "$CSV2")
        if [ -n "${r30:-}" ]; then
            printf "%s,%s,%s,%s,%s,FROM_T3\n" "$T4_IND" "$BASE_DEPTH" "$BASE_READS" \
                "$(stat -c%s "$ARCH_DIR/$T4_IND.capsule" 2>/dev/null || echo 0)" "$r30" >> "$CSV4"
        fi
        for DEPTH in ${T4_DEPTHS:-10 15 20}; do
            mark "T2.2 $T4_IND: ${DEPTH}x"
            step "T2.2  ${DEPTH}x  (subsample -> compress -> call from archive)"
            local frac sub arc4 d4 t0d
            t0d=$(date +%s)
            frac=$(awk -v d="$DEPTH" -v b="$BASE_DEPTH" 'BEGIN{printf "%.4f", d/b}')
            sub="$WD/${T4_IND}_${DEPTH}x.fq"; arc4="$WD/${T4_IND}_${DEPTH}x.capsule"
            # fixed seed: the sweep must be reproducible run to run
            seqtk sample -s11 "$t4fq" "$frac" > "$sub" 2>"$WD/seqtk.log" \
                || { err "T2.2 ${DEPTH}x: seqtk failed"; printf "%s,%s,,,,,,,,,,,SUBSAMPLE_FAILED\n" "$T4_IND" "$DEPTH" >> "$CSV4"; continue; }
            local NR4; NR4=$(( $(wc -l < "$sub") / 4 ))
            /usr/bin/time -v env CAPS_SPANS=1 CAPS_NAMES=1 CAPS_QUAL=1 INPUT="$sub" ARCHIVE="$arc4" BEST="$BEST" \
                bash "$HERE/scripts/encode_adaptive.sh" >/dev/null 2>"$WD/t4c.txt"
            if [ ! -s "$arc4" ]; then
                err "T2.2 ${DEPTH}x: no archive"; debug_dump "T2.2 ${DEPTH}x encode" "${arc4}.log"
                printf "%s,%s,%s,,,,,,,,,,ENCODE_FAILED\n" "$T4_IND" "$DEPTH" "$NR4" >> "$CSV4"
                rm -f "$sub"; continue; fi
            d4="$OUT_DIR/c2_${T4_IND}_${DEPTH}x"
            bash "$HERE/scripts/run_fullchr20_archive_capsule.sh" \
                 "$DEC" "$HERE/scripts" "$REFS/chr20.fa" "$arc4" "$T4_IND" "$d4" > "$d4.log" 2>&1
            local SL4; SL4=$(grep -aE '^SNV ' "$d4.log" | tail -1)
            local tv4 w4 h4; tv4=$(parse_time_v "$d4.log"); w4=${tv4% *}; h4=${tv4#* }
            if [ -n "${SL4:-}" ]; then
                ok "T2.2 ${DEPTH}x  F1=$(echo "$SL4"|grep -oP 'F1=\K[0-9.]+')  ($NR4 reads, $(mbs $(stat -c%s "$arc4")))"
                checkpoint "$T4_IND T2.2 ${DEPTH}x done -- F1=$(echo "$SL4"|grep -oP 'F1=\K[0-9.]+')"
                printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,DONE\n" "$T4_IND" "$DEPTH" "$NR4" \
                    "$(stat -c%s "$arc4")" \
                    "$(echo "$SL4"|grep -oP 'TP=\K[0-9]+')" "$(echo "$SL4"|grep -oP 'FP=\K[0-9]+')" \
                    "$(echo "$SL4"|grep -oP 'FN=\K[0-9]+')" "$(echo "$SL4"|grep -oP ' P=\K[0-9.]+')" \
                    "$(echo "$SL4"|grep -oP ' R=\K[0-9.]+')" "$(echo "$SL4"|grep -oP 'F1=\K[0-9.]+')" \
                    "$w4" "$h4" >> "$CSV4"
            else
                err "T2.2 ${DEPTH}x: no SNV line"; debug_dump "T2.2 ${DEPTH}x" "$d4.log"
                printf "%s,%s,%s,%s,,,,,,,,,NO_SNV_LINE\n" "$T4_IND" "$DEPTH" "$NR4" "$(stat -c%s "$arc4")" >> "$CSV4"
            fi
            rm -f "$sub" "$arc4"          # the transient is the whole point of deleting it
            inf "T2.2 ${DEPTH}x took $(( $(date +%s) - t0d ))s"
        done
        banner "T2.2 COMPLETE  ->  $CSV4"
    fi

    # ═══ T2.4 — MULTI-ALLELIC, and T2.5 — TETRAPLOID ═══════════════════════
    # Both are locked parts of Claim 2 (docs/CLAIM2_TABLES_AND_INDEL_SCAN.md)
    # and both had working runners that this benchmark never called, so a full
    # sweep produced 6 of the 8 tables and looked complete.
    #
    # They are REGION benchmarks, not whole-chr20, by design: multi-allelic
    # truth sites are rare (HG002 has 952 on all of chr20) and the tetraploid
    # construction concatenates two individuals' reads, so both are scoped to a
    # window and that window is recorded in the table rather than implied.
    # Both drive the ENCODER, not the decoder.
    if case ",$PHASES," in *,2,*) true;; *) false;; esac; then
      if [ -x "$BEST" ]; then
        # ---- T2.4 multi-allelic (one diploid individual's own GT=1/2 sites) ----
        mark "T2.4 multi-allelic"
        banner "T2.4 — multi-allelic sites recovered"
        local MA_IND="${T24_IND:-HG002}" MA_REG="${T24_REGION:-20:1000000-6000000}" mad
        mad="$OUT_DIR/c2_T2.4_${MA_IND}"
        if /usr/bin/time -v bash "$HERE/scripts/run_multiallelic_bench_capsule.sh" \
              "$BEST" "$HERE/scripts" "$REFS/chr20.fa" "$MA_IND" "$MA_REG" "$mad" \
              > "$mad.log" 2>&1; then
            local NT CH DH tvm wm hm
            NT=$(grep -aoP 'truth multi-allelic sites: \K[0-9]+' "$mad.log" | tail -1)
            CH=$(grep -aoP 'CAPSULE\s+sites with a call at that position: \K[0-9]+' "$mad.log" | tail -1)
            DH=$(grep -aoP 'DiscoSNP\+\+ sites with a call at that position: \K[0-9]+' "$mad.log" | tail -1)
            # STRICT is the metric T2.4 claims; anycall is a weaker proxy kept
            # for continuity. They disagree in DIRECTION on real data
            # (anycall 6 vs 15 = loss; strict 5 vs 0 = win), so recording only
            # one of them would be a choice, not a measurement.
            local CS DS
            CS=$(grep -aoP 'CAPSULE\s+both-allele sites: \K[0-9]+' "$mad.log" | tail -1)
            DS=$(grep -aoP 'DiscoSNP\+\+ both-allele sites: \K[0-9]+' "$mad.log" | tail -1)
            tvm=$(parse_time_v "$mad.log"); wm=${tvm% *}; hm=${tvm#* }
            if [ -n "${NT:-}" ] && [ -n "${CH:-}" ]; then
                ok "T2.4  strict both-allele: CAPSULE ${CS:-?}/$NT  DiscoSNP++ ${DS:-NOT_RUN}/$NT   (anycall proxy $CH vs ${DH:-NOT_RUN})"
                checkpoint "T2.4 multi-allelic done -- $CH/$NT vs ${DH:-0}/$NT"
                printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,DONE\n" "$MA_IND" "$MA_REG" "$NT" \
                    "$CH" "${DH:-NOT_RUN}" "${CS:-}" "${DS:-NOT_RUN}" \
                    "$(awk -v a=${CS:-0} -v n=${NT:-1} 'BEGIN{printf "%.3f",(n?a/n:0)}')" \
                    "$(awk -v a=${DS:-0} -v n=${NT:-1} 'BEGIN{printf "%.3f",(n?a/n:0)}')" \
                    "$wm" "$hm" >> "$CSV6"
            else err "T2.4: could not parse counts"; debug_dump "T2.4" "$mad.log"
                 printf "%s,%s,,,,,,,,,,PARSE_FAILED\n" "$MA_IND" "$MA_REG" >> "$CSV6"; fi
        else err "T2.4 runner failed"; debug_dump "T2.4" "$mad.log"
             printf "%s,%s,,,,,,,,,,FAILED\n" "$MA_IND" "$MA_REG" >> "$CSV6"; fi

        # ---- T2.5 tetraploid (two real diploids concatenated, Cooke 2022) ----
        mark "T2.5 tetraploid"
        banner "T2.5 — tetraploid SNV + indel"
        local TA="${T25_A:-HG003}" TB="${T25_B:-HG004}" TP4="${T25_PLOIDY:-4}"
        local TREG="${T25_REGION:-20:3000000-3400000}" ted
        ted="$OUT_DIR/c2_T2.5_${TA}_${TB}"
        if /usr/bin/time -v bash "$HERE/scripts/run_tetraploid_bench_capsule.sh" \
              "$BEST" "$HERE/scripts" "$REFS/chr20.fa" "$TA" "$TB" "$TP4" "$TREG" "$ted" \
              > "$ted.log" 2>&1; then
            # P= and R= MUST be anchored. The score line is
            #   CAPSULE_SNV TP=430 FP=24 FN=143 P=0.947 R=0.750 F1=0.837
            # so a bare `P=\K` also matches the P inside TP= and FP=, returning
            # THREE values; printf then reused its format string and emitted
            # three mangled lines per tool instead of one row. Caught by running
            # it -- the schema check passes either way, because the format
            # string is correct; it is the ARGUMENTS that multiplied.
            local tvt wt ht wrote=0
            tvt=$(parse_time_v "$ted.log"); wt=${tvt% *}; ht=${tvt#* }
            while read -r nm rest; do
                [ -n "${nm:-}" ] || continue
                local tool cls
                case "$nm" in CAPSULE_*) tool=CAPSULE;; DiscoSNP*) tool="DiscoSNP++";; *) continue;; esac
                case "$nm" in *_SNV) cls=SNV;; *_INDEL) cls=INDEL;; *) continue;; esac
                printf "%s+%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,DONE\n" "$TA" "$TB" "$TP4" "$TREG" \
                    "$tool" "$cls" \
                    "$(echo "$rest"|grep -oP 'TP=\K[0-9]+')" "$(echo "$rest"|grep -oP 'FP=\K[0-9]+')" \
                    "$(echo "$rest"|grep -oP 'FN=\K[0-9]+')" "$(echo "$rest"|grep -oP '(?<![A-Z])P=\K[0-9.]+')" \
                    "$(echo "$rest"|grep -oP '(?<![A-Z0-9])R=\K[0-9.]+')" "$(echo "$rest"|grep -oP 'F1=\K[0-9.]+')" \
                    "$wt" "$ht" >> "$CSV7"
                wrote=$((wrote+1))
            done < <(grep -aE '^(CAPSULE|DiscoSNP)[A-Za-z+_]*_(SNV|INDEL) ' "$ted.log")
            if [ "$wrote" -gt 0 ]; then ok "T2.5  $wrote rows written"; checkpoint "T2.5 tetraploid done -- $wrote rows"
            else err "T2.5: no score lines parsed"; debug_dump "T2.5" "$ted.log"
                 printf "%s+%s,%s,%s,,,,,,,,,,,NO_SCORE_LINES\n" "$TA" "$TB" "$TP4" "$TREG" >> "$CSV7"; fi
        else err "T2.5 runner failed"; debug_dump "T2.5" "$ted.log"
             printf "%s+%s,%s,%s,,,,,,,,,,,FAILED\n" "$TA" "$TB" "$TP4" "$TREG" >> "$CSV7"; fi
      else
        err "encoder $BEST missing -- T2.4 and T2.5 skipped"
      fi
    fi

    checkpoint "PHASE 2 COMPLETE"
    banner "PHASE 2 COMPLETE  ->  $CSV2"
}

# ═══════════════════════════════════════════════════════════════════════════
#  PHASE 3 — CLAIM 3
# ═══════════════════════════════════════════════════════════════════════════
run_phase3(){
    banner "PHASE 3 — CLAIM 3 (archive analysis): reads the archives phase 1 kept"
    say "  T3.1 export   vs SPAdes                — 6 datasets (one per kingdom)"
    say "  T3.2 coverage vs bwa+samtools+mosdepth — same 6"
    say "  T3.3 query    — every archive present (partial random access; see note)"
    local SPADES="$HOME/SPAdes-4.0.0-Linux/bin/spades.py"
    local i=0 n t0 arc d t_a t_b sp ref src

    # ---- T6c: cheap, every archive ----
    say ""; inf "T3.3 — query, over every retained archive"
    for arc in "$ARCH_DIR"/*.capsule; do
        [ -s "$arc" ] || continue
        local ds; ds=$(basename "$arc" .capsule)
        if /usr/bin/time -v "$DEC" query "$arc" "$WD/q.fa" 0-100000 >/dev/null 2>"$WD/q.log"; then
            read -r sp qram <<< "$(parse_time_v "$WD/q.log")"
            local QB QR; QB=$(stat -c%s "$WD/q.fa" 2>/dev/null || echo 0); QR=$(grep -c . "$WD/q.fa" 2>/dev/null || echo 0)
            ok "T3.3 query  $ds  ${sp}s  RAM=$(ramg $qram)  $(mbs $QB)  $QR rows"
            printf "T3.3,%s,query,%s,%s,%s,%s,none,,,,DONE,partial decode by pseudogenome coordinate. SPRING --decompress-range addresses by READ INDEX and measured 93-99% of a FULL decode (4.93s vs 5.28s on E.coli) vs ours at 8-17%; CRAM/genocat need a reference\n" \
                "$ds" "$sp" "$qram" "$QB" "$QR" >> "$CSV3"
        else
            err "T3.3 query failed on $ds"; debug_dump "$ds query" "$WD/q.log"
            printf "T3.3,%s,query,,,,,none,,,,FAILED,\n" "$ds" >> "$CSV3"
        fi
        rm -f "$WD/q.fa"
    done

    # ---- T6a + T6b: the six ----
    if [ -n "${SANITY_ONLY:-}" ]; then
        local keep=""
        for _e in $C3_BASE; do [ -s "$ARCH_DIR/${_e%%:*}.capsule" ] && keep="$keep $_e"; done
        C3_BASE="$(echo $keep | cut -d' ' -f1)"
        [ -z "$C3_BASE" ] && { err "SANITY_ONLY: no archive matches a Claim 3 baseline -- T6a/T6b not rehearsed"; }
        say "  SANITY_ONLY: baselines restricted to:${C3_BASE:- none}"
    fi
    n=$(echo $C3_BASE | wc -w)
    for entry in $C3_BASE; do
        i=$((i+1)); t0=$(date +%s)
        local DS="${entry%%:*}" REFNAME="${entry##*:}"
        say ""; say "──────────────────────────────────────────────────────────────────────"
        inf "PHASE 3  [$i/$n]  $DS   (baseline reference: $REFNAME)"
        say "──────────────────────────────────────────────────────────────────────"
        arc="$ARCH_DIR/$DS.capsule"
        [ -s "$arc" ] || { err "$DS: archive missing (phase 1 must run first) — skipping"
                           printf "T3.1,%s,export,,,,,SPAdes,,,,NO_ARCHIVE,\n" "$DS" >> "$CSV3"; continue; }
        if [ "$REFNAME" = chr20 ]; then ref="$REFS/chr20.fa"; else ref="$REFS/c3_${REFNAME}.fa"; fi
        src="$DATA_DIR/${DS}_1.fq"; [ -s "$src" ] || src="$DATA_DIR/${DS}_pooled.fq"

        # T3.1 export vs SPAdes
        mark "P3 $DS: export"
        step "T3.1  our export"
        d="$OUT_DIR/c3_$DS"; mkdir -p "$d"
        /usr/bin/time -v "$DEC" export "$arc" "$d/contigs.fa" >/dev/null 2>"$d/export.log"
        local OURS_EXP EXP_RAM EXP_B EXP_R
        read -r OURS_EXP EXP_RAM <<< "$(parse_time_v "$d/export.log")"
        EXP_B=$(stat -c%s "$d/contigs.fa" 2>/dev/null || echo 0)
        EXP_R=$(grep -c '^>' "$d/contigs.fa" 2>/dev/null || echo 0)
        if [ -s "$d/contigs.fa" ]; then ok "our export  ${OURS_EXP}s  RAM=$(ramg $EXP_RAM)  $(mbs $EXP_B)  $EXP_R contigs"
        else err "our export produced nothing for $DS"; debug_dump "$DS export" "$d/export.log"; fi

        step "T3.1  SPAdes de-novo (the slow baseline — minutes to hours)"
        if [ -x "$SPADES" ] && [ -s "$src" ]; then
            /usr/bin/time -v python3 "$SPADES" -s "$src" -o "$d/spades" -t "$NPROC" -m $(( $(free -g | awk '/^Mem:/{print $2}') - 8 )) \
                > "$d/spades.log" 2>&1
            local RC=$?
            if [ $RC -eq 0 ] && [ -s "$d/spades/contigs.fasta" ]; then
                local SP_T SP_RAM; read -r SP_T SP_RAM <<< "$(parse_time_v "$d/spades.log")"
                ok "SPAdes      ${SP_T}s   ->  speedup $(awk -v s=$SP_T -v o=$OURS_EXP 'BEGIN{printf "%.1fx",s/o}')"
                checkpoint "$DS  T3.1 export done -- ours ${OURS_EXP}s vs SPAdes ${SP_T}s"
                printf "T3.1,%s,export,%s,%s,%s,%s,SPAdes,%s,%s,%s,DONE,\n" "$DS" "$OURS_EXP" "$EXP_RAM" \
                    "$EXP_B" "$EXP_R" "$SP_T" "$SP_RAM" \
                    "$(awk -v s=$SP_T -v o=$OURS_EXP 'BEGIN{printf "%.1f",s/o}')" >> "$CSV3"
            else
                err "SPAdes did not complete on $DS (rc=$RC)"; debug_dump "$DS SPAdes" "$d/spades.log"
                printf "T3.1,%s,export,%s,%s,%s,%s,SPAdes,,,,BASELINE_DNF,SPAdes did not complete (rc=%s)\n" "$DS" "$OURS_EXP" "$EXP_RAM" "$EXP_B" "$EXP_R" "$RC" >> "$CSV3"
            fi
            rm -rf "$d/spades/K"* "$d/spades/tmp" 2>/dev/null
        else
            err "SPAdes or input missing for $DS"
            printf "T3.1,%s,export,%s,%s,%s,%s,SPAdes,,,,BASELINE_MISSING,\n" "$DS" "$OURS_EXP" "$EXP_RAM" "$EXP_B" "$EXP_R" >> "$CSV3"
        fi

        # T3.2 coverage vs bwa+samtools+mosdepth
        mark "P3 $DS: coverage"
        step "T3.2  our coverage"
        /usr/bin/time -v "$DEC" coverage "$arc" "$d/coverage.tsv" >/dev/null 2>"$d/coverage.log"
        local OURS_COV COV_RAM NROW COV_B
        read -r OURS_COV COV_RAM <<< "$(parse_time_v "$d/coverage.log")"
        NROW=$(wc -l < "$d/coverage.tsv" 2>/dev/null || echo 0)
        COV_B=$(stat -c%s "$d/coverage.tsv" 2>/dev/null || echo 0)
        ok "our coverage ${OURS_COV}s  RAM=$(ramg $COV_RAM)  $(mbs $COV_B)  ${NROW} rows"

        step "T3.2  bwa + samtools sort + mosdepth (the conventional route)"
        if [ -s "$ref.bwt" ] && [ -s "$src" ] && command -v mosdepth >/dev/null; then
            # Timed as ONE pipeline under time -v: the conventional route is
            # align+sort+index+depth, and quoting only one of those would flatter us.
            /usr/bin/time -v bash -c "bwa mem -t $NPROC '$ref' '$src' \
              | samtools sort -@ 4 -o '$d/aln.bam' - \
              && samtools index '$d/aln.bam' \
              && mosdepth -t 4 '$d/md' '$d/aln.bam'" > "$d/bwa.log" 2>&1
            local RC=$?
            if [ $RC -eq 0 ]; then
                local CV_T CV_RAM; read -r CV_T CV_RAM <<< "$(parse_time_v "$d/bwa.log")"
                ok "bwa+mosdepth ${CV_T}s  ->  speedup $(awk -v s=$CV_T -v o=$OURS_COV 'BEGIN{printf "%.1fx",s/o}')"
                checkpoint "$DS  T3.2 coverage done -- ours ${OURS_COV}s vs bwa+mosdepth ${CV_T}s"
                printf "T3.2,%s,coverage,%s,%s,%s,%s,bwa+samtools+mosdepth,%s,%s,%s,DONE,\n" "$DS" "$OURS_COV" "$COV_RAM" \
                    "$COV_B" "$NROW" "$CV_T" "$CV_RAM" \
                    "$(awk -v s=$CV_T -v o=$OURS_COV 'BEGIN{printf "%.1f",s/o}')" >> "$CSV3"
            else
                err "bwa/mosdepth baseline failed on $DS"; debug_dump "$DS bwa" "$d/bwa.log"
                printf "T3.2,%s,coverage,%s,%s,%s,%s,bwa+samtools+mosdepth,,,,BASELINE_FAILED,\n" "$DS" "$OURS_COV" "$COV_RAM" "$COV_B" "$NROW" >> "$CSV3"
            fi
            rm -f "$d/aln.bam" "$d/aln.bam.bai"    # BAMs are large, the timing is what we keep
        else
            err "bwa index / input / mosdepth missing for $DS — baseline skipped"
            printf "T3.2,%s,coverage,%s,%s,%s,%s,bwa+samtools+mosdepth,,,,BASELINE_MISSING,\n" "$DS" "$OURS_COV" "$COV_RAM" "$COV_B" "$NROW" >> "$CSV3"
        fi
        inf "[$i/$n] $DS done in $(( $(date +%s) - t0 ))s   (elapsed $(_el))"
    done
    # ═══ T3.4 — LOCUS RETRIEVAL FIDELITY ═══════════════════════════════════
    # T3.1-T3.3 measure what retrieval COSTS. This measures whether what comes
    # back is USABLE. A het locus is not one place in a compression-optimal
    # pseudogenome -- it is N parallel contigs -- so a COORDINATE returns one
    # haplotype with the variation gone, while content addressing resolves all
    # representatives and returns a real pileup. Both modes are run on the same
    # archive and the same GIAB sites, so nothing else can explain the gap.
    # Correctness only: no timing is recorded and the queries run in parallel.
    if [ -s "$ARCH_DIR/HG002.capsule" ] && [ -s "$REFS/chr20.fa" ]; then
        mark "T3.4 locus fidelity"
        banner "T3.4 — locus retrieval fidelity (coordinate vs content)"
        local t34d="$OUT_DIR/c3_T3.4"
        if bash "$HERE/scripts/run_locus_fidelity.sh" "$DEC" "$ARCH_DIR/HG002.capsule" \
              "$REFS/chr20.fa" HG002 "$t34d" "${T34_N:-100}" > "$t34d.log" 2>&1; then
            local row; row=$(grep -a '^T34,' "$t34d.log" | tail -1)
            if [ -n "${row:-}" ]; then
                printf "%s,DONE\n" "${row#T34,}" >> "$CSV8"
                ok "T3.4  $(grep -a 'content recovers both' "$t34d.log" | tail -1)"
                checkpoint "T3.4 locus fidelity done"
            else err "T3.4: no result row"; debug_dump "T3.4" "$t34d.log"
                 printf "HG002,,,,,,,,,,NO_RESULT\n" >> "$CSV8"; fi
        else err "T3.4 runner failed"; debug_dump "T3.4" "$t34d.log"
             printf "HG002,,,,,,,,,,FAILED\n" >> "$CSV8"; fi
    else
        inf "SKIP T3.4: needs the HG002 archive and chr20.fa"
    fi

    checkpoint "PHASE 3 COMPLETE"
    banner "PHASE 3 COMPLETE  ->  $CSV3"
}

# ═══════════════════════════════════════════════════════════════════════════
case ",$PHASES," in *,1,*) run_phase1;; esac
case ",$PHASES," in *,2,*) run_phase2;; esac
case ",$PHASES," in *,3,*) run_phase3;; esac
stop_heartbeat

banner "BENCHMARK 1 COMPLETE"
say "finished  : $(date '+%Y-%m-%d %H:%M:%S')"
say "total time: $(_el)"
say "datasets  : $N_OK ok, $N_FAIL failed${FAILED_LIST:+ ($FAILED_LIST)}"
say ""
say "  T1.1/T1.2 (Claim 1): $CSV1"
say "  T2.1  (Claim 2): $CSV2"
say "  T2.2  (Claim 2): $CSV4"
say "  T2.3  (Claim 2): $CSV5"
say "  T2.4  (Claim 2): $CSV6"
say "  T2.5  (Claim 2): $CSV7"
say "  T3.4  (Claim 3): $CSV8"
say "  T6    (Claim 3): $CSV3"
say "  archives kept  : $ARCH_DIR  ($(du -sh "$ARCH_DIR" 2>/dev/null | cut -f1))"
say "  full log       : $LOG"
say ""
say "── CLAIM 1 ──"; column -s, -t "$CSV1" 2>/dev/null | head -30 | tee -a "$LOG"
say ""; say "── CLAIM 2 ──"; column -s, -t "$CSV2" 2>/dev/null | tee -a "$LOG"
say ""; say "── CLAIM 2 — T2.2 coverage sweep ──"; column -s, -t "$CSV4" 2>/dev/null | tee -a "$LOG"
say ""; say "── CLAIM 2 — T2.3 het-indel ──";      column -s, -t "$CSV5" 2>/dev/null | tee -a "$LOG"
say ""; say "── CLAIM 2 — T2.4 multi-allelic ──";  column -s, -t "$CSV6" 2>/dev/null | tee -a "$LOG"
say ""; say "── CLAIM 2 — T2.5 tetraploid ──";     column -s, -t "$CSV7" 2>/dev/null | tee -a "$LOG"
say ""; say "── CLAIM 3 ──"; column -s, -t "$CSV3" 2>/dev/null | head -30 | tee -a "$LOG"
say ""; say "── CLAIM 3 — T3.4 locus retrieval fidelity ──"; column -s, -t "$CSV8" 2>/dev/null | tee -a "$LOG"
rm -rf "$WD"
[ "$N_FAIL" -eq 0 ] && exit 0 || exit 1
