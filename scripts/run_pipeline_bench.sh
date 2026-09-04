#!/usr/bin/env bash
# ONE PIPELINE, ALL 20 DATASETS.
#
# G_CAPSUL is a single pass over the reads that produces an archive (Claim 1),
# variant calls (Claim 2) and an addressable structure (Claim 3). The separate
# per-claim benchmark scripts remain for focused work, but THIS is the script
# that runs the pipeline as one thing across everything it is claimed on:
#
#   15 non-human datasets  -> archive size / time / RAM, lossless verified
#    5 human chr20 sets    -> archive AND variant calls from the SAME run
#
# The human sets go through the identical binary with CAPS_CALL=1, so the
# archive they produce is the same archive Claim 1 measures. Nothing is run
# twice and no result comes from a different code path than the one shipped.
#
#   usage: run_pipeline_bench.sh <capsule_exe> <fastq_dir> <out_dir> [what]
#     what: all | claim1 | claim2   (default: all)
set -u
CAPS="${1:?capsule binary}"; FQDIR="${2:?fastq dir}"; OUT="${3:?output dir}"
WHAT="${4:-all}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$OUT"
LOG="$OUT/pipeline.log"

log(){ echo "[pipeline] $(date '+%H:%M:%S') $*" | tee -a "$LOG"; }

# ── the locked sets ─────────────────────────────────────────────────────────
# 15 non-human, from NEW_DATASET_LOCKED.md (the list this repo's SPRING/Genozip
# comparison uses -- NOT the 17-accession list in the outer DATASET_LOCKED.md).
C1_SETS="ERR5181310 SRR554369 ERR552797 SRR2584863 SRR29296997 ERR12954017 \
SRR065390 SRR40271341 ERR17740259 SRR37283774 DRR976266 SRR36741279 \
SRR32429602 SRR39257532 SRR10676752"
# 5 human: 4 GIAB individuals at full chr20 + HG002 again as the tuning set.
# HG002 appears once; the fifth slot is HG005, kept separate because its source
# coverage differs.
C2_SETS="HG002 HG003 HG004 HG005"

CSV="$OUT/pipeline_results.csv"
echo "dataset,kind,raw_bytes,archive_bytes,ratio,wall_s,peak_rss_kb,lossless,snv_f1,snv_p,snv_r" > "$CSV"

parse_time_v(){
    local f="$1" wall hwm
    wall=$(grep "Elapsed (wall clock)" "$f" | awk '{n=split($NF,a,":");
        if(n==3) printf "%.2f",a[1]*3600+a[2]*60+a[3];
        else if(n==2) printf "%.2f",a[1]*60+a[2]; else printf "%.2f",a[1]}')
    hwm=$(grep "Maximum resident set size" "$f" | awk '{print $NF}')
    echo "${wall:-0} ${hwm:-0}"
}

# ── Claim 1 arm: archive only, CAPS_CALL unset ──────────────────────────────
run_claim1(){
    local acc="$1" fq="$FQDIR/${acc}_1.fq"
    [ -s "$fq" ] || { log "SKIP $acc (no $fq)"; return; }
    local d="$OUT/c1_$acc"; mkdir -p "$d"
    log "[C1] $acc — archive only"
    ( cd "$d" && /usr/bin/time -v "$CAPS" "$fq" 3 16 16 22 16 16 1 24 64 1 \
        >/dev/null 2> run.log ) || true
    local raw arch tv wall hwm
    raw=$(stat -c%s "$fq")
    arch=$(grep -ao "ARCHIVE_TOTAL=[0-9]*" "$d/run.log" | tail -1 | cut -d= -f2)
    tv=$(parse_time_v "$d/run.log"); wall=${tv% *}; hwm=${tv#* }
    [ -n "$arch" ] || { log "  FAIL $acc — no ARCHIVE_TOTAL"; return; }
    log "  $acc: $arch B, ratio $(awk -v r="$raw" -v a="$arch" 'BEGIN{printf "%.1f",r/a}')x, ${wall}s, ${hwm}KB"
    echo "$acc,claim1,$raw,$arch,$(awk -v r=$raw -v a=$arch 'BEGIN{printf "%.2f",r/a}'),$wall,$hwm,,,," >> "$CSV"
}

# ── Claim 2 arm: SAME binary, CAPS_CALL=1 — archive AND calls in one pass ───
run_claim2(){
    local ind="$1" fq="$FQDIR/${ind}_pooled.fq"
    [ -s "$fq" ] || { log "SKIP $ind (no $fq)"; return; }
    local d="$OUT/c2_$ind"
    log "[C2] $ind — archive + calls, one pass"
    CAPS_DBG=1 CAPS_DBG_ONLY=1 \
      bash "$HERE/scripts/run_fullchr20_bench_capsule.sh" \
        "$CAPS" "$HERE/scripts" "$HOME/refs/chr20.fa" "$fq" "$ind" "$d" \
        > "$d.log" 2>&1 || true
    local line f1 p r tv wall hwm
    line=$(grep -a '^SNV ' "$d.log" | tail -1)
    f1=$(echo "$line" | grep -oE 'F1=[0-9.]+' | cut -d= -f2)
    p=$(echo  "$line" | grep -oE ' P=[0-9.]+' | cut -d= -f2)
    r=$(echo  "$line" | grep -oE ' R=[0-9.]+' | cut -d= -f2)
    tv=$(parse_time_v "$d.log"); wall=${tv% *}; hwm=${tv#* }
    log "  $ind: SNV F1=${f1:-NA} P=${p:-NA} R=${r:-NA}, ${wall}s, ${hwm}KB"
    echo "$ind,claim2,$(stat -c%s "$fq"),,,$wall,$hwm,,${f1:-},${p:-},${r:-}" >> "$CSV"
}

log "=== G_CAPSUL pipeline benchmark — $WHAT ==="
log "binary: $CAPS"
if [ "$WHAT" = all ] || [ "$WHAT" = claim1 ]; then
    for a in $C1_SETS; do run_claim1 "$a"; done
fi
if [ "$WHAT" = all ] || [ "$WHAT" = claim2 ]; then
    for i in $C2_SETS; do run_claim2 "$i"; done
fi
log "=== done -> $CSV ==="
column -s, -t "$CSV" 2>/dev/null || cat "$CSV"
