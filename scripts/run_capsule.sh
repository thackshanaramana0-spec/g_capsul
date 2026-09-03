#!/bin/bash
# CAPSULE — one entry point, three claims. Written for a reviewer who wants
# to type one thing per claim, not read the docs first.
#
#   bash scripts/run_capsule.sh <1|2|3> [phase]
#
# Claims are INDEPENDENT, not sequential -- there is no requirement to run
# 1 before 3. Each claim builds whatever it needs from scratch (its own
# archive, its own binaries). Running `run_capsule.sh 3` alone is a valid,
# complete thing to do. Each invocation prints which locked dataset it used
# and where that path comes from, so nothing is silently assumed.
#
# No dataset path is hardcoded in this file or the phase scripts it calls --
# every one of them reads from scripts/capsule_config.sh. To point at a
# different data location, either edit that one file, or override for one
# run:
#   CAPSULE_DATA_DIR=/mnt/other/fastq bash scripts/run_capsule.sh 1
#
# Phases, per claim (default phase is marked *):
#   1  verify    * -- one locked dataset, encode+decode+diff, real lossless
#                     proof (fast, default -- always runnable with no setup
#                     beyond the dataset itself)
#      compress    -- verify, plus keeps the archive for inspection
#      sweep       -- the real 14-dataset SPRING/Genozip comparison. NOT YET
#                     wrapped as one command (docs/INDUSTRIAL_CHECKLIST_CLAIM1.md
#                     names this gap) -- this phase prints where to find the
#                     manual steps rather than pretending to automate them
#   2  test      * -- synthetic regression test, no downloads needed (fast, default)
#      giab       -- real GIAB het-SNV+indel benchmark (needs chr20.fa + truth VCFs)
#      window     -- single-window benchmark for one individual
#   3  test      * -- synthetic regression test, no downloads needed (fast, default)
#      full       -- the real export/coverage/query benchmark (builds SPAdes if absent)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$HERE/scripts/capsule_config.sh"

CLAIM="${1:?Usage: $0 <1|2|3> [phase]}"
PHASE="${2:-}"
mkdir -p "$CAPSULE_BIN_DIR" "$CAPSULE_OUT_DIR"

log() { echo "[capsule] $*"; }

case "$CLAIM" in

1)
    PHASE="${PHASE:-verify}"
    log "Claim 1 (COMPACT), phase: $PHASE"
    log "Dataset dir: \$CAPSULE_DATA_DIR = $CAPSULE_DATA_DIR (edit scripts/capsule_config.sh to change)"
    case "$PHASE" in
        verify)
            FQ="${3:-$CAPSULE_DATA_DIR/SRR2584863_1.fq}"
            [ -s "$FQ" ] || { log "FAIL: $FQ not found -- see docs/SERVER_SETUP_AND_DOWNLOADS.md sec 2"; exit 1; }
            bash "$HERE/scripts/verify_lossless.sh" "$FQ"
            ;;
        compress)
            FQ="${3:-$CAPSULE_DATA_DIR/SRR2584863_1.fq}"
            [ -s "$FQ" ] || { log "FAIL: $FQ not found -- see docs/SERVER_SETUP_AND_DOWNLOADS.md sec 2"; exit 1; }
            BEST="$CAPSULE_BIN_DIR/best106"; [ -x "$BEST" ] || bash "$HERE/scripts/build106.sh" "$BEST" >/dev/null
            INPUT="$FQ" ARCHIVE="$CAPSULE_OUT_DIR/out.capsule" BEST="$BEST" bash "$HERE/scripts/encode_adaptive.sh"
            log "Archive: $CAPSULE_OUT_DIR/out.capsule ($(stat -c%s "$CAPSULE_OUT_DIR/out.capsule") bytes)"
            ;;
        sweep)
            log "This is the real 14-dataset SPRING/Genozip comparison -- slow, real data required."
            log "Not yet wrapped as a single script (docs/INDUSTRIAL_CHECKLIST_CLAIM1.md notes this gap)."
            log "Run manually per docs/SOTA_COMPARISON.md's regeneration section, or use 'compress'/'verify' for a fast single-dataset check."
            exit 2
            ;;
        *) log "Unknown Claim 1 phase: $PHASE (try: compress, verify, sweep)"; exit 2 ;;
    esac
    ;;

2)
    PHASE="${PHASE:-test}"
    log "Claim 2 (FAITHFUL), phase: $PHASE"
    case "$PHASE" in
        test)
            bash "$HERE/scripts/test_claim2.sh"
            ;;
        giab)
            log "Reference dir: \$CAPSULE_REFS_DIR = $CAPSULE_REFS_DIR"
            REF="${3:-$CAPSULE_REFS_DIR/chr20.fa}"
            FQ="${4:-$CAPSULE_DATA_DIR/HG002_pooled.fq}"
            [ -s "$REF" ] || { log "FAIL: $REF not found -- see docs/SERVER_SETUP_AND_DOWNLOADS.md sec 3"; exit 1; }
            [ -s "$FQ" ]  || { log "FAIL: $FQ not found -- see docs/SERVER_SETUP_AND_DOWNLOADS.md sec 3"; exit 1; }
            BEST="$CAPSULE_BIN_DIR/best106"; [ -x "$BEST" ] || bash "$HERE/scripts/build106.sh" "$BEST" >/dev/null
            bash "$HERE/scripts/run_giab_indel_capsule.sh" "$BEST" "$HERE/scripts" "$FQ" "$REF"
            ;;
        window)
            IND="${3:-HG002}"; WIN="${4:-r2}"
            REF="$CAPSULE_REFS_DIR/chr20.fa"
            [ -s "$REF" ] || { log "FAIL: $REF not found -- see docs/SERVER_SETUP_AND_DOWNLOADS.md sec 3"; exit 1; }
            BEST="$CAPSULE_BIN_DIR/best106"; [ -x "$BEST" ] || bash "$HERE/scripts/build106.sh" "$BEST" >/dev/null
            bash "$HERE/scripts/run_window_bench_capsule.sh" "$BEST" "$HERE/scripts" "$REF" "$IND" "$WIN"
            ;;
        *) log "Unknown Claim 2 phase: $PHASE (try: test, giab, window)"; exit 2 ;;
    esac
    ;;

3)
    PHASE="${PHASE:-test}"
    log "Claim 3 (ADDRESSABLE), phase: $PHASE"
    case "$PHASE" in
        test)
            bash "$HERE/scripts/test_claim3.sh"
            ;;
        full)
            log "Dataset dir: \$CAPSULE_DATA_DIR = $CAPSULE_DATA_DIR"
            bash "$HERE/scripts/run_claim3.sh" "$CAPSULE_DATA_DIR" "$CAPSULE_OUT_DIR/claim3"
            ;;
        *) log "Unknown Claim 3 phase: $PHASE (try: test, full)"; exit 2 ;;
    esac
    ;;

*)
    echo "Usage: $0 <1|2|3> [phase]"
    echo "  1: verify | compress | sweep  (default: verify)"
    echo "  2: test | giab | window       (default: test)"
    echo "  3: test | full                (default: test)"
    exit 2
    ;;
esac
