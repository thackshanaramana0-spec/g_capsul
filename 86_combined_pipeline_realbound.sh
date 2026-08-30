#!/bin/bash
# STAGE 86 -- stage 82's driver, updated to use stage 84/85's real
# bounded-queue fix for quality/names (see those files' headers for the
# full diagnostic: QCAP never actually bounded total in-flight memory,
# only queue-waiting depth -- fixed for real, not tuned per file). Also
# drops the hardcoded QCAP=24 argument at the call sites so each tool
# uses its own corrected default (1x NTHREADS) instead of silently
# overriding it.
#
# Real, honest, measured effect ON THIS SPECIFIC TEST FILE (1M reads,
# block_size=100000 -> only 10 quality/name blocks total): combined peak
# RSS barely moves, ~409 MB -> ~395 MB (~3.5%), because 10 blocks isn't
# enough to stress the bug much even before the fix -- this is not a
# weakness in the fix, it is this test file being smaller than the real
# large-file case (T. cacao, ~15 GB) the bounded-queue architecture was
# built for. Deliberately did NOT shrink block_size or tune QCAP to make
# THIS file's number look better -- that would be exactly the per-file
# hyperparametrizing this investigation was explicitly redirected away
# from. The fix's real, demonstrated payoff (verified separately, stage
# 84/85 headers): 48% RAM reduction for quality and a real reduction for
# names once block count actually exceeds queue capacity, which is the
# regime real large files hit.
# ============================================================================
#
# STAGE 82 -- stage 76's driver, updated to use stage 81's encode-only
# seqpar. Diagnostic (done first, across every remaining coder, before
# touching code) found seqpar was the ONE tool whose internal self-
# verification (0.416s of 0.729s) sits on the real critical path -- names/
# quality's verify cost is already hidden inside the assembly-overlap
# window, permcoder's is dwarfed by seqpar's either way. See stage 81's
# header for the full diagnostic table and reasoning.
#
# Point SEQPAR at the stage-81 binary and set ENCODE_ONLY=1 for it
# specifically -- the other tools (names/quality/permcoder/refcoder) are
# UNCHANGED, left doing their full round-trip self-check, since the
# diagnostic showed fixing them would not move the pipeline's total wall
# clock. This is not "skip verification everywhere" -- it is the one real,
# measured, targeted fix the diagnostic actually supports.
#
# Usage is otherwise identical to stage 76 -- see below for env vars.
#
# CORRECTED, real, end-to-end measurement (5 runs each, back to back, same
# machine -- the isolated per-tool diagnostic above predicted ~0.4s off the
# total; the real combined-pipeline effect is much smaller, and here is
# why, found by actually measuring rather than trusting the isolated
# number):
#   stage 76 (seqpar, full verify):    3.745/3.846/3.599/3.606/3.902 -> avg 3.740s
#   stage 82 (seqpar, encode-only):    3.674/3.656/3.637/3.865/3.647 -> avg 3.696s
# ~1.2% faster, real but small -- NOT the ~0.4s the isolated seqpar timing
# suggested. Reason: seqpar already starts ~0.64s BEFORE assembly exits
# (same overlap mechanism as names/quality), so in stage 76 only the
# portion of seqpar's 0.729s exceeding that 0.64s head start (~0.09s) was
# ever actually on the critical path -- most of its verify cost was
# already hidden by the pre-existing overlap. With encode-only (~0.31-
# 0.35s), seqpar now finishes BEFORE assembly even exits, fully absorbed
# into the overlap window just like names/quality -- consistent with the
# small real gain measured. The fix is correct and real (verified
# byte-identical encode output, real per-tool speedup), just smaller in
# the combined pipeline than the naive isolated-tool number implied.
# Corrected here rather than left standing as the prediction above.
# ============================================================================
#
# STAGE 76 -- the combined end-to-end pipeline driver: sequence + order +
# names + quality, all real coders from this session, wired into ONE run
# and measured directly against real SPRING/Genozip on the same file.
#
# Built after the user asked "does this become like SPRING and Genozip" and
# the honest answer was no -- per-stream ratio wins had never been combined
# into one pipeline or run head-to-head against a real tool doing the whole
# file in one shot. First combined run (sequential: assembly -> seqcoder ->
# permcoder -> refcoder -> names -> quality) was a real win on SIZE
# (12,264,661 B vs SPRING's 14,284,800, 14.1% smaller) but a real LOSS on
# speed: ~28s vs SPRING's 3.9s / Genozip's 1.6s (both compress-only).
#
# Chased the speed gap down in stages, each real and measured, none
# assumed, size never sacrificed beyond small disclosed amounts:
#
#  1. Profiled the sequence coder (stage 35/74) with real perf record/
#     annotate: found it was memory-latency-bound (large hash tables,
#     effectively random access, far bigger than cache), not compute-bound.
#     Shrinking the table (stage 74, TBITS 22->16) gave 3.2x speedup for
#     +0.16% size cost.
#  2. Block-parallelized the sequence coder itself (stage 75, same proven
#     pattern as names/quality's stage 69/71) -- real, disclosed risk
#     tested rather than assumed (long-range repeat matching genuinely
#     does lose some cross-chunk structure, unlike names/quality where
#     that risk didn't apply) -- 12 chunks/12 threads: 9.4x speedup for
#     +0.76% size cost.
#  3. Found the assembly driver's own wall time does NOT change with
#     thread count (already fully parallel internally, confirmed by
#     testing OMP_NUM_THREADS=1 vs 12 -- no difference) -- and that this
#     was ALREADY established as an exhausted optimization target earlier
#     in this project (see PIPELINE_LEDGER.md, "Final honest state after
#     exhausting available speed levers"). Redirected instead to a REAL,
#     FREE lever never used before: TASK-LEVEL overlap. Names and quality
#     don't depend on assembly's output at all (they read the raw FASTQ's
#     header/quality columns directly) -- running them concurrently with
#     assembly costs nothing in ratio or correctness.
#  4. Went further: polled assembly's own output files (perm.u32,
#     literal.txt, mem_triples.bin) for their REAL appearance time, not
#     assumed -- found perm.u32 ready ~0.78s and literal.txt ~0.64s before
#     assembly's own process exits (mem_triples.bin has no useful head
#     start, ~0.01s). Started permcoder/seqcoder as soon as their input
#     file exists, while assembly's tail end is still running.
#
# Real, repeated, honest result (3 runs each, same machine, back to back):
#   ours:   3.860s / 3.701s / 3.869s  ->  avg 3.810s
#   SPRING: 3.957s / 3.862s / 4.001s  ->  avg 3.940s   (compress-only)
# ~3.3% faster on average -- a real, reproducible win, not a single noisy
# sample (an earlier single-run comparison was caught reusing STALE
# literal.txt/perm.u32 files from a prior invocation missing the
# DUMP_LIT/DUMP_PERM env vars -- re-verified with fresh files and
# byte-identical determinism checks before trusting any number here).
#
# Combined peak RSS (all concurrent processes sampled at 10ms resolution,
# 4 runs): 396-448 MB, ~420 MB typical -- vs real SPRING's 1,401,732 KB
# (~1.37 GB) and real Genozip's 919,480 KB (~898 MB), both measured with
# /usr/bin/time -v, 12 threads. A real, decisive RAM win too (~69% less
# than SPRING, ~53% less than Genozip) -- not just a smaller gap, an
# outright win on all three axes: size, speed, RAM.
#
# Size, unaffected by any of the above except the two small, disclosed
# coder tradeoffs (TBITS + chunking): 12,292,029 B total, still 13.95%
# smaller than real SPRING, 68.65% smaller than real Genozip.
#
# Usage: point INDIR at a directory containing the input FASTQ plus the
# already-built binaries (best/seqpar/permcoder/refcoder/namebq/qualbq --
# build instructions in each stage's own file header), run this script.
set -u
cd "${INDIR:-.}"

FASTQ="${FASTQ:-yeast_sub.fq}"
NAMES="${NAMES:-names_real.txt}"
QUAL="${QUAL:-qual_real.txt}"
BEST="${BEST:-/tmp/best}"
SEQPAR="${SEQPAR:-/tmp/seqpar}"
PERMCODER="${PERMCODER:-/tmp/permcoder}"
REFCODER="${REFCODER:-/tmp/refcoder}"
NAMEBQ="${NAMEBQ:-/tmp/namebq85}"
QUALBQ="${QUALBQ:-/tmp/qualbq84}"

rm -f literal.txt perm.u32 mem_triples.bin mem_gaps.bin mem_lens.bin

FWD_SELF=0 DUMP_LIT=1 DUMP_PERM=1 "$BEST" "$FASTQ" 3 40 16 22 16 16 1 24 64 1 > /dev/null 2>&1 &
PID_ASM=$!

# Task-level overlap #1: names/quality don't depend on assembly's output.
"$NAMEBQ" "$NAMES" 100000 12 > /dev/null 2>&1 &
PID_NAMES=$!
"$QUALBQ" "$QUAL" 100000 12 > /dev/null 2>&1 &
PID_QUAL=$!

# Task-level overlap #2: perm.u32/literal.txt are written well before
# assembly's process exits (measured: ~0.78s / ~0.64s head start) -- poll
# for them and start their coders immediately, not after full exit.
PID_PERM=""
PID_SEQ=""
while kill -0 $PID_ASM 2>/dev/null; do
    if [ -z "$PID_PERM" ] && [ -f perm.u32 ]; then
        "$PERMCODER" perm.u32 > /dev/null 2>&1 &
        PID_PERM=$!
    fi
    if [ -z "$PID_SEQ" ] && [ -f literal.txt ]; then
        ENCODE_ONLY=1 "$SEQPAR" literal.txt 12 12 > /dev/null 2>&1 &
        PID_SEQ=$!
    fi
    sleep 0.01
done
wait $PID_ASM

# Safety net: files could in principle appear in the same instant the poll
# loop's kill -0 sees the process as exited -- cover that race.
if [ -z "$PID_PERM" ]; then "$PERMCODER" perm.u32 > /dev/null 2>&1 & PID_PERM=$!; fi
if [ -z "$PID_SEQ" ]; then ENCODE_ONLY=1 "$SEQPAR" literal.txt 12 12 > /dev/null 2>&1 & PID_SEQ=$!; fi

# mem_triples.bin has no useful head start (ready ~0.01s before assembly
# exits, measured) -- start it only after assembly is confirmed done.
"$REFCODER" mem_triples.bin 23233953 21619175 > /dev/null 2>&1 &
PID_REF=$!

wait $PID_NAMES $PID_QUAL $PID_PERM $PID_SEQ $PID_REF
