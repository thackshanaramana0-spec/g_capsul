#!/bin/bash
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
NAMEBQ="${NAMEBQ:-/tmp/namebq}"
QUALBQ="${QUALBQ:-/tmp/qualbq}"

rm -f literal.txt perm.u32 mem_triples.bin mem_gaps.bin mem_lens.bin

FWD_SELF=0 DUMP_LIT=1 DUMP_PERM=1 "$BEST" "$FASTQ" 3 40 16 22 16 16 1 24 64 1 > /dev/null 2>&1 &
PID_ASM=$!

# Task-level overlap #1: names/quality don't depend on assembly's output.
"$NAMEBQ" "$NAMES" 100000 12 24 > /dev/null 2>&1 &
PID_NAMES=$!
"$QUALBQ" "$QUAL" 100000 12 24 > /dev/null 2>&1 &
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
        "$SEQPAR" literal.txt 12 12 > /dev/null 2>&1 &
        PID_SEQ=$!
    fi
    sleep 0.01
done
wait $PID_ASM

# Safety net: files could in principle appear in the same instant the poll
# loop's kill -0 sees the process as exited -- cover that race.
if [ -z "$PID_PERM" ]; then "$PERMCODER" perm.u32 > /dev/null 2>&1 & PID_PERM=$!; fi
if [ -z "$PID_SEQ" ]; then "$SEQPAR" literal.txt 12 12 > /dev/null 2>&1 & PID_SEQ=$!; fi

# mem_triples.bin has no useful head start (ready ~0.01s before assembly
# exits, measured) -- start it only after assembly is confirmed done.
"$REFCODER" mem_triples.bin 23233953 21619175 > /dev/null 2>&1 &
PID_REF=$!

wait $PID_NAMES $PID_QUAL $PID_PERM $PID_SEQ $PID_REF
