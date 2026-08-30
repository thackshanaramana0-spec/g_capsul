#!/bin/bash
# STAGE 89 -- combined pipeline, back to SEQUENTIAL execution (revert of
# the stage 76+ task-level overlap trick, per direct user instruction
# after stage 87 uncovered a real race in that overlap mechanism: it read
# literal.txt/perm.u32 while `best` might still be mid-write, confirmed to
# fail 3/3 times under direct testing, and shown to have silently inflated
# the earlier "beats SPRING" speed/RAM numbers -- seqpar finishing fast
# because it sometimes read truncated data, not because it was actually
# faster).
#
# Safe by construction, not by patching around a race: assembly runs to
# full completion first, THEN each downstream tool runs, one at a time.
# Nothing ever reads a file while another process might still be writing
# it, so there is no race to close in the first place.
#
# Keeps every other real, verified improvement from today's work: stage
# 87's atomically-safe writes (good defensive practice regardless), stage
# 74/75's sequence-coder speed fix (TBITS + block-parallelism), stage
# 84/85's real bounded-queue RAM fix for names/quality. Only the risky
# overlap SCHEDULING is reverted, not the underlying per-tool work.
#
# Usage: same env vars as stage 88.
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

FWD_SELF=0 DUMP_LIT=1 DUMP_PERM=1 "$BEST" "$FASTQ" 3 40 16 22 16 16 1 24 64 1 > /dev/null 2>&1
ENCODE_ONLY=1 "$SEQPAR" literal.txt 12 12 > /dev/null 2>&1
"$PERMCODER" perm.u32 > /dev/null 2>&1
"$REFCODER" mem_triples.bin 23233953 21619175 > /dev/null 2>&1
"$NAMEBQ" "$NAMES" 100000 12 > /dev/null 2>&1
"$QUALBQ" "$QUAL" 100000 12 > /dev/null 2>&1
