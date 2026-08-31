#!/bin/bash
# STAGE 91 -- resolves the reorder-vs-no-reorder loss found on SARS-CoV-2
# (SCOPE_AND_PLAN.md item 6/7). Root cause: the pipeline always reorders and
# always pays the full log2(n!) permutation cost, even when reordering
# doesn't reduce sequence size enough to justify it (confirmed on
# SARS-CoV-2: unreordered raw-sequence encode is 69,993 B vs the reorder
# path's 227,766 B seq+perm+mem total -- 3.25x smaller without reordering).
# Confirmed the fix must be a real MEASUREMENT, not a per-file-tuned
# threshold, by testing the opposite case: on P. aeruginosa, raw-unreordered
# encode is 2,394,279 B vs the reorder path's 1,900,348 B -- reordering
# genuinely wins there, so a fixed "always skip reorder" or "always reorder"
# rule would be wrong in one direction or the other. This decides by running
# BOTH and keeping the actual smaller total -- no threshold, no per-file
# constant, generalizes by construction.
#
# No-reorder path is simpler than the reorder path: it turns out our own
# seqpar coder's context model already exploits cross-read redundancy on its
# own (0.025 bits/base on SARS-CoV-2's raw, unreordered sequence) -- no
# separate MEM-matching or permutation stream is needed. The raw per-read
# text (one read per line, original order, newlines included) round-trips
# through seqpar exactly (verified: newlines are part of the coded stream,
# so decode reconstructs read boundaries with no extra index needed).
#
# Usage: same env vars as stage 89 (sequential pipeline), point BEST at a
# stage-90 (or later) binary -- the read-length fix is required for correct
# behavior on datasets with >255bp reads regardless of this stage's change.
set -u
cd "${INDIR:-.}"

FASTQ="${FASTQ:-yeast_sub.fq}"
NAMES="${NAMES:-names_real.txt}"
QUAL="${QUAL:-qual_real.txt}"
BEST="${BEST:-/tmp/best90}"
SEQPAR="${SEQPAR:-/tmp/seqpar}"
PERMCODER="${PERMCODER:-/tmp/permcoder}"
REFCODER="${REFCODER:-/tmp/refcoder}"
NAMEBQ="${NAMEBQ:-/tmp/namebq85}"
QUALBQ="${QUALBQ:-/tmp/qualbq84}"

rm -f literal.txt perm.u32 mem_triples.bin mem_gaps.bin mem_lens.bin literal_noreorder.txt

# ── candidate A: reorder path (unchanged) ───────────────────────────────
FWD_SELF=0 DUMP_LIT=1 DUMP_PERM=1 "$BEST" "$FASTQ" 3 40 16 22 16 16 1 24 64 1 > asm_stdout.log 2>&1
PG_LEN=$(grep -oP 'PG_LEN \K[0-9]+' asm_stdout.log)
MAIN_END=$(grep -oP 'MEM main   : \K[0-9]+' asm_stdout.log)

SEQ_A=$("$SEQPAR" literal.txt 12 12 2>&1 | grep -oP 'coded=\K[0-9]+')
PERM_A=$("$PERMCODER" perm.u32 2>&1 | grep -oP 'this coder     = \K[0-9]+')
REF_A=$([ -s mem_triples.bin ] && "$REFCODER" mem_triples.bin "$PG_LEN" "$MAIN_END" 2>&1 | grep -oP 'this coder     \K[0-9]+' || echo 0)
TOTAL_A=$((SEQ_A + PERM_A + REF_A))

# ── candidate B: no-reorder path (raw sequence, original order) ────────
awk 'NR%4==2' "$FASTQ" > literal_noreorder.txt
SEQ_B=$("$SEQPAR" literal_noreorder.txt 12 12 2>&1 | grep -oP 'coded=\K[0-9]+')
TOTAL_B=$SEQ_B

# ── decide by actual measurement, not a threshold ───────────────────────
if [ "$TOTAL_B" -lt "$TOTAL_A" ]; then
    MODE="noreorder"
    SEQ_TOTAL=$TOTAL_B
    rm -f literal.txt perm.u32 mem_triples.bin   # unused streams for this file
else
    MODE="reorder"
    SEQ_TOTAL=$TOTAL_A
    rm -f literal_noreorder.txt
fi

"$NAMEBQ" "$NAMES" 100000 12 > names.log 2>&1
"$QUALBQ" "$QUAL" 100000 12 > qual.log 2>&1
NAMES_B=$(grep -oP 'coded=\K[0-9]+' names.log)
QUAL_B=$(grep -oP 'coded=\K[0-9]+' qual.log)

TOTAL=$((SEQ_TOTAL + NAMES_B + QUAL_B))
echo "MODE=$MODE  seq_total=$SEQ_TOTAL (reorder would be $TOTAL_A, no-reorder would be $TOTAL_B)  names=$NAMES_B  qual=$QUAL_B  GRAND_TOTAL=$TOTAL"
