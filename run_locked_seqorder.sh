#!/bin/bash
# Single-command runner for the locked sequence+order scope
# (see LOCKED_SEQORDER_SCOPE.md for what/why). Usage:
#   INPUT=your.fq bash run_locked_seqorder.sh
set -euo pipefail
INPUT="${INPUT:?set INPUT=path/to.fq}"
BEST="${BEST:-/tmp/best100}"
cd "$(dirname "$INPUT")" 2>/dev/null || true
WORK="$(mktemp -d)"
cd "$WORK"
cp "$OLDPWD/$(basename "$INPUT")" ./in.fq 2>/dev/null || cp "$INPUT" ./in.fq

rm -f literal.txt perm.u32 mem_triples.bin pos_delta.bin pos_strand.bin \
      mm_ref.bin mm_obs.bin mm_pos.bin mm_ctx3.bin n_reads.txt n_indices.bin

ASM_OUT=$(DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 "$BEST" in.fq 3 40 16 22 16 16 1 24 64 1 2>&1)
PG_LEN=$(echo "$ASM_OUT" | grep -oP 'PG_LEN \K[0-9]+')
MAIN_END=$(echo "$ASM_OUT" | grep -oP 'MEM main   : \K[0-9]+')

L1=$(ENCODE_ONLY=1 /tmp/seqpar literal.txt 12 12 2>&1 | grep -oP 'coded=\K[0-9]+')
L2=$(/tmp/permcoder perm.u32 2>&1 | grep -oP 'this coder     = \K[0-9]+')
L3=$([ -s mem_triples.bin ] && /tmp/refcoder mem_triples.bin "$PG_LEN" "$MAIN_END" 2>&1 | grep -oP 'this coder     \K[0-9]+' || echo 0)
L4a=$(xz -9 -c pos_delta.bin 2>/dev/null | wc -c)
L4b=$(xz -9 -c pos_strand.bin 2>/dev/null | wc -c)
L5=$([ -s mm_ref.bin ] && /tmp/mmcoder mm_ref.bin mm_obs.bin 2>&1 | grep -oP 'coded=\K[0-9]+' || echo 0)
L6a=$([ -s n_reads.txt ] && /tmp/seqpar n_reads.txt 1 1 2>&1 | grep -oP 'coded=\K[0-9]+' || echo 0)
L6b=$([ -s n_indices.bin ] && xz -9 -c n_indices.bin 2>/dev/null | wc -c || echo 0)

TOTAL=$((L1+L2+L3+L4a+L4b+L5+L6a+L6b))
echo "LAYER1_sequence=$L1 LAYER2_order=$L2 LAYER3_memref=$L3 LAYER4_pos=$((L4a+L4b)) LAYER5_mismatch=$L5 LAYER6_nreads=$((L6a+L6b))"
echo "LOCKED_SEQORDER_TOTAL=$TOTAL"

rm -rf "$WORK"
