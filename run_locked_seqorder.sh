#!/bin/bash
# Single-command runner for the locked sequence+order scope
# (see LOCKED_SEQORDER_SCOPE.md for what/why). Usage:
#   INPUT=your.fq bash run_locked_seqorder.sh
#
# STAGE 100 REWRITE: matches the simplified, byte-identical-verified encoder
# (100_locked_seqorder.cpp). Dropped perm.u32/orig2uid_ranks.bin/uidorder.bin
# (the whole rank/bucket indirection scheme, replaced by PgRC2's real,
# simpler direct-array design) and mem_triples.bin now carries a real is_rc
# flag (13 bytes/record, not 12) needed to fix a genuine reverse-complement
# reconstruction bug found while building the first full-pipeline decoder.
# Confirmed 0-diff byte-identical against the true original file on E. coli
# with this exact set of layers -- see decode_locked_seqorder.py.
set -euo pipefail
INPUT="${INPUT:?set INPUT=path/to.fq}"
BEST="${BEST:-/tmp/best100}"
cd "$(dirname "$INPUT")" 2>/dev/null || true
WORK="$(mktemp -d)"
cd "$WORK"
cp "$OLDPWD/$(basename "$INPUT")" ./in.fq 2>/dev/null || cp "$INPUT" ./in.fq

rm -f literal.txt mem_triples.bin pos_abs.bin pos_strand.bin read_lengths.bin \
      orig2uid.bin mm_ref.bin mm_obs.bin mm_pos.bin mm_ctx3.bin mm_count_per_read.bin \
      n_reads.txt n_indices.bin

ASM_OUT=$(DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 "$BEST" in.fq 3 40 16 22 16 16 1 24 64 1 2>&1)
PG_LEN=$(echo "$ASM_OUT" | grep -oP 'PG_LEN \K[0-9]+')
MAIN_END=$(echo "$ASM_OUT" | grep -oP 'MEM main   : \K[0-9]+')

L1=$(ENCODE_ONLY=1 /tmp/seqpar literal.txt 12 12 2>&1 | grep -oP 'coded=\K[0-9]+')
L3=$([ -s mem_triples.bin ] && /tmp/refcoder mem_triples.bin "$PG_LEN" "$MAIN_END" 2>&1 | grep -oP 'this coder     \K[0-9]+' || echo 0)
L4a=$(xz -9 -c pos_abs.bin 2>/dev/null | wc -c)
L4b=$(xz -9 -c pos_strand.bin 2>/dev/null | wc -c)
L5=$([ -s mm_ref.bin ] && /tmp/mmcoder mm_ref.bin mm_obs.bin 2>&1 | grep -oP 'coded=\K[0-9]+' || echo 0)
L6a=$([ -s n_reads.txt ] && /tmp/seqpar n_reads.txt 1 1 2>&1 | grep -oP 'coded=\K[0-9]+' || echo 0)
L6b=$([ -s n_indices.bin ] && xz -9 -c n_indices.bin 2>/dev/null | wc -c || echo 0)
L7=$([ -s read_lengths.bin ] && xz -9 -c read_lengths.bin 2>/dev/null | wc -c || echo 0)
L8=$([ -s orig2uid.bin ] && xz -9 -c orig2uid.bin 2>/dev/null | wc -c || echo 0)

TOTAL=$((L1+L3+L4a+L4b+L5+L6a+L6b+L7+L8))
echo "LAYER1_sequence=$L1 LAYER3_memref=$L3 LAYER4_pos=$((L4a+L4b)) LAYER5_mismatch=$L5 LAYER6_nreads=$((L6a+L6b)) LAYER7_lengths=$L7 LAYER8_orig2uid=$L8"
echo "LOCKED_SEQORDER_TOTAL=$TOTAL"

rm -rf "$WORK"
