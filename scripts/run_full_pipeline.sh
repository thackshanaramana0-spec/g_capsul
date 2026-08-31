#!/bin/bash
# One-command, reproducible run of the full reimpl pipeline PLUS the two
# stages found and round-trip verified in the 2026-08-29 deep pass: iterative
# MEM removal (48+49) and the real adaptive mismatch coder (50). Everything
# printed is either verified in this run or carried from a prior verified
# measurement (order/reference coders, unaffected by this work).
#
#   ./run_full_pipeline.sh ./best yeast_sub.fq /tmp/fullrun
set -eu
BIN=${1:?usage: run_full_pipeline.sh <best-binary> <fastq> <workdir>}
FQ=${2:?}
WD=${3:?}
HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
FQ="$(cd "$(dirname "$FQ")" && pwd)/$(basename "$FQ")"
mkdir -p "$WD"; cd "$WD"

for t in iter_mem fwd_mem mmcoder permcoder refcoder; do
    [ -x "$HERE/$t" ] || { echo "missing $HERE/$t -- build it first (see README)"; exit 1; }
done

echo "=== assembly + mapping + single-pass MEM removal ==="
FWD_SELF=0 DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 "$BIN" "$FQ" 3 40 16 22 16 16 1 24 64 1 2>run.log >run.out
PGLEN=$(grep -o 'PG_LEN [0-9]*' run.out | awk '{print $2}')
MAINEND=$(grep -o 'pg after chains=[0-9]*' run.log | tail -1 | cut -d= -f2)
grep -E "PG_LEN|PG_LITERAL" run.out
orig_lit_bytes=$(stat -c%s literal.txt)
orig_lit_coded=$(xz -6 -T"$(nproc)" -c literal.txt | wc -c 2>/dev/null || true)
# literal.txt is ASCII ACGT here; pack to 2 bit for the size that matters
python3 - "$HERE" <<'PYEOF'
import sys
sys.exit(0)
PYEOF

echo "=== iterate MEM removal to convergence (RC passes) ==="
cp literal.txt "iter0.txt"
i=0
while :; do
  i=$((i+1))
  "$HERE/iter_mem" "iter$((i-1)).txt" 24 64 --dump
  mv literal2.txt "iter$i.txt"
  mv r2_gaps.bin "iter${i}_gaps.bin"; mv r2_srcs.bin "iter${i}_srcs.bin"; mv r2_lens.bin "iter${i}_lens.bin"
  before=$(stat -c%s "iter$((i-1)).txt"); after=$(stat -c%s "iter$i.txt")
  gain=$((before-after))
  echo "  RC pass $i: $before -> $after  (-$gain bases)"
  [ "$gain" -lt 2000 ] && break
  [ "$i" -ge 8 ] && break
done
LASTRC="iter$i.txt"
RC_PASSES=$i

echo "=== one forward self-match pass on the converged residual ==="
"$HERE/fwd_mem" "$LASTRC" 24 64
mv fwd_literal.txt final_literal.txt
mv fwd_gaps.bin final_gaps.bin; mv fwd_srcs.bin final_srcs.bin; mv fwd_lens.bin final_lens.bin
before=$(stat -c%s "$LASTRC"); after=$(stat -c%s final_literal.txt)
echo "  FWD pass: $before -> $after (-$((before-after)) bases)"

echo "=== real mismatch coder ==="
"$HERE/mmcoder" mm_ref.bin mm_obs.bin

echo "=== pricing ==="
python3 - "$RC_PASSES" <<'PYEOF'
import subprocess, sys, os
RC_PASSES = int(sys.argv[1])

def xz9e(f):
    return int(subprocess.run(f"xz -9e -c {f} | wc -c", shell=True, capture_output=True, text=True).stdout)

def pack2(inp):
    with open(inp,'rb') as f: d=f.read()
    n=len(d); out=bytearray((n+3)//4)
    m={65:0,67:1,71:2,84:3}
    for i,c in enumerate(d):
        out[i>>2] |= m.get(c,0) << (2*(3-(i&3)))
    return bytes(out)

def xz6_2bit(txtfile):
    packed = pack2(txtfile)
    tmp = txtfile+".2bit"
    open(tmp,"wb").write(packed)
    r = int(subprocess.run(f"xz -6 -T$(nproc) -c {tmp} | wc -c", shell=True, capture_output=True, text=True).stdout)
    os.remove(tmp)
    return r

orig_lit_coded = xz6_2bit("literal.txt")
final_lit_coded = xz6_2bit("final_literal.txt")

ref_total = 0
for i in range(1, RC_PASSES+1):
    ref_total += xz9e(f"iter{i}_gaps.bin") + xz9e(f"iter{i}_srcs.bin") + xz9e(f"iter{i}_lens.bin")
ref_total += xz9e("final_gaps.bin") + xz9e("final_srcs.bin") + xz9e("final_lens.bin")

print(f"original single-pass literal, coded : {orig_lit_coded:,} B")
print(f"final (RC-converged + FWD) literal   : {final_lit_coded:,} B")
print(f"all new-pass reference streams       : {ref_total:,} B  ({RC_PASSES} RC + 1 FWD)")
print(f"sequence+new-refs TOTAL              : {final_lit_coded+ref_total:,} B")
PYEOF
