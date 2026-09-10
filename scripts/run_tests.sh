#!/bin/bash
# ============================================================================
# G_CAPSUL test suite. Self-contained: generates its own inputs, builds, and
# asserts. No external data, no network, no locked datasets required.
#
#   bash scripts/run_tests.sh          # all tests
#
# Exit 0 = every test passed. Exit 1 = at least one failed, with which.
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
PASS=0; FAIL=0
ok(){ printf "  PASS  %s\n" "$1"; PASS=$((PASS+1)); }
no(){ printf "  FAIL  %s -- %s\n" "$1" "${2:-}"; FAIL=$((FAIL+1)); }

echo "G_CAPSUL test suite"
echo "==================="

echo "[build]"
bash "$HERE/scripts/build106.sh"     "$W/enc" >"$W/b1.log" 2>&1 && ok "encoder builds" || { no "encoder builds" "see $W/b1.log"; echo "cannot continue"; exit 1; }
bash "$HERE/scripts/build_decode.sh" "$W/dec" >"$W/b2.log" 2>&1 && ok "decoder builds" || { no "decoder builds" "see $W/b2.log"; echo "cannot continue"; exit 1; }

# --- deterministic synthetic FASTQ -----------------------------------------
gen(){ # gen <reads> <len> <file> [seed]
  python3 - "$1" "$2" "$3" "${4:-1}" <<'PY'
import random,sys
n,L,path,seed=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3],int(sys.argv[4])
random.seed(seed)
G=''.join(random.choice('ACGT') for _ in range(max(L*4, n*L//8 + L)))
with open(path,'w') as f:
    for i in range(n):
        p=random.randrange(0,max(1,len(G)-L)); s=list(G[p:p+L])
        if random.random()<0.3: s[random.randrange(len(s))]=random.choice('ACGT')
        f.write(f"@read{i} synthetic\n{''.join(s)}\n+\n{'I'*len(s)}\n")
PY
}
enc(){ env CAPS_CALL=1 CAPS_NAMES=1 CAPS_QUAL=1 DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 \
        ARCHIVE="$2" "$W/enc" "$1" 3 16 16 22 16 16 1 24 64 1 >"$2.log" 2>&1; }

echo "[roundtrip] lossless on synthetic inputs"
for spec in "200 100" "1000 150" "5000 75"; do
  set -- $spec; N=$1; L=$2
  gen "$N" "$L" "$W/r.fq" "$N"
  if enc "$W/r.fq" "$W/r.arc"; then
    rm -rf "$W/o"; mkdir -p "$W/o"
    "$W/dec" "$W/r.arc" "$W/o" "$W/o/reads.seq" >/dev/null 2>&1
    awk 'NR%4==2' "$W/r.fq" > "$W/orig.seq"
    cmp -s "$W/orig.seq" "$W/o/reads.seq" && ok "roundtrip ${N}x${L}bp lossless" \
      || no "roundtrip ${N}x${L}bp lossless" "sequence differs"
  else no "roundtrip ${N}x${L}bp" "encode failed"; fi
done

echo "[edge] degenerate and out-of-scope inputs must be REFUSED, never crash or silently corrupt"
# Empty input used to succeed silently, producing a structurally valid archive
# containing nothing. That is now a defect, not a feature: the encoder must
# REFUSE (nonzero exit, no archive written, no crash) rather than report
# success on data it did nothing with.
: > "$W/empty.fq"
if enc "$W/empty.fq" "$W/e0.arc"; then
  no "empty input refused" "encoder reported success on an empty file"
elif [ -s "$W/e0.arc" ]; then
  no "empty input refused" "exited nonzero but still wrote an archive"
else
  ok "empty input refused (no archive, no crash)"
fi
# Any read past the fixed per-read stack buffer size (MAX_READ_LEN in
# stages/106_inprocess.cpp) must be refused loudly and specifically -- this is
# the boundary between short-read (supported) and long-read (Nanopore/PacBio,
# unsupported) input. Repeated several times: the failure mode this guards
# against was a std::thread destructor race that only reproduced outside a
# debugger, so a single pass proves nothing.
python3 -c "
import random; random.seed(11)
g=''.join(random.choice('ACGT') for _ in range(5000))
with open('$W/long.fq','w') as f:
    for i in range(20):
        f.write(f'@long{i}\n{g}\n+\n{\"I\"*len(g)}\n')"
LONG_OK=1
for rep in 1 2 3 4 5; do
  rm -f "$W/elong.arc"
  env CAPS_NAMES=1 CAPS_QUAL=1 CAPS_CALL=1 DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 \
      ARCHIVE="$W/elong.arc" "$W/enc" "$W/long.fq" 3 16 16 22 16 16 1 24 64 1 \
      >"$W/elong.log" 2>&1
  RC=$?
  if [ "$RC" -eq 0 ] || [ -s "$W/elong.arc" ]; then LONG_OK=0; break; fi
  if [ "$RC" -gt 128 ]; then LONG_OK=0; break; fi        # crashed (signal), not refused
  grep -q "structural limit" "$W/elong.log" || { LONG_OK=0; break; }
done
[ "$LONG_OK" -eq 1 ] && ok "oversize (long-read-shaped) input refused, 5/5 reps, no crash" \
                     || no "oversize input refused" "see $W/elong.log"
printf "@r\nACGT\n+\nIIII\n" > "$W/t.fq";              enc "$W/t.fq" "$W/e1.arc"    && ok "4bp single read"     || no "4bp single read"
# A header long enough to overrun nmc's per-token model arrays (MAXTOK=1024)
# used to silently truncate its token stream and decode to something other
# than what was encoded -- names, not sequence, and only under CAPS_NAMES=1.
# Same contract as the oversize-read check: refuse loudly, never crash, never
# silently corrupt.
python3 -c "
hdr='@'+'.'.join(str(i) for i in range(2000))
with open('$W/bighdr.fq','w') as f:
    for i in range(5): f.write(f'{hdr}_{i}\nACGTACGTAC\n+\nIIIIIIIIII\n')"
rm -f "$W/bh.arc"
env CAPS_NAMES=1 DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 ARCHIVE="$W/bh.arc" \
    "$W/enc" "$W/bighdr.fq" 3 16 16 22 16 16 1 24 64 1 >"$W/bh.log" 2>&1
BHRC=$?
if [ "$BHRC" -eq 0 ] || [ -s "$W/bh.arc" ]; then
  no "oversize header (CAPS_NAMES) refused" "encoder reported success"
elif [ "$BHRC" -gt 128 ]; then
  no "oversize header (CAPS_NAMES) refused" "crashed (signal), not refused"
elif grep -q "tokenizer bound" "$W/bh.log"; then
  ok "oversize header (CAPS_NAMES) refused, no crash"
else
  no "oversize header (CAPS_NAMES) refused" "wrong reason: see $W/bh.log"
fi
python3 -c "
open('$W/n.fq','w').write(''.join(f'@r{i}\n{\"N\"*100}\n+\n{\"I\"*100}\n' for i in range(50)))"
if enc "$W/n.fq" "$W/e2.arc"; then
  rm -rf "$W/on"; mkdir -p "$W/on"
  "$W/dec" "$W/e2.arc" "$W/on" "$W/on/reads.seq" >/dev/null 2>&1
  awk 'NR%4==2' "$W/n.fq" > "$W/on.orig"
  cmp -s "$W/on.orig" "$W/on/reads.seq" && ok "all-N reads lossless" || no "all-N reads lossless" "N handling regressed"
else no "all-N reads lossless" "encode failed"; fi

echo "[robustness] malformed archives must be REJECTED, not misread"
head -c 4096 /dev/urandom > "$W/junk"
"$W/dec" "$W/junk" "$W/oj" >/dev/null 2>&1; [ $? -ne 0 ] && ok "random data rejected" || no "random data rejected" "exit 0 on garbage"
printf 'CAPSULE\0\x63\x00\x00\x00xx' > "$W/badver"
"$W/dec" "$W/badver" "$W/ov" >/dev/null 2>&1; [ $? -ne 0 ] && ok "bad version rejected" || no "bad version rejected"
gen 500 100 "$W/tr.fq" 7; enc "$W/tr.fq" "$W/tr.arc"
head -c $(( $(stat -c%s "$W/tr.arc") / 2 )) "$W/tr.arc" > "$W/trunc.arc"
"$W/dec" "$W/trunc.arc" "$W/ot" >/dev/null 2>&1; [ $? -ne 0 ] && ok "truncated archive rejected" || no "truncated archive rejected"

echo "[claim3] addressable operations, exact identities"
gen 3000 120 "$W/c3.fq" 11
if enc "$W/c3.fq" "$W/c3.arc"; then
  "$W/dec" export   "$W/c3.arc" "$W/x.fa"  >"$W/x.log" 2>&1
  PG=$(grep -oE 'pg rebuilt: [0-9]+' "$W/x.log" | grep -oE '[0-9]+$' | head -1)
  B=$(grep -v '^>' "$W/x.fa" | tr -d '\n' | wc -c)
  [ -n "$PG" ] && [ "$B" = "$PG" ] && ok "export emits exactly PG_LEN bases" || no "export emits exactly PG_LEN bases" "got $B want ${PG:-?}"
  "$W/dec" coverage "$W/c3.arc" "$W/c.tsv" >/dev/null 2>&1
  COV=$(awk 'NR>1{t+=($3-$2)*$4} END{print t+0}' "$W/c.tsv")
  LEN=$(awk 'NR%4==2{t+=length($0)} END{print t+0}' "$W/c3.fq")
  [ "$COV" = "$LEN" ] && ok "coverage identity exact" || no "coverage identity exact" "$COV vs $LEN"
  "$W/dec" query "$W/c3.arc" "$W/q.fa" 0-5000 >/dev/null 2>&1
  OUT=$(grep -c '^>' "$W/q.fa" 2>/dev/null || echo 0)
  BAD=$(grep '^>' "$W/q.fa" 2>/dev/null | grep -oE 'pos=[0-9]+' | grep -oE '[0-9]+' | awk '$1>=5000' | wc -l)
  [ "$OUT" -gt 0 ] && [ "$BAD" = "0" ] && ok "query returns only overlapping reads ($OUT)" \
    || no "query returns only overlapping reads" "$OUT returned, $BAD out of range"
else no "claim3 operations" "encode failed"; fi

echo "[determinism] same input must give the same archive"
gen 800 110 "$W/d.fq" 3
enc "$W/d.fq" "$W/d1.arc"; enc "$W/d.fq" "$W/d2.arc"
cmp -s "$W/d1.arc" "$W/d2.arc" && ok "encoder is deterministic" || no "encoder is deterministic" "two runs differ"

echo "==================="
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ] || exit 1
