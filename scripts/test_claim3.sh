#!/bin/bash
# Claim 3 synthetic regression test — deterministic, no external data needed.
#
# This does NOT hand-predict the pseudogenome's exact bytes (the overlap
# chaining + MEM matching that builds it is too complex to hand-verify for
# arbitrary input). Instead it checks INVARIANTS that must hold regardless of
# assembly details, chosen specifically because each one is exactly what a
# past real bug violated (see docs/CLAIM3_LOCKED.md sec 3):
#
#   export   : output is non-empty and pure ACGT (no assembly artifact bytes)
#   coverage : total covered-base-units (sum of depth*width) equals the sum
#              of every ORIGINAL read's length -- this is precisely the
#              invariant the 20% duplicate-read undercount (bug 1) violated.
#              A future regression of that bug fails this test immediately.
#   query    : querying the full pseudogenome range returns exactly one
#              record per UNIQUE read (not per original read) -- verifies the
#              dedup semantics documented in capsule_decode.cpp are what
#              actually happens, and that every read's placement overlaps
#              the range it was returned for (no off-by-one leakage).
#   roundtrip: the underlying archive is still lossless (delegates to the
#              existing scripts/verify_lossless.sh so this test cannot regress
#              silently if compression correctness breaks first).
#
# Usage: bash scripts/test_claim3.sh [workdir]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
W="${1:-$(mktemp -d)}"; mkdir -p "$W"
PASS=0; FAIL=0
check() { if [ "$1" = "1" ]; then echo "  [PASS] $2"; PASS=$((PASS+1)); else echo "  [FAIL] $2"; FAIL=$((FAIL+1)); fi; }

echo "== Claim 3 synthetic regression test =="
echo "workdir: $W"

# ── Generate a deterministic synthetic FASTQ ─────────────────────────────────
# 300 reads of length 60, tiled with overlap across a fixed 3000bp pseudo-
# genome (seeded PRNG, same every run), plus 40 EXACT duplicates injected --
# duplicates are the specific case bug 1 silently dropped.
python3 - "$W/synth.fq" <<'PYEOF'
import random, sys
random.seed(1234)
GLEN=3000
genome=''.join(random.choice('ACGT') for _ in range(GLEN))
RLEN=60
reads=[]
pos=0
while pos+RLEN<=GLEN:
    reads.append(genome[pos:pos+RLEN])
    pos+=17  # overlapping tiling, non-multiple-of-RLEN stride
# inject 40 exact duplicates of already-emitted reads
for i in range(40):
    reads.append(reads[i*3 % len(reads)])
with open(sys.argv[1],'w') as f:
    for i,r in enumerate(reads):
        f.write(f"@r{i}\n{r}\n+\n{'I'*RLEN}\n")
print(f"synthetic reads: {len(reads)}, read length: {RLEN}", file=sys.stderr)
PYEOF

N_READS=$(( $(wc -l < "$W/synth.fq") / 4 ))
SUM_LEN=$(awk 'NR%4==2{s+=length($0)} END{print s}' "$W/synth.fq")
echo "-> $N_READS reads, total bases $SUM_LEN"

# ── Build + encode + decode ──────────────────────────────────────────────────
BEST="$W/best106"; [ -x "$BEST" ] || bash "$HERE/scripts/build106.sh" "$BEST" >/dev/null
DEC="$W/capsule_decode"; [ -x "$DEC" ] || bash "$HERE/scripts/build_decode.sh" "$DEC" >/dev/null
ARC="$W/synth.capsule"
INPUT="$W/synth.fq" ARCHIVE="$ARC" BEST="$BEST" bash "$HERE/scripts/encode_adaptive.sh" >"$W/encode.log" 2>&1
# encode_adaptive.sh's own stdout only has ARCHIVE_TOTAL; the per-candidate
# "unique=" line is in the archive's own .log sibling file (see the
# DUMP_LIT/DUMP_PERM/DUMP_MM invocation inside encode_adaptive.sh).
UNIQUE=$(grep -oP 'unique=\K[0-9]+' "${ARC}.log" | tail -1)
echo "-> encoder reports unique=$UNIQUE"

# ── export: non-empty, pure ACGT ─────────────────────────────────────────────
"$DEC" export "$ARC" "$W/export.fa" 2>"$W/export.log"
BAD_CHARS=$(grep -v '^>' "$W/export.fa" | tr -d 'ACGTN\n' | wc -c)
EXPORT_NONEMPTY=$([ -s "$W/export.fa" ] && echo 1 || echo 0)
check "$EXPORT_NONEMPTY" "export produced non-empty output"
check "$([ "$BAD_CHARS" -eq 0 ] && echo 1 || echo 0)" "export output is pure ACGT/N (0 stray bytes, found $BAD_CHARS)"

# ── coverage: total covered-base-units == sum of ALL original read lengths ──
"$DEC" coverage "$ARC" "$W/coverage.tsv" 2>"$W/coverage.log"
COV_UNITS=$(awk -F'\t' 'NR>1{u+=($3-$2)*$4} END{print u}' "$W/coverage.tsv")
check "$([ "$COV_UNITS" = "$SUM_LEN" ] && echo 1 || echo 0)" \
    "coverage total ($COV_UNITS) == sum of all $N_READS original read lengths ($SUM_LEN) -- this is exactly the invariant bug 1 (20% undercount) violated"

# ── query: full-range query returns exactly one record per UNIQUE read ──────
PGLEN=$(awk -F'\t' 'NR>1{if($3>m)m=$3} END{print m}' "$W/coverage.tsv")
"$DEC" query "$ARC" "$W/query_full.fa" "0-$PGLEN" 2>"$W/query.log"
N_QUERY=$(grep -c '^>' "$W/query_full.fa")
check "$([ "$N_QUERY" = "$UNIQUE" ] && echo 1 || echo 0)" \
    "full-range query returned $N_QUERY records == encoder's unique count $UNIQUE (dedup semantics correct)"

# every returned read's stated position must fall inside [0,PGLEN)
BAD_POS=$(grep '^>' "$W/query_full.fa" | grep -oP 'pos=\K[0-9]+' | awk -v p="$PGLEN" '$1>=p{c++} END{print c+0}')
check "$([ "$BAD_POS" -eq 0 ] && echo 1 || echo 0)" "no query record has a position >= PGLEN ($PGLEN)"

# ── roundtrip: the archive is still lossless ─────────────────────────────────
if bash "$HERE/scripts/verify_lossless.sh" "$W/synth.fq" "$W/vl" > "$W/vl.log" 2>&1; then
    check 1 "underlying archive round-trips LOSSLESS (scripts/verify_lossless.sh)"
else
    check 0 "underlying archive round-trips LOSSLESS (scripts/verify_lossless.sh) -- see $W/vl.log"
fi

echo ""
echo "== $PASS passed, $FAIL failed =="
[ "$FAIL" -eq 0 ] || exit 1
