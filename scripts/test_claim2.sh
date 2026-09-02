#!/bin/bash
# Claim 2 synthetic regression test — deterministic, no GIAB download needed.
#
# Before this, include/caps_caller.h (2132 lines) had ZERO automated tests.
# Every existing verification was either a real-data benchmark (needs GIAB
# downloads, hours to set up) or hand-inspection of one window's VCF. This
# gives a fast (<1 min), self-contained sanity check that would catch a
# caller regression immediately: a synthetic diploid genome with KNOWN het
# SNVs and clean (non-homopolymer) het indels at KNOWN positions, reads
# simulated from both haplotypes, run through the exact same
# contig-call -> bwa-align -> lift_vcf.py pipeline the real GIAB benchmark
# uses (scripts/run_giab_indel_capsule.sh), then checked for recall against
# the truth this script itself generated -- not hand-computed pg bytes.
#
# This is NOT a replacement for the real GIAB benchmark (real human
# homopolymer/STR indel behavior cannot be synthesized honestly -- see
# docs/HET_INDEL_SOTA.md). It exists to catch "the caller crashed" or "the
# caller stopped finding anything" regressions cheaply, on every change,
# without needing real data on hand.
#
# Usage: bash scripts/test_claim2.sh [workdir]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
W="${1:-$(mktemp -d)}"; mkdir -p "$W"
PASS=0; FAIL=0
check() { if [ "$1" = "1" ]; then echo "  [PASS] $2"; PASS=$((PASS+1)); else echo "  [FAIL] $2"; FAIL=$((FAIL+1)); fi; }

command -v bwa >/dev/null || { echo "SKIP: bwa not installed, cannot lift contig calls to genome coords"; exit 0; }

echo "== Claim 2 synthetic regression test =="
echo "workdir: $W"

# ── Generate a synthetic diploid genome with KNOWN het variants ─────────────
# 5000bp haploid reference; hap1 = reference; hap2 = reference + 15 het SNVs
# + 5 clean (non-homopolymer) het indels, all spaced >=150bp apart so no two
# variants share a read and no indel sits in a repeat run. 30x coverage,
# 150bp reads, sampled from BOTH haplotypes into one FASTQ -- exactly the
# "one diploid sample, no reference" scenario this caller targets.
python3 - "$W" <<'PYEOF'
import random, sys, os
random.seed(4242)
W = sys.argv[1]
GLEN = 5000
ref = list(''.join(random.choice('ACGT') for _ in range(GLEN)))

# Decide variant positions/kinds in REFERENCE coordinates first (all reads
# against the immutable `ref` list), then apply them to a hap2 COPY in
# HIGH-TO-LOW position order so an earlier (lower-position) edit is never
# shifted by a later (higher-position) one that hasn't been applied yet --
# indexing a shrinking list left-to-right silently drifts every downstream
# truth position by the cumulative prior deletion size (caught by this
# test's own recall check failing at 1/10 on the first run: the caller was
# actually finding the variants, just not where truth.txt said they were).
plan = []  # list of (pos, kind) in reference coordinates
pos = 150
while pos < GLEN - 150:
    if len(plan) % 3 != 0 and sum(1 for _,k in plan if k=='SNV') < 15:
        plan.append((pos, 'SNV'))
    elif sum(1 for _,k in plan if k=='INDEL') < 5:
        if ref[pos] != ref[pos+1] and ref[pos+1] != ref[pos+2]:
            plan.append((pos, 'INDEL'))
    pos += 170

hap2 = ref[:]
truth_snv, truth_indel = [], []
for p, kind in sorted(plan, reverse=True):  # high-to-low: stable indices
    if kind == 'SNV':
        alt = random.choice([b for b in 'ACGT' if b != ref[p]])
        hap2[p] = alt
        truth_snv.append(p)
    else:
        del hap2[p:p+2]
        truth_indel.append(p)

refseq = ''.join(ref)
hap2seq = ''.join(hap2)
with open(f"{W}/ref.fa","w") as f:
    f.write(">synth\n")
    for i in range(0,len(refseq),60): f.write(refseq[i:i+60]+"\n")

RLEN=150; DEPTH=30
def emit_reads(seq, tag, fh):
    n_reads = int(len(seq)*DEPTH/RLEN)
    for i in range(n_reads):
        p = random.randint(0, max(0,len(seq)-RLEN))
        r = seq[p:p+RLEN]
        if len(r) < RLEN: continue
        fh.write(f"@{tag}_{i}\n{r}\n+\n{'I'*len(r)}\n")

with open(f"{W}/synth.fq","w") as fh:
    emit_reads(refseq, "h1", fh)
    emit_reads(hap2seq, "h2", fh)

with open(f"{W}/truth.txt","w") as f:
    for p in truth_snv: f.write(f"SNV\t{p+1}\n")     # 1-based, matches VCF convention
    for p in truth_indel: f.write(f"INDEL\t{p+1}\n")

print(f"truth: {len(truth_snv)} het SNVs, {len(truth_indel)} het indels", file=sys.stderr)
PYEOF

N_READS=$(( $(wc -l < "$W/synth.fq") / 4 ))
echo "-> $N_READS reads simulated from both haplotypes"

# ── Build + run the caller (same invocation as run_giab_indel_capsule.sh) ──
BEST="$W/best106"; [ -x "$BEST" ] || bash "$HERE/scripts/build106.sh" "$BEST" >/dev/null
cd "$W"
export CAPS_CALL=1 CALL_VCF="$W/calls.vcf" CAPS_DUMP_CONTIGS="$W/contigs.tsv"
"$BEST" "$W/synth.fq" 3 16 16 22 16 16 1 24 64 1 > /dev/null 2> "$W/capsule_call.log"
cp "$W/contigs.tsv" "$W/contigs.fa"

check "$([ -s "$W/calls.vcf" ] && echo 1 || echo 0)" "caller produced a non-empty calls.vcf"

# ── Lift contig-coordinate calls to genome coordinates ──────────────────────
bwa index "$W/ref.fa" >/dev/null 2>&1
bwa mem -t2 "$W/ref.fa" "$W/contigs.fa" 2>/dev/null > "$W/c2r.sam"
python3 "$HERE/scripts/lift_vcf.py" "$W/calls.vcf" "$W/c2r.sam" "$W/ref.fa" synth "$W/lifted.vcf" "$W/contigs.fa" 2>"$W/lift.log" || true

check "$([ -s "$W/lifted.vcf" ] && echo 1 || echo 0)" "lift_vcf.py produced a non-empty genome-coordinate VCF"

# ── Compare lifted calls against the truth this script generated ────────────
( python3 - "$W" <<'PYEOF'
import sys
W = sys.argv[1]
truth_snv, truth_indel = set(), set()
for line in open(f"{W}/truth.txt"):
    kind, pos = line.split()
    (truth_snv if kind=="SNV" else truth_indel).add(int(pos))

called = set()
for line in open(f"{W}/lifted.vcf"):
    if line.startswith('#'): continue
    f = line.split('\t')
    called.add(int(f[1]))

TOL = 2  # +/-2bp tolerance for indel coordinate conventions
def recall(truth):
    found = sum(1 for t in truth if any(abs(t-c) <= TOL for c in called))
    return found, len(truth)

snv_found, snv_total = recall(truth_snv)
indel_found, indel_total = recall(truth_indel)
print(f"RESULT_SNV {snv_found}/{snv_total}")
print(f"RESULT_INDEL {indel_found}/{indel_total}")
PYEOF
) > "$W/result.txt"
cat "$W/result.txt"

SNV_FOUND=$(grep RESULT_SNV "$W/result.txt" | awk '{split($2,a,"/"); print a[1]}')
SNV_TOTAL=$(grep RESULT_SNV "$W/result.txt" | awk '{split($2,a,"/"); print a[2]}')
INDEL_FOUND=$(grep RESULT_INDEL "$W/result.txt" | awk '{split($2,a,"/"); print a[1]}')
INDEL_TOTAL=$(grep RESULT_INDEL "$W/result.txt" | awk '{split($2,a,"/"); print a[2]}')

# Thresholds are deliberately loose (this is a crash/regression detector, not
# a precision benchmark -- real accuracy numbers come from the GIAB
# benchmark only, see docs/HET_INDEL_SOTA.md). A healthy caller should find
# nearly all of these clean, well-separated, non-homopolymer variants; the
# floor exists to catch "the caller stopped finding anything" outright.
check "$([ "$SNV_TOTAL" -gt 0 ] && [ "$SNV_FOUND" -ge $((SNV_TOTAL*7/10)) ] && echo 1 || echo 0)" \
    "recovered >=70% of clean synthetic het SNVs ($SNV_FOUND/$SNV_TOTAL)"
check "$([ "$INDEL_TOTAL" -gt 0 ] && [ "$INDEL_FOUND" -ge $((INDEL_TOTAL*4/10)) ] && echo 1 || echo 0)" \
    "recovered >=40% of clean synthetic het indels ($INDEL_FOUND/$INDEL_TOTAL) -- low floor is deliberate, see docs/HET_INDEL_SOTA.md for why indel recall is the harder side of this caller"

echo ""
echo "== $PASS passed, $FAIL failed =="
[ "$FAIL" -eq 0 ] || exit 1
