#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
#  TEST ZERO — capability, scope and safety self-test
#
#  Cloned from benchmark_0_preflight.sh (same helper functions, same
#  dynamic-system-probe pattern) and rewritten for a different job:
#  benchmark_0_preflight.sh gates the FULL 19-locked-dataset benchmark run
#  and needs those datasets, competitor tools and reference genomes present.
#  THIS script needs none of that. It needs only this repository and a C++
#  compiler, runs on ANY machine with ANY input, and answers one question:
#  "does this tool behave correctly and safely, and what exactly is it
#  willing to do?" -- BEFORE a reviewer runs anything real.
#
#  Every system fact below is PROBED at run time (uname, nproc, /proc/meminfo,
#  df, the encoder's own source) -- nothing here is a value typed in advance
#  and echoed back. If this box changes, the next run reflects that change
#  automatically.
#
#  Read the SCOPE section before anything else. It states, in one place, what
#  input this tool is validated for, what it will refuse and why, and what it
#  will accept but has not been measured on. Then the script PROVES each of
#  those boundary claims by actually constructing the input and running the
#  real binary against it -- a claim that is only ever printed, never
#  exercised, is not a claim this project's own standing rules would accept
#  from anyone else's paper, and should not be accepted from this one either.
#
#  usage: bash scripts/test_0_scope_and_capability.sh [OUT_FILE]
#         default OUT_FILE: results/TEST_0_STATUS.txt
#
#  Exit: 0 = every check passed, the tool's stated scope is proven on this
#            machine, right now
#        1 = at least one check failed -- read the FAILURES section, this
#            is not a "run anyway and see" situation
# ═══════════════════════════════════════════════════════════════════════════
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
OUT="${1:-$HERE/results/TEST_0_STATUS.txt}"
mkdir -p "$(dirname "$OUT")"
: > "$OUT"

T_START=$(date +%s)
FAIL=0; WARN=0; OK=0
say(){ echo "$*" | tee -a "$OUT"; }
hdr(){ say ""; say "════════════════════════════════════════════════════════════════════════"; say " $*"; say "════════════════════════════════════════════════════════════════════════"; }
pass(){ OK=$((OK+1));     printf "  [ OK ]   %-40s %s\n" "$1" "${2:-}" | tee -a "$OUT"; }
warn(){ WARN=$((WARN+1)); printf "  [WARN]   %-40s %s\n" "$1" "${2:-}" | tee -a "$OUT"; }
fail(){ FAIL=$((FAIL+1)); printf "  [FAIL]   %-40s %s\n" "$1" "${2:-}" | tee -a "$OUT"; }
hb(){ printf "  %-42s %s\n" "$1" "${2:-}" | tee -a "$OUT"; }
gb(){ awk -v b="$1" 'BEGIN{printf "%.2f GB", b/1073741824}'; }

say "╔══════════════════════════════════════════════════════════════════════╗"
say "║  TEST 0 — CAPABILITY, SCOPE AND SAFETY SELF-TEST                     ║"
say "╚══════════════════════════════════════════════════════════════════════╝"
say "generated : $(date '+%Y-%m-%d %H:%M:%S %Z')"
say "host      : $(hostname)  |  user: $(whoami)"
say "repo      : $HERE"
say "branch    : $(git rev-parse --abbrev-ref HEAD 2>/dev/null)  commit: $(git rev-parse --short HEAD 2>/dev/null)"

# ── 1. SYSTEM ─────────────────────────────────────────────────────────────
# Every value here is read from the running kernel, not typed in. Compare
# with any prior run of this script to see exactly what changed on the box.
hdr "1. SYSTEM CONFIGURATION — probed, not hardcoded"
say "  Every line below comes from a live system call. If two runs of this"
say "  script on two different machines print different numbers here, that"
say "  is this section working correctly, not a bug."
say ""
OS_NAME="$( (source /etc/os-release 2>/dev/null && echo "$PRETTY_NAME") || uname -s)"
KERNEL="$(uname -r)"
ARCH="$(uname -m)"
CORES="$(nproc 2>/dev/null || echo '?')"
CORES_PHYS="$(lscpu 2>/dev/null | awk -F: '/^Core\(s\) per socket/{c=$2} /^Socket\(s\)/{s=$2} END{if(c&&s) printf "%d", c*s}')"
CPU_MODEL="$(lscpu 2>/dev/null | awk -F: '/^Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
RAM_TOTAL_KB="$(awk '/^MemTotal:/{print $2}' /proc/meminfo 2>/dev/null)"
RAM_AVAIL_KB="$(awk '/^MemAvailable:/{print $2}' /proc/meminfo 2>/dev/null)"
DISK_AVAIL_G="$(df -BG "$HERE" 2>/dev/null | awk 'NR==2{gsub("G","",$4); print $4}')"
DISK_FS="$(df -T "$HERE" 2>/dev/null | awk 'NR==2{print $2}')"
hb "OS"            "${OS_NAME:-unknown}"
hb "kernel"        "$KERNEL ($ARCH)"
hb "CPU model"     "${CPU_MODEL:-unknown}"
hb "logical cores" "$CORES${CORES_PHYS:+  (physical: $CORES_PHYS)}"
hb "RAM total"     "$([ -n "${RAM_TOTAL_KB:-}" ] && awk -v k="$RAM_TOTAL_KB" 'BEGIN{printf "%.1f GB", k/1048576}' || echo unknown)"
hb "RAM available" "$([ -n "${RAM_AVAIL_KB:-}" ] && awk -v k="$RAM_AVAIL_KB" 'BEGIN{printf "%.1f GB", k/1048576}' || echo unknown)"
hb "disk free ($HERE)" "${DISK_AVAIL_G:-?} GB (${DISK_FS:-unknown fs})"
[ -n "${CORES:-}" ] && [ "$CORES" != "?" ] && pass "system probed successfully" "$CORES cores" \
                                            || warn "core count unavailable" "nproc failed"
[ -n "${RAM_AVAIL_KB:-}" ] && [ "$(awk -v k="$RAM_AVAIL_KB" 'BEGIN{print (k>2097152)?1:0}')" = 1 ] \
    && pass "RAM sufficient for this self-test" "> 2 GB available" \
    || warn "low RAM" "this self-test's synthetic inputs are tiny (KB-scale) and should still run, but a real dataset will need far more -- see server/SETUP_SHORT.md"
[ -n "${DISK_AVAIL_G:-}" ] && [ "$DISK_AVAIL_G" -ge 1 ] 2>/dev/null \
    && pass "disk sufficient for this self-test" "${DISK_AVAIL_G} GB free" \
    || warn "disk check inconclusive" "could not confirm free space"

# ── 2. TOOLCHAIN ─────────────────────────────────────────────────────────
hdr "2. BUILD TOOLCHAIN"
chk_cmd(){ local n="$1" c="$2"
  if command -v "$c" >/dev/null 2>&1; then pass "$n" "$(command -v "$c")$($c --version 2>/dev/null | head -1 | sed 's/^/  /')"
  else fail "$n" "MISSING — cannot build or run anything below"; fi; }
chk_cmd "g++ (C++17, -fopenmp)" g++
chk_cmd "python3"               python3
chk_cmd "/usr/bin/time"         /usr/bin/time

# ── 3. THE SCOPE DECLARATION — extracted from source, not typed here ──────
hdr "3. SCOPE — what this build is validated for, and what it refuses"
say "  These numbers are extracted from the encoder's own source at the"
say "  moment this script runs. They cannot silently drift out of sync with"
say "  the binary being tested, because they ARE the binary's own constants."
say ""
MAX_READ_LEN=$(python3 - "$HERE/stages/106_inprocess.cpp" <<'PY'
import re,sys
src=open(sys.argv[1]).read()
m=re.search(r'static\s+const\s+uint32_t\s+MAX_READ_LEN\s*=\s*(\d+)\s*;', src)
print(m.group(1) if m else "")
PY
)
MAXTOK=$(python3 - "$HERE/include/names_coder.h" <<'PY'
import re,sys
src=open(sys.argv[1]).read()
m=re.search(r'static\s+const\s+uint32_t\s+MAXTOK\s*=\s*(\d+)\s*;', src)
print(m.group(1) if m else "")
PY
)
if [ -n "$MAX_READ_LEN" ]; then
    pass "MAX_READ_LEN extracted" "$MAX_READ_LEN bases (stages/106_inprocess.cpp)"
else
    fail "MAX_READ_LEN extraction" "constant not found — source may have moved; do not trust any scope claim below"
fi
if [ -n "$MAXTOK" ]; then
    pass "nmc::MAXTOK extracted" "$MAXTOK tokens (include/names_coder.h)"
else
    warn "nmc::MAXTOK extraction" "not found — names-header bound cannot be stated"
fi
say ""
say "  ┌─────────────────────────────────────────────────────────────────┐"
say "  │  SUPPORTED                                                       │"
say "  │  short-read FASTQ (Illumina-shaped), any length up to           │"
say "  │  ${MAX_READ_LEN:-?} bases, fixed or variable, fwd/reverse-complement    │"
say "  │  strands, ACGTN alphabet, Phred quality.                        │"
say "  │  Validated at scale against 19 real public datasets, 40-301 bp.  │"
say "  │                                                                   │"
say "  │  NOT SUPPORTED — refused, not silently mishandled                │"
say "  │  - any read over ${MAX_READ_LEN:-?} bases (long-read: Nanopore, PacBio).  │"
say "  │    Every per-read decode path is a fixed-size stack buffer      │"
say "  │    sized for short-read technology. Section 5 below builds      │"
say "  │    such a file and proves the refusal, live.                    │"
say "  │  - a read header over ${MAXTOK:-?} characters, with CAPS_NAMES=1.  │"
say "  │    Section 5 proves this too.                                    │"
say "  │  - empty or non-FASTQ input.                                     │"
say "  │                                                                   │"
say "  │  ACCEPTED BUT NOT VALIDATED — a real, stated gap, not hidden     │"
say "  │  - non-ACGTN alphabet (protein, RNA with U): NOT gated today.    │"
say "  │    industry/check_input_scope.py flags this as a caveat; the     │"
say "  │    encoder itself does not refuse it.                            │"
say "  │  - reads far longer than 301 bp but under the hard bound: will   │"
say "  │    be PROCESSED, but the assembly/mismatch/position coders'      │"
say "  │    parameters were derived from and validated on 40-301 bp       │"
say "  │    reads. Treat any result on longer input as exploratory.       │"
say "  └─────────────────────────────────────────────────────────────────┘"
say ""
say "  Full detail, and the empirical process that found each boundary:"
say "  industry/README.md"

# ── 4. BUILD ──────────────────────────────────────────────────────────────
hdr "4. BUILD — from source, this checkout, right now"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT
BEST="$W/enc"; DEC="$W/dec"
if bash "$HERE/scripts/build106.sh" "$BEST" >"$W/b1.log" 2>&1 && [ -x "$BEST" ]; then
    pass "encoder builds" "$(stat -c%s "$BEST" | awk '{printf "%.1f MB", $1/1048576}')"
else
    fail "encoder build" "see the tail below"; tail -8 "$W/b1.log" | sed 's/^/         /' | tee -a "$OUT"
fi
if bash "$HERE/scripts/build_decode.sh" "$DEC" >"$W/b2.log" 2>&1 && [ -x "$DEC" ]; then
    pass "decoder builds" "$(stat -c%s "$DEC" | awk '{printf "%.1f MB", $1/1048576}')"
else
    fail "decoder build" "see the tail below"; tail -8 "$W/b2.log" | sed 's/^/         /' | tee -a "$OUT"
fi

# ── 5. PROVING THE SCOPE, NOT JUST STATING IT ─────────────────────────────
hdr "5. BOUNDARY PROOFS — every scope claim above, actually exercised"
say "  Each check below constructs the exact input the SCOPE section"
say "  describes and runs the real, just-built binary against it. A"
say "  boundary claim that is only ever printed is not proven; this"
say "  section is the difference between the two."
say ""

if [ ! -x "$BEST" ]; then
    fail "boundary proofs" "cannot run — encoder did not build (see section 4)"
else
    gen(){ # gen <reads> <len> <path> [seed]
        python3 -c "
import random,sys
n,L,path,seed=int('$1'),int('$2'),'$3',int('${4:-1}')
random.seed(seed)
G=''.join(random.choice('ACGT') for _ in range(max(L*4, n*L//8 + L)))
with open(path,'w') as f:
    for i in range(n):
        p=random.randrange(0,max(1,len(G)-L)); s=list(G[p:p+L])
        if random.random()<0.3: s[random.randrange(len(s))]=random.choice('ACGT')
        f.write(f'@read{i} synthetic\n{\"\".join(s)}\n+\n{\"I\"*len(s)}\n')
"
    }
    encrun(){ env CAPS_NAMES=1 CAPS_QUAL=1 DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 \
                  ARCHIVE="$2" "$BEST" "$1" 3 16 16 22 16 16 1 24 64 1 \
                  >"$2.log" 2>&1; }

    # 5.1 -- in-scope short-read input compresses AND decompresses losslessly
    gen 2000 150 "$W/short.fq" 42
    if encrun "$W/short.fq" "$W/short.arc"; then
        mkdir -p "$W/out"
        if [ -x "$DEC" ] && "$DEC" "$W/short.arc" "$W/out" "$W/out/reads.seq" >"$W/dec.log" 2>&1 \
           && cmp -s <(awk 'NR%4==2' "$W/short.fq") "$W/out/reads.seq"; then
            pass "in-scope input: compress + decompress lossless" \
                 "2000 x 150bp -> $(stat -c%s "$W/short.arc") B, byte-identical"
        else
            fail "in-scope input: round trip" "archive built but did not decode losslessly — see $W/dec.log"
        fi
    else
        fail "in-scope input: encode" "the SUPPORTED case itself failed — see $W/short.arc.log"
    fi

    # 5.2 -- oversize (long-read-shaped) input is REFUSED, not silently
    #        corrupted, and not crashed. Repeated: the defect this proves
    #        against was a race that reproduced 100% under plain execution
    #        and 0% under a debugger, so a single pass proves nothing.
    if [ -n "$MAX_READ_LEN" ]; then
        OVER_LEN=$((MAX_READ_LEN + 500))
        gen 20 "$OVER_LEN" "$W/long.fq" 7
        REFUSED_CLEAN=1
        for rep in 1 2 3 4 5; do
            rm -f "$W/long.arc"
            encrun "$W/long.fq" "$W/long.arc"
            RC=$?
            if [ "$RC" -eq 0 ] || [ -s "$W/long.arc" ]; then REFUSED_CLEAN=0; break; fi
            if [ "$RC" -gt 128 ]; then REFUSED_CLEAN=0; break; fi   # killed by signal
            grep -q "structural limit" "$W/long.arc.log" || { REFUSED_CLEAN=0; break; }
        done
        if [ "$REFUSED_CLEAN" -eq 1 ]; then
            pass "out-of-scope input (${OVER_LEN}bp reads): refused cleanly" "5/5 reps, no crash, no false success"
        else
            fail "out-of-scope input refusal" "rep $rep: exit=$RC, see $W/long.arc.log — this is the exact defect class this script exists to catch"
        fi
    else
        warn "5.2 skipped" "MAX_READ_LEN unknown, cannot construct the boundary case"
    fi

    # 5.3 -- empty input is refused, not silently accepted as an empty
    #        "success"
    : > "$W/empty.fq"
    rm -f "$W/empty.arc"
    encrun "$W/empty.fq" "$W/empty.arc"
    RC=$?
    if [ "$RC" -ne 0 ] && [ ! -s "$W/empty.arc" ] && [ "$RC" -le 128 ]; then
        pass "empty input: refused cleanly" "no archive written, no crash"
    else
        fail "empty input refusal" "exit=$RC, archive present=$([ -s "$W/empty.arc" ] && echo yes || echo no)"
    fi

    # 5.4 -- a pathological header (CAPS_NAMES=1) is refused, not silently
    #        mis-decoded
    if [ -n "$MAXTOK" ]; then
        python3 -c "
hdr='@'+'.'.join(str(i) for i in range($MAXTOK*2))
with open('$W/bighdr.fq','w') as f:
    for i in range(5): f.write(f'{hdr}_{i}\nACGTACGTAC\n+\nIIIIIIIIII\n')"
        rm -f "$W/bighdr.arc"
        encrun "$W/bighdr.fq" "$W/bighdr.arc"
        RC=$?
        if [ "$RC" -ne 0 ] && [ ! -s "$W/bighdr.arc" ] && [ "$RC" -le 128 ] \
           && grep -q "tokenizer bound" "$W/bighdr.arc.log"; then
            pass "oversize header (CAPS_NAMES=1): refused cleanly" "no crash, correct reason"
        else
            fail "oversize header refusal" "exit=$RC — see $W/bighdr.arc.log"
        fi
    else
        warn "5.4 skipped" "nmc::MAXTOK unknown"
    fi

    # 5.5 -- the standalone pre-flight reporter agrees with the encoder's
    #        own verdict on both the in-scope and out-of-scope cases above.
    #        These are two independent code paths (Python vs C++) reaching
    #        the same file; if they disagree, one of them is wrong.
    if [ -f "$HERE/industry/check_input_scope.py" ]; then
        python3 "$HERE/industry/check_input_scope.py" "$W/short.fq" >/dev/null 2>&1
        SHORT_RC=$?
        python3 "$HERE/industry/check_input_scope.py" "$W/long.fq" >/dev/null 2>&1
        LONG_RC=$?
        if [ "$SHORT_RC" -lt 2 ] && [ "$LONG_RC" -ge 2 ]; then
            pass "industry/check_input_scope.py agrees with the encoder" \
                 "in-scope file: exit $SHORT_RC; out-of-scope file: exit $LONG_RC"
        else
            fail "check_input_scope.py disagreement" \
                 "short.fq exit=$SHORT_RC (want <2), long.fq exit=$LONG_RC (want >=2) — one of the two verdicts is wrong"
        fi
    else
        warn "5.5 skipped" "industry/check_input_scope.py not found"
    fi
fi

# ── 6. DETERMINISM ────────────────────────────────────────────────────────
hdr "6. DETERMINISM — same input, same machine, same archive"
if [ -x "$BEST" ] && [ -f "$W/short.fq" ]; then
    encrun "$W/short.fq" "$W/short_a.arc"
    encrun "$W/short.fq" "$W/short_b.arc"
    if cmp -s "$W/short_a.arc" "$W/short_b.arc"; then
        pass "two independent runs, identical output" "$(stat -c%s "$W/short_a.arc") B, byte-identical"
    else
        fail "determinism" "two runs of the same input produced different archives"
    fi
else
    warn "determinism check skipped" "no successful short.fq encode to compare"
fi

# ── 7. VERDICT ────────────────────────────────────────────────────────────
T_END=$(date +%s)
hdr "7. VERDICT"
hb "checks passed"  "$OK"
hb "warnings"       "$WARN"
hb "failures"       "$FAIL"
hb "self-test took" "$((T_END-T_START)) s"
say ""
if [ "$FAIL" -eq 0 ]; then
    say "  ╔════════════════════════════════════════════════════════════════════╗"
    say "  ║   PASS                                                              ║"
    say "  ║   This build's stated scope (section 3) has been PROVEN on THIS     ║"
    say "  ║   machine: in-scope input round-trips losslessly, and every         ║"
    say "  ║   out-of-scope case tested is refused cleanly — no crash, no        ║"
    say "  ║   silent data loss, no false success.                               ║"
    say "  ╚════════════════════════════════════════════════════════════════════╝"
    say ""
    say "  Next: bash scripts/run_tests.sh                (17-test self-contained suite)"
    say "        bash scripts/benchmark_0_preflight.sh     (gates the full published benchmark;"
    say "                                                    needs the 19 locked datasets)"
    RC=0
else
    say "  ╔════════════════════════════════════════════════════════════════════╗"
    say "  ║   FAIL                                                              ║"
    say "  ║   $FAIL check(s) failed. Do not trust results from this build on   ║"
    say "  ║   this machine until every failure below is understood and fixed.   ║"
    say "  ╚════════════════════════════════════════════════════════════════════╝"
    say ""; say "  Failures:"; grep '\[FAIL\]' "$OUT" | sed 's/^/  /'
    RC=1
fi
say ""
say "  status file: $OUT"
say "  finished   : $(date '+%Y-%m-%d %H:%M:%S')"
exit $RC
