#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
#  BENCHMARK ZERO — pre-run dependency and readiness check
#
#  Runs FIRST. Changes NOTHING. Verifies every tool, dataset, reference and
#  truth file the full benchmark needs, checks the code still reproduces its
#  validated numbers, then writes ONE status file with a GO / NO-GO verdict.
#
#  The status file is emptied and rewritten on every run, so what you are
#  reading is always the current state -- never a stale mix of two runs.
#
#  usage: bash scripts/benchmark_0_preflight.sh [OUT_FILE]
#         default OUT_FILE: results/BENCHMARK_0_STATUS.txt
#
#  Exit: 0 = GO, 1 = NO-GO (a required item is missing/broken)
# ═══════════════════════════════════════════════════════════════════════════
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
OUT="${1:-$HERE/results/BENCHMARK_0_STATUS.txt}"
mkdir -p "$(dirname "$OUT")"
: > "$OUT"                      # fresh file, every run

T_START=$(date +%s)
FAIL=0; WARN=0; OK=0
say(){ echo "$*" | tee -a "$OUT"; }
hdr(){ say ""; say "════════════════════════════════════════════════════════════════════════"; say " $*"; say "════════════════════════════════════════════════════════════════════════"; }
pass(){ OK=$((OK+1));     printf "  [ OK ]   %-34s %s\n" "$1" "${2:-}" | tee -a "$OUT"; }
warn(){ WARN=$((WARN+1)); printf "  [WARN]   %-34s %s\n" "$1" "${2:-}" | tee -a "$OUT"; }
fail(){ FAIL=$((FAIL+1)); printf "  [FAIL]   %-34s %s\n" "$1" "${2:-}" | tee -a "$OUT"; }
hb(){ printf "  %-36s %s\n" "$1" "${2:-}" | tee -a "$OUT"; }
gb(){ awk -v b="$1" 'BEGIN{printf "%.2f GB", b/1073741824}'; }

say "╔══════════════════════════════════════════════════════════════════════╗"
say "║  BENCHMARK 0 — PRE-RUN STATUS                                        ║"
say "╚══════════════════════════════════════════════════════════════════════╝"
say "generated : $(date '+%Y-%m-%d %H:%M:%S %Z')"
say "host      : $(hostname)  |  user: $(whoami)"
say "repo      : $HERE"
say "branch    : $(git rev-parse --abbrev-ref HEAD 2>/dev/null)  commit: $(git rev-parse --short HEAD 2>/dev/null)"
say "tree      : $(if [ -z "$(git status --porcelain 2>/dev/null)" ]; then echo CLEAN; else echo 'DIRTY (uncommitted changes present)'; fi)"

# ── 1. SYSTEM ─────────────────────────────────────────────────────────────
hdr "1. SYSTEM RESOURCES"
CORES=$(nproc); RAM_GB=$(free -g | awk '/^Mem:/{print $2}'); DISK_AVAIL_G=$(df -BG / | awk 'NR==2{gsub("G","",$4); print $4}')
hb "cores"        "$CORES"
hb "RAM total"    "${RAM_GB} GB"
hb "RAM available" "$(free -g | awk '/^Mem:/{print $7}') GB"
hb "disk free"    "${DISK_AVAIL_G} GB"
[ "$CORES" -ge 8 ]        && pass "cores >= 8"            "$CORES"        || warn "cores < 8"  "$CORES (slower, not fatal)"
[ "$RAM_GB" -ge 32 ]      && pass "RAM >= 32 GB"          "${RAM_GB} GB"  || fail "RAM < 32 GB" "${RAM_GB} GB"
[ "$DISK_AVAIL_G" -ge 60 ] && pass "disk >= 60 GB free"   "${DISK_AVAIL_G} GB" || fail "disk < 60 GB free" "${DISK_AVAIL_G} GB — peak transient is ~35 GB"

# ── 2. TOOLS ──────────────────────────────────────────────────────────────
hdr "2. TOOLS — state of the art we benchmark against"
chk_cmd(){ local n="$1" c="$2" req="$3" v
  if command -v "$c" >/dev/null 2>&1; then v=$(command -v "$c"); pass "$n" "$v"
  else [ "$req" = req ] && fail "$n" "MISSING — required" || warn "$n" "MISSING — optional"; fi; }
chk_path(){ local n="$1" p="$2" req="$3"
  if [ -e "$p" ]; then pass "$n" "$p"
  else [ "$req" = req ] && fail "$n" "MISSING: $p" || warn "$n" "MISSING: $p"; fi; }

say "  -- Claim 1 (compression) --"
chk_cmd "SPRING"          spring     req
chk_cmd "Genozip"         genozip    req
chk_cmd "Genounzip"       genounzip  req
say "  -- Claim 2 (variant calling) --"
chk_path "DiscoSNP++"     "$HOME/DiscoSnp/run_discoSnp++.sh" req
chk_cmd  "rtg (vcfeval)"  rtg        req
if [ -d "$HOME/miniconda3/envs/kmer2snp_r" ]; then pass "Kmer2SNP (conda env)" "$HOME/miniconda3/envs/kmer2snp_r"
else warn "Kmer2SNP (conda env)" "MISSING — Claim 2 runs with DiscoSNP++ only, T3 loses its 3rd arm"; fi
say "  -- Claim 3 (archive analysis) --"
chk_path "SPAdes"         "$HOME/SPAdes-4.0.0-Linux/bin/spades.py" req
chk_cmd  "bwa"            bwa        req
chk_cmd  "samtools"       samtools   req
chk_cmd  "mosdepth"       mosdepth   req
say "  -- build toolchain --"
chk_cmd "g++"             g++        req
chk_cmd "/usr/bin/time"   /usr/bin/time req

# ── 3. BINARIES ───────────────────────────────────────────────────────────
hdr "3. OUR BINARIES"
source "$HERE/scripts/capsule_config.sh" 2>/dev/null || true
BIN_DIR="${CAPSULE_BIN_DIR:-/tmp/capsule_bin}"; mkdir -p "$BIN_DIR"
BEST="$BIN_DIR/best106"; DEC="$BIN_DIR/capsule_decode"
hb "CAPSULE_BIN_DIR" "$BIN_DIR"
say ""
say "  building encoder (best106) ..."
if bash "$HERE/scripts/build106.sh" "$BEST" >/tmp/_b0_enc.log 2>&1 && [ -x "$BEST" ]; then
    pass "encoder builds"   "$BEST ($(stat -c%s "$BEST" | awk '{printf "%.1f MB", $1/1048576}'))"
else fail "encoder build"   "see /tmp/_b0_enc.log"; tail -5 /tmp/_b0_enc.log | sed 's/^/         /' | tee -a "$OUT"; fi
say "  building decoder (capsule_decode) ..."
if bash "$HERE/scripts/build_decode.sh" "$DEC" >/tmp/_b0_dec.log 2>&1 && [ -x "$DEC" ]; then
    pass "decoder builds"   "$DEC ($(stat -c%s "$DEC" | awk '{printf "%.1f MB", $1/1048576}'))"
else fail "decoder build"   "see /tmp/_b0_dec.log"; tail -5 /tmp/_b0_dec.log | sed 's/^/         /' | tee -a "$OUT"; fi

# ── 4. DATASETS ───────────────────────────────────────────────────────────
hdr "4. DATASETS — 19 total (15 non-human + 4 GIAB human)"
DATA_DIR="${CAPSULE_DATA_DIR:-/data/fastq}"
hb "data dir" "$DATA_DIR"
say ""
printf "  %-4s %-16s %-26s %12s  %s\n" "#" "accession" "organism" "size" "path" | tee -a "$OUT"
say "  ----------------------------------------------------------------------------------------"
TOTAL_BYTES=0; DS_OK=0; DS_MISS=0
check_ds(){ local i="$1" acc="$2" org="$3" f="$DATA_DIR/$4"
    if [ -s "$f" ]; then local s; s=$(stat -c%s "$f"); TOTAL_BYTES=$((TOTAL_BYTES+s)); DS_OK=$((DS_OK+1))
        printf "  %-4s %-16s %-26s %12s  %s\n" "$i" "$acc" "$org" "$(gb "$s")" "$f" | tee -a "$OUT"
    else DS_MISS=$((DS_MISS+1)); FAIL=$((FAIL+1))
        printf "  %-4s %-16s %-26s %12s  %s\n" "$i" "$acc" "$org" "MISSING" "$f" | tee -a "$OUT"; fi; }
check_ds 1  SRR2584863  "E. coli B REL606"       SRR2584863_1.fq
check_ds 2  ERR552797   "M. tuberculosis H37Rv"  ERR552797_1.fq
check_ds 3  SRR554369   "P. aeruginosa PAO1"     SRR554369_1.fq
check_ds 4  ERR5181310  "SARS-CoV-2"             ERR5181310_1.fq
check_ds 5  ERR17740259 "S. aureus"              ERR17740259_1.fq
check_ds 6  DRR976266   "S. cerevisiae"          DRR976266_1.fq
check_ds 7  SRR36741279 "Leishmania major"       SRR36741279_1.fq
check_ds 8  SRR37283774 "P. falciparum"          SRR37283774_1.fq
check_ds 9  SRR32429602 "HCMV"                   SRR32429602_1.fq
check_ds 10 SRR39257532 "Aspergillus fumigatus"  SRR39257532_1.fq
check_ds 11 SRR29296997 "Halobacterium salinarum" SRR29296997_1.fq
check_ds 12 ERR12954017 "Sulfolobus acidocaldarius" ERR12954017_1.fq
check_ds 13 SRR40271341 "Helicobacter pylori"    SRR40271341_1.fq
check_ds 14 SRR065390   "C. elegans N2"          SRR065390_1.fq
check_ds 15 SRR10676752 "Utricularia gibba"      SRR10676752_1.fq
check_ds 16 HG002       "GIAB Ashkenazi son"     HG002_pooled.fq
check_ds 17 HG003       "GIAB Ashkenazi father"  HG003_pooled.fq
check_ds 18 HG004       "GIAB Ashkenazi mother"  HG004_pooled.fq
check_ds 19 HG005       "GIAB Han Chinese son"   HG005_pooled.fq
say "  ----------------------------------------------------------------------------------------"
hb "datasets present" "$DS_OK / 19"
hb "datasets missing" "$DS_MISS"
hb "total input size" "$(gb $TOTAL_BYTES)"
[ "$DS_MISS" -eq 0 ] && pass "all 19 datasets present" || fail "$DS_MISS dataset(s) missing" "cannot run the locked benchmark"

# ── 5. REFERENCES + TRUTH ─────────────────────────────────────────────────
hdr "5. REFERENCES AND TRUTH SETS"
REFS="${CAPSULE_REFS_DIR:-$HOME/refs}"; TRUTH="${CAPSULE_TRUTH_DIR:-$HOME/giab_truth}"
hb "refs dir" "$REFS"; hb "truth dir" "$TRUTH"; say ""
say "  -- Claim 2 --"
for f in chr20.fa chr20.sdf chr1.fa chr1.sdf; do
    if [ -e "$REFS/$f" ]; then pass "$f" "$(du -sh "$REFS/$f" 2>/dev/null | cut -f1)"; else
        case "$f" in chr20.*) fail "$f" "MISSING — required for Claim 2";; *) warn "$f" "MISSING — chr1 axis only";; esac; fi
done
NT=$(ls "$TRUTH"/*.vcf.gz 2>/dev/null | wc -l)
[ "$NT" -ge 4 ] && pass "GIAB truth VCFs" "$NT found" || fail "GIAB truth VCFs" "only $NT found, need 4 (HG002-HG005)"
say "  -- Claim 3 baselines (6 datasets, one per kingdom) --"
C3_REFS="ecoli sarscov2 halobacterium pfalciparum scerevisiae"
for r in $C3_REFS; do
    if [ -s "$REFS/c3_$r.fa" ]; then
        if [ -s "$REFS/c3_$r.fa.bwt" ]; then pass "c3_$r.fa (+bwa index)" "$(du -sh "$REFS/c3_$r.fa" | cut -f1)"
        else warn "c3_$r.fa" "present, NOT bwa-indexed — benchmark_1 will index it"; fi
    else warn "c3_$r.fa" "MISSING — benchmark_1 will fetch it (T6b coverage row)"; fi
done
[ -s "$REFS/chr20.fa.bwt" ] && pass "chr20 bwa index" "present" || warn "chr20 bwa index" "missing — needed for HG002 T6b row"

# ── 6. CODE HEALTH ────────────────────────────────────────────────────────
hdr "6. CODE HEALTH — does it still reproduce its validated numbers?"
say "  Gate: the k-mer multiset AND the SNV F1 must both be unchanged."
say "  k-mer identity is the real gate -- a change can pass an F1 check while"
say "  silently altering kc (it has happened twice in this project)."
say ""
KC_EXPECT=1063607; F1_EXPECT=0.886
SAN_DIR=/tmp/_b0_sanity; rm -rf "$SAN_DIR"; mkdir -p "$SAN_DIR"
if [ -x "$BEST" ] && [ -s "$REFS/chr20.fa" ]; then
    say "  running sanity window (HG002 r2, ~1 min) ..."
    T0=$(date +%s)
    env CAPS_DBG=1 CAPS_DBG_ONLY=1 \
        bash "$HERE/scripts/run_window_bench_capsule.sh" "$BEST" "$HERE/scripts" \
             "$REFS/chr20.fa" HG002 r2 "$SAN_DIR" > "$SAN_DIR/run.log" 2>&1
    T1=$(date +%s)
    KC_GOT=$(grep -oP 'kc nodes=\K[0-9]+' "$SAN_DIR/capsule.log" 2>/dev/null | head -1)
    F1_GOT=$(grep -aE '^SNV ' "$SAN_DIR/run.log" 2>/dev/null | grep -oP 'F1=\K[0-9.]+' | head -1)
    P_GOT=$(grep -aE '^SNV ' "$SAN_DIR/run.log" 2>/dev/null | grep -oP ' P=\K[0-9.]+' | head -1)
    R_GOT=$(grep -aE '^SNV ' "$SAN_DIR/run.log" 2>/dev/null | grep -oP ' R=\K[0-9.]+' | head -1)
    hb "sanity wall time" "$((T1-T0)) s"
    [ "${KC_GOT:-0}" = "$KC_EXPECT" ] && pass "k-mer multiset identity" "kc=$KC_GOT (expected $KC_EXPECT)" \
                                      || fail "k-mer multiset identity" "kc=${KC_GOT:-NONE}, expected $KC_EXPECT"
    [ "${F1_GOT:-0}" = "$F1_EXPECT" ] && pass "SNV F1 unchanged" "F1=$F1_GOT P=$P_GOT R=$R_GOT" \
                                      || fail "SNV F1 changed" "F1=${F1_GOT:-NONE}, expected $F1_EXPECT"
else fail "sanity run" "encoder or chr20.fa unavailable — cannot verify code health"; fi

# ── 7. CONFIG DRIFT ───────────────────────────────────────────────────────
hdr "7. CONFIGURATION — the exact values that will run"
say "  Every value below is a DEFAULT compiled into the binary. If one has"
say "  drifted from the validated value, every published number is at risk."
say ""
H="$HERE/include/caps_caller.h"
cfg(){ printf "  %-26s %-14s %-10s %s\n" "$1" "$2" "$3" "$4" | tee -a "$OUT"; }
cfg "PARAMETER" "SHIPPED" "VALIDATED" "BASIS"
say "  ----------------------------------------------------------------------------------------"
# These declarations span several lines (env override on one, default on the
# next), so a line-oriented grep cannot see the default and reports "?" --
# which is a FALSE drift warning, and a warning nobody can act on is worse
# than none. Extraction reads the whole file and takes the default after the
# last ':' of the declaration.
chk_default(){ local name="$1" want="$2" basis="$3" got
    got=$(python3 - "$H" "$name" <<'PYEOF'
import re,sys
src=open(sys.argv[1]).read(); name=sys.argv[2]
# Anchor on the DECLARATION, not the first textual match: these names also
# appear inside comment tables (e.g. "HALFW=31   63 bp   TP=336"), and matching
# a comment produced a bogus "?" and a false drift warning.
m=re.search(r'const\s+\w+(?:\s+\w+)?\s+'+re.escape(name)+r'\s*=\s*(.*?);', src, re.S)
if m:
    tail=m.group(1)
    d=re.findall(r':\s*([0-9]+)\s*u?\s*$', tail.strip(), re.S)
    if not d: d=re.findall(r':\s*([0-9]+)', tail)
    print(d[-1] if d else "?")
else: print("?")
PYEOF
)
    cfg "$name" "${got:-?}" "$want" "$basis"
    [ "${got:-x}" = "$want" ] || { warn "config drift: $name" "shipped=${got:-?} validated=$want"; }; }
chk_default "MINC"    2  "structural (drop singletons)"
chk_default "MINQ"    20 "phred convention"
chk_default "MAXPOLY" 1  "swept: beats DiscoSNP++ P=3 held-out"
chk_default "HALFW"   31 "held-out rejected the derived 44"
say "  ----------------------------------------------------------------------------------------"
for v in CAPS_DBG_IBFS CAPS_KC_FREQMIN CAPS_DBG_STR CAPS_DBG_SB CAPS_DBG_INDEL; do
    if [ -n "${!v:-}" ]; then warn "env override active: $v=${!v}" "this is NOT the validated configuration"
    else hb "env $v" "unset (correct)"; fi
done
say ""
say "  L3 memory ceiling  : 60% of measured MemAvailable, spill auto-enables"
say "  L15 parallelism    : traversal parallel (was serial under CAPS_DBG_SB)"
say "  L2 freq minimizers : OFF (measured +14% volume — rejected)"
say "  spill format       : superkmer (key-only spill silently drops 30% of k-mers)"

# ── 8. PLAN ───────────────────────────────────────────────────────────────
hdr "8. WHAT BENCHMARK 1 WILL RUN"
say "  PHASE 1 — Claim 1: all 19 datasets, ONE AT A TIME"
say "      per dataset: encode(8-candidate adaptive sweep, concurrent) -> archive KEPT"
say "                   decode -> lossless compare -> scratch deleted"
say "                   SPRING compress+decompress, Genozip compress+decompress"
say "      tables: T1 archive size | T2 wall time + peak RAM (all 3 tools)"
say "      est: ~1-2 h (encoder is 5.4x faster as of 9e39e23)  archives kept: ~7 GB"
say ""
say "  PHASE 2 — Claim 2: 4 GIAB human sets only"
say "      per set: our caller (compress+call, one pass) -> DiscoSNP++ -> Kmer2SNP"
say "               -> rtg vcfeval against GIAB truth"
say "      tables: T3 het-SNV F1 (3-way) | T4 coverage sweep 10/15/20/30x"
say "      est: 3-4 h"
say ""
say "  PHASE 3 — Claim 3: reads the archives Phase 1 kept"
say "      T6a export   vs SPAdes            — 6 datasets (one per kingdom)"
say "      T6b coverage vs bwa+samtools+mosdepth — same 6"
say "      T6c query    — ALL 19 (no competitor exists)"
say "      est: 2-2.5 h"
say ""
say "  TOTAL ESTIMATE: 6-9 h sequential. Peak transient disk ~35 GB, peak RAM ~20 GB."
say "  NOTE: timings are ESTIMATES from two anchors (E. coli encode 32.46 s,"
say "        HG002 chr20 25 min). Phase 1's first dataset replaces them with"
say "        a measured rate."

# ── 9. VERDICT ────────────────────────────────────────────────────────────
T_END=$(date +%s)
hdr "9. VERDICT"
hb "checks passed"  "$OK"
hb "warnings"       "$WARN"
hb "failures"       "$FAIL"
hb "preflight took" "$((T_END-T_START)) s"
say ""
if [ "$FAIL" -eq 0 ]; then
    say "  ╔════════════════════════════════════════════════════════════════════╗"
    say "  ║   VERDICT: GO                                                      ║"
    say "  ║   Every required tool, dataset, reference and truth set is present.║"
    say "  ║   The code reproduces its validated k-mer set and F1.              ║"
    say "  ╚════════════════════════════════════════════════════════════════════╝"
    [ "$WARN" -gt 0 ] && { say ""; say "  $WARN warning(s) above — the run proceeds, with that scope reduced."; }
    say ""
    say "  START WITH:  bash scripts/benchmark_1_run.sh"
    say "  (add SANITY_ONLY=1 to run just the first dataset and get a measured rate)"
    RC=0
else
    say "  ╔════════════════════════════════════════════════════════════════════╗"
    say "  ║   VERDICT: NO-GO                                                   ║"
    say "  ║   $FAIL required item(s) failed. Fix these, then re-run preflight.  ║"
    say "  ╚════════════════════════════════════════════════════════════════════╝"
    say ""; say "  Failures:"; grep '\[FAIL\]' "$OUT" | sed 's/^/  /'
    RC=1
fi
say ""
say "  status file: $OUT"
say "  finished   : $(date '+%Y-%m-%d %H:%M:%S')"
exit $RC
