#!/usr/bin/env bash
# T1 (archive size) + T2 (compress/decompress time + peak RAM), whole-file
# (sequence + names + quality + line 3), CAPSULE vs SPRING vs Genozip, on the
# 14 locked datasets. Adapted from the outer ARCS project's proven
# benchmark/run_block1.sh, using THIS repo's own binaries
# (best106/capsule_decode) in place of the outer `arcs` binary.
#
#   usage: run_claim1_bench.sh [DATA_DIR] [OUT_DIR]
#   DATA_DIR default: $CAPSULE_DATA_DIR (scripts/capsule_config.sh)
#   OUT_DIR  default: ./results/claim1_bench
#
# Tools expected on PATH: spring, genozip, genounzip
# NEVER run this alongside another timed benchmark job -- contaminates
# wall-clock and RSS measurements (CLAUDE.md rule 3).
#
# LOSSLESS CHECK, stated precisely (this is the one thing to get right):
# the FASTQ '+' line is spec-redundant (the record ID already appears on the
# '@' line), so per this project's own established convention
# (/root/arcs-clean/CLAUDE.md rule 5: "the lossless check normalizes the '+'
# line on both sides"), this script does the same: both the original file and
# the reassembled decode output have their 3rd line replaced with a bare '+'
# before comparison. This is NOT a workaround for a defect in CAPSULE's own
# round trip -- it is the same fair standard SPRING is held to, since SPRING
# strips '+'-line IDs entirely and reports LOSSY on a literal byte diff.
#
# ONE REAL GAP FOUND while writing this, disclosed rather than hidden: the
# archive DOES store enough information to reconstruct the ORIGINAL exact
# '+' line (a `line3_mode` byte -- 0=bare '+', 1='+'+full header repeat,
# 2=arbitrary/must-store -- set in include/names_coder.h's encode path), but
# `capsule_decode.cpp` never reads this byte back out or uses it when writing
# the decoder's own outputs. So an EXACT (non-normalized) full-file byte
# reassembly is not actually achievable via decoder output alone today for
# any dataset using line3_mode=1 or 2 -- confirmed by inspection: SRR2584863,
# ERR552797, and ERR5181310 (3 of the 14 locked datasets) all use the
# '+'+header-repeat convention, not bare '+'. This means the wording in
# README.md/docs/TECHNICAL_ARCHITECTURE.md ("verified byte-identical, same
# MD5, to the original input") is broader than what the current decoder can
# support for most real datasets -- true only for bare-'+' files. The
# NORMALIZED check this script performs is still a complete, honest,
# spec-compliant lossless proof; it is not the same claim as an exact MD5
# match of the whole original file, and those two claims should not be
# conflated in the paper.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$HERE/scripts/capsule_config.sh"
DATA_DIR="${1:-$CAPSULE_DATA_DIR}"
OUT_DIR="${2:-$HERE/results/claim1_bench}"
mkdir -p "$OUT_DIR"
WD="$CAPSULE_BIN_DIR/claim1_wd"; mkdir -p "$WD"
BEST="$CAPSULE_BIN_DIR/best106"
DEC="$CAPSULE_BIN_DIR/capsule_decode"
NPROC=$(nproc 2>/dev/null || echo 4)

DATASETS="ERR5181310 SRR554369 ERR552797 SRR2584863 SRR29296997 ERR12954017 \
SRR40402583 SRR40271341 ERR17740259 SRR37283774 DRR976266 SRR36741279 \
SRR32429602 SRR39257532"

_PHASE=""
phase()  { _PHASE="$1"; echo ""; echo "[Phase $_PHASE] $2"; }
pdone()  { echo "[Phase $_PHASE] DONE — $*"; }
pfail()  { echo "[Phase $_PHASE] FAIL — $*"; exit 1; }
pskip()  { echo "[Phase $_PHASE] SKIP — $*"; }
pinfo()  { echo "[Phase $_PHASE]   -> $*"; }

parse_time_v() {
    local logf="$1" wall vmhwm
    wall=$(grep "Elapsed (wall clock)" "$logf" | awk '{
        n=split($NF,a,":");
        if(n==3) printf "%.2f", a[1]*3600+a[2]*60+a[3];
        else if(n==2) printf "%.2f", a[1]*60+a[2];
        else printf "%.2f", a[1];
    }')
    vmhwm=$(grep "Maximum resident set size" "$logf" | awk '{print $NF}')
    echo "${wall:-0} ${vmhwm:-0}"
}

# Normalized whole-file lossless compare: reassemble CAPSULE's 3 decoder
# outputs (sequence/.names/.qual, one read per line each -- names file has NO
# '@' prefix, decode_to_file() is called with at=false, see
# capsule_decode.cpp:604) into the same tab-joined 4-field-per-record shape
# `paste - - - -` produces from a real 4-line FASTQ, with a literal '+' as
# field 3 (this project's established redundant-'+' normalization, not a
# workaround -- see this script's own header comment). Uses `paste`'s default
# tab join across exactly 3 files, not a fragile multi-delimiter interleave.
# $1=original.fq  $2=decoded_seq  $3=decoded_names  $4=decoded_qual
losscmp_capsule() {
    paste <(sed 's/^/@/' "$3") "$2" "$4" \
        | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,"+",$3}' | sort > /tmp/_dec_$$
    paste - - - - < "$1" | awk 'BEGIN{FS=OFS="\t"}{$3="+";print}' | sort > /tmp/_orig_$$
    local res=LOSSY
    cmp -s /tmp/_orig_$$ /tmp/_dec_$$ && res=LOSSLESS
    rm -f /tmp/_orig_$$ /tmp/_dec_$$
    echo "$res"
}
# Same normalization, for SPRING/Genozip's own decompressed output (already a
# real 4-line file — just normalize '+' on both sides).
losscmp_plain() {
    paste - - - - < "$1" | tr -d '\r' | awk 'BEGIN{FS=OFS="\t"}{$3="+";print}' | sort > /tmp/_o_$$
    paste - - - - < "$2" | tr -d '\r' | awk 'BEGIN{FS=OFS="\t"}{$3="+";print}' | sort > /tmp/_d_$$
    local res=LOSSY
    cmp -s /tmp/_o_$$ /tmp/_d_$$ && res=LOSSLESS
    rm -f /tmp/_o_$$ /tmp/_d_$$
    echo "$res"
}

phase "0" "Prerequisites"
[ -x "$BEST" ] || { pinfo "building encoder"; bash "$HERE/scripts/build106.sh" "$BEST" >/dev/null; }
[ -x "$DEC" ]  || { pinfo "building decoder"; bash "$HERE/scripts/build_decode.sh" "$DEC" >/dev/null; }
for t in spring genozip genounzip; do
    command -v "$t" >/dev/null || pfail "$t not on PATH -- see docs/SERVER_SETUP_AND_DOWNLOADS.md sec 5"
done
pdone "encoder/decoder built, spring/genozip present"

RESULTS_CSV="$OUT_DIR/t1_t2_results.csv"
echo "dataset,tool,archive_bytes,compress_s,decompress_s,peak_ram_kb,lossless" > "$RESULTS_CSV"

DS_IDX=0
for DS in $DATASETS; do
    DS_IDX=$((DS_IDX+1))
    SRC="$DATA_DIR/${DS}_1.fq"
    [ -s "$SRC" ] || { pskip "$DS: $SRC not found"; continue; }
    IN="$WD/$DS.fq"; cp "$SRC" "$IN"
    RAW=$(stat -c %s "$IN")
    phase "${DS_IDX}.0" "$DS ($RAW B)"

    # ---- CAPSULE: full-file scope (names + quality), fixed candidate ------
    phase "${DS_IDX}.1" "[$DS] CAPSULE compress"
    A="$WD/$DS.capsule"
    TF="$WD/t_caps_c_$$"
    /usr/bin/time -v env CAPS_NAMES=1 CAPS_QUAL=1 INPUT="$IN" ARCHIVE="$A" BEST="$BEST" \
        bash "$HERE/scripts/encode_adaptive.sh" >/dev/null 2>"$TF"
    read -r CWALL CRAMKB <<< "$(parse_time_v "$TF")"
    ARCH=$(stat -c %s "$A")
    pinfo "compress: archive=${ARCH}B time=${CWALL}s ram=${CRAMKB}KB"

    phase "${DS_IDX}.2" "[$DS] CAPSULE decompress"
    OUTDIR="$WD/${DS}_dec"; mkdir -p "$OUTDIR"
    TF2="$WD/t_caps_d_$$"
    /usr/bin/time -v "$DEC" "$A" "$OUTDIR" "$OUTDIR/reads.seq" >/dev/null 2>"$TF2"
    read -r DWALL DRAMKB <<< "$(parse_time_v "$TF2")"
    LL=$(losscmp_capsule "$IN" "$OUTDIR/reads.seq" "$OUTDIR/reads.seq.names" "$OUTDIR/reads.seq.qual")
    rm -f "$TF" "$TF2" "$A"; rm -rf "$OUTDIR"
    printf "%s,CAPSULE,%s,%.2f,%.2f,%s,%s\n" "$DS" "$ARCH" "$CWALL" "$DWALL" "$CRAMKB" "$LL" >> "$RESULTS_CSV"
    pdone "CAPSULE — archive=${ARCH}B ctime=${CWALL}s dtime=${DWALL}s $LL"
    [ "$LL" = "LOSSLESS" ] || pfail "$DS CAPSULE LOSSY -- halt per CLAUDE.md rule 5, do not proceed"

    # ---- SPRING (needs -g on decompress, per CLAUDE.md) -------------------
    phase "${DS_IDX}.3" "[$DS] SPRING"
    A="$WD/$DS.spring"; OUTFQ="$WD/$DS.spring.dec"
    TF="$WD/t_spr_c_$$"
    /usr/bin/time -v spring -c -i "$IN" -o "$A" -t "$NPROC" -g 2>"$TF"
    read -r CWALL CRAMKB <<< "$(parse_time_v "$TF")"
    ARCH=$(stat -c %s "$A")
    TF2="$WD/t_spr_d_$$"
    /usr/bin/time -v spring -d -i "$A" -o "$OUTFQ" -t "$NPROC" -g 2>"$TF2"
    read -r DWALL DRAMKB <<< "$(parse_time_v "$TF2")"
    LL=$(losscmp_plain "$IN" "$OUTFQ")
    rm -f "$TF" "$TF2" "$A" "$OUTFQ"
    printf "%s,SPRING,%s,%.2f,%.2f,%s,%s\n" "$DS" "$ARCH" "$CWALL" "$DWALL" "$CRAMKB" "$LL" >> "$RESULTS_CSV"
    pdone "SPRING — archive=${ARCH}B ctime=${CWALL}s dtime=${DWALL}s $LL"

    # ---- Genozip ------------------------------------------------------------
    phase "${DS_IDX}.4" "[$DS] Genozip"
    A="$WD/$DS.genozip"; OUTFQ="$WD/$DS.gz.dec"
    TF="$WD/t_gz_c_$$"
    /usr/bin/time -v genozip --force -o "$A" "$IN" 2>"$TF"
    read -r CWALL CRAMKB <<< "$(parse_time_v "$TF")"
    ARCH=$(stat -c %s "$A")
    TF2="$WD/t_gz_d_$$"
    /usr/bin/time -v genounzip --force -o "$OUTFQ" "$A" 2>"$TF2"
    read -r DWALL DRAMKB <<< "$(parse_time_v "$TF2")"
    LL=$(losscmp_plain "$IN" "$OUTFQ")
    rm -f "$TF" "$TF2" "$A" "$OUTFQ"
    printf "%s,Genozip,%s,%.2f,%.2f,%s,%s\n" "$DS" "$ARCH" "$CWALL" "$DWALL" "$CRAMKB" "$LL" >> "$RESULTS_CSV"
    pdone "Genozip — archive=${ARCH}B ctime=${CWALL}s dtime=${DWALL}s $LL"

    rm -f "$IN"
    phase "${DS_IDX}.5" "[$DS] complete"
done

echo ""
echo "=== SUMMARY ==="
cat "$RESULTS_CSV"
echo ""
echo "Results: $RESULTS_CSV"
echo "Cross-validate: SPRING bpb for SRR554369 (expect 0.2416) and SRR870667"
echo "(expect 1.2621, not in the 14-set) per CLAUDE.md's cross-validation rule."
