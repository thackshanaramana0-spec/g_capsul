#!/bin/bash
# ============================================================================
# G_CAPSUL — ONE end-to-end benchmark. Compress, then serve every claim from
# the archive alone, reporting SIZE / SPEED / RAM / DISK per stage.
#
#   bash scripts/benchmark_final.sh <reads.fq> [outdir]
#
# There is nothing else to run and nothing to remember: this builds its own
# binaries, encodes, verifies losslessness, and exercises Claim 1 (size),
# Claim 2 (calling from the archive) and Claim 3 (export/coverage/query).
# Every number it prints is measured in this run.
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN="${1:?usage: benchmark_final.sh <reads.fq> [outdir]}"
W="${2:-$(mktemp -d)}"; mkdir -p "$W"
[ -r "$IN" ] || { echo "FATAL: cannot read $IN" >&2; exit 1; }

# --- preflight: disk. The k-mer counter spills, and a full disk CORRUPTS
# --- rather than fails (see ARCHIVE_CALLER_SPEED_RAM.md section 6).
NEED_GB=40
FREE_GB=$(df -BG --output=avail "$W" | tail -1 | tr -dc 0-9)
if [ "${FREE_GB:-0}" -lt "$NEED_GB" ]; then
  echo "FATAL: ${FREE_GB}GB free at $W, need >=${NEED_GB}GB." >&2
  echo "       A full disk silently truncates the k-mer spill and produces" >&2
  echo "       WRONG results rather than an error. Refusing to run." >&2
  exit 1
fi

R="$W/BENCHMARK.txt"; : > "$R"
say(){ echo "$*" | tee -a "$R"; }
IN_BYTES=$(stat -c%s "$IN"); IN_READS=$(( $(wc -l < "$IN") / 4 ))

say "=============================================================="
say " G_CAPSUL end-to-end benchmark"
say " input   : $IN"
say "           $(numfmt --to=iec $IN_BYTES 2>/dev/null || echo $IN_BYTES)B, $IN_READS reads"
say " host    : $(nproc) vCPU, $(free -g | awk '/^Mem:/{print $2}')GB RAM, ${FREE_GB}GB free disk"
say " started : $(date -Is)"
say "=============================================================="

T0=$(date +%s.%N)
# ---------------------------------------------------------------- build ----
say ""
say "[0/5] build"
B0=$(date +%s.%N)
bash "$HERE/scripts/build106.sh"      "$W/enc" >/dev/null 2>&1 || { say "  FATAL: encoder build failed"; exit 1; }
bash "$HERE/scripts/build_decode.sh"  "$W/dec" >/dev/null 2>&1 || { say "  FATAL: decoder build failed"; exit 1; }
say "  encoder + decoder built in $(echo "$(date +%s.%N) $B0" | awk '{printf "%.1f",$1-$2}')s"

# ------------------------------------------------------------- compress ----
# DUMP_PERM=1 is REQUIRED: pos_abs/pos_strand/read_lengths live behind it, and
# an archive without them cannot drive `call`. See the memory note on this flag.
say ""
say "[1/5] COMPRESS  (Claim 1: size)"
/usr/bin/time -f "%e %M" -o "$W/t.enc" env \
  CAPS_CALL=1 CAPS_NAMES=1 CAPS_QUAL=1 DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 \
  ARCHIVE="$W/a.capsule" "$W/enc" "$IN" 3 16 16 22 16 16 1 24 64 1 \
  >/dev/null 2>"$W/enc.log"
ARC=$(stat -c%s "$W/a.capsule" 2>/dev/null || echo 0)
[ "$ARC" -gt 0 ] || { say "  FATAL: no archive produced (see $W/enc.log)"; exit 1; }
ENC_S=$(awk '{print $1}' "$W/t.enc"); ENC_M=$(awk '{printf "%.0f",$2/1024}' "$W/t.enc")
say "  archive        : $ARC B   ratio $(echo "$IN_BYTES $ARC" | awk '{printf "%.2fx", $1/$2}')   $(echo "$IN_READS $ARC" | awk '{printf "%.3f bits/base-ish (%.1f B/read)", $2*8/($1*150), $2/$1}')"
say "  time / peak RAM: ${ENC_S}s / ${ENC_M} MB"
say "  per-stream (top 6):"
grep -oE '\[archive\] [a-z_0-9]+ +[0-9]+ B' "$W/enc.log" | sort -k3 -rn | head -6 \
  | awk -v t="$ARC" '{printf "     %-14s %10s B  %5.1f%%\n",$2,$3,100*$3/t}' | tee -a "$R"

# ------------------------------------------------------------- lossless ----
say ""
say "[2/5] VERIFY    (round trip from the archive alone)"
V0=$(date +%s.%N); mkdir -p "$W/rt"
"$W/dec" "$W/a.capsule" "$W/rt" "$W/rt/reads.seq" >"$W/dec.log" 2>&1
awk 'NR%4==2' "$IN" > "$W/orig.seq"
if cmp -s "$W/orig.seq" "$W/rt/reads.seq"; then LL="LOSSLESS"; else LL="*** LOSSY ***"; fi
say "  decompress     : $(echo "$(date +%s.%N) $V0" | awk '{printf "%.1f",$1-$2}')s   $(wc -l < "$W/rt/reads.seq") reads"
say "  sequence       : $LL"

# ------------------------------------------------- Claim 2: call from archive
say ""
say "[3/5] CALL      (Claim 2: variants from the archive, no FASTQ)"
/usr/bin/time -f "%e %M" -o "$W/t.call" env CAPS_CALL_INDELS=1 \
  "$W/dec" call "$W/a.capsule" "$W/calls.vcf" "$W/callwk" >/dev/null 2>"$W/call.log"
say "  time / peak RAM: $(awk '{print $1}' "$W/t.call")s / $(awk '{printf "%.0f",$2/1024}' "$W/t.call") MB"
say "  VCF records    : $(grep -vc '^#' "$W/calls.vcf" 2>/dev/null || echo 0)"
say "  placements     : $(grep -oE '[0-9]+/[0-9]+ read placements rebuilt' "$W/call.log" | head -1)"
say "  stage flow:"
grep -hE 'call-timing|CAPS-CALL-TIMING' "$W/call.log" \
  | sed -E 's/^\s*//; s/\[[a-z-]+\] //' | awk '{printf "     %s\n",$0}' | head -8 | tee -a "$R"
rm -rf "$W/callwk"

# ------------------------------------------- Claim 3: addressable operations
say ""
say "[4/5] SERVE     (Claim 3: addressable, all from the archive)"
for OP in export coverage query; do
  ARG=""; [ "$OP" = query ] && ARG="0-100000"
  /usr/bin/time -f "%e %M" -o "$W/t.$OP" "$W/dec" "$OP" "$W/a.capsule" "$W/out.$OP" $ARG \
    >"$W/$OP.log" 2>&1
  say "  $(printf '%-9s' $OP) $(awk '{print $1}' "$W/t.$OP")s / $(awk '{printf "%.0f",$2/1024}' "$W/t.$OP") MB   out $(stat -c%s "$W/out.$OP") B   pg_rebuilds=$(grep -c 'pg rebuilt' "$W/$OP.log")"
done
# correctness of coverage against an identity that a single dropped read breaks
COV=$(awk 'NR>1{t+=($3-$2)*$4} END{print t+0}' "$W/out.coverage")
LEN=$(awk 'NR%4==2{t+=length($0)} END{print t+0}' "$IN")
say "  coverage identity: $COV covered bases vs $LEN read bases -> $([ "$COV" = "$LEN" ] && echo EXACT || echo "DIFF $((LEN-COV))")"

# ------------------------------------------------------------------ disk ----
say ""
say "[5/5] DISK"
say "  input          : $IN_BYTES B"
say "  archive        : $ARC B"
say "  peak scratch   : $(du -sb "$W" | cut -f1) B  (includes the verification decompress)"
say "  archive alone is what ships; everything else in $W is reproducible."

TOT=$(echo "$(date +%s.%N) $T0" | awk '{printf "%.1f",$1-$2}')
say ""
say "=============================================================="
say " TOTAL wall: ${TOT}s   |   archive $ARC B   |   sequence $LL"
say " full report: $R"
say "=============================================================="
