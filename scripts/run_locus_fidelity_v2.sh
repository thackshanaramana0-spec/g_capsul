#!/usr/bin/env bash
# T3.4-v2 — LOCUS RETRIEVAL FIDELITY, decomposed into its two failure modes.
#
#   usage: run_locus_fidelity_v2.sh <capsule_decode> <archive> <ref.fa> [IND] [outdir] [N]
#
# WHY THIS EXISTS, AND WHY IT IS A SEPARATE SCRIPT.
# `run_locus_fidelity.sh` produced the published 345/400 and is left untouched
# so that figure stays exactly reproducible. This script asks a different
# question of the same archive and the same sites: WHERE do the residual
# one-allele sites actually lose the alternate allele?
#
# Two candidate stages, and they are independent:
#
#   SEARCH   the probe is matched against the pseudogenome CONSENSUS. At a het
#            locus the two haplotypes are separate contigs, so a reference
#            probe differs from the alternate contig wherever a SECOND variant
#            falls inside the 40 bp window. Exact matching then resolves only
#            the haplotype that agrees with the probe. (CAPS_QUERY_MM=k makes
#            the search tolerate k mismatches, pigeonhole seed-and-verify --
#            the same pattern the encoder already uses to place reads.)
#
#   SCORING  the original scorer requires each returned read to CONTAIN the
#            reference probe (`s.find(probe)`). A read carrying the variant may
#            also differ from the probe elsewhere, so containment discards
#            exactly the reads that matter. The COORDINATE arm never had this
#            restriction -- it reads the base at `off = v - p` from the read's
#            stored placement -- so the two arms were being scored by different
#            rules on the same archive. Offset scoring removes that asymmetry.
#
# Both are reported, in a 2x2, so the effect of each is separable and this
# experiment can DISPROVE the hypothesis rather than only confirm it. The
# prediction under test, from a 400-site simulation of the algorithm:
#
#     exact  + containment   ~= the published number
#     tolerant + containment ~= unchanged      (search alone does nothing)
#     exact  + offset        ~= unchanged      (scoring alone does nothing)
#     tolerant + offset      >> both           (synergistic, not additive)
#
# If the real archive shows the two corrections are additive, or that either
# alone recovers the residue, the simulation's mechanism is wrong and this
# table says so.
set -u
DEC="${1:?usage: run_locus_fidelity_v2.sh <capsule_decode> <archive> <ref.fa> [IND] [outdir] [N]}"
ARC="${2:?}"; REF="${3:?}"
IND="${4:-HG002}"; OUT="${5:-/tmp/t34v2_$IND}"; N="${6:-100}"
TRUTH="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
CHROM="${T34_CHROM:-20}"; LO="${T34_LO:-3000000}"; HI="${T34_HI:-3600000}"
MM="${T34_MM:-2}"                  # mismatch tolerance for the tolerant arm
PAR="${T34_PAR:-6}"

abspath(){ case "$1" in /*) printf '%s' "$1";; *) printf '%s/%s' "$(pwd)" "$1";; esac; }
# Resolve the scorer BEFORE the cd below -- $0 may be relative, and this script
# cd's into its output directory, which silently broke the lookup otherwise.
SCORER="$(cd "$(dirname "$(abspath "$0")")" && pwd)/score_locus_fidelity_v2.py"
DEC="$(command -v "$DEC" 2>/dev/null || abspath "$DEC")"
ARC="$(abspath "$ARC")"; REF="$(abspath "$REF")"
case "$OUT" in /*) ;; *) OUT="$(abspath "$OUT")";; esac

for t in tabix samtools python3; do command -v $t >/dev/null || { echo "MISSING: $t" >&2; exit 1; }; done
[ -x "$DEC" ] || { echo "MISSING/not executable: $DEC" >&2; exit 1; }
[ -s "$REF" ] || { echo "MISSING ref: $REF" >&2; exit 1; }
[ -s "$ARC" ] || { echo "MISSING archive: $ARC" >&2; exit 1; }
[ -s "$TRUTH" ] || { echo "MISSING truth: $TRUTH" >&2; exit 1; }
mkdir -p "$OUT/q"; cd "$OUT"

echo "=== T3.4-v2 locus fidelity, decomposed — $IND chr$CHROM:$LO-$HI (mm=$MM) ==="
[ -s pg.fa ] || "$DEC" export "$ARC" pg.fa >/dev/null 2>&1
[ -s pg.fa ] || { echo "export failed" >&2; exit 1; }
[ -s ref.txt ] || samtools faidx "$REF" $CHROM:$((LO-1000))-$((HI+1000)) 2>/dev/null \
    | tail -n +2 | tr -d '\n' | tr 'acgt' 'ACGT' > ref.txt
tabix "$TRUTH" $CHROM:$LO-$HI 2>/dev/null \
  | awk -F'\t' '$10~/^0[|\/]1|^1[|\/]0/ && length($4)==1 && length($5)==1 {print $2"\t"$4"\t"$5}' > hets.tsv
echo "  GIAB het SNVs in window: $(wc -l < hets.tsv)"

# Identical probe construction to v1: 40 bp of REFERENCE ending 5 bp before the
# variant, so the probe can never contain the variant and cannot rig the result.
python3 - "$N" "$((LO-1000))" > probes.tsv <<'PY'
import sys
N=int(sys.argv[1]); REF0=int(sys.argv[2])
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg.fa") if l[0]!='>')
n=0
for line in open("hets.tsv"):
    pos,r,a=line.split(); pos=int(pos)
    lo=pos-45-REF0; hi=pos-5-REF0
    if lo<0 or hi>len(ref): continue
    p=ref[lo:hi]
    if len(p)!=40 or set(p)-set("ACGT"): continue
    at=pg.find(p)
    if at<0: continue
    print(f"{pos}\t{r}\t{a}\t{p}\t{at+45}")
    n+=1
    if n>=N: break
PY
echo "  probes usable: $(wc -l < probes.tsv)"

# Three queries per site: sequence-exact, sequence-tolerant, coordinate.
# stderr is KEPT for the sequence arms -- it carries the resolved occurrences,
# which is what makes offset scoring possible without re-deriving them.
cat probes.tsv | xargs -P "$PAR" -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r POS R A P V <<< "{}"
  CAPS_QUERY_MM=0   '"$DEC"' query '"$ARC"' q/s0_$POS.fa "$P" >/dev/null 2>q/s0_$POS.occ
  CAPS_QUERY_MM='"$MM"' '"$DEC"' query '"$ARC"' q/s1_$POS.fa "$P" >/dev/null 2>q/s1_$POS.occ
  '"$DEC"' query '"$ARC"' q/c_$POS.fa $((V-100))-$((V+100)) >/dev/null 2>&1'

[ -s "$SCORER" ] || { echo "MISSING scorer: $SCORER" >&2; exit 1; }
python3 "$SCORER" "$IND" "$MM"
