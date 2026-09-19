#!/usr/bin/env bash
# T3.4-v2 on a REAL archive, in ~/t34real.
#
# Same three queries per site and the SAME scorer as
# scripts/score_locus_fidelity_v2.py -- the scorer is the artifact under test
# and is already unit-tested, so it is used unmodified.
#
# Differences from scripts/run_locus_fidelity_v2.sh, which targets the
# benchmark box: sites come from the locally streamed GIAB window rather than a
# local truth VCF, and the reference slice is already on disk as ref.txt.
set -u
cd ~/t34real
DEC=~/bin/capsule_decode
ARC=real.capsule
N="${N:-400}"
MM="${MM:-2}"
PAR="${PAR:-6}"
REF0=2999001

rm -rf q && mkdir -p q

# Probe construction, identical to v1/v2: 40 bp of REFERENCE ending 5 bp before
# the variant, so the probe can never contain the variant and cannot rig the
# result.
python3 - "$N" "$REF0" > probes.tsv <<'PY'
import sys
N=int(sys.argv[1]); REF0=int(sys.argv[2])
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg.fa") if l[0]!='>')
n=0
for line in open("sites.tsv"):
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
echo "probes usable: $(wc -l < probes.tsv)"

export CAPS_SPANS=1
cat probes.tsv | xargs -P "$PAR" -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r POS R A P V <<< "{}"
  CAPS_QUERY_MM=0   '"$DEC"' query '"$ARC"' q/s0_$POS.fa "$P" >/dev/null 2>q/s0_$POS.occ
  CAPS_QUERY_MM='"$MM"' '"$DEC"' query '"$ARC"' q/s1_$POS.fa "$P" >/dev/null 2>q/s1_$POS.occ
  '"$DEC"' query '"$ARC"' q/c_$POS.fa $((V-100))-$((V+100)) >/dev/null 2>&1'

python3 ~/gc/scripts/score_locus_fidelity_v2.py HG002_chr20_real "$MM"
