#!/usr/bin/env bash
# NEGATIVE CONTROL for the T3.4-v2 real run.
#
# The positive result (400/400) only shows the mm+offset configuration finds
# two alleles at sites that HAVE two alleles. That is worthless on its own
# unless the same configuration REFUSES to find two alleles where there is only
# one. Otherwise "both" would just mean "we looked".
#
# So: 400 HOMOZYGOUS-REFERENCE positions -- positions with no GIAB variant
# within +/-200 bp, and identical on both simulated haplotypes by construction.
# A correct method scores these as "one" (reference only). Any "both" here is a
# false heterozygous call, and the ALT count is the false-positive rate.
#
# Scored with the SAME unmodified scorer.
set -u
cd ~/t34real
DEC=~/bin/capsule_decode
ARC=real.capsule
N="${N:-400}"
MM="${MM:-2}"
PAR="${PAR:-6}"
REF0=2999001

rm -rf qneg && mkdir -p qneg

python3 - "$N" "$REF0" > probes_neg.tsv <<'PY'
import sys, random
N=int(sys.argv[1]); REF0=int(sys.argv[2])
ref=open("ref.txt").read().strip()
hap1=open("hap1.txt").read().strip()
hap2=open("hap2.txt").read().strip()
pg="".join(l.strip() for l in open("pg.fa") if l[0]!='>')
var=set()
for line in open("giab_win.vcf"):
    if line.startswith("#"): continue
    f=line.split("\t"); var.add(int(f[1]))
def near(p):
    return any((p+d) in var for d in range(-200,201))
rng=random.Random(4242)
cands=list(range(REF0+1000, REF0+len(ref)-1000))
rng.shuffle(cands)
n=0
for pos in cands:
    if near(pos): continue
    i=pos-REF0
    # must be truly homozygous in the simulated sample
    if hap1[i]!=hap2[i]: continue
    r=ref[i]
    if r not in "ACGT": continue
    # a fake "alt" that is simply the other 3 bases -- pick one deterministically
    a="ACGT"[("ACGT".index(r)+1)%4]
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
echo "negative-control sites: $(wc -l < probes_neg.tsv)"

export CAPS_SPANS=1
cat probes_neg.tsv | xargs -P "$PAR" -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r POS R A P V <<< "{}"
  CAPS_QUERY_MM=0   '"$DEC"' query '"$ARC"' qneg/s0_$POS.fa "$P" >/dev/null 2>qneg/s0_$POS.occ
  CAPS_QUERY_MM='"$MM"' '"$DEC"' query '"$ARC"' qneg/s1_$POS.fa "$P" >/dev/null 2>qneg/s1_$POS.occ
  '"$DEC"' query '"$ARC"' qneg/c_$POS.fa $((V-100))-$((V+100)) >/dev/null 2>&1'

# the scorer reads probes.tsv and q/ -- point it at the negative set
rm -rf q_pos_backup && mv q q_pos_backup && mv probes.tsv probes_pos_backup.tsv
mv qneg q && cp probes_neg.tsv probes.tsv
python3 ~/gc/scripts/score_locus_fidelity_v2.py NEGATIVE_CONTROL "$MM"
mv q qneg && mv q_pos_backup q && mv probes_pos_backup.tsv probes.tsv
