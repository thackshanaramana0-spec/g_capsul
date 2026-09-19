#!/usr/bin/env bash
# BILATERAL ANCHORING on real HG002 reads.
#
# Measured cause of the nine residual sites: the probe is 40 bp of reference
# ending 5 bp BEFORE the variant, so it anchors upstream only. The pseudogenome
# is built by greedy overlap chaining, so a contig's neighbourhood is an
# artefact of which reads overlapped rather than the genomic neighbourhood. When
# chaining breaks the upstream side, no mismatch tolerance can recover it --
# upstream anchoring succeeds at 0 of 6 failing sites while downstream succeeds
# at 5 of 6.
#
# Fix: query with BOTH a 40 bp window ending 5 bp before the variant and a 40 bp
# window starting 5 bp after it, and take the union of occurrences. Neither
# window contains the variant, so neither can rig the result, and the downstream
# window is exactly as reference-derived as the upstream one always was.
#
# Both arms are run over the het set AND the homozygous control, because a union
# of two searches admits more spurious matches and that cost has to be measured
# rather than assumed away.
set -u
cd ~/t34real
export CAPS_SPANS=1
DEC=~/bin/capsule_decode
ARC=realreads.capsule
MM="${MM:-2}"
N="${N:-400}"
REF0=2999001

build(){   # $1 = het|neg -> pos, ref, alt, upprobe, downprobe
python3 - "$N" "$REF0" "$1" <<'PY'
import sys, random
N=int(sys.argv[1]); REF0=int(sys.argv[2]); MODE=sys.argv[3]
PL=40; GAP=5
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')
def emit(pos,r,a):
    i=pos-REF0
    up   = ref[i-GAP-PL : i-GAP]
    down = ref[i+1+GAP : i+1+GAP+PL]
    if len(up)!=PL or len(down)!=PL: return None
    if set(up)-set("ACGT") or set(down)-set("ACGT"): return None
    # keep the SAME admission rule as the published experiment: the site is
    # usable if the upstream probe anchors somewhere in the pg. Admitting sites
    # on the strength of the new probe would change the denominator.
    if pg.find(up)<0: return None
    return f"{pos}\t{r}\t{a}\t{up}\t{down}"
n=0
if MODE=="het":
    for line in open("sites.tsv"):
        pos,r,a=line.split()
        s=emit(int(pos),r,a)
        if s: print(s); n+=1
        if n>=N: break
else:
    var=set()
    for line in open("giab_win.vcf"):
        if line.startswith("#"): continue
        var.add(int(line.split("\t")[1]))
    def near(p): return any((p+d) in var for d in range(-200,201))
    rng=random.Random(4242)
    c=list(range(REF0+1000, REF0+len(ref)-1000)); rng.shuffle(c)
    for pos in c:
        if near(pos): continue
        r=ref[pos-REF0]
        if r not in "ACGT": continue
        a="ACGT"[("ACGT".index(r)+1)%4]
        s=emit(pos,r,a)
        if s: print(s); n+=1
        if n>=N: break
PY
}

query(){   # $1 probes  $2 outdir
  rm -rf "$2"; mkdir -p "$2"
  cat "$1" | xargs -P 6 -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r POS R A U D <<< "{}"
    CAPS_QUERY_MM='"$MM"' '"$DEC"' query '"$ARC"' '"$2"'/up_$POS.fa   "$U" >/dev/null 2>'"$2"'/up_$POS.occ
    CAPS_QUERY_MM='"$MM"' '"$DEC"' query '"$ARC"' '"$2"'/down_$POS.fa "$D" >/dev/null 2>'"$2"'/down_$POS.occ'
}

build het > bl_het.tsv
build neg > bl_neg.tsv
echo "het sites: $(wc -l < bl_het.tsv)   homozygous controls: $(wc -l < bl_neg.tsv)"
query bl_het.tsv qbl_het
query bl_neg.tsv qbl_neg
python3 ~/gc/scripts/score_bilateral.py
