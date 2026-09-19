#!/usr/bin/env bash
# T3.4-v2 on REAL HG002 READS.
#
# Everything is real here: the reads come from the GIAB 300x Illumina BAM,
# range-fetched for chr20:2,999,001-3,601,000 and subsampled to ~30x, so they
# carry real error profiles, real coverage bias and real duplicates. The
# variants are real GIAB truth, the reference is real GRCh37, and the encoder,
# archive, sidecar, query and scorer are the shipped code.
#
# Positive set : 400 real GIAB heterozygous SNVs in the window.
# Negative set : 400 positions with NO GIAB variant within +/-200 bp. With real
#                reads there are no simulated haplotypes to consult, so GIAB's
#                own benchmark calls are the only statement of homozygosity --
#                which is the correct authority anyway.
set -u
cd ~/t34real
DEC=~/bin/capsule_decode
ARC=realreads.capsule
PG=pg_real.fa
N="${N:-400}"
MM="${MM:-2}"
PAR="${PAR:-6}"
REF0=2999001

export CAPS_SPANS=1

build_probes(){   # $1 = mode (het|neg), writes to stdout
python3 - "$N" "$REF0" "$1" "$PG" <<'PY'
import sys, random
N=int(sys.argv[1]); REF0=int(sys.argv[2]); MODE=sys.argv[3]; PGF=sys.argv[4]
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open(PGF) if l[0]!='>')
def emit(pos,r,a):
    lo=pos-45-REF0; hi=pos-5-REF0
    if lo<0 or hi>len(ref): return None
    p=ref[lo:hi]
    if len(p)!=40 or set(p)-set("ACGT"): return None
    at=pg.find(p)
    if at<0: return None
    return f"{pos}\t{r}\t{a}\t{p}\t{at+45}"
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
    cands=list(range(REF0+1000, REF0+len(ref)-1000)); rng.shuffle(cands)
    for pos in cands:
        if near(pos): continue
        r=ref[pos-REF0]
        if r not in "ACGT": continue
        a="ACGT"[("ACGT".index(r)+1)%4]      # a base that is NOT the reference
        s=emit(pos,r,a)
        if s: print(s); n+=1
        if n>=N: break
PY
}

run_set(){        # $1 = probes file, $2 = query dir
  rm -rf "$2" && mkdir -p "$2"
  cat "$1" | xargs -P "$PAR" -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r POS R A P V <<< "{}"
    CAPS_QUERY_MM=0   '"$DEC"' query '"$ARC"' '"$2"'/s0_$POS.fa "$P" >/dev/null 2>'"$2"'/s0_$POS.occ
    CAPS_QUERY_MM='"$MM"' '"$DEC"' query '"$ARC"' '"$2"'/s1_$POS.fa "$P" >/dev/null 2>'"$2"'/s1_$POS.occ
    '"$DEC"' query '"$ARC"' '"$2"'/c_$POS.fa $((V-100))-$((V+100)) >/dev/null 2>&1'
}

build_probes het > probes_real.tsv
build_probes neg > probes_real_neg.tsv
echo "het probes: $(wc -l < probes_real.tsv)   negative probes: $(wc -l < probes_real_neg.tsv)"

run_set probes_real.tsv     qreal
run_set probes_real_neg.tsv qrealneg

SC=~/gc/scripts/score_locus_fidelity_v2.py
cp probes_real.tsv probes.tsv;     rm -rf q && cp -r qreal q
python3 "$SC" HG002_REAL_READS "$MM"
cp probes_real_neg.tsv probes.tsv; rm -rf q && cp -r qrealneg q
python3 "$SC" REAL_NEGATIVE_CONTROL "$MM"
