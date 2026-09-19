#!/usr/bin/env bash
# PROBE LENGTH vs MISMATCH TOLERANCE, on real HG002 reads.
#
# WHY. Raising CAPS_QUERY_MM above 2 changed nothing on a 40 bp probe, and that
# is not the hypothesis failing -- it is the seed floor. Pigeonhole search needs
# k+1 seeds of at least MINSEED=12 bp, so a probe of length P supports at most
#     k_max = floor(P/12) - 1
# and a 40 bp probe is therefore capped at k=2 no matter what k is requested.
# The nine residual sites need MORE than two mismatches, so the only way to
# reach them is a LONGER probe.
#
# This is not free, and the sweep is built to expose the cost rather than hide
# it. A longer probe also spans more sequence, so it meets more variants and
# more sequencing error. Whether it wins is an empirical question about how
# k_max ~ P/12 races against expected mismatches ~ P x density. Both sets are
# swept at every length, and the number that matters is the SEPARATION between
# real het sites and homozygous controls.
set -u
cd ~/t34real
export CAPS_SPANS=1
DEC=~/bin/capsule_decode
ARC=realreads.capsule
PG=pg_real.fa
SC=~/gc/scripts/score_locus_fidelity_v2.py
N="${N:-400}"
REF0=2999001

build(){   # $1 mode(het|neg)  $2 probe length -> stdout
python3 - "$N" "$REF0" "$1" "$PG" "$2" <<'PY'
import sys, random
N=int(sys.argv[1]); REF0=int(sys.argv[2]); MODE=sys.argv[3]
PGF=sys.argv[4]; P=int(sys.argv[5])
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open(PGF) if l[0]!='>')
def emit(pos,r,a):
    lo=pos-5-P-REF0; hi=pos-5-REF0          # probe ends 5 bp before the variant
    if lo<0 or hi>len(ref): return None
    p=ref[lo:hi]
    if len(p)!=P or set(p)-set("ACGT"): return None
    at=pg.find(p)
    if at<0: return None
    return f"{pos}\t{r}\t{a}\t{p}\t{at+P+5}"
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

query(){   # $1 probes  $2 outdir  $3 k
  rm -rf "$2"; mkdir -p "$2"
  cat "$1" | xargs -P 6 -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r POS R A P V <<< "{}"
    CAPS_QUERY_MM=0 '"$DEC"' query '"$ARC"' '"$2"'/s0_$POS.fa "$P" >/dev/null 2>'"$2"'/s0_$POS.occ
    CAPS_QUERY_MM='"$3"' '"$DEC"' query '"$ARC"' '"$2"'/s1_$POS.fa "$P" >/dev/null 2>'"$2"'/s1_$POS.occ
    '"$DEC"' query '"$ARC"' '"$2"'/c_$POS.fa $((V-100))-$((V+100)) >/dev/null 2>&1'
}

score(){   # $1 probes  $2 qdir -> "both" of the mm+offset arm (field 16)
  cp "$1" probes.tsv; rm -rf q; cp -r "$2" q
  python3 "$SC" X 2 | grep '^T34V2' | cut -d, -f16
}

printf '%6s %5s %8s %14s %22s %12s\n' probe k sites "het both" "homozygous false" separation
for P in 40 60 80 100 120; do
  K=$(( P/12 - 1 ))
  build het "$P" > ph_$P.tsv
  build neg "$P" > pn_$P.tsv
  NH=$(wc -l < ph_$P.tsv)
  query ph_$P.tsv "qp$P" "$K"
  query pn_$P.tsv "qpn$P" "$K"
  H=$(score ph_$P.tsv "qp$P")
  G=$(score pn_$P.tsv "qpn$P")
  printf '%6s %5s %8s %14s %22s %12s\n' "$P" "$K" "$NH" "$H" "$G" "$((H-G))"
done
