#!/usr/bin/env bash
# Stage 2 of a held-out window: everything after the reads are fetched.
#
#   usage: finish_window.sh <workdir> <ref0> [N]
#
# Split from the fetch because the GIAB server throttles range requests and the
# download is by far the longest and most fragile step. If it is interrupted the
# partial win.bam is kept and only this stage is re-run.
#
# Waits for any in-flight fetch to exit before starting, so it can be chained
# directly onto a running download.
set -eu
WD="${1:?usage: finish_window.sh <workdir> <ref0> [N]}"
REF0="${2:?1-based coordinate of ref.txt[0]}"
N="${3:-400}"
MM="${MM:-2}"; PAR="${PAR:-3}"; PL=40; GAP=5
BIN="${BIN:-$HOME/bin}"

while pgrep -f "samtools view -b -F 0x900" >/dev/null 2>&1; do sleep 10; done
cd "$WD"
[ -s win.bam ] || { echo "no win.bam"; exit 1; }
echo "  win.bam: $(stat -c%s win.bam) B"

[ -s reads.fq ] || samtools fastq -n win.bam > reads.fq 2>/dev/null
echo "  reads: $(( $(wc -l < reads.fq) / 4 ))"

awk -F'\t' '$7=="PASS" && length($4)==1 && length($5)==1 && $5~/^[ACGT]$/ {
  split($10,g,":"); if(g[1]=="0|1"||g[1]=="1|0"||g[1]=="0/1"||g[1]=="1/0") print $2"\t"$4"\t"$5 }' \
  giab_win.vcf > sites.tsv
echo "  het SNVs: $(wc -l < sites.tsv)"

export CAPS_SPANS=1 DUMP_PERM=1 DUMP_MM=1 DUMP_LIT=1
[ -s win.capsule ] || { "$BIN/best106" reads.fq 3 16 16 22 16 16 1 24 64 1 > enc.log 2>&1
                        cp out.arcs2 win.capsule; }
echo "  archive: $(stat -c%s win.capsule) B"

rm -rf dec; mkdir -p dec
"$BIN/mp_decode" win.capsule dec dec/reads.out >/dev/null 2>&1
awk 'NR%4==2' reads.fq > orig.txt
cmp -s orig.txt dec/reads.out && echo "  LOSSLESS: yes" || { echo "  LOSSLESS: NO"; exit 1; }
rm -rf dec orig.txt

[ -s pg.fa ] || "$BIN/mp_decode" export win.capsule pg.fa >/dev/null 2>&1
[ -s win.capsule.qidx ] || CAPS_PILEUP=1 "$BIN/mp_decode" index win.capsule win.capsule.qidx 2>&1 | grep '\[index\]'

build(){
python3 - "$N" "$REF0" "$1" "$PL" "$GAP" <<'PY'
import sys, random
N=int(sys.argv[1]); REF0=int(sys.argv[2]); MODE=sys.argv[3]
PL=int(sys.argv[4]); GAP=int(sys.argv[5])
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg.fa") if l[0]!='>')
def emit(pos,r,a):
    i=pos-REF0
    up=ref[i-GAP-PL:i-GAP]; down=ref[i+1+GAP:i+1+GAP+PL]
    if len(up)!=PL or len(down)!=PL: return None
    if set(up)-set("ACGT") or set(down)-set("ACGT"): return None
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
    for l in open("giab_win.vcf"):
        if not l.startswith("#"): var.add(int(l.split("\t")[1]))
    rng=random.Random(4242)
    c=list(range(REF0+1000, REF0+len(ref)-1000)); rng.shuffle(c)
    for pos in c:
        if any((pos+d) in var for d in range(-200,201)): continue
        r=ref[pos-REF0]
        if r not in "ACGT": continue
        a="ACGT"[("ACGT".index(r)+1)%4]
        s=emit(pos,r,a)
        if s: print(s); n+=1
        if n>=N: break
PY
}

query(){
  rm -rf "$2"; mkdir -p "$2"
  cat "$1" | xargs -P "$PAR" -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r POS R A U D <<< "{}"
    CAPS_QUERY_MM='"$MM"' '"$BIN"'/mp_decode query win.capsule '"$2"'/up_$POS.fa   "$U" >/dev/null 2>'"$2"'/up_$POS.occ
    CAPS_QUERY_MM='"$MM"' '"$BIN"'/mp_decode query win.capsule '"$2"'/down_$POS.fa "$D" >/dev/null 2>'"$2"'/down_$POS.occ'
}

build het > bl_het.tsv
build neg > bl_neg.tsv
echo "  het probes: $(wc -l < bl_het.tsv)   controls: $(wc -l < bl_neg.tsv)"
query bl_het.tsv qbl_het
query bl_neg.tsv qbl_neg
python3 "$HOME/gc/scripts/score_bilateral.py"
