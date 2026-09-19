#!/usr/bin/env bash
# END-TO-END LOCUS RETRIEVAL EXPERIMENT FOR ANY chr20 WINDOW.
#
#   usage: run_window.sh <start> <end> <workdir> [N_sites]
#
# Parameterised deliberately. The bilateral-anchoring mechanism was DISCOVERED
# by inspecting the failures in one window, which is a form of selection, so the
# only honest test of it is a window that had no part in the discovery. This
# script makes a held-out replicate a one-line command rather than a bespoke
# effort, and it changes nothing between windows -- same probes, same tolerance,
# same scorer, no per-window constants.
#
# Everything real: reference from UCSC, variants and phasing from the GIAB
# benchmark VCF, reads range-fetched from the GIAB 300x Illumina BAM and
# subsampled. The reference slice is verified base-for-base against the VCF REF
# column before anything downstream runs.
set -eu
LO="${1:?usage: run_window.sh <start> <end> <workdir> [N]}"
HI="${2:?}"; WD="${3:?}"; N="${4:-400}"
MM="${MM:-2}"; PAR="${PAR:-6}"; PL=40; GAP=5
BIN="${BIN:-$HOME/bin}"
PAD=1000
REF0=$((LO-PAD+1))

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
BAM=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam

mkdir -p "$WD"; cd "$WD"
echo "=== window chr20:$LO-$HI -> $WD ==="

[ -s giab_win.vcf ] || tabix "$GIAB" 20:$LO-$HI > giab_win.vcf
echo "  GIAB records: $(wc -l < giab_win.vcf)"

[ -s ref.txt ] || { curl -sS "https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr20;start=$((LO-PAD));end=$((HI+PAD))" -o seq.json
  python3 -c "
import json;open('ref.txt','w').write(json.load(open('seq.json'))['dna'].upper())"; }

# GATE: the reference slice and coordinate frame must agree with GIAB, or every
# number downstream is silently wrong.
python3 - "$REF0" <<'PY'
import sys
REF0=int(sys.argv[1]); ref=open("ref.txt").read().strip()
ok=bad=0
for l in open("giab_win.vcf"):
    if l.startswith("#"): continue
    f=l.split("\t"); p=int(f[1]); r=f[3]
    ok += ref[p-REF0:p-REF0+len(r)]==r
    bad += ref[p-REF0:p-REF0+len(r)]!=r
print(f"  reference check: {ok} match, {bad} mismatch")
assert bad==0 and ok>0, "reference/coordinate frame does not agree with GIAB"
PY

# real reads for exactly this window, 300x -> ~30x
[ -s win.bam ] || samtools view -b -F 0x900 -s 42.10 "$BAM" 20:$((LO-PAD+1))-$((HI+PAD)) > win.bam
[ -s reads.fq ] || samtools fastq -n win.bam > reads.fq 2>/dev/null
echo "  reads: $(( $(wc -l < reads.fq) / 4 ))"

# real GIAB het SNVs
awk -F'\t' '$7=="PASS" && length($4)==1 && length($5)==1 && $5~/^[ACGT]$/ {
  split($10,g,":"); if(g[1]=="0|1"||g[1]=="1|0"||g[1]=="0/1"||g[1]=="1/0") print $2"\t"$4"\t"$5 }' \
  giab_win.vcf > sites.tsv
echo "  het SNVs: $(wc -l < sites.tsv)"

export CAPS_SPANS=1 DUMP_PERM=1 DUMP_MM=1 DUMP_LIT=1
if [ ! -s win.capsule ]; then
  "$BIN/best106" reads.fq 3 16 16 22 16 16 1 24 64 1 > enc.log 2>&1
  cp out.arcs2 win.capsule
fi
echo "  archive: $(stat -c%s win.capsule) B"

# GATE: lossless, or the archive is not a valid product to query.
rm -rf dec && mkdir -p dec
"$BIN/mp_decode" win.capsule dec dec/reads.out >/dev/null 2>&1
awk 'NR%4==2' reads.fq > orig.txt
cmp -s orig.txt dec/reads.out && echo "  LOSSLESS: yes" || { echo "  LOSSLESS: NO"; exit 1; }

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
    if pg.find(up)<0: return None          # same admission rule as the published run
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

query(){   # $1 probes  $2 outdir -- ONE invocation carries both probes
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
