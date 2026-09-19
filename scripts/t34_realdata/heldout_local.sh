#!/usr/bin/env bash
# TWO HELD-OUT TESTS THAT NEED NO NEW DOWNLOAD.
#
# Bilateral anchoring was discovered by inspecting the nine failures among the
# first 400 het sites of chr20:3.0-3.6 Mb, on one archive. That is selection, so
# the result has to be re-tested on data that had no part in the discovery.
#
#   A. DIFFERENT ASSEMBLY, same window. The real reads are re-subsampled from
#      the same 300x BAM with a different seed and depth, then re-encoded. This
#      produces a DIFFERENT pseudogenome -- different chaining, different
#      contigs, different neighbourhoods. Since the failure mode is an artefact
#      of chaining, a different assembly is the sharpest available retest: if
#      bilateral only worked because of the particular contigs in archive one,
#      it will not work here.
#
#   B. HELD-OUT SITES, same archive. The published run used the first 400 het
#      SNVs in the window. Any remaining sites were never looked at.
#
# Nothing is tuned between runs. Same probes, same tolerance, same scorer.
set -eu
SRC=~/t34real
BIN="${BIN:-$HOME/bin}"
MM="${MM:-2}"; PAR="${PAR:-6}"; PL=40; GAP=5; REF0=2999001
SEED="${SEED:-99}"; FRAC="${FRAC:-07}"      # samtools -s SEED.FRAC

echo "############ A. DIFFERENT ASSEMBLY (seed $SEED, fraction .$FRAC) ############"
A=~/hold_assembly
rm -rf "$A"; mkdir -p "$A"
cp "$SRC"/ref.txt "$SRC"/giab_win.vcf "$SRC"/sites.tsv "$A"/
cd "$A"
samtools view -b -s "$SEED.$FRAC" "$SRC/win30x.bam" > sub.bam
samtools fastq -n sub.bam > reads.fq 2>/dev/null
echo "  reads: $(( $(wc -l < reads.fq) / 4 ))"

export CAPS_SPANS=1 DUMP_PERM=1 DUMP_MM=1 DUMP_LIT=1
"$BIN/best106" reads.fq 3 16 16 22 16 16 1 24 64 1 > enc.log 2>&1
cp out.arcs2 win.capsule
rm -rf dec && mkdir -p dec
"$BIN/mp_decode" win.capsule dec dec/reads.out >/dev/null 2>&1
awk 'NR%4==2' reads.fq > orig.txt
cmp -s orig.txt dec/reads.out && echo "  LOSSLESS: yes" || { echo "  LOSSLESS: NO"; exit 1; }
"$BIN/mp_decode" export win.capsule pg.fa >/dev/null 2>&1
CAPS_PILEUP=1 "$BIN/mp_decode" index win.capsule win.capsule.qidx 2>&1 | grep '\[index\]'
echo "  pseudogenome: $(grep -v '>' pg.fa | tr -d '\n' | wc -c) bp   (archive one: 2153559 bp)"

build(){   # $1 het|neg   $2 skip   $3 take
python3 - "$REF0" "$1" "$PL" "$GAP" "$2" "$3" <<'PY'
import sys, random
REF0=int(sys.argv[1]); MODE=sys.argv[2]; PL=int(sys.argv[3]); GAP=int(sys.argv[4])
SKIP=int(sys.argv[5]); TAKE=int(sys.argv[6])
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg.fa") if l[0]!='>')
def emit(pos,r,a):
    i=pos-REF0
    up=ref[i-GAP-PL:i-GAP]; down=ref[i+1+GAP:i+1+GAP+PL]
    if len(up)!=PL or len(down)!=PL: return None
    if set(up)-set("ACGT") or set(down)-set("ACGT"): return None
    if pg.find(up)<0: return None
    return f"{pos}\t{r}\t{a}\t{up}\t{down}"
rows=[]
if MODE=="het":
    for line in open("sites.tsv"):
        pos,r,a=line.split()
        s=emit(int(pos),r,a)
        if s: rows.append(s)
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
        if s: rows.append(s)
        if len(rows)>=SKIP+TAKE: break
for s in rows[SKIP:SKIP+TAKE]: print(s)
PY
}

query(){
  rm -rf "$2"; mkdir -p "$2"
  cat "$1" | xargs -P "$PAR" -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r POS R A U D <<< "{}"
    CAPS_QUERY_MM='"$MM"' '"$BIN"'/mp_decode query win.capsule '"$2"'/up_$POS.fa   "$U" >/dev/null 2>'"$2"'/up_$POS.occ
    CAPS_QUERY_MM='"$MM"' '"$BIN"'/mp_decode query win.capsule '"$2"'/down_$POS.fa "$D" >/dev/null 2>'"$2"'/down_$POS.occ'
}

build het 0 400 > bl_het.tsv
build neg 0 400 > bl_neg.tsv
echo "  het probes: $(wc -l < bl_het.tsv)  controls: $(wc -l < bl_neg.tsv)"
query bl_het.tsv qbl_het
query bl_neg.tsv qbl_neg
python3 "$HOME/gc/scripts/score_bilateral.py"

echo "############ B. HELD-OUT SITES 401+ on the ORIGINAL archive ############"
cd "$SRC"
cp "$SRC"/pg_real.fa /tmp/pgb.fa
B=~/hold_sites
rm -rf "$B"; mkdir -p "$B"
cp ref.txt giab_win.vcf sites.tsv "$B"/
cp pg_real.fa "$B"/pg.fa
cp realreads.capsule "$B"/win.capsule
cp realreads.capsule.qidx "$B"/win.capsule.qidx
cd "$B"
build het 400 200 > bl_het.tsv
build neg 400 200 > bl_neg.tsv
echo "  held-out het probes: $(wc -l < bl_het.tsv)  controls: $(wc -l < bl_neg.tsv)"
if [ "$(wc -l < bl_het.tsv)" -gt 0 ]; then
  query bl_het.tsv qbl_het
  query bl_neg.tsv qbl_neg
  python3 "$HOME/gc/scripts/score_bilateral.py"
else
  echo "  no held-out sites remain in this window"
fi
