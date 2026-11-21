#!/usr/bin/env bash
# T3.4 / T3.5 — EXACT-MATCH RECALL AND HETEROZYGOUS LOCUS RETRIEVAL, one individual,
# one chr20 window, end to end from public GIAB data.
#
#   usage: run_locus_retrieval_capsule.sh <individual> <start> <end> <workdir> [N_sites]
#     individual : HG002 | HG003 | HG004 | HG005
#     start,end  : GRCh37 chr20 coordinates, e.g. 3000000 3600000
#     workdir    : created if missing; every intermediate file lands here
#     N_sites    : maximum het SNV sites to test (default 400)
#
# Output: <workdir>/locus_retrieval_<IND>_<start>_<end>.csv, one row holding both
# tables' values for this individual. benchmark_1_run.sh (PHASES=4) splits these rows
# into results/claim3/claim3_T3.4_exact_match_recall.csv (Table 5) and
# results/claim3/claim3_T3.5_locus_retrieval.csv (Table 6).
#
# Everything is real data: reference slice from UCSC (hg19), variants from the GIAB
# v4.2.1 benchmark VCF, reads range-fetched from the GIAB 300x Illumina BAM and
# subsampled to ~30x (samtools -s 42.10, a fixed seed). The reference slice is
# verified base-for-base against the VCF REF column, and the archive is verified
# byte-identical to its input, before any query result is allowed to count.
#
# TWO CONFIGURATIONS ARE MEASURED ON THE SAME ARCHIVE AND THE SAME SITES:
#   native      the archive, no index file            CAPS_QUERY_MM=$MM_NATIVE
#   completion  archive + FM-index, built from it     CAPS_QUERY_MM=$MM_XMI CAPS_QUERY_CONTAIN=1
# "The archive" means exactly that: the query is pointed at a sidecar path
# that does not exist, so it decodes the archive itself and emits what the
# archive alone yields. The FM-index step is the one post-processing command
# (`index` with CAPS_FMI=1): it writes the decoded sidecar (.qidx) and the
# FM-index over the reads the pseudogenome cannot represent (.qidx.fmi,
# include/caps_fmi.h). COMPLETION=xmi re-runs the older k-mer completion index
# (.qidx.xmi) in its place, for comparison.
# "Sites" is the number of het SNV sites that pass the admission rule below; it is
# the same set for both configurations and does not depend on the completion index.
#
# ADMISSION RULE (identical for every individual and window, nothing tuned): a het
# site is counted only if both 40-bp flanking reference probes are clean A/C/G/T and
# the upstream probe occurs in the pseudogenome. Sites are taken in VCF order.
#
# REQUIREMENTS (none are installed by this script):
#   built binaries best106 and capsule_decode in $CAPSULE_BIN_DIR
#       (cmake -S . -B build && cmake --build build; or the Dockerfile)
#   samtools, tabix, curl, python3, awk, xargs; internet access to NCBI and UCSC
# Paths come from scripts/benchmark_scripts/capsule_config.sh (env vars override).
set -euo pipefail

IND="${1:?usage: run_locus_retrieval_capsule.sh <HG002|HG003|HG004|HG005> <start> <end> <workdir> [N]}"
LO="${2:?}"; HI="${3:?}"; WD="${4:?}"; N="${5:-400}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../benchmark_scripts/capsule_config.sh
source "$HERE/../benchmark_scripts/capsule_config.sh"
BIN="$CAPSULE_BIN_DIR"
ENC="$BIN/best106"; DEC="$BIN/capsule_decode"
for t in "$ENC" "$DEC"; do
  [ -x "$t" ] || { echo "missing binary: $t (build first, see header)" >&2; exit 1; }
done
for t in samtools tabix curl python3 awk xargs; do
  command -v "$t" >/dev/null 2>&1 || { echo "missing tool on PATH: $t" >&2; exit 1; }
done

MM_NATIVE="${MM_NATIVE:-2}"; MM_XMI="${MM_XMI:-3}"; PAR="${PAR:-6}"
COMPLETION="${COMPLETION:-fmi}"
case "$COMPLETION" in fmi|xmi) ;; *) echo "COMPLETION must be fmi or xmi" >&2; exit 1;; esac
PL=40; GAP=5; PAD=1000
REF0=$((LO-PAD+1))

# ── public GIAB data (same BAM URLs as run_window_bench_capsule.sh) ──────────
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
case "$IND" in
  HG002) T=AshkenazimTrio/HG002_NA24385_son
         BAM="$G/data/$T/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam";;
  HG003) T=AshkenazimTrio/HG003_NA24149_father
         BAM="$G/data/$T/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam";;
  HG004) T=AshkenazimTrio/HG004_NA24143_mother
         BAM="$G/data/$T/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam";;
  HG005) T=ChineseTrio/HG005_NA24631_son
         BAM="$G/data/$T/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.hs37d5.300x.bam";;
  *) echo "unknown individual: $IND" >&2; exit 1;;
esac
VCF="$G/release/$T/NISTv4.2.1/GRCh37/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"

mkdir -p "$WD"; cd "$WD"
echo "=== $IND chr20:$LO-$HI -> $WD ==="

[ -s giab_win.vcf ] || tabix "$VCF" 20:$LO-$HI > giab_win.vcf
echo "  GIAB records: $(wc -l < giab_win.vcf)"

[ -s ref.txt ] || { curl -sS "https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr20;start=$((LO-PAD));end=$((HI+PAD))" -o seq.json
  python3 -c "
import json;open('ref.txt','w').write(json.load(open('seq.json'))['dna'].upper())"; }

# GATE: reference slice and coordinate frame must agree with GIAB.
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

# real reads for exactly this window: 300x -> ~30x
[ -s win.bam ] || samtools view -b -F 0x900 -s 42.10 "$BAM" 20:$((LO-PAD+1))-$((HI+PAD)) > win.bam
[ -s reads.fq ] || samtools fastq -n win.bam > reads.fq 2>/dev/null
NREADS=$(( $(wc -l < reads.fq) / 4 ))
echo "  reads: $NREADS"

# real GIAB het SNVs
awk -F'\t' '$7=="PASS" && length($4)==1 && length($5)==1 && $5~/^[ACGT]$/ {
  split($10,g,":"); if(g[1]=="0|1"||g[1]=="1|0"||g[1]=="0/1"||g[1]=="1/0") print $2"\t"$4"\t"$5 }' \
  giab_win.vcf > sites.tsv
echo "  het SNVs: $(wc -l < sites.tsv)"

# encode (ARCHIVE names the output; the three DUMP_* flags are required for the
# decoder's plain-reconstruction path, see the Dockerfile comment)
if [ ! -s win.capsule ]; then
  ARCHIVE=win.capsule CAPS_SPANS=1 DUMP_PERM=1 DUMP_MM=1 DUMP_LIT=1 \
    "$ENC" reads.fq 3 16 16 22 16 16 1 24 64 1 > enc.log 2>&1
  [ -s win.capsule ] || { [ -s out.arcs2 ] && cp out.arcs2 win.capsule; }
fi
[ -s win.capsule ] || { echo "encode produced no archive, see $WD/enc.log" >&2; exit 1; }
ARCBYTES=$(stat -c%s win.capsule)
echo "  archive: $ARCBYTES B"

# GATE: lossless, or the archive is not a valid product to query.
rm -rf dec && mkdir -p dec
"$DEC" win.capsule dec dec/reads.out >/dev/null 2>&1
awk 'NR%4==2' reads.fq > orig.txt
cmp -s orig.txt dec/reads.out && echo "  LOSSLESS: yes" || { echo "  LOSSLESS: NO" >&2; exit 1; }

[ -s pg.fa ] || "$DEC" export win.capsule pg.fa >/dev/null 2>&1
PGBP=$(grep -v '^>' pg.fa | tr -d '\n' | wc -c)
# sidecar (.qidx) plus the completion index, built together in one pass
if [ "$COMPLETION" = fmi ]; then
  [ -s win.capsule.qidx.fmi ] || CAPS_PILEUP=1 CAPS_FMI=1 "$DEC" index win.capsule win.capsule.qidx 2>&1 | grep -E '\[(index|fmi)\]' || true
else
  [ -s win.capsule.qidx.xmi ] || CAPS_PILEUP=1 CAPS_XMI=1 "$DEC" index win.capsule win.capsule.qidx 2>&1 | grep -E '\[(index|xmi)\]' || true
fi
[ -s "win.capsule.qidx.$COMPLETION" ] || { echo "completion index (.qidx.$COMPLETION) was not built" >&2; exit 1; }

build(){   # $1 = het | neg  -> probe table on stdout
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
    if pg.find(up)<0: return None          # admission rule
    return f"{pos}\t{r}\t{a}\t{up}\t{down}"
n=0
if MODE=="het":
    for line in open("sites.tsv"):
        pos,r,a=line.split()
        s=emit(int(pos),r,a)
        if s: print(s); n+=1
        if n>=N: break
else:   # homozygous negative controls: positions >200 bp from any GIAB variant
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

query(){   # $1 probes  $2 outdir  $3 mismatch tolerance  $4 0 = archive | 1 = archive + FM-index
  rm -rf "$2"; mkdir -p "$2"
  export QD="$2" QMM="$3" QCT="$4" QDEC="$DEC" QARC="$PWD/win.capsule" QSEL="$COMPLETION"
  export QNOIDX="$PWD/.no_index_for_archive_queries"
  rm -f "$QNOIDX"
  xargs -P "$PAR" -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r POS R A U D <<< "{}"
    if [ "$QCT" = 1 ]; then export CAPS_QUERY_CONTAIN=1 CAPS_QUERY_INDEX="$QSEL"
    else export CAPS_QIDX="$QNOIDX"; fi
    CAPS_QUERY_MM="$QMM" "$QDEC" query "$QARC" "$QD/up_$POS.fa"   "$U" >/dev/null 2>"$QD/up_$POS.occ"
    CAPS_QUERY_MM="$QMM" "$QDEC" query "$QARC" "$QD/down_$POS.fa" "$D" >/dev/null 2>"$QD/down_$POS.occ"' < "$1"
}

build het > bl_het.tsv
build neg > bl_neg.tsv
NSITES=$(wc -l < bl_het.tsv)
echo "  admitted het sites: $NSITES   negative controls: $(wc -l < bl_neg.tsv)"

query bl_het.tsv q_het_native "$MM_NATIVE" 0
query bl_neg.tsv q_neg_native "$MM_NATIVE" 0
query bl_het.tsv q_het_comp    "$MM_XMI"    1
query bl_neg.tsv q_neg_comp    "$MM_XMI"    1

# csv field helpers: score_locus_retrieval prints CSV,label,n,<up b,o,n>,<down b,o,n>,<bilateral b,o,n>
HETN=$(python3 "$HERE/score_locus_retrieval.py" "$IND-native" bl_het.tsv q_het_native bl_neg.tsv q_neg_native | tee score_native.txt | grep '^CSV,')
HETX=$(python3 "$HERE/score_locus_retrieval.py" "$IND-completion" bl_het.tsv q_het_comp bl_neg.tsv q_neg_comp | tee score_completion.txt | grep '^CSV,')
cat score_native.txt score_completion.txt
# bilateral "both" count is field 10; line 1 = het panel, line 2 = negative control
nat_het=$(echo "$HETN" | sed -n 1p | cut -d, -f10);  nat_fp=$(echo "$HETN" | sed -n 2p | cut -d, -f10)
xmi_het=$(echo "$HETX" | sed -n 1p | cut -d, -f10);  xmi_fp=$(echo "$HETX" | sed -n 2p | cut -d, -f10)

rec(){  # $1 queries dir -> exact-match recall of the het probes against brute-force truth
  DECODED=dec/reads.out QDIR="$1" PROBES=bl_het.tsv python3 "$HERE/score_exact_match_recall.py" "$NSITES" | tee "exact_$1.txt" | grep '^CSV,' | cut -d, -f4
}
rec_nat=$(rec q_het_native); rec_xmi=$(rec q_het_comp)

OUTCSV="locus_retrieval_${IND}_${LO}_${HI}.csv"
{
  echo "individual,region,reads,pseudogenome_bp,archive_bytes,sites,locus_archive_both,locus_archive_fm_index_both,negctl_fp_archive,negctl_fp_archive_fm_index,exact_recall_archive,exact_recall_archive_fm_index"
  echo "$IND,chr20:$LO-$HI,$NREADS,$PGBP,$ARCBYTES,$NSITES,$nat_het,$xmi_het,$nat_fp,$xmi_fp,$rec_nat,$rec_xmi"
} > "$OUTCSV"
echo "=== result: $WD/$OUTCSV ==="; cat "$OUTCSV"
