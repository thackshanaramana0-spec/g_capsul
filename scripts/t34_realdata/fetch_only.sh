#!/usr/bin/env bash
# Stage 1 only: fetch the held-out window's reads. Split out because the full
# pipeline has twice been killed by Windows-side memory pressure during the
# download, and a partial win.bam is worth keeping rather than restarting.
#
# Window is 200 kb rather than 600 kb: the NCBI server began throttling after
# repeated range requests, and a third of the data is a third of the transfer.
# ~160 het SNVs is still a meaningful independent replicate.
set -u
LO=4000000; HI=4200000; PAD=1000
WD=~/hold3
BAM=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz

mkdir -p "$WD"; cd "$WD"
[ -s giab_win.vcf ] || tabix "$GIAB" 20:$LO-$HI > giab_win.vcf
echo "GIAB records: $(wc -l < giab_win.vcf)"

[ -s ref.txt ] || { curl -sS "https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr20;start=$((LO-PAD));end=$((HI+PAD))" -o seq.json
  python3 -c "import json;open('ref.txt','w').write(json.load(open('seq.json'))['dna'].upper())"; }
echo "reference bp: $(wc -c < ref.txt)"

if [ ! -s win.bam ]; then
  samtools view -b -F 0x900 -s 42.10 "$BAM" 20:$((LO-PAD+1))-$((HI+PAD)) > win.bam
fi
echo "win.bam bytes: $(stat -c%s win.bam)"
samtools index win.bam 2>/dev/null || true
echo "reads: $(samtools view -c win.bam)"
