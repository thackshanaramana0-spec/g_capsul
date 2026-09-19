#!/usr/bin/env bash
# Cost of becoming locus-addressable, both sides. bwa is used for the baseline
# because the published Claim 3 comparison already uses bwa.
#
# Every stage is checked for a non-empty result before its timing is believed --
# the first attempt at this silently reported a 0.23 s "alignment" that had
# produced nothing, because the binary was not on PATH.
set -u
cd ~/t34real
T=/usr/bin/time
ok(){ [ -s "$1" ] || { echo "FAILED: $1 is empty or missing"; exit 1; }; }
secs(){ grep "Elapsed" "$1" | sed 's/.*: //'; }
rss(){  grep "Maximum resident" "$1" | sed 's/.*: //'; }

echo "############ 1. build: from reads to locus-addressable ############"
echo
echo "--- G_CAPSUL: sidecar from the capsule, NO reference ---"
rm -f realreads.capsule.qidx
CAPS_PILEUP=1 $T -v ~/bin/capsule_decode index realreads.capsule realreads.capsule.qidx \
  >/dev/null 2> sc_index.time
ok realreads.capsule.qidx
printf '  wall %-10s peak RSS %s kB\n' "$(secs sc_index.time)" "$(rss sc_index.time)"
printf '  archive %s B   sidecar %s B\n' \
  "$(stat -c%s realreads.capsule)" "$(stat -c%s realreads.capsule.qidx)"

echo
echo "--- baseline: reference + bwa index + bwa mem + sort + index ---"
python3 -c "
s=open('ref.txt').read().strip()
w=open('ref.fa','w'); w.write('>20_window\n')
[w.write(s[i:i+60]+'\n') for i in range(0,len(s),60)]"
ok ref.fa
printf '  reference the baseline REQUIRES: %s B\n' "$(stat -c%s ref.fa)"

$T -v bwa index ref.fa > /dev/null 2> bwaidx.time
ok ref.fa.bwt
$T -v bash -c 'bwa mem -t 6 ref.fa real_reads.fq > aln.sam' 2> bwamem.time
ok aln.sam
$T -v samtools sort -@ 4 -o aln.bam aln.sam >/dev/null 2> sort.time
ok aln.bam
$T -v samtools index aln.bam 2> idx.time
ok aln.bam.bai
for s in bwaidx bwamem sort idx; do
  printf '  %-8s wall %-10s peak RSS %s kB\n' "$s" "$(secs $s.time)" "$(rss $s.time)"
done
printf '  bam %s B   bai %s B   bwa index %s B\n' \
  "$(stat -c%s aln.bam)" "$(stat -c%s aln.bam.bai)" \
  "$(cat ref.fa.amb ref.fa.ann ref.fa.bwt ref.fa.pac ref.fa.sa | wc -c)"

echo
echo "############ 2. per-query latency, 50 loci ############"
head -50 probes_real.tsv > q50.tsv
rm -rf qbench && mkdir -p qbench
$T -v bash -c 'while IFS=$'"'"'\t'"'"' read -r POS R A P V; do
    CAPS_QUERY_MM=2 ~/bin/capsule_decode query realreads.capsule qbench/$POS.fa "$P" >/dev/null 2>&1
  done < q50.tsv' 2> qcaps.time
n=$(ls qbench | wc -l)
printf '  G_CAPSUL by SEQUENCE   wall %-10s peak RSS %s kB   (%s answers)\n' \
  "$(secs qcaps.time)" "$(rss qcaps.time)" "$n"

$T -v bash -c 'while IFS=$'"'"'\t'"'"' read -r POS R A P V; do
    samtools view -c aln.bam 20_window:$((POS-2999001-100))-$((POS-2999001+100))
  done < q50.tsv > bamcounts.txt 2>/dev/null' 2> qbam.time
ok bamcounts.txt
printf '  samtools by COORDINATE wall %-10s peak RSS %s kB   (%s answers, %s reads)\n' \
  "$(secs qbam.time)" "$(rss qbam.time)" "$(wc -l < bamcounts.txt)" \
  "$(awk '{s+=$1} END{print s}' bamcounts.txt)"
