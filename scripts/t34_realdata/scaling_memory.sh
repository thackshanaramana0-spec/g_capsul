#!/usr/bin/env bash
# HOW DOES THE INDEX-BUILD MEMORY MARGIN SCALE?
#
# A single-point ratio ("6.7x lighter than bwa") is the wrong claim shape. bwa's
# peak RSS is dominated by the FM-index of the REFERENCE, which does not grow
# with the number of reads. Ours is dominated by the pseudogenome and the
# placement arrays, which DO. So the ratio must shrink with input size, and
# quoting one point either flatters us or will be shown to.
#
# This measures both sides at four input sizes from the same real reads, so the
# trend is characterised rather than a single number defended. Runs strictly
# sequentially -- never two timed jobs at once (standing rule 4).
set -u
cd ~/t34real
BIN="${BIN:-$HOME/bin}"
T=/usr/bin/time
rss(){ grep "Maximum resident" "$1" | sed 's/.*: //'; }
secs(){ grep "Elapsed" "$1" | sed 's/.*: //'; }

# the reference bwa needs and we do not
if [ ! -s ref.fa ]; then
  python3 -c "
s=open('ref.txt').read().strip()
w=open('ref.fa','w'); w.write('>20_window\n')
[w.write(s[i:i+60]+'\n') for i in range(0,len(s),60)]"
fi
[ -s ref.fa.bwt ] || bwa index ref.fa >/dev/null 2>&1

printf '%8s %10s %14s %14s %10s\n' frac reads ours_peak_kB bwa_peak_kB ratio
for FR in 25 50 75 100; do
  if [ "$FR" = "100" ]; then cp -f real_reads.fq s.fq
  else samtools view -b -s "3.$FR" win30x.bam 2>/dev/null | samtools fastq -n - > s.fq 2>/dev/null
  fi
  N=$(( $(wc -l < s.fq) / 4 ))

  export CAPS_SPANS=1 DUMP_PERM=1 DUMP_MM=1 DUMP_LIT=1
  rm -f out.arcs2 s.capsule
  "$BIN/best106" s.fq 3 16 16 22 16 16 1 24 64 1 >/dev/null 2>&1
  cp out.arcs2 s.capsule
  CAPS_PILEUP=1 $T -v "$BIN/mp_decode" index s.capsule s.qidx >/dev/null 2>ours.time
  OURS=$(rss ours.time)

  $T -v bash -c "bwa mem -t 6 ref.fa s.fq > /dev/null" 2>bwa.time
  BWA=$(rss bwa.time)

  printf '%8s %10s %14s %14s %10s\n' "$FR%" "$N" "$OURS" "$BWA" \
    "$(echo "$BWA $OURS" | awk '{printf "%.1fx", $1/$2}')"
done
rm -f s.fq s.capsule s.qidx
echo
echo "READ THE TREND, NOT THE LAST ROW. If the ratio falls as reads grow, the"
echo "margin is a property of this input size and must not be quoted at WGS scale."
