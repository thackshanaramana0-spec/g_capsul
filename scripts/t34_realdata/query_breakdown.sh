#!/usr/bin/env bash
# Separate the two things bundled into our 59 ms/query figure.
#
# Every `capsule_decode query` is a fresh process that loads the sidecar before
# it does any work, so a per-query number taken from 50 separate invocations
# charges the load 50 times. samtools pays a much smaller fixed cost because a
# BAM+BAI is designed to be seeked into.
#
# Running the COORDINATE arm and the SEQUENCE arm over the same 50 loci, in the
# same way, isolates the parts: the coordinate arm pays load + seek, the
# sequence arm pays load + seek + search. The difference is what content
# addressing actually costs, and the coordinate arm is the like-for-like
# comparison against samtools view.
set -u
cd ~/t34real
T=/usr/bin/time
secs(){ grep "Elapsed" "$1" | sed 's/.*: //'; }
head -50 probes_real.tsv > q50.tsv

echo "--- ours, COORDINATE arm (load + seek) ---"
rm -rf qc && mkdir -p qc
$T -v bash -c 'while IFS=$'"'"'\t'"'"' read -r POS R A P V; do
    ~/bin/capsule_decode query realreads.capsule qc/$POS.fa $((V-100))-$((V+100)) >/dev/null 2>&1
  done < q50.tsv' 2> qc.time
echo "    wall $(secs qc.time)   answers $(ls qc | wc -l)"

echo "--- ours, SEQUENCE arm exact (load + seek + exact search) ---"
rm -rf qs0 && mkdir -p qs0
$T -v bash -c 'while IFS=$'"'"'\t'"'"' read -r POS R A P V; do
    CAPS_QUERY_MM=0 ~/bin/capsule_decode query realreads.capsule qs0/$POS.fa "$P" >/dev/null 2>&1
  done < q50.tsv' 2> qs0.time
echo "    wall $(secs qs0.time)   answers $(ls qs0 | wc -l)"

echo "--- ours, SEQUENCE arm mm=2 (load + seek + tolerant search) ---"
rm -rf qs2 && mkdir -p qs2
$T -v bash -c 'while IFS=$'"'"'\t'"'"' read -r POS R A P V; do
    CAPS_QUERY_MM=2 ~/bin/capsule_decode query realreads.capsule qs2/$POS.fa "$P" >/dev/null 2>&1
  done < q50.tsv' 2> qs2.time
echo "    wall $(secs qs2.time)   answers $(ls qs2 | wc -l)"

echo "--- baseline, samtools view by coordinate ---"
$T -v bash -c 'while IFS=$'"'"'\t'"'"' read -r POS R A P V; do
    samtools view -c aln.bam 20_window:$((POS-2999001-100))-$((POS-2999001+100))
  done < q50.tsv > /dev/null 2>&1' 2> qb.time
echo "    wall $(secs qb.time)"

echo
echo "--- fixed cost per invocation (usage print, no archive touched) ---"
$T -v bash -c 'for i in $(seq 50); do ~/bin/capsule_decode >/dev/null 2>&1; done' 2> qz.time
echo "    wall $(secs qz.time)"
