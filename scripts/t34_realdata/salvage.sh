#!/usr/bin/env bash
# Can the partial download be used?
#
# `samtools view -b <region>` writes reads in coordinate order, so a truncated
# output is not corrupt data -- it is COMPLETE coverage of a shorter window,
# ending wherever the transfer stopped. That is still a legitimate held-out
# window, and it needs no further download.
#
# This finds the usable extent, trimming a margin off the end so the last reads
# are not partially covered.
set -u
cd ~/hold3
echo "win.bam bytes: $(stat -c%s win.bam)"

# tolerate a truncated final BGZF block
samtools view win.bam 2>/dev/null | awk '{print $4}' > pos.txt || true
n=$(wc -l < pos.txt)
echo "reads readable: $n"
if [ "$n" -lt 1000 ]; then echo "too few reads to use"; exit 1; fi
lo=$(head -1 pos.txt); hi=$(tail -1 pos.txt)
echo "coordinate span: $lo .. $hi"
echo "usable window (trimming 5 kb off the end): $lo .. $((hi-5000))"
echo "span bp: $(( hi-5000-lo ))"
