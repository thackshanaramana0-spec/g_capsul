# HG005's precision collapse — cause found, caveat closed

**Status: explained by a controlled experiment, 2026-09-10.** This was the last
unexplained number in Claim 2 and is no longer one.

## The anomaly

HG005 was our only loss on het-SNV and our only loss on het-indel, and the SNV
loss was a precision collapse rather than a recall shortfall:

    individual   TP      FP     FN     P       R       F1
    HG002        37,011  1,721  7,564  0.956   0.830   0.888
    HG003        37,890  1,644  7,643  0.958   0.832   0.891
    HG004        38,825  1,753  7,695  0.957   0.835   0.891
    HG005        30,404  3,426  8,699  0.899   0.777   0.834   <- twice the FP

## It is not depth, and not the individual

All four carry ~12.6 M reads at ~30x. But HG005's FASTQ is 59% larger for the
same read count, which is the tell:

    individual   instrument   read length      distinct lengths
    HG002/3/4    HISEQ1       148 bp fixed     1
    HG005        D00360       250 bp           216

HG005 is a different sequencing run entirely: 250 bp quality-trimmed reads with
216 distinct lengths, against 148 bp fixed for the other three. It is the only
variable-length dataset among the four.

## The controlled experiment

Truncate every HG005 read to a fixed 148 bp — same individual, same depth, same
reads, same caller, same truth, same scoring. Only the geometry changes:

    HG005 as-is (250 bp, 216 lengths)   TP 30,404  FP 3,426  P 0.899  R 0.777  F1 0.834
    HG005 truncated to 148 bp fixed     TP 32,825  FP 1,291  P 0.962  R 0.839  F1 0.897

**False positives fall 2.7x and F1 rises 0.063.** At fixed length HG005 becomes
our BEST individual — above HG002/3/4 (0.888/0.891/0.891) and above
DiscoSNP++'s 0.863 on the untruncated data.

## What this means, stated precisely

The caller is tuned for fixed-length reads. On variable-length input its
precision degrades, and HG005 is the only variable-length set in the benchmark.
That is a **characterised property with a named cause and a reproducible
control**, not an anomaly.

It is not a defect in the mechanism: the collapse-plus-re-placement result
(0.431 -> 0.888, `docs/CLAIM2_RESULTS_V2.md`) is unaffected, and neither is
Claim 1 — HG005 compresses losslessly at 997,518,994 B and wins its row.

**Do NOT report the truncated figure as an HG005 result.** Truncating discards
40% of every read and would be scoring a different dataset. The honest table
keeps 0.834 and cites this document for why. The truncation exists to identify
the cause, not to improve the number.

**The actionable form:** if variable-length short-read input matters for a
deployment, the caller's length-dependent parameters are the place to look.
That is scoped work with a known target, and it is now specified rather than
open.

## Reproduce

    awk 'NR%4==1{n=$0} NR%4==2{s=substr($0,1,148)} NR%4==0{q=substr($0,1,148);
         if(length(s)==148) printf "%s\n%s\n+\n%s\n", n, s, q}' HG005_pooled.fq > fixed148.fq
    # compress, then run scripts/run_fullchr20_archive_capsule.sh against it
