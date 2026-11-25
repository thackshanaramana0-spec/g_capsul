# T1.2: compression and decompression wall time

This file records the wall-time comparison of G_CAPSUL, SPRING and Genozip on the 19 datasets of
Claim 1 (see [`overview_claim1.md`](overview_claim1.md)). G_CAPSUL is not the fastest tool on any
dataset, for either compression or decompression. The explanation is the cost of the structure
that the archive retains, described below. The numbers are in
`results/claim1/claim1_T1.1_T1.2.csv` (columns `compress_s` and `decompress_s`).

## What was measured

Compression and decompression wall time for each tool on the same 19 datasets, in the same run on
the same machine, one timed job at a time. This is the protocol stated in the manuscript table
caption.

## Result

The per-dataset comparison of `compress_s` and `decompress_s` across the three tools, recomputed
from the CSV, gives:

| | G_CAPSUL fastest |
|---|---|
| Compression, out of 19 datasets | 0 |
| Decompression, out of 19 datasets | 0 |

SPRING and Genozip share the fastest position between them. The result is a plain loss on this
axis and should be described that way.

The ERR5181310 row, which agrees with the manuscript table, shows the pattern. Compression took
5.76 s for G_CAPSUL, 3.67 s for SPRING (fastest) and 4.54 s for Genozip. Decompression took 2.70 s
for G_CAPSUL, 6.44 s for SPRING and 0.51 s for Genozip (fastest).

## Why G_CAPSUL is slower

The archive keeps a placement for every read, and building it takes work that a compressor
without this structure does not do:

- assembling the pseudogenome by overlap chaining (section 3.3 of
  [`../technical_architecture.md`](../technical_architecture.md)),
- mapping the remaining reads and finding MEM self-matches (sections 3.4 and 3.7),
- trying several encodings of each stream and keeping the smallest (section 3.9),
- searching a small grid of parameter settings per dataset (section 3.10).

This cost buys the smaller archive of T1.1 and the placements that Claims 2 and 3 reuse. The
sequence-only comparison with PgRC2 shows the same trade-off, with G_CAPSUL about 1.7 times slower
and 2.5 times heavier in memory at worst (see
[`mechanism_insight_claim1.md`](mechanism_insight_claim1.md)). The full-FASTQ comparison on 19
datasets is T1.2.

## A measured saving behind the shipped time

The encoder evaluates a grid of `MAXMAP` and `MINOV` settings and keeps the smallest archive.
`MINOV` is used in the round-2 sweep and `MAXMAP` only later, at pigeonhole mapping. Candidates
that share a `MINOV` therefore computed a bit-identical round-2 assembly, up to four times over.
The production logs confirmed this, with `round2:` printed as two distinct results across the 8
candidates and not as 8. On a full HG002 run the repeated work cost about 2,360 s of CPU, and
because the candidates ran concurrently at one thread each, it sat on the wall-clock critical
path. The encoder now computes round 2 once per `MINOV` value and forks the candidates from it
(the "candidates share a prefix" design in `src/encoder.cpp`, described in
[`code_mapping_claim1.md`](code_mapping_claim1.md)).

## Scope of this table

T1.2 measures compression and decompression against SPRING and Genozip. It is a different
measurement from the speed tables of Claim 3, where export is 129 to 784 times faster than SPAdes
(T3.1) and coverage is 16 to 54 times faster than bwa with samtools and mosdepth (T3.2). Those
compare operations performed on the already-built archive with other tools, so a slow T1.2 and
fast T3.1 and T3.2 do not contradict each other.
