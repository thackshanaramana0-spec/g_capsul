# T1.1: archive size on 19 datasets

This file records the T1.1 result of Claim 1 (see [`overview_claim1.md`](overview_claim1.md)). It
gives the aggregate figures recomputed from the raw CSV and explains what the LOSSLESS column of
that file does and does not establish.

## What it measures

The size of the lossless whole-FASTQ archive (sequence, names, quality, line 3 and read order,
which together reconstruct the original file) for G_CAPSUL, SPRING and Genozip on 19 real
datasets from human, bacteria, archaea, protista, fungi, virus and plants.

## Result

The aggregate was recomputed by summing `archive_bytes` over the 19 datasets for each tool in
`results/claim1/claim1_T1.1_T1.2.csv`:

| | Against SPRING | Against Genozip |
|---|---|---|
| Datasets where G_CAPSUL is smaller | 19 of 19 | 19 of 19 |
| Aggregate margin | -6.03% | -43.26% |
| Archives verified lossless | 57 of 57 (19 datasets, 3 tools) | |

These match the headline figures of the development notes. For ERR5181310, the raw file of
471,183,048 bytes gives 8,679,796 bytes for G_CAPSUL (8.68 MB), 9,984,000 for SPRING (9.98 MB) and
9,222,378 for Genozip (9.22 MB). This agrees with the manuscript row `8.68 | 9.98 | 9.22`. All 19
rows are in the CSV and in the manuscript table.

The two margins should be cited separately and not averaged. SPRING is the closer comparator
(-6.03%), and the gap to Genozip (-43.26%) is much larger.

## What the LOSSLESS column establishes

The manuscript states that every archive was independently decoded and checked against its input
by comparing the sorted read names, sequences and quality strings. All 57 archives, including the
SPRING and Genozip outputs, were checked this way and none is assumed correct on the tool's own
word.

The comparison is on sorted records, so it does not by itself test that the original read order
and line 3 are reproduced. The exact byte-for-byte round trip, which does, is
`scripts/build_scripts/verify_roundtrip.sh`. It is described in the
[README](../../README.md) and was run on simulated inputs (including a shuffled-order file) and on
SRR29296997. The benchmark makes only the sorted-record claim for the other datasets.

The discipline behind the column has a history. The project once shipped an incomplete verification
path and had to correct four silent data-loss bugs that were found only by decoding archives and
not by trusting dump-level checks (see [`code_mapping_claim1.md`](code_mapping_claim1.md)).
