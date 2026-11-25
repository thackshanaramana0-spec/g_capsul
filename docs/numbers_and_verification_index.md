# Numbers and verification index

This file lists the tables behind the manuscript, the result file each number comes from, and how
each was checked. Cite a number from here together with its result file.

## How numbers were checked

| Status | Meaning |
|---|---|
| Re-derived | Recomputed independently from raw data. |
| Cross-checked | Compared with its result file, row by row. |
| Cited | Taken from a result file produced by the benchmark scripts and not recomputed separately. |

The result files in [`../results/`](../results/) are authoritative. Where this file and a result file
disagree, the result file is correct. Scripts that produce each file are listed in
[`../README.md`](../README.md#reproducing-the-results).

---

# Compression

## Archive size, all 19 datasets (MB)

Status: re-derived for the aggregate and cross-checked row by row. Source: `results/claim1/claim1_T1.1_T1.2.csv`.

| Dataset | Raw MB | CAPSULE MB | SPRING MB | Genozip MB |
|---|---|---|---|---|
| DRR976266 | 2248.04 | 56.54 | 61.76 | 201.14 |
| ERR12954017 | 333.14 | 15.16 | 17.00 | 39.43 |
| ERR17740259 | 1336.96 | 83.20 | 92.40 | 161.81 |
| ERR5181310 | 471.18 | 8.68 | 9.98 | 9.22 |
| ERR552797 | 431.08 | 46.97 | 52.13 | 82.80 |
| HG002 | 4279.20 | 573.77 | 598.19 | 951.00 |
| HG003 | 4279.57 | 567.12 | 591.14 | 944.22 |
| HG004 | 4279.33 | 606.19 | 634.75 | 1041.66 |
| HG005 | 6809.02 | 997.52 | 1081.38 | 1694.51 |
| SRR065390 | 11282.99 | 887.41 | 942.64 | 1597.08 |
| SRR10676752 | 15372.50 | 1650.03 | 1718.55 | 2954.90 |
| SRR2584863 | 695.16 | 68.43 | 74.09 | 114.59 |
| SRR29296997 | 210.23 | 15.66 | 17.54 | 27.23 |
| SRR32429602 | 1868.92 | 55.12 | 57.55 | 107.18 |
| SRR36741279 | 1661.16 | 106.25 | 117.62 | 194.66 |
| SRR37283774 | 1019.82 | 67.20 | 71.34 | 94.99 |
| SRR39257532 | 1761.85 | 108.27 | 153.86 | 227.02 |
| SRR40271341 | 290.89 | 40.55 | 45.78 | 62.26 |
| SRR554369 | 456.44 | 57.29 | 59.20 | 88.81 |

G_CAPSUL is smaller than SPRING on all 19 datasets (aggregate -6.03%) and smaller than Genozip on
all 19 (aggregate -43.26%). All 57 archives (19 datasets, 3 tools) are marked lossless in the CSV.
This is a sorted-record comparison, see [`limitations_and_scope.md`](limitations_and_scope.md).

## Wall time (s), compression and decompression

Status: re-derived. Source: same CSV.

| Dataset | CAPS compress | SPRING compress | Genozip compress | CAPS decompress | SPRING decompress | Genozip decompress | Fastest compress | Fastest decompress |
|---|---|---|---|---|---|---|---|---|
| DRR976266 | 38.30 | 11.66 | 4.93 | 10.73 | 9.71 | 1.98 | Genozip | Genozip |
| ERR12954017 | 6.28 | 5.55 | 1.50 | 1.91 | 6.72 | 0.57 | Genozip | Genozip |
| ERR17740259 | 23.55 | 10.84 | 3.43 | 8.40 | 5.72 | 1.46 | Genozip | Genozip |
| ERR5181310 | 5.76 | 3.67 | 4.54 | 2.70 | 6.44 | 0.51 | SPRING | Genozip |
| ERR552797 | 8.44 | 8.99 | 3.14 | 4.05 | 8.94 | 1.34 | Genozip | Genozip |
| HG002 | 210.49 | 54.00 | 25.59 | 47.24 | 26.05 | 10.75 | Genozip | Genozip |
| HG003 | 207.33 | 53.05 | 25.61 | 46.87 | 26.48 | 10.79 | Genozip | Genozip |
| HG004 | 218.24 | 57.99 | 8.99 | 50.02 | 26.33 | 3.50 | Genozip | Genozip |
| HG005 | 579.51 | 133.24 | 21.09 | 80.99 | 48.05 | 9.99 | Genozip | Genozip |
| SRR065390 | 367.25 | 112.94 | 25.93 | 82.72 | 41.24 | 12.11 | Genozip | Genozip |
| SRR10676752 | 945.38 | 388.19 | 26.31 | 223.56 | 88.63 | 15.10 | Genozip | Genozip |
| SRR2584863 | 12.38 | 9.79 | 3.93 | 5.81 | 5.30 | 1.74 | Genozip | Genozip |
| SRR29296997 | 4.53 | 4.72 | 1.56 | 1.49 | 4.19 | 0.73 | Genozip | Genozip |
| SRR32429602 | 43.03 | 11.56 | 4.24 | 11.64 | 9.62 | 1.46 | Genozip | Genozip |
| SRR36741279 | 43.10 | 15.02 | 6.13 | 12.90 | 6.87 | 2.58 | Genozip | Genozip |
| SRR37283774 | 26.65 | 8.75 | 4.18 | 7.46 | 3.66 | 1.64 | Genozip | Genozip |
| SRR39257532 | 56.94 | 34.49 | 4.59 | 13.72 | 10.55 | 1.97 | Genozip | Genozip |
| SRR40271341 | 6.44 | 10.87 | 1.35 | 3.38 | 9.26 | 0.50 | Genozip | Genozip |
| SRR554369 | 11.68 | 5.57 | 3.47 | 4.66 | 3.41 | 1.65 | Genozip | Genozip |

G_CAPSUL is the fastest tool in none of these 38 cells (19 compression and 19 decompression).
Genozip is fastest in 37 and SPRING in 1 (ERR5181310 compression). The manuscript's bold marking
agrees cell for cell.

---

# Variant calling

## Heterozygous SNVs, four individuals

Status: cross-checked against the manuscript. Source: `results/claim2/claim2_T2.1_snv.csv`.

| Individual | Tool | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|
| HG002 | CAPSULE | 37,011 | 1,721 | 7,564 | 0.956 | 0.830 | **0.888** |
| HG002 | DiscoSNP++ | 34,036 | 1,735 | 10,560 | 0.952 | 0.763 | 0.847 |
| HG002 | Kmer2SNP | 13,565 | 289 | 31,010 | 0.979 | 0.304 | 0.464 |
| HG003 | CAPSULE | 37,890 | 1,644 | 7,643 | 0.958 | 0.832 | **0.891** |
| HG003 | DiscoSNP++ | 34,811 | 1,633 | 10,744 | 0.955 | 0.764 | 0.849 |
| HG003 | Kmer2SNP | 13,682 | 303 | 31,851 | 0.978 | 0.300 | 0.460 |
| HG004 | CAPSULE | 38,825 | 1,753 | 7,695 | 0.957 | 0.835 | **0.891** |
| HG004 | DiscoSNP++ | 35,706 | 1,617 | 10,837 | 0.957 | 0.767 | 0.851 |
| HG004 | Kmer2SNP | 14,335 | 258 | 32,185 | 0.982 | 0.308 | 0.469 |
| HG005 | CAPSULE | 30,404 | 3,426 | 8,699 | 0.899 | 0.777 | 0.834 |
| HG005 | DiscoSNP++ | 30,726 | 1,331 | 8,396 | 0.959 | 0.785 | **0.863** |
| HG005 | Kmer2SNP | 13,294 | 220 | 25,809 | 0.984 | 0.340 | 0.505 |

Mean F1 is 0.876 for G_CAPSUL, 0.853 for DiscoSNP++ and 0.475 for Kmer2SNP. G_CAPSUL is ahead on
three individuals and DiscoSNP++ on HG005.

## Heterozygous SNV F1 against depth, HG002

Status: cross-checked. Source: `results/claim2/claim2_T2.2_coverage_sweep.csv`.

| Depth | Reads | Archive bytes | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|---|
| 10x | 4,201,737 | 203,851,561 | 15,308 | 2,933 | 29,267 | 0.839 | 0.343 | 0.487 |
| 15x | 6,301,720 | 296,083,050 | 25,322 | 2,059 | 19,253 | 0.925 | 0.568 | 0.704 |
| 20x | 8,402,407 | 388,576,775 | 30,792 | 1,903 | 13,783 | 0.942 | 0.691 | 0.797 |
| 30x | 12,604,917 | 573,767,964 | 37,011 | 1,721 | 7,564 | 0.956 | 0.830 | 0.888 |

The 30x row is identical to the HG002 row above and was carried over, not re-run.

## Heterozygous indels, four individuals

Status: cross-checked. Source: `results/claim2/claim2_T2.3_indel.csv`.

| Individual | Tool | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|
| HG002 | CAPSULE | 3,871 | 824 | 3,913 | 0.825 | 0.497 | **0.620** |
| HG002 | DiscoSNP++ | 3,254 | 273 | 4,527 | 0.923 | 0.418 | 0.576 |
| HG003 | CAPSULE | 3,916 | 732 | 3,727 | 0.843 | 0.512 | **0.637** |
| HG003 | DiscoSNP++ | 3,293 | 289 | 4,346 | 0.919 | 0.431 | 0.587 |
| HG004 | CAPSULE | 4,029 | 759 | 3,926 | 0.842 | 0.506 | **0.632** |
| HG004 | DiscoSNP++ | 3,484 | 276 | 4,465 | 0.927 | 0.438 | 0.595 |
| HG005 | CAPSULE | 2,483 | 582 | 2,827 | 0.810 | 0.467 | 0.593 |
| HG005 | DiscoSNP++ | 2,400 | 223 | 2,907 | 0.915 | 0.452 | **0.605** |

Kmer2SNP is marked `NOT_APPLICABLE_SNP_ONLY` for all four individuals because it has no indel
model. Mean F1 is 0.621 for G_CAPSUL and 0.591 for DiscoSNP++. DiscoSNP++ is ahead on HG005.

## Multi-allelic SNV sites, complete census of chromosome 20

Status: cross-checked. Source: `results/claim2/claim2_T2.4_multiallelic.csv`.

| Individual | Scope | Sites | CAPSULE both alleles | CAPSULE rate | DiscoSNP++ multi-allelic records emitted | DiscoSNP++ total records |
|---|---|---|---|---|---|---|
| HG002 | chr20 COMPLETE CENSUS | 26 | **17** | 65.4% | 0 | 3,989 |

DiscoSNP++ emitted no multi-allelic record anywhere in its output for this individual, not only at
these 26 sites.

## Tetraploid calling, HG003 and HG004 reads, chr20:3,000,000-3,400,000

Status: cross-checked. Source: `results/claim2/claim2_T2.5_tetraploid.csv`.

| Class | Tool | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|
| SNV | CAPSULE | 472 | 7 | 101 | 0.985 | 0.824 | **0.897** |
| SNV | DiscoSNP++ | 373 | 8 | 200 | 0.979 | 0.651 | 0.782 |
| INDEL | CAPSULE | 51 | 8 | 61 | 0.864 | 0.455 | **0.597** |
| INDEL | DiscoSNP++ | 44 | 3 | 68 | 0.936 | 0.393 | 0.553 |

G_CAPSUL is ahead in both classes.

---

# Addressability

## Export speed, 6 datasets

Status: cited. Source: `results/claim3/claim3_T3.1_T3.2_T3.3.csv`.

| Dataset | CAPSULE (s) | SPAdes (s) | Speedup |
|---|---|---|---|
| ERR5181310 | 0.40 | 51.73 | 129.3x |
| SRR2584863 | 0.37 | 205.04 | 554.2x |
| SRR29296997 | 0.13 | 83.95 | 645.8x |
| SRR37283774 | 1.49 | 380.88 | 255.6x |
| DRR976266 | 0.70 | 548.82 | 784.0x |
| HG002 | 5.67 | 2,678.83 | 472.5x |

## Export correctness, SRR2584863 (E. coli)

Status: re-derived. The one dataset on which both tools can be compared fairly.

| Metric | CAPSULE (final, shipping config) | SPAdes | Winner |
|---|---|---|---|
| Genome fraction | 98.698% | 98.346% | CAPSULE |
| Indels /100kbp | 0.25 | 0.51 | CAPSULE |
| Mismatches /100kbp | 19.40 | 2.75 | SPAdes |
| Duplication ratio | 2.021 | 1.000 | SPAdes |
| N50 | 2,054 | 139,621 | SPAdes |
| Unaligned length | 1,895,813 | 9,582 | SPAdes |
| Misassemblies | 7 | 1 | SPAdes |

## Coverage speed, 6 datasets

Status: cited. Source: `results/claim3/claim3_T3.1_T3.2_T3.3.csv`.

| Dataset | CAPSULE (s) | bwa+samtools+mosdepth (s) | Speedup |
|---|---|---|---|
| ERR5181310 | 0.15 | 4.15 | 27.7x |
| SRR2584863 | 0.33 | 7.57 | 22.9x |
| SRR29296997 | 0.11 | 2.53 | 23.0x |
| SRR37283774 | 0.56 | 15.98 | 28.5x |
| DRR976266 | 1.49 | 23.98 | 16.1x |
| HG002 | 3.04 | 164.90 | 54.2x |

## Range query timing, all 19 datasets

Status: cited. Source: same CSV. There is no comparator.

| Dataset | Time (s) | Rows |
|---|---|---|
| SRR29296997 | 0.15 | 18,002 |
| ERR12954017 | 0.20 | 15,744 |
| SRR40271341 | 0.20 | 15,082 |
| ERR552797 | 0.34 | 10,334 |
| SRR2584863 | 0.46 | 23,604 |
| ERR5181310 | 0.50 | 30,510 |
| ERR17740259 | 0.57 | 36,328 |
| SRR554369 | 0.65 | 25,470 |
| DRR976266 | 1.57 | 24,348 |
| SRR37283774 | 1.95 | 18,486 |
| SRR36741279 | 3.14 | 16,816 |
| SRR39257532 | 3.77 | 7,184 |
| SRR32429602 | 3.79 | 40,964 |
| HG003 | 6.17 | 16,598 |
| HG002 | 6.27 | 16,506 |
| HG004 | 6.77 | 15,844 |
| HG005 | 10.15 | 11,356 |
| SRR065390 | 12.11 | 30,458 |
| SRR10676752 | 118.48 | 4,990 |

## Exact-match recall, four individuals

Status: re-derived. Sources: `docs/claim3/t34_exact_match_claim3.md` and
`results/claim3/claim3_T3.4_exact_match_recall.csv`.

| Individual | Reads with a probe | Recall, archive | Recall, archive + FM-index |
|---|---:|---:|---:|
| HG002 | 12,059 | 0.8536 | **1.0000** |
| HG003 | 9,454 | 0.8457 | **1.0000** |
| HG004 | 12,455 | 0.8099 | **1.0000** |
| HG005 | 10,121 | 0.6047 | **1.0000** |

Measured with `scripts/claim3_scripts/run_locus_retrieval_capsule.sh`, using the same probes and
scorer in every column. A query that also has the decoded read deviations gives 0.9999, 0.9998,
1.0000 and 0.9990. The FM-index gives 1.0000 at mismatch tolerance 0 and 3. Earlier archive-only
values (0.9619, 0.9612, 0.8563) could not be reproduced from the repository and are replaced.

## Locus retrieval, four individuals, one locus

Status: re-derived. Sources: `docs/claim3/t35_locus_retrieval_claim3.md` and
`results/claim3/claim3_T3.5_locus_retrieval.csv`.

| Individual | Locus | Sites | Position query (tolerance 2) | Archive + FM-index | Negative-control hits, position query / FM-index |
|---|---|---|---|---|---|
| HG002 | chr20:3.0-3.6Mb | 400 | 400/400 | 400/400 | 29 / 48 of 400 |
| HG003 | chr20:3.0-3.6Mb | 335 | 335/335 | 335/335 | 35 / 52 of 400 |
| HG004 | chr20:3.0-3.6Mb | 400 | 399/400 | **400/400** (miss fixed) | 29 / 49 of 400 |
| HG005 | chr20:3.0-3.6Mb | 317 | 317/317 | 317/317 | 66 / 87 of 400 |
| *HG005, additional check (not in the aggregate, earlier k-mer index)* | chr20:4.0-4.6Mb | 400 | 400/400 | 400/400 | 63 of 400 |

Over the four chr20:3.0-3.6 Mb rows, 1,452 of 1,452 sites (100%) are recovered with the FM-index
and 1,451 of 1,452 (99.93%) by position alone. The extra HG005 window is not part of this total.

---

# Cross-cutting figures

| Number | Value | Status | Source |
|---|---|---|---|
| PgRC2 sequence-only margin | +1.88%, 6 wins and 1 loss | Cited | `docs/claim1/mechanism_insight_claim1.md` |
| PgRC2 speed and memory against G_CAPSUL | about 1.7 times slower and 2.5 times heavier | Cited | same |
| Ablation, collapse and re-placement | F1 0.431 to 0.888 (HG002, chr20) | Cited | manuscript Results, `docs/claim2/` |
| `contig_spans` archive cost | 0 bytes (18,282,397 bytes with and without) | Re-derived | `docs/claim3/t31_export_claim3.md` |
| FM-index size against archive | 1.9 to 3.5 times (0.9 to 2.6 MB per individual), against 7 to 13 times for the k-mer index it replaces | Re-derived | `docs/claim3/t34_exact_match_claim3.md` |

## Reproducibility across cores

Every number was produced on one 12-core machine with 82 GB of memory. Re-running the same encode
on E. coli (SRR2584863) with 2, 3, 7 and 12 cores produced the same archive each time (68,429,027
bytes), and only the wall time changed (36.17 s down to 12.38 s). The encoder chooses how many
adaptive candidates run at once from the free memory, so machine state affects the schedule and
never the selection: every candidate is always evaluated, and running them one at a time or four
at a time gives identical archives. Decoded reads were MD5-identical at 2, 7 and 12 cores and
matched the original FASTQ.

---

# Verification notes

The tables were verified in two tiers. Ten tables were checked in total. Four were re-executed
against an independent measurement (multi-allelic, tetraploid, the HG002 row of exact-match
recall, and the HG002 rows of SNV and indel calling). The rest were checked by parsing the run log
and comparing it with the CSV. The second tier is weaker, because the log line and the CSV row are
written from the same variables in the same run.

**Defects found and corrected**

1. **Multi-allelic recovery was first reported as 21 of 26.** The denominator was re-derived from
   the GIAB HG002 truth set (952 multi-allelic truth sites, 26 of them SNV-only), which matched.
   Three re-runs, on identical reads, at the preceding commit and with the original working
   directory traced, all gave 17 of 26. The 21 came from a stale 5 Mb window result carried
   forward as the full-chromosome number. It was corrected to 17 of 26 everywhere.
2. **Exact-match and locus retrieval need `CAPS_PILEUP=1` when the query sidecar is built.**
   Without it `query` returns the consensus instead of the reads and the coordinate-only arm scores
   zero by construction. The requirement is documented in
   `docs/claim3/t34_t35_multi_individual_claim3.md`.
3. **Two ignore patterns hid evidence files.** `*.log` hid reproduction logs and `*.vcf` hid the
   truth file used for the multi-allelic correction. Both were fixed with targeted exceptions.
4. **The file-integrity manifest included untracked files.** A fresh clone failed with eight
   missing-file errors. The manifest was regenerated from `git ls-files`.

**Limits of the verification**

- No second physical machine was available. A fresh clone with a restricted core count was used
  instead, and it passed.
- The exact-match rows for HG003, HG004 and HG005 were not re-executed independently of the
  HG002 run.
- The verification was carried out by one person. A second reviewer and a re-execution of the
  log-checked tier remain the most useful next steps.

One incorrect number was found in ten tables. That does not show that there was only one, and the
check was not exhaustive.

**Canonical version for citation:** `v1.2.0`. The verification above was carried out at
`v1.0.2-capsule`. Later tags changed documentation, and the exact-match and locus-retrieval
results were re-derived on the current source.

---

# Superseded figures

Older documents in `docs/extras/` may quote the values on the left. Use the replacement.

| Superseded value | Measure | Replacement | Reason |
|---|---|---|---|
| 80.8% genome fraction, 90% unaligned | Export correctness | 98.698% genome fraction | Measured on an export that concatenated all contigs into two records |
| 345 of 400 | Locus retrieval | 1,452 of 1,452 | Measured with an earlier query that had a strand error and covered one individual |
| 0 of 400 | Coordinate-only retrieval | a small fraction of sites | The query emitted the consensus instead of the reads |
| 21 of 26 | Multi-allelic recovery | 17 of 26 | No supporting raw log, see the verification notes |
| 5 of 111 and 11 of 18 | Multi-allelic recovery | 17 of 26 | One used a check that cannot score indel sites, and one was not reproducible |
| 0.637 and 0.666 | Heterozygous indel F1 | 0.621 | Measured on 8 chromosome 20 windows |
| 0.890 | Heterozygous SNV F1 | 0.876 | Measured before the 19-dataset sweep |
| 0.836 | Tetraploid SNV F1 | 0.897 | Measured before the final sweep |
| 555 to 656 times | Export speed-up | 129 to 784 times | Measured before the final sweep |
| 1.0000 recall without an index | Exact-match recall | 0.60 to 0.85 archive only, 1.0000 with the FM-index | The earlier figure held only for HG002 |
