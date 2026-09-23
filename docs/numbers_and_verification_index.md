---
Date: 2026-09-19
Title: Master Numbers and Verification Index — Every Table, Every Row,
  First to Last
Purpose: The complete, literal, row-by-row record of every table in the
  paper (not just headline aggregates), each row's exact source, and its
  verification status this session. The deepest single reference in this
  folder tree.
When to refer to this file: Before citing ANY number from this paper
  anywhere; when a reviewer or collaborator asks for the full data behind
  a headline figure; the single place to find a specific dataset's row in
  any table without opening a CSV.
Keywords: numbers index, master table, full data, every row, verification
  status, source files, re-derived, cross-checked, cited
---

# Master numbers and verification index — full depth

**Verification status legend**: **RE-DERIVED** = independently recomputed
from raw data this session. **CROSS-CHECKED** = manuscript number compared
directly against its CSV source, row by row, this session. **CITED** =
read from a prior, documented source, not independently re-run this
session.

---

# CLAIM 1 — COMPACT

## T1.1 — Archive size, all 19 datasets (MB). Status: RE-DERIVED (aggregate) + spot cross-checked (rows). Source: `results/claim1/claim1_T1.1_T1.2.csv`.

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

**CAPSULE is smaller than SPRING on all 19, and smaller than Genozip on all
19 — verified by inspection of every row above, not just the aggregate.**
Aggregate: 19/19 vs SPRING (−6.03%), 19/19 vs Genozip (−43.26%). Every one
of the 57 archives (19×3 tools) is `LOSSLESS` per the CSV's own column.

## T1.2 — Wall time (s), compress and decompress, same 19 datasets. Status: RE-DERIVED. Source: same CSV.

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

**CAPSULE is the fastest tool in exactly 0 of these 38 cells** (19 compress
+ 19 decompress). Genozip wins 37/38; SPRING wins 1/38 (ERR5181310
compress). Confirmed against the manuscript's own bold-marking, which
agrees cell for cell.

---

# CLAIM 2 — FAITHFUL

## T2.1 — Het-SNV calling, full counts, 4 individuals. Status: CROSS-CHECKED against manuscript. Source: `results/claim2/claim2_T2.1_snv.csv`.

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

**Mean F1: CAPSULE 0.876, DiscoSNP++ 0.853, Kmer2SNP 0.475.** CAPSULE wins
3/4 individuals; DiscoSNP++ wins HG005.

## T2.2 — Het-SNV F1 vs sequencing depth, HG002. Status: CROSS-CHECKED. Source: `results/claim2/claim2_T2.2_coverage_sweep.csv`.

| Depth | Reads | Archive bytes | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|---|
| 10x | 4,201,737 | 203,851,561 | 15,308 | 2,933 | 29,267 | 0.839 | 0.343 | 0.487 |
| 15x | 6,301,720 | 296,083,050 | 25,322 | 2,059 | 19,253 | 0.925 | 0.568 | 0.704 |
| 20x | 8,402,407 | 388,576,775 | 30,792 | 1,903 | 13,783 | 0.942 | 0.691 | 0.797 |
| 30x | 12,604,917 | 573,767,964 | 37,011 | 1,721 | 7,564 | 0.956 | 0.830 | 0.888 |

30x row is identical to T2.1's HG002/CAPSULE row (carried, not re-run, per
the manuscript's own stated methodology).

## T2.3 — Het-indel calling, full counts, 4 individuals. Status: CROSS-CHECKED. Source: `results/claim2/claim2_T2.3_indel.csv`.

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

Kmer2SNP: `NOT_APPLICABLE_SNP_ONLY` on all four individuals — no indel
model, by construction. **Mean F1: CAPSULE 0.621, DiscoSNP++ 0.591.**
CAPSULE wins 3/4; DiscoSNP++ wins HG005 (same individual as T2.1's
exception).

## T2.4 — Multi-allelic SNV sites, complete chr20 census. Status: CROSS-CHECKED. Source: `results/claim2/claim2_T2.4_multiallelic.csv`.

| Individual | Scope | Sites | CAPSULE both alleles | CAPSULE rate | DiscoSNP++ multi-allelic records emitted | DiscoSNP++ total records |
|---|---|---|---|---|---|---|
| HG002 | chr20 COMPLETE CENSUS | 26 | **17** | 65.4% | 0 | 3,989 |

Status field in CSV: `CORRECTED_20260910`. DiscoSNP++ emitted zero
multi-allelic records across its ENTIRE output on this individual, not
just these 26 sites.

## T2.5 — Tetraploid calling, HG003+HG004 real reads, chr20:3,000,000-3,400,000. Status: CROSS-CHECKED. Source: `results/claim2/claim2_T2.5_tetraploid.csv`.

| Class | Tool | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|
| SNV | CAPSULE | 472 | 7 | 101 | 0.985 | 0.824 | **0.897** |
| SNV | DiscoSNP++ | 373 | 8 | 200 | 0.979 | 0.651 | 0.782 |
| INDEL | CAPSULE | 51 | 8 | 61 | 0.864 | 0.455 | **0.597** |
| INDEL | DiscoSNP++ | 44 | 3 | 68 | 0.936 | 0.393 | 0.553 |

Both classes: CAPSULE wins.

---

# CLAIM 3 — ADDRESSABLE

## T3.1 — Export speed, 6 datasets. Status: CITED. Source: `results/claim3/claim3_T3.1_T3.2_T3.3.csv`.

| Dataset | CAPSULE (s) | SPAdes (s) | Speedup |
|---|---|---|---|
| ERR5181310 | 0.40 | 51.73 | 129.3x |
| SRR2584863 | 0.37 | 205.04 | 554.2x |
| SRR29296997 | 0.13 | 83.95 | 645.8x |
| SRR37283774 | 1.49 | 380.88 | 255.6x |
| DRR976266 | 0.70 | 548.82 | 784.0x |
| HG002 | 5.67 | 2,678.83 | 472.5x |

## T3.1 — Correctness (E. coli/SRR2584863 only, the one dataset both tools can be fairly compared on). Status: RE-DERIVED this session. Source: `docs/T3.1_CORRECTNESS_FINAL_20260919.md`.

| Metric | CAPSULE (final, shipping config) | SPAdes | Winner |
|---|---|---|---|
| Genome fraction | 98.698% | 98.346% | CAPSULE |
| Indels /100kbp | 0.25 | 0.51 | CAPSULE |
| Mismatches /100kbp | 19.40 | 2.75 | SPAdes |
| Duplication ratio | 2.021 | 1.000 | SPAdes |
| N50 | 2,054 | 139,621 | SPAdes |
| Unaligned length | 1,895,813 | 9,582 | SPAdes |
| Misassemblies | 7 | 1 | SPAdes |

## T3.2 — Coverage speed, same 6 datasets. Status: CITED. Source: `results/claim3/claim3_T3.1_T3.2_T3.3.csv`.

| Dataset | CAPSULE (s) | bwa+samtools+mosdepth (s) | Speedup |
|---|---|---|---|
| ERR5181310 | 0.15 | 4.15 | 27.7x |
| SRR2584863 | 0.33 | 7.57 | 22.9x |
| SRR29296997 | 0.11 | 2.53 | 23.0x |
| SRR37283774 | 0.56 | 15.98 | 28.5x |
| DRR976266 | 1.49 | 23.98 | 16.1x |
| HG002 | 3.04 | 164.90 | 54.2x |

## T3.3 — Range query timing, all 19 datasets, no comparator. Status: CITED. Source: `results/claim3/claim3_T3.1_T3.2_T3.3.csv`.

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

## T3.4 — Exact-match recall, 4 individuals. Status: RE-DERIVED this session. Source: `docs/T34_T35_MULTI_INDIVIDUAL_20260919.md`.

| Individual | Recall, plain query (no the completion index) | Recall, with the completion index |
|---|---|---|
| HG002 | 1.0000 | 1.0000 (already at ceiling, not re-tested with the completion index) |
| HG003 | 0.9619 | **1.0000** |
| HG004 | 0.9612 | **1.0000** |
| HG005 | 0.8563 | **1.0000** |

## T3.5 — Locus retrieval, 4 individuals, 2 loci. Status: RE-DERIVED this session. Source: `docs/T34_T35_MULTI_INDIVIDUAL_20260919.md`.

| Individual | Locus | Sites | Native (no the completion index) | With the completion index | Negative-control FP (with the completion index) |
|---|---|---|---|---|---|
| HG002 | chr20:3.0-3.6Mb | 400 | 400/400 | 400/400 | 42/400 |
| HG003 | chr20:3.0-3.6Mb | 335 | 335/335 | 335/335 | 52/400 |
| HG004 | chr20:3.0-3.6Mb | 400 | 399/400 | **400/400** (miss fixed) | 49/400 |
| HG005 | chr20:3.0-3.6Mb | 317 | 317/317 | 317/317 | 87/400 |
| HG005 | chr20:4.0-4.6Mb (new locus) | 400 | 400/400 | 400/400 | 63/400 |

**Aggregate: 1,452/1,452 (100%) with the completion index; 1,451/1,452 (99.93%) native.**

---

# Cross-cutting figures

| Number | Value | Status | Source |
|---|---|---|---|
| PgRC2 sequence-only margin | +1.88%, 6 wins/1 loss | CITED | `CLAUDE.md`, `docs/CLAIM1_FINAL_VERDICT.md` §42, `docs/FINAL_HEADROOM.md` |
| PgRC2 speed/RAM vs G_CAPSUL | ~1.7x slower, ~2.5x heavier | CITED | same |
| Ablation: allele-split collapse effect | F1 0.431 → 0.888 | CITED, not re-run this session | manuscript Results, ~line 306 |
| `contig_spans` archive cost | 0 bytes (18,282,397 B identical with/without) | RE-DERIVED this session | `docs/T3.1_CORRECTNESS_FINAL_20260919.md` |
| the completion index index size vs archive | 7-13x | RE-DERIVED this session | `docs/T34_T35_MULTI_INDIVIDUAL_20260919.md` |

## Cross-machine / cross-core reproducibility — verified, not assumed

Every number in the paper was produced on one 12-core, 82 GB machine.
**Confirmed** (`docs/REPRODUCIBILITY.md`, tested 2026-09-10): re-running the
identical encode on E. coli (SRR2584863) at 2, 3, 7, and 12 cores produces
a **byte-identical archive** at every core count (68,429,027 B, unchanged);
only wall time varies (36.17 s → 12.38 s, 2.9× range). The encoder derives
how many adaptive candidates run concurrently from measured free RAM, so
machine state enters the *schedule* but never the *selection* — every
candidate is always evaluated regardless of concurrency, confirmed
identical archives whether candidates run sequentially (`K=1`) or in
parallel (`K=4`). Decoded reads are MD5-identical at 2, 7, and 12 cores and
lossless against the original FASTQ. **Answer: yes for output, no for wall
time — which is the correct, expected behavior**, not a caveat to qualify
away.

## The project's own forensic audit — what it checked, what it found, what it admits it cannot establish

Verified against root-level `AUDIT.md` (2026-09-10, tag `v1.0.2-capsule`,
the canonical frozen state — read in full this pass, not previously in
`refer_paper_docs/`). This is the single most rigorous self-verification
document in the project and belongs here as the methodological backbone
for every RE-DERIVED/CROSS-CHECKED status above.

**Two-tier verification, and why the distinction matters.** Ten tables
were checked; only four were re-executed against an independent
measurement (T2.4, T2.5, T3.4's HG002 row, T2.1/T2.3's HG002 rows) — the
rest were checked by parsing the run log and diffing against the CSV
("transcription-checked"). The audit is explicit that transcription
checking is structurally weaker: the log line and the CSV row are written
from the *same shell variables in the same execution*, so a wrong
underlying measurement would produce a matching wrong number in both
places and the check would report clean. **The one real defect this audit
found (T2.4) was in the re-executed tier and could not have been caught by
transcription checking** — stated as the honest measure of what the
weaker tier is actually worth, not asserted as a flaw glossed over.

**Four real defects found, each with its own forensic trail:**

1. **T2.4's 21/26 had no supporting raw log anywhere on the server.**
   Investigated by re-deriving the denominator independently from the raw
   GIAB HG002 BAM (952 multi-allelic truth sites, 26 SNV-only — matched
   the published denominator exactly), then re-running the caller three
   separate times under three different hypotheses to rule each out in
   turn: suspected nondeterminism (re-ran on identical already-streamed
   reads — 17/26 again, every intermediate count matching); suspected an
   intervening encoder fix (built at the prior commit via `git worktree`,
   re-ran — 17/26 a third time, the fix confirmed inert); finally traced
   the original working directory and found its `regions.bed` spans a
   5 Mb window, not the full chromosome, and its truth file has 111 sites
   (7 SNV-only) — matching this project's own *earlier, already-superseded*
   window figures. **21/26 was a stale window result mistakenly carried
   forward as if it were the full-chromosome number.** Corrected to 17/26,
   propagated to all 25 files that cited it.
2. **T3.4 needed an undocumented required flag, and finding it took three
   false failures.** Reproducing T3.4 gave `coordinate = 0` three times,
   matching the *old*, already-acknowledged-wrong figure — each time a
   confident wrong diagnosis. The actual cause: `CAPS_PILEUP=1` must be
   passed to `capsule_decode index` before building the query sidecar;
   without it the sidecar carries placements but not each read's
   deviations, so `query` emits the consensus rather than the reads and
   the coordinate arm scores 0 *by construction*, not by measurement. This
   flag appeared in **zero** documentation files at the time — a reviewer
   following the documented steps would have hit the identical wall and
   drawn the identical wrong conclusion, three times, exactly as happened
   here. Now documented in `REPRODUCE_EVERYTHING.md` and `RESULT_CODE.md`.
3. **`.gitignore` was silently excluding the audit's own evidence — twice.**
   First `*.log` excluded preserved reproduction logs already cited as
   "in the repository" when they were not; the identical failure mode was
   then found again for `*.vcf`, which had excluded the T2.4 correction's
   own truth file. Both fixed with targeted exceptions.
4. **The file-integrity manifest itself was generated wrong**, found only
   by actually performing a fresh clone (the exact scenario the manifest
   exists to support) — it failed with 8 `No such file` errors because the
   original manifest was built with a plain `find` that pulled in a
   concurrent agent's untracked files. Regenerated from `git ls-files` only.

**What the audit explicitly does NOT claim to have established** — stated
plainly rather than omitted: no second physical machine exists to test on
(a fresh clone + restricted core count was the closest substitute, and it
passed, but the audit itself calls this "not the same as a genuinely
independent machine"); T3.4's HG003/HG004/HG005 rows were never re-run,
only HG002 — the other three are "confirmed by inheritance," not
execution; no independent human reviewer exists — single session, single
party throughout.

**The self-assessed detection rate, stated as a number, not a feeling**:
"one bad number in ten tables was found; the honest inference is not
'there was exactly one' but 'the process that produced one can produce
another, and the audit that caught it is not exhaustive.'" The document's
own stated highest-value next step is re-executing the transcription-only
tier and getting a second reviewer — not more documentation.

**Canonical version, for citation**: `v1.0.2-capsule`. Five tags were cut
the same day, each named as if final; four of the five are byte-identical
in code and results to their predecessor, differing only in documentation.
The last commit that changed any code or result was `7056c4f` (the T2.4
correction) — everything after it, including the tag currently in use, is
prose-only.

---

# How to use this file

Find the table, find the row, cite it. Every number above was traced to a
specific CSV cell or a specific dated doc this session — nothing here was
typed from memory. If a status says CITED rather than RE-DERIVED or
CROSS-CHECKED, it means this session read and trusted a prior session's
verified work rather than re-running it — say so explicitly if a reviewer
or collaborator asks whether something was freshly verified.
