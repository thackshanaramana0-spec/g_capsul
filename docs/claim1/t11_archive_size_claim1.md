---
Date: 2026-09-19
Title: T1.1 — Lossless Whole-FASTQ Archive Size, 19 Datasets
Purpose: Complete, verified record of T1.1's numbers, independently
  recomputed from raw data this session, plus the losslessness
  verification methodology behind the "LOSSLESS" column.
When to refer to this file: Writing/checking T1.1's table or text; citing
  the 19/19, −6.03%, −43.26% aggregate figures; explaining how
  losslessness was actually verified (not just claimed).
Keywords: T1.1, archive size, SPRING, Genozip, lossless, 19 datasets,
  aggregate margin, byte-identical
---

# T1.1 — Archive Size

## What it measures

Lossless whole-FASTQ archive size (sequence + names + quality + line-3 mode
+ read order — everything needed to reconstruct the exact original file),
across 19 real datasets spanning seven kingdoms (human, bacteria, archaea,
protista, fungi, virus, plantae), vs SPRING and Genozip.

## Results — independently recomputed from `results/claim1/claim1_T1.1_T1.2.csv`, not just cited

**Aggregate (recomputed this session by summing all 19 datasets'
`archive_bytes` per tool from the raw CSV and deriving the margin
directly)**:

| | vs SPRING | vs Genozip |
|---|---|---|
| Datasets won | **19 / 19** | **19 / 19** |
| Aggregate margin | **−6.03%** | **−43.26%** |
| Losslessness | **57 / 57** archives (19 datasets × 3 tools) verified LOSSLESS |

Exact match to `BTR_NOTES.md`'s headline figure — confirmed independently, not
assumed correct because it was previously stated.

**Sample row, cross-checked byte-for-byte against the manuscript table**
(ERR5181310): raw 471,183,048 B → CAPSULE 8,679,796 B (8.68 MB) vs SPRING
9,984,000 B (9.98 MB) vs Genozip 9,222,378 B (9.22 MB). Manuscript's T1.1
row for this dataset: `8.68 | 9.98 | 9.22` — exact match.

Full 19-row table: `results/claim1/claim1_T1.1_T1.2.csv`, `archive_bytes`
column per tool, or the manuscript's own T1.1 table (lines 364-393).

## What "LOSSLESS" actually means here, verified — not just a column label

The manuscript states (verified against its own text): *"Every archive was
independently decoded and verified byte-identical to the source FASTQ
before being counted."* This is a real, methodologically load-bearing
claim: the CSV's `lossless` column reads `LOSSLESS` for all 57 rows (19
datasets × 3 tools), meaning **every single archive**, including SPRING's
and Genozip's own outputs, was actually decoded and diffed against the
original — not assumed correct because the tool claims lossless
compression. This matters because Claim 1's own project history (see
`code_mapping_claim1.md` and `BTR_NOTES.md` §6.1/6.3) records that this
project once shipped incomplete archive-decode verification and had to
correct four silent data-loss bugs found only by actually decoding
archives rather than trusting dump-level checks — the discipline behind
this column is itself a real, previously-tested lesson, not a default
assumption.

## Why 19/19 with a mixed margin matters more than a single average

Stating "19/19 wins" alongside the two separate margins (−6.03% vs SPRING,
a competitive general compressor; −43.26% vs Genozip, a much larger gap) is
more honest and more informative than a single blended percentage across
both comparators — the two margins tell different stories (a close win
against the stronger comparator, a decisive one against the other) and
should be kept separate in any citation, not averaged.
