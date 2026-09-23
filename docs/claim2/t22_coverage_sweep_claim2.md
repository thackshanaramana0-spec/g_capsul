---
Date: 2026-09-19
Title: T2.2 — Het-SNV F1 as a Function of Sequencing Depth (HG002)
Purpose: Complete, verified record of T2.2's numbers and why it is a
  sensitivity curve, not a win/loss table.
When to refer to this file: Writing/checking T2.2's table or text;
  answering "why does T2.2 have no comparator column."
Keywords: T2.2, coverage sweep, depth, sensitivity, HG002, seqtk
---

# T2.2 — Coverage Sweep

## What it measures

Het-SNV F1 on HG002 as sequencing depth varies (10x, 15x, 20x, 30x), each
depth an **independent full re-compress and re-call** from reads
downsampled with `seqtk sample` — depth is the only variable that changes
row to row. The 30x row is not re-run; it is the identical measurement
carried from T2.1 (same individual, same method), not duplicated work.

## Results — verified against `results/claim2/claim2_T2.2_coverage_sweep.csv`

| Depth | Reads | Archive bytes | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|---|
| 10x | 4,201,737 | 203,851,561 | 15,308 | 2,933 | 29,267 | 0.839 | 0.343 | 0.487 |
| 15x | 6,301,720 | 296,083,050 | 25,322 | 2,059 | 19,253 | 0.925 | 0.568 | 0.704 |
| 20x | 8,402,407 | 388,576,775 | 30,792 | 1,903 | 13,783 | 0.942 | 0.691 | 0.797 |
| 30x | 12,604,917 | 573,767,964 | 37,011 | 1,721 | 7,564 | 0.956 | 0.830 | **0.888** |

## Why this table has no comparator column, by design

This is a **sensitivity/degradation curve** for one method, not a
head-to-head comparison. Its purpose is to show *how* performance depends on
input coverage — where recall is the bottleneck at low depth (0.343 at 10x,
climbing to 0.830 at 30x, while precision is already reasonably high even
at 10x, 0.839) — not to win or lose against another tool at each depth
point. Framing it as a "T2.2 vs DiscoSNP++ at four depths" table would
require running the same sweep for DiscoSNP++, which was not done and is
not claimed.

## What the shape of the curve says, honestly

Recall is the limiting factor at low coverage, not precision — consistent
with fewer reads meaning fewer overlapping-read opportunities to detect an
alternate allele in the first place (candidate discovery is coverage-bound,
not a specificity problem). This is the same "once an alternative sequence
structure has not been proposed, additional filtering cannot recover the
missed event" limitation the manuscript states for indels in T2.3 — the
coverage sweep is direct, measured evidence of that same limitation's shape
for SNVs.
