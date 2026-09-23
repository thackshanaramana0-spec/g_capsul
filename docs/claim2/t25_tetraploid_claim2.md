---
Date: 2026-09-19
Title: T2.5 — Tetraploid Calling, Real HG003+HG004 Mix vs DiscoSNP++
Purpose: Complete, verified record of T2.5's numbers, why nothing at the
  sequence level is simulated, and what this table generalizes about the
  method (ploidy, not just diploid biology).
When to refer to this file: Writing/checking T2.5's table or text; citing
  tetraploid F1 numbers; explaining the truth-set construction.
Keywords: T2.5, tetraploid, ploidy, CAPS_PLOIDY, HG003, HG004, real reads,
  not simulated, Cooke et al.
---

# T2.5 — Tetraploid Calling

## What it measures

SNV and indel calling at ploidy 4, constructed by concatenating **real**
HG003 and HG004 reads (method following Cooke et al. 2022) over
chr20:3,000,000-3,400,000 (400kb), truth = the union of both individuals'
own real GIAB calls over that region. **Nothing at the sequence level is
simulated — only the ploidy relationship is constructed** (manuscript's own
wording, verified accurate: the reads themselves are real sequencing data
from two real people; what's artificial is only the act of treating their
combined reads as one 4-ploid sample).

## Results — verified against `results/claim2/claim2_T2.5_tetraploid.csv`, exact match to manuscript

| Class | Ours F1 | DiscoSNP++ F1 |
|---|---|---|
| SNV | **0.897** | 0.782 |
| Indel | **0.597** | 0.553 |

Full counts (both classes, same CSV): SNV TP 472/FP 7/FN 101 (ours) vs TP
373/FP 8/FN 200 (DiscoSNP++) — DiscoSNP++'s much higher FN count (200 vs
101) is the main driver of its lower recall (0.651 vs our 0.824) and hence
lower F1, at comparable precision (0.979 vs 0.985). Indel: TP 51/FP 8/FN 61
(ours) vs TP 44/FP 3/FN 68 (DiscoSNP++) — same pattern, DiscoSNP++ higher
precision, ours higher recall.

## Why this table matters beyond "another win"

T2.1-T2.4 all operate at the standard diploid (ploidy 2) setting. T2.5 is
the one test of whether the method's ploidy parameter is a real, working
mechanism rather than a diploid-only implementation with a cosmetic
`CAPS_PLOIDY` flag. **Both wins here (SNV and indel) demonstrate the
multi-allelic reconciliation machinery generalizes along the ploidy axis,
not just at the fixed diploid case** — directly relevant to T2.4's
multi-allelic result, since a tetraploid sample is structurally similar to
a multi-allelic site: more than two haplotype copies must be reconciled at
once. See `mechanism_insight_claim2.md`.

## Scope, stated honestly

This is a **400kb region**, not a full chromosome — smaller in scope than
T2.1/T2.3's full-chr20 tests. The manuscript states this precisely via the
explicit coordinate range in its caption; do not imply full-chromosome
scale for this specific table when citing it. It remains a real-data test
(no simulated sequence), which is the more important scope claim to
preserve.
