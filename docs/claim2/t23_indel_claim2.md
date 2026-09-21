---
Date: 2026-09-19
Title: T2.3 — Heterozygous Indel Calling vs DiscoSNP++
Purpose: Complete, verified record of T2.3's numbers, the CAPS_CALL_INDELS
  flag dependency, and the withdrawn-numbers history this table replaced.
When to refer to this file: Writing/checking T2.3's table or text; citing
  any het-indel F1 number; verifying the CAPS_CALL_INDELS=1 requirement.
Keywords: T2.3, het-indel, DiscoSNP++, Kmer2SNP, CAPS_CALL_INDELS,
  gapped_indel_scan, withdrawn numbers
---

# T2.3 — Heterozygous Indel Calling

## What it measures

Same scoring pipeline as T2.1, same four GIAB individuals, called from the
archive, but for **indels** instead of SNVs. Kmer2SNP is correctly marked
`NOT_APPLICABLE_SNP_ONLY` in the CSV — it has no indel model by
construction, not a missing measurement.

## Results — verified against `benchmark/results/claim2_T2.3_indel.csv`, exact match to manuscript lines 582-610

| Individual | Ours F1 | DiscoSNP++ F1 |
|---|---|---|
| HG002 | **0.620** | 0.576 |
| HG003 | **0.637** | 0.587 |
| HG004 | **0.632** | 0.595 |
| HG005 | 0.593 | **0.605** |
| **Mean** | **0.621** | 0.591 |

Every TP/FP/FN count and every derived rate cross-checked line by line
against the manuscript this session — exact match, no discrepancy.

## The pattern, and what it means (verified against the manuscript's own stated interpretation)

DiscoSNP++ maintains higher precision throughout (0.915-0.927 vs our
0.810-0.843); ours recovers more true positives on the three winning
individuals. **The manuscript's own stated explanation, checked and
correct**: "the remaining indel limitation lies largely in candidate
discovery... once an alternative sequence structure has not been proposed,
additional filtering cannot recover the missed event." This is consistent
with the mechanism file's discussion of `gapped_indel_scan` — indel
detection depends on finding a structural gap in the alignment scan, and a
missed gap cannot be recovered downstream by any scoring adjustment.

## Same honest exception as T2.1 — HG005

DiscoSNP++ wins F1 on HG005 again (0.605 vs 0.593), consistent with T2.1's
pattern on the same individual. Not hidden, not averaged away.

## Critical flag dependency — every number above requires `CAPS_CALL_INDELS=1`

Verified in `PIPELINE.md`'s flags table (committed this session, but the
flag itself and its necessity predate tonight): without
`CAPS_CALL_INDELS=1`, `capsule_decode call` silently takes the graph-only
SNV path instead of erroring — producing a different, narrower
configuration, not a failure, and definitely not the numbers above. This is
the single most important flag for reproducing T2.3 (and T2.4, T2.5 — see
those files).

## Withdrawn numbers — history, kept per this project's own retraction rule

`docs/CLAIM2_FINAL_VERDICT.md` records that earlier het-indel figures
(0.637, then 0.666, both measured on 8 chr20 WINDOWS, not full chromosomes)
are **superseded** by the final 0.621 figure above (full chr20, 4
individuals, from the 2026-09-09/10 sweep). The 0.621 figure in this file
and in the manuscript is the final, correct one — do not cite 0.637 or
0.666 anywhere.
