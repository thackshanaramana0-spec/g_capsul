---
Date: 2026-09-19
Title: Claim 2 (FAITHFUL) — Overview, What It Is And Is Not
Purpose: Entry point for understanding what Claim 2 claims, what it does
  not, and how T2.1-T2.5 relate, before reading any single table's detail
  file.
When to refer to this file: Before writing/reviewing any Claim 2 paper
  section; answering "is Claim 2 locked"; deciding what belongs in main
  tables vs discussion.
Keywords: FAITHFUL, Claim 2, variant calling, het-SNV, het-indel,
  multi-allelic, tetraploid, DiscoSNP++, Kmer2SNP, F1, overview, scope
---

# Claim 2 — FAITHFUL. What it is, what it is not.

## What Claim 2 actually claims, in one sentence

**The same pseudogenome structure retained for compression (Claim 1) still
contains enough biological signal to call heterozygous variants —
reference-free, from the archive alone — without returning to the original
FASTQ or building a second, independent representation.**

## What Claim 2 is NOT

- **Not a claim to beat every variant caller in existence.** The comparator
  set is deliberately narrow and justified: per `refer_paper_docs/HET_INDEL_SOTA.md`'s
  literature survey, DiscoSNP++ is the only other general-purpose
  reference-free caller that handles heterozygous indels from a single
  diploid sample with no reference at all. Everything else surveyed
  (DeepVariant, GATK, Strelka2, ska lo, eBWT2SNP) either requires a
  reference or doesn't call indels. Kmer2SNP is SNV-only by construction —
  included in T2.1, correctly marked NOT_APPLICABLE in T2.3 (indel), not a
  missing measurement.
- **Not uniformly a win.** HG005 is the honest exception across T2.1 and
  T2.3 — DiscoSNP++ wins F1 on that individual in both tables. Stated
  plainly in both the CSV and the manuscript text, not hidden.
- **Not simulated data anywhere.** T2.5 (tetraploid) is real HG003+HG004
  reads concatenated — "nothing at the sequence level is simulated, only
  the ploidy relationship is constructed" (manuscript's own wording,
  verified against the CSV).

## The five tables, one line each

| Table | Question | Comparators | Result |
|---|---|---|---|
| T2.1 | Het-SNV calling, 4 individuals, full chr20 | DiscoSNP++, Kmer2SNP | WIN on 3/4 individuals, mean F1 0.876 vs 0.853/0.475 |
| T2.2 | F1 as a function of sequencing depth (HG002) | none (sensitivity curve) | Not a win/loss table — shows where the method degrades |
| T2.3 | Het-indel calling, same 4 individuals | DiscoSNP++ (Kmer2SNP N/A) | WIN on 3/4, mean F1 0.621 vs 0.591 |
| T2.4 | Multi-allelic SNV sites, complete chr20 census | DiscoSNP++ | WIN, 17/26 vs 0/26 — structural, not marginal |
| T2.5 | Tetraploid SNV+indel, real HG003+HG004 mix | DiscoSNP++ | WIN both classes, SNV 0.897 vs 0.782, indel 0.597 vs 0.553 |

## Verification status of this overview and its 7 companion files

Every number above and in the detail files was cross-checked this session
directly against `results/claim2/claim2_T2.*.csv` **and** the live
manuscript tables in `capsul_paper/capsul_manuscript.tex` (lines 499-660+),
number by number, not sampled. Zero discrepancies found — Claim 2's
manuscript tables are current and correct, unlike Claim 3's T3.1/T3.4
tables which needed correction this session (see `refer_paper_docs/claim3/`).

## How the tables relate to each other

1. **T2.1 is the baseline capability**: does archive-native calling work at
   all, on the simplest variant class (biallelic het-SNV).
2. **T2.2 is a robustness check on T2.1**, not a separate claim — same
   individual, same method, varying only depth, showing the method
   degrades gracefully rather than failing a cliff.
3. **T2.3 repeats T2.1's structure for a harder variant class** (indels),
   using the same scoring pipeline, same comparator, same four individuals
   — a controlled generalization, not a new method.
4. **T2.4 is the sharpest evidence for the paper's mechanism** (see
   `mechanism_insight_claim2.md`): it tests a case ordinary biallelic
   calling cannot even represent — a locus with more than one alternate
   allele. This is the SAME reconciliation mechanism Claim 3/T3.5 uses for
   locus retrieval, applied here to variant calling instead.
5. **T2.5 is a generalization test along a different axis** (ploidy),
   confirming the mechanism isn't specific to diploid biology.

## File index for this folder

- `t21_snv_claim2.md`
- `t22_coverage_sweep_claim2.md`
- `t23_indel_claim2.md`
- `t24_multiallelic_claim2.md`
- `t25_tetraploid_claim2.md`
- `architecture_vs_discosnp_claim2.md` — layer-by-layer implementation
  comparison against DiscoSNP++ (18-row correspondence table: shared
  textbook algorithm vs. their design vs. ours), including the
  "assembly is the byproduct, not the graph" self-correction and one
  refuted headroom candidate (free read threading).
- `mechanism_insight_claim2.md` — how Claim 2 connects to the paper's
  central insight (the same allele-splitting/reconciliation mechanism as
  Claim 3/T3.5), and what a fair-comparison honest reading of DiscoSNP++
  and Kmer2SNP's own capabilities actually supports.
- `code_mapping_claim2.md` — exact function/line references in
  `include/caps_caller.h` and `src/decoder.cpp` for every table,
  and the exact source file for every number cited.

## Source-of-truth rule for this whole folder

Where any file here disagrees with `results/`, the result file
wins. This folder is navigation/synthesis, not a replacement source of
truth.
