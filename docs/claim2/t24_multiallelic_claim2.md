---
Date: 2026-09-19
Title: T2.4 — Multi-Allelic SNV Sites, Complete chr20 Census
Purpose: Complete, verified record of T2.4's numbers, why it is a complete
  census (not a sample), the withdrawn-numbers history (21/26, 11/18,
  5/111), and why this is the sharpest evidence of the paper's mechanism.
When to refer to this file: Writing/checking T2.4's table or text; citing
  17/26 anywhere; explaining why DiscoSNP++ structurally scores 0.
Keywords: T2.4, multi-allelic, complete census, 17/26, CAPS_PLOIDY,
  withdrawn numbers, structural not marginal
---

# T2.4 — Multi-Allelic Sites

## What it measures

Sites where a chromosome-20 locus has **more than one** alternate allele —
a case ordinary biallelic calling cannot represent at all, not just detects
poorly. **26 is stated as the total population of such sites under this
scoring convention, not a sample** — verified against the manuscript
caption and the CSV's own `scope` field: "chr20 COMPLETE CENSUS."

## Results — verified against `benchmark/results/claim2_T2.4_multiallelic.csv`, exact match to manuscript

| Tool | Recovered | Fraction |
|---|---|---|
| G_CAPSUL | **17 / 26** | **65.4%** |
| DiscoSNP++ | 0 / 26 | 0.0% (structural) |

CSV status field: `CORRECTED_20260910`. `disco_multiallelic_records_emitted`
= 0 out of `disco_total_records` = 3,989 — DiscoSNP++ did not emit a single
multi-allelic record across its entire output on this individual, not just
missed these specific 26 sites.

## Why DiscoSNP++'s 0 is structural, not a tuning gap — do not soften this

Verified against the manuscript's own stated reasoning, and consistent with
DiscoSNP++'s documented design (bubble-based, biallelic by construction):
DiscoSNP++'s output format has no record type that represents more than one
ALT allele at one locus. This is not "DiscoSNP++ missed these 26 sites" in
the sense of a detectable-but-undetected event — it is "DiscoSNP++'s output
schema cannot express the answer," the same category of limitation as
asking a tool with no indel model to call an indel. Stating "0/26" without
this context risks implying a fixable miss; the correct framing is a
categorical capability difference.

## Withdrawn numbers — three of them, and why each was wrong

Per `CLAUDE.md` and `docs/CLAIM2_FINAL_VERDICT.md`, ALL of the following are
withdrawn and must not be cited:

- **"5/111"**: scored 104 indel-bearing sites with a single-base check that
  cannot evaluate them — a measurement-method defect, not a caller result.
- **"11/18"**: unreproducible on re-derivation.
- **"21/26"**: published once, had no supporting raw log, and could not be
  reproduced — three independent re-derivation runs all gave 17/26 instead,
  with identical intermediate counts every time. **17/26 is the verified
  answer.**

This history is worth keeping in mind for discussion/methods sections: it
demonstrates the project's actual practice of re-deriving and correcting
its own published numbers when they could not be reproduced, not just
accepting a favorable-looking figure.

## Ploidy mechanism — how the caller represents this at all

Multi-allelic emission is native VCF output via `CAPS_PLOIDY=k`
(default 2), not a post-hoc merge of separate biallelic calls — confirmed
in `include/caps_caller.h` by the presence of explicit ploidy-gated
rejection tracking (`n_multi_overploidy`, logged per-run) alongside
successful multi-allelic emission counts. See `code_mapping_claim2.md` for
exact line references.

## Connection to the paper's central mechanism — the important one

**This table is the sharpest evidence, alongside T3.5, that allele
separation at compression time is a real, general phenomenon, not specific
to heterozygous SNV calling.** A multi-allelic site is exactly the case
where more than two versions of one locus must be reconciled, not just two
— the same reconciliation logic that recovers both alleles of a simple het
site is what allows a *third* or *fourth* allele to also be represented
here. Full argument: `mechanism_insight_claim2.md`.
