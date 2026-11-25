# T2.1 — Heterozygous SNV Calling

> Record of the T2.1 numbers and comparators, including the one dataset where G_CAPSUL does not win (HG005).

## What it measures

Heterozygous SNV calls made directly from the archive (no FASTQ re-read, no
reference genome during candidate discovery — reference coordinates are
introduced only after calling, to score against GIAB truth), across four
GIAB individuals, full chr20, ~30x, against DiscoSNP++ and Kmer2SNP. Scored
with `rtg vcfeval --squash-ploidy` against GIAB v4.2.1 confident regions.

## Results — verified against `results/claim2/claim2_T2.1_snv.csv`, exact match to manuscript

| Individual | Ours F1 | DiscoSNP++ F1 | Kmer2SNP F1 |
|---|---|---|---|
| HG002 | **0.888** | 0.847 | 0.464 |
| HG003 | **0.891** | 0.849 | 0.460 |
| HG004 | **0.891** | 0.851 | 0.469 |
| HG005 | 0.834 | **0.863** | 0.505 |
| **Mean** | **0.876** | 0.853 | 0.475 |

Precision/recall underlying the F1s (all from the same CSV): ours is
recall-favoring relative to DiscoSNP++ on the three winning individuals
(e.g. HG002: ours 0.956P/0.830R vs DiscoSNP++ 0.952P/0.763R) — we recover
more true positives at a small precision cost. Kmer2SNP is the opposite
extreme: very high precision (0.978-0.984) but very low recall (0.30-0.34),
consistent with a SNV-only, conservative-candidate-generation design.

## The exception — HG005

DiscoSNP++ wins F1 on HG005 (0.863 vs 0.834) — the only individual where
this happens in T2.1. Precision is the driver: ours drops to 0.899 on
HG005 (vs 0.956-0.958 on the other three) while recall stays comparable
(0.777 vs 0.830-0.835). The CSV and the manuscript both report it, and the row is not averaged away.

## Why DiscoSNP++ and Kmer2SNP are the right comparators, not a weak set

According to the heterozygous-indel literature survey:
DiscoSNP++ is the only other general-purpose reference-free caller handling
heterozygous variants from a single diploid sample with zero reference.
Everything else considered (DeepVariant, GATK, Strelka2, ska lo, eBWT2SNP)
requires a reference genome, which would make the comparison apples-to-
oranges against a reference-free method. Kmer2SNP is included because it is
also reference-free and SNV-capable, giving a second, independent
reference-free baseline — even though its very low recall means it is not
competitive, including it avoids selecting only the stronger
comparator alone.

## What this table does and does not establish

Establishes: archive-native calling recovers real heterozygous SNVs at
competitive-to-superior F1 against the field's own reference-free
comparators, on real GIAB individuals, full chromosome scale — not a
window or a subsample.

Does not establish: performance at lower coverage (see
`t22_coverage_sweep_claim2.md`), on indels (see `t23_indel_claim2.md`), or
why archive-native calling works mechanistically (see
`mechanism_insight_claim2.md`).
