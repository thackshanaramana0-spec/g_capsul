# Claim 2 (FAITHFUL): reference-free variant calling from the archive

This folder documents Claim 2 of the manuscript. The pseudogenome structure that the archive keeps
for compression (Claim 1) still holds enough signal to call heterozygous variants, from the archive
alone, without a reference and without going back to the FASTQ or building a second representation.
This file states what the claim covers and does not cover and how its five tables (T2.1 to T2.5,
the result tables of Claim 2) relate. The result files in `results/claim2/` are authoritative.

## What Claim 2 does not claim

- **It does not claim to beat every variant caller.** The comparator set is narrow on purpose.
  According to the heterozygous-indel literature survey in the development notes, DiscoSNP++ is the
  only other general-purpose, reference-free caller that handles heterozygous indels from a
  single diploid sample. The other tools surveyed (DeepVariant, GATK, Strelka2, ska lo and
  eBWT2SNP) either need a reference or do not call indels. Kmer2SNP is SNV-only by construction. It
  appears in T2.1, and in T2.3 it is marked `NOT_APPLICABLE` because it has no indel model, which
  means no measurement is missing.
- **It is not a win on every individual.** HG005 is the exception in T2.1 and T2.3, where DiscoSNP++
  has the higher F1. The CSV and the manuscript text both say so.
- **It does not use simulated data.** T2.5 (tetraploid) concatenates real HG003 and HG004 reads.
  In the manuscript's words, nothing at the sequence level is simulated, and only the ploidy
  relationship is constructed.

## The five tables, one line each

| Table | Question | Comparators | Result |
|---|---|---|---|
| T2.1 | Het-SNV calling, 4 individuals, full chr20 | DiscoSNP++, Kmer2SNP | WIN on 3/4 individuals, mean F1 0.876 vs 0.853/0.475 |
| T2.2 | F1 as a function of sequencing depth (HG002) | none (sensitivity curve) | Not a win/loss table — shows where the method degrades |
| T2.3 | Het-indel calling, same 4 individuals | DiscoSNP++ (Kmer2SNP N/A) | WIN on 3/4, mean F1 0.621 vs 0.591 |
| T2.4 | Multi-allelic SNV sites, complete chr20 census | DiscoSNP++ | WIN, 17/26 vs 0/26 — structural, not marginal |
| T2.5 | Tetraploid SNV+indel, real HG003+HG004 mix | DiscoSNP++ | WIN both classes, SNV 0.897 vs 0.782, indel 0.597 vs 0.553 |

## Verification of the figures

Every number in this folder was compared, one by one, with `results/claim2/claim2_T2.*.csv` and
with the manuscript tables. No discrepancy was found.

## How the tables relate to each other

1. **T2.1 is the baseline capability**: does archive-native calling work at
   all, on the simplest variant class (biallelic het-SNV).
2. **T2.2 is a robustness check on T2.1**, not a separate claim — same
   individual, same method, varying only depth, showing the method
   degrades gracefully rather than failing a cliff.
3. **T2.3 repeats T2.1's structure for a harder variant class** (indels),
   using the same scoring pipeline, same comparator, same four individuals
   — a controlled generalization, not a new method.
4. **T2.4 is the clearest evidence for the paper's mechanism** (see
   `mechanism_insight_claim2.md`): it tests a case ordinary biallelic
   calling cannot even represent — a locus with more than one alternate
   allele. This is the SAME reconciliation mechanism Claim 3/T3.5 uses for
   locus retrieval, applied here to variant calling instead.
5. **T2.5 is a generalization test along a different axis** (ploidy),
   confirming the mechanism isn't specific to diploid biology.

## File index for this folder

- [`t21_snv_claim2.md`](t21_snv_claim2.md): heterozygous SNVs.
- [`t22_coverage_sweep_claim2.md`](t22_coverage_sweep_claim2.md): F1 against depth.
- [`t23_indel_claim2.md`](t23_indel_claim2.md): heterozygous indels.
- [`t24_multiallelic_claim2.md`](t24_multiallelic_claim2.md): multi-allelic sites.
- [`t25_tetraploid_claim2.md`](t25_tetraploid_claim2.md): tetraploid calling.
- [`architecture_vs_discosnp_claim2.md`](architecture_vs_discosnp_claim2.md): layer-by-layer
  comparison with DiscoSNP++ (an 18-row table of the shared textbook algorithm, DiscoSNP++'s design
  and G_CAPSUL's), the correction that assembly is a byproduct and not the graph, and one refuted
  speed-up (free read threading).
- [`mechanism_insight_claim2.md`](mechanism_insight_claim2.md): how Claim 2 rests on the same
  allele-splitting and reconciliation mechanism as Claim 3 (T3.5), and what DiscoSNP++ and Kmer2SNP
  can and cannot do.
- [`code_mapping_claim2.md`](code_mapping_claim2.md): the functions in `include/caps_caller.h` and
  `src/decoder.cpp` behind each table, and the source file of every number.

## Source of truth

Where a file in this folder disagrees with `results/`, the result file is correct. This folder
explains the results and does not replace them.
