# Table numbering — claim-based scheme

Renumbered 2026-09-09 so every table names the claim it belongs to. The old
flat numbering (T1..T6 plus T5.2/T5.3) appears throughout older docs and
results; this is the mapping, kept so those remain readable.

| old | new | content | produced by |
|---|---|---|---|
| T1 | **T1.1** | archive size, 19 datasets x 3 tools | `claim1_T1.1_T1.2.csv` |
| T2 | **T1.2** | compress/decompress wall + peak RAM | same CSV |
| T3 | **T2.1** | het-SNV F1 — ours vs DiscoSNP++ vs Kmer2SNP | `claim2_T2.1_snv.csv` |
| T4 | **T2.2** | coverage sweep 10/15/20/30x (HG002) | `claim2_T2.2_coverage_sweep.csv` |
| T5 | **T2.3** | het-indel F1 — ours vs DiscoSNP++ | `claim2_T2.3_indel.csv` |
| T5.2 | **T2.4** | multi-allelic sites recovered | `claim2_T2.4_multiallelic.csv` |
| T5.3 | **T2.5** | tetraploid SNV + indel (Cooke et al. 2022) | `claim2_T2.5_tetraploid.csv` |
| T6a | **T3.1** | export vs SPAdes | `claim3_T3.1_T3.2_T3.3.csv` |
| T6b | **T3.2** | coverage vs bwa+samtools+mosdepth | same CSV |
| T6c | **T3.3** | query (no competitor exists) | same CSV |

**Caution when reading old documents:** under the OLD scheme `T3` meant
het-SNV; under the NEW scheme `T3.x` means Claim 3 (addressability). A bare
"T3" in a pre-2026-09-09 document is het-SNV, i.e. today's T2.1.

Every table carries wall time, peak RAM and an output size or count, so each
claim can be audited on all three axes. Claim 3's rows are tagged with their
sub-table id in the first column.
