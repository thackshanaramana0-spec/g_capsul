# Code-to-table mapping (Claim 2)

This file maps the Claim 2 operations (T2.1 to T2.5) to the code that performs them and lists the
source file of every number cited in this folder. The line numbers are from `include/caps_caller.h`
(5,970 lines) and `src/decoder.cpp`. They drift as the files change, so search for the function
name. The file maps entry points and structure, and it does not claim a line-by-line read of the
caller. The caller is described in [`../technical_architecture.md`](../technical_architecture.md),
section 7.

## How to read this file

Each result table of the manuscript is produced by a specific part of the program. This file links the two in three steps:

1. **Summary table.** For each result it gives the operation, the code that performs it, and the result file that holds the numbers. Start here.
2. **Walkthrough.** The sections after the summary describe each piece of code in a few sentences, so that you can find the function and understand what it does.
3. **Number sources.** The last two sections list where every number in this folder comes from, and which numbers were recomputed and which are cited from earlier work.

Names in backticks are functions, variables or environment variables. A line number such as `decoder.cpp:528` is the line in the current source. Line numbers move when the code is edited, so use the function name to search for the code if a line no longer matches.

## Summary

| Result | What the code does | Where | Result file |
|---|---|---|---|
| T2.1 SNVs | Reconciles the contigs, then finds bubbles that mark heterozygous SNVs | `collapse_contigs`, `build_substrate`, `extract_snv_bubble`, `run_variant_call` in `include/caps_caller.h` | `results/claim2/claim2_T2.1_snv.csv` |
| T2.2 depth sweep | The same code on HG002 reads downsampled to 10x, 15x, 20x and 30x | the same functions | `results/claim2/claim2_T2.2_coverage_sweep.csv` |
| T2.3 indels | Adds the read-level gapped scan and the indel bubbles, switched on by `CAPS_CALL_INDELS=1` | `gapped_indel_scan`, `extract_bubble`, `scan_pair` | `results/claim2/claim2_T2.3_indel.csv` |
| T2.4 multi-allelic sites | Emits more than one alternate allele at a site, within the configured ploidy | the `MAF_K` and `n_multi_overploidy` logic in `run_variant_call` | `results/claim2/claim2_T2.4_multiallelic.csv` |
| T2.5 tetraploid | The same caller with ploidy 4, set by `CAPS_PLOIDY=4` | the ploidy handling in `run_variant_call` | `results/claim2/claim2_T2.5_tetraploid.csv` |

## Walkthrough

### Entry point

- `src/decoder.cpp:2512` dispatches `capsule_decode call <in.capsule> <out.vcf> [workdir]` to
  `capsule_call_from_archive`, defined at `decoder.cpp:2233`.
- The call reads a stored `.capsule` file. It does not use in-memory state from compression, which
  confirms at the code level that calling is archive-native and does not read the FASTQ again.

### Core functions in `include/caps_caller.h`

- **`extract_bubble` (line 258) and `extract_snv_bubble` (line 427)** find divergent regions between
  two contig sequences (bubbles), the signal that a heterozygous site produces when its two alleles
  sit on different pseudogenome segments.
- **`scan_pair` (line 373)** runs bubble extraction across a pair of aligned sequences.
- **`build_substrate` (line 791)** is the reconciliation step. It builds the calling substrate from
  the contig sequences and `CallData`. It is very probably the function that the manuscript's
  ablation toggles (F1 0.431 without collapse and re-placement, 0.888 with both). This is inferred
  from names and from the commit history around `1fec5e0` ("Close the indel_pass accounting:
  build_substrate is 54% and is not waste"), and the ablation was not re-run.
- **`collapse_contigs` (line 625)** is the contig-merging step, referred to as "collapse" in the
  ablation.
- **`gapped_indel_scan` (line 707)** is the read-level gapped alignment scan that produces the indel
  evidence behind T2.3.
- **`run_variant_call` (line 1188)** is the dominant function of the file. It orchestrates the
  pipeline: SNV pileup substrate, indel and bubble substrate, and ploidy-aware emission.
- **Ploidy.** `CAPS_PLOIDY` is read at line 1262 (`if (v >= 2 && v <= 4) PLOIDY = v`) and used in
  the multi-allelic frequency threshold at line 2262 (`MAF_K = MAF * 2.0 / PLOIDY`). Candidates
  that exceed the configured ploidy are counted in `n_multi_overploidy` (from line 3438). This is
  the mechanism behind T2.4 (ploidy 2, multi-allelic within that ploidy) and T2.5 (ploidy 4, set by
  `CAPS_PLOIDY=4`).

### Flags the tables depend on

- `CAPS_CALL_INDELS=1` is required for T2.3, T2.4 and T2.5. Without it, `capsule_decode call`
  takes the narrower graph-only SNV path without an error. The comment at `decoder.cpp:2241`
  documents it.
- `CAPS_PLOIDY=N` is required for T2.5 (`N=4`) and defaults to 2.

The full list of environment variables is in [`../configuration_reference.md`](../configuration_reference.md).

## Number sources

| Number | Source |
|---|---|
| T2.1 (het-SNV, all rows) | `results/claim2/claim2_T2.1_snv.csv` |
| T2.2 (coverage sweep, all rows) | `results/claim2/claim2_T2.2_coverage_sweep.csv` |
| T2.3 (het-indel, all rows) | `results/claim2/claim2_T2.3_indel.csv` |
| T2.4 (multi-allelic, 17/26) | `results/claim2/claim2_T2.4_multiallelic.csv` |
| T2.5 (tetraploid, all rows) | `results/claim2/claim2_T2.5_tetraploid.csv` |
| Withdrawn multi-allelic numbers (5/111, 11/18, 21/26) and why | An internal reproduction log that is not part of this repository. The current figure (17 of 26) is in `results/claim2/claim2_T2.4_multiallelic.csv` |
| Withdrawn het-indel numbers (0.637, 0.666) and why | the development notes (the Claim 2 final-verdict note, marked superseded) |
| Claim 2 development history | the first Claim 2 results note, the second Claim 2 results note, the Claim 2 tables and indel-scan note, the het-indel re-scan note |
| Literature survey justifying DiscoSNP++/Kmer2SNP as comparators | the het-indel literature survey |
| Ablation numbers (F1 0.431 vs 0.888) | Manuscript body text. The underlying run is recorded in the development notes |

## What was re-derived and what is cited

The manuscript tables T2.1, T2.3, T2.4 and T2.5 were compared with their CSV files and match
exactly. The mechanics of `CAPS_PLOIDY` and `CAPS_CALL_INDELS` were confirmed by reading the code.

The original run log of the F1 0.431 to 0.888 ablation is cited and was not re-run. Whether
`build_substrate` and `collapse_contigs` are precisely the functions the ablation toggled is
inferred from names and commit history. This distinction should be kept when the ablation numbers
are cited.
