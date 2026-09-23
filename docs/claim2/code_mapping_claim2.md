---
Date: 2026-09-19
Title: Claim 2 — Code-to-Table Mapping and Result File Index
Purpose: Exact function/line references in include/caps_caller.h and
  src/decoder.cpp for every Claim 2 operation, and the exact
  file path for every number cited anywhere in this folder. Traced
  directly from the actual repo this session.
When to refer to this file: Verifying a claim in another file of this
  folder against the real code; before citing any number, to find its
  exact source file.
Keywords: code mapping, caps_caller.h, capsule_call_from_archive,
  build_substrate, collapse_contigs, extract_snv_bubble,
  gapped_indel_scan, CAPS_PLOIDY, benchmark results, file index
---

# Code-to-table mapping

Line numbers from `include/caps_caller.h` (6,908 lines) and
`src/decoder.cpp`, as committed. This file maps *structure and
entry points*, verified directly by grep and targeted reads this session —
it is not a claim to have read all 6,908 lines of `caps_caller.h` (not
feasible or necessary for this purpose). If a specific line number below
has shifted, re-verify with `grep -n` before trusting it; the described
mechanism should still hold even if a line moves.

## Entry point

- `decoder.cpp:2136-2137` — dispatch: `capsule_decode call <in.capsule>
  <out.vcf> [workdir]` routes to `capsule_call_from_archive`, defined at
  `decoder.cpp:1857`.
- This is the "Claim 2 from a STORED archive" path (comment at line 2135) —
  calling reads a `.capsule` file, not live in-memory state from
  compression, confirming the "archive-native, no FASTQ re-read" claim at
  the code level, not just in prose.

## Core mechanism functions, `include/caps_caller.h`

- **`extract_bubble` / `extract_snv_bubble`** (lines ~195, ~401) — find
  divergent regions between two contig sequences ("bubbles"), the raw
  signal that a het site produces when its two alleles are separated onto
  different pseudogenomic segments.
- **`scan_pair`** (line ~347) — runs bubble extraction across a pair of
  aligned sequences.
- **`build_substrate`** (line 816) — the reconciliation mechanism. Builds
  the calling substrate from raw contig sequences plus `CallData`; this is
  the function whose presence/absence is almost certainly what the
  manuscript's ablation (F1 0.431 without collapse/re-placement, 0.888
  with both) is measuring — verify against `docs/HET_INDEL_SOTA.md` or the
  git history around commit `301d1a5` ("Close the indel_pass accounting:
  build_substrate is 54% and is not waste") before stating this as fully
  confirmed rather than strongly inferred.
- **`collapse_contigs`** (line 632) — the contig-merging step referenced by
  name in the ablation's "collapse" terminology.
- **`gapped_indel_scan`** (line 731) — read-level gapped alignment scan
  producing indel evidence; the function behind T2.3's indel calls.
- **`run_variant_call`** (line 1260 onward — the dominant function in the
  file, spanning roughly 5,600 of the file's 6,908 lines) — orchestrates
  the full pipeline: SNV pileup substrate, indel/bubble substrate,
  ploidy-aware emission.
- **Ploidy handling**: `CAPS_PLOIDY` read at `caps_caller.h:1347`
  (`if (v >= 2 && v <= 4) PLOIDY = v`), used in multi-allelic-frequency
  logic at line 2512 (`MAF_K = MAF * 2.0 / PLOIDY`), and in rejection
  accounting for multi-allelic candidates that exceed the configured ploidy
  (`n_multi_overploidy`, logged around line 4325). This is the mechanism
  behind both T2.4 (ploidy=2, multi-allelic within that ploidy) and T2.5
  (ploidy=4, via `CAPS_PLOIDY=4`).

## Flag dependency, confirmed against `PIPELINE.md`

- `CAPS_CALL_INDELS=1` — required for T2.3/T2.4/T2.5. Without it,
  `capsule_decode call` silently takes the graph-only SNV path (a
  different, narrower configuration, not an error) — documented in
  `PIPELINE.md`'s flags table (committed prior to this session, re-verified
  present and correct this session).
- `CAPS_PLOIDY=N` — required for T2.5 specifically (`N=4`); defaults to 2.

---

# Result file index — exact source for every number cited in this folder

| Number / claim | File |
|---|---|
| T2.1 (het-SNV, all rows) | `results/claim2/claim2_T2.1_snv.csv` |
| T2.2 (coverage sweep, all rows) | `results/claim2/claim2_T2.2_coverage_sweep.csv` |
| T2.3 (het-indel, all rows) | `results/claim2/claim2_T2.3_indel.csv` |
| T2.4 (multi-allelic, 17/26) | `results/claim2/claim2_T2.4_multiallelic.csv` |
| T2.5 (tetraploid, all rows) | `results/claim2/claim2_T2.5_tetraploid.csv` |
| Withdrawn multi-allelic numbers (5/111, 11/18, 21/26) and why | Internal reproduction log, not included in this public repo; the surviving, correct number (17/26) is the one reported in this folder and in `results/claim2/claim2_T2.4_multiallelic.csv` |
| Withdrawn het-indel numbers (0.637, 0.666) and why | `docs/CLAIM2_FINAL_VERDICT.md` (marked SUPERSEDED at top) |
| Full session-by-session history of Claim 2's development | `docs/CLAIM2_RESULTS.md`, `docs/CLAIM2_RESULTS_V2.md`, `docs/CLAIM2_TABLES_AND_INDEL_SCAN.md`, `docs/HET_INDEL_FRESH_SCAN.md` |
| Literature survey justifying DiscoSNP++/Kmer2SNP as comparators | `docs/HET_INDEL_SOTA.md` |
| Ablation numbers (F1 0.431 vs 0.888) | Manuscript body text (verify against `docs/CLAIM2_DATA_AND_TOOLS.md` or `docs/INDEL_LOSS_SKELETAL.md` for the underlying run before citing in a new context) |

## What was verified this session vs cited from earlier work

**Verified this session, directly, number by number**: T2.1, T2.3, T2.4,
T2.5's manuscript tables cross-checked against their CSVs — exact match,
zero discrepancies. `CAPS_PLOIDY` and `CAPS_CALL_INDELS` flag mechanics
confirmed by reading the actual code, not assumed from documentation.

**Not independently re-verified this session** (cited from prior,
documented work): the exact F1 0.431→0.888 ablation numbers' original run
log, and whether `build_substrate`/`collapse_contigs` are precisely the
functions that ablation toggled — strongly inferred from naming and commit
history, not confirmed by re-running the ablation this session. Flag this
distinction if the ablation numbers are cited in new analysis.
