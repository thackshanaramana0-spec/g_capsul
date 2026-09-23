---
Title: G_CAPSUL — Master Repository Map
Purpose: One file connecting every claim (COMPACT / FAITHFUL / ADDRESSABLE) to
  its figure, its implementing code, its result data, and its supporting
  documentation, using real, verified paths as they exist in this repository.
Scope: This file describes the PUBLIC repository (this `docs/` tree and its
  siblings at repo root). The manuscript itself, full development history,
  and internal research notes live outside this repo, in the project's
  private working tree, and are not published here.
Verification: Every path below was checked against the live repository at
  the time this file was written. Where a cited number disagrees with the
  CSV it names, the CSV wins.
---

# Master repository map

This repository has two layers. The **root** is the clean, reproducible
view: source code, raw result data, and curated documentation. Nothing here
is draft material or intermediate history.

```
g_capsul/
├── src/            the encoder and decoder, and the headers they share
├── include/         the same headers, duplicated at root for #include paths
├── scripts/         the minimal build / encode / verify closure
├── thirdparty/       vendored dependencies (never modified, never claimed as ours)
├── results/          raw CSVs and plots behind every number in the paper
├── docs/             this folder — the curated, cross-checked reference tree
└── Figures/          the four manuscript figures, Fig1.tif–Fig4.tif
```

Each of the three claims below is a row in three tables: what the claim
says, what proves it, and where every piece of that proof physically lives.

---

## Claim 1 — COMPACT (lossless compression)

**What it claims:** the archive can store a complete FASTQ file (sequence,
identifiers, quality, separator line) losslessly, competitively with SPRING
and Genozip, and losslessly against PgRC2 on the sequence stream alone.

| Piece | Path |
|---|---|
| Figure | `Figures/Fig1.tif` — architecture overview (encoder stages → archive → four read paths) |
| Figure | `Figures/Fig2.tif` — archive size as % of raw input, 19 datasets, G_CAPSUL vs SPRING vs Genozip |
| Encoder | `src/encoder.cpp` |
| Decoder | `src/decoder.cpp` |
| Sequence coder | `include/seqpar_core.h` — the DNA context-mixing coder shared by encoder and decoder |
| Stream selector / general coders | `include/coders_inproc.h` (selector, transforms), `include/coders_pgrc.h` (PPMd7 / FSE / range coder) |
| Identifiers | `include/names_coder.h` |
| Quality scores | `include/quality_coder.h` (wraps vendored `thirdparty/htscodecs/fqzcomp_qual.c`) |
| Build | `scripts/build106.sh` → binary `best106` (must include `-fopenmp`, see the script's own header comment) |
| Encode / verify | `scripts/encode_adaptive.sh`, `scripts/verify_lossless.sh` |
| Result data | `results/claim1/claim1_T1.1_T1.2.csv` — one row per (dataset, tool), archive size + compress/decompress time + peak RAM, all 19 datasets × 3 tools |
| Result plot source | `results/plots/claim1/make_claim1_bar_vertical.py` → `claim1_margin_vertical.{pdf,png}` (the `.tif` copy is `Figures/Fig2.tif`) |
| Docs — overview | `docs/claim1/overview_claim1.md` |
| Docs — per-table detail | `docs/claim1/t11_archive_size_claim1.md` (T1.1/T1.3), `docs/claim1/t12_wall_time_claim1.md` (T1.2) |
| Docs — names/quality detail | `docs/claim1/names_and_quality_claim1.md` |
| Docs — mechanism and PgRC2 comparison | `docs/claim1/mechanism_insight_claim1.md` |
| Docs — exact function/line references | `docs/claim1/code_mapping_claim1.md` |

**Note on `src/`:** the folder holds exactly two `.cpp` files, `encoder.cpp`
and `decoder.cpp`, plus their own `include/` copy of the 7 shared headers.
Older internal documentation (and the git history) refers to these same two
files by their pre-restructure names, `106_inprocess.cpp` and
`capsule_decode.cpp` — they are the identical files, renamed when the
repository was split into a public root and a private working tree. Every
reference to the old names inside `docs/` has been corrected to `src/encoder.cpp`
/ `src/decoder.cpp`.

---

## Claim 2 — FAITHFUL (reference-free variant calling)

**What it claims:** heterozygous SNVs, indels, multi-allelic sites, and
tetraploid genotypes can be called directly from the archive, without a
reference genome, by first reconciling the fragmented allelic structure a
size-minimizing compressor produces.

| Piece | Path |
|---|---|
| Figure | `Figures/Fig3.tif` — reconciliation (aggressive/mild contig collapse) → candidate generation (SNV/indel paths) |
| Caller | `include/caps_caller.h` — reconciliation, coverage-adaptive thresholds, SNV/indel/multi-allelic/tetraploid candidate generation and genotyping |
| Decoder entry point | `src/decoder.cpp` — `call` mode dispatches into `caps_caller.h` |
| Result data | `results/claim2/claim2_T2.1_snv.csv` (het-SNV), `claim2_T2.2_coverage_sweep.csv` (T2.2), `claim2_T2.3_indel.csv` (het-indel), `claim2_T2.4_multiallelic.csv` (T2.4), `claim2_T2.5_tetraploid.csv` (T2.5) |
| Result plot source | `results/plots/claim2/make_claim2_t21_bar.py`, `make_claim2_t22_line.py` → `claim2_t21_snv_f1.{pdf,png,tif}`, `claim2_t22_coverage_sweep.{pdf,png,tif}` |
| Docs — overview | `docs/claim2/overview_claim2.md` |
| Docs — per-table detail | `docs/claim2/t21_snv_claim2.md`, `t22_coverage_sweep_claim2.md`, `t23_indel_claim2.md`, `t24_multiallelic_claim2.md`, `t25_tetraploid_claim2.md` |
| Docs — mechanism (allele-splitting, the ablation) | `docs/claim2/mechanism_insight_claim2.md` |
| Docs — DiscoSNP++/Kmer2SNP architecture comparison | `docs/claim2/architecture_vs_discosnp_claim2.md` |
| Docs — exact function/line references | `docs/claim2/code_mapping_claim2.md` |

---

## Claim 3 — ADDRESSABLE (export, coverage, query, locus retrieval)

**What it claims:** the same archive supports pseudogenome export, per-base
coverage, coordinate-range query, exact-match retrieval (via an optional
completion index), and locus retrieval at heterozygous sites (via bilateral
upstream/downstream querying) — all without a second assembly or a
reference genome.

| Piece | Path |
|---|---|
| Figure | `Figures/Fig4.tif` — decoded representation → coverage / export / search (coordinate, sequence, completion-index) / heterozygous locus retrieval |
| Decoder — export / coverage / query | `src/decoder.cpp` — `export`, `coverage`, `query`, and `index` modes |
| Contig boundaries at encode time | `src/encoder.cpp` — `contig_spans` stream, gated by `CAPS_SPANS` |
| Result data | `results/claim3/claim3_T3.1_T3.2_T3.3.csv` (export/coverage/range-query timing, all 19 datasets), `results/claim3/claim3_T3.4_locus_fidelity.csv` (original single-individual locus-fidelity trail; superseded numbers now cited from `docs/claim3/t35_locus_retrieval_claim3.md`, see below) |
| Docs — overview | `docs/claim3/overview_claim3.md` |
| Docs — per-table detail | `docs/claim3/t31_export_claim3.md`, `t32_coverage_claim3.md`, `t33_range_query_claim3.md`, `t34_exact_match_claim3.md`, `t35_locus_retrieval_claim3.md` |
| Docs — mechanism, novelty vs PgRC2/BEETL | `docs/claim3/mechanism_insight_claim3.md` |
| Docs — exact function/line references | `docs/claim3/code_mapping_claim3.md` |

**Note on T3.4/T3.5 numbering:** the internal docs in this folder split
locus addressability into T3.4 (exact-match recall, compared against
BWT-family tools) and T3.5 (locus retrieval, both alleles). The manuscript's
own table numbering is independent of this folder's numbering — it runs
T3.1 through T3.5 in its own sequence, matching the same underlying five
operations. `docs/numbers_and_verification_index.md` is the canonical
source for the final, multi-individual T3.4/T3.5 numbers (1,452/1,452 with
the completion index, 1,451/1,452 native).

---

## Cross-claim documentation (applies to more than one claim)

| File | What it covers |
|---|---|
| `docs/central_insight_and_mechanism.md` | The one finding both Claim 2 and Claim 3 rest on: a size-minimizing compressor places the two alleles of a heterozygous site on disconnected pseudogenomic segments, which is why both calling and locus retrieval require reconciliation |
| `docs/technical_architecture.md` | Full architecture description independent of any single claim |
| `docs/paper_overview.md` | How the paper's sections map onto the three claims |
| `docs/numbers_and_verification_index.md` | The single canonical source for every final, locked number cited anywhere in this folder |
| `docs/novelty_and_prior_art.md` | What is and is not novel relative to PgRC2, PgRC, Quip, BEETL-fastq, CIndex, sFASTQ, BFQzip, DiscoSNP++, Kmer2SNP |
| `docs/honest_limitations_and_scope.md` | Every limitation and scope boundary that must appear in the paper's Discussion |
| `docs/coverage_manifest.md` | Audit trail proving this folder's completeness against the project's full internal documentation |
| `docs/extras/` | Background engineering history (PgRC2 reimplementation notes, layer-by-layer profiling, SOTA comparison tables, server setup, locked dataset list) — lower-priority context, not required to verify any single claim |

---

## Vendored dependencies (`thirdparty/`)

| Package | License | Used for |
|---|---|---|
| `thirdparty/ppmd/` | Public domain (LZMA SDK) | PPMd7 entropy coding, one of the general-purpose stream codecs `coders_pgrc.h` selects between |
| `thirdparty/fse/` | BSD (Yann Collet) | FSE / Huffman entropy coding, another selectable stream codec |
| `thirdparty/htscodecs/` | BSD 3-clause (Genome Research Ltd / James Bonfield) | `fqzcomp_qual.c` — the vendored quality-score codec wrapped by `include/quality_coder.h` |

None of the three is modified from upstream. Both `scripts/build106.sh` and
`scripts/build_decode.sh` compile the vendored C sources with `gcc` and the
project's own C++ with `g++`, then link them together (see either script's
header comment for why the two compilers cannot be merged into one
invocation).

---

## What is not in this repository, and where it actually lives

This is deliberate, not an oversight — the root of this repository is the
reproducible code+data+docs view, not the full development history.

- **The manuscript itself** (`capsul_manuscript.tex` and its siblings) is not
  published here.
- **Full development history** — `stages/` (the 96-file experimental
  progression), `server/` (machine setup), `benchmark/documentation/`
  (dated reproduction logs), and superseded internal docs — is not
  published here. Where a file in this folder's own history once cited one
  of those paths, the citation has been corrected to either point at the
  equivalent public path (for code and result CSVs) or reworded to say
  plainly that the source is internal and not published, rather than left
  pointing at a path that does not resolve in this repository.

---

## How to reproduce a number end to end

1. Pick a claim above and its result CSV.
2. Find the exact command that produced that CSV's row: `scripts/build106.sh`
   to build, `scripts/encode_adaptive.sh` to encode, then the relevant mode
   of the built decoder (`call`, `export`, `coverage`, `query`, `index`) to
   produce the downstream number.
3. Cross-check the row against the claim's own `docs/claimN/tXX_*.md` file,
   which states where that exact number is verified against the manuscript.
4. If a number in a `docs/` file and its result CSV ever disagree, the CSV
   is authoritative — this rule applies throughout this folder.
