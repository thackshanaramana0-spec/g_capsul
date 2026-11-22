# G_CAPSUL: Genomic, Compact, Addressable, Pseudogenome-Structured, Unified Lossless

[![CI](https://github.com/thackshanaramana0-spec/g_capsul/actions/workflows/ci.yml/badge.svg)](https://github.com/thackshanaramana0-spec/g_capsul/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Docker-lightgrey.svg)](#installation)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Lossless](https://img.shields.io/badge/lossless-byte--exact-brightgreen.svg)](#results)

---

## Overview

G_CAPSUL is a **lossless, reference-free FASTQ compressor** designed to keep the structure
discovered during compression useful after the archive is written. Following the **ARCH
design principle, Assemble, Retain, Compress, and Harness**, it builds a pseudogenome from
the reads and retains their placements, orientations, lengths, and local sequence
differences as a reusable representation. When compression separates allelic evidence
across pseudogenomic regions, G_CAPSUL reconciles those relationships during analysis so
reads can again be interpreted at their shared biological locus. The resulting archive
supports exact FASTQ reconstruction, reference-free variant discovery, coverage, and locus
retrieval, extending the compressed representation from a storage endpoint into a substrate
for downstream genomic analysis.

---

## The Problem

Lossless FASTQ compressors such as SPRING and Genozip focus on compact storage and exact
reconstruction, while the structure inferred during compression is not retained as a
general-purpose representation for downstream analysis. Reference-free variant callers such
as DiscoSNP++ and Kmer2SNP therefore rebuild related structure from the reads again. G_CAPSUL
targets this gap by retaining the compression-derived representation so that downstream
analysis can begin from the archive itself rather than rebuilding the underlying organization
from raw FASTQ, extending compression from a storage endpoint into a reusable substrate for
genomic analysis.

## Proposed Solution

G_CAPSUL follows the **ARCH design principle: Assemble, Retain, Compress, and Harness**.
During compression, it assembles a **pseudogenome** from the reads using greedy suffix-prefix
overlap chaining, pigeonhole mapping for the remaining reads, and self-matching to remove
residual redundancy. The resulting structure is retained inside the archive and reused for
three purposes:

1. **Compression (COMPACT)**: the pseudogenome, together with per-read placements, mismatches,
   and original read order, is entropy-coded into the smallest lossless archive among the
   three tools compared.
2. **Reference-free variant calling (FAITHFUL)**: the retained read organization is reused
   for variant analysis through aggressive collapse for SNV pileup and milder collapse for
   bubble and indel extraction, allowing heterozygous SNV, indel, and multi-allelic calling
   without constructing a second assembly.
3. **Archive-native addressability (ADDRESSABLE)**: the pseudogenome and retained placement
   information support `export`, `coverage`, and `query`, allowing sequence export, depth
   calculation, and coordinate-range read retrieval to begin from the archive rather than
   rebuilding the genomic organization from the original FASTQ.

## Key Features

- **Lossless reconstruction**: sequence, read order, identifiers, quality scores, and line-3 mode are preserved for supported input (see [Input requirements](#input-requirements)). Every reported archive was decoded and its records compared with the original input.
- **Smallest complete lossless archive** across 19 real datasets totaling 59.1 GB of FASTQ, with **19/19 wins against SPRING and 19/19 against Genozip**. Across all three tools, **57/57 archives reconstructed to byte-identical input**.
- **Reference-free variant calling from the archive** without returning to the original FASTQ, using a reference genome, or building a second assembly. Across four GIAB individuals, heterozygous-SNV mean F1 reached **0.876**, compared with 0.853 for DiscoSNP++ and 0.475 for Kmer2SNP.
- **Multi-allelic and higher-ploidy support**: G_CAPSUL recovered **17/26** multi-allelic sites in the complete chromosome-20 census, while DiscoSNP++ produced no complete multi-allelic records. On a tetraploid dataset formed from real HG003 and HG004 reads, G_CAPSUL achieved SNV F1 **0.897** versus 0.782.
- **Archive-level addressability**: pseudogenome export was **129–784×** faster than SPAdes, and per-base coverage was **16–54×** faster than BWA+mosdepth. The `query` operation also recovers both allelic read sets at heterozygous loci even when compression separates them across different pseudogenome coordinates, succeeding at **1,452 of 1,452** tested sites across four individuals at one chr20 locus. An FM-index built from the archive after compression returns every read containing a query sequence (exact-match recall **1.0000** on all four), without changing the archive.

---

## Contents

- [Repository layout](#repository-layout)
- [Results](#results)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Try it on your own FASTQ](#try-it-on-your-own-fastq)
- [Reproduce the paper tables](#reproduce-the-paper-tables)
- [Configuration](#configuration)
- [Design](#design)
- [Archive format](#archive-format)
- [Limitations](#limitations)
- [Known issues](#known-issues)
- [Citation](#citation)
- [License](#license)

---

## Repository layout

```
src/, include/    the shipped encoder and decoder: assembly, mapping, stream coding,
                  the FAITHFUL variant caller, and export/coverage/query
scripts/          build and test scripts, plus the benchmark harness behind every
                  number in results/
tests/            FM-index self-test and end-to-end test (run by ctest)
thirdparty/       vendored PPMd7, FSE/Huf0, htscodecs/fqzcomp, each under its own license
results/          the reported numbers: per-claim CSVs plus generated plots
docs/             reference documentation: claim1/, claim2/, claim3/, the mechanism,
                  novelty, limitations and configuration documents, and docs/extras/
Figures/          the manuscript figures
```

| you want | go to |
|---|---|
| **the master cross-reference: claim → figure → code → result → doc** | **[`docs/REPO_MAP.md`](docs/REPO_MAP.md)** |
| what's *not* claimed, explicitly | [`docs/limitations_and_scope.md`](docs/limitations_and_scope.md) |
| what's novel vs. prior art | [`docs/novelty_and_prior_art.md`](docs/novelty_and_prior_art.md) |
| the single canonical source for every locked number | [`docs/numbers_and_verification_index.md`](docs/numbers_and_verification_index.md) |

The manuscript itself and the full development history (experimental stages, machine
setup, dated reproduction logs) are not published in this repository. `docs/REPO_MAP.md`
explains that split in full. **Where a document and a result CSV disagree, the CSV wins.**

---

## Results

### Compression (COMPACT)

**Manuscript Table 1.** Whole-file lossless archive (sequence + names + quality + line 3) on the 19
locked datasets, 59.1 GB of FASTQ spanning bacteria, archaea, viruses, fungi, protists,
plants, animals and human. Every archive was decoded **before** being counted, and its
(name, sequence, quality) records were compared with the input's. All 57/57 archives across
the three tools passed. That benchmark comparison sorts the records, so it does not test read
order or line 3; ordered byte-exact reconstruction, line 3 included, is shown by the procedure
in [Try it on your own FASTQ](#try-it-on-your-own-fastq).

Percentage is how much smaller G\_CAPSUL's archive is than the competitor's, per row
(negative = smaller). All sizes in MB.

| kingdom | dataset | original | G_CAPSUL | SPRING | Genozip | vs SPRING | vs Genozip |
|---|---|---:|---:|---:|---:|---:|---:|
| Virus | ERR5181310 | 471.18 | **8.68** | 9.98 | 9.22 | −13.06% | −5.88% |
| | SRR32429602 | 1868.92 | **55.12** | 57.55 | 107.18 | −4.22% | −48.58% |
| Bacteria | SRR554369 | 456.44 | **57.29** | 59.20 | 88.81 | −3.23% | −35.50% |
| | ERR552797 | 431.08 | **46.97** | 52.13 | 82.80 | −9.90% | −43.27% |
| | SRR2584863 | 695.16 | **68.43** | 74.09 | 114.59 | −7.64% | −40.28% |
| | SRR40271341 | 290.89 | **40.55** | 45.78 | 62.26 | −11.43% | −34.87% |
| | ERR17740259 | 1336.96 | **83.20** | 92.40 | 161.81 | −9.96% | −48.58% |
| Archaea | SRR29296997 | 210.23 | **15.66** | 17.54 | 27.23 | −10.71% | −42.49% |
| | ERR12954017 | 333.14 | **15.16** | 17.00 | 39.43 | −10.79% | −61.54% |
| Protista | SRR37283774 | 1019.82 | **67.20** | 71.34 | 94.99 | −5.80% | −29.25% |
| | SRR36741279 | 1661.16 | **106.25** | 117.62 | 194.66 | −9.67% | −45.42% |
| Fungi | DRR976266 | 2248.04 | **56.54** | 61.76 | 201.14 | −8.45% | −71.89% |
| | SRR39257532 | 1761.85 | **108.27** | 153.86 | 227.02 | −29.63% | −52.31% |
| Plantae | SRR10676752 | 15372.50 | **1650.03** | 1718.55 | 2954.90 | −3.99% | −44.16% |
| Animalia | SRR065390 | 11282.99 | **887.41** | 942.64 | 1597.08 | −5.86% | −44.44% |
| | HG002 | 4279.20 | **573.77** | 598.19 | 951.00 | −4.08% | −39.67% |
| | HG003 | 4279.57 | **567.12** | 591.14 | 944.22 | −4.06% | −39.94% |
| | HG004 | 4279.33 | **606.19** | 634.75 | 1041.66 | −4.50% | −41.81% |
| | HG005 | 6809.02 | **997.52** | 1081.38 | 1694.51 | −7.76% | −41.13% |
| | **aggregate** | **59087.49** | **6011.37** | 6396.90 | 10594.52 | **−6.03%** | **−43.26%** |

**19/19 wins vs SPRING and 19/19 vs Genozip.**

Raw CSV: [`results/claim1/claim1_T1.1_T1.2.csv`](results/claim1/claim1_T1.1_T1.2.csv).
Full breakdown, including timing and peak memory: [`docs/claim1/`](docs/claim1/).

Separately, sequence-only content against PgRC2's own binary (the closest architectural
relative, GPL-3, run from its own source, never vendored): **+1.88% aggregate** on the
seven datasets both tools could process. That is a narrower measurement (PgRC2 stores no
names, no quality, and no line 3) and is deliberately not in the CSV above. The comparison
is limited to seven datasets, not the full non-human collection, because PgRC2 cannot
process the rest at all: attempting to run it on datasets with variable-length reads either
produces an explicit refusal or crashes, on six of the fourteen originally tested
non-human datasets. G_CAPSUL processes every one of these without modification, since
variable-length reads are an ordinary case throughout, not a special path. Full breakdown:
[`docs/claim1/mechanism_insight_claim1.md`](docs/claim1/mechanism_insight_claim1.md).

### Reference-free variant calling (FAITHFUL)

**Manuscript Tables 2 and 3, and the multi-allelic and tetraploid results.** Real GIAB HG002–HG005 chr20 at 30×, scored by third-party `rtg vcfeval`
against GIAB v4.2.1 inside the confident regions. **Every G_CAPSUL number here is called
from the compressed archive**: the archive is written with `CAPS_CALL=1`, and the calls come
from `CAPS_CALL_INDELS=1 capsule_decode call archive out.vcf workdir` (the full caller, as
run by `scripts/claim3_scripts/run_fullchr20_archive_capsule.sh`), with no FASTQ, no
reference and no separate assembly graph. Without `CAPS_CALL_INDELS=1`, `call` runs a faster
graph-only caller that emits SNVs only; the reported numbers do not come from that mode. Competitors get the identical reads, truth,
regions, normalisation and scorer, only the caller differs.

| comparison | G_CAPSUL | DiscoSNP++ | Kmer2SNP |
|---|---:|---:|---:|
| het-SNV F1 (mean of 4 individuals) | **0.876** | 0.853 | 0.475 |
| het-indel F1 (mean of 4) | **0.621** | 0.591 | not applicable (SNP-only by construction) |
| multi-allelic sites recovered | **17 / 26** | **0 / 26** | not applicable |
| tetraploid SNV F1 | **0.897** | 0.782 | not applicable |
| tetraploid indel F1 | **0.597** | 0.553 | not applicable |

Per individual, het-SNV: 0.888 / 0.891 / 0.891 / 0.834, **3 wins and 1 loss**, and the
loss is reported per individual. HG005 is the only variable-length dataset in the
set (250 bp quality-trimmed, 216 distinct read lengths, against 148 bp fixed elsewhere).
Truncating the same reads to 148 bp in a controlled diagnostic experiment raised its SNV F1
from 0.834 to 0.897, isolating read-length distribution as the cause. The published number
is the untruncated 0.834.
[`docs/claim2/t21_snv_claim2.md`](docs/claim2/t21_snv_claim2.md).

The multi-allelic row is a **complete census of chr20**, not a sample: 26 is every
SNV-only multi-allelic site there is. DiscoSNP++'s 0 is structural: across its entire
3,989-record output it emits no record with more than one ALT allele.

**Why reconciliation is necessary.** Lossless reconstruction alone does not guarantee that
the organization created for compression remains suitable for biological analysis. At a
heterozygous locus, reads supporting different alleles can be placed in separate,
internally-consistent pseudogenomic regions. Every read is still preserved exactly, but the
relationship identifying those reads as alternative forms of the same biological locus is no
longer explicit. G_CAPSUL therefore reconciles the separated regions before interpreting
their variant evidence. The full chr20 ablation isolates this effect:

| | F1 | Precision | Recall |
|---|---:|---:|---:|
| neither | 0.431 | 0.967 | 0.278 |
| re-placement only | 0.426 | 0.962 | 0.274 |
| collapse only | 0.648 | 0.963 | 0.488 |
| **both** | **0.888** | 0.956 | **0.830** |

The effect is synergistic, not additive: without collapse, precision remains above 0.96
while recall falls to about 0.27, consistent with variant evidence being separated rather
than incorrectly interpreted. Combining collapse with re-placement restores much of that
missing evidence, raising recall to 0.830 and F1 to 0.888. The correction operates during
analysis and does not alter the stored compressed representation, so it introduces **no
archive-size penalty**.

Raw CSVs: [`results/claim2/`](results/claim2/). Survey of why these two competitors
are the applicable ones, and the full mechanism behind the ablation:
[`docs/claim2/overview_claim2.md`](docs/claim2/overview_claim2.md),
[`docs/claim2/mechanism_insight_claim2.md`](docs/claim2/mechanism_insight_claim2.md).

### Archive-native addressability (ADDRESSABLE)

**Manuscript Table 4 (export, coverage) and the query operation.** G_CAPSUL reuses the representation already retained in the archive for sequence
export, per-base coverage, and locus retrieval rather than rebuilding the same organization
for each downstream task.

| operation | comparison | result |
|---|---|---:|
| export (pseudogenome as FASTA) | SPAdes v4.0.0 | **129–784× faster** |
| coverage (per-base depth) | bwa + samtools + mosdepth, timed as one pipeline | **16–54× faster** |
| query (reads at a locus) | no equivalent speed comparison | **locus-retrieval claim** |

**Read the export ratio with its correctness caveat.** Exported per-contig (using the
`contig_spans` stream, free at encode time), our export reaches **98.7% genome fraction**
against SPAdes's 98.3% and a **lower indel rate**, while SPAdes still wins on mismatch
rate, duplication ratio, and N50. These gaps trace to repeat resolution and
read error correction, which are the algorithmic content of a dedicated assembler and are not present
in a pseudogenome built to minimize compressed size. The speed win is unconditional. The
correctness comparison is mixed and is reported metric by metric. Coverage carries no such caveat: the depth numbers are
exact, not estimated.

**`query` is not a speed claim.** `genocat --head=100` extracts in 0.19 s
against our 0.46 s without an index, dropping to ≈0.03 s once the post-processing index step keeps
the decoded pseudogenome and placements. What `query` does that nothing else can is resolve
a **locus**, not just a coordinate or a matching string.

#### Locus retrieval

A biological locus is not necessarily one coordinate in a compression-derived pseudogenome.
At heterozygous sites, reads carrying different alleles can be placed on disconnected
pseudogenomic regions, so querying a single stored coordinate can recover only part of the
biological evidence. G_CAPSUL reconciles these separated regions and retrieves the read sets
belonging to the same locus:

| individual | locus | sites | both alleles recovered (archive + FM-index) |
|---|---|---:|---:|
| HG002 | chr20:3.0–3.6 Mb | 400 | 400/400 |
| HG003 | chr20:3.0–3.6 Mb | 335 | 335/335 |
| HG004 | chr20:3.0–3.6 Mb | 400 | 400/400 |
| HG005 | chr20:3.0–3.6 Mb | 317 | 317/317 |
| **all** | chr20:3.0–3.6 Mb | **1,452** | **1,452/1,452 (100%)** |

The 1,452 sites are four individuals at one locus. The count differs between individuals because a site is tested only if both 40-bp flanks are unambiguous and the upstream flank occurs in the pseudogenome, the same rule everywhere. A separate HG005 run at chr20:4.0–4.6 Mb (400/400) is recorded in [`docs/claim3/t35_locus_retrieval_claim3.md`](docs/claim3/t35_locus_retrieval_claim3.md) and is not part of this total. This table and the exact-match recall table are phase 4 of `scripts/benchmark_scripts/benchmark_1_run.sh`, which runs it by default (see [Reproduce the paper tables](#reproduce-the-paper-tables)). `PHASES=4` runs that phase alone, and needs `samtools`, `tabix`, `curl` and `python3`. One individual alone: `bash scripts/claim3_scripts/run_locus_retrieval_capsule.sh HG003 3000000 3600000 <workdir>`. Results land in `claim3_T3.4_exact_match_recall.csv` and `claim3_T3.5_locus_retrieval.csv`, and committed reference copies are in [`results/claim3/`](results/claim3/).

A heterozygous locus is **not one place** in a compression-optimal pseudogenome: reads
carrying the two alleles are routinely placed on different, disconnected pseudogenomic
segments, occasionally many megabases apart. `query` resolves this by probing each site
from both upstream and downstream and combining the two retrieved read sets. Position
alone leaves one HG004 site unresolved: its alternate-allele reads differ from the
pseudogenome consensus by an indel rather than a substitution, so they were placed
elsewhere at compression time, beyond a positional query's reach at any tolerance. The
**FM-index**, built after compression from the placements and deviations the archive
already retains, matches by read content instead of position and recovers that site, so
both alleles are found at all 1,452 sites. The cost is disclosed: tolerant locus search
also returns reads at homozygous negative-control sites (48, 52, 49 and 87 of 400 for
HG002 to HG005).

Same exact-match question against BWT-family tools (BEETL, CIndex). Searching the archive
gives recall 0.8536, 0.8457, 0.8099 and 0.6047 for HG002 to HG005, because a read that
differs from the pseudogenome is not a substring of it. Compression has already split the
reads in two: those that match the pseudogenome exactly, which the pseudogenome search
reaches, and those that deviate from it. The FM-index covers only the second group (14–37%
of reads here), so pseudogenome search plus the FM-index returns every read containing the
probe **by construction**: recall **1.0000** on all four individuals, at mismatch tolerance
0, the completeness a BWT index gives. It is 0.9–2.6 MB per individual, 3.5–3.8× smaller
than the k-mer completion index it replaces (`include/caps_fmi.h`; tests in `tests/`).

This is distinct from exact-sequence search: finding a matching sequence does not
necessarily identify other disconnected pseudogenomic regions representing the same
biological locus. The ADDRESSABLE claim is therefore about **biological locus retrieval
from the retained archive**, not query speed.

Full traceability: [`docs/claim3/`](docs/claim3/): `t31_export_claim3.md` through
`t35_locus_retrieval_claim3.md`, plus `mechanism_insight_claim3.md` for the full argument
and `code_mapping_claim3.md` for exact function and line references.

---

## Quick Start

```bash
git clone https://github.com/thackshanaramana0-spec/g_capsul.git
cd g_capsul
scripts/build_scripts/build106.sh /tmp/best106             # encoder (links -fopenmp)
scripts/build_scripts/build_decode.sh /tmp/capsule_decode  # decoder and the export/coverage/query/call/index operations
BEST=/tmp/best106 scripts/build_scripts/verify_lossless.sh path/to/reads.fq   # quick sequence round trip
```

`verify_lossless.sh` encodes the file, rebuilds the sequence column from the encoder's raw
stream dumps with `scripts/build_scripts/decode_105.py`, and compares it with the input's
sequence column. It checks the assembly, mapping and read order, not the archive's entropy
layer, names or quality. The full ordered, byte-for-byte check, through the archive itself, is
`BEST=/tmp/best106 DEC=/tmp/capsule_decode scripts/build_scripts/verify_roundtrip.sh path/to/reads.fq`
(prints `BYTE-IDENTICAL`, or `REFUSED` with the reason for unsupported input); the individual
steps are in [Try it on your own FASTQ](#try-it-on-your-own-fastq).

---

## Installation

**Requirements.** Linux (tested: Ubuntu 24.04 in CI, Ubuntu 26.04 under WSL2) or macOS with
Homebrew GCC (CI). A C++17 compiler with OpenMP (GCC; Apple Clang lacks `-fopenmp`),
`liblzma` headers, and `python3` for the scripts and tests. On Ubuntu/Debian:

```bash
sudo apt-get install build-essential cmake liblzma-dev python3
```

**Build with the scripts** (Quick Start above), or **with CMake**, which also builds the
FM-index tests:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"                 # build/best106, build/capsule_decode
(cd build && ctest --output-on-failure)          # fmi_selftest, fmi_end_to_end
```

Both builds default to a portable `-march=x86-64-v2` on x86-64 (`CAPS_BUILD_NATIVE=1` for the
scripts, `-DCAPS_BUILD_NATIVE=ON` for CMake, selects `-march=native`). Archives written by one
build decode with the other. `make` and `make test` are thin wrappers over the CMake commands
above.

**Docker**: build the image from source with `docker build -t g_capsul .` (the same two-stage
`Dockerfile` CI uses), or pull a published image:

```bash
docker pull ghcr.io/thackshanaramana0-spec/g_capsul:1.2.0    # or :latest
mkdir -p outdir
docker run --rm -v "$PWD":/data ghcr.io/thackshanaramana0-spec/g_capsul:1.2.0 \
    capsule_decode /data/out.capsule /data/outdir /data/outdir/reads.fq
```

The image ships `best106` (compress) and `capsule_decode` (decompress, plus its
`export`/`coverage`/`query`/`call`/`index` subcommands) and sets `DUMP_PERM=1 DUMP_MM=1`, which
a usable archive requires (see [Configuration](#configuration)). CI builds the image and
checks a round trip on every push ([`.github/workflows/ci.yml`](.github/workflows/ci.yml),
[`.github/workflows/publish-image.yml`](.github/workflows/publish-image.yml)).

The comparison tools used in the benchmarks (SPRING, Genozip, PgRC2, DiscoSNP++, Kmer2SNP,
SPAdes, bwa, samtools, mosdepth, rtg-tools) are installed separately; exact commands are in
[`docs/extras/SERVER_SETUP_AND_DOWNLOADS.md`](docs/extras/SERVER_SETUP_AND_DOWNLOADS.md).

---

## Usage

```bash
# Compress the FULL FASTQ (sequence + names + quality + line 3)
CAPS_NAMES=1 CAPS_QUAL=1 INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 \
    bash scripts/build_scripts/encode_adaptive.sh

# Compress sequence + read order only (the default: NOT a full FASTQ)
INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 bash scripts/build_scripts/encode_adaptive.sh

# Decompress. The output directory MUST exist (see Known issues).
mkdir -p outdir
/tmp/capsule_decode out.capsule outdir outdir/reads.fq
#   outdir/reads.fq        one sequence per line, original order
#   outdir/reads.fq.names  one header per line, without '@' (archives written with CAPS_NAMES=1)
#   outdir/reads.fq.qual   one quality line per read      (archives written with CAPS_QUAL=1)

# Addressable operations: each reads the same archive independently
/tmp/capsule_decode export   out.capsule contigs.fa
/tmp/capsule_decode coverage out.capsule coverage.tsv     # region, start, end, depth (run-length intervals)
/tmp/capsule_decode query    out.capsule region.fa 0-100000
/tmp/capsule_decode index    out.capsule out.qidx         # optional post-processing index for repeated queries

# Reference-free variant calling FROM the archive (archive written with CAPS_CALL=1)
CAPS_CALL=1 CAPS_NAMES=1 CAPS_QUAL=1 INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 \
    bash scripts/build_scripts/encode_adaptive.sh
CAPS_CALL_INDELS=1 /tmp/capsule_decode call out.capsule calls.vcf callwork   # SNVs and indels, the reported mode
/tmp/capsule_decode call out.capsule calls.vcf callwork                      # graph-only, SNVs only

# Variant calling fused with compression (single pass; CAPS_PLOIDY=4 for a tetraploid sample)
CAPS_CALL=1 CALL_VCF=calls.vcf /tmp/best106 reads.fq 3 16 16 22 16 16 1 24 64 1
```

Calls are in the coordinates of the archive's contigs. The benchmarks lift them to genome
coordinates with `scripts/benchmark_scripts/lift_vcf.py` (bwa alignment of the contigs) before
scoring with `rtg vcfeval`; see `scripts/claim3_scripts/run_fullchr20_archive_capsule.sh`.

Exact function and line references for every operation, organized by claim:
[`docs/claim1/code_mapping_claim1.md`](docs/claim1/code_mapping_claim1.md),
[`docs/claim2/code_mapping_claim2.md`](docs/claim2/code_mapping_claim2.md),
[`docs/claim3/code_mapping_claim3.md`](docs/claim3/code_mapping_claim3.md).

---

## Try it on your own FASTQ

This compresses a FASTQ through the archive, decodes it, rebuilds the 4-line FASTQ from the
decoder output alone, and compares it with the original byte for byte, read order and line 3
included. Run from the repository root after the [Quick Start](#quick-start) builds.

```bash
IN=reads.fq                 # your single-end FASTQ (uncompressed)
W=capsul_try; mkdir -p "$W/dec"

# 1. Optional input pre-check (the encoder and scope checker refuse the same cases, see
#    Input requirements). Each count must be 0.
grep -c $'\r' "$IN"                                                    # CRLF line endings
awk 'NR%4==2 && /[^ACGTN]/ {b++} END {print b+0}' "$IN"               # bases outside upper-case ACGTN
awk 'NR%4==2 && length($0)>256 && index(substr($0,257),"N") {b++} END {print b+0}' "$IN"   # N after base 256
awk 'NR%4==1 {h=substr($0,2)} NR%4==3 {if ($0=="+") p++; else if ($0=="+"h) c++; else o++}
     END {print (o>0 || (p>0 && c>0)) ? 1 : 0}' "$IN"                 # mixed or other line 3

# 2. Compress (full FASTQ). Unsupported input stops here with an error and exit code 1 or 2.
CAPS_NAMES=1 CAPS_QUAL=1 INPUT="$IN" ARCHIVE="$W/reads.capsule" \
    BEST=/tmp/best106 bash scripts/build_scripts/encode_adaptive.sh

# 3. Decompress
/tmp/capsule_decode "$W/reads.capsule" "$W/dec" "$W/dec/reads.fq"

# 4. Rebuild the FASTQ. Line 3 is stored as a mode: a bare '+' for every read, or '+'
#    followed by the read's header. Use the first command for '+', the second for '+header'.
paste -d'\n' <(sed 's/^/@/' "$W/dec/reads.fq.names") "$W/dec/reads.fq" \
             <(sed 's/.*/+/' "$W/dec/reads.fq") "$W/dec/reads.fq.qual" > "$W/rebuilt.fq"
# paste -d'\n' <(sed 's/^/@/' "$W/dec/reads.fq.names") "$W/dec/reads.fq" \
#              <(sed 's/^/+/' "$W/dec/reads.fq.names") "$W/dec/reads.fq.qual" > "$W/rebuilt.fq"

# 5. Compare
cmp "$IN" "$W/rebuilt.fq" && echo "BYTE-IDENTICAL"
ls -l "$IN" "$W/reads.capsule"
```

Expected: the four counts in step 1 print `0`, step 2 prints `ARCHIVE_TOTAL=<bytes>`, and
step 5 prints `BYTE-IDENTICAL`. In the fresh-environment test of this repository this was
obtained on simulated diploid reads (24,000 × 150 bp), variable-length reads with N
(60–150 bp, `+header` line 3), 200–600 bp reads, 1–30 bp reads, heavy duplication, all-N
reads and a single-read file. Encoder peak memory on those inputs was about 0.6 GB; the
largest benchmark dataset (15.4 GB of FASTQ) peaked at 11.4 GB.

### Input requirements

The archive is lossless for single-end FASTQ in which:

- sequences use upper-case `A`, `C`, `G`, `T`, `N` only;
- lines end in LF, not CRLF;
- every record has four lines, and each quality line has the length of its sequence;
- line 3 is a bare `+` in every record, or `+` followed by that record's header in every record;
- in a read longer than 256 bp, no `N` occurs after the first 256 bases, and no read has more
  than 255 `N` (an N's offset is stored in one byte);
- reads are at most 1023 bp.

The encoder checks all of this before it writes anything and refuses a file that breaks any
condition, printing the first offending record and exiting with a nonzero status (2 from the
encoder, 1 from `encode_adaptive.sh`, which prints the encoder's `FATAL` line). The refused
cases are: reads over 1023 bp, quality lines of a different length from the sequence,
truncated or malformed records, CRLF line endings, any sequence byte other than `A`, `C`, `G`,
`T`, `N` (lower case and IUPAC codes included), an `N` that violates the offset rule above,
mixed line-3 content (checked when names are stored, `CAPS_NAMES=1`), and files with no valid
record. `scripts/benchmark_scripts/check_input_scope.py` reports the same verdicts without
encoding. Convert CRLF with `tr -d '\r'`. Paired-end data must be supplied as separate
single-end files; mate relationships are not kept. Line 3 is stored as a mode, not per read,
which is why mixed content is refused. Without `CAPS_NAMES=1 CAPS_QUAL=1` the archive holds
sequences only and is not a full FASTQ.

---

## Reproduce the paper tables

Two scripts drive every manuscript table. Both expect the benchmark datasets, the comparison
tools and the GIAB truth files described in
[`docs/extras/SERVER_SETUP_AND_DOWNLOADS.md`](docs/extras/SERVER_SETUP_AND_DOWNLOADS.md) and
[`docs/extras/NEW_DATASET_LOCKED.md`](docs/extras/NEW_DATASET_LOCKED.md) (accessions, sizes,
and the chr20 GIAB HG002–HG005 inputs).

```bash
bash scripts/benchmark_scripts/benchmark_0_preflight.sh   # tools, data, network, binaries, FM-index self-test -> GO / NO-GO
bash scripts/benchmark_scripts/benchmark_1_run.sh         # all phases
PHASES=1,3 bash scripts/benchmark_scripts/benchmark_1_run.sh   # selected phases only
```

| manuscript | benchmark_1 phase | result file |
|---|---|---|
| Table 1, archive size | 1 | `results/claim1/claim1_T1.1_T1.2.csv` |
| Table 2, heterozygous SNV | 2 | `results/claim2/claim2_T2.1_snv.csv` |
| Table 3, heterozygous indel | 2 | `results/claim2/claim2_T2.3_indel.csv` |
| text, multi-allelic (17/26) | 2 | `results/claim2/claim2_T2.4_multiallelic.csv` |
| text, tetraploid | 2 | `results/claim2/claim2_T2.5_tetraploid.csv` |
| Table 4, export and coverage | 3 | `results/claim3/claim3_T3.1_T3.2_T3.3.csv` |
| Table 5, exact-match recall | 4 | `results/claim3/claim3_T3.4_exact_match_recall.csv` |
| Table 6, locus retrieval | 4 | `results/claim3/claim3_T3.5_locus_retrieval.csv` |

**Resources.** The reported runs used a 12-core Linux server. Phase 1 compresses 59.1 GB of
FASTQ (G_CAPSUL alone: about 47 minutes of compression in total, at most 16 minutes and
11.4 GB peak memory for one dataset), and runs SPRING and Genozip on the same files. Phase 2
calls variants on chr20 for four individuals (G_CAPSUL peak 11.4–19.9 GB per individual) and
runs DiscoSNP++ and Kmer2SNP. Plan for at least 32 GB of RAM and several hundred GB of disk
for phases 1–3.

**Phase 4 one individual at a time.** Phase 4 streams GIAB alignments for a 600 kb chr20
window and needs `samtools`, `tabix`, `curl` and `python3`. Run as one job on a machine with
about 7 GB of memory, the combined phase was stopped for lack of memory. One individual per
run completes on that machine, and is how the committed Table 5 and 6 rows were produced:

```bash
for IND in HG002 HG003 HG004 HG005; do
    bash scripts/claim3_scripts/run_locus_retrieval_capsule.sh "$IND" 3000000 3600000 "work_$IND"
done
```

**What is measured where.** Every number in `results/` comes from these scripts; where a
document and a CSV disagree, the CSV is authoritative. Each CSV row records its status and
source. The benchmark lossless check (`losscmp` in `benchmark_1_run.sh`) compares sorted
(name, sequence, quality) records, as described under [Results](#results).

---

## Configuration

Exactly what each configuration puts in the archive and gives back on decode, including the
important point that the default is sequence-only, not a FASTQ:

| you set | what's IN the archive | what decoding gives you |
|---|---|---|
| *(nothing, the default)* | sequence + read order only | sequence only: no names, no quality, no `+` line |
| `CAPS_NAMES=1` | + names/read-ID column | `<outreads>.names` in addition |
| `CAPS_QUAL=1` | + quality scores | `<outreads>.qual` in addition |
| `CAPS_NAMES=1 CAPS_QUAL=1` | full FASTQ content | sequence, names and quality files, from which the 4-line FASTQ is rebuilt byte for byte ([Try it on your own FASTQ](#try-it-on-your-own-fastq)) |
| `CAPS_CALL=1` | + `contig_spans` stream (as `CAPS_SPANS=1`), and writes a VCF as a side effect | `$CALL_VCF` gets heterozygous SNV/indel calls, combinable with any of the above |

`encode_adaptive.sh` and the Docker image also set `DUMP_PERM=1 DUMP_MM=1`. When `best106` is
run directly these two must be set, or the archive lacks the read-order and mismatch streams and
decodes to no reads. Every other switch and every positional argument is listed, with its
default, in [`docs/configuration_reference.md`](docs/configuration_reference.md).

---

## Design

![G_CAPSUL architecture: FASTQ in, one encoder assembling and compressing a pseudogenome archive, five read paths served from that same archive](Figures/Fig1.png)

G_CAPSUL follows the **ARCH design principle: Assemble, Retain, Compress, and Harness**.
Reads are organized once into a reference-free pseudogenome. Their placements, orientations,
lengths, and local sequence differences are retained in the archive rather than discarded
after compression. The same representation is then reused for exact reconstruction,
reference-free variant discovery, coverage, and locus-aware retrieval. When compression
separates biologically related allelic structures, G_CAPSUL reconciles them during
downstream analysis rather than constructing a second assembly.

### Algorithmic contributions

- **Multi-region pseudogenome construction.** Exact suffix-prefix chaining builds the primary
  pseudogenomic region. Reads that cannot be incorporated are placed against existing
  sequence when possible, while the remaining reads form a secondary pseudogenomic region.
  Repeated sequence within the resulting representation is subsequently reduced through
  maximal exact matching.
- **Heterozygous structure reconciliation and read re-placement.** Related pseudogenomic
  regions are identified from sequence overlap and local duplication, then brought into a
  shared analytical frame. Reads are placed again against the reconciled structure so
  evidence originally separated across pseudogenomic regions can contribute to the same
  biological locus.
- **Dual-substrate variant calling.** Two reconciled representations are maintained at
  different collapse levels. The more aggressive representation supports SNV pileup, while
  the milder representation preserves alternative sequence structure needed for indel
  discovery.
- **Multi-channel indel discovery.** Bubble traversal, gapped sequence comparison,
  cross-contig analysis, and an eBWT2SNP-inspired positional-clustering channel expose indel
  candidates that no single local representation captures reliably.
- **Placement-derived coverage.** Per-base coverage is calculated directly from retained read
  positions and lengths without decoding individual read sequences or performing a new
  alignment.
- **Biological locus retrieval across pseudogenomic addresses.** Queries can recover reads
  belonging to the same biological locus even when alternative alleles were placed at
  different pseudogenomic coordinates during compression. Bilateral flanking queries and
  retained placement information reconnect these separated representations.

---

## Archive format

A `.capsule` file is a small fixed header followed by named streams (all integers
little-endian):

```
"CAPSULE\0"  8 B    magic
version      u16    2
pg_len       u64    pseudogenome length in bases
main_end     u64    boundary between the main region and the second region
minmem       u32    minimum MEM length
n_streams    u16
n_streams x  [name_len u8][name][coded_len u64][payload]
payload      [method u8][raw_len u64][coded bytes]   (method selected per stream)
```

The decoder finds streams by name, so an added stream does not break older archives. Up to 29
streams hold the pseudogenome (literal bases and MEM references), read placements and strands,
read lengths, original read order, mismatches, N positions, and, when enabled, names
(`names_*`), quality and `contig_spans`. Each stream is documented in the table at the top of
[`src/encoder.cpp`](src/encoder.cpp), and the layout and pipeline in
[`docs/technical_architecture.md`](docs/technical_architecture.md). The optional `index` output
and the FM-index are separate files beside the archive and are not part of the reported archive
size.

**Tests.** `ctest` in the CMake build runs the FM-index self-test against brute force and an
end-to-end test on synthetic archives (variable lengths, contained reads, duplicates, N,
repeats, high error). CI additionally runs the lossless round trip and the regression tests
`scripts/claim2_scripts/test_claim2.sh` and `scripts/claim3_scripts/test_claim3.sh`.
[`docs/REPO_MAP.md`](docs/REPO_MAP.md) links each claim to its benchmark, result files,
implementation, figures and documentation.

---

## Limitations

- **Current scope is short-read sequencing.** The evaluation uses unpaired read sets and
  does not retain paired-end mate relationships. Long reads, large structural variants,
  population-scale inference, and general naturally polyploid analysis remain outside the
  current claims.
- **Downstream biological validation is not yet whole-genome.** Variant calling and locus
  retrieval are evaluated on chromosome 20 from four GIAB individuals at approximately 30×
  coverage. Performance should therefore not be extrapolated to complete human genomes or
  to non-human variant calling.
- **The current implementation trades computation for reusable structure.** Pseudogenome
  construction performs substantially more structural work than conventional FASTQ
  compression, making compression slower and, on larger datasets, more memory-intensive.
  Several construction stages are still serial or only partly parallelized.
- **Whole-genome calling is limited by the present implementation.** The current k-mer spill
  representation creates large intermediate-storage requirements at larger scales, making
  whole-genome variant calling impractical in its present form.
- **Variant calling remains sensitive to read geometry and candidate generation.** HG005
  shows that heterogeneous read lengths can reduce calling accuracy, while some indels
  remain inaccessible when the required alternative sequence structure is not exposed
  during candidate generation.
- **Exact-match and locus retrieval use one optional post-processing step.** It builds an
  FM-index over the reads the pseudogenome cannot represent and keeps the decoded
  pseudogenome for repeated queries. Together these are larger than the archive, are not
  part of the reported archive size, and grow with the share of reads that deviate from
  the pseudogenome. Tolerant locus search also returns reads at homozygous
  negative-control sites.
- **The input format is narrow and strictly enforced.** The sequence coder assumes upper-case
  ACGTN, LF line endings and a single line-3 convention, so files outside that are refused
  rather than converted. See [Input requirements](#input-requirements).

Full, current status for each claim: [`docs/claim1/overview_claim1.md`](docs/claim1/overview_claim1.md),
[`docs/claim2/overview_claim2.md`](docs/claim2/overview_claim2.md),
[`docs/claim3/overview_claim3.md`](docs/claim3/overview_claim3.md), and the single
cross-claim source for every scope boundary that must be disclosed:
[`docs/limitations_and_scope.md`](docs/limitations_and_scope.md).

---

## Known issues

- **Unsupported input is refused, not converted.** CRLF line endings, lower-case or IUPAC
  bases, an `N` after base 256 of a long read and mixed line-3 content stop the encoder with
  an error (see [Input requirements](#input-requirements)).
- **Reported archives were checked with a sorted comparison.** The benchmark's lossless check
  (`losscmp` in `benchmark_1_run.sh`) compares sorted (name, sequence, quality) records, so it
  confirms record content but not read order or line 3 for the 19 published datasets. The
  ordered byte-for-byte check, `scripts/build_scripts/verify_roundtrip.sh`, was run on the
  test inputs listed in [Try it on your own FASTQ](#try-it-on-your-own-fastq), not on the 19
  datasets.
- **Decoder output directory.** `capsule_decode <archive> <outdir> [reads]` requires `<outdir>` to
  exist and exits with status 2 if it does not. The `export`, `coverage`, `query`, `index` and
  `call` subcommands write a single file whose parent directory must exist.

---

## Citation

If you use G_CAPSUL in your research, please cite:

> Thackshanaramana B (2026). *G_CAPSUL: a unified pseudogenome for lossless FASTQ
> compression, reference-free variant calling, and archive-native addressability.*
> Manuscript in preparation.

To cite the exact code version used, see [`CITATION.cff`](CITATION.cff) (GitHub's own "Cite
this repository" button reads this automatically) or reference a tagged
[release](https://github.com/thackshanaramana0-spec/g_capsul/releases) directly. Current
release: [v1.2.0](https://github.com/thackshanaramana0-spec/g_capsul/releases/tag/v1.2.0).

**Author:** Thackshanaramana B, SRM Institute of Science and Technology.

---

## License

**Apache License 2.0**: see [`LICENSE`](LICENSE) for the full text.

Vendored third-party code keeps its own license, all of them Apache-2.0-compatible.
Full table and terms in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
