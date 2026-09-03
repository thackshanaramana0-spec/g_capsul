# CAPSULE: Compact, Addressable, Pseudogenome-Structured, Unified Lossless Compressor

[![CI](https://github.com/thackshanaramana0-spec/capsule/actions/workflows/ci.yml/badge.svg)](https://github.com/thackshanaramana0-spec/capsule/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-undecided-lightgrey.svg)](docs/INDUSTRIAL_CHECKLIST_OVERALL.md)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](#build)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Lossless](https://img.shields.io/badge/lossless-byte--exact-brightgreen.svg)](#results)

> Compress Illumina short reads to a smaller-than-SPRING-and-Genozip lossless archive, call heterozygous SNVs, indels and multi-allelic sites with no reference genome, and export the assembly, per-base depth, or a coordinate range straight from the archive — one pseudogenome, three capabilities.

---

## The Problem

Lossless FASTQ compressors (SPRING, Genozip) treat compression as a storage-only problem and
discard the internal assembly structure they build to get there. Reference-free variant
callers (DiscoSNP++, Kmer2SNP) reprocess the same reads from scratch to build their own
assembly. Neither exposes what a de-novo assembly already contains once built: the sample's
own approximate genome, and where every read sits on it. No tool this project found combines
smallest-lossless-archive, reference-free variant calling, and direct archive-level
export/coverage/coordinate-query in one pass over the reads.

## Proposed Solution

CAPSULE builds one structure — a **pseudogenome** assembled directly from the reads at
compress time by greedy suffix-prefix overlap chaining, pigeonhole mapping of the
remainder, and self-matching to remove residual redundancy — and reuses it for three
purposes instead of discarding it after compression:

1. **COMPACT** — the pseudogenome plus per-read placement, mismatches, and read order are
   entropy-coded into the smallest lossless archive of the three tools compared.
2. **FAITHFUL** — the same read placements support a dual-substrate variant caller
   (aggressive collapse for SNV pileup, mild collapse for bubble/indel extraction), so
   heterozygous SNV, indel, and multi-allelic calling are a side effect of compressing,
   not a second pass.
3. **ADDRESSABLE** — `export`, `coverage`, and `query` read the pseudogenome and placement
   index directly from the archive: an assembly, a depth profile, or a coordinate-range
   read set, without decompressing the file or running an aligner.

## Key Features

- **Smallest lossless archive** on 14 real datasets: beats SPRING by 11.68% and Genozip by 48.93% aggregate, 14/14 wins against both
- **Reference-free variant calling** as a side effect of compression: het-SNV F1 0.890, het-indel F1 0.666, both beating the strongest applicable competitor
- **Multi-allelic and tetraploid calling**: native multi-allelic VCF output (a capability DiscoSNP++ structurally lacks), and real-data tetraploid SNV/indel wins built from real GIAB samples, not synthetic data
- **Archive-native addressability**: pseudogenome export 254–656× faster than SPAdes, per-base coverage 23–33× faster than bwa+mosdepth, coordinate-range query at 132× output selectivity
- **Byte-exact lossless**: sequence, read order, names, quality and line-3 mode all verified to reconstruct the original file, same MD5

---

## Contents

- [Results](#results)
- [Quick test](#quick-test)
- [Build](#build)
- [Usage](#usage)
- [Design](#design)
- [Repository layout](#repository-layout)
- [Reproducing benchmarks](#reproducing-benchmarks)
- [What's open, honestly](#whats-open-honestly)
- [Citation](#citation)
- [Author](#author)
- [License](#license)

---

## Results

### Compression (COMPACT)

**Table 1.** Whole-file lossless archive (sequence + names + quality + line 3) on 14 real
public datasets spanning bacteria, viruses, fungi, protists, and one human virus, every
archive decoded back to byte-identical input before being counted.

| dataset | organism | CAPSULE | SPRING | Genozip |
|---|---|---:|---:|---:|
| SRR2584863 | E. coli B REL606 | 68,677,977 | 74,045,440 | 119,618,198 |
| ERR552797 | M. tuberculosis H37Rv | 46,964,185 | 52,101,120 | 82,799,524 |
| SRR554369 | P. aeruginosa PAO1 | 57,320,645 | 59,166,720 | 88,813,037 |
| ERR5181310 | SARS-CoV-2 | 8,686,774 | 9,390,080 | 9,218,949 |
| ERR17740259 | S. aureus | 83,421,422 | 92,303,360 | 164,126,889 |
| DRR976266 | S. cerevisiae | 56,684,465 | 61,716,480 | 201,140,420 |
| SRR36741279 | L. major | 106,342,010 | 117,565,440 | 194,786,333 |
| SRR37283774 | P. falciparum | 67,382,606 | 71,168,000 | 94,989,854 |
| SRR32429602 | HCMV | 55,694,667 | 57,344,000 | 106,983,461 |
| SRR39257532 | A. fumigatus | 108,558,425 | 153,589,760 | 227,157,899 |
| SRR29296997 | H. salinarum | 15,654,489 | 17,500,160 | 27,233,618 |
| ERR12954017 | S. acidocaldarius | 15,159,145 | 16,967,680 | 39,545,150 |
| SRR40271341 | H. pylori | 40,620,675 | 45,742,080 | 62,259,448 |
| SRR40402583 | C. jejuni | 9,040,464 | 9,543,680 | 30,863,017 |

**14/14 wins vs SPRING (+11.68% aggregate). 14/14 wins vs Genozip (+48.93% aggregate).**
Raw CSV: [`results/phase_a/allphases_14dataset.csv`](results/phase_a/allphases_14dataset.csv).
One locked dataset (Utricularia gibba, `SRR10676752`) is not yet run — see
[What's open, honestly](#whats-open-honestly).

Separately, sequence-only content against PgRC2's own binary (the closest architectural
relative, GPL-3, run from its own source): **+1.88% aggregate, 6 wins, 1 loss** (S.
acidocaldarius, −0.83%) across 7 datasets. Full breakdown:
[`docs/CLAIM1_FINAL_VERDICT.md`](docs/CLAIM1_FINAL_VERDICT.md).

### Reference-free variant calling (FAITHFUL)

**Table 2.** Real GIAB HG002–HG005 chr20 data, 8 independent chr20-window evaluations plus
a real tetraploid construction, scored by third-party `rtg vcfeval`.

| comparison | CAPSULE | DiscoSNP++ | Kmer2SNP |
|---|---:|---:|---:|
| het-SNV F1 | **0.890** | 0.874 | 0.464 |
| het-indel F1 | **0.666** | 0.639 | not applicable (SNP-only by construction) |
| multi-allelic sites recovered | **11/18** | 0/18 | not applicable |
| tetraploid SNV F1 | **0.836** | 0.782 | not applicable |
| tetraploid indel F1 | **0.567** | 0.553 | not applicable |

DiscoSNP++ and Kmer2SNP are the only two general-purpose reference-free callers found
applicable to this exact task (single sample, no reference; literature survey in
[`docs/HET_INDEL_SOTA.md`](docs/HET_INDEL_SOTA.md)). The tetraploid rows are built from
real HG003+HG004 reads mixed following Cooke, Wedge & Lunter, *Genome Research* 2022 —
nothing simulated. Full tables and the two measurement bugs whose fixes flipped het-indel
and tetraploid-indel from documented losses to wins:
[`docs/CLAIM2_TABLES_AND_INDEL_SCAN.md`](docs/CLAIM2_TABLES_AND_INDEL_SCAN.md).

### Archive-native addressability (ADDRESSABLE)

**Table 3.** `export`/`coverage`/`query` served directly from the archive vs the
conventional pipeline that would otherwise compute the same thing, on real data.

| operation | vs | speedup |
|---|---|---:|
| export (pseudogenome as FASTA) | SPAdes v4.0.0 | **555–656×** |
| coverage (per-base depth) | bwa + samtools + mosdepth | **23–33×** |
| query (coordinate-range reads) | full decompression | 1.62× time, **132× fewer reads / 112× fewer bytes** returned |

Full numbers, exact commands, and the two real bugs found and fixed while verifying them:
[`docs/CLAIM3_LOCKED.md`](docs/CLAIM3_LOCKED.md).

---

## Quick test

```bash
bash scripts/run_capsule.sh 1   # COMPACT    — verify losslessness on the locked E. coli dataset
bash scripts/run_capsule.sh 2   # FAITHFUL   — fast synthetic caller regression test, no downloads needed
bash scripts/run_capsule.sh 3   # ADDRESSABLE — fast synthetic decoder regression test, no downloads needed
```

Each claim is independent — there is no requirement to run them in order. Every command
builds whatever binaries it needs on first use and prints which dataset it used. No dataset
path is hardcoded: [`scripts/capsule_config.sh`](scripts/capsule_config.sh) is the single
file to edit if your data moves, or override for one run:

```bash
CAPSULE_DATA_DIR=/mnt/other/fastq bash scripts/run_capsule.sh 1
```

---

## Build

```bash
git clone https://github.com/thackshanaramana0-spec/capsule.git
cd capsule
scripts/build106.sh /tmp/best106            # encoder (must link -fopenmp, see the script's own note)
scripts/build_decode.sh /tmp/capsule_decode # decoder + all three Claim 3 operations
```

**Dependencies:** C++17 compiler, `liblzma-dev`. PPMd7, FSE/Huf0, and htscodecs/fqzcomp are
vendored under `thirdparty/`, each with its own license included — nothing else to install
for the core binaries. Benchmark/comparison tools (SPRING, Genozip, PgRC2, DiscoSNP++,
Kmer2SNP, MEGAHIT, SPAdes, bwa, samtools, mosdepth, rtg-tools) are separate, with exact,
verified install commands in
[`docs/SERVER_SETUP_AND_DOWNLOADS.md`](docs/SERVER_SETUP_AND_DOWNLOADS.md).

```bash
# Ubuntu / Debian
sudo apt-get install build-essential liblzma-dev
```

---

## Usage

```bash
# Compress (sequence + read order only — the default, NOT a full FASTQ)
INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 bash scripts/encode_adaptive.sh

# Compress the FULL FASTQ (sequence + names + quality + line 3)
CAPS_NAMES=1 CAPS_QUAL=1 INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 \
    bash scripts/encode_adaptive.sh

# Decompress (byte-exact round trip)
/tmp/capsule_decode out.capsule outdir outdir/reads.fq

# Reference-free variant calling, fused with compression (single pass)
CAPS_CALL=1 CALL_VCF=calls.vcf /tmp/best106 reads.fq 3 16 16 22 16 16 1 24 64 1

# Multi-allelic / polyploid calling
CAPS_CALL=1 CAPS_PLOIDY=4 CALL_VCF=calls.vcf /tmp/best106 reads.fq 3 16 16 22 16 16 1 24 64 1

# Addressable archive operations — each reads the SAME archive independently
/tmp/capsule_decode export   out.capsule contigs.fa
/tmp/capsule_decode coverage out.capsule coverage.tsv
/tmp/capsule_decode query    out.capsule region.fa 0-100000
```

Exactly what each configuration puts in the archive and gives back on decode — including
the important point that the default is sequence-only, not a FASTQ — is spelled out in
full in [What you get, by configuration](#what-you-get-by-configuration) below. Every
command CAPSULE supports, organized by claim: [`docs/COMMANDS_REFERENCE.md`](docs/COMMANDS_REFERENCE.md).

### What you get, by configuration

| you set | what's IN the archive | what decoding gives you |
|---|---|---|
| *(nothing — the default)* | sequence + read order only | sequence only — no names, no quality, no `+` line |
| `CAPS_NAMES=1` | + names/read-ID column | `<outreads>.names` in addition |
| `CAPS_QUAL=1` | + quality scores | `<outreads>.qual` in addition |
| `CAPS_NAMES=1 CAPS_QUAL=1` | full FASTQ content | a full 4-line FASTQ, verified byte-identical (same MD5) to the original |
| `CAPS_CALL=1` | no archive change — writes a VCF as a side effect | `$CALL_VCF` gets heterozygous SNV/indel calls, combinable with any of the above |

---

## Design

### Algorithmic contributions

- **Multi-region pseudogenome assembly.** Greedy exact suffix-prefix chaining builds a main
  region from well-tiling reads; leftovers are pigeonhole-mapped or assembled into a second
  region; both regions are self-matched to remove residual redundancy.
- **Dual-substrate variant calling.** The same contig set is rebuilt at two collapse
  aggressiveness levels — one tuned for reliable minor-allele-fraction estimation in SNV
  pileup, one left closer to its pre-collapse form to preserve the two-path bubble structure
  indel calling needs.
- **Positional-clustering indel channel.** An eBWT2SNP-inspired channel clusters reads by a
  right-context anchor unique across the contig set, catching indels neither pileup nor
  bubble extraction reaches alone.
- **Hoisted, index-only coverage.** `coverage` needs only the pseudogenome's length and every
  read's placement, not its content — so it runs entirely before the pseudogenome is
  rebuilt, skipping the cost export and query both pay.
- **Adversarial correctness discipline.** Four silent data-loss bugs (mismatch positions
  above 256bp, orphaned unique-read desync, reverse-complement contained-read indexing, an
  FSE-RLE decode defect) were found by actually decoding archives and diffing against the
  original file, not by trusting a passing size table — documented in
  [`CLAUDE.md`](CLAUDE.md) §6.3.

---

## Repository layout

```
stages/106_inprocess.cpp     the shipped encoder — assembly, mapping, stream coding,
                              and (gated on env vars) names/quality/calling
stages/capsule_decode.cpp    the shipped decoder — full round trip plus export/coverage/query
stages/01...105               the full experimental progression, one file per decision
include/caps_caller.h        the Claim 2 variant caller
include/*_coder.h            stream-specific coders (names, quality, sequence, generic)
scripts/run_capsule.sh       single entry point for all three claims
scripts/capsule_config.sh    the one file to edit if dataset paths move
scripts/test_claim2.sh       synthetic caller regression test
scripts/test_claim3.sh       synthetic decoder regression test
scripts/run_claim3.sh        one-command real export/coverage/query benchmark
thirdparty/                  PPMd7 (public domain), FSE/Huf0 (BSD), htscodecs/fqzcomp (BSD)
docs/                        architecture, per-claim results, checklists, refuted ideas
                              — start at docs/INDEX.md, 68 files, navigable by topic
results/                     raw measurement CSVs, including reverted work
DATASET_LOCKED.md            the locked accessions — do not substitute without documenting why
```

---

## Reproducing benchmarks

Full datasets are not stored in this repository. Exact, verified download commands for
every dataset and every comparison tool — including the S3-mirror trick for large SRA
accessions and the chr20-only streaming method for GIAB BAMs — are in
[`docs/SERVER_SETUP_AND_DOWNLOADS.md`](docs/SERVER_SETUP_AND_DOWNLOADS.md).

```bash
bash scripts/run_capsule.sh 2 giab                 # real GIAB het-SNV+indel benchmark
bash scripts/run_capsule.sh 3 full                 # real export/coverage/query benchmark
```

---

## What's open, honestly

- **Utricularia gibba (`SRR10676752`)**, the 15th locked dataset, has not been run — Claim 1
  is 14/14, not yet 15/15.
- **Claim 2 is validated on chr20 windows and a tetraploid construction, not the full 30×
  individual** the project's own spec commits to.
- **`export`/`coverage`/`query` exist only in this repository**, not in the outer ARCS
  binary.
- **No CI, no top-level license file yet** — both named, neither silently assumed.

Full, current status for each claim: [`docs/CLAIM1_FINAL_VERDICT.md`](docs/CLAIM1_FINAL_VERDICT.md),
[`docs/CLAIM2_FINAL_VERDICT.md`](docs/CLAIM2_FINAL_VERDICT.md),
[`docs/CLAIM3_LOCKED.md`](docs/CLAIM3_LOCKED.md).

---

## Citation

If you use CAPSULE in your research, please cite:

> Thackshanaramana B (2026). *CAPSULE: a unified pseudogenome for lossless FASTQ
> compression, reference-free variant calling, and archive-native addressability.*
> Manuscript in preparation.

---

## Author

**Thackshanaramana B**
SRM Institute of Science and Technology

---

## License

Not yet decided for this repository's own code — see
[`docs/INDUSTRIAL_CHECKLIST_OVERALL.md`](docs/INDUSTRIAL_CHECKLIST_OVERALL.md).
Vendored third-party code keeps its own license: PPMd7 (public domain,
`thirdparty/ppmd/`), FSE/Huf0 (BSD, [`thirdparty/fse/LICENSE`](thirdparty/fse/LICENSE)),
htscodecs/fqzcomp (BSD 3-clause, [`thirdparty/htscodecs/LICENSE.md`](thirdparty/htscodecs/LICENSE.md)).
