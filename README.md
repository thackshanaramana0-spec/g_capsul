# G_CAPSUL: Genomic, Compact, Addressable, Pseudogenome-Structured, Unified Lossless

> **Status:** frozen 2026-09-10 at tag `v1.0.1-capsule` — code and results final.
> Authoritative numbers live in `benchmark/results/`; verification status in
> [`AUDIT.md`](AUDIT.md). Where this file and a result file disagree, the result file wins.

[![CI](https://github.com/thackshanaramana0-spec/g_capsul/actions/workflows/ci.yml/badge.svg)](https://github.com/thackshanaramana0-spec/g_capsul/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
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

G_CAPSUL builds one structure — a **pseudogenome** assembled directly from the reads at
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

- **Smallest lossless archive** on 19 real datasets, 55.0 GB of FASTQ: **19/19 wins against SPRING (−6.03% aggregate) and 19/19 against Genozip (−43.26%)**, with **57/57 archives decoded back to byte-identical input**
- **Reference-free variant calling from the archive itself** — no FASTQ, no reference, no second assembly pass: het-SNV mean F1 **0.876** vs DiscoSNP++ 0.853 and Kmer2SNP 0.475, across four GIAB individuals
- **Multi-allelic and tetraploid calling**: **17/26** multi-allelic sites recovered — a *complete* chr20 census — where DiscoSNP++ emits **zero** multi-allelic records in 3,989; and tetraploid SNV F1 **0.897** vs 0.782 from real HG003+HG004 reads, nothing simulated
- **Archive-native addressability**: pseudogenome export **129–784×** faster than SPAdes and per-base coverage **16–54×** faster than bwa+mosdepth — plus `query`, which resolves a heterozygous locus that **no coordinate can name** (see below)
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

## Start here if you are a new reviewer — run this first

    bash scripts/test_0_scope_and_capability.sh

Self-contained: needs nothing but this checkout and a C++ compiler, runs on
any machine, touches no locked dataset. It probes YOUR machine (no value in
its output is typed in advance), states this build's exact scope — what read
lengths, technologies and alphabets it supports and why, extracted live from
the encoder's own source constants — and then PROVES every boundary claim by
constructing the input and running the real binary: a short-read file
compresses and decompresses losslessly, and long-read-shaped input, an
oversize header, and empty input are each refused cleanly (no crash, no
silent data loss, no false success). If it exits 0, that scope is proven on
your machine right now. Full detail: `industry/README.md`.

## Where everything is

| you want | go to |
|---|---|
| **whether the numbers were audited** | **[`AUDIT.md`](AUDIT.md)** — every table's verification status, one correction, and what the audit did NOT establish |
| the results, as executed | [`benchmark/results/`](benchmark/results/) — 8 CSVs, one run |
| **where any number came from** | [`benchmark/documentation/RESULT_CODE.md`](benchmark/documentation/RESULT_CODE.md) |
| to re-run it all yourself | [`benchmark/documentation/REPRODUCE_EVERYTHING.md`](benchmark/documentation/REPRODUCE_EVERYTHING.md) |
| to rebuild the machine | [`server/`](server/) — setup, tool versions, `verify_environment.sh` |
| the methods, for the paper | [`paper/METHODS.md`](paper/METHODS.md) |
| which scripts / run dirs are live | [`benchmark/documentation/SCRIPT_INVENTORY.md`](benchmark/documentation/SCRIPT_INVENTORY.md), [`RESULTS_INVENTORY.md`](benchmark/documentation/RESULTS_INVENTORY.md) |

Documents outside those folders are the historical record, and several are marked
superseded in place. **Where a document and a result file disagree, the result file wins.**

---

## Results

### Compression (COMPACT)

**Table 1.** Whole-file lossless archive (sequence + names + quality + line 3) on the 19
locked datasets — 55.0 GB of FASTQ spanning bacteria, archaea, viruses, fungi, protists,
plants, animals and human — every archive decoded back to byte-identical input **before**
being counted.

| dataset | G_CAPSUL | SPRING | Genozip | lossless |
|---|---:|---:|---:|:---:|
| ERR5181310 | **8,679,796** | 9,984,000 | 9,222,378 | 3/3 |
| SRR554369 | **57,286,360** | 59,197,440 | 88,813,291 | 3/3 |
| ERR552797 | **46,970,874** | 52,131,840 | 82,799,532 | 3/3 |
| SRR2584863 | **68,429,027** | 74,086,400 | 114,589,268 | 3/3 |
| SRR29296997 | **15,663,251** | 17,541,120 | 27,233,698 | 3/3 |
| ERR12954017 | **15,164,609** | 16,998,400 | 39,429,441 | 3/3 |
| SRR065390 | **887,410,845** | 942,643,200 | 1,597,083,547 | 3/3 |
| SRR40271341 | **40,551,005** | 45,783,040 | 62,260,895 | 3/3 |
| ERR17740259 | **83,196,587** | 92,395,520 | 161,813,831 | 3/3 |
| SRR37283774 | **67,201,215** | 71,342,080 | 94,989,797 | 3/3 |
| DRR976266 | **56,541,178** | 61,757,440 | 201,135,432 | 3/3 |
| SRR36741279 | **106,248,319** | 117,616,640 | 194,658,228 | 3/3 |
| SRR32429602 | **55,118,753** | 57,548,800 | 107,183,883 | 3/3 |
| SRR39257532 | **108,270,384** | 153,856,000 | 227,016,274 | 3/3 |
| SRR10676752 | **1,650,034,057** | 1,718,548,480 | 2,954,900,932 | 3/3 |
| HG002 | **573,767,964** | 598,190,080 | 950,997,605 | 3/3 |
| HG003 | **567,123,999** | 591,144,960 | 944,223,940 | 3/3 |
| HG004 | **606,192,930** | 634,746,880 | 1,041,658,192 | 3/3 |
| HG005 | **997,518,994** | 1,081,384,960 | 1,694,512,440 | 3/3 |
| **aggregate** | **6,011,370,147** | 6,396,897,280 | 10,594,522,604 | **57/57** |

**19/19 wins vs SPRING (−6.03% aggregate). 19/19 vs Genozip (−43.26%). 57/57 LOSSLESS.**

One run, start to finish: `benchmark_1_run.sh` at commit `21ee619`, 5 h 07 m, 0 failures.
Raw CSV: [`benchmark/results/claim1_T1.1_T1.2.csv`](benchmark/results/claim1_T1.1_T1.2.csv).
Traced to the code that produced it in
[`benchmark/documentation/RESULT_CODE.md`](benchmark/documentation/RESULT_CODE.md).

> **Superseded.** This table previously reported *14 datasets, +11.68% vs SPRING and
> +48.93% vs Genozip*. Those numbers are withdrawn: they came from a run that
> `CLAUDE.md` §6.3 records as **VOID for 4 of the 14 datasets**, whose archives could not
> reproduce their input because of four silent data-loss bugs. All four are fixed, and the
> table above is the post-fix sweep with losslessness verified per archive. The margin is
> smaller and it is real.

Separately, sequence-only content against PgRC2's own binary (the closest architectural
relative, GPL-3, run from its own source): **+1.88% aggregate**. That is a different
measurement — PgRC2 stores no names, no quality and no line 3 — and is deliberately not
in the CSV above. Breakdown:
[`docs/CLAIM1_FINAL_VERDICT.md`](docs/CLAIM1_FINAL_VERDICT.md).

### Reference-free variant calling (FAITHFUL)

**Table 2.** Real GIAB HG002–HG005 chr20 at 30×, scored by third-party `rtg vcfeval`
against GIAB v4.2.1 inside the confident regions. **Every G_CAPSUL number here is called
from the compressed archive** — `capsule_decode call archive out.vcf`, with no FASTQ, no
reference and no separate assembly graph. Competitors get the identical reads, truth,
regions, normalisation and scorer; only the caller differs.

| comparison | G_CAPSUL | DiscoSNP++ | Kmer2SNP |
|---|---:|---:|---:|
| het-SNV F1 (mean of 4 individuals) | **0.876** | 0.853 | 0.475 |
| het-indel F1 (mean of 4) | **0.621** | 0.591 | not applicable (SNP-only by construction) |
| multi-allelic sites recovered | **17 / 26** | **0 / 26** | not applicable |
| tetraploid SNV F1 | **0.897** | 0.782 | not applicable |
| tetraploid indel F1 | **0.597** | 0.553 | not applicable |

Per individual, het-SNV: 0.888 / 0.891 / 0.891 / 0.834 — **3 wins and 1 loss**, and the
loss is stated rather than averaged away. HG005 is the only variable-length dataset in the
set (250 bp quality-trimmed, 216 distinct read lengths, against 148 bp fixed elsewhere) and
the caller is tuned for fixed-length reads; at fixed length the same reads give 0.897, our
best. The published number is the untruncated 0.834.
[`docs/HG005_EXPLAINED.md`](docs/HG005_EXPLAINED.md).

The multi-allelic row is a **complete census of chr20**, not a sample: 26 is every
SNV-only multi-allelic site there is. DiscoSNP++'s 0 is structural — across its entire
3,989-record output it emits no record with more than one ALT allele.

**Why this works at all** — the finding the project is built on: *the representation that
compresses a heterozygous site best is the one that conceals it.* Optimal compression puts
each allele on its own internally-consistent contig, so ref and alt reads never share a
coordinate and the variant is **not in the data structure**. Ablation, full chr20:

    neither             F1 0.431   P 0.967   R 0.278
    re-placement only      0.426     0.962     0.274   <- worse than neither
    collapse only          0.648     0.963     0.488
    both                   0.888     0.956     0.830

Synergy, not additivity, and the error *shape* confirms it: without collapse, precision
holds at 0.96 while recall collapses to 0.27. The caller is not mistaken, it is **blind** —
exactly what "the alt reads are on another contig" predicts. **And the correction costs
nothing in compression ratio.**

Raw CSVs: [`benchmark/results/`](benchmark/results/). Survey of why these two competitors
are the applicable ones: [`docs/HET_INDEL_SOTA.md`](docs/HET_INDEL_SOTA.md).

### Archive-native addressability (ADDRESSABLE)

**Table 3.** `export`/`coverage`/`query` served directly from the archive, against the
conventional pipeline that would otherwise compute the same thing.

| operation | vs | speedup |
|---|---|---:|
| export (pseudogenome as FASTA) | SPAdes v4.0.0 | **129–784×** |
| coverage (per-base depth) | bwa + samtools + mosdepth, timed as one pipeline | **16–54×** |
| query (reads at a locus) | — | see below; **not a speed claim** |

**Read the export ratio with its caveat**, which is why the CSV carries output bytes and
row counts for both sides: our export emits 2 records (the pseudogenome), SPAdes emits tens
of thousands of biological contigs. It measures time-to-a-reference-free-coordinate-system
from an archive that had to exist anyway — **not "the same output, faster."** T3.2 and T3.3
carry no such caveat.

**`query` is not a speed claim, and we say so.** `genocat --head=100` extracts in 0.19 s
against our 0.46 s. What `query` does that nothing else can is resolve a **locus**:

| addressing mode, same archive, same 400 GIAB het sites | both alleles returned |
|---|---:|
| by coordinate | 81 / 400 |
| by content (ours) | **345 / 400** |

A heterozygous locus is **not one place** in a compression-optimal pseudogenome — it is N
parallel places (median 4: two haplotypes × two strands), measured up to 18.6 Mb apart. So
the archive's own coordinate system cannot name it, and **adding a coordinate API would not
fix that** — the 81/400 measures precisely that failure. Content addressing resolves every
parallel representative at once. Same mechanism as Claim 2, one layer out.

Cost of all this: the `contig_spans` stream, 232,509 B — **0.041%** of a 573 MB archive.

Full traceability: [`benchmark/documentation/RESULT_CODE.md`](benchmark/documentation/RESULT_CODE.md).
Mechanism: [`docs/CLAIM3_MECHANISM.md`](docs/CLAIM3_MECHANISM.md).

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
git clone https://github.com/thackshanaramana0-spec/g_capsul.git
cd g_capsul
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
command G_CAPSUL supports, organized by claim: [`docs/COMMANDS_REFERENCE.md`](docs/COMMANDS_REFERENCE.md).

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
benchmark/results/           ** THE CITABLE NUMBERS ** — 8 CSVs, one run, 19 datasets
results/                     working output of 12 runs, incl. reverted work — NOT citable
                              (see benchmark/documentation/RESULTS_INVENTORY.md)
NEW_DATASET_LOCKED.md        ** THIS repo's locked set ** — 15 non-human + 4 GIAB human = 19
DATASET_LOCKED.md            the OUTER ARCS project's 17-accession list — does NOT govern this repo
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

The four items that used to sit here — the unrun 15th dataset, Claim 2 on windows rather
than a full individual, no CI, no license — are all closed. `SRR10676752` is in the sweep,
Claim 2 runs on full chr20 at 30× for four individuals, and both `LICENSE` and
`.github/workflows/ci.yml` exist. These are what is genuinely open:

- **Human validation is chr20 only, at 30×.** Claims 2 and 3 are measured on four GIAB
  individuals across the whole chromosome, not on windows — but it is one chromosome,
  because the archives are chr20. Nothing here is evidence about whole-genome behaviour.
- **No non-human diploid variant validation.** The 15 non-human datasets carry Claim 1
  only; there is no comparable truth set for them.
- **HG005 loses, and we publish the loss.** The caller is tuned for fixed-length reads;
  HG005 is the only variable-length dataset. Cause identified by controlled experiment
  ([`docs/HG005_EXPLAINED.md`](docs/HG005_EXPLAINED.md)), not fixed.
- **Compression is slower than SPRING** — 2.99× on compress, 1.79× on decompress. We are
  faster on 3 of 19 and lighter on 9 of 19. Assembly costs time; that is the trade.
- **`query` is not a speed win** against `genocat --head` (0.46 s vs 0.19 s), and Table 3
  says so rather than leading with the indexed 0.03 s.
- **The export ratio is not "the same output, faster."** See the caveat under Table 3.

Full, current status for each claim: [`docs/CLAIM1_FINAL_VERDICT.md`](docs/CLAIM1_FINAL_VERDICT.md),
[`docs/CLAIM2_FINAL_VERDICT.md`](docs/CLAIM2_FINAL_VERDICT.md),
[`docs/CLAIM3_LOCKED.md`](docs/CLAIM3_LOCKED.md).

---

## Citation

If you use G_CAPSUL in your research, please cite:

> Thackshanaramana B (2026). *G_CAPSUL: a unified pseudogenome for lossless FASTQ
> compression, reference-free variant calling, and archive-native addressability.*
> Manuscript in preparation.

---

## Author

**Thackshanaramana B**
SRM Institute of Science and Technology

---

## License

**MIT** — see [`LICENSE`](LICENSE) for the full text, which also documents every
vendored dependency's own terms.

Vendored third-party code keeps its own license, all of them MIT-compatible:

| component | license | where |
|---|---|---|
| PPMd7 (LZMA SDK) | public domain | [`thirdparty/ppmd/LICENSE`](thirdparty/ppmd/LICENSE) |
| FSE / Huff0 | BSD | [`thirdparty/fse/LICENSE`](thirdparty/fse/LICENSE) |
| htscodecs / fqzcomp | BSD 3-clause | [`thirdparty/htscodecs/LICENSE.md`](thirdparty/htscodecs/LICENSE.md) |
| liblzma | public domain | runtime dependency, linked with `-llzma` |

**PgRC2 (GPL-3) is deliberately NOT vendored.** It is used only as an external
comparison binary, cloned separately; no PgRC2 source is included in or linked
into this repository, so its GPL-3 terms do not attach here. That separation
was a design decision, not an accident — see `docs/REIMPL_NOTES.md`, which
records that PgRC2's assembler was reimplemented from the algorithm rather than
copied, precisely so this stays true.
