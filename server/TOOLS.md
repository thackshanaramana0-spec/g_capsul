# TOOLS — exact versions, provenance, and what each is for

Versions are what was **actually installed on the machine that produced every
number in `../benchmark/results/`**, read off the binaries rather than from a
requirements file.

## Platform

    OS        Ubuntu 24.04
    CPU       AMD EPYC 9555, 12 vCPU exposed
    caches    L1d 768 KiB, L2 6 MiB, L3 192 MiB
    RAM       82 GB
    disk      233 GB SSD

The L3 size is worth recording: at 192 MiB the 12.6 MB `readMM` array is
cache-resident, so profiling conclusions about cache misses on this box do not
transfer to a machine with a small L3.

## Toolchain

| tool | version | source |
|---|---|---|
| gcc / g++ | 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | apt |
| cmake | 3.28.3 | apt |
| python3 | 3.12.3 | apt |

GCC 13 is the reason SPRING needs a `<cstdint>` patch — see `SETUP_FULL.md` §3.5.

## Compression baselines (Claim 1)

| tool | version | source | role |
|---|---|---|---|
| **SPRING** | built from `github.com/shubhamchandak94/SPRING` @ master | source + patch | primary FASTQ competitor |
| **Genozip** | 15.0.87, `distribution=github` | source, `nasm` needed | second FASTQ competitor |
| **PgRC2** | 2.0 (Kowalski & Grabowski, 2024) | `arcs-clean/method_c`, GPL-3, cloned not vendored | DNA-stream reference point |

**Genozip needs an interactive activation** and a network-reachable
`/dev/stdout`; its archive drifts ~37 B between runs. **SPRING's `-g` flag
means "input is gzipped"** and must not be passed to our plain `.fq` inputs.
**PgRC2 is not a FASTQ compressor** — its output is bare DNA, no names, no
quality; it cannot round-trip a FASTQ file. All three facts are load-bearing
for how the results are read.

## Variant-calling baselines (Claim 2)

| tool | version | source | role |
|---|---|---|---|
| **DiscoSNP++** | `/root/DiscoSnp` | source | reference-free SNV/indel competitor |
| **Kmer2SNP** | `/root/Kmer2SNP` + conda env `kmer2snp_r` | source | third T2.1 arm |
| **KMC** | 3.2.4 (2024-02-09) | apt | k-mer counting for Kmer2SNP, replacing its broken DSK wrapper |
| **rtg-tools** | RTG Tools 3.13 | GitHub release | `vcfeval` — the scorer for every F1 in Claim 2 |

DiscoSNP++ must be **reachable on PATH**, is always called with `-G <ref.fa>`,
and its documented POS off-by-one is corrected in our harness. Kmer2SNP must be
run with the conda env's python and be fed `k_K_pair.snp` — see `SETUP_FULL.md`
§3.8-3.9 for the three traps.

## Claim 3 baselines

| tool | version | source | role |
|---|---|---|---|
| **SPAdes** | 4.0.0 | `/root/SPAdes-4.0.0-Linux` | T3.1, de novo assembly baseline |
| **bwa** | 0.7.17-r1188 | apt | T3.2, alignment half of the conventional route |
| **samtools** | 1.19.2 | apt | T3.2, sort + index |
| **mosdepth** | 0.3.6 | apt | T3.2, per-base depth |

T3.2's baseline is timed as **one pipeline** (align → sort → index → depth),
because that is the conventional route; timing only part of it would flatter us.

## Supporting

| tool | version | role |
|---|---|---|
| bcftools | 1.19 | VCF normalisation before scoring |
| tabix | 1.19 (htslib) | truth-set region extraction |
| seqtk | apt | T2.2's coverage-sweep subsampling (fixed seed) |
| SRA toolkit | current | `prefetch`, `fasterq-dump` — **symlink, never copy** |
| AWS CLI | v2, official installer | S3 open-data mirror for the two large SRA files |

## Our own binaries

| binary | built by | what it is |
|---|---|---|
| `best106` | `scripts/build106.sh` | the encoder (stage 106, `stages/106_inprocess.cpp`) |
| `capsule_decode` | `scripts/build_decode.sh` | decoder + `export` / `coverage` / `query` / `call` / `index` |
| `arcs` | `cmake --build build` in ARCS | the outer project's binary |

**Always build with the scripts.** `build106.sh` exists because the command
once lived only in shell history, which is how `-fopenmp` went missing —
silently discarding the parallel sweep and costing 45% of wall time.

Optional build flags, both off by default and neither shipped:

    -DMAPFUNNEL    mapping-funnel instrumentation (CAPS_MAPDBG=1 at run time).
                   Compile-time because a runtime branch in that loop measured
                   ~1.8% -- an instrument must not cost what it measures.

## Reproducibility of the tools themselves

Verified on this machine and recorded in
`../benchmark/documentation/REPRODUCIBILITY.md`: our archive and our VCFs are
**byte-identical across 2/3/7/12 cores** and across candidate-concurrency and
memory-budget settings. Only wall time and peak RSS vary. Genozip's archive is
the one exception that varies by design.
