# Limitations and scope

This document lists what the evaluation covers, what it does not cover, and the costs that come
with the design. It applies to the results in the manuscript and in `results/`. For the input that
the program accepts, see [`PROJECT_SCOPE.md`](PROJECT_SCOPE.md).

## 1. What was evaluated

| Area | Evaluated on | Not evaluated |
|---|---|---|
| Compression | 19 public short-read datasets (59.1 GB) from 15 non-human and 4 human sources, single-end | Long reads, paired-end mate relationships, whole human genomes |
| Variant calling | Chromosome 20 of GIAB HG002 to HG005 at about 30× coverage, and a 400 kb tetraploid mixture of HG003 and HG004 reads | Whole genomes, non-human organisms, structural variants, population-scale inference |
| Export | 6 datasets that SPAdes can also process | The 13 other datasets (export speed was measured only where SPAdes could be run) |
| Coverage | 6 datasets that bwa can also process | As above |
| Range query timing | All 19 datasets, no comparator | |
| Locus retrieval | Four GIAB individuals, one 600 kb window of chromosome 20 (3.0 to 3.6 Mb), plus one extra HG005 window | Other loci, other organisms, whole-genome retrieval |

Ground truth for variant calling and locus retrieval exists only for the human GIAB individuals,
because no comparable validated heterozygous truth set exists for other organisms. Restricting
scoring to chromosome 20 and to the GIAB confident regions follows common practice. The
tetraploid table covers a 400 kb region and should not be read as chromosome-scale.

## 2. Costs of the design

| Area | Cost | Detail |
|---|---|---|
| Compression speed | Slower than both comparators | G_CAPSUL is the fastest tool in none of the 38 compression and decompression cells. Genozip is fastest in 37 and SPRING in 1. This is the cost of the structure that variant calling and retrieval reuse. |
| Memory | Higher than conventional compressors on larger datasets | Peak memory grows with input size, roughly 1.2 GB plus 0.64 GB per GB of raw FASTQ across the 19 datasets (R² 0.97), and the ratio to input falls from 4.3 to 0.7. Peak memory in the result files is for the compression phase only. |
| Low coverage | Loses to SPRING | At about 2.6× coverage, which is outside the 19-dataset comparison, G_CAPSUL is about 3.97% larger than SPRING. Most reads then fail to chain and are handled by region-to-region matching, which costs about 21 bits per mismatch against about 2 bits for read-to-read matching. This is not a tunable parameter. |
| Retrieval index | Disk space | The optional FM-index is 0.9 to 2.6 MB per individual. The archives of the 600 kb retrieval window are 0.48 to 0.75 MB, so the index is 1.9 to 3.5 times as large as the archive it serves. It grows with the share of reads that deviate from the pseudogenome. It is built after the archive and is not part of the reported archive size. |
| Locus retrieval specificity | Extra reads at control sites | Adding the FM-index raises hits at homozygous control sites from 29, 35, 29 and 66 to 48, 52, 49 and 87 of 400 (HG002 to HG005). Completeness is obtained at the cost of specificity. |
| Parallelism | Partly serial | The scope check and the pseudogenome rebuild are serial. Chaining, mapping and stream coding use threads and forked candidate encodings, and each candidate is mostly serial. |

## 3. Limits on individual results

- **Variant calling on HG005.** DiscoSNP++ scores higher on HG005 for both SNVs (0.863 against
  0.834) and indels (0.605 against 0.593). HG005 is the only variable-length dataset (250 bp,
  quality-trimmed, 216 distinct read lengths). Truncating its reads to 148 bp in a diagnostic run
  raised the G_CAPSUL SNV F1 to 0.897, which points at read-length distribution as the cause. The
  reported value is the untruncated one.
- **Multi-allelic sites.** The 17 of 26 result is a complete census of chromosome 20 for one
  individual (HG002), not a sample. Nine sites were not recovered.
- **Export correctness.** G_CAPSUL is better than SPAdes on genome fraction (98.7% against 98.3%)
  and indel rate, and worse on mismatch rate, duplication ratio, N50 and misassemblies. Closing the
  gap requires repeat resolution and error correction, which a pseudogenome built to minimize size
  does not do. Four attempts to close it (region splitting, a parameter sweep, contig
  deduplication and consensus polishing) were tried and did not help.
- **ERR5181310 export.** This amplicon dataset has more than 1000-fold duplicated coverage of a
  29,903 bp genome. The exported pseudogenome is 24,072,995 bp, about 800 times the genome,
  although duplicate collapse already reduces the reads from 1,825,142 to 152,085. The working
  explanation is that near-duplicate reads overlap the chain tail by chance, so the chain does
  not terminate. This was not confirmed by tracing the code. The dataset is reported for speed
  and is not scored for correctness against SPAdes.
- **Query cost.** `query` reads the whole archive and has no windowed decode. In a self-referential
  compressor the matcher optimizes match length and ignores match distance. On a measured E. coli
  archive 44% of references point to a copy about 15 MB further away than an identical nearer
  copy that the sampled index did not offer. A 100 kb window therefore has a large transitive
  closure (179,322 bytes, 0.67% of the pseudogenome). A nearest-copy tie-break was implemented
  and rejected because ties are only 0.035% of candidates. Increasing the index sampling density
  is the remaining option and trades against mapping speed. The decoded-index cache makes repeated
  queries about 17 times faster and does not remove the cost.
- **Baselines.** No speed or index-size comparison against BEETL-fastq, CIndex or BFQzip was run.
  See [`novelty_and_prior_art.md`](novelty_and_prior_art.md).

## 4. Validation of reconstruction

The benchmark decoded every archive before counting it and compared the sorted (name, sequence,
quality) records with the input. All 57 archives (19 datasets, 3 tools) passed. This check does not
test read order or line 3. The ordered byte-for-byte check,
`scripts/build_scripts/verify_roundtrip.sh`, was run on simulated test inputs covering variable
lengths, N bases, duplication, short and long reads and a single-read file, and on benchmark
dataset SRR29296997. It was not run on the other 18 datasets.

## 5. Statistical methodology

No confidence intervals or significance tests are reported, because comparable compression and
reference-free calling comparisons do not report them for these comparison types. Sensitivity to
the input is addressed with sweeps: the MAXMAP and MINOV search for compression, the coverage
sweep (10×, 15×, 20× and 30×) and the duplicate-fraction plateau for variant calling, and the
per-dataset speed-up range for export and coverage. Each of the four GIAB individuals is a
separate result and the table keeps them separate.

## 6. Implementation notes

- **`.sites` pileup strand.** The per-position tally written by `index` ignores read strand, while
  `query` corrects for it. A fix was written and passed its checks, and was not adopted, because
  it would change the published native-pileup figures and those were not re-measured. The tally
  site carries a comment describing this.
- **Consensus polishing from the pileup.** Polishing an exported sequence with the archive's own
  pileup makes correctness worse (11.86 to 407 mismatches per 100 kb), because the pileup reflects
  compression-driven placement and not homology. The negative result is recorded as a code comment.
- **Mapping ceiling.** Raising `MAXMAP` to improve an export metric would cost about 11.6% of the
  compression ratio, so it was not adopted.
- **`CAPS_PILEUP=1`.** The query sidecar needs it so that reads carry their deviations. Without it
  `query` returns the consensus and not the reads.

## 7. Platform and tooling limits

Native Windows is not supported, because the encoder and decoder use POSIX process calls. Use
WSL or Docker. macOS needs Homebrew GCC for OpenMP. See [`PROJECT_SCOPE.md`](PROJECT_SCOPE.md).
