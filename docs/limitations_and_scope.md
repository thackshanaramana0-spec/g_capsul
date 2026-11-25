# Limitations and scope

This document states what the evaluation covers, where the current implementation falls short,
and which directions follow from those shortfalls. It is based on the Discussion of the
manuscript and on the result files. For the input the program accepts, see
[`PROJECT_SCOPE.md`](PROJECT_SCOPE.md), and for the numbers, see
[`results_tables.md`](results_tables.md).

## 1. Scope of the evaluation

| Area | Evaluated on | Outside the evaluation |
|---|---|---|
| Compression | 19 public short-read datasets (59.1 GB), single-end, from 15 non-human and 4 human sources | Long reads, paired-end mate relationships, whole human genomes |
| Variant calling | Chromosome 20 of GIAB HG002 to HG005 at about 30× coverage, and a 400 kb tetraploid mixture of HG003 and HG004 reads | Whole genomes, non-human organisms, structural variants, population-scale inference, naturally polyploid organisms |
| Export | 6 datasets that SPAdes can also process | The other 13 datasets |
| Coverage | 6 datasets that bwa can also process | The other 13 datasets |
| Range query | All 19 datasets, timed without a comparator | |
| Locus retrieval | Four GIAB individuals, chr20:3.0 to 3.6 Mb, plus one further HG005 window | Other loci, other organisms, whole-genome retrieval |

Ground truth for variant calling and locus retrieval exists only for the human GIAB individuals,
because no comparable validated heterozygous truth set exists for other organisms. Restricting
scoring to chromosome 20 and to the GIAB high-confidence regions follows common practice. The
tetraploid table covers a 400 kb region and does not establish chromosome-scale performance. It
shows that the ploidy-aware machinery works beyond the diploid case, and it does not show general
performance in naturally polyploid organisms.

Reference-based workflows, such as GATK and DeepVariant, remain appropriate when a high-quality
reference and the additional information those workflows use are available. G_CAPSUL addresses a
different question: how much useful organization can be recovered from the reads themselves,
kept through reference-free compression and carried into later analysis.

## 2. Computational cost

Constructing a pseudogenome requires the compressor to infer a substantial amount of organization
before the archive can be written. Reads have to be associated with the emerging sequence
structure, their residual differences have to be represented, and the relationships needed for
reconstruction and later analysis have to be kept.

| Aspect | Observation |
|---|---|
| Speed | G_CAPSUL is the fastest tool in none of the 38 compression and decompression cells measured. Genozip is fastest in 37 and SPRING in 1. |
| Memory | Peak memory grows with input size, roughly 1.2 GB plus 0.64 GB per GB of raw FASTQ across the 19 datasets (R² 0.97). The ratio to input falls from 4.3 to 0.7. The memory in the result files is the compression phase only. |
| Parallelism | The scope check and the pseudogenome rebuild are serial, chaining is partly parallel, and each forked candidate encoding is mostly serial. Mapping, MEM search and stream coding use threads. |
| Low coverage | At 2.6× coverage, outside the 19-dataset comparison, G_CAPSUL is about 3.97% larger than SPRING. Too few reads chain to support region-scale organization, so more sequence is represented by region-to-region matching. That costs about 21 bits per mismatch, against about 2 bits for SPRING's read-to-read matching. The cause is architectural and not a tunable parameter. |

These costs come from the current implementation. They do not follow from any requirement that a
reusable archive use the same algorithms. Parallel construction and more compact temporary
structures are direct ways to reduce them without changing the archive model.

## 3. Variant calling

| Limitation | Detail |
|---|---|
| Whole-genome calling | The k-mer spill representation needs large intermediate storage at larger scales, so whole-genome calling is not practical in the current form. |
| Read geometry | DiscoSNP++ scores higher on HG005 for SNVs (0.863 against 0.834) and for indels (0.605 against 0.593). HG005 is the only variable-length dataset, with 250 bp quality-trimmed reads of 216 distinct lengths. Truncating its reads to 148 bp in a diagnostic run raised the SNV F1 to 0.897, which isolates read-length distribution as the cause. The reported value is the untruncated one. |
| Indel candidates | An indel is found only if its alternative sequence appears in the candidate space before scoring. Heterozygous indel calling was once a closed loss under the original candidate generator, and no threshold or filter applied to that output could recover it. A different candidate channel, positional clustering, closed the gap. Some remaining failures still occur before scoring. |
| Multi-allelic sites | The result of 17 of 26 recovered sites is a complete census of chromosome 20 for one individual. Nine sites were not recovered. |
| Dependence on compression | The caller is influenced both by how compression organizes the sequence and by the characteristics of the reads that organization is inferred from. |

## 4. Export, query and retrieval

| Limitation | Detail |
|---|---|
| Export correctness | G_CAPSUL is better than SPAdes on genome fraction (98.7% against 98.3%) and indel rate and worse on mismatch rate, duplication ratio, N50 and misassemblies. Closing the gap needs repeat resolution and error correction, which a pseudogenome built to minimize size does not do. Four attempts (region splitting, a parameter sweep, contig deduplication, consensus polishing) did not close it. |
| Amplicon data | ERR5181310 has more than 1000-fold duplicated coverage of a 29,903 bp genome. Its exported pseudogenome is 24,072,995 bp, about 800 times the genome, although duplicate collapse already reduces the reads from 1,825,142 to 152,085. The working explanation is that near-duplicate reads overlap the chain tail by chance, so chains do not terminate. This was not confirmed by tracing the code. The dataset is reported for speed and is not scored for correctness. |
| Query is not an indexed database | The base archive was not designed as a random-access database. `query` reads the whole archive and has no windowed decode. In a self-referential compressor the matcher optimizes match length and ignores match distance. On a measured E. coli archive 44% of references point to a copy about 15 MB further away than an identical nearer copy that the sampled index did not offer, so a 100 kb window has a large transitive closure (179,322 bytes, 0.67% of the pseudogenome). A nearest-copy tie-break was implemented and dropped, because ties are only 0.035% of candidates. |
| Index cost | Unlike BEETL-fastq, which builds its index during compression, G_CAPSUL needs a post-processing step. It keeps the decoded pseudogenome for repeated queries and builds an FM-index over the reads that deviate from it. The decoded sidecar is 3.0 to 3.4 times the archive, and the FM-index is 0.9 to 2.6 MB per individual, which is 1.9 to 3.5 times the 0.48 to 0.75 MB archives of the retrieval window. Both grow with the share of reads that deviate from the pseudogenome and are not part of the reported archive size. The decoded-index cache makes repeated queries about 17 times faster. |
| Specificity | Adding the FM-index raises hits at homozygous control sites from 29, 35, 29 and 66 to 48, 52, 49 and 87 of 400 (HG002 to HG005). |
| Missing comparisons | No speed or index-size comparison against BEETL-fastq, CIndex or BFQzip was run. See [`novelty_and_prior_art.md`](novelty_and_prior_art.md). |

## 5. Reconstruction and validation

The benchmark decoded every archive before counting it and compared the sorted (name, sequence,
quality) records with the input. All 57 archives (19 datasets, 3 tools) passed. This check does not
test read order or line 3. The ordered byte-for-byte check,
`scripts/build_scripts/verify_roundtrip.sh`, was run on simulated test inputs covering variable
lengths, N bases, duplicates, short and long reads, a single-read file and a shuffled read order,
and on benchmark dataset SRR29296997. It was not run on the other 18 datasets.

No confidence intervals or significance tests are reported, because comparable compression and
reference-free calling comparisons do not report them for these comparison types. Sensitivity to the
input is addressed with sweeps: the `MAXMAP` and `MINOV` search for compression, the 10×, 15×, 20×
and 30× coverage sweep and the duplicate-fraction plateau for calling, and the per-dataset speed-up
range for export and coverage. Each of the four GIAB individuals is a separate result.

## 6. Implementation notes

- **`.sites` pileup strand.** The per-position tally written by `index` ignores read strand, while
  `query` corrects for it. A fix passed its checks and was not adopted, because it would change the
  published native-pileup figures and those were not re-measured. A comment at the tally site
  describes this.
- **Consensus polishing from the pileup.** Polishing an exported sequence with the archive's own
  pileup makes correctness worse (11.86 to 407 mismatches per 100 kb), because the pileup reflects
  compression-driven placement and not homology. The negative result is recorded as a code comment.
- **Mapping ceiling.** Raising `MAXMAP` to improve an export metric would cost about 11.6% of the
  compression ratio, so it was not adopted.
- **`CAPS_PILEUP=1`.** The query sidecar needs it so that reads carry their deviations. Without it
  `query` returns the consensus and not the reads.
- **DiscoSNP++ output.** Its VCF positions are shifted by one against its own description, found by
  reading its source. All comparisons use positions corrected for this. Uncorrected, its measured
  accuracy would be close to zero.
- **Platform.** Native Windows is not supported, because the encoder and decoder use POSIX process
  calls. Use WSL or Docker. macOS needs Homebrew GCC for OpenMP.

## 7. Extensions that follow from these limits

| Direction | Addresses |
|---|---|
| Parallelize pseudogenome construction and reduce working memory | The largest runtime and memory cost (section 2) |
| Replace the k-mer spill representation | Whole-genome calling (section 3) |
| Expose indel candidates earlier in the pipeline | Indel failures that occur before scoring (section 3) |
| Retain paired-end mate relationships in the pseudogenome | The current loss of pairing information, which would give compression, calling and retrieval another source of structure |
| Build the retrieval index during compression | The post-processing step and its storage (section 4) |
| Multi-sample analysis | The archive already keeps within-sample relationships that conventional workflows rebuild |
| Longer reads, other sequencing technologies and general ploidy | These need representations and validation suited to those settings and not extrapolation from short-read results |

## 8. Verification and corrections

The result tables were verified in two tiers. Ten tables were checked in total. Four were
re-executed against an independent measurement (multi-allelic, tetraploid, the HG002 row of
exact-match recall, and the HG002 rows of SNV and indel calling). The rest were checked by parsing the
run log and comparing it with the result file. The second tier is weaker, because the log line and
the file row are written from the same variables in the same run.

**Defects found and corrected**

1. Multi-allelic recovery was first reported as 21 of 26. The denominator was re-derived from the
   GIAB HG002 truth set (952 multi-allelic truth sites, 26 of them SNV-only) and matched. Three
   re-runs, on identical reads, at the preceding commit and with the original working directory
   traced, all gave 17 of 26. The 21 came from a stale 5 Mb window result carried forward as the
   full-chromosome number. It was corrected to 17 of 26 everywhere.
2. Exact-match and locus retrieval need `CAPS_PILEUP=1` when the query sidecar is built. Without it
   `query` returns the consensus and the coordinate-only arm scores zero by construction.
3. Two ignore patterns hid evidence files. `*.log` hid reproduction logs and `*.vcf` hid the truth
   file used for the multi-allelic correction. Both were fixed with targeted exceptions.
4. The file-integrity manifest included untracked files, so a fresh clone failed with eight
   missing-file errors. It was regenerated from `git ls-files`.

**Limits of the verification.** No second physical machine was available. A fresh clone with a
restricted core count was used and it passed. The exact-match rows for HG003, HG004 and HG005 were
not re-executed independently of the HG002 run. One person carried out the verification, and a second
reviewer is the most useful next step. One incorrect number was found in ten tables, which does not
show that there was only one.

**Reproducibility across cores.** Every number was produced on one 12-core machine with 82 GB of
memory. Re-encoding E. coli (SRR2584863) on 2, 3, 7 and 12 cores gave the same archive each time
(68,429,027 bytes), and only the wall time changed (36.17 s down to 12.38 s). Decoded reads were
MD5-identical at 2, 7 and 12 cores.

**Superseded figures.** Older documents in `docs/extras/` may quote the values on the left.

| Superseded value | Measure | Replacement | Reason |
|---|---|---|---|
| 80.8% genome fraction, 90% unaligned | Export correctness | 98.698% genome fraction | Measured on an export that joined all contigs into two records |
| 345 of 400 | Locus retrieval | 1,452 of 1,452 | Measured with an earlier query that had a strand error and covered one individual |
| 0 of 400 | Coordinate-only retrieval | a small fraction of sites | The query emitted the consensus instead of the reads |
| 21 of 26 | Multi-allelic recovery | 17 of 26 | No supporting raw log |
| 5 of 111 and 11 of 18 | Multi-allelic recovery | 17 of 26 | One used a check that cannot score indel sites, and one was not reproducible |
| 0.637 and 0.666 | Heterozygous indel F1 | 0.621 | Measured on 8 chromosome 20 windows |
| 0.890 | Heterozygous SNV F1 | 0.876 | Measured before the 19-dataset sweep |
| 0.836 | Tetraploid SNV F1 | 0.897 | Measured before the final sweep |
| 555 to 656 times | Export speed-up | 129 to 784 times | Measured before the final sweep |
| 1.0000 recall without an index | Exact-match recall | 0.60 to 0.85 archive only, 1.0000 with the FM-index | The earlier figure held only for HG002 |
