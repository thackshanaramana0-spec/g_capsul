# Related work and positioning

This document places G_CAPSUL next to existing tools, one capability at a time. For each area it
states what the other tools do, how G_CAPSUL differs, and how that statement was established.
Reference details are in the manuscript bibliography.

**Evidence labels used in the tables**

| Label | Meaning |
|---|---|
| Measured | Compared on the same data in this repository, with the result in `results/`. |
| Source read | Established by reading the other tool's source code or primary paper. |
| Documented | Taken from the tool's documentation or paper and not checked in its source. |
| Inferred | An argument from documented design or observed output, not directly verified. |

## 1. Summary

| Capability | Closest existing tools | What G_CAPSUL adds | Main evidence |
|---|---|---|---|
| Lossless FASTQ compression | SPRING, Genozip, PgRC2 | A retained placement structure in the archive, at a size smaller than SPRING and Genozip on all 19 datasets | Measured |
| Variant calling without a reference | DiscoSNP++, Kmer2SNP | Calling from the compressed archive, without the original reads | Measured |
| Sequence export | SPAdes | Export from the archive, much faster, with different correctness trade-offs | Measured |
| Per-base coverage | bwa, samtools, mosdepth | Coverage from retained placements, without alignment | Measured |
| Retrieval by content | BEETL-fastq, CIndex | An optional index over the reads that deviate from the pseudogenome | Measured and Source read |
| Retrieval by locus | none identified | Recovery of both allelic read sets at a locus | Measured |

The combination of the capabilities in one archive is the contribution. Several individual
capabilities have established tools, and the tables below say where G_CAPSUL is better, worse or
only different.

## 2. FASTQ compression

| Tool | Approach | Relation to G_CAPSUL | Evidence |
|---|---|---|---|
| SPRING (Chandak et al., 2019) | Reordering and hash-based read matching | Primary comparator for whole-FASTQ size. G_CAPSUL is smaller on 19 of 19 datasets (aggregate −6.03%). | Measured |
| Genozip (Lan et al., 2021) | Universal genomic data compressor | Second comparator. G_CAPSUL is smaller on 19 of 19 datasets (aggregate −43.26%). | Measured |
| PgRC and PgRC2 (Kowalski and Grabowski, 2020, 2025) | Pseudogenome-based read compression | Closest architecture. PgRC2 stores sequence only (no names, quality or line 3), so the comparison is restricted to sequence and to 7 datasets it can process: +1.88% aggregate, 6 wins and 1 loss, with G_CAPSUL about 1.7 times slower and up to 2.5 times heavier. | Measured |
| CRAM 3.1 (Bonfield, 2022) and reference-based compressors (Fritz et al., 2011) | Alignment against a reference genome | Need a reference. G_CAPSUL needs none. | Documented |
| fqzcomp (Bonfield and Mahoney, 2013) | Context-mixing coding of FASTQ fields | The quality coder is vendored from it (see `THIRD_PARTY_NOTICES.md`), so quality compression is not claimed as new. | Source read |
| DSRC 2, FaStore, SCALCE, minicom, NanoSpring | Other FASTQ and read compressors | Cited as context. Not benchmarked here. | Documented |

Wall time and memory are the cost of the added structure. G_CAPSUL is the fastest tool in none of
the 38 compression and decompression cells measured, and uses more memory than SPRING and Genozip
on larger datasets. See [`limitations_and_scope.md`](limitations_and_scope.md).

## 3. Reference-free variant calling

| Tool | Approach | Relation to G_CAPSUL | Evidence |
|---|---|---|---|
| DiscoSNP++ (Peterlongo et al., 2017) | Bubbles in a de Bruijn graph built from raw reads | Primary comparator. It is the other general-purpose reference-free caller that handles heterozygous indels from one diploid sample. Mean het-SNV F1 is 0.853 against 0.876 for G_CAPSUL, and DiscoSNP++ is ahead on HG005. | Measured |
| Kmer2SNP (Li et al., 2020) | k-mer matching for SNPs | Second baseline. SNV only by design. Mean het-SNV F1 0.475. | Measured |
| Cortex (Iqbal et al., 2012) | Coloured de Bruijn graphs | Needs several samples or a reference to genotype. Cited as background. | Documented |
| eBWT2SNP (Prezza et al., 2019) and BWT-based discovery (Prezza et al., 2020) | Positional clustering on the extended BWT | The positional-clustering idea inspired one indel channel in the G_CAPSUL caller. The precision reported for eBWT2SNP comes from simulated chromosomes. On the one real-data chromosome in that paper, DiscoSNP++ was more precise (74.57% against 66.62%). | Source read |

Two further points about DiscoSNP++:

- It emitted no multi-allelic record across its entire output on the multi-allelic census (0 of
  3,989 records, and 0 of 26 multi-allelic sites). The explanation, that its bubble output is
  biallelic by design, is **inferred** from its documented behaviour and was not verified in its
  source.
- Its VCF positions are shifted by one relative to its own description. All comparisons use
  positions corrected for this shift. Without the correction its accuracy would be measured as
  close to zero.

## 4. Assembly and export

| Tool | Approach | Relation to G_CAPSUL | Evidence |
|---|---|---|---|
| SPAdes (Bankevich et al., 2012) | De Bruijn assembly with repeat resolution and error correction | Comparator for sequence export. Export is 129 to 784 times faster. Genome fraction is 98.7% against 98.3% and the indel rate is lower. SPAdes is better on mismatch rate, duplication ratio, N50 and misassemblies. | Measured |
| Eulerian-path assembly (Pevzner et al., 2001) | Graph formulation of assembly | Background for de Bruijn assembly. | Documented |
| Purge Haplotigs (Roach et al., 2018) and Redundans (Pryszcz and Gabaldón, 2016) | Remove duplicated heterozygous contigs from an assembly | They show that heterozygous regions often assemble into separate contigs, which is the phenomenon described in [`central_insight_and_mechanism.md`](central_insight_and_mechanism.md). They remove the duplicate to obtain one representative sequence. A lossless archive cannot delete a contig, because that would delete an allele, so G_CAPSUL keeps both and reconciles them by content. | Documented |

SPAdes is the right comparator for 5 of the 6 export datasets. The sixth, ERR5181310, is amplicon
data with more than 1000-fold duplicated coverage, which is outside SPAdes's own documented scope
and where the exported pseudogenome is far larger than the 29,903 bp genome. That dataset is
reported for speed only. See [`limitations_and_scope.md`](limitations_and_scope.md).

## 5. Coverage

| Tool | Approach | Relation to G_CAPSUL | Evidence |
|---|---|---|---|
| bwa-mem (Li, 2013), samtools (Li et al., 2009) and mosdepth (Pedersen and Quinlan, 2018) | Align, sort, then count depth | The standard pipeline for the same output. G_CAPSUL reads depth from retained read positions and lengths, so no alignment is needed. It is 16 to 54 times faster, with exact depth values. | Measured |

## 6. Retrieval by content and by locus

| Tool | Approach | Relation to G_CAPSUL | Evidence |
|---|---|---|---|
| BEETL (Cox et al., 2012) and BEETL-fastq (Janin et al., 2014) | BWT-based compressed archive of reads | Search by sequence content is already available from BWT indexes, so it is not new. A BWT stores reads as independent strings and has no coordinate system. | Source read |
| CIndex (Huo et al., 2022) | Compressed indexes for FASTQ retrieval | Same class. Content search from a compressed index. | Documented |
| FM-index (Ferragina and Manzini, 2000) | Compressed full-text index | G_CAPSUL builds one after compression, over the reads that differ from the pseudogenome. | Measured |
| Population BWT (Dolle et al., 2017) | Compressed structure over many samples | Targeted queries at pre-specified sites, and quality scores are not kept. | Documented |
| Compressive genomics (Loh et al., 2012) | Search on compressed sequence databases | Same general idea of computing on compressed data. It accelerates search across many reference databases, whereas G_CAPSUL serves variant evidence from one sample's reads. | Documented |

What G_CAPSUL adds in this area is retrieval by **locus** rather than by exact match. A query
returns the reads that cover a site, including reads whose sequence differs from the pseudogenome
and therefore cannot be found by exact search. Measured: archive-only search finds 0.60 to 0.85 of
reads containing a probe, and the archive plus the FM-index finds all of them (1.0000). Recovery of
both alleles at 1,452 of 1,452 sites, with a cost in extra reads at homozygous control sites,
is reported in [`claim3/t35_locus_retrieval_claim3.md`](claim3/t35_locus_retrieval_claim3.md).

The reason a BWT archive does not offer locus retrieval is its data model, not the lack of an
index. A BWT sorts suffixes and never computes where a read sits. Position is a free byproduct of
overlap chaining in G_CAPSUL, so the FM-index here was a small addition. For BEETL this reasoning
was checked against how the method works, and it was not tested by implementing it.

No speed or index-size benchmark against BEETL-fastq or CIndex was run. This is the most useful
missing experiment.

## 7. Lossy variant-preserving compression

| Tool | Approach | Relation to G_CAPSUL |
|---|---|---|
| BFQzip (Guerrini, Louza and Rosone, arXiv:2304.08534) | Modifies bases and quality scores using the extended BWT so that an external caller still finds the variants | Poses a similar question, one structure for compression and variant evidence. It is lossy and prepares data for someone else's caller. G_CAPSUL is lossless and performs the calling itself. |

BFQzip was not benchmarked here, and the two serve different goals.

## 8. Combined capabilities

Whether any tool offers all three capabilities from one representation was checked against the
candidate architectures:

| Tool | Lossless compression | Reference-free calling | Locus retrieval | One representation |
|---|---|---|---|---|
| SPRING, PgRC, PgRC2, Genozip, NanoSpring | yes | no | no | not applicable |
| DiscoSNP++, Kmer2SNP | no (raw FASTQ input) | yes | no | not applicable |
| CRAM | yes | needs an external reference | needs an external reference | no |
| BEETL-fastq | yes | no, hands off to external alignment and an external reference | content match only | no |
| Population BWT | no, quality scores are not kept | at pre-specified sites only | content (k-mer) only | no |
| **G_CAPSUL** | yes | yes, from the archive alone | yes, from the archive alone | yes |

Evidence for the last row: a source scan of the encoder (`src/encoder.cpp`), the caller
(`include/caps_caller.h`) and the decoder (`src/decoder.cpp`) finds no calls to `system()`,
`popen()` or the `exec` family. The encoder uses `fork()` only to run its own code for parallel
candidate encoding. No reference genome file is opened by the core pipeline. The entropy coders
(LZMA, PPMd7, FSE/Huf0) are compiled in, and the vendored `fqzcomp` code is used for quality
scores only, never by the caller or by export, coverage and query. These claims about other tools
come from their papers, so they carry the evidence label Documented.

## 9. Statements that are not made

- No speed or size comparison against BEETL-fastq, CIndex or BFQzip.
- No claim that DiscoSNP++'s multi-allelic behaviour was verified in its source.
- No claim that PgRC2's pseudogenome splits alleles on real data. That is a prediction from its
  design, not a measurement. Its decoder was read and computes each read's position and strand
  during decompression, but exposes no query interface. Whether the position data persists or is
  transient was not determined.
- No claim that G_CAPSUL is the first reference-free FASTQ compressor. The claim made is
  narrower: no tool found combines lossless compression, reference-free variant calling and
  reference-free locus retrieval in one archive, to the best of the authors' knowledge.
