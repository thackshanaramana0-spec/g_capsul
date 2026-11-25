# Benchmark datasets: 15 organisms and 4 human individuals

This file lists the 19 datasets of the benchmark with their accessions, organisms and coverage,
and states which datasets belong to which comparison. The 19 datasets give 57 archives (three
tools each), and all 57 were verified lossless. The accessions are public (NCBI SRA, EBI ENA and
DDBJ), and the download steps are in [`SERVER_SETUP_AND_DOWNLOADS.md`](SERVER_SETUP_AND_DOWNLOADS.md).
The archive sizes are in `results/claim1/claim1_T1.1_T1.2.csv`, and the result file is
authoritative where this file differs.

## How the 19 divide

| Group | Datasets | Used for |
|---|---|---|
| 15 organisms | The table below | The SPRING and Genozip comparison (T1.1, T1.2) and Claim 3 |
| 4 human individuals | HG002, HG003, HG004, HG005 (GIAB, chromosome 20, 30x) | Claim 2 (variant calling), and also carried in T1.1 and T1.2 because the same archive serves both claims |

The four human sets are not part of the 15-organism comparison and should not be counted into it.
The published Claim 1 table therefore has 19 rows.

## The 15 organisms

Coverage is for the first mate file (`_1`) only, which is the file the benchmark uses.

| # | Accession | Organism | Kingdom | Coverage |
|---|-----------|----------|---------|----------|
| 1 | SRR2584863 | E. coli B REL606 | Bacteria | about 59x |
| 2 | ERR552797 | M. tuberculosis H37Rv | Bacteria | about 49x |
| 3 | SRR554369 | P. aeruginosa PAO1 | Bacteria | about 24x |
| 4 | ERR5181310 | SARS-CoV-2 | Virus | about 500x |
| 5 | ERR17740259 | S. aureus | Bacteria | about 158x |
| 6 | DRR976266 | S. cerevisiae | Fungi | about 63x |
| 7 | SRR36741279 | Leishmania major | Protista | about 22x |
| 8 | SRR37283774 | P. falciparum | Protista | about 12.6x |
| 9 | SRR32429602 | Human betaherpesvirus 5 (HCMV) | Virus | not recorded |
| 10 | SRR39257532 | Aspergillus fumigatus | Fungi | about 20x |
| 11 | SRR29296997 | Halobacterium salinarum | Archaea | about 54x |
| 12 | ERR12954017 | Sulfolobus acidocaldarius | Archaea | about 116x |
| 13 | SRR40271341 | Helicobacter pylori | Bacteria | about 139x |
| 14 | SRR065390 | C. elegans N2 WGS | Animalia | about 34x |
| 15 | SRR10676752 | Utricularia gibba | Plantae | about 65x |

By kingdom, the set has 5 bacteria, 2 viruses, 2 fungi, 2 protists, 2 archaea, 1 plant and 1
animal.

Two entries need a note:

- **Utricularia gibba (SRR10676752)** has a genome of 100,688,548 bp (assembly U_gibba_v2), the
  smallest known sequenced plant genome, about 4.3 times smaller than T. cacao. This makes a
  full-coverage plant run possible at a moderate size. The run is Illumina MiSeq, paired, WGS, from
  BioProject PRJNA595351 (study SRP237333). It has 26,658,699 spots with mates of about 246 bp,
  which gives about 65x for `_1`. The `_1` file is about 3.5 GB.
- **C. elegans N2 (SRR065390)** has 33,808,546 reads of 100 bp (11,282,985,734 bytes for `_1`) and
  about 34x coverage. It is a large single-pass input, so memory and time are proportionally higher
  than for the small datasets.

### Datasets considered and not used

Drosophila melanogaster (SRR40104876 and the low-coverage SRR40104920) was dropped, and Utricularia
gibba took its place so that a plant is represented. Campylobacter jejuni (SRR40402583) was dropped
to make room for C. elegans, because bacteria were over-represented (6 of 15) and no animal was
included. Its earlier result row is superseded by the 19-dataset benchmark.

### Low-coverage control

Arabidopsis thaliana ERR17716639 (2.8x) is a deliberate low-coverage control outside the 15. On it
G_CAPSUL's archive is 2.74% larger than SPRING's, and this is reported as a limitation of the
method at low coverage.

## The four human individuals

All four were standardized to 30x on chromosome 20, so that the SNV and indel comparisons are not
confounded by the source coverage of each GIAB individual (60x to 300x). They were streamed from the
GIAB S3 whole-genome BAMs with samtools and downsampled.

| # | File | Individual | Population | Depth | Additional role |
|---|---|---|---|---|---|
| 16 | `HG002_pooled.fq` | NA24385 | Ashkenazi son | 30x chr20 | Sole subject of the T2.2 coverage sweep (10x, 15x, 20x, 30x) |
| 17 | `HG003_pooled.fq` | NA24149 | Ashkenazi father | 30x chr20 | One half of the synthetic tetraploid of T2.5 |
| 18 | `HG004_pooled.fq` | NA24143 | Ashkenazi mother | 30x chr20 | The other half of the T2.5 tetraploid |
| 19 | `HG005_pooled.fq` | NA24631 | Han Chinese son | 30x chr20 | The one individual on which G_CAPSUL loses the SNV and indel comparison, see [`../limitations_and_scope.md`](../limitations_and_scope.md) |

All four appear in T2.1 (heterozygous SNVs), T2.3 (heterozygous indels), T3.4 (exact-match
recall) and T3.5 (locus retrieval). T2.4 (multi-allelic sites) uses HG002 only.

HG005 is the only variable-length set in the benchmark, with 250 bp trimmed reads and 216 distinct
lengths, against a fixed 148 bp for the other three. This is the identified cause of its loss in
T2.1 and T2.3. The loss is reported and the reads were not truncated to match the others.

**Truth sets and reference.** The truth sets are the GIAB v4.2.1 VCF and high-confidence BED of each
individual, for chromosome 20 only. The GRCh37 chromosome 20 sequence (`chr20.fa`, with the
chromosome named "20") was used to build probes and for liftover. The caller uses none of these.
G_CAPSUL calls from the archive without a reference, and the truth sets and the reference serve
only to score the calls.

## Where the files go

`scripts/benchmark_scripts/benchmark_1_run.sh` reads the FASTQ files from `$CAPSULE_DATA_DIR`,
which defaults to `/data/fastq/`, as `<accession>_1.fq`. Substituting a dataset invalidates the
comparison with the published tables.

## Results for the two datasets added last

These two were added after the first benchmark run and were run in full. The values are from
`results/claim1/claim1_T1.1_T1.2.csv`, in bytes.

| Accession | G_CAPSUL | SPRING | Genozip | Lossless |
|---|---|---|---|---|
| SRR10676752 (U. gibba) | 1,650,034,057 | 1,718,548,480 | 2,954,900,932 | 3 of 3 |
| SRR065390 (C. elegans) | 887,410,845 | 942,643,200 | 1,597,083,547 | 3 of 3 |
