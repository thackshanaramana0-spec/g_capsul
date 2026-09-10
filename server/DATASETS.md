# DATASETS — the locked set, sizes, and where each file goes

**19 datasets, 55.0 GB of FASTQ.** Plus references and truth sets (~2.9 GB).

The authoritative list is `NEW_DATASET_LOCKED.md` in this repository — **not**
the 17-accession list in ARCS's `DATASET_LOCKED.md`, which is an earlier lock.
One swap was made 2026-09-02 (Drosophila out, *Utricularia gibba* in); that
file records why.

## Claim 1 — all 19 (compression)

| accession | organism | bytes | GB |
|---|---|---:|---:|
| `SRR2584863` | E.coli B REL606 | 695,163,748 | 0.65 |
| `ERR552797` | M.tuberculosis H37Rv | 431,078,636 | 0.40 |
| `SRR554369` | P.aeruginosa PAO1 | 456,443,722 | 0.43 |
| `ERR5181310` | SARS-CoV-2 | 471,183,048 | 0.44 |
| `ERR17740259` | S.aureus | 1,336,957,140 | 1.25 |
| `DRR976266` | S.cerevisiae | 2,248,042,490 | 2.09 |
| `SRR36741279` | Leishmania major | 1,661,157,034 | 1.55 |
| `SRR37283774` | P.falciparum | 1,019,818,902 | 0.95 |
| `SRR32429602` | HCMV | 1,868,917,402 | 1.74 |
| `SRR39257532` | Aspergillus fumigatus | 1,761,853,698 | 1.64 |
| `SRR29296997` | Halobacterium salinarum | 210,225,246 | 0.20 |
| `ERR12954017` | Sulfolobus acidocaldarius | 333,144,048 | 0.31 |
| `SRR40271341` | Helicobacter pylori | 290,894,264 | 0.27 |
| `SRR065390` | C.elegans N2 | 11,282,985,734 | 10.51 |
| `SRR10676752` | Utricularia gibba | 15,372,502,354 | 14.32 |
| `HG002` | GIAB Ashkenazi son | 4,279,197,941 | 3.99 |
| `HG003` | GIAB Ashkenazi father | 4,279,569,332 | 3.99 |
| `HG004` | GIAB Ashkenazi mother | 4,279,327,517 | 3.99 |
| `HG005` | GIAB Han Chinese son | 6,809,023,884 | 6.34 |
| **total** | | **59,087,486,140** | **55.0** |

Expected location: `/data/fastq/<ACC>_1.fq`, and `/data/fastq/<IND>_pooled.fq`
for the four GIAB individuals.

## Claim 2 — the four GIAB individuals only

`HG002` `HG003` `HG004` `HG005`, chr20 at a standardised 30x. Sourced from the
GIAB S3 WGS BAMs (chr20 streamed via samtools) and downsampled, because the S3
source coverage varies per individual (60-300x) and standardising removes that
confound.

**HG005 is not like the other three, and it matters.** It is a different
sequencing run: instrument `D00360`, **250 bp quality-trimmed reads with 216
distinct lengths**, against `HISEQ1` and **148 bp fixed** for HG002/3/4. That
single fact explains its lower variant-calling precision — see
`../benchmark/documentation/HG005_EXPLAINED.md`. Do not treat the four as
interchangeable.

The other 15 are bacterial, archaeal, viral or fungal — **haploid**. Claim 2
and T3.4 are meaningless on them (no heterozygous sites), which is a scope
statement, not a gap.

## References and truth sets

| what | path | note |
|---|---|---|
| GRCh37 chr20 | `~/refs/chr20.fa` | contig must be named `20`, not `chr20` |
| chr20 SDF | `~/refs/chr20.sdf` | for `rtg vcfeval`; built with `rtg format` |
| BWA index | `~/refs/chr20.fa.bwt` etc. | for T3.2's baseline |
| GIAB truth | `~/giab_truth/<IND>_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` | + `.tbi` |
| GIAB regions | `~/giab_truth/<IND>_..._benchmark_noinconsistent.bed` | falls back to `..._benchmark.bed` |

## Fetching

    bash /root/arcs-clean/benchmark/download.sh /data/fastq 2>&1 | tee /data/download.log
    tail -30 /data/download.log     # last line must contain DOWNLOAD COMPLETE

**The two big SRA files should come from the S3 open-data mirror**, which is
far faster than NCBI prefetch:

    aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR065390/SRR065390 \
        /data/fastq/prefetch/SRR065390/SRR065390.sra
    aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR10676752/SRR10676752 \
        /data/fastq/prefetch/SRR10676752/SRR10676752.sra

then re-run the download script — it skips prefetch and goes straight to
`fasterq-dump`.

## Verifying what you fetched

    bash scripts/benchmark_0_preflight.sh    # section 4 lists all 19 with sizes

Sizes must match the table above. A short file is usually an interrupted
`fasterq-dump`; delete and refetch rather than resuming.

## Banned accessions — never substitute these in

`SRR390728` `SRR988075` `SRR327342` `SRR1663585` `SRR1296601` `ERR015526`
`SRR1294122` `ERR174310` `SRR16357346` `SRR1945765`

They appear in the literature but are excluded here (wrong platform, wrong
scale, or already used as a tuning set). Substituting a dataset invalidates the
comparison against the locked numbers.
