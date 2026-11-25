# Server setup and dataset downloads

This runbook lists the tools and data that the benchmark needs and the commands that install or
fetch each of them, so that a new machine can reach the same state. The commands come from the
scripts used for the reported runs and were checked against the server that produced them. Version
numbers are the ones installed there. Where something could not be verified, the text says so.
References to "the development notes" and to the per-tool benchmark notes mean the author's
internal working record, which is not distributed with the repository.

The accession list is in [`DATASETS.md`](DATASETS.md). The sections are independent of each other.

---

## 0. Base system tools

```bash
apt-get update
apt-get install -y build-essential cmake liblzma-dev bwa samtools tabix \
    mosdepth bcftools megahit python3 python3-pip curl unzip git
```

Verified exact versions currently on this server (`dpkg -l`):

| package | version |
|---|---|
| bwa | 0.7.17-7 |
| samtools | 1.19.2-1build2 |
| tabix | 1.19+ds-1.1build3 |
| mosdepth | 0.3.6+ds-1 |
| bcftools | 1.19-1build2 |
| megahit | 1.2.9-5 |
| liblzma-dev | 5.6.1+really5.4.5-1ubuntu0.3 |

## 1. SRA download tooling (Claim 1 datasets)

```bash
# SRA Toolkit (prefetch, fasterq-dump) -- verified 3.4.1 on this server
# Install via the official SRA Toolkit release tarball for your platform:
# https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
# (already present at /usr/local/bin/prefetch, /usr/local/bin/fasterq-dump)

# AWS CLI v2 -- verified aws-cli/2.36.30 on this server, no credentials
# needed (every S3 access in this project uses --no-sign-request against
# public buckets)
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip && sudo ./aws/install
```

## 2. Claim 1: SRA datasets

**The exact accessions** are listed in [`DATASETS.md`](DATASETS.md) (15 organisms and 4 human
individuals). The set must not be changed without recording the reason, because a substitution
invalidates the comparison with the published tables.

**The actual, working download command** (from the dataset download script,
run for every accession in the list):

```bash
DATA_DIR=/data/fastq
prefetch --output-directory "$DATA_DIR/prefetch" "$ACCESSION"
fasterq-dump --threads "$(nproc)" \
    -O "$DATA_DIR" "$DATA_DIR/prefetch/$ACCESSION/$ACCESSION.sra"
# Default mode. For paired runs it writes ${ACCESSION}_1.fastq and ${ACCESSION}_2.fastq, and
# writes reads that have no mate to ${ACCESSION}.fastq. Keep ${ACCESSION}_1.fastq as _1.fq.
# Do not add --split-files: it keeps reads without a mate in the first-mate file, which then
# no longer matches the recorded sizes (ERR5181310: 472,726,004 bytes with --split-files,
# 471,183,048 bytes in default mode, which is the recorded size).
```

**For the two large datasets (C. elegans SRR065390 ~22GB, T. cacao
SRR870667 ~15GB SRA), `prefetch` is too slow — use the S3 Open Data mirror
instead**, as recorded in the development notes:

```bash
aws s3 cp --no-sign-request \
    s3://sra-pub-run-odp/sra/SRR065390/SRR065390 \
    /data/fastq/prefetch/SRR065390/SRR065390.sra
aws s3 cp --no-sign-request \
    s3://sra-pub-run-odp/sra/SRR870667/SRR870667 \
    /data/fastq/prefetch/SRR870667/SRR870667.sra
# then run fasterq-dump on the downloaded .sra exactly as above --
# it skips prefetch automatically when the .sra file already exists
```

**This is the actual, measured-faster method for any large SRA accession**,
not a fallback: `s3://sra-pub-run-odp/sra/<ACCESSION>/<ACCESSION>` is a
public, unsigned S3 object for every SRA run, and `aws s3 cp` saturates
bandwidth far better than `prefetch`'s own transport for multi-GB files.
**Rule of thumb used on this project: any accession whose `.sra` is
expected above ~5GB, use the S3 mirror first, not `prefetch`.**

**Fetching the whole list.** Loop the commands above over the accessions of `DATASETS.md`, using
the S3 mirror for any accession whose `.sra` is above about 5 GB. Check that every
`${ACCESSION}_1.fastq` exists and has the expected size before starting the benchmark, then rename
it to `${ACCESSION}_1.fq` in `$CAPSULE_DATA_DIR` (default `/data/fastq`).

## 3. Claim 2: GIAB HG002 to HG005, chromosome 20 at 30x

**This is a targeted S3 stream and downsample, not a bulk download**, because
the full WGS BAMs are hundreds of GB per individual and only chr20 at 30×
is needed. This is the download step used for the reported runs:

```bash
export HTS_S3_ADDRESS_STYLE=path
CHR20_LEN=63025520; TARGET_COV=30; READ_LEN=150
N_NEEDED=$(( TARGET_COV * CHR20_LEN / READ_LEN ))    # ~12.6M reads

# 1. Auto-discover the largest (highest-coverage) WGS BAM for the individual
BASE="s3://giab/data/AshkenazimTrio/HG002_NA24385_son"   # swap per individual, see table below
BAM_KEY=$(aws s3 ls --no-sign-request "$BASE/" --recursive \
    | grep '\.bam$' | grep -v '\.bai' | grep -iv 'indel\|sv\|snv\|phased' \
    | sort -k3 -n -r | head -1 | awk '{print $4}')
S3_BAM="s3://giab/${BAM_KEY}"

# 2. Stream ONLY chr20 reads (random-access via the BAM index, not a full download)
samtools view -b "$S3_BAM" chr20 \
    | samtools sort -n -@ "$(nproc)" \
    | samtools fastq -o HG002_chr20_full.fq -

# 3. Downsample to exactly 30x (seqtk, or the python fallback in download.sh
#    if seqtk isn't installed)
FRAC=$(awk "BEGIN{printf \"%.6f\", $N_NEEDED / $(wc -l < HG002_chr20_full.fq | awk '{print $1/4}')}")
seqtk sample -s 42 HG002_chr20_full.fq "$FRAC" > HG002_pooled.fq
```

| individual | S3 base path |
|---|---|
| HG002 | `s3://giab/data/AshkenazimTrio/HG002_NA24385_son` |
| HG003 | `s3://giab/data/AshkenazimTrio/HG003_NA24149_father` |
| HG004 | `s3://giab/data/AshkenazimTrio/HG004_NA24143_mother` |
| HG005 | `s3://giab/data/ChineseTrio/HG005_NA24631_son` |

Source GIAB S3 coverage varies per individual (60-300×) — this is exactly
why the downsample-to-30× step exists: it removes that confound so that F1
comparisons across individuals do not depend on which individual was sequenced deeper.

**chr20 reference and GIAB truth VCFs** :

```bash
# GRCh37/hg19 chr20, sequence renamed from ">chr20" to ">20" to match GIAB's convention
curl -L -o chr20.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz"
zcat chr20.fa.gz | sed 's/^>chr/>/' > ~/refs/chr20.fa

# GIAB v4.2.1 truth VCF + confident BED, per individual (~700MB VCF each)
FTP=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
curl -L -o HG002_truth.vcf.gz \
    "$FTP/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
tabix -p vcf HG002_truth.vcf.gz
curl -L -o HG002_confident.bed \
    "$FTP/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
# HG003/HG004/HG005: same pattern, swap the individual's directory name
# (HG005 uses ChineseTrio/HG005_NA24631_son and its BED has no "_noinconsistent" suffix)
```

**rtg vcfeval** (the third-party scorer used for every F1 number in Claim 2):

```bash
# Installed on this server at /root/rtg-tools-3.13/, from rtg-tools.zip
# Official releases: https://github.com/RealTimeGenomics/rtg-tools/releases
curl -L -o rtg-tools.zip \
    "https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.13/rtg-tools-3.13-linux-x64.zip"
unzip rtg-tools.zip
```

## 4. Claim 3: E. coli export, coverage and query

Uses the same `SRR2584863` accession as Claim 1 (dataset 1 of `DATASETS.md`), so no separate
download is needed. `scripts/benchmark_scripts/benchmark_1_run.sh` builds the other binaries it
needs, and SPAdes if it is absent (see section 5).

## 5. Comparison tools: versions and install methods

| tool | version (verified) | install method (verified against this server) |
|---|---|---|
| **SPRING** | (no `--version` flag; boost-based build) | `git clone https://github.com/shubhamchandak94/SPRING.git` at `/root/SPRING`, built per its own README (CMake + boost). |
| **Genozip** | 15.0.87 (`genozip --version`) | Source tarball `genozip-src.tar.gz` (7.6MB, present at `/root/genozip-src.tar.gz`) extracted to `/root/genozip-genozip-15.0.87/` and built there. Official source: `https://github.com/genozip/genozip` (download a release tarball, e.g. `https://github.com/genozip/genozip/archive/refs/tags/genozip-15.0.87.tar.gz`, then `make`). |
| **PgRC2** | git clone, GPL-3, never vendored | `git clone https://github.com/kowallus/PgRC.git` at `method_c/` in the working directory (renamed on clone). Build per its own CMake instructions. |
| **DiscoSNP++** | v2.6.2-12 (per the DiscoSNP++ comparison note) | `git clone https://github.com/GATB/DiscoSnp.git` at `/root/DiscoSnp`. Build with its own `./INSTALL` or CMake script; needs `-G <ref.fa>` at run time (mandatory, not optional — see the development notes' command reference). |
| **Kmer2SNP** | — | `git clone https://github.com/yanboANU/Kmer2SNP.git` at `/root/Kmer2SNP`. **Real compatibility fixes required** for this 2020-era code on a modern Python (documented in full in the Kmer2SNP benchmark note): `time.clock()` → `time.perf_counter()` (Python 3.8 removed the former), `networkx.connected_component_subgraphs` → a documented generator-expression equivalent (removed in networkx 2.4+), and `findGSE` (an R package) bypassed via Kmer2SNP's own documented `--c1/--c2/--r` flags with data-derived values from a k-mer histogram, rather than installing R. |
| **KMC** | 3.2.4 | `apt-get install kmc`. The k-mer counter used for Kmer2SNP in the reported results, in place of DSK (see `scripts/claim2_scripts/run_kmer2snp.sh`). |
| **DSK** | v2.3.3 | Not used in the reported results. Kmer2SNP's own driver script hardcodes a `/path2dsk/` placeholder, so KMC replaces it. |
| **MEGAHIT** | 1.2.9-5 | `apt-get install megahit` — a real Ubuntu package, no manual build needed. |
| **SPAdes** | v4.0.0 | **Prebuilt binary, not apt** (no current Ubuntu package tracks the latest release): `curl -sL https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz -o spades.tar.gz && tar xzf spades.tar.gz` — installed this way at `~/SPAdes-4.0.0-Linux/`, verified working on the E. coli assembly of T3.1 (239 contigs). |
| **bwa / samtools / mosdepth / tabix / bcftools** | see §0 table | plain `apt-get install`, no manual build for any of these. |

## 6. Mistakes to avoid

1. **Do not use `prefetch` alone for any SRA accession likely to exceed
   ~5GB.** It works but is far slower than the S3 mirror (§2) for large
   files — this was a measured, not assumed, speed difference.
2. **Do not install R/`findGSE` to run Kmer2SNP.** The tool's own
   `--c1/--c2/--r` flags accept the same inputs `findGSE` would have
   computed — use those instead of adding a whole R toolchain for one
   dependency (the Kmer2SNP benchmark note).
3. **Do not run DiscoSNP++ without `-G <ref.fa>`.** Without it, its output
   cannot be lifted to genome coordinates for scoring — this is mandatory
   for every benchmark in this project, not an optional accuracy flag.
4. **Do not skip `-g` on SPRING decompress.** `spring -d` without `-g`
   produces the wrong output format — the development notes' own command reference
   flags this explicitly because it was a real, previously-hit mistake.
5. **Do not download full WGS BAMs for GIAB individuals.** Stream chr20
   only via `samtools view -b <S3 BAM> chr20` (§3) — a full WGS BAM per
   individual is hundreds of GB and none of it outside chr20 is used.
6. **Do not build the encoder without `-fopenmp`.** `scripts/build_scripts/build106.sh`'s
   own header comment documents this exact mistake: without the flag, the
   single `#pragma omp parallel for` in the largest stage is silently
   discarded and it runs serial — a real, previously-hit, -45%-wall-time
   mistake, which is why the build command is now a checked-in script
   rather than something retyped from memory.
