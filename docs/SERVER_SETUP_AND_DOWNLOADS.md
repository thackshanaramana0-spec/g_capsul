# Server setup — every dataset and tool, exactly as installed, verified

Written 2026-09-03. Every command below is either (a) copied from a
checked-in script (`benchmark/download.sh` in the outer `/root/arcs-clean`
repo) and re-verified against this server's actual state, or (b)
reconstructed from a verified git remote / installed version / package
manager record on this exact server — never guessed. Where something could
not be independently verified, it says so explicitly rather than presenting
a guess as fact.

**Purpose:** if this server is destroyed, this file plus the two repos
(`arcs-clean`, `c_star_pg_advance`) is everything a new server needs to
reach the same state — same datasets, same tool versions, same install
method. Follow it top to bottom; each phase is independent.

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

## 2. Claim 1 — 10-17 locked SRA datasets

**Source of truth for the exact accessions: `/root/arcs-clean/DATASET_LOCKED.md`**
(10 primary + 7 extended = 17) and, for the G_CAPSUL-sandbox-specific
15-dataset SPRING/Genozip comparison, `/root/arcs-clean/c_star_pg_advance/NEW_DATASET_LOCKED.md`
(one swap: Drosophila out, Utricularia gibba in — documented there with the
exact reasoning, not silently substituted). **Never change these lists
without documenting why**, per this project's own standing rule.

**The actual, working download command** (from `benchmark/download.sh`,
run for every accession in the locked list):

```bash
DATA_DIR=/data/fastq
prefetch --output-directory "$DATA_DIR/prefetch" "$ACCESSION"
fasterq-dump --split-files --threads "$(nproc)" \
    -O "$DATA_DIR" "$DATA_DIR/prefetch/$ACCESSION/$ACCESSION.sra"
# fasterq-dump emits ${ACCESSION}.fastq for single-end; rename to _1.fq
```

**For the two large datasets (C. elegans SRR065390 ~22GB, T. cacao
SRR870667 ~15GB SRA), `prefetch` is too slow — use the S3 Open Data mirror
instead**, exactly as documented in `CLAUDE.md`:

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

**Run the whole locked list in one command** (already scripted, not manual):

```bash
bash benchmark/download.sh /data/fastq 2>&1 | tee /data/download.log
tail -30 /data/download.log   # must end with "DOWNLOAD COMPLETE"
```

**Utricularia gibba (Plantae, `SRR10676752`)** — the one dataset in the
15-set list not yet downloaded as of this writing. Per `NEW_DATASET_LOCKED.md`,
its `.sra` (7.3 GB) was already staged at `/tmp/newdl/SRR10676752/` via the
S3-mirror method above; only the `fasterq-dump` conversion step remains.

## 3. Claim 2 — GIAB HG002-HG005 chr20 at standardized 30×

**Not a bulk file download — a targeted S3 stream + downsample**, because
the full WGS BAMs are hundreds of GB per individual and only chr20 at 30×
is needed. This is `benchmark/download.sh` Phase 3, reproduced here exactly:

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
why the downsample-to-30× step exists: it removes that confound so T3 F1
comparisons across individuals are apples-to-apples, not an accident of
whichever individual happened to be sequenced deeper.

**chr20 reference and GIAB truth VCFs** (also Phase 5 of `download.sh`):

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

## 4. Claim 3 — E. coli export/coverage/query baseline

Uses the same `SRR2584863` accession already fetched for Claim 1 (locked
accession #1) — no separate download needed. `scripts/run_claim3.sh` in
`c_star_pg_advance` builds everything else it needs (its own binaries,
SPAdes if absent) automatically — see §5 below for the SPAdes part
specifically.

## 5. SOTA and allied tools — exact install method, verified per tool

| tool | version (verified) | install method (verified against this server) |
|---|---|---|
| **SPRING** | (no `--version` flag; boost-based build) | `git clone https://github.com/shubhamchandak94/SPRING.git` at `/root/SPRING`, built per its own README (CMake + boost). |
| **Genozip** | 15.0.87 (`genozip --version`) | Source tarball `genozip-src.tar.gz` (7.6MB, present at `/root/genozip-src.tar.gz`) extracted to `/root/genozip-genozip-15.0.87/` and built there. Official source: `https://github.com/genozip/genozip` (download a release tarball, e.g. `https://github.com/genozip/genozip/archive/refs/tags/genozip-15.0.87.tar.gz`, then `make`). |
| **PgRC2** | git clone, GPL-3, never vendored | `git clone https://github.com/kowallus/PgRC.git` at `/root/arcs-clean/method_c` (renamed on clone). Build per its own CMake instructions. |
| **DiscoSNP++** | v2.6.2-12 (per `docs/HOW_DISCOSNP_WINS.md`) | `git clone https://github.com/GATB/DiscoSnp.git` at `/root/DiscoSnp`. Build with its own `./INSTALL` or CMake script; needs `-G <ref.fa>` at run time (mandatory, not optional — see `CLAUDE.md`'s command reference). |
| **Kmer2SNP** | — | `git clone https://github.com/yanboANU/Kmer2SNP.git` at `/root/Kmer2SNP`. **Real compatibility fixes required** for this 2020-era code on a modern Python (documented in full in `docs/KMER2SNP_BENCHMARK.md`): `time.clock()` → `time.perf_counter()` (Python 3.8 removed the former), `networkx.connected_component_subgraphs` → a documented generator-expression equivalent (removed in networkx 2.4+), and `findGSE` (an R package) bypassed via Kmer2SNP's own documented `--c1/--c2/--r` flags with data-derived values from a k-mer histogram, rather than installing R. |
| **DSK** | v2.3.3 | `git clone https://github.com/GATB/dsk.git` at `/root/dsk`. Kmer2SNP's driver script hardcodes a `/path2dsk/` placeholder — must be pointed at the real built binary after cloning. |
| **MEGAHIT** | 1.2.9-5 | `apt-get install megahit` — a real Ubuntu package, no manual build needed. |
| **SPAdes** | v4.0.0 | **Prebuilt binary, not apt** (no current Ubuntu package tracks the latest release): `curl -sL https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz -o spades.tar.gz && tar xzf spades.tar.gz` — installed this way in this session at `~/SPAdes-4.0.0-Linux/`, verified working (real 239-contig E. coli assembly, `docs/CLAIM3_LOCKED.md` §6.2). |
| **bwa / samtools / mosdepth / tabix / bcftools** | see §0 table | plain `apt-get install`, no manual build for any of these. |

## 6. What NOT to do — real errors this project hit, so a fresh server doesn't repeat them

1. **Do not use `prefetch` alone for any SRA accession likely to exceed
   ~5GB.** It works but is far slower than the S3 mirror (§2) for large
   files — this was a measured, not assumed, speed difference.
2. **Do not install R/`findGSE` to run Kmer2SNP.** The tool's own
   `--c1/--c2/--r` flags accept the same inputs `findGSE` would have
   computed — use those instead of adding a whole R toolchain for one
   dependency (`docs/KMER2SNP_BENCHMARK.md`).
3. **Do not run DiscoSNP++ without `-G <ref.fa>`.** Without it, its output
   cannot be lifted to genome coordinates for scoring — this is mandatory
   for every benchmark in this project, not an optional accuracy flag.
4. **Do not skip `-g` on SPRING decompress.** `spring -d` without `-g`
   produces the wrong output format — `CLAUDE.md`'s own command reference
   flags this explicitly because it was a real, previously-hit mistake.
5. **Do not download full WGS BAMs for GIAB individuals.** Stream chr20
   only via `samtools view -b <S3 BAM> chr20` (§3) — a full WGS BAM per
   individual is hundreds of GB and none of it outside chr20 is used.
6. **Do not build the encoder without `-fopenmp`.** `scripts/build106.sh`'s
   own header comment documents this exact mistake: without the flag, the
   single `#pragma omp parallel for` in the largest stage is silently
   discarded and it runs serial — a real, previously-hit, -45%-wall-time
   mistake, which is why the build command is now a checked-in script
   rather than something retyped from memory.
