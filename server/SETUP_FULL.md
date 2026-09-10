# SETUP_FULL — every step, and every failure that has actually happened

`SETUP_SHORT.md` is the happy path. This is the same sequence with the reasons,
plus a catalogue of real failures with their exact error text so you can search
for the string you are staring at.

Merged from `arcs-clean/NEW_SERVER_MANDATORY.md` (the original handover, whose
tool-install gotchas are preserved verbatim below because they are still true)
and everything learned in the sessions since — Kmer2SNP, KMC, the Genozip
`/dev/stdout` failure, the 19-dataset lock, the Capsule build, and the
reproducibility work.

---

## 1. tmux before anything

    tmux new -s work
    # detach: Ctrl-b d      reattach: tmux attach -t work

Downloads take hours and the benchmark takes five. A dropped SSH session
without tmux loses the run and, worse, can leave a half-written archive that a
later step will happily read.

## 2. Repositories

    git clone https://github.com/thackshanaramana0-spec/ARCS.git /root/arcs-clean
    git clone https://github.com/thackshanaramana0-spec/g_capsul.git \
              /root/arcs-clean/c_star_pg_advance

They are separate repositories that nest; ARCS `.gitignore:73` excludes
`c_star_pg_advance/`. Committing in one never touches the other.

**Branches.** The repo carries two branch names, `c_star_pg_advance` and
`gpt2026`, and **they point at the same commit** — the second was created
mid-project and the work simply continued on it. Take either; `git log` is
identical. (The untracked `gpt2026_*` source files in the working tree are a
concurrent agent's experiment, are on no branch, and are not part of this
result.)

## 3. Toolchain

### 3.1 apt

    apt-get update -y
    apt-get install -y cmake build-essential git curl wget unzip nasm \
        samtools bcftools tabix bwa mosdepth seqtk kmc default-jre

`nasm` is needed by Genozip's build. `kmc` is needed by Kmer2SNP and is the
easiest to forget — without it Claim 2's third arm records `NOT_AVAILABLE` and
T2.1 quietly becomes a two-tool table.

### 3.2 `/dev/stdout` — the failure that blanks a whole column

    [ -e /dev/stdout ] || ln -sfn /proc/self/fd/1 /dev/stdout

**Symptom:** `genozip: LICENSE ERROR in zfile_compress_genozip_header:83 ...
Failed to upload a telemetry record to the Genozip server.` and **no archive**,
after it has already compressed the entire input.

**Cause**, from Genozip's own source (`arch.c:552`, version 15.0.87):

    return !system("which curl > /dev/null 2>&1") && file_exists("/dev/stdout");

On this box `/dev/stdout` had gone missing while `/dev/stdin` and `/dev/stderr`
survived. Both `curl_available()` and `wget_available()` therefore returned
false, Genozip concluded it had no way to upload the telemetry its Student
licence requires, and refused to write the archive header. `strace` showed
**zero AF_INET syscalls** — it never even attempted the network.

**Two things that make this nasty:** `command -v genozip` still passes, so a
presence check does not catch it; and Genozip only attempts the upload above
roughly 1 MB of input, so a small smoke test passes on a machine where every
real dataset fails. `benchmark_0_preflight.sh` therefore asserts `/dev/stdout`
exists *and* round-trips a >1 MB probe.

### 3.3 AWS CLI

`apt-get install -y awscli` has **no installation candidate** on this server's
repo configuration. Use the official installer:

    cd /root
    curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip
    unzip -q awscliv2.zip && ./aws/install && aws --version

### 3.4 SRA Toolkit — symlink, do not copy

    cd /root
    wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzf sratoolkit.current-ubuntu64.tar.gz
    ln -sf /root/sratoolkit.*/bin/prefetch     /usr/local/bin/prefetch
    ln -sf /root/sratoolkit.*/bin/fasterq-dump /usr/local/bin/fasterq-dump
    prefetch --version && fasterq-dump --version

**Symptom if you copy instead:** `failed to exec prefetch: No such file or
directory`.

**Cause:** they are symlink chains (`prefetch → prefetch.3 → prefetch.3.4.1 →
sratools.3.4.1`) that expect a sibling `<tool>-orig.3.4.1` **in the same
directory** at runtime. Flattening breaks the sibling lookup.

### 3.5 SPRING — required source patch under GCC 13

    cd /root && git clone https://github.com/shubhamchandak94/SPRING.git
    cd SPRING && git submodule update --init --recursive
    bash /root/arcs-clean/c_star_pg_advance/server/patch_spring.sh
    mkdir -p build && cd build && cmake .. && make -j$(nproc)
    cp spring /usr/local/bin/

**Symptom without the patch:** `'uint8_t' does not name a type`, cascading into
dozens of "no member named" errors.

**Cause:** GCC 13 no longer transitively includes `<cstdint>`. The patch script
inserts `#include <cstdint>` before the first include of every **compiled**
source/header. `old_src/` is not built and does not need it.
`-Walloc-size-larger-than=` warnings during the build are benign.

**A separate SPRING trap, at RUN time, not build time:** `-g` means *"the input
is gzipped"*. Our inputs are plain `.fq`.

    spring -c -i in.fq -o out.spring -t $(nproc)      # correct
    spring -d -i out.spring -o dec.fq -t $(nproc)     # correct

Passing `-g` on decompress writes **gzip bytes into a `.fq`**, which then reads
as garbage and reports a false LOSSY. That mistake has been made twice in this
project, and once put a false `LOSSY` into a results table.

### 3.6 Genozip — build, then activate interactively

    cd /root
    curl -sL https://github.com/divonlan/genozip/archive/refs/tags/genozip-15.0.87.tar.gz \
         -o genozip-src.tar.gz
    tar -xzf genozip-src.tar.gz && cd genozip-genozip-15.0.87
    make -j$(nproc) && cp genozip genounzip genocat genols /usr/local/bin/
    genozip --activate

No prebuilt Linux binaries exist on GitHub releases — source only.

**Activation is mandatory and interactive** (license key → email confirmation →
T&C acceptance → emailed 6-digit code). It cannot be reliably driven through an
agent sandbox; a prior attempt via a Python `pty` wrapper was blocked partway
through by a safety classifier. **Run it yourself in a terminal.** A free
Genozip Student licence (genozip.com/student) is valid about a year.

**Genozip's archive is not byte-reproducible** — it embeds run metadata and
drifts about 37 B (0.005%) between runs of the same input. Never gate on
byte-identity of a Genozip archive.

### 3.7 rtg-tools

    curl -s https://api.github.com/repos/RealTimeGenomics/rtg-tools/releases/latest \
      | grep browser_download_url

The `.../releases/latest/download/rtg-tools-non-commercial-linux-x64.zip` URL
**404s** — that asset name is not on current releases. Take the
`*-linux-x64.zip` asset from the API listing (3.13 at time of writing).

    cd /root && curl -sL <that-url> -o rtg-tools.zip && unzip -q rtg-tools.zip
    cp -r rtg-tools-3.13 /usr/local/rtg-tools
    ln -sf /usr/local/rtg-tools/rtg /usr/local/bin/rtg && rtg version

### 3.8 DiscoSNP++ — must be on PATH, not merely present

    export PATH="$HOME/DiscoSnp:$PATH"

`run_discoSnp++.sh` resolves helpers **through PATH**. A present-but-unreachable
install fails later as *"no SNV line"*, which reads like a scoring bug rather
than a missing tool. `benchmark_1_run.sh` exports the PATH itself; the
preflight checks reachability rather than file existence.

`-G <ref.fa>` is mandatory in every invocation. The POS off-by-one in its
output (`$2-1`) is corrected in our harness — without that correction it scores
F1 0.004 instead of 0.85, which looks like the tool being broken.

### 3.9 Kmer2SNP — bypass its own dependency wrappers

Kmer2SNP is at `/root/Kmer2SNP` with a conda env at
`/root/miniconda3/envs/kmer2snp_r`. Two traps, both hit:

**Trap 1 — its DSK wrapper is hardcoded to a path that does not exist.**
`script/run_dsk.sh` calls `/tmp/dsk-v2.3.3-bin-Linux/bin/dsk`. It fails silently
(`/path2dsk/dsk: not found`), leaves `hete.para` empty, and `kmer2snp.py` then
dies in findGSE parsing with `IndexError: list index out of range` — which
reads like a Kmer2SNP bug rather than a missing dependency.

**It never needs DSK or findGSE.** Both are only called to FILL IN arguments
(`kmer2snp.py:59-70`). Supplying `--t1` (the k-mer counts) and `--c1/--c2/--r`
(the coverage band) bypasses both. Our `scripts/run_kmer2snp.sh` counts with
KMC and derives the band from the histogram.

**Trap 2 — `networkx` lives in the conda env, not system python3.**
`python3 kmer2snp.py` dies on `ModuleNotFoundError: No module named 'networkx'`.
Use `/root/miniconda3/envs/kmer2snp_r/bin/python`.

**Trap 3 — the pairs file.** Kmer2SNP writes both `k_K_pair.snp` (SNP pairs)
and `k_K_pair.non` (NON-SNP pairs). Feeding `.non` yields **zero** calls,
because every entry differs at more than one base. That would publish
`Kmer2SNP F1 = 0.000` as if measured.

### 3.10 SPAdes and PgRC2

    # SPAdes 4.0.0 -> /root/SPAdes-4.0.0-Linux (Claim 3's T3.1 baseline)
    # PgRC2       -> /root/arcs-clean/method_c, cmake build -> build/PgRC

PgRC2 is GPL-3 and is **cloned separately, never vendored** into our source.
Note for interpreting results: `PgRC -h` offers compress, `-d` decompress and
tuning flags only — it has no query/export interface, and its decompressed
output is **bare DNA** (no names, no `+`, no quality). It is a DNA-stream
compressor, not a FASTQ compressor.

## 4. Data

See `DATASETS.md` for the full manifest, sizes and verification. Summary:
**19 datasets, ~55 GB of FASTQ**, plus `~/refs/chr20.fa` + `.sdf` and
`~/giab_truth/` (GIAB v4.2.1, ~2.8 GB).

The two large SRA files should come from the S3 open-data mirror, which is far
faster than NCBI prefetch:

    aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR065390/SRR065390 \
        /data/fastq/prefetch/SRR065390/SRR065390.sra
    aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR10676752/SRR10676752 \
        /data/fastq/prefetch/SRR10676752/SRR10676752.sra

then re-run the download script, which skips prefetch and proceeds to
`fasterq-dump`.

## 5. Build

    cd /root/arcs-clean/c_star_pg_advance
    scripts/build106.sh      /tmp/capsule_bin/best106
    scripts/build_decode.sh  /tmp/capsule_bin/capsule_decode

**`build106.sh` exists because the build command once lived only in shell
history — which is exactly how `-fopenmp` went missing.** Without that flag the
`#pragma omp parallel for` in the sweep is silently discarded and the largest
stage runs serial; linking it correctly was worth −45% wall time. Never build
by hand.

For ARCS (the outer repo):

    cd /root/arcs-clean && cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel $(nproc) && ctest --test-dir build   # 7/7

## 6. Disk, RAM and the ENOSPC trap

Peak transient during the sweep is about **36 GB** on top of the data; the
preflight refuses to start below **60 GB free**.

**This is not caution.** A full disk once silently truncated the caller's
k-mer spill file, which dropped k-mers and produced a **fake F1 regression**
(0.1722 against a true 0.2613) with no error anywhere. The distinct-k-mer count
came out 72,485,292 instead of 84,181,752. Every harness now asserts free space
*and* re-derives that k-mer count before trusting a result.

Watch the largest dataset: `SRR10676752` (15.4 GB) took free space to **14 GB**
mid-run during its decode. It passed LOSSLESS, but it is the tightest point.

## 7. Verify before running anything long

    scripts/run_tests.sh                       # 15/15 self-contained assertions
    bash scripts/benchmark_0_preflight.sh      # must print VERDICT: GO

`benchmark_0` rebuilds both binaries, checks every tool/dataset/reference/truth
set, verifies the archive-calling path end to end on a synthetic input,
re-derives the validated k-mer count and F1, prints the projected runtime per
phase, and gives a GO / NO-GO. **A NO-GO names the missing item.** Do not
proceed past it.

Then one dataset end to end:

    CLAIMS=1 bash scripts/sanity_archive_one.sh SRR2584863     # ~45 s

Expected: archive exactly **68,429,027 B**, ratio 9.844%, and LOSSLESS from all
three tools. That number is machine-independent — verified byte-identical
across 2/3/7/12 cores and every candidate-concurrency setting.

## 8. Catalogue of failures that have actually happened

| symptom | cause | fix |
|---|---|---|
| `failed to exec prefetch` | SRA binaries copied, not symlinked | §3.4 |
| `'uint8_t' does not name a type` | GCC 13 + unpatched SPRING | §3.5 |
| SPRING reports LOSSY | `-g` passed on decompress | §3.5 |
| Genozip: no archive, LICENSE ERROR | `/dev/stdout` missing | §3.2 |
| Genozip works on small input, fails on real data | telemetry only above ~1 MB | §3.2 |
| Kmer2SNP `IndexError: list index out of range` | its DSK wrapper's hardcoded path | §3.9 |
| Kmer2SNP `ModuleNotFoundError: networkx` | system python3 instead of the conda env | §3.9 |
| Kmer2SNP returns 0 calls | fed `k_K_pair.non` instead of `.snp` | §3.9 |
| DiscoSNP++ "no SNV line" | not reachable on PATH | §3.8 |
| DiscoSNP++ F1 = 0.004 | POS off-by-one uncorrected | §3.8 |
| rtg-tools download 404 | asset name not on current releases | §3.7 |
| Fake F1 regression, no error | disk full, truncated k-mer spill | §6 |
| Encoder much slower than expected | `-fopenmp` missing; built by hand | §5 |
| Compression 17× slower than SPRING | `CAPS_CALL=1` runs the caller inline; use `CAPS_SPANS=1` | benchmark docs |
| `ARCHIVE LACKS contig_spans` | archive built without `CAPS_SPANS=1` | benchmark docs |
| Archive smaller than expected, `reads written: 0` | encoder invoked bare, without `DUMP_PERM=1` | §9 |

## 9. Never invoke the encoder directly

    # WRONG -- produces an incomplete archive that fails later as a "code bug"
    /tmp/capsule_bin/best106 reads.fq ...

    # RIGHT
    INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/capsule_bin/best106 \
        bash scripts/encode_adaptive.sh

The encoder needs a specific environment and argument vector. Run bare, it
writes an archive missing `pos_abs`/`pos_strand`/`read_lengths` (they are gated
on `DUMP_PERM=1`), which then decodes with `reads written: 0` rather than
erroring. The tell is that the archive is a few MB *smaller* than a known-good
one for the same input.

**Five separate "the tool is broken" diagnoses in this project were all wrong
invocations.** Check the command before blaming the code.
