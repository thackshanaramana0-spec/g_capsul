# SETUP_SHORT — fresh server to first result

Copy-paste path. Ubuntu 24.04, root, ~250 GB disk, ~32 GB RAM.
If any step fails, stop and open `SETUP_FULL.md` at the matching section —
every failure listed there has already happened here once.

---

## 0. tmux first

    tmux new -s work

Everything below outlives a dropped SSH session only if it is inside tmux.
Downloads run for hours; the benchmark runs for five.

## 1. Repos

    git clone https://github.com/thackshanaramana0-spec/ARCS.git /root/arcs-clean
    git clone https://github.com/thackshanaramana0-spec/g_capsul.git \
              /root/arcs-clean/c_star_pg_advance
    cd /root/arcs-clean/c_star_pg_advance

## 2. apt packages

    apt-get update -y
    apt-get install -y cmake build-essential git curl wget unzip nasm \
        samtools bcftools tabix bwa mosdepth seqtk kmc default-jre

`kmc` is required by Kmer2SNP (Claim 2's third arm) and is easy to miss —
without it that arm silently records NOT_AVAILABLE.

## 3. `/dev/stdout` must exist

    [ -e /dev/stdout ] || ln -sfn /proc/self/fd/1 /dev/stdout

Genozip tests for curl by calling `file_exists("/dev/stdout")`. If it is
missing, Genozip compresses your whole file and then **exits 1 with no
archive**, blanking an entire competitor column. This actually happened.

## 4. AWS CLI (apt has no candidate)

    cd /root && curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip
    unzip -q awscliv2.zip && ./aws/install && aws --version

## 5. SRA toolkit — SYMLINK, never copy

    cd /root
    wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzf sratoolkit.current-ubuntu64.tar.gz
    ln -sf /root/sratoolkit.*/bin/prefetch     /usr/local/bin/prefetch
    ln -sf /root/sratoolkit.*/bin/fasterq-dump /usr/local/bin/fasterq-dump
    prefetch --version && fasterq-dump --version

Copying the binary breaks it: they are symlink chains that resolve a sibling
`*-orig` binary at runtime.

## 6. SPRING (needs a GCC 13 patch)

    cd /root && git clone https://github.com/shubhamchandak94/SPRING.git
    cd SPRING && git submodule update --init --recursive
    # apply the <cstdint> patch -- see SETUP_FULL 3.6 for the exact file list
    bash /root/arcs-clean/c_star_pg_advance/server/patch_spring.sh
    mkdir -p build && cd build && cmake .. && make -j$(nproc)
    cp spring /usr/local/bin/ && spring --help

## 7. Genozip (build, then activate INTERACTIVELY)

    cd /root
    curl -sL https://github.com/divonlan/genozip/archive/refs/tags/genozip-15.0.87.tar.gz \
         -o genozip-src.tar.gz
    tar -xzf genozip-src.tar.gz && cd genozip-genozip-15.0.87
    make -j$(nproc) && cp genozip genounzip genocat genols /usr/local/bin/
    genozip --activate            # INTERACTIVE. Do this yourself in a terminal.

Then verify with a probe **larger than 1 MB** — Genozip only attempts its
telemetry upload above roughly that size, so a small test passes on a machine
where every real dataset fails:

    head -12000 /data/fastq/SRR2584863_1.fq > /tmp/p.fq
    genozip --force -o /tmp/p.genozip /tmp/p.fq && echo GENOZIP OK

## 8. rtg-tools (Claim 2 scoring)

    curl -s https://api.github.com/repos/RealTimeGenomics/rtg-tools/releases/latest \
      | grep browser_download_url        # take the *-linux-x64.zip asset
    cd /root && curl -sL <that-url> -o rtg-tools.zip && unzip -q rtg-tools.zip
    cp -r rtg-tools-3.13 /usr/local/rtg-tools
    ln -sf /usr/local/rtg-tools/rtg /usr/local/bin/rtg && rtg version

## 9. Competitor callers and assembler

    # DiscoSNP++  -> /root/DiscoSnp        (must be on PATH at run time)
    # Kmer2SNP    -> /root/Kmer2SNP        (+ conda env with networkx)
    # SPAdes      -> /root/SPAdes-4.0.0-Linux
    # PgRC2       -> /root/arcs-clean/method_c   (cmake build)
    export PATH="$HOME/DiscoSnp:$PATH"     # benchmark_1 does this itself

See `TOOLS.md` for the exact provenance of each.

## 10. Data (~55 GB)

    bash server/fetch_data.sh              # or benchmark/download.sh in ARCS

`SRR065390` (11 GB) and `SRR10676752` (15 GB) are the slow ones. Pull them from
the S3 open-data mirror rather than NCBI prefetch:

    aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR065390/SRR065390 \
        /data/fastq/prefetch/SRR065390/SRR065390.sra

Also required for Claim 2: `~/refs/chr20.fa` (GRCh37, contig named `20`),
`~/refs/chr20.sdf`, and `~/giab_truth/` (GIAB v4.2.1 VCF + BED for HG002-HG005).

## 11. Build and verify

    cd /root/arcs-clean/c_star_pg_advance
    scripts/build106.sh      /tmp/capsule_bin/best106          # encoder
    scripts/build_decode.sh  /tmp/capsule_bin/capsule_decode   # decoder
    scripts/run_tests.sh                                       # must be 15/15
    bash scripts/benchmark_0_preflight.sh                      # must print GO

`benchmark_0` re-builds both binaries itself, checks every tool, dataset,
reference and truth set, and re-derives the validated k-mer count. **If it does
not print `VERDICT: GO`, do not run the benchmark** — fix what it names.

## 12. Run

    bash scripts/benchmark_1_run.sh          # full sweep, ~5 h, 19 datasets
    # or, one dataset end to end first:
    CLAIMS=1 bash scripts/sanity_archive_one.sh SRR2584863   # ~45 s

The sanity run must produce an archive of exactly **68,429,027 B** and report
LOSSLESS for all three tools. If the size differs on any machine, that is a
real regression, not an environment difference — see
`../benchmark/documentation/REPRODUCIBILITY.md`.
