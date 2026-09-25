# G_CAPSUL — reference-free FASTQ compressor, variant caller, and
# locus-addressable archive. Two-stage build: the first stage builds
# best106/capsule_decode with the full toolchain, the second stage ships
# only the built binaries and their runtime library, so the final image is
# a few MB, not a full compiler toolchain.
#
# Build:
#   docker build -t g_capsul .
#
# Run (single command, compress a FASTQ mounted from the host -- sequence +
# read order only, matching the README's own default/documented behavior):
#   docker run --rm -v "$PWD":/data -e ARCHIVE=/data/out.capsule g_capsul best106 /data/reads.fq
#
# For a full byte-identical round trip (sequence + names + quality + line 3),
# set CAPS_NAMES/CAPS_QUAL, exactly as the README's "full FASTQ" example does:
#   docker run --rm -v "$PWD":/data -e CAPS_NAMES=1 -e CAPS_QUAL=1 -e ARCHIVE=/data/out.capsule \
#     g_capsul best106 /data/reads.fq
#
# Decompress (writes out.capsule.names/.qual alongside it when those streams exist):
#   docker run --rm -v "$PWD":/data g_capsul capsule_decode /data/out.capsule /data/outdir /data/outdir/reads.fq
#
# Any other mode (export/coverage/query/call/index) works the same way --
# see docs/PROJECT_SCOPE.md for every command.
#
# Verified this exact round trip end-to-end before shipping (and re-verified
# on every push -- see .github/workflows/ci.yml's docker job): compressed and
# decompressed a real FASTQ through these two commands, reconstructed the
# 4-line records from the separate seq/names/qual outputs the same way
# scripts/benchmark_1_run.sh's own losscmp does, and confirmed byte-identical.

FROM ubuntu:24.04 AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake liblzma-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY CMakeLists.txt ./
COPY src/ src/
COPY include/ include/
COPY thirdparty/ thirdparty/
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build -j"$(nproc)"

FROM ubuntu:24.04
RUN apt-get update && apt-get install -y --no-install-recommends \
    liblzma5 libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/build/best106 /usr/local/bin/best106
COPY --from=build /src/build/capsule_decode /usr/local/bin/capsule_decode
COPY THIRD_PARTY_NOTICES.md /usr/local/share/g_capsul/THIRD_PARTY_NOTICES.md
COPY LICENSE /usr/local/share/g_capsul/LICENSE
# These three are set by scripts/encode_adaptive.sh (the documented entry
# point) on every real invocation, including every one behind the published
# results -- without them, best106 silently skips the streams the decoder's
# plain-reconstruction path needs (orig2uid, read_lengths), producing a
# complete-looking archive that decodes to zero reads with no warning.
# Baked in here as container-wide defaults so a direct `best106 <file>`
# single command is correct by construction, not dependent on the caller
# knowing this. Still overridable (docker run -e DUMP_LIT=0 ...) if ever
# needed.
ENV DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1
WORKDIR /data
ENTRYPOINT []
CMD ["best106", "--version"]
