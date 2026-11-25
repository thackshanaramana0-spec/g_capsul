# Capabilities and input contract

This document describes what the program does and exactly which input it accepts, refuses or
leaves unsupported. Every limit below was taken from the source (`src/encoder.cpp`,
`src/decoder.cpp`, `include/`). If this document and the source ever disagree, the source is
correct. Results and their limits are in [`limitations_and_scope.md`](limitations_and_scope.md).

## 1. What the program does

G_CAPSUL compresses short-read FASTQ into a lossless archive. While compressing, it assembles the
reads into a pseudogenome and keeps that structure in the archive. The archive can then be queried
for variants, coverage and sequence without a reference genome and without decompressing it
completely first.

## 2. Capabilities

| Capability | What it provides | How to use it |
|---|---|---|
| Lossless compression | Byte-identical reconstruction of the sequence, and optionally read order, read names, quality scores and the line-3 field, from the archive alone | `encode_adaptive.sh`, `capsule_decode <archive> <outdir> [reads]` |
| Read names | Names stored and restored | `CAPS_NAMES=1` |
| Quality scores | Stored with the vendored fqzcomp coder | `CAPS_QUAL=1` |
| Variant calling | Heterozygous SNVs, indels and multi-allelic sites, ploidy 2 to 4, with no reference genome | `CAPS_CALL=1` at compression and `capsule_decode call`, with `CAPS_PLOIDY` |
| Export | The pseudogenome as FASTA | `capsule_decode export` |
| Coverage | Per-base depth from the retained placements | `capsule_decode coverage` |
| Locus retrieval | Every read covering a coordinate range, probed upstream and downstream | `capsule_decode query` |
| Content retrieval | Every read containing a query sequence, through an optional FM-index | `CAPS_FMI=1` at `index` time and `CAPS_QUERY_CONTAIN=1` at query time |
| Parameter selection | Candidate settings are tried and the smallest archive is kept | `scripts/build_scripts/encode_adaptive.sh` |
| Input pre-check | A report on whether a file is supported, with numbers taken from the encoder source | `scripts/benchmark_scripts/check_input_scope.py` |

The names, quality and line-3 options are independent switches that cost nothing when unset. The
default archive holds the sequence and read order only. See the Configuration table in the README.

## 3. Platforms

| Platform | Status |
|---|---|
| Linux | Supported. Built and tested in CI. |
| macOS | Supported through CMake with Homebrew GCC, because Apple Clang has no OpenMP. Built and tested in CI. |
| Docker | `ghcr.io/thackshanaramana0-spec/g_capsul:1.2.0`, publicly pullable. CI rebuilds the image and runs a compress, decompress and `cmp` round trip on every push, and again against the published image on every release. |
| Windows | Not supported natively, because the encoder and decoder use POSIX process calls (`fork`, `waitpid`). Use WSL or Docker. |

Specific versions are installable from the tagged
[releases](https://github.com/thackshanaramana0-spec/g_capsul/releases). Each release keeps
`CAPS_VERSION` in the source and `CITATION.cff` in step.

## 4. Input contract

| Property | Accepted | Behaviour outside it |
|---|---|---|
| File format | Plain-text FASTQ with four lines per record | Malformed or truncated records are refused |
| Compression | Uncompressed only, no `.gz` or `.bz2` | Not detected specially, and usually ends in the "no usable records" refusal |
| Read length | 1 to 1023 bases (`MAX_READ_LEN`) | Refused above 1023 |
| Read count | Up to 536,870,911 reads (2²⁹ − 1, the limit of the 29-bit read id in the MEM stage) | Refused above this |
| Sequence alphabet | Upper-case `A`, `C`, `G`, `T`, `N` | Any other byte is refused. `N` is substituted internally and its positions are stored separately. |
| N bases | No `N` after base 256 of a read, and at most 255 `N` per read | Refused |
| Header length | Up to 1023 characters when `CAPS_NAMES=1` | Refused above this. Unbounded when names are not stored. |
| Line endings | LF | CRLF is refused |
| Quality line | Same length as the sequence | Any mismatch is refused |
| Line 3 | A bare `+` in every record, or `+` followed by the header in every record | Mixed content is refused when names are stored |
| Read layout | Single-end, fixed or variable length within one file | Mate relationships are not kept |
| Validated read lengths | 40 to 301 bp | Reads up to 1023 bp are processed, but compression ratio and calling accuracy are not validated beyond 301 bp |

## 5. Input that is refused

The encoder prints a `FATAL:` message naming the first offending record and exits without writing
an archive. Refused cases:

- reads longer than 1023 bases, which includes Oxford Nanopore and PacBio data (this is a fixed
  per-read buffer size and not a policy threshold)
- more than 536,870,911 usable records
- a read header of 1023 or more characters when `CAPS_NAMES=1`
- a record whose quality line differs in length from its sequence line
- an empty file, or a file with no usable record
- CRLF line endings
- any sequence byte other than upper-case `A`, `C`, `G`, `T`, `N` (including lower case and IUPAC
  codes), an `N` after base 256, or more than 255 `N` in one read
- truncated or malformed records
- mixed line-3 content when `CAPS_NAMES=1`
- no arguments, or an input file that does not exist or cannot be opened

## 6. Input that is not supported and not specially detected

Most of these are still stopped by the general record and alphabet checks.

- **Gzip-compressed FASTQ.** The file must be decompressed first.
- **FASTA, SAM, BAM, CRAM and other formats.** They are not FASTQ and are not detected as a special
  case.
- **Protein sequence and RNA with `U`.** Refused by the alphabet check.
- **Paired-end data.** There is no R1 and R2 linking or interleaving. Each input file is treated
  as an independent stream of reads, so paired data can be encoded as two files, but the format
  does not record which reads are mates.
- **Numeric command-line parameters.** `MAXMM`, `MINOV` and similar positional arguments are
  parsed with `atoi` and are not range-checked. An invalid value is not refused, and the behaviour
  is not validated.
- **An unrecognized decoder mode.** Anything other than `export`, `coverage`, `query`, `index` or
  `call` is treated as a file path and fails with a generic "bad archive" error. The exit status is
  nonzero and no output is written, but the message does not name the mistake.
- **Whole human genomes.** Full 3.1 Gb datasets were outside the evaluation, which covers
  chromosome 20 regions and the 19-dataset collection. Beyond the read-count ceiling above, this is
  a statement about validation and not a code-level block.
