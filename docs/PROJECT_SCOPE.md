# G_CAPSUL — project scope: what it does, what it can do, and its input contract

This file is the code-level scope reference: what the tool actually does,
what it can be used for, and exactly what input it will and will not accept.
It is not about benchmark results or paper claims — for that, see
[`honest_limitations_and_scope.md`](honest_limitations_and_scope.md). Every
number here is taken directly from the current source
(`src/encoder.cpp`, `src/decoder.cpp`, `include/`), not from memory or the
paper — where the two could ever disagree, the source wins.

---

## 1. What it does, in one sentence

G_CAPSUL compresses short-read FASTQ into a compact, lossless archive by
assembling the reads into a de novo "pseudogenome" as it compresses, then
keeps that assembly around after compression so the archive itself can be
queried — for variants, coverage, and sequence — without a reference genome
and without fully decompressing first.

## 2. The three capabilities

- **COMPACT** — lossless FASTQ compression. Byte-exact round-trip of
  sequence, and optionally read order, names, quality, and the FASTQ line-3
  field, all reconstructible from the archive alone.
- **FAITHFUL** — reference-free variant calling directly from the archive
  (heterozygous SNVs, indels, multi-allelic sites), with no FASTQ
  decompression and no external reference genome.
- **ADDRESSABLE** — direct queries against the archive: export the
  assembled sequence, compute per-base depth/coverage, or retrieve every
  read covering an arbitrary coordinate range — all without a separate
  decompress-then-search step.

## 3. What it can do

- Losslessly compress single-end short-read FASTQ, reconstructing the
  original sequence byte-for-byte.
- Optionally preserve exact read order, read names/IDs (`CAPS_NAMES=1`),
  quality scores (`CAPS_QUAL=1`, vendored fqzcomp/htscodecs), and the FASTQ
  line-3 field, each independently toggleable and gated to zero cost when
  unset.
- Call heterozygous SNVs, indels, and multi-allelic variants with no
  reference genome (`CAPS_CALL=1`), for ploidy 2 through 4
  (`CAPS_PLOIDY`, default 2/diploid).
- Export the assembled pseudogenome as FASTA directly from the archive.
- Compute per-base sequencing depth/coverage directly from the archive.
- Retrieve every read covering a given coordinate range of the
  pseudogenome via bilateral (upstream + downstream) probing, with an
  optional post-archive completion index (`CAPS_QUERY_CONTAIN=1`) that
  trades archive-adjacent disk space and a higher homozygous
  false-positive rate for complete recall.
- Adaptively select its own compression parameters by sweeping candidate
  settings and keeping the smallest resulting archive
  (`scripts/encode_adaptive.sh`), rather than requiring the caller to tune
  them.
- Assess a FASTQ file's fitness for this pipeline before committing to a
  full encode, via a standalone pre-flight report
  (`scripts/check_input_scope.py`) that derives its numbers from the live
  encoder source rather than a hardcoded copy.
- Refuse cleanly and loudly on out-of-scope or malformed input (see §5),
  with a specific, actionable message rather than a crash or a silently
  wrong/empty archive.
- Build and run on Linux (native) and macOS (via CMake, real Homebrew GCC
  rather than AppleClang) — verified in CI on both. Native Windows is not
  currently supported (see §5).

## 4. Input this tool requires and accepts

| Property | Accepted range | Enforced? |
|---|---|---|
| File format | Plain-text FASTQ, standard 4 lines/record (`@header`, sequence, `+`/description, quality) | Malformed records refused |
| Compression | Uncompressed only — no `.gz`/`.bz2`/etc. | Not applicable, see §5 |
| Read length | 1–1023 bases per read (`MAX_READ_LEN`, `src/encoder.cpp`) | Hard refusal above this |
| Read count | Up to 536,870,911 reads per input (2²⁹−1, the MEM stage's 29-bit packed read-id ceiling) | Hard refusal above this |
| Sequence alphabet | Uppercase `A`/`C`/`G`/`T`; `N` is substituted internally and tracked separately | See §5 for what is *not* validated |
| Header length | Up to 1023 characters if `CAPS_NAMES=1` (`nmc::MAXTOK`); unbounded if unset | Hard refusal above this, `CAPS_NAMES=1` only |
| Seq/quality parity | Quality line length must equal sequence line length, per record | Hard refusal on any mismatch |
| Read shape | Single-end; fixed or variable length within one file | No built-in paired-end/mate awareness (see §5) |
| Validated read-length range | 40–301 bp (this project's locked benchmark datasets) | Reads up to 1023 bp are processed but unvalidated for compression ratio/calling accuracy past 301 bp |

## 5. What this tool explicitly cannot, or should not, be given

**Hard-refused (the encoder will exit with a `FATAL:` message, not process the file):**
- Long-read data (Oxford Nanopore, PacBio, or anything else over 1023 bases
  per read) — this is a fixed per-read stack buffer size, not a policy
  threshold.
- A file with more than 536,870,911 usable records.
- A read header longer than 1023 characters, if `CAPS_NAMES=1` is set.
- Any record whose quality-line length doesn't match its sequence-line
  length.
- An empty file, or a file with zero usable (well-formed) FASTQ records.
- Running the encoder with no arguments, or pointing it at a file that
  doesn't exist or can't be opened.

**Not supported, and not always cleanly rejected — real, disclosed gaps, not silent guarantees:**
- Gzip-compressed FASTQ (`.fastq.gz`) — must be decompressed to plain text
  first; the encoder does not detect or reject this specially, it will
  simply fail to parse valid FASTQ records from it and most likely hit the
  "zero usable records" refusal.
- FASTA (no quality line), SAM/BAM/CRAM, or any other non-FASTQ format —
  same as above, not specially detected.
- Protein sequence, RNA (with `U`), IUPAC ambiguity codes
  (`R`/`Y`/`S`/`W`/etc.), or lowercase bases (`a`/`c`/`g`/`t`) — the base
  packer only recognizes uppercase `A`/`C`/`G`/`T` as valid; anything else
  is treated the same as an encoding edge case internally rather than
  being explicitly detected and refused at the pre-flight stage. Do not
  trust a resulting archive on such input without decoding and diffing it
  against the source.
- Paired-end/mate-pair reads — there is no R1/R2 linking or interleaving
  awareness anywhere in the format; each input file is treated as one
  independent stream of reads. Paired data can be given to the tool (e.g.
  concatenated, or as two separate encodes), but the format itself does
  not know or care that two reads are mates.
- Numeric CLI tuning parameters (`MAXMM`, `MINOV`, and similar positional
  arguments) are parsed with `atoi` and are not range-checked — an invalid
  value (negative, non-numeric) will not be refused, and the resulting
  behavior is not validated for out-of-range inputs.
- An unrecognized mode string passed to the decoder (anything other than
  `export`, `coverage`, `query`, or `call`) does not produce an "unknown
  command" message — it falls through to the plain-decode path, which then
  fails with a generic "bad archive" error because the mode string is
  treated as a file path. The failure is still a clean refusal (confirmed:
  a nonzero exit, no crash, no wrong output), just not a message that names
  the actual mistake.
- Native Windows — the encoder and decoder use POSIX-only APIs
  (`fork()`/`waitpid()`, `sys/wait.h`) for their adaptive parameter-sweep
  and huge-page-hint code paths; these do not exist on native Windows.
  Linux and macOS (including Apple Silicon) are supported and CI-verified;
  Windows would require either WSL or a source-level port of that specific
  process-management code, neither of which is done today.
- Human whole-genome-scale datasets (full ~3.1 Gb) are outside this
  project's validated scope by an explicit prior decision — chr20-region
  and 19-dataset-sweep scale inputs are what has actually been measured.
  This is a validation-scope statement, not a hard code-level block beyond
  the read-count ceiling above.

---

*Every hard limit and every disclosed gap in this file was verified against
the current source during this session's own edge-case audit and
SPRING/PgRC2 comparison, and re-confirmed with a real adversarial test pass
against the built binaries (zero-arg, missing/empty/garbage/truncated input,
corrupted archives, malformed queries) — not asserted from documentation
alone.*
