# c_star_pg_advance — CAPSULE

An independent, from-scratch pseudogenome-based FASTQ archive, developed as
a sandbox alongside ARCS. **As of 2026-09-03 this covers all three claims,
updated from this README's earlier, now-stale scope statement** (see
`docs/TECHNICAL_ARCHITECTURE.md` §0 for the full architecture):

- **COMPACT** — sequence, read order, names, and quality, all independently
  toggleable (§ "What you get, by configuration" below).
- **FAITHFUL** — reference-free heterozygous SNV/indel calling as an
  in-process side effect of compression (`include/caps_caller.h`).
- **ADDRESSABLE** — `export`/`coverage`/`query` served directly from the
  archive, no reference genome, no full decompression
  (`stages/capsule_decode.cpp`).

**Start here for a 3-line quick check of all three:**
```bash
bash scripts/run_capsule.sh 1   # COMPACT
bash scripts/run_capsule.sh 2   # FAITHFUL
bash scripts/run_capsule.sh 3   # ADDRESSABLE
```
Full command reference: `docs/COMMANDS_REFERENCE.md`. Fresh-server
bootstrap (every dataset and tool, exact verified commands):
`docs/SERVER_SETUP_AND_DOWNLOADS.md`.

## What you get, by configuration

**Every feature below is additive and independently OFF by default.**
Turning one off does not affect the others' correctness — it only removes
that column from the output. This matters because "compress a FASTQ" is
not one fixed operation here; it's a sequence-only archive unless you ask
for more.

| you set | what's IN the archive | what decoding gives you | what's NOT there if you don't set it |
|---|---|---|---|
| *(nothing — the default)* | sequence + read order only | `capsule_decode <archive> <outdir> <outdir>/reads.seq` → one sequence per line, in original file order | no read names, no quality scores, no `+` line — this is NOT a FASTQ, it's the sequence column only |
| `CAPS_NAMES=1` | + names/read-ID column | decode also writes `<outreads>.names` | quality still absent unless also set |
| `CAPS_QUAL=1` | + quality scores | decode also writes `<outreads>.qual` | names still absent unless also set |
| `CAPS_NAMES=1 CAPS_QUAL=1` | full FASTQ content | decode writes sequence + `.names` + `.qual`; a full 4-line FASTQ is reassembled from these three plus the recovered line-3 mode, verified byte-identical (same MD5) to the original input | nothing — this is the complete round trip |
| `CAPS_CALL=1` | *(no archive change — writes a VCF as a side effect)* | `$CALL_VCF` gets heterozygous SNV/indel calls in contig coordinates | does not affect what the archive itself contains; combinable with any of the above with no interaction (verified this session, see `docs/FINAL_ALGORITHMIC_SCAN.md`) |

**Claim 3's three operations are independent reads of the SAME archive —
running one does not run or require the others:**

| command | what you get | what you do NOT get |
|---|---|---|
| `capsule_decode export <arc> out.fa` | the assembled pseudogenome as FASTA — an assembly, not individual reads | no per-read output, no depth information, no coordinate lookup |
| `capsule_decode coverage <arc> out.tsv` | a `#region start end depth` TSV — per-base depth only | no sequence content anywhere in the output |
| `capsule_decode query <arc> out.fa START-END` | FASTA records for reads overlapping that coordinate range only | reads outside the range; also does NOT give you `.names`/`.qual` for those reads — query returns sequence only |
| `capsule_decode <arc> <outdir> <outdir>/reads.fq` (no mode keyword) | the full round trip — every original read, in order, full sequence | this is the only one of the four that reconstructs actual per-original-read output; the other three are all reads FROM the archive's own internal state (pg content or placement index), not from replaying every read |

**In short: pick exactly what you need.** Want to check compression ratio
only → default, no flags. Want the full file back → both `CAPS_*` flags on
at encode time. Want just the assembly → `export`, nothing else runs. Want
depth only → `coverage`, which is deliberately cheaper because it skips
rebuilding the assembly entirely (`TECHNICAL_ARCHITECTURE.md` §10.1).

## Where things are

| path | contents |
|---|---|
| `stages/` | the full 01 → 106 progression, one file per experiment |
| `stages/106_inprocess.cpp` | **the shipped encoder** — single process, in-memory streams |
| `include/` | `coders_inproc.h`, `coders_pgrc.h`, `seqpar_core.h` |
| `scripts/` | `build106.sh`, `encode_adaptive.sh`, `decode_105.py` |
| `thirdparty/` | PPMd7 (LZMA SDK, public domain), FSE/Huf0 (Yann Collet, BSD) |
| `docs/` | analysis and plans, including the layer-by-layer breakdown |
| `results/phase_a/` | raw measurement output, including reverted work |
| `DATASET_LOCKED.md` | the locked accessions — do not substitute |

## Build and run

```bash
scripts/build106.sh /tmp/best106
INPUT=reads.fq ARCHIVE=out.arc BEST=/tmp/best106 bash scripts/encode_adaptive.sh
```

## Verify losslessness

One command: encodes, decodes the streams back, and compares against the
original file's sequence column. Not a coder-level round trip -- the reference
is the input FASTQ itself.

```bash
scripts/verify_lossless.sh reads.fq
# LOSSLESS  archive=2618332 bytes  input=reads.fq
```

Confirmed on H. salinarum (2,618,332), S. acidocaldarius (3,192,790),
E. coli (7,965,683) and SARS-CoV-2 (838,654) -- the last being variable-length,
which PgRC2 declines outright.

## Current standing

Size, against PgRC2's own binary on the same inputs — **+3.19% aggregate, 7 of
8 wins**, every archive verified byte-identical after decoding back to the
original FASTQ:

| dataset | kingdom | ours | PgRC2 | |
|---|---|---|---|---|
| H. salinarum | Archaea | 2,618,332 | 3,050,477 | +14.17% |
| E. coli | Bacteria | 7,965,683 | 8,864,420 | +10.14% |
| P. falciparum | Protista | 16,321,262 | 17,219,695 | +5.22% |
| S. aureus | Bacteria | 13,241,852 | 13,595,003 | +2.60% |
| L. major | Protista | 27,561,086 | 28,272,652 | +2.52% |
| P. aeruginosa | Bacteria | 8,969,147 | 9,043,181 | +0.82% |
| S. cerevisiae | Fungi | 21,787,190 | 21,841,286 | +0.25% |
| S. acidocaldarius | Archaea | 3,192,790 | 3,114,782 | **−2.50%** |

Two datasets PgRC2 cannot archive at all, measured against SPRING and Genozip
instead:

| dataset | ours | SPRING | Genozip | PgRC2 |
|---|---|---|---|---|
| SARS-CoV-2 | 838,654 | 1,556,480 | 976,083 | refuses: variable-length reads |
| A. fumigatus | 34,526,106 | 74,813,440 | 147,388,530 | crashes at constant 150 bp |

Variable-length support and >255 bp are documented limits in PgRC2's README, so
those are scope, not defects. A. fumigatus is different: it satisfies every
constraint they state and still aborts, while a 100k-read subsample succeeds.

## Speed and memory

Phase A was four restructurings held to one rule — **no output byte may change**.
Three passed, one was reverted.

| | before | after |
|---|---|---|
| 7-file wall time | 243.78 s | **130.60 s** (−46.4%) |
| worst-case peak RSS | 1388 MB | **930 MB** (−33.0%) |
| adaptive path (E. coli) | 41.14 s | **32.87 s** (−21.1%) |

Still behind PgRC2 on both: roughly 1.7× slower and 2.5× heavier at worst.
The remaining speed gap is the coding phase, whose wall time cannot fall below
its longest single job — ours is 4.76 s against their 352 ms, because we run 13
coarse jobs where they run 82 fine ones, and because `best_encode` probes seven
coders where they estimate.

## What is deliberately recorded as failed

`results/phase_a/04_gate_A2_REVERTED.txt` keeps a change that made things worse
(+3.2%), because the reason matters: it estimated coder cost from source bytes,
and bytes do not predict time — `literal` is the largest stream and takes
0.43 s, `orig2uid` is smaller and takes 4.76 s. Several conclusions in `docs/`
are likewise marked retracted rather than deleted.

## Third-party

PPMd7 from the LZMA SDK (public domain) and FSE/Huf0 by Yann Collet (BSD) are
vendored under `thirdparty/`. PgRC2 itself is **not** vendored; it is cloned
separately for benchmarking:

```bash
git clone https://github.com/kowallus/PgRC.git
```
