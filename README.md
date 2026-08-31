# c_star_pg_advance

An independent, from-scratch reimplementation of PgRC2's pseudogenome-based read
compressor, developed as a sandbox alongside ARCS. Scope is the **sequence and
read-order** layers only: assembly, mapping, self-matching, and the streams a
decoder needs to rebuild every read in its original position. Name and quality
coders are present in the stage history (61-75, 84-85, 92) but are not part of
the current scope.

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
