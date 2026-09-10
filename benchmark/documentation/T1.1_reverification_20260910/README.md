# T1.1 independently re-verified, 2026-09-10 — E. coli, byte-exact

`paper/PAPER.md` tells a reviewer to run this first:

    CLAIMS=1 bash benchmark/scripts/sanity_archive_one.sh SRR2584863
    -> archive 68,429,027 B exactly, 3/3 LOSSLESS

This directory is the raw evidence that it does. The run was a **clean rebuild
from current source**, independent of the 2026-09-09 sweep that produced the
published tables, and it reproduced the published `archive_bytes` exactly:

    dataset     tool     raw_bytes  archive_bytes  ratio_pct  compress_s  decompress_s  peak_ram_kb  lossless
    SRR2584863  CAPSULE  695163748  68429027       9.8436     11.68       5.82          2014084      LOSSLESS
    SRR2584863  SPRING   695163748  74076160       10.6559    9.67        5.14          2184896      LOSSLESS
    SRR2584863  Genozip  695163748  119674636      17.2153    1.89        0.83          1063900      LOSSLESS

All three tools decode-and-compare clean.

## The competitor columns did not reproduce, and ours did

Against the published table:

| tool | published 2026-09-09 | this re-run 2026-09-10 | delta |
|---|---|---|---|
| **CAPSULE** | 68,429,027 | 68,429,027 | **byte-exact** |
| SPRING | 74,086,400 | 74,076,160 | −0.014% |
| Genozip | 114,589,268 | 119,674,636 | +4.4% |

The Genozip gap was checked rather than assumed. Three further runs of the same
command on the same file gave 119,615,499 / 119,662,147 / 119,631,835 — a 0.04%
spread, so the 4.4% is not run-to-run noise. Thread count is not the cause
either (`-@ 12/8/4/2` → 119.618/119.661/119.639/119.639 MB). `--vblock` is:
16/64/128/512 give 120.0/112.4/99.8/86.2 MB on this file. Genozip's default
VBlock is selected from machine conditions, so its default output is not a
fixed target across runs on differently-loaded machines.

**This is not documented by Genozip.** Its paper mentions `--vblock` only as a
compression/speed tradeoff and makes no statement about reproducible output;
its online compression guide says the default is "selected based on
characteristics of the data" and likewise says nothing about reproducibility.

Both competitor figures were produced with default flags, exactly as ours were
with zero flags — default against default, which is the comparison the paper
describes. The published Genozip number is the one that favours Genozip, so
Claim 1's margin is understated, not inflated. See `paper/LIMITATIONS.md` §10. The CAPSULE check decodes the
**archive**, not the encoder's intermediate dumps — the distinction that hid
four silent data-loss bugs earlier in this project (CLAUDE.md §6.1).

| file | what it is |
|---|---|
| `sanity.log` | the full run, timestamped |
| `claim1_T1.1_T1.2.csv` | the row it produced — compare against `benchmark/results/claim1_T1.1_T1.2.csv` |
| `_manifest.tsv` | input file identity (size + checksum) |
| `SRR2584863.capsule.log` | the encoder's own stage-by-stage log |
| `encode.stdout` | encoder stdout |

## Why this directory exists

The working directory this came from (`results/sanity_SRR2584863_20260910_022835/`)
is **834 MB of archives and decoded FASTQ and is gitignored** — it does not ship
in a clone. `RESULTS_INVENTORY.md` called it "the useful one" while it was
invisible to anyone who cloned the repo. The text evidence, 20 KB of it, is
copied here so the claim travels with the repository.
