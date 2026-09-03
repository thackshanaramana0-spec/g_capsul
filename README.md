# c_star_pg_advance — CAPSULE

An independent, from-scratch pseudogenome-based FASTQ archive, developed as
a sandbox alongside ARCS. One archive format, evaluated as three claims:

- **COMPACT** — competitive lossless compression (sequence, read order,
  names, quality — all independently toggleable).
- **FAITHFUL** — reference-free heterozygous SNV/indel calling as an
  in-process side effect of compression, no separate assembler or
  alignment step.
- **ADDRESSABLE** — `export`/`coverage`/`query` served directly from the
  archive: no reference genome, no full decompression, because the
  compressor already built the assembly and placement index these
  operations need.

Full architecture: `docs/TECHNICAL_ARCHITECTURE.md`. Authoritative status
for each claim, including what's locked and what's still open:
`docs/CLAIM1_FINAL_VERDICT.md`, `docs/CLAIM2_FINAL_VERDICT.md`,
`docs/CLAIM3_FINAL_VERDICT.md`.

## Quick start

```bash
bash scripts/run_capsule.sh 1   # COMPACT    — verify losslessness on the locked E. coli dataset
bash scripts/run_capsule.sh 2   # FAITHFUL   — fast synthetic caller regression test
bash scripts/run_capsule.sh 3   # ADDRESSABLE — fast synthetic decoder regression test
```

Each claim is independent — there's no need to run them in order. Every
command builds whatever binaries it needs on first use. Full command
reference, including the real (slower, real-data) benchmarks behind each
claim's headline numbers: `docs/COMMANDS_REFERENCE.md`. Setting up a fresh
server from nothing (every dataset and tool, exact verified commands, not
guessed): `docs/SERVER_SETUP_AND_DOWNLOADS.md`.

## What you get, by configuration

**Every feature below is additive and OFF by default.** Turning one off
does not affect the others' correctness — it only removes that column
from the output. "Compress a FASTQ" is not one fixed operation here: the
default is a sequence-only archive, not a full FASTQ, unless you ask for
more.

| you set | what's IN the archive | what decoding gives you | what's NOT there if you don't set it |
|---|---|---|---|
| *(nothing — the default)* | sequence + read order only | `capsule_decode <archive> <outdir> <outdir>/reads.seq` → one sequence per line, in original file order | no read names, no quality scores, no `+` line — this is NOT a FASTQ, it's the sequence column only |
| `CAPS_NAMES=1` | + names/read-ID column | decode also writes `<outreads>.names` | quality still absent unless also set |
| `CAPS_QUAL=1` | + quality scores | decode also writes `<outreads>.qual` | names still absent unless also set |
| `CAPS_NAMES=1 CAPS_QUAL=1` | full FASTQ content | sequence + `.names` + `.qual`; a full 4-line FASTQ is reassembled from these plus the recovered line-3 mode, verified byte-identical (same MD5) to the original input | nothing — this is the complete round trip |
| `CAPS_CALL=1` | *(no archive change — writes a VCF as a side effect)* | `$CALL_VCF` gets heterozygous SNV/indel calls in contig coordinates | does not affect archive contents; combinable with any of the above (verified, `docs/FINAL_ALGORITHMIC_SCAN.md`) |

Claim 3's three operations are independent reads of the **same** archive —
running one does not run or require the others:

| command | what you get | what you do NOT get |
|---|---|---|
| `capsule_decode export <arc> out.fa` | the assembled pseudogenome as FASTA — an assembly, not individual reads | no per-read output, no depth, no coordinate lookup |
| `capsule_decode coverage <arc> out.tsv` | a `#region start end depth` TSV — per-base depth only, computed without rebuilding the assembly | no sequence content anywhere in the output |
| `capsule_decode query <arc> out.fa START-END` | FASTA for reads overlapping that coordinate range only | reads outside the range; also no `.names`/`.qual` — query returns sequence only |
| `capsule_decode <arc> <outdir> <outdir>/reads.fq` (no mode keyword) | the full round trip — every original read, full sequence, in order | this is the only mode that reconstructs actual per-original-read output; the other three read the archive's internal state directly |

## Where things are

| path | contents |
|---|---|
| `stages/106_inprocess.cpp` | **the shipped encoder** — assembly, mapping, stream coding, and (gated on env vars) names/quality/calling, all in one process |
| `stages/capsule_decode.cpp` | **the shipped decoder** — full round trip plus the three Claim 3 modes |
| `include/caps_caller.h` | the Claim 2 variant caller (dual-substrate SNV pileup + bubble/indel extraction + positional clustering) |
| `include/coders_inproc.h`, `coders_pgrc.h`, `seqpar_core.h`, `names_coder.h`, `quality_coder.h` | the stream-specific coders |
| `stages/01…106` | the full experimental progression — one file per decision, several later contradicted by measurement and kept, not deleted |
| `scripts/run_capsule.sh` | single entry point for all three claims (§ Quick start) |
| `scripts/capsule_config.sh` | the one file to edit if dataset paths move — nothing else hardcodes a path |
| `scripts/build106.sh`, `build_decode.sh` | build the two binaries |
| `scripts/encode_adaptive.sh`, `verify_lossless.sh` | Claim 1's compress-and-verify path |
| `scripts/run_giab_indel_capsule.sh`, `run_window_bench_capsule.sh`, `run_polyploid_bench_capsule.sh` | Claim 2's real-GIAB benchmarks |
| `scripts/run_claim3.sh` | Claim 3's one-command export/coverage/query benchmark |
| `scripts/test_claim2.sh`, `test_claim3.sh` | fast synthetic regression tests (no downloads needed) |
| `thirdparty/` | PPMd7 (LZMA SDK, public domain), FSE/Huf0 (Yann Collet, BSD — `thirdparty/fse/LICENSE`), htscodecs/fqzcomp (BSD 3-clause) |
| `docs/` | architecture, per-claim results, checklists, and this project's own record of every refuted idea |
| `results/phase_a/` | raw measurement CSVs, including reverted work |
| `DATASET_LOCKED.md`, `NEW_DATASET_LOCKED.md` | the locked accessions — do not substitute without documenting why |

## Current standing

### COMPACT

Whole-file archive (sequence + names + quality + line 3) against SPRING and
Genozip on the same 14 real datasets, every archive decoded back to
byte-identical input before being counted:

**14/14 wins vs SPRING, +11.68% aggregate. 14/14 wins vs Genozip, +48.93% aggregate.**

| dataset | organism | ours | SPRING | Genozip |
|---|---|---:|---:|---:|
| SRR2584863 | E. coli B REL606 | 68,677,977 | 74,045,440 | 119,618,198 |
| ERR552797 | M. tuberculosis H37Rv | 46,964,185 | 52,101,120 | 82,799,524 |
| SRR554369 | P. aeruginosa PAO1 | 57,320,645 | 59,166,720 | 88,813,037 |
| ERR5181310 | SARS-CoV-2 | 8,686,774 | 9,390,080 | 9,218,949 |
| ERR17740259 | S. aureus | 83,421,422 | 92,303,360 | 164,126,889 |
| DRR976266 | S. cerevisiae | 56,684,465 | 61,716,480 | 201,140,420 |
| SRR36741279 | L. major | 106,342,010 | 117,565,440 | 194,786,333 |
| SRR37283774 | P. falciparum | 67,382,606 | 71,168,000 | 94,989,854 |
| SRR32429602 | HCMV | 55,694,667 | 57,344,000 | 106,983,461 |
| SRR39257532 | A. fumigatus | 108,558,425 | 153,589,760 | 227,157,899 |
| SRR29296997 | H. salinarum | 15,654,489 | 17,500,160 | 27,233,618 |
| ERR12954017 | S. acidocaldarius | 15,159,145 | 16,967,680 | 39,545,150 |
| SRR40271341 | H. pylori | 40,620,675 | 45,742,080 | 62,259,448 |
| SRR40402583 | C. jejuni | 9,040,464 | 9,543,680 | 30,863,017 |

Regenerate this table: `results/phase_a/allphases_14dataset.csv`,
column order documented in `docs/SOTA_COMPARISON.md`. **Not included above:
Utricularia gibba (`SRR10676752`), the 15th locked dataset, not yet run —
`docs/CLAIM1_FINAL_VERDICT.md` names this as the one open item.**

Separately, sequence-only content against PgRC2's own binary (the direct
architectural relative, GPL-3, run from its own source, never vendored),
7 datasets, every archive decoded back to byte-identical: **+1.88%
aggregate, 6 wins, 1 loss** (S. acidocaldarius, −0.83%). Speed/RAM on the
same 7 files: compress 106.2 s vs PgRC2's 68.3 s (1.6× slower), 1032 MB vs
371 MB (2.8× heavier); decompress 9.0 s vs 4.0 s (2.2× slower), 806 MB
peak. Full breakdown, including per-stream entropy-bound verification:
`docs/CLAIM1_FINAL_VERDICT.md`.

### FAITHFUL

Against DiscoSNP++ and Kmer2SNP — the only two reference-free callers
applicable to this exact task (single diploid sample, no reference
genome; literature survey in `docs/HET_INDEL_SOTA.md`) — on real GIAB
HG002-HG005 chr20 data, scored by third-party `rtg vcfeval`:

| table | result |
|---|---|
| het-SNV F1 | **WIN** — 0.890 (ours) vs 0.874 (DiscoSNP++) vs 0.464 (Kmer2SNP) |
| multi-allelic sites recovered | **WIN** — 11/18 vs 0/18 (DiscoSNP++ structurally cannot emit true multi-allelic records) |
| het-indel F1 | **loss** — 0.637 vs 0.663, kept on the record rather than dropped; root cause traced to a structural property of the candidate generator, not tuning (`docs/CLAIM2_TABLES_AND_INDEL_SCAN.md`) |

**Open item, not hidden:** every number above is measured on a chr20
window (~75K reads), not the full 30× individual (~12.6M reads) the
project's own spec commits to. `docs/CLAIM2_FINAL_VERDICT.md` names this as
the one blocker to calling this claim locked.

### ADDRESSABLE

Per-operation, against the conventional pipeline that would otherwise
compute the same thing, on real data (E. coli for export/query, real GIAB
windows for coverage):

| operation | vs | speedup | spec target |
|---|---|---|---|
| export | SPAdes v4.0.0 (spec-exact) | **555–656×** | ≥40× |
| coverage | bwa + samtools + mosdepth | **23–33×** | 2–5× |
| query | full decompression | 1.62× in time, **132× fewer reads / 112× fewer bytes** returned | not spec'd — selectivity is the real advantage, not raw speed |

Full numbers, exact commands, and the two real bugs found and fixed while
verifying them: `docs/CLAIM3_LOCKED.md`.

## What's deliberately recorded as failed or open

This project keeps refuted ideas and open gaps on the record rather than
deleting or hiding them:

- **`results/phase_a/04_gate_A2_REVERTED.txt`** — a real change that made
  size worse (+3.2%) because it estimated coder cost from source bytes,
  and bytes don't predict coding time.
- **Four indel-precision filter attempts** (Claim 2) — each measured, each
  neutral or negative, kept behind flags in `caps_caller.h` rather than
  deleted (`docs/HOW_DISCOSNP_WINS.md` §4).
- **A query-speedup figure corrected mid-project** — an earlier "3.3×" was
  found to be measured against the wrong baseline; corrected to the real
  1.62× rather than left standing (`docs/CLAIM3_LOCKED.md` §6.4).
- **A cross-claim bug found by testing the product as a whole**, not one
  claim at a time — names/quality were silently dropped when compressing
  through the candidate-sweep script with a relative input path. Found,
  fixed, verified byte-identical against every existing locked result.
  Full writeup: `docs/FINAL_ALGORITHMIC_SCAN.md`.
- **Open, product-wide**: no CI, no top-level project license, and Claim
  2's full-scale run — see `docs/INDUSTRIAL_CHECKLIST_OVERALL.md` and
  `docs/RESEARCH_CHECKLIST_OVERALL.md` for the complete, current list.

## Third-party

Vendored under `thirdparty/`, each with its license included: PPMd7 (LZMA
SDK, public domain), FSE/Huf0 (Yann Collet, BSD — `thirdparty/fse/LICENSE`),
htscodecs/fqzcomp (BSD 3-clause, `thirdparty/htscodecs/LICENSE.md`).

**Not vendored, cloned/installed separately for benchmarking** (exact,
verified commands for every one of these: `docs/SERVER_SETUP_AND_DOWNLOADS.md`):
PgRC2 (GPL-3), SPRING, Genozip, DiscoSNP++, Kmer2SNP, DSK, MEGAHIT, SPAdes,
bwa, samtools, mosdepth, rtg-tools.

This repository's own code does not yet have a top-level LICENSE file —
an open decision, not an oversight (`docs/INDUSTRIAL_CHECKLIST_OVERALL.md`).
