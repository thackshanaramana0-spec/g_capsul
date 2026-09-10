# RESULT_CODE — every reported number, traced to code and to an output file

The purpose of this file is to make it impossible to quote a number that did
not come out of an executed script. For each table: what went in, what it was
compared against, how it was configured, what was measured, **which script
produced it**, and **which file holds the raw output**.

Paths are relative to the repository root (`c_star_pg_advance/`).

**Source of truth.** Where a document and a result file disagree, the result
file wins. Where a result file and the script that would regenerate it
disagree, re-run the script. Nothing here was typed in by hand from memory.

**Provenance of the sweep.** All 19-dataset numbers come from one run:
`benchmark_1_run.sh` at commit `21ee619`, started 2026-09-09 02:46:07, finished
07:53:21 (5 h 07 m), 19 datasets, 0 failures. It was executed from a **clean
git worktree at HEAD**, not the working tree, because the working tree carried
a concurrent agent's uncommitted encoder edits and `benchmark_0` rebuilds the
encoder from source — running in place would have benchmarked that work instead
of the baseline. The binaries were verified free of those symbols before launch.
Raw log: `benchmark/results/run.log.txt`. Checkpoints:
`benchmark/results/PROGRESS.txt`. Preflight verdict:
`benchmark/results/benchmark_0_status.txt`.

**Code changes after the sweep, and whether they affect any published
number.** Three commits landed after `21ee619` and before this document's
current revision, each fixing a defect in how the encoder handles input
*outside* the locked 19-dataset set (an unconditional refusal on reads past
1023 bases, a `std::thread` teardown crash exposed while fixing that, and a
names-tokenizer bound on pathological headers — see `industry/README.md` for
the full account and `paper/LIMITATIONS.md` §0):

    b28dd40  the three defect fixes
    478a241  scope self-test + industry/ safety layer
    7410a6a  paper documentation of the scope decision

**None of them changes a single published number.** Verified directly, not
assumed: the canonical checkpoint
(`CLAIMS=1 bash benchmark/scripts/sanity_archive_one.sh SRR2584863`) reproduces
**68,429,027 B, byte-identical**, both before and after all three commits, and
the fixes only ever engage on input that fails a condition none of the 19
locked datasets meet (every read in every locked dataset is well under the
structural bound; every header is far short of the tokenizer bound). The
manifest below fixes the exact file states this applies to.

**Independent cross-check of every table against the raw log (2026-09-10).**
Each published CSV was re-derived from `benchmark/results/run.log.txt` by a
parser written against the log format, then diffed against the CSV. Results:

| table | rows re-derived | matched |
|---|---|---|
| T1.1 + T1.2 | 57 | **57** — including every LOSSLESS verdict |
| T2.1 het-SNV | 12 | **12** |
| T2.2 coverage | 3 | **3** (30× is carried from T2.1, not re-run) |
| T2.3 het-indel | 8 | **8** |
| T3.1 export | 6 | **6** — output bytes, contig count and speedup |
| T3.2 coverage | 6 | **6** |
| T3.3 query | 19 | **19** |
| T2.4 multi-allelic | — | **no raw log existed; re-executed instead, see below** |
| T2.5 tetraploid | — | **re-executed instead: all 16 values confirmed exactly** |
| T3.4 locus fidelity | — | **not in this log; re-executed instead, see below** |

**T2.5 had no independent evidence in the log** — the sweep's end-of-run
section prints it via `column -s, -t "$CSV7"`, i.e. it pretty-prints the CSV it
is supposed to corroborate, which is circular. It was therefore **re-executed
from scratch on 2026-09-10, and all 16 published values reproduce exactly**
(every TP, FP, FN, precision, recall and F1 for both tools and both variant
classes), along with the intermediate counts. Raw log and full comparison:
`benchmark/documentation/T2.5_reproduction_20260910/`.

**File integrity.** `benchmark/documentation/MANIFEST.sha256` lists a SHA-256
for every file a published number depends on — every results CSV and log,
every paper document, and the encoder/decoder source itself. Verify with
`sha256sum -c benchmark/documentation/MANIFEST.sha256` from the repository
root; a mismatch means something changed after the manifest was generated and
should be investigated before trusting any number here.

---

## Claim 1 — COMPACT

### T1.1 archive size · T1.2 wall time + peak RAM

| | |
|---|---|
| **input** | all 19 locked datasets, 55.0 GB FASTQ (`server/DATASETS.md`) |
| **compared with** | SPRING (source build), Genozip 15.0.87 |
| **configuration** | `CAPS_SPANS=1 CAPS_NAMES=1 CAPS_QUAL=1`, adaptive candidate sweep. **Not `CAPS_CALL=1`** — that also runs the caller inline at ~20x the cost |
| **metric** | archive bytes; wall seconds and peak RSS from `/usr/bin/time -v`; lossless verdict by decoding the ARCHIVE and comparing the sequence column |
| **RAM/time method** | one timed job at a time, idle box; `parse_time_v` takes the LAST `time -v` block and the MAX resident size |
| **script** | `benchmark/scripts/benchmark_1_run.sh` → `run_phase1` |
| **encoder invoked via** | `benchmark/scripts/encode_adaptive.sh` (never the binary directly) |
| **output** | `benchmark/results/claim1_T1.1_T1.2.csv` (57 rows = 19 × 3 tools) |

**Reported:** 19/19 wins vs SPRING, 19/19 vs Genozip, **57/57 rows LOSSLESS**.
Aggregate ours **6,011,370,147 B** vs SPRING 6,396,897,280 (**−6.03%**) vs
Genozip 10,594,522,604 (**−43.26%**). Time 2.99× SPRING compress, 1.79×
decompress; faster than SPRING on 3 of 19 and lighter on 9 of 19.

**Independently re-verified 2026-09-10** by a clean rebuild + one-dataset run:
E. coli archive **68,429,027 B**, byte-exact, 3/3 LOSSLESS.
Script: `benchmark/scripts/sanity_archive_one.sh` with `CLAIMS=1`.

**+1.88% vs PgRC2** is a separate, earlier measurement on the DNA stream only
(PgRC2 stores no names, no quality, no line 3 — verified by running it: a
440,048 B FASTQ in, 151,000 B of bare DNA out). It is not part of this sweep
and is not in the CSV.

---

## Claim 2 — FAITHFUL. Every number is called FROM THE ARCHIVE.

Common to T2.1–T2.3: reads are never re-read. The caller runs
`capsule_decode call <archive>` with no FASTQ and no reference. Scoring is
`rtg vcfeval --squash-ploidy` against GIAB v4.2.1 inside the confident-region
BED, het-restricted, with normalisation before classification. DiscoSNP++ gets
the identical truth, regions, normalisation and scorer — only the caller
differs.

### T2.1 het-SNV F1 (three arms, four individuals)

| | |
|---|---|
| **input** | HG002-HG005 chr20 @ 30x, from the archives Phase 1 kept |
| **compared with** | DiscoSNP++ (`-k 31 -c 3 -D 100 -P 3 -b 0 -G ref`, POS off-by-one corrected), Kmer2SNP (KMC counts, band derived from the histogram) |
| **metric** | TP/FP/FN, precision, recall, F1 from `rtg vcfeval` |
| **scripts** | ours `run_fullchr20_archive_capsule.sh`; DiscoSNP++ `run_fullchr20_bench_disco.sh`; Kmer2SNP `run_kmer2snp.sh` |
| **output** | `benchmark/results/claim2_T2.1_snv.csv` |

**Reported:** mean F1 **0.876** (ours) / 0.853 (DiscoSNP++) / 0.475 (Kmer2SNP).
Per individual 0.888 / 0.891 / 0.891 / 0.834. **3 wins, 1 loss.**
HG005's loss is explained in `HG005_EXPLAINED.md` — it is the only
variable-length dataset.

### T2.2 coverage sweep

| | |
|---|---|
| **input** | HG002 subsampled with `seqtk sample -s11` at fraction depth/30 |
| **configuration** | each depth is a full re-compress and re-call, so depth is the only variable. 30× is **carried from T2.1**, not re-run |
| **script** | `benchmark_1_run.sh`, the T2.2 block |
| **output** | `benchmark/results/claim2_T2.2_coverage_sweep.csv` |

**Reported:** 10× **0.487**, 15× **0.704**, 20× **0.797**, 30× **0.888** — monotone.

### T2.3 het-indel F1

Same scripts and scoring as T2.1; the INDEL line was already being produced in
the same pass and is now captured. **Output:** `claim2_T2.3_indel.csv`.

**Reported:** ours 0.620 / 0.637 / 0.632 / 0.593 vs DiscoSNP++ 0.576 / 0.587 /
0.595 / 0.605. **3 wins, 1 loss** (HG005, same cause). Kmer2SNP records
`NOT_APPLICABLE_SNP_ONLY` — it has no indel model, so the blank is a property
of the tool.

### T2.4 multi-allelic sites recovered

| | |
|---|---|
| **input** | HG002; every `GT=1/2` site on chr20 whose two ALT alleles are both single-base |
| **compared with** | DiscoSNP++ on the identical reads |
| **script** | `benchmark/scripts/run_multiallelic_bench_capsule.sh` |
| **output** | `benchmark/results/claim2_T2.4_multiallelic.csv` |

**Reported:** **17/26 (65.4%)** — and **26 is the complete census** for chr20,
not a sample. DiscoSNP++ emits **0 multi-allelic records in 3,989 total
records**: a structural inability, not a miss.

> **CORRECTED 2026-09-10 — this table previously reported 21/26 (80.8%). That
> figure was withdrawn during a reproducibility audit: no raw log anywhere on
> this server supported a caller run against the full chromosome that
> recovered 21 sites, and three independent re-derivations — two different
> encoder builds (before and after commit `b28dd40`), one exact repeat on
> identical input — all gave **17/26**, with identical intermediate SNV/indel/
> contig counts every time. The denominator (26 SNV-only sites out of 952
> total multi-allelic sites) is independently confirmed and unchanged; only
> the recovery count was wrong. Full investigation, all three raw logs, and
> what the original evidence actually contained:
> `benchmark/documentation/T2.4_reproduction_20260910/README.md`.

### T2.5 tetraploid SNV + indel

| | |
|---|---|
| **input** | HG003+HG004 real reads concatenated, ploidy 4, truth = union of their GIAB calls (Cooke, Wedge & Lunter, *Genome Research* 2022) |
| **script** | `benchmark/scripts/run_tetraploid_bench_capsule.sh` |
| **output** | `benchmark/results/claim2_T2.5_tetraploid.csv` |

**Reported:** SNV **0.897** vs 0.782; INDEL **0.597** vs 0.553.

---

## Claim 3 — ADDRESSABLE

### T3.1 export · T3.2 coverage · T3.3 query

| | |
|---|---|
| **input** | the archives Phase 1 kept; 6 datasets spanning kingdoms for T3.1/T3.2, all 19 for T3.3 |
| **compared with** | T3.1 SPAdes 4.0.0; T3.2 bwa + samtools sort + index + mosdepth **timed as one pipeline** |
| **metric** | wall + peak RSS via `/usr/bin/time -v`, plus output bytes and row count for both sides |
| **script** | `benchmark_1_run.sh`, `run_phase3` |
| **output** | `benchmark/results/claim3_T3.1_T3.2_T3.3.csv` (31 rows, tagged by sub-table) |

**Reported:** export **129–784×** vs SPAdes; coverage **16–54×** vs
bwa+mosdepth; query on all 19.

**Read T3.1 with its caveat, which is why `output_bytes`/`rows` are in the
table:** our export emits **2 records** (the pseudogenome), SPAdes emits tens of
thousands of biological contigs. The ratio measures time-to-a-reference-free-
coordinate-system from an archive that had to exist anyway — **not "the same
output, faster"**. T3.2 and T3.3 carry no such caveat.

**T3.3 is not a speed win.** `genocat --head=100` extracts in 0.19 s against our
0.46 s. With the optional sidecar index we reach 0.03 s, but Genozip is the
honest comparator and Claim 3 does not lead with speed.

### T3.4 locus retrieval fidelity — the mechanism table

| | |
|---|---|
| **input** | HG002-HG005, 100 GIAB het SNV sites each, chr20:3.0-3.6 Mb |
| **compared with** | **an internal control**: coordinate addressing vs content addressing on the SAME archive and the SAME sites |
| **why internal** | no other tool can produce a row. SPRING addresses by read index; `genocat --regions` is refused on FASTQ; PgRC2/NanoSpring expose no read-out; BEETL returns reads *containing* a string, not covering a locus; CRAM needs an external reference |
| **probe design** | 40 bp of REFERENCE ending 6 bp **before** the variant, so it never contains the variant — a probe taken from one haplotype's own sequence could only match that haplotype and would rig the result |
| **metric** | does the returned read set contain BOTH truth alleles |
| **script** | `benchmark/scripts/run_locus_fidelity.sh` |
| **REQUIRED SETUP** | the sidecar must exist and must have been built **with `CAPS_PILEUP=1`**: `CAPS_PILEUP=1 capsule_decode index <archive> <archive>.qidx`. Without that flag the sidecar omits per-read deviations, `query` emits the consensus instead of the reads, and the coordinate arm scores **0**, not 18 — reproduced three times on 2026-09-10 before the cause was found. `query` finds `<archive>.qidx` automatically. |
| **output** | `benchmark/results/claim3_T3.4_locus_fidelity.csv` |

**Reported:** coordinate **81/400 (20.2%)**, content **345/400 (86.2%)**,
pooled allele balance **1.00**.

> **Correction, and it must not be lost.** An earlier version of this table
> reported coordinate at **0/400**. That was an artifact of our own query
> emitting the CONSENSUS rather than the reads, so no allele difference could
> ever appear and coordinate addressing looked absolutely incapable. The query
> now applies each read's deviations (commit `d58fc23`), and the honest figure
> is a **4.3× gap, not an infinite one**. Any document still saying 0/400 is
> superseded.

**INDEPENDENTLY REPRODUCED 2026-09-10.** The HG002 row was re-derived from a
freshly built archive (573,767,964 B) with a freshly built decoder and matches
the published row in **every field**: `100,18,82,0,85,15,0,2771,3023`. Full
account, including the undocumented `CAPS_PILEUP=1` flag that must be set when
building the sidecar or the coordinate arm silently scores 0 instead of 18:
`benchmark/documentation/T3.4_reproduction_20260910/README.md`. HG003/HG004/
HG005 were not re-run.

**Window independence** (4 further windows on HG002 at 10/25/40/55 Mb, 60 sites
each) was measured under the *earlier* consensus-emitting query and gave
210/240 content vs 0/240 coordinate. Those coordinate figures carry the same
artifact and are retained only as evidence that the **content** result does not
depend on one window.

---

## What is deliberately NOT in the results

- **PgRC2 comparison** (+1.88%) — earlier, DNA-stream only, not in this sweep.
- **The ablation** (0.431 / 0.426 / 0.648 / 0.888) — a mechanism experiment, not
  a benchmark table. Script: the `CAPS_NO_COLLAPSE` / `CAPS_NO_REMAP` toggles;
  documented in `CLAIM2_RESULTS_V2.md`.
- **HG005 truncated to fixed length** (F1 0.897) — a diagnostic that identified
  a cause. **It is not an HG005 result** and must never be reported as one:
  truncating discards 40% of every read and scores a different dataset.
- **The sidecar-accelerated query** (0.03 s) — an optional index, not the
  archive. The archive-alone figure is 0.46 s.
