# benchmark/ — the results, the code that made them, and the chain between

Three folders, one rule: **no number appears here that did not come out of a
script in `scripts/` and land in a file in `results/`.**
`documentation/RESULT_CODE.md` is that chain, table by table.

    scripts/        the authoritative harness (symlinks into ../scripts,
                    so the working code stays where its internal paths expect)
    documentation/  claims, methodology, traceability, reproducibility
    results/        the CSVs, the raw log, the preflight verdict

## Start here

| you want | read |
|---|---|
| the results | `results/*.csv`, or the summary below |
| where a number came from | **`documentation/RESULT_CODE.md`** |
| to re-run everything yourself | **`documentation/REPRODUCE_EVERYTHING.md`** |
| whether it reproduces on another machine | `documentation/REPRODUCIBILITY.md` |
| what the claims actually are | `documentation/CLAIMS_FINAL.md` |
| which table is which | `documentation/TABLE_MAP.md` |

## The eight tables

| table | claim | result | file |
|---|---|---|---|
| T1.1 | COMPACT | **19/19 vs SPRING, 19/19 vs Genozip**, −6.03% / −43.26% | `claim1_T1.1_T1.2.csv` |
| T1.2 | COMPACT | 2.99× SPRING compress; faster on 3/19, lighter on 9/19 | same |
| T2.1 | FAITHFUL | het-SNV mean F1 **0.876** / 0.853 / 0.475 — 3 wins, 1 loss | `claim2_T2.1_snv.csv` |
| T2.2 | FAITHFUL | 10× 0.487 → 30× 0.888, monotone | `claim2_T2.2_coverage_sweep.csv` |
| T2.3 | FAITHFUL | het-indel, 3 wins 1 loss vs DiscoSNP++ | `claim2_T2.3_indel.csv` |
| T2.4 | FAITHFUL | **21/26** (complete chr20 census); DiscoSNP++ 0 of 3,989 records | `claim2_T2.4_multiallelic.csv` |
| T2.5 | FAITHFUL | tetraploid SNV **0.897**/0.782, INDEL **0.597**/0.553 | `claim2_T2.5_tetraploid.csv` |
| T3.1/3.2/3.3 | ADDRESSABLE | export 129–784×, coverage 16–54×, query ×19 | `claim3_T3.1_T3.2_T3.3.csv` |
| T3.4 | ADDRESSABLE | coordinate **81/400**, content **345/400** | `claim3_T3.4_locus_fidelity.csv` |

**57/57 Claim 1 rows LOSSLESS.** All Claim 2 numbers are called **from the
archive** — no FASTQ, no reference.

## How to run it

    bash scripts/benchmark_0_preflight.sh          # must print VERDICT: GO
    CLAIMS=1 bash scripts/sanity_archive_one.sh SRR2584863   # 45 s rehearsal
    bash scripts/benchmark_1_run.sh                # the full sweep, ~5 h

`benchmark_0` is a gate, not a formality: it rebuilds the binaries, checks
every tool/dataset/reference/truth set, exercises the archive-calling path, and
re-derives the validated k-mer count. `benchmark_1` halts the entire run on any
LOSSY archive.

## Methodology, in one place

- **One timed job at a time**, idle machine. Two concurrent jobs contaminate
  both wall time and peak RAM.
- **Timing and RAM** from `/usr/bin/time -v`. `parse_time_v` takes the **last**
  `Elapsed` block and the **max** resident size — a log can hold more than one
  block, and taking them all once glued two numbers into `wall=78.3384.97s`.
- **Lossless** is checked by decoding the **archive**, not the encoder's
  intermediate dumps: the dumps bypass the entropy layer, where one whole class
  of bug lived.
- **Competitors get identical everything** — same reads, same truth VCF, same
  confident-region BED, same normalisation order, same `rtg vcfeval` call. Only
  the caller differs.
- **T3.4 is the deliberate exception**: it measures correctness, not cost, so
  it records no timing and its queries run in parallel. Stated rather than left
  as an inconsistency.
- **A competitor producing nothing fails the run.** A blank column reads as
  "we did not check", which is worse than a tool error.

## Superseded material — kept, not deleted

Earlier benchmark scripts remain in `../scripts/` and are **not** symlinked
here: `run_claim1_bench.sh`, `run_claim3.sh`, `benchmark_final.sh`,
`benchmark_sanity_one.sh`, `run_window_bench_capsule.sh`,
`run_polyploid_bench_capsule.sh`, `run_fullchr20_bench_capsule.sh` (the FASTQ
path, superseded by the archive path), `sim_indel_bench.py`, and others. They
are the evidence for decisions that were later reversed and are kept for that
reason.

Two specific retractions worth knowing before reading any older document:

- **`docs/PHASE2B_RESULT.md` is VOID for 4 of 14 datasets.** Its numbers
  predate four silent data-loss bugs found and fixed the same day.
- **T3.4's "0/400"** appears in documents written before 2026-09-10. The query
  was emitting the consensus rather than the reads; the corrected figure is
  **81/400**. Any document still saying 0/400 is superseded.
