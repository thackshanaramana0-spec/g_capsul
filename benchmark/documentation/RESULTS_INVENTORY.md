# RESULTS_INVENTORY — which run directory is the real one

`results/` at the repo root holds twelve run directories. Exactly **one** is the
source of the published tables. The others are aborted runs, rehearsals and
diagnostics, kept because the project's rule is to archive uncertain evidence
rather than delete it — but they are **not citable**, and several contain
numbers that were later corrected.

## The authoritative run

    results/FULL_SWEEP_20260909/     8 CSVs, 125 data rows, run.log.txt, PROGRESS.txt

`benchmark_1_run.sh` at commit `21ee619`, started 2026-09-09 02:46:07, finished
07:53:21 (5 h 07 m), 19 datasets, 0 failures, executed from a clean git worktree
at HEAD.

**`benchmark/results/` is a byte-identical copy of it** — verified with `cmp` on
all 11 files (8 CSVs + the three logs), 11/11 identical. Nothing was edited on
the way into `benchmark/`. Cite either path; they are the same bytes.

Two tables in it were **re-measured after the sweep and written back**, and the
copies agree:

| table | why re-measured | value now in both copies |
|---|---|---|
| T2.4 | the windowed 5/111 scored 104 indel-bearing sites a single-base check cannot evaluate | **21/26**, complete chr20 census |
| T3.4 | our query emitted the consensus, so no allele difference could appear | coordinate **81/400**, content **345/400** |

## Everything else, and what it actually is

| directory | what it is | why it is not citable |
|---|---|---|
| `benchmark_1_20260908_212344` | first launch | 3 rows, aborted |
| `benchmark_1_20260908_213517` | second launch | 6 rows, aborted (35 MB of transient output) |
| `benchmark_1_20260908_213722` | third launch | 7 rows, aborted |
| `benchmark_1_20260909_000229` | fourth launch | 6 rows, aborted |
| `benchmark_1_20260909_011600` | last pre-sweep rehearsal | 24 rows, partial — this is the run that exposed the `parse_time_v` bug (`wall=78.3384.97s`) |
| `sanity_HG002_*` (3 dirs) | single-dataset rehearsals | one dataset, plus multi-GB working files |
| `sanity_SRR2584863_20260908_192453` | E. coli rehearsal | superseded by the one below |
| `sanity_SRR2584863_20260910_022835` | **the independent re-verification** | not a table, but the useful one: a clean rebuild reproduced the E. coli archive at **68,429,027 B byte-exact**, 3/3 LOSSLESS |
| `claim2_timing` | a timing note | diagnostic |

## The one number to check on any machine

    CLAIMS=1 bash benchmark/scripts/sanity_archive_one.sh SRR2584863
    ->  archive 68,429,027 B exactly, ratio 9.844%, lossless 3/3

Verified byte-identical across 2/3/7/12 cores and every candidate-concurrency
setting. **If this differs on any machine of any size, it is a real regression,
not an environment difference.**
