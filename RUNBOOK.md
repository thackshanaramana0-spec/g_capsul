# RUNBOOK — what to run, in order

Four scripts, one order, no guessing. Each is idempotent and each refuses
rather than guesses when a precondition is missing.

| # | command | what it does | how long |
|---|---|---|---|
| 0 | `make test` | self-contained correctness suite. No datasets, no network. **Run this first, always.** | ~3 min |
| 1 | `bash scripts/benchmark_0_preflight.sh` | verifies every tool, dataset, reference and truth file the full run needs. Changes nothing. | ~2 min |
| 2 | `bash scripts/benchmark_sanity_one.sh` | ONE dataset, all three claims, end to end. The rehearsal. | ~15 min |
| 3 | `bash scripts/benchmark_1_run.sh` | the full run: Claims 1, 2, 3 across the locked set, **including SPRING, Genozip and DiscoSNP** | hours |

`scripts/benchmark_final.sh <reads.fq>` is a fifth, lighter entry point: one
dataset, our tool only, no competitors, reporting size/speed/RAM/disk per stage.
Use it to characterise a single input quickly. It does NOT replace step 3 --
the comparison tables come from `benchmark_1_run.sh`.

## Rules these scripts obey, and you must too

1. **One timed job at a time.** Concurrency corrupts wall time and peak RAM.
   Every script here serialises; if you run two yourself, both sets of numbers
   are void.
2. **Disk is a correctness dependency, not a convenience.** The k-mer counter
   spills; below ~40 GB free the spill truncates and results are WRONG rather
   than failed. Step 1 checks this. So does `benchmark_final.sh`.
3. **An archive that must serve `call` needs `CAPS_CALL=1` and `DUMP_PERM=1`
   at compress time.** Without them the tool refuses cleanly and names the
   missing stream. See `docs/FORMAT.md`.

## Interpreting results

- The SNV+indel caller is deterministic: two runs give byte-identical VCFs, so
  regressions can be gated on `cmp`.
- The graph-only path is NOT: it renumbers bubble contigs, so compare it after
  lifting to genome coordinates, never byte-wise. See
  `docs/CLAIM3_POSTSPLIT_AUDIT.md` for the measurement.
