# G_CAPSUL — shipping readiness

> **START HERE INSTEAD, unless you specifically want the edge-case and robustness checks.**
> `README.md` for what the project is and what it measured; `CLAUDE.md` for the
> rules and the refuted-ideas list; `benchmark/documentation/RESULT_CODE.md` for
> where any number came from; `benchmark/documentation/REPRODUCE_EVERYTHING.md`
> to re-run it. **This document predates the final 2026-09-09 sweep, so where it
> disagrees with a file in `benchmark/results/`, the result file wins.**

What was verified, how, and what a user needs to know. Everything here is
measured in this repository, not projected.

---

## 1. One command, end to end

```bash
bash scripts/benchmark_final.sh <reads.fq> [outdir]
```

Builds both binaries, compresses, verifies the round trip, calls variants from
the archive, exercises all three addressable operations, and reports SIZE,
SPEED, RAM and DISK per stage. There is nothing else to run and no flags to
remember. It refuses to start below 40 GB free disk, because a full disk
silently truncates the k-mer spill and produces WRONG results rather than an
error.

**Measured run** (SRR29296997, 201 MB, 460,501 reads, 12 vCPU / 82 GB):

| stage | time | peak RAM | output |
|---|---|---|---|
| build | 38.5 s | -- | two binaries |
| compress | 26.7 s | 957 MB | **15,492,420 B, 13.57x** |
| verify round trip | 1.8 s | -- | **LOSSLESS**, 460,501 reads |
| call from archive | 6.7 s | 811 MB | 7 records, **460,501/460,501 placements** |
| export | 0.20 s | 49 MB | 18,302,206 B |
| coverage | 0.12 s | 97 MB | identity **EXACT** |
| query | 0.20 s | 55 MB | 1,390,256 B |
| **total** | **74.8 s** | | |

Disk: input 210 MB, archive 15 MB, peak scratch 317 MB (the scratch is the
verification decompress and is reproducible; only the archive ships).

---

## 2. Correctness — what is guaranteed and how it is checked

| property | check | result |
|---|---|---|
| sequence losslessness | decompress and `cmp` against the input's sequence column | LOSSLESS |
| coverage exactness | covered bases must equal the sum of read lengths -- one dropped read breaks it | EXACT (69,535,651) |
| export exactness | emitted bases must equal PG_LEN | EXACT (18,002,109) |
| query exactness | returned reads must equal an independent `pos_abs`+`read_lengths` scan | EXACT (7,785, none outside range) |
| placement completeness | every read restored from the archive | 460,501/460,501 |
| encoder determinism | same input, same archive bytes | byte-identical across builds |

## 3. Edge cases — measured, not assumed

| input | encode | round trip | export / coverage / query | call |
|---|---|---|---|---|
| empty file | ok, 361 B archive | -- | -- | -- |
| 1 read | ok | **LOSSLESS** | ok | clean refusal |
| 4-base read | ok | **LOSSLESS** | ok | clean refusal |
| 10 reads | ok | **LOSSLESS** | ok | clean refusal |
| all-N reads | ok | **LOSSLESS** | ok | clean refusal |

Nothing crashes, hangs, or produces a silently wrong answer. The `call`
refusals are correct and actionable -- those archives were built without
`CAPS_CALL=1`, so `contig_spans` is genuinely absent, and the tool says exactly
that and how to fix it. **All-N losslessness matters specifically**: N-read
handling was a historical data-loss bug in this project.

## 4. Known operating limits

- **Disk**: the k-mer counter spills. Below ~40 GB free the spill truncates and
  results are WRONG, not failed. The benchmark refuses to start; any other
  harness must check this itself.
- **`call` requires `CAPS_CALL=1` at compress time** (writes `contig_spans`) and
  `DUMP_PERM=1` (writes `pos_abs`/`pos_strand`/`read_lengths`). Without them the
  tool refuses cleanly rather than guessing.
- **Determinism differs by path**: the SNV+indel caller is deterministic and can
  be gated on byte-identity. The graph-only path renumbers bubble contigs
  non-deterministically, so it must be compared after lifting to genome
  coordinates, never byte-wise.

## 5. Verification gaps, stated plainly

- Encoder speed wins measured on **3 of 15** locked datasets; size win on
  **4 of 15**. Mechanism proven in both cases, aggregate across the locked set
  not yet run. This is measurement, not engineering.
- `scripts/verify_lossless.sh` reports LOSSY for SRR40271341 and ERR552797 --
  **for the baseline encoder too**. It tests the DUMPED streams via
  `decode_105.py`, not the archive, and selects its own parameters
  (`MAXMAP=60 MINOV=105` vs production's `11/16`). The ARCHIVE round trip is
  LOSSLESS for both at production parameters. That script needs its own
  investigation; it is not a defect in the shipping path.
