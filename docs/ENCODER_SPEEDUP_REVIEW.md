# Encoder speedup — correctness review

Reviewed by reasoning over the code, not by running more cases. Every entry
path and every shared object is enumerated below, because two of these changes
are concurrency changes and a passing test does not prove a race absent — one
race in this work passed twice and then hung for 13 minutes.

Verified afterwards: byte-identical archives on 5 datasets spanning 100, 151,
221, 300 and 301 bp reads (including variable-length), and 5 consecutive runs
with no hang and identical output.

---

## A. Parallel trials inside `encode_block` (quality_coder.h)

`encode_block` tries every (offset, strategy) pair and keeps the smallest:
2 offsets x 4 strategies = 8 full compressions of the same block. Those 8 are
now concurrent.

| concern | why it is safe |
|---|---|
| race on `res` / `ok` | each iteration writes only index `t`; the outer vectors are sized before the region and never resized inside it |
| race on `mylen` / `myflags` | declared **inside** the loop body, so they are per-iteration, not shared |
| `fqz_compress` mutates its slice | it writes `s->flags[rec]` (fqzcomp_qual.c:655) and `s->len[i]` (:790). Each trial owns its own arrays. `fqz_slice` has exactly three fields — `num_records`, `len`, `flags` — and all three are initialised |
| race on `shifted` | fully written before the parallel region, read-only inside |
| **determinism** | the selection scan runs `t` ascending, and `t = ci*4 + strat`, so it visits (offset, strategy) in exactly the serial order and uses the same strictly-less test. Ties resolve to the lowest (ci, strat), as before. Completion order cannot affect the result |
| `ncand == 1` | `ntrial = 4`, `ci = t/4 = 0` throughout, `shifted` has one element, `cands[0]` — consistent |
| `ncand == 2` | `ntrial = 8`, `ci = 0` for `t<4` and `1` for `t>=4` |
| `qmin_out` | set to `cands[t/4]`, which is the serial loop's `off` for that trial |
| `fqz_compress` returns NULL | `ok[t]` stays 0 and the trial is skipped, as in the serial version |
| OpenMP compiled out | the pragma is ignored, the loop runs serially, the result is unchanged |
| nested parallelism | `encode_block` is reached from `encode_from_fastq`, which is never called from inside a parallel region |
| memory | `shifted` grows from one copy to `ncand` copies (one extra block), and up to 8 compressed outputs are held. Bounded, and see B: only ONE process now does this |

**The bug this review would have caught, and did not until measured:** the first
version copied the `fqz_slice` struct but not the arrays it points at, so eight
threads wrote the same `len`/`flags`. The archive grew 1.55% and stopped
matching serial output. Copying a struct that holds pointers copies the
pointers, not the storage.

---

## B. Hoisting names and quality above the candidate fork (106_inprocess.cpp)

Neither column depends on MAXMAP or MINOV, so an 8-point sweep encoded the
identical quality column eight times. They are now computed once by the parent
and inherited through copy-on-write.

**Every entry path:**

| path | `gs_child` | `child` | hoist runs? | outcome |
|---|---|---|---|---|
| `CANDIDATES` unset | — | — | no — the enclosing block (1197) is never entered | each process encodes its own; unchanged behaviour |
| GSEARCH parent | false | false | no — it `return 0`s at 1275, before the hoist at 1320 | unchanged |
| GSEARCH child | **true** | false | no — guarded by `!gs_child` | encodes its own |
| candidate parent | false | false | **yes**, then tears the OpenMP runtime down at 1342 | encodes once |
| candidate child | false | **true** | no — forked at 1418, after the hoist | reuses the parent's result via the `g_*_done` guard |

| concern | why it is safe |
|---|---|
| **OpenMP-then-fork deadlock** | this is the one that bit. The teardown at 1211 runs BEFORE the hoist, and the hoist's trial search re-creates the pool, so the fork at 1418 would inherit broken libgomp locks. A second teardown at 1342 closes it. The GSEARCH fork at 1243 precedes the hoist, so it is unaffected |
| `CAPS_NAMES` set, `CAPS_QUAL` not | `g_NM_done` true, `g_QL_done` false; the teardown fires on the `||`; children reuse names and encode quality themselves |
| `CAPS_QUAL` set, `CAPS_NAMES` not | symmetric |
| neither set | no hoist, no OpenMP region, and the teardown is correctly skipped |
| teardown with no live threads | `omp_pause_resource_all` is a no-op in that case |
| memory | quality is now encoded once rather than once per candidate, so peak is LOWER than before, not higher |
| output | names and quality are functions of the input file alone, so every candidate produced identical bytes for them; computing once cannot change the archive |

---

## What is proven vs what is measured

| change | basis |
|---|---|
| concurrent candidates | **proven.** `speedup = (N·S+a)/(S+a) >= 1` for all S,a >= 0, so it can never be slower. Children write distinct files and the parent selects by size, so scheduling cannot change which archive wins |
| 8-point grid | **proven.** The 8-point grid contains the 4-point grid and the smallest archive is kept, so `min(8) <= min(4)` by construction |
| quality hoist + parallel trials | **measured, not proven.** The case analysis above covers every path and shared object, and 5 datasets are byte-identical, but this rests on review plus measurement rather than a structural guarantee |

The distinction matters: the first two cannot regress on a dataset nobody has
run. The third is well-reviewed and well-tested, and that is a weaker claim.

**Backstop:** `benchmark_1_run.sh` decodes every archive and compares it against
the original, and a LOSSY result halts the entire run. A defect on an untested
dataset stops the benchmark rather than silently corrupting a published number.

**Not verified:** the 4-14 GB datasets (each old-binary comparison run costs
10+ minutes) and reads containing N. Neither mechanism depends on input size or
on N handling, but that is an argument, not a measurement.
