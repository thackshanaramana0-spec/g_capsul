# Making the full CAPS_CALL path competitive

## The actual scoreboard (full chr20, same box)

| | time | RAM | SNV F1 |
|---|---|---|---|
| DiscoSNP++ | 75.7 s | 3.29 GB | 0.847 |
| **ours, Option 0 (SNV only)** | **111.5 s** | **6.21 GB** | **0.8766** |
| **ours, full CAPS_CALL (SNV+indel)** | **~1900 s** | **~14.2 GB** | 0.8766 |

Option 0 is 1.47x time / 1.89x RAM off DiscoSNP -- close, and it wins on F1.
**The full path is 25x time / 4.3x RAM off, and that is where indels live.**

(RAM figures for this path must be sampled, not read off a stage boundary.
Observed during one T5 run: 10.5 GB entering `indel_pass`, 14.2 GB mid-stage,
**22.96 GB** later in the same stage. `indel_pass` alone adds 12+ GB. The
identified structures -- the second seed index (~2 GB) and `kidx` (1.45 GB) --
do NOT account for that, so something else in the stage is allocating heavily
and has yet to be found. Do not quote a single number for this path without
saying when it was sampled.)
So the competitive problem is not Option 0; it is `indel_pass`.

## Where the ~1900 s goes

| stage | time | share |
|---|---|---|
| ridx_build (substrate 1) | 147.2 s | 7.7% |
| kc_H_build | 22.2 s | 1.2% |
| parallel_loop | 491.4 s | 25.7% |
| filter_snv_emit | 0.9 s | 0.0% |
| **indel_pass (incl. substrate 2)** | **1250.0 s** | **65.4%** |

`indel_pass` and `parallel_loop` are 91% of the run.

## Finding 1 — `indel_pass` is ENTIRELY SERIAL

There is not a single `#pragma omp` in its ~1,100 lines. 1,250 s of
single-threaded work on a 12-core box, while `parallel_loop` -- the stage right
before it, over the same contigs -- is parallelised.

Three loops over `cdb.contigs` and one over all reads dominate it:

1. **`cov` allocation** (`:4484`) — each `ci` touches only `cov[ci]`. Trivially
   parallel.
2. **read-coverage accumulation** (`:4486`) — `++cov[cid][p]`, a scatter. Each
   read belongs to exactly ONE `cid`, so partitioning BY CONTIG (not by read)
   is race-free without atomics.
3. **`kidx` build** (`:4510`) — ~65M entries x 24 B = **~1.5 GB**, built by
   serial `push_back`. Parallelisable with per-thread buffers concatenated in
   contig order.
4. **`stable_sort` of kidx** — serial. `__gnu_parallel::stable_sort` is a drop-in
   that keeps the stability the code depends on.

**CRITICAL CONSTRAINT.** The comment at `:4518` records that a plain `sort`
changed which pair is recorded first and therefore moved downstream anchor
positions, DP and AF -- caught by byte-comparison, not by reasoning. So any
parallel version must preserve *exactly* the insertion order (contig 0,1,2...,
each in increasing position) and be verified by byte-comparing the VCF, never
by inspection.

Expected: if `indel_pass` parallelises like `parallel_loop` did (~8x on 12
cores), 1250 s -> ~160 s, i.e. **~1900 s -> ~800 s**.

## Finding 2 — the seed index costs 4 GB, and half of it is waste

Measured: DBG_ONLY peaks at 6.21 GB, full CAPS_CALL at 10.25 GB. The 4.04 GB
gap is `build_substrate`'s seed index, built TWICE (two substrates at different
`dup_frac`, so they are genuinely different and cannot be shared).

Inside each build:

    flat (key, ci<<32|pos)   16 B x ~65M positions = 0.97 GB
    skey + scid + spos       16 B/entry, reserved at flat.size() = 0.97 GB
    both live simultaneously                        = 1.94 GB
    x2 substrates                                   = 3.87 GB   (observed 4.04)
    kidx: 24 B x ~65M positions                     = 1.45 GB   (indel_pass only)

The arithmetic predicted the gap to within 4% BEFORE any change, which is what
makes this a diagnosis rather than a guess.

**Fix applied:** the 64-hits-per-key cap means the outputs are usually far
smaller than `flat`, but they were reserved at `flat.size()` regardless. Now
counted first and reserved exactly, plus `flat.shrink_to_fit()` after the build
(`reserve(contigs*128)` over-allocates). Output is unchanged by construction --
same entries, same order, only the allocation size differs.

## What is NOT worth doing

**Streaming `kc` from disk, GATB-style.** DiscoSNP holds 3.29 GB because GATB
streams a 0.99 GB HDF5 graph rather than keeping it resident. Our traversal does
**1.1 billion** probes into `kc`; at ~80 ns in RAM that is bounded, at ~10 us
per NVMe seek it is hours. GATB gets away with it because its structure is a
Bloom-filter cascade answering *membership* without touching disk in the common
case -- a different data structure, not ours paged out. We need exact counts
(`MINC`, `H`, per-path support), so a probabilistic index changes what the
caller can ask.

**Also note our `kc` is legitimately 2.33x DiscoSNP's** (140.7M vs 60.3M k-mers)
because `MINC=2` keeps what their `-c 3` discards -- that is where the recall
advantage comes from. 6.21 GB is close to the floor for exact counts at that
sensitivity.

## Order of work

1. `indel_pass` parallelisation (the 65%), byte-verified against the current VCF
2. seed-index exact reserve (applied, untested -- ~1 GB)
3. re-measure both paths and update the scoreboard


---

# indel_pass, fully profiled (2026-09-06, 4M-read subset)

Nine timers plus a sum check that prints `accounted X of stage` and says
explicitly that a mismatch means the timers are wrong. That check is the single
most useful thing added, because THREE separate targets in this session were
identified by reading code and were wrong:

| section | cost | share |
|---|---|---|
| anchor scan + merge | 39.0 s | 6% |
| pcluster: build 3 hash maps | 121.9 s | 20% |
| pcluster: 2 erase passes | 17.8 s | 3% |
| pcluster: scan + emit | 251.3 s | 42% |
| im emit loop | **0.02 s** | 0% |
| setup / medcov / scan_pair / matching | 1.9 s | 0% |
| **UNACCOUNTED** | **167.9 s** | **28%** |
| indel_pass total | 599.8 s | |

**The `im` emit loop is 0.02 s.** A fix for it (`kcount_at`, removing a
per-probe `substr` allocation) was written and ready to ship as "the answer"
before the timer showed it costs 20 milliseconds.

## pcluster index: flat sorted array instead of 3 hash maps -- FAST, BUT FAILS THE GATE

    build index       121.9 s -> 0.67 s
    erase passes       17.8 s -> 7.56 s (sort+select)
    139.7 s -> 8.2 s = 17x
    indel_pass        599.8 s -> 443.7 s (-26%)
    peak RAM         14.75 GB -> 13.00 GB
    pcluster indels       300 -> 296     <-- FAILS
    VCF content                DIFFERS   <-- FAILS

**CORRECTION — it SHIPS. The "gate failure" was my own measurement error.**

An in-process check (`CAPS_PCLUSTER_VERIFY=1`) runs BOTH selections on the SAME
contigs in the SAME process:

    [PCLUSTER-VERIFY] original=21,605,670  flat=21,605,670
                      only_orig=0  only_flat=0  value_diff=0

Identical anchor count, identical keys, identical values. The selections are
equivalent.

**AND THAT EXPLANATION WAS ALSO WRONG — retracted a second time.** Repeating
both builds shows they are each DETERMINISTIC and disagree by exactly 4:

    original build:  300, 300
    flat build:      296, 296

So it is not run-to-run nondeterminism. The anchor SET is identical (proven
above); the difference is downstream.

**Root cause, found:** `for (auto& kv : pkidx)` iterates an `unordered_map`,
whose order depends on INSERTION HISTORY -- and the flat build inserts sorted by
k-mer where the original inserted by contig. `ploc` records one location per
event with first-writer-wins, so anchor order decides which location an event
gets, and 4 events resolve differently.

**Fix: sort the anchors by (contig, pos, key) before the loop.** The result then
depends only on the data, not on container bucket layout. That is a correctness
improvement in its own right -- the ORIGINAL was quietly depending on an
implementation detail -- and it makes both builds agree.

Two hypotheses for the 4 lost indels were tested and REFUTED:
1. *Selection predicate differs.* Replayed both rules on 400k synthetic
   entries: same 72 keys kept, same values, ZERO disagreement.
2. *`pack25` skips create zero-filled holes.* Both the count pass and the fill
   pass use the same guard, so counts match fills exactly.

Remaining explanation: the flat version keeps **21,605,670 unique forward
anchors**, where the original's map was seeded at 2M and shrunk by two erase
passes. More surviving anchors changes which reads hit the `vec.size() < 200`
cap in `rc_reads` downstream, and 4 indels fall out differently. That is a
behavioural difference, not a predicate bug -- and it means the ORIGINAL was
silently capacity-limited.

## The real remaining target: 340M hash probes

`pcluster: scan + emit` (251 s) is dominated by one line:

    if (!pkidx.count(cn)) continue;     // ~340M probes, one per read k-mer

4M reads x ~85 k-mers each, each an `unordered_map` probe. The insertions after
it ARE filtered (`provably dead otherwise`), so the lookups themselves are the
cost. Keeping `pkidx` as a sorted array with a Bloom/bitset pre-filter would
reject the ~99% that miss without touching the map.

Plus **167.9 s still unaccounted** after nine timers, sitting between
`pcluster: scan + emit` and the stage end. Not yet read.

## Honest projection, corrected twice

I projected "522 -> ~80 s" before the timers split the stage; that was wrong by
5x. With the index fix alone it is ~460 s. The 419 s of scan + unaccounted tail
(70% of the stage) is where a drastic cut has to come from, and neither is
built.
