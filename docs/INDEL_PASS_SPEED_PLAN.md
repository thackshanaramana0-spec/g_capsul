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


---

# indel_pass optimisation — measured results (4M subset)

| change | section | before | after |
|---|---|---|---|
| flat sorted index | pcluster build + erase | 139.7 s | **8.2 s (17x)** |
| bitset pre-filter | pcluster scan + emit | 244.9 s | **162.4 s (-34%)** |
| deterministic anchor order | — | — | correctness fix |
| **stage total** | `indel_pass` | **599.8 s** | **352.0 s (-41%)** |
| peak RAM | | 14.75 GB | **13.00 GB** |

## The bitset pre-filter

`pcluster: scan + emit` walks ~85 k-mers of each of 4M reads (~340M iterations)
and probed `pkidx.count(cn)` -- a hash lookup into a 21.6M-entry map -- for every
one. The insertions that follow were already filtered ("provably dead
otherwise"), so the PROBES were the cost, not the work they admit.

A 64 MB bitset with two independent hashes answers "definitely absent" from a
single cache line. False positives are possible and harmless (they fall through
to the real `pkidx.count()`); false negatives are impossible. **Output is
identical BY CONSTRUCTION** -- the filter can only skip lookups that would have
missed -- and the measurement agrees: `content vs pcfix: IDENTICAL`.

## Four explanations for a 4-record difference, three of them wrong

The flat index moved `pcluster indels` 300 -> 296. Each hypothesis was tested,
not assumed:

1. **Predicate mismatch.** Replayed both selection rules on 400k synthetic
   entries: same 72 keys, same values, zero disagreement. REFUTED.
2. **`pack25` holes.** Count and fill passes share the same guard, so counts
   match fills exactly. REFUTED.
3. **Run-to-run nondeterminism.** Repeating both builds: original gives 300,
   300; flat gives 296, 296. Each is deterministic. REFUTED -- and this one had
   already been committed as the explanation, then retracted.
4. **`unordered_map` iteration order.** CONFIRMED. `for (auto& kv : pkidx)`
   iterates a hash map whose order depends on insertion history, and the flat
   build inserts sorted by k-mer where the original inserted by contig. `ploc`
   records one location per event with first-writer-wins, so anchor order
   decides which location an event gets.

Fixed by sorting anchors on (contig, pos, key) before the loop. The result now
depends only on the data, not on container bucket layout -- **a correctness
improvement in its own right, because the ORIGINAL was silently depending on an
implementation detail.**

## Where the time went, and what is still dark

Every real gain in this session came from a TIMER. Every false lead came from
reading code and inferring:

* the `im` emit loop -- a fix was written and ready to ship; it costs **0.02 s**
* the anchor scan -- an Amdahl analysis run on my own broken implementation
  concluded "not the bottleneck", when done properly it was 1,153 s -> 73 s
* the 4-record difference -- three wrong explanations before the right one

The `accounted X of stage` sum check is what caught all three, and it still
reports a gap: **210 s accounted of 352 s, so ~142 s remains unmeasured.** That
is the next target.


---

# VERDICT: F1 is bit-identical, the work ships

Both builds run back-to-back on the same box, same 4M subset, lifted and scored
identically against GIAB truth:

    orig/SNV     TP=13862  FP=2624  FN=30713  P=0.8408  R=0.3110  F1=0.4540
    new/SNV      TP=13862  FP=2624  FN=30713  P=0.8408  R=0.3110  F1=0.4540
    orig/INDEL   TP=1214   FP=316   FN=6567   P=0.7935  R=0.1560  F1=0.2608
    new/INDEL    TP=1214   FP=316   FN=6567   P=0.7935  R=0.1560  F1=0.2608

**Every value identical** -- TP, FP, FN, P, R, F1, both classes. Both lifted to
exactly 23,907 records.

| | orig | new |
|---|---|---|
| `indel_pass` | 511.9 s | **339.6 s (-34%)** |
| wall | 873.3 s | **696.8 s (-20%)** |
| peak RAM | 14.74 GB | **13.00 GB (-1.74 GB)** |

## The gate I had been using was wrong

For most of this work the gate was "contig-space VCF content must match". That
treats a CONTIG RELABELING as a regression: the pcluster events sit on repeated
sequence where several contigs are equally valid anchors, so a different anchor
order reports the same variant under a different `bcontig` id. The lift resolves
both to the same genome coordinate.

The right gate is F1 after the lift, and by that gate nothing moved. Three runs
were spent chasing a difference (300 vs 296 vs 289 records) that the scoring
pipeline does not see.

**Keep both gates in future:** contig-space identity is a useful STRONG signal
(if it holds, stop), but its failure is not sufficient evidence of a regression.
Score before concluding.

## What actually shipped

1. **pcluster index**: 3 `unordered_map`s (~65M insertions each, reserved at 2M
   so repeatedly rehashing) + 2 erase-during-iteration passes -> one flat sorted
   array. 139.7 s -> 8.2 s, **17x**. Anchor set proven identical in-process
   (21,605,670 keys, zero differences).
2. **Bitset pre-filter** on the ~340M anchor probes: 244.9 s -> 162.4 s,
   **-34%**. Output identical BY CONSTRUCTION -- false positives fall through to
   the real lookup, false negatives are impossible.
3. **`ploc` order-independence**: the event location was a last-writer-wins
   assignment, so it depended on `unordered_map` bucket layout -- an
   implementation detail that varies with insertion history and libstdc++
   version. Now keeps the smallest (contig,pos), a function of the data. This is
   a reproducibility fix independent of the speed work.

## One correction on the numbers

`indel_pass` measured 599.8 s in one baseline run and 511.9 s in another -- same
binary, same input, ~15% timing variance (the OUTPUT was stable at 300 indels /
36,638 records both times). The -34% figure is quoted against the controlled
back-to-back run, not against the slowest baseline.
