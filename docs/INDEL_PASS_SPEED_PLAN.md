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
