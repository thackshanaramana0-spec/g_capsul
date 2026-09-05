# Where the encoder's time actually goes — five hypotheses, one win

Written after instrumenting rather than reasoning. Three earlier attempts at
this stage were wrong because they were argued from timers and comments instead
of measured, so this file records the measurements, the refutations, and the
one change that survived.

## The measured breakdown (3M human reads, 12 threads, uninstrumented)

| stage | time | share |
|---|---|---|
| MEM matching (`allrefs built`) | 22.1 s | 34% |
| pigeonhole mapping | 15.5 s | 24% |
| stream coding pool | 8.3 s | 13% |
| round 1 (division) | 7.0 s | 11% |
| round 2 (assembly) | 2.9 s | 4% |

Scaling to the full 12.6M-read file: round 2 15.7 s, pigeonhole 128 s, MEM
37 s. Pigeonhole scales worst (8.1x for 4.2x more data); MEM is sublinear.

## The probe funnels

Round 2, 3M human reads:

    probes 94,268,287  seed_hit 16.1%  cand_examined 371,742,373  links 1.89%

Pigeonhole, same input:

    seeds 198,569,468 -> maybe 40.3% -> find 16.4%
    candidates 1,931,909,447 (59.5 per find)
    99.5% reach the full read-vs-text comparison
    accepted ~2,000,000 (0.1%)

Both stages are dominated by a CANDIDATE WALK, not by failed lookups. The whole
sweep runs at IPC 0.62 with 3.9 billion cache misses: memory-latency bound.

## What worked

**Sequential discriminator in round 2 (shipped).** `pext[]` holds the 8 bases
following each index entry's seed, in index order, so the walk reads it off the
cache line already carrying `pent[q]`. A candidate whose bases [32,40) differ
cannot pass `rcmp`, so it is rejected without touching the read store. Applied
only when `L >= SW+8`, since below that the extra bases lie outside the overlap.

    round 2   3.23 s -> 2.90 s
    total     65.13 s -> 61.26 s   (-6%)
    archive   byte-identical

## What did not, and why

| idea | result | reason |
|---|---|---|
| "round 2 is serial" | false | it is already parallel; I read the timer's line range, not the function's |
| "rolling seed saves 32x" | false | `rseed` is already O(1) -- one `w32` load and a shift |
| Bloom filter on seed misses | byte-identical, **no gain** | the 84% of failed lookups were never the cost; the candidate walk is |
| reorder mismatch chunks to exit sooner | **17% WORSE** (99.8 s -> 117.0 s) | sequential prefetch is worth more than an earlier `mm>lim` exit; the modulo and scattered access cost more than they save |
| full-seed check in pigeonhole | byte-identical, **no gain** | the index key is only the low 32 bits of a 32-base seed (16 bases), and restoring the full check rejects candidates that the mismatch loop was already discarding after one 32-base chunk |
| candidate concurrency K | already optimal | K=8/4/2 on 6M human reads: 564 / 549 / 596 s, identical archives. The "8 candidates thrash the cache" hypothesis does not survive the total. |

## The measurement error worth recording

**Atomic counters in a 1.93-billion-iteration loop inflated pigeonhole 6x** --
99.8 s instrumented against 15.5 s clean. Every timing taken from an
instrumented build was wrong, and the "pigeonhole is the bottleneck" reading
came from those numbers. The COUNTS are valid; the TIMES were not. Any future
funnel work here should count on a sampled subset, or count in thread-local
scalars merged once, never with atomics in the hot path.

## Standing conclusion

After five hypotheses, the encoder's hot stages are close to memory-bandwidth
bound and already parallel. The remaining large win is not a micro-optimisation
but an algorithmic change: replacing the level-sweep overlap search with an
FM-index (the all-pairs suffix-prefix problem, as SGA solves it). That changes
which chaining decisions are made, and therefore every Claim 1 archive size, so
it is a post-paper project rather than a tuning pass.
