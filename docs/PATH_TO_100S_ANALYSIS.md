# Can the archive path reach 100 s? Analysis and answer

**Answer: no, not without changing what is computed. The structural floor is
~200 s, and the run is already at 207.6 s.**

## Where the 207.6 s is, fully accounted

| item | time | share |
|---|---|---|
| `collapse_contigs` x2 | 51.8 s | 25% |
| pcluster scan | 47.6 s | 23% |
| read re-placement x2 | 37.9 s | 18% |
| decode reads+qual | 18.4 s | 9% |
| emit + XSNV | 12.7 s | 6% |
| kc_H_build | 11.8 s | 6% |
| seed index x2 | 4.2 s | 2% |
| rest | ~23 s | 11% |

## `collapse_contigs` is at its algorithmic floor -- measured, not argued

A standalone benchmark of the exact operation mix:

    insert 140M into unordered_set<uint64_t> : 52.3 s
    query   28M                              :  2.1 s
    TOTAL                                    : 54.4 s

`collapse_contigs` measures **51.8 s**. So the stage is essentially 100% hash
INSERTION, at the floor of what `unordered_set` does, and queries are only 4% of
it.

It also already inserts the minimum: only KEPT contigs' k-mers enter the set
(42% of contigs on call 1). 42% of 140M is 59M inserts, which at the benchmarked
rate is ~22 s against a measured 24.25 s -- consistent within 10%.

### Two fixes were built, measured, and REVERTED

| attempt | rationale | result |
|---|---|---|
| cache the canonical k-mers across both calls | the two calls compute identical k-mers | **5 s gained for +1.7 GB, wall got WORSE (210.8 vs 207.6)** -- proving the k-mer arithmetic was never the cost |
| Bloom pre-filter on `claimed` | the same lever took pcluster 112 -> 47.6 s | **223.4 s, WORSE** -- it guards 28M queries but costs 140M adds |

Both failed for the same reason, which the benchmark then explained: I was
optimising the 4% (queries) and adding overhead to the 96% (inserts).

### Why nothing else applies

* **Sorted array** -- `claimed` is built INCREMENTALLY (each contig queries the
  set as it stands, then adds to it), so there is a sequential dependency and no
  batch to sort up front. This is why the flat-array trick that gave 17x
  elsewhere cannot be used here.
* **Parallelism** -- same dependency: contig *i*'s decision depends on every
  prior kept contig.

The 140M inserts are not overhead; they ARE the algorithm.

## What reaching 100 s would require

Removing ~108 s. The only item large enough is the **second `build_substrate`**
(collapse 27.5 + re-place 21.6 + index 2.7 + downstream = ~64 s).

It exists because the SNV pileup wants aggressive collapse (`dup_frac` 0.45) and
the bubble/indel passes want mild collapse (0.92) -- the code records F1
measurements behind that split. Removing it changes indel results, which the
brief explicitly forbids.

Everything else is genuine computation on 4M reads and ~140M k-mers: entropy
decode, k-mer counting, and read re-placement that is already parallel.

## What WAS achieved

    671.0 s   FASTQ path, unoptimised (session start)
    432.0 s   archive path, unoptimised
    207.6 s   archive path, six structural fixes

**3.2x, F1 bit-identical throughout** (SNV 0.4539, INDEL 0.2613 -- every TP, FP
and FN value unchanged).
