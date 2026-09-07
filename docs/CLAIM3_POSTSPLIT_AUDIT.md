# Claim 3 — re-verified after the `pos_abs` region split, and why `query` is bounded

The `pos_abs` split (see `CLAIM1_STREAM_ANALYSIS.md` §2) changed the stream that
**two of Claim 3's three operations read**. Claim 3's own history contains a
20% coverage undercount that passed every existing test, so wiring inspection
was not treated as sufficient. All three operations were re-measured against
independent ground truth.

## 1. Correctness after the split — exact, not approximate

| operation | check | result |
|---|---|---|
| `export` | bases emitted vs `PG_LEN` from `pg_params.txt` | 18,002,109 = 18,002,109 **exact** |
| `coverage` | covered bases vs `sum(read_lengths)` over all reads | 69,535,651 = 69,535,651, **difference 0**, and 0 bases legitimately clipped at the pg boundary |
| `query` | returned reads vs an independent scan of `pos_abs`+`read_lengths` | 7,785 = 7,785 **exact**, and 0 returned reads fail to overlap the requested range |

The `coverage` identity is the strongest of the three: it re-proves the historical
undercount bug stays fixed, because a single dropped duplicate would break the
equality.

**One harness error worth recording:** the first check reported `query` returning
0 reads. `query` emits FASTA (`>`), and the checker counted FASTQ (`@`). The
operation was correct throughout. This is the fifth harness-vs-code confusion in
this line of work -- always confirm the checker before believing a failure.

## 2. The architecture claim, confirmed by measurement

| operation | peak RSS | pg rebuilds |
|---|---|---|
| `coverage` | **4 MB** | **0** |
| `export` | 36 MB | 1 |
| `query` | 36 MB | 1 |

`coverage` runs in 9x less memory than the other two *because* it is hoisted
above the pseudogenome rebuild -- it needs only `pos_abs` and `read_lengths`.
That is the "the compressor already stored the structure this operation needs"
claim made concrete rather than asserted.

## 3. Why `query` is 1.62x and cannot simply be made faster

`query` returns 1.69% of reads (7,785 of 460,501) for a 100 kb window, and needs
only `pg[qa..qb]` -- 100,200 of 18,002,109 bytes, **0.56%** of the pseudogenome.
It nevertheless rebuilds 100% of it. That looks like an obvious ~178x
inefficiency, and the natural fix is a range-limited rebuild.

**It is not available, and the reason is structural.** The pseudogenome is built
by SELF-REFERENTIAL copies:

```cpp
memcpy(&pg[d], &pg[src[i]], L);        // source is an EARLIER pg region
```

so the content at `d` depends on `src[i]`, which may itself have been produced by
an earlier reference, transitively. Measured over all 175,401 references on
SRR29296997:

    backward distance dst-src:  median 687,834,479   p90 3,351,183,573
    references reaching back more than 100 kb:  104,750/105,057 = 99.7%

The dependency graph is global. `pg[qa..qb]` cannot be materialised without
resolving references scattered across the whole pseudogenome, so no windowed
rebuild exists without changing the ARCHIVE FORMAT -- e.g. periodic
self-contained restart points, which would cost compression ratio, the axis
Claim 1 defends.

So `query`'s modest 1.62x is a **property of the compression design**, not an
unoptimised code path. It should be reported that way: `query` earns its place
on SELECTIVITY (59x fewer reads returned here, 132x on E. coli) and on being
reference-free coordinate retrieval that no competitor offers -- not on speed.

## 4. Status

Claim 3 is intact after the split. All three operations verified exact against
independent ground truth, the architectural distinction between them is visible
in the measurements, and `query`'s known weak result now has a measured
structural explanation rather than an open question.
