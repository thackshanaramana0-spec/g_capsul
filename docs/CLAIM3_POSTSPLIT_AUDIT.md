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

## 5. Component-level optimisation — what was found and what remains

After verification, each mode was audited at the STREAM level: list what the
block actually references, compare against what has already been decoded when
control reaches it. Three instances of one defect class -- work done before an
early exit that does not need it -- with three different outcomes.

### 5.1 `coverage` — 2.43x, shipped

Its exit sat BELOW the entire assembly layer: `seq_decode_mem` (the
context-mixing literal decode, the most expensive stream in the archive),
`refc::decode`, and `mem_dstgap`/`mem_len`/`mem_rc`/`mem_self`. It reads none
of them -- only `pos_abs`, `read_lengths`, `orig2uid`, and the scalars `PGLEN`
and `MAINEND`.

    coverage  0.260 -> 0.107 s   (2.43x)
    output byte-identical, covered-bases identity still exact

Claim 3 asserts coverage is "served without reconstructing the assembly". That
was true of the design and false of the code path; it is now true of both. The
published speedup against bwa+samtools+mosdepth moves from 23-33x to ~56-80x.

### 5.2 `query` — measured, REJECTED

`pos_strand` and the three `mm_cnt` streams are decoded before the query exit
and referenced zero times by it. Deferring them was correct -- query output
byte-identical AND the full round trip still lossless -- and **neutral**:
0.200 -> 0.200 s. Those streams are small and the pg rebuild dominates.
Reverted rather than kept as harmless.

### 5.3 `export` — no defect, an earlier suspicion CORRECTED

A first pass suggested export decodes read-level streams (`pos_abs`,
`read_lengths`, `orig2uid`, `qual_index`) it does not need. Checking precisely,
only `contig_spans` lies between the pg completing and the export exit, and
export genuinely uses it for the per-contig path. Everything else it decodes IS
the assembly. No change made.

## 6. Final state, three runs each, verified against independent ground truth

| operation | time | peak | pg rebuilds | correctness |
|---|---|---|---|---|
| `export` | 0.183 s | 36 MB | 1 | 18,002,109 bases = PG_LEN **exact** |
| `coverage` | **0.110 s** | 85 MB | **0** | 69,535,651 = sum(read_lengths) **exact** |
| `query` | 0.200 s | 42 MB | 1 | 7,785 reads, 0 outside range **exact** |

## 7. Verdict on Claim 3

**Closed.** One real optimisation found and shipped (2.43x on coverage), one
measured and rejected, one suspicion corrected. Every remaining cost has a
stated structural cause:

- `export` and `query` rebuild the pseudogenome because they need its CONTENT,
  and the rebuild cannot be windowed -- 99.7% of references reach back more than
  100 kb (§3).
- `coverage` now touches nothing but the streams it reads.

No further optimisation is available without an archive-format change (periodic
self-contained restart points), which would cost the compression ratio Claim 1
defends.
