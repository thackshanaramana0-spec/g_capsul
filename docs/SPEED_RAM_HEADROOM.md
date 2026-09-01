# Speed and RAM headroom — measured, and only what generalizes

Baseline: the verified build. 7/7 datasets compress to one archive and
decompress from that archive alone to byte-identical reads.

    compress    112.8 s total, 965 MB peak   = 2.55x slower, 2.63x heavier than PgRC2
    decompress   21.3 s total, 296 MB peak   = 5.3x faster than compressing

## Where the wall clock goes

Phase timing (CAPS_PHASE=1), P. falciparum and L. major:

| phase | P. falciparum | L. major |
|---|---|---|
| assembly + mapping + positions | 17.63 s (66%) | 23.35 s (48%) |
| mismatches + literal | 3.55 s (13%) | 8.75 s (18%) |
| coding | 5.36 s (20%) | 16.97 s (35%) |

## Finding 1 — the coding phase costs exactly ONE job. Verified 7/7.

| dataset | pool wall | longest job | equal? | speedup on 12 cores |
|---|---|---|---|---|
| halo | 1.03 | pos_abs 1.03 | YES | 3.13x |
| ecoli | 2.82 | pos_abs 2.82 | YES | 3.22x |
| sulfo | 1.93 | mm_pos 1.93 | YES | 2.50x |
| paeru | 2.92 | pos_abs 2.92 | YES | 2.68x |
| saureus | 4.07 | pos_abs 4.07 | YES | 4.30x |
| pfalc | 5.10 | pos_abs 5.10 | YES | 3.08x |
| lmajor | 16.55 | pos_abs 16.55 | YES | 1.77x |

On every dataset the pool wall equals the single longest job to within 0.06 s.
Every other job hides completely behind the biggest one, so **adding cores
cannot help** -- speedup is 1.77-4.30x on 12 cores and falls as the largest
stream grows. The bound is structural, not a scheduling accident: we parallelise
ACROSS streams and never WITHIN one.

The mechanism to fix it already exists in the codebase. `literal` is coded as
SEQT chunks with a chunk table (added so the stream could be decoded at all).
The largest stream has no such split.

### The cost is real and does NOT generalise

Measured by coding each dataset's actual pos_abs whole versus in k chunks:

| chunks | L. major | P. falciparum | H. salinarum |
|---|---|---|---|
| 1 | +0.000% (14.85 s) | +0.000% (4.63 s) | +0.000% (0.67 s) |
| 2 | +0.071% (5.55 s) | +0.383% (1.95 s) | +0.021% (0.36 s) |
| 4 | +0.095% (3.04 s) | +0.604% (0.90 s) | +0.056% (0.20 s) |
| 8 | +0.170% (1.34 s) | +0.658% (0.50 s) | +0.133% (0.12 s) |
| 12 | +0.225% (0.89 s) | +0.705% (0.35 s) | +0.159% (0.09 s) |

(time shown is the longest chunk, which is what the pool wall becomes)

P. falciparum pays **6x** what L. major pays for the same split, because the
context reset costs more where the positions carry long-range structure. Net
effect on the published margins:

| dataset | pos_abs share | size cost at k=4 | margin now | margin after |
|---|---|---|---|---|
| L. major | 55.6% | 14,716 B | +1.46% | +1.41% |
| P. falciparum | 30.2% | 31,247 B | +0.59% | **+0.41%** |
| H. salinarum | 45.4% | 710 B | +8.58% | +8.55% |

So this is a **speed/size tradeoff, not free headroom**. It buys roughly 25% of
L. major's wall clock for 0.05% of its archive, but costs P. falciparum nearly a
third of its winning margin. Any split factor must be derived from a measured
property of the stream, not fixed -- and on H. salinarum chunking is *slower*
sequentially (0.67 -> 0.82 s), so a size threshold is part of the rule.

Note also that chunking is faster even without parallelism on the large file
(14.85 -> 10.42 s sequential at k=4): LZMA is superlinear in block size.

## Finding 2 — peak RAM has two different owners

| | phase 1 peak | coding adds | total |
|---|---|---|---|
| L. major | 852 MB | +85 MB | 937 MB |
| P. falciparum | 349 MB | +434 MB | 783 MB |

peak = max(assembly floor, concurrent coder demand), and which one wins depends
on the dataset. This is why bounding coder concurrency looked like a 33% win on
P. falciparum (766 -> 511 MB, free) and was worth only 11% on L. major while
costing 5 s. **That was a knob, it did not generalise, and it has been removed.**

The coder side is LZMA encoder state, which is ~11.5x the dictionary, and the
dictionary is bounded by the stream size. So coder demand scales with the
largest streams and assembly demand scales with read count.

Assembly floor on L. major, 4,739,289 reads, accounted:

    rpk 2-bit packed reads (5 u64/read)   190 MB
    woff u64 per read                      38 MB
    pent u64 per read                      38 MB
    ppos u64 per read                      38 MB
    nxt/prv/ovl/ch_h/ch_t u32 x n          95 MB
    ment + mtab (MEM index, measured)     100 MB
    rlen/prc/readMM                        19 MB
    ------------------------------------------
    accounted                             518 MB
    measured                              852 MB
    unaccounted                           334 MB   <- not yet identified

Reads are already 2-bit packed (rbase decodes via rseed), so the obvious
saving is not available. `woff` as u64 per read is 38 MB where u32 would do,
but that is 4% of peak, not a structural win.

## IMPLEMENTED: chunking the largest stream (method 7)

Built and verified 7/7 lossless through the real archive.

| | size | wall | peak RAM |
|---|---|---|---|
| baseline | 81,581,661 | 112.77 s | 965 MB |
| chunked, concurrency 2 | 81,634,777 (+0.065%) | 100.27 s (-11.1%) | 997 MB (+3.3%) |

Margin vs PgRC2 +1.90% -> +1.83%. Speed 2.55x -> 2.26x slower.
L. major alone: 47.60 -> 37.91 s (-20.4%).

The design separates the SPLIT from the CONCURRENCY, because they have
different costs:

- Splitting is nearly free and helps on its own -- LZMA is superlinear in block
  size, so coding L. major's pos_abs as 4 pieces takes 10.42 s sequentially
  against 14.85 s whole, with no extra memory.
- Concurrency is what costs memory: each chunk in flight holds its own LZMA
  encoder state, ~11.5x its dictionary.

Measured on L. major, varying concurrency at a fixed split:

| concurrent chunks | wall | peak RAM |
|---|---|---|
| 1 | 41.69 s | 941 MB |
| 2 | 38.31 s | 989 MB |
| 4 | 35.56 s | 1118 MB |
| 12 | 36.30 s | 1441 MB |

Twelve is both SLOWER and 49% heavier than four: the gain saturates well before
the core count because the job stops being the pool's bottleneck, while memory
keeps climbing. Running all chunks at once cost +44.9% peak RAM across the suite
for only 2.9 points more speed than two does. The default is 2.

The split factor is derived from stream size, never fixed: below CHUNK_MIN (4 MB)
the stream is coded whole, which is why H. salinarum and S. acidocaldarius are
byte-identical to baseline (+0.000%) -- their streams are too small for splitting
to pay, and on H. salinarum chunking is slower sequentially.

## Honest ranking

1. **Assembly is 48-66% of the wall clock and has not been attacked at all.**
   Every speed finding above concerns the coding phase, which is 20-35%. This is
   the larger target and the least explored.
2. **Chunk the largest stream** -- removes a structural Amdahl bound proven on
   7/7 datasets, but costs real size and the cost varies 6x across datasets.
   Needs a rule keyed on stream size and measured context loss.
3. **334 MB of L. major's assembly peak is unaccounted for.** Identifying it is
   a prerequisite to any RAM claim, not a result in itself.

## What is NOT headroom

- **Bounding coder thread count.** Measured; a knob; does not generalise.
- **2-bit packing the reads.** Already done.
- **More cores for coding.** Provably useless while one job bounds the pool.
