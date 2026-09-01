# Final headroom, all three axes, measured

## Where we stand (7/7 lossless, both tools measured fresh on an idle machine)

| axis | ours | PgRC2 | |
|---|---|---|---|
| size | 81,631,156 | 83,192,412 | **+1.88%** |
| compress | 106.2 s | 68.3 s | 1.6x slower |
| compress RAM | 1032 MB | 371 MB | 2.8x |
| decompress | 40.1 s | 4.0 s | **10.0x slower** |
| decompress RAM | 1563 MB | - | higher than our own compress |

## SIZE: closed. No headroom.

Every stream measured against a bound:

| stream | verdict |
|---|---|
| read order | 0.996x its information bound |
| mismatch positions | below order-0 AND count-conditioned entropy |
| read strands | 16,314 B vs an order-0 bound of 16,116; orders 1-8 give 16,095 |
| references | per-match within 11% of PgRC2, cheaper on lengths |
| mismatch symbols | we are 40% ahead |
| MINMEM / MAXMAP | both verified interior optima on the one losing dataset |

Two structural attempts this session, both refuted by measurement: the
second-region self-match (correct now, but LZMA already had the redundancy at
0.24 bits/base against 21.5 bits per reference) and PgRC2's both-side-overlap
admission rule (worse on 3/3). Roughly 44% of every archive is the cost of
preserving read order, which is irreducible while we keep that promise.

## DECOMPRESS: this is where the headroom is, and it is large.

Split of the 40.1 s:

| | time | share | peak RAM |
|---|---|---|---|
| capsule_d (C++) | 7.95 s | 20% | 305 MB |
| **decode_105.py (Python)** | **32.19 s** | **80%** | **1563 MB** |

**80% of decompression is a Python script**, and it sets the memory peak --
1563 MB, higher than our compress peak of 1032 MB. It iterates every read in
the file in the interpreter.

Ported to C++, decompress lands near 7.9 s against PgRC2's 4.0 s: from 10.0x to
roughly 2x. Nothing else on any axis is worth this much.

Already taken this session: `seq_decode_mem` decoded its chunks sequentially
even though the chunk table needed to parallelise them already existed (it was
added so the stream could be decoded at all). Decoding them in parallel took
L. major's C++ half from 11.10 s to 3.33 s, byte-identical, +30 MB.

## COMPRESS: some, in the assembly.

| phase | share |
|---|---|
| greedy sweep | 43% |
| mismatches + literal | 23% |
| coding | 20% |
| MEM + positions | 15% |

The coding phase was Amdahl-bound on its single largest job on all 7 datasets;
chunking that stream cut total compress 11% for 0.065% size. The sweep, at 43%,
has never been attacked for speed -- only for correctness.

## RAM: bounded by per-read state.

corr(read count, peak RAM) = **+0.983**, far above corr(total bases) = +0.883.
At 100 bp the per-read arrays are rpk 32 B + woff 8 + rlen 2 +
nxt/prv/ovl/ch_h/ch_t 20 + ppos 8 + pent 8 + prc/readMM 2 = 80 B/read, so
L. major's 4,739,289 reads cost 379 MB before any index or coder. woff and ppos
are u64 where u32 would serve, which is ~38 MB, or 4% -- real but not
structural.

## The order of work, by value

1. **Port decode_105.py to C++.** 80% of decompress and the memory peak. Takes
   the worst axis from 10x to ~2x.
2. **The greedy sweep**, 43% of compress, never optimised for speed.
3. **Nothing on size.** It is measurably exhausted within this scope.
