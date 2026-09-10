# Making `query` cost track the answer, not the archive

**Status:** measured and planned, not implemented. The measurements below
overturn a conclusion this project had recorded as settled.

## 1. Where the time actually goes

`CAPS_QTIME=1`, E. coli (SRR2584863), 100 kb range, 26.96 MB pseudogenome:

    literal decode (seq_decode_mem)   0.242 s   55%
    pg rebuild (replay references)    0.129 s   29%
    reference streams decode          0.044 s   10%
    placement streams unpack          0.007 s    2%
    scan all reads + emit the answer  0.006 s    1.3%
                                      -------
    total                             0.47 s

**98.6% of a query is spent building a 27 MB pseudogenome in order to hand back
1.6 MB.** The answer itself is already free. Every optimisation must therefore
attack pg materialisation; nothing else is worth touching.

## 2. The blocker was mis-diagnosed

This project recorded that windowed rebuild is impossible because "99.7% of
references reach back >100kb, so no windowed rebuild is possible"
([[reimpl notes]], and repeated in the caller docs).

The reach-back figure is correct — measured on 240,762 references:

    median 13,428,096   p90 21,676,383   p99 25,537,480   max 26,941,182
    reaching back <1kb: 0.1%     <100kb: 0.2%

But reach-back distance is not the quantity that decides windowed rebuild.
**Closure size is.** A reference reaches 13 MB back yet carries only ~85 bytes
(20.58 MB of reference coverage over 240,762 refs). Measured closure, taking
the transitive dependency of a window to fixpoint:

    window   median closure   % of pg   intervals   depth
      1 kb          2,453 B     0.01%          29       8
     10 kb         19,182 B     0.07%         157       8
    100 kb        179,322 B     0.67%         181       8
      1 MB      1,728,600 B     6.41%      10,633      10

**Closure is ~1.7-1.8x the window and converges in ~8 rounds.** A 100 kb query
depends on 0.67% of the pseudogenome, not 100%. The conclusion "no windowed
rebuild is possible" confused *far* with *large* and should be retracted.

Reproduce: `CAPS_DUMP_REFS=refs.tsv capsule_decode query <archive> /dev/null 0-1000`
then the closure script in this document's commit.

## 3. What has to change, in order of the time it buys

### 3.1 Make `literal` block-addressable — buys 55%

`stages/capsule_decode.cpp`: `seq_decode_mem(S["literal"]...)` decodes the
whole literal stream before anything else can happen. literal is 23.7% of pg
(6.38 MB here) and is decoded in full for a query that needs 0.67% of it.

**Encoder** (`stages/106_inprocess.cpp`, literal emission): code literal in
independently-decodable blocks of B bases, and emit a new `literal_index`
stream of `(pg_offset_of_block_start, byte_offset_in_coded_stream)` per block.
**Decoder**: decode only the blocks intersecting the closure.

B is the whole tradeoff. Smaller B = finer granularity but more context resets
= worse ratio. Start at B = 256 KB (25 blocks here, index ~200 B) and measure
the ratio cost on the 4-file regression subset before going smaller. **If the
ratio cost exceeds ~0.05% this is not worth shipping** — Claim 1 is the
stronger claim and must not be spent on Claim 3.

### 3.2 Make references seekable by dst — buys 10%

`dst` is already monotonically ascending (it is delta-coded: `dst[i] = prev +
gaps[i]`), so no re-sorting is needed. Today all 240,762 references are decoded
to answer any query.

Add `mem_ckpt`: every 1024th reference, store `(dst, byte offset into each of
mem_dstgap / mem_len / mem_rc / mem_triples)`. ~240 entries, a few KB. The
decoder binary-searches to the first checkpoint at or before the window and
decodes forward. Closure resolution needs the same seek for each of its
intervals, which is why the checkpoint stride wants to be small-ish.

### 3.3 Replace the full replay with closure materialisation — buys 29%

Today the decoder replays all 240,762 references into a 27 MB buffer. Instead:

    closure(W):
      need = {W}; have = {}
      until need is empty (~8 rounds):
        for each interval I in need:
          for each reference R overlapping I:            # via 3.2
            map I's part of R back to R's source range
            add it to need unless already in have
        move need into have
      materialise: decode the literal blocks (via 3.1) covering `have`,
                   then replay references over `have` in dst order

Materialise into a sparse structure -- a sorted vector of
`(pg_start, len, offset_into_scratch)` -- and serve read emission from that
instead of a flat `pg` buffer. The read-emit loop changes from `pg.data()+a` to
a lookup in that structure; everything else in the query path is unchanged.

### 3.4 A read-position index — measure before building

The scan over all reads costs 0.006 s here (12.6 M reads on HG002 will be
worse, but placement unpack was only 0.007 s there too). Bucket read ids by
pg position (e.g. 64 kb buckets) only if measurement says it matters. **Do not
build this on the assumption that O(n) is bad; it is currently 1.3% of query.**

## 4. Expected result and how it will be judged

If closure work is proportional to closure size, a 100 kb query should touch
~0.7% of what it touches today:

    today          0.47 s
    projected      5-20 ms          (Genozip --head=100 is 0.19 s)

That would make query dominant on **time** as well as on capability, and would
change the claim from "we are the only one that can be asked" to "we are the
only one that can be asked, and we are also the fastest."

**Gates, in order:**
1. Full decode stays byte-identical (`cmp` on decoded reads, all 4 regression
   datasets). The full path must not regress to buy a query path.
2. Archive size regression < 0.05% aggregate, or the change is reverted. Claim 1
   outranks Claim 3.
3. Query output byte-identical to today's for the same range, on E. coli and
   HG002.
4. Only then, the timing claim.

## 5. Why this is worth doing

Query is currently our weakest table: Genozip extracts 100 reads in 0.19 s
against our 0.46 s, so the honest framing today is "slower, but the only one
that can be asked a locus question". Closing the time gap removes the only
concession in Claim 3 — and it does so by making the cost model right
(proportional to the answer), not by micro-optimising a constant.
