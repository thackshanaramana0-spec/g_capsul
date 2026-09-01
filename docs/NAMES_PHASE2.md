# Phase 2 — names / read-ID column, wired into the CAPSULE archive

## What shipped

`include/names_coder.h` (namespace `nmc`, following the project's per-coder
convention: `mmc`, `refc`, `pgc`), wired into `106_inprocess.cpp` as three
additive archive streams and inverted in `capsule_decode.cpp`.

Gated behind `CAPS_NAMES=1`. With it unset the archive is **byte-identical**
to the pre-change binary (verified by `cmp` against a build of `HEAD`), so the
locked Phase-1 result is provably untouched.

## Result — exact, both competitors, same protocol

Names cost is isolated by differencing against a constant-name control file
(identical SEQ/QUAL, names replaced by `@r`), so no per-field stats rounding
enters the comparison:

| dataset | SPRING | Genozip | ours | vs SPRING | vs Genozip |
|---|---|---|---|---|---|
| DRR976266, 1 M names | 2.181 | 1.906 | **1.779** | **-18.4%** | **-6.7%** |
| SRR2584863 (E. coli), 500 K | 5.693 | 2.708 | **2.710** | **-52.4%** | +0.075% |

We beat SPRING decisively on both and Genozip on one; on E. coli we are
1,011 B behind Genozip on a 1.35 MB stream. Disclosed, not rounded away.

## Both verification gates, actually met this time

1. **Lossless from the real archive.** `capsule_decode` runs as a separate
   process with only the `.capsule` file and reproduces both the sequence
   column and the name column byte-for-byte. The earlier stage-86 round trip
   decoded its own in-memory buffers -- structurally the same shortcut as
   decoding the dumped intermediate streams, which is what hid the
   incomplete-`mem_triples` archive.
2. **Accounting exact.** Parsing the container and summing header +
   per-stream overhead + payloads gives 4,852,586 B against 4,852,586 B on
   disk, zero bytes unattributed.

## Four things a mechanical port of stage 86 would have got wrong

1. **Stage 86's dictionary header was never serialized** -- `entries*4`
   arithmetic printed with a `~=`. Real serialization (sorted, delta-varint,
   then through `best_encode` like every other stream) is 18,912 B where the
   estimate said 25,356 B.
2. **Stage 86 has no block index at all**; its round trip re-reads the
   original file to learn each block's name count. A decoder holding only the
   archive cannot. That cost -- a real stream -- was simply missing from its
   1.790 B/name. This is the third instance of this project's recurring bug
   class (mm_pos/mm_count omitted from totals; mem_triples' dst/len/rc dropped
   as "diagnostics only").
3. **Feeding the coder headers captured during the main FASTQ parse** would
   hold every header resident (~5 GB on T. cacao) and defeat the bounded queue
   the whole stage-70+85 merge exists to provide. Both passes stream the ID
   column straight from the file instead; peak RSS 118 MB on 5 M names.
4. Decode streams block-by-block to a `FILE*` rather than returning every
   name in a vector, keeping the decoder's memory bounded too.

## Measured and REJECTED (do not retry)

- **Bit-length-class + mantissa coding of the zigzag delta**, replacing the
  two adaptive byte models: worse on both files (+0.081 and +0.008 B/value on
  E. coli, +0.092 / +0.147 on DRR). The two byte models already capture the
  structure.
- **Conditioning the coordinate models on tile change** (reads come off a
  flowcell in raster order, so coordinates reset at a tile boundary): real
  signal, negligible size -- 155 B on E. coli, 572 B on DRR, against a
  8,973 B coordinate gap. Refuted before implementing.
- **File-constant token elision.** 23 of 27 tokens on E. coli and 24 of 28 on
  DRR are byte-identical in every read, and each still spent a `token_type`
  symbol per read. Implemented and kept (it is strictly correct and free),
  but it is worth only **971 B**, not the ~43 KB predicted: the adaptive model
  had already driven those ID_MATCH symbols to ~0.0007 bits. The estimate that
  motivated it assumed 0.03 bit/symbol and was wrong by 40x.
  The mechanism generalizes past a fixed token count -- an index is elided
  only if it is invariant AND present in every read, so `ID_END` still
  terminates variable-length headers. E. coli needs this: it carries two
  flowcells (`244:H73TDADXX` 27 tokens, `245:H73R4ADXX` 29 tokens) and a
  fixed-count rule silently disabled the whole mechanism there.

## Where the remaining 1,011 B on E. coli is

Decomposed by freezing the X/Y coordinates and re-measuring both tools
against the same control:

| | Genozip | ours | delta |
|---|---|---|---|
| non-coordinate fields | 8,431 | **2,222** | **-6,209** |
| X + Y coordinates | **1,345,325** | 1,354,298 | +8,973 |

We win the structural part by a wide margin and lose the coordinates. The
cause is **not** the representation -- the zigzag-delta stream's own order-0
entropy is 2.659 B/name, below Genozip's achieved 2.691 -- it is per-block
adaptive model warm-up, paid once per block. Block size sweep, E. coli:
1,356,520 B at 100 K/block, 1,354,790 at 250 K, 1,354,147 as a single block.
Genozip's own `--STATS` reports 18.0 MB vblocks against our 6.5 MB, i.e. it
pays the warm-up far less often.

Default block is therefore 250 K names (~16 MB, matching Genozip's vblock
granularity). Both datasets improve; peak RSS stays at 118 MB. The remaining
slack is the warm-up that a single block would remove, which is not available
without giving up parallelism and bounded memory.

**Open, honest:** priming each block's models from pass-1 statistics would
remove that warm-up while keeping blocks independent -- pass 1 already reads
the whole file. Not built: the priming histograms have to travel in the
header, and on a 500 K-name file that header would cost more than the
2,373 B it saves. It pays off as block count grows, so it is the right next
step for large files specifically, and it should be judged there.
