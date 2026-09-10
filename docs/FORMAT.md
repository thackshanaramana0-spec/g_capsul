# The G_CAPSUL archive format (version 2)

> **Status:** frozen 2026-09-10 at tag `v1.0.1-capsule` — code and results final.
> Authoritative numbers live in `../benchmark/results/`; verification status in
> [`../AUDIT.md`](../AUDIT.md). Where this file and a result file disagree, the result file wins.

Read directly from `read_capsule()` in `stages/capsule_decode.cpp` and the
emission sites in `stages/106_inprocess.cpp`. All integers are little-endian.
An implementation that follows this document can read a `.capsule` without
reference to our source.

## 1. Container

```
offset  size  field
     0     8  magic       "CAPSULE\0"
     8     2  version     uint16, currently 2. A reader MUST reject any other.
    10     8  pg_len      uint64, pseudogenome length in bases
    18     8  main_end    uint64, end of the MAIN region; [main_end, pg_len) is
                          the SECOND region
    26     4  minmem      uint32, minimum MEM length used at encode time
    30     2  n_streams   uint16
    32     -  streams     n_streams records, in the order below
```

Each stream record:

```
  1  name_len   uint8
  n  name       name_len bytes, ASCII, not NUL-terminated
  8  coded_len  uint64
  m  payload    coded_len bytes, entropy-coded (see section 3)
```

A reader must treat unknown stream names as skippable: the count and lengths
make the container self-delimiting. A file whose stream payloads do not all
read completely is malformed and MUST be rejected -- the reference decoder
returns non-zero, and the test suite asserts this for random data, a wrong
version, and truncation.

## 2. Streams, and who needs them

| stream | holds | required by |
|---|---|---|
| `literal` | pseudogenome literal bases, 2-bit coded | decompress, export, query, call |
| `mem_triples`, `mem_dstgap`, `mem_len`, `mem_rc`, `mem_self` | the self-referential copies that rebuild the pseudogenome | decompress, export, query, call |
| `mem_extmm_cnt`, `mem_extmm_pos`, `mem_extmm_obs` | mismatches applied during pg reconstruction | decompress, export, query, call |
| `pos_abs` | MAIN-region read positions (see section 4) | decompress, coverage, query, call |
| `pos_sec` | SECOND-region positions, zigzag-varint deltas | as above |
| `pos_region` | per-read region bitmap, prefixed with an exact read count | as above |
| `pos_strand` | one bit per read | decompress |
| `read_lengths` | uint16 per ORIGINAL read | decompress, coverage, query, call |
| `orig2uid_flags`, `orig2uid_vals` | duplicate-read aliasing | decompress, coverage, query, call |
| `mm_sym`, `mm_pos`, `mm_cnt`, `mm_cnt_flags`, `mm_cnt_vals` | per-read mismatches against the pseudogenome | decompress |
| `n_pos`, `n_indices`, `n_cnt` | positions of non-ACGT bases | decompress |
| `contig_spans` | contig boundaries within the pseudogenome | export (per-contig), call |
| `names_body`, `names_dict`, `names_index` | read identifiers | decompress with names |
| `qual_body`, `qual_index` | quality strings | decompress with quality |

**Two streams are written only when asked for at compress time.** An archive
without them is still valid and still decompresses losslessly; the operations
that need them refuse with a message naming the missing stream:

- `contig_spans` requires `CAPS_CALL=1`
- `pos_abs`/`pos_strand`/`read_lengths` require `DUMP_PERM=1`

## 3. Stream coding

Each payload is independently entropy-coded. The encoder selects per stream
from LZMA, PPMd, FSE, a range coder, and byte-plane variants, choosing the
smallest on a 40% probe of the stream (`probe_frac()` in
`include/coders_inproc.h`). Large streams are split into chunks of about 2 MB
with a chunk table, so the selection is per chunk. The chosen method is
recorded in the payload, so a decoder does not need to know how the choice was
made.

`literal` is coded by a context-mixing model (`include/seqpar_core.h`): ten
hashed context orders, logistic mixing, an APM/SSE stage, and a match model,
coded as two binary decisions per base.

`qual_body` is coded by fqzcomp (`thirdparty/htscodecs`), which is itself
searched over four strategies per block.

## 4. Why positions are split across three streams

`pos_abs` used to hold every read's position. Measured, it carried two
populations with an 18x entropy gap: MAIN-region positions are effectively
random (delta entropy 18.45 b), while SECOND-region reads are APPENDED
consecutively, so 93.7% of their consecutive deltas are exactly the read
length (delta entropy 1.03 b). One stream forced the model to straddle both.

Splitting them costs one bitmap and saves 0.52% to 6.18% of a sequence-only
archive, tracking the second-region fraction. A reader reconstructs the full
position array with `caps_join_positions()`: walk the bitmap, take main
positions in order from `pos_abs`, and accumulate zigzag-varint deltas from
`pos_sec` starting at `main_end`.

**Backward compatibility:** an archive without `pos_region` is pre-split, and
`pos_abs` is then complete on its own. The reference decoder detects this and
takes the old path unchanged.

## 5. What the format does not do

- No index for random access by read id. `query` is by pseudogenome
  coordinate.
- No windowed pseudogenome reconstruction. The pg is built by self-referential
  copies and 99.7% of references reach back more than 100 kb, so a coordinate
  range cannot be materialised without resolving the whole chain. Windowing
  would need periodic self-contained restart points, at a cost in ratio.
- No checksum over the payload. Truncation and corruption are caught
  structurally (stream lengths must all read completely), not cryptographically.
