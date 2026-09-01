# Stream-by-stream against PgRC2

PgRC2's archive is fully attributable: every stream goes through `writeCompressed`
or `CompressionJob::writeCompressedCollectiveParallel`, and CompressionJob
carries a `label`. Instrumenting those two functions in a copy of their tree
prints raw and coded bytes per stream, and the total accounts for **99.9%** of
their archive (3,114,267 of 3,116,030).

Measured on S. acidocaldarius, the one dataset where we lose.

| component | ours | PgRC2 | delta |
|---|---|---|---|
| read order | 1,414,637 | 1,398,530 | +16,107 |
| mismatch positions | 731,919 | 721,906 | +10,013 |
| pseudogenome literal | 597,736 | 595,838 | +1,898 |
| mismatch symbols | 124,280 | 208,234 | **-83,954** |
| mismatch counts | 159,665 | 163,068 | **-3,403** |
| **references** | **94,935** | **26,388** | **+68,547** |
| strand / misc | 17,184 | 303 | +16,881 |
| **archive** | **3,140,692** | **3,116,030** | **+24,662** |

## The finding: our references cost 3.6x theirs

    ours    mem_triples 57,262 + mem_dstgap 11,929 + mem_len 23,185 + mem_rc 2,559 = 94,935
    PgRC2   offsets 11,133 + lengths 5,316 + RC info 9,939                          = 26,388

This is the single largest component of the deficit and it is larger than the
deficit itself: without it we would win this dataset by ~44 KB. Every other
stream is within 2% or we are ahead.

Note what this overturns. mm_pos was the obvious suspect -- it is 23.4% of our
archive here. But PgRC2 spends 23.2% on the same thing (721,906 vs our 731,919,
a 1.4% difference), so mismatch positions are not where this is lost. Both tools
pay about the same, and we already code that stream below its order-0 and
count-conditioned entropy.

Where we are genuinely ahead: **mismatch symbols, 124,280 against 208,234** --
we spend 40% less than they do, because of the exclusion-model coder.

## Why their references are cheaper

Their offsets cost 11,133 B against our 57,262 + 11,929 = 69,191 for the same
role. Two structural differences, both already suspected and now quantified:

1. **No destination stream.** A MATCH_MARK inside their literal codebook says
   "a match starts here", so destinations cost nothing beyond the marks. Our
   mem_dstgap is a separate stream because our literal carries no in-band
   marker.
2. **Fewer, longer matches.** Their lengths stream is 5,316 B against our 23,185,
   which is consistent with substantially fewer references, not merely better
   coding of the same number.

## Read order: they are 1.1% ahead, and it is not a coding gap

Their order stream is a single block: raw 2,067,616 B = 516,904 reads x 4 B
exactly, coded to 1,398,530 = 21.65 bits/read. Ours is 1,414,637 = 21.89
bits/read. We measured our pos_abs at 0.996x its information bound, so this
1.1% is not slack in the coder -- it is a slightly different decomposition, and
it is the smallest of the three losses.

## Target order, by size

1. **References, +68,547.** Larger than the whole deficit. An in-band match
   marker removes mem_dstgap outright; the length stream ratio (4.4x) points at
   match count, not coding.
2. **Mismatch positions, +10,013.** Only 1.4% apart. Little to win.
3. **Read order, +16,107.** At 0.996x its bound; the gap is structural.

## What we should NOT target

- **mm_pos coding.** Both tools spend ~23% here and we are within 1.4%. The
  transpose (now implemented) bought 0.5% of the stream, and the earlier -5.4%
  estimate was an order-0 artefact.
- **Mismatch symbols.** We are already 40% ahead.
