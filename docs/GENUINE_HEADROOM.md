# Genuine headroom — measured against a verified baseline

Every number here is taken from an archive that has been decompressed back to
byte-identical reads with the archive as the ONLY input (strace-confirmed: the
decoder opens the .capsule and nothing else).

## Verified baseline, 7/7 LOSSLESS end to end

| dataset | archive | PgRC2 | margin |
|---|---|---|---|
| H. salinarum | 2,788,602 | 3,050,477 | +8.58% |
| E. coli | 8,199,823 | 8,864,420 | +7.50% |
| L. major | 27,861,244 | 28,272,652 | +1.46% |
| P. aeruginosa | 8,961,752 | 9,043,181 | +0.90% |
| S. aureus | 13,506,922 | 13,595,003 | +0.65% |
| P. falciparum | 17,118,878 | 17,219,695 | +0.59% |
| S. acidocaldarius | 3,144,440 | 3,114,782 | **-0.95%** |
| **aggregate** | **81,581,661** | **83,160,210** | **+1.90%** |

## Where the bytes are

| stream | S. acidocaldarius | H. salinarum |
|---|---|---|
| pos_abs | **45.0%** | **45.4%** |
| mm_pos | 23.4% | 4.8% |
| literal | 19.0% | 22.7% |
| mem_triples | 1.8% | 15.2% |

## 1. pos_abs is CLOSED — 45% of the archive has nothing left

Measured on S. acidocaldarius: 516,904 positions coded in 1,414,637 B = 21.89
bits each, against a uniform floor of log2(6,757,000) = 22.69 bits. We are
already BELOW the naive bound.

Decomposed properly, the information content is:

    sorted gaps (WHERE the reads sit)     287,433 B   gap entropy 4.45 bits
    permutation (WHICH order they go in) 1,133,108 B   log2(516,904!)
    total bound                          1,420,541 B
    we pay                               1,414,637 B   ratio 0.996

**There is no coding headroom in the largest stream.** Any further gain must
change what is stored, not how.

## 2. The real structural opportunity: read order is 36% of the archive

The permutation term is 80% of pos_abs and 36% of the entire archive. It is the
price of preserving the original read ORDER, and it is irreducible while we
promise to.

PgRC2 makes this optional -- `-o` preserves order, and without it they reorder
freely. All our comparisons use their `-o` mode, so this is apples to apples
today. But it means an order-free mode would remove ~1.13 MB from S.
acidocaldarius alone. Both tools would shrink; the comparison would need
re-running in that mode on both sides.

This is the single largest quantity in the archive and it is a PRODUCT decision,
not a coding one.

## 3. Second-region self-match -- quantified, implemented, currently LOSSY

The second region is matched CROSS-wise against the main pg only, never against
itself. On SARS-CoV-2 that removes 2.2% of a 14.7 MB region while 74.3% of its
32-mers are repeat occurrences. Two implementations exist; both give
SARS-CoV-2 -94,279 to -108,432 and E. coli +61,792 to +95,373, discriminated by
how much the cross-pass already removed (97.8% / 39.2% / 12.7% unremoved).
Reverted because the archive it produces is LOSSY. See
docs/SECOND_REGION_SELF_MATCH.md.

## 4. mm_pos on S. acidocaldarius -- 23.4%, never analysed at the current optimum

735,667 B, the second largest stream on that file and 5x its share on
H. salinarum. An earlier bucketing attempt was measured at 100k reads and lost
at full scale; it has not been revisited since the coder work or since the
optimum moved to MAXMAP=49. Unknown, not closed.

## 5. mem_triples on H. salinarum -- 15.2%, 8x its share on S. acidocaldarius

423,445 B against 57,262. The destination stream (mem_dstgap) exists only
because our literal carries no in-band match marker; PgRC2 pays nothing for
destinations because a MATCH_MARK rides inside their literal codebook. Adopting
that removes a stream outright.

## What is NOT headroom, measured

- **pos_abs coding**: 0.996x its bound. Closed.
- **The 3-way pg split**: we already beat their main region (2,368,666 vs
  2,450,403 surviving on S. acidocaldarius). Four structural variants tested and
  refuted -- see docs/DO_WE_NEED_THEIR_3WAY.md.
- **MINOV / MAXMAP tuning**: golden-section already finds the true optimum, and
  the objective is unimodal in MAXMAP but NOT in MINOV.
- **Chain matching order**: 96.8% of tails have zero candidates at a level, so
  there is no contention to resolve.

## Honest ranking

1. **Order-free mode** -- 36% of the archive, but changes the product's promise
   and requires re-benchmarking both tools.
2. **Second-region self-match** -- large on amplicon-shaped data, implementation
   exists, needs one correctness bug found.
3. **mm_pos** -- 23% of the worst file, genuinely unexamined.
4. **In-band match marker** -- removes the mem_dstgap stream, 15.2% of
   H. salinarum's archive is mem_triples+mem_len.
