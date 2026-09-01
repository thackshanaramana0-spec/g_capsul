# Do we need to beat PgRC2's 3-way split? Measured answer: NO.

Four levers tested. One won and shipped; three were refuted by measurement.

## The question

PgRC2 splits its pseudogenome three ways -- hqPg (both-side-overlapped reads,
assembled), lqPg (the rest, mapped onto hqPg), nPg (N-reads) -- each
self-matched separately. We have main pg + second region. Is the third split
the reason we lose S. acidocaldarius by 2.46%?

## Answer: we already have the structural equivalent, and we beat them on it

Round 1 (division) is their hq/lq split: it labels reads overlapped on both
sides, and round 2 assembles only those. Measured on S. acidocaldarius, main
region surviving: **ours 2,368,666 vs theirs 2,450,403 -- we win by 81,737.**

## Refutation 1: per-region MINMEM

97% of our MEM references come from the second region (E. coli: 16,197 main vs
479,148 second, against PgRC2's ~4-5K total). That looked like a reference
explosion from applying one acceptance length to two statistically different
regions -- the same error B3 fixed on orig2uid.

Implemented and swept. **Refuted, monotonically:**

| MINMEM2 | E. coli size | delta |
|---|---|---|
| 24 (= main) | 7,855,707 | 0 (byte-identical baseline) |
| 32 | 7,886,359 | +30,652 |
| 45 | 8,024,521 | +168,814 |
| 128 | 9,865,746 | +2,010,039 |

Those 431K references are not waste. Each removes enough bases to pay for
itself several times over. Reverted.

## The real structure: positions dominate everything

At the true optima, coded stream sizes:

| | pos_abs | literal | mem_triples | mm_pos |
|---|---|---|---|---|
| S. acidocaldarius | **45.8%** | 19.2% | 3.2% | 22.7% |
| E. coli | **46.2%** | 21.2% | 10.9% | - |
| H. salinarum | **48.5%** | 25.8% | 19.0% | 2.8% |

Nearly half of every archive is read positions, at log2(pg_span) bits each. Our
pre-removal pg is 1.5-2.4x theirs (E. coli 34.9M vs 22.6M; S. acidocaldarius
9.8M vs ~3.8M), because our second region is huge (E. coli 25.2M against their
lqPg of 1.5M). That span penalty is ~87K bytes on S. acidocaldarius, which is
essentially its entire 76,559 byte loss.

## Refutation 2: splitting the position stream by region

86% of E. coli reads sit in the main 9.7M and pay log2(34.9M) bits. Splitting
by region should save 170,358 B against the uniform bound.

**Refuted: we already pay less than the split bound.**

    uniform bound            3,870,083 B
    region-split bound       3,699,725 B
    what we actually pay     3,629,158 B   (0.94x uniform)

The byte-plane + LZMA coder already extracts more than the structural
decomposition would. Not implemented.

## Refutation 3: MINOV, the pg-span lever

Span drives position cost, so a tighter pg should pay. Swept on both:

    S. acidocaldarius: span 11.03M -> 8.96M (-19%) moves pos_abs only -25,952 B
                       (-1.8%), because log2 is slow; TOTAL is worst outside 87
    E. coli:           16 is optimal; 24 +2,170, 32 +11,793, 40 +6,590

MINOV is already at its optimum on both, and the curve is **not unimodal**
(E. coli 32 is worse than 40), so no search method is licensed for it either.

## What DID win: golden-section on MAXMAP

total_size(MAXMAP) IS unimodal on every dataset swept, which licenses
golden-section. 7 files, all smaller, -154,223 B. Chosen ceilings 18, 24, 9, 45,
8, 19, 24 -- none of them L/13 or L/5, the only two values the old grid could
reach. Aggregate vs PgRC2 +4.32% -> **+4.36%**.

## Conclusion

The S. acidocaldarius loss is not a missing third split. It is a genuine
tradeoff difference:

- **They** build a tight pre-removal pg and remove little (36%), so positions
  are cheap and literal is moderate.
- **We** build a loose pg and remove aggressively (74%), so literal is cheaper
  (we win the main region by 81,737) but every position pays ~1 extra bit.

On 9 of 10 datasets our tradeoff wins. On S. acidocaldarius -- low repeat, poor
second-region collapse -- theirs does, by 2.46%.

Copying their split would not fix it, and the two structural decompositions that
would have targeted it are both already beaten by our existing coder.

## What remains unexamined for size
- `mm_pos`, 22.7% of the S. acidocaldarius archive, never analysed at the
  current optimum (an earlier bucketing attempt lost at real scale).
- `orig2uid_vals`, 10% of E. coli after B3 split it; the values half is still
  786,787 B.
