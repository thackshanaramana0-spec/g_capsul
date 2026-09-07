# Claim 1: stream-level analysis

`Claim 1 LOCKED` meant its NUMBER was defended (+1.88% vs PgRC2), not that its
components had been understood. Sub-dividing one stream found -0.52% to -6.18%
on four datasets. This file records the method and what it has and has not
reached, so the next stream gets the same treatment rather than a fresh guess.

## 1. The composition (SRR29296997, 460,501 reads)

Before the split, 2,756,247 B:

| stream | coded B | share | bits/unit | vs order-0 |
|---|---|---|---|---|
| pos_abs | 1,264,855 | 45.9% | 21.97 b/pos | 0.82x |
| literal | 654,371 | 23.7% | 1.663 b/base | 0.86x |
| mem_triples | 448,192 | 16.3% | | 0.32x |
| mem_len | 115,462 | 4.2% | | |
| mm_pos | 93,233 | 3.4% | | 0.75x |
| everything else | 180,134 | 6.5% | | |

## 2. The method that worked: look for a MIXTURE

A stream that carries two populations with different statistics always measures
"near its bound", because the bound is computed on the mixture. The test is to
split by a structural property and compare entropies.

`pos_abs` carries positions from two pseudogenome regions:

    main-region    371,009 positions   delta entropy 18.45 b   (random)
    second-region   89,492 positions   delta entropy  1.03 b   (93.7% of
                                        consecutive deltas are exactly the
                                        read length -- appended, not placed)

An 18x entropy gap inside one stream. Splitting it (main raw uint32,
second-region zigzag-varint deltas, plus a region bitmap) gives:

    SRR29296997  -6.18%   (19.4% of reads in the second region)
    SRR554369    -1.94%   ( 6.4%)
    SRR40271341  -1.81%   ( 9.1%)
    ERR552797    -0.52%   ( 2.8%)

**The gain tracks the second-region fraction.** That dose-response is the
evidence the mechanism is real rather than a per-dataset fluke.

## 3. Where the same test says NO

Applying it honestly matters more than applying it everywhere.

- **`literal`**: its two regions measure H0 1.917 vs 1.946 b/base -- effectively
  the same population. Splitting LOSES shared context: 713,472 B combined
  against 720,924 B split. Not taken.
- **`mem_triples`**: already region-aware. `refc::bound_sel` bounds each source
  by MAINEND depending on whether the reference is a second-region self-match,
  which is why it sits at 0.32x order-0. The region idea was already applied
  here -- it had simply never been applied to `pos_abs`.
- **`pos_region`** (the new flag): order-0 through order-8 entropy are all
  0.710 b/read, so there is no run structure. Conditioning on mismatch count
  reveals that second-region reads ALWAYS have zero mismatches (P(second)=0.000
  for every mm>=1 read, they are stored verbatim), but mm=0 dominates so it
  saves only 2.5 KB. Not worth the coupling.

## 4. What is genuinely at its bound

- **`pos_abs` main region**: 371,009 positions over a 4,702,229 universe.
  Set-of-positions bound 231 KB + permutation bound 772 KB = 1,003 KB; we ship
  1,010 KB = **1.007x**. The permutation is 77% of it and measures as random
  (51.8% ascending adjacent pairs, 50% = random), so it is irreducible. Note
  the repo had already measured that 96.8% of chain tails have ZERO candidates,
  so there is no tie-breaking freedom to correlate placement with file order.
- **`literal`**: 1.663 b/base, and our coder beats `xz -9e` on it by 8%
  (654,371 vs 713,472 B).

## 5. Open, in order of size

1. `literal` at 25% -- a DNA-specific context model (order-16 mixing) is the
   only untried lever; general-purpose coders are already being beaten.
2. `mem_len` 4.5% and `mm_pos` 3.6% -- not yet sub-divided at all.
3. The dedup threshold is a hardcoded 0.15, which violates the project's own
   rule against fitted constants. The coding-cost break-even is
   `p * log2(PGLEN/N) > H(p)`, which for typical PGLEN/N ~ 10 puts p* near
   0.25, not 0.15 -- so datasets with 15-25% duplicates may be turning dedup ON
   where it costs. No locked dataset is currently known to sit in that band;
   worth checking before relying on it.
