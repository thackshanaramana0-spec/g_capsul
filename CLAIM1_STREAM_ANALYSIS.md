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

## 6. Verification, and a pre-existing issue found on the way

**The split is LOSSLESS through the archive on every dataset tested**, at the
production parameter set (`3 16 16 22 16 16 1 24 64 1`), comparing decoded
reads against the original FASTQ's sequence column:

    SRR29296997   460,501 reads   LOSSLESS
    SRR40271341   386,537 reads   LOSSLESS
    ERR552797     590,692 reads   LOSSLESS

Plus, on SRR29296997: decoded positions byte-identical to the pre-change
archive, and every other stream byte-identical.

### 6.1 A pre-existing LOSSY result that is NOT this change

`scripts/verify_lossless.sh` reports **LOSSY** for SRR40271341 and ERR552797 --
and it reports LOSSY for the **BASELINE encoder too**, at the same parameters:

    SRR40271341   baseline LOSSY 3,814,744 B   |  split LOSSY 3,775,454 B
    ERR552797     baseline LOSSY 4,717,476 B   |  split LOSSY 4,707,212 B

Identical verdict on both binaries, so the split did not cause it. Two things
distinguish that check from the archive round trip above, and both matter:

1. **It does not test the archive.** Per this repo's own CLAUDE.md section 6.1,
   `verify_lossless.sh` decodes the DUMPED intermediate streams through
   `decode_105.py`, not the container. It exercises the algorithm but not the
   entropy layer.
2. **It selects its own parameters.** It derives a candidate grid from read
   length and chose `MAXMAP=60 MINOV=105` and `MAXMAP=60 MINOV=16` here, where
   the production encode uses `MAXMAP=11 MINOV=16`. So it is testing a
   configuration the production path does not use.

Read counts match exactly (386,537 / 590,692), so no reads are lost -- content
differs. That is the signature of the silent data-loss class recorded in
CLAUDE.md section 6.3.

**This needs its own investigation and is independent of the pos_abs work.**
Either `decode_105.py` is stale relative to the encoder for those parameters,
or a real defect exists at `MINOV=105`. Do not treat it as cleared just because
the archive path round-trips; and do not treat the archive path as suspect just
because this script fails -- they are different code.

## 7. `literal` (25%): investigated to the end, and it is DONE

The stated lever was "a DNA-specific context model". **That premise was wrong:
one already exists.** `run_chunk` in `include/seqpar_core.h` is a
context-mixing coder of the GeCo/XM class -- ten hashed context orders
(`ORD = {1,2,3,4,6,8,11,14,18,22}`), logistic mixing with weights learned per
(node, history) context, an APM/SSE refinement stage, a match model for
LZ-style repeats, and binary arithmetic coding at 2 bits/base.

### 7.1 It beats every general-purpose coder available

On the 3,148,818-base literal, ours codes 654,371 B (1.663 b/base):

| coder | bytes | vs ours |
|---|---|---|
| **ours (CM)** | **654,371** | -- |
| xz -9e | 713,472 | +9.0% |
| zstd -19 --ultra --long | 727,204 | +11.1% |
| bzip2 -9 | 751,246 | +14.8% |

### 7.2 Table size: refuted twice, at two scales

`TBITS=16` gives each order 65,536 slots while `ORD` reaches 22, so the high
orders collide heavily. That looked like an obvious defect. It is not worth
fixing:

| TBITS | archive (3.1 Mbase literal) | delta | enc peak |
|---|---|---|---|
| 16 | 2,585,780 | -- | 318 MB |
| 18 | 2,584,917 | -863 | 332 MB |
| 20 | 2,583,757 | -2,023 | 422 MB |
| 22 | 2,582,909 | -2,871 (**0.11%**) | **795 MB** |

The first explanation was that the model is DATA-starved rather than
table-starved -- 3.1 Mbases cannot populate order-11+ statistics. That makes a
falsifiable prediction: the gain should grow with literal size. **It was tested
at 10x scale and the prediction FAILED:**

| TBITS | archive (chr20 4M subset) | literal | enc s | peak |
|---|---|---|---|---|
| 16 | 34,775,005 | 15,625,547 | 52.69 | 1,265 MB |
| 22 | 34,730,681 | 15,581,223 | 65.99 | 1,660 MB |

**-0.13% for +25% encode time and +31% RAM.** Flat against the 0.11% at small
scale. So the high orders carry little information on this data at any size
tested, and `TBITS=16` is the correct setting on the three-axis trade. Reverted.

### 7.3 Chunk count: already at its knee

The chunk count is a COMPRESSION parameter, not only a threading one, because
the model tables are allocated per chunk and restart with each. This was
already measured and is documented in the source: E. coli literal at 12 chunks
1,809,073 B, 4 chunks 1,804,223 B, 1 chunk 1,801,710 B. Four costs 2.5 KB
against the single-chunk optimum while staying parallel. `SEQT` tunes it.

### 7.4 Verdict

`literal` is at its practical bound. The model is the right class, its table
size is correct on the three-axis trade at two scales, and its chunking is at
the documented knee. No further work here without a fundamentally different
model, and general-purpose coders are already 9-15% behind.

## 8. The reference streams (6.8%): measured, rejected

`mem_len`, `mem_dstgap` and `mem_rc` are coded with plain `best_encode`, with
none of the region awareness `mem_triples` gets from `refc::bound_sel` -- so
they looked like the same opportunity as `pos_abs`. Instrumented before
building anything:

    main refs    n=  3,916   H(len) = 9.863 b   P(rc) = 0.760
    second refs  n=157,993   H(len) = 6.741 b   P(rc) = 0.383
    combined                 H(len) = 6.910 b

The distributions genuinely differ -- but the population is **97.6% second-
region**, so separating them recovers almost nothing:

    potential mem_len saving: 1,892 B = 1.4% of the stream = 0.07% of archive
    potential mem_rc saving:  ~200 B

**Rejected.** A format change is not worth 0.07%.

This is the third time the mixture test has said no (`literal`, `mem_triples`
already region-aware, and now the reference streams), against one time it said
yes (`pos_abs`, an 18x entropy gap across a 19%/81% split). The test is
discriminating on two things at once: the populations must differ AND be
reasonably balanced. `pos_abs` had both; these have only the first.

## 9. Claim 1: where it stands

Four streams covering **82% of the archive** have now been examined at the
component level:

| stream | share | verdict |
|---|---|---|
| `pos_abs` | 46% -> 40% | **WON** -0.52% to -6.18%, verified lossless |
| `literal` | 25% | at bound; CM beats general coders 9-15%, table size refuted at two scales |
| `mem_triples` | 17% | already region-aware, 0.32x order-0 |
| `mem_len`/`dstgap`/`rc` | 6.8% | measured, 0.07% available, rejected |

Untested remainder is ~10% (`mm_pos` 3.6%, `mm_sym` 2.0%, `mm_cnt` 2.3%,
`pos_region` 1.6%, `pos_sec` 0.7%). Even a uniform 20% win across all of it
would be ~2% of the archive, so the realistic remaining upside under this
decomposition is small.

**The one structural lever left is not a coder change: shrink the pseudogenome
itself.** CLAUDE.md already names second-region self-match as "the single
largest quantified opportunity", implemented twice, not shipped because the
archive comes out LOSSY (the RC path with a non-zero SRCBASE is the documented
first place to look). Its measured effect is dataset-dependent -- SARS -94,279
to -108,432 B, but E. coli +61,792 to +95,373 B -- so it needs the mismatch
economics reworked as well as the bug fixed. That is the next real piece of
work on Claim 1, and it is a larger one than anything in this file.

## 10. `literal`: five knobs swept, all at optimum -- and why nothing shipped

Continuing past the first "at bound" claim, every tunable in the model was
swept. All on SRR29296997, archive 2,585,780 B; the encoder is deterministic so
these sizes are exact, not noisy.

| knob | variants tried | outcome |
|---|---|---|
| table size `TBITS` | 16/18/20/22, at two data scales | 16 optimal; 22 buys 0.11% (small literal) and 0.13% (10x scale) for +477 MB and +25% time |
| chunk count | already measured in source | 4 is the knee, 2.5 KB from the single-chunk optimum |
| context orders `ORD` | lean6, nolow, dense14, hi_ext | baseline best; every variant +134 to +5,210 B |
| mixer context `MCB` | 4/6/8 | 4 best; +268 and +1,636 B -- finer selection dilutes the weights |
| SSE context `APMB` | 6/8/10/11/12/14/16 | 12 is the argmax, -1,423 B, 4/4 datasets positive |

### 10.1 The one that "won" was NOT taken

`APMB=12` is a **fitted constant** -- the argmax of a sweep. This project's own
standing rule (CLAUDE.md rule 1) forbids exactly that, and four datasets
agreeing does not turn an argmax into a formula. At a mean of 0.036% it would
not justify a format change on its own either. **Reverted, not shipped.**

The finding underneath it is worth keeping even though the constant is not: the
SSE stage behaves differently from every other component because it refines an
ALREADY-FORMED probability through a 33-bucket curve, so each of its contexts
needs far less evidence than a predictor does. That is why it had headroom where
the mixer and the context tables were saturated -- and it means the earlier
generalisation "this model is evidence-limited everywhere" was too broad. If an
SSE size is ever DERIVED from a measured input property, that would be
shippable. The argmax is not.

### 10.2 What the sweeps were for

Diagnosis, not a fix. Five independent sweeps landing on the shipped
configuration is what turns "literal is at bound" from an assertion into a
measurement. **`literal` is closed by evidence, with nothing shipped from it.**

The structural win in Claim 1 remains the `pos_abs` region split (section 2):
not a tuned constant, but a decomposition keyed on a measured property of each
read, separating two populations with an 18x entropy gap.
