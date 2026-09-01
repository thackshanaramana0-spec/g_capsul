# Where the size work stands, and what actually comes next

## Everything in the current scope is at its limit

Checked against bounds, not intuition:

| | verdict |
|---|---|
| read order | 0.996x its information bound |
| mismatch positions | below order-0 AND count-conditioned entropy |
| read strands | 16,314 B against an order-0 bound of 16,116; orders 1-8 give 16,095 |
| references | per-match within 11% of PgRC2 and cheaper on lengths |
| mismatch symbols | we are 40% AHEAD |
| MINMEM | swept, genuine interior optimum (24; worse at 20 and at 28) |
| MAXMAP | swept, genuine interior optimum (49; worse at 40 and at 60) |

Both parameters on our only losing dataset are at verified optima, and no coder
has slack. **The sequence+order scope is closed.**

## What predicts our margin

| dataset | second-region share | reads mapped | margin |
|---|---|---|---|
| H. salinarum | 70.5% | 21.9% | +8.58% |
| E. coli | 67.2% | 15.0% | +7.48% |
| L. major | 12.5% | 20.3% | +1.37% |
| P. aeruginosa | 27.0% | 27.4% | +0.89% |
| S. aureus | 65.8% | 20.5% | +0.64% |
| P. falciparum | 62.6% | 15.1% | +0.43% |
| **S. acidocaldarius** | **8.9%** | **51.3%** | **-0.83%** |

corr(second-region share, margin) = **+0.607**; corr(mapped fraction, margin) =
**-0.445**, over 7 datasets.

This **inverts** the thesis in DO_WE_NEED_THEIR_3WAY.md, which held that our
oversized second region was the problem ("that span penalty is ~87K bytes on
S. acidocaldarius, essentially its entire loss"). Measured now, a LARGER second
region goes with a BETTER margin. That analysis predates the coder work that
brought positions to 0.996x their bound, which is what removed the span penalty.

S. acidocaldarius is the extreme on both axes: it is the dataset where our
mapper places the most reads (51.3%, against 15-27% everywhere else), and
mapping is what generates the 983,565 mismatches that make its mm_pos 23.4% of
the archive. PgRC2 spends 23.2% on the same thing, so neither tool handles that
shape well; theirs simply handles it 0.83% better.

## The recommendation

**Stop optimizing sequence size.** Six of seven datasets win, aggregate +1.84%,
every stream is at its bound, and the one loss is at a verified parameter
optimum on the one data shape that suits their construction.

The real remaining work is scope, not compression:

1. **Names and quality are not in the archive at all.** They are roughly two
   thirds of a FASTQ. Earlier sessions measured them as 12.2% smaller than
   SPRING and 6.4% smaller than Genozip -- but that was measured by SUMMING
   CODED STREAM SIZES, the same method that hid 2.29 MB of missing reference
   data for a whole session. Those numbers are unverified.
2. **Verify them through capsule_roundtrip.sh, then integrate.** The container
   already matches streams by name, so adding `names` and `quality` is additive
   and needs no format break. The gate proves them the way it now proves
   sequence.
3. **Breadth.** Three locked accessions and T. cacao (Plantae, the only
   uncovered kingdom) have never completed.

A full-FASTQ tool that is verified end to end is worth far more than another
0.1% on the sequence stream, and the sequence stream has nothing left to give.
