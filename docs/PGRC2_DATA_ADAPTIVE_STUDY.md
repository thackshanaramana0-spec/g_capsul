# What PgRC2 actually does with the data, read from source and tested

## Their chain (pgrc-encoder.cpp:108-150)

    1. runQualityBasedDivision()          split reads by QUALITY SCORES
    2. runPgGeneratorBasedReadsDivision() split reads by OVERLAP STRUCTURE
    3. runHQPgGeneration()                assemble the pg from the kept reads
    4. runMappingLQReadsOnHQPg()          map everything else onto it

## Which stage actually matters, measured by disabling each

| flag | meaning | their archive on S. acidocaldarius |
|---|---|---|
| default | | 3,116,017 |
| `-q 1000` | quality division OFF | 3,116,012 (**5 bytes**) |
| `-g 0` | generator division OFF | **6,588,254 (2.11x worse)** |
| `-M 8` | min chars per mismatch 8 | 3,126,456 |
| `-p 100` | RC repeat length 100 | 3,122,302 |

**Quality scores contribute essentially nothing.** The whole benefit is the
overlap-structure division. This matters for scope: they are not exploiting
quality information that our sequence-only pipeline cannot see.

## The mechanism

`getBothSidesOverlappedReads` (AbstractOverlapPseudoGenomeGenerator.cpp:64-85):

    if (prevOverlap[i] && hasSuccessor(i)) continue;              // both sides -> keep
    if (hasSuccessor(i) && overlap[i] == readLength(i)) continue; // contained -> keep
    if (prevOverlap[i] == readLength(i)) continue;                // contained -> keep
    res[i-1] = false;                                             // else -> LQ, map it

A read enters the pseudogenome only if it is overlapped on BOTH sides. Chain
heads and tails do not; they are mapped instead. We emit whole chains,
endpoints included.

Its stop condition (GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp):

    overlapIterations = maxReadLength * overlappedReadsCountStopCoef   // 0.65

so the division sweep runs down to 0.35 x readLength. Our swept MINOV lands on
that same value on the long-read datasets -- H. salinarum 52 against 52.8,
S. aureus 51 against 51.8, S. acidocaldarius 87 against 87.8 -- and deviates
sharply on the 100-150 bp ones (16 against 35-52.5), where we win by the most.

## Implemented their admission rule. REFUTED, 3/3.

| dataset | ours | with both-side admission | delta |
|---|---|---|---|
| S. acidocaldarius | 3,140,692 | 3,143,537 | +2,845 |
| H. salinarum | 2,788,617 | 2,834,363 | +45,746 |
| E. coli | 8,200,969 | 8,279,566 | +78,597 |

On S. acidocaldarius it does shrink the pseudogenome by 17% (6,157,270 ->
5,111,602 bases), but the reads it evicts then have to be mapped, and their
mismatches cost more than the smaller span saves. Kept behind BOTHSIDE, off by
default, all outputs still lossless.

## A correction

While arguing for this I claimed their pseudogenome was 2.65x smaller than ours
(2,546,054 against 6,757,250). **That was wrong.** Their 2,546,054 is "Joined
mapped sequences", their literal AFTER MEM removal; our comparable figure is
2,500,501, which is SMALLER. I had compared our pre-removal pseudogenome against
their post-removal literal.

## The repetitiveness hypothesis: not supported

Sampled 32-mers over 120,000 reads per dataset:

| dataset | 32-mer repeat rate | top-1k share | margin |
|---|---|---|---|
| P. falciparum | 44.15% | 10.10% | +0.43% |
| S. acidocaldarius | 28.63% | 0.23% | **-0.83%** |
| S. aureus | 27.02% | 7.88% | +0.64% |
| H. salinarum | 21.58% | 2.21% | +8.58% |
| E. coli | 13.38% | 2.89% | +7.48% |
| L. major | 6.14% | 1.71% | +1.37% |
| P. aeruginosa | 5.06% | 0.35% | +0.89% |

corr(repeat rate, margin) = **-0.244**. P. falciparum is the most repetitive
dataset and we win it; H. salinarum is middling and we win it biggest. Repeat
content does not predict where we lose.

What does separate S. acidocaldarius is that our mapper places **51.3%** of its
reads, against 15-27% everywhere else, and each placed read averages 3.94
mismatches on 251 bp. That is what makes its mm_pos 23.4% of the archive -- and
PgRC2 spends 23.2% on exactly the same thing, so the shape defeats both tools.
