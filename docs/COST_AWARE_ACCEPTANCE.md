# Cost-aware match acceptance (COST_GATE)

## What it replaces

MEM matches were accepted purely by length: `len >= MINMEM`. MINMEM is a
length PROXY for "is this reference worth its own storage cost" -- it has no
idea that a reference costs a source varint + a length varint + a gap varint
+ 1 RC bit regardless of how long the match is, while literal costs a
near-fixed ~2 bits/base. A match that just clears MINMEM can still lose.

## The gate

Computed from the SAME byte-cost formulas the real coders use, not a guessed
constant:

    src_bits  = log2(dst)                    -- upper bound on refc's own
                                                 bound (true bound is <= dst
                                                 in every case, self or cross,
                                                 so this can only OVER-count,
                                                 never wrongly accept a bad
                                                 match)
    len_bits  = 8 * vint_bytes(len - MINMEM)  -- exact, same vint() as the
                                                 real mem_len stream
    gap_bits  = 8 * vint_bytes(dst - lastend) -- exact, same vint() as the
                                                 real mem_dstgap stream
    rc_bits   = 1

    lit_bits_saved = len * 2.0   -- the information-theoretic MAXIMUM for 4
                                    symbols. A match this rejects loses even
                                    under the most generous possible
                                    assumption about literal's actual cost.

    accept = (src_bits + len_bits + gap_bits + rc_bits) < lit_bits_saved

Off by default (COST_GATE env var). Byte-identical to the base build when off
-- verified directly.

## Measured, isolated from the (already-shipped, off-by-default) mismatch
## tolerance work, so the two effects are not conflated

**Drosophila SRR40104920, 2.6x coverage** (MEM_MAXMM forced to 0 on both
sides, so this is purely the acceptance-criterion effect):

    without gate: 36,982,418 B, 640,685 matches
    with gate:    36,971,536 B, 633,958 matches   (-10,882 B, -0.029%)

6,727 matches correctly identified and rejected as unprofitable. Lossless,
verified through the real archive.

**All 7 locked datasets**, lossless on every one:

| dataset | delta |
|---|---|
| halo | +0.0075% |
| ecoli | +0.0027% |
| sulfo | +0.0000% |
| paeru | +0.0014% |
| saureus | +0.0015% |
| pfalc | -0.0125% |
| lmajor | -0.0084% |
| **aggregate** | **-0.0046%** |

## Honest assessment

At normal coverage (25-59x) MINMEM=24 was already a verified interior
optimum -- almost every match that clears it is already profitable, so the
gate has almost nothing to reject. The result is noise-level (all 7 under
0.01%, both directions), consistent with greedy-parse interaction effects
(rejecting one match shifts where the NEXT seed search starts) rather than a
real gain or loss.

At low coverage, where MANY more marginal matches exist, the gate finds real
work to do -- but the absolute size (10,882 B) is nowhere near enough to
close SPRING's 6.5% (~2.4 MB) advantage on the same file. This is a genuine,
principled, ALGORITHMIC improvement -- the first mechanism change this
session with a positive rather than negative effect at low coverage -- but it
is a partial step, not a solution. Combined with the earlier finding (three
independent architectures -- overlap-graph, hash-reorder, de Bruijn graph --
all hit the same coverage floor, and our own literal/pos_abs streams are
independently confirmed at their entropy bounds), the honest conclusion
stands: there is no larger lever left to find in this scope.
