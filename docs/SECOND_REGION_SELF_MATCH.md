# The second region is never self-matched — a real, quantified, unclaimed win

## The gap

`run(Q, qlen, CROSS, ...)` matches the second region against the MAIN pg only.
There is no self-match, so the second region can never reference itself. Where
the main pg is large this is invisible; where it is small, nothing is removed.

Measured on SARS-CoV-2:

    MEM second: 14,714,404 -> 14,392,090   (2.2% removed, against a 1,110-byte main pg)
    32-mers in that region: 1,839,106 sampled, 473,440 distinct
    repeat occurrences: 1,365,666 = 74.3%

**74.3% of the region is repeat content and 2.2% of it is being removed.** The
redundancy is there; it is unreachable by construction.

## Why SARS-CoV-2 is shaped this way

Amplicon sequencing. Reads start at ~98 fixed primer positions, so after dedup
(92.3%) the surviving unique reads are error-variants at the SAME coordinates --
aligned, not staggered. Suffix-prefix chaining needs staggering, so round 1 finds
73 both-side-overlapped reads out of 72,244 (0.1%), the pseudogenome reaches
2,639 bases against a 29,900-base genome, and 85% of reads (61,120) are appended
raw. The algorithm is correct; the data shape defeats it.

Floor estimate: a 30 kb genome plus 72,244 positions at log2(30,000) = 14.9 bits
is of order 250-350 KB. We pay 841 KB. Roughly 2.5-3x headroom.

## What was built, and what it measured

Two implementations, both correct in structure:

1. **Whole-pg self-match** -- treat main+second as one SELF_FWD text.
   E. coli +95,373. Merging imposes SELF_RC's halved cap (capL=(qlen-qp-s)/2) on
   cross-matches that were previously uncapped.
2. **Third pass** -- keep both cross-passes, add SELF_FWD/SELF_RC over the second
   region with SRCBASE=main_pg_end. Required making the seed index rebuildable
   and the source text a parameter (it was hardcoded to pg.data()/main_pg_end).

Results, and a clean discriminator:

| | unremoved after cross | self-pass |
|---|---|---|
| SARS-CoV-2 | 97.8% | **-94,279** |
| S. acidocaldarius | 39.2% | -171 |
| E. coli | 12.7% | **+61,792** |

Where the cross-pass already removed nearly everything, the extra references cost
more than the short matches they buy. Where it removed nothing, the win is large.

## Why it is NOT shipped

**The archive it produces is LOSSY.** Verified: SARS-CoV-2 decodes to the right
number of reads with the wrong content. SELF_FWD guarantees s<qp within the
region, so src<dst holds globally once both are offset by main_pg_end -- but
something in the reconstruction does not honour references whose source lies in
the second region. A flag that emits corrupt archives is worse than no flag, so
the whole change was reverted rather than gated.

## What would be needed

1. Find why second-region-sourced references do not reconstruct. The decoder
   applies references in dst order; the RC path with a non-zero SRCBASE is the
   first place to look.
2. Only then decide when to run it. The discriminator above is clean and is
   measurable at the right time -- after the cross-pass, before this one -- but
   turning it into a cutoff is a threshold, and a cost-based rule needs a model
   of self-redundancy that has not been validated. Three cost models were
   attempted earlier this session and all three failed.

The opportunity is real and quantified. The implementation is not correct yet.
