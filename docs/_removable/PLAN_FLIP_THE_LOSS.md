# Plan: flip S. acidocaldarius by fixing chain quality

Everything below is measured. No step is a parameter change.

## Exact diagnosis (verified against PgRC2's own log, not inferred)

    PgRC2   HQ 241,774 reads -> pseudogenome 5,106,918   (21.1 bases/read)
    OURS    chained 233,124  -> pseudogenome 9,134,100   (39.2 bases/read)

Near-identical read counts; our pg is **1.86x longer**. Per-stream, that one
fact produces the whole 76,559 B loss:

| stream | ours | PgRC2 | delta |
|---|---|---|---|
| references | 100,709 | 11,130 | **+89,579** |
| positions | 1,459,981 | 1,398,530 | **+61,451** |
| literal/pg | 612,611 | 595,837 | +16,774 |
| RC info | 15,320 | 9,939 | +5,381 |
| mismatch positions | 723,746 | 721,906 | +1,840 |
| mismatch counts | ~158,233 | 163,068 | -4,835 |
| lengths | 11 | 5,315 | -5,304 |
| mismatch symbols | 119,903 | 208,234 | **-88,331** |
| | | | **+76,559** |

Positions cost log2(span) each; a longer pg also leaves more redundancy, needing
52,100 references against their ~6,700 (our matches average 131 bases, theirs
389). Our mismatch coding already beats theirs by 88,331 -- that advantage is
real and must not be traded away.

## Where our pg length goes

    [ovl] 219,688 links, mean overlap 89.6%
      90-100%  144,895 links (66%)  ~13 bases each  ->  1.88M
      <90%      74,793 links (34%)  38-163 each     ->  4.5M
      chain heads 13,436           251 each         ->  3.37M  (37% of the pg)

Two costs, both structural: links formed at poor overlap, and reads that never
link at all (each chain head pays a full 251 bases).

## Why this is NOT fixable by MINOV

MINOV was swept: 40 -> +29,708, 60 -> +15,115, 87 -> best, 113 -> +423,
150 -> +33,536. It is already optimal, and the curve is not unimodal, so no
search applies. Lowering it buys links at 87-base overlap that contribute 164
bases each -- worse than the head they save. The lever is link QUALITY at a
fixed floor, not the floor.

## The algorithmic change: minimum-degree-first matching

At each overlap length L the sweep solves a bipartite matching problem: open
tails on one side, open heads on the other, an edge where they overlap by
exactly L. We currently resolve it by **first-fit in array order** -- tail i
takes the first available head it finds.

First-fit is the worst greedy order for matching. A tail with many options can
consume the only head available to some other tail, leaving that tail unlinked;
it then either links at a lower L (paying 38-163 bases) or becomes a chain head
(paying 251).

The standard fix is to process vertices in increasing degree: a tail with one
candidate takes it before a tail with eight, which still has alternatives. This
is textbook greedy matching, provably no worse than arbitrary order in the
number matched, and typically far better.

**We already compute the degree.** The parallel candidate phase fills `ccnt[i]`,
the number of surviving candidates for tail i. The serial commit then walks
`i = 0..w` in array order. The change is to walk it in increasing `ccnt` order.

Cost: one counting sort over w entries per level, O(w) with 9 buckets since
`ccnt <= CCAP = 8`. No extra memory beyond a bucket table, no extra passes over
the reads, no new parameter.

## Why this preserves what we win

The change alters WHICH tail claims a contested head, not which links are legal.
The overlap floor, the mismatch ceiling, the coders and the mismatch streams are
untouched -- so the -88,331 mismatch-symbol advantage is unaffected by
construction. Output changes, so it is gated on size like any Phase B step.

## Expected effect, from the measured decomposition

If it converts a meaningful share of the 74,793 sub-90% links into 90-100%
links and reduces the 13,436 chain heads:

- pg length falls; positions save log2(old/new) x 516,904 / 8, at the **1.34x
  response measured** on the MINOV sweep (the coder responds slightly more than
  the log2 model, verified over five span points)
- references fall roughly in proportion to removed volume at 131 bases/match

A pg of 9.13M -> ~7M would give roughly 39,000 B on positions and ~29,000 B on
references. That is ~68,000 against a 76,559 loss: close, not certain.

## Verification, in order

1. Byte-identity is NOT expected (output changes by design). Instead:
   `[ovl]` histogram must show mass moving INTO the 90-100% bucket and chain
   count falling. If it does not, the mechanism is wrong -- stop.
2. pg length must fall on all 7 files.
3. Size gate: no file may regress; aggregate must improve.
4. Lossless on at least 3 files.

If step 1 shows movement but step 3 regresses, the mechanism works and the cost
is elsewhere -- report and revert, do not patch with a threshold.

## Rejected alternatives, with reasons

- Lower MINOV: measured optimal, curve not unimodal.
- Per-region MINMEM: refuted, monotonically worse.
- Splitting positions by region: refuted, we already code below that bound.
- Eliminating the second region: only 6% of span on this file, worth 7,752 B.
- Copying their 3-way split: we already beat them on the main region
  (2,368,666 vs 2,450,403 surviving).

---

# RESULT: min-degree-first matching REFUTED

Implemented and measured. Output byte-identical -- links 219,688, mean overlap
0.896, chains 13,436, pg 9,134,100, archive 3,191,341, all unchanged.

The reordering was verified to be happening (516,888 of 516,904 positions
reordered at L=250). It changed nothing because the premise was wrong:

    ccnt histogram at L=250, w=516,904 tails
      0 candidates : 500,145  (96.8%)
      1 candidate  :  13,330  ( 2.6%)
      2+           :   3,429  ( 0.66%)

**There is no contention to resolve.** 97% of open tails have no candidate at a
given overlap level, and under 1% have more than one. Matching order is
irrelevant when nearly every edge is uncontested.

## What this reveals instead

The assembly is **candidate-starved, not contended**. A tail links at the
highest L where a partner exists at all, and for 97% of tails no partner exists
at that L. So the lever is candidate GENERATION -- finding overlaps that exist
but are not detected -- not the order in which found candidates are consumed.

Three things bound generation, none yet measured:
1. the round-2 seed is SW=16 exact bases at one offset; a single sequencing
   error inside that window hides the overlap entirely
2. `rcmp` then requires the FULL L-base overlap to match exactly, so one error
   anywhere in a 225-base overlap kills the link
3. CCAP=8 is not binding here (97% have zero, not eight)

At ~1% error rate a 225-base exact overlap survives with probability ~0.99^450,
so long overlaps are detectable only between error-free read pairs. PgRC2's
`compareSuffixWithPrefix` is also exact, so this is not by itself the
difference -- but it bounds how much any matching-order change could ever have
achieved, which is zero, as measured.

## Status of the diagnosis

The diagnosis stands and is verified: our pg is 9,134,100 against their
5,106,918 for near-identical read counts, and that drives both the reference
(+89,579) and position (+61,451) losses. What is NOT yet identified is the
mechanism by which their chaining reaches 21.1 bases/read where ours reaches
39.2. Four hypotheses have now been tested and refuted; this one is the fourth.
