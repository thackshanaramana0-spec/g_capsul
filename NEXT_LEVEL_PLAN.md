# Where the next real gain is

Written 2026-08-29, after the three-axis lock (size 0.92x, speed parity, RSS
0.89x). This is the plan for what comes after, and it opens by withdrawing a
conclusion this file's predecessors stated twice.

## Withdrawn: "their assembler is better than ours"

The claim was that PgRC2 tiles the genome better -- mean overlap 133.7 against our
119.9, their pg 21.10M against our 21.62M -- and that closing it needed a better
layout algorithm. **That was a read-count error.** Their stage log, on this exact
input:

    HQ reads count: 999340
    Found 148065 duplicates          <- full set: 851,275 unique, same as ours
    HQ reads count: 708819
    Found 148065 duplicates          <- every duplicate landed in HQ
    LQ reads count: 290521
    Target pseudogenome length: 21104500

708,819 - 148,065 = **560,754 unique reads** build their main pg, and
560,754 + 290,521 = 851,275 -- our unique count exactly, which is what confirms
the reading.

| | reads in main pg | main pg | B/read |
|---|---|---|---|
| ours | 586,345 | 21,619,175 | **36.87** |
| PgRC2 | 560,754 | 21,104,500 | 37.64 |

**Our layout is 2% tighter per read than theirs.** Their pg is smaller because it
holds 25,591 fewer reads, not because it packs them better. The greedy is not the
problem and never was.

What differs is the **division policy**: they push 290,521 reads onto the mapping
path, we push 264,930. A mapped read costs mismatch symbols plus a position; an
assembled read costs ~36.9 bases of literal. Which is cheaper decides the optimum,
and right now we cannot answer that, for the reason in item 1.

## Also refuted, cheaply: cycle handling

PgRC2 does **not** avoid cycles during the sweep. `findOverlappingReads` calls
`overlapSortedReadsAndMergeSortSuffixes<false>` -- `avoidCyclesMode = false` --
so it links freely and calls `removeCyclesAndPrepareComponents()` afterwards to
break each cycle. Ours refuses a link that would close a cycle (`ch_h[a]==b`) and
falls through to a worse candidate. That looked like a real algorithmic edge for
them.

Instrumented, over the whole run:

    round1 refusals=9  reads-blocked=8  mean L at first block=69.4  downgraded=1  left-unlinked=8
    round2 refusals=9  reads-blocked=8  mean L at first block=87.0  downgraded=0  left-unlinked=9

Nine refusals. Worth ~500 bases. Not a lever. (Their cycle-breaking loop also
shadows `minOverlap` and `minOverlapIdx` with fresh declarations inside the `if`,
so the intended "break at the weakest link" never updates and it always breaks at
the entry point -- a real bug in PgRC2, and immaterial for the same reason.)

## 1. The unmeasured stream that blocks everything else

**This is the same failure that cost us twice already** -- relaxed acceptance
hiding cost in the mismatch stream (stage 27), aggressive MEM removal hiding cost
in the reference stream (stage 32). An optimisation looks free because its price
lands somewhere nobody counts.

We map 252,865 reads into the pseudogenome and **store nothing for their
positions.** PgRC2 pays 683,370 B of reads-list offsets and 176,241 B of
mismatch counts for exactly this. Our four-stream total of 5,815,545 against
their 6,352,312 compares the four streams both sides produce, which is honest as
far as it goes -- but it means the cost of moving a read from the assembled path
to the mapped path is currently **zero in our accounting and positive in theirs.**

Until that stream exists, every division-threshold experiment will report that
mapping more reads is free, which is false.

Rough economics, to show why the sign is genuinely unknown:

  - a read assembled into the pg: ~36.9 bases of literal, ~47% of which the MEM
    stage removes, at 1.906 bits/base -> **~4.6 B**
  - a read mapped instead: 4.94 mismatches at ~1.55 bits -> ~0.96 B, **plus** a
    pg offset, which at their rate is ~2.7 B -> **~3.7 B**

Those are close enough that the answer depends on the actual coded position
stream, not on an estimate. Build it, code it, round-trip it, then decide.

**Deliverable:** emit (pg offset, orientation) per mapped read, code it, and add
it to both sides of the table. Expect our total to rise by ~500-600 KB and
theirs by 859,611 B; the lead should survive and may widen, since we map fewer
reads. Report the corrected table whatever it says.

## 2. The division threshold, swept with the full cost model

Only after item 1. `R1MINOV` (round-1 classifier floor, currently 16, split from
`SEEDW` at stage 17) is the admission knob: raise it and fewer reads qualify as
both-sides-overlapped, so the main pg shrinks and the mapping path grows. It has
never been swept against a total that prices both sides, because one side had no
price.

Sweep R1MINOV over roughly 16, 24, 32, 40, 48, and at each point report
literal + references + order + mismatches + **positions**, plus peak RSS. The
objective is the sum; any single component moving is not evidence.

The prior from the numbers above is that we are admitting too many reads and the
optimum is stricter than 16 -- but that prior comes from PgRC2's choice, and
their cost model is not ours (they pay more per assembled read and we pay less,
which pushes our optimum the other way). Measure it.

## 3. Maximum-weight cycle cover instead of greedy

This is the genuinely different algorithm, and it is third because items 1 and 2
may absorb the available gain first.

Greedy overlap chaining -- what we do, what PgRC2 does, what SPRING and every
other reordering compressor does -- is the classic approximation to shortest
common superstring. Its proven guarantee is at most (13 + sqrt 57)/6 ~ 3.425
(Kaplan et al.; recent refinements), conjectured 2. The better algorithms in that
literature do not improve the greedy; they replace it:

  - build the overlap graph
  - compute a **maximum-weight cycle cover**, which is exactly a maximum-weight
    bipartite perfect matching (left = each read as a predecessor slot, right =
    each read as a successor slot, weight = overlap length) -- the assignment
    problem, solvable *exactly* in polynomial time, unlike the Hamiltonian path
    we actually want
  - break each cycle at its minimum-weight edge to get a path cover
  - concatenate

That route underpins the 2.475-approximations (Mucha; Paluch et al.) against
greedy's 3.425. Greedy is a 1/2-approximation to the matching; the cycle cover is
optimal. The gap between them is the headroom, and it is not a constant factor we
can bound in advance for real read data -- it has to be measured.

Feasibility is the real question, and it is better than it looks:

  - the overlap graph is **sparse** -- each read has a handful of candidates, and
    our flat prefix index already enumerates them (stage 43)
  - exact Hungarian is O(n^3) and hopeless at n = 851,275, but
    **(1-epsilon)-approximate maximum weight matching runs in near-linear time**
    (Duan-Pettie 2014, O(m log(1/eps)/eps)), which at m ~ a few million edges is
    seconds, not hours
  - it parallelises, and it is a *global* optimum where greedy is a sequence of
    local first-fits, so it should be robust where greedy is order-dependent

**Nobody in the FASTQ compression literature does this.** PgRC, PgRC2, SPRING,
Minicom, zDUR (BMC Bioinformatics 2026, clustering plus similarity trees) all use
greedy or clustering. A matching-based layout would be a genuine methodological
contribution rather than a tuning of theirs -- which is the thing this project has
been missing, since the current 8.45% lead rests almost entirely on the order
coder (Lehmer + Fenwick at log2(n!)), a component anyone could port into PgRC2 in
an afternoon.

**Staging, so it can be killed early and cheaply:**

  1. Dump the candidate overlap graph the current sweep already computes (edges
     with overlap >= MINOV, weight = overlap). Report edge count and degree
     distribution. If the graph is nearly a matching already, greedy is near
     optimal and this stops here.
  2. Compute the greedy matching weight and an upper bound on the optimal
     matching (LP relaxation, or a cheap dual). **Decision gate:** if greedy is
     within ~1% of the bound, the ceiling is too low to justify the work, and
     that is a publishable negative result on its own.
  3. Only if the gap is real, implement approximate max-weight matching, break
     cycles at minimum weight, and re-run the whole pipeline. Every downstream
     verified result (permutation round trip, MEM counts, coder round trips) must
     be re-checked, not assumed.

## What is deliberately not on this list

**The order stream.** It is 40% of the archive at 2,309,967 B, and it is already
at log2(999340!) ~ 2,309,375 -- the information-theoretic floor for a uniform
permutation. There is no headroom unless the input order carries structure
(paired-end mates, coordinate-sorted input), and SPRING already exploits pairing,
so that would be catch-up rather than a lead. Worth one measurement on a locked
accession before believing either way; not worth a project.

**Generalisability, which is a gap not a lever.** Every number in this progression
and in `BEST_METHOD_C_REIMPLEMENTATION.md` comes from one dataset, `yeast_sub.fq`.
`ecoli_sub.fq` and `pf_sub.fq` exist and have never been run end to end against
PgRC2. That is a prerequisite for any claim in a paper and it is cheap. Do it
before item 3, not after.

## Sources

  - Kaplan & Shafrir, and the recent line on greedy SCS bounds:
    https://arxiv.org/pdf/2111.03968 , https://arxiv.org/pdf/2407.20422
  - Cycle-cover route to 2.475: Mucha (SODA 2013), Paluch et al.
  - SCS with reverse complements: https://arxiv.org/pdf/2603.26176
  - Duan & Pettie, linear-time approximate maximum weight matching
  - PgRC2: https://academic.oup.com/bioinformatics/article/41/3/btaf101/8051895
  - zDUR: https://link.springer.com/article/10.1186/s12859-025-06364-1
