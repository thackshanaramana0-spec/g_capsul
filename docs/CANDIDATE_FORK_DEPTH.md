# The candidate fork was at the right place and the wrong depth

## The finding

`encode_adaptive.sh` searches an 8-point grid: **4 MAXMAP x 2 MINOV**. The
encoder forks after round 1 so all eight candidates share load, the seed index
and round 1 (~24.8% of runtime).

But the two parameters are consumed at very different points:

| parameter | first read | source |
|---|---|---|
| `MINOV` | round-2 sweep | `106_inprocess.cpp:1682` |
| `MAXMAP` | pigeonhole mapping | `106_inprocess.cpp:1987` |

Round 2 and chain emission sit **between** them. So round 2 is a function of
`(input, MINOV)` alone, and the four candidates sharing a MINOV compute a
bit-identical assembly — three of which are discarded.

## The evidence — from the production log, not from reading the code

Across the eight candidates of a full HG002 run, `round2:` prints exactly two
distinct results:

```
probes=165361129 links=10254332   x4   (MINOV=16)
probes=173265221 links=10369827   x4   (MINOV=51)
```

and the per-candidate stage times were:

```
round 2:      794.06  786.40  785.36  781.53      <- one MINOV group
              124.72  124.59  115.42  114.00      <- the other
pigeonhole:   471.98 452.29 448.60 440.33 405.24 397.44 395.74 387.32   <- all 8 differ
emit chains:   25.19  16.35  14.83   7.79  7.75  7.57  7.54  7.47
```

Three redundant copies of a 785 s assembly and three of a 120 s one: about
**2,360 s of CPU recomputing an assembly that was already in hand**, and because
the candidates run concurrently at one thread each, the round-2 wall time is the
*single-threaded* cost of the slowest group.

## The change

Two-level fork.

* **Level 1** (after round 1, where the fork already was) forks one child per
  **distinct MINOV**. That child runs round 2 and chain emission once, with a
  real share of the cores instead of one thread.
* **Level 2** (immediately after `lap("emit chains")`, the last point before
  MAXMAP goes live) forks one grandchild per **MAXMAP** in the group.

Per-candidate state — `MAXMAP`, `ARCHIVE`, the `.cand<N>.d` scratch dir, the
caller's VCF/contig paths and the `chdir` — moves from level 1 to level 2. It
was verified that nothing writes to the cwd between the two fork points, which
is what makes deferring the `chdir` safe.

## Why this cannot change the output

A grandchild's state at the level-2 fork is exactly the state the corresponding
candidate had at the same point, because round 2 is a deterministic function of
(shared prefix, MINOV). It is also thread-count independent — already
established separately, since K=8/4/2 (1, 3 and 6 threads per candidate)
produced byte-identical archives.

The parent's selection rule is untouched: every candidate still writes its own
`.cand<N>` file and the smallest wins.

## Memory is held where it was

`K` is the number of candidate-sized resident images the box was measured to
afford. With two levels the concurrent images are (groups x members), so the
budget is split — `Kg = min(groups, K)` and `conc = K/Kg` — keeping the product
at `K`. At `K=1` (a memory-starved box) the groups run serially and round 2 is
still computed once per MINOV rather than once per candidate, so the win
degrades gracefully instead of inverting.

## Verification

| dataset | shipped | two-level | verdict |
|---|---|---|---|
| H. salinarum SRR29296997 | 15,650,029 | 15,650,029 | BYTE-IDENTICAL |
| S. acidocaldarius ERR12954017 | 15,159,145 | 15,159,145 | BYTE-IDENTICAL |
| E. coli SRR2584863 | 68,669,925 | 68,669,925 | BYTE-IDENTICAL |

Same winning candidate in every case. Sub-GB inputs gain little because round 2
is a small share of their runtime; the win scales with round-2 dominance.

## RETRACTED: "mapping is genuinely per-candidate"

This document originally ended with the claim below, written the same day:

> Pigeonhole mapping is genuinely per-candidate: all eight times differ, because
> MAXMAP is the acceptance ceiling and the mapping loop is greedy -- a read
> accepted under a looser ceiling is removed from the pool and changes what every
> later read sees. The placements are not separable into "search once, filter per
> candidate".

**That is wrong, and it was refuted by the log I already had.** The mapping loop
is not greedy over a shrinking pool: the pseudogenome is fixed for the whole
stage, so each read's placement is independent of every other read's. And the
search keeps the *minimum*-mismatch placement --
`lim = (cur==255) ? MAXMAP : cur-1`, updated only on strict improvement -- so
MAXMAP is the initial cap that decides whether a read is placed at all, never
which placement wins.

The proof was sitting in the production log. The `[MM] hist:` lines for the four
candidates of one MINOV group are bit-identical in every bin below the smaller
ceiling:

```
MAXMAP=7    0:207172 1:928067 2:239064 3:117518 4:76881 5:57665 6:45931 7:38044   8+:0
MAXMAP=11   ...identical bins 0-7...                       8:32688 9:28422 10:25114 11:22563
MAXMAP=18   ...identical bins 0-11...                                    >=12:109346
MAXMAP=29   ...identical bins 0-11...                                    >=12:195906
```

and the cumulative sums reproduce every candidate's mapped count exactly:

| ceiling | cumulative sum | reported `mapped=` |
|---|---|---|
| 7 | 1,710,342 | 1,710,342 |
| 11 | 1,710,342 + 108,787 = 1,819,129 | 1,819,129 |
| 18 | 1,819,129 + 109,346 = 1,928,475 | 1,928,475 |
| 29 | 1,819,129 + 195,906 = 2,015,035 | 2,015,035 |

A larger ceiling only **adds** reads in higher bins. So one search at the group's
ceiling answers every member by thresholding `readMM`, and the level-2 fork moves
again -- from after `emit chains` to after `lap("pigeonhole mapping")`.

Un-placing a read restores exactly the never-searched state, because the arrays
the search writes are initialised to precisely those values: `readMM` 255,
`matched` 0, `prc.assign(n,0)`, `ppos.assign(n,UINT64_MAX)` -- and the
second-region chain emission sets `ppos` for every admitted read afterwards.
`MAXMAP` itself is never read after the mapping stage.

| dataset | committed (round 2 shared) | + mapping shared | verdict |
|---|---|---|---|
| H. salinarum | 15,650,029 / 7.16 s | 15,650,029 / 6.86 s | BYTE-IDENTICAL |
| S. acidocaldarius | 15,159,145 / 10.52 s | 15,159,145 / 9.22 s | BYTE-IDENTICAL |
| E. coli | 68,669,925 / 21.59 s | 68,669,925 / 19.51 s | BYTE-IDENTICAL |
| **P. falciparum** | 67,382,606 / **94.81 s** | 67,382,606 / **64.54 s** | BYTE-IDENTICAL, **-32%** |

**The lesson, which is the same one this project keeps paying for.** I asserted a
structural property of the code ("greedy over a shrinking pool") from a plausible
reading rather than from the loop, and used it to close off the largest remaining
opportunity. The refuting evidence was already in a log on disk. Read the loop,
then check the histogram, before writing "not separable".

## Result at full scale

| | total | round 2 | pigeonhole |
|---|---|---|---|
| original | 914.80 s | 8 runs, max 478.30 s | 8 runs, 182-218 s |
| + huge pages | 885.14 s | 8 runs, max 478.30 s | 8 runs, 182-218 s |
| + round 2 shared | 812.79 s | **2 runs**, 62.4 / 66.9 s | 8 runs, 189-565 s |
| + mapping shared | **513.16 s** | **2 runs**, 62.3 / 66.2 s | **2 runs**, 112.1 / 112.9 s |

**1.78x faster, archive byte-identical (574,014,786 B), peak RSS +0.6%.**

Correctness is exact, not approximate: all eight `mapped=` / `appended=` counts
reproduce the baseline read for read.

```
leftovers=2112037  mapped=1660493 1767491 1873738 1954531
leftovers=2187048  mapped=1710342 1819129 1928475 2015035
```

## What is left

The shared prefix is now 235 s of the 513 s (load 14, round 1 35, round 2 66,
emit chains 7, mapping 113). The remaining ~278 s is the genuinely per-candidate
tail: the second-region sweep, MEM matching (102-114 s per candidate) and stream
coding, run for 8 candidates at concurrency 4. MEM matching is the largest piece
and is genuinely per-candidate -- it runs over a pseudogenome whose second region
depends on which reads this candidate appended, and those counts differ
(157,506 to 476,706). That is the next place to look, and it should be looked at
by reading the loop and checking a histogram, not by asserting a property.
