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

## What this does NOT fix

Pigeonhole mapping is genuinely per-candidate: all eight times differ, because
MAXMAP is the acceptance ceiling and the mapping loop is greedy — a read
accepted under a looser ceiling is removed from the pool and changes what every
later read sees. The placements are not separable into "search once, filter per
candidate". After this change pigeonhole becomes the dominant stage.
