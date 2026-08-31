# Phase A — output-identical speed and RAM work

Every step here is a restructuring of *when and where* work happens. None may
change a single output byte. The gate after each step is:

  1. all 7 archives `cmp`-identical to /tmp/cmp3/perlayer_omp/*.arc
  2. lossless decode vs the original FASTQ sequence column

A step that fails either gate is reverted, not argued with. This is the same
standard that killed the uncapped parallel-probe change (+54-86% RAM for 6-11%
speed) and the architecture rewrite (0% size, 0% speed, +54 MB RAM).

## What is NOT in this plan, and why

- `OMP_WAIT_POLICY=passive`. Measured 9% faster on 19% less CPU, and it is a
  parameter: it tunes how threads wait instead of removing the wait. Kept as
  evidence that barrier overhead is real (see A1), not as a fix.
- Any per-dataset or per-kingdom switch.
- Assembly. Measured 1.64 s against PgRC2's 2.75 s on E. coli -- we are ahead
  there and it is not a target.
- Optimum Search Schemes / bidirectional FM-index. Provably beats pigeonhole
  for k mismatches, but needs a larger index; RAM is the column we lose, so
  this is the wrong direction.

## Measured baseline (OpenMP build, 7 files, single encode)

| layer | seconds | share |
|---|---|---|
| stream coding | 51.42 | 38.4% |
| pigeonhole mapping | 20.60 | 15.4% |
| pg MEM matching | 19.10 | 14.3% |
| round 1 (division) | 16.33 | 12.2% |
| load+filter+dedup | 16.02 | 12.0% |
| round 2 | 5.50 | 4.1% |
| emit chains | 4.09 | 3.1% |
| prefix seed index | 0.88 | 0.7% |
| **total** | **133.94** | unattributed 0.00-0.04 s/file |

Independent confirmation (perf, cycles by DSO, E. coli): our code 45.68%,
libgomp 21.29%, kernel 15.97%, liblzma 15.44%. PgRC2 on the same input:
own code 72.05%, libgomp 24.84%, kernel **1.72%**. Our kernel time is 9x
theirs -- repeated team creation and std::thread pools, not compression.

## A1 — hoist the sweep's parallel region

`sweep()` (line 582) runs a level loop, 134 levels in round 1 and 110 in round
2, and enters `#pragma omp parallel for` *inside* the loop. With three sweep()
calls that is ~600 team creations, each a fork/join with a barrier, for work
that now totals 21.8 s across 7 files.

Restructure to one team per sweep() call:

    #pragma omp parallel
    { for(L...){
        #pragma omp single      // compact tails, resize cand, assign ccnt
        #pragma omp for         // the per-tail search (unchanged body)
        #pragma omp single      // the serial commit, still in i order
    } }

Correctness argument, unchanged from the existing code: iteration i writes only
cand[i*CCAP..], ccnt[i], and seed[a]/ok[a] for its own unique tail; everything
else is read-only; the commit stays serial and in original order.

Two mechanical details:
- `break` cannot leave a structured block, so `tails.empty()` sets a shared
  flag inside the single and the loop breaks after the barrier.
- `num_threads(T)` with T=(w<4096)?1:NT cannot vary per level once the team is
  hoisted. Dropping it is semantically free -- the guard existed only to avoid
  fork cost, which is exactly what this step removes. Writes stay disjoint.

## A2 — compress-and-release in the coding pool

All 13 jobs currently hold their full compressed output in `results[]` plus
their coder work buffers, while every source stream is still alive. Measured
cost: +163/+189/+93/+86/+145 MB on ecoli/paeru/halo/sulfo/pfalc. On saureus and
lmajor it adds 0, because pigeonhole already pushed the peak past it (see A4).

Write each result into the archive and free it as soon as its turn comes, and
bound how many coders may hold work buffers at once.

## A3 — share the adaptive prefix

CORRECTION to an earlier claim: the shared prefix is NOT 61%. `MINOV` (argv[3])
first takes effect at line 751, before round 2; round 1 uses R1MINOV/SW, fixed
across all four candidates. So:

- shared by all 4 candidates: load + seed index + round 1 = 33.23 s (24.8%)
- shared by the 2 MAXMAP values at each MINOV: round 2 + chains = 9.59 s (7.2%)
- per candidate: mapping + MEM + coding = 91.12 s (68%)

4 x 133.94 = 535.8 s today; with two-level sharing 33.23 + 2(9.59) + 4(91.12)
= 416.9 s. That is 4x -> ~3.1x, a 22% saving. Real, but smaller than coding.

## A4 — invert pigeonhole mapping to streaming

The RAM cause on deep-coverage data. We index the pg and accumulate every
candidate hit into hit[]/hmm[] before choosing, so memory tracks coverage
depth: +877 MB on saureus, +715 MB on lmajor.

PgRC2 does the same pigeonhole filter with O(reads) memory: `addPattern()`
indexes the reads, `iterateOver()` streams the pg once, and each hit is
verified on arrival with only best-so-far kept in two flat arrays
(ReadsMatchers.cpp, executeMatching). No candidate list exists.

Adopt that structure. Same filter, same acceptance rule, same best-match
tie-break -- so the output must not move. This is the single largest RAM win
available and it also removes the allocation traffic.

## Order and expected effect

| step | targets | expected | risk |
|---|---|---|---|
| A1 | ~600 team creations, 21.3% libgomp + 16.0% kernel | speed | low, contained |
| A2 | coder buffers | -86..-189 MB, shallow files | low |
| A3 | the 4x adaptive multiplier | 4x -> ~3.1x | medium, main() refactor |
| A4 | candidate materialisation | 1357 -> ~650 MB, big files | high, largest change |

Coding's Amdahl floor (longest job 4.76 s vs PgRC2's 352 ms, 13 jobs vs their
82) is the largest single cost at 38.4%, but every fix for it either changes
the archive (block splitting) or costs RAM (uncapped probe parallelism), so it
belongs to Phase B and is gated on size, not attempted here.
