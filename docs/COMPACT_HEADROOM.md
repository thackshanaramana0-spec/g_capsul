# CAPSULE COMPACT — where the seconds actually are

Method: `perf record -F 299` over a full HG002 encode with a `-g` build
(633,156 samples), plus the encoder's own per-stage timers, plus purpose-built
funnel instrumentation. Machine: AMD EPYC 9555, 12 vCPU, L1d 768 KiB,
L2 6 MiB, **L3 192 MiB**. HG002 chr20, 12.6 M reads, 3.99 GB.

## 1. Wall-clock decomposition (baseline 240.73 s)

| stage | s | share | threads |
|---|---|---|---|
| pigeonhole mapping | 64.4 | 27% | 6 |
| pg MEM self-match | 47.4 | 20% | 6 |
| round 1 (division) | 43.2 | 18% | 12 |
| round 2 (assembly) | 27.2 | 11% | 6 |
| load+filter+dedup | 22.5 | 9% | 12 |
| coding pool (bounded by `pos_abs` 19.5 s) | 19.5 | 8% | 6 |
| emit chains | 6.4 | 3% | 6 |

**Only 6 of 12 cores run the post-fork 70% of the pipeline.** The adaptive
sweep forks `conc=2 x threads=6/6`. On this input the two candidates differ
only in MINOV, and the winner beats the loser by **99,117 B — 0.0173%**. That
is the price of halving the machine for 165 s of the 237 s run. It is a real
size/speed trade and is the user's call, not a defect.

## 2. Fixed: the MEM matcher had no candidate discriminator (−12.3%)

The hottest single line in the encoder was `extendTol`'s compare loop
(:3591, 9.09% of all samples) with its `lookup` at 3.16%. Every candidate --
up to `MAXCAND=64` per query position -- went straight into `extendTol`, a
**random access into a 152 MB ASCII text**. The false candidates were paying
for the cache miss, not the compare.

Round 2 had already solved this exact pattern with `pext` and never carried it
across. `sext[]` does, with a proof of output-preservation (see the commit).

    baseline  240.73 s   573,767,964 B
    sext      211.04 s   573,767,964 B   BYTE-IDENTICAL, RAM -27 MB
    MEM stage 47.44 s -> 21.80 s (-54%)

Verified byte-identical on 6 datasets spanning both tolerance regimes
(`MEM_MAXMM=0` and `=2`), 148/301/variable read lengths. The gain scales with
pseudogenome size, which is the mechanism's own prediction: L. major (1.66 GB)
MEM 10.78 -> 7.09 s (-34%); E. coli and SARS-CoV-2 unchanged, because their pg
fits in L3 and there is no miss to remove.

## 3. Measured, not yet fixed: mapping is repeat-bound

`CAPS_MAPDBG=1` on a `-DMAPFUNNEL` build bins each probe's bucket walk:

    probes=43,289,862  cands=3,663,240,336  (84.6 per probe)

| bucket size | probes | **work** |
|---|---|---|
| <= 256 | 89.0% | 45.0% |
| 257-512 | 2.08% | 17.85% |
| **513-1024 (at the cap)** | **3.05%** | **36.89%** |

**3% of probes cause 37% of the work; 5% cause 55%.** The walk is dominated by
a handful of repeat-derived keys sitting at the `[seedcap]` ceiling of 1024.

**Why the `sext` fix does not port here.** A window discriminator can only
reject when the window's mismatches exceed the tolerance. Mapping accepts up to
`MAXMAP=29` mismatches, and an 8-base window holds at most 8 -- it can never
prove rejection. To reject at 29 you need a >29-base window, where a random
candidate averages ~22 mismatches and so is almost never rejected. The filter
is not merely weak here, it is structurally inapplicable.

**What the literature says to do instead.** This is the classic
high-multiplicity-seed problem. The current answer -- cap the bucket at 1024 --
drops 126,723 entries (0.828%) across 76 keys **in index order, i.e.
arbitrarily**. LAST's adaptive seeds (Kiełbasa et al., *Genome Research* 2011)
lengthen the seed at over-frequent keys until multiplicity falls below a
threshold, which discards an entry **only when it actually disagrees**.
Adaptive extension therefore dominates capping on both axes: strictly less
arbitrary sensitivity loss, and it removes the 37%-of-work tail rather than
truncating it.

**Why it is not done here.** It is not output-preserving: it changes which
candidates are examined, so its gate is archive size across the locked set, not
byte-identity. That is a larger experiment than this pass, and it must be
measured against size, not asserted. Recorded as the next target with the
measurement that justifies it.

## 4. What is NOT the bottleneck (so nobody re-derives it)

- `readMM` (12.6 MB) is **not** a cache problem: L3 is 192 MiB and holds it.
  The 3.9 billion misses come from `rpk`, the 480 MB packed read store.
- Quality coding is ~9.6% of samples and is vendored fqzcomp, already measured
  against its own bound elsewhere in this project.
- The mapping verification is already 2-bit packed, popcount-based and
  software-prefetched. Its cost is candidate COUNT, not per-candidate work.
