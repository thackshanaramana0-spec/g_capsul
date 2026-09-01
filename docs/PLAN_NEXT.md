# What is open, and the order to do it in

Written from measurement, not recall. Nothing below has been started.

## Where we actually stand

| | ours | PgRC2 | |
|---|---|---|---|
| size, 7 files | 79,693,000 | 83,160,210 | **+4.17%** |
| wall, 7 files | 113.6 s | 42.01 s | 2.70x slower |
| peak RSS, worst file | ~930 MB (L. major) | 368 MB | 2.5x heavier |

Session movement: wall 243.78 s -> 113.6 s (-53%) with sizes improving, not
traded. 7/7 lossless against the original FASTQ. Repo clean, nothing unpushed.

## Facts that set the order

1. **Every result today was validated on inputs <= 1.7 GB.** C. elegans (11 GB)
   and T. cacao (22 GB) are on disk and have never been run with ANY of this
   work -- not OpenMP, not A1/A3/A4, not B1/B2/B3.
2. **A4 introduced an O(reads x threads) allocation.** The per-thread `slot`
   array is 4 B x reads x 12. Measured 0.21 GB on L. major; projected 1.37 GB
   on C. elegans and 2.74 GB on T. cacao. 82 GB RAM available, so it should
   fit -- but it is a scale behaviour nobody has observed.
3. **The probe still disagrees after B3.** Re-measured: u32 streams pick
   method 6 / 6 / 6 / 1 across E. coli, H. salinarum, S. acidocaldarius; byte
   streams pick 0, 2, 3, 4. Deleting the probe today would cost size. So the
   route to removing it is more decomposition, not a shorter candidate list.
4. **Three streams are still undecomposed**, with measured floors:
   - `read_lengths` 3,106,518 B, one distinct value (150), entropy **0**
   - `mm_count_per_read` 2,471,292 B, 85.79% zeros, floor 145,332 B
   - `pos_strand` 154,456 B, 87.14% zeros, floor 19,231 B

## The order

### 1. Scale validation FIRST  (correctness, blocks everything)
Run yeast (2.1 GB) -> C. elegans (11 GB) -> T. cacao (22 GB), one at a time,
recording wall, peak RSS and archive size, and verifying lossless.

Rationale: a failure at 22 GB invalidates more than any optimisation below adds.
Watch specifically for the A4 slot arrays and for disk (53 GB free against a
22 GB input plus dumps). This is also the only way to refresh the 8-file
aggregate, which is stale since the coder work.

### 2. `read_lengths` -- constant-stream elimination
Entropy 0, one distinct value. The exact representation is `(count, value)`,
about 3 bytes. Removes a 3.1 MB allocate-and-code per file. Same proven method
as B3; near-zero risk. Note it must stay correct for genuinely variable-length
input (SARS-CoV-2), so the constant case is detected, not assumed.

### 3. `pos_strand` and `mm_count` -- finish the decomposition
`pos_strand` is 87% zeros and unsplit. `mm_count` is split already, but codes
BOTH the split and the flat form and keeps the smaller -- three encodes where
one would do, now that the split is known to win.

### 4. Re-test the probe
After 2 and 3, re-run the pick log. If each remaining stream is homogeneous
enough that one coder dominates, the 7-way probe becomes deletable -- that is
the ~16.7%-of-profile item and most of the remaining page faults. If it still
disagrees, the answer is more decomposition, and this records which stream is
still mixed.

### 5. `pos_abs` -- scheduling, not modelling
It is the new pool floor at 2.65 s but sits at 0.94x its set-plus-permutation
bound, so there is no size in it. Only touch after 4, because if the probe goes
away the floor moves again.

### 6. Close the verification gap
There is **no archive decoder**. Losslessness is proven by decoding the raw
stream dumps and diffing against the original FASTQ -- which is a real round
trip -- but nothing reads the archive back. Every method header is
self-describing; no inverse is written. This should be closed before any
external claim about the format.

### 7. Refresh the external comparisons
SPRING and Genozip were only ever run on the virus and fungi files. The 7 core
datasets have no SPRING/Genozip numbers. Also unexplained: Genozip produced
162 MB on yeast where SPRING produced 24 MB -- a 6.6x gap between two mature
tools is an anomaly to diagnose before it appears in any table.

### 8. S. acidocaldarius, the last loss (-2.50%)
Cause is known and structural: their 3-way pg split (hqPg / lqPg / nPg, each
self-matched separately) needs 16-72K of references where our single pg needs
99-473K. This is a project, not a tweak, and it is last because it is the only
item whose payoff is smaller than its risk.

## Not in this plan, deliberately
- Shrinking the probe's candidate list. Methods 2, 3, 4 and 5 have never won on
  a u32 stream in ten observations, so trimming would likely be free -- but it
  treats the symptom. Decomposition is the cause-level fix.
- Any per-dataset or per-kingdom switch.
- Assembly. Measured 1.64 s against their 2.75 s; we are ahead there.
