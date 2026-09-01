# Coverage-aware placement tolerance (MAXMAP)

## What was wrong

MAXMAP is the mismatch ceiling for accepting a read PLACEMENT via pigeonhole
mapping -- our cheap, read-relative matching path, where a mismatch position is
an offset WITHIN a read (small number, coded via mm_pos/mm_sym) rather than an
offset into a 159 MB pseudogenome.

It was set to `Lmax/13` (~11-19), chosen by a minimax-regret sweep over the 7
locked datasets. Every one of those measures leftover_frac 0.222-0.518 -- all
normal coverage. The constant was never swept against low-coverage data, where
the ALTERNATIVE to accepting a mismatch-heavy placement is not "slightly more
literal": it is falling through to the region-scale MEM path, whose mismatches
cost ~21 bits each against ~2 bits on the read-relative path.

SPRING's equivalent (`THRESH_ENCODER`) is 24. Ours was 11.

## Measured

Sweeping MAXMAP on Drosophila SRR40104920 (2.6x coverage, leftover_frac 0.812),
everything else fixed:

| MAXMAP | archive | reads placed |
|---|---|---|
| 11 (old default) | 36,982,418 | 158,351 |
| 16 | 36,910,105 | 182,024 |
| 20 | 36,878,860 | 199,432 |
| 24 (SPRING's value) | 36,856,348 | 217,697 |
| **30** | **36,836,970** | **250,841** |
| 39 (our break-even ESTIMATE) | 36,873,742 | 312,464 |
| 50 | 37,050,398 (worse than baseline) | 398,515 |

A clean interior optimum at 30 -- worse on both sides, so this is a real
optimum, not "more is always better". Note the measured optimum (30) sits well
below our own documented break-even estimate of ~0.26*L = 39, so the estimate
overshoots and the sweep was necessary.

## The gate

Same measured signal and same T0/T1 margin as the (rejected) MEM_MAXMM work:

    lf = leftovers / n
    T0 = 0.60   -- above the locked datasets' measured max of 0.518, with margin
    T1 = 0.85
    div = MAXMAP_DIV - (MAXMAP_DIV - 5) * clamp((lf-T0)/(T1-T0))
    MAXMAP = max(Lmax / div, existing Lmax/MAXMAP_DIV)

Never lowers MAXMAP, only raises it, and only above T0. Verified: halo
(lf=0.385) and E. coli (lf=0.298) produce MAXMAP=11 and byte-identical archives
against the pre-change build; Drosophila (lf=0.812) ramps to 25.

Honest limit: the ramp's SHAPE (linear in leftover_frac) is principled, but its
endpoint (div=5) is calibrated on ONE low-coverage dataset and is not yet
cross-validated on a second.

## Also in this commit: extension mismatch tolerance turned OFF by default

Built, debugged (it exposed a real decoder bug -- batch decode-then-apply reads
a stale byte when a later match's source overlaps an earlier match's
destination, fixed with a streaming decoder) and verified lossless, then
measured honestly: 36,982,418 -> 37,473,181 B on the dataset it was designed
for, a 490,763 B LOSS. It removed 11 M literal bases worth 212,041 B while its
references cost 369,424 B. It had been left ON by default; that was costing real
bytes on exactly the data it was meant to help. Code stays (it is correct, and
the decoder bug it found was real); enable with MEM_MAXMM_OVERRIDE only to
reproduce the measurement.

## Result

Drosophila 2.6x: 36,982,418 -> 36,851,932 B (-130,486 B, -0.353%), lossless.
All 7 locked datasets byte-identical and lossless.
