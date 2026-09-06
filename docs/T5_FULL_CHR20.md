# T5 (het-indel) at FULL chr20 — measured 2026-09-06

The published T5 (0.666 vs DiscoSNP++ 0.639) is a **400 kb window**. Every other
Claim 2 table except T3 shares that caveat. This is T5 repeated at **12.6M reads
/ whole chr20**, same box, same scoring order (normalise -> classify -> dedup,
`rtg vcfeval`, het-restricted, inside the GIAB confident BED), with DiscoSNP++'s
documented POS off-by-one correction applied as everywhere else in the project.

## Result — WIN, and the margin WIDENS at scale

| | G_CAPSUL | DiscoSNP++ | margin |
|---|---|---|---|
| **het-indel F1** | **0.5916** | 0.4190 | **+0.173** |
| precision | 0.8391 | 0.6887 | |
| recall | 0.4569 | 0.3011 | |
| TP | 3,556 | 2,343 | |
| FP | 682 | 1,059 | |
| FN | 4,226 | 5,438 | |

Truth: 8,882 het-indels on chr20 inside the confident regions.

**We win on precision AND recall** — more true positives with fewer false
positives, not a trade.

## Both tools degrade at scale; the competitor degrades more

| | window | full chr20 |
|---|---|---|
| G_CAPSUL | 0.666 | 0.5916 (-0.074) |
| DiscoSNP++ | 0.639 | 0.4190 (-0.220) |
| **margin** | **+0.027** | **+0.173** |

The published +0.027 was the narrowest margin in Claim 2 and the table most
exposed to a scale challenge. At full scale it is **6x wider**. The honest
framing is that the window UNDERSTATED this result.

**The absolute numbers must be restated at full scale, not the window values.**
0.666 does not hold at 12.6M reads.

## Cost — and this is the problem

| | time | peak RAM |
|---|---|---|
| G_CAPSUL full `CAPS_CALL` | 1,980 s wall (1,815 s caller) | **22.96 GB** |
| DiscoSNP++ | 75.7 s | 3.29 GB |

    ridx_build       147.2 s
    kc_H_build        22.2 s
    parallel_loop    491.4 s
    indel_pass     1,153.8 s   <- 64% of the caller
    TOTAL          1,815.4 s

**26x the time and 7x the RAM.** Indels require the full path: `CAPS_DBG_ONLY`
skips `indel_pass` entirely, so Option 0's 111 s / 6.21 GB archive path cannot
produce indels at all.

RSS must be SAMPLED for this path, not read at stage boundaries: 10.5 GB
entering `indel_pass`, 22.96 GB mid-stage, 17.2 GB at its end. The stage peaks
~12.5 GB above entry and releases ~5.7 GB before finishing, so the stage-end log
line understates the peak by 5.7 GB.

## Note on the SNV number from this run

This run also scored SNV at **0.8488** (P 0.9475, R 0.7687). That is NOT
comparable to Option 0's 0.8766: this is the full path's pileup-based SNV
channel on a differently-collapsed substrate, while Option 0 is the DBG_ONLY
bubble channel. Consistent with the A-vs-B result measured the same day
(bubbles 0.8766 > pileup 0.8072), the bubble path is the better SNV caller.

## What to fix next, in order

1. **`indel_pass` is 64% of the caller and had ZERO `#pragma omp`** in ~1,100
   lines -- single-threaded on 12 cores. Parallelisation is written (per-thread
   maps over independent kidx runs, deterministic merge, `CAPS_INDEL_SERIAL=1`
   restores serial) and awaits byte-comparison against this baseline.
2. **The ~9 GB transient inside `indel_pass`** is not explained by the
   structures identified so far (second seed index ~2 GB, `kidx` 1.45 GB).
   Prime suspect: `std::string Aalt = opp ? rc_str(cdb.contigs[ac]) : ...`
   constructs a FULL CONTIG COPY per pair inside an O(occ^2) loop. To be
   confirmed by sampling, not assumed.
3. Recall (0.457) is the weaker half. An indel needs a shared 25-mer anchor
   between two contigs; with 1,065,546 contigs at full scale those pairings get
   rarer. Testable against the FN set.
