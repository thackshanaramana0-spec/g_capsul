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

---

# Why the remaining low-coverage gap does not close (both candidates measured)

After the MAXMAP ramp, Drosophila 2.6x sits at 36,851,932 B against SPRING's
34,723,840 -- a 2,169,648 B gap. Stream-level accounting locates it exactly:

    SEQUENCE content   ours 29,193,716   SPRING 29,891,230    -697,514  (we WIN)
    METADATA/refs      ours  7,623,810   SPRING  4,756,648  +2,867,162  (we lose)

We already beat SPRING on the actual DNA. The whole gap is metadata, and within
it, one block SPRING has no equivalent for: region-scale reference streams
(mem_triples + mem_dstgap + mem_len + mem_rc) = 3,091,068 B, which is 1.42x the
gap. Two ways to remove it were identified and both were MEASURED, not assumed.

## Candidate 1: convert region-scale matches to read-scale placements. DEAD.

Read-scale placement (pigeonhole + mm_pos/mm_sym) costs ~2 bits per mismatch
against ~21 for a region reference, so converting looked worth ~1.23 MB.

It cannot be done, for a structural reason: **the average MEM match in the
second region is 44.9 bases** (31,974,409 bases over 711,422 matches) against a
150 bp read. These are mid-read PARTIAL OVERLAPS, not duplicate reads. Read-scale
placement requires a read to match another read over its full length; chaining
requires suffix-prefix adjacency. Neither can express "my middle 45 bases match
your middle 45 bases" -- only region coordinates can. The matches are the wrong
SHAPE to convert, not merely expensive.

## Candidate 2: drop MEM on the second region, let a BWT compressor find the
## fragment redundancy implicitly (what SPRING does with read_unaligned + bsc).

Measured directly on the extracted raw second region (143,457,338 bases):

    bsc on raw, no MEM at all          31,641,218 B   (1.7645 bits/base)
    ours (literal share + refs)        30,066,062 B

**Our explicit MEM approach beats it by 1,575,156 B.** The 3.09 MB of reference
metadata is load-bearing: removing it costs more than it saves.

## Conclusion

The region-scale metadata cannot be re-expressed more cheaply (candidate 1) and
cannot be removed (candidate 2). Both were tested rather than argued. At 2.6x
coverage the gap to SPRING is not closable by any mechanism identified here.

This is a bounded, specific claim: at NORMAL coverage -- every dataset in the
locked suite -- we beat SPRING by 21-37% and Genozip by 27%, and Genozip loses
to us by 27% even at 2.6x. The unclosed case is one regime, on a dataset outside
the locked suite, where PgRC2 crashes outright and SPRING was never benchmarked
by its own authors.

---

# LOWCOV: the chaining floor, and why it is opt-in rather than a default

SEEDW (round-1 chain seed width) was 16, with MINOV floored at 16 because a seed
shorter than the overlap cannot find it. Both were set for normal coverage.
Swept together on Drosophila SRR40104920 (2.6x), everything else fixed:

| SEEDW/MINOV | archive | leftover_frac | second region |
|---|---|---|---|
| 16/16 | 36,851,932 | 0.812 | 134,787,263 |
| 12/12 | 36,783,282 | 0.810 | 133,729,577 |
| 10/10 | 36,558,158 | 0.802 | 133,000,521 |
| **8/8** | **36,103,845** | **0.785** | **128,485,041** |

Monotonic to the floor (8 is the hard limit -- "below this a seed is noise").
-2.38% against the tuned baseline, and it works at the SOURCE: more reads chain,
so the expensive second region shrinks by 6.3 M bases, rather than re-encoding
the result.

## Why it is NOT a default

On the three datasets checked it wins on H. salinarum (-1.33%) and Drosophila
(-2.38%) but LOSES on E. coli (+0.86%). Coverage does not explain the split --
H. salinarum is 25.8x and E. coli 50.6x, both normal. No validated property
separates them, and picking a threshold that happens to suit the 7 locked files
is precisely the overfit this project forbids: an eighth dataset could land the
wrong side of it.

So LOWCOV changes no default. With it off, all 7 locked datasets are
BYTE-IDENTICAL, which means an unseen dataset cannot regress -- a guarantee, not
a hope. With it on, Drosophila 2.6x is 36,103,845 B and verified LOSSLESS.

## Automatic detection was attempted and REJECTED

Gating this automatically needs a coverage signal available BEFORE round 1,
since that is the stage SEEDW controls. Two cheap estimators were built and
measured:

- 25-mer frequency histogram mode: returns 2 for EVERY dataset regardless of
  true depth (a fixed-size read sample has low sample-coverage for all of them).
- unique-kmer fraction: ranks P. aeruginosa at 26x (0.929) ABOVE Drosophila at
  2.6x (0.919), because it conflates genome size with depth.

Neither discriminates. The only signal that does separate the regimes is
leftover_frac (locked 0.222-0.518, Drosophila 0.812), but it is known only
AFTER round 1, so using it would mean re-running chaining. That restructuring is
not done, and the mode stays explicit rather than guessed.

## Low-coverage literature: the survey is now complete

Four tools, four architectures, no low-coverage analysis in any of them:

- SPRING / HARC (hash reorder): benchmarks 25-100x; singletons still 20-40% of
  archive after their realignment rescue.
- PgRC / PgRC2 (overlap graph): datasets 13-373x; crashes outright on ours.
- Leon (probabilistic de Bruijn graph): "degrades similarly to other de novo
  methods -- graph fragmentation is inherent to k-mer-level approaches".
- mstcom (Hamming-shifting MST): no coverage statistics, no low-redundancy
  analysis, no singleton fraction reported, and explicitly assumes a connected
  graph.

Every one silently avoids the regime.

## Position after this work

Drosophila 2.6x: 36,982,418 -> 36,103,845 B. Gap to SPRING 6.5% -> 3.97%.
The residual 1,380,005 B is architectural: SPRING matches read-to-read (~2 bits
per mismatch), we match region-to-region (~21 bits). Both routes around that are
measured and refuted above.
