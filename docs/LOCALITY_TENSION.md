# Claim 3's structural mechanism: the matcher optimises length and is blind to distance

**Status: measured, not implemented.** The measurements are reproducible from
`refs.tsv` (dump with `CAPS_DUMP_REFS=refs.tsv capsule_decode query <arc> /dev/null 0-1000`).

## The problem with the sidecar

`capsule_decode index` makes query 17x faster (0.51 s -> 0.03 s) by caching the
pseudogenome. It works, and it costs the archive nothing. But it is a **cache,
not a mechanism** — it sidesteps the reason query is slow instead of removing
it. An extension feature, not a contribution.

## The tension, stated

**Compression works by creating dependencies. Random access requires
independence.** In a self-referential compressor every byte of the
pseudogenome may depend on a byte far behind it, so answering a local question
requires a global decode. That is why *every* member of this family — PgRC,
NanoSpring, SPRING — offers compress and decompress and nothing in between.
It is not an interface omission. The representation genuinely does not support
a local answer.

This is the same SHAPE as Claim 2's finding, one layer down:

| | optimised objective | property silently destroyed |
|---|---|---|
| Claim 2 | bits per het site | the variant is no longer visible |
| Claim 3 | match length | locality is no longer available |

## The tension is FALSE, and here is the measurement

Two facts, both from the 240,762 references of an E. coli archive:

**1. Dependencies are long but SPARSE.** Reach-back distance is genuinely huge
(median 13,428,096; only 0.2% under 100 kb). But each reference carries ~85
bytes, so the transitive closure of a window is tiny:

    window   median closure   % of pg   depth
      1 kb          2,453 B     0.01%       8
     10 kb         19,182 B     0.07%       8
    100 kb        179,322 B     0.67%       8
      1 MB      1,728,600 B     6.41%      10

Closure is **~1.7x the window**, not ~1x the archive. The recorded conclusion
"99.7% of references reach back >100kb, so no windowed rebuild is possible"
confused *far* with *large* and is retracted.

**2. The distance is mostly GRATUITOUS.** Sampling 300 references and searching
the pseudogenome for a nearer exact occurrence of the same string:

    a NEARER exact copy exists within 200 kb : 132  (44.0%)
    the chosen source was already the nearest :   0  ( 0.0%)
    no copy within 200 kb                     : 168  (56.0%)
    median distance saved where one exists    : 14,933,040 bytes

**The matcher never picks the nearest.** It optimises match LENGTH and is
blind to distance, so 44% of references point ~15 MB further away than they
need to, at an identical string.

## Why fixing it is free — and probably cheaper

An exact match is an exact match. Referencing a nearer identical occurrence
reproduces the same bytes with the same length and the same (zero) mismatches.
Match quality is untouched. Only the coded source value changes.

And that value is our documented weak spot. `include/coders_inproc.h:292`
records sources as the single biggest loss against PgRC2:

    sources   ours 25.38 bits/match   PgRC2 23.83   +58,181 B

because `src` is coded as an absolute position in `[0, dst)`. A *near* source
delta-codes small. Modelling 44% of references relocated within 200 kb:

    today                      23.99 bits/match   721,913 B
    locality-aware             21.18 bits/match   637,456 B
    saving                      84,456 B  (0.12% of archive)
    PgRC2's figure             23.83 bits/match

So the change plausibly **closes the source-stream gap to PgRC2 and overtakes
it**, while collapsing the closure enough to make windowed access work without
any sidecar.

## What to implement

1. **Nearest-among-equals in the matcher** (`stages/106_inprocess.cpp`, the MEM
   candidate loop). It already enumerates candidates per seed and keeps the
   longest. Add: among candidates achieving the maximum length, keep the one
   with the smallest `dst - src`. This is a tie-break, not a new search, so it
   costs the loop nothing — the candidates are already in hand.
   **Note:** this is output-CHANGING, so the gate is archive size across the
   locked set, not byte-identity.
2. **Delta-code `src` as `dst - src`** in `refc` (`include/coders_inproc.h:292`)
   so the shorter distances actually pay off. Decoder already knows `dst`.
3. **Then windowed rebuild** per `docs/QUERY_WINDOWED_PLAN.md`, which becomes
   cheap once the closure is local, and retires the sidecar.

## Gates

- Archive size must not regress on any locked dataset; aggregate must improve.
  (Expected +0.12%, but the model is idealised — measure, do not assume.)
- Lossless round-trip re-verified by decoding the ARCHIVE on all 4 regression
  datasets.
- Closure re-measured after the change: the claim is that a 100 kb window's
  closure drops well below 0.67%.
- Only then is the "no sidecar needed" claim earned.

## The one-sentence version

> A self-referential compressor's matcher optimises match length and ignores
> match distance; 44% of its references therefore point ~15 MB further than an
> identical nearer copy. That gratuitous distance is what destroys random
> access — and removing it is free in match quality and cheaper to encode.
