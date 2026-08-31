# Variable-length reads: where the field is, and what the design has to do

Measured 2026-08-29. Inputs: `fx100k.fq` (100k yeast reads, all 150 bp) and
`vl100k.fq` (the same reads truncated to 100-150 bp, `150-((i-1)%51)`).

## Where the field actually is

| tool | variable-length short reads | how |
|------|------------------------------|-----|
| **PgRC2** | **refuses** — `Unsupported variable length reads.` | hard check before any work; README states constant-length, <=255 bp |
| **SPRING** | full support, no degradation | same stream layout as fixed-length (`read_pos`, `read_noise`, `read_noisepos`, `read_rev`, `read_unaligned`) — verified by listing both archives |
| **Genozip** | full support | general-purpose per-field contexts |
| **ARCS** | full support, lossless | `pg_readlen` per-read length stream; quality path auto-switches to adaptive-CM when `var_len` is detected (`encoder.cpp:2694-2699`) |

**Correction to a common assumption:** SPRING's documented degradation ("`-l`
applies general-purpose BSC, gains may be lower") is its **long-read** mode. It
does not apply to variable-length short reads, where SPRING runs its full
assembly pipeline. Do not claim SPRING falls back here — it does not.

Measured on `vl100k.fq` (26.05 MB, 100k reads of 100-150 bp):

| tool | archive | compress | peak RSS | lossless |
|------|---------|----------|----------|----------|
| **ARCS** | **2,085,416** | 4.47 s | 582 MB | yes (verified round trip) |
| SPRING | 2,146,013 (summed tar members) | 1.97 s | 482 MB | — |
| Genozip | 2,935,349 | 0.65 s | 184 MB | yes (self-verified) |
| PgRC2 | — | — | — | **refuses to run** |

ARCS is 2.8% smaller than SPRING and 29% smaller than Genozip, at 12.49x. So
ARCS's variable-length support is real and competitive, not a checkbox. The
strong claim is only against PgRC2, which cannot process this input at all.

## What variable length actually breaks

Two effects, both measured on the same 100k reads:

| | fixed (150 bp) | variable (100-150 bp) |
|---|----------------|------------------------|
| unique after dedup | 93,350 (**6.65% collapsed**) | 99,978 (**0.03% collapsed**) |
| reads that are a prefix of a longer read | **0** | **7,371 (7.4%)** |

**These are the same phenomenon.** At fixed length, one read containing another
is only possible when they are equal, so exact-duplicate dedup catches all of
it. At variable length the two separate: a 3'-trimmed read is a strict prefix of
a less-trimmed read from the same start position, and exact-hash dedup misses
every one. The 6.6 percentage points that vanish from dedup reappear almost
exactly as the 7.4% now contained.

The current pipeline handles neither case:
- dedup is exact-hash equality (`fnv` + `==`), so a prefix is a different read;
- the sweep finds only suffix-prefix overlaps, never "B is a substring of A".

A contained read therefore cannot extend anything, fails to chain, and falls
through to the pigeonhole mapper. It is not lost — mapping places it at zero
mismatches, costing a position rather than a full read — but it has already
wasted a full sweep's work at every L level on the way there, and while it sits
in `admit` it can be selected as a chain partner ahead of a read that would
have extended further.

## Design

**Containment removal, before the sweep.** This is the textbook answer, not an
invention: Myers' string graph (2005) and SGA (Simpson & Durbin, 2012) both
remove contained reads before graph construction, for exactly this reason -- a
contained read adds no information the container does not already carry.

The insight that makes it cheap here: **containment subsumes dedup.** At fixed
length the two coincide, which is why exact dedup has been sufficient. Replacing
dedup with containment removal is therefore not an added stage but a
generalisation of one that already exists, and it costs nothing on fixed-length
input, where it degenerates to exactly the current behaviour.

Two tiers, in the order their measured value justifies:

**Tier 1 — prefix containment, via the sort that dedup already needs.**
Sort reads lexicographically; a strict prefix sorts immediately before its
extensions, so one linear scan of the sorted array finds every prefix
containment. Record `(container_id, offset 0, len)` and drop the read from
`admit`. O(n log n), no extra index, no new data structure.

This is the dominant case and not by accident: real variable-length data comes
from 3' quality trimming, so two reads differ in length only when they start at
the same genome position on the same strand -- which is precisely prefix
containment. **Covers the measured 7,371 of 7,371.**

**Tier 2 — interior containment (A occurs inside B, not at its start).**
Requires probing A's prefix seed against an index of *all* positions of all
reads, not just read starts. That index is ~128M entries on this input against
the current ~851k, so it is a real cost for a case Tier 1's measurement suggests
is rare. **Do not build it until measured.** The honest check is to run Tier 1
and count how many reads still fail to chain yet map at zero mismatches --
those are the interior-contained ones. If that count is small, Tier 2 is dead.

## What this does to the three axes

- **Size:** small and positive. Contained reads already cost only a position via
  mapping, so the gain is second-order -- removing them frees chain partners for
  reads that extend further, which shortens the pg.
- **Speed:** the clear win. 7.4% fewer reads through the sweep, and the sweep is
  now 2.2 s of a 7.3 s run.
- **Memory:** 7.4% fewer reads held.

**Predicted, not measured.** Tier 1 must be implemented and run before any of
these numbers are quoted.

## What NOT to do

- Do not claim SPRING degrades on variable-length short reads. It does not; only
  its long-read `-l` mode does.
- Do not build Tier 2 before measuring what Tier 1 leaves behind.
- Do not model variable length by truncating reads from the 3' end and then
  reading assembly quality off the result. That was tried: it drove
  both-sides-overlapped from 69.8% to 0.4% and produced 954k survivors, because
  truncation removes exactly the suffixes the overlap search needs. It measures
  the test, not the code. A fair test needs reads whose overlap structure
  survives -- real trimmed data, or truncation applied to a minority of reads.
