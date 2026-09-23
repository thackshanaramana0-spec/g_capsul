---
Date: 2026-09-19
Title: T3.4 — Exact-Match Retrieval (Recall vs BWT-Family Tools)
Purpose: Complete, verified record of what T3.4 measures, the completion index
  mechanism and why it is a completeness proof (not a lucky statistic), the
  full multi-individual generalization result, and its real, disclosed
  costs.
When to refer to this file: Writing/checking T3.4's table or text; deciding
  how to frame the comparison to BEETL/CIndex/sFASTQ (capability + cost,
  never head-to-head parity framing); citing any exact-match recall number.
Keywords: T3.4, exact match, recall, completion index, k-mer index, BWT, BEETL, CIndex,
  sFASTQ, find_occurrences, completeness proof, CAPS_QUERY_CONTAIN
---

# T3.4 — Exact-Match Retrieval

## What it measures

Given a query DNA string (a "probe"), return every read that literally
contains it — the same question a BWT/FM-index answers by construction with
recall 1.0.

## Two-stage mechanism: plain search, then completion

**Stage 1 — plain pg-search (`find_occurrences`)**: searches the 2-bit
packed pseudogenome for the probe, tolerating up to `k` mismatches via
pigeonhole seed-and-extend, where `k_max = floor(probe_len / MINSEED) - 1`.
`MINSEED` is *derived from the haystack size*
(`clamp(ceil(log4(|pg|))-2, 8, 16)`), not a fitted constant — this matters
for the generalization argument below.

**Why plain search alone cannot reach recall 1.0000 in general**: a read is
stored as `pg[position] + its own deviations`. If a read's deviations fall
inside the probe window, the read's actual letters differ from the
consensus there — searching the *consensus* cannot find that read, even
though it genuinely contains the probe. Only reads carrying deviations can
be missed this way (~14% of reads, since a deviation-free read is
byte-identical to a pg substring and is always findable).

**Stage 2 — the completion index completion** (`CAPS_QUERY_CONTAIN=1`): a k-mer index
built at `capsule_decode index` time, over exactly the deviation-carrying
reads, indexing their *actual reconstructed sequence* rather than the
consensus. This is not approximate matching — it is a genuine **completeness
proof**: a probe of length `P >= k+stride-1` (27bp here) is guaranteed to
contain at least one stride-aligned indexed k-mer, so any read genuinely
containing the probe is guaranteed to surface as a *candidate*. Every
candidate is then verified by reconstructing the read in full (deviations,
strand, N all applied) and doing a literal substring test — confirmed
directly in the code this session (`decoder.cpp` ~line 705-720). A
genuine match cannot be missed, on any dataset, given the code is correct —
which is a stronger guarantee than "we tested N samples and it worked."

The result returned is the **union** of both stages, so one call answers
both "which reads cover this locus" and "which reads contain this probe."

## Results — locked, multi-individual, generalized

Source: `docs/T34_T35_MULTI_INDIVIDUAL_20260919.md`.

**First pass, plain query alone (no the completion index)** — this revealed a real
generalization gap:

| Individual | Recall, plain query |
|---|---|
| HG002 | 1.0000 |
| HG003 | 0.9619 |
| HG004 | 0.9612 |
| HG005 | 0.8563 |

HG002's clean 1.0000 was **luck specific to that test set** — no read in it
happened to exceed the tolerance ceiling. HG003/4/5 were not so lucky.

**With the completion index completion (`CAPS_QUERY_CONTAIN=1`) — the general, correct
configuration:**

| Individual | Recall |
|---|---|
| HG002 | 1.0000 |
| HG003 | 1.0000 |
| HG004 | 1.0000 |
| HG005 | 1.0000 |

Confirmed by directly ruling out confounds first (binary choice, mismatch
tolerance setting each moved recall by ~0.0001 — not the explanation).

## The real, disclosed cost of the completion index — do not omit

1. **Index size: 7-13x the archive itself.** HG005: archive 748,892 B,
   the completion index 9,925,044 B. This is a real number to put in a cost table, not
   hide. It is optional and built post-archive (`CAPS_XMI=1` at `index`
   time), so it never affects Claim 1's compression numbers.
2. **False positives on homozygous negative controls rise ~35-75%
   relative** when the completion index is engaged (HG002 31→42, HG003 38→52, HG004
   28→49, HG005 71→87 — measured in the T3.5 bilateral test, which shares
   the same underlying mechanism). Completeness goes up; specificity goes
   down. Real trade-off, not a free win.
3. **Minimum probe length 27bp.** Shorter probes are refused loudly rather
   than silently degraded.

## How to frame the comparison to BEETL/CIndex/sFASTQ — do this, not that

**Do**: "The archive answers both retrieval types — reads covering a locus
(native) and reads containing a string (recall 1.0000, via an index derived
from the archive after compression)." State capability and the real
the completion index cost. This matches the reciprocity argument the paper already makes
(see `mechanism_insight_claim3.md`).

**Do not**: claim head-to-head parity with BWT-family tools as if in
competition. That framing invites "why not just use a BWT then," which has
a real answer (see file 06) but is better made in discussion, not implied
by table framing. **No speed comparison and no index-size comparison against
BEETL/CIndex/sFASTQ has been measured.** Both are explicitly future work,
not claimed.

## What generalizes and why — not "we tested enough samples"

The completeness proof (stage 2 above) is architectural: any dataset, any
individual, given the code is correct, a genuinely-containing read cannot
be missed once the probe clears the length bound. The 4-individual result
above is confirmation of that proof holding in practice, not the source of
confidence in it — this distinction matters when writing the generalization
claim in the paper (see `mechanism_insight_claim3.md` for the same
point made about T3.5, where the underlying argument is geometric rather
than a completeness proof, and is honestly weaker in that specific sense).
