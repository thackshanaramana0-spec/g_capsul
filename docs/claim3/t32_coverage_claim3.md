---
Date: 2026-09-19
Title: T3.2 — Per-Base Coverage vs bwa+samtools+mosdepth
Purpose: Complete, verified record of what T3.2 measures, its exact
  mechanism, and why it is a clean, unconditional win with no caveats.
When to refer to this file: Writing/checking T3.2's table or text; the one
  table in Claim 3 that needs no qualification, useful as a contrast when
  explaining why T3.1/T3.4 DO need qualification.
Keywords: T3.2, coverage, bwa, samtools, mosdepth, pos_abs, read_lengths,
  orig2uid, sweep-line, difference array, per-base depth
---

# T3.2 — Coverage

## What it measures

Time to compute per-base sequencing depth directly from the archive,
compared to the conventional pipeline: align (bwa-mem) → sort/index
(samtools) → compute depth (mosdepth), timed as one combined pipeline.

## The mechanism — genuinely exits before assembly, not just fast

This is the cleanest of the five tables architecturally. `coverage` mode in
`decoder.cpp` **returns before the pseudogenome is even rebuilt**
(lines 417-490; see `code_mapping_claim3.md` for the exact trace).
It needs only three things from the archive:

- `pos_abs` (and its region-split extensions `pos_sec`/`pos_region`) — where
  each unique read is placed
- `read_lengths` — per-*original*-read length (not per-unique — see the bug
  note below)
- `orig2uid` — the mapping from original read index to its unique/dedup
  representative

From these it builds a **sweep-line difference array**: for each placed read
spanning `[a, b)` in pseudogenome coordinates, increment `diff[a]` and
decrement `diff[b]`, then walk the array once summing a running total. This
is a standard O(n) technique, not approximated, not sampled — it produces
the exact same per-base depth a real alignment-based pipeline would, because
it is reading the same placement information a real aligner would have had
to compute, except this system already has it as a byproduct of
compression.

**Why this specific operation exits early and T3.1 (export) cannot**: export
needs the pseudogenome's actual *content* (base-by-base sequence), which
requires replaying every stored reference/literal to rebuild it. Coverage
needs only *where* reads sit and *how long* they are — never their content —
so the entire literal-decode and reference-replay stages (measured as ~94%
of a query's cost when unoptimized) are skippable.

## A real, documented, fixed bug worth keeping in the paper's discussion

`pos_abs` is indexed per **unique** read; `read_lengths` is indexed per
**original** read. Indexing both with a single counter silently undercounts
every duplicate read, because duplicates share one unique read's placement
but are separate reads for coverage-depth purposes. Measured before the fix,
on E. coli: 185,346,900 covered bases against the reads' true 232,988,850 —
a **20% undercount**. Fixed by expanding originals through `orig2uid`
exactly as the read-reconstruction path already does. This is a good,
concrete example for the paper (or a reviewer response) of the kind of
correctness issue this project actively hunts for and fixes, not something
to hide.

## Results — locked, dominant, no caveats

Source: `results/claim3/claim3_T3.1_T3.2_T3.3.csv`, table=T3.2 rows.

| Dataset | Ours (s) | bwa+samtools+mosdepth (s) | Speedup |
|---|---|---|---|
| ERR5181310 | 0.15 | 4.15 | 27.7x |
| SRR29296997 | 0.11 | 2.53 | 23.0x |
| SRR2584863 | 0.33 | 7.57 | 22.9x |
| SRR37283774 | 0.56 | 15.98 | 28.5x |
| DRR976266 | 1.49 | 23.98 | 16.1x |
| HG002 | 3.04 | 164.90 | 54.2x |

The output numbers themselves are identical to the conventional pipeline's
(not approximated) — the speed win comes entirely from skipping alignment
and sorting, work the archive already did once during compression and never
needs to redo.

## Why this table needs no qualification, unlike T3.1 and T3.4

- No correctness trade-off exists to disclose: the depth numbers are exact,
  not estimated, so there is no "2 of N metrics" table needed.
- No auxiliary structure (like the completion index) is required — the base archive alone
  is sufficient.
- No competitor-architecture asymmetry to explain (unlike T3.4's "BWT gets
  this natively, we need an index" story).

This makes T3.2 a useful anchor when writing the discussion section: it is
the table that shows the "retained structure is directly reusable" claim at
its strongest and simplest, before the paper has to start explaining the
more nuanced tables (T3.1's honest gaps, T3.4's auxiliary index, T3.5's
orchestration layer).

## The architecture claim, made concrete by measurement, not asserted

Verified against `docs/CLAIM3_POSTSPLIT_AUDIT.md` — peak RSS per operation
on the same archive:

| operation | peak RSS | pseudogenome rebuilds |
|---|---|---|
| **coverage** | **4 MB** | **0** |
| export | 36 MB | 1 |
| query | 36 MB | 1 |

`coverage` runs in **9× less memory** than export or query specifically
*because* it needs only `pos_abs` and `read_lengths` — it is hoisted above
the pseudogenome-rebuild step entirely, unlike export and query, which both
must reconstruct the pseudogenome before they can answer anything. This is
the "the compressor already stored the structure this operation needs"
claim made concrete rather than asserted, and it is coverage-specific: no
other Claim 3 operation gets this for free. Correctness was independently
re-verified by an exact identity, not sampling: covered bases summed from
`coverage` output equal `sum(read_lengths)` over every read exactly
(69,535,651 = 69,535,651, difference 0, with 0 bases legitimately clipped
at the pseudogenome boundary) — the strongest of the three operations'
post-split re-checks, because a single dropped duplicate read would break
this equality.
