---
Date: 2026-09-19
Title: T3.3 — Coordinate-Range Retrieval, All 19 Datasets
Purpose: Complete, verified record of what T3.3 measures, why it has no
  baseline column, and its precise code-level relationship to T3.4 (they
  are the same function).
When to refer to this file: Writing/checking T3.3's table or text;
  explaining why no competitor comparison exists for this table; explaining
  the T3.3/T3.4 code relationship when someone asks "why does T3.4 need an
  index but T3.3 doesn't."
Keywords: T3.3, range query, coordinate query, genocat, no competitor,
  19 datasets, rr vector, sidecar, qidx
---

# T3.3 — Coordinate-Range Query

## What it measures

Given a pseudogenome coordinate range `[start, end)`, return every read
placed anywhere within it. Run across **all 19 locked datasets** (not just
the 6-dataset T3.1/T3.2 subset — see `overview_claim3.md` and
`mechanism_insight_claim3.md` for exactly why the dataset counts
differ between tables).

## The mechanism — the simplest of the query-family operations

Inside `decoder.cpp`'s `mode=="query"` handler, a variable `rr`
(a vector of start/end pairs) gets filled, and everything downstream —
loading the sidecar, applying deviations/strand/N, emitting FASTA — is
shared with T3.4. For T3.3 specifically, `rr` is filled by **directly
parsing the input string as two integers**:

```cpp
const char* d=strchr(modearg.c_str(),'-');
if(!d){ fprintf(stderr,"query needs START-END or a DNA sequence\n"); return 2; }
rr.push_back({strtoull(modearg.c_str(),nullptr,10), strtoull(d+1,nullptr,10)});
```

No search happens. The coordinates ARE the input; the archive doesn't have
to find anything, only look up what it already knows about that span from
the placement stream. This is why T3.3 is described as "native" — there is
no search cost, no tolerance, no false-positive risk, none of the
complications T3.4 has to manage.

Exact code path and how it diverges from T3.4: `code_mapping_claim3.md`.

## Why there is no baseline column — precise, not hand-waved

The published note in `results/claim3/claim3_T3.1_T3.2_T3.3.csv` says:
*"no competitor exists for coordinate-range retrieval."* The precise reason,
worth stating carefully in the paper (this session verified it against real
literature, see below): it is **not** that no other tool does range
queries — CRAM/BAM indexed retrieval and Genozip's `genocat` both do range
retrieval by coordinate. The reason is that **coordinates only mean
something relative to a shared coordinate system**, and this archive's
coordinate system (pseudogenome position) is internal to itself — nobody
else's coordinates refer to the same numbering. Comparing "our range query
time" to "their range query time" would not be comparing the same operation
on the same input; it would be comparing two operations that happen to
share a name. This asymmetry is a real, documented finding in bioinformatics
literature (assembly evaluation without a shared reference is explicitly
called "ambiguous" in the field — see PAQman/Merqury discussion cited in
this session's research), not an excuse invented for this paper.

## Results — locked, all 19 datasets, timing only

Source: `results/claim3/claim3_T3.1_T3.2_T3.3.csv`, table=T3.3 rows (19
entries, `ours_s` column, `0.15s` to `118.48s`, no baseline columns
populated by design).

Representative range: `SRR29296997` at 0.15s up to `SRR10676752` at 118.48s
(the latter is a large plant genome dataset with a correspondingly large
pseudogenome to search through without the optional sidecar index).

## The critical code-level fact: T3.3 and T3.4 are NOT two mechanisms

This is worth repeating precisely because it is easy to misstate: **T3.3 and
T3.4 are the same `mode=="query"` code path in `decoder.cpp`.** The
only difference is what populates `rr`:

- T3.3: `rr` = literal parsed integers.
- T3.4: `rr` = output of `find_occurrences()` (pg-search), optionally
  unioned with the completion index k-mer completion.

Once `rr` is populated, the exact same emit loop serves both — deviations
applied, strand corrected, N restored, FASTA written. There is no separate
"T3.3 engine" and "T3.4 engine." This matters for how the paper's
architecture diagram should be drawn (see `docs/CLAIM3_ARCHITECTURE_DIAGRAM.md`)
and for correctly answering "why doesn't T3.3 need the completion index but T3.4 does" —
the answer is that the completion index only helps when `rr` comes from a *search*
(T3.4), because search has a tolerance ceiling that literal-coordinate
lookup (T3.3) never encounters in the first place.
