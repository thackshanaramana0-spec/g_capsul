---
Date: 2026-11-21
Title: T3.4 — Exact-Match Retrieval (Recall vs BWT-Family Tools)
Purpose: Complete, verified record of what T3.4 measures, how the FM-index
  completes exact-match search and why that is complete by construction (not
  a lucky statistic), the four-individual result, and its real, disclosed
  costs.
When to refer to this file: Writing/checking T3.4's table or text; deciding
  how to frame the comparison to BEETL/CIndex/sFASTQ (capability + cost,
  never head-to-head parity framing); citing any exact-match recall number.
Keywords: T3.4, exact match, recall, FM-index, BWT, BEETL, CIndex, sFASTQ,
  find_occurrences, completeness by construction, CAPS_QUERY_CONTAIN, caps_fmi.h
---

# T3.4 — Exact-Match Retrieval

## What it measures

Given a query DNA string (a "probe"), return every read that literally
contains it — the same question a BWT/FM-index answers by construction with
recall 1.0.

## Mechanism: search the pseudogenome, then the FM-index

Compression has already split the reads into two groups, and exact-match
retrieval uses that split.

**Group 1 — reads that match the pseudogenome exactly.** A read with no
deviations and no N is emitted byte-for-byte as a slice of the pseudogenome,
so any string it contains also occurs in the pseudogenome inside that read's
span. The pseudogenome search (`find_occurrences`, `decoder.cpp`) finds it.
This search is exact: recall for this group does not depend on a mismatch
tolerance.

**Group 2 — reads that deviate from the pseudogenome.** A read stored as
`pg[position] + its own deviations` (or carrying N, which the pseudogenome
never holds) can contain the probe where the consensus does not, so searching
the consensus cannot find it. These are 14–37% of the reads in the four
individuals here. They are indexed in an **FM-index** (`include/caps_fmi.h`),
built after compression from the placements and deviations the archive
already retains, as a separate file `<archive>.qidx.fmi`. Each read is
indexed exactly as `query` emits it (pseudogenome slice, deviations applied
at strand-corrected offsets, N restored). A read the archive reproduces
exactly but that `query` cannot emit as itself — several original reads
sharing one placement with different lengths, N patterns or deviations,
which occurs with variable-length input — is stored verbatim in the same
file and indexed like any other.

**Why this is complete by construction.** Every read is in group 1 or group
2. Group 1 is reached by an exact search of the pseudogenome, and group 2 by
an exact FM-index lookup. No tolerance, seed length or minimum probe length
enters either lookup, so a read that contains the probe cannot be missed, on
any data. The result is the union of both, so one call answers both "which
reads cover this locus" and "which reads contain this probe".

**Index.** The index is a BWT over the group-2 reads, built with SA-IS in
linear time (Nong, Zhang and Chan 2009), stored as 2 bits per base with rank
checkpoints every 512 rows, and searched by backward search. The index grows
with the number of reads that deviate from the pseudogenome, not with the
number of reads.

## Results — four individuals, chr20:3.0-3.6 Mb

Reads containing a probe, found by scanning every decoded read:

| Individual | Reads with a probe | Recall, archive | Recall, archive + FM-index |
|---|---:|---:|---:|
| HG002 | 12,059 | 0.8536 | **1.0000** |
| HG003 | 9,454 | 0.8457 | **1.0000** |
| HG004 | 12,455 | 0.8099 | **1.0000** |
| HG005 | 10,121 | 0.6047 | **1.0000** |

"Archive" queries the archive alone: each query decodes the pseudogenome and
searches it (mismatch tolerance 2). The archive stores each read as a
pseudogenome position plus its own deviations, and without the decoded
read-level information `query` returns the consensus under each read, not the
read itself. A read that carries deviations is therefore returned as a
different sequence, which does not count as recovering it. With the
FM-index, recall is 1.0000 at mismatch tolerance 0 and at tolerance 3 alike. The result files are
`results/claim3/claim3_T3.4_exact_match_recall.csv` and are produced by phase
4 of `scripts/benchmark_scripts/benchmark_1_run.sh`.

For reference, a query that also has the decoded read deviations, strand and
N positions (tolerance 2, no FM-index) reaches 0.9999, 0.9998, 1.0000 and
0.9990, because it returns each read as the read. What it still misses is a
read whose own deviations spell the probe where the consensus does not, and
those few reads, 1, 2, 0 and 10 for HG002 to HG005, are what the FM-index
adds.

## Tests

- `tests/fmi_selftest.cpp`: the suffix-array construction against a naive
  sort on 415 texts (single-symbol runs, periodic strings, tiny alphabets),
  and every FM-index lookup against a brute-force scan on 15,256 queries over
  63 read sets, including reads with N, reads longer than the locate sample
  interval, empty reads, and a save/load round trip.
- `tests/fmi_end_to_end.py`: eight synthetic archives built to break the
  index (variable lengths, contained reads, duplicates, N, repeats, high
  error). Each is compressed, checked lossless, indexed, and queried with
  probes of every length from 4 to 80 bp, on both strands, including probes
  with N and random strings. Every read containing a probe is returned, zero
  missed in all eight.
- Both run under `ctest`, and the self-test also runs inside
  `benchmark_0_preflight.sh`, so a broken index is a NO-GO before any table
  is produced.

## The real, disclosed cost — do not omit

1. **Index size.** 0.9 to 2.6 MB per individual: 928,092 B for HG002 and
   2,640,919 B for HG005, against archives of 480,627 B and 748,892 B. It is
   optional and built after compression (`CAPS_FMI=1` at `index` time), so it
   never affects Claim 1's compression numbers. It is 3.5 to 3.8 times
   smaller than the k-mer completion index it replaces (3.3 to 9.9 MB).
2. **Size scales with the share of reads that deviate from the pseudogenome.**
   Noisier data, or data where compression absorbs fewer reads, gives a
   larger index. This has not been measured beyond these four individuals and
   the synthetic archives.
3. **Build and query time.** The FM-index adds 0.2 to 0.6 s to the `index`
   step here. 800 queries take 7 to 12.6 s with the FM-index, against 13 to
   18 s with the k-mer index, at 8 queries in parallel. A first query
   without the decoded pseudogenome (`archive` above) decodes the archive
   each time: 90 to 208 s for 800 queries.
4. **False positives on homozygous negative controls.** Tolerant locus
   search also returns reads at homozygous negative-control sites (48, 52, 49
   and 87 of 400 for HG002 to HG005; see T3.5). That belongs to the locus
   query's mismatch tolerance, not to exact-match search, which is exact.

## How to frame the comparison to BEETL/CIndex/sFASTQ — do this, not that

**Do**: "The archive answers both retrieval types — reads covering a locus
(by position) and reads containing a string (recall 1.0000, via an FM-index
built from the archive after compression)." State capability and the real
cost. A BWT index answers the exact-match question by construction because
it is built from every read at compression time. Here the pseudogenome has
already absorbed the reads that match it, so the FM-index covers only the
rest.

**Do not**: claim head-to-head parity with BWT-family tools as if in
competition, or call the FM-index free. **No speed comparison and no
index-size comparison against BEETL/CIndex/sFASTQ has been measured.** Both
are explicitly future work, not claimed.

## What generalizes and why — not "we tested enough samples"

The completeness argument above is architectural: any dataset, any
individual, a read containing the probe is either a pseudogenome substring or
in the FM-index. The four-individual result is confirmation that it holds in
practice, and the synthetic suite probes the cases real data rarely produces
(contained reads, duplicates, N, repeats, 5% error). See
`mechanism_insight_claim3.md` for the same point made about T3.5, where the
argument is geometric and honestly weaker in that specific sense.
