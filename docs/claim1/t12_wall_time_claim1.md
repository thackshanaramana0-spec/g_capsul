---
Date: 2026-09-19
Title: T1.2 — Wall Time vs SPRING and Genozip, 19 Datasets
Purpose: Complete, verified, honest record of T1.2's numbers — CAPSULE is
  never the fastest tool, on either axis, on any dataset. Stated plainly,
  not softened.
When to refer to this file: Writing/checking T1.2's table or text; citing
  any speed number for the base archive (not to be confused with Claim 3's
  export/coverage/query speeds, which ARE dominant wins); explaining the
  size-vs-speed trade-off.
Keywords: T1.2, wall time, compress, decompress, SPRING, Genozip, trade-off,
  honest loss, not mixed
---

# T1.2 — Wall Time

## What it measures

Compression and decompression wall time, same 19 datasets, same run, same
machine, one timed job at a time (manuscript's own stated protocol —
verified present in the caption).

## Results — a real, disclosed loss, not a mixed result

**Independently computed this session, directly from the raw CSV**
(`results/claim1/claim1_T1.1_T1.2.csv`), by comparing `compress_s` and
`decompress_s` per dataset across all three tools:

| | CAPSULE fastest |
|---|---|
| Compress, out of 19 datasets | **0** |
| Decompress, out of 19 datasets | **0** |

**CAPSULE is not the fastest tool at either compression or decompression on
any of the 19 datasets.** SPRING and Genozip trade the "fastest" position
between themselves row by row; CAPSULE never takes it. This should be
stated in the paper exactly this plainly — do not describe this as "mixed"
or "competitive," which would overclaim relative to the actual 0/19, 0/19
result.

**Sample row** (ERR5181310, cross-checked against the manuscript):
compress — CAPSULE 5.76s, SPRING **3.67s** (fastest), Genozip 4.54s;
decompress — CAPSULE 2.70s, SPRING 6.44s, Genozip **0.51s** (fastest).
Consistent with the aggregate 0/19 pattern.

## Why this is the honest cost of Claim 1's own design, not a flaw to hide

This is the direct, real price of the additional structure Claim 1 retains
and Claims 2/3 reuse: building a pseudogenome via overlap chaining,
tracking placements/deviations for every read, and running the full stream
selector (which measures multiple candidate encoders per stream and keeps
the smallest — see `code_mapping_claim1.md`) all cost real time that a
simpler compressor without this structure does not spend. **This is the
same trade-off already stated candidly elsewhere in the project's own
documentation** (`PIPELINE.md`/session history: "~1.7x slower, ~2.5x
heavier at worst" against PgRC2, on the earlier, narrower sequence-only
comparison — see `mechanism_insight_claim1.md`). T1.2 is the full-FASTQ,
19-dataset version of that same honest disclosure.

## Precisely where the RAM cost comes from, verified file-and-line

Verified against `docs/HEADROOM_ANALYSIS.md`. The dominant RAM item is a
single structure: `src/encoder.cpp`'s duplicate-detection map
(`std::unordered_map<uint64_t, std::vector<uint32_t>>`) costs **~80 B per
distinct read** (hash node + a separately heap-allocated inner vector).
PgRC2's equivalent — a flat `vector<uint32_t>` of sorted read indices,
compared adjacently for duplicates — costs **~4 B per read**, because the
*same* sorted array then also serves as the overlap-search index: **one
structure, two jobs**, where this project builds a hash map for dedup *and*
a separate prefix index for overlap. Measured headroom from this
difference alone: 90 MB on E. coli, 344 MB on L. major, 428 MB on yeast,
scaling to 2.2 GB projected on C. elegans and 4.4 GB on T. cacao — this is
the concrete, file-and-line explanation for the "~2.5x heavier" RAM
disclosure above, not a vague resource-usage gap. **Speed, by contrast, is
not where PgRC2 is ahead in assembly/comparison itself**: this project's
`rcmp` compares 32 bases per 64-bit operation at arbitrary offset, against
PgRC2's `compareSuffixWithPrefix`, which falls back to one base per
iteration through a lookup table whenever the suffix is not byte-aligned —
measured 1.64 s vs 2.75 s in this project's favor for that specific
operation. The overall slowdown is not from losing at the comparison
primitive; it is paid elsewhere (RAM-driven cache pressure, and the
decompress-side cost below).

## Most of the disclosed decompress slowdown is a Python/C++ artifact, not an algorithmic cost

Verified against `docs/FINAL_HEADROOM.md`. Of a measured 40.1 s total
decompress (against PgRC2's 4.0 s — the ~10x figure), **80% (32.19 s) is a
Python script** (`scripts/decode_105.py`) iterating every read in the
interpreter, and it also sets the peak memory (1,563 MB — higher than this
project's own *compress* peak of 1,032 MB). The C++ portion
(`capsule_decode`) alone takes 7.95 s. Ported fully to C++, decompress
would land near 7.9 s against PgRC2's 4.0 s — **roughly 2x, not 10x**.
This is disclosed precisely so a reader does not mistake an
implementation-language artifact (a reference decoder script never
optimized because the archive's own C++ decoder is what ships and is
measured for correctness) for an algorithmic limitation of the compression
architecture itself. **Compression-side speed, separately, has real,
structural headroom of its own**: the greedy assembly sweep is 43% of
compress time and "has never been attacked for speed — only for
correctness," per the same source, an explicitly open item.

## One real engineering fix behind the shipped speed number, worth keeping visible

Verified against `docs/CANDIDATE_FORK_DEPTH.md`: the encoder searches a grid
of `MAXMAP`×`MINOV` candidates and keeps whichever produces the smallest
archive (a genuine measured-optimum search, not a heuristic — see
`docs/REPRODUCIBILITY.md`'s discipline). A real inefficiency was found and
fixed: `MINOV` is consumed at the round-2 sweep, `MAXMAP` only later at
pigeonhole mapping, so candidates sharing the same `MINOV` were computing a
**bit-identical round-2 assembly** up to four times over — confirmed
directly from production logs, not inferred (`round2:` printed exactly two
distinct results across all 8 candidates, not 8). Measured cost: ~2,360 s of
CPU on a full HG002 run recomputing an assembly already in hand, and because
candidates ran concurrently at one thread each, the wasted work sat directly
on the wall-clock critical path. Fixed by sharing the round-2 result across
same-`MINOV` candidates rather than reforking per candidate. This is the
origin of `encoder.cpp`'s "candidates share a prefix" design
(`code_mapping_claim1.md`), not an incidental detail — it is a real,
measured, non-trivial fraction of what makes T1.2's wall-time cost as low as
it is, and belongs in Methods as evidence the speed trade-off was minimized,
not merely accepted.

## What NOT to write in the paper

- Do not write "G_CAPSUL is competitive on speed" or imply any row wins —
  verified 0/19 on both axes.
- Do not conflate this with Claim 3's speed tables (T3.1 export 129-784x
  faster, T3.2 coverage 16-54x faster) — those measure a *different*
  operation (recovering already-computed structure post-archive) against a
  *different* comparator (SPAdes, bwa+samtools+mosdepth), not
  compression/decompression against SPRING/Genozip. A reader confusing
  "G_CAPSUL is slow at T1.2" with "G_CAPSUL is slow at T3.1/T3.2" would be
  factually wrong in the opposite direction — T3.1/T3.2 are dominant wins.
