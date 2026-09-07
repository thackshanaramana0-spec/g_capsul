# Archive-path caller: 128.5 s -> 71.9 s, 13.3 GB -> 10.5 GB

Target was 100 s. Final is **71.9 s**. Every number below is `capsule_decode
call` on the 4M-read chr20 subset archive (`sub/a.capsule`), full caller
(`CAPS_CALL_INDELS=1`), one job at a time.

## The gate every change passed

- SNV  F1 **0.4539**  TP 13858 / FP 2624 / FN 30717
- INDEL F1 **0.2613**  TP 1217 / FP 318 / FN 6564
- VCF 36,647 records; kc distinct 84,181,752; anchors 21,605,670
- placements 1,689,820 / 2,070,932 then 917,051 / 3,031,190

Not one of these moved across the whole session.

## What actually cost the time

| lever | stage effect | wall |
|---|---|---|
| hash-bucket directory over the seed index | 52 probes/lookup -> 2 | 128.5 -> 126.9 |
| exact branch-and-bound in the alignment scan | place 33.3 -> 11.3 s | 126.9 -> 105.1 |
| rolling canonical 25-mer at four walks | collapse 13.9 -> 7.3, read scan 12.8 -> 6.7 | 105.1 -> 93.1 |
| parallel in-memory quality-as-text decode | quality 11.0 -> 5.2 s | 93.1 -> 86.5 |
| contigs handed back from step 1 | export 5.6 -> 0.0 s | 86.5 -> 81.1 |
| parallel bin-wise merge for kc | merge 8.2 -> 0.7 s | 81.1 -> 73.4 |
| one-shot anchor directory (drops the Bloom) | read scan 6.6 -> 4.8 s | 75.1 -> 73.1 |
| word-wise skip in the alignment scan | place 11.3 -> 10.9 s | 71.8 -> 71.0 |

RAM, separately:

| lever | peak |
|---|---|
| claimed set sized from what goes in (4.29 GB -> 1.02 GB table) | no peak change |
| 16-byte pcluster anchor record | 13.30 -> 11.32 GB |
| 16-byte kidx record | 11.30 -> 10.51 GB |

## Two recurring shapes

**A container/search structure, not an algorithm.** Five of the eight wins are
the same move: something was looked up through a structure whose cost is
dependent memory probes (hash-map nodes, a binary search over tens of millions
of keys, a Bloom filter in front of a binary search, a priority queue fed by
12-byte freads). Replacing the structure while keeping the exact same contents
and ordering is bit-identical by construction. The largest single win in this
whole set -- branch-and-bound -- was NOT of this shape, which is why guessing
"it must be cache misses" was wrong twice.

**Work recomputed because two callers each did it from scratch.** The
pseudogenome was rebuilt twice per run (once to reconstruct reads, once to
export contigs) and the quality column was decoded to a file only to be read
straight back. Both were pure duplication visible in the log as a repeated
line.

## Measurement notes for whoever runs this next

- **Run-to-run noise is +/-1.5 s** (four consecutive baseline runs spanned
  71.32-74.11 s). A single before/after pair cannot resolve anything smaller.
  Every result above comes from two interleaved A/B/A/B runs, and stage timers
  are trusted over wall clock for sub-second effects.
- **A full disk corrupts this pipeline rather than failing it.** The k-mer
  counter spills to disk; a truncated spill silently drops k-mers. One run
  reported kc distinct 72,485,292 instead of 84,181,752, and its INDEL F1 came
  out 0.1722 -- which looked exactly like a bad code change and was not one.
  The harness now refuses to start below 40 GB free and asserts the kc identity
  after every run. Do not remove those checks.
- The `[INDEL-PROF] === accounted Xs of stage ===` line must match the
  `indel_pass` total. It has caught wrongly-placed timers four times.
