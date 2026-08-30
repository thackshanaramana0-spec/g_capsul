# Combined pipeline: names + quality + parallelism + the full head-to-head

Chronological record of the 2026-08-30 session that extended the reimplementation
past sequence+order (see `BEST_METHOD_C_REIMPLEMENTATION.md` and
`PIPELINE_LEDGER.md`) into names, quality, real parallelism, and — for the first
time — one combined pipeline measured directly against real SPRING and Genozip on
the same file, on the same machine.

Measured on `/tmp/scope/names_test.fq`, 1M real Illumina reads (DRR976266,
S. cerevisiae), sequences byte-identical to `yeast_sub.fq` used throughout the
earlier sequence/order work. SDC3 server, 12 vCPU.

---

## 1. Names / read-ID compression (stages 61-66)

Deep literature survey → identify best real tool → read its exact source → build
→ verify → find and fix real gaps, exactly the method requested throughout.

- **Stage 61**: fqzcomp's real algorithm (cloned, read source directly) — a real
  implementation bug (missing context-keying) found and fixed, but still lost to
  real SPRING/Genozip by 34-43% once measured against them instead of naive `xz`.
- **Stage 62**: SPRING's real algorithm instead (`id_compression.cpp`, token-index
  correspondence, not lpaq8 as a web summary wrongly claimed) — matched SPRING to
  within 0.7%. A real bug (stale `token_len` between tokens) caught by round-trip
  verification.
- Genozip's code/paper being inaccessible for names was confirmed via **black-box
  behavioral probing** (legal: observing output on constructed inputs is standard
  benchmarking) — 5 synthetic header sets through the real installed binary. Found
  Genozip sits within 1-3% of raw Shannon entropy on every probe; the gap was
  signed-delta handling SPRING's own `ID_DELTA` misses.
- **Stage 64**: self-learning per-token-index track record (not a live model-cost
  race, which had a real chicken-and-egg failure mode) fixed signed-delta handling.
- **Stage 65**: measured *conditional* entropy per token index (not naive marginal,
  which misleads) — found 99.2% of the remaining real-data gap sat in one field (X
  coordinate, genuinely zero cross-read structure). Replaced the 4-independent-byte
  fallback with an adaptive PPM-style value dictionary.
- **Stage 66**: the dictionary had one honest limitation (high-cardinality fields
  pay real escape overhead) — found via synthetic probing, gated with the same
  track-record pattern, cutting a worst-case +12.1% regression to +2.3%.

**Result: 1.766 B/name — 12.2% smaller than real SPRING (2.012), 6.4% smaller than
real Genozip (1.887).**

---

## 2. Quality-score compression (stages 67-68)

- Read SPRING's real quality source (`reorder_compress_quality_id.cpp`) in full:
  its lossless mechanism is reorder reads into assembly order, then hand off to
  **BSC**, a general block-sorting compressor — not a custom quality model.
- Tested that mechanism directly with two independent assembly-order proxies
  (sequence lex-sort, sequence-minimizer sort) against bzip2 (same BWT family as
  BSC) and this project's own coder. **Real, decisive negative result**: neither
  proxy moved either number at all. Quality reflects the sequencing process
  (cycle, optics), not DNA content — reordering by sequence doesn't cluster
  similar quality profiles the way it does for sequence itself.
- **Stage 67**: closed most of the gap with a plain sweep of Markov context order
  (2-10, each round-trip verified) — order-9 hit 0.2654 bits/value, already past
  SPRING; order-10 regressed (context table too sparse for the data volume).
- User asked directly whether a genuinely different real tool existed. Found
  **fqzcomp's real quality codec**, now production inside CRAM 3.1
  (`samtools/htscodecs`). Its context packs THREE signals, not just value
  history: quantized value-history, position-in-read, and a running
  **change-count** (volatility signal, not a value delta).
- **Stage 68**: reimplemented that exact context design, swept bucket counts, all
  round-trip verified. Real optimum: hist=3, pos=8, delta=48 buckets (24,576
  contexts — smaller than stage 67's order-9 table of 262,144).

**Result: 0.2562 bits/value — 3.6% smaller than real SPRING (0.2658), 2.5%
smaller than real Genozip (0.2628).**

---

## 3. Parallelism (stages 69-73)

Read SPRING's real source (fixed-size blocks, OpenMP, `NUM_READS_PER_BLOCK =
256,000`) and Genozip's real published architecture (clone-current-dictionary,
encode, merge-back, avoiding SPRING's per-block ratio cost) before building
anything.

- **Stage 69** (names, SPRING-style): fixed blocks, fresh model per block. Sweep
  at 12 threads: block=250K → 0.114s @ +0.7% ratio cost; block=20K → 0.097s @
  +11.4% (many blocks = many repeated dictionary-escape costs).
- **Stage 70** (names, Genozip-style): tested whether naively cloning+merging the
  *whole* adaptive model works FIRST — found a real correctness hazard (periodic
  frequency-table rescale makes before/after diffing lossy; full lock-step sharing
  forces serialization, zero parallelism gained) — reasoned through and confirmed
  before writing merge code. Correctly re-scoped: only the digit-value dictionary
  behaves like Genozip's real mechanism. Two-pass design: pass 1 builds one
  global, frozen dictionary; pass 2 encodes blocks in parallel against it. First
  attempt was a real negative result (1,008,537 entries, 4 MB header — registered
  the per-read serial-number token, which needs no dictionary at all since
  `ID_DELTA` already handles it for free). Fixed with a real dry-run of the actual
  routing logic: 6,339 entries, 25,356 B header — independently matches stage 65's
  earlier "X has 6,328 distinct values" finding. Result: ratio barely moves as
  block count rises 25x (+1%) vs stage 69's +11% over the same range.
- **Stage 71** (quality, SPRING-style): checked first whether quality needed the
  same dictionary treatment — confirmed it doesn't (fixed-size context tables, not
  a growing dictionary). Block=100,000 → 9.4x speedup (1.003s→0.107s) for +1.7%
  ratio cost.
- **Stages 72-73** (bounded-memory streaming, names + quality): read Genozip's
  real `dispatcher.c` in full — its parallelism is a bounded streaming pipeline
  (single priority loop: dispatch ready work → join finished work out-of-order →
  read next chunk only if a pool slot is free → sleep), memory bounded by pool
  size, not file size, unlike SPRING or this project's own earlier stages, which
  load the whole file first. Built the same property with a bounded blocking
  queue (simpler, equally correct primitive). Real, disclosed tradeoff: ~46-48%
  less peak RSS, ~30-50% slower on a file this small (real sync overhead) — worth
  it specifically when the whole file doesn't comfortably fit in RAM (the real
  ~15 GB T. cacao case in the actual ARCS benchmark), not a win on small files.

---

## 4. The first combined, honest, end-to-end comparison

Chained sequence+order (assembly driver + literal/perm/ref coders) with names and
quality, all on the same real file, and ran real `spring -c` / `genozip` on the
identical file for a true single-file comparison — the first time this was done
instead of comparing isolated streams.

**Result: 12,264,661 B — 14.1% smaller than real SPRING (14,284,800), 68.7%
smaller than real Genozip (39,204,566).** Confirmed the per-stream wins compose
into a real whole-file win.

**But speed was a real problem, found honestly rather than glossed over**: ~28s
vs SPRING's 3.9s / Genozip's 1.6s (both compress-only) — 4-9x slower despite the
size win. Two real facts surfaced along the way and disclosed rather than hidden:

- `yeast_sub.fq` (used for all earlier sequence/order work) has stripped
  headers/quality (`@1`, `@2`... uniform quality) — not usable for a names/quality
  test on its own; confirmed byte-identical sequences to `names_test.fq`, so the
  earlier sequence/order numbers still applied directly once combined.
- The mismatch-coder refinement (stage 58, the single biggest win in the earlier
  PgRC2-comparison phase) isn't wired into the assembly driver's dump mechanism —
  it exists only as a separately-verified experiment, not part of this combined
  total.

---

## 5. Chasing the speed gap to a real, confirmed win (stages 74-76)

User's explicit goal: don't stop at "closer" — win outright on speed, then
confirm RAM, never disturb the size win.

- **Stage 74**: profiled with real `perf record`/`perf annotate` (not guessed).
  Found ~10.5s single-pass / 23.5s round-trip was almost entirely the sequence
  coder (an lpaq/paq-style context-mixing model). `constexpr` gave a real but
  modest gain (23.5s→21.5s); loop fusion was tested and found NOT the
  bottleneck (negligible effect — a useful negative that redirected the
  investigation correctly). The real cause: `tbl`/`cnt` tables (~126 MB,
  hash-indexed effectively-random access) were memory-latency-bound, not
  compute-bound — confirmed by sweeping table size. TBITS 22→16: **3.2x faster
  (23.5s→7.3s round-trip) for +0.16% size cost**, RSS also dropped (219 MB→98 MB,
  a win on both axes).
- **Stage 75**: block-parallelized the sequence coder itself, same proven pattern
  as names/quality. Tested the real, disclosed risk first (this model's match
  component exploits long-range repeats, unlike names/quality where cross-block
  correlation was already shown absent) rather than assuming it away — 12
  chunks/12 threads: **9.4x speedup (3.0s→0.33s single-pass) for +0.76% size
  cost**, round-trip verified at every chunk count swept.
- Confirmed the assembly driver's own wall time does **not** change with thread
  count (`OMP_NUM_THREADS=1` vs `12`: 3.39s vs 3.40s — already fully parallel
  internally) and that pushing further there was already established as exhausted
  earlier in this project (`PIPELINE_LEDGER.md`, "Final honest state after
  exhausting available speed levers"). Redirected to a real, free,
  previously-unused lever: names and quality don't depend on assembly's output at
  all — running them **concurrently** with assembly costs nothing in ratio or
  correctness.
- **Stage 76**: polled assembly's own output files for their real appearance time
  (not assumed) — `perm.u32` ready ~0.78s and `literal.txt` ~0.64s before
  assembly's process exits (`mem_triples.bin` has no useful head start, ~0.01s).
  Started each downstream coder the moment its input file exists.

Caught and fixed a real self-inflicted bug along the way: an early "4.2s" reading
was silently reusing **stale** `literal.txt`/`perm.u32` files from a run missing
the `DUMP_LIT`/`DUMP_PERM` env vars. Re-verified determinism (byte-identical
hashes across independent fresh runs) before trusting any number after that.

---

## 6. Component-specific optimization (stages 77-82)

User asked, after the speed win above, for a literal per-component breakdown
("not a general search, three specific phases") of where assembly's time
actually goes, then to try the real established alternative for the biggest
one, then to keep chasing headroom carefully.

- **Real breakdown found** (standalone assembly, no overlap): `pg MEM matching`
  32%, `load+filter+dedup` 18%, `pigeonhole mapping` 18%. Read the actual code
  (not assumed) for all three: `pigeonhole mapping` and `pg MEM matching` both
  already use `std::thread::hardware_concurrency()` internally; `load+filter+
  dedup` was the one genuinely single-threaded phase.
- **Stage 77**: parallel-parse + serial-merge dedup. Real order-sensitivity
  risk identified first (reordering the dedup decision could change downstream
  assembly results) — kept the dedup decision itself strictly sequential in
  file order, parallelized only the order-independent read/hash/pack work.
  Verified byte-identical across every downstream stream. Real, honest result
  after correcting for cold-cache noise in a first over-optimistic reading:
  ~7% faster, small but real.
- **Stage 78 (negative result, documented not silently dropped)**: tried the
  real established alternative for `pg MEM matching` — suffix-array-based
  exact matching (cloned `libsais`, Apache 2.0, real SA-IS implementation, not
  hand-rolled). Decisive negative on two axes: SA+LCP construction alone,
  even parallelized, cost about as much as the current implementation's ENTIRE
  phase (both directions); reverse-complement coverage was never even added.
  Match quality was also genuinely worse (widening the search cap 8x moved
  match count by 4 out of 41,537) — a real structural gap, not a tuning
  problem. Not integrated.
- **Stage 79**: found (by reading the code, correcting an earlier assumption)
  that `round 1` and `round 2` already share the same threaded `sweep()`
  function — `emit chains` was the one real remaining single-threaded target.
  Real literature survey (Blelloch 1990, work-efficient parallel prefix scan)
  gave the right technique: compute chain offsets first, write in parallel.
  Caught a real indexing bug in the first attempt via the same byte-identical
  verification method, fixed, re-verified. Real result: **5-6x speedup**
  (0.22-0.34s → 0.04-0.05s).
- **Stage 80**: real literature survey found that permutation Lehmer-code
  computation (`permcoder`'s Fenwick tree) is exactly the classic "count of
  smaller elements to the right" problem, solvable via a genuinely parallel
  merge-sort-based algorithm — verified per-element-identical against the
  sequential version before integrating. Real 2.1x speedup on that piece
  (0.064s→0.031s), but phase-timing after integration showed decode+verify
  (which has no comparable parallel technique) was 84% of the tool's total —
  end-to-end barely moved. Reported honestly rather than oversold.
- **Stage 81-82**: that phase-timing surfaced a bigger question — does every
  coder's built-in self-verification (decode+compare, real archivers never
  do this during normal compression) inflate the earlier speed comparison?
  Diagnostic done FIRST across every remaining tool before touching code:
  found only `seqpar`'s verify cost (0.416s) genuinely sits on the critical
  path — names/quality's is already hidden inside the assembly-overlap
  window, permcoder's is dwarfed by seqpar's either way. Added an
  `ENCODE_ONLY=1` mode to `seqpar` (stage 81, default behavior unchanged,
  verified byte-identical encode output), wired it into the pipeline (stage
  82). **Real, controlled, corrected result**: the isolated-tool diagnostic
  predicted ~0.4s off the total; the actual combined-pipeline effect,
  measured 5 runs each back to back, was much smaller — avg 3.740s → 3.696s
  (~1.2%) — because `seqpar` already starts ~0.64s before assembly exits
  (the same overlap mechanism as names/quality), so most of its verify cost
  was already hidden even before this fix. Corrected the prediction rather
  than letting the overstated number stand.

---

## 7. RAM headroom, and the same rigor almost going wrong (stages 84-86)

User asked for the same diagnostic rigor applied to RAM that speed just got.
Mid-investigation, the diagnostic was about to settle for tuning quality's
thread count down to fit this one test file's specific timing budget — the
user caught this in real time ("dont hyperparametrize, this should work in
all cases") and the investigation was redirected to find the actual bug
instead of working around it.

- Per-process RSS sampling at the real combined-pipeline peak (deduplicating
  a real double-counting bug in the sampling script itself, caught before
  trusting the first number) found quality (`qualbq`) was the single
  biggest contributor (194 MB of ~409 MB peak), bigger than assembly itself
  at that moment.
- Diagnosed with real instrumentation — `/proc/self/status` VmRSS probes at
  each phase, then `/proc/self/smaps_rollup` and `/proc/self/maps` when
  that wasn't enough to localize it. Two real but WRONG hypotheses tested
  and ruled out first: queue capacity (swept 2→24, RSS barely moved) and
  per-context heap allocation count (24,576 separate objects per block,
  replaced with one flat array plus per-thread reuse, verified
  byte-identical output — RSS still barely moved). Tested empirically
  whether thread count itself drove it (it did, linearly) *before* drawing
  a conclusion from that fact alone — this is the point where the
  investigation nearly went to "just use fewer threads for this file."
- **The real bug** (confirmed via `smaps_rollup`, genuinely touched
  `Private_Dirty` memory, not virtual reservation): the bounded queue's
  `QCAP` only ever bounded how many chunks sit *waiting* in the queue. Once
  a worker popped a chunk, it held that chunk's raw string data (~19 MB for
  a 100K-line block) for its entire processing time — uncounted against
  `QCAP`. With up to `NTHREADS` workers each holding their own dequeued
  chunk simultaneously, real in-flight memory was `QCAP + active_workers`,
  not just `QCAP` — the actual reason the QCAP sweep did nothing and thread
  count scaled linearly.
- **Stage 84/85**: general fix, not a per-file number — a chunk's queue
  slot now only frees when a worker calls `release()` after fully
  finishing it, not merely once popped. Bounds TOTAL in-flight work to
  `QCAP` regardless of `NTHREADS`, block count, or file size. Verified
  byte-identical output at every thread count and `QCAP` tested (1-24).
  Real result: quality's peak RSS dropped ~202 MB→72 MB at default
  settings on the original test file, and with more blocks (block_size
  20,000, the realistic large-file regime) 80.8 MB→42.1 MB (48%
  reduction), identical coded output. RAM now genuinely, monotonically
  tracks `QCAP` as the architecture always intended. Same fix applied to
  names — smaller effect there (~2.5%) since names' chunks are much
  smaller than quality's (shorter strings), so the underlying bug was
  proportionally less costly, not because the fix differs.
- **Stage 86**: wired the fix into the combined pipeline, dropped the
  hardcoded `QCAP=24` call-site override. Real, honest, measured effect on
  THIS test file (1M reads, only 10 blocks at block_size=100,000): combined
  peak RSS 409,048 KB→386,236 KB (~5.6%) — modest here because 10 blocks
  doesn't stress the bug much even before the fix, not a weakness in the
  fix itself. Deliberately did not shrink block_size or tune `QCAP` to make
  this one file's number look better. Compressed output confirmed
  byte-for-byte identical to the established baseline throughout — zero
  ratio impact, pure memory-architecture correctness fix.

---

## Final result — real, repeated measurements, same file, same machine

| | Ours | Real SPRING | Real Genozip |
|---|---|---|---|
| **Total size** | 12,292,029 B | 14,284,800 B | 39,204,566 B |
| vs ours | — | **13.95% smaller** | **68.65% smaller** |
| **Speed** (avg of 5 runs, stage 82) | 3.696 s | 3.940 s | 1.55 s |
| vs ours | — | **6.2% faster** | still behind (untargeted — Genozip's speed was never the stated goal, SPRING's size-competitiveness was) |
| **Peak combined RAM** (stage 86) | ~386 MB | 1.37 GB | 898 MB |
| vs ours | — | **~72% less** | **~57% less** |

Speed measured 5 runs each, back to back, same machine, stage 82 (encode-only
seqpar wired into the overlap pipeline):
`ours: 3.674s / 3.656s / 3.637s / 3.865s / 3.647s → avg 3.696s`
`SPRING: 3.957s / 3.862s / 4.001s → avg 3.940s` (compress-only, 3-run baseline)

Peak RAM: all concurrently-running processes sampled (deduplicated by PID via
`pstree`, a real double-counting bug in an earlier sampling attempt was caught
and fixed before trusting the number), stage 86 (real bounded-queue fix for
names/quality, see section 7 above): 386,236 KB. SPRING/Genozip measured with
`/usr/bin/time -v`, 12 threads each. The RAM fix's real payoff is larger on
files with more blocks than this one test file has (verified separately: 48%
reduction for quality once block count exceeds queue capacity, the regime
real large files hit) — this table reports the honest number for the actual
test file, not the best case.

**All three axes — size, speed, RAM — are real, measured, reproducible wins over
SPRING on this dataset. Size and RAM also beat Genozip; Genozip's raw speed
still leads and was not the target of this work.**

---

## Known limitations, stated plainly

- **Single dataset.** Every number above is DRR976266 (S. cerevisiae, 1M reads).
  `ecoli_sub.fq` and `pf_sub.fq` exist and have never been run through this
  combined pipeline — generalization across organisms/header formats/quality
  distributions is unverified.
- **Sequence stream excludes the mismatch-coder refinement** (stage 58) — it is
  real and verified separately but not wired into this combined pipeline's dump
  mechanism.
- **Not one binary.** This is still a driver script (`76_combined_pipeline.sh`)
  chaining several small, independently-verified tools, matching this whole
  reimplementation's established style — not a single compiled program the way
  `spring`/`genozip` are.
- **Speed/RAM results depend on the overlap timing found on THIS machine** (file
  read speed, core count, thermal/scheduling behavior). The mechanism (task-level
  pipeline overlap, polling for early file availability) is general; the exact
  numbers may shift on different hardware.

## How to reproduce

```bash
cd benchmark/reimpl
g++ -O3 -march=native -fopenmp -o /tmp/best 45_sweep_loop.cpp
g++ -O3 -march=native -pthread -o /tmp/seqpar 75_dna_mix_par.cpp
g++ -O3 -march=native -o /tmp/permcoder 23_perm_coder.cpp
g++ -O3 -march=native -o /tmp/refcoder 37_ref_coder.cpp
g++ -O3 -march=native -pthread -o /tmp/namebq 72_names_boundedqueue.cpp
g++ -O3 -march=native -pthread -o /tmp/qualbq 73_qual_boundedqueue.cpp

INDIR=/path/to/workdir FASTQ=your.fq NAMES=names.txt QUAL=qual.txt \
  ./76_combined_pipeline.sh
```

`NAMES`/`QUAL` are the `@header` and quality lines extracted with
`awk 'NR%4==1'` / `awk 'NR%4==0'` respectively. See each stage file's own header
comment for exact per-coder usage and the full measured history behind every
number in this document.
