# Caller RAM and speed — skeletal analysis

Written 2026-09-03 while the first full-chr20 HG002 run was in flight, from
the code plus that run's own phase instrumentation. Every structure size below
is derived from **measured** full-scale counts (12,604,917 reads; 451,578 →
233,369 contigs), not from a window extrapolated by hand.

Scope: `include/caps_caller.h` only. The compressor is not the problem and is
not discussed except to bound it — see §0.

---

## 0. The compressor is not the cost centre. Measured, this run.

| stage | time | RSS | peak |
|---|---|---|---|
| load+filter+dedup | 16.6 s | 706 MB | 1,430 MB |
| prefix seed index | 1.0 s | 942 MB | 1,430 MB |
| round 1 (division) | 63.6 s | 1,399 MB | 1,783 MB |
| round 2 (assembly) | 25.6 s | 1,399 MB | 1,783 MB |
| emit chains | 7.1 s | 1,672 MB | 1,783 MB |
| pigeonhole mapping | 102.7 s | 1,498 MB | **2,353 MB** |
| **total** | **216.7 s** | | **2.35 GB** |

The entire assembly the caller consumes costs 3.6 minutes and 2.35 GB. Every
RAM and speed problem in this document is the caller's, and the "one pipeline"
reuse saves 216.7 s — real, but 5–9% of a 40–70 min caller run, not the
headline. This is stated here so no future doc reads the 2.35 GB and the
30 GB as one number.

---

## 1. Parallelism: one pragma in 2,535 lines

    $ grep -n "pragma omp" include/caps_caller.h
    1100:    #pragma omp parallel for schedule(dynamic, 64)

That is the whole of it. On a 12-vCPU server, five of the caller's six phases
are single-threaded:

| phase | line range | time (full chr20) | threaded? |
|---|---|---|---|
| entry | 751 | 0.0 s | — |
| `ridx_build` | 751–900 | **738.0 s** | **no** |
| `kc_H_build` | 900–1090 | 166.7 s | **no** |
| `parallel_loop` | 1090–1214 | — | yes (the one pragma) |
| `filter_snv_emit` | 1214–1353 | — | **no** |
| `indel_pass` | 1353–2529 | historically 57% of total | **no** |

### 1.1 The `ridx_build` label is misleading and cost us the diagnosis

`ridx` is disabled by default (`WANT_RSUB = getenv("CAPS_FORCE_READSUB")`,
line 799). So the 738 s under that label contains **no ridx work at all** — it
is entirely the `build_substrate(seqs, cd_in)` call at line 756. The phase
marker sits after both, so the name attributes the time to the dead branch.
Rename before optimizing; a mislabelled profile is how a 12-minute serial loop
hides in plain sight.

### 1.2 `build_substrate` step 3 is embarrassingly parallel and runs twice

The 738 s is the read-placement loop (line 634, `for (size_t o = 0; o < n; ++o)`).
Work per read:

    2 strands x 8 seed offsets x <=64 index hits x 2 implied orientations
      x ~148-base mismatch scan

Iteration `o` writes **only** `S.read_cid[o]`, `S.read_pos[o]`, `S.read_rc[o]`,
`S.read_clip[o]` — disjoint per iteration. `idx` and `S.contigs` are read-only
once step 2 finishes. The only shared mutable state is the `placed` /
`improved` counters, which are a textbook reduction.

**It is called twice**: once at line 756 (pileup substrate, `dup=0.45`) and
again at line 1397 inside `indel_pass` (bubble substrate, `dup=0.92`). So this
single serial loop is both the largest measured phase *and* a large share of
the phase that has always been the caller's most expensive. Fixing it once
fixes both.

### 1.3 A prior negative result that does NOT apply here

The comment at line 716 records: *"OpenMP over contigs measured SLOWER at both
400kb and 5Mb scale (+11%, +13%)."* That was parallelising the **pileup over
contigs** — few items, wildly unequal sizes, poor load balance. Parallelising
**over reads** in `build_substrate` is a different loop shape: 12.6M uniform
items, `schedule(dynamic)` handles the tail. The earlier result is not
evidence against this one, and it is also not evidence *for* it. Measure it.

### 1.4 `kc_H_build`'s k-way merge is O(output x runs), and output is huge

The chunked k-mer counter (line ~945) is right in structure — bounding the
transient to one batch is exactly the fix the earlier reverted sorted-array
attempt needed. Its **merge**, however, carries a comment that was true at
window scale and is false at full scale:

> *"Runs are few (total k-mer occurrences / BATCH_KMERS -- roughly 47 at full
> chr20 scale), so a linear scan per output element is simple and correct; a
> heap would be faster but is not needed at this run count."*

The run count is indeed small. The **output** is not: distinct canonical
31-mers over 12.6M reads is on the order of 100–400M *(est., inflated by error
k-mers)*. The merge does **two** linear passes over all runs per output
element:

    ~400M outputs x 47 runs x 2 passes  =  ~37 BILLION comparisons

That is the 166.7 s. A `priority_queue` over run heads makes it
O(output x log runs) — roughly 400M x 5.5 — an **~8-17x reduction in
comparisons** for a mechanical change with no effect on the result (the merge
is a pure function of the runs). The comment's reasoning should be corrected
in place rather than deleted: it was sound for the scale it was written at.

Two further points in the same block:

- **`kc.push_back` has no `reserve`.** It grows to ~400M x 12 B ≈ 4.8 GB by
  doubling, so at the final reallocation it transiently holds old + new ≈
  **14.4 GB**.

  **The obvious fix is wrong and is recorded here so it is not tried again.**
  Reserving the summed run lengths *is* an exact upper bound, but it is a
  terrible one: each of ~47 runs carries its own copy of every common key, so
  the sum is bounded by 47 x 32M ≈ 1.5B entries ≈ **18 GB reserved up front** —
  worse than the doubling spike it was meant to remove. An upper bound is only
  useful as a reservation when it is *tight*, and this one is off by ~4x.

  The fix that does work is a **two-pass merge**: run the heap merge once
  counting distinct keys only (no `push_back`), `reserve` that exact count,
  then run it again to fill. That is 2x the merge cost — but since the heap
  already makes the merge ~8-17x cheaper, two heap passes still beat one
  linear-scan pass by ~4-8x *and* allocate exactly once, with no spike.
- **Batch collection and `std::sort(batch)` are serial.** Per-thread batches
  plus a parallel sort would thread the collect side too.

### 1.5 A gratuitous allocation, 25.2M times per call

    std::string r = strand ? rc_str(raw) : raw;      // line 645

Strand 0 **copies the read for no reason**. Two heap allocations per read x
12.6M reads x 2 calls = ~50M allocations that do nothing. Use a `const&` for
strand 0 and a reusable (thread-local, once parallel) buffer for the reverse
complement.

---

## 2. Where the RAM actually is

Sizes below use the run's own counts. Contig total ≈ 233,369 x ~335 bp
≈ 78 Mbp, so distinct contig 25-mers ≈ 70M (marked *est.* where it matters).

### 2.1 `rc_reads` — the largest structure, and most of it is never read

Line 2115, inside the pcluster channel, which is **ON by default**
(`if (!std::getenv("CAPS_NO_PCLUSTER"))`, line 2037):

```cpp
std::unordered_map<uint64_t, std::vector<std::pair<uint32_t,uint32_t>>> rc_reads;
for (uint32_t i = 0; i < seqs.size(); ++i)
    for (size_t j = LW; j + AK <= q.size(); ++j)      // AK=25, LW=40, stride 1
        rc_reads[cn].push_back({i, j});
```

With L=148, AK=25, LW=40 that is 84 positions per read:

    12,604,917 reads x 84 = 1.059 BILLION insertions

Distinct keys are genomic 25-mers **plus every error k-mer** (a single
sequencing error spawns up to 25 novel k-mers), plausibly 100–400M. Per key an
`unordered_map` node costs ~40 B + allocator overhead + a bucket slot + a
separate heap buffer for the vector — call it ~96 B before payload.

**Estimated 15–30 GB.**

**The decisive observation** — its only consumer is line 2141:

```cpp
for (auto& kv : pkidx) {                 // contig-unique anchors only
    auto rit = rc_reads.find(kv.first);  // <-- the ONLY lookup
    if (rit == rc_reads.end()) continue;
```

`rc_reads` is **only ever queried with keys that are in `pkidx`.** Every entry
whose key is absent from `pkidx` — which is the overwhelming majority, since
`pkidx` is filtered down to contig-*unique* anchors and error k-mers are never
contig k-mers at all — is built, hashed, stored, and never read once.

Filtering at insert time (`if (!pkidx.count(cn)) continue;`) is **byte-identical
by construction**: a loop that only reads keys in `pkidx` cannot observe the
absence of keys not in `pkidx`. This is not a heuristic or an accuracy
trade — it is deleting provably dead data.

Cost note: it requires `pkidx` to be built before `rc_reads`, which it already
is (line 2093 vs 2115).

### 2.2 `pkidx` / `pcount` / `porient` — three parallel maps over the same keys

Lines 2093–2099 build three separate `unordered_map`s keyed on the same ~70M
contig 25-mers, then erase from `pkidx` using the other two:

```cpp
std::unordered_map<uint64_t, std::pair<uint32_t,uint32_t>> pkidx;
std::unordered_map<uint64_t, uint32_t> pcount;
std::unordered_map<uint64_t, uint8_t>  porient;
```

Three node-based maps at ~70M keys each: **estimated 8–12 GB**, all alive
simultaneously, to compute a filter that a single sorted array would give for
free (a run of equal keys *is* its count, and the orientation rides along).

### 2.3 `idx` in `build_substrate` — ~6.7 GB, built serially

Line 620, `unordered_map<uint64_t, vector<pair<uint32_t,uint32_t>>>` over all
contig 25-mers, values capped at 64. ~70M keys x ~96 B ≈ **6.7 GB** *(est.)*,
and it is rebuilt on the second `build_substrate` call too.

### 2.4 Already fixed, kept for contrast

- `ridx` — eliminated (default off) after byte-identical verification on all
  5 windows. Was the largest single structure.
- `kidx` — already the flat sorted-array pattern, 16 B/entry, ~1.1 GB at 70M.
  **This is the template the structures above should follow.**
- `kc` — chunked counting; observed this run spiking to 25 GB and releasing
  back to 14.95 GB, i.e. behaving exactly as designed.
- `cov` / `pcov` / `lcov` — ~156 MB each. Not worth touching.

---

## 3. Ordered plan

Ranked by (value / risk), not by size. Nothing here is applied yet; the
full-chr20 run must finish first so there is a real before-number to gate
against — and per the standing project rule each change is measured
individually, not as a batch.

| # | change | expected | risk | validation gate |
|---|---|---|---|---|
| 1 | `#pragma omp parallel for schedule(dynamic,256) reduction(+:placed,improved)` on `build_substrate` step 3 | 738 s → ~90–120 s, **x2 sites** | LOW | byte-identical VCF (iterations are disjoint + deterministic) |
| 2 | drop the strand-0 string copy; thread-local rc buffer | ~50M allocations removed | LOW | byte-identical VCF |
| 3 | filter `rc_reads` inserts by `pkidx` membership | **15–30 GB → ~6–10 GB** | LOW | byte-identical VCF (provably dead data) |
| 4 | `pkidx`/`pcount`/`porient` → one flat sorted array | ~10 GB → ~1.6 GB | MED | F1-neutral on 5 windows |
| 5 | `idx` → flat sorted array (the `kidx` pattern) | 6.7 GB → ~1.1 GB, x2 sites | MED | F1-neutral on 5 windows |
| 6 | parallel fill + parallel sort for the flat arrays | serial 70M inserts → threaded | LOW after 4–5 | byte-identical vs its own serial version |
| 7 | `kc` k-way merge: linear scan → `priority_queue` | ~37B → ~2.2B comparisons (166.7 s phase) | LOW | byte-identical (merge is a pure function of the runs) |
| 8 | two-pass heap merge: count distinct, `reserve` exactly, then fill | removes a ~14.4 GB doubling spike; still ~4–8x faster than today | LOW | byte-identical |

**Why 4 and 5 are MED, not LOW.** Converting `kidx` from `unordered_map` to a
sorted array already taught this lesson once: `unordered_map` iteration order
is arbitrary, and several downstream structures are *first-write-wins*, so the
converted version produced a **different but not worse** call set. It had to be
validated on 5-window mean F1 rather than byte-identity. Expect the same here,
and do not report a byte-identity gate that cannot hold.

**Projected combined effect:** peak RAM ~32 GB → **~12–15 GB**; caller wall
time 40–70 min → **~15–20 min**. Both are projections from structure sizes and
core count, and both are to be replaced by measurements. They are written down
now precisely so that a later measurement can contradict them on the record.

---

## 4. What this does NOT fix

DiscoSNP++ runs comparable jobs in ≤6 GB because GATB partitions k-mer
counting **to disk** with minimizer-based buckets. Even at ~12–15 GB we would
still be 2–2.5x its footprint. Closing that gap is an architectural change
(disk-partitioned counting), not a tuning exercise, and nothing above attempts
it. If the full-scale run confirms we lose the RAM axis to DiscoSNP++, that
belongs in the paper as a stated limitation next to the accuracy result — the
162x compression number is a different measurement and does not offset it.
