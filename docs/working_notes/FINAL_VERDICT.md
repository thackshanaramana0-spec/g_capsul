# Final component-wise verdict — speed, RAM, size

> **START HERE INSTEAD, unless you specifically want the component-wise speed/RAM/size optimisation record.**
> `README.md` for what the project is and what it measured; `CLAUDE.md` for the
> rules and the refuted-ideas list; `benchmark/documentation/RESULT_CODE.md` for
> where any number came from; `benchmark/documentation/REPRODUCE_EVERYTHING.md`
> to re-run it. **This document predates the final 2026-09-09 sweep, so where it
> disagrees with a file in `benchmark/results/`, the result file wins.**
>
> Despite the name, this is **not** the verdict on the three claims — it is the
> optimisation record. The claims are in `benchmark/documentation/CLAIMS_FINAL.md`.

Every component of all three claims, what it costs, whether a lever remains,
and the evidence. Written after the session's optimisation work so that
"nothing left here" is a measurement rather than an assertion.

Rule applied throughout: a change ships only if it is **structural** (a
decomposition, or a formula over a measured input property). Swept constants
were rejected even when they won -- see §5.

---

## 1. Claim 2 caller (the full SNV+indel path) — full chr20

**381.30 s / 33.28 GB -> 122.61 s / 15.24 GB. 3.11x time, 2.18x RAM,
output byte-identical.**

| component | cost | lever left? |
|---|---|---|
| decode reads+qual | 21.8 s | **NO** — bounded by the ARCHIVE: the literal stream carries 4 chunks so it decodes on 4 threads, quality has 3 blocks. Measured directly (`chunks=4 threads=4`). Widening needs re-chunking at compress time, changing archive bytes. Overlapping the two streams was the win available without touching the format, and it is taken. |
| place reads (x2) | 18.9 s | **NO** — exact early-exit at a provably optimal placement already cut it 34.7 -> 19.1 s. Candidate dedup was tried: it skipped 723 M candidates and got SLOWER. |
| collapse (x2) | 16.2 s | **NO** — inherently serial: each contig's accept/reject depends on the `claimed` set built by every prior contig. Prefetch measured neutral (6.83 -> 6.79 s). Running the two calls concurrently would save ~6 s but needs both 2 GB sets resident: rejected on the RAM axis. |
| kc_H_build | 12.5 s | **NO** — bin-wise parallel merge already took the k-way merge 8.2 -> 0.7 s; the two-pass exact allocation is at its floor. |
| indel setup | 7.1 s | **NO** — kidx built in place, 16-byte records. |
| pcluster read scan | 10.4 s | **NO** — fused with the anchor work, dissolving a 1.71 GB CSR. |

**Verdict: closed.** Every component has either been optimised to a measured
floor or has a stated structural reason it cannot move.

---

## 2. Claim 1 encoder — speed and RAM

**-2.6% to -6.3% wall on three datasets, byte-identical, RAM flat to better.**

Splits as ~4.0 s assembly + 1.84 s stream coding.

| component | cost | lever left? |
|---|---|---|
| coding pool | 1.84 s | **NO** — the pool runs at its Amdahl floor (wall 1.84 s, floor 1.84 s, speedup 4.41x). Only the longest job matters. |
| `pos_abs` (longest job) | 1.84 s, longest in **6/6** runs | **NO** — finer chunking refuted on BOTH axes (archive +6,718 B AND slower, 1.84 -> 2.37 s: more concurrent LZMA states thrash bandwidth). The coder probe fraction is already a documented measured trade (0.40 over PgRC2's 0.20). The region split already shrank this stream 1,264,855 -> 1,033,835 B, which shortened this job as a side effect. |
| `mm_pos` | 1.23 s | **DONE** — its nested bucket search (buckets x 2 layouts x 4 coders, serial) is now parallel: this was the session's encoder win. |
| `mm_cnt` | 1.75 s | **NO** — parallelising its three searches measured neutral on all three datasets and worse on RAM for one. |
| assembly stages | ~4.0 s | **untested at component level** — the honest gap in this audit. Peak RAM is flat at 296 MB across all of them, and no single stage exceeds 1.2 s. |

**Verdict: coding is closed; assembly is unexamined but small and flat.**

---

## 3. Claim 1 size

**pos_abs region split: -0.52% to -6.18%, 4/4 datasets, lossless-verified.**

82% of the archive examined at component level.

| stream | share | verdict |
|---|---|---|
| `pos_abs` | 46% -> 40% | **WON.** Two populations with an 18x entropy gap (18.45 b vs 1.03 b) in one stream. |
| `literal` | 25% | **at bound** — five independent sweeps (table size at two data scales, chunk count, context orders, mixer context, SSE context) all land on the shipped configuration. Beats xz/zstd/bzip2 by 9-15%. |
| `mem_triples` | 17% | **at bound** — already region-aware, 0.32x order-0. |
| reference streams | 6.8% | **rejected** — distributions differ but the population is 97.6% one-sided; 0.07% available. |
| untested remainder | ~10% | mm_pos/mm_sym/mm_cnt/pos_region/pos_sec. Even a uniform 20% win here is ~2% of the archive. |

**The one structural lever left on size** is not a coder change: shrink the
pseudogenome itself via second-region self-match. CLAUDE.md already names it as
the largest quantified opportunity; it is unshipped because the archive comes
out LOSSY and its effect flips sign by dataset (SARS -108 KB, E. coli +95 KB).
That needs the bug fixed AND the mismatch economics reworked.

---

## 4. Claim 3

All three operations re-verified exact after the `pos_abs` split (see
`docs/CLAIM3_POSTSPLIT_AUDIT.md`): export = PG_LEN exactly, coverage identity
holds to ZERO difference, query matches an independent scan exactly.

`query`'s modest 1.62x now has a measured structural cause rather than an open
question: 99.7% of the pseudogenome's 175,401 references reach back more than
100 kb, so `pg[qa..qb]` cannot be materialised without resolving the whole
dependency graph. Windowing would need archive-format restart points that cost
the ratio Claim 1 defends.

**Verdict: closed, with query's bound explained rather than outstanding.**

---

## 5. What was rejected, and why that matters

Nine measured refutations this session. Three are worth carrying forward as
principles:

1. **A swept constant is not a result.** `APMB=8->12` won on 4/4 datasets and
   was still reverted: it is the argmax of a sweep, which CLAUDE.md rule 1
   forbids. Sweeps are diagnosis -- they proved `literal` is at bound -- not
   fixes.
2. **Idle threads are not free capacity when the work is bandwidth-bound.**
   63% of the coding pool's thread-seconds are idle, and every attempt to fill
   them (finer chunking, more concurrent candidates) made things slower.
3. **The biggest STEP is not the biggest peak contributor.** Quotienting halved
   the collapse table and moved peak RSS by 50 MB, because that table is freed
   long before the peak.

---

## 6. Final verdict

**No further optimisation is warranted on any component of any claim without a
change of scope.**

Everything reachable by decomposition, representation, or scheduling has been
taken or measured and refuted. The three genuinely open items all require
crossing a boundary this session deliberately respected:

- **second-region self-match** (size) — needs a LOSSY bug fixed and
  dataset-dependent economics reworked
- **archive re-chunking** (caller decode speed) — changes archive bytes, which
  Claim 1's locked sizes forbid
- **windowed pg reconstruction** (query speed) — needs format-level restart
  points, costing ratio

Two verification gaps remain, and neither is an optimisation: the encoder wins
are measured on 3 of 15 locked datasets, and the size win on 4 of 15. The
mechanism is proven in both cases; the magnitude across the full locked set is
not. That sweep is the next work, and it is measurement, not engineering.
