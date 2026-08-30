# Layer-by-layer analysis: ours vs PgRC2, real profiling, no changes

Pure analysis, no code touched. Real numbers from real runs, both sides, same
file (E. coli, `ecoli_test.fq`, 1,000,000 reads, 150bp, byte-identical
verified on our side — see `LOCKED_SEQORDER_SCOPE.md` for that proof).

## 1. Our pipeline, per-stage real profiling

### 1a. Assembly (encoder) — internal phase timing, from real `/usr/bin/time -v`

| Phase | Time | Cumulative RSS |
|---|---|---|
| load+filter+dedup | 0.67s | 49 MB |
| prefix seed index | 0.06s | 64 MB |
| round 1 (division) | 0.49s | 92 MB |
| round 2 (assembly) | 0.18s | 111 MB |
| emit chains | 0.23s | 121 MB |
| pigeonhole mapping | 0.47s | 151 MB |
| pg MEM matching | 0.58s | 135 MB |
| **Total assembly** | **2.68s** | **peak 166 MB** |

Real finding: assembly's own peak RAM (166 MB) is NOT where the pipeline's
total peak RAM comes from — see 1b.

### 1b. Per-coder layer timing/RAM (isolated, real `/usr/bin/time -v` per tool)

| Layer | Tool | Time | Peak RSS | Output size |
|---|---|---|---|---|
| L1 (sequence) | seqpar | 0.22s | **820 MB** | 1,209,399 B |
| L3 (MEM refs) | refcoder | <0.01s | 4.5 MB | 137,065 B |
| L4 (pos+strand) | xz -9 (×2) | 0.70s | 96 MB | 526,888 B |
| L5 (mismatches) | mmcoder | 0.11s | 11 MB | 727,297 B |
| L6a (N-reads seq) | seqpar | 0.07s | 71 MB | 23,478 B |
| L8 (orig2uid) | xz -9 | 0.69s | 104 MB | 890,332 B |

**Major, previously-unstated finding: `seqpar` (L1) alone uses 820 MB RSS —
essentially the ENTIRE pipeline's reported peak RAM (819,572 KB from the
end-to-end wrapper run).** Assembly (166 MB) is not the RAM bottleneck this
project's own documentation has repeatedly implied — the sequence *coder*
is. This reframes the whole "assembly costs RAM" narrative that's been
assumed all session; it was never actually verified layer-by-layer until now.

Real question this opens, not yet answered: WHY does seqpar need 820 MB to
encode 5 MB of literal text? That ratio (~164x) is large enough to be a real
target — worth checking its context-table sizing (same class of question
already root-caused once this session for the quality coder's alphabet-size
issue) before assuming it's inherent.

**Checked for generalization immediately, not left as a single-dataset
guess (per the standing rule):** re-ran seqpar on P. aeruginosa's literal.txt
(7,012,675 B, ~1.4x larger than E. coli's 5,039,422 B). RAM: 821,536 KB —
essentially IDENTICAL to E. coli's 819,612 KB despite the different input
size. **This confirms the 820 MB is a FIXED allocation, not proportional to
input** — almost certainly a fixed-size context table sized for the DNA
alphabet's full context space (4 bases × some fixed context depth), not
scaled to the actual data volume. This is the exact same bug shape as the
quality coder's pre-fix behavior (fixed, oversized context table regardless
of alphabet/data size, found and fixed earlier this session via an
alphabet-adaptive switch) — real, strong reason to believe the same fix
pattern applies here, though NOT attempted or confirmed this pass (analysis
only, per instruction).

## 2. PgRC2's real per-stream breakdown (from ITS OWN stderr, same file)

Extracted directly from a live compress run — not estimated, not the paper's
numbers, this session's own measurement on the identical E. coli file.

| PgRC2 stream | Real reported cost | Time |
|---|---|---|
| **Order (permutation)** | 4,000,000 → **2,915,270 B** (LZMA) | 221ms |
| Sequence mapping compound (var-len + LZMA) | 6,323,005 → 1,312,634 B | 135ms |
| Good/Bad/N sequence mapping offsets+lengths | ~72,858 B combined | ~11ms |
| Mismatched symbols (PPMd ord=5) | 3,744,373 → 740,304 B | 50ms |
| Mismatch counts (zero flags + values) | ~214,129 B combined | ~3ms |
| ~50 small division-related range/PPMd streams (period 1-50) | sums to roughly 600-900 KB (not fully reconciled — see open question below) | <10ms each |
| Reverse-complement info | 981,389 → 23,524 B | 1ms |

**Major finding: PgRC2's OWN order/permutation stream is its single most
expensive stream — 2,915,270 B, bigger than its entire sequence stream
(1,312,634 B).** This directly contradicts an assumption this session
implicitly carried ("PgRC2's order is cheap, ours is the expensive part") —
it isn't. Order/read-list-restoration is expensive for BOTH tools; PgRC2
just pays it as one big LZMA'd raw permutation, while we pay it as two
smaller, more targeted arrays (L4 position 526,888 B + L8 orig2uid 890,332 B
= 1,417,220 B combined) — **a real, already-confirmed advantage on this
specific axis** (1.4M vs 2.9M), not something needing further work.

**Open question, not resolved in this pass:** the ~50 "range_coder/ppmd,
period=1..50" lines in PgRC2's log don't map cleanly to a single named
stream in the stderr output — they're almost certainly per-read-length-group
division/quality-adjacent bookkeeping (their tool interleaves quality-based
read division logic even in sequence-only SE mode, per `prepareChainData()`
read earlier this session), but confirming exactly what they encode would
need reading `runQualityBasedDivision()`'s output-writing code, not done
this pass — flagged, not guessed at.

## 3. Direct line-by-line code correspondence (both sides read this session)

| Our code | PgRC2 equivalent | What's the same | What's different |
|---|---|---|---|
| `100_locked_seqorder.cpp:1043-1074` (position/strand/length dump) | `SeparatedPseudoGenomePersistence.cpp:446` (`compressReadsPgPositions`, `singleFileMode` branch) | Both store position per read, absolute (not delta, in our current reverted version) | We index once per UNIQUE read + correlate via `orig2uid`; they store once per ORIGINAL read directly (no unique-read dedup for position at all) |
| `100_locked_seqorder.cpp:1377-1405` (trim-on-overlap match application) | `SimplePgMatcher.cpp:99-131` | Structurally identical algorithm now (sort by dst, trim overlap, mark-based literal emission) — confirmed by code, not just claimed, this session | Their marks are literal characters embedded in a string; ours are implicit via sorted triples + a separately-dumped literal.txt |
| `37_ref_coder.cpp` (range-coded src position) | Their sequence-mapping offset streams (LZMA, `Good/Bad sequence mapping - offsets`) | Both code a src/offset value per match | We use a range coder bounded by `log2(pg_len)`; they use generic LZMA over raw offsets — ours is a tighter, purpose-built floor (confirmed: `bounded floor 137,058 B` vs `this coder 137,065 B`, 99.995% of floor) |
| `50_mismatch_coder_real.cpp` (mmcoder) | `Mismatched symbols codes` (PPMd ord=5) | Both code substituted bases at mismatch positions | Untested this pass whether our context (ref+pos+prev-base) beats their PPMd-5 pound-for-pound — real open comparison, not done here |

## 3b. Follow-up: 4-dataset regression profiling (E. coli, yeast win; P.
aeruginosa, P. falciparum lose) — real per-layer fractions, one real
reversal found

| Dataset | L1(seq) | L4(pos) | L5(mm) | L8(orig2uid) | Result |
|---|---|---|---|---|---|
| E. coli (WIN) | 20.9% | 48.3% | 12.6% | 15.4% | +5.8% |
| yeast (WIN) | 42.6% | 41.9% | 3.5% | 9.0% | +2.0% |
| P. aeruginosa, small slice (LOSE) | 83.2% | 13.0% | 1.1% | 1.4% | -5.9% |
| P. falciparum (LOSE) | 43.4% | 27.2% | 6.1% | 20.6% | -2.7% |

Initial hypothesis from this table: losses correlate with L1 (sequence)
dominance. **Checked directly with real bits/base across all 4 — refuted.**
Sequence-coder efficiency is nearly identical everywhere (1.72-1.94
bits/base, P. aeruginosa is not a real outlier: 1.9425 vs E. coli's 1.9200).
The coder itself is not the differentiator.

**Re-tested P. aeruginosa at proper coverage** (the same coverage-mismatch
fix that flipped it to a win vs SPRING/Genozip earlier this session,
`pa_test_big.fq`, 1.4M reads) — **real, important finding: this does NOT
flip it against PgRC2.** Still loses, slightly worse (-6.1% vs -5.9% on the
small slice). Coverage-mismatch was a real bug for the SPRING/Genozip
comparison specifically; it does not explain or fix the PgRC2 comparison —
two genuinely different mechanisms, confirmed by testing rather than
assumed to be the same bug.

At the larger, coverage-matched scale, P. aeruginosa's layer-fraction
profile actually SHIFTS to look like a winning dataset (L1 29.1%, L4 51.5%
— now L4-dominant like E. coli/yeast) **yet still loses overall.** This
directly refutes the fraction-based hypothesis as causal — it was likely
confounded by dataset size (small slices amortize fixed-overhead layers
less), not a real signal of where the loss lives.

**The real, precise mechanism, found by comparing PgRC2's own real stderr
output for this exact file:** their order/permutation stream —
`5,600,000 → 4,188,183 B` (raw LZMA over the permutation, no dedup/rank
scheme at all) — is CHEAPER than our combined position+orig2uid cost on
this file (L4 4,393,572 + L8 715,096 = 5,108,668 B). **This is the opposite
of what E. coli showed** (there, ours beat theirs: 1.4M vs 2.9M). Real
generalization lesson: our position+correlation-array approach is not
universally better than PgRC2's simple raw-permutation+LZMA — it depends on
the file's actual duplication/order structure, and P. aeruginosa's
structure favors their simpler approach. This is exactly the kind of
single-dataset overfit risk the standing "algorithmic, not per-dataset"
rule exists to catch — found here specifically because a second dataset
was checked, not assumed from the first.

**Not yet resolved:** why P. aeruginosa's permutation structure favors raw
LZMA over our scheme specifically — real next question, not answered this
pass (analysis only, per instruction).

## 4. What this pass did NOT do (explicitly, per instruction — analysis only)

- Did not modify any code.
- Did not fix the 820 MB seqpar RAM question — flagged as the single
  highest-value real lead found this pass, not yet investigated further.
- Did not fully decode PgRC2's ~50 period-N streams.
- Did not run this same profiling on the other 5 datasets (P. aeruginosa,
  S. aureus, L. major, P. falciparum, yeast) — this is one dataset, and
  per the standing generalization rule, no conclusion here should be
  treated as proven until checked on more than one file.
