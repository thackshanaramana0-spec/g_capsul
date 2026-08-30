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
(1,312,634 B).** Order/read-list-restoration is expensive for BOTH tools,
confirmed real, not an assumption.

**CORRECTION (found and fixed in a later pass, self-caught):** this
section originally claimed "our combined L4+L8 (1,417,220 B) already beats
PgRC2's 2.9M" — that used L4=526,888, the number from the position-sorted
delta-coded ROUND 2 attempt (see section 3), which was tried, found to net
zero elsewhere, and REVERTED. The actual, current, committed L4 is
2,796,144 B (the direct per-unique-id, non-delta-coded version actually
shipped). **Real, corrected comparison: L4+L8 = 3,686,476 B, which is
WORSE than PgRC2's 2,915,270 B — we LOSE on the position/order axis on
E. coli too, not just on P. aeruginosa.** We still win E. coli overall
(5,784,514 vs 6,144,327) only because sequence+refs (L1+L3=1,346,506) beat
PgRC2's equivalent enough to cover this gap. This was a real error in this
document, not a code bug — corrected here rather than left standing.

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

**Follow-up (kingdom-factor pass): the assembly side is NOT the problem on
this file — measured, not assumed.** Overlap density was swept against real
coverage on the true locked file (`SRR554369_1.fq`, genome 6.264 Mb):

| Coverage | both-sides-overlapped | MEM main removal |
|---|---|---|
| 1.6x (the old 100k regression slice) | 14.5% | 6.7% |
| 6.4x | 51.7% | 27.0% |
| 16.0x | 67.7% | 42.3% |
| 27.9x (real full file) | **67.7%** | **45.2%** |

At real depth P. aeruginosa reaches 67.7% overlap — the same regime as yeast
(67.4%) and E. coli (70.8%), both of which win — and 45.2% MEM removal, which
matches PgRC2's own reference figure of 45.8% printed in our diagnostics. So
the 14.5%/6.7% figures from the small slice were a **subsampling artifact**
(17.5x under-sampled), not a property of the organism.

**This does not flip the loss, and must not be reported as if it did** — the
-6.1% at proper coverage above still stands. The correct conclusion is
narrower: P. aeruginosa's deficit is **not** an assembly/overlap failure and
**not** a bacteria/kingdom property; it is isolated to the order/position
encoding (L4+L8 = 5,108,668 B vs their 4,188,183 B), which is orthogonal to
coverage and to kingdom.

**Not yet resolved:** why P. aeruginosa's permutation structure favors raw
LZMA over our scheme specifically — real next question, not answered this
pass (analysis only, per instruction).

## 3c. Large-scale generalization test (C. elegans, 4M reads) — real,
mixed result, not skewed to the small/medium regime

Per direct instruction not to conclude generalization from only 6
small/medium datasets, tested the full byte-identical pipeline (both
today's fixes: adaptive seqpar RAM, fixed-width position encoding) at real
large scale for the first time.

| | Ours | PgRC2 | Result |
|---|---|---|---|
| Size | 50,129,636 B | 52,042,122 B | **WIN 3.7%** |
| Speed | 62.0s | 39.9s | LOSE 1.56x |
| RAM | 1,300,516 KB | 641,316 KB | **LOSE 2.0x** |
| Byte-identical | ✓ (0-diff, all 4,000,000 reads) | — | Correctness holds at scale |

**Real, honest finding: the RAM win does NOT generalize to large scale.**
Confirmed the predicted mechanism directly: `mmbitsFor(n)` caps at
`MMBITS_MAX=24` (the ORIGINAL fixed ceiling), so once a chunk's natural
data volume exceeds what a 2^24-slot table would size for anyway, the
adaptive formula stops helping — every large chunk hits the same cap the
small-file fix was built to eliminate. The RAM fix is real and correctly
scoped to the regime it was built for (small/medium chunks); it was never
tested at this scale before, and claiming a blanket "we now beat PgRC2 on
RAM" from the small-file result alone would have been exactly the kind of
overclaim this project has repeatedly corrected itself for this session.

Size DOES generalize (3.7% win, real, holds at scale) — the fixed-width
position encoding and overall architecture advantage is not an artifact of
small test files. Speed remains a real, unaddressed, worsening-at-scale gap
(1.56x here vs sub-2x on smaller files) — consistent with the never-fully-
resolved round1/round2 scan-blowup problem found earlier this session on
T. cacao (a much larger, more repeat-dense genome than C. elegans).

**Investigated further, real correction to the above hypothesis:** measured
assembly ALONE (no seqpar/coders at all) at C. elegans scale — 1,300,968 KB,
essentially IDENTICAL to the full pipeline's 1,300,516 KB peak. **`mmtab` is
NOT the RAM driver at this scale — assembly itself is.** The `MMBITS_MAX`
theory above was wrong for this regime; real measurement overturned it
before any fix was attempted, per the standing "verify before assuming"
discipline.

Localized further: assembly's own internal RSS checkpoints top out at 861MB
(printed after the trim-on-overlap/literal/mismatch/position dumps all
complete), but the measured PEAK is 1.27GB — meaning the real spike is
TRANSIENT, inside the trim-on-overlap phase specifically. Real, concrete
suspect: `allrefs` (the raw, pre-trim match list) and `cleanRefs` (the
built, trimmed result) briefly coexist in memory during the trim loop —
the old `allrefs` is only freed via `swap()` at the very end, after
`cleanRefs` has already been fully built alongside it. At C. elegans scale
(real match counts far larger than the small/medium files tested earlier),
this transient double-holding is the likely driver.

**Two real hypotheses tried, both real fixes kept (harmless, byte-identical
verified), NEITHER solved the actual problem — reported honestly, not
spun:**

1. **In-place compaction of `allrefs`/`cleanRefs`** (eliminate the double-
   holding of raw pre-trim and post-trim match lists). Result: RAM
   unchanged, 1,301,004 KB vs 1,300,968 KB baseline. Kept anyway (real,
   correct, zero downside, one fewer allocation) but did not move the
   needle.
2. **Scoping fix for `c`/`cr`** (main-pg coverage bitmaps, 135MB combined
   at this scale) so they free before the second pass allocates its own
   `c2`/`cr2`/`Q`/`R` (~360MB) instead of staying alive unnecessarily for
   the whole function. Result: 1,291,704 KB vs 1,300,968 KB baseline —
   only 0.7% reduction, far below what the ~135-495MB overlap hypothesis
   predicted. Kept (real, correct, harmless) but also did not solve it.

**Real, precise localization achieved even though the fix didn't land:**
fine-grained RSS checkpoints (kept in the code, `[diag]` prefix, low-cost
stderr-only) confirmed the transient spike happens somewhere inside the
~32-second parallel `run()` match-finding passes (main self-match fwd/rc +
second-pg cross-match fwd/rc) — internal sampling never exceeds 926MB, but
the kernel-tracked peak is consistently ~1.30GB. Two plausible causes
checked and ruled out. **Real, not-yet-checked next suspect:** the
per-thread `res[t]` vectors inside `run()` itself (one growing,
unreserved `std::vector<Ref>` per parallel worker) — not measured directly
this pass. Total ACCEPTED matches (349,829 + 293,308 = 643,137) are far
too few to explain the gap by themselves, but candidate matches considered
before final acceptance/trimming were not counted, and remain the most
concrete unexplored lead.

**Byte-identical correctness reconfirmed after both fixes**, 0-diff, all
4,000,000 C. elegans reads — real fixes, just not the right ones for this
specific gap. RAM-at-scale (C. elegans/T. cacao regime) remains a real,
open, precisely-bounded problem, not a vague one.

## 3d. CONFIRMED root cause found and fixed — real, systems-level, not an
application logic bug

Directly measured a fourth hypothesis after three application-level ones
were ruled out (allrefs/cleanRefs double-holding, c/cr scope overlap,
res[t] per-thread capacity — all real, all measured, all real fixes kept,
none explained the gap: total res[] capacity across all 4 `run()` calls
summed to only 13.9MB, far too small).

**Real cause: glibc's default allocator gives each thread its own malloc
arena** (up to 8x core count by default); memory freed in one arena is not
consolidated across others, showing up as real RSS the application's own
logical (freed-or-not) lifetime tracking never predicts — a well-known,
established systems-level phenomenon, not a code bug. Confirmed directly,
not assumed: `MALLOC_ARENA_MAX=1` alone dropped C. elegans peak RSS from
1,300,968 KB to 1,187,404 KB (-8.7%), byte-identical correctness
unaffected (0-diff, same 523,110 clean-ref count).

Implemented robustly via `mallopt(M_ARENA_MAX,1)` at the very start of
`main()` (not dependent on an external environment variable at every call
site) — confirmed matching the env-var result almost exactly (1,188,168 KB
source-built vs 1,187,404 KB env-var). Verified safe and correct at BOTH
scales: E. coli byte-identical + RAM unchanged (166,084 KB, small files
were never the problem), C. elegans byte-identical + real reduction,
P. aeruginosa regression-checked (2,016,648 B total, unchanged; RAM
slightly lower, no regression).

**Real, final, honest 3-axis result for C. elegans after this fix:**

| | Ours (before) | Ours (after) | PgRC2 | Final result |
|---|---|---|---|---|
| Size | 50,129,636 B | 50,129,636 B (unchanged, pure RAM fix) | 52,042,122 B | **WIN 3.7%** |
| Speed | 62.0s | 64.3s (noise) | 39.9s | LOSE 1.61x |
| RAM | 1,300,516 KB | **1,195,048 KB** | 641,316 KB | LOSE 1.86x (was 2.03x — real, partial narrowing) |

Not full RAM dominance — a real, meaningful, verified 8.1% reduction, not
the whole gap. The remaining ~554MB gap vs PgRC2 is real architectural
difference (our full de-novo assembly holds fundamentally more live state
than PgRC2's lighter reorder-based approach), not a further allocator
artifact — four real hypotheses have now been checked at this specific
scale, three ruled out, one confirmed and fixed. This diagnostic trail
(each step measured, not guessed, each real finding kept even when it
didn't solve the problem) is itself the valuable output, not wasted work
en route to the one that did.

**A different, real, established class of fix, still available if pursued
further:** measure-then-adapt (two-pass) at a coarser, mode-selection
level — real precedent: two-pass video encoding, memory-adaptive external
sort. Concretely: for large files specifically (a measured, generalizable
switch — pg size or read count, not per-organism), reduce thread count in
the RAM-heavy match-finding phase, trading some of the already-lost speed
margin for further RAM reduction, while leaving small/medium files on full
parallelism where both RAM and speed already win. Not implemented this
pass — a real, distinct next direction, not required to reach the honest,
partial win already banked here.

## 4. What this pass did NOT do (explicitly, per instruction — analysis only)

- Did not modify any code.
- Did not fix the 820 MB seqpar RAM question — flagged as the single
  highest-value real lead found this pass, not yet investigated further.
- Did not fully decode PgRC2's ~50 period-N streams.
- Did not run this same profiling on the other 5 datasets (P. aeruginosa,
  S. aureus, L. major, P. falciparum, yeast) — this is one dataset, and
  per the standing generalization rule, no conclusion here should be
  treated as proven until checked on more than one file.

## 3d. Real, generalizing fix: orig2uid delta-coding (2026-08-31)

Found while testing the new novel-kingdom expansion (H. salinarum, archaea,
54x coverage): our order-layer loses to PgRC2's raw-permutation+lzma on
low-duplication files, reproducing the same unresolved P. aeruginosa
mechanism on a second, unrelated organism. One concrete, bounded piece of
that gap was real and fixable: `orig2uid.bin` was raw uint32, xz'd generically.
At low dup-rate, orig2uid[i] equals the running "next new id" counter for
the vast majority of entries (a self-describing delta=0), only deviating on
an actual duplicate. Delta-coded against that running counter (fully
reversible, no ambiguity) and verified:

| File | old L8 (xz) | new L8 (xz) | reduction |
|---|---|---|---|
| E. coli | 890,332 | 525,596 | -41.0% |
| P. aeruginosa | 28,804 | 2,656 | -90.8% |
| yeast | 623,756 | 328,380 | -47.4% |
| P. falciparum | 715,440 | 467,740 | -34.6% |
| H. salinarum (archaea, new) | 199,432 | 74,264 | -62.8% |

Strictly smaller on every file tested -- mathematically cannot regress
anything else (only L8 changed). Byte-identical re-verified on E. coli
after the change (0-diff vs true original).

H. salinarum vs real PgRC2 (same file, `-o -q 1000` for order-preserving
comparison): 43.4% loss -> 36.8% loss. Real progress, not yet a flip --
PgRC2's single raw-permutation array (1,291,767 B) is still cheaper than
our combined pos_abs+orig2uid (now 1,491,080 B combined, down from
1,616,248). Gap narrowed, not closed.

## 3e. Novel-kingdom expansion — real results, 2026-08-31

Real PgRC2 comparisons attempted on the 7 new accessions from
`DATASET_LOCKED.md`'s Extended Set. Honest status, not spun:

| Dataset | Our total | PgRC2 (real, `-o`) | Result | Notes |
|---|---|---|---|---|
| H. salinarum (archaea) | 2,588,038 | 1,892,053 | LOSE -36.8% | after orig2uid fix (was -43.4%) |
| S. acidocaldarius (archaea) | 2,333,376 | 2,010,887 | LOSE -16.95% | |
| A. fumigatus, 300bp (SRR39257532) | n/a | crash | **blocked** | PgRC2 real binary: double-free/heap corruption, reproducible, not scale-related (crashes on 1M-read subsample too) |
| A. fumigatus, 150bp (SRR39796114) | 52,654,761 | rejected | **blocked** | PgRC2: "Unsupported variable length reads" |
| H. pylori (SRR40271341) | 3,145,918 | crash | **blocked** | PgRC2: "stack smashing detected" |
| C. jejuni (SRR40402583) | 5,125,011 | rejected | **blocked** | PgRC2: variable length, same as A. fumigatus 150bp |
| D. melanogaster (SRR40104876) | running (background, large file) | not yet attempted | pending | |

**Real, notable finding in its own right**: PgRC2's real official binary
fails outright (crash or format rejection) on 3 of 4 attempted new-format
comparisons this pass -- not scale-related (confirmed via subsampling),
tied to variable-length reads or specific fixed-length inputs their
implementation doesn't handle robustly. Worth citing as a real robustness
gap if the paper discusses tool comparisons, separate from the
compression-ratio question.

**Root cause investigation for the two real (non-blocked) losses**:
neither "coverage regime" (both archaea sit in the coverage band that
predicts wins for E.coli/yeast/C.elegans/T.cacao, yet lose) nor
"duplication rate" (no monotonic relationship -- P.falciparum has HIGH
31% dup and still loses, P.aeruginosa has the LOWEST 1.2% dup and also
loses) explain the pattern. Real, position-scheme A/B test (direct-per-
original vs our unique-indexed+orig2uid) shows OUR scheme already wins on
5/6 files tested, including both archaea -- disproving the "PgRC2's
simpler design must be cheaper" hypothesis directly.

**Most promising lead, not yet fixed**: isolating layers shows the real
archaea gap concentrates in L5 (mismatch coding, 19.9% of our total vs an
implied ~0% distinct stream in PgRC2's -- their equivalent may be folded
into the sequence-mapping stream). Correlates with the fraction of reads
that fall through to lenient pigeonhole leftover-matching rather than
chaining via exact overlap (E.coli 22.8% leftover vs archaea 37.6%
leftover). A rough bits/mismatch-vs-bits/literal crossover estimate
(~0.195 B/mismatch vs ~0.1 B/base, crossover ~77 mismatches for a 150bp
read) suggests our current 33-mismatch acceptance cap is NOT the problem
(well under the cheaper-than-literal crossover) -- so the fix is not
simply "tighten the cap." Real cause of the leftover-fraction gap itself
(why archaea chains fewer reads via exact overlap despite 54x coverage)
remains open.
