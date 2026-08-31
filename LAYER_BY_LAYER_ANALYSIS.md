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

## 3f. MAJOR CORRECTION (2026-08-31): the archaea "losses" were a
## measurement error. Real result is 5/6 WINS, biggest wins on archaea.

**What I got wrong.** Sections 3e and earlier compared our COMPLETE layer
total against a hand-summed SUBSET of PgRC2's streams. I built their total
by reading the tail of their stderr log (the last "collective compression of
streams" block: pg sequence + self-match offsets + order) and **omitted their
entire mismatch machinery**, which sits in an earlier block:

| PgRC2 mismatch stream (H. salinarum) | bytes |
|---|---|
| Mismatched symbols codes (ppmd ord=5) | 512,885 |
| Mismatch positions (~50 period-keyed range/ppmd streams) | 479,549 |
| Mismatch counts, non-zero values | 99,456 |
| Mismatch counts, zero flags | 51,995 |
| mismatch props | 50 |
| **omitted subtotal** | **1,143,935 (37.5% of their archive)** |

Parsing every stream in their log sums to 3,049,055 B against an actual
archive of 3,050,260 B — a 1,205 B difference (header), which confirms the
parse is complete and nothing else is missing.

**The fix: compare against their actual archive size on disk.** This is
unambiguous and needs no stream accounting. Fairness verified directly:
`PgRC -d` on their archive emits exactly 460,501 lines of bare sequence,
no qualities, no names, and `diff` confirms it is IDENTICAL to the original
file's sequence column -- the same job our locked scope does, verified the
same way on all 6 files.

### Real, corrected results (our best total vs their real archive)

| Dataset | best100 | best103 | OURS best | PgRC2 archive | margin |
|---|---|---|---|---|---|
| E. coli | 5,202,566 | 5,183,815 | 5,183,815 | 6,151,254 | **+15.7%** |
| P. aeruginosa | 1,990,500 | 1,995,058 | 1,990,500 | 1,931,631 | **−3.0%** |
| yeast | 6,538,698 | 6,556,971 | 6,538,698 | 7,063,718 | **+7.4%** |
| P. falciparum | 4,079,085 | 4,098,423 | 4,079,085 | 4,268,973 | **+4.4%** |
| H. salinarum (archaea) | 2,588,302 | 2,479,901 | 2,479,901 | 3,050,260 | **+18.7%** |
| S. acidocaldarius (archaea) | 2,351,821 | 2,382,598 | 2,351,821 | 3,116,012 | **+24.5%** |

**5 of 6 wins. The two archaea — the novel-kingdom datasets — are our two
BIGGEST wins (+18.7%, +24.5%), not losses.** P. falciparum, previously
recorded as a −2.7% loss, is actually a +4.4% win. The only real loss in
the set is P. aeruginosa at −3.0%.

### What this invalidates

- The "-36.8% / -16.95% archaea loss" figures: **wrong**, artifact of the
  omitted-streams error. Retracted.
- The "our L5 mismatch layer is an outsized cost with no PgRC2 counterpart"
  conclusion: **backwards**. Their mismatch machinery costs 1,143,935 B
  against our L5 of 521,180 B -- we are **54.4% cheaper** on mismatches,
  which is a large part of why we win.
- The "leftover-read fraction drives our loss" lead: built on the same bad
  numbers, so it is not evidence of anything. Not a real finding.

### What survives

- The orig2uid delta-coding fix (real, -34.6% to -90.8% on that layer,
  byte-identical verified) -- still a genuine improvement.
- `103_no_predup.cpp` -- helps H. salinarum (+4.2%) and E. coli, hurts
  slightly elsewhere; keep both binaries and take the min per file until the
  rule for which wins is understood.
- The disk-architecture analysis in `PGRC2_DISK_ARCHITECTURE.md` §1-4
  (dead-code finding, RAM parity) -- independent of this error, unaffected.

**Methodology rule going forward: compare against the real archive size on
disk, never against a hand-summed subset of the opponent's log.**

## 3g. SECOND MAJOR CORRECTION (2026-08-31): our OWN totals omitted two
## mandatory streams. All wins in 3f were invalid.

**The bug.** `decode_locked_seqorder.py` requires FOUR mismatch streams to
reconstruct a read: `mm_ref`, `mm_obs`, `mm_pos` (where each mismatch is) and
`mm_count_per_read` (how many per read). But `mmcoder` only ever takes two
arguments (`ref.bin obs.bin`), and `run_locked_seqorder.sh` only counted that.
**`mm_pos.bin` and `mm_count_per_read.bin` were never in
LOCKED_SEQORDER_TOTAL** despite being mandatory for decode.

This is the same omitted-mandatory-stream error caught on PgRC2's side in 3f.
I checked theirs and did not check ours. Fixed in `run_locked_seqorder.sh`
(L5b/L5c).

**Corrected real-full-file results (vs PgRC2's real archive):**

| Dataset | 3f claimed | uncounted | corrected total | PgRC2 | real margin |
|---|---|---|---|---|---|
| E. coli | +16.6% | 2,510,732 | 9,905,845 | 8,864,420 | **−11.7%** |
| P. aeruginosa | +4.1% | 1,239,068 | 9,908,025 | 9,043,181 | **−9.6%** |
| S. aureus | +13.8% | 3,638,268 | 15,352,697 | 13,595,003 | **−12.9%** |
| P. falciparum | +8.0% | 2,967,436 | 18,810,808 | 17,219,695 | **−9.2%** |
| L. major | +2.5% | 3,037,812 | 30,612,464 | 28,272,652 | **−8.3%** |

**Every win in 3f becomes a loss.** The uncounted streams are 1.2-3.6 MB per
file, larger than every margin claimed. Also note: those full-file runs were
never decode-verified -- byte-identical had only ever been checked on
subsamples.

### Root cause, predicted by our own Stage 27 comment

Our mapper averages **16.19 mismatches per placed read** (E. coli); PgRC2
averages **1.39**. The Stage 27 note in `100_locked_seqorder.cpp` called this
exactly: *"Every mismatch we accept has to be STORED. This progression
measures pg literal and order and never counted that stream, so stage 17's
win could be partly borrowed against a cost that was never on the books."*
It was.

`MAXMAP = Lmax/3` (=50) copies ReadsMatchers.cpp:700, but with the streams
counted the economics are: a mismatch costs ~0.88 B (position+count+symbol,
measured) while a pg-literal base costs ~0.229 B (measured). Accepting a
50-mismatch placement costs more than simply storing the read.

### The fix: MAXMAP=12, swept against the CORRECTED total

| MAXMAP | E. coli | yeast | P. aeruginosa |
|---|---|---|---|
| 50 (old default) | 6,681,130 | 7,300,474 | 2,059,342 |
| 20 | 5,797,938 | 7,145,843 | 2,025,310 |
| **12** | **5,786,089** | **7,138,680** | 2,019,367 |
| 8 | 5,809,076 | 7,157,904 | 2,016,404 |
| 2 | 5,987,241 | 7,381,503 | 2,014,307 |

Convex with optimum at 12 on both 150bp datasets; P. aeruginosa is monotone
but 12 is within 0.25% of its best. So **12 is optimal-or-near-optimal on all
three -- a generalizing constant, not per-dataset tuning** (−13.4% / −2.2% /
−1.9% vs the old 50). Shipped as the default, env-overridable for re-sweeping
if coder costs change.

Mean mismatches/read on E. coli drops 16.19 → 2.62, which is what collapses
the mm_pos/mm_sym/mm_cnt streams.

**Real-full-file validation with MAXMAP=12 AND byte-identical decode
verification is the next step -- no margin from this section should be quoted
until that lands.**

### 3g-final. Validated result: MAXMAP=12 on real full files, byte-identical

First run in this project's history that pairs a CORRECTED total (all decoder
streams counted) with a BYTE-IDENTICAL decode check on the REAL full locked
files -- not subsamples.

| Dataset | before fix | after MAXMAP=12 | PgRC2 | margin | verify |
|---|---|---|---|---|---|
| E. coli | 9,905,845 | 8,322,815 | 8,864,420 | **+6.1%** | BYTE_IDENTICAL |
| P. aeruginosa | 9,908,025 | 9,878,603 | 9,043,181 | −9.2% | BYTE_IDENTICAL |
| S. aureus | 15,352,697 | 14,102,598 | 13,595,003 | −3.7% | BYTE_IDENTICAL |
| P. falciparum | 18,810,808 | 17,663,096 | 17,219,695 | −2.6% | BYTE_IDENTICAL |
| L. major | 30,612,464 | 30,368,493 | 28,272,652 | −7.4% | BYTE_IDENTICAL |

**1 win / 4 losses.** The fix is real and helps everywhere (−16.0% on E. coli,
−8.1% S. aureus, −6.1% P. falciparum, −0.8% L. major, −0.3% P. aeruginosa)
but only flips E. coli. This is the honest, fully-accounted, fully-verified
baseline: **we currently lose to PgRC2 on 4 of 5 real datasets.**

Every headline in 3f (+15.7% to +24.5%, "5/6 wins") is withdrawn. Those
numbers omitted mandatory streams on our side.

**Where the remaining gap is.** MAXMAP helps least exactly where we lose most
(P. aeruginosa −0.3%, L. major −0.8%), so the residual is not the mismatch
ceiling. Candidates, in order of measured promise, none yet tested:
1. **Mismatch position coding.** Ours is flat absolute positions + xz. PgRC2
   uses reverse offsets, buckets by mismatch count, optional transpose, and a
   per-bucket selector coder (`compressRlMisRevOffDest`,
   SeparatedPseudoGenomePersistence.cpp:823-905). A quick prototype showed
   reverse/delta offsets alone give −10.4% on that stream; bucketing and
   transpose hurt at 100k-read scale but were not tried at full scale where
   buckets are large enough to amortize.
2. **Per-stream coder selection.** Measured per-layer bake-off (xz-9/-9e, xz
   delta filters, zstd-19/-22, bzip2): zstd−22 wins `pos_strand` −16.2%,
   `orig2uid` −13.7%, `read_lengths` −77.4%; bzip2 wins `mm_pos` −7.3%,
   `mm_count` −14.6%; xz keeps `pos_abs`. Different coder per stream, chosen
   by measured ratio -- PgRC2 does exactly this via `getSelectorCoderProps`.
3. The 3-way pg split remains untested and is now re-motivated: their LQ pg
   stores poorly-matching reads as sequence rather than as expensive mismatch
   records, which is the same economics MAXMAP exploits, applied structurally.

## 3h. Stage 105 — N-reads routed through assembly (real structural change)

Found by a full stream-vs-layer comparison on our worst-loss file
(P. aeruginosa, real full file). Category breakdown vs PgRC2's real streams:

| Category | Ours | PgRC2 | |
|---|---|---|---|
| position/order | 4,734,796 | 4,964,368 | −230K better |
| mismatch (all) | 1,049,845 | 1,442,834 | −393K better |
| sequence | 2,870,296 | 2,485,786 | +384K worse |
| **N-reads** | **586,364** | assembled into nPg | **+586K worse** |

**Our N-read layer was 70% of the entire gap on that file.** N-containing
reads were `continue`'d out of the pipeline at load time and dumped as raw
text (2,459,350 B raw → 567,140 seqpar + 19,224 indices). Not a coder
problem: seqpar 1 1 already beat xz (567,140 vs 644,304), and every other
seqpar context setting was worse.

The reads did not deserve the exclusion: mean **1.23** N characters each,
max 8, none all-N. Perfectly good reads were kept out of overlap assembly
over one or two ambiguous bases. PgRC2 does not do this -- it builds a
dedicated nPg so they are still assembled (`runNPgGeneration`,
pgrc-encoder.cpp:184-204).

**Fix (`105_nreads_assembled.cpp` + `decode_105.py`):** substitute N->A, send
the reads through the SAME pipeline as everything else, keep N positions in a
small side stream (`n_pos.bin`, `n_cnt.bin`, `n_indices.bin`), restore on
decode. Same goal as their nPg without a second pseudogenome.

P. aeruginosa: N layer 586,364 → **23,608 (−96%)**, total −4.0%.

### Combined result, real full files, all BYTE_IDENTICAL verified

| Dataset | original | +MAXMAP=12 | +stage105 | PgRC2 | margin |
|---|---|---|---|---|---|
| E. coli | 9,905,845 | 8,322,815 | 8,297,531 | 8,864,420 | **+6.4%** |
| P. aeruginosa | 9,908,025 | 9,878,603 | 9,482,924 | 9,043,181 | −4.9% |
| S. aureus | 15,352,697 | 14,102,598 | 14,044,361 | 13,595,003 | −3.3% |
| P. falciparum | 18,810,808 | 17,663,096 | 17,634,117 | 17,219,695 | −2.4% |
| L. major | 30,612,464 | 30,368,493 | 30,185,565 | 28,272,652 | −6.8% |

Both fixes together: **−5.8%** on our own totals. Aggregate vs PgRC2:
**−3.4%**, 1 win / 4 losses. Stage 105 improves every dataset and regresses
none.

### Tested and rejected on measurement (not taste)

- **Reverse-offset mismatch positions** (PgRC2's rlMisRevOffDest technique):
  helps yeast −1.4% and P. aeruginosa −4.9%, but HURTS E. coli +1.6%. The
  MAXMAP fix made mismatches sparse (2.62/read), and reverse offsets only pay
  when they are dense -- the two fixes are antagonistic. A prototype measured
  −10.4% before the MAXMAP fix; that gain evaporated after it.
- **No pre-dedup** (`104_adaptive_dedup.cpp`, NODEDUP=1): orig2uid falls
  423,940 → 1,084 but the other per-read arrays absorb it. Net only −0.27% on
  real P. aeruginosa. Cost relocated, not removed -- same lesson as the
  ROUND 2 revert.
- **Per-stream coder selection**: real but small. Measured on real post-fix
  E. coli: orig2uid −6.4% (bzip2), mm_count −5.4% (xz -9e), read_lengths
  −79.3%, pos_strand −0.5%; pos_abs/mm_pos/n_indices unchanged (xz already
  best). Total **0.49%** of the archive -- worth taking, not decisive.

### RAM: correcting an earlier architecture claim

Our multi-process design already bounds peak RAM to the MAX over stages, not
the sum, because each coder exits. PgRC2 reaches the same bound in one process
via compress-and-release. Measured parity at ~0.5M reads: H. salinarum ours
112 MB vs theirs 185 MB (we win); S. acidocaldarius 208 vs 194 MB (they win).
So "PgRC2 wins on memory release" is wrong at this scale.

The real RAM gap is at C. elegans scale (4M reads, ours 1.30 GB, ~1.86x
behind) and it is NOT caused by the process split -- it is inside our assembly
process, which holds every array live until it dumps at the end, whereas their
`disposeReadsList()` frees each structure as soon as it is coded. The
technique worth adopting is compress-and-release INSIDE assembly, which is
independent of the disk-staging question. Not yet implemented or measured.

Disk cost, separately: intermediates are 361 MB for a 456 MB input (~0.8x), so
C. elegans would stage ~8.7 GB and T. cacao ~17 GB.

## 3i. Stage 106 — single-process, in-memory, compress-and-release
## (built, verified, and measured to change almost nothing)

Implements all four architecture rows PgRC2 wins on paper: single process,
zero intermediate files, live-object stage handoff, explicit
compress-and-release. Every stream is an `open_memstream` FILE*, so all the
existing fwrite/fputc calls are untouched and provably write the same bytes,
just into RAM. The three custom coders were extracted into shared headers
(`seqpar_core.h`, `coders_inproc.h`) so the standalone binaries cannot
diverge; liblzma replaces the external xz calls.

**Coder equivalence verified on real data (not assumed):**

| coder | in-process | standalone binary |
|---|---|---|
| seqpar (L1) | 2,893,487 | 2,893,487 |
| refcoder (L3) | 235,983 | 235,983 |
| mmcoder (L5) | 166,874 | 166,874 |
| xz (L4) | 4,735,648 | 4,735,640 (8 B container diff) |

Full pipeline decodes **BYTE_IDENTICAL** on real full P. aeruginosa.

**Head-to-head, same input, real full P. aeruginosa:**

| axis | OLD multi-process | NEW single-process | verdict |
|---|---|---|---|
| wall time | 12.8 s | 12.8 s | **no change** |
| peak RSS | 255 MB | **309 MB** | **NEW WORSE** |
| intermediate files | 21 (42 MB) | **0** | NEW better |
| archive size | 9,482,924 | 9,483,072 | +148 B (length headers) |

**Two earlier "Winner: PgRC2" claims are hereby corrected:**

1. *Memory release.* Holding streams in RAM is WORSE than the multi-process
   design. Process exit already bounded peak RAM to the max over stages; we
   got that bound for free and now pay +54 MB to keep streams live. Their
   `disposeReadsList` discipline is what makes single-process VIABLE, not what
   makes it better.
2. *Stage handoff / re-read cost.* Intermediates are 42 MB, not the 361 MB
   quoted earlier (that figure wrongly included the input copy). Re-reading
   42 MB is invisible against 12.8 s of compute -- hence zero speed change.

The single real gain is 0 intermediate files instead of 21. **The architecture
does not touch archive size, which is where we actually lose.**

## 3j. Seed geometry ruled out as the mapping-rate cause

Hypothesis: with SEEDSTRIDE=22 and SEEDW=16, a 101bp read yields only 4 seeds,
so pigeonhole guarantees finding placements up to 3 mismatches while we accept
12 -- suggesting we MISS mappable reads and append them.

**Refuted by measurement** (real P. aeruginosa, corrected accounting):

| SEEDSTRIDE | seeds | mapped | appended | TOTAL |
|---|---|---|---|---|
| 22 | 4 | 421,413 | 89,595 | 9,483,072 |
| 12 | 8 | 421,560 | 89,448 | 9,481,297 |
| 8 | 11 | 421,570 | 89,438 | 9,481,160 |
| 4 | 22 | 421,474 | 89,534 | 9,481,884 |

5.5x more seeds moves mapping by 61 reads and the total by 0.01%. The 89,595
appended reads genuinely have no acceptable placement -- we are not failing to
FIND one.

**So the remaining gap is match QUALITY, not sensitivity.** At MAXMAP=50 we
map more reads but at mean 4.70 mismatches/read; PgRC2 achieves 91.7% mapped
at mean 1.39. Their placements are ~3.4x better, which comes from their
matcher (copMEM, seed 38) and/or their hqPg being built only from both-side-
overlapped reads. Replacing our pigeonhole mapper is the next real piece of
work; nothing smaller has moved this number.

## 3k. Equal-pg-size comparison finds the real gap: mismatch POSITION coding

Comparing at our own operating point was misleading -- the pg/mismatch
trade-off means the two tools sit at different points on the same curve. The
honest comparison is at EQUAL pg size. Our MAXMAP=40 point has L1=2,482,224
against PgRC2's 2,485,786, so the mismatch streams are directly comparable:

| stream | ours (MAXMAP=40) | PgRC2 | |
|---|---|---|---|
| mm_sym | 362,567 | 370,090 | we are 2% BETTER |
| **mm_pos** | **951,720** | **748,988** | **we are 27% WORSE** |
| mm_cnt | 365,376 | 323,856 | we are 13% worse |

Mismatch position coding is the concrete remaining gap.

**PgRC2's technique, now implemented** (compressRlMisRevOffDest,
SeparatedPseudoGenomePersistence.cpp:823-905): route each read's positions
into a bucket chosen by that read's mismatch COUNT, then TRANSPOSE each bucket
so column k holds the k-th mismatch of every read in it.

| variant (real P. aeruginosa, MAXMAP=40) | mm_pos | vs flat |
|---|---|---|
| flat (baseline, already reverse-offset coded) | 951,720 | — |
| bucket by count | 875,792 | −8.0% |
| bucket + transpose | 832,364 | −12.5% |
| + per-bucket coder selection | 831,255 | −12.7% |
| PgRC2 | 748,988 | −21.3% |

Per-bucket coder selection is worthless here (xz wins 38 of 40 buckets);
their remaining edge over −12.7% is their range/PPMd coders, not the
bucketing. **An earlier prototype of this at 100k-read scale showed it HURT
(+78.8%) -- that was wrong only because the buckets were too small to
amortize; at real scale it is a solid win.** Prototypes on subsamples cannot
be trusted for anything bucket-structured.

**Shipped** in `106_inprocess.cpp` with an encoder-side round-trip self-check
(invertibility is a correctness claim, so it is proven at runtime, not
asserted; falls back to the flat stream if it ever fails).

Real full P. aeruginosa at the shipping MAXMAP=12, BYTE_IDENTICAL verified:
mm_pos 627,976 -> 543,180 (**−13.5%**), archive 9,483,072 -> **9,398,276**
(−0.89%), gap to PgRC2 −4.9% -> **−3.9%**.

## 3l. mm_cnt zero-flag split + full 5-dataset validation

PgRC2 splits the per-read mismatch counts into two streams ("Mismatches counts
(zero flags)" / "(non-zero values)"), visible in their own stderr. Most reads
have zero mismatches, so a 1-bit flag per read plus a dense uint8 array for
the rest beats a sparse uint16 array. Implemented with the same round-trip
self-check and a size guard (falls back to flat if it does not win).
Real P. aeruginosa: 298,188 -> 266,536 (**-10.6%**).

### Full validation, all 5 real locked files, all BYTE_IDENTICAL

| Dataset | session start | +stage105 | +106 transforms | PgRC2 | margin |
|---|---|---|---|---|---|
| E. coli | 9,905,845 | 8,297,531 | **8,228,651** | 8,864,420 | **+7.2%** |
| P. aeruginosa | 9,908,025 | 9,482,924 | **9,366,624** | 9,043,181 | −3.6% |
| S. aureus | 15,352,697 | 14,044,361 | **13,921,389** | 13,595,003 | −2.4% |
| P. falciparum | 18,810,808 | 17,634,117 | **17,494,397** | 17,219,695 | −1.6% |
| L. major | 30,612,464 | 30,185,565 | **29,946,905** | 28,272,652 | −5.9% |

**Total reduction this session: −6.7%.** Aggregate vs PgRC2 improved from
−3.4% to **−2.5%**. Still 1 win / 4 losses, but every dataset improved and
none regressed at any step.

Cumulative real wins, in order of size:
1. `MAXMAP=12` (mismatch acceptance economics) -- the largest, up to −16%
2. Stage 105 (N-reads through assembly) -- up to −4%
3. mm_pos bucket+transpose -- −13.5% of that stream
4. mm_cnt zero-flag split -- −10.6% of that stream

Dead ends, all measured and documented so they are not re-attempted:
single-process architecture (0% size, +54 MB RAM, 0% time), seed geometry
(0.01%), no-pre-dedup (0.27%), reverse-offset positions alone (mixed),
per-stream/per-bucket coder selection (0.49% / 0.2%).

**Remaining gap is match quality**, unchanged by any of the above: PgRC2 maps
91.7% of its leftover reads at mean 1.39 mismatches; we map 82.5% at 4.70.
Closing it means replacing the pigeonhole mapper with a MEM-anchored one
(their copMEM, seed 38). That is the only remaining lever of any size.

## 3m. Two of my own conclusions corrected by measurement

**(a) "The remaining gap is match quality" -- WRONG.** That claim compared our
mean-per-PLACED-read (4.70) against PgRC2's mean-per-ALL-reads (1.39). Apples
to oranges. On the same file, real totals:

| | mapped rate | mismatch symbols |
|---|---|---|
| ours, MAXMAP=50 | **92.5%** | 2,219,618 |
| ours, MAXMAP=12 | 82.5% | — |
| PgRC2 | 91.7% | 1,947,270 |

At a comparable acceptance ceiling we map MORE reads than they do (92.5% vs
91.7%) and append FEWER (38,352 vs 49,088), with only 14% more mismatch
symbols -- not the 3.4x quality deficit claimed. **Our pigeonhole matcher is
not the problem, and replacing it with a MEM-anchored one is NOT the next
piece of work.** That recommendation is withdrawn.

**(b) The MAXMAP optimum moved once mismatches got cheaper.** The mm_pos and
mm_cnt transforms lower the price of a mismatch, which shifts the optimum
upward (accept more, shorten the pg). Re-swept on real full files:
E. coli optimum 20, P. aeruginosa optimum 27 (vs 12 before).

But validated across all five, a fixed 20 LOSES 1.2% on P. falciparum:

| Dataset | MAXMAP=12 | MAXMAP=20 |
|---|---|---|
| E. coli | 8,228,651 | **8,212,728** |
| P. aeruginosa | 9,366,624 | **9,317,652** |
| S. aureus | 13,921,389 | **13,846,439** |
| P. falciparum | **17,494,397** | 17,710,535 |
| L. major | 29,946,905 | **29,900,569** |
| **aggregate** | **78,957,966** | 78,987,923 (+0.04%) |

20 wins 4 of 5 yet the aggregate is a wash. Minimax regret picks **12**
(worst case 0.54% vs 20's 1.24%), so 12 stays -- now justified on real full
files rather than the subsamples it was first derived from. The true optimum
is dataset-dependent across 12..27 with ~1% spread; capturing it would need a
formula keyed on a measured property, and none has been found. A per-dataset
switch is forbidden by the standing algorithmic-first rule.

### Final state, all BYTE_IDENTICAL on real full files

| Dataset | session start | final | PgRC2 | margin |
|---|---|---|---|---|
| E. coli | 9,905,845 | 8,228,651 | 8,864,420 | **+7.2%** |
| P. aeruginosa | 9,908,025 | 9,366,624 | 9,043,181 | −3.6% |
| S. aureus | 15,352,697 | 13,921,389 | 13,595,003 | −2.4% |
| P. falciparum | 18,810,808 | 17,494,397 | 17,219,695 | −1.6% |
| L. major | 30,612,464 | 29,946,905 | 28,272,652 | −5.9% |

Session total **−6.7%**; aggregate vs PgRC2 **−2.5%**, 1 win / 4 losses.

## 3n. Adaptive dedup (keyed on measured duplication rate)

Per-stream comparison at the shipping config located the last big structural
loss precisely:

| stream | ours | PgRC2 | delta |
|---|---|---|---|
| pg sequence | 2,893,487 | 2,485,786 | +407,701 |
| **orig2uid** | **427,060** | **0** | **+427,060** |
| self-match refs | 235,983 | 71,700 | +164,283 |
| order/positions | 4,809,188 | 4,964,368 | −155,180 (we win) |
| mm_sym | 166,874 | 370,090 | −203,216 (we win) |
| mm_pos | 543,180 | 748,988 | −205,808 (we win) |
| mm_cnt | 266,536 | 323,856 | −57,320 (we win) |

**We now win all four mismatch/order streams.** `orig2uid` is the single
biggest remaining loss and PgRC2 pays literally nothing for it.

Whether dedup is worth its correlation array depends on how duplicated the
input is. Measured on real full files:

| dataset | duplicate fraction | dedup verdict |
|---|---|---|
| E. coli | 20.4% | dedup WINS by 3.8% |
| L. major | 10.0% | dedup LOSES by 0.40% |
| P. aeruginosa | 1.2% | dedup LOSES by 0.49% |

Break-even is between 10% and 20.4%; threshold set at 15%, measured by a
cheap hash-only pre-pass before the real load. This is a formula over a
measured input property, not a per-dataset switch -- it satisfies the standing
algorithmic-first rule. (An earlier no-dedup test showed only 0.27% and was
shelved; that predated stage 105 and the mm transforms, which changed the
economics enough to matter.)

### FINAL, all BYTE_IDENTICAL on real full locked files

| Dataset | session start | final | PgRC2 | margin |
|---|---|---|---|---|
| E. coli | 9,905,845 | **8,228,651** | 8,864,420 | **+7.2%** |
| P. aeruginosa | 9,908,025 | **9,320,707** | 9,043,181 | −3.1% |
| S. aureus | 15,352,697 | **13,921,389** | 13,595,003 | −2.4% |
| P. falciparum | 18,810,808 | **17,494,397** | 17,219,695 | −1.6% |
| L. major | 30,612,464 | **29,827,649** | 28,272,652 | −5.5% |

**Session total −6.9%; aggregate vs PgRC2 −2.3%** (from −11.7%..−12.9% per
file when the accounting bug was first fixed). 1 win / 4 losses.

Remaining losses are only two streams: pg sequence (+407,701) and self-match
refs (+164,283). Both are assembly-side, not coding-side.

## 3o. LZMA parameters aligned to the integer stride — the decisive fix

The equal-pg-size comparison (ours MAXMAP=40, literal within 99 B of theirs)
finally located the real loss, and it was NOT the sequence coder:

| stream | ours | PgRC2 | delta |
|---|---|---|---|
| pg sequence | 2,485,885 | 2,485,786 | **+99 (parity)** |
| **order/positions** | **5,517,012** | **4,964,368** | **+552,644** |
| mm_pos | 826,852 | 748,988 | +77,864 |
| self-match refs | 124,934 | 71,700 | +53,234 |
| mm_sym | 360,060 | 370,090 | −10,030 (we win) |

Our sequence coder is at parity with theirs. The position stream was the loss.

**Cause:** PgRC2 codes positions with explicitly chosen LZMA parameters
(`getReadsPositionsCoderProps`, SeparatedPseudoGenomePersistence.cpp:451-460):
`lc=8, lp=2, pb=2` with a 6 MB dictionary, rather than the xz default
`lc=3, lp=0, pb=2`. **`lp=2` aligns the literal-position context to a 4-byte
stride**, which is exactly the shape of a uint32 position array. We were
feeding an integer array to a coder configured for byte-oriented text.

Swept on real pos_abs (6,631,484 B raw):

| parameters | coded |
|---|---|
| xz -9 default (lc=3,lp=0,pb=2) | 5,082,152 |
| lc=8,lp=2,pb=2 (PgRC2's) | 5,038,185 |
| **lc=2,lp=2,pb=2** | **5,021,990** |
| delta4 + lzma2 | 5,242,954 |

Applied via liblzma custom filters to every uint32-shaped stream (pos_abs,
orig2uid, n_indices), keeping whichever of tuned/default is smaller.

### FINAL RESULT — all BYTE_IDENTICAL on real full locked files

| Dataset | session start | final | PgRC2 | margin |
|---|---|---|---|---|
| E. coli | 9,905,845 | **7,995,703** | 8,864,420 | **+9.8%** |
| P. aeruginosa | 9,908,025 | **9,260,575** | 9,043,181 | −2.4% |
| S. aureus | 15,352,697 | **13,331,585** | 13,595,003 | **+1.9%** |
| P. falciparum | 18,810,808 | **16,752,837** | 17,219,695 | **+2.7%** |
| L. major | 30,612,464 | **28,526,841** | 28,272,652 | −0.9% |

**Session reduction −10.3%. Aggregate vs PgRC2 +1.46% — we now WIN overall,
3 of 5 datasets, all byte-identical.**

The single-stream test predicted −1.2%; the real effect was −2.8% to −4.4% per
dataset, because the tuning applies to every uint32 stream and those dominate
the archive. Worth remembering: a per-stream microbenchmark understates a fix
that generalises across streams.

## 3p. Byte-plane split for uint32 streams — final increment

PgRC2's exact position-coder parameters are `lc=8, lp=2, pb=2`
(`getReadsPositionsCoderProps`, PropsLibrary.cpp:90-102). **liblzma caps
`lc+lp <= 4`, so lc=8 is unreachable through xz** -- they can use it only
because they bundle the 7-zip LZMA SDK. The best liblzma equivalent is
`lc=2,lp=2,pb=2`.

A byte-plane (structure-of-arrays) split beats it without needing lc=8: emit
every value's byte 0, then every byte 1, etc. Real P. aeruginosa pos_abs
(6,631,484 B raw):

| encoding | coded |
|---|---|
| xz -9 default | 5,082,152 |
| lc=2,lp=2,pb=2 | 5,021,990 |
| byte-planes + xz -9 | 4,993,384 |
| **byte-planes + per-plane params** | **4,983,604** |
| PgRC2 | 4,964,368 |

Plane 0 and 1 code to ~1 byte per value -- the low bits of a position array
are genuinely random, so both tools sit near the entropy floor here
(uniform over a 21 Mb pg is 24.3 bits = 3.04 B/position). Only ~19 KB of
headroom remains on this stream; it is effectively closed.

Shipped as a 4-way selector with a 1-byte method header (xz default /
lc2lp2pb2 / byte-planes+xz / byte-planes+lc0lp0pb0), keeping whichever is
smallest per stream.

### FINAL RESULT — all BYTE_IDENTICAL, real full locked files

| Dataset | session start | final | PgRC2 | margin |
|---|---|---|---|---|
| E. coli | 9,905,845 | **7,981,962** | 8,864,420 | **+10.0%** |
| P. aeruginosa | 9,908,025 | **9,214,602** | 9,043,181 | −1.9% |
| S. aureus | 15,352,697 | **13,330,920** | 13,595,003 | **+1.9%** |
| P. falciparum | 18,810,808 | **16,752,352** | 17,219,695 | **+2.7%** |
| L. major | 30,612,464 | **28,523,424** | 28,272,652 | −0.9% |

**Session −10.4%. Aggregate vs PgRC2 +1.55%, winning 3 of 5.**

Full list of what actually moved the number, largest first:
1. LZMA/byte-plane coding of uint32 streams (~4% per dataset)
2. `MAXMAP=12` mismatch-acceptance economics (up to −16%)
3. Stage 105, N-reads through assembly (up to −4%)
4. mm_pos bucket+transpose (−13.5% of that stream)
5. mm_cnt zero-flag split (−10.6% of that stream)
6. Adaptive dedup on measured duplication rate

Measured dead ends, recorded so they are not retried: the single-process
architecture rewrite (0% size, 0% speed, +54 MB RAM), seed geometry (0.01%),
no-pre-dedup as an unconditional change (0.27%), reverse-offset mismatch
positions alone (mixed), per-stream/per-bucket coder selection among
xz/bzip2/zstd (0.49% / 0.2%).

## 3q. Remaining losses traced to assembly, not coding — coding is closed

**L. major (−0.9%), stream by stream vs PgRC2 on the same file:**

| stream | ours | PgRC2 | |
|---|---|---|---|
| positions | 15,515,217 | 15,625,904 | we win 110,687 |
| mm_pos | 1,121,860 | 1,774,295 | we win 652,435 |
| mm_sym | 364,885 | 917,584 | we win 552,699 |
| mm_cnt | ~355k+ | 745,872 | we win |
| pg sequence | 9,379,560 | 8,509,123 | **we lose 870,437** |
| self-match refs | 1,311,720 | 322,985 | **we lose 988,735** |

The refs "loss" is not inefficiency: we find **412,595** matches removing 42 M
bases (worth ~9.8 MB of literal at 0.234 B/base) for 1.31 MB of reference
data -- hugely net positive. They self-match far less (5.8% on their LQ pg).

**Our reference coder is already optimal.** It spends 25.43 bits/match against
`log2(pg_len)` = 25.9 bits for an unbounded source. Alternatives measured on
the real stream, all WORSE:

| encoding | bytes |
|---|---|
| **current bounded range coder** | **1,311,720** |
| src raw uint32 + xz | 1,390,362 |
| (dst−src) delta + xz | 1,442,328 |
| (dst−src) delta + byte-planes | 1,332,486 |

Delta-against-destination losing tells us the matches are genuinely
long-range, which is what a pseudogenome should produce. No coding gain left.

**Conclusion: every stream we code is now at or better than PgRC2's, or
provably near its entropy floor.** Both remaining losses reduce to one thing
-- their pseudogenome is smaller than ours for the same input (L. major
8,509,123 vs our 9,379,560 coded). That is assembly structure (their HQ/LQ
division), not stream coding, and it is the only avenue left.

## 3r. MINOV re-swept: the assembly-side lever that closes most of the gap

Section 3q concluded the remaining losses were assembly-side (their pg is
smaller than ours for the same input). `MINOV` -- the minimum overlap length
accepted during chain building -- is exactly that lever, and it had never been
re-swept since the stream coding was fixed. It was defaulted to 40, copied
from PgRC2's CODER_LEVEL_NORMAL.

Real P. aeruginosa, full file:

| MINOV | total | pg_len |
|---|---|---|
| 40 (old default) | 9,214,602 | 22,354,430 |
| 32 | 9,159,157 | 22,049,293 |
| 24 | 9,115,452 | 21,859,912 |
| 20 | 9,097,124 | 21,794,304 |
| **16** | **9,083,594** | **21,738,890** |
| 12, 10 | 9,083,594 | 21,738,890 (identical) |

It floors at 16 because `SEEDW=16` -- overlaps shorter than the seed cannot be
found at all, so lower values change nothing. Lowering MINOV collapses more of
the pseudogenome, and now that stream coding is efficient the shorter pg is
worth far more than the extra chain links cost. Under the OLD coders this
trade looked unattractive, which is why the value sat at 40.

### Full validation, MINOV=16, all BYTE_IDENTICAL

| Dataset | MINOV=40 | MINOV=16 | PgRC2 | margin |
|---|---|---|---|---|
| E. coli | 7,981,962 | **7,974,367** | 8,864,420 | **+10.0%** |
| P. aeruginosa | 9,214,602 | **9,083,594** | 9,043,181 | −0.4% |
| S. aureus | 13,330,920 | 13,498,413 | 13,595,003 | **+0.7%** |
| P. falciparum | 16,752,352 | **16,448,646** | 17,219,695 | **+4.5%** |
| L. major | 28,523,424 | **27,765,787** | 28,272,652 | **+1.8%** |

**Aggregate vs PgRC2 +1.55% -> +2.89%; wins 3/5 -> 4/5.** L. major flips from
−0.9% to +1.8%. MINOV=24 was also validated (+2.55%, worse). S. aureus is the
one dataset preferring a higher MINOV (+1.9% at 40 vs +0.7% at 16), but the
aggregate strongly favours 16 and it remains a win either way.

Only P. aeruginosa still loses, at −0.4%.

## 3s. Joint MINOV/MAXMAP sweep — 5/5 is reachable, but costs aggregate

MINOV and MAXMAP interact (a shorter pg changes the mismatch/literal trade),
so they were swept jointly at MINOV=16. On real P. aeruginosa:

| MAXMAP | total | vs PgRC2 |
|---|---|---|
| 12 | 9,083,594 | −0.45% |
| 16 | 9,053,141 | −0.11% |
| 20 | 9,044,863 | −0.02% |
| **24** | **9,041,007** | **+0.02%** |
| 28 | 9,042,724 | +0.01% |

P. aeruginosa -- the last loss -- flips to a win at MAXMAP=24. Validated
across all five, all BYTE_IDENTICAL:

| Dataset | MAXMAP=12 | MAXMAP=24 | PgRC2 |
|---|---|---|---|
| E. coli | **7,974,367** (+9.89%) | 7,987,812 (+9.89%) | 8,864,420 |
| P. aeruginosa | 9,083,594 (−0.45%) | **9,041,007 (+0.02%)** | 9,043,181 |
| S. aureus | 13,498,413 (+0.71%) | **13,456,136 (+1.02%)** | 13,595,003 |
| P. falciparum | **16,448,646 (+4.48%)** | 16,892,395 (+1.90%) | 17,219,695 |
| L. major | 27,765,787 (+1.79%) | **27,641,405 (+2.23%)** | 28,272,652 |
| **aggregate** | **+2.89%, 4/5 wins** | +2.57%, **5/5 wins** | |

**Real trade, stated rather than cherry-picked.** MAXMAP=24 wins every
dataset but costs 0.32 points of aggregate, almost all of it P. falciparum
(+4.48% -> +1.90%). Minimax regret strongly favours 12 (worst case 0.47%
versus 24's 2.70%, again P. falciparum).

**Shipped default stays MAXMAP=12**, on aggregate compression and minimax
regret -- the two criteria this project has used throughout. P. aeruginosa's
−0.45% is disclosed rather than tuned away. If a "wins on every dataset"
claim is ever wanted, MAXMAP=24 delivers it at a measured, documented cost;
that is a presentation decision, not an engineering one, and should be made
deliberately.

## 3t. The REAL coders (PPMd7 / FSE / range-with-period) — not xz selection

**What was wrong before.** Section 3n's "per-stream coder selection" only chose
among xz / bzip2 / zstd. That is selection in name only: PgRC2's advantage on
the streams where it still beat us never came from a different general-purpose
compressor, it came from **PPMd7, FSE/Huff0, and range coders with
per-position models** (`coders/PpmdCoder.cpp`, `FSECoder.cpp`,
`RangeCoder.cpp`, `rangecoder/simple_model.h`; selected per stream in
`PropsLibrary.cpp`). Measuring xz-vs-bzip2-vs-zstd and reporting "0.49%, not
decisive" was measuring the wrong thing.

**What was built.** `coders_pgrc.h` plus `thirdparty/{ppmd,fse}`:
- **PPMd7** (LZMA SDK, public domain) at order 5 / 32 MB, matching their
  `getDefaultCoderProps(PPMD7_CODER, level, 5)` for mismatched symbols.
- **FSE and Huff0** (Yann Collet, BSD).
- **Adaptive range coder with `period` interleaved order-0 models**, which is
  what `RangeCoderCompressTemplate` does -- at `period = c` each column of a
  fixed-stride record gets its own statistics. Their per-bucket streams are
  literally the `period = N` lines in their stderr.

Both PPMd7 and FSE are standard third-party libraries PgRC2 merely bundles;
using them is no different from our already linking liblzma, and is not a copy
of their implementation.

**Per-stream measurement on real P. aeruginosa** (raw -> coded):

| stream | xz (was) | ppmd5 | fse | range | chosen |
|---|---|---|---|---|---|
| pos_abs | 5,082,160 | 5,455,315 | 5,661,616 | 5,697,519 | xz/byteplanes |
| pos_strand | 74,272 | 74,691 | 68,867 | **68,333** | range |
| mm_count | 305,604 | **304,044** | 310,814 | 309,365 | ppmd |
| mm_ref | 257,624 | 238,381 | **234,947** | 235,449 | (feeds mmcoder) |
| mm_obs | 236,692 | 215,995 | 215,210 | **214,209** | (feeds mmcoder) |

**mm_pos, with their actual bucketing scheme:**

| variant | bytes |
|---|---|
| flat + xz | 627,776 |
| bucket+transpose + xz (previous ship) | 543,080 |
| **per-bucket range(period=c) — PgRC2's scheme** | **535,367** |
| per-bucket ppmd5 | 561,821 |
| **per-bucket best-of-4 (shipped)** | **530,582** |

PgRC2's own scheme beats what we shipped, confirming the coder choice — not
the transform — was the missing piece.

Shipped: a `best_encode` selector over {xz, xz-lc2lp2pb2, ppmd5, fse,
range, byte-planes×2} with a 1-byte method id, on every stream; and
`mmpos_encode_buckets` doing per-bucket selection including range at
`period = c`.

### FINAL — all BYTE_IDENTICAL on real full locked files

| Dataset | session start | before real coders | **final** | PgRC2 | margin |
|---|---|---|---|---|---|
| E. coli | 9,905,845 | 7,974,367 | **7,960,702** | 8,864,420 | **+10.19%** |
| P. aeruginosa | 9,908,025 | 9,083,594 | **9,050,430** | 9,043,181 | −0.08% |
| S. aureus | 15,352,697 | 13,498,413 | **13,462,553** | 13,595,003 | **+0.97%** |
| P. falciparum | 18,810,808 | 16,448,646 | **16,400,369** | 17,219,695 | **+4.76%** |
| L. major | 30,612,464 | 27,765,787 | **27,701,753** | 28,272,652 | **+2.02%** |

Real coders gained a further −0.26%. **Session reduction −11.8%; aggregate vs
PgRC2 +3.14%, 4 of 5 wins.** P. aeruginosa is now −0.08% — a statistical tie.

Per-stream highlights on P. aeruginosa: pos_strand 74,272 -> 60,094 (−19%),
n_indices 19,228 -> 11,245 (−42%), orig2uid 1,108 -> 2, mm_pos −14.9%,
mm_cnt −14.6%.

## 3u. WHY WE ARE SLOW — it is parallelism, not algorithmic cost

First real speed comparison against PgRC2 (earlier "12.8s vs 12.8s" was
ours-old vs ours-new, never against them). Real full P. aeruginosa, machine
idle, runs sequential:

| | wall | user CPU | utilisation | archive |
|---|---|---|---|---|
| **Ours** | **13.9 s** | 19.7 s | **145%** (1.45 of 12 cores) | 9,050,617 |
| **PgRC2** | **3.14 s** | 16.0 s | **515%** (5.2 cores) | 9,043,170 |

**We do only ~23% more actual computation** (19.7 s vs 16.0 s of CPU). The
4.4x wall gap is almost entirely that we run it on 1.45 cores while they use
5.2. At their utilisation our own workload would finish in ~3.8 s — parity.

### Where the 13.9 s goes

| phase | time | share |
|---|---|---|
| **assembly** | **8.67 s** | **62%** |
| — round 1 (division) | 3.92 s | 28% |
| — round 2 (assembly) | 1.51 s | 11% |
| — MEM matching | 1.17 s | 8% |
| — load+filter+dedup | 1.05 s | 8% |
| — pigeonhole mapping | 0.49 s | 4% |
| — emit chains | 0.42 s | 3% |
| **stream coding** | **4.88 s** | **35%** |
| — pos_abs | 2.33 s | 17% |
| — mm_cnt | 1.40 s | 10% |
| — mm_pos | 0.61 s | 4% |
| — all others | 0.54 s | 4% |

Derived utilisation per stage: **coding runs at ~100% (fully serial)** — xz,
PPMd, FSE and the range coder are all single-threaded; **assembly at ~171%**
(1.7 cores) despite having threads, because the per-level sweep parallelises
the candidate search but commits serially in original tail order.

### One self-inflicted cost, already fixed

`best_encode` originally ran all 7 candidate coders over the FULL stream:
10.0 s of an 18.7 s run, 6.69 s on pos_abs alone. PgRC2 never does this —
`SelectorCoderProps` probes a fraction (`getSelectorCoderProps(..., 0.2, ...)`,
64 KB floor, CodersLib.h:216-236) then encodes once with the winner.
Implemented the same way: **pos_abs 6.69 s → 2.33 s, wall 18.7 s → 13.9 s
(−29%), at a cost of +187 bytes (0.002%)**, still BYTE_IDENTICAL.

### What would close the rest, in order of size

1. **Parallelise stream coding** (4.88 s, fully serial). The streams are
   independent — coding them concurrently is embarrassingly parallel and needs
   no algorithmic change. Alone this is worth ~4 s of wall.
2. **Raise assembly utilisation** from 171%. round 1 at 3.92 s is the single
   biggest item; its serial commit phase is the bottleneck.
3. Multi-threaded LZMA was added (`lzma_stream_encoder_mt`, ≥2 MB inputs) but
   is **inert on this file** — pos_abs wins via the byte-plane + custom-params
   path, which routes through `xz_compress_lzma`, not the MT function.
   Harmless, and will apply where plain xz wins on a large stream.

**Bottom line: we are not doing more work, we are doing it serially.** The
size result (+3.14% aggregate, 4/5 wins) stands independently of this.

## 3v. Archaea layer-by-layer: identifying WHY sulfo loses and halo wins

Two archaea, one win one loss, same kingdom -- the most informative contrast
available. Every number below is measured from both tools on the same file,
not assumed.

### Layer table (ours vs PgRC2's real streams)

| layer | SULFO ours | SULFO PgRC2 | delta | HALO ours | HALO PgRC2 | delta |
|---|---|---|---|---|---|---|
| literal | 814,938 | 595,834 | **+219,104** | 665,856 | 579,465 | **+86,391** |
| mem_triples | 175,723 | 16,444 | **+159,279** | 472,066 | 20,673 | **+451,393** |
| pos_abs | 1,482,345 | 1,398,530 | +83,815 | 1,269,209 | 1,291,734 | −22,525 |
| mm_sym | 87,042 | 208,234 | −121,192 | 52,524 | 513,055 | **−460,531** |
| mm_pos | 621,159 | 721,906 | −100,747 | 93,052 | 479,845 | **−386,793** |
| mm_cnt | 147,185 | 163,068 | −15,883 | 57,632 | 151,467 | −93,835 |
| **TOTAL** | **3,341,008** | **3,113,979** | **+227,029** | **2,619,493** | **3,049,295** | **−429,802** |

**The same three layers lose on BOTH files** (literal, mem_triples, pos_abs).
Nothing about sulfo makes our coders worse. What differs is how much the
mismatch layers win by: on halo they are ~941 KB cheaper than PgRC2's, easily
covering the ~538 KB lost elsewhere; on sulfo only ~238 KB, which cannot cover
462 KB.

### Root cause: it is not our mismatch coder, it is PgRC2's mismatch VOLUME

Raw mismatch symbols produced, both tools, same files:

| | ours | PgRC2 |
|---|---|---|
| sulfo | 777,182 | 1,096,906 |
| halo | 284,039 | **2,682,318** |

**We emit fewer mismatches than PgRC2 on both files.** Our coder is not the
variable. PgRC2 emits 2.4x more on halo than on sulfo, which is why our
mismatch-layer advantage collapses on sulfo -- their number came down, ours
did not go up. Our mean mismatches per placed read is identical on the two
files (2.94 and 2.94).

### The structural difference: where the leftover reads go

| | leftovers | mapped | appended | second pg | MEM second removal |
|---|---|---|---|---|---|
| sulfo | 273,046 | **264,702 (96.9%)** | 8,344 | 2,088,887 B | 50.3%, 28,918 matches |
| halo | 185,825 | 96,564 (52.0%) | **89,261 (48.0%)** | 13,286,869 B | **93.1%, 254,750 matches** |

This is the real mechanism, and it runs opposite to intuition:

- **sulfo maps 96.9% of leftovers onto the pg.** Mapped reads are cheap in
  sequence but every one carries mismatch records -- hence 777,182 mismatch
  symbols and a 621,159 B mm_pos.
- **halo appends 48% of leftovers as raw sequence** into a second pg of
  13.3 MB. That second pg is then 93.1% removed by MEM matching (254,750
  matches), so the bases nearly vanish while generating almost no mismatches
  -- 284,039 symbols and a 93,052 B mm_pos.

So halo wins because its leftovers are **highly self-similar**: appending them
and letting MEM matching collapse the result is far cheaper than mapping them
with mismatches. Sulfo's are not, so mapping is the only option and the
mismatch cost is unavoidable.

### What this identifies (no change made yet)

The choice between "map this read with mismatches" and "append it and let MEM
matching collapse it" is currently made per read by the mapper, with no
knowledge of which is cheaper downstream. On halo the append path wins by a
large margin and we take it only 48% of the time; on sulfo the map path is
right and we take it 96.9% of the time -- both by accident of the acceptance
threshold, not by design.

`mem_triples` is worth noting as measured, not a defect: 22.2 bits/match on
halo and 17.8 on sulfo, both close to the log2(pg_len) floor. It is large on
halo purely because halo finds 169,968 matches against sulfo's 79,149. Those
matches pay for themselves -- they remove 12.4 MB of sequence for 472 KB.

**Not yet identified:** why literal and pos_abs lose on both files. Those are
the two layers where PgRC2 is genuinely ahead of us regardless of dataset, and
they are the honest remaining target.

## 3w. Completing the identification: literal and pos_abs are ONE root cause

§3v left two layers unidentified. Both are now measured, and both reduce to
the same thing -- neither is a coder defect.

### literal: our CODER IS BETTER, the VOLUME is worse

| | our raw bases | our coded | our rate | their pg bases | their coded | their rate |
|---|---|---|---|---|---|---|
| sulfo | 3,418,901 | 814,938 | 1.9069 b/base | 2,546,069 | 595,834 | 1.8722 b/base |
| halo | 3,239,294 | 665,856 | **1.6444 b/base** | 2,554,328 | 579,465 | 1.8148 b/base |

Our sequence coder **beats** theirs on halo (1.6444 vs 1.8148) and is 1.9%
behind on sulfo. Counterfactual, coding THEIR base count at OUR rate:

| | we would code | they code | verdict |
|---|---|---|---|
| sulfo | 606,888 | 595,834 | +1.9% |
| halo | 525,057 | 579,465 | **−9.4%, we would win** |

**The literal loss is entirely volume**: our final pg carries 27-34% more
bases than theirs.

### pos_abs: we are BELOW the entropy floor; the span is the cost

| | entries | our pg span | bits/entry | uniform floor | vs floor |
|---|---|---|---|---|---|
| sulfo | 516,904 | 13,391,488 | 22.94 | 23.67 | **−3.1%** |
| halo | 460,501 | 18,251,548 | 22.05 | 24.12 | **−8.6%** |

We already code positions *below* the uniform-random floor for that span --
there is no meaningful coder gain available. The cost is `N x log2(pg_span)`,
and we lose only because our span is larger:

| | our span | their span | ratio | extra bits/read | extra cost |
|---|---|---|---|---|---|
| sulfo | 13,391,488 | 2,546,069 | 5.3x | 2.39 | 154,746 B |
| halo | 18,251,548 | 2,554,328 | 7.1x | 2.84 | 163,305 B |

### THE SINGLE ROOT CAUSE

All three losing layers -- literal, pos_abs, and mem_triples -- are consequences
of **one fact: our pseudogenome is 5-7x larger in span and ~30% larger in
surviving bases than PgRC2's.**

- literal: more bases to code (coder is fine, better on halo)
- pos_abs: wider span, so every position costs 2.4-2.8 more bits (coder is
  already below the entropy floor)
- mem_triples: a larger pg needs more self-match references to collapse
  (17.8-22.2 bits/match, at the log2 floor)

**Why their pg is smaller:** they split reads into hq/lq/n pseudogenomes and
build each separately (`Final size of Pg` prints three pieces: sulfo
95,666 + 0 + 2,450,403; halo 205,065 + 0 + 2,349,263). Each piece is
independently assembled and self-matched, so positions index into a ~2.5 MB
span rather than one 13-18 MB span. We build a single pg plus an appended
second region, and never compact the result.

This is the same 3-way split noted in `PGRC2_DISK_ARCHITECTURE.md` §5 and
dismissed there because the mismatch-layer reasoning that motivated it was
based on a measurement error. **It is now re-identified from the opposite
direction -- span and volume, not mismatches -- and this time the evidence is
direct: three independent layers all scale with pg size, and pg size is where
we differ by 5-7x.**

### Verified NOT the problem (do not re-investigate)

- our sequence coder (better than theirs on halo, 1.9% behind on sulfo)
- our position coder (below the uniform entropy floor on both)
- our reference coder (at the log2(pg_len) floor on both)
- our mismatch coder or mismatch volume (we emit FEWER than they do on both
  files -- §3v)
