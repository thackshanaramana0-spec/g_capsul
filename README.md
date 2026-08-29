# Independent reimplementation of PgRC2's assembler

Written from the algorithm, not copied. PgRC2's source was read closely to
understand what each stage does — and reading it is what found several of the
results below — but every line here is our own. `method_c/` is the unmodified
upstream clone used as the measurement target; it is GPL-3 and is never
vendored (see `.gitignore` and `../PGRC2_COMPARISON.md`).

Purpose: settle whether PgRC2's speed, memory and assembly quality are
reachable from textbook methods alone, and extract what transfers to ARCS.
These are throwaway measurement instruments, kept for the numbers in them.

Target (yeast_sub.fq, 1M reads, 851,275 unique after N-filter and dedup):
PgRC2 = 22.60M pseudogenome, 4.0 s, 232 MB.

**CORRECTION (2026-08-28): 22.60M is their pg BEFORE their MEM-removal stage.**
Read off PgRC2's own stdout on this exact input:

    Target pseudogenome length: 21104500          <- HQ pg
    Destination pseudogenome length: 1504035      <- LQ pg
      Final size of Pg: 1230006   (removed 278457, 18.5%)
    Destination pseudogenome length: 2100         <- N pg
      Final size of Pg: 1849      (removed 255, 12.1%)
    Final size of Pg: 11463320    (removed 9679949, 45.8%)   <- HQ self-match

21,104,500 + 1,504,035 + 2,100 = 22,610,635, which is where "22.60M" came from.
Their **final literal is 12,695,175** — they remove 43.8% of the pseudogenome in
stage 7 (`SimplePgMatcher::matchPgsInPg`), a stage stages 01-15 never had at
all. So every "N% short of PgRC2" figure below compares two pre-removal
pseudogenomes. That comparison is internally consistent, but it concealed that
an entire stage was missing rather than merely underperforming. Stage 16 adds
it.

## Results

| # | file | what it adds | pg | notes |
|---|------|--------------|-----|-------|
| 01 | greedy_exact | textbook greedy SCS, cycle avoidance | 49.84M | 40 s; matches their filter+dedup to the read (851,275) |
| 02 | minimizer_tolerant | minimizer-seeded tolerant overlap | 38.88M | tolerance replaces the quality split |
| 03 | pigeonhole_mapping | map leftovers into the pg (q-gram lemma) | 36.92M | + in-place RC rescan |
| 04 | rc_in_greedy | reverse complement during growth | 38.54M | **worse** — refuted |
| 05 | single_pass_overlap | packed 32-base seeds, rolling, outward walk | 34.46M | **40 s → 6 s**; the O(n·L²) fix |
| 06 | multi_candidate | keep backup partners | 35.46M | **worse** — weak merges cost bytes |
| 08 | true_sweep | full descending-length sweep | 36.34M | **worse** — same reason |
| 09 | two_round_division | round 1 labels, round 2 builds | 34.16M | chains 29.90M → 25.13M |
| 10 | scaled_seed_width | part width = readLen/(mm+1) | 31.59M | fixed 32 capped tolerance at 3 |
| 11 | second_pseudogenome | assemble survivors; admit[] bugfix | 29.06M | main pg **20.15M < their 21.10M** |
| 12 | classifier_sweep | round 1 links at widest range | 28.91M | best config |
| 14 | relaxed_division | admit non-singletons (86.7%) | 34.56M | **worse** — strictness is the point |
| 15 | parameterised_seed | sweep seed width | **27.49M** | optimum at 16; 12 and 10 regress |
| 16 | pg_mem_matching | MEM removal within and between pgs | **15.60M literal** | the missing stage 7; see below |
| 17 | relaxed_mapping | acceptance decoupled from seed geometry | **13.36M literal** | 5.3% short of PgRC2; see below |
| 18 | parallel_mem | CSR index, threaded MEM search, early frees | **13.36M literal** | same bytes, 14.3 s / 388 MB; see below |
| 19 | parallel_sweep | threaded sweep search, serial commit | **13.36M literal** | same bytes, 7.3 s / 389 MB; see below |
| 20 | dense_seed_map | sliding seeds at any offset, threaded scan, CSR index | **13.11M literal** | **+3.3% of PgRC2**, 5.8 s / 422 MB; see below |
| 21 | packed_reads | reads 2 bits/base; seed retuned on copMEM's lemma | **12.93M literal** | **+1.8% of PgRC2**, 5.3 s / **214 MB (ahead)**; see below |
| 22 | order_info | track pg position of every original read | (permutation) | their stage 6, never previously attempted |
| 23 | perm_coder | Lehmer + Fenwick + range coder | **2,309,967 B** | **-19.0% vs PgRC2's 2,852,758**; 0.022% off the floor |
| 24 | coprime_scan | sample the query side too | -- | **refuted** — verification-bound, not probe-bound |
| 25 | packed_verify | packed-vs-packed XOR + popcount | -- | mapping 0.95 -> 0.84 s |
| 26 | probe_filter | 2 MB cache-resident presence filter | -- | mapping 0.78 s, total 4.53 s |
| 27 | mismatch_cost | count what first-acceptable placement costs | -- | 10.74 mismatches/read; a ~260 KB stream nobody counted |
| 28 | best_match | keep the fewest-mismatch placement | -- | 4.94/read; total **below** theirs |
| 29 | dna_cm | single order-k context model | 2.0206 bpb | **refuted** — dilution; worse with more context |
| 30 | dna_mix | multi-order logistic mixing + SSE | **1.9174 bpb** | **beats PgRC2's 1.9261**; round trip verified |
| 31 | rc_probe | measure RC-overlap headroom before building it | -- | 9,188 links available; >=367 KB, exceeds the deficit |
| 32 | ref_cost | price the MEM reference streams; sweep MINMEM | -- | **+122,700 B uncounted**; MINMEM tuning worth ~5 KB |
| 33 | maximal_mem | extend matches BACKWARD; MINMEM/MAXCAND swept | -- | **-77 KB**; literal now within 1,774 bases of theirs |
| 34 | lazy_parse | one-position lookahead in the MEM parse | -- | -3,634 B; greedy parsing was not the cause |
| 35 | match_model | CM + match model, tested as a REPLACEMENT for MEM removal | 1.8487 bpb | **refuted as a replacement**; kept as a coder, -12 KB |
| 36 | rc_only_selfmatch | drop forward self-matching, as PgRC2 does | -- | **-36 KB**; forward self-matches were a net loss |
| 37 | ref_coder | sources range-coded bounded by their own destination | -- | **-12.5 KB**, 23.68 bits/match, beats their 23.83 |

## Stage 16 — the missing MEM-removal stage

Adds PgRC2's stage 7 (`SimplePgMatcher::matchPgsInPg`): match the survivor pg
into the main pg, and the main pg into itself, reverse complement enabled,
minimum match 45 (their `CODER_LEVEL_NORMAL` default), replacing each maximal
exact match with an (offset, length) reference.

Not copied: they use CopMEM with sparse sampling, which can miss matches. This
reuses the packed-32-base seed already in the file and samples the SOURCE every
`45-32+1 = 14` positions, which is exact — a match of length >= 45 spans >= 14
consecutive candidate seed starts, so it must cover a sampled one.

Measured on yeast_sub.fq, `./s16 reads.fq 3 40 16`:

| | ours | PgRC2 | |
|---|------|-------|---|
| main pg before removal | 21,619,175 | 21,104,500 | |
| main pg after removal | **11,632,777** (46.2% removed) | 11,463,320 (45.8%) | +1.5% |
| survivor pg after removal | **3,971,773** (58.1% removed) | 1,231,855 | **+222%** |
| final literal | **15,604,550** | 12,695,175 | **+22.9%** |
| time | 26.2 s | 4.0 s | 6.5x |
| peak RSS | 508 MB | 232 MB | 2.2x |

**The 6.5x is not a like-for-like number.** PgRC2's 4.0 s is `-t 12`; every
prototype here is single-threaded. Measured on the same input:

| | time | peak RSS |
|---|------|----------|
| PgRC2 `-t 12` | 4.0 s | 232 MB |
| **PgRC2 `-t 1`** | **13.68 s** | 208 MB |
| ours (serial) | 26.19 s | 508 MB |

Against their serial build we are **1.9x slower**, not 6.5x. PgRC2 v2's headline
speed gain over v1 was parallelisation (their paper: 8-9x compress, ~5x of it
from threading), so most of the remaining distance on this axis is a thing this
progression has not attempted rather than a thing it does worse.

At the untuned seed width 32 the main pg lands at **11,384,786 — 0.7% SMALLER
than theirs** — while total literal is worse (15.79M). So the main assembler is
at or past parity after removal; the seed width trades main-pg quality against
survivor count.

**94% of the remaining gap is the survivor pg** (2.74M of 2.91M). The cause is
upstream, in the division: we admit 593,968 both-side-overlapped reads where
they admit 708,819, leaving us 66,775 survivors against their 11,298.

Note stage 14 refuted relaxing the division — but that was measured BEFORE this
stage existed. The README's own finding says "weak merges cost bytes" inverts
once reads pass a tiling filter. MEM removal is exactly such a filter, so stage
14's conclusion needs re-testing on top of stage 16, not assumed.

## Stage 17 — acceptance is not seed geometry

Stage 16 left 94% of the gap in the survivor pseudogenome. The cause was not
assembly. It was one line of ours tying two unrelated things together:

    const uint32_t NPARTS = MAXMM + 1;      // pigeonhole: 4 parts -> 3 mismatches

That makes the seed layout decide acceptance, so a read is rejected past 3
mismatches because that is all the pigeonhole lemma can certify at that part
width. PgRC2 keeps them separate (`ReadsMatchers.cpp:700`):

    uint8_t maxMismatches = readLength / minCharsPerMismatch;   // 150/3 = 50

Their seed (`readsExactMatchingChars` = 38) only proposes candidates; acceptance
is a separate and far looser test at up to **fifty** mismatches per 150-base
read. Ours stays at a 32-base seed because that is what a uint64 holds exactly;
a shorter seed proposes more candidates, so it is not a handicap.

| | mapped | survivors | survivor pg |
|---|--------|-----------|-------------|
| stage 16 | 198,155 (74.8%) | 66,775 | 9,482,167 |
| **stage 17** | **249,634 (94.2%)** | **15,296** | **2,095,440** |
| PgRC2 | (96.1%) | 11,298 | 1,504,035 |

Result, `./s17 reads.fq 3 40 16`:

| | ours | PgRC2 | |
|---|------|-------|---|
| main pg after removal | 11,632,777 | 11,463,320 | +1.5% |
| survivor pg after removal | 1,731,204 (17.4% removed) | 1,231,855 (18.5%) | +40.5% |
| **final literal** | **13,363,981** | **12,695,175** | **+5.3%** |
| time | 23.6 s | 13.68 s (`-t 1`) | 1.7x |
| peak RSS | 509 MB | 208 MB | 2.4x |

From 22.9% short to 5.3% short on one parameter, at no time cost (mapping 2.95 s
vs 3.15 s). Our MEM removal rate on the survivor pg (17.4%) now lands on theirs
(18.5%), as the main pg's already did, which is further evidence the matcher is
behaving like theirs rather than accidentally.

Remaining gap is 668,806 B: 499,349 of it (75%) still the survivor pg, 169,457
(25%) the main pg.

Not copied: they re-scan and keep the BEST match per read (the "better-matches"
column in their log, `ReadsMatchers.cpp:315-328`). That shrinks their mismatch
stream, which this progression does not measure, so first-acceptable wins here.

## Stage 18 — the MEM stage parallelised, and its memory cut

Stage 17 left the MEM pass as the single largest time block (9.32 s of 23.6 s)
and the largest allocation. Three changes, none of which touch what gets
matched:

**Flat CSR index.** `unordered_map<uint64_t,vector<uint32_t>>` allocated a
separate heap block per distinct key -- ~1.5M of them -- plus a node header
each. Sorting (key,pos) pairs once and pointing an open-addressing table at the
first occurrence of each key gives identical answers from two flat arrays: after
the sort every occurrence of a key is adjacent, so a lookup is one probe and a
forward scan. Same exact-size-CSR fix already carried into ARCS from here.

**Threaded search.** The greedy parse is serial by nature -- it skips ahead by
each accepted match -- so what is parallelised is the search, exactly as ARCS's
`repeat_elim` does: every query only reads the shared, already-built index.
Threads collect their matches locally and the consumed bitmap is unioned after
the join, so no two threads ever write the same byte. Slicing the destination
can only differ from one serial parse at the T-1 slice boundaries, bounded by
~45 bytes each.

**Early frees.** Nothing after the second pseudogenome reads the reads, the
prefix index, or the sweep's per-read arrays -- the MEM stage works purely on
`pg`. 851,275 std::strings at ~150 bases plus a 32-byte header apiece is
~168 MB, and the prefix index ~80 MB. Released before the MEM index is built so
the two never coexist.

| | stage 17 | **stage 18** | |
|---|----------|----------|---|
| MEM stage time | 9.32 s | **1.01 s** | 9.2x |
| total time | 23.63 s | **14.33 s** | 1.65x |
| peak RSS | 509 MB | **388 MB** | -24% |
| final literal | 13,363,981 | **13,363,981** | **identical** |

The literal is byte-identical, so the boundary effect predicted above cost
nothing measurable: main-pg matches went 52,820 -> 52,824 and removed exactly
the same 9,986,398 bytes.

### Where the three axes now stand

Stage 18's MEM stage calls `std::thread::hardware_concurrency()`, which is 12
on this machine, so **14.33 s is already a 12-thread number** and must be
compared against PgRC2 `-t 12`:

| threads | ours | PgRC2 |
|---------|------|-------|
| 1 | 23.63 s (stage 17) | 13.68 s |
| **12** | **14.33 s** | **4.0 s** |

| | ours | PgRC2 |
|---|------|-------|
| final literal | 13,363,981 | 12,695,175 (+5.3%) |
| time @12 threads | 14.33 s | 4.0 s (**3.6x**) |
| peak RSS | 388 MB | 232 MB (1.9x) |

**At equal thread count we are 3.6x slower.** An earlier draft of this section
claimed "parity" by comparing our 12-thread run against their 1-thread run --
the same apples-to-oranges error this file corrects elsewhere. Only the MEM
stage is threaded here; round 1 (6.4 s), round 2 (3.1 s) and mapping (3.0 s) are
serial, 12.5 s of the 14.33 s. PgRC2 threads all of them.

Peak RSS is now set by the assembly phase, not the MEM phase, so further memory
work has to go there: `reads` as 851,275 separate std::strings is ~168 MB and
the prefix index ~80 MB, both live simultaneously during the sweep.

## Stage 19 — the sweep parallelised without changing a byte

After stage 18 the sweep was 9.5 s of a 14.3 s run. It splits the same way the
MEM stage did: everything expensive (seed roll, hash probe, L-byte memcmp
against each candidate) reads only state fixed for the whole of one L level,
while the two conditions that change *within* a level -- `prv[b]!=NONE` and the
cycle check `ch_h[a]==b` -- are a couple of array loads.

So the search runs on 12 threads and the commit stays serial in the original
tail order. The parallel phase records every candidate passing the static tests
in bucket order; the serial phase applies the dynamic tests to that list and
takes the same first survivor the fully serial loop would have.

**The candidate cap is the part that needed care.** Retaining 8 candidates per
tail is not free: a tail whose first 8 are all claimed would, serially, have
gone on to a 9th. Measured, that cost 14 of 673,334 links and 1,288 bytes of
literal -- 0.0096%, small, but a real regression and not what "identical" means.
The fix is a serial re-probe for exactly the tails whose list filled to the cap
and whose every entry turned out taken. It costs 0.11 s and restores the fully
serial result.

| | serial (17) | +MEM par (18) | **+sweep par (19)** |
|---|-------------|---------------|-----------------|
| round 1 | 6.41 s | 6.41 s | **1.45 s** |
| round 2 | 3.08 s | 3.08 s | **0.73 s** |
| MEM stage | 9.32 s | 1.01 s | 1.06 s |
| mapping | 2.95 s | 2.95 s | 2.95 s |
| **total** | 23.63 s | 14.33 s | **7.28 s** |
| peak RSS | 509 MB | 388 MB | 389 MB |
| final literal | 13,363,981 | 13,363,981 | **13,363,981** |

Links 673,334, both-sides-overlapped 593,968, MEM main 11,632,777 over 52,824
matches -- every one identical to the serial run.

### Three axes, both at 12 threads

| | stage 19 | **stage 20** | PgRC2 | |
|---|----------|----------|-------|---|
| final literal | 13,363,981 | **13,113,591** | 12,695,175 | **+3.3%** |
| time @12 threads | 7.28 s | **5.83 s** | 4.0 s | **1.46x** |
| peak RSS | 389 MB | 418 MB | 232 MB | 1.80x |

Size and speed both closed materially (5.3% -> 3.3%, 1.82x -> 1.46x). Memory
went the wrong way by 29 MB and the profile below explains why -- the denser
seed index costs real memory, and the read storage underneath it was always the
binding constraint.

## Stage 20 — the mapper was looking in four places

Seeds were indexed at fixed multiples of 32 (offsets 0, 32, 64, 96). Two defects
at a 150-base read: a read whose only clean 32-mer starts at offset 17 is
invisible, and bases 128-149 are in no seed at all.

PgRC2 does not do this. At `CODER_LEVEL_NORMAL` the matcher is mode `'c'`,
`CopMEMReadsApproxMatcher` (`ReadsMatchers.cpp:715-735`), which finds maximal
exact matches at ANY offset -- their 38 is a minimum MEM length, not a position.
Worth noting their pigeonhole guarantee is only
`targetMismatches = readLength/exactMatchingChars - 1 = 150/38-1 = 2` against a
`maxMismatches` of 50, so past 2 mismatches their seed is a heuristic filter
too. The difference was never the guarantee -- it was that their filter looks
everywhere and ours looked in four places.

Sliding the seed by 8 gives 15 seeds per read instead of 4. The pg is still
swept once with one probe per position, so only the candidate lists lengthen.

Also threaded the pg scan (disjoint slices, local collection, union on join --
each slice starts SEEDW-1 early so seeds spanning a boundary are not lost).

| | stage 19 | **stage 20** | |
|---|----------|----------|---|
| reads mapped | 249,634 (94.2%) | **251,555 (95.0%)** | PgRC2 96.1% |
| survivors | 15,296 | **13,375** | PgRC2 11,298 |
| survivor pg after MEM | 1,731,204 | **1,480,814** | PgRC2 1,231,855 |
| mapping time | 2.95 s | **1.06 s** | |
| total time | 7.28 s | **5.83 s** | |
| **final literal** | 13,363,981 | **13,113,591** | **+3.3%** (was +5.3%) |

### The index structure mattered more than the parallelism

First cut used `unordered_map<uint64,vector<pair<uint32,uint8>>>` for the denser
index: ~4M entries over ~3.5M distinct keys, each key paying a node header plus
a separately-allocated vector. **Peak RSS went 389 MB -> 673 MB and mapping only
reached 2.12 s.** Replacing it with the same flat CSR + open addressing used for
the MEM index in stage 18 -- sort `(key,rid,part)` once, point a table at each
key's first occurrence -- gave 422 MB and 1.06 s from identical output.

So the dense-seed idea would have looked like a memory regression and a weak
speedup if the container had been left alone. Same lesson as stage 18: on these
sizes the allocation strategy dominates the algorithm.

### Memory: measured, and the remaining gap is structural

Peak RSS profiled per stage (stride 8):

    load+filter+dedup      rss=224 MB      <-- already 96% of PgRC2's total peak
    prefix seed index      rss=238 MB
    round 1 (division)     rss=267 MB
    round 2 (assembly)     rss=285 MB
    emit chains            rss=305 MB
    pigeonhole mapping     rss=377 MB
    pg MEM matching        rss=353 MB
    PEAK                   418 MB

**224 MB is reached before any index is built.** That is read storage:
851,275 `std::string`s at 150 bases each pay a 32-byte header plus allocator
rounding, ~192 B per 150 B of sequence.

Two things were tried against the index side and neither is the answer. Freeing
the prefix index across the mapping stage (it is not consulted there, and is
next needed by a survivor sweep over ~13k reads rather than 851k) moved peak
422 -> 418 MB, because the peak is not there. It was kept anyway -- correct, and
it lowers the mapping-stage plateau -- but it is not a memory fix.

**The real fix is 2-bit packing, which is exactly what PgRC2 does:**
`SymbolsPackingFacility` stores reads packed and compares packed-against-
unpacked directly (`countSequenceMismatchesVsUnpacked`,
`compareSequenceWithUnpacked`). At 4 bases/byte, 851,275 x 150 bases is 32 MB
against our ~163 MB. That single change accounts for essentially the whole
remaining memory gap.

Not implemented here: it touches every comparison in the file (the sweep's
memcmp, the mapper's mismatch loop, the MEM extension), so it is real surgery on
a prototype rather than a local change. Identified and measured, not done.

## Stage 21 — reads at 2 bits per base

The profile in stage 20 showed peak RSS reaching 224 MB at load, before any
index existed -- 96% of PgRC2's entire 232 MB peak, spent purely on holding the
reads. 851,275 `std::string`s pay a 32-byte header plus allocator rounding on
every 150 bases: ~192 B to store 150 B of sequence.

PgRC2 does not do that. `SymbolsPackingFacility` keeps reads packed and compares
packed against unpacked in place (`countSequenceMismatchesVsUnpacked`,
`compareSequenceWithUnpacked`). Same idea here: one flat `vector<uint64_t>` at
32 bases per word, each read starting on a word boundary (wastes under a word
per read, makes a prefix load aligned), plus a word offset and a `uint8_t`
length. Safe because N-containing reads are already filtered, so the alphabet
really is 2-bit.

Comparisons work on the packed form rather than unpacking to compare:
`rcmp(a,off,b,L)` walks both reads 32 bases at a time with a masked tail, and
`rseed(i,o,w)` extracts a seed from an arbitrary base offset by shifting across
the word boundary. Only the mapper's mismatch loop unpacks, into a 256-byte
stack buffer, because it compares against the unpacked pg.

| | stage 20 | **stage 21** | |
|---|----------|----------|---|
| rss at load | 224 MB | **115 MB** | -49% |
| **peak RSS** | 418 MB | **280 MB** | **-33%** |
| round 1 | 1.52 s | 2.42 s | slower |
| total time | 5.83 s | 6.33 s | +0.5 s |
| final literal | 13,113,591 | **13,113,591** | **identical** |

Output is byte-identical and so are the intermediate counts (links 673,334,
both-sides 593,968, MEM main 11,632,777 over 52,824 matches), which is the check
that the packed comparators agree with the memcmp they replaced.

**The 0.5 s is real and worth naming.** `rcmp` does ~5 shifted word loads where
`memcmp` did one SIMD-friendly call over contiguous bytes, and the sweep runs it
20M times. Trading 138 MB for 0.5 s is the right side of that trade here, but it
is a trade, not a free win.

### Three axes after stage 21

| | ours | PgRC2 | |
|---|------|-------|---|
| final literal | 13,113,591 | 12,695,175 | **+3.3%** |
| time @12 threads | 6.33 s | 4.0 s | **1.58x** |
| peak RSS | 280 MB | 232 MB | **1.21x** |

...and after retuning the seed on copMEM's insight (below), at seed 16 /
stride 22:

| | ours | PgRC2 | |
|---|------|-------|---|
| final literal | **12,926,925** | 12,695,175 | **+1.83%** |
| time @12 threads | **5.33 s** | 4.0 s | **1.33x** |
| peak RSS | **214 MB** | 232 MB | **0.92x -- ahead** |

or at seed 16 / stride 12 if size matters more than memory: **12,899,224
(+1.61%)**, 5.54 s, 264 MB.

## copMEM's actual trick, and what it buys

Read `CopMEMMatcher::calcCoprimes` (`matching/copmem/CopMEMMatcher.cpp:111`)
and `initParams` (:69):

    tempVar = L - K + 1              // K-mer start positions inside an L-match
    k1 = sqrt(tempVar)+1; k2 = k1-1  // consecutive => coprime, k1*k2 <= tempVar

They sample the **reference every k1** and the **query every k2**, coprime. Any
match of length >= L spans `tempVar` consecutive K-mer starts, and because
gcd(k1,k2)=1 the Chinese remainder theorem puts a position inside any k1*k2
window where both samplings coincide. So no match of length >= L is missed,
at a fraction of the probes a dense scan would need.

Also `ReadsMatchers.cpp:424` -- they index the **pseudogenome** and stream each
read against it. We index the reads and stream the pg. Both are valid
orientations of the same lemma.

Ours is the degenerate case k1 = stride, k2 = 1. What matters is the resulting
**sensitivity floor = SEEDW + STRIDE - 1**, the shortest exact read/pg stretch
guaranteed to be found:

    ours (stage 20/21):  32 + 8  - 1 = 39
    copMEM at L=38:      K=28, k1*k2=10  ->  28 + 10 - 1 = 37

That reframes the tuning. Shortening the STRIDE improves the floor but inflates
the index. Shortening the SEED improves the floor *and shrinks the index*,
because fewer seeds per read are indexed at a wider stride. Measured:

| SEEDW | stride | floor | survivors | final literal | time | peak RSS |
|-------|--------|-------|-----------|---------------|------|----------|
| 32 | 4 | 35 | 13,113 | 13,076,725 | 5.98 s | **431 MB** |
| 32 | 8 | 39 | 13,375 | 13,113,591 | 5.38 s | 280 MB |
| 28 | 10 | 37 | 12,831 | 13,035,869 | 5.21 s | 274 MB |
| 24 | 14 | 37 | 12,479 | 12,985,655 | 5.26 s | 253 MB |
| 20 | 18 | 37 | 12,222 | 12,948,614 | 5.20 s | 245 MB |
| 24 | 10 | 33 | 12,362 | 12,968,854 | 5.27 s | 275 MB |
| 20 | 8 | 27 | 11,980 | 12,914,345 | 6.01 s | 337 MB |
| 18 | 6 | 23 | 11,824 | **12,893,559** | 6.12 s | 383 MB |
| **16** | **12** | **27** | **11,867** | **12,899,224** | 5.54 s | 264 MB |
| **16** | **22** | **37** | 12,065 | 12,926,925 | 5.33 s | **214 MB** |
| 16 | 26 | 41 | 12,199 | 12,945,818 | **5.08 s** | 212 MB |

Stride 4 is the trap: it buys 262 reads for 150 MB. Seed 20 at stride 18 beats
it on every axis -- 890 fewer survivors, 128 KB smaller, 186 MB lighter, faster.

Two points are worth naming. **Seed 16 / stride 12** is the best size available
at sane memory: it is within 6 KB of the smallest archive in the table while
using 119 MB less than it and running faster. **Seed 16 / stride 22** trades
28 KB of size for another 50 MB, and is the first configuration whose peak RSS
falls *below PgRC2's*.

Note the last three rows share a floor of 37 yet keep improving. The floor is a
worst-case guarantee, not a description of typical behaviour: a shorter seed is
also more likely to land inside a clean stretch of a read that carries
mismatches, so it wins beyond what the guarantee alone predicts.

## Stages 22-23 — OrderInfo, their stage 6, and the biggest win yet

Stage 6 was the one PgRC2 stage this progression had never attempted, and it is
the gate for losslessness: without it an archive cannot restore the original
file order. Stage 22 tracks the pg position of every original read (chained,
mapped, and survivor, with the original->unique map kept across dedup) and emits
the inverse permutation. Stage 23 codes it.

**What they do.** `SeparatedPseudoGenomePersistence.cpp:226-233`, single-end
with complete order info:

    for (i...) rev[orgIdxs[i]] = i;
    writeCompressed(pgrcOut, (char*)rev.data(), rev.size()*sizeof(uint32), LZMA);

A raw uint32 array handed to LZMA. No permutation-specific coding at all. From
their own log on yeast_sub.fq:

    lzma ... compressed 4000000 bytes to 2852758 bytes (ratio 0.713)

**A false start worth recording.** The first estimate of their order cost came
from subtracting archives, 7,063,459 - 4,898,620 = 2,164,839 B. That is *below*
log2(999340!) = 2,309,466 B, which is impossible for a permutation, and the
impossibility is what exposed the error: `-o` is `PGRC_ORD_SE_MODE` and default
is `PGRC_SE_MODE`, two pipelines differing in more than order, so the delta
attributes other changes to order. Their real cost is the 2,852,758 B their log
prints. Had the subtraction not produced an impossible number it would have
been believed, and the conclusion drawn from it -- that their permutation is
structured and a uniform coder would lose -- was exactly backwards.

**Why a general compressor cannot win here.** LZMA does not know the array is a
permutation. Once you do know every value occurs exactly once, element i only
has to be identified among the n-i values still unused, so the content is
sum log2(i) = log2(n!) = 18.49 bits/read, against the 22.82 they spend. Our own
permutation is random-like (47.4% descents), and xz on it gives 2,605,236 B --
better than their LZMA but still 13% above the floor, which is the point: no
general-purpose coder reaches it.

**Stage 23** maps the permutation to its Lehmer code (each element replaced by
its rank among values not yet used) and codes digit i as a uniform symbol in
[0, n-i), ranks coming from a Fenwick tree over still-available values. Textbook
throughout -- factorial number system, Fenwick tree, LZMA-style carry-cached
range coder. Their single-end path has no permutation coder to copy.

| | bytes | bits/read | |
|---|-------|-----------|---|
| raw uint32 | 3,997,360 | 32.00 | |
| PgRC2 (LZMA) | 2,852,758 | 22.82 | |
| xz on the same array | 2,605,236 | 20.85 | |
| **stage 23** | **2,309,967** | **18.49** | **-19.0% vs PgRC2** |
| log2(n!) floor | 2,309,466 | 18.49 | overhead **0.022%** |

**501 bytes off the floor, and the decoder was run**: it rebuilds the
permutation from the coded stream and every element is compared before any size
is reported. The input is also checked to be a genuine permutation first, since
the Lehmer mapping is meaningless otherwise.

For scale: this saves 542,791 B, which is more than twice the entire remaining
sequence-coding gap (231,750 B) that stages 16-21 were fighting over. Order
information was the largest single lever on the board and it sat untouched
because it was the one stage nobody had built.

## Final check on the two axes we still lose

### Speed: 1.46x, not 1.33x, and mapping is the whole story

Their real total on this input at `-t 12` is **3.66 s**, from their own log --
not the 4.0 s quoted from the older comparison doc. Every "1.33x" above should
read **1.46x**. Per stage:

| stage | ours | PgRC2 | |
|-------|------|-------|---|
| load / filter / dedup | 0.57 s | 0.46 s | 1.24x |
| round 1 (division) | 1.32 s | 1.21 s | **1.09x, parity** |
| round 2 (assembly) | 0.74 s | 0.56 s | 1.33x |
| emit chains | 0.22 s | -- | |
| **pigeonhole mapping** | **1.06 s** | **0.41 s** | **2.6x** |
| pg MEM matching | 0.83 s | ~0.5 s | 1.7x |

Mapping is the loss. Their paper notes the RC-match search "is currently
serial", so our stage 7 is already parallel where theirs is not; the deficit is
concentrated in stage 4.

**The obvious fix does not work.** Stage 20 took copMEM's lemma for sensitivity
(sampling the reads every k1) but left the query side at k2 = 1, probing all
21.6M pg positions. Their k1*k2 samples BOTH sides, which is why their matching
is 0.41 s. Sampling our query too, holding the floor constant:

| k1 | k2 | floor | mapping | literal | peak RSS |
|----|----|-------|---------|---------|----------|
| **22** | **1** | 37 | **0.95 s** | 12,926,925 | **239 MB** |
| 11 | 2 | 37 | 1.09 s | 12,919,342 | 288 MB |
| 7 | 3 | 36 | 1.36 s | 12,916,013 | 344 MB |

Monotonically **worse**. Halving the probes doubles the read index, so each
surviving probe carries twice the candidates -- and our verification (a full
read against the pg at up to 50 mismatches) costs far more than the hash probe.
The stage is verification-bound, not probe-bound, so trading probes for
candidates loses. copMEM's balance is the mirror image: they index the huge pg
and query small reads, so probes dominate for them and sampling both sides pays.

k2 = 1 is therefore the right operating point *for this orientation*, and the
lemma does not transfer to the speed axis the way it transferred to sensitivity.
Flipping speed would mean flipping the orientation -- indexing the pg and
streaming reads, as they do -- which is a rewrite of stage 4, not a parameter.

### Sequence: +101,610 B coded, and most of it is the stand-in coder

| | ours | PgRC2 |
|---|------|-------|
| literal | 12,926,925 | 12,694,903 |
| coded | 3,158,084 (2-bit + xz -9e) | 3,056,474 (VarLenDNA + LZMA) |
| bits/base | 1.955 | 1.926 |

232,022 more literal bases, but literal compresses ~4:1, so the algorithmic part
of the gap is ~56 KB coded; the other ~45 KB is our stand-in coder being weaker
than their tuned one. Neither is a missing stage -- it is mean overlap, already
diagnosed as a layout problem.

### Speed, chased further: 1.46x -> 1.24x, and the bottleneck was misdiagnosed twice

Stage 24 refuted "fewer probes". Two more attempts located the real cost.

**Stage 25 -- packed verify.** Stage 21 packed the reads but then, for every
candidate, unpacked all 150 bases into a stack buffer and compared byte by byte:
~300 operations per candidate, run ~26M times. PgRC2 never unpacks
(`SymbolsPackingFacility.cpp:344-362`) -- it packs the pattern a word at a time
and compares whole packed words, falling back to per-symbol only when a word
differs. Packing the text once per direction goes further: both sides packed
makes a 150-base comparison 5 loads, 5 XORs and 5 popcounts. **Mapping
0.95 -> 0.84 s.** Real, but only 12% -- so verification was not dominant either.

**Stage 26 -- the probe filter.** 43M probes (21.6M pg positions x 2 directions)
against ~24 MB of CSR arrays plus a 16 MB open-addressing table is past L3, so
essentially every probe is a memory access. That is also the retrospective
explanation for stage 24: it halved the probes but doubled the table, and the
table was the problem. Almost every probe misses, so a one-bit-per-slot presence
filter -- 2 MB, cache-resident -- rejects misses without touching the 40 MB
structure. A false positive costs one wasted lookup and cannot change the result.
**Mapping 0.84 -> 0.78 s, total 5.33 -> 4.53 s.**

| | mapping | total | vs PgRC2 3.66 s |
|---|---------|-------|------------------|
| stage 21/22 | 0.95 s | 5.33 s | 1.46x |
| stage 25 (packed verify) | 0.84 s | 5.20 s | 1.42x |
| **stage 26 (+ probe filter)** | **0.78 s** | **4.53-4.75 s** | **1.24x** |

Output byte-identical throughout: 12,926,925 literal, 252,865 mapped, 12,065
survivors. Still behind, so speed does not flip -- but the residual is now 0.9 s,
and the honest read is that most of it went to a cache effect rather than to
anything algorithmic.

**Memory moved the wrong way and it is not free.** Peak RSS is 238 MB against
214 MB at stage 21, because stage 22's order tracking adds `ppos` (n x 8 B) and
`orig2uid` (~4 MB), and stage 25/26 add the packed text and the filter. Against
PgRC2's 232 MB that is 1.03x -- parity, not the win stage 21 had. Adding
OrderInfo costs memory, and the earlier "0.92x, ahead" number was measured
before OrderInfo existed. Both numbers are real; they are just not the same
configuration.

### Whole-architecture check: a stream that was never on the books

Everything above measures pg literal plus order. A real archive also stores,
for every read placed by the mapper, where it sits and how it differs. PgRC2's
log prices that stream:

    Mismatched symbols codes ... 1,369,413 bytes to 265,900   (ppmd)  = 1.55 bits/mismatch
    Mismatches counts (zero flags)  988,702 ->  93,290
    Mismatches counts (non-zero)    214,603 ->  82,951
    Reads list offsets              988,702 -> 683,370

~1.15 MB this progression never counted. And stage 17 had bought its survivor
win precisely by relaxing acceptance to 50 mismatches while taking the FIRST
placement that qualified, with the note that best-match "only shrinks their
mismatch stream, which this progression does not measure". That was wrong: it
shrinks OURS too, and the stream is real whether or not it is measured.

**Stage 27 measured it.** First-acceptable placement: 252,865 reads placed with
**2,715,908 mismatches, 10.74 per read**, against their 1.39. 30.6% of
placements carried 12 or more. At their coded rate that is ~526 KB against their
266 KB -- a ~260 KB liability, more than half the order-coding win, hidden
purely by not being on the books.

**Stage 28 fixes it.** PgRC2 re-scans and keeps the best placement
(`ReadsMatchers.cpp:315-328`). What makes that cheap is that the current best
doubles as the early-exit bound -- a candidate is abandoned as soon as it is
worse than what the read already holds, so better placements also mean less
work per candidate.

| | placed | total mismatches | mean | zero-mm |
|---|--------|------------------|------|---------|
| stage 27 (first acceptable) | 252,865 | 2,715,908 | 10.74 | 9.9% |
| **stage 28 (best match)** | 252,865 | **1,250,111** | **4.94** | **18.2%** |
| PgRC2 | (988,702 in reads list) | 1,369,413 | 1.39 | -- |

Literal, placements and survivors byte-identical: 12,926,925 / 252,865 / 12,065.

**The right comparison is total, not per-read.** Our chained reads carry zero
mismatches by construction, so across the whole read set ours is 1,250,111
against their 1,369,413 -- **8.7% below**, roughly 242 KB against their 266 KB.
The liability becomes a small asset.

Cost is 0.44 s (4.53 -> 4.97 s), and it is not the rescan: skipping reads
already placed at or below 1, 2 or 4 mismatches changes neither the time
(4.97 / 5.16 / 5.14 s) nor the mismatch total (within 0.3%). The cost is the
extra candidate verification itself.

### Three axes, whole architecture, everything measured

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded | 3,158,084 | 3,056,474 | +101,610 |
| order | **2,309,967** | 2,852,758 | **-542,791** |
| mismatch symbols | **~242,209** | 265,900 | **-23,691** |
| **sum** | **5,710,260** | 6,175,132 | **-464,872, 7.5% ahead** |

With stage 30's coder replacing 2-bit + xz:

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded | **3,098,294** | 3,056,474 | +41,820 |
| order | **2,309,967** | 2,852,758 | -542,791 |
| mismatch symbols | **~242,209** | 265,900 | -23,691 |
| **sum** | **5,650,470** | 6,175,132 | **-524,662, 8.50% ahead** |

| | ours | PgRC2 | |
|---|------|-------|---|
| time @12 threads | 4.97 s | 3.66 s | 1.36x |
| peak RSS | 238 MB | 232 MB | 1.03x |

Still not a whole-archive comparison -- their reads-list offsets (683,370 B) and
mismatch-count streams (176,241 B) have no counterpart here, since this
progression computes positions and discards them. But the largest omitted stream
has now been measured rather than assumed away, and it moved the answer.

### The sequence gap, taken apart: two refutations and the real cause

+101,610 B splits almost evenly and the halves have different causes:

    232,022 extra literal bases      -> 55,863 B at their rate   (assembly)
    1.9544 vs 1.9261 bits/base       -> 45,747 B                 (coder)

**Refuted 1: an order-k context model.** The reasoning was that MEM removal has
already stripped every exact repeat of 45+ bases, so the literal is near
repeat-free and LZMA must be paying for match machinery it cannot use, while a
context model would spend everything on predicting the next base. Stage 29
implements it (order-k adaptive counts, same range coder as stage 23):

| k | bits/base |
|---|-----------|
| 8 | 2.0206 |
| 10 | 2.5925 |
| 11 | 2.8329 |
| 12 | 2.5984 |

Worse than 2-bit + xz (1.9544) at every order, and worse the more context it is
given. 12.9M bases cannot populate 4^10 contexts -- at k=10 that is ~3
observations per counter, so the model starves. Single-order CM is the wrong
tool here and the "repeat-free" premise was wrong too: LZMA is still finding
real value.

**Refuted 2: LZMA literal-context tuning.** Their log shows `lc=8, lp=0, pb=0`
for the pg stream, and lc=8 means eight bits of previous-byte context for
literals -- on packed DNA, an order-4-bases model inside the entropy coder. xz
caps lc at 4, which looked like a handicap worth most of the deficit. Measured
on our 2-bit stream:

    lc=3,pb=2 (default)  3,157,542   1.9541 bpb
    lc=4,pb=0            3,155,659   1.9529 bpb
    lc=0,lp=2,pb=0       3,158,366   1.9546 bpb

2,425 B. Not the cause either.

**The real cause: byte alignment, and it is visible in one comparison.**

| | pre-entropy | after LZMA | LZMA's gain |
|---|-------------|------------|-------------|
| PgRC2 VarLen | 3,670,728 B (**2.313** bpb) | 3,056,474 (1.926) | **16.7%** |
| ours 2-bit | 3,231,732 B (**2.000** bpb) | 3,158,084 (1.954) | **2.3%** |

Their intermediate stream is *larger* and still finishes *smaller*. Packing four
bases into a byte at fixed boundaries means a repeat that does not begin at a
multiple of four becomes an entirely different byte string, so LZMA sees almost
nothing -- 2.3%. Their 1-4 base phrases keep boundaries tied to content, and
LZMA still works -- 16.7%. The dictionary compounds it by grouping phrases on
their LAST base (`VarLenDNACoder.cpp:200-219`: all phrases ending A, then C, G,
T), so a code byte's high bits carry the terminal base and the stream has
further structure to match on.

Confirmed at the other extreme: raw one-base-per-byte ACGT, which has perfect
alignment but 8 bits/base to start, gives 2.012 bpb with xz -9e. So 1 base/byte
= 2.012, 4 bases/byte = 1.954, their 1-4 variable = 1.926. Theirs is the
optimum of that trade, not an accident.

**Conclusion:** the coder half of the gap is real, is worth ~45 KB, and is not
reachable with off-the-shelf tools or a naive context model. Closing it means
implementing alignment-preserving phrase coding -- adopting their design
knowingly, which is a different thing from the transfers elsewhere in this file
where ours turned out stronger. The other half, 55,863 B of extra literal, is
assembly quality and still needs a better layout.

### Beating the phrase coder instead of copying it

The choice was: implement their VarLenDNACoder, which caps at 1.9261 bits/base
and is their design, or find something better. The literature says better exists
-- DNA-COMPACT reports **1.838 bits/base on yeast**, and GeCo3 gets its gains the
same way: several finite context models of different orders combined by logistic
mixing, rather than one model or a dictionary.

That is exactly the fix for what killed stage 29. A starved high order
contributes little weight instead of dominating, low orders carry the
prediction, and the mixer learns the balance per position instead of it being
fixed.

Stage 30, standard lpaq-style and textbook throughout -- their path has no
context model to copy: each base as two binary decisions; one hashed table of
12-bit predictions per order (1,2,3,4,6,8,11,14,18,22); combination in the
logistic domain, `p = squash(sum w_i * stretch(p_i))`; mixer weights selected by
(node, order-2 context) and trained online; an SSE/APM stage refining against a
low-order context.

| coder | bits/base | on our 12,926,925 bases |
|-------|-----------|--------------------------|
| single order-8 (stage 29) | 2.0206 | 3,264,940 |
| 2-bit + xz -9e | 1.9544 | 3,158,084 |
| PgRC2 VarLen + LZMA | **1.9261** | (3,112,319 equivalent) |
| **stage 30 mixing** | **1.9174** | **3,098,294** |

**Our coder is now better than theirs**: on the same literal it would save
14,025 B. The 45,747 B coder deficit is gone and is a small surplus.

**Two bugs, and they mattered more than the architecture.** The first version
scored 2.0730 -- worse than a single order. Weights initialised at 1<<14 with
ten models put the dot product near 5000, far past squash's +-2047 domain, so
the mixer saturated on the first symbol; `w = 65536/NM` makes the initial mix a
plain average. And the counter update was `delta*2/rate` with rate starting at
2, a full jump to 0 or 4095 on one observation. Fixing those two moved
2.0730 -> 1.9174. The architecture was never the problem.

Round trip is verified, not assumed: the decoder reruns the identical model with
the same update order and rebuilds all 12,926,925 bases before any size is
reported, on the standard set in stage 23 that a number which has not survived a
decode is not evidence.

**The sequence axis is still +41,820 B behind** -- but the cause is now entirely
the 232,022 extra literal bases, i.e. assembly, not coding. Every byte of the
coder deficit has been recovered.

### The last deficit is assembly, and RC is worth re-testing

After stage 30 the sequence gap is +41,820 B and it is entirely the 232,022
extra literal bases -- coding is settled. Mean overlap is the cause:

    admitted reads 593,968 x 150 = 89,095,200 bases
    pg after chains                21,619,175
    over 546,899 links           = 123.4 bases/link

against the ~133.7 this file records for them.

Their merge loop is not the difference. `GreedySwiping...cpp:196-204` scans
forward through every read whose prefix matches the current suffix, skipping
self-links and cycles until one is acceptable -- structurally the same as our
capped candidate list plus the stage 19 re-probe. Candidate coverage is equal.

The one clearly-identified untested lever is **reverse complement in the
greedy**. This file refutes it twice (stage 04, and again in a second
architecture) -- but both refutations predate MEM removal, exactly the situation
that made stage 14's refutation stale and that this file already flags as
needing a re-test.

Stage 31 measures the headroom instead of assuming either way. Index every
read's RC prefix, then for each chain TAIL after round 2 ask whether its suffix
matches the RC prefix of some chain HEAD:

    chain tails after round 2      47,069
      would link via RC             9,188  (19.5%)
      forward seed existed anyway  37,086

A tail joining a head drops that head's contribution from a full read to
`rlen - L`, so 9,188 links are worth at least 40 bases each at the MINOV floor:
**>=367,520 bases, already more than the 232,022-base deficit**, and ~1.13M at
the 123-base mean.

**This is an upper bound, and the two prior refutations are the reason to say so
loudly.** A head can be claimed once, so some of the 9,188 conflict; taking an
RC link may displace a better forward link; and RC merges put sequence into the
pg in an orientation the later mapping and MEM stages then see differently,
which is the most likely explanation for why RC measured *worse* before. Nothing
here says the realised gain is positive -- only that the raw opportunity is
larger than the deficit, which it was not obvious it would be.

**Cost of finding out.** Orientation has to be carried through chain emission,
the mapping stage, MEM removal and the order permutation, and each of those
currently has a verified byte-identical or round-trip-verified result that the
change would put at risk. That is a real piece of work, not a parameter, and it
is the only thing left on the board.

## Stage 32 — the MEM references were never counted either

`PG_LITERAL` is `lit_main + lit_second`: surviving bases only. The
`(offset, length)` pairs for 59,013 accepted matches cost nothing in it. That is
the same failure stage 27 caught with mismatches -- an aggressive stage looking
free because its cost lands in a stream nobody measures -- and it had been sitting
under every sequence number in this file.

Emitted and priced (destination gaps, source offsets, lengths as varints, xz -9e):

| stream | raw | coded |
|--------|-----|-------|
| gaps | 74,740 | 54,576 |
| sources | 223,528 | **182,168** |
| lengths | 79,213 | 63,136 |
| **total** | 377,481 | **299,880** |
| PgRC2's equivalent | | **177,180** |

**We are +122,700 B worse on this stream**, having 59,013 matches against their
~38,762 (implied by their 155,048 B raw offset stream).

**The source stream cannot be coded better.** 59,013 offsets into a 23.2 Mbase
pg is 59,013 x log2(23,233,953)/8 = 180,494 B, and we spend 182,168 -- within 1%
of the floor. A zigzag delta in destination order was tried and is *worse*
(187,904), because self-match sources have no locality: a repeat points anywhere.
The only way to shrink it is fewer matches.

### MINMEM swept, with the full objective

`MINMEM = 45` was PgRC2's `CODER_LEVEL_NORMAL` default (`pgrc-params.h:145`),
copied and never swept. With references finally priced the trade is computable,
and the memory objection dissolves: `STEP = MINMEM - SEED + 1`, so shrinking the
seed alongside the threshold holds the index size constant.

| MINMEM | literal | coded seq | refs | **total** | peak RSS |
|--------|---------|-----------|------|-----------|----------|
| 45 | 12,953,039 | 3,104,520 | 304,780 | 3,409,300 | 239 MB |
| **36** | 12,787,587 | 3,064,865 | 339,392 | **3,404,257** | 239 MB |
| 30 | 12,657,987 | 3,033,803 | 375,788 | 3,409,591 | 239 MB |
| 24 | 12,515,916 | 2,999,752 | 423,824 | 3,423,576 | 239 MB |

**Worth 5,043 B.** The curve is almost flat and turns over at 36. A simple model
said each match costs ~5.07 B and saves `L * 1.9174/8`, implying break-even near
L = 21 and a large win from lowering the threshold. That was wrong: a source
offset costs ~24.5 bits regardless of match length, so short matches pay nearly
the same reference price for much less literal. Measuring both sides is what
made the flatness visible -- optimising literal alone would have read MINMEM 24
as a 105 KB win when it is a 14 KB loss.

Note the MINMEM=45 row here is 12,953,039 rather than the 12,926,925 recorded
elsewhere: the parameterised seed gives 31/step 15 where the original used
32/step 14. A 26,114-base difference from sampling density alone, worth
remembering before reading small literal differences as algorithmic.

### Corrected position

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded (MINMEM 36) | 3,064,865 | 3,056,474 | +8,391 |
| MEM references | 339,392 | 177,180 | **+162,212** |
| order | **2,309,967** | 2,852,758 | -542,791 |
| mismatch symbols | **~242,209** | 265,900 | -23,691 |
| **sum** | **5,956,433** | 6,352,312 | **-395,879, 6.23% ahead** |

The lead is real but **6.23%, not the 8.50% reported before this stage**. Every
earlier summary in this file that quotes 8.50% or 7.5% predates the reference
stream being counted and is superseded by this table.

## Stage 33 — our matches were not maximal

Stage 32 left the reference stream +162,212 B behind and I attributed the rest to
assembly quality. Both readings were wrong, and one number said so: our matches
averaged **179 bases against copMEM's 250**.

copMEM reports *maximal* exact matches. Ours extended **forward only** from a
seed hit. The source is sampled every `STEP`, so a true maximal match at source
`ss` is first seen at the first sampled position `p >= ss` inside it -- and
extending forward from there captures a SUFFIX, abandoning the leading `p-ss`
bases (about STEP/2 each) to literal. Same match count, shorter matches, more
literal, for 58,782 matches.

Extending backward before accepting costs nothing in references -- the match
count is unchanged, so the streams are unchanged, and the recovered bases come
straight off the literal. Care needed on three bounds: no overlap with the
previous match, no run past the source start, and for SELF_FWD the source must
still end at or before the destination. SELF_RC needs no extra cap, since its
condition `src+L <= qlen-qp-L` is unchanged on both sides by backward extension.

| at MINMEM 45 | literal | refs | total |
|---|---------|------|-------|
| forward only (stage 32) | 12,953,039 | 304,780 | 3,409,300 |
| **maximal (stage 33)** | **12,696,677** | 300,212 | **3,343,288** |

**Literal is now within 1,774 bases of their 12,694,903.** The "232,022-base
assembly deficit" diagnosed over the previous several stages was almost entirely
this defect, not mean overlap. That earlier conclusion is withdrawn.

### MINMEM re-swept on maximal matches

| MINMEM | matches | literal | refs | **total** |
|--------|---------|---------|------|-----------|
| **36** | 66,831 | 12,517,972 | 332,012 | **3,332,257** |
| 45 | 58,782 | 12,696,677 | 300,212 | 3,343,288 |
| 60 | 46,041 | 13,144,569 | 244,152 | 3,394,577 |
| 80 | 35,704 | 13,640,589 | 194,756 | 3,464,064 |

### Two coding ideas refuted, and one of them twice

**MAXCAND is not binding.** Candidates are scanned in index order and the longest
kept, so a cap can hide the best match. At 64 and 512 the output is
byte-identical -- the lists are simply shorter than 64.

**The source stream is at its floor.** 66,831 offsets into a 23.2 Mbase pg has a
~205 KB floor and we spend 212,880 absolute. Two alternatives tried:

    consecutive-source delta   worse  (sources have no locality with each other)
    repeat distance dst-src    210,600 vs 212,880, 1.1%

The second was the better-motivated idea -- every self-match points backwards, so
`dst-src` is an LZ77 match distance, and distances are normally far more skewed
than absolute positions. It barely helps here because a pseudogenome interleaves
reads from across the genome, so a repeat's source is not near its destination.
Kept anyway on keep-the-smaller. The only real lever on this stream remains
fewer matches, which is what the MINMEM sweep optimises.

### Position after stage 33

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded | 3,000,245 | 3,056,474 | **-56,229** |
| MEM references | 332,012 | 177,180 | +154,832 |
| order | **2,309,967** | 2,852,758 | -542,791 |
| mismatch symbols | **~242,209** | 265,900 | -23,691 |
| **sum** | **5,884,433** | 6,352,312 | **-467,879, 7.37% ahead** |

Sequence itself is now **ahead by 56,229 B**. The single remaining loss is the
reference stream, and it is a match-count difference (66,831 against their
~43,190) whose coding is already at the entropy floor.

### copMEM2 is not the lever, and here is why

PgRC2 ships copMEM (`matching/copmem/`), and copMEM2 has since been published
(Bioinformatics 39(5), btad313). Its abstract is explicit about what it improves:
multithreading, a carefully built predecessor query structure, sort-procedure
selection, and a memory-frugal mode -- **execution speed and RAM**. Both tools
find the *same* MEMs; a MEM is defined by the strings, not the finder. Upgrading
versions cannot change their match count or their reference stream, so it cannot
be the reason we spend 332 KB where they spend 177 KB.

That left match SELECTION, which is a real algorithmic difference: we take the
longest match at each position, left to right, which is greedy parsing, and the
classic LZ result is that greedy is not optimal -- a shorter match now can enable
a much longer one a base later, and each fragment costs a reference.

Stage 34 adds lazy matching, the cheap form of lookahead used in DEFLATE: having
found a match at `qp`, look at `qp+1`, and if that is strictly longer leave `qp`
as a literal.

| | matches | literal | refs | total |
|---|---------|---------|------|-------|
| greedy | 66,831 | 12,517,972 | 332,012 | 3,332,257 |
| **lazy** | **65,994** | 12,517,415 | **328,512** | **3,328,623** |

**Worth 3,634 B.** Real but small, so greedy parsing was not the cause either.

### What the match-count difference actually is

At an equal MINMEM of 45: ours 58,782 matches removing 10,537,276 bases (179
each), theirs ~43,190 removing 9,958,661 (231 each). We find 36% more matches and
remove 5.8% more sequence. So we are not failing to find long matches -- we are
additionally finding shorter ones they do not, on a pseudogenome that is not the
same string as theirs (23.2M against 21.1M, built by a different assembler).

Four things have now been refuted on this stream, and the elimination is the
result:

- **MAXCAND** -- byte-identical at 64 and 512, the cap never binds
- **consecutive-source delta** -- worse than absolute; sources have no mutual locality
- **repeat distance `dst-src`** -- 1.1%, because a pseudogenome interleaves reads
  from across the genome, so a repeat's source is not near its destination
- **lazy parsing** -- 3,634 B

and the source stream is within 1% of `n*log2(pg_len)`, its entropy floor. The
remaining +95 KB is a property of which matches exist in our pg versus theirs,
not of how we find, select, or code them.

### Position after stage 34

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded | **3,000,111** | 3,056,474 | **-56,363** |
| MEM references | 328,512 | 177,180 | +151,332 |
| order | **2,309,967** | 2,852,758 | -542,791 |
| mismatch symbols | **~242,209** | 265,900 | -23,691 |
| **sum** | **5,880,799** | 6,352,312 | **-471,513, 7.42% ahead** |

## Stage 35 — could a match model delete the reference stream entirely?

Four attempts to shrink the reference stream failed and the structural check
found zero mergeable matches, so the question became whether the stream needs to
exist at all. MEM removal is only there because an order-22 context model cannot
see a repeat 10 Mbases back. A MATCH MODEL can: hash a recent k-mer to where it
last occurred and, while the match holds, predict the base that followed. Long
repeats would then cost almost nothing AND cost no references, because nothing is
removed and nothing has to be pointed at.

Tested on the pseudogenome BEFORE removal (23,233,953 bases), against
literal-coded + references = 3,332,257 B:

    bases 23,233,953   coded 5,369,074 B   1.8487 bits/base   round trip VERIFIED

**Refuted, and not narrowly.** The per-base rate is genuinely better than our
1.9174 -- and lands on DNA-COMPACT's published 1.838 for yeast, so the coder is
now at literature-competitive rates. But there are 1.86x more bases to code, and
5,369,074 against 3,332,257 is not close. **MEM removal is strongly justified and
the reference stream is worth paying for.**

Why the match model cannot substitute: it follows one pointer, anchored on the
most recent occurrence of a 24-mer, and drops on a single mismatch. In a
pseudogenome the other copy of a repeat is usually not the most recent
occurrence, so the pointer is wrong more often than not. Explicit MEM references
find the right copy by construction. That is the argument for the architecture,
and it is now measured rather than assumed.

### Kept as a coder, though

The same model applied to the literal it was meant to replace:

| coder | bits/base | on 12,926,925 bases |
|-------|-----------|---------------------|
| stage 30 (mixing, no match model) | 1.9174 | 3,098,294 |
| **stage 35 (+ match model)** | **1.9098** | **3,086,012** |

12,282 B, round trip verified. Measured on the stage-22 literal dump; on the
current best literal (12,517,415 at MINMEM 36) the same rate is worth ~11,500 B.

### Position

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded | **~2,988,300** | 3,056,474 | **-68,174** |
| MEM references | 328,512 | 177,180 | +151,332 |
| order | **2,309,967** | 2,852,758 | -542,791 |
| mismatch symbols | **~242,209** | 265,900 | -23,691 |
| **sum** | **~5,869,000** | 6,352,312 | **-483,300, 7.61% ahead** |

## Stage 36 — we were doing a whole matching pass they never do

The reference stream sat at +151,332 B through five refuted attempts to shrink
it. Every one of those attacked how references are found, selected or coded.
None asked why there are 65,994 of them against their 43,190. Reading their
stage 7 line by line answers it:

`SimplePgMatcher::exactMatchPg` (`SimplePgMatcher.cpp:32-42`):

    if (revComplMatching) {
        if (destPgIsSrcPg) {
            string queryPg = reverseComplement(destPg);
            matcher->matchTexts(textMatches, queryPg, ...);
        } else { reverseComplementInPlace(destPg); ... }
    } else
        matcher->matchTexts(textMatches, destPg, ...);

`revComplMatching` is `true` for every call in `matchPgsInPg`, so the `else`
branch never runs there. For the self-match they build
`reverseComplement(destPg)` and match *that*. Their parameter is even named
`minimalReverseComplementedRepeatLength` (the `-p` flag, `PgRC.cpp:197`).

**Their stage 7 is reverse-complement only. They never search for forward
repeats inside the pseudogenome.** We always did both, and that entire extra
pass is where the match-count difference comes from.

It is a net loss. A match only pays when `L * bits_per_base / 8` beats the ~40
bits its reference costs, and forward self-matches in a pseudogenome are mostly
short:

| | matches | literal | seq | refs | **total** |
|---|---------|---------|-----|------|-----------|
| forward + RC | 65,994 | 12,517,415 | 2,988,220 | 328,512 | 3,316,732 |
| **RC only** | **49,766** | 12,730,875 | 3,039,178 | **255,780** | **3,294,958** |

The forward pass removed 213,460 more literal bases -- worth ~50,958 B once
coded -- while costing 72,732 B of extra references.

### MINMEM re-swept, because its optimum moved

The threshold had been tuned with forward matching on, so it was optimising a
different trade:

| MINMEM | matches | literal | refs | **total** | vs PgRC2 |
|--------|---------|---------|------|-----------|----------|
| 45 | 45,166 | 12,896,442 | 236,060 | 3,314,763 | +81,109 |
| 36 | 49,766 | 12,730,875 | 255,780 | 3,294,958 | +61,304 |
| 30 | 53,650 | 12,621,364 | 272,728 | 3,285,763 | +52,109 |
| **24** | 58,908 | 12,506,313 | 294,772 | **3,280,342** | **+46,688** |
| 22 | 60,666 | 12,476,703 | 302,304 | 3,280,805 | +47,151 |

**Together: 3,316,732 -> 3,280,342, a 36,390 B gain**, and the sequence+reference
deficit falls from +83,078 to **+46,688**. 4.40 s, 239 MB.

The lesson is the one this file keeps relearning: five attempts failed because
they all assumed the work was necessary and only argued about how to pay for it.
The question worth asking was whether to do it at all, and their source answered
it in five lines.

### Position after stage 36

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded | **2,985,570** | 3,056,474 | **-70,904** |
| MEM references | 294,772 | 177,180 | +117,592 |
| order | **2,309,967** | 2,852,758 | -542,791 |
| mismatch symbols | **~242,209** | 265,900 | -23,691 |
| **sum** | **5,832,518** | 6,352,312 | **-519,794, 8.18% ahead** |

## Stage 37 — profiling the reference stream by component, not in aggregate

Every previous attempt treated "MEM references" as one number. Split by
component at MINMEM 24, RC-only:

| component | ours | bits/match | PgRC2 | bits/match | delta |
|-----------|------|------------|-------|------------|-------|
| sources | 186,852 | **25.38** | 128,671 | **23.83** | +58,181 |
| lengths | 63,988 | 8.69 | 48,509 | 8.99 | +15,479 |
| dest gaps | 43,932 | 5.97 | 0 | -- | +43,932 |
| TOTAL | 294,772 | 40.03 | 177,180 | 32.82 | +117,592 |

**The gap stream is not a loss, and an earlier note in this file was wrong to
imply it.** Their destination positions are in-band -- a '%' MATCH_MARK sits in
the literal and rides inside the VarLenDNA codebook, which carries phrases like
"T%", "AT%", "GG%". The fair unit is literal plus positions: ours
2,985,570 + 43,932 = 3,029,502 against their 3,056,474, so we are **26,972
ahead** there. Comparing the gap stream against a zero was double-counting.

**Sources are the loss, and a varint was the reason.** log2(23,233,953) = 24.47
bits is the naive floor; we spent 25.38 because a varint near 23M needs four
bytes. Theirs sits at 23.83, *below* their own naive floor, so LZMA was finding
structure we were not.

The structure: in a self-match the source always precedes the destination, and
matches are stored in destination order, so src is uniform in [0, dst) rather
than [0, pg_len). Range-coding it against that bound needs no extra information
-- the decoder knows dst before it reads src.

    naive floor    180,183 B   24.47 bits/match
    bounded floor  174,339 B   23.68
    this coder     174,346 B   23.68   round trip VERIFIED
    varint + xz    186,852 B   25.38   <- replaced
    PgRC2          128,671 B   23.83

**12,506 B**, and the rate now beats theirs. Seven bytes off the bounded floor.

### Two real bugs this uncovered

Building it failed twice, and both failures were defects in the pipeline rather
than in the coder.

**RC destinations were never converted back.** An RC pass reports its
destination in reverse-complement coordinates while the source stays in forward
ones, so the two were not comparable. PgRC2 converts in
`correctDestPositionDueToRevComplMatching` (`SimplePgMatcher.cpp:58-60`); we
stored the raw RC position.

**Cross-match destinations were relative.** The survivor cross-match is handed
`pg.data()+main_pg_end`, so its `qp` counts from the survivor pg, not the
pseudogenome. Stored raw, a cross-match at 0 claimed to be the very start of the
pg. Together these made 12,003 of 58,908 rows violate `src < dst`; after both
fixes, **0 of 58,908**.

Neither changed the literal (12,506,313 throughout), so no size measurement in
this file was wrong because of them -- but any real decoder would have produced
garbage, and only building a coder that depends on the invariant exposed them.

### Position after stage 37

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded | **2,985,570** | 3,056,474 | -70,904 |
| MEM references | 282,266 | 177,180 | +105,086 |
| order | **2,309,967** | 2,852,758 | -542,791 |
| mismatch symbols | **~242,209** | 265,900 | -23,691 |
| **sum** | **5,820,012** | 6,352,312 | **-532,300, 8.38% ahead** |

## Stage 37b — the reference stream is now closed, all three components

After the bounded source coder, MINMEM was re-swept because references got
cheaper (25.38 -> 23.68 bits/match) and a cheaper reference makes more matches
profitable, so the optimum should move down:

| MINMEM | matches | literal | seq | refs (gaps/src/len) | **total** | vs PgRC2 |
|--------|---------|---------|-----|---------------------|-----------|----------|
| 30 | 53,650 | 12,621,364 | 3,013,035 | 260,901 (44,268 / 158,713 / 57,920) | 3,273,936 | +40,282 |
| 26 | 56,988 | 12,544,994 | 2,994,804 | 275,431 (45,760 / 168,639 / 61,032) | 3,270,235 | +36,581 |
| **24** | 58,908 | 12,506,313 | 2,985,570 | 283,686 (46,592 / 174,346 / 62,748) | **3,269,256** | **+35,602** |
| 22 | 60,666 | 12,476,703 | 2,978,501 | 291,382 (47,456 / 179,530 / 64,396) | 3,269,883 | +36,229 |

It did not move -- 24 remains optimal, and every row round-trip verified.

### All three components are now at their floors

| stream | we spend | bits/match | order-0 floor | headroom |
|--------|----------|------------|---------------|----------|
| sources | 174,346 | 23.68 | 174,339 (bounded) | **7 B** |
| lengths | 62,748 | 8.52 | 63,418 (8.61 bits) | **-670 B, below it** |
| gaps | 46,592 | 6.33 | 46,219 (6.28 bits) | **373 B** |

Lengths coming in *below* their order-0 entropy is not an error: varint plus xz
exploits correlation between neighbouring lengths that an order-0 measure cannot
see. Either way there is nothing left to extract by coding.

**PgRC2's own lengths cost 8.99 bits/match against our 8.52, and their sources
23.83 against our 23.68.** We are better on every per-match rate. The entire
remaining +35,602 is match count -- 58,908 against their 43,190 -- and the count
is optimal for our pseudogenome, which the sweep confirms twice under two
different reference prices.

### Confirmation from their source that the bounded coder is the right idea

`resolveMappingCollisionsInTheSameText` (`SimplePgMatcher.cpp:160-165`):

    if (match.posSrcText > match.posDestText) { swap(posSrcText, posDestText); }

They explicitly force `src < dst` by swapping, which is exactly the invariant
stage 37's coder exploits. They then `sort` and `unique` the match list
(`:94-95`) before trimming overlaps. So the structure LZMA was finding in their
offsets is the structure we now code directly.

### Where this leaves the sequence axis

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| literal, coded | **2,985,570** | 3,056,474 | **-70,904** |
| references | 283,686 | 177,180 | +106,506 |
| **sequence total** | **3,269,256** | 3,233,654 | **+35,602** |

Every coding avenue on this stream is now measured and closed: MAXCAND does not
bind, consecutive-source delta is worse, repeat distance gives 1.1%, lazy parsing
3,634 B, a match model cannot replace removal, the varint is fixed, and all three
components sit at their floors. What remains is that our pseudogenome contains
more short repeats than theirs, which is the assembler.

## Stages 40-41 — three targeted memory cuts measured zero, and then 45 MB

Peak RSS sat at 233 MB against PgRC2's 232 MB all session, and three separate
reductions moved it by nothing at all:

| change | what it removed | peak RSS effect |
|--------|-----------------|-----------------|
| `ppos` uint64 -> uint32 (stage 38) | 3.4 MB of array | **0** |
| parallel chunked loader (stage 38) | -- | **+5 MB, worse** |
| packed mapping index (stage 40) | 13 B -> 8 B/entry, and a 30 MB temp | **0** |

Stage 40 is worth describing because it *should* have worked. The index was
`mkey`(8) + `mrid`(4) + `mpart`(1) built through a `vector<E>` of 16 B/entry that
coexists with all three during the copy -- 29 B/entry, ~54 MB at the build
moment. Packing it into one sorted `uint64` of `(key32 << 32) | (rid << 3) |
part` gives 8 B/entry with no temp and no copy: 24 MB -> 14.1 MB, verified by the
index's own printout. Sorting the whole word sorts by key because the key is in
the high bits, so the run-detection and first-occurrence logic are untouched. The
key truncates to 32 bits, exact at the default SEEDW 16 and merely a partial key
above it -- partial is safe, since a collision only produces a candidate and every
candidate is verified against the text.

**It changed peak RSS by zero.** Three real reductions, no effect, which stops
being a coincidence and starts being a measurement: **the peak is not the live
set.** Summing what is actually live during mapping gives ~155 MB against a
measured 232 MB.

### The 77 MB was glibc, and ARCS already knew

`src/vodbg_pg.cpp` says it outright: *"malloc_trim is what makes the release
visible: glibc keeps freed arenas by default, so an earlier measurement recovered
only 52 MB of a 312 MB free until the arena was returned explicitly."*

The sweep's per-level candidate array and the 851k-entry prefix map are both
freed before mapping. Both stay resident. One `malloc_trim(0)` after the free:

| | peak RSS | mapping stage RSS |
|---|----------|-------------------|
| without | 238,828 / 238,908 / 238,932 KB | 232 MB |
| **with** | **193,648 / 193,732 KB** | **186 MB** |

**45 MB, reproducible across reps, output byte-identical** (`PG_LITERAL
12,506,313`).

### Memory is now won

| | ours | PgRC2 | |
|---|------|-------|---|
| peak RSS | **189 MB** | 232 MB | **0.81x -- 19% ahead** |

The lesson is the one this file keeps paying for: three plausible fixes were
built and measured before the *fourth* question -- "is the thing I am shrinking
even at the peak?" -- got asked. The answer was in this repository's own source
the whole time.

## A measurement harness, a trustworthy baseline, and a plan killed by arithmetic

### bench.sh

Wall clock is unusable here -- another session runs concurrently and the same
binary measured 1.32 s and 1.91 s for the same stage. `bench.sh <reps> <bin>
<args>` fixes what can be fixed:

- **min over repetitions**, because contention only ever ADDS time, so the
  minimum is the least contaminated estimate while the mean is dragged by
  whatever else was running;
- **per stage**, so a change is judged on its own timer rather than drowned in
  the other six;
- **spread reported**, so a result is believed only when the spread is small
  against the effect claimed;
- **output equality checked across reps**, since timings of runs that computed
  different things are not comparable.

### Baseline, min of 5, outputs identical

| stage | min | spread | PgRC2 |
|-------|-----|--------|-------|
| load+filter+dedup | 0.56 | 0.03 | 0.47 |
| prefix seed index | 0.16 | 0.03 | -- |
| round 1 (division) | 1.50 | **0.25** | 1.19 |
| round 2 (assembly) | 0.71 | **0.13** | 0.60 |
| emit chains | 0.23 | 0.02 | ~0 |
| **pigeonhole mapping** | **0.83** | **0.02** | **0.40** |
| pg MEM matching | 0.67 | 0.03 | ~0.56 |
| **total** | **4.66** | | **3.70** |

Mapping is both the largest fixable delta and the most measurable one, so it was
the obvious next target. Round 1's 0.25 spread, by contrast, is comparable to its
own 0.31 deficit -- that stage cannot support a claim on this machine.

### The orientation flip, refuted before writing it

The hypothesis: PgRC2 indexes the pseudogenome and streams reads through it while
we do the reverse, and stage 24 already showed that orientation is why coprime
sampling failed for us. Flipping looked structural and worth days.

Modelled first. The sensitivity floor requires `SEEDW + k1*k2 - 1 >= 37`, so with
SEEDW 16 the sampling product is capped at 22, and the two orientations trade
index size against probe count:

    current: index reads/22 = 1.85M entries, scan pg twice = 43.2M probes = 45.1M ops
    flipped, unconstrained best (k1=2, k2=11):            18.0M ops -- 2.5x fewer

That 2.5x is real and irrelevant: k1=2 means 10.8M index entries, **141 MB for
the index alone**, and peak RSS is already 233 MB against their 232 MB. There is
no memory to spend.

At **equal index size** the flip inverts:

    flipped, k1=12 k2=1: index 1.80M entries, probes 79.5M = 81.3M ops
    -> 80% WORSE, and the floor drops to 27 rather than 37

The asymmetry is simply that the pg is 21.6M positions while the reads are 264,930
x 150 = 39.7M bases. The volumes are comparable, so neither orientation is
inherently cheaper -- what differs is which side you can afford to sample densely,
and that is decided by the memory budget, which is spent.

**PgRC2 affords their side of it.** copMEM at their parameters gives k1=5, so
their pg index is ~4.3M entries against our 1.85M, inside the same 232 MB. Their
mapping advantage is a denser index, not a better orientation, and the way to
match it is a more compact index entry -- ours is 13 bytes (key 8 + rid 4 + part
1) where a position-only entry is 4 -- not a rewrite of the stage.

Recorded rather than built: a day of work avoided by ten lines of arithmetic
against a constraint that was already measured.

## Stage 39 — PgRC2's sort-merge sweep, built and measured

The one place the evidence said we were algorithmically behind: their division
runs 1186 ms with a **parallel sort and a serial merge**
(`GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp:6-10`, pstld), while ours
took 1.3-1.9 s fully threaded. Beating us with less parallelism means the
algorithm differs, and no amount of threading closes it. So it was rebuilt.

Their `overlapSortedReadsAndMergeSortSuffixes` is a two-list merge:
- reads sorted by full sequence are already sorted by *every* prefix length, so
  the prefix side is built once and never rebuilt;
- as L drops by one the suffix key loses its **first** base, so the suffix side
  is re-sorted in O(n) rather than O(n log n);
- two sorted lists then merge in a linear pass.

**The algorithm is validated.** Against the hash sweep, on identical input:

| | hash sweep | merge sweep |
|---|-----------|-------------|
| round 1 links | 673,334 | **673,334** identical |
| round 2 links | 546,899 | **547,569** (+670) |
| final literal | 12,506,313 | **12,506,278** (-35 B) |
| comparisons | 20,879,865 | **1,909,923 (11x fewer)** |
| round 1 + 2 wall | **2.15 s** | 6.81 s |

Same-or-better output for a fraction of the comparisons -- and 3.2x slower.

### Three bugs, each found by measuring rather than reasoning

**The re-sort was a no-op.** Bucketing by the dropped base and *concatenating*
A/C/G/T rebuilds sorted-by-(dropped_base, rest), which is the OLD key. Every
level after the first merged an order that had never changed. It showed as 41,564
links at L=149 and almost none across the 130 levels below. The buckets must be
**merged**, not concatenated -- which is exactly what PgRC2 keeps in `ssiOrder` +
`updateSuffixQueue`, machinery that reads like bookkeeping and is in fact the
algorithm.

**Allocation churn was not the cost.** Three pairwise merges with vector locals
looked like the obvious waste; replacing them with one preallocated 4-way merge
moved 10.81 -> 10.52 s. Nothing.

**Random access was the cost.** Comparing through the packed reads means two
loads into a 34 MB array per comparison, ~320M cache misses. Giving each read a
32-base key in a compact array -- prefix key fixed, suffix key rolled one base
per level -- halved it, 10.52 -> 5.31 s.

### Why it still loses, and it is not the algorithm

Parallelising the merge search (stage 19's split: search wide, commit ordered,
slice boundaries found by binary search since the walk is monotone) gave only
5.31 -> 4.32 s. **1.2x from 12 threads.**

The reason is structural: the search is not where the time goes. The per-level
bucket, 4-way merge, compaction and key roll are each O(n) and each inherently
sequential, across 134 levels. Parallelising the one part that parallelises
cannot help when the serial remainder dominates.

**So the hash sweep is the right choice on this hardware**, despite doing 11x
more comparisons: it is embarrassingly parallel and this machine has 12 cores.
PgRC2 reaches 1.79 s with the merge because their per-level maintenance is far
tighter than ours -- but they have no parallel alternative to compare against,
and we do.

Kept behind `MERGE_SWEEP=1`, not made default. The measurement that matters is
the comparison count and the output equality, both of which are deterministic
and therefore trustworthy on a machine where wall clock is not.

## Stage 38 — two speed/memory attempts, both refuted, neither kept

The checkpoint after stage 37 was size won, memory 1.3 MB behind, speed ~1.2x
behind. Two changes were built to close the last two axes. Both failed, and the
failures are worth more than the code was.

**Parallel chunked load.** Measurement said this stage is CPU-bound, not
I/O-bound: the input is fully page-cached (`cat` 0.01 s, `wc -l` 0.02 s on
311 MB) while the serial getline loop took 0.59 s, so ~0.57 s was parsing,
hashing and packing. PgRC2 has no OpenMP anywhere in reads loading either -- only
`CodersLib`, `ReadsMatchers`, `CopMEMMatcher` and the decoder -- so this was ours
to win, and their equivalent stage is 464 ms serial.

It worked as designed (load 0.59 -> 0.51 s, output byte-identical) and still
lost. A first version read the whole 311 MB file at once and took peak RSS from
239 MB to **563 MB** -- 325 MB for 0.10 s. Chunking to a 16 MB window fixed the
blow-up but still measured **~5 MB above baseline**, consistently across reps,
because the window plus per-batch packed arrays outweigh what was saved. Memory
was the axis in play, so paying 5 MB for a time gain smaller than this machine's
noise is the wrong trade.

**Narrowing `ppos` from uint64 to uint32.** Positions index a 23M-base
pseudogenome, so 32 bits is ample and this halves a 6.8 MB array. Measured over
three interleaved reps:

    s36 (uint64)  238,776 / 238,668 / 238,752 KB
    s38 (uint32)  238,800 / 238,692 / 238,892 KB

**No saving. Identical within noise.** Peak RSS occurs during the mapping stage
(232 MB), where the seed index and packed text dominate; `ppos` is allocated by
then but is not what sets the peak. Halving a 6.8 MB array cannot move a peak
that other structures determine. The reasoning was arithmetic on an array size
without checking where the peak actually is.

Neither file is kept. Both were verified byte-identical (`PG_LITERAL 12,506,313`,
851,275 unique) so no size result is affected.

### A measurement-methodology failure worth recording

The first A/B run gave s36 4.56 s and s38 5.29 s, which looked like a clean
refutation. It was not measurable at all: `ps` showed a concurrent `arcs` job at
256% CPU from another session, violating CLAUDE.md rule 3 ("never run two timed
benchmark jobs concurrently"). The same baseline binary then measured 4.56 s and
6.17 s, and round 1 read 1.32 s in one run and 1.91 s in another.

Re-run on an idle machine, interleaved, the timings still overlap:
s36 5.46-6.09 s against s38 5.91-6.02 s, ordering flipping between pairs. **This
machine cannot resolve a 0.2 s effect on a 5 s run.** Memory is stable to ~0.1%
and is the only axis these prototypes can measure reliably at this scale, which
is why the memory numbers above are trusted and the timing ones are not.

Any future speed claim here needs either a much larger input, or a stage-level
timer around the specific change rather than end-to-end wall clock.

### Verdict

Neither loss flips. **But the total does**, because stage 6 is worth more than
both:

| stream | ours | PgRC2 | |
|--------|------|-------|---|
| sequence, coded | 3,158,084 | 3,056,474 | +101,610 |
| order | **2,309,967** | 2,852,758 | **-542,791** |
| **combined** | **5,468,051** | **5,909,232** | **-441,181, 7.5% ahead** |

Final three-axis position, all at 12 threads, order info included on both sides:

| | ours | PgRC2 | |
|---|------|-------|---|
| sequence + order, coded | **5,468,051** | 5,909,232 | **7.5% ahead** |
| time | 4.53-4.75 s | 3.66 s | 1.24x behind |
| peak RSS | 238 MB | 232 MB | 1.03x, parity |

Scope, stated plainly: these are the two largest streams, not whole archives.
Their `-o` archive is 7,063,459 B and these two are 84% of it; the remaining
~1.15 MB is mismatch/offset/length streams this progression computes and
discards. We do not have an archive, so the claim is "ahead on the two streams
both tools produce and that were measured", not "smaller archive".

## Correction: the division was never behind

Several sections above cite a read deficit -- "we admit 593,968 both-side
overlapped reads where they admit 708,819" -- and treat it as a defect to close.
**That comparison was wrong and the deficit is not real.**

Their log reads `HQ reads count: 999340`, then `Found 708819 both-side
overlapped reads`. 999,340 is 1,000,000 minus 660 N-reads: **they do not
deduplicate.** `initAndFindDuplicates`
(`GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp:110-118`) links each
duplicate to the next with `setReadSuccessor(..., maxReadLength)` -- a
full-length overlap, costing zero pg bytes but counting as an overlap.

So their rate is 708,819 / 999,340 = **70.9%**; ours is 593,968 / 851,275 =
**69.8%**. The divisions are equivalent. The 115k gap was an artifact of
comparing a deduplicated population against a non-deduplicated one.

This also corrects the per-read figure below. Their 708,819 carries the 148,065
duplicates our dedup removed, and over-represents them, since a duplicate is
always both-side-overlapped. Their HQ set holds at most ~603,700 distinct reads:

    theirs:  <=603,700 unique -> 21,104,500 pg  =  >=34.96 B/read
    ours:     593,968 unique -> 21,619,175 pg  =    36.40 B/read

**~4% worse per read, not 17%.** The 17% figure divided their pg by their
duplicate-inflated read count.

Which settles solvability. The remaining 1.83% is not a missing stage, a
mistuned parameter, or a weak division -- all three have been checked and
closed. It is ~4% of mean overlap inside a greedy layout already doing the right
thing, against a competitor doing the same thing slightly better. Closing it
needs a better layout algorithm, and there is no cheap lever left.

## Is the remaining gap solvable? A refuted shortcut, and where it actually lives

At seed 16 / stride 22 the gap decomposes as:

| | ours | PgRC2 | gap | share |
|---|------|-------|-----|-------|
| main pg after MEM | 11,632,777 | 11,463,320 | +169,457 | **73%** |
| survivor pg after MEM | 1,294,148 | 1,231,855 | +62,293 | 27% |

**The gap has inverted.** At stage 17 it was 94% survivor; it is now 73% main
pg -- the one thing stages 16-21 never touched.

Sweeping the round-1 classifier gives a clean monotone trade:

| SW | main pg after MEM | survivors | survivor pg | total |
|----|-------------------|-----------|-------------|-------|
| 32 | **11,384,786** (78,534 BELOW theirs) | 14,972 | 1,567,625 | 12,952,411 |
| 24 | 11,511,100 | 13,510 | 1,431,590 | 12,942,690 |
| 16 | 11,632,777 | 12,065 | **1,294,148** | **12,926,925** |

The tempting arithmetic: SW=32's main pg plus SW=16's survivor pg is
**12,678,934**, below PgRC2's 12,695,175. Both halves of a win exist, at
different settings.

**Tested, and it does not work.** The hypothesis was that `SW` conflated two
roles -- the rolling seed width and round 1's classification floor -- in the
same way stage 17 found acceptance conflated with seed geometry. Splitting them
(`R1MINOV` separate from `SEEDW`) reproduces the old numbers exactly: R1MINOV=32
with seed 16 gives main 11,384,786, survivors 14,972, literal 12,952,411 --
identical to SW=32. The split is a genuine no-op.

The reason is conservation, and it is obvious in hindsight: **the reads a
stricter classifier excludes from the main pg are exactly the reads that become
survivors.** Tightening one pool fills the other. The two best halves are two
different partitions of the same 851,275 reads and cannot be held at once.

So the remaining 1.83% is not reachable by tuning. Where it actually is:

- Our MEM removal is already **better** than theirs on both pools (main 46.2% vs
  45.8%, survivor 19.9% vs 18.5%).
- Our per-read packing is worse. At SW=32 we admit ~594k reads and build a
  20,665,755 pg (34.8 B/read); they admit 708,819 and build 21,104,500
  (29.8 B/read). They fit **17% more sequence per byte**.

That is mean overlap, measured earlier at 119.9 against their 133.7 (ideal
tiling at this coverage ~136). Closing it needs a better overlap search -- more
candidates considered per read, or a non-greedy layout -- not another parameter.
Real work, and honestly assessed as the only remaining lever on size.

### Bug found while building it

The first version extended every candidate match before checking that the
source preceded the destination. In a self-match the seed at position `qp` also
occurs at `qp` itself, so the extension ran to the end of the text — O(n) work
at ~2M sampled positions. It had not finished after 5 minutes. PgRC2 documents
the fix as capping usable length at `i - src`; capping before extending rather
than rejecting after brings the whole MEM stage to 12.0 s.

Best: **27.49M, 11 s, 421 MB** against their 22.60M, 4.0 s, 232 MB.

## Stage coverage against PgRC2's own 7 stages

PgRC2 names its stages in its `-B`/`-E` help text (`method_c/PgRC.cpp:215`):
`1:QualDivision 2:PgGenDivision 3:Pg(HQ) 4:ReadsMatching 5:Pg(LQ&N)
6:OrderInfo 7:PgSequences`. Mapping this progression onto them:

| # | their stage | ours | status |
|---|-------------|------|--------|
| 1 | QualDivision | 02, tolerance replaces it | mastered (theirs is inactive by default) |
| 2 | PgGenDivision | 09 -> 12, structural, from round-1 chain topology | mastered -- the hardest one |
| 3 | Pg(HQ) | main pg 20.15M vs their 21.10M | **exceeded** |
| 4 | ReadsMatching | 03, q-gram/pigeonhole lemma | mastered |
| 5 | Pg(LQ&N) | 11, second pseudogenome | **partial -- this is the entire gap** |
| 6 | OrderInfo | -- | not attempted |
| 7 | PgSequences | -- | not attempted |

Mastered through stage 4, exceeded at stage 3, stalled at stage 5. Stages 6-7
were never attempted by design: this progression measures pg length in bases,
so it has no serialization stage at all.

The arithmetic agrees that stage 5 is the whole story. Splitting both totals
into main pg + remainder:

| | main pg | remainder | total |
|---|---------|-----------|-------|
| PgRC2 | 21.10M | 1.50M | 22.60M |
| ours (15) | 20.15M | 7.34M | 27.49M |
| delta | **-0.95M (we win)** | **+5.84M (we lose)** | +4.89M |

The remainder gap (5.84M) is LARGER than the total gap (4.89M): 100% of the
deficit is in the survivor tail, and stage 5 is the survivor tail. Nothing
upstream of it has anything left to win.

Note 5.84M / 141k unplaced reads is ~41 bytes/read, but a raw survivor is ~150
bytes and a chained one ~30 -- so those reads are NOT sitting raw. Stage 11 is
already chaining most of them; it is simply about half as effective per read as
the main pg is. Anyone resuming this should test that before designing a fix.

## Why this stopped here

Measured 2026-08-28 on the same yeast_sub.fq, ARCS's shipping Method B
assembler (`ARCS_VODBG_TIMING=1`, log in `../scope/results/`):

    [HQ-TRIAL] selected frac=0.45 pg_len=16345321 (3 trials)

| assembler | pg_len |
|-----------|--------|
| ARCS Method B (shipping default) | **16.35M** |
| PgRC2 | 22.60M |
| this progression, stage 15 | 27.49M |
| ARCS Method A | 30.15M |

Method B is 27.7% shorter than PgRC2's pg and 40% shorter than stage 15's.
Closing stage 5 here would be catching up to a target ARCS had already passed.
So this stops at stage 15 deliberately -- not abandoned, finished. What it was
built to settle, it settled, and the four transfers below are already in ARCS.

Caveat for anyone quoting the table above: pg_len alone is not a fair
scoreboard. ARCS buys a shorter pg by tolerating mismatches and paying for them
in aux_blob (687 KB here), which PgRC2's pg figure does not carry. The
defensible cross-tool number is total_seq = 4,495,972 B against PgRC2's
4,898,620 B whole archive (see ../PGRC2_COMPARISON.md).

## What this settled

Their speed and memory ARE reachable from textbook parts — same class on both,
with none of their code. Assembly quality was not fully reached (21.6% short),
though the main pseudogenome alone beats theirs; the deficit is that ours holds
141k fewer reads, so more survivors fail to place.

Every gain above came from fixing a defect in this code, not from adding
something PgRC2 has and we lacked:
  - a fixed 32-base seed silently ceilings mismatch tolerance at 3
  - a missing admit[] check emitted excluded reads twice (5 MB)
  - the survivor pseudogenome stage was absent entirely
  - the round-1 classifier was needlessly conservative
  - O(n·L²) re-hashing, not the algorithm, was the 40 s
  - MINOV had been swept in a different architecture and the conclusion carried over

## What transferred to ARCS

  - packed 32-base seeds as integer keys (hash_apsp.cpp)
  - index one side, stream the other
  - outward offset walk: first verified hit is the longest overlap
  - exact-size CSR table, which also fixes silent k-mer loss on bucket overflow

## Findings worth keeping

  - PgRC2's quality-based HQ/LQ split is INACTIVE by default; the division that
    does the work is structural, computed from round-1 chain topology
  - "weak merges cost bytes" is CONDITIONAL: true before the division, inverted
    after it. Once reads pass a tiling filter a marginal merge converts a
    150-byte survivor into a ~30-byte chain member
  - their read-to-pg matching is the q-gram/pigeonhole lemma, textbook
  - reverse complement in the greedy is refuted twice, in two architectures

Build any of these standalone:
    g++ -O3 -march=native -o scs <file>.cpp
    ./scs reads.fq [maxmm] [minov] [seedwidth]
