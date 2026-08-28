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
