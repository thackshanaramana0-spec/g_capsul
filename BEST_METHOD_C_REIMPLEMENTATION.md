# Best case: the method-c reimplementation

The canonical configuration and its measured result.

The wins in this progression are spread across stages 16-45, each file being the
previous one plus a single change. Nothing recorded which combination is *the*
result, so this does: one configuration, one run, every number from that run.

Measured 2026-08-29 on `yeast_sub.fq` (1M reads of *S. cerevisiae*, 851,275
unique after N-filter and dedup), SDC3 server, 12 vCPU.

## The configuration

**Pipeline: `45_sweep_loop.cpp`** — the newest file, and therefore the one
carrying every accepted change: the sweep level-loop rework (45), the containment
guard (44), the flat prefix index (43), prefix-containment removal (42), the arena
release (41), the packed mapping index (40) and RC-only self-matching (36), which
in turn inherit maximal MEM extension (33), lazy parsing (34), best-match
placement (28), the probe filter (26), packed reads (21), and the rest.

    g++ -O3 -march=native -fopenmp -o best 45_sweep_loop.cpp
    FWD_SELF=0 ./best yeast_sub.fq 3 40 16 22 16 16 1 24 64 1

Stages 43-45 changed only how the work is done, never what is computed: every
size-bearing stream is **byte-identical** to stage 41's (`literal.txt`,
`perm.u32`, `mem_triples.bin`, `mem_gaps.bin`, `mem_lens.bin`, `mem_srcs.bin`,
all compared with `cmp`), so the size table below is untouched by them and only
the time and memory rows moved.

Positional arguments: `maxmm=3 minov=40 SW=16 stride=22 SEEDW=16 R1MINOV=16 K2=1
MINMEM=24 MAXCAND=64 LAZY=1`, with `FWD_SELF=0` selecting reverse-complement-only
self-matching. Add `DUMP_LIT=1 DUMP_PERM=1` to emit the stream dumps the coders
below consume.

**Stream coders**, each a separate file so it can be verified alone:

| stream | coder | verification |
|--------|-------|--------------|
| sequence | `35_match_model.cpp` — multi-order CM + match model | round trip on all 12,506,313 bases |
| order | `23_perm_coder.cpp` — Lehmer + Fenwick + range coder | round trip on all 999,340 elements |
| MEM sources | `37_ref_coder.cpp` — range-coded, bounded by destination | round trip on all 58,908 |
| MEM gaps, lengths | varint + `xz -9e` | -- |

## The result

| stream | ours | PgRC2 | delta |
|--------|------|-------|-------|
| sequence (literal, coded) | **2,979,683** | 3,056,474 | **-76,791** |
| MEM references | 283,686 | 177,180 | +106,506 |
| order (permutation) | **2,309,967** | 2,852,758 | **-542,791** |
| mismatch symbols (estimated) | **242,209** | 265,900 | **-23,691** |
| **TOTAL** | **5,815,545** | 6,352,312 | **-536,767, 8.45% ahead** |

References break down as gaps 46,592 + sources 174,346 + lengths 62,748.

The sequence row is the context-mixing coder. Substituting 2-bit + `xz -6 -T12`
gives 3,057,532 B, a total of 5,893,394 and a 7.22% lead, at 1/44th of the coding
time -- see the time discussion below for why that is the configuration to quote
when speed is also being claimed.

**The time axis has to be stated carefully, because the obvious comparison is not
a fair one.** PgRC2's wall clock is its whole binary: assemble *and* entropy-code
the archive. Ours, above, is the assembler alone -- the stream coders are separate
programs and their cost is not in it. Timing our assembler against their whole
pipeline flatters us, and it hid a large cost for most of this session.

Coding cost, measured on the streams this run produces:

| stream | coder | time | output |
|--------|-------|------|--------|
| sequence | `35_match_model.cpp`, CM + match model | **19.01 s** | 2,979,683 B |
| sequence | 2-bit + `xz -6 -T12` | **0.43 s** | 3,057,532 B |
| order | `23_perm_coder.cpp` | 0.18 s | 2,309,967 B |
| MEM refs | `37_ref_coder.cpp` + xz | 0.06 s | 283,686 B |

The context model buys 77,865 B for 18.6 s. That is a bad trade at the pipeline
level and it was invisible while only the assembler was being timed. So there are
two honest configurations, and both are reported:

| axis | ours, size config | ours, speed config | PgRC2 |
|------|-------------------|--------------------|-------|
| size (four streams) | **5,815,545** -- 0.92x | **5,893,394** -- 0.93x | 6,352,312 |
| wall, end to end | 21.9 s -- 6.4x, **lost** | **3.34 s** -- 0.97x | 3.43 s |
| peak RSS | **219 MB** -- 0.93x | **209.6 MB** -- 0.89x | 234.7 MB |

**Size is won decisively either way** -- 7.2% ahead even with the fast coder.
**Memory is won either way.** **Speed is parity, not a win:** 3.34 s against
3.43 s is a 3% margin, inside the spread of these measurements. The defensible
statement is that speed stopped being a deficit -- it was 1.27x behind and is now
level -- not that this is faster than PgRC2.

The speed configuration runs the four coders concurrently (they are independent);
the sequence coder dominates at 0.43 s, so the whole coding step is 0.42 s
measured, on top of the 2.92 s assembler.

Assembler stage alone, which is what stages 43-45 moved:

| axis | ours | PgRC2 (whole binary) |
|------|------|----------------------|
| assembler wall @12 threads | **2.92 s** | 3.43 s |

That row is a stage result and is labelled as one. It is not a like-for-like
comparison and must not be quoted as the headline.

Time and RSS are the minimum of three runs of each binary, taken back to back on
an otherwise idle machine (`load average 0.11`) so the two sides face the same
conditions; PgRC2 is `method_c/build/PgRC -i yeast_sub.fq -t 12` under
`/usr/bin/time -v`. Earlier records in this file quoted ~4.7 s against 3.70 s --
both sides were measured while a second session held the machine, and both
numbers moved when it was quiet. The comparison, not either figure alone, is what
carries over.

## Per-unit rates — ahead on every one

| | ours | PgRC2 |
|---|------|-------|
| sequence | **1.9060 bits/base** | 1.9261 |
| order | **18.49 bits/read** | 22.82 |
| MEM sources | **23.68 bits/match** | 23.83 |
| MEM lengths | **8.52 bits/match** | 8.99 |

The one deficit is not a rate. It is match count: 58,908 against their 43,190, on
a pseudogenome of 23.2M against their 21.1M. The MINMEM sweep confirms that count
is optimal *for our pseudogenome* under two different reference prices -- fewer
matches costs more literal than it saves -- so closing it requires a different
assembler, not better coding.

## Speed: where the 1.40 s came from

This section is about the **assembler stage**, not the pipeline; see the time
discussion above for why those are different numbers. Stages 43-45 are pure
engineering -- the output is byte-identical throughout, so nothing here trades
size for time. Per-stage minimum over 5 reps, `bench.sh`,
both binaries measured in the same quiet window:

| stage | 41 (previous record) | 45 (now) | |
|-------|---------------------|----------|---|
| load + filter + dedup | 0.54 | 0.60 | +0.06 |
| prefix seed index | 0.15 | **0.06** | -0.09 |
| round 1 (division) | 1.19 | **0.55** | **-0.64** |
| round 2 (assembly) | 0.58 | **0.27** | **-0.31** |
| emit chains | 0.22 | 0.22 | -- |
| pigeonhole mapping | 0.74 | **0.59** | -0.15 |
| pg MEM matching | 0.61 | **0.59** | -0.02 |
| **wall clock, min of 3** | **4.32 s** | **2.92 s** | **-1.40 s** |

Three findings, in the order they were worth:

**The level loop, not the search (45).** The sweep runs once per overlap level --
134 levels in round 1, 110 in round 2 -- so anything paid per level is paid 244
times. Three such costs were sitting there, none of them the actual overlap
search. `cand.assign(w*CCAP, NONE)` memset 27 MB at the widest level to write a
sentinel that is never read, because the serial commit reads `cand[i*CCAP+c]`
only for `c < ccnt[i]`. Twelve `std::thread`s were created and joined per level,
2,928 spawns in a run under four seconds. And the per-thread probe counter was a
`vector<size_t>`, putting twelve counters in two cache lines and incrementing
them 20.9M times from twelve cores, so every increment stole the line back.
Growing the buffer once, letting OpenMP hold one pool across all levels, and
accumulating the counter in a register: **round 1 and round 2 together fell from
1.77 s to 0.82 s**, with `links`, `both-sides-overlapped` and `probes` all
unchanged to the digit.

**The prefix index was still a hash map (43).** Stages 18 and 40 had replaced
`unordered_map<K, vector<V>>` elsewhere for large wins; this one survived in the
place it is hit hardest, 20.9M probes across the two rounds. 851k distinct keys
each carried a node header plus a separately allocated vector. Flat instead:
`(key32 << 32) | rid` in one sorted array with an open-addressing table pointing
at each key's first entry. Sorting the whole word sorts by key, because the key
is in the high bits, and rid ascends within a key -- so candidates are visited in
the same order the map's insertion-ordered vector gave, and the output is
unchanged. The 32-bit key is exact at SW 16 and partial above it, which is safe:
a collision only proposes a candidate, and `rcmp` verifies every candidate
against the full L bases anyway.

**A pass proving a negative (44).** Stage 42's containment scan is a full
lexicographic sort of all 851k reads. At constant read length one read can only
contain another by equalling it, which dedup has already handled, so on
fixed-length input the pass is guaranteed to find nothing -- 0.31 s, 8% of the
run. It is now guarded by a length-variation check and runs only when lengths
actually vary. Verified on both paths: fixed-length output still `PG_LITERAL
12,506,313`, variable-length still folds 69,466 contained reads.

## Memory: a deliberate 20 MB

Peak RSS moved 189 MB -> 209.6 MB across these stages, and that was the intended
trade rather than a regression. The flat prefix index is ~15 MB of the increase
and it is what bought most of round 1; the containment machinery is the rest.
Against PgRC2's 234.7 MB the margin narrows from 0.81x to 0.89x, which is the
axis being spent, and it buys 1.40 s on the axis being won. Nothing here is
close to the 3-4x kind of trade that would make the memory claim vanish.

## Variable-length reads

PgRC2 refuses them outright -- `Unsupported variable length reads.` -- so this is
capability it does not have rather than a margin over it. Verified on a realistic
input: 30% of reads trimmed 1-20 bases from the 3' end, as quality trimming
produces, giving 21 distinct lengths from 130 to 150.

| | fixed 150 bp | variable, no containment removal | variable, with it (42) |
|---|--------------|----------------------------------|------------------------|
| both-side overlapped | 69.8% | 40.0% | 40.8% |
| survivors | 12,065 | 57,230 | **53,808** |
| pg length | 23,233,953 | 24,066,661 | **23,972,320** |
| literal | 12,506,313 | 15,564,505 | **15,309,682** |

**Fixed-length output is byte-identical with and without stage 42**
(`PG_LITERAL 12,506,313`), which is the property that matters: containment
subsumes dedup, so at constant read length one read contains another only when
they are equal and the pass degenerates to the existing behaviour. None of the
three-axis results above are affected.

On variable-length input it folds 69,466 contained reads into their containers
and recovers 254,823 bases of literal. The remaining gap to fixed-length is not
this: it is that trimming genuinely destroys 3' overlap, which no amount of
bookkeeping restores. Tier 2 of `VARIABLE_LENGTH_DESIGN.md` (interior
containment, where a read sits inside another rather than at its start) is still
unimplemented and is the next thing to try if variable-length quality matters.

## Scope: what this reimplementation does NOT do

Stated plainly so the three-axis result is not read as more than it is.

**It is not a compressor.** No container format, no quality stream, no read
names, and no decoder for the pipeline as a whole -- though every individual
stream coder round-trips, which is what makes the byte counts meaningful. It
measures assembly and stream coding, nothing else.

**Claim 2 (FAITHFUL -- variant calling) is not implemented here.** The assembler
emits the (contig, position, strand) triples a caller needs -- that is the
`CallData` contract ARCS uses -- but no calling, no F1, no GIAB comparison exists
in this progression. PgRC2 has no variant caller either, so nothing is being
conceded; it is simply out of scope for a size/speed/memory comparison.

**Claim 3 (ADDRESSABLE -- export, coverage, query) is not implemented here.**
Same reasoning: it is an ARCS capability with no PgRC2 counterpart, so it cannot
appear in a head-to-head.

Both are ARCS work, not method-c work, and neither can change the numbers above.
They are recorded here so a reader knows the boundary of what was measured.

## What is and is not claimed

**Measured, deterministic, reproducible:** every size number above, and the peak
RSS. All four coders round-trip.

**Estimated:** the mismatch stream, priced at PgRC2's observed 1.55 bits per
mismatch applied to our 1,250,111 mismatches. Everything else is measured.

**Not a whole-archive comparison.** These are the four streams both tools
produce. PgRC2's `-o` archive is 7,063,459 B and these four are ~90% of it; their
reads-list offsets (683,370 B) and mismatch-count streams (176,241 B) have no
counterpart here, because this progression computes read positions and discards
them. It is not a compressor: there is no container, no quality stream, no names,
no decoder for the pipeline as a whole.

**Time is the most fragile number.** The same binary has measured 4.56 s and
6.17 s when a second session held the machine. Two rules make it usable: measure
the two tools back to back in the same quiet window, and take the minimum of
several runs, since contention only ever adds. `bench.sh` does the per-stage
version of this. Note that its "sum of mins" is optimistic -- the minima do not
co-occur in one run -- so the headline figure above is wall-clock minimum, which
for this pipeline matches the stage sum within 0.01 s (there is no untimed tail).

## Reproducing

    cd benchmark/reimpl
    g++ -O3 -march=native -fopenmp -o best 45_sweep_loop.cpp
    g++ -O3 -march=native -o seqcoder 35_match_model.cpp
    g++ -O3 -march=native -o permcoder 23_perm_coder.cpp
    g++ -O3 -march=native -o refcoder  37_ref_coder.cpp

    FWD_SELF=0 DUMP_LIT=1 DUMP_PERM=1 ./best yeast_sub.fq 3 40 16 22 16 16 1 24 64 1
    ./seqcoder  literal.txt                              # 2,979,683
    ./permcoder perm.u32                                 # 2,309,967
    ./refcoder  mem_triples.bin 23233953 21619175        # 174,346
    xz -9e -c mem_gaps.bin | wc -c                       #    46,592
    xz -9e -c mem_lens.bin | wc -c                       #    62,748

    ./bench.sh 5 ./best yeast_sub.fq 3 40 16 22 16 16 1 24 64 1   # per-stage timing
    /usr/bin/time -v ./best yeast_sub.fq 3 40 16 22 16 16 1 24 64 1   # wall + peak RSS

The input is built by `../scope/make_scope_inputs.sh` from
`/data/fastq/DRR976266_1.fq`; PgRC2's numbers come from `method_c/build/PgRC -i
yeast_sub.fq -t 12 out.pgrc` and its own stdout.
