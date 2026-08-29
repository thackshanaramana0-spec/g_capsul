# Best case: the method-c reimplementation

The canonical configuration and its measured result.

The wins in this progression are spread across stages 16-41, each file being the
previous one plus a single change. Nothing recorded which combination is *the*
result, so this does: one configuration, one run, every number from that run.

Measured 2026-08-29 on `yeast_sub.fq` (1M reads of *S. cerevisiae*, 851,275
unique after N-filter and dedup), SDC3 server, 12 vCPU.

## The configuration

**Pipeline: `42_containment.cpp`** — the newest file, and therefore the one
carrying every accepted change: prefix-containment removal (42), the arena
release (41), the packed mapping index (40) and RC-only self-matching (36), which
in turn inherit maximal MEM extension (33), lazy parsing (34), best-match
placement (28), the probe filter (26), packed reads (21), and the rest.

    g++ -O3 -march=native -pthread -o best 42_containment.cpp
    FWD_SELF=0 ./best yeast_sub.fq 3 40 16 22 16 16 1 24 64 1

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

| axis | ours | PgRC2 | |
|------|------|-------|---|
| size (four streams) | **5,815,545** | 6,352,312 | **0.92x -- ahead** |
| peak RSS | **189 MB** | 232 MB | **0.81x -- ahead** |
| time @12 threads | ~4.7 s | 3.70 s | 1.27x -- behind |

Two of three axes won.

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

**Time is the weakest number.** Another session runs on this machine and the same
binary has measured 4.56 s and 6.17 s. Use `bench.sh` (min over reps, per stage)
for anything that depends on it; a single wall-clock figure here is indicative
only.

## Reproducing

    cd benchmark/reimpl
    g++ -O3 -march=native -pthread -o best 42_containment.cpp
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

The input is built by `../scope/make_scope_inputs.sh` from
`/data/fastq/DRR976266_1.fq`; PgRC2's numbers come from `method_c/build/PgRC -i
yeast_sub.fq -t 12 out.pgrc` and its own stdout.
