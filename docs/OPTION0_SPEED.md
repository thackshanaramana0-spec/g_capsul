# Option 0 speed work — calling from a stored archive

## Clean baseline (single job, nothing else on the box)

    wall            172.34 s      peak 5.40 GB
      decode+export+pack   ~72 s
      caller              100.23 s
        kc_H_build          23.45 s
        traversal           62.19 s

For comparison, A (call at compress time) has a caller of 105.98 s, so
**Option 0's caller is already slightly FASTER than A's.** There is no caller
regression on the archive path.

## A measurement error to not repeat

An earlier baseline of 254.92 s was **contaminated**: `option0.sh` ran its own
call step and `option0b.sh` waited only for `best106` to exit, not for that
call, so two full decode+call jobs ran concurrently on 12 cores (their timing
files finished 2 seconds apart). CLAUDE.md rule 3 exists for exactly this.

Everything derived from that number was wrong, including a "1.7x caller
regression" and an "82 s huge-page win". Both were artifacts of contention.
**Do not cite 254.92 s.**

## Lever 1 — huge pages for `kc`

`kc` is ~2.1 GB (140,719,632 x 16 B) and the bubble traversal BINARY-SEARCHES
it, so that walk is dTLB-bound. This box runs
`/sys/kernel/mm/transparent_hugepage/enabled = [madvise]`, meaning a region gets
huge pages only if it asks — and nothing asked.

Two traps hit while implementing it, both already in this repo's history:

1. **madvise after fill is inert.** The hint must land between `reserve()` and
   the fill, on untouched memory.
2. **There are THREE `kc.reserve()` sites** (in-RAM, spill/superkmer, k-way
   merge) and at full scale the SPILL path runs. Hinting one site did nothing
   and printed nothing, which is how the first attempt looked like an 82 s win
   when it had in fact never executed.

Now hinted at all three, and it prints on FAILURE as well as success, so it can
never be silently inert again.

## Lever 2 — the decode round-trip

`capsule_call_from_archive` currently does:

    decode  -> builds `flat`, every read contiguous, IN MEMORY, filled by 12 threads
            -> fwrite 1.86 GB to disk
    call    -> std::getline it back, 12.6M times, into a temporary std::string
            -> capspack::pack_seq that temporary into a second std::string

That is 2-bit -> text -> disk -> text -> 2-bit, twice over (sequence and
quality), and ~25M heap allocations. Range overloads
(`pack_seq(const char*, size_t)`, `pack_qual(const char*, size_t, int)`) are in
`include/caps_pack.h` so the temporary can be dropped; the `std::string`
entry points forward to them, so there is still exactly one definition of the
format.

## Results so far (each solo, each output-verified)

    172.34 s   clean baseline
    162.58 s   + kc huge pages          traversal 62.19 -> 55.69 s (-10.5%)
    158.42 s   + in-memory read handoff ~0 s; keeps a 1.86 GB write saving
    147.27 s   + parallel packing       14.56 -> 1.41 s (10.3x)

Phase budget at 147.27 s:

    decode reads+qual   48.12 s     (reads 8.50 + QUALITY 35.25)
    export pg            6.45 s
    read back + pack     1.41 s
    caller              90.72 s     (traversal ~53 + kc_H ~23)

## Is the quality filter load-bearing? YES — measured

The coherence pass skips sub-Q20 bases when counting per-path read support, and
that is the ONLY thing it does that `kc` cannot. If it were not load-bearing the
35.25 s quality decode could be skipped outright. Measured with
`CAPS_CALL_NOQUAL=1`:

| | with quality | without |
|---|---|---|
| wall | 147.27 s | 121.73 s |
| TP | 36,855 | 36,866 |
| FP | 2,656 | **9,067** |
| precision | 0.9328 | 0.8026 |
| **F1** | **0.8766** | **0.8146** |

**Quality contributes precision only, never sensitivity** — recall is unchanged
(within 11 calls) while false positives TRIPLE. It removes 6,411 false calls and
finds nothing new. Not tradeable for 25 s; `CAPS_CALL_NOQUAL` stays off.

So the quality decode is structurally required, and the remaining lever is to
make it FASTER rather than skip it (below).

## Lever 3 — do not decode the reads at all (analysed, not built)

The read set's k-mers ARE the pseudogenome's k-mers, with multiplicity given by
read coverage, plus corrections where reads differ from pg. Measured on this
input:

    reads with mismatches   1,767,182 of 12,604,917   (86% are EXACT pg substrings)
    total mismatches        3,801,365                 (2.15 per mismatched read)

So a pg-derived build costs ~65 M pg k-mers + 3.8 M x 31 = 118 M corrections =
**~183 M operations against 1.86 G from scanning reads — roughly 10x less.**

**REFUTED — the decode cannot be eliminated.** The bubble coherence pass
(`DBG-COH`, `sup1`/`sup2`) builds a probe map of bubble k-mers and then SCANS
THE READS, using their placements and quality bitmaps to tell a real
heterozygous site from a repeat collapse. `caps_caller.h:3049` states why this
cannot move to k-mer space:

> "k-mer counts cannot see this distinction by construction: the collapsed
> copies are identical in k-mer space, which is why they collapsed. Placement
> is the only signal that separates them"

So a pg-derived build would still produce a correct `kc` 10x cheaper, but the
reads are needed afterwards regardless. Removing the decode would change
results, which is out of scope. **Lever 3 is closed**: the remaining wins are
in HOW the reads are decoded and handed over (Lever 2), not in avoiding it.

Prior art (searched): existing k-mer work compresses k-mer SETS, and read
compressors are noted to "require the substantial overhead of running a k-mer
counter as part of decompression" — i.e. this overhead is acknowledged in the
literature and not solved for calling.

### Correction to the REASON this was refuted

The first justification given here was wrong twice, and the conclusion survived
only by luck:

- "placement is the only signal that separates a het site from a repeat" —
  placement dispersion is INSTRUMENTATION (`EXCL`, i.e. `CAPS_DBG_EXCL`, is off
  by default). It gates nothing.
- "k-mer counts cannot substitute" — the code says the OPPOSITE: "a canonical
  k-mer's count in `kc` IS the number of reads containing it, so walking the
  path and taking the minimum gives per-path read support directly".

The ACTUAL blocker is quality: the coherence scan skips bases below Q20
(`qs` bitmap test), and `kc` carries no quality. That is the one thing the read
scan provides which k-mer counts cannot — and the measurement above shows it is
worth 0.062 F1, so it cannot be dropped.

Independently, DiscoSNP++ runs `kissreads2` to "map back the reads ... in order
to determine the read coverage per allele", so read-scanning is inherent to this
method family, not a quirk of this implementation.

## Lever 4 — parallel quality decode straight to bitmaps

`qlc::decode_to_file` is a SERIAL loop over blocks that are fully independent:
each carries its own length, count and qmin, is `fqz_decompress`-ed on its own,
and nothing crosses a block boundary. 35.25 s on one core of twelve.

It also decodes to the wrong representation: the caller reduces every quality
string to ONE BIT per base, so materialising 1.86 GB of text, writing it to disk
and reading it back to extract those bits is pure overhead.
`qlc::decode_to_bitmaps` walks the index once to find each block's offset and
first read, then decompresses blocks concurrently, writing packed bitmaps
directly into disjoint ranges (no synchronisation).

**MEASURED: 35.25 s -> 5.64 s, 6.25x.** F1 unchanged at 0.8766
(TP 36,855 / FP 2,655 / FN 7,720 -- the published figures exactly).

## FINAL: 172.34 s -> 117.51 s (-31.8%), output identical throughout

| step | wall | change |
|---|---|---|
| clean baseline | 172.34 s | |
| + kc huge pages | 162.58 s | traversal 62.19 -> 55.69 s |
| + in-memory handoff | 158.42 s | ~0 s (keeps: no 1.86 GB write) |
| + parallel packing | 147.27 s | 14.56 -> 1.41 s (10.3x) |
| **+ parallel quality bitmaps** | **117.51 s** | **35.25 -> 5.64 s (6.25x)** |

Final phase budget:

    decode reads+qual   18.14 s   (reads 8.20 + quality 5.64)
    export pg            6.41 s
    read back + pack     1.12 s
    caller              91.33 s   <- now 78% of the run
                       --------
                       117.51 s   peak 5.95 GB

Non-caller overhead fell from ~72 s to ~26 s. Everything further is inside the
caller itself (traversal ~53 s, kc_H_build ~23 s), which is shared with
configuration A and therefore governed by Claim 2's own tuning, not by the
archive path.

## Deep pass — the caller itself

With the non-caller phases down to ~26 s, the remaining budget was the caller:
traversal ~53 s + kc_H_build ~23 s.

### Lever 5 — interleave the four successor probes (KEPT)

`succs()` does FOUR independent `kc_find` binary searches, each ~27 hops over a
2.1 GB table. Run back to back, every hop's cache miss is waited out in full;
the four keys are known before any search starts, so the four latencies can
overlap instead.

Checked it was LATENCY-bound before writing it (this project has been burned by
a filter that won in a latency-bound loop and lost in a bandwidth-bound one):
a model of the access pattern predicts ~113 s of pure latency against 53.5 s
observed -- same order, i.e. latency-dominated -- and a standalone benchmark on
a table of identical shape measured 2.123 s sequential vs 1.393 s interleaved
(1.52x).

**MEASURED IN SITU: traversal 53.49 -> 48.21 s (-9.9%), output identical.**
Well below the 1.52x microbenchmark, because the real loop's `cnt < MINC`
early-exit already skips most nodes and the surviving probes have better
locality than synthetic random keys.

### Levers measured and REJECTED

- **Split `KC{kmer,cnt}` into parallel key/count arrays.** The search touches
  only keys, so the working set halves 2.10 -> 1.05 GB. Benchmarked **1.09x**:
  at 27 hops the first ~15 are cold regardless of element size, because
  consecutive midpoints are gigabytes apart. Not worth a ~30-call-site refactor
  for speed. Remains a 0.52 GB RAM option if RAM ever binds.
- **Drop `cnt==1` singletons from kc.** Audited all six consumers and every one
  already requires `cnt >= 2` (`cnt_max` can't be a singleton; the histogram is
  gated on `>= 2`; `kcount()` feeds a `> 1.5*H` test; `PLMIN = max(2, H/4)`; the
  dBG probes use MINC=2), so it was provably output-identical and predicted
  140.7M -> ~56M entries. **MEASURED: a NO-OP -- kc stayed at exactly
  140,719,632 nodes.** The "60.4% singletons" figure in the MINC comment is a
  WINDOW-scale measurement; at full chr20 the superkmer path has already summed
  duplicate keys before any filter sees them. Reverted. Lesson recorded at the
  site: measure the fraction AT THE SCALE YOU RUN.
- **Bypassing the k-mer spill** (`CAPS_MAXRAM_MB`). Not pursued. The absolute
  5000 MB ceiling is what makes full chr20 run at ~6 GB instead of 28.3 GB;
  the default stays untouched.

## FINAL, ceiling and spill unchanged

    172.34 s / 5.40 GB  ->  111.53 s / 6.21 GB
    F1 0.8766, TP 36,855 / FN 7,720 -- identical to published
    (FP 2,654 vs 2,655: the caller's known bubble-ordering nondeterminism)

    decode reads+qual   17.9 s   (reads 8.2 + quality 5.5)
    export pg            6.8 s
    read back + pack     1.1 s
    caller              84-90 s  (traversal ~48, kc_H ~23)

**-35% wall.** Non-caller overhead fell from ~72 s to ~26 s. What remains is the
caller, shared with configuration A.

## What this pass got wrong, recorded so it is not repeated

1. **A contaminated baseline.** 254.92 s came from two of my own jobs running
   concurrently (CLAUDE.md rule 3). Everything derived from it -- a "1.7x caller
   regression", an "82 s huge-page win" -- was an artifact. Option 0's caller is
   in fact FASTER than A's (91 s vs 105.98 s).
2. **A silently inert optimisation.** The huge-page hint went to one of three
   `kc.reserve()` sites; at full scale the SPILL path allocates. It now hints at
   all three and logs on FAILURE, so it cannot go quiet again.
3. **A stale marker across sessions.** A follow-up script waited on a `_DONE`
   marker in a filename that already existed from an EARLIER session, so it
   fired immediately and ran a second job alongside the first -- the same
   contention error as (1), from reusing filenames. Use fresh names per session.
4. **An unmeasured assumption.** Removing a 1.86 GB file write and 12.6M
   getline calls bought ZERO -- the file was page-cache resident and the real
   cost was 12.6M serial string allocations. The phase timers that revealed this
   are now permanent, so phase costs need never be inferred by subtraction.
