# Handoff prompt — CAPSULE / G_CAPSUL, COMPACT axis

Copy everything below the line into the receiving model, together with the repo
and the exported session transcript.

---

## Who you are and what this is

You are a research engineer taking over an existing, working FASTQ compressor
called **CAPSULE** (repo name `g_capsul`, branch `c_star_pg_advance`). It
compresses losslessly, and it also calls variants and answers coordinate
queries directly from its own archive. Target venue: Nature Methods.

The lineage is **PgRC2** (Kowalski & Grabowski) — pseudogenome-based read
compression. Method, in order: greedy suffix-prefix overlap chaining builds a
pseudogenome from well-tiling reads → remaining reads are pigeonhole-mapped
onto it → unmapped reads are appended and assembled as a second region → the
pseudogenome is self-matched to remove redundancy → everything is emitted as
separate streams and entropy-coded.

**Your mandate is the COMPACT axis: archive size, wall time, peak RAM.**

### The ambition, stated plainly

PgRC2 beats SPRING and Genozip on **all three axes simultaneously** — smaller
*and* faster *and* lighter. That is the existence proof that a three-way win is
achievable, and it is the standing instruction of this project: **never
conclude "you cannot win on all three axes."** If you find yourself about to
write that sentence, you have stopped looking.

### What you may and may not change

- **You MAY restructure substantially.** Replace algorithms, replace data
  structures, redesign representations, fuse or delete stages, change the
  dataflow, make things streaming or incremental. Structural change is wanted.
- **You MAY NOT rewrite from scratch.** This code encodes years of measured
  decisions, many of them counter-intuitive and several of them recorded as
  refutations of the obvious approach. Throwing it away discards the
  measurements, not just the code.
- The PgRC2-style architecture is the **starting point, not a constraint**. If
  a component is genuinely limiting, say so and name what should replace it.
  But the burden is a measurement, not an intuition.

---

## THE MOST IMPORTANT INSTRUCTION: do not inherit the conclusions

You are being given a long transcript of the previous engineer's work. It is
full of statements like *"this stream is at its entropy bound"*, *"this is
structurally impossible"*, *"the curve is not unimodal so no search applies"*,
*"architecture is not the gap."*

**Treat every negative conclusion in that transcript as an untested hypothesis,
not as a finding.** Use the transcript for its *measurements, instrumentation,
invocation recipes and traps* — those are gold and will save you days. Do not
use it for its verdicts.

This is not politeness. It is the empirical record of that very project:

- Its own standing note says *"every floor claim in this project was wrong
  within the hour"* — "100 s impossible" became 55 s; "kc is at its floor"
  had 4 bytes of padding in it.
- On 2026-09-08 the author profiled, declared the MEM matcher's cost
  intrinsic, and then found a **54% reduction** in the same stage an hour
  later — a discriminator the codebase had *already implemented elsewhere and
  never carried across*.
- A code comment asserting a parameter was `0` "on every locked dataset" was
  **false for the human dataset** (it is 2). Building on that comment instead
  of the runtime value would have produced a silently wrong archive.

So: **verify parameter values at runtime, not from comments. Re-measure before
believing any bound.** When you disagree with the transcript, the transcript is
the thing that needs re-testing.

Equally: do not invert this into contrarianism. Several refutations in there
are real and were paid for in machine-days (see "already tested and REFUTED" in
`CLAUDE.md`). Re-run them cheaply before redoing them expensively.

---

## Orient yourself in the code

Roughly 15,800 lines of C++ across seven files. Read in this order.

| file | lines | what it is |
|---|---|---|
| `stages/106_inprocess.cpp` | 4,863 | **the encoder.** The whole pipeline. Start here. |
| `stages/capsule_decode.cpp` | 1,243 | the decoder + `export` / `coverage` / `query` / `call` |
| `include/coders_inproc.h` | 931 | stream coders, codec selector, transforms |
| `include/seqpar_core.h` | 357 | the DNA coder, shared by both paths so they cannot diverge |
| `include/names_coder.h` | 954 | read-name coder |
| `include/quality_coder.h` | 561 | quality wrapper around vendored fqzcomp |
| `include/caps_caller.h` | 6,908 | the variant caller (**not your axis** — do not touch) |

### Stage map of the encoder — real line numbers, use these as entry points

    :946   lap("load+filter+dedup")        parse, filter, optional dedup
    :1027  lap("prefix seed index")        builds pent/ptab/pext
    :1551  lap("round 1 (division)")
    :2005  lap("round 2 (assembly)")       descending-length overlap sweep
    :2083  lap("emit chains")
    :2852  lap("pigeonhole mapping")       maps leftover reads onto the pg
    :3504  [MEM] MINMEM/seed/step          pg self-match setup
    :4095  "both run() passes done"        end of MEM matching
    :4150  lap("pg MEM matching")
    :4565..:4812                            stream prep, then the coding pool

The encoder prints per-stage timings, stream sizes and funnel counters to
stderr on every run. **Read that log before profiling anything** — it already
answers many questions.

### Build and run

    scripts/build106.sh /tmp/best106            # encoder  (must keep -fopenmp)
    scripts/build_decode.sh /tmp/capsule_decode # decoder
    scripts/run_tests.sh                        # 15 self-contained assertions

    INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 \
        bash scripts/encode_adaptive.sh         # the ONLY correct way to encode

**Never invoke the encoder binary directly.** It needs a specific env and
argument vector; running it bare produces an incomplete archive that fails
later in a way that looks like a code bug. Five separate "the tool is broken"
diagnoses in this project were all wrong invocations.

Benchmarks (one dataset, all tools, same methodology as the paper):

    bash scripts/benchmark_0_preflight.sh                    # must say GO
    CLAIMS=1 bash scripts/sanity_archive_one.sh SRR2584863   # COMPACT only
    bash scripts/benchmark_1_run.sh                          # full sweep, 6-9 h

Data lives in `/data/fastq`. The locked dataset list is
`NEW_DATASET_LOCKED.md` (**not** the 17-accession list in `DATASET_LOCKED.md`).
Do not substitute datasets.

---

## The measurements you are inheriting (evidence, not verdicts)

Machine: AMD EPYC 9555, 12 vCPU, L1d 768 KiB, L2 6 MiB, **L3 192 MiB**.
HG002 chr20, 12.6 M reads, 3.99 GB input. `perf record -F 299`, 633 k samples,
`-g` build. Full detail in `docs/COMPACT_HEADROOM.md`.

Wall-clock decomposition (baseline 240.73 s, now 211.04 s):

| stage | s | share | threads |
|---|---|---|---|
| pigeonhole mapping | 64.4 | 27% | 6 |
| pg MEM self-match | 47.4 → **21.8** | 20% → 9% | 6 |
| round 1 (division) | 43.2 | 18% | 12 |
| round 2 (assembly) | 27.2 | 11% | 6 |
| load+filter+dedup | 22.5 | 9% | 12 |
| coding pool (longest job `pos_abs` 19.5 s) | 19.5 | 8% | 6 |
| emit chains | 6.4 | 3% | 6 |

Head-to-head, all LOSSLESS, one tool at a time on an idle box:

| dataset | ours | SPRING | Genozip | time vs SPRING |
|---|---|---|---|---|
| E. coli 0.70 GB | 68,429,027 B | 74,086,400 | 119,614,774 | 1.20x |
| L. major 1.66 GB | 106,248,319 B | 117,616,640 | — | 2.99x |
| HG002 4.00 GB | 573,767,964 B | 598,169,600 | 950,996,735 | 3.90x |

Our throughput **degrades with genome size and repeat content** (60 → 38 → 19
MB/s) while SPRING's stays roughly flat (~70–113 MB/s). On E. coli we are also
*lighter* than SPRING (2.05 vs 2.20 GB peak). The size margin is **thinnest
exactly where we are slowest** (human: −4.08%).

---

## Open levers, with the numbers behind them

Ranked by expected impact × confidence ÷ risk. These are leads, not orders —
re-measure and re-rank.

1. **The candidate sweep forks 2 candidates at 6 threads each**, so ~70% of the
   run uses half the machine. On HG002 the winning candidate beat the loser by
   **99,117 B = 0.0173%**. This is a size/speed policy decision that has never
   been priced properly. Ask whether the choice can be *predicted from a
   measured property of the input* instead of searched by running the pipeline
   twice.
2. **Mapping is repeat-bound.** Instrumented (`-DMAPFUNNEL`, then
   `CAPS_MAPDBG=1`): 43.3 M probes walk **3.66 billion** candidates, 84.6 per
   probe. Buckets at the `seedcap` ceiling of 1024 are **3.05% of probes and
   36.89% of the work**; 5% of probes cause 55%. The current answer truncates
   over-frequent keys at 1024, discarding 126,723 entries (0.828%) **in index
   order, i.e. arbitrarily**. LAST-style adaptive seeds (Kiełbasa et al.,
   *Genome Research* 2011) discard an entry only when it actually disagrees.
   Note: **not output-preserving** — its gate is archive size across the locked
   set, not byte-identity.
3. **Round 1 (43 s, 18%) and load (22 s, 9%) have never been decomposed at line
   level.** Genuinely unexplored. Nobody has looked.
4. **Two documented, unclosed SIZE gaps**, both written in the code as open:
   the seedcap's arbitrary 0.828% loss (above), and MEM matches averaging
   **179 bases against copMEM's 250** (`106_inprocess.cpp:3507`; the seedcap is at `:2514`). copMEM is by
   PgRC2's own authors and is the direct comparable.
5. **Cross-layer opportunities.** The biggest win of the last session was
   noticing that three sibling seed-and-verify matchers (round 2, mapping, MEM)
   had drifted to three different levels of care, and one of them was missing a
   discriminator the other already had. Look for that shape: the same idea
   applied unevenly, the same array built twice, the same text represented two
   different ways. Prefer one change that removes several costs at once over
   isolated 2–3% wins.

---

## Hard gates — non-negotiable

The encoder and decoder are **coupled**. Change a stream in
`106_inprocess.cpp` without matching `capsule_decode.cpp` and you get an
archive that decodes to the wrong data *without erroring*. Four silent
data-loss bugs have already reached "verified" state this way.

- **Output-preserving change** → the archive must be **byte-identical**
  (`cmp`) on E. coli plus at least one large repetitive dataset.
- **Output-changing change** → must not regress any locked dataset, must
  improve the aggregate, and must be **re-verified LOSSLESS by decoding the
  archive** (not the intermediate dumps — the dumps bypass the entropy layer,
  which is where one whole class of bug lived).
- **Any `LOSSY` result halts everything.** It is a bug, not a tradeoff.
- `scripts/run_tests.sh` must stay 15/15.
- **Never run two timed jobs concurrently** — it contaminates wall time and
  peak RAM. Check the box is idle before every measurement.
- Assert ≥40 GB free disk before any A/B. A full disk once produced a *fake*
  regression by silently truncating a spill file.

## Standing rules of this project

1. **A fix must be a formula over a measured property of the input** — never a
   fitted constant, never a per-dataset special case. A parameter sweep is
   diagnosis, not a fix. If you cannot state the change as a formula or as a
   decomposition of a structure, it does not ship.
2. **A change that fails its gate is reverted and recorded, not tuned until it
   passes.**
3. **Retractions stay in the docs, marked in place, never deleted.**
4. **An instrument must not cost what it measures.** A runtime counter in the
   hottest loop measured ~1.8% — the same order as the effects being hunted.
   Compile-time gate it.
5. **Verify the gate itself.** A comparison has twice reported "DIFFERS" for
   byte-identical output because it was comparing files that did not exist.
   When a gate fires, confirm the files exist and are non-empty before acting.

## Method

    measure -> explain -> hypothesise -> smallest meaningful change ->
    benchmark against an untouched baseline, same data, same threads ->
    inspect -> keep or revert -> re-profile -> attack the new largest cost

Profile before guessing. Decompose a stage's *parts*, not its total — a "16.3 s
stage" here hid a 0.10 s loop that got optimised twice by mistake. When you hit
a dominant cost, search the literature for **that specific operation** (MEM
finding, all-pairs suffix-prefix, approximate matching under high repeat
content), not for "FASTQ compression optimisation" in general.

## What to deliver

For every finding, in this form:

> This exact stage/function/data structure does **X**. It costs **Y**, because
> **Z** (with the measurement). The literature/competitor uses principle
> **A**. In CAPSULE that applies specifically at **file:line**, changing
> **B → C**. Expected effect: __. Second-order consequences elsewhere: __.
> Gate: __.

Not "consider multithreading" or "use a better compressor."

If, after genuinely exhausting the alternatives, a cost is irreducible —
**prove it from measurements**, state what would have to change for it to move,
and go attack the next one.
