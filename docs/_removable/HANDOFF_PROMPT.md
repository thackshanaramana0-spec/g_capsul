# Handoff prompt — CAPSULE / G_CAPSUL, COMPACT axis

Copy everything below the line into the receiving model, together with the repo
and the exported session transcript.

**Glance over the transcript at
<https://github.com/thackshanaramana0-spec/g_capsul/blob/c_star_pg_advance/2026-09-08-194956-this-session-is-being-continued-from-a-previous-c.txt>
for context only — it is a working Claude session, not a final or authoritative
document.** Skim it to understand how the system was built, measured and
debugged; refer back to it for specifics, but never treat anything in it as
settled.

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

### The goal

1. **The win must be STRUCTURAL and it must GENERALIZE.** A different
   algorithm, representation, data structure or dataflow — whose gain holds
   across the locked datasets and grows with input size, not one that appears
   on a single favourable file. Tuning a constant is allowed as *part* of a
   change, never as the change itself: a fitted constant is not a result.
2. **Win all three axes at once** — smaller archive *and* less wall time *and*
   less peak RAM than SPRING, Genozip and PgRC2.
3. **The insight must hold.** From the `.capsule` archive **alone** — no FASTQ,
   no reference — Claim 2 (variant calling) and Claim 3 (export, coverage,
   query) must still work. A size or speed win that breaks either is a loss,
   not a trade-off.
4. **A margin worth a Nature Methods paper** — large, reproducible, and
   explained by its mechanism.

### How to work

5. **Additive only; destroy nothing.** Your own branch `gpt2026`, namespaced
   symbols, substantial replacements behind flags that default OFF, and the
   existing path byte-identical while your flag is off. **Never delete, rewrite
   or force-push anything on GitHub** — no history rewriting, no removing
   existing code, streams, scripts or docs. Everything you add must be
   deletable in one clean diff.
6. **Higher-order, out-of-the-box, critical and precise engineering.** Profile
   one level deeper than the last profile. Ask what is computed twice, what is
   represented in the wrong form, which stage could be fused, moved, made
   streaming or deleted outright. Prefer one change that removes several costs
   at once over isolated 2-3% wins.
7. **All three claims already work and run in parallel. Your depth goes to
   Claim 1 (COMPACT).** Do not touch the caller.
8. **Do not fear the previous engineer's dead ends.** That transcript is full
   of refutations, "at its bound" verdicts and abandoned directions. Re-measure
   them; its own author overturned his "this cost is intrinsic" verdict with a
   54% win an hour later. Take its measurements, instrumentation and traps.
   Leave its conclusions.

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

## Orient yourself — by reading the code, not by being told

**Scan the repository yourself and build your own map before you change
anything.** You are deliberately not being given a file-by-file guide or a list
of line numbers. Two reasons, both practical:

- Line numbers move. Several in the transcript are already stale after a single
  session's edits. Anchor on quoted comment text or symbol names, never on
  `:NNNN`.
- A map handed to you imports the previous engineer's idea of what matters, and
  that is exactly the bias you are here to avoid. The largest win of the last
  session came from noticing something the map would not have shown: three
  sibling matchers had drifted to three different levels of care.

Read the encoder end to end before forming hypotheses. Understand every
representation, buffer, transformation, codec invocation, allocation, copy,
thread boundary and temporary — and **why each exists and what it costs**. The
encoder prints per-stage timings, stream sizes and funnel counters to stderr on
every run; read that log early, it answers a lot.

### Operational facts you cannot infer by reading

These are invocation contracts, not code understanding. Getting them wrong
produces failures that look like bugs in the code:

    scripts/build106.sh /tmp/best106            # encoder (the -fopenmp matters)
    scripts/build_decode.sh /tmp/capsule_decode # decoder
    scripts/run_tests.sh                        # must stay 15/15

    INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 \
        bash scripts/encode_adaptive.sh         # the ONLY correct way to encode

**Never invoke the encoder binary bare.** It needs a specific environment and
argument vector; run directly it silently produces an incomplete archive. Five
separate "the tool is broken" diagnoses in this project were all wrong
invocations, not bugs.

    bash scripts/benchmark_0_preflight.sh                    # must print GO
    CLAIMS=1 bash scripts/sanity_archive_one.sh SRR2584863   # COMPACT only
    bash scripts/benchmark_1_run.sh                          # full sweep, 6-9 h

Data is in `/data/fastq`. The locked dataset list is `NEW_DATASET_LOCKED.md`
(**not** the 17-accession list in `DATASET_LOCKED.md`). Do not substitute
datasets.

---

## Isolation: how to work without endangering what exists

Everything you do is **additive and reversible**. Concretely:

1. **Work on your own branch.** Branch from `c_star_pg_advance` and name it
   `gpt2026`. Never commit to `c_star_pg_advance` or `main`. Tag a backup ref
   before any history-altering operation.

       git checkout -b gpt2026

2. **Namespace everything new.** New symbols, structs, files, env vars and
   build flags carry a `gpt2026` marker — `gpt2026_sext`, `GPT2026_ADAPTIVE`,
   `docs/GPT2026_FINDINGS.md`. Anyone reading a diff must be able to see at a
   glance what is yours and delete it cleanly.
3. **Prefer adding a path over editing one in place.** Where a replacement is
   substantial, put it behind a flag (default OFF, so the shipped behaviour is
   untouched), prove it, and only then propose making it the default. The
   existing path must remain runnable and must keep producing its current
   output byte-for-byte while your flag is off.
4. **Do not delete or rename existing code, streams, scripts or docs.** The
   retractions and refuted experiments in this repo are evidence, not clutter;
   its own rules say retractions stay marked in place, never removed.

## Do not break the pseudogenome insight, or the other two claims

This is the constraint most likely to be violated by a purely COMPACT-minded
optimisation, so read it twice.

The archive is **not just a compressed file**. The whole thesis of this project
is that the pseudogenome built during compression *is* a reference-free
assembly and coordinate system, so the same archive also serves:

- **Claim 2** — calling variants directly from the archive, no FASTQ, no
  reference (`capsule_decode call`).
- **Claim 3** — `export`, `coverage` and `query` off the archive alone.

Those paths consume named streams from the container. Verified consumers
include `contig_spans`, `pos_abs`, `pos_sec`, `pos_region`, `pos_strand`,
`orig2uid_flags`, `orig2uid_vals`, plus the read lengths and the assembly
streams (`literal`, `mem_triples` and companions) that rebuild the
pseudogenome. **Derive the authoritative list yourself from
`stages/capsule_decode.cpp`** — do not trust this paragraph as complete.

Therefore:

- A change that drops, renames, reshapes or reorders a stream may shrink the
  archive and **silently break Claim 2 or Claim 3**. That is a regression even
  if every COMPACT number improves.
- A change that destroys the *addressability* of the pseudogenome — its
  coordinate meaning, the read placements, the mapping from archive
  coordinates back to reads — breaks the central claim of the paper. Compression
  ratio does not buy that back.
- Whatever you do to the assembly, **it must remain an assembly**. If your
  change makes the pseudogenome cheaper to produce but no longer a usable
  coordinate system, it is not a win; it converts this project into an
  ordinary compressor and throws away its reason to exist.

**Gate for this:** after any change that touches streams, the container, the
pseudogenome construction or the read placements, run the archive through the
other claims and confirm they still work:

    CLAIMS=1   bash scripts/sanity_archive_one.sh SRR2584863   # COMPACT
    CLAIMS=13  bash scripts/sanity_archive_one.sh SRR2584863   # + export/coverage/query
    bash scripts/sanity_archive_one.sh HG002                   # all three claims

The last one must still produce a populated Claim 2 table (het-SNV F1 ~0.888
from the archive) and a Claim 3 table. Cost of addressability, for scale: the
`contig_spans` stream that makes archive-path calling possible is 232,509 B —
**0.041%** of a 573 MB archive. Do not "optimise" it away.

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
   **179 bases against copMEM's 250** — grep the encoder for the comment
   `copMEM's 250`, and for `[seedcap]` for the other. copMEM is by
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
