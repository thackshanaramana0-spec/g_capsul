# Tagging chaining's own branch points — full engineering plan

Read-only analysis, 2026-09-04. No code changed. This turns last night's vague
idea ("transform the trunk into a graph") into a specific, gated, measurable
plan, using the code that actually exists — `stages/106_inprocess.cpp`'s
overlap sweep and `include/caps_caller.h`'s substrate/bubble machinery.

**The one-paragraph version:** the branch information the caller currently
pays 738 s + tens of GB to *rediscover* is computed once already, for free,
inside the sweep that builds the pseudogenome — it is simply thrown away at
commit time. Capturing it costs O(1) extra work at the exact moment it already
exists, and — because Claim 2 runs in-process rather than from a saved
archive — **it does not have to touch the archive at all**, which changes
last night's "cost to Claim 1" framing more than the idea itself does.

---

## 1. Where the decision actually happens — exact lines

`stages/106_inprocess.cpp`, the overlap sweep (`sweep()`, used for both round 1
and round 2; round 2 is the one that matters, see §2):

- **Line 887-909**, parallel phase, per open tail `a` at overlap level `L`:
  finds every read `b` whose first `L` bases exactly match `a`'s last `L`
  bases (`rcmp`), and records **up to `CCAP=8` of them** into
  `cand[i*CCAP+0 .. i*CCAP+ccnt[i]-1]`, in the *order the hash bucket produced
  them*.
- **Line 934-949**, serial commit, same order: walks `cand[i*CCAP+c]` for
  `c = 0, 1, 2, ...` and takes the **first** `b` that is not already claimed
  (`prv[b]==NONE`) and would not close a 2-cycle. That is the entire
  decision rule: **first untaken candidate wins**, not best-scored.
- **Line 823-826** (`open_tails` filter): once `nxt[a]` is set, `a` leaves
  `open_tails` **permanently** — it is never reconsidered at a shorter `L`.
  So the commit above happens **exactly once per read**, not once per level;
  there is one and only one moment where `a`'s candidate list existed and a
  choice was made.

**The exact quantity that defines an ambiguous decision:** `ccnt[i] >= 2` at
the level where `nxt[a]` is actually set (not at every level `a` is probed —
most levels a tail is probed with `ccnt[i]==0` because no candidate passed;
the interesting case is when the winning level had a genuine second option).

## 2. Round 1 vs round 2 — which one to hook

Round 1 (line 998, `sweep_minov=R1MINOV`) only **labels** reads
(`both_sides`/`admit`) — it decides which reads are well-tiled enough to enter
round 2, then its `nxt`/`prv`/`ovl` arrays are **fully rebuilt** in round 2
(line 1229 `lap("round 2 (assembly)")` runs `sweep()` again from scratch over
`admit`-filtered reads). Round 1's candidate ties are therefore not
meaningful — they get thrown away regardless.

**Round 2's commit is the one to hook.** Its `nxt`/`ovl`/`ppos` arrays are what
the chain-emission loop (line 1231-1277) walks to build `pg` and
`g_contig_spans` — the exact structure `CAPS_CALL` receives. A tie captured
here is a tie in the structure Claim 2 actually consumes.

## 3. What the tie actually means, biologically — and the filter it needs

`rcmp(a, off, b1, L)` and `rcmp(a, off, b2, L)` both passing means `b1` and
`b2` **share an identical `L`-base prefix** (the overlap with `a`'s suffix).
They are not required to agree beyond that shared window. So a genuine tie is:
*two reads that agree on the whole anchor length but can differ immediately
past it* — which is structurally identical to what `caps_caller.h`'s
`extract_bubble` (line 153) already calls a bubble: shared opening anchor,
divergence after it, required closing anchor before being trusted
(line 1616). **Chaining has already found the opening anchor and enumerated
both sides of it as a side effect of assembly**, before the caller does any
work at all.

This is exactly the same phenomenon last night's `PG_AS_GRAPH_REFUTED.md`
measured as branch degree, seen from the other side: `PG_AS_GRAPH_REFUTED.md`
found that finished contigs have out-degree 1 because chaining *resolves*
these ties — this plan proposes recording the tie **before** it gets resolved
away.

**The filter that finding demands here too, not optionally:** `ccnt[i]==2` is
the het-bubble shape. `ccnt[i]>=3` is very likely a repeat or paralog family —
the same class `PG_AS_GRAPH_REFUTED.md §3` found enriched 2.9x in false
positives. **Any capture mechanism must gate on `ccnt[i]==2` specifically**,
not "ccnt[i]>=2", or it inherits the same FP class the graph idea couldn't
filter its way out of.

## 4. Where the false-positive risk actually comes from, and how to bound it

Two sources of a *spurious* tie, both addressable without new machinery:

1. **Sequencing errors at short `L`.** `rcmp` is an exact match, so a
   collision at `L` near `Lmax` (long, specific context) is astronomically
   unlikely from error alone. As the sweep descends toward `sweep_minov`/`SW`,
   exact-match collisions by chance become more plausible. **Bound: record
   the commit `L` alongside the tie, and require `L` above some fraction of
   `Lmax` (a swept, formula-derived threshold per standing rule 1 — not a
   fitted constant) before trusting it as a het signal**, mirroring how
   `MEM_MAXMM` and the duplicate-collapse threshold are already each a
   measured function of a property of the run, not a hardcoded number.
2. **Repeat families**, handled by the `ccnt[i]==2` gate in §3.

Neither requires new code to *detect* — both are properties of data already
in hand (`ccnt[i]`, the commit `L`) at the exact moment of capture.

## 5. What to record, and the key realization about where it lives

Minimal record at commit time (line ~941, inside the existing `if` branch):

    struct BranchTag { uint32_t read_a; uint32_t chosen_b; uint32_t rejected_b; uint32_t commit_L; };

No sequence data needs to be copied — `seqs[rejected_b]` is already resident
in memory for the whole run. This is a list of **four integers per tie**,
appended only when `ccnt[i]==2` at commit (§3's gate already excludes the
`>=3` case, so no list is even needed for those).

**Translating to contig coordinates happens for free in the emit-chains
walk.** Line 1263-1269 already computes `ppos[cur]` for every committed read
as it walks the chain to build `pg`. A second, tiny pass over the recorded
`BranchTag` list — after `pg` and `g_contig_spans` exist — looks up
`ppos[read_a]` (already computed, O(1)) to get the exact `(contig_id,
pg_offset)` anchor position. No new placement work; it reuses output that
already exists.

### 5.1 The realization that changes the "cost to Claim 1" framing

This entire mechanism is naturally expressed as **in-memory state gated on
`CAPS_CALL`**, exactly like `g_contig_spans` itself (`static
std::vector<...> g_contig_spans;` gated `if(CAPS_CALL) g_contig_spans.push_
back(...)` throughout the file). `CAPS_CALL` is off for every Claim 1 run —
it is only ever set for the HG002-005 diploid runs, per the architecture this
session already established (Claim 2 calls at compress time, off in-process
state; it does not read a saved archive).

**So the honest engineering answer is not "a few hundred KB added to the
archive" — it is closer to zero archive bytes**, because nothing here needs to
survive past the in-process run that already builds the VCF. The earlier
estimate (§6 below) is the answer to a *different*, harder question — "what if
a future decode-only Claim 2 needs this from a saved archive" — which is not
today's design and should not be assumed into this plan.

## 6. If it ever DID need to be archived (kept for completeness, not proposed)

Only relevant if a future decode-from-archive Claim 2 is built (explicitly out
of scope today). Basis: real ties should track true heterozygosity rate, which
`HETSCAN` already measures per-run (`pair_frac`, e.g. 0.0264 on HG002,
scale-invariant, confirmed at both window and full-chr20 scale tonight) — not
the volume of raw reads. At chr20 scale that is on the order of 45,000-55,000
loci (matching the truth het-SNV+indel count order of magnitude), at an
estimated 3-6 bytes/entry before entropy coding: **roughly 150-350 KB**,
against a 26.4 MB archive (~0.5-1.5%). This number is unchanged from last
night's estimate; it is retained here only as the answer to the harder
question, not as this plan's actual cost.

## 7. What this could replace or cheapen in the caller, and why it matters more than a new accuracy channel

`caps_caller.h` currently spends most of its cost **rediscovering** exactly
this information from scratch:

- `build_substrate`'s second call (`dup=0.92`, line 1397) exists specifically
  to keep haplotypes apart so bubbles are findable — 738 s, called because the
  caller has no other way to know where the two haplotype paths are.
- `pkidx`/`pcount`/`porient` (line 2093) and `rc_reads` (line 2147, the
  15-30 GB structure fixed last night) exist to **re-derive anchor
  uniqueness and candidate read pairs** that chaining already computed and
  discarded.

If the branch-tag list from round 2 already names `(anchor contig position,
chosen read, rejected read)` for every genuine tie, the caller's bubble search
does not need to rebuild `rc_reads` or scan for anchors at all for those
loci — it can seed directly from the tag list and only fall back to the
existing full search for loci chaining never flagged (recall on non-tagged
loci is preserved by keeping the existing path as a fallback, not a
replacement).

**This is the connection worth stating plainly:** the RAM/speed skeletal
analysis (`CALLER_RAM_SPEED_SKELETAL.md`) found the caller's cost is mostly
*re-deriving* structure. The graph-refutation finding
(`PG_AS_GRAPH_REFUTED.md`) found *why* re-deriving it after the fact is hard
(the pg is already linearized). This plan is the synthesis: **capture the
structure before linearization destroys it, instead of trying to recover it
afterward.** If it works, it is a candidate fix for both the precision gap
*and* a chunk of the RAM/time gap simultaneously — which no single mechanism
tried this session (six for het-indel, five for RAM/speed) targeted at once.

## 8. Risks, named specifically, each with how it would be closed

| # | risk | how it is closed before any code ships |
|---|---|---|
| 1 | Ties may be dominated by repeats, not het sites | `ccnt[i]==2` gate (§3); measure the actual ratio before trusting the mechanism (§9 step 1) |
| 2 | Ties may be dominated by sequencing error at short `L` | commit-`L` threshold, swept as a formula not a constant (§4) |
| 3 | The rejected read `b` may never be placed anywhere else, so its "other haplotype" sequence is an isolated fragment with no contig context | acceptable for a first version: the caller only needs `seqs[b]` itself, not its placement, to extract the divergent bases relative to `a`'s continuation — same information `extract_bubble` needs today |
| 4 | Overlap with the existing bubble channel — could double-count or conflict | tag list is a *seed set* for search, not a replacement; existing dedup-by-locus logic (`pkey` in `caps_caller.h`, keyed on local contig sequence) already collapses duplicate evidence for the same event |
| 5 | Might only fire where the existing channel already succeeds, buying nothing | directly testable (§9 step 2) before writing any caller integration |
| 6 | Second-region sweep (leftovers, line ~1761) and `BOTHSIDE` variant path both call chain-building logic too | scope v1 to the main round-2 sweep only; extend after it is validated there, per standing rule 2 (every change gated individually) |
| 7 | Any capture inside the parallel loop (line 887-909) must not introduce a race | commit happens in the `#pragma omp single` serial block (line 934), already single-threaded — safe to append to a plain `std::vector` with no synchronization needed |
| 8 | Changing what is captured must not change what is computed | capture is pure observation of `cand`/`ccnt`/`nxt`/`prv` — reads existing values, writes to a new side vector only; zero risk of altering `pg`, `g_contig_spans`, or archive bytes when `CAPS_CALL` is unset (guarded exactly like every other `CAPS_CALL` branch already in this file) |

## 9. Exact next steps, in order — building toward a win, not gating one

**Reframe, stated explicitly because it changes how the steps below are read:**
this is not a go/no-go filter to decide whether the idea deserves to exist. The
idea is adopted — chaining computes this for free and nothing else in the
codebase gets it this cheaply. The two measurements below are **diagnostic,
not evaluative**: they tell us which of two concrete integration paths to build
first and how to tune the filter in §3-4, not whether to build at all. Either
outcome is a build instruction, not a stop condition:

1. **Instrument the round-2 commit site**: record, for every `nxt[a]`
   assignment, `ccnt[i]` and the commit `L`. Run on HG002 full chr20 (data and
   binary already on disk). Read the `ccnt==2` vs `>=3` split and the `L`
   distribution — this sets the exact threshold in §3-4's filter (a swept
   formula, not a guess) so the capture is precise from the first version
   rather than tuned after the fact.
2. **Check overlap with the existing bubble channel** on the same run. This
   picks which of two wins to build first, not whether to build:
   - **High overlap** → build the §7 integration: seed the caller's bubble
     search from the tag list instead of rebuilding `rc_reads`/`pkidx` from
     scratch. This is the RAM/time win — cutting into the 34.94 GB / 58 min
     the caller currently spends re-deriving what chaining already knew.
   - **Low overlap** → build it as an independent third calling channel
     first. This is the recall win — real het/indel signal the existing
     channel structurally cannot see, on top of the 0.592 we already lead
     DiscoSNP++ with.
   Both are real wins available from the same captured data; the measurement
   picks the order, not the outcome.
3. Implement the capture (§5, an `O(1)` append at an already-serial commit
   site, `CAPS_CALL`-gated exactly like `g_contig_spans`), then the chosen
   integration from step 2, gated byte-identical on the `CAPS_CALL`-off path
   (nothing here touches it) and F1-measured against tonight's full-chr20
   result, per standing rule 2.

Nothing in this document has been implemented yet. Step 1-2 are the fastest
path to building the right version first — not a checkpoint that could end
with "not worth it."
