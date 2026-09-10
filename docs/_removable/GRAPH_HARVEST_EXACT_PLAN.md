# Graph harvest — exact, line-by-line implementation plan

Written 2026-09-04. Read-only analysis; **nothing in here is implemented yet**.
Every line number is from the current working tree (`stages/106_inprocess.cpp`,
`include/caps_caller.h` as of this session's uncommitted state).

Supersedes the vaguer `CHAINING_BRANCH_TAG_PLAN.md` — that document framed this
as a "side-list of tags." That framing was too small. The correct statement is
below.

---

## 0. The thesis, stated exactly

**The overlap graph already exists inside the sweep. It is computed, used once,
and discarded — every level, for every read.**

`stages/106_inprocess.cpp:906`:

    cand[i*CCAP+c]=b;

For open tail `a` at overlap level `L`, `cand[i*CCAP + 0 .. i*CCAP+ccnt[i]-1]`
is **exactly the out-edge set of node `a` in an overlap graph**: every read `b`
whose first `L` bases are verified equal (`rcmp`, line 903) to `a`'s last `L`
bases. Chaining then keeps the first untaken entry (line 941-942) and drops the
rest, permanently (`open_tails` filter, line 841-843, removes `a` once `nxt[a]`
is set).

So we are **not building a graph.** We are stopping the deletion of one.

This is the precise complement of `docs/PG_AS_GRAPH_REFUTED.md`: that document
measured 97.7% out-degree-1 on *finished contigs* and correctly concluded a
graph cannot be recovered after the fact. The out-degree is 1 there **because
this commit already collapsed it**. Harvesting at line 941 reads the degree
before the collapse.

---

## PHASE 1 — Harvest (encoder side)

### 1.1 Declaration

**File:** `stages/106_inprocess.cpp`
**Insert after line 787** (`std::vector<uint32_t> cand; std::vector<uint8_t> ccnt;`)

    // Harvested overlap-graph branch points (CAPS_CALL only; see Phase 1.3).
    struct GEdge { uint32_t a, b_taken, b_alt, L; };
    std::vector<GEdge> g_branch;

`g_branch` must be declared **outside** `sweep()` (which is a lambda starting at
line 762) so it survives the call, and **before** line 819's level loop.

### 1.2 Capture site — the exact insertion

**File:** `stages/106_inprocess.cpp`
**Location:** inside the `#pragma omp single` commit block, lines 934-949.

Current code (lines 936-947):

    for(size_t i=0;i<w;++i){
        const uint32_t a=tails[i];
        bool done=false;
        for(uint8_t c=0;c<ccnt[i];++c){
            const uint32_t b=cand[i*CCAP+c];
            if(prv[b]!=NONE||ch_h[a]==b) continue;       // taken, or would cycle
            nxt[a]=b; prv[b]=a; ovl[a]=L;
            uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h; ++links;
            done=true; break;
        }

**The edit:** capture the *rejected but still-viable* alternative at the moment
the winner is committed. Immediately before `done=true; break;`:

    if(CAPS_CALL && ccnt[i]==2){
        const uint32_t other = cand[i*CCAP + (c==0?1:0)];
        if(other!=b && prv[other]==NONE && ch_h[a]!=other)
            g_branch.push_back({a,b,other,L});
    }

**Why every clause is there:**

| clause | reason |
|---|---|
| `CAPS_CALL &&` | Claim 1 runs must be bit-identical. Same gate pattern as `g_contig_spans` (lines 1269, 1277, 1765). |
| `ccnt[i]==2` | Exactly-2 is the het-bubble shape. `>=3` is a repeat/paralog family — the FP class `PG_AS_GRAPH_REFUTED.md §3` measured at 2.9x enrichment. **This gate is not optional.** |
| `prv[other]==NONE` | The alternative must be genuinely unclaimed. If it is already in another chain, this is not a branch — it is two chains merging. |
| `ch_h[a]!=other` | Same cycle guard the commit itself uses (line 941). |
| placed inside `omp single` | Line 934 is already serial. **No synchronization needed, no race possible.** Placing it in the parallel block at line 906 would require a critical section. |

### 1.3 Cost, exactly

`sizeof(GEdge)` = 16 B. Growth is bounded by one push per committed link, and
only for `ccnt==2` commits, so `g_branch.size() <= links`.

- If ties track true heterozygosity (~55K loci on chr20): **~900 KB**
- Absolute worst case (every one of 12.6M commits is a 2-way tie): 200 MB

**Archive cost: ZERO bytes.** `CAPS_CALL` is off for all 14 Claim 1 datasets,
and even when on, this never reaches the container writer — Claim 2 consumes it
in the same process. This is the correction to the earlier "150-350 KB archive
cost" estimate: that answered a different question (a hypothetical
decode-from-archive Claim 2), which is not the architecture.

### 1.4 Round 1 vs round 2

`sweep()` is called at line 997 (round 1) and line 1228 (round 2). **Round 1's
result is discarded** — line 764-765 refills `nxt`/`prv`/`ovl` on the next
call, and round 1 exists only to compute `admit` (line 990-994).

**Therefore:** `g_branch.clear()` must be the first statement of `sweep()`
(insert at line 763, beside `links=0; probes=0;`). Only round 2's harvest
survives, which is correct: round 2's `nxt`/`ovl` are what build `pg` and
`g_contig_spans` (lines 1273-1277), the exact structure the caller receives.

### 1.5 Coordinate translation — free, reusing existing output

**File:** `stages/106_inprocess.cpp`, **insert after line 1296** (`lap("emit chains")`),
i.e. after `pg`, `ppos` and `g_contig_spans` are all fully built.

    if(CAPS_CALL){
        // ppos[a] is the pg offset of read a, already computed by the
        // chain-emission walk (lines 1274-1276). Binary-search g_contig_spans
        // exactly as find_cid_mm() already does at line 2822.
        for(auto& e : g_branch){
            if(ppos[e.a]==UINT64_MAX) continue;   // a wasn't emitted into pg
            /* resolve (contig_id, offset_in_contig) from ppos[e.a] */
        }
    }

No new placement work: `ppos[]` is populated at lines 1274-1276 during the walk
that builds `pg`, and the binary search over `g_contig_spans` is the same
2-line idiom already used at line 2822.

---

## PHASE 2 — Consume (caller side)

### 2.1 What the caller receives

`include/caps_caller.h`'s `CallData` (used at `run_variant_call`, line 733)
gains one field, populated exactly like `cd.contigs` is at
`106_inprocess.cpp:1782`:

    std::vector<GraphBranch> branches;   // (contig_id, contig_pos, read_a, read_taken, read_alt, L)

### 2.2 The bubble walk — reusing what already exists

**This is the part that does NOT need new algorithm work.** For each branch:

- **Path A** = `seqs[e.a]` continued by the chain: `nxt[e.b_taken]`, etc. But
  the simpler and sufficient form is that path A is **already in the contig** —
  `cd.contigs[cid]` from position `contig_pos` onward IS path A, because that
  is precisely the chain the commit built.
- **Path B** = `seqs[e.b_alt]`, the rejected read, whose first `L` bases align
  to `a`'s suffix by construction (`rcmp` guaranteed it at line 903).

Then call the **existing, tuned, already-shipping** extractor:

    Bubble bub = extract_bubble(contigA, posA, seqB, posB, MAXINDEL, FLANK);

`extract_bubble` is `include/caps_caller.h:153`. It already handles:
- generic SNV/indel divergence + re-convergence (`flank_match`, line 126)
- the **homopolymer run-length branch** (lines 176-205) written specifically
  because every missed 1 bp indel on HG002 r2 was inside a 7-8 base homopolymer

**Nothing about the bubble algorithm is new work.** The novelty is exclusively
*where the two paths come from* — harvested free from assembly, rather than
reconstructed at a cost of 738 s + tens of GB.

### 2.3 What this can replace, with measured costs

From this session's own full-chr20 phase profile:

| structure / phase | current cost | why it exists | replaceable? |
|---|---|---|---|
| `build_substrate` #2 (`caps_caller.h:1397`, `dup=0.92`) | **738 s** | keep haplotypes apart so bubbles are findable | **yes** — branches give haplotype pairs directly |
| `rc_reads` (`caps_caller.h:2147`) | **15-30 GB** | find reads sharing a right-context anchor | **yes** — that is what a branch IS |
| `pkidx`/`pcount`/`porient` (line 2093) | ~8-12 GB | find contig-unique anchors | **yes** — same |
| `kidx` bubble scan (line ~1443) | in `indel_pass`'s 2091 s | pair contigs by shared k-mer | **partially** — keep as fallback for loci with no harvested branch |
| pileup path (`parallel_loop`, 480 s) | 480 s | SNV calling by depth | **keep** — different evidence class, and it is where SNV 0.849 comes from |

### 2.4 Fallback policy — non-negotiable

The harvested set is a **seed set, not a replacement**, in v1. Loci with no
harvested branch still go through the existing `kidx` path. Recall cannot
regress by construction; only cost falls. Removing the fallback is a separate,
separately-gated decision after v1 is measured.

---

## PHASE 3 — Order of work, with gates

| # | step | gate |
|---|---|---|
| 3.1 | Phase 1.1-1.4 harvest + a `fprintf` of `g_branch.size()`, `ccnt` histogram, `L` distribution. **No consumption.** | Claim 1 path byte-identical (E. coli, `CAPS_CALL` unset) |
| 3.2 | Run on HG002 full chr20 (binary + data on disk) | reports: how many 2-way ties, at what `L` |
| 3.3 | Phase 1.5 coordinate translation + dump branches to TSV | spot-check N branches against known GIAB het sites |
| 3.4 | Phase 2.1-2.2, emitting into the VCF as a new channel **alongside** existing ones | F1 must not regress vs tonight's 0.849/0.592 on full chr20 |
| 3.5 | Phase 2.3: disable the structures the branches replace, one at a time | each individually gated on F1 + measured RAM/time |

Standing rule 2 applies throughout: output-preserving changes gated
`cmp`-identical; output-changing changes must not regress and must improve the
aggregate. Standing rule 3: a failed gate is reverted **and recorded**, not
tuned until it passes.

---

## PHASE 4 — Realistic expectations

Same table maintained all session. **Left column is measured tonight. Right
column is a projection and is explicitly labelled as one.**

### 4.1 Cost build-up — where every second and every GB goes

| stage | cost | why |
|---|---|---|
| trunk (chaining + pg build) | **3.6 min, 2.35 GB** | this is Claim 1. We pay it regardless. |
| graph harvest | **~0** | already computed inside the sweep (line 906); we stop deleting it |
| bubble walk (`extract_bubble`) | **seconds, ~1 GB** | a traversal over harvested pairs, not a rebuild |
| **total** | **~4 min, ~3.5 GB** | |

Both figures are *additive over a measured floor*, which is why they are tight
rather than a range: the trunk numbers are measured from tonight's run, and the
only additions are a 16 B/edge vector and a walk over it.

### 4.2 The three axes

| axis | now (measured, full chr20 HG002) | DiscoSNP++ (measured, same data) | after graph harvest (PROJECTED) |
|---|---|---|---|
| SNV F1 | **0.849** | 0.847 | 0.845 - 0.855 |
| indel F1 | **0.592** | 0.576 | 0.60 - 0.65 |
| peak RAM | **34.94 GB** | **3.45 GB** | **~3.5 GB** |
| wall time | **62.5 min** | **1.3 min** | **~4 min** |
| archive | **26.4 MB** | none | 26.4 MB (unchanged) |

**Basis for each projection, so it can be checked rather than trusted:**

- **RAM ~3.5 GB:** trunk floor 2.35 GB (measured tonight) + `g_branch`
  (~900 KB) + the bubble walk's working set (~1 GB). `rc_reads` (15-30 GB),
  the `pkidx` trio (8-12 GB) and `build_substrate` #2 are *deleted*, not
  shrunk — §2.3. This lands at DiscoSNP++'s 3.45 GB, not near it.
- **Time ~4 min:** trunk 3.6 min (measured) + a walk over ~55K harvested
  branches (seconds). `build_substrate` #2 (738 s) and the bulk of
  `indel_pass` (2091 s) are deleted. `parallel_loop` (480 s) is retained only
  if the pileup path stays — see the note below.
- **indel 0.60-0.65:** current 0.592 is recall-limited (R 0.457, P 0.839).
  Harvested branches are read-level evidence at loci the contig-pair scan may
  never pair, so recall should rise. Bounded conservatively because six
  het-indel mechanisms this session did *not* beat the shipped channel.
- **SNV ~flat:** SNV comes from the pileup path, which this does not touch.

**One dependency that the ~4 min figure rests on, stated so it is not
buried:** `parallel_loop` (480 s = 8 min) is the pileup, and it is where SNV
0.849 comes from. Keeping it makes the total ~12 min, not ~4. The ~4 min figure
assumes the graph walk also supplies SNVs (bubbles yield SNVs and indels
alike — it is how DiscoSNP++ gets both from one traversal) and the pileup is
retired. **That is exactly what step 3.5 gates.** If pileup SNV accuracy proves
irreplaceable, the honest number is ~12 min at ~3.5 GB, which is still 5x
faster and 10x lighter than tonight.

**The honest limit, stated plainly:** at ~3.5 GB we reach **parity** with
DiscoSNP++ on RAM. On speed we land at ~4 min against their 1.3 min — still
~3x behind, and that entire residual is the trunk (3.6 min), which exists
because we also produce a 26 MB lossless archive and an addressable structure
(Claim 3). DiscoSNP++ produces neither: its graph is built and discarded.

So the defensible claim after this work is: **equal or better accuracy, equal
RAM, ~3x slower — while additionally emitting a lossless archive the
competitor cannot produce at all.** That is a strong position and it is what
the numbers support. "We beat them on all three axes" is not, and must not be
written: the trunk cost is real and does not go away.
