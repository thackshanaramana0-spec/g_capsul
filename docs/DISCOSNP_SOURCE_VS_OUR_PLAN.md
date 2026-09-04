# Reading DiscoSNP++'s real source against our graph-harvest plan

2026-09-04. Read-only. Source: `/root/DiscoSnp/tools/kissnp2/src/Bubble.cpp`
(1074 lines) and GATB's `Graph.cpp`/`Graph.hpp`.

**Why this document exists:** `GRAPH_HARVEST_EXACT_PLAN.md` was written from an
*assumed* model of DiscoSNP++'s algorithm ("walk both branches until they
reconverge"). Reading the actual code shows that model was wrong in one
load-bearing way. The plan's arithmetic survives; **one of its core claims does
not**, and the correction changes what we should build.

---

## 1. What their algorithm actually is, line by line

### 1.1 Seeding — branching nodes only

`Bubble.cpp:300-347`, `BubbleFinder::start(Bubble&, const BranchingNode&)`:

    GraphVector<Node> successors = graph.successors((Node&)node);
    if(successors.size()<2) return;                 // line 309
    for (i...) for (j=i+1...) {
        bubble.begin[0] = successors[i];
        bubble.begin[1] = successors[j];
        start_snp_prediction();                     // line 337
        start_indel_prediction();                   // line 345
    }

They enumerate **branching nodes** and take **every unordered pair of
successors** as a candidate bubble opening. Lines 317-322 deduplicate: if a
successor has multiple predecessors, only the smallest-kmer predecessor is
allowed to open the bubble.

**This is exactly the structure we harvest** — `cand[i*CCAP+...]` at
`106_inprocess.cpp:906` is a successor set, and `ccnt[i]==2` is
`successors.size()==2`. Our harvest maps onto their seeding step correctly.

### 1.2 Extension — LOCKSTEP, and this is what the plan got wrong

`Bubble.cpp:509-527`, inside `expand()`:

    successors = graph.successors (node1, node2);   // TWO-node successors
    if (successors.size() != 1) break;
    if (successors[0].first == successors[0].second) break;   // ended bubble
    local_extended_string1 += ascii(graph.getNT(successors[0].first,  sizeKmer-1));
    local_extended_string2 += ascii(graph.getNT(successors[0].second, sizeKmer-1));

`graph.successors(node1, node2)` is `getNodesCouple` (`Graph.hpp:705`), and its
implementation (`Graph.cpp:1602-1640`) loops `nt = 0..3` and emits a pair
**only when BOTH `forward1` and `forward2` extended by the SAME nucleotide `nt`
exist in the graph**:

    Type forward1 = ((graine1 << 2) + nt) & mask;
    Type forward2 = ((graine2 << 2) + nt) & mask;
    if (data.contains(forward1) && data.contains(forward2)) { ...emit pair... }

**The two paths advance in lockstep, one base at a time, always by the same
nucleotide.** The bubble closes when the two nodes become equal
(`successors[0].first == successors[0].second`, line 517).

### 1.3 Indels are a breadth-first search, not a walk

`Bubble.cpp:200-270`, `start_indel_prediction()`: it pushes onto
`breadth_first_queue`, extends ONE path while holding the other fixed, and at
each depth checks whether the extended path's last character equals the other
path's `end_insertion` character — then calls `expand()` to try to close.
Bounded by `max_indel_size` and `max_recursion_depth` (queue size).

So indel detection is: *hold one path, BFS the other up to `max_indel_size`,
attempt closure at each depth.* Not a symmetric walk.

---

## 2. The assumption that was wrong, stated plainly

`GRAPH_HARVEST_EXACT_PLAN.md` §2.2 claimed:

> "Path A is already in the contig ... Path B = `seqs[e.b_alt]`, the rejected
> read ... Then call the existing `extract_bubble`. **Nothing about the bubble
> algorithm is new work.**"

**That is wrong, and the reason matters.**

Their extension is a lockstep walk **through a graph**, where at every step
they ask "does the SAME next base exist on both paths?" — a graph membership
query (`data.contains`). A bubble closes when both paths arrive at the
*identical k-mer node*.

Our harvested path B is **a single read** (`seqs[e.b_alt]`, ~148 bp), not a
graph path. We cannot ask "does this k-mer exist on path B" beyond that read's
end, because there is no graph to query — chaining consumed the alternative
into some other chain or dropped it. `extract_bubble` compensates by comparing
two *strings* and requiring a byte-identical `FLANK` re-convergence
(`caps_caller.h:126`).

**Consequence:** we do not have their closure test. We have a string-suffix
approximation of it, which is precisely the thing this project already measured
as failing — `caps_caller.h:132-146` records the tolerant-flank experiment
(`flank_match_tol`) as **net negative** (r2 indel F1 0.602→0.595), with the
verdict *"the bubble GEOMETRY is not what limits us."*

So harvesting branch pairs and feeding them to `extract_bubble` does **not**
reproduce DiscoSNP++'s mechanism. It reproduces the *seeding* and keeps our
existing, weaker closure.

## 3. What this does and does not invalidate

**Still true, and independently valuable:**

- The harvest itself (Phase 1) is correct, nearly free, and archive-free. The
  successor set at line 906 genuinely is graph adjacency.
- The seeding correspondence is exact: their `successors.size()>=2` is our
  `ccnt[i]>=2`, and their smallest-predecessor dedup (lines 317-322) has a
  direct analogue.
- Everything the harvest could *replace* (§2.3 of the plan: `rc_reads`,
  `pkidx` trio, `build_substrate` #2) is still replaceable, because those
  structures exist to *find candidate pairs* — exactly what the harvest
  supplies. **The RAM and time argument survives intact.**

**No longer true:**

- "Nothing about the bubble algorithm is new work." Reproducing their closure
  requires a k-mer membership structure to walk against, which we would have to
  build — and building one is the "second project" this whole line of work
  exists to avoid.
- The projected indel gain (0.60-0.65) was predicated on getting their
  mechanism. With our existing closure test, the honest expectation is **indel
  F1 approximately unchanged**; the harvest changes *which loci we test*, not
  *how well we resolve one*.

## 4. The corrected position — and it is still worth building

Reading their source separates the plan into two independent claims, one of
which is much stronger than the other:

| claim | status after reading their code |
|---|---|
| **A. Harvest replaces the caller's candidate-finding machinery** | **INTACT.** `rc_reads`/`pkidx`/`build_substrate`#2 exist to find pairs; the harvest hands them over free. This is the RAM/time win and it does not depend on closure at all. |
| **B. Harvest gives us DiscoSNP++'s accuracy mechanism** | **REFUTED.** Their accuracy comes from lockstep graph extension + node-identity closure. We have neither, and acquiring them means building a k-mer graph. |

**So the build order in `GRAPH_HARVEST_EXACT_PLAN.md` §3 should be re-read as
pursuing claim A only**, with accuracy held flat as a *gate* rather than
projected as a *gain*.

### Revised expectations

| axis | now (measured) | DiscoSNP++ (measured) | after harvest (claim A only) |
|---|---|---|---|
| SNV F1 | **0.849** | 0.847 | 0.845 - 0.855 (held flat, gated) |
| indel F1 | **0.592** | 0.576 | 0.58 - 0.60 (held flat, gated) |
| peak RAM | **34.94 GB** | **3.45 GB** | ~3.5 GB |
| wall time | **62.5 min** | **1.3 min** | ~4 min (~12 min if pileup retained) |
| archive | **26.4 MB** | none | 26.4 MB |

The RAM and speed numbers are unchanged from the previous document — they never
depended on the closure mechanism. **The accuracy column is now a gate, not a
projection**, which is the honest form: we already lead indel 0.592 vs 0.576 and
tie SNV, so holding flat while collapsing 34.94 GB → ~3.5 GB and 62.5 min →
~4 min is a large, real win on its own.

## 5. What a genuine accuracy gain would require, priced honestly

To get their closure test we need k-mer membership over the reads — `contains(kmer)`.
We **already build exactly this** and throw it away: `kc` in
`caps_caller.h` (`kc_H_build`, 166.7 s) is a sorted `(kmer, count)` array over
all canonical 31-mers of all reads, with `kc_find` as its binary search. That
is a membership oracle.

So a lockstep walk in the style of `Bubble.cpp:509-527` is implementable
against `kc` without building any new index: extend both paths by the same
base, test membership of both extensions in `kc`, stop when the two extended
strings' terminal k-mers are equal. This is a real, bounded piece of work —
**not** a second assembler, and not the "combining two projects" the project
rejected.

**But it must be measured, not assumed**, and it is a *separate* change from
claim A. Sequencing: land claim A first (RAM/time, accuracy gated flat), then
evaluate the lockstep closure against `kc` as its own gated change with its own
before/after F1. Bundling them would make it impossible to tell which one moved
the numbers — the exact failure mode standing rule 2 exists to prevent.
