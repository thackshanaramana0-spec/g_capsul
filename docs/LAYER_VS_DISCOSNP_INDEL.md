# Indel calling, layer by layer, against DiscoSNP++ source

Written after reading `~/DiscoSnp/tools/kissnp2/src/Bubble.cpp` rather than
reasoning about it. The question this answers: **our SNV precision is 0.950 and
our indel precision is 0.658, against their indel precision of 0.936. The
designs are broadly the same shape, so where exactly do they diverge?**

Nine earlier attempts failed because they all added filters *downstream*. This
compares the two pipelines layer by layer and asks which layer differs.

## The comparison

| # | layer | DiscoSNP++ (kissnp2) | ours (Method B) | verdict |
|---|---|---|---|---|
| 1 | k-mer set | GATB dBG, cascading Bloom | `kc` hash table, built by the compressor | equivalent; ours is free |
| 2 | bubble start | branching node, >=2 successors | same, both orientations | same |
| 3 | SNV closure | lockstep `expand`, `nextNode1==nextNode2` | lockstep pairwise walk | same — and we win here (0.879 vs 0.847) |
| 4 | **indel proposal** | `start_indel_prediction`: BFS extends ONE path until its last base equals the other path's start base, bounded by `max_indel_size` (D=100) | `find_sb` (Onodera superbubble) over unitig traces | **DIFFERENT** |
| 5 | **indel closure geometry** | `expand` walks both paths in LOCKSTEP after the inserted block, so the alleles are identical apart from ONE contiguous insertion — enforced by the traversal | two unitig traces that reach a common node; nothing requires them to agree in between | **DIFFERENT — this is the gap** |
| 6 | branching during extension | `checkBranching`, `authorised_branching=0` (default `b=0`): reject if either path node branches | `branch_walk` stops at the first junction, so traces are non-branching unitigs | already equivalent |
| 7 | indel positional ambiguity | `checkRepeatSize`: `size_repeat = k-2-min(ext1,ext2)`, reject if > 20 | absent | ported; measured (below) |
| 8 | low complexity | `checkLowComplexity` (DUST over trinucleotides) — **disabled by default**, `l="-l"` | absent | not a gap: off in their default runs |
| 9 | read coherence | kissreads2: per-path mapping, quality, per-position coverage | containment + quality bitmap from `kc` | ours is weaker (per-position test is inert by construction) |
| 10 | multi-polymorphism | `P=3` | `MAXPOLY=1` | theirs is more permissive; we still win SNVs |

Layers 1-3 are the same and we win the SNV comparison outright, which is the
evidence that the shared parts are sound. Everything specific to indels sits in
layers 4, 5 and 7.

## Layer 7 — ported, and it does NOT explain the gap

`checkRepeatSize` (Bubble.cpp:952) measures how far the bubble closure was
truncated by a repeat: a clean event leaves the shorter path an internal
extension of `k-1`, and an ambiguous one closes early.

The port is correct — instrumenting `min_ext` on our own bubbles reproduces
their predicted geometry, a clear mode at 29-31 against their expected
`k-1 = 30`:

    min_ext histogram: 2:2 3:2 7:2 9:2 10:1 12:2 13:2 14:4 15:9 16:1 17:7
                       18:6 19:2 20:8 21:8 22:4 23:6 24:2 25:2 26:4 27:2
                       28:6 29:10 30:4 31:5 32:2 33:2

But sweeping it is a **monotone loss** — precision barely moves while recall
collapses, so it removes true and false positives at nearly the same rate:

| `max_indel_ambiguity` | requires `min_ext >=` | TP | FP | P | R | F1 |
|---|---|---|---|---|---|---|
| off (baseline) | — | 43 | 16 | 0.729 | 0.512 | **0.601** |
| 20 (their default) | 9 | 41 | 15 | 0.732 | 0.488 | 0.586 |
| 15 | 14 | 39 | 15 | 0.722 | 0.464 | 0.565 |
| 12 | 17 | 33 | 13 | 0.717 | 0.393 | 0.508 |
| 10 | 19 | 28 | 10 | 0.737 | 0.333 | 0.459 |
| 8 | 21 | 24 | 7 | 0.774 | 0.286 | 0.417 |
| 5 | 24 | 16 | 6 | 0.727 | 0.191 | 0.302 |

SNV F1 is 0.886 at **every** setting — their code returns early for
equal-length paths and so does ours, so this is provably SNV-safe.

Kept behind `CAPS_DBG_NOAMB` (default: the gate is on, at their published 20).
It is their gate at their setting, so it is not a value fitted here.

## Layer 5 — the actual gap

Their `expand` closes only when both paths reach the same node while walking in
LOCKSTEP, so the two alleles are identical everywhere except one contiguous
inserted block. **The geometry is guaranteed by the traversal**; a non-indel
shape cannot be emitted.

Ours reads allele strings off two unitig traces that reach a common node via
`find_sb`. Nothing requires the traces to agree in between, so two paths
differing at several positions are emitted as one multi-base REF/ALT — which
vcfeval scores as an INDEL. This repository had already recorded the symptom,
in a comment at the emission site, without connecting it to a cause:

> "Two same-length paths can still differ at several positions (e.g. ACGT vs
> TGCA); those are emitted as a multi-base REF/ALT, which vcfeval scores as an
> INDEL."

The test is exact and parameter-free. `s_long` is `s_short` with one contiguous
block inserted iff the common prefix and common suffix together already cover
`s_short`:

    LCP + LCS >= |s_short|

Anything that fails it differs by more than a single indel and is a traversal
artefact. Note what this does NOT ask: it says nothing about the sequence being
repetitive, which is the question already measured not to separate (true
positives are MORE homopolymeric than false ones, 0.531 vs 0.427 at n=5,001).

Measured on the tuning window — the first change in this effort that removes
false positives and no true ones:

| | TP | FP | FN | P | R | F1 |
|---|---|---|---|---|---|---|
| baseline | 43 | 16 | 41 | 0.729 | 0.512 | 0.601 |
| + shape test | 43 | **14** | 41 | **0.754** | 0.512 | **0.610** |

SNV unchanged at 0.886. Gated by `CAPS_DBG_NOSHAPE` for A/B.

## Why the earlier nine attempts could not have worked

Every one of them was a filter over the emitted records — length bounds, STR
context, coverage gates, closure ordering. The defect is in **how the records
are generated** (layer 5), not in how they are scored. A filter cannot recover
a bubble whose two alleles were never required to be one event.
