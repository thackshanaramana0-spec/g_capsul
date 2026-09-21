---
Date: 2026-09-19
Title: The Central Insight and Mechanism — How Claims 1, 2, 3 Are One Argument
Purpose: Synthesizes the three per-claim mechanism files into the single
  cross-cutting argument the paper actually makes, anchored directly in
  the manuscript's own Introduction and Conclusion text.
When to refer to this file: Writing the Discussion section; answering "what
  is this paper's one-sentence contribution"; checking whether a new result
  belongs to the paper's thesis or is a tangent; deciding how to frame ARCH.
Keywords: insight, mechanism, allele splitting, reconciliation, substrate,
  build_substrate, collapse, ARCH, cross-claim, unifying argument
---

# The central insight and mechanism

## The insight, in the paper's own words (verified, not paraphrased from memory)

**Introduction** (line ~342): *"ADDRESSABLE asks whether that same
organization can also provide a usable positional basis for operations
that would otherwise require a new assembly, alignment, or coordinate-
bearing representation."*

**Conclusion** (line ~1303-1304): *"The broader implication is that a
compression-derived representation can remain useful as a computational
substrate rather than serving only as an intermediate for reconstruction...
G_CAPSUL therefore provides a practical basis for a broader class of
archive-native genomic analyses in which storage and computation reuse the
same underlying organization rather than reconstructing it independently
for each task."*

**This is the paper's actual thesis, and it is bigger than genomics**: a
compressor's inferred structure need not be a disposable intermediate. It
can be the basis for the analysis that follows, without rebuilding a
second representation per task.

## The specific mechanism that proves it, stated once, used three times

**The mechanism** (verified against the manuscript's own Results section,
line ~305-306): a size-minimizing, reference-free compressor is
structurally pushed to place the two alleles of a heterozygous site on
*different* pseudogenomic segments once their divergence passes some
threshold — one location plus a growing mismatch list costs more bits than
two small, separate placements. This is not an implementation quirk; it
follows from the three preconditions any such compressor meets (builds a
shared reference, encodes position+mismatches, minimizes total size — see
`refer_paper_docs/claim3/mechanism_insight_claim3.md` Layer 2 for the full,
freestanding argument).

**The same mechanism, measured three separate times, three separate ways —
this is the actual structure of the paper's evidence, not three
coincidences:**

1. **Claim 2 (calling)**: F1 0.431 → 0.888 in ablation when the split
   alleles are reconciled before calling (manuscript Results, ~line 306).
   T2.4 (multi-allelic) and T2.5 (tetraploid) are the same reconciliation
   generalized to N alleles and to non-diploid ploidy respectively — not
   two more separate wins, the same mechanism at higher order (see
   `refer_paper_docs/claim2/mechanism_insight_claim2.md`).
2. **Claim 3 (retrieval)**: a single pseudogenome coordinate reaches both
   alleles of a het site only a small fraction of the time; bilateral
   anchoring (T3.5) — reconciling both split addresses — recovers
   1,451/1,452 (99.93%) natively, 100% with the completion index. Same mechanism, applied
   to "which reads cover this locus" instead of "what variant is here"
   (see `refer_paper_docs/claim3/mechanism_insight_claim3.md`).
3. **Claim 1 (compression)**: this is where the mechanism is *caused*, not
   fixed. The same overlap-chaining process that produces allele-splitting
   is what makes T1.1's compression ratio possible in the first place, and
   the placement/deviation data it necessarily computes is what Claims 2
   and 3 reuse for free (`CAPS_SPANS`/`contig_spans` measured at zero
   archive-byte cost — see `refer_paper_docs/claim1/mechanism_insight_claim1.md`).

**Additional refutation-based evidence for the same tension, not a fourth
pillar**: a tempting shortcut (build a de Bruijn graph directly over the
pseudogenome, to run DiscoSNP++'s bubble caller on it) was measured and
refuted — 97.7% of anchors have out-degree 1, because greedy overlap
chaining *linearizes the data on purpose*, which is exactly what makes the
compression ratio good and exactly what destroys the branch structure a
bubble caller needs. See
`refer_paper_docs/claim2/architecture_vs_discosnp_claim2.md`. This is not
counted as a fourth instance of the mechanism above — it is the same
tension observed from the failure side, confirming by direct measurement
rather than assumption that Claim 1's compression and graph-branch
structure trade off against each other.

## Why this makes the three claims one argument, not three separate results

A reader could mistake COMPACT/FAITHFUL/ADDRESSABLE for three unrelated
"we are also good at X" claims. The verified mechanism above is what
prevents that reading: **Claim 1 causes the phenomenon that Claim 2 and
Claim 3 each independently discover and independently fix, using
structurally the same reconciliation idea, in two different domains
(calling and retrieval).** That two separately-built systems (a variant
caller and a query engine) converge on the same fix for the same
underlying cause is itself evidence the mechanism is real and general, not
tuned to make one benchmark look good.

## Why Claim 3's operations are not five bolted-on features — the architectural version of the same point

Verified in code this session (`refer_paper_docs/claim3/code_mapping_claim3.md`):
T3.3 and T3.4 are literally the same function in `capsule_decode.cpp`; T3.5
is not a separate code path at all, only an orchestration layer calling
T3.4's primitive twice. This is the structural signature of one retained
data structure being reused in increasingly demanding ways, not five
independently-engineered features — the code-level evidence for the same
claim the Conclusion states in prose.

## The honest asymmetry between Claim 2's and Claim 3's versions of the mechanism

Claim 2's reconciliation (`build_substrate`/`collapse_contigs`) happens
*before* scoring, baked into the caller's substrate construction. Claim
3's reconciliation (bilateral anchoring) happens *after* the archive is
built, at query time, via two calls to an unmodified primitive. **This
difference is itself worth stating in the paper**: it demonstrates the same
underlying fix can be applied either at analysis-construction time or
after the fact, on a finished archive — reinforcing the "substrate,
reusable in more than one way" framing rather than implying the fix had to
be anticipated and hard-coded once.

## The one-sentence version of this entire file, from the project's own prior verified work

Found this session in `refer_paper_docs/CLAIM3_LOCKED.md` (verified 2026-09-09, still
true independent of that document's otherwise-superseded numbers):

*"Everyone else's archive stores what the reads SAY. Ours stores where
they SIT. The compressor had to assemble in order to compress, so the
coordinate system is a byproduct rather than an addition — which is why
locus retrieval, per-base coverage and reference-free variant calling all
come out of the same structure."*

This is the whole paper's argument in two sentences, and it is worth
having ready as a direct quote for a talk, an abstract, or a reviewer
response.

## Where ARCH fits, precisely, if adopted

ARCH is a proposed name for the general pattern this file describes
(construct once, retain, harness for compression and for downstream
computation) — not for the allele-splitting mechanism specifically, which
is the *evidence* for the pattern, not the pattern itself. Conflating the
two in the paper's wording would be an error: name ARCH for the general
substrate-reuse claim (paired with G_CAPSUL's introduction, per
`paper_overview.md`), and describe allele-splitting/reconciliation as the
mechanism that demonstrates it, not as ARCH itself.
