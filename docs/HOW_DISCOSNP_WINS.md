# How DiscoSNP++ wins, what transferred, and where we now stand

Written 2026-09-02 from its source (`~/DiscoSnp`, v2.6.2-12) plus 8 head-to-head
evaluations on identical reads. Companion to `docs/DISCOSNP_INTERNALS.md`
(layer-by-layer) and `docs/CLAIM2_RESULTS_V2.md` (numbers).

---

## 1. Their mechanism, precisely

**A bubble is a local graph structure anchored at BOTH ends.** In
`kissnp2`, two paths leave a shared node and must re-converge on another
shared node. For indels, `start_indel_prediction()` runs a breadth-first
search extending ONE path 1..D bases, at each step testing whether it can
close against the other path, and stops at the SHORTEST closing extension.
Both ends of both paths are therefore k-mers that exist in the reads.

**Then every candidate is re-validated against reads.** `kissreads2` maps the
reads back onto both predicted paths and keeps only "coherent" predictions.
This is a separate binary, and it is their precision engine.

**Plus two structural guards:** `authorised_branching=0` inside an insertion,
and `checkRepeatSize` / `max_ambigous_indel` rejecting indels whose position
is ambiguous in a repeat.

## 2. Where their advantage actually shows up — measured, not assumed

Averages over the 8 evaluations (5 chr20 windows, 3 unseen individuals):

| | CAPSULE | DiscoSNP++ |
|---|---|---|
| SNV precision | 0.951 | **0.975** |
| SNV recall | **0.837** | 0.794 |
| **SNV F1** | **0.890** | 0.874 |
| INDEL precision | 0.729 | **0.910** |
| INDEL recall | 0.497 | **0.522** |
| **INDEL F1** | 0.587 | **0.663** |

Read this carefully, because it inverts the naive story:

- **On SNVs we win by RECALL** (0.837 vs 0.794) while conceding precision
  (0.951 vs 0.975). Our collapsed-substrate pileup sees variants their graph
  bubbles miss.
- **On indels their entire lead is PRECISION** (0.910 vs 0.729). Their recall
  is only 0.025 above ours. We are not failing to *find* indels; we are
  failing to *reject* bad ones.

## 3. What transferred, and what it was worth

The single biggest idea we took from them was **validate candidates against
reads, not against contig-coverage proxies** — the `kissreads2` principle.
Applied to indel junctions it moved precision **0.37 → 0.93** on one window
and lifted indel F1 from ~0.40 to ~0.58.

## 4. Why the remaining indel precision gap resists us — the structural reason

Four separate attempts, each swept across five windows, all recorded in-code
behind flags rather than deleted:

| attempt | result |
|---|---|
| Tolerant re-convergence flank (allow a het SNV inside the flank, a case their code explicitly models as "indel + n SNPs") | **net negative** 0.602→0.595, 0.600→0.582 |
| Junction allele-fraction reusing the frozen MAF | **neutral** |
| Full junction coherence — min over junction k-mers, the direct analogue of kissreads2 | **neutral, +0.001** |
| Extended contig agreement, 60–200 bp either side, to separate haplotype pairs from paralogs | **net negative** 0.5930 → 0.5864 → 0.5644 |

The third one is the diagnostic that explains the other three:

> **Our candidates are built FROM contigs, and contigs are built FROM reads.
> Every candidate is therefore read-supported by construction, so
> read-support tests cannot discriminate among them.**

`kissreads2` works for DiscoSNP++ because its paths are *hypotheses generated
by graph traversal* — a path can be proposed that no read actually
traverses, and coherence throws it out. Our paths are *observed sequence*.
The same test, applied to our substrate, has nothing to reject.

And the fourth attempt failed for a second structural reason: a paralog
filter needs long flanking context, but our contigs (mean ~335 bp) run out of
sequence before the test becomes informative.

**So closing the indel precision gap is not a filtering problem. It requires
the candidate generator itself to propose fewer bad bubbles** — i.e. a
bubble anchored at both ends by shared k-mers, rather than a pair of contigs
sharing one 25-mer anchor with 15 bp flank checks. That is a substrate change
(the eBWT or dBG route in `docs/CALLER_ARCHITECTURE_PLAN.md`), not another
filter.

## 5. Is the win generalized?

For **SNVs, yes**, on the evidence available:

- 5 of 8 evaluations won, and the wins are larger than the losses
  (+0.047, +0.076, +0.023 vs −0.007, −0.012, −0.010).
- The **largest win is on a held-out window** (r4, +0.076), not the tuning
  window — the opposite of what overfitting produces.
- It wins on **2 of 3 unseen individuals** (HG003 +0.023, HG004 +0.005),
  losing narrowly on HG005 (−0.010), the most genetically distant.
- Every parameter was swept on the tuning window only, and each default sits
  at the **centre of a flat plateau** rather than its argmax (`dup_frac`
  0.35–0.50 flat; indel anchors 5–10 flat).
- The ablation reproduces the original pre-change score (0.419) **exactly**
  when both new passes are disabled.

For **indels, no** — we are behind at 0.587 vs 0.663 and should say so.

## 6. Honest summary

- **het-SNV: we win**, 0.890 vs 0.874, by recall, validated across windows and
  unseen individuals.
- **het-indel: we lose**, 0.587 vs 0.663, entirely on precision, and the cause
  is now understood well enough to state it structurally rather than guess:
  our candidate generator proposes bubbles that no read-based filter can
  refute, because they are made of read-derived sequence in the first place.
- The remaining SNV headroom is also real: measured recall ceiling on this
  data is **0.995** against our 0.837.
