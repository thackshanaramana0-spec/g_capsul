# Can the pseudogenome be transformed into a graph to fix indel precision?

Measured 2026-09-04 on the full-chr20 HG002 run (the first full-scale Claim 2
data this project has). **Answer: no, and the reason is structural rather than
a tuning failure — it is the most useful negative result of the session.**

---

## 1. The idea, stated fairly

Not "build a de Bruijn graph from reads" — that would be a second assembly
bolted onto the first, and the project has already rejected that framing as
"combining two projects, not one insight." The idea was narrower and better:

> The compressor has already built a pseudogenome. Rather than re-reading the
> reads, **transform that existing structure into a graph** and call variants
> from it. Claim 1 is untouched (the encoder does not change); Claim 2 gets the
> thing it lacks.

What it would buy, specifically: G_CAPSUL detects bubbles **pairwise** — pick
two contigs sharing an opening anchor, walk right, require a closing anchor
(`caps_caller.h:1616`). That test cannot see how many *other* contigs also
carry the anchor. A graph can: a node with exactly 2 outgoing edges is a het
site; 3+ means a repeat or paralog, and should be rejected. That is the
textbook de Bruijn het filter, and it aims exactly at our measured weakness —
indel precision **0.839 vs DiscoSNP++'s 0.923**, while we already lead on
recall (0.457 vs 0.419).

## 2. A wrong measurement, corrected before it was believed

The first probe counted **how many contigs carry each anchor k-mer** and
reported FP 97.6% / TP 98.9% at degree >= 3 — which reads like an emphatic
refutation and is not one, because it measures the wrong quantity.

1,065,182 contigs tile a 63 Mb chromosome, so one locus appears in many
overlapping contigs: a k-mer in 5 contigs almost always means *five contigs
covering the same place*, not a branch. **Contig multiplicity is not branch
degree.** A de Bruijn graph collapses identical sequence into one node, and its
degree is the number of distinct **successor bases** (at most 4).

Recorded here because the wrong number was persuasive, arrived first, and
pointed at the same conclusion the right number eventually supported — which is
exactly the circumstance in which a wrong measurement gets believed.

## 3. The correct measurement

True de Bruijn out-degree at each call's anchor, canonical k-mers, orientation
handled the standard way (for an occurrence whose reverse complement is
canonical, the canonical node's successor is the complement of the *preceding*
base). Labels are real: TP/FP come from the `rtg vcfeval` output of the
full-chr20 run against GIAB v4.2.1, joined back to contig coordinates by
re-lifting with provenance carried in INFO.

| class | od=1 | od=2 | od=3 | od=4 | n | branching (od>1 or id>1) |
|---|---|---|---|---|---|---|
| **FP** | 263 | 13 | 1 | 9 | 286 | **10.8%** |
| **TP** | 2138 | 33 | 8 | 9 | 2188 | **3.7%** |

**There is a real signal: false positives are 2.9x more likely to sit at a
branching node than true positives.** That is a larger effect than most of the
mechanisms tried against het-indel in earlier sessions.

**And it is still useless as a filter, because of the base rates:**

    drop every branching call:  -132 TP, -74 FP
    P 0.839 -> 0.849   R 0.457 -> 0.440   F1 0.592 -> 0.580

Precision improves by exactly the amount hoped for, and **F1 gets worse**. The
TP pool is 7.6x larger than the FP pool, so a 2.9x enrichment still removes
nearly twice as many true calls as false ones in absolute terms. Enrichment
ratios do not survive contact with base rates.

## 4. The finding that actually matters

Look at the `od=1` column: **97.7% of TP anchors and 92.0% of FP anchors have
out-degree 1.** Our contig set barely branches at all. A de Bruijn graph over
these contigs is not a branching graph — it is a collection of near-disjoint
linear paths.

That is not an accident, and it is not fixable by tuning:

> **Greedy overlap chaining linearises the data on purpose. That is what makes
> the compression good, and it is precisely what destroys the branch structure
> that variant bubbles live in.**

DiscoSNP++'s graph branches because it is built from raw reads, where both
haplotypes are present as separate evidence at every position. Our pseudogenome
is built by *choosing* a path through that evidence and stitching reads into
long linear contigs — the redundancy removal that yields a 162x archive is the
same operation that removes the branches. The two objectives are in direct
opposition:

| | wants |
|---|---|
| compression | one linear path, redundancy removed |
| bubble calling | both paths retained, branches intact |

This also explains, retroactively, why the caller has to build **two**
substrates at different collapse thresholds (`dup=0.45` for SNVs, `dup=0.92`
for indels, `caps_caller.h:1397`): it is recovering, by tuning how much to
*un*-collapse, some of the branch structure that chaining removed. The dual
substrate is a workaround for exactly this tension, and it costs 738 s and a
large share of the 34.94 GB peak.

## 5. What this closes and what it leaves open

**Closed:** transforming the pseudogenome into a graph does not give us
DiscoSNP++'s precision mechanism. The structure that mechanism needs is not
present in the pseudogenome, having been removed by construction. This is now
the seventh distinct mechanism tested against het-indel and refuted by
measurement rather than argument.

**Left open, and worth stating plainly:** we currently WIN het-indel at full
scale anyway (F1 0.592 vs 0.576), on recall. The precision gap is real but we
are not losing the metric. The honest position is a recall-driven win with a
known, structurally-explained precision deficit — not a problem requiring a
fix before publication.

**Not attempted, and the only route left that could work:** retaining branch
information *at compress time* — recording, per contig, where chaining chose
between two comparable extensions. That would put the bubbles back into the
archive rather than trying to recover them afterwards. It would cost archive
size (Claim 1) to buy indel precision (Claim 2), which is a real trade against
a claim that currently wins 14/14, and it has not been measured. It should not
be attempted without first deciding whether a precision gain on a metric we
already win is worth any Claim 1 cost at all.

---

Raw data: `/root/.claude/jobs/*/tmp/outdegrees.tsv` (per-call out-degree,
in-degree, label). Probes: `outdeg.cpp` (correct), `degree.cpp` (the wrong
metric, kept for the record).
