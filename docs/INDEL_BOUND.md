# Why the indel table does not flip — three measurements and one bound

> ## OVERTAKEN 2026-09-09 — het-indel is now a WIN, and this file is the reason why
>
> This document concluded that het-indel was a closed loss and should be
> withdrawn from the claim. **That conclusion is no longer true**, and the
> executed sweep is the authority:
>
> | | this file (earlier caller) | **final sweep, T2.3** |
> |---|---|---|
> | HG002 TP | 3,290 | **3,871** (+581) |
> | HG002 FP | 1,711 | **824** (-887) |
> | HG002 F1 | 0.516 | **0.620** |
> | vs DiscoSNP++ | 0.576 — a loss | 0.576 — **a win** |
> | across 4 individuals | — | **3 wins, 1 loss** (0.620/0.637/0.632/0.593 vs 0.576/0.587/0.595/0.605) |
>
> `benchmark/results/claim2/claim2_T2.3_indel.csv`, traced in
> `benchmark/documentation/RESULT_CODE.md`.
>
> **The analysis below was not wrong — it was RIGHT, and it named the way out.**
> Its bound `F1 = 2TP/(TP+FP+truth)` still holds exactly: at today's numbers it
> gives 0.620, which is what the CSV records. What it proved was that no FILTER
> could flip the result at TP = 3,290, and it said so explicitly: *"closing the
> gap requires a different candidate generator, not another filter on top of the
> current one."* A different generator is what changed. Both halves of the
> arithmetic moved — more true indels found AND fewer false ones — which is
> exactly what a filter cannot do and a generator can.
>
> What still stands, unchanged: the representational limit (a length change
> inside a homopolymer creates no second k-mer path, so no bubble exists for us
> or for any de Bruijn caller), the finding that true positives are MORE
> homopolymeric than false ones (so context filters remove TPs first), and the
> list of nine attempts that failed. Do not retry those.

Status: **closed**. Indels are a measured loss (F1 0.516 vs DiscoSNP++ 0.576 at
full chr20) and this file records why no further attempt is warranted, so the
next person does not spend another nine tries on it.

## The bound

Truth carries 7,768 indels in the confident region, so `FN = 7768 - TP` and

    F1 = 2*TP / (TP + FP + 7768)

At our TP = 3,290, beating 0.576 requires **FP < 366** — removing **79%** of
1,711 false positives with zero true-positive loss. Perfect precision
(FP = 0, TP unchanged) yields only **0.595**. The entire filtering family is
therefore bounded to a margin no filter has ever produced here.

Raising recall instead needs TP > 3,833 (+543, recall 0.424 -> 0.494) with no
new FPs. The next section is why that is unavailable.

## 1. Recall is capped by representation, not by thresholds

Every missed indel was traced into `kc`. At the missed sites the branching node
has coverage 12-30 but exactly ONE successor, and that holds with **no coverage
floor at all** — raw successor counts return `[0,0,21,0]`. The alternate
haplotype's k-mers are absent, not filtered.

66% of missed indels are homopolymer-length changes; most of the rest are
tandem-repeat expansions. Inserting a base into a run of that base does not
create a second path. Probing reads directly: for a `G->GT` call inside a T-run
the REF and ALT probes are the SAME STRING
(`TCTGGGTTTTTGTTTTTCGGGTTTTTTTTTTT`), each present in 6 reads. No bubble exists
— for us or for any de Bruijn caller, DiscoSNP++ included.

## 2. Reference context cannot fix precision

Measured over all 5,001 scored indels at full chr20 (not a window):

| class | n | mean run length | frac run>=4 | frac run>=6 | mean tandem score |
|---|---|---|---|---|---|
| TP | 3,290 | 5.72 | 0.531 | 0.419 | 0.447 |
| FP | 1,711 | 5.23 | 0.427 | 0.344 | 0.450 |
| FN | 4,478 | 7.66 | 0.555 | 0.458 | 0.583 |

True positives are **more** homopolymeric than false positives. Any
homopolymer/STR filter removes TPs preferentially — which is exactly what the
window-scale test measured (dropped 9 TPs to remove 7 FPs, F1 0.601 -> 0.535,
`CAPS_DBG_STR`). This confirms it at n = 5,001 and closes the whole family.

## 3. The normalisation artefact is real but small

False positives that sit near a false negative — i.e. the event was found but
scored at the wrong place or shape:

| window | FPs near an FN | share |
|---|---|---|
| exact position | 103 | 6.0% |
| +-2 bp | 210 | 12.3% |
| +-20 bp, same net length change | 230 | 13.4% |

Repairing the exact-position cases perfectly (103 FP -> TP, 103 FN -> TP) gives
F1 0.531. The +-2 bp band gives 0.548. Neither reaches 0.576, and both are
upper bounds that assume a flawless repair.

## What would actually be required

Length-aware evidence — pileup depth across the run — which resolves run length
by counting reads at each length rather than by finding a second path. That is
a different mechanism from bubble calling, not a parameter on this one, and it
is not attempted here.

## Consequence for the paper

Scope the variant-calling claim to SNVs, where the win is validated on four
independent axes (individuals incl. one held out, coverage 10-30x, ploidy 2 and
4, chromosomes 20 and 1). DiscoSNP++'s own paper leads with SNP precision and
recall. Report the indel result honestly with the mechanism above; it is a
property of de Bruijn variant calling, not of this implementation.

## Attempts, all measured, none retained

Four precision filters (two catastrophic at F1 0.035), a junction probe, four
recall gates, shortest-closure ordering, homopolymer/STR filtering. Every one
operated downstream of a signal that was never in the graph.
