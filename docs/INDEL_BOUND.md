# Why the indel table does not flip — three measurements and one bound

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
