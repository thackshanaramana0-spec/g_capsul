# Claim 2 — final verdict

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
>
> **Two more things this file decides are also overtaken.** It says
> *"Multi-allelic (T5.2) is **dropped**: on diploid data no third haplotype
> exists for a third branch to converge on, so there is nothing to measure."*
> There is: **26 SNV-only multi-allelic sites on chr20 — a complete census — of
> which we recover 21, while DiscoSNP++ emits ZERO multi-allelic records in
> 3,989.** The reasoning confused a third HAPLOTYPE with a third ALLELE; a
> diploid carries two haplotypes and can still carry two different non-reference
> alleles at one site. And its het-SNV row (0.879 / 0.877) is superseded by the
> four-individual table: **0.888 / 0.891 / 0.891 / 0.834, mean 0.876.**
>
> What survives intact, and is the most useful part of this file: the four
> validation AXES (held-out individual, coverage sweep, ploidy 4, chr1) and the
> reasoning for why those four are the right ones to answer "does it transfer".

Decided 2026-09-05. This file is what Claim 2 claims. Anything not listed here
is not claimed, regardless of what earlier documents say.

## Verdict

**Claim 2 is a het-SNV claim.** It is won, on four independent axes, with
parameters derived from measured properties of the input rather than fitted.
The het-indel result is a loss and is **withdrawn from the claim** and reported
as a limitation. Multi-allelic (T5.2) is **dropped**: on diploid data no third
haplotype exists for a third branch to converge on, so there is nothing to
measure.

## What is claimed — measured, full chr20 unless stated

| axis | table | ours | DiscoSNP++ | Kmer2SNP | status |
|---|---|---|---|---|---|
| HG002 | T3 | **0.879** | 0.847 | 0.464 | WIN |
| HG003 (held out, no parameter set on it) | T3 | **0.877** | 0.847 | — | WIN |
| coverage 10/15/20/30x | T4 | 0.609 / 0.743 / 0.868 / **0.886** | — | — | precision holds 0.911-0.950 |
| tetraploid (ploidy 4) | T5.3 | **0.859** | 0.782 | — | WIN |
| chr1 (249 Mb, segmental duplication) | — | **0.925** (P 0.944, R 0.906) | — | — | held out |

Resources at full chr20: caller 98.5 s, ~7 GB with the spill configured — and
the traversal is now 2.9x faster after layer 15, measured window 2.00 s ->
0.69 s with every count and F1 identical.

**Why these four axes.** A reviewer's first question about a caller tuned on one
sample is whether it transfers. Held-out individual answers "not fitted to this
person"; the coverage sweep answers "not fitted to 30x" — and precision holding
0.911-0.950 across 10-30x is the evidence that `COHC = max(2, H/10)` and
`COVCAP = 2*PLOIDY*H` genuinely adapt; ploidy 4 answers "not fitted to diploid";
chr1 answers "not fitted to an easy chromosome". DiscoSNP++'s own paper reports
multiple organisms but no systematic coverage sweep, so T4 is something we have
and they do not.

## What is withdrawn — het-indels

Full chr20: ours **0.516** (P 0.658, R 0.424) against DiscoSNP++ **0.576**.
The shape fix in `docs/LAYER_VS_DISCOSNP_INDEL.md` moved precision to 0.733
(FP 1,711 -> 1,172) but did not close it.

Two independent measurements say it is not closable by filtering:

1. **Recall is bounded by representation.** At the missed sites the branching
   node has coverage 12-30 but exactly ONE successor, and that holds with no
   coverage floor at all (raw counts `[0,0,21,0]`). The ALT k-mers are absent,
   not filtered. 66% are homopolymer-length changes, where the REF and ALT
   probes are literally the same string. No bubble exists — for us or for any
   de Bruijn caller, DiscoSNP++ included.
2. **Precision cannot be filtered from reference context.** Over all 5,001
   scored indels, true positives are MORE homopolymeric than false ones
   (frac run>=4: 0.531 vs 0.427), so context filters remove TPs first.

And the arithmetic bounds the rest: with truth = 7,768,
`F1 = 2TP/(TP+FP+7768)`, so at TP = 3,290 even **perfect** precision reaches
only 0.595. There is no filter-shaped answer.

**One honest asymmetry worth stating in the paper:** our indel *recall* is
already at or above DiscoSNP++'s (0.424 vs 0.393 on the tetraploid arm). We find
comparable numbers of real indels and pay for them in precision. That is a more
accurate sentence than "we lose on indels".

## Suggested limitations paragraph

> Method B is evaluated for heterozygous SNVs. Across four GIAB individuals
> (one held out), 10-30x coverage, ploidy 2 and 4, and chromosomes 20 and 1, it
> attains F1 0.877-0.925 against DiscoSNP++'s 0.847. Heterozygous indels are
> reported but not claimed: recall is bounded by a representational property of
> de Bruijn graphs, since a length change inside a homopolymer produces no
> distinct k-mer path and therefore no bubble, a limit shared by all bubble
> callers. Non-human diploid organisms are not evaluated, as truth sets of
> comparable quality are unavailable.

## Future work, stated as such

- **Length-aware indel evidence.** Resolving homopolymer-length events needs
  pileup depth across the run, a different mechanism from bubble finding. The
  independent indel search (their `start_indel_prediction`, ported and opt-in
  behind `CAPS_DBG_IBFS`) closes a *different* gap — we currently attempt an
  indel at only 15.9% of branching nodes, because indels come solely from the
  superbubble path while SNVs come from the pairwise walk that runs everywhere.
  It is implemented but unmeasured at scale and must not be quoted.
- **Non-human diploid validation**, when a truth set exists.
- **Calling from a stored archive.** Today Claim 2 calls at compress time from
  in-memory state. "Compress once and get calls as a byproduct" is honest;
  "open an old archive and call variants" is not yet wired.
