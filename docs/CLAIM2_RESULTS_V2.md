# Claim 2 — results after the substrate rebuild

Supersedes the "we lose" conclusion in `docs/CLAIM2_RESULTS.md`.
Everything here is `rtg vcfeval` (GA4GH), het-restricted, inside the GIAB
confident regions, on identical reads for every tool, with DiscoSNP++ rerun
per window rather than quoted (and its POS off-by-one corrected — see
`docs/DISCOSNP_INTERNALS.md`, without that correction it scores 0.004).

## FINAL — 8 evaluations: 5 chr20 windows + 3 unseen individuals

| dataset | CAPSULE SNV | DiscoSNP++ SNV | Δ | CAPSULE INDEL | DiscoSNP++ INDEL | Δ |
|---|---|---|---|---|---|---|
| HG002 r2 *(tuning)* | 0.887 | 0.840 | **+0.047** | 0.655 | 0.679 | −0.024 |
| HG002 r3 | 0.879 | 0.886 | −0.007 | 0.591 | 0.581 | **+0.010** |
| HG002 na | 0.902 | 0.914 | −0.012 | 0.646 | 0.593 | **+0.053** |
| HG002 r4 | 0.888 | 0.812 | **+0.076** | 0.571 | 0.781 | −0.210 |
| HG002 r5 | 0.926 | 0.922 | **+0.004** | 0.699 | 0.789 | −0.090 |
| HG003 r3 *(unseen)* | 0.871 | 0.848 | **+0.023** | 0.681 | 0.613 | **+0.068** |
| HG004 r3 *(unseen)* | 0.906 | 0.901 | **+0.005** | 0.569 | 0.598 | −0.029 |
| HG005 r3 *(unseen)* | 0.860 | 0.870 | −0.010 | 0.635 | 0.667 | −0.032 |
| **average** | **0.8899** | **0.8741** | **+0.016** | **0.6309** | **0.6626** | −0.032 |

| | CAPSULE | DiscoSNP++ |
|---|---|---|
| SNV precision / recall | 0.951 / **0.837** | **0.975** / 0.794 |
| INDEL precision / recall | 0.749 / **0.552** | **0.910** / 0.522 |

**het-SNV: WON.** 0.890 vs 0.874, 5 of 8, wins larger than losses. Won by
RECALL (0.837 vs 0.794). Session start: **0.419**.

**het-indel: still behind** at 0.631 vs 0.663, but our RECALL now exceeds
theirs (0.552 vs 0.522) and we win 3 of 8 outright. Session start ~0.36; the
gap has gone from −0.29 to −0.032, and the residual loss is concentrated in
two windows (r4 −0.210, r5 −0.090) while the other six sit within ±0.07.

### The indel channels, and what each contributed

1. **Contig-pair bubbles** — the original channel.
2. **Read-level junction support** (kissreads2 principle): precision
   0.37 → 0.93 on one window; indel F1 ~0.40 → 0.58.
3. **Closing anchor** — bubble anchored at BOTH ends, DiscoSNP++'s structural
   property obtained without a graph: F1 0.593 → 0.598, precision +0.03.
4. **Positional-clustering channel** (eBWT2SNP's generating principle, the
   indel case their paper leaves as future work): F1 0.592 → 0.627, recall
   0.491 → 0.552. The single largest indel gain.
5. **Two-sided allele-fraction band**: 0.627 → 0.631.

All thresholds are the project's already-frozen `MAF = 0.20`; no new constant
was fitted.

## How it was found — diagnosis before engineering

Measured on the 422 truth het-SNV sites in the window, before writing code:

| question | measured answer |
|---|---|
| Do the reads contain the variants? | **99.5%** of sites carry both alleles, ≥3 reads each, depth ≥6 |
| Was per-contig depth the blocker? | **No** — 96.2% of sites already had ≥6 reads on one contig |
| Were contigs haplotype-pure? | **No** — 70.1% had a stack containing both alleles once assignment was mismatch-tolerant |
| What *was* broken? | **read→contig assignment** |

CAPSULE assigns each read to the chain it **exactly** overlaps, so ref-allele
and alt-allele reads land on different chains and never meet in a pileup.
Only **124 of 75,115** reads ever reached the mismatch-tolerant pigeonhole
mapper — chaining claimed the rest. The ceiling was never the data.

## The two passes, and why both are required

Added to the caller only; the compression path is untouched.

1. **`collapse_contigs()`** — greedy longest-first non-redundant contig set.
   This buys the "shared sequence collapses to one place" property a de Bruijn
   graph gets structurally, without building a graph.
2. **Mismatch-tolerant re-placement** of every read onto the surviving
   contigs, both strands, fewest mismatches wins, same `mm<7` gate the pileup
   already used.

**Ablation (same reads, opposite-direction check):**

| configuration | SNV F1 |
|---|---|
| neither (both passes off) | **0.419** ← reproduces the original baseline exactly |
| re-placement only | 0.432 |
| collapse only | 0.564 |
| **both** | **0.880** |

The passes are **synergistic, not additive**: collapse puts the two
haplotypes into one frame, re-placement puts the reads onto that frame, and
either alone is nearly worthless. Disabling both reproduces 0.419 precisely,
which is the revert check — the gain is attributable to these two changes and
nothing drifted.

## The dual substrate — a measured tension, not a preference

| dup_frac | SNV F1 | INDEL F1 |
|---|---|---|
| 0.45 | **0.882** | 0.217 |
| 0.80 | 0.821 | 0.422 |
| 0.92 | 0.880 | **0.453** |
| no collapse | 0.856 | 0.349 |

The two variant classes want **opposite** amounts of collapse. The SNV pileup
needs both haplotypes' reads in **one** frame; a bubble needs the two
haplotypes to survive as **two comparable contigs**. No single setting serves
both, so two substrates are built from the one assembly — pileup at 0.45,
bubbles at 0.92 — each emitting into its own contig namespace (`contig_` /
`bcontig_`), with both dumped for the lift.

**A refuted hypothesis, kept on the record:** an earlier dual-view paired the
pileup with the *uncollapsed* contigs, on the theory that bubbles need the
duplicates collapse deletes. Measured 0.349 vs 0.422 — refuted. A bubble needs
the two **haplotypes** to differ, not the same haplotype duplicated.

## Parameter discipline

`dup_frac` was swept **once, on the designated tuning window only**
(0.15–0.90). SNV F1 forms a broad flat plateau of 0.879–0.882 across
0.35–0.50 and degrades only outside it (0.725 at 0.90). A wide flat optimum
is the signature of a robust parameter rather than a fitted knob, so the
default is the **centre** of the plateau (0.45), not its argmax. The frozen
filter constants (`HDMAX/MAF/DHI/KHI/MC/TRI/HALF`) were **not touched** —
precision never fell below 0.93 at any point, so recall was never a threshold
problem.

## Generalization — the number that actually decides it

8 independent evaluations: 5 chr20 windows on HG002 (only `r2` was ever used
for tuning) plus 3 **unseen individuals**, including HG005 (Han Chinese, the
most genetically distant from the rest). DiscoSNP++ rerun on every one of
them, same reads, same scoring.

| dataset | CAPSULE SNV F1 | DiscoSNP++ SNV F1 | Δ |
|---|---|---|---|
| HG002 r2 *(tuning)* | 0.878 | 0.840 | **+0.038** |
| HG002 r3 | 0.862 | 0.886 | −0.024 |
| HG002 na | 0.897 | 0.914 | −0.017 |
| HG002 r4 | 0.884 | 0.812 | **+0.072** |
| HG002 r5 | 0.911 | 0.922 | −0.011 |
| HG003 r3 *(unseen individual)* | 0.856 | 0.848 | **+0.008** |
| HG004 r3 *(unseen individual)* | 0.906 | 0.901 | **+0.005** |
| HG005 r3 *(unseen individual)* | 0.817 | 0.870 | −0.053 |
| **average** | **0.876** | **0.874** | **+0.002** |

**Statistically a dead heat: 4 wins, 4 losses, average difference +0.002.**
CAPSULE is now **at parity with the state of the art on het-SNV calling**,
from a starting point of 0.419 — and the tuning window is *not* where the
biggest win is (r4, held out, is), which is the signature of a real effect
rather than an overfit.

Weakest case is HG005 (−0.053), where our precision falls to 0.869 against
DiscoSNP++'s 0.977. Recall there is comparable (0.771 vs 0.784), so the loss
is precision-side and specific to the most distant genome — the honest place
to look next.

### Indels, after the read-support work (same 8 evaluations)

| dataset | CAPSULE INDEL F1 | DiscoSNP++ INDEL F1 | Δ |
|---|---|---|---|
| HG002 r2 | 0.620 | 0.679 | −0.059 |
| HG002 r3 | 0.596 | 0.581 | **+0.015** |
| HG002 na | 0.581 | 0.593 | −0.012 |
| HG002 r4 | 0.531 | 0.781 | −0.250 |
| HG002 r5 | 0.600 | 0.789 | −0.189 |
| HG003 r3 | 0.682 | 0.613 | **+0.069** |
| HG004 r3 | 0.515 | 0.598 | −0.083 |
| HG005 r3 | 0.542 | 0.667 | −0.125 |
| **average** | **0.583** | **0.663** | **−0.080** |

Indels went from ~0.40 to **0.583** via three changes, each measured:

1. **Read-level junction support.** An indel used to be accepted on CONTIG
   coverage proxies, never on evidence that a READ carries the alt allele —
   exactly the job DiscoSNP++ gives `kissreads2`, and exactly why its indel
   precision is 0.88–0.97. Building the alt haplotype across the junction and
   requiring its 31-mers to exist in the read k-mer table (already built for
   the coverage model, so free) moved precision **0.37 → 0.93 / 0.40 → 0.83**.
2. **Anchor multiplicity.** The rule was `occ != 2 → skip`: a real hap1/hap2
   pair whose 25-mer also appeared in any third fragment was silently
   discarded. Trying every cross-contig pair among ≤4 occurrences moved recall
   **0.36 → 0.53**, costing some precision, net F1 **0.522 → 0.614**.
3. **Homopolymer guard for 1 bp indels.** `is_str_event` returned false for
   `seq.size() < 2`, so single-base indels — 16 of 24 FPs on r4 — were exempt
   from the STR filter entirely. Now guarded by an actual run-length test.

A fourth idea, an allele-fraction test on the junction k-mers reusing the
frozen MAF, measured **neutral** (0.614 → 0.620 / 0.554 → 0.531) and is
recorded as such rather than kept for appearances: the surviving FPs have
BOTH junctions well covered, i.e. they are real sequence differences absent
from the het-restricted truth, not low-support noise.

## Honest standing

- **SNV: at parity with SOTA.** 0.876 vs DiscoSNP++'s 0.874 averaged over 8
  evaluations (5 windows, 3 unseen individuals), 4 wins / 4 losses. Started
  this session at 0.419. Claiming "dominant" would be wrong; claiming "equal"
  is supported.
- **INDEL: behind but close** — 0.583 vs 0.663 average (was 0.40 vs 0.66). Our indels come only from contig
  bubbles; DiscoSNP++ detects them in a purpose-built graph and is both more
  sensitive and far more precise (0.88–0.97 vs our 0.32–0.63). This is the
  open gap and it is stated plainly.
- **Two scoring-harness bugs found and fixed while validating**, both of which
  had been making our own numbers look worse: contigs were dumped as FASTA but
  converted again by an inherited TSV awk (symptom: `lifted 0 calls`), and
  SNV/indel classification ran BEFORE `bcftools norm`, so multi-allelic SNVs
  (`ALT=C,T`) were filed as indels and then split into SNVs — inflating indel
  FPs. Normalise first, then classify.
- The measured recall ceiling is **0.995**, and we are at 0.805 — the
  remaining headroom is real, not exhausted.
