# Claim 2 — results after the substrate rebuild

Supersedes the "we lose" conclusion in `docs/CLAIM2_RESULTS.md` for SNVs.
Everything here is `rtg vcfeval` (GA4GH), het-restricted, inside the GIAB
confident regions, on identical reads for every tool, with DiscoSNP++ rerun
per window rather than quoted (and its POS off-by-one corrected — see
`docs/DISCOSNP_INTERNALS.md`, without that correction it scores 0.004).

## Headline — HG002 r2 (the designated tuning window)

| tool | SNV P | SNV R | **SNV F1** | INDEL F1 |
|---|---|---|---|---|
| **CAPSULE** | 0.953–0.976 | 0.805 | **0.873–0.882** | 0.453 |
| DiscoSNP++ | 0.971 | 0.740 | **0.840** | 0.679 |

**SNV: we now beat DiscoSNP++** (0.88 vs 0.84), driven by higher recall at
comparable precision. **INDEL: we still lose** (0.45 vs 0.68).

Progression on this window: **0.419 → 0.541 → 0.821 → 0.882**.

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

## Honest standing

- **SNV: dominant** on the tuning window (0.88 vs 0.84). Held-out and
  cross-individual validation is the number that decides it — see the
  generalization section once complete.
- **INDEL: behind** (0.45 vs 0.68). Our indels come only from contig bubbles;
  DiscoSNP++ detects them in a purpose-built graph and is more sensitive.
  This is the open gap and it is not hidden.
- The measured recall ceiling is **0.995**, and we are at 0.805 — the
  remaining headroom is real, not exhausted.
