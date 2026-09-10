# Multi-allelic calling (real data) and polyploid calling (synthetic only)

> **SUPERSEDED 2026-09-10 — numbers only.** This document predates the final
> 19-dataset sweep of 2026-09-09, and its multi-allelic rows say 11/18 vs 0/18. Every one of those figures has been
> replaced:
>
> | | superseded | **final** |
> |---|---|---|
> | Claim 1 | 14 datasets, various margins | **19/19 vs SPRING (−6.03%), 19/19 vs Genozip (−43.26%), 57/57 LOSSLESS** |
> | het-SNV F1 | 0.890 | **0.876** (mean of 4 individuals, full chr20, called from the archive) |
> | het-indel F1 | 0.637 / 0.666 | **0.621** |
> | multi-allelic | 11/18 | **21/26**, a complete chr20 census; DiscoSNP++ 0 of 3,989 records |
> | tetraploid SNV | 0.836 | **0.897** |
> | export vs SPAdes | 555–656x | **129–784x** |
> | T3.4 coordinate | 0/400 | **81/400** (the 0 was our query emitting consensus, not reads) |
>
> Citable source: `benchmark/results/*.csv`, traced in
> `benchmark/documentation/RESULT_CODE.md`. **Where this document and a result
> file disagree, the result file wins.**
>
> Kept because its *reasoning* is still the record of how the conclusion was
> reached — this repo does not delete superseded work, it marks it.

Written 2026-09-02, corrected 2026-09-03.

## 0. THE TWO CLASSES ARE NOT THE SAME — an earlier version of this document
## conflated them, and the claim must be stated separately for each

* **Multi-allelic** = ONE SITE carrying two or more alternate alleles
  (`REF=C ALT=A,G`). An ordinary **diploid** sample has these wherever the
  genotype is `1/2`: both haplotypes non-reference and different from each
  other. HG002 has **952** on chr20 alone.
* **Polyploid** = an ORGANISM with three or more haplotypes (triploid wheat,
  tetraploid potato). A different thing entirely.

They share a code path (`CAPS_PLOIDY=k` admits up to k co-occurring alleles),
which is why they were run together — but they are different claims and only
one of them has real-data evidence.

| class | data available | result | claimable? |
|---|---|---|---|
| **multi-allelic** | **REAL** (GIAB GT=1/2 sites) | **G_CAPSUL 11/18 vs DiscoSNP++ 0/18** | **YES — this is the claim** |
| polyploid | **synthetic only** | G_CAPSUL 0.981 vs DiscoSNP++ 0.994 | **NO — we lose, and it is not real data** |

**Why polyploid has no real-data evidence.** Humans are diploid, so GIAB
provides none. Pooling individuals does not create it — measured: pooling all
four GIAB samples yields **1 multi-allelic site in 763**, because human SNPs
are essentially always biallelic in the population (§1). Genuinely polyploid
organisms are sequenced (wheat, potato, strawberry), but **no GIAB-quality
truth set exists for any of them**, so accuracy cannot be scored. The only
polyploid test available is therefore synthetic — and we lose it.

**Conclusion: drop the polyploid claim. Claim multi-allelic calling on real
data instead.** It is the capability the polyploid code path actually
provides, it is measured on real reads against real truth, and DiscoSNP++
cannot do it at all.

---
## 1. Real human data CANNOT test polyploid calling — measured

The obvious real-data design is to pool GIAB individuals: N diploid people
→ 2N haplotypes, truth = the union of their alleles. It does not work, and the
reason is biological rather than technical.

Counting multi-allelic sites (>1 distinct ALT at one position) in the union
truth over chr20:3.0–3.4M:

| pool | SNV sites | **multi-allelic** |
|---|---|---|
| HG002+HG003 | 603 | **1 (0.2%)** |
| HG002+HG005 | 704 | **1 (0.1%)** |
| HG002+HG004+HG005 | 761 | **1 (0.1%)** |
| HG002+HG003+HG004+HG005 | 763 | **1 (0.1%)** |

**Human SNPs are almost always biallelic in the population**: at a given
position essentially every individual carrying a variant carries the *same*
alternate base. Adding more individuals adds more *sites*, not more *alleles
per site*. Pooling four unrelated-to-related individuals still yields exactly
one multi-allelic site.

**Consequence.** The "real-data pooled polyploid" run reported earlier
(G_CAPSUL 0.709 vs DiscoSNP++ 0.710) was **not testing polyploid calling at
all** — with 1 multi-allelic site in 603, both tools were being scored on
ordinary biallelic het-SNV calling, and the tie says nothing about the
polyploid capability. That result is withdrawn as a polyploid measurement.

## 1b. REAL-DATA multi-allelic testing IS possible — via GT=1/2 sites

§1 shows pooling cannot create multi-allelic sites. But a **single diploid
sample already contains them**: wherever the genotype is `1/2`, both haplotypes
carry a non-reference allele AND the two differ from each other, so the VCF
record is genuinely multi-allelic (`REF=C ALT=A,G`).

Counted in GIAB truth on chr20:

| sample | GT=1/2 sites on chr20 |
|---|---|
| HG002 | **952** |
| HG003 | **940** |

These are real reads, real variants, real multi-allelic truth, and a single
sample — no pooling, no simulation. They exercise exactly the capability under
test: a caller that can only ever emit ONE ALT per site gets every one of them
half-wrong.

Benchmark built on this basis (`~/pld/realmulti`): HG002, chr20:1–6 Mb, 30×,
**971,250 real reads, 5,216 truth het sites of which 111 are multi-allelic**.
This supersedes both the pooled design (§1, no multi-allelic sites) and the
synthetic design (§5, a generator whose truth had to be repaired) as the
primary evidence for this class.

## 2. What a legitimate polyploid benchmark requires

Multi-allelic sites must be common, which means one of:

1. **A genuinely polyploid organism** (tetraploid potato, hexaploid wheat).
   No GIAB-quality truth set exists for these, so accuracy cannot be scored.
2. **Pooled distinct strains/species**, where the pooled genomes differ from
   each other at many positions. Truth derivable from the reference genomes.
3. **Synthetic k-ploid data** — planted multi-allelic sites, exact truth.

Only (3) gives exact truth today, which is why it is used here and why the
outer ARCS project used it as well. It must be **labelled synthetic** in the
paper, and the human-pooling result above should be cited as the reason a real
human polyploid benchmark is not possible rather than merely absent.

## 3. Scoring rule for this class

DiscoSNP++ emits **separate biallelic records** at a multi-allelic site, plus
alt-vs-alt records whose REF does not match the reference genome. Two rules
follow:

* **Symmetric hygiene**: drop records whose REF disagrees with the reference,
  for BOTH tools (G_CAPSUL has none to drop). This is legitimate normalisation.
* **Do NOT merge one tool's records into multi-allelic form**: that is a
  transformation its own pipeline never performs, applied to only one side, and
  it inflates its score from 0.800 to 1.000. An earlier version of these
  results did exactly that and had to be retracted.

Under symmetric hygiene on synthetic triploid data: **G_CAPSUL 0.886 vs
DiscoSNP++ 0.800**.

## 4. The capability difference, stated separately from accuracy

G_CAPSUL emits **multi-allelic VCF records natively** (`CAPS_PLOIDY=k` admits up
to k co-occurring alleles). DiscoSNP++ does not, and its own paper does not
claim polyploid calling — it presents the tool for diploid and haploid data and
reports biallelic variants only. This is a capability difference, not merely an
accuracy margin, and it is the more defensible claim of the two.

---

## 5. FINAL polyploid result, after fixing the generator and the scoring

Two defects had to be removed before any number here meant anything.

**Defect 1 — the truth set was broken (fixed, commit `d92a5b1`).**
`sim_polyploid.py` mutated the reference haplotype at multi-allelic sites, so
**22 of 80 truth records (28%), all of them multi-allelic, had a REF that
disagreed with the `ref.fa` they were generated against**. A caller reporting
the correct genome base could never match them. After the fix: 0 disagreements.

**Defect 2 — the scoring convention decided the winner.** Measured three ways
on the same corrected data:

| convention | G_CAPSUL | DiscoSNP++ | who "wins" |
|---|---|---|---|
| raw, no normalisation | 0.886 | 0.800 | G_CAPSUL |
| split multi-allelics in truth + both calls | 0.658 | 0.664 | tie |
| **join to multi-allelic in truth + both calls** | **0.981** | **0.994** | **DiscoSNP++** |

The split convention is invalid for this class: decomposing a truth record
`A>C,G  GT 1/2` into two `0/1` records at one position is genotypically
impossible, and rtg then credits one and rejects the other **for both tools** —
which is exactly what the per-site inspection shows (we call both `A>C` and
`A>G`, one scores TP, the other FP).

**The joined convention is the correct one**, and under it both tools are
near-perfect with DiscoSNP++ narrowly ahead: **G_CAPSUL 0.981 vs 0.994**. Our
gap is 2 missed sites of 80; precision is equal (0.987 vs 0.988).

### Verdict

**G_CAPSUL does not win polyploid.** Every earlier margin in its favour came
from a broken truth set or a favourable representation. Stated honestly:

* **Accuracy: near-parity, DiscoSNP++ marginally ahead** (0.981 vs 0.994) on
  synthetic triploid data.
* **Capability: ours is native.** G_CAPSUL emits multi-allelic VCF records
  directly; DiscoSNP++ emits separate biallelic rows that must be joined before
  they can be scored as multi-allelic — `bcftools norm -m +any` failed on its
  output entirely, and the join had to be written by hand. Its own paper does
  not claim polyploid support.
* **Real human data cannot test this class at all** (§1): 1 multi-allelic site
  per 763, so the capability is unexercised there.

The defensible paper claim is the capability, not an accuracy win.

---

## 6. REAL-DATA RESULT — multi-allelic calling, HG002 chr20:1–6 Mb

The first polyploid/multi-allelic measurement in this project made on **real
reads with real truth**. 971,250 real Illumina reads at 30×, GIAB v4.2.1 truth
restricted to the confident regions, 4,376 truth SNV sites after joining, both
tools normalised identically (SNV-only, REF must match the genome, records
joined to multi-allelic form, DiscoSNP++'s POS off-by-one corrected).

### Overall SNV accuracy

| | TP | FP | FN | P | R | **F1** |
|---|---|---|---|---|---|---|
| **G_CAPSUL** | 3316 | 192 | 930 | 0.945 | **0.781** | **0.855** |
| DiscoSNP++ | 3247 | 112 | 999 | **0.967** | 0.765 | 0.854 |

G_CAPSUL wins on recall and on F1 (marginally), DiscoSNP++ on precision. This is
**12× more truth sites than a single 400 kb window**, so it is the most
statistically substantial head-to-head in the project.

### Multi-allelic sites — the capability under test

| | sites called of 7 | with BOTH alleles correct |
|---|---|---|
| **G_CAPSUL** | **5** | **5** |
| DiscoSNP++ | **0** | **0** |

Examples: truth `A>G,T` → G_CAPSUL `G,T`; truth `G>C,T` → G_CAPSUL `C,T`.
DiscoSNP++ returns nothing at any of them.

Record counts confirm the structural difference: G_CAPSUL emitted **8
multi-allelic records**, DiscoSNP++ emitted **0**. It cannot represent two
alternate alleles at one position, so every GT=1/2 site is out of its reach.

### Verdict for this class

**On real data G_CAPSUL wins multi-allelic calling outright (5/7 vs 0/7) while
matching DiscoSNP++ on overall SNV F1 (0.855 vs 0.854).** Unlike the synthetic
triploid comparison — where a broken generator and the choice of representation
decided the result — this is real reads, real truth, identical normalisation,
and the margin is a capability the competitor does not have rather than a
threshold artefact.

## 7. INDEPENDENT VALIDATION — chr20:6–26 Mb (non-overlapping, 4× larger)

3,931,113 real reads, 16,188 truth SNV sites after joining.

| | P | R | F1 | multi-allelic sites correct | multi-allelic records emitted |
|---|---|---|---|---|---|
| **G_CAPSUL** | 0.941 | **0.808** | 0.870 | **6 / 11** | **29** |
| DiscoSNP++ | **0.956** | 0.802 | **0.873** | **0 / 11** | **0** |

**The capability result replicates exactly**: DiscoSNP++ recovers none of the
multi-allelic sites in either region and emits zero multi-allelic records in
25 Mb of real sequence. Combined over both regions: **G_CAPSUL 11/18, DiscoSNP++
0/18.**

Overall SNV F1 is a tie at scale (0.855 vs 0.854 on 5 Mb; 0.870 vs 0.873 on
20 Mb) — G_CAPSUL consistently higher recall, DiscoSNP++ consistently higher
precision.

**Final position for this class: the win is the capability, replicated on two
independent real-data regions, with overall SNV accuracy at parity.**
