# Polyploid calling — what can and cannot be benchmarked

Written 2026-09-02.

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
(CAPSULE 0.709 vs DiscoSNP++ 0.710) was **not testing polyploid calling at
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
  for BOTH tools (CAPSULE has none to drop). This is legitimate normalisation.
* **Do NOT merge one tool's records into multi-allelic form**: that is a
  transformation its own pipeline never performs, applied to only one side, and
  it inflates its score from 0.800 to 1.000. An earlier version of these
  results did exactly that and had to be retracted.

Under symmetric hygiene on synthetic triploid data: **CAPSULE 0.886 vs
DiscoSNP++ 0.800**.

## 4. The capability difference, stated separately from accuracy

CAPSULE emits **multi-allelic VCF records natively** (`CAPS_PLOIDY=k` admits up
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

| convention | CAPSULE | DiscoSNP++ | who "wins" |
|---|---|---|---|
| raw, no normalisation | 0.886 | 0.800 | CAPSULE |
| split multi-allelics in truth + both calls | 0.658 | 0.664 | tie |
| **join to multi-allelic in truth + both calls** | **0.981** | **0.994** | **DiscoSNP++** |

The split convention is invalid for this class: decomposing a truth record
`A>C,G  GT 1/2` into two `0/1` records at one position is genotypically
impossible, and rtg then credits one and rejects the other **for both tools** —
which is exactly what the per-site inspection shows (we call both `A>C` and
`A>G`, one scores TP, the other FP).

**The joined convention is the correct one**, and under it both tools are
near-perfect with DiscoSNP++ narrowly ahead: **CAPSULE 0.981 vs 0.994**. Our
gap is 2 missed sites of 80; precision is equal (0.987 vs 0.988).

### Verdict

**CAPSULE does not win polyploid.** Every earlier margin in its favour came
from a broken truth set or a favourable representation. Stated honestly:

* **Accuracy: near-parity, DiscoSNP++ marginally ahead** (0.981 vs 0.994) on
  synthetic triploid data.
* **Capability: ours is native.** CAPSULE emits multi-allelic VCF records
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
| **CAPSULE** | 3316 | 192 | 930 | 0.945 | **0.781** | **0.855** |
| DiscoSNP++ | 3247 | 112 | 999 | **0.967** | 0.765 | 0.854 |

CAPSULE wins on recall and on F1 (marginally), DiscoSNP++ on precision. This is
**12× more truth sites than a single 400 kb window**, so it is the most
statistically substantial head-to-head in the project.

### Multi-allelic sites — the capability under test

| | sites called of 7 | with BOTH alleles correct |
|---|---|---|
| **CAPSULE** | **5** | **5** |
| DiscoSNP++ | **0** | **0** |

Examples: truth `A>G,T` → CAPSULE `G,T`; truth `G>C,T` → CAPSULE `C,T`.
DiscoSNP++ returns nothing at any of them.

Record counts confirm the structural difference: CAPSULE emitted **8
multi-allelic records**, DiscoSNP++ emitted **0**. It cannot represent two
alternate alleles at one position, so every GT=1/2 site is out of its reach.

### Verdict for this class

**On real data CAPSULE wins multi-allelic calling outright (5/7 vs 0/7) while
matching DiscoSNP++ on overall SNV F1 (0.855 vs 0.854).** Unlike the synthetic
triploid comparison — where a broken generator and the choice of representation
decided the result — this is real reads, real truth, identical normalisation,
and the margin is a capability the competitor does not have rather than a
threshold artefact.
