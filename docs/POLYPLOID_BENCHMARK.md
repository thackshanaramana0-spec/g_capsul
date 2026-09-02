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
