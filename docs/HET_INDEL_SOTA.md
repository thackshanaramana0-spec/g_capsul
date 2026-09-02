# Reference-free het-indel calling — who the SOTA actually is

Literature survey 2026-09-02, done to establish what our 0.631 should be
compared against. Result: **the competitive field for this exact task is very
small, and our 0.631 vs DiscoSNP++'s 0.663 is a comparison against the
strongest applicable tool.**

## 1. The task, stated precisely

Heterozygous indels, from **ONE diploid sample**, with **no reference genome**.
Each qualifier removes most of the literature.

## 2. Who does NOT qualify, and why

| tool / class | reported indel performance | why it does not apply |
|---|---|---|
| **DeepVariant, DRAGEN, Clair3, Octopus, Strelka2, GATK, FreeBayes** | indel F1 **0.89–0.97** (DRAGEN 96.99%, VarSome/Sentieon 89–93% on GIAB WES) | **reference-based.** They align to a known genome. Not comparable — different problem. |
| **ska lo** (MBE 2025, newest reference-free tool) | 97–98% insertion sensitivity | **haploid bacteria, between-sample.** Its own paper: "identify within-strain variants in pathogen WGS data" by *comparing samples*. It does not call heterozygous variants from a single diploid sample. |
| **eBWT2SNP** | 99.13% precision | **SNPs only.** Its paper leaves indel typing explicitly unimplemented ("extract the left-context ... perform a local alignment" — future work). |
| **Kmer2SNP** | — | **SNPs only**, by construction (heterozygous k-mer pair matching). |
| **SKA (base)** | — | ignores indels entirely; `ska lo` exists to add them. |

## 3. Who DOES qualify

**DiscoSNP++ is effectively the only general-purpose reference-free caller
that reports heterozygous indels from a single diploid sample.** That is why
this project benchmarks against it, and it means our 0.631 vs its 0.663 is a
comparison against the applicable state of the art, not against a weak baseline.

Their own published indel numbers are worth recording: on **simulated** human
chr1, 72% recall / 71.2% precision unfiltered (F1 ≈ 0.716), or 99.25%
precision at 58% recall when filtered (F1 ≈ 0.732). We measure it at **0.663
on real GIAB data** — the usual simulated-to-real drop, and the reason our own
numbers must never be compared to anyone's published simulated figures.

## 4. What this means for the claim

* **0.631 vs 0.663 is a 4.8% relative gap against the only applicable
  competitor**, on real data, scored identically. That is a competitive
  result, not a weak one.
* Both tools sit far below reference-based callers (0.89–0.97). That gap is
  the price of needing no reference, and it is a property of the problem, not
  of either implementation.
* The class we lose — 1 bp indels in 7–8 bp homopolymers — is the class GIAB's
  own stratification resource (Nat Commun 2024) identifies as the hardest
  standard stratum, scoring worse than low-mappability or GC-extreme regions.
* Our complementarity result stands: on HG002 r2 we find **38 true indels to
  DiscoSNP++'s 26** with fewer false positives, and 14 of our 21 exclusive
  finds are ≥2 bp, while all 9 of theirs are 1 bp homopolymer events.

**Honest framing for the paper:** "competitive with the only applicable
reference-free het-indel caller, with complementary strengths by event length"
— not "state of the art", and not a loss to be buried.
