# Design principle and mechanism

This document explains why one archive can serve compression, variant calling and locus retrieval.
It states the design principle (ARCH), the mechanism that makes reconciliation necessary, and how
each capability applies it. The per-claim documents give the measurements.

## The design principle

A reference-free compressor has to work out how the reads fit together before it can store them
compactly. Most compressors discard that knowledge once the archive is written. G_CAPSUL keeps it.
The design principle is ARCH:

| Step | What happens |
|---|---|
| **A**ssemble | Reads are chained and mapped into a pseudogenome, a single sequence built from the reads themselves. |
| **R**etain | The placement, orientation and length of every read, and its differences from the pseudogenome, are stored in the archive. |
| **C**ompress | The pseudogenome and the retained information are entropy-coded as separate streams. |
| **H**arness | The retained structure is reused for analysis, so that downstream tasks do not rebuild it from the FASTQ. |

Because the coordinate system is produced while compressing, using it costs no extra storage.
Measured on the benchmark, recording the contig boundaries needed for per-contig export adds 0 bytes
to the archive (18,282,397 bytes with and without them, E. coli SRR2584863).

## The mechanism: separated alleles

A size-minimizing compressor that builds a shared reference and codes each read as a position plus
mismatches is pushed to separate the two alleles of a heterozygous site. Beyond some divergence it is
cheaper to give each allele its own placement than to keep one placement with a growing list of
mismatches. This follows from the design goal, not from a flaw in one implementation. The full
argument is in [`claim3/mechanism_insight_claim3.md`](claim3/mechanism_insight_claim3.md).

The consequence is that reads are preserved exactly, but the fact that two groups of reads are
alternative forms of one biological locus is no longer visible in the coordinates. Counting
variant evidence or fetching the reads at a locus from the stored coordinates alone therefore
misses part of it.

G_CAPSUL reconciles the separated regions during analysis. Reconciliation changes nothing in the
stored archive, so it adds no archive-size cost.

## How each capability applies the mechanism

| Capability | What goes wrong without reconciliation | What reconciliation does | Evidence |
|---|---|---|---|
| Variant calling | Allele evidence is split across regions, so the caller sees few supporting reads at a het site. | Related regions are collapsed into a shared frame and reads are placed again against it. | Chromosome 20 ablation on HG002: F1 0.431 with neither step and 0.888 with both. Without the collapse, precision stays at 0.96 while recall falls to 0.28. |
| Multi-allelic and polyploid calling | The same split, with more than two alleles. | The same reconciliation, generalised to more alleles and to higher ploidy. | 17 of 26 multi-allelic sites on chromosome 20 and tetraploid SNV F1 of 0.897. |
| Locus retrieval | A single stored coordinate reaches only one allele's reads at most sites. | Probes are anchored upstream and downstream, and the two read sets are combined. | 1,451 of 1,452 sites by position, and 1,452 of 1,452 with the FM-index. |

Claim 1 is where the separation arises. The same overlap chaining that gives the compression ratio
is what places divergent alleles apart, and the placement and deviation data it computes is what
Claims 2 and 3 reuse.

## A related observation: why the caller does not use a de Bruijn graph

Building a de Bruijn graph directly over the pseudogenome would let a bubble-based caller such as
DiscoSNP++ run on it. This was measured and does not work. Greedy overlap chaining linearizes the
data, so 97.7% of true-positive anchors have an out-degree of 1 and the branching structure a
bubble caller needs is gone. The result is detailed in
[`claim2/architecture_vs_discosnp_claim2.md`](claim2/architecture_vs_discosnp_claim2.md). It shows
the same tension from the other side: the structure that compresses well is not the structure that
exposes variants, which is why reconciliation is needed.

## One archive, several operations

The retrieval operations are not independent features. Range retrieval and exact-match retrieval
share one primitive in `src/decoder.cpp`, and locus retrieval calls that primitive twice (upstream
and downstream) and merges the results. See
[`claim3/code_mapping_claim3.md`](claim3/code_mapping_claim3.md) for the function-level mapping.

Reconciliation also happens at different times. For variant calling it is part of building the
caller's substrate, before any call is scored. For locus retrieval it is applied at query time to
a finished archive. The same correction therefore works either at analysis construction or after
the fact.

## What the mechanism does not claim

- It does not claim that other compressors split alleles in the same way. That is a prediction
  from the design goal, and it was not measured on other tools.
- It does not claim to replace dedicated assemblers or callers across all conditions. The
  limitations of each capability are in [`limitations_and_scope.md`](limitations_and_scope.md).
