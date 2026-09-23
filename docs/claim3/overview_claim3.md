---
Date: 2026-09-19
Title: Claim 3 (ADDRESSABLE) — Overview, What It Is And Is Not
Purpose: Entry point for anyone (human or Claude) needing to understand what
  Claim 3 claims, what it does not claim, and how the five sub-tables relate
  to each other before reading any single table's detail file.
When to refer to this file: Before writing/reviewing any Claim 3 paper
  section, before answering "is Claim 3 locked," before deciding what belongs
  in the main tables vs discussion, or when confused about how T3.1-T3.5
  relate to each other.
Keywords: ADDRESSABLE, Claim 3, G_CAPSUL, pseudogenome, export, coverage,
  query, locus retrieval, exact match, coordinate-range, overview, scope
---

# Claim 3 — ADDRESSABLE. What it is, what it is not.

## What Claim 3 actually claims, in one sentence

**The same archive built for compression (Claims 1/2) can also serve as the
basis for five distinct downstream operations — export, coverage, range
query, exact-match retrieval, and locus retrieval — without rebuilding a
separate representation for each one.**

This is an architectural claim, not a "we are better at assembly / better at
search than tool X" claim. Every sub-table exists to demonstrate one specific
operation working *from the same retained archive*, not to win a competition
against a specialized tool built to do only that one thing.

## What Claim 3 is NOT

- **Not a claim that the pseudogenome is a correct standalone genome
  assembly.** T3.1's correctness numbers are honestly mixed (2 of 5-7
  QUAST metrics won, see `04_T3.1_EXPORT.md`). This was never the point —
  see `mechanism_insight_claim3.md` for why "is it a good assembly"
  is the wrong question to ask of this archive.
- **Not a claim to beat BEETL/CIndex/sFASTQ at exact-match search on their
  own terms.** T3.4 reaches the same recall (1.0000) they get by
  construction, using an auxiliary index (the completion index) built from the archive
  after compression. No speed or index-size comparison against them has
  been measured or is claimed. See `01_T3.4_EXACT_MATCH.md`.
- **Not five independent mechanisms.** Architecturally there are four real
  code paths, not five (T3.3 and T3.4 are the *same* function in
  `decoder.cpp`, differing only in how one variable gets filled; T3.5
  is not a distinct C++ path at all, it is an orchestration pattern built by
  calling T3.4's primitive twice). See `code_mapping_claim3.md` and
  the architecture diagram at `docs/CLAIM3_ARCHITECTURE_DIAGRAM.md`.

## The five tables, one line each

| Table | Question it answers | Competitor | Status |
|---|---|---|---|
| T3.1 | How fast can the pseudogenome be recovered as a standalone FASTA, and how correct is it as a genome? | SPAdes (de novo assembler) | Speed: dominant (129-784x). Correctness: 2 of ~7 QUAST metrics won (genome fraction, indels), rest disclosed as SPAdes's. |
| T3.2 | How fast can per-base coverage be computed from the archive? | bwa+samtools+mosdepth | Dominant win, 16-54x, identical numbers (not approximated). |
| T3.3 | Can an arbitrary coordinate range be retrieved directly? | None exists | Native, all 19 datasets, no baseline column applies. |
| T3.4 | Can every read containing a query string be retrieved (exact match)? | BWT-family (BEETL, CIndex, sFASTQ) | Recall 1.0000, 4 individuals, via completion index (an index derived from the archive, built after compression). |
| T3.5 | Can both alleles of a heterozygous site be recovered even though compression may have split them onto different pseudogenome coordinates? | None exists | 99.93% native (1,451/1,452), 100% with the completion index, 4 individuals, 2 loci. This is the paper's actual novel claim. |

## Why T3.1 and T3.2 matter beyond the benchmark numbers — the real applications

Verified against `docs/CLAIM3_APPLICATIONS.md`. Claim 3 is fundamentally an
**amortization argument**: the pseudogenome is built once, during
compression, so analyses that normally need a separate tool run are served
directly from stream decodes instead —
`export: literal+mem_triples -> pseudogenome (no assembly)`,
`coverage: pos_abs+read_lengths -> per-position depth (no alignment)`,
`query: pos_abs+pg -> reads in a range (no full decode)`. Each stops as
soon as its own required streams are decoded; none pays for full
reconstruction. The real-world motivation is not abstract: **T3.1 (export)**
targets settings where reference-based workflows fail outright — hospital
outbreak investigation (assembling isolates, clustering, SNP distances,
transmission chains), novel-pathogen detection, and metagenomics, all cases
where no trustworthy reference genome is guaranteed to exist for an
arbitrary clinical isolate. **T3.2 (coverage)** targets sequencing QC and
dosage decisions — clinical labs routinely certify "≥30x over 95% of
target" before a run is reportable, and coverage is also the direct input
to copy-number-variant detection. Scope is stated precisely, not implied
beyond what is true: T3.2's coverage is per-contig in *pseudogenome* space
— it answers uniformity and relative depth, sufficient for the QC
application, but not locus-level CNV without a coordinate mapping this
project does not have.

## How the tables relate to each other — the one thing to hold onto

Do not think of T3.1-T3.5 as five separate feature checkboxes. Think of it as:

1. **T3.1 and T3.2 are "can you get back what's already stored, fast."**
   Both exit the assembly-rebuild layer early or entirely; neither needs to
   search anything.
2. **T3.3 is "can you address an arbitrary span."** The simplest form of the
   archive's internal coordinate system working as advertised.
3. **T3.4 is "can you also do what a totally different architecture (BWT)
   does natively."** Answering yes, at a real cost (the completion index), demonstrates
   the archive is not *limited* to positional queries — it can be extended.
4. **T3.5 is "can you use the coordinate system to answer a question T3.3's
   plain version structurally cannot."** At heterozygous sites, compression
   splits one biological locus across multiple disconnected pseudogenome
   addresses — the same phenomenon that drives Claim 2's F1 improvement
   (0.431 → 0.888 when alleles are reconciled instead of read as one
   coordinate). T3.5 is that same insight, applied to retrieval instead of
   calling. **This is the sentence that ties Claim 3 back to the paper's
   actual thesis** — expanded fully in `mechanism_insight_claim3.md`.

## File index for this folder

- `t31_export_claim3.md`
- `t32_coverage_claim3.md`
- `t33_range_query_claim3.md`
- `t34_exact_match_claim3.md`
- `t35_locus_retrieval_claim3.md`
- `mechanism_insight_claim3.md` — how Claim 3 connects to the paper's
  central insight; honest, verified (not guessed) comparison to PgRC2 and
  BEETL on whether they could replicate this with a small addition.
- `code_mapping_claim3.md` — exact function/line references in
  `src/decoder.cpp` for every operation above, and exactly which
  benchmark/results files and docs contain the numbers cited throughout this
  folder.

## Source-of-truth rule for this whole folder

Where any file in this folder disagrees with `results/`, the
result file wins. Where it disagrees with `docs/T3.1_CORRECTNESS_FINAL_20260919.md`
or `docs/T34_T35_MULTI_INDIVIDUAL_20260919.md`, those (dated, more recent,
more thoroughly gated) win. This folder is a *navigation and synthesis*
layer over those sources, not a replacement source of truth.
