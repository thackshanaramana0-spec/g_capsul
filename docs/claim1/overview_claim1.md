---
Date: 2026-09-19
Title: Claim 1 (COMPACT) — Overview, What It Is And Is Not
Purpose: Entry point for understanding what Claim 1 claims, what it does
  not, and how T1.1/T1.2 relate, before reading either table's detail file.
When to refer to this file: Before writing/reviewing any Claim 1 paper
  section; answering "is Claim 1 locked"; understanding why Claim 1 is the
  foundation every other claim's cost gets measured against.
Keywords: COMPACT, Claim 1, compression, lossless, SPRING, Genozip,
  archive size, wall time, 19 datasets, overview, scope, foundation
---

# Claim 1 — COMPACT. What it is, what it is not.

## What Claim 1 actually claims, in one sentence

**G_CAPSUL produces smaller lossless FASTQ archives than SPRING and Genozip
across 19 real datasets spanning seven kingdoms, while also retaining
enough structure to support the downstream operations Claims 2 and 3
measure — and this is stated as the load-bearing constraint on everything
else in the paper, not just the first result.**

## Why this claim is the foundation, not just "first"

Every other claim in this paper is evaluated against a compressed
representation that this claim's own discipline protects. The project's
own standing rule, verified in this session's Claim 3 work: *"Claim 1 is
the stronger claim and must not be spent buying Claim 3 a constant"*
(directly quoted from code comments found this session, `decoder.cpp`
line ~511). Concretely, this was tested and enforced, not just stated —
when this session's own investigation drifted toward retuning the encoder
to chase a Claim 3 (T3.1) metric, it was correctly stopped because doing so
would have cost ~11.6% archive size. That the project has a real, working
discipline of protecting Claim 1's numbers from being spent elsewhere is
itself evidence worth citing.

## What Claim 1 is NOT

- **Not a claim to beat every compressor that exists.** SPRING and Genozip
  are the two general-purpose, widely-used FASTQ compressors chosen as
  comparators — not a hand-picked weak set. PgRC2 is compared separately
  and explicitly flagged as not directly comparable (see below).
- **Not measured on a convenient subset.** 19 datasets, seven kingdoms,
  chosen to include human, bacteria, archaea, protista, fungi, virus, and
  plantae — the manuscript states this explicitly as a deliberate test of
  whether the method holds under "substantially different sequencing
  conditions," not cherry-picked for favorable margins.
- **Not the same comparison as the PgRC2 numbers cited elsewhere in the
  project's history.** PgRC2 stores sequence only — no read names, no
  quality scores, no line-3 mode byte. Comparing it to G_CAPSUL's
  whole-FASTQ archive size would be comparing different objects. Where
  PgRC2 is cited, it must be on the sequence-stream-only subset, separately
  stated — see `mechanism_insight_claim1.md`.

## The two tables, one line each

| Table | Question | Comparators | Result |
|---|---|---|---|
| T1.1 | Lossless whole-FASTQ archive size, 19 datasets | SPRING, Genozip | 19/19 wins vs both. Aggregate: −6.03% vs SPRING, −43.26% vs Genozip. |
| T1.2 | Wall time (compress + decompress), same 19 datasets | SPRING, Genozip | **CAPSULE is 0/19 fastest at compress and 0/19 fastest at decompress** — verified by independent recomputation from the raw CSV this session, not assumed. This is a real, disclosed trade-off (smaller archives, slower construction), not a mixed result. |

## Verification status of this overview and its 4 companion files

The T1.1 aggregate figures above were **independently recomputed from the
raw CSV this session** (`results/claim1/claim1_T1.1_T1.2.csv`), not just
cited from `BTR_NOTES.md` — summed all 19 datasets' archive bytes for each
tool and derived the percentage margins from scratch. Result: exact match
to `BTR_NOTES.md`'s stated headline (19/19 vs SPRING −6.03%, 19/19 vs Genozip
−43.26%, 57/57 LOSSLESS). One individual row (ERR5181310) additionally
cross-checked byte-for-byte against the manuscript table. This is the
strongest verification standard applied across any of the three claim
folders in `refer_paper_docs/` — re-derivation, not just cross-checking.

## How the tables relate to each other

1. **T1.1 is the primary claim — size.** This is what every downstream
   claim's cost gets measured against (see the `CAPS_SPANS`/`contig_spans`
   free-byproduct discussion in `mechanism_insight_claim1.md`, which
   directly connects to Claim 3's T3.1 mechanism).
2. **T1.2 is the honest secondary axis — speed, and it is a real loss, not
   a mixed result.** CAPSULE is fastest on 0 of 19 datasets for compress and
   0 of 19 for decompress — SPRING and Genozip trade the "fastest" bold
   mark between themselves, CAPSULE never takes it. This is the direct,
   disclosed cost of the extra structure Claim 1 retains and Claims 2/3
   reuse — smaller, richer archives cost more time to build and unpack.
   State this plainly as a real trade-off, not softened as "mixed" or
   implied to be occasionally competitive.

## File index for this folder

- `t11_archive_size_claim1.md`
- `t12_wall_time_claim1.md`
- `names_and_quality_claim1.md` — the two non-DNA-sequence columns: the
  names tokenizer (incl. this project's own `ID_ZDELTA`/`ID_SEQLEN`
  additions), why quality is vendored (fqzcomp) rather than reimplemented,
  and the one-byte line-3 encoding.
- `mechanism_insight_claim1.md` — how Claim 1's own compression mechanism
  (overlap chaining, position+deviation encoding) is the structural cause of
  the allele-splitting phenomenon Claims 2 and 3 build on; honest PgRC2
  sequence-only comparison.
- `code_mapping_claim1.md` — exact stage/function references in
  `src/encoder.cpp` and the shared coder headers, and exact source
  file for every number cited.

## Source-of-truth rule for this whole folder

Where any file here disagrees with `results/`, the result file
wins. This folder is navigation/synthesis, not a replacement source of
truth.
