# Claim 3 (ADDRESSABLE): operations served from the archive

This folder documents Claim 3 of the manuscript. The archive built for compression (Claims 1 and 2)
can serve five further operations, which are export, coverage, range query, exact-match retrieval
and locus retrieval, without a separate representation for each. This file states what the claim
covers and does not cover and how its five tables (T3.1 to T3.5) relate. The result files in
`results/claim3/` are authoritative.

The claim concerns the architecture. Each table shows one operation working from the same retained
archive, and none is meant to beat a tool built for that single operation.

## What Claim 3 does not claim

- **It does not claim that the pseudogenome is a correct standalone assembly.** The T3.1
  correctness numbers are mixed (G_CAPSUL is ahead on 2 of about 7 QUAST metrics, see
  [`t31_export_claim3.md`](t31_export_claim3.md)). That was not the purpose, and
  [`mechanism_insight_claim3.md`](mechanism_insight_claim3.md) explains why "is it a good
  assembly" is the wrong question for this archive.
- **It does not claim to beat BEETL, CIndex or sFASTQ at exact-match search.** T3.4 reaches the
  recall of 1.0000 that they have by construction, using an FM-index over only the reads that
  deviate from the pseudogenome, built from the archive after compression. No speed or index-size
  comparison with them was measured, and none is claimed (see
  [`t34_exact_match_claim3.md`](t34_exact_match_claim3.md)).
- **It does not describe five independent mechanisms.** There are four code paths. T3.3 and T3.4
  are the same function in `decoder.cpp`, differing in how one variable is filled, and T3.5 is not
  a separate C++ path. It is an orchestration pattern that calls the T3.4 primitive twice (see
  [`code_mapping_claim3.md`](code_mapping_claim3.md)).

## The five tables, one line each

| Table | Question it answers | Competitor | Status |
|---|---|---|---|
| T3.1 | How fast can the pseudogenome be recovered as a standalone FASTA, and how correct is it as a genome? | SPAdes (de novo assembler) | Speed: 129 to 784 times faster. Correctness: ahead on 2 of about 7 QUAST metrics (genome fraction, indels), and SPAdes is ahead on the rest. |
| T3.2 | How fast can per-base coverage be computed from the archive? | bwa+samtools+mosdepth | 16 to 54 times faster, with identical numbers (not approximated). |
| T3.3 | Can an arbitrary coordinate range be retrieved directly? | None exists | Native, all 19 datasets, no baseline column applies. |
| T3.4 | Can every read containing a query string be retrieved (exact match)? | BWT-family (BEETL, CIndex, sFASTQ) | Recall 1.0000, 4 individuals, via an FM-index derived from the archive and built after compression. |
| T3.5 | Can both alleles of a heterozygous site be recovered even though compression may have split them onto different pseudogenome coordinates? | None exists | 1,452/1,452 (100%) with the FM-index (1,451/1,452, 99.93%, by position alone), 4 individuals, 1 locus. This is the paper's central claim. |

## Why T3.1 and T3.2 matter beyond the benchmark numbers

Claim 3 is an amortization argument. The pseudogenome is built once, during compression, so analyses
that normally need a separate tool run are served by decoding only the streams they need:

- export: `literal` and `mem_triples` give the pseudogenome, with no assembly,
- coverage: `pos_abs` and `read_lengths` give per-position depth, with no alignment,
- query: `pos_abs` and the pseudogenome give the reads in a range, with no full decode.

Each stops as soon as its own streams are decoded. **T3.1 (export)** is relevant where reference-based
workflows fail, such as outbreak investigation (assembling isolates, clustering, SNP distances and
transmission chains), novel-pathogen detection and metagenomics, where no trustworthy reference
may exist for an isolate. **T3.2 (coverage)** is relevant to sequencing quality control and dosage
decisions, since clinical labs certify a minimum depth over a target before a run is reportable,
and since coverage is the input to copy-number detection. T3.2's coverage is per contig in
pseudogenome space. It answers uniformity and relative depth, which suffices for quality control,
and it does not give locus-level copy number without a coordinate mapping that G_CAPSUL does not
have.

## How the tables relate

The five tables are one progression and not five separate features:

1. **T3.1 and T3.2 are "can you get back what's already stored, fast."**
   Both exit the assembly-rebuild layer early or entirely; neither needs to
   search anything.
2. **T3.3 is "can you address an arbitrary span."** The simplest form of the
   archive's internal coordinate system working as advertised.
3. **T3.4 asks whether the archive can also do what a BWT architecture does natively.** The
   answer is yes, at a real cost (the FM-index), so the archive is not limited to positional
   queries and can be extended.
4. **T3.5 is "can you use the coordinate system to answer a question T3.3's
   plain version structurally cannot."** At heterozygous sites, compression
   splits one biological locus across multiple disconnected pseudogenome
   addresses — the same phenomenon that drives Claim 2's F1 improvement
   (0.431 to 0.888 when alleles are reconciled instead of read as one
   coordinate). T3.5 applies the same insight to retrieval instead of
   calling, which ties Claim 3 to the central argument of the paper. It is expanded in
   [`mechanism_insight_claim3.md`](mechanism_insight_claim3.md).

## File index for this folder

- [`t31_export_claim3.md`](t31_export_claim3.md): export against SPAdes.
- [`t32_coverage_claim3.md`](t32_coverage_claim3.md): coverage against bwa, samtools and mosdepth.
- [`t33_range_query_claim3.md`](t33_range_query_claim3.md): coordinate-range query.
- [`t34_exact_match_claim3.md`](t34_exact_match_claim3.md): exact-match retrieval.
- [`t35_locus_retrieval_claim3.md`](t35_locus_retrieval_claim3.md): locus retrieval.
- [`t34_t35_multi_individual_claim3.md`](t34_t35_multi_individual_claim3.md): the replication of
  T3.4 and T3.5 on HG003, HG004 and HG005.
- [`mechanism_insight_claim3.md`](mechanism_insight_claim3.md): how Claim 3 connects to the central
  insight, and whether PgRC2 or BEETL could reproduce it with a small addition.
- [`code_mapping_claim3.md`](code_mapping_claim3.md): the functions in `src/decoder.cpp` behind each
  operation, and the source file of every number.

## Source of truth

Where a file in this folder disagrees with `results/`, the result file is correct. Where it
disagrees with [`t31_export_claim3.md`](t31_export_claim3.md) or
[`t34_t35_multi_individual_claim3.md`](t34_t35_multi_individual_claim3.md), those files are
correct, since they are the more recent and more thoroughly checked. This folder explains the
results and does not replace them.
