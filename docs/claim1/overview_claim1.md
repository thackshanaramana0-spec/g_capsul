# Claim 1 (COMPACT): lossless archive size and wall time

This folder documents Claim 1 of the manuscript: G_CAPSUL produces smaller lossless FASTQ
archives than SPRING and Genozip on 19 real datasets, and it keeps enough structure in the
archive for the operations measured in Claims 2 and 3. This file states what the claim covers and
does not cover and how its two tables (T1.1 and T1.2, the first two numbered result tables of the
manuscript) relate. Tables are numbered by claim, so T1.x belong to Claim 1, T2.x to Claim 2 and
T3.x to Claim 3. The sources cited as "the development notes" are the author's internal working
record, which is not distributed with the repository. The result files in `results/` are
authoritative.

## The claim in one sentence

G_CAPSUL produces smaller lossless FASTQ archives than SPRING and Genozip across 19 real datasets
from seven kingdoms, and it retains the structure that Claims 2 and 3 use.

## Why Claim 1 constrains the others

The other claims are evaluated on the archive that this claim measures, so the archive size is not
traded away for them. A comment in `decoder.cpp` records the rule: "Claim 1 is the stronger claim
and must not be spent buying Claim 3 a constant". It was applied once in practice. A retuning of
the encoder that would have improved a T3.1 metric was stopped because it would have cost about
11.6% of the archive size.

## What Claim 1 does not claim

- **It does not claim to beat every compressor.** SPRING and Genozip are the two general-purpose
  FASTQ compressors used as comparators. PgRC2 is compared separately, on sequence only, because
  it is not directly comparable.
- **It is not measured on a convenient subset.** The 19 datasets cover human, bacteria, archaea,
  protista, fungi, virus and plants. The manuscript states that they were chosen to test the method
  under substantially different sequencing conditions, and not for favorable margins.
- **It does not mix whole-FASTQ and sequence-only numbers.** PgRC2 stores sequence only, with no
  read names, quality scores or line 3. Comparing it with G_CAPSUL's whole-FASTQ archive size would
  compare different objects, so PgRC2 appears only in the sequence-stream comparison in
  [`mechanism_insight_claim1.md`](mechanism_insight_claim1.md).

## The two tables, one line each

| Table | Question | Comparators | Result |
|---|---|---|---|
| T1.1 | Lossless whole-FASTQ archive size, 19 datasets | SPRING, Genozip | 19/19 wins vs both. Aggregate: −6.03% vs SPRING, −43.26% vs Genozip. |
| T1.2 | Wall time (compress and decompress), same 19 datasets | SPRING, Genozip | G_CAPSUL is the fastest tool on 0 of 19 datasets for compression and on 0 of 19 for decompression, recomputed from the raw CSV. This is a disclosed trade-off of smaller archives against slower construction. |

## Verification of the figures

The T1.1 aggregates were recomputed from the raw CSV (`results/claim1/claim1_T1.1_T1.2.csv`) by
summing the archive bytes of all 19 datasets for each tool. The result matches the development
notes: 19/19 against SPRING (-6.03%), 19/19 against Genozip (-43.26%) and 57/57 archives
lossless. The ERR5181310 row was also compared with the manuscript table value for value.

## How the two tables relate

1. **T1.1 (size) is the primary result.** The cost of every later claim is measured against it. The
   `CAPS_SPANS` and `contig_spans` discussion in
   [`mechanism_insight_claim1.md`](mechanism_insight_claim1.md) connects it to T3.1.
2. **T1.2 (time) is the secondary axis, and G_CAPSUL loses it.** SPRING and Genozip share the
   fastest time between them and G_CAPSUL never takes it. This is the cost of the structure
   that Claim 1 retains and that Claims 2 and 3 reuse.

## File index for this folder

- [`t11_archive_size_claim1.md`](t11_archive_size_claim1.md): archive size, 19 datasets.
- [`t12_wall_time_claim1.md`](t12_wall_time_claim1.md): wall time, 19 datasets.
- [`names_and_quality_claim1.md`](names_and_quality_claim1.md): the two non-sequence columns, which
  are the names tokenizer (including the `ID_ZDELTA` and `ID_SEQLEN` additions), the vendored
  fqzcomp quality coder, and the one-byte line-3 encoding.
- [`mechanism_insight_claim1.md`](mechanism_insight_claim1.md): variable-length support and the
  stream-by-stream comparison with PgRC2.
- [`code_mapping_claim1.md`](code_mapping_claim1.md): the stage and function in `src/encoder.cpp`
  and the shared coder headers behind each number.

## Source of truth

Where a file in this folder disagrees with `results/`, the result file is correct. This folder
explains the results and does not replace them.
