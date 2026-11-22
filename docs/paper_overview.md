# Paper overview

> Entry point: how the paper maps onto the three claims and onto the documents in this folder.

## Names

**G_CAPSUL** (Genomic, Compact, Addressable, Pseudogenome-Structured, Unified Lossless) is the
name of the tool and the subject of the paper.

**ARCH** (Assemble, Retain, Compress, Harness) is the name of the design principle the tool
follows. During compression the reads are assembled into a pseudogenome, the resulting structure
is retained in the archive, and it is reused for downstream analysis. The allele-splitting
mechanism described in `central_insight_and_mechanism.md` is evidence for the principle and is not
the principle itself.

## The three claims

| Claim | Name | Question | Folder |
|---|---|---|---|
| 1 | COMPACT | Is the archive smaller than general-purpose FASTQ compressors, and is it lossless? | `docs/claim1/` |
| 2 | FAITHFUL | Does the retained structure still support reference-free variant calling? | `docs/claim2/` |
| 3 | ADDRESSABLE | Can the same archive serve export, coverage, range query, exact match and locus retrieval? | `docs/claim3/` |

The three claims are not independent results. Claim 1 builds the structure. Claims 2 and 3 are two
separately measured demonstrations that the structure can be reused. The full argument is in
`central_insight_and_mechanism.md`.

## Which documents are authoritative

For every number, the CSV files in `results/` are the authoritative source. Where a document in
`docs/` disagrees with a result file, the result file wins. The folders `docs/claim1/`,
`docs/claim2/` and `docs/claim3/` are the current explanation of each claim. Some documents refer
by title to the project's internal working notes. Those notes are not part of this repository and
nothing in `docs/` depends on reading them.

## File index

- `technical_architecture.md`: the container format and the assembly pipeline (chaining,
  pigeonhole mapping, the MAXMAP ramp, the second region, MEM self-match). How the archive is
  built, independent of any one claim's results.
- `central_insight_and_mechanism.md`: the mechanism shared by all three claims.
- `novelty_and_prior_art.md`: comparison with PgRC2, BEETL, CIndex, sFASTQ, DiscoSNP++, Kmer2SNP,
  SPAdes and bwa.
- `limitations_and_scope.md`: trade-offs, withdrawn numbers and scope limits in one place.
- `numbers_and_verification_index.md`: index of every headline number with its source and
  verification status.
- `configuration_reference.md`: every environment variable and positional argument.
- `claim1/`, `claim2/`, `claim3/`: per-claim detail.
