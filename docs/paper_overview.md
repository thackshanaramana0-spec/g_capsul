# Documentation guide

This folder explains how the G_CAPSUL results were produced and where each statement in the
manuscript can be checked. It complements the README, which covers installation and usage, and it
does not repeat the manuscript or its supplementary material.

## What G_CAPSUL is

G_CAPSUL is a lossless, reference-free compressor for short-read FASTQ files. It assembles the
reads into a pseudogenome while compressing and keeps the read placements and local sequence
differences in the archive. The same retained structure serves three purposes, which the
documentation groups as three claims.

| Claim | Question it answers | Detail |
|---|---|---|
| 1. Compression | Is the archive smaller than that of general-purpose FASTQ compressors, and is it lossless? | [`claim1/`](claim1/) |
| 2. Variant calling | Can heterozygous variants be called from the archive without the original reads or a reference? | [`claim2/`](claim2/) |
| 3. Addressability | Can the archive export sequence, report coverage and retrieve reads by locus or by content? | [`claim3/`](claim3/) |

The claims are connected. Claim 1 builds the structure, and Claims 2 and 3 are two separate
demonstrations that the structure can be reused. The argument is in
[`central_insight_and_mechanism.md`](central_insight_and_mechanism.md).

The name G_CAPSUL stands for Genomic, Compact, Addressable, Pseudogenome-Structured, Unified
Lossless. ARCH (Assemble, Retain, Compress, Harness) names the design principle behind it.

## Which document to read

| If you want | Read |
|---|---|
| to see how each claim connects to its figure, code, results and documents | [`REPO_MAP.md`](REPO_MAP.md) |
| to understand the idea and the mechanism behind the three capabilities | [`central_insight_and_mechanism.md`](central_insight_and_mechanism.md) |
| to see how the archive and the pseudogenome are built | [`technical_architecture.md`](technical_architecture.md) |
| to see how G_CAPSUL relates to other tools, area by area | [`novelty_and_prior_art.md`](novelty_and_prior_art.md) |
| to look up a reported number, its source file and how it was checked | [`numbers_and_verification_index.md`](numbers_and_verification_index.md) |
| to know what is supported, what is refused and what was not evaluated | [`PROJECT_SCOPE.md`](PROJECT_SCOPE.md) and [`limitations_and_scope.md`](limitations_and_scope.md) |
| to find every environment variable and positional argument | [`configuration_reference.md`](configuration_reference.md) |
| background engineering notes and the dataset and server setup | [`extras/`](extras/) |

## Authoritative sources

Numbers come from the CSV files in [`../results/`](../results/). If a document and a result file
disagree, the result file is correct. Each result file records the script that produced it, and
[`numbers_and_verification_index.md`](numbers_and_verification_index.md) lists them together with
how each number was checked.
