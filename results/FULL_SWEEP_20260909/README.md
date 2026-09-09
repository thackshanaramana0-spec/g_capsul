# Full benchmark sweep — 2026-09-09

`benchmark_1_run.sh` at commit `21ee619`, run from a CLEAN git worktree at HEAD
(`/tmp/bench_head`) because the main tree had uncommitted third-party work in
`stages/106_inprocess.cpp`; building there would have measured that instead of
the baseline. Binaries verified free of those symbols before launch.

    started 02:46:07, finished 07:53:21, total 05:07:14
    19 datasets, 0 failures, 57/57 Claim 1 rows LOSSLESS
    machine idle at launch (load 0.55), 67 GB disk free, one timed job at a time

## Claim 1 — COMPACT (T1.1 size, T1.2 time + RAM)

**19/19 wins vs SPRING. 19/19 wins vs Genozip.**

    aggregate   ours    6,011,370,147 B
                SPRING  6,396,897,280 B   -6.03%
                Genozip 10,594,522,604 B  -43.26%

Range vs SPRING -3.23% (P. aeruginosa) to -29.63% (A. fumigatus); vs Genozip
-5.88% to -71.89% (S. cerevisiae). The four human sets are our thinnest
margins, -4.06% to -7.76%.

## Claim 2 — FAITHFUL (all five tables, ALL FROM THE ARCHIVE)

Every Claim 2 number comes from `capsule_decode call <archive>`: no FASTQ, no
reference. T2.4 and T2.5 were converted to this path on 2026-09-09; they had
been driving the encoder over the FASTQ.

T2.1 het-SNV F1        ours / DiscoSNP++ / Kmer2SNP
    HG002  0.888 / 0.847 / 0.464
    HG003  0.891 / 0.849 / 0.460
    HG004  0.891 / 0.851 / 0.469
    HG005  0.834 / 0.863 / 0.505     <- our loss
    mean   0.876 / 0.853 / 0.475     3 wins / 1 loss

T2.2 coverage sweep (HG002)   10x 0.487 | 15x 0.704 | 20x 0.797 | 30x 0.888
    Monotone, no cliff. 30x is carried from T2.1 rather than re-run.

T2.3 het-indel F1      ours / DiscoSNP++
    HG002  0.620 / 0.576      HG003  0.637 / 0.587
    HG004  0.632 / 0.595      HG005  0.593 / 0.605   <- our loss
    3 wins / 1 loss. Kmer2SNP is NOT_APPLICABLE_SNP_ONLY -- it has no indel
    model, so the blank is a property of the tool, not a gap here.

T2.4 multi-allelic (HG002 chr20:1-6Mb, 111 truth sites)
    strict, both ALT alleles recovered:  ours 5/111   DiscoSNP++ 0/111
    proxy, any call at the position:     ours 14/111  DiscoSNP++ 15/111
    Both recorded because they disagree in DIRECTION. The strict metric is
    what the claim asserts; DiscoSNP++'s 0 is the documented structural result.

T2.5 tetraploid (HG003+HG004, ploidy 4, real reads, Cooke et al. 2022)
    SNV    ours 0.897 / DiscoSNP++ 0.782
    INDEL  ours 0.597 / DiscoSNP++ 0.553

## Claim 3 — ADDRESSABLE (T3.1 export, T3.2 coverage, T3.3 query)

    dataset        T3.1 export vs SPAdes     T3.2 coverage vs bwa+mosdepth
    ERR5181310     0.40s vs 51.73s   129.3x  0.15s vs 4.15s    27.7x
    SRR2584863     0.37s vs 205.04s  554.2x  0.33s vs 7.57s    22.9x
    SRR29296997    0.13s vs 83.95s   645.8x  0.11s vs 2.53s    23.0x
    SRR37283774    1.49s vs 380.88s  255.6x  0.56s vs 15.98s   28.5x
    DRR976266      0.70s vs 548.82s  784.0x  1.49s vs 23.98s   16.1x
    HG002          5.67s vs 2678.83s 472.5x  3.04s vs 164.90s  54.2x

T3.3 query ran on ALL 19 archives. **"No competitor exists" is FALSE and was
removed** -- this project's own docs/CLAIM3_APPLICATIONS.md had already recorded
that, and the scripts still carried the claim. Verified against the literature
2026-09-09:

  SPRING        --decompress-range start end   -> by READ INDEX (1..N in file
                order, or reordered order if -r was used). Our own baseline.
  BEETL-fastq   sFASTQ searchable archive      -> by SEQUENCE (BWT substring)
  CIndex        compressed FASTQ indexes       -> by record identifier
  genocat --regions / CRAM+samtools            -> by REFERENCE coordinate,
                which requires aligning to an external reference first

So partial random access into compressed reads is well established. The
defensible statement is about the COORDINATE SYSTEM, not the capability:
SPRING returns "reads 1000-2000", an arbitrary file slice with no biological
meaning; BEETL returns "reads containing string S"; CRAM returns "reads at
chr20:1-100000" but only after alignment to a reference. We return reads at a
coordinate range of an assembly THE COMPRESSOR ITSELF BUILT, with no reference
anywhere in the pipeline. That is the claim; "no competitor" is not.

### T3.3 head-to-head, MEASURED 2026-09-09

An earlier note here said T3.3 had no competitor number. It does. SPRING has
`--decompress-range`, so the fair question is what a small slice COSTS:

    dataset       ours query   SPRING --decompress-range 1 100   SPRING full
    ERR5181310      0.50 s              5.90 s                     5.94 s
    SRR2584863      0.46 s              4.93 s                     5.28 s

SPRING's range decode costs 93-99% of decoding the ENTIRE archive: it inflates
everything and then slices. Ours costs 8-17% of our own full decode, because
`query` decodes the assembly layer (literal + mem_triples -> pseudogenome) plus
placements and stops -- it never touches quality (the largest stream), never
touches names, and never reconstructs reads. See stages/capsule_decode.cpp: the
quality thread only starts when a quality output is requested, and the
query/coverage branch returns before names.

So there are TWO distinct claims and both survive:

  1. COST. We serve a slice without a full decode; SPRING measurably cannot.
     ~11x on these two datasets.
  2. ADDRESSING. SPRING takes read INDICES -- "reads 1000-2000", an arbitrary
     file slice. It cannot answer "what is at this locus" at any cost, because
     a FASTQ archive has no loci. CRAM/genocat --regions can, but only after
     aligning to an external reference. We take a pseudogenome coordinate, from
     an assembly the compressor built, with no reference anywhere.

What is still NOT done: the same head-to-head against `genocat --regions` and
indexed CRAM, which do have coordinates and would need the reference-build cost
counted honestly on their side.

**Read T3.1 carefully.** `output_bytes`/`rows` are in the table for this
reason: our export emits 2 records (the pseudogenome), SPAdes emits tens of
thousands of biological contigs. The ratio measures time-to-a-reference-free
coordinate-system from an archive that had to exist anyway, NOT "same output,
faster". T3.2 and T3.3 carry no such caveat -- they produce the same object as
the conventional route.

## Honest notes

- HG005 loses on BOTH T2.1 and T2.3, and it is the same cause on the SNV side:
  precision falls to 0.899 with FP 3,426 against ~1,700 elsewhere. HG005 is the
  6.34 GB Han Chinese son, deeper than the standardised 30x. Not investigated.
- T2.4's 5/111 does not match the locked doc's 11/18. Different denominator,
  and that original analysis was ad-hoc and never saved. Do not quote T2.4 in
  the paper until it is reconciled.
- Phase 1 took 2h35m against a 1h21m projection (1.9x). The projection
  under-modelled competitor decompress on multi-GB inputs; results unaffected.

### T3.3 competitor numbers — measured 2026-09-09, E. coli (SRR2584863)

Retrieving a small subset from an already-compressed archive:

    tool     command                        wall     its own full decode   addresses by
    ours     query <range or sequence>      0.46 s   5.81 s                pseudogenome coordinate
    Genozip  genocat --head=100             0.19 s   0.85 s                first-N only
    SPRING   --decompress-range 1 100       4.93 s   5.28 s                read index

Read this honestly, in three parts:

1. GENOZIP IS FASTER THAN US at raw extraction, 0.19 s vs 0.46 s. Do not lead
   with speed anywhere in Claim 3.
2. COST vs POSITION. Ours is flat: 0.46-0.50 s for ranges at offset 0, 1 Mb,
   3 Mb and 26.76 Mb of a 26.96 Mb pseudogenome. Genozip's is not -- --head=N
   is sequential from the start and climbs 0.19 s -> 0.83 s as N goes
   100 -> 2,000,000, i.e. to the full-decode cost. SPRING's range decode is
   93-99% of a full decode at any N.
3. WHAT CAN BE ASKED. Genozip's --regions (real locus retrieval) is refused on
   a FASTQ archive -- "not supported for this file because it was not indexed
   during compression" -- because FASTQ carries no coordinates. SPRING takes
   read indices, which are file positions, not genome positions. Only ours
   accepts a locus, and only ours accepts a bare DNA sequence.

NOT RUN: BEETL-fastq, the one tool that also searches a compressed archive by
sequence. The argument against it is structural rather than timed -- it is a
BWT text index and returns reads CONTAINING the query, where we return reads
COVERING the locus (median 38% of ours do not contain the probe, IQR 26-46%,
n=50). That argument follows from what the index computes and does not need a
timing run, but the head-to-head has not been done and should be expected as a
reviewer request.
