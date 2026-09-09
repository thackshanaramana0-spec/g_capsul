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

T3.3 query ran on ALL 19 archives; no competitor exists for coordinate-range
retrieval from a compressed archive.

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
