# All eight tables — final, with the file that holds each

Full sweep 2026-09-09 (commit `21ee619`, clean worktree, idle box, one timed
job at a time), plus T3.4 added 2026-09-10. Reproducibility across machines is
verified in `docs/REPRODUCIBILITY.md`.

| table | claim | what it measures | file |
|---|---|---|---|
| T1.1 | COMPACT | archive size, 19 datasets x 3 tools | `claim1_T1.1_T1.2.csv` |
| T1.2 | COMPACT | compress/decompress wall + peak RAM | same file |
| T2.1 | FAITHFUL | het-SNV F1, 3 arms, 4 individuals | `claim2_T2.1_snv.csv` |
| T2.2 | FAITHFUL | coverage sweep 10/15/20/30x | `claim2_T2.2_coverage_sweep.csv` |
| T2.3 | FAITHFUL | het-indel F1 | `claim2_T2.3_indel.csv` |
| T2.4 | FAITHFUL | multi-allelic sites recovered | `claim2_T2.4_multiallelic.csv` |
| T2.5 | FAITHFUL | tetraploid SNV + indel | `claim2_T2.5_tetraploid.csv` |
| T3.1 | ADDRESSABLE | export vs SPAdes | `claim3_T3.1_T3.2_T3.3.csv` |
| T3.2 | ADDRESSABLE | coverage vs bwa+samtools+mosdepth | same file |
| T3.3 | ADDRESSABLE | query cost | same file |
| T3.4 | ADDRESSABLE | **locus retrieval FIDELITY** | `claim3_T3.4_locus_fidelity.csv` |

## Headline numbers

**T1.1 — 19/19 wins vs SPRING, 19/19 vs Genozip, 57/57 rows LOSSLESS.**

    ours    6,011,370,147 B
    SPRING  6,396,897,280 B   -6.03%
    Genozip 10,594,522,604 B  -43.26%

**T1.2** — 2.99x SPRING on compress, 1.79x on decompress in aggregate. Faster
than SPRING on 3 of 19 and lighter on 9 of 19, so it is not a uniform trade.
The two >10 GB sets dominate: without them we are ~2.3x.

**T2.1 het-SNV F1** (mean 0.876 / 0.853 / 0.475) — 3 wins, 1 loss (HG005).

    HG002 0.888 | HG003 0.891 | HG004 0.891 | HG005 0.834
    DiscoSNP++    0.847 / 0.849 / 0.851 / 0.863
    Kmer2SNP      0.464 / 0.460 / 0.469 / 0.505

**T2.2** — 10x 0.487, 15x 0.704, 20x 0.797, 30x 0.888. Monotone, no cliff.

**T2.3 het-indel F1** — 0.620 / 0.637 / 0.632 / 0.593 against DiscoSNP++
0.576 / 0.587 / 0.595 / 0.605. 3 wins, 1 loss. Kmer2SNP has no indel model.

**T2.4 — RECONCILED 2026-09-10.** Strict both-allele **5/7 vs DiscoSNP++ 0/7**
on the SNV-only subset. The earlier 5/111 was scoring 104 indel-bearing sites
that a single-base comparison cannot evaluate; the earlier "11/18" is withdrawn
because 111 is the reproducible denominator (confirmed by
docs/POLYPLOID_BENCHMARK.md's own setup: 971,250 reads, 111 multi-allelic
sites, both reproduced exactly). Quotable now, with its scope: 7 sites, a
capability demonstration rather than a rate.

**T2.5 tetraploid** — SNV 0.897 vs 0.782, INDEL 0.597 vs 0.553.

**T3.1/T3.2/T3.3** — export 129-784x vs SPAdes, coverage 16-54x vs
bwa+mosdepth, query on all 19 archives. Read T3.1 with its caveat: our export
emits 2 records (the pseudogenome), SPAdes emits tens of thousands of contigs.
The ratio is time-to-a-reference-free-coordinate-system, not "same output".

**T3.4 — the mechanism table.**

    individual   sites   coordinate both   content both   balance
    HG002          100                 0             81      1.06
    HG003          100                 0             88      2.43
    HG004          100                 0             78      0.92
    HG005          100                 0             92      0.84
    ALL            400                 0            339      0.97

Plus 4 further windows on HG002 (10/25/40/55 Mb): 210/240 by content, 0/240 by
coordinate. **640 sites, 4 individuals, 5 windows: coordinate 0, content 549.**

## Every number that is NOT clean

- **HG005's losses are explained**, see `docs/HG005_EXPLAINED.md`: it is the
  only variable-length set (250 bp, 216 distinct lengths vs 148 bp fixed).
  Truncated to fixed length the same reads give F1 0.897 (FP 3,426 -> 1,291),
  our best individual. The table keeps the untruncated 0.834 deliberately.
- **T3.4's 12-15% residue** is explained (mismatch-encoded variants, invisible
  to a consensus-emitting query) but not recovered.
- **T3.4 is chr20 only, at 30x.** Not replicated on another chromosome.
- **Genozip's archive drifts ~37 B run to run** (0.005%); it embeds run
  metadata. Never gate on its byte-identity.
- **T3.3 is not a speed win**: `genocat --head=100` extracts in 0.19 s against
  our 0.46 s. With the optional sidecar index we reach 0.03 s, but Genozip is
  the honest comparator and Claim 3 does not lead with speed.
