# PAPER — the whole submission, in one place

**G_CAPSUL / CAPSULE** — a retained assembly for unified lossless FASTQ
compression, reference-free variant calling, and reference-free locus retrieval.

This folder is the paper's source material, condensed from 114 working documents
in `../docs/`. **Where any file here and a result file disagree, the result file
wins** (`../benchmark/results/`), and where a result file and the script that
regenerates it disagree, re-run the script.

| read this | for |
|---|---|
| **[METHODS.md](METHODS.md)** | what the system does, as implemented |
| **[ARCHITECTURE.md](ARCHITECTURE.md)** | the same thing in diagrams — pipeline, streams, both mechanisms |
| **[RESULTS.md](RESULTS.md)** | every table, as executed |
| **[DISCUSSION.md](DISCUSSION.md)** | the field, what is ours, what is not, and six refuted hypotheses |
| **[LIMITATIONS.md](LIMITATIONS.md)** | what this work does not establish |
| [NOVELTY_FINAL.md](NOVELTY_FINAL.md) | the novelty argument with prior art read as primary sources |
| [CLAIMS_FINAL.md](CLAIMS_FINAL.md) | the three claims, stated to be quoted |
| [CLAIM3_MECHANISM.md](CLAIM3_MECHANISM.md) | the locus-fragmentation measurement in full |
| [FORMAT.md](FORMAT.md) · [PIPELINE.md](PIPELINE.md) | archive format · command reference |
| [ARCHITECTURE_VS_DISCOSNP.md](ARCHITECTURE_VS_DISCOSNP.md) | layer-by-layer against the closest caller |
| [COMPRESSION_DERIVED_CALLING.md](COMPRESSION_DERIVED_CALLING.md) | the channels that come free, and why they still lose |

---

## Abstract, in the form the numbers support

Lossless FASTQ compressors that reach the best ratios do so by assembling the
reads — and then discard the assembly. NanoSpring states this of its own design:
the contigs are *"strictly a compression intermediate ... not preserved or made
available for downstream genomic analysis."* We show that the discarded
structure is not merely unexploited but **deformed by the compression objective
itself**: the representation that compresses a heterozygous site best is the one
that conceals it, because two internally-consistent contigs compress better than
one contig plus a column of disagreements. The ref-allele and alt-allele reads
are then never placed at the same coordinate, and the variant is absent from the
data structure — no read-out layer can recover it.

We identify this deformation, measure it, and correct it inside the caller
alone, leaving the compression path untouched. On 19 datasets (55.0 GB) the
archive is smaller than SPRING on 19/19 and Genozip on 19/19 (−6.03% / −43.26%
aggregate), with 57/57 archives verified lossless by decode-and-diff. From that
archive — no FASTQ, no reference, no second assembly pass — heterozygous SNV
calling reaches mean F1 0.876 against DiscoSNP++ 0.853 and Kmer2SNP 0.475 across
four GIAB individuals, and recovers 21 of the 26 multi-allelic sites on chr20
where DiscoSNP++ emits zero such records in 3,989. The same deformation appears
one layer out: a heterozygous locus is not one place in a compression-optimal
pseudogenome but N parallel places (median 4, up to 18.6 Mb apart), so the
archive's own coordinate system cannot name a locus — content addressing
resolves 345 of 400 sites against 81 by coordinate. The correction costs nothing
in compression ratio, and addressability costs 0.041% of the archive.

---

## The three claims, and what each is worth

### Claim 1 — COMPACT

**19/19 vs SPRING (−6.03%), 19/19 vs Genozip (−43.26%), 57/57 LOSSLESS.**
Separately, +1.88% smaller than PgRC2 on the DNA stream — which is a different
measurement, because PgRC2 stores no names, no `+` line and no quality and
cannot reproduce a FASTQ.

*Honest position:* incremental, and we say so. A strong engineering result in a
mature field, not a new idea. Its job is to earn the right to make Claims 2 and
3 — without it, "our archive is also analysable" is a consolation prize.

### Claim 2 — FAITHFUL

Called **from the archive**. het-SNV mean F1 **0.876** / 0.853 / 0.475;
het-indel **0.621** / 0.591; multi-allelic **21/26** vs **0/26**; tetraploid SNV
**0.897** vs 0.782, from real HG003+HG004 reads. Coverage 10× → 30× monotone.

The mechanism is the paper's centre, and the ablation shows it is a mechanism
rather than a fit: F1 0.431 → 0.888, with **re-placement alone worse than doing
nothing** (0.426), and precision *holding* at 0.96 while recall collapses to
0.27 without collapse — the caller is blind, not mistaken.

*Honest position:* competitiveness per class and unification overall, not
per-class dominance. 3 wins and 1 loss on SNV; 3 and 1 on indel.

### Claim 3 — ADDRESSABLE

`export` 129–784× vs SPAdes, `coverage` 16–54× vs bwa+mosdepth, and the
capability nothing else has: retrieving the reads at a locus from a
reference-free archive. **81/400 by coordinate against 345/400 by content.**

*Honest position:* not a speed claim — `genocat --head=100` is faster than our
query. Not first to search a compressed archive — BEETL, 2014. What is ours is
the mechanism: adding a coordinate API would not work, and 81/400 measures why.

---

## What a reviewer should check first

1. **`../benchmark/documentation/RESULT_CODE.md`** — every number traced to the
   script that made it and the file that holds it.
2. **The one-command reproduction:**
   `CLAIMS=1 bash benchmark/scripts/sanity_archive_one.sh SRR2584863`
   → archive **68,429,027 B exactly**, 3/3 LOSSLESS. Verified byte-identical
   across 2/3/7/12 cores. A different size on any machine is a real regression.
3. **`LIMITATIONS.md`** — including the two corrections we made against
   ourselves (T3.4's 0/400 → 81/400, and het-indel's withdrawal → 3-win).

## Provenance

One run produced every table: `benchmark_1_run.sh` at commit `21ee619`,
2026-09-09, 5 h 07 m, 19 datasets, 0 failures, executed from a clean git
worktree at HEAD. `../benchmark/results/` is a `cmp`-verified byte-identical
copy of it.
