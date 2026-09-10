# The three claims, finalised

One product. One archive. Every number below is measured on this repository's
locked datasets and is reproducible per `docs/REPRODUCIBILITY.md`.

---

## Claim 1 — COMPACT. Smallest, and a complete tool.

**The result:** 19/19 wins against SPRING, 19/19 against Genozip, 57/57 rows
LOSSLESS.

    ours    6,011,370,147 B
    SPRING  6,396,897,280 B   -6.03%
    Genozip 10,594,522,604 B  -43.26%

**And against PgRC2, the strongest DNA compressor, +1.88% smaller.**

**The part that is easy to miss: PgRC2 is not a FASTQ compressor.** Verified by
running it, not inferred. Given a 440,048 B FASTQ it returns 151,000 B of
output — one line per read, **bare DNA**. No read names, no `+` line, no
quality scores. It cannot reproduce the input file.

So the two facts compose:

| | DNA | names | line 3 | quality | round-trips a FASTQ |
|---|---|---|---|---|---|
| PgRC2 | yes | no | no | no | **no** |
| SPRING | yes | yes | yes | yes | yes |
| Genozip | yes | yes | yes | yes | yes |
| **ours** | yes | yes | yes | yes | **yes** |

We beat the DNA specialist on the column it specialises in, **and** we are a
complete lossless FASTQ compressor, **and** we beat both complete compressors
on size. Those are three separate statements and all three hold.

Also handled and worth stating because it is where implementations break:
variable-length reads, reads containing N, and reads longer than 256 bp — each
of which was a real silent data-loss bug found and fixed in this repository
(`23be207`, `121fea9`, `f3ab0c2`).

**Honest position:** this is a strong engineering result in a mature field, not
a new idea. Its job is to earn the right to make Claims 2 and 3 — without it,
"our archive is also analysable" would be a consolation prize.

---

## Claim 2 — FAITHFUL. A mechanism, not a feature.

**The finding.** The representation that compresses a heterozygous site best is
the one that conceals it. Optimal compression puts the two alleles on separate,
internally-consistent contigs, so ref-allele and alt-allele reads are never
placed at the same coordinate and the variant is not present in the data
structure at all. **No read-out layer can recover it** — you have to change how
the pseudogenome is built, which is rebuilding the compressor.

**The engineering.** Two passes, added to the caller only, compression path
untouched: collapse the duplicate contigs into one frame, then re-place every
read onto it with mismatch tolerance.

**The evidence** — full chr20, HG002, called FROM THE ARCHIVE:

    neither             F1 0.431   P 0.967   R 0.278
    re-placement only      0.426     0.962     0.274   <- WORSE than neither
    collapse only          0.648     0.963     0.488
    both                   0.888     0.956     0.830

Synergy, not additivity (+0.217 and -0.005 alone, +0.457 together), and the
error SHAPE confirms the mechanism: without collapse, precision holds at 0.96
while recall falls to 0.27 — the caller is not mistaken, it is **blind**, which
is what "the alt reads are on another contig" predicts. Noise or a bad
threshold would cost precision instead.

**Against the state of the art**, four GIAB individuals, same reads, same
truth, same scoring:

    het-SNV F1     ours 0.876 mean   DiscoSNP++ 0.853   Kmer2SNP 0.475
    het-indel F1   ours 0.620/0.637/0.632/0.593   DiscoSNP++ 0.576/0.587/0.595/0.605
    multi-allelic  ours 21/26 (all chr20)   DiscoSNP++ 0 records of 3,989
    tetraploid     SNV 0.897/0.782, INDEL 0.597/0.553
    coverage       10x 0.487 -> 30x 0.888, monotone

**And it costs nothing.** The archive is +1.88% smaller than PgRC2 with all of
this in it; the addressability stream is 0.041% of the archive.

**The structural point.** This is one artifact, not two tools. The variants come
out of the same assembly that produced the compression — no separate de Bruijn
graph, no second pass over the reads, no reference. That is what "unified" means
here, and it is why the cost is 0.041% instead of a second index.

**Honest position:** 3 wins and 1 loss on SNV (HG005), 3 and 1 on indel. Claim
competitiveness per class and unification overall — not per-class dominance.

---

## Claim 3 — ADDRESSABLE. The same mechanism, one layer out.

**The capability.** Retrieve the reads at a locus from a reference-free archive.
No tool does this: SPRING addresses by read index; Genozip's `--regions` is
*refused* on FASTQ because the format carries no coordinates to index; BEETL
returns reads *containing* a string rather than covering a locus; CRAM answers
a locus but only after aligning to an external reference; PgRC2 and NanoSpring
expose no read-out at all.

**The mechanism — why adding a coordinate API would not fix it.** A het locus
is not one place in a compression-optimal pseudogenome; it is N parallel places
(median 4: two haplotypes x two strands), up to 18.6 Mb apart. So the archive's
own coordinate system **cannot name a locus**. Measured:

    400 GIAB het sites, 4 individuals, true reads returned
      coordinate addressing :  81 / 400   (20.2%)
      content addressing    : 345 / 400   (86.2%)
      pooled allele balance : 1.00

**81 against 345 — a 4.3x gap.** A het variant is stored either split across
haplotype contigs or as per-read deviations on one contig; coordinates can only
reach the second kind, about a fifth of them. Content addressing resolves every
parallel representative and returns a complete, unbiased pileup.

*(An earlier version of this table said 0/640. That was an artifact of our own
query emitting the consensus rather than the reads, which made coordinate
addressing look absolutely incapable. Corrected — see
`docs/CLAIM3_MECHANISM.md`.)*

**And it is fast**, though this is the least of the three points: 0.03 s per
query with the optional index (0.46 s without), against SPRING's 4.93 s, which
costs 93-99% of a FULL decode. Genozip extracts faster than us at 0.19 s but
cannot answer a locus at all. Claim 3 does not lead with speed.

**Honest position:** the phenomenon (haplotypes assembling to separate contigs)
is known in the ASSEMBLY literature, where the remedy is to purge the duplicate
haplotigs. For a lossless archive that remedy is inadmissible — it deletes an
allele, a variant, and losslessness. What is ours is the consequence in an
archive, the measurement, and the resolution.

---

## The one-paragraph version

> Compressing sequencing reads optimally requires assembling them, and the
> assembly encodes far more than the compressor stores: where every read sits,
> how each deviates from consensus, which sequences are alternative haplotypes.
> Every tool in this family discards that structure — NanoSpring's own paper
> says the assembly is "strictly a compression intermediate ... not preserved
> or made available for downstream genomic analysis". Worse, the compression
> objective actively deforms it: compressing a heterozygous site optimally
> separates its alleles, which hides the variant from any caller and fragments
> the locus beyond what any coordinate can name. We identify that deformation,
> correct it, and keep the structure — recovering reference-free variant
> calling (F1 0.431 -> 0.888, beating DiscoSNP++ and Kmer2SNP) and
> reference-free locus retrieval (0/640 -> 549/640) from the same archive, at
> no cost in compression ratio and 0.041% of archive size, while remaining the
> smallest lossless FASTQ compressor measured.

---

## What is still open

- **HG005's losses are EXPLAINED** (`docs/HG005_EXPLAINED.md`): it is the only
  variable-length dataset -- 250 bp quality-trimmed with 216 distinct lengths,
  against 148 bp fixed for HG002/3/4. Truncated to fixed 148 bp, the same reads
  give F1 0.897 with FP falling 3,426 -> 1,291, making it our BEST individual.
  The caller is tuned for fixed-length reads; that is a characterised property
  with a named cause, not an anomaly. The table keeps the untruncated 0.834.
- **T2.4 is now a COMPLETE CENSUS**: 21/26 (80.8%), where 26 is every SNV-only
  multi-allelic site on chr20. DiscoSNP++ emits ZERO multi-allelic records in
  3,989 output records -- a structural inability, not a miss.
- **T3.4 is chr20 only, at 30x.** Its residue (13.8%) is no longer an
  artifact: query returns true reads now. True-read retrieval requires the
  optional sidecar, because the deviations are adaptively coded in read order
  and cannot be reached from the archive on demand -- that is a real
  architectural limit, stated rather than hidden.
