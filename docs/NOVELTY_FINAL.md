# What is novel, finalised — with the prior art actually read

Written 2026-09-09 after a literature pass, not from intuition. Each claim
below states what is ours, what is not, and the measurement that supports it.
Where prior art overlaps, it is named.

---

## The field, as it actually stands

**1. Compressive genomics — the paradigm exists.**
Loh, Baym & Berger, *Nature Biotechnology* 2012, established computing directly
on compressed data (CaBLAST/CaBLAT). This is our paradigm and must be cited as
such. **But it was instantiated for genome and protein DATABASES to accelerate
SEARCH**, exploiting redundancy across many similar genomes. It was not
instantiated for read archives, and not for variant evidence.

**2. Assembly-based read compressors build an assembly and throw it away.**
This is not our characterisation — it is theirs. NanoSpring (*Scientific
Reports* 2023), on its own design: the contigs and consensus sequences are
"internal to the compression process and are discarded after encoding", and
"the assembly is strictly a compression intermediate ... not preserved or made
available for downstream genomic analysis". Its stated future work is
compression ratio and adding quality/identifier streams — exposing the assembly
is not contemplated. PgRC2's binary confirms the same by its interface: `PgRC
-h` offers compress, `-d` decompress, and nothing else.

**3. Searchable archives exist, and return a different object.**
BEETL-fastq (2014) and the population BWT (*Genome Research* 2017) offer
"a highly compressed, searchable, scalable archival format". They are BWT text
indexes: they return reads **containing** a query string. Measured against
ours, that is a different object — median 38% of the reads we return at a locus
do not contain the probe (n=50). They also compress far worse, because an
FM-index is not a compressor.

**4. Haplotype fragmentation is known — in assembly, with the opposite remedy.**
Purge Haplotigs (2018), Redundans (2016), and PLOS Comp Biol 2020 all describe
heterozygous regions assembling into separate contigs. **Their remedy is to
purge the redundant haplotigs** and keep a pseudo-haploid reference. For a
lossless archive that is inadmissible: deleting a haplotig deletes an allele,
a variant, and losslessness.

---

## Claim 1 — COMPACT

**Result.** 19/19 wins vs SPRING, 19/19 vs Genozip, 57/57 rows LOSSLESS.
Aggregate −6.03% vs SPRING, −43.26% vs Genozip. +1.88% vs PgRC2.

**Novelty:** incremental and honestly so. This is a strong engineering result in
a mature field, not a new idea. It earns the right to make Claims 2 and 3 —
without it, "our archive is also analysable" would be a consolation prize.

---

## Claim 2 — FAITHFUL. **The mechanism, and the strongest result.**

**Finding.** The representation that compresses a heterozygous site best is the
one that conceals it. Optimal compression puts the two alleles on separate,
internally-consistent contigs; the ref-allele and alt-allele reads are then
never placed at the same coordinate, so the variant is not present in the data
structure at all. No read-out layer can recover it.

**Evidence — full chr20, from the archive, ablation:**

    neither             F1 0.431   P 0.967   R 0.278
    re-placement only      0.426     0.962     0.274   <- WORSE than neither
    collapse only          0.648     0.963     0.488
    both                   0.888     0.956     0.830

Synergy, not additivity: +0.217 and −0.005 alone, +0.457 together. And the
error SHAPE confirms the mechanism rather than merely fitting it — without
collapse, precision holds at 0.96 while recall collapses to 0.27. The caller is
not mistaken, it is **blind**, which is what "the alt reads are on another
contig" predicts. Noise or a bad threshold would cost precision instead.

**And it costs nothing in ratio** (+1.88% vs PgRC2 stands).

**Novelty:** we found no prior statement of this tension in the compression
literature. It is a mechanism, it is measured, and the resolution is free. This
is the paper's centre.

---

## Claim 3 — ADDRESSABLE. **The same mechanism, one layer out.**

**Finding.** A het locus is not one place in a compression-optimal
pseudogenome; it is N parallel places (median 4: two haplotypes × two strands),
up to 18.6 Mb apart. **The archive's own coordinate system therefore cannot
name a locus.**

**Evidence — 100 GIAB het SNVs, both modes, same archive, same sites (T3.4):**

                              both alleles   one allele
    coordinate addressing                0          100
    content addressing (ours)           81           19
    allele balance: 3,072 REF / 3,263 ALT   (ratio 1.06)

**Novelty, stated precisely.** Three parts, in increasing order of strength:

1. *Capability.* No tool retrieves the reads at a locus from a reference-free
   read archive. SPRING addresses by read index; Genozip's `--regions` is
   refused on FASTQ for want of coordinates; BEETL returns reads containing a
   string; CRAM needs an external reference; PgRC2/NanoSpring expose nothing.
   Being first at a capability is the weakest of the three.
2. *Instantiation of an established paradigm in a domain where it was not
   attempted.* Compressive genomics said compute on the compressed form; the
   read-compression field built assemblies and discarded them, by its own
   account. We keep the assembly and compute on it. That is a bridge between
   two literatures that had not been connected.
3. *Mechanism — the load-bearing part.* **Adding a coordinate API would not
   work.** 0/100 measures exactly that. The obvious engineering fix fails for a
   structural reason, and the reason is the same one Claim 2 identified. That
   converts "nobody has done it" into "here is why it cannot be done that way,
   and here is what does work".

**What is NOT claimed.** Not speed — `genocat --head=100` extracts in 0.19 s
against our 0.46 s (0.03 s with the optional index, but Genozip is the honest
comparator and we do not lead with this). Not first to search a compressed
archive — BEETL, 2014. Not that the haplotype phenomenon is new — it is known
in assembly. Not that competitors are incapable in principle.

---

## The unifying statement

> A compressor computes far more structure about sequencing data than it
> stores: an assembly, per-read placements, per-read deviations from consensus.
> It then keeps only the projection decompression requires, and the projections
> analysis requires are discarded — not because they are expensive, but because
> the compression objective never asks for them. Worse, the objective actively
> deforms them: compressing a heterozygous site optimally separates its alleles,
> which both hides the variant from a caller and fragments the locus beyond what
> any coordinate can name.
>
> We identify that deformation, measure it, and correct it — recovering
> variant calling (F1 0.431 → 0.888) and locus retrieval (0/100 → 81/100) from
> the same archive, at no cost in compression ratio and 0.041% of archive size.

---

## Honest weaknesses, for the discussion section

- Claim 3's validation is one individual, one chromosome, a 600 kb window,
  100 sites. HG003-HG005 would close this cheaply.

**UPDATED 2026-09-10 — the five weaknesses listed here have since been closed
or corrected. Current status:**

- **T3.4 scope: CLOSED.** Now 4 individuals x 100 sites, plus 4 further windows
  on HG002 (10/25/40/55 Mb). Still chr20 only, at 30x -- that remains the one
  real scope limit, because the archives are chr20.
- **The one-allele residue: EXPLAINED, and partly closed.** It was 19/100 and
  uncharacterised. Instrumented: every such site resolved MULTIPLE loci that
  agreed, none resolved a single locus -- so it is not a missing haplotype. It
  is the mechanism's own complement: variants stored as per-read deviations
  rather than as separate contigs. `query` now applies those deviations and
  returns true reads, and the residue is 13.8%.
- **T2.4: RECONCILED and quotable.** 21/26 (80.8%), where 26 is every SNV-only
  multi-allelic site on chr20 -- a complete census, not a sample. DiscoSNP++
  emits ZERO multi-allelic records in 3,989 output records: structural, not a
  miss. The old "5/111" scored 104 indel-bearing sites that a single-base check
  cannot evaluate; the "11/18" denominator is unreproducible and is withdrawn.
- **HG005: EXPLAINED.** Not a precision collapse of unknown origin. HG005 is
  the only variable-length dataset -- 250 bp quality-trimmed, 216 distinct
  lengths, against 148 bp fixed for HG002/3/4. Truncated to fixed 148 bp the
  same reads give F1 0.897 with FP falling 3,426 -> 1,291, making it our best
  individual. The caller is tuned for fixed-length reads. See
  `HG005_EXPLAINED.md`. The published table keeps the untruncated 0.834.
- **Probe recall and repeat loci: MEASURED.** Probe length works 20-150 bp on
  both E. coli and HG002; long probes are the most SPECIFIC (148 bp resolves to
  2 loci, 12 bp to 2,029). Repeat probes are not capped -- the worst case
  measured costs 3.9 MB and 15.7 s, and truncating would turn a correct answer
  into an arbitrary one -- but above 32 loci the query says so.

**AND ONE CORRECTION AGAINST US.** T3.4's headline was published as coordinate
0/400. That was an artifact of our own query emitting the CONSENSUS rather than
the reads, so no allele difference could appear and coordinate addressing
looked absolutely incapable. With true reads the honest figure is **81/400
against 345/400 -- a 4.3x gap, not an infinite one.**
