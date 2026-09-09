# Claim 3's mechanism: a heterozygous locus is not one place

**Discovered and verified 2026-09-09 on HG002 chr20 against GIAB v4.2.1 truth.**
This supersedes the earlier framing of Claim 3 as speed or as "nobody else
offers the interface".

## The finding

In a pseudogenome built to minimise bits, **a heterozygous locus is not one
place — it is N parallel places.** The compressor gives each haplotype its own
consensus contig, because two internally-consistent contigs compress better
than one contig plus a column of disagreements. Each is then stored again in
reverse complement. Measured across 250 GIAB het SNVs on chr20:

    both alleles present in the pseudogenome : 200 / 250   (80.0%)
    only one allele found                    :  42 / 250   (16.8%)
    probe not located                        :   8 / 250
    parallel loci per site: median 4   mean 15.58   max 1007

Worked example, GIAB `20:3001343 C>T` (het). A 40 bp probe ending 6 bp before
the variant occurs four times in the pseudogenome, and the consensus base at
the variant offset differs between them:

    pg 101,247,805 (+)   C    <- REF haplotype
    pg 119,833,369 (+)   T    <- ALT haplotype
    pg   4,237,310 (-)   G    -> complement C   REF
    pg 136,671,198 (-)   A    -> complement T   ALT

The two haplotypes sit **18.6 Mb apart** in pseudogenome space.

## Why this makes coordinate addressing structurally insufficient

A coordinate names one of those places. It therefore returns one haplotype,
and the reads it returns agree with each other perfectly — a pileup with the
variation removed. Measured, same 12 sites, coordinate query at +/-100 bp
around one locus:

    sites returning BOTH alleles :  0 / 12
    sites returning ONE allele   : 12 / 12

This is not a bug in our query and it is not an implementation limit. It is
forced by the compression objective, and it is the **same mechanism Claim 2
identified**: the representation that compresses a het site best is the one
that separates its alleles. Claim 2 meets it at call time (collapse the
duplicate contigs); Claim 3 meets it at query time.

**A corollary worth stating**: the mismatch stream is therefore almost pure
sequencing noise. Real variants become separate contigs and produce NO
mismatch, because reads match their own haplotype's consensus exactly.
Measured on HG002: 7,466,871 deviations over 6,430,543 positions, only 15,515
(position,base) pairs recurring 5x or more — **98.6% of the stream is error,
1.4% is signal**. Any design that looks for variants in the mismatch stream is
looking in the wrong place.

## The engineered resolution: content addressing

Asking by SEQUENCE resolves every parallel representative of a locus at once —
both haplotypes, both strands — because they are parallel precisely in the
sense of sharing content. The query then returns their union: a real pileup.

**FULL VALIDATION — all four GIAB individuals, 100 het SNVs each, both query
modes on the identical archive and identical sites (T3.4):**

    individual   sites   coordinate both   content both   allele balance
    HG002          100                 0             81         1.06
    HG003          100                 0             88         2.43
    HG004          100                 0             78         0.92
    HG005          100                 0             92         0.84
    ALL            400                 0            339         0.97

    coordinate addressing : 0/400   (0.0%)
    content addressing    : 339/400 (84.8%)
    overall allele balance: 25,161 REF / 24,489 ALT  (ratio 0.97)

**0 of 400 against 339 of 400**, replicated across four unrelated individuals
(three Ashkenazi trio members and one Han Chinese), and the pooled allele
balance is 0.97 -- the returned evidence is unbiased. The archives were rebuilt
from scratch for this run and reproduce the sweep's sizes byte for byte
(HG003 567,123,999; HG004 606,192,930; HG005 997,518,994), so the result does
not depend on a particular build.

Per-individual detail for HG002, from the earlier n=100 run:

                              BOTH alleles   one allele   neither
    COORDINATE query                     0          100         0
    SEQUENCE  query (ours)              81           19         0

    aggregate allele balance: 3,072 REF / 3,263 ALT  (ratio 1.06)

**0/100 against 81/100.** The 1.06 ratio matters as much as the headline: the
returned evidence is unbiased, which is what makes it a pileup rather than a
read dump. An earlier n=20 run gave 16/18 by sequence and 0/12 by coordinate;
the effect does not weaken with scale.

Earlier 20-site detail, kept for the per-site view:

    sites where reads were returned      : 18 / 20
    sites returning BOTH alleles         : 16 / 18   (88.9%)

    20:3001343  C>T   33 reads   14 REF   6 ALT
    20:3008909  C>T   39 reads   14 REF   9 ALT
    20:3013651  A>G   20 reads    2 REF  13 ALT
    20:3034916  T>C   84 reads   13 REF  45 ALT
    ...

Against the coordinate control's 0/12, this is the result: **content addressing
recovers the pileup that coordinate addressing structurally cannot.**

## PRIOR ART — the phenomenon is known, in a different field

Checked 2026-09-09, after the measurements above. **Heterozygous regions
assembling into separate contigs is well documented in the DE NOVO ASSEMBLY
literature** and must not be presented as a new observation:

- Purge Haplotigs (BMC Bioinformatics 2018) — identifies "suspected duplicated
  sequences assumed to be allelic haplotypes".
- Redundans (NAR 2016) — an assembly pipeline for highly heterozygous genomes.
- Regional sequence expansion or collapse in heterozygous genome assemblies
  (PLOS Comp Biol 2020) — heterozygous regions "duplicate existing sequences
  and lead to higher rates of fragmentation as they are resolved into separate
  contigs".

The correct statement of what is ours is therefore NARROWER than "we discovered
haplotype fragmentation". It is:

1. **The setting.** In assembly, haplotigs are a QUALITY defect: the assembly
   is redundant and fragmented. In a compression archive the assembly is a
   byproduct and the archive is the product, so the same phenomenon is an
   ADDRESSABILITY defect: the archive's own coordinate system cannot express
   "this locus". That consequence is not described anywhere we could find.
2. **The measurement in that setting** (250 GIAB het sites, median 4 parallel
   loci, 0/12 coordinate queries recovering both alleles).
3. **The resolution, which is the OPPOSITE of the field's.** The assembly
   field's remedy is to PURGE the redundant haplotigs and keep a pseudo-haploid
   reference. For an archive that is exactly wrong: deleting a haplotig deletes
   an allele, and with it the lossless property and the variant. We keep every
   representative and reach them all by content instead.

Point 3 is the load-bearing one. The known remedy destroys what an archive
needs, so importing it would be a correctness bug, not a fix.

## What is claimed, and what is not

**Claimed.** A compression-optimal read archive fragments each heterozygous
locus across multiple contigs, so its own coordinate system cannot express
"this locus"; content addressing is required, and with it the archive serves a
complete reference-free pileup. Discovered by measurement, mechanism explained,
resolution implemented and validated against an external truth set.

**Not claimed.**
- Not speed. `genocat --head=100` extracts faster than we do (0.19 s vs 0.46 s
  without the index). Speed is in `docs/QUERY_WINDOWED_PLAN.md` and is a
  separate, lesser result.
- Not "first to search a compressed archive". BEETL-fastq did sequence search
  in 2014. The difference is what comes back: a BWT index returns reads
  CONTAINING the query, we return reads COVERING the locus (median 38% of ours
  do not contain the probe), and critically we return them from ALL parallel
  representatives of the locus.
- Not that competitors could not implement this. The point is that the
  fragmentation is a consequence of compressing well, so any tool in this
  family has it, and none has recognised or addressed it. PgRC2 stores per-read
  positions and exposes no read-out at all (`PgRC -h`: compress, `-d`, nothing
  else).

## Honest limits

- 2 of 20 sites returned no reads: the probe was not located in the
  pseudogenome. 8 of 250 in the wider scan. Recall of the probe search is not
  100% and has not been characterised.
- Mean parallel loci is 15.58 against a median of 4, with a maximum of 1007:
  repeat-rich loci return large read sets. A repeat-aware cap is not
  implemented.
- Allele balance is not uniform (e.g. 20:3013473 gave 20 REF / 1 ALT). The
  query returns evidence; it does not genotype. Genotyping is Claim 2.
- Replicated on 4 individuals x 100 sites, but all within the same chr20
  600 kb window. Not replicated on other chromosomes or regions.
- 61 of 400 sites returned only one allele (8-22% per individual). Cause not
  characterised: the probe may land in a repeat, or one haplotype may be
  absent from the assembly at that locus.
- Per-individual allele balance varies (0.84 to 2.43) even though the pooled
  figure is 0.97. The query returns evidence, not genotypes.

## The claim, finalised

> **No existing tool can retrieve the reads at a locus from a reference-free
> read archive**, and the reason is structural rather than a missing feature.
>
> SPRING addresses by read index and Genozip by first-N; neither is a locus,
> and `genocat --regions` is refused outright on FASTQ because the format
> carries no coordinates to index. BEETL-fastq searches by sequence but returns
> reads CONTAINING the query, not reads COVERING the locus (median 38% of ours
> do not contain the probe). CRAM answers a locus but only after aligning to an
> external reference. PgRC2 and NanoSpring build the coordinates and expose no
> read-out at all.
>
> **And the obvious fix does not work.** Adding a coordinate API to an
> assembly-based compressor returns ONE haplotype — measured, 0 of 100 GIAB het
> sites recovered both alleles — because compressing a het locus optimally
> splits its alleles onto separate contigs. The archive's own coordinate system
> therefore cannot name a locus.
>
> Content addressing resolves every parallel representative at once and returns
> a complete, balanced pileup: 81 of 100 sites, allele ratio 1.06. This is the
> capability, and the mechanism is why nobody has it.

## Reproduce

    capsule_decode export <archive> pg.fa                       # the pseudogenome
    capsule_decode query  <archive> out.fa <40bp probe>         # content address
    capsule_decode query  <archive> out.fa <START-END>          # coordinate control
