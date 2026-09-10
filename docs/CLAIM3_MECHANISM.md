# Claim 3's mechanism: a heterozygous locus is not one place

> **Status:** frozen 2026-09-10 at tag `v1.0-capsule` — code and results final.
> Authoritative numbers live in `../benchmark/results/`; verification status in
> [`../AUDIT.md`](../AUDIT.md). Where this file and a result file disagree, the result file wins.


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

**CORRECTED 2026-09-10 — read this before quoting any earlier figure.**

An earlier version of this table reported coordinate addressing at 0/400. That
was an ARTIFACT OF OUR OWN QUERY, not a property of coordinate addressing.
`query` emitted the CONSENSUS at each read's position rather than the read, so
every returned read agreed with every other and no allele difference could ever
appear — which made coordinate addressing look absolutely incapable and
flattered the comparison.

`query` now applies each read's own deviations and returns the true reads
(measured on the same locus: 121 overlapping read pairs went from 0 disagreeing
bases to 142). Re-measured on that basis:

    individual   sites   coordinate both   content both   allele balance
    HG002          100                18             85         1.09
    HG003          100                16             88         2.37
    HG004          100                23             79         0.91
    HG005          100                24             93         0.86
    ALL            400          81 (20.2%)    345 (86.2%)        1.00

**81/400 against 345/400 — a 4.3x gap, not an infinite one.** The pooled allele
balance is 1.00.

**The mechanism is unchanged and is what the 20.2% measures.** A het variant is
stored one of two ways: split across separate haplotype contigs, or kept as
per-read deviations on a single contig. Coordinate addressing can only ever
reach the second kind, and that is about a fifth of them. Content addressing
reaches both, because it resolves every parallel representative of the locus.
The claim is therefore that coordinates reach a MINORITY of het loci by
construction — which is still decisive, and is now defensible against a
reviewer who runs our own tool and gets 20% rather than 0%.

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

## The three limits, now closed by measurement

**1. Was it one window?** No. Four further windows spanning chr20, 60 sites
each, HG002:

    window (chr20)   content both   coordinate both
    10.0 - 10.6 Mb        55 / 60             0 / 60
    25.0 - 25.6 Mb        50 / 60             0 / 60
    40.0 - 40.6 Mb        56 / 60             0 / 60
    55.0 - 55.6 Mb        49 / 60             0 / 60
    ALL                  210 / 240 (87.5%)    0 / 240 (0.0%)

Consistent with the original 3.0-3.6 Mb window. Combined with the
four-individual run this is **640 het sites across 4 individuals and 5 windows,
content 549 (85.8%)**.

> **The coordinate column of the window sweep (0/240) is NOT the final figure.**
> Those four windows were measured under the earlier consensus-emitting query,
> where no allele difference could appear by construction, so coordinate
> addressing was guaranteed a 0 and the number measures our bug rather than the
> mechanism. They were not re-run. **The citable coordinate figure is the
> four-individual table: 81/400 against content 345/400** — a 4.3x gap, not an
> infinite one. The window sweep is retained only as evidence that the CONTENT
> result does not depend on one window, which is what it was run to show.

**2. Why did some sites return one allele?** Not because a haplotype is
missing. Instrumented: of all 30 one-allele sites across the window sweep,
**0 resolved a single locus and 30 resolved MULTIPLE loci that agreed**. The
parallel contigs exist; they simply carry the same base at that position.

That is the expected complement of the mechanism. A het variant is stored one
of two ways: split into separate contigs (recovered by content addressing) or
kept as a per-read MISMATCH on a single contig.

> **UPDATE 2026-09-10 — the paragraph that follows has been partly overtaken by
> a fix.** It argued that mismatch-encoded variants "cannot be made visible
> cheaply". That was true of the query as it then stood, and it is why the
> coordinate arm read 0/400. `query` now applies each read's deviations at emit
> and returns TRUE READS rather than the consensus (commit `d58fc23`), so those
> variants ARE visible, the coordinate arm is 81/400, and the content residue
> fell to 13.8%. The cost argument below still explains why the residue is
> bounded rather than zero — the ordering constraint on `mm_sym` is real — but
> "invisible" is now too strong.

Query emitted consensus, so
mismatch-encoded variants were invisible to it — and they cannot be recovered
at arbitrary offsets cheaply, because `mm_sym` is coded with an adaptive model
in READ order, so
read k's deviations require decoding all k-1 before it. **The 12-15% residue is
therefore explained and bounded, not unknown**: it is exactly the variants the
compressor chose to encode as deviations rather than as structure.

**3. Why does allele balance vary per individual (0.84-2.43)?** Because the
query returns EVIDENCE, not genotypes: it reports every read at every parallel
representative, and those representatives can have unequal read depth. Pooled
over 400 sites the ratio is 0.97, which is the meaningful figure. Genotyping is
Claim 2's job, and Claim 2 does it from the same archive.

## Remaining limits

- 2 of 20 sites returned no reads: the probe was not located in the
  pseudogenome. 8 of 250 in the wider scan. Recall of the probe search is not
  100% and has not been characterised.
- Mean parallel loci is 15.58 against a median of 4, with a maximum of 1007:
  repeat-rich loci return large read sets. A repeat-aware cap is not
  implemented.
- Allele balance is not uniform (e.g. 20:3013473 gave 20 REF / 1 ALT). The
  query returns evidence; it does not genotype. Genotyping is Claim 2.
- 4 individuals x 100 sites plus 4 further windows on HG002 (640 sites, 5
  windows) — but all on chr20, because the archives are chr20 at 30x. Not
  replicated on another chromosome or another coverage depth.
- The 12-15% residue is explained (mismatch-encoded variants, see above) but
  NOT recovered. Recovering it needs the mismatch stream to be positionally
  addressable, which its adaptive read-order coding prevents.

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
