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

Same 20 GIAB het sites, real tool, end to end:

    sites where reads were returned      : 18 / 20
    sites returning BOTH alleles         : 16 / 18   (88.9%)

    20:3001343  C>T   33 reads   14 REF   6 ALT
    20:3008909  C>T   39 reads   14 REF   9 ALT
    20:3013651  A>G   20 reads    2 REF  13 ALT
    20:3034916  T>C   84 reads   13 REF  45 ALT
    ...

Against the coordinate control's 0/12, this is the result: **content addressing
recovers the pileup that coordinate addressing structurally cannot.**

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
- All measurements are HG002 chr20. Not replicated on other individuals or
  chromosomes.

## Reproduce

    capsule_decode export <archive> pg.fa                       # the pseudogenome
    capsule_decode query  <archive> out.fa <40bp probe>         # content address
    capsule_decode query  <archive> out.fa <START-END>          # coordinate control
