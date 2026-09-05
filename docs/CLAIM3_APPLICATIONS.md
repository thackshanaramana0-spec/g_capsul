# Claim 3: what each table is for, and how strong it is

Claim 3 is an AMORTIZATION argument. The pseudogenome is built and stored during
compression, so three analyses that normally need separate tools are served from
the archive as stream decodes. `stages/capsule_decode.cpp:160` states the
mechanism exactly:

    export   : literal + mem_triples  -> the pseudogenome      (no assembly)
    coverage : pos_abs + read_lengths -> per-position depth    (no alignment)
    query    : pos_abs + pg           -> reads in a range      (no full decode)

Each stops as soon as its streams are decoded. None pays for full reconstruction.
Every number must be reported beside the one-time compression cost (~200 s on
HG002), or the comparison is dishonest.

## T6a -- export: reference-free assembly

**Application.** Analysis when no trustworthy reference exists. Reference-based
workflows rely on "public repositories having closely-related high quality
references, which is not guaranteed for arbitrary clinical isolates", and
surveillance "requires a priori knowledge of the pathogen of interest [and]
cannot detect the emergence of completely novel pathogens". Concretely: hospital
outbreak investigation (assemble isolates, cluster, SNP distances, transmission
chains), novel-pathogen detection, and metagenomics where most organisms have no
reference at all.

**Competitors.** SPAdes 4.0.0 AND MEGAHIT v1.2.9. MEGAHIT is mandatory -- it is
~6x faster than SPAdes, so quoting SPAdes alone inflates the margin ~6x. Both
installed. QUAST is required so the exported FASTA is shown to be the same KIND
of object (N50, genome fraction, misassemblies).

## T6b -- coverage: QC and dosage

**Application.** Deciding whether a run is usable, and detecting dosage change.
Coverage is used "to assess the quality of sequencing data (e.g. percentage of
genome with >=30x read depth) or to identify genomic regions with insufficient
reads for reliable variant calling", and is the input to copy-number variant
detection. Clinical labs certify ">=30x over 95% of target" before reporting.

**Competitors.** bwa-mem + samtools + mosdepth. bwa-mem2 is 50-100% faster and is
NOT installed -- install it or state that the margin roughly halves.

**Scope.** Per-contig in pseudogenome space: answers uniformity and relative
depth, not "coverage at chr20:31,000,000". Sufficient for the QC application, not
for locus-level CNV without a coordinate mapping we do not have.

## T6c -- query: partial decode

**Application.** Retrieving a subset from cold storage without inflating the
archive: re-analysis of one region, checking a variant, sharing a slice.

**This is a REAL capability, not a wrapper.** `query` decodes `pos_abs` + `pg`
and stops -- it never materialises the reads. That is partial random access into
a compressed archive.

**Competitors.** `genocat --regions` (Genozip 15.0.87, installed) and indexed
CRAM via `samtools view` both do regional extraction, so the earlier claim that
"no competitor exists" is FALSE and must be deleted. The defensible comparison is
head-to-head on wall time and bytes read for an equivalent slice. The real
difference is the coordinate system: ours is pseudogenome space, theirs is
reference space -- which is a property of being reference-free, not a defect, but
it does mean the user needs a pg offset rather than a chromosome locus.

## What to fix before publication

1. Add MEGAHIT to T6a; add QUAST validation.
2. Install bwa-mem2 for T6b, or disclose the halved margin.
3. Delete "no competitor exists" from T6c; benchmark against genocat --regions.
4. Report the one-time compression cost in every table.
