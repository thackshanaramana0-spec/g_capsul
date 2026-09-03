#!/usr/bin/env python3
"""
Tetraploid truth, v2 -- drops NOTHING.

v1 keyed on genome position and skipped any site where the two individuals'
REF strings disagreed. Audit showed that was 631 sites, 100% of them indels,
all of them STR/repeat loci where the two samples merely used different-length
representations of the same event (e.g. GAT->G vs GATATAT->G). Skipping them
systematically biased the indel truth toward easy indels -- exactly the wrong
direction for an honest benchmark.

v2: normalise each individual first (bcftools norm -f -m -any: left-align and
split multi-allelics), then take the union of ALLELE records (pos, ref, alt).
Records that are genuinely the same event unify after normalisation; records
that remain distinct are kept as separate VCF rows at the same position, which
is legal VCF and which rtg vcfeval handles correctly because it compares
haplotypes, not literal rows.

Genotype: emitted as 0/1 (allele present). --squash-ploidy scores allele
presence, not per-copy dosage, which is the only thing a reference-free caller
can be held to; the 4-copy dosage information is therefore not needed and is
not faked.
"""
import gzip, sys

def load(path):
    """(pos, ref, alt) -> True for every non-ref allele carried by this sample"""
    out = set()
    for line in gzip.open(path, 'rt'):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if f[0] != '20':
            continue
        pos, ref, alt = int(f[1]), f[3], f[4]
        if alt in ('.', ''):
            continue
        gt = f[9].split(':')[0].replace('|', '/')
        # post-normalisation each record is biallelic, so any non-zero call
        # means this sample carries this ALT
        if any(g not in ('0', '.') for g in gt.split('/')):
            out.add((pos, ref, alt))
    return out

a = load(sys.argv[1])
b = load(sys.argv[2])
union = sorted(a | b)

with open(sys.argv[3], 'w') as o:
    o.write('##fileformat=VCFv4.2\n')
    o.write('##contig=<ID=20,length=63025520>\n')
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
    for pos, ref, alt in union:
        o.write(f"20\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/1\n")

print(f"v2 truth alleles: {len(union)}  (HG003 {len(a)}, HG004 {len(b)}, "
      f"shared {len(a & b)}, dropped 0)", file=sys.stderr)
