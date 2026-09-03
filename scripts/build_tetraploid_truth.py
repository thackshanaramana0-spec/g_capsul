#!/usr/bin/env python3
"""
Build a TETRAPLOID truth VCF from two real diploid GIAB samples, following
Cooke, Wedge & Lunter, "Benchmarking small-variant genotyping in polyploids",
Genome Research 2022 (PMC8805713): the polyploid sample is the union of the
parents' real variants, and its reads are the parents' real reads concatenated.

Nothing here is simulated: every allele comes from a real GIAB v4.2.1 call on a
real individual. The only artificial step is treating two real diploids as one
4-copy organism, which is exactly what the cited paper does.

Genotype construction: HG003 contributes 2 chromosome copies, HG004 contributes
2. At each union site, the tetraploid genotype is the multiset of the two
samples' alleles (HG003's two + HG004's two), with 0 = reference for a sample
that has no call there.

Output is written with --squash-ploidy-style scoring in mind: we emit the set of
non-reference ALT alleles present, which is what a reference-free caller can
actually be held to (it recovers alleles, not per-copy dosage).
"""
import gzip, sys

def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def load(path):
    """pos -> (ref, [alt alleles present in this sample's genotype])"""
    out = {}
    for line in opener(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if f[0] != '20':
            continue
        pos, ref, alt = int(f[1]), f[3], f[4]
        if alt == '.' or not alt:
            continue
        gt = f[9].split(':')[0].replace('|', '/')
        alts = alt.split(',')
        present = []
        for g in gt.split('/'):
            if g in ('.', '0'):
                continue
            try:
                i = int(g)
            except ValueError:
                continue
            if 1 <= i <= len(alts):
                present.append(alts[i - 1])
        if present:
            out[pos] = (ref, present)
    return out

a = load(sys.argv[1])   # HG003
b = load(sys.argv[2])   # HG004
outp = sys.argv[3]

positions = sorted(set(a) | set(b))
written = skipped = 0
with open(outp, 'w') as o:
    o.write('##fileformat=VCFv4.2\n')
    o.write('##contig=<ID=20,length=63025520>\n')
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
    for p in positions:
        ra, aa = a.get(p, (None, []))
        rb, ab = b.get(p, (None, []))
        ref = ra if ra is not None else rb
        # if the two samples disagree on the REF string at this position the
        # site is not cleanly mergeable (different indel representations);
        # skip rather than guess -- this is a truth set, it must not be wrong.
        if ra is not None and rb is not None and ra != rb:
            skipped += 1
            continue
        alts = []
        for al in aa + ab:
            if al != ref and al not in alts:
                alts.append(al)
        if not alts:
            continue
        # tetraploid genotype: 4 copies, HG003's 2 then HG004's 2.
        # a copy is 0 when that sample carries reference there.
        def code(sample_alts):
            g = []
            for al in sample_alts:
                g.append(str(alts.index(al) + 1) if al in alts else '0')
            while len(g) < 2:
                g.append('0')
            return g[:2]
        gt = '/'.join(code(aa) + code(ab))
        o.write(f"20\t{p}\t.\t{ref}\t{','.join(alts)}\t50\tPASS\t.\tGT\t{gt}\n")
        written += 1

print(f"tetraploid truth sites written: {written}, skipped (REF mismatch): {skipped}",
      file=sys.stderr)
