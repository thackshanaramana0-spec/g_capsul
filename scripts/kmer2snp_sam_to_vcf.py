import re,sys
# Convert Kmer2SNP SNP-kmer-pair output to a genome-coordinate VCF so it can be
# scored by the SAME rtg vcfeval pipeline as ARCS (apples-to-apples). Each Kmer2SNP
# line is two k-mers differing at one base (the SNP); we align kmer1 to the reference
# (SAM) to place the SNP, take the genome REF base, and emit the alt allele.
#   usage: kmer2snp_to_vcf.py <pairs.snp> <kmer1.sam> <ref.fa> <chrom> <out.vcf>
PAIRS,SAM,REF,CHROM,OUT=sys.argv[1:6]
comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(comp.get(c,'N') for c in reversed(s))

# reference sequence for CHROM
seq=[];cur=None
for line in open(REF):
    if line[0]=='>': cur=line[1:].split()[0]; continue
    if cur==CHROM: seq.append(line.strip())
refseq=''.join(seq)
def gref(p1): return refseq[p1-1].upper() if 1<=p1<=len(refseq) else 'N'

# read pairs, index by kmer1 (the query we align). Store the two alleles + diff offset.
pairs={}   # kmer1 -> (offset, allele1, allele2)
for line in open(PAIRS):
    p=line.split()
    if len(p)<2: continue
    k1,k2=p[0],p[1]
    if len(k1)!=len(k2): continue
    diffs=[i for i in range(len(k1)) if k1[i]!=k2[i]]
    if len(diffs)!=1: continue                      # exactly one SNP
    d=diffs[0]
    pairs[k1]=(d,k1[d],k2[d])

# map each aligned kmer1 to genome, emit VCF
rows=[]
for line in open(SAM):
    if line[0]=='@': continue
    f=line.split('\t')
    qn=f[0];flag=int(f[1]);mapq=int(f[4]);cig=f[5];pos=int(f[3]);rname=f[2]
    if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20: continue
    if rname!=CHROM: continue
    seqf=f[9]                                        # aligned query seq (may be RC of kmer1)
    rev=bool(flag&0x10)
    # recover original kmer1 to find which key this alignment belongs to
    k1 = rc(seqf) if rev else seqf
    if k1 not in pairs: continue
    d,a1,a2=pairs[k1]
    # walk CIGAR: query index d -> genome position
    qtarget=d
    refidx=pos-1; readidx=0; gpos=None
    # if reverse-aligned, the SAM query is rc(k1); the base at k1[d] sits at query index (len-1-d)
    qidx = (len(k1)-1-d) if rev else d
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            if readidx<=qidx<readidx+n: gpos=refidx+(qidx-readidx)+1; break
            refidx+=n;readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
    if gpos is None: continue
    gR=gref(gpos)
    if gR=='N': continue
    # alleles in genome frame
    g1 = comp[a1] if rev else a1
    g2 = comp[a2] if rev else a2
    if   g1==gR and g2!=gR: alt=g2
    elif g2==gR and g1!=gR: alt=g1
    elif g1!=gR and g2!=gR: alt=(g1 if g1==g2 else g1+","+g2)
    else: continue
    rows.append((gpos,gR,alt))

rows.sort()
with open(OUT,'w') as o:
    o.write("##fileformat=VCFv4.2\n##contig=<ID=%s>\n"%CHROM)
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    seen=set()
    for gpos,gR,alt in rows:
        if gpos in seen: continue
        seen.add(gpos)
        gt="1/2" if "," in alt else "0/1"
        o.write("%s\t%d\t.\t%s\t%s\t30\tPASS\t.\tGT\t%s\n"%(CHROM,gpos,gR,alt,gt))
print("Kmer2SNP: %d SNP calls -> %s"%(len(seen),OUT))
