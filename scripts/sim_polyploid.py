#!/usr/bin/env python3
# Synthetic POLYPLOID genome for validating multi-allelic reference-free SNV calling.
# k haplotypes share a common backbone (hap0 = reference); at planted sites 2..k of the
# haplotypes carry different bases, giving biallelic AND multiallelic (triallelic for k=3)
# heterozygous SNVs. Reads are drawn evenly from all k haplotypes. Emits ref.fa (hap0),
# reads.fq, truth.vcf (het multiallelic), regions.bed. Deterministic.
#   usage: sim_polyploid.py <outdir> [k] [genome_len] [cov] [seed]
import random, sys, os
out = sys.argv[1]
K   = int(sys.argv[2]) if len(sys.argv) > 2 else 3
L   = int(sys.argv[3]) if len(sys.argv) > 3 else 40000
cov = int(sys.argv[4]) if len(sys.argv) > 4 else 60
seed= int(sys.argv[5]) if len(sys.argv) > 5 else 42
random.seed(seed); os.makedirs(out, exist_ok=True)
B="ACGT"; comp={'A':'T','C':'G','G':'C','T':'A'}; rl,err=150,0.002; CHROM="chrPLD"
def rc(s): return ''.join(comp[c] for c in reversed(s))

hap0=[random.choice(B) for _ in range(L)]
haps=[list(hap0) for _ in range(K)]

sites=list(range(400,L-400,300)); random.shuffle(sites)
n_multi, n_bi = 40, 40                                  # triallelic + biallelic het sites
chosen=sorted(sites[:n_multi+n_bi]); random.shuffle(chosen)
multi=sorted(chosen[:n_multi]); bi=sorted(chosen[n_multi:])
truth=[]
for p in multi:                                        # all K haplotypes differ (up to 4)
    alleles=random.sample(B,K)                          # K distinct bases
    # BUG FIX 2026-09-02: the old guard only fired when hap0's base was ABSENT
    # from `alleles`. When it was present but NOT at index 0, haps[0][p] was set
    # to alleles[0] -- a different base -- so ref.fa (written from haps[0])
    # diverged from truth's REF (taken from the untouched hap0 list). Measured:
    # 22 of 80 truth records (28%), ALL of them multi-allelic sites, had a REF
    # that disagreed with the reference genome they were generated against, so
    # a caller reporting the CORRECT genome base could never match them.
    # Guarantee alleles[0] == hap0[p] in both cases.
    if hap0[p] in alleles:
        _i = alleles.index(hap0[p]); alleles[0], alleles[_i] = alleles[_i], alleles[0]
    else:
        alleles[0] = hap0[p]
    for h in range(K): haps[h][p]=alleles[h]
    alts=sorted(set(alleles)-{hap0[p]})
    if alts: truth.append((p+1,hap0[p],",".join(alts)))
for p in bi:                                            # two allele classes among K haps
    alt=random.choice([b for b in B if b!=hap0[p]])
    for h in range(K):
        if h%2==1: haps[h][p]=alt
    truth.append((p+1,hap0[p],alt))

haps_s=[''.join(h) for h in haps]
def emit(fh,hap,n,tag):
    for i in range(n):
        s=random.randint(0,len(hap)-rl); r=list(hap[s:s+rl])
        for j in range(rl):
            if random.random()<err: r[j]=random.choice([b for b in B if b!=r[j]])
        r=''.join(r)
        if random.random()<0.5: r=rc(r)
        fh.write(f"@{tag}_{i}\n{r}\n+\n{'I'*rl}\n")
per=(L*cov)//rl//K
with open(f"{out}/reads.fq","w") as fh:
    for h in range(K): emit(fh,haps_s[h],per,f"h{h}")
with open(f"{out}/ref.fa","w") as fa:
    fa.write(f">{CHROM}\n")
    for i in range(0,len(haps_s[0]),70): fa.write(haps_s[0][i:i+70]+"\n")
with open(f"{out}/regions.bed","w") as bd: bd.write(f"{CHROM}\t0\t{L}\n")
with open(f"{out}/truth.vcf","w") as tv:
    tv.write("##fileformat=VCFv4.2\n"+f"##contig=<ID={CHROM},length={L}>\n")
    tv.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
    tv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    for pos,r,a in sorted(truth):
        na=a.count(",")+1
        gt="/".join(str(x) for x in range(1,na+1)) if na>1 else "0/1"   # multiallelic GT
        tv.write(f"{CHROM}\t{pos}\t.\t{r}\t{a}\t30\tPASS\t.\tGT\t{gt}\n")
print(f"K={K} genome={L} truth: {len([t for t in truth if ',' in t[2]])} multiallelic + "
      f"{len([t for t in truth if ',' not in t[2]])} biallelic het SNV -> {out}/")
