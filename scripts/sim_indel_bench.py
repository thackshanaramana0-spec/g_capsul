#!/usr/bin/env python3
# Large synthetic diploid benchmark for reference-free indel calling, scored by the SAME
# gold-standard engine (rtg vcfeval) as the GIAB het-SNV numbers. hapA is the reference;
# hapB carries known het SNVs + het indels of varied length. Emits: ref.fa (hapA),
# reads.fq (50/50 diploid, RC, small error), truth.vcf (het SNVs+indels, reference coord),
# regions.bed (whole chrom). Deterministic. Mutations are well-separated and placed in
# non-repetitive local context so truth/call representations are canonical for vcfeval.
#   usage: sim_indel_bench.py <outdir> [genome_len] [cov] [seed]
import random, sys, os
out = sys.argv[1]
L   = int(sys.argv[2]) if len(sys.argv) > 2 else 60000
cov = int(sys.argv[3]) if len(sys.argv) > 3 else 50
seed= int(sys.argv[4]) if len(sys.argv) > 4 else 12345
random.seed(seed)
os.makedirs(out, exist_ok=True)
B = "ACGT"
comp = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(comp[c] for c in reversed(s))
rl, err = 150, 0.002
CHROM = "chrSIM"

hapA = [random.choice(B) for _ in range(L)]

# choose well-separated mutation sites (>=300 bp apart, away from the ends)
sites = list(range(500, L - 500, 350))
random.shuffle(sites)
n_snv, n_del, n_ins = 45, 25, 25
chosen = sorted(sites[:n_snv + n_del + n_ins])
random.shuffle(chosen)
snv_sites = sorted(chosen[:n_snv])
del_sites = sorted(chosen[n_snv:n_snv + n_del])
ins_sites = sorted(chosen[n_snv + n_del:])

def uniq_seq(k, left, right):
    # a k-mer differing from the flanking bases (avoids left-align drift); for k>=2 also
    # avoid a homopolymer run so truth/call representations stay canonical.
    for _ in range(1000):
        s = ''.join(random.choice(B) for _ in range(k))
        if k >= 2 and len(set(s)) == 1: continue
        if s[0] == right or s[-1] == left: continue
        return s
    return s

muts = {}   # pos -> ('snv',alt) | ('del',dlen) | ('ins',seq)
for p in snv_sites:
    muts[p] = ('snv', random.choice([b for b in B if b != hapA[p]]))
for p in del_sites:
    dlen = random.choice([1, 1, 2, 3, 3, 4, 6, 8, 10, 12])
    muts[p] = ('del', dlen)
for p in ins_sites:
    ilen = random.choice([1, 1, 2, 3, 3, 4, 6, 8, 10, 12])
    muts[p] = ('ins', uniq_seq(ilen, hapA[p], hapA[p + 1]))

# build hapB (reference-coordinate walk) + truth records
hapA_s = ''.join(hapA)
hapB = []
truth = []
cur = 0
for p in sorted(muts):
    kind = muts[p]
    hapB.append(hapA_s[cur:p])
    if kind[0] == 'snv':
        hapB.append(kind[1]); cur = p + 1
        truth.append((p + 1, hapA[p], kind[1], "0/1"))
    elif kind[0] == 'del':
        d = kind[1]
        # anchor = hapA[p-1]; delete hapA[p..p+d-1]
        hapB.append(hapA_s[p:p])            # nothing (deleted)
        truth.append((p, hapA_s[p-1:p+d], hapA_s[p-1], "0/1"))
        cur = p + d
    else:  # ins after p
        s = kind[1]
        hapB.append(hapA_s[p:p+1]); hapB.append(s)
        truth.append((p + 1, hapA_s[p], hapA_s[p] + s, "0/1"))
        cur = p + 1
hapB.append(hapA_s[cur:])
hapB_s = ''.join(hapB)

def emit(fh, hap, n, tag):
    for i in range(n):
        st = random.randint(0, len(hap) - rl)
        r = list(hap[st:st+rl])
        for j in range(rl):
            if random.random() < err: r[j] = random.choice([b for b in B if b != r[j]])
        r = ''.join(r)
        if random.random() < 0.5: r = rc(r)
        fh.write(f"@{tag}{i}\n{r}\n+\n{'I'*rl}\n")

nreads = (L * cov) // rl
with open(f"{out}/reads.fq", "w") as fh:
    emit(fh, hapA_s, nreads // 2, "A")
    emit(fh, hapB_s, nreads // 2, "B")

with open(f"{out}/ref.fa", "w") as fa:
    fa.write(f">{CHROM}\n")
    for i in range(0, len(hapA_s), 70): fa.write(hapA_s[i:i+70] + "\n")

with open(f"{out}/regions.bed", "w") as bd:
    bd.write(f"{CHROM}\t0\t{L}\n")

with open(f"{out}/truth.vcf", "w") as tv:
    tv.write("##fileformat=VCFv4.2\n")
    tv.write(f"##contig=<ID={CHROM},length={L}>\n")
    tv.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
    tv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    for pos, r, a, gt in sorted(truth):
        tv.write(f"{CHROM}\t{pos}\t.\t{r}\t{a}\t30\tPASS\t.\tGT\t{gt}\n")

n_ti = sum(1 for _,r,a,_ in truth if len(r)!=len(a))
print(f"genome={L} reads={nreads} truth: {n_snv} SNV + {len(truth)-n_snv} indel "
      f"({n_del} del + {n_ins} ins) -> {out}/")
