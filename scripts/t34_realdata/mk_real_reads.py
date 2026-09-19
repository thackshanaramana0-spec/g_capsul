#!/usr/bin/env python3
"""
Build a REAL diploid test input for T3.4-v2 from REAL data:

  reference   GRCh37 chr20:2,999,001-3,601,000, fetched from UCSC and already
              verified base-for-base against the GIAB REF column (870/870).
  variants    the GIAB HG002 v4.2.1 benchmark VCF, streamed for that window.
              Phase is taken from GIAB's own GT field -- haplotypes are not
              invented here.

Only the READS are simulated, because HG002's actual FASTQ is ~200 GB and this
box has 7 GB of RAM. Everything the experiment is actually testing -- the
encoder's pseudogenome construction, the archive, the sidecar, the query search
and the scorer -- runs for real on the output.

Writes: hap1.txt hap2.txt reads.fq sites.tsv
"""
import random, sys

REF0 = 2999001                      # 1-based coord of ref[0]
L    = 148                          # read length
COV  = 30                           # total diploid coverage
ERR  = 0.002                        # per-base substitution rate
SEED = 20260915

COMP = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))

def main():
    rng = random.Random(SEED)
    ref = open("ref.txt").read().strip().upper()
    h1 = list(ref)
    h2 = list(ref)

    n_het = n_hom = n_skip = 0
    hets = []                       # (pos, ref, alt) real GIAB het SNVs
    for line in open("giab_win.vcf"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        pos, r, a, filt = int(f[1]), f[3], f[4], f[6]
        if filt != "PASS":
            n_skip += 1; continue
        if len(r) != 1 or len(a) != 1 or a not in "ACGT" or r not in "ACGT":
            n_skip += 1; continue    # indels/MNPs are not modelled here
        gt = f[9].split(":")[0]
        i = pos - REF0
        if i < 0 or i >= len(ref):
            n_skip += 1; continue
        if gt in ("1|1", "1/1"):
            h1[i] = a; h2[i] = a; n_hom += 1
        elif gt in ("0|1", "0/1"):
            h2[i] = a; n_het += 1; hets.append((pos, r, a))
        elif gt in ("1|0", "1/0"):
            h1[i] = a; n_het += 1; hets.append((pos, r, a))
        else:
            n_skip += 1

    hap1 = "".join(h1); hap2 = "".join(h2)
    open("hap1.txt","w").write(hap1)
    open("hap2.txt","w").write(hap2)
    with open("sites.tsv","w") as fh:
        for pos,r,a in hets:
            fh.write(f"{pos}\t{r}\t{a}\n")

    # reads: COV/2 per haplotype, uniform starts, random strand, substitution noise
    nper = (len(ref) * (COV // 2)) // L
    nreads = 0
    with open("reads.fq","w") as fq:
        for hname, hap in (("h1", hap1), ("h2", hap2)):
            for _ in range(nper):
                st = rng.randrange(0, len(hap) - L)
                s = list(hap[st:st+L])
                if 'N' in s:
                    continue
                for j in range(L):
                    if rng.random() < ERR:
                        s[j] = rng.choice([c for c in "ACGT" if c != s[j]])
                seq = "".join(s)
                if rng.random() < 0.5:
                    seq = rc(seq)
                nreads += 1
                fq.write(f"@r{nreads}_{hname}_{st}\n{seq}\n+\n{'I'*L}\n")

    print(f"window            : chr20:{REF0}-{REF0+len(ref)-1}  ({len(ref)} bp)")
    print(f"GIAB het SNVs     : {n_het}")
    print(f"GIAB hom-alt SNVs : {n_hom}")
    print(f"skipped (indel/etc): {n_skip}")
    print(f"haplotype diffs   : {sum(1 for a,b in zip(hap1,hap2) if a!=b)}")
    print(f"reads written     : {nreads}  ({L} bp, ~{nreads*L/len(ref):.1f}x)")

if __name__ == "__main__":
    main()
