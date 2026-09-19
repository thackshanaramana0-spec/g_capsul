#!/usr/bin/env python3
"""
Are the ALT alleles of the nine residual sites present in the INPUT READS?

Counts alleles directly from win30x.bam -- the reads fed to the encoder, before
any compression -- via samtools mpileup. Without a reference FASTA mpileup
emits literal bases rather than '.'/',', so bases are counted by letter.

If ALT is present in the input but the query returns none, the loss is ours.
If ALT is absent from the input, no query could have returned it and the site
is a property of the 30x subsample.

A control group of sites that SUCCEEDED is measured the same way, so a broken
measurement shows up as both groups reading zero.
"""
import re, subprocess, sys, collections

BAM="win30x.bam"
FAIL=[3036066,3053969,3059523,3163883,3172647,3172774,3325265,3345618,3349190]

def giab():
    d={}
    for line in open("giab_win.vcf"):
        if line.startswith("#"): continue
        f=line.split("\t")
        d[int(f[1])]=(f[3],f[4])
    return d

def pileup(pos):
    out=subprocess.run(["samtools","mpileup","-r",f"20:{pos}-{pos}",
                        "-d","100000","-Q","0","-B",BAM],
                       capture_output=True,text=True).stdout.strip()
    if not out: return None
    f=out.split("\t")
    if len(f)<5: return None
    bases=f[4]
    # strip indel runs and read-start markers so only aligned bases remain
    bases=re.sub(r'\^.','',bases).replace('$','')
    bases=re.sub(r'[+-](\d+)', lambda m:' '*0, bases)
    cnt=collections.Counter(b.upper() for b in bases if b.upper() in "ACGT")
    return int(f[3]), cnt

def show(title, positions, V):
    print(f"=== {title} ===")
    print(f"{'pos':>10}{'ref':>5}{'alt':>5}{'REF':>7}{'ALT':>7}{'depth':>7}  verdict")
    for pos in positions:
        if pos not in V:
            print(f"{pos:>10}  not in GIAB window"); continue
        r,a=V[pos]
        got=pileup(pos)
        if got is None:
            print(f"{pos:>10}{r:>5}{a:>5}{'-':>7}{'-':>7}{'-':>7}  no pileup"); continue
        dp,cnt=got
        nr,na=cnt.get(r,0),cnt.get(a,0)
        if na==0:   v="ALT ABSENT from input"
        elif na<=2: v="ALT near-absent (<=2)"
        else:       v="ALT PRESENT -- lost downstream"
        print(f"{pos:>10}{r:>5}{a:>5}{nr:>7}{na:>7}{dp:>7}  {v}")
    print()

V=giab()
show("the 9 residual sites", FAIL, V)

succ=[]
for line in open("probes_real.tsv"):
    p=int(line.split("\t")[0])
    if p not in FAIL: succ.append(p)
show("control: 9 sites that SUCCEEDED", succ[:9], V)
