#!/usr/bin/env python3
"""
Threshold sensitivity for the T3.4-v2 real run.

The published "both / one / neither" rule counts a site as recovering both
alleles if at least ONE read supports each. At 30x with a 0.2% substitution
rate a single miscalled base satisfies that, which is why the homozygous
negative control shows any "both" at all.

This asks the question that actually matters: does the 400/400 result survive a
stricter rule, and do the false positives die under the same rule? If the
positive set holds while the negative set collapses, the recovery is real
signal and not a loosened metric. If both move together, it is metric-gaming
and the claim fails.

Run from ~/t34real. Reads the already-written q/ (het sites) and qneg/
(homozygous sites) outputs -- no re-querying, so both sets are scored by
exactly the same code on exactly the same files as the headline run.
"""
import re, sys, os

COMP = {'A':'T','C':'G','G':'C','T':'A'}

def read_fa(path):
    out=[]; p=None
    try:
        for l in open(path):
            if l.startswith('>'):
                m=re.search(r'pos=(\d+)', l); p=int(m.group(1)) if m else None
            elif p is not None:
                s=l.strip()
                if s: out.append((p,s))
    except FileNotFoundError:
        pass
    return out

def read_occ(path):
    out=[]
    try:
        for l in open(path):
            m=re.match(r'\[query\] occ (\d+) (\d+) ([+-])', l)
            if m: out.append((int(m.group(1)), int(m.group(2)), m.group(3)))
    except FileNotFoundError:
        pass
    return out

def score_offset(reads, occs, ref, alt):
    vs=[(at+45,'+') if st=='+' else (at-6,'-') for at,_e,st in occs]
    nr=na=0
    for p,s in reads:
        for v,st in vs:
            off=v-p
            if 0<=off<len(s):
                b=s[off]
                if st=='-': b=COMP.get(b,'N')
                if b==ref: nr+=1
                elif b==alt: na+=1
                break
    return nr,na

def collect(probes, qdir):
    res=[]
    for line in open(probes):
        pos,ref,alt,probe,v = line.rstrip("\n").split("\t")
        r1=read_fa(f"{qdir}/s1_{pos}.fa"); o1=read_occ(f"{qdir}/s1_{pos}.occ")
        res.append(score_offset(r1,o1,ref,alt))
    return res

def main():
    pos_set = collect("probes.tsv",     "q")
    neg_set = collect("probes_neg.tsv", "qneg")
    print(f"het sites: {len(pos_set)}   homozygous sites: {len(neg_set)}")
    print()
    print("rule: a site recovers BOTH alleles if minor-allele support meets the threshold")
    print()
    print(f"{'threshold':<34}{'het both':>10}{'hom both':>10}{'separation':>13}")
    rules = []
    for k in (1,2,3,4,5):
        rules.append((f"minor allele >= {k} read(s)", lambda nr,na,k=k: min(nr,na)>=k))
    for f in (0.05,0.10,0.15,0.20):
        rules.append((f"minor allele fraction >= {f:.2f}",
                      lambda nr,na,f=f: (nr+na)>0 and min(nr,na)/(nr+na)>=f))
    rules.append(("min 2 reads AND fraction >= 0.10",
                  lambda nr,na: min(nr,na)>=2 and (nr+na)>0 and min(nr,na)/(nr+na)>=0.10))
    for label,rule in rules:
        hp=sum(1 for nr,na in pos_set if rule(nr,na))
        hn=sum(1 for nr,na in neg_set if rule(nr,na))
        print(f"{label:<34}{hp:>10}{hn:>10}{hp-hn:>13}")
    print()
    # distribution of minor-allele fraction, the thing that actually separates
    def frac(s):
        return sorted(min(nr,na)/max(nr+na,1) for nr,na in s)
    fp=frac(pos_set); fn=frac(neg_set)
    def pct(a,q): return a[int(q*(len(a)-1))]
    print("minor-allele fraction distribution")
    print(f"{'':<12}{'p5':>8}{'p25':>8}{'median':>8}{'p75':>8}{'p95':>8}")
    for nm,a in (("het",fp),("homozygous",fn)):
        print(f"{nm:<12}"+ "".join(f"{pct(a,q):>8.3f}" for q in (0.05,0.25,0.5,0.75,0.95)))

if __name__=="__main__":
    main()
