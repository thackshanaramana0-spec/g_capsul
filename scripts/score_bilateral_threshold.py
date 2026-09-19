#!/usr/bin/env python3
"""
Does bilateral 400/400 survive a stricter support rule, and do its extra false
positives die under the same rule?

The "both alleles" rule counts a site on the strength of a SINGLE supporting
read, so at 30x any one miscalled base satisfies it. Bilateral anchoring unions
two searches and therefore admits more spurious matches -- false positives rise
from 21 to 31 of 400. If 400/400 only holds because the rule is permissive, the
homozygous control will hold with it and the result is metric-gaming.

Reuses score_bilateral's own functions, so the positive and negative sets are
scored by exactly the same code that produced the headline.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_bilateral import read_fa, read_occ, variants_from, score

def collect(probes, qdir):
    out=[]
    for line in open(probes):
        pos,ref,alt,up,down = line.rstrip("\n").split("\t")
        ru=read_fa(f"{qdir}/up_{pos}.fa");   ou=read_occ(f"{qdir}/up_{pos}.occ")
        rd=read_fa(f"{qdir}/down_{pos}.fa"); od=read_occ(f"{qdir}/down_{pos}.occ")
        out.append(score(ru+rd, variants_from(ou,'up')+variants_from(od,'down'),
                         ref, alt))
    return out

def fr(s): return sorted(min(a,b)/max(a+b,1) for a,b in s)
def pct(a,q): return a[int(q*(len(a)-1))]

def main():
    het = collect("bl_het.tsv", "qbl_het")
    neg = collect("bl_neg.tsv", "qbl_neg")
    print(f"het sites: {len(het)}   homozygous sites: {len(neg)}")
    print()
    print(f"{'rule':<34}{'het both':>10}{'hom both':>10}{'separation':>13}")
    rules=[]
    for k in (1,2,3,4,5):
        rules.append((f"minor allele >= {k} read(s)", lambda a,b,k=k: min(a,b)>=k))
    for f in (0.05,0.10,0.15,0.20):
        rules.append((f"minor allele fraction >= {f:.2f}",
                      lambda a,b,f=f: (a+b)>0 and min(a,b)/(a+b)>=f))
    rules.append(("min 2 reads AND fraction >= 0.10",
                  lambda a,b: min(a,b)>=2 and (a+b)>0 and min(a,b)/(a+b)>=0.10))
    for label,rule in rules:
        hp=sum(1 for a,b in het if rule(a,b))
        hn=sum(1 for a,b in neg if rule(a,b))
        print(f"{label:<34}{hp:>10}{hn:>10}{hp-hn:>13}")
    print()
    fh, fn = fr(het), fr(neg)
    print("minor-allele fraction distribution")
    print(f"{'':<12}{'p5':>8}{'p25':>8}{'median':>8}{'p75':>8}{'p95':>8}")
    for nm,a in (("het",fh),("homozygous",fn)):
        print(f"{nm:<12}"+"".join(f"{pct(a,q):>8.3f}" for q in (0.05,0.25,0.5,0.75,0.95)))

if __name__ == "__main__":
    main()
