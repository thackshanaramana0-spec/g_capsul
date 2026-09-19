#!/usr/bin/env python3
"""
The 34 reads that contain a probe but are not returned -- why?

Tolerance is ruled out: recall is identical at MINSEED 12/8/6/5 (k_max 2/4/5/7).
The remaining candidate explanation is structural -- the read sits somewhere the
CONSENSUS does not carry the probe at all, at any tolerance, because the read's
own deviations are what spell the probe.

Test: for each missed read, find where its sequence occurs in the pseudogenome
(it may not occur at all, since the pg is a consensus), and report the best
mismatch distance between the probe and the pg at the read's own span.

  pg has the probe within k   -> we should have found it, this is a real defect
  pg is far from the probe    -> architecture, not a defect
  read not placeable on pg    -> the read is stored via deviations, as expected
"""
import re, sys
COMP={'A':'T','C':'G','G':'C','T':'A'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))
def canon(s):
    r=rc(s); return s if s<=r else r

pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')

def read_fa(path):
    out=set(); p=None
    try:
        for l in open(path):
            if l.startswith('>'):
                m=re.search(r'pos=(\d+)',l); p=int(m.group(1)) if m else None
            elif p is not None:
                s=l.strip()
                if s: out.add(canon(s))
    except FileNotFoundError: pass
    return out

reads=[l.strip() for l in open("decr/reads.out") if l.strip()]

def best_dist(hay_sub, pat):
    """min mismatches of pat against any offset in hay_sub"""
    best=len(pat)
    for i in range(0, max(1,len(hay_sub)-len(pat)+1)):
        d=sum(1 for a,b in zip(hay_sub[i:i+len(pat)],pat) if a!=b)
        if d<best: best=d
    return best

sites=[]
for line in open("bl_het.tsv"):
    pos,ref,alt,up,down=line.rstrip("\n").split("\t")
    sites.append((pos,up,down))
    if len(sites)>=50: break

cat={"probe_in_pg_near":0,"probe_far_from_pg":0,"read_not_in_pg":0}
examples=[]
for pos,up,down in sites:
    pats=[up,rc(up),down,rc(down)]
    truth={}
    for s in reads:
        for q in pats:
            if q in s: truth[canon(s)]=s; break
    got=read_fa(f"qbl_het/up_{pos}.fa")|read_fa(f"qbl_het/down_{pos}.fa")
    for c,orig in truth.items():
        if c in got: continue
        # where does this read sit on the pg?
        at=pg.find(orig)
        if at<0: at=pg.find(rc(orig))
        if at<0:
            cat["read_not_in_pg"]+=1
            if len(examples)<3: examples.append((pos,"read sequence absent from pg consensus"))
            continue
        span=pg[max(0,at-60):at+len(orig)+60]
        d=min(best_dist(span,up),best_dist(span,down))
        if d<=2:
            cat["probe_in_pg_near"]+=1
            if len(examples)<3: examples.append((pos,f"probe is {d} mismatches from pg here -- SHOULD have been found"))
        else:
            cat["probe_far_from_pg"]+=1
            if len(examples)<3: examples.append((pos,f"probe is {d} mismatches from pg at the read's span"))

print("why the missed reads were missed:")
for k,v in cat.items(): print(f"  {k:<22} {v}")
print()
for e in examples: print(" ",e)
