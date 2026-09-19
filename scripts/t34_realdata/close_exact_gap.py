#!/usr/bin/env python3
"""
CAN THE EXACT-MATCH GAP BE CLOSED, AND AT WHAT COST?

Established: all 34 reads missed by exact-match retrieval have literal sequences
that do not appear in the pseudogenome. A read whose sequence is absent from the
consensus MUST carry deviations -- if it had none it would BE a pg substring.

So every possible miss is confined to the deviation-carrying subset, which the
sidecar already enumerates (16,415 of 117,565 reads here, 14%). Scanning only
that subset should recover the full 1.000 recall at 14% of the cost of scanning
everything.

This tests the claim instead of asserting it: it takes the reads the query
returned, adds a scan restricted to deviation-carrying reads, and re-measures
recall against brute force over ALL reads.
"""
import re, sys, os, time

COMP={'A':'T','C':'G','G':'C','T':'A'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))
def canon(s):
    r=rc(s); return s if s<=r else r

NS=int(sys.argv[1]) if len(sys.argv)>1 else 50
pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')
reads=[l.strip() for l in open("decr/reads.out") if l.strip()]

# the deviation-carrying subset: a read absent from the consensus as a literal
# substring is exactly a read that deviates from it
t0=time.time()
dev=[s for s in reads if pg.find(s)<0 and pg.find(rc(s))<0]
t_part=time.time()-t0
print(f"all reads                : {len(reads)}")
print(f"deviation-carrying subset: {len(dev)}  ({100.0*len(dev)/len(reads):.1f}%)  [{t_part:.1f}s to derive]")

def read_fa(path):
    out=set(); p=None
    try:
        for l in open(path):
            if l.startswith('>'):
                m=re.search(r'pos=(\d+)',l); p=int(m.group(1)) if m else None
            elif p is not None:
                x=l.strip()
                if x: out.add(canon(x))
    except FileNotFoundError: pass
    return out

sites=[]
for line in open("bl_het.tsv"):
    pos,ref,alt,up,down=line.rstrip("\n").split("\t")
    sites.append((pos,up,down))
    if len(sites)>=NS: break

tt=tf_q=tf_aug=0
t_scan=0.0
for pos,up,down in sites:
    pats=[up,rc(up),down,rc(down)]
    truth=set()
    for s in reads:
        for q in pats:
            if q in s: truth.add(canon(s)); break
    got=read_fa(f"qbl_het/up_{pos}.fa")|read_fa(f"qbl_het/down_{pos}.fa")
    t0=time.time()
    extra=set()
    for s in dev:
        for q in pats:
            if q in s: extra.add(canon(s)); break
    t_scan+=time.time()-t0
    tt+=len(truth); tf_q+=len(truth&got); tf_aug+=len(truth&(got|extra))

print()
print(f"reads containing a probe (brute force over ALL reads): {tt}")
print(f"  query alone                : {tf_q}   recall {tf_q/max(tt,1):.4f}")
print(f"  query + deviation-subset scan: {tf_aug}   recall {tf_aug/max(tt,1):.4f}")
print(f"  added scan cost: {t_scan/len(sites)*1000:.1f} ms/site over {len(dev)} reads")
