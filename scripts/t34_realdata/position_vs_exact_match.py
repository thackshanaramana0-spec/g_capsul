#!/usr/bin/env python3
"""
LOCALITY vs CONTAINMENT — the distinction that survives the prior art.

BWT/FM-index archives (BEETL-fastq, CIndex, sFASTQ) retrieve reads that CONTAIN
the query string exactly. They store each read as an independent string, so
there is no notion of one read lying near another and no internal coordinate
system.

A pseudogenome archive has one. Reads are PLACED on a shared string, so once a
probe anchors, every read whose placement overlaps that span can be returned --
including reads that do not contain the probe at all, because they carry a
sequencing error inside it or only partially overlap it.

This measures the gap directly, on the real HG002 archive:

  returned   reads our query returns for a locus (placement overlap)
  contain    how many of those actually contain the probe exactly, which is the
             most a containment-based retrieval could return for the same query

The difference is the evidence an exact-match query cannot reach. It matters for
genotyping because those reads carry alleles too, and dropping them biases the
pileup toward the allele the probe was written from.

Run from ~/t34real after bilateral_run.sh.
"""
import re, statistics

COMP={'A':'T','C':'G','G':'C','T':'A'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))

def read_fa(path):
    out=[]; p=None
    try:
        for l in open(path):
            if l.startswith('>'):
                m=re.search(r'pos=(\d+)', l); p=int(m.group(1)) if m else None
            elif p is not None:
                s=l.strip()
                if s: out.append((p,s))
    except FileNotFoundError: pass
    return out

tot_ret=tot_con=0
per_site=[]
sites=0
for line in open("bl_het.tsv"):
    pos,ref,alt,up,down = line.rstrip("\n").split("\t")
    reads={}
    for arm,probe in (("up",up),("down",down)):
        for p,s in read_fa(f"qbl_het/{arm}_{pos}.fa"):
            reads[(p,s)]=None
    if not reads: continue
    sites+=1
    ret=len(reads)
    rcu, rcd = rc(up), rc(down)
    con=sum(1 for (p,s) in reads
            if up in s or rcu in s or down in s or rcd in s)
    tot_ret+=ret; tot_con+=con
    per_site.append((ret,con))

print(f"sites: {sites}")
print(f"reads returned by POSITION (placement overlap)  : {tot_ret}")
print(f"  of those, containing a probe exactly          : {tot_con}")
print(f"  reachable ONLY by position                    : {tot_ret-tot_con}"
      f"  ({100.0*(tot_ret-tot_con)/max(tot_ret,1):.1f}%)")
print()
fr=[c/max(r,1) for r,c in per_site]
fr.sort()
def pct(a,q): return a[int(q*(len(a)-1))]
print("per-site fraction of returned reads that contain a probe")
print(f"  p5 {pct(fr,0.05):.3f}   p25 {pct(fr,0.25):.3f}   median {pct(fr,0.5):.3f}"
      f"   p75 {pct(fr,0.75):.3f}   p95 {pct(fr,0.95):.3f}")
print()
print("An exact-match archive can return at best the 'contain' column for the")
print("same query. The remainder is evidence only an internal coordinate")
print("system can reach.")
