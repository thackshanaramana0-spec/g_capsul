#!/usr/bin/env python3
"""
CAN THIS ARCHIVE ALSO ANSWER THE EXACT-MATCH QUESTION, COMPLETELY?

Position-based retrieval is native here. Exact-match retrieval -- "every read
that CONTAINS this string", which is what a BWT/FM-index guarantees -- can be
approximated by filtering the returned reads for containment.

But there is a real gap to test. The search runs over the PSEUDOGENOME, which is
a consensus. A read carries its own deviations from that consensus, so a read
could contain the query string in a place where the consensus does not. A BWT
indexes the reads themselves and cannot miss such a read. We might.

This measures the recall directly. Ground truth is obtained by brute-force
scanning every decoded read for the probe, which is exactly what a complete
exact-match index would return. That is compared against what our query returns.

  truth      reads containing the probe, found by scanning all decoded reads
  found      of those, how many our query actually returned
  recall     found / truth

Run from ~/t34real. Uses decr/reads.out (the lossless decode) as the read set.
"""
import re, sys, os

COMP={'A':'T','C':'G','G':'C','T':'A'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))

NSITES=int(sys.argv[1]) if len(sys.argv)>1 else 50
QDIR  =os.environ.get("QDIR","qbl_het")
PROBES=os.environ.get("PROBES","bl_het.tsv")

def canon(s):
    """A read may be returned in pseudogenome orientation, i.e. reverse
    complemented relative to how it was stored in the FASTQ. Comparing raw
    strings would count those as misses and fake a recall shortfall, so both
    sides are canonicalised to min(seq, revcomp(seq))."""
    r=rc(s)
    return s if s<=r else r

def read_fa(path):
    out=set(); p=None
    try:
        for l in open(path):
            if l.startswith('>'):
                m=re.search(r'pos=(\d+)', l); p=int(m.group(1)) if m else None
            elif p is not None:
                s=l.strip()
                if s: out.add(canon(s))
    except FileNotFoundError: pass
    return out

reads=[l.strip() for l in open("decr/reads.out") if l.strip()]
print(f"reads in archive: {len(reads)}")

sites=[]
for line in open(PROBES):
    pos,ref,alt,up,down = line.rstrip("\n").split("\t")
    sites.append((pos,up,down))
    if len(sites)>=NSITES: break
print(f"sites tested: {len(sites)}")

tot_truth=tot_found=0
missed_sites=0
for pos,up,down in sites:
    pats=[up,rc(up),down,rc(down)]
    truth=set()
    for s in reads:
        for q in pats:
            if q in s:
                truth.add(canon(s)); break
    got = read_fa(f"{QDIR}/up_{pos}.fa") | read_fa(f"{QDIR}/down_{pos}.fa")
    found = len(truth & got)
    tot_truth += len(truth); tot_found += found
    if found < len(truth): missed_sites += 1

print()
print(f"reads containing a probe, by brute force : {tot_truth}")
print(f"  of those, returned by our query        : {tot_found}")
print(f"  MISSED                                 : {tot_truth-tot_found}")
print(f"  recall                                 : {tot_found/max(tot_truth,1):.4f}")
print(f"sites with at least one miss             : {missed_sites} of {len(sites)}")
print()
print("A BWT/FM-index answers this question with recall 1.0 by construction.")
print("Anything below that is the honest cost of searching a consensus rather")
print("than the reads themselves.")
