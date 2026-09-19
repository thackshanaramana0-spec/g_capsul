#!/usr/bin/env python3
"""
THE TALLY CARRIES COUNTS, AND THE FIRST PASS IGNORED THEM.

consensus_genotype.py accepted an allele if it appeared AT ALL -- either as the
consensus base at an anchored position, or as any deviation recorded in the
sites tally. Presence only. That reached 400/400 but at 39 false positives on
the homozygous control, against 31 for the read-decoding path.

But the tally stores (position, base) -> COUNT. A single deviation at a position
is overwhelmingly a sequencing error. Requiring a minimum count before believing
a deviation should cut false positives without touching real het sites, where
the alternate allele is carried by roughly half the reads.

Consensus bases are always trusted -- they are the assembled agreement of many
reads, not a single observation. Only tally-derived alleles are thresholded.

This sweeps that threshold on both sets. If the fast path can reach the read
path's false-positive rate while staying at 400/400, it is strictly better:
19.6x faster AND no worse.

Run from ~/t34real.
"""
import re, sys

COMP={'A':'T','C':'G','G':'C','T':'A'}
PL=40; GAP=5
SITES = sys.argv[1] if len(sys.argv)>1 else "full.qidx.sites"
pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')

def load_sites(path):
    d=open(path,'rb').read(); n=int.from_bytes(d[8:16],'little')
    out={}; i=16; prev=0
    for _ in range(n):
        shift=0; delta=0
        while i<len(d):
            b=d[i]; i+=1
            delta |= (b & 0x7F) << shift
            if not (b & 0x80): break
            shift+=7
        base=d[i]; cnt=d[i+1]; i+=2
        pos=prev+delta; prev=pos
        out.setdefault(pos,{})[chr(base)]=cnt
    return out

def read_occ(p):
    o=[]
    try:
        for l in open(p):
            m=re.match(r'\[query\] occ (\d+) (\d+) ([+-])',l)
            if m: o.append((int(m.group(1)),int(m.group(2)),m.group(3)))
    except FileNotFoundError: pass
    return o

def vf(occs,which):
    out=[]
    for at,end,st in occs:
        v = (at+PL+GAP if st=='+' else at-GAP-1) if which=='up' else (at-GAP-1 if st=='+' else at+PL+GAP)
        out.append((v,st))
    return out

tally=load_sites(SITES)

def both_at(probes,qdir,minc):
    n=0
    for line in open(probes):
        pos,ref,alt,up,down = line.rstrip("\n").split("\t")
        vs = vf(read_occ(f"{qdir}/up_{pos}.occ"),'up') + vf(read_occ(f"{qdir}/down_{pos}.occ"),'down')
        seen=set()
        for v,st in vs:
            if not (0 <= v < len(pg)): continue
            b=pg[v]
            seen.add(COMP.get(b,'N') if st=='-' else b)     # consensus always trusted
            for db,cnt in tally.get(v,{}).items():
                if cnt >= minc:                              # deviations thresholded
                    seen.add(COMP.get(db,'N') if st=='-' else db)
        if ref in seen and alt in seen: n+=1
    return n

print(f"{'min deviation count':<22}{'het both':>10}{'hom false':>11}{'separation':>13}")
best=None
for minc in (1,2,3,4,5,6,8,10):
    h=both_at("bl_het.tsv","qbl_het",minc)
    g=both_at("bl_neg.tsv","qbl_neg",minc)
    tag=""
    if h==400 and (best is None or g<best[1]):
        best=(minc,g); tag="  <- 400/400 at fewest false positives"
    print(f"{minc:<22}{h:>10}{g:>11}{h-g:>13}{tag}")
print()
if best: print(f"best operating point: deviations need >= {best[0]} reads -> 400/400 with {best[1]} false positives")
print("read-decoding path for comparison: 400/400 with 31 false positives")
