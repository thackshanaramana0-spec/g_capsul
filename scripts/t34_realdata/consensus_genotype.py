#!/usr/bin/env python3
"""
GENOTYPING WITHOUT DECODING READS -- exploiting two layers the query already
builds but never uses together.

THE OBSERVATION. "A heterozygous locus is not one place" cuts both ways. If the
two haplotypes sit on different contigs, then the pseudogenome CONSENSUS at
those two places already disagrees -- one carries the reference base, the other
carries the alternate. Both alleles are therefore visible in the consensus
itself, at the positions bilateral anchoring resolves, with no read decoding at
all.

THE SECOND LAYER. `capsule_decode index` already tallies
(pg position, observed base) -> count over every read and writes sites above a
threshold to <sidecar>.sites. That is a complete allele-resolved deviation
pileup, and `query` has never read it. It supplies per-allele read counts to go
with the consensus bases.

Together: allele identity from the consensus, allele support from the tally.

This script checks whether that reproduces the read-derived answer, because a
faster path that is wrong is worthless. Compared against the bilateral
read-based result on the same 400 sites.

Run from ~/t34real. Needs pg_real.fa, bl_het.tsv, qbl_het/, and a .sites file.
"""
import re, sys, os

COMP={'A':'T','C':'G','G':'C','T':'A'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))
PL=40; GAP=5
SITES = sys.argv[1] if len(sys.argv)>1 else "/tmp/full.qidx.sites"

pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')

def load_sites(path):
    """<magic:8><nsite:8> then per site: varint delta position, base byte, count byte"""
    d=open(path,'rb').read()
    if len(d)<16: return {}
    n=int.from_bytes(d[8:16],'little')
    out={}; i=16; prev=0
    for _ in range(n):
        shift=0; delta=0
        while i<len(d):
            b=d[i]; i+=1
            delta |= (b & 0x7F) << shift
            if not (b & 0x80): break
            shift+=7
        if i+1>=len(d)+1: break
        base=d[i]; cnt=d[i+1]; i+=2
        pos=prev+delta; prev=pos
        out.setdefault(pos,{})[chr(base) if isinstance(base,int) else base]=cnt
    return out

tally=load_sites(SITES)
print(f"sites tally entries: {len(tally)}")

def read_occ(path):
    out=[]
    try:
        for l in open(path):
            m=re.match(r'\[query\] occ (\d+) (\d+) ([+-])', l)
            if m: out.append((int(m.group(1)),int(m.group(2)),m.group(3)))
    except FileNotFoundError: pass
    return out

def variants_from(occs, which):
    out=[]
    for at,end,st in occs:
        if which=='up': v = at+PL+GAP if st=='+' else at-GAP-1
        else:           v = at-GAP-1  if st=='+' else at+PL+GAP
        out.append((v,st,at,end))
    return out

def read_fa(path):
    out=[]; p=None
    try:
        for l in open(path):
            if l.startswith('>'):
                m=re.search(r'pos=(\d+)',l); p=int(m.group(1)) if m else None
            elif p is not None:
                s=l.strip()
                if s: out.append((p,s))
    except FileNotFoundError: pass
    return out

def score_reads(reads, occvs, ref, alt):
    """best-overlap assignment, identical in form to score_bilateral"""
    nr=na=0
    for p,s in reads:
        best=None; bo=0
        for v,st,at,end in occvs:
            ov = min(p+len(s), end) - max(p, at)
            if ov > bo and 0 <= v-p < len(s):
                bo=ov; best=(v,st)
        if best is None:
            for v,st,at,end in occvs:
                if 0 <= v-p < len(s): best=(v,st); break
        if best is None: continue
        v,st=best; b=s[v-p]
        if st=='-': b=COMP.get(b,'N')
        if b==ref: nr+=1
        elif b==alt: na+=1
    return nr,na

def run(probes, qdir, label):
    both_c=one_c=none_c=0; agree=disagree=0
    for line in open(probes):
        pos,ref,alt,up,down = line.rstrip("\n").split("\t")
        ou=read_occ(f"{qdir}/up_{pos}.occ"); od=read_occ(f"{qdir}/down_{pos}.occ")
        vs = variants_from(ou,'up') + variants_from(od,'down')

        # ---- CONSENSUS + TALLY PATH: no reads touched ----
        seen=set()
        for v,st,at,end in vs:
            if 0 <= v < len(pg):
                b=pg[v]
                if st=='-': b=COMP.get(b,'N')
                seen.add(b)
                for db,cnt in tally.get(v,{}).items():
                    bb=COMP.get(db,'N') if st=='-' else db
                    seen.add(bb)
        c_both = (ref in seen) and (alt in seen)
        if c_both: both_c+=1
        elif (ref in seen) or (alt in seen): one_c+=1
        else: none_c+=1

        # ---- READ PATH ----
        reads = read_fa(f"{qdir}/up_{pos}.fa") + read_fa(f"{qdir}/down_{pos}.fa")
        nr,na = score_reads(reads, vs, ref, alt)
        if (nr>0 and na>0) == c_both: agree+=1
        else: disagree+=1
    print(f"=== {label} ===")
    print(f"  consensus+tally : both {both_c}  one {one_c}  neither {none_c}")
    print(f"  agreement with read path: {agree}/{agree+disagree}")
    print()
    return both_c

h=run("bl_het.tsv","qbl_het","REAL HET SITES")
g=run("bl_neg.tsv","qbl_neg","HOMOZYGOUS NEGATIVE CONTROL")
print(f"separation: {h} - {g} = {h-g}")
