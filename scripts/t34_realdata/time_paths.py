#!/usr/bin/env python3
"""Scoring-phase cost: read decoding vs consensus+tally lookup, 400 sites.

The search phase is shared by both paths, so this times only what differs --
parsing and scoring returned reads, against looking up the consensus base and
the deviation tally. Reported per site and in total.
"""
import re, time

COMP={'A':'T','C':'G','G':'C','T':'A'}
PL=40; GAP=5
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
        out.append((v,st,at,end))
    return out

sites=[l.rstrip("\n").split("\t") for l in open("bl_het.tsv")]

t0=time.time(); tally=load_sites("full.qidx.sites"); t_load=time.time()-t0

# --- consensus + tally ---
t0=time.time()
for pos,ref,alt,up,down in sites:
    vs=vf(read_occ(f"qbl_het/up_{pos}.occ"),'up')+vf(read_occ(f"qbl_het/down_{pos}.occ"),'down')
    seen=set()
    for v,st,at,end in vs:
        if 0<=v<len(pg):
            b=pg[v]
            if st=='-': b=COMP.get(b,'N')
            seen.add(b)
            for db in tally.get(v,{}):
                seen.add(COMP.get(db,'N') if st=='-' else db)
    _=(ref in seen and alt in seen)
t_cons=time.time()-t0

# --- read decoding ---
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

t0=time.time()
for pos,ref,alt,up,down in sites:
    vs=vf(read_occ(f"qbl_het/up_{pos}.occ"),'up')+vf(read_occ(f"qbl_het/down_{pos}.occ"),'down')
    reads=read_fa(f"qbl_het/up_{pos}.fa")+read_fa(f"qbl_het/down_{pos}.fa")
    nr=na=0
    for p,s in reads:
        best=None; bo=0
        for v,st,at,end in vs:
            ov=min(p+len(s),end)-max(p,at)
            if ov>bo and 0<=v-p<len(s): bo=ov; best=(v,st)
        if best is None: continue
        v,st=best; b=s[v-p]
        if st=='-': b=COMP.get(b,'N')
        if b==ref: nr+=1
        elif b==alt: na+=1
t_read=time.time()-t0

n=len(sites)
print(f"sites: {n}")
print(f"tally load (once)      : {t_load*1000:7.1f} ms   [124 KB]")
print(f"consensus+tally scoring: {t_cons*1000:7.1f} ms   ({t_cons/n*1000:.3f} ms/site)")
print(f"read-decoding scoring  : {t_read*1000:7.1f} ms   ({t_read/n*1000:.3f} ms/site)")
print(f"speedup on the scoring phase: {t_read/max(t_cons,1e-9):.1f}x")
