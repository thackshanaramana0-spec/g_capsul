#!/usr/bin/env python3
"""
Diagnose the sites the mm+offset arm still misses on REAL reads.

The simulated run reached 400/400 and the real run reaches 391/400. The nine
differing sites are the only place where real data disagrees with the model, so
they are the only place a NEW mechanism can be hiding. This characterises them
against the sites that succeeded, instead of assuming they are noise.

Run from ~/t34real after run_realreads_t34.sh.
"""
import re, collections

COMP={'A':'T','C':'G','G':'C','T':'A'}

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

def read_occ(path):
    out=[]
    try:
        for l in open(path):
            m=re.match(r'\[query\] occ (\d+) (\d+) ([+-])', l)
            if m: out.append((int(m.group(1)),int(m.group(2)),m.group(3)))
    except FileNotFoundError: pass
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

ref_seq=open("ref.txt").read().strip(); REF0=2999001

# GIAB variant density around each site, to test the "second variant in window"
# story directly on the failures.
var={}
for line in open("giab_win.vcf"):
    if line.startswith("#"): continue
    f=line.split("\t")
    var[int(f[1])]=(f[3],f[4],f[9].split(":")[0])

fail=[]; ok=[]
for line in open("probes_real.tsv"):
    pos,r,a,probe,v = line.rstrip("\n").split("\t"); pos=int(pos)
    r0=read_fa(f"qreal/s0_{pos}.fa"); o0=read_occ(f"qreal/s0_{pos}.occ")
    r1=read_fa(f"qreal/s1_{pos}.fa"); o1=read_occ(f"qreal/s1_{pos}.occ")
    nr,na = score_offset(r1,o1,r,a)
    rec = dict(pos=pos, ref=r, alt=a, nr=nr, na=na,
               occ0=len(o0), occ1=len(o1),
               reads0=len(r0), reads1=len(r1),
               nearby=sum(1 for d in range(-60,61) if (pos+d) in var and d!=0))
    (fail if (nr==0 or na==0) else ok).append(rec)

print(f"succeeded: {len(ok)}   failed: {len(fail)}")
print()
print("FAILING SITES")
print(f"{'pos':>10}{'ref':>4}{'alt':>4}{'REF rds':>9}{'ALT rds':>9}"
      f"{'occ k=0':>9}{'occ k=2':>9}{'reads':>8}{'GIAB vars +/-60':>17}")
for d in sorted(fail, key=lambda x:x['pos']):
    print(f"{d['pos']:>10}{d['ref']:>4}{d['alt']:>4}{d['nr']:>9}{d['na']:>9}"
          f"{d['occ0']:>9}{d['occ1']:>9}{d['reads1']:>8}{d['nearby']:>17}")

def avg(xs,k): return sum(x[k] for x in xs)/max(len(xs),1)
print()
print(f"{'':<22}{'failed':>10}{'succeeded':>12}")
for k,label in (("occ0","occurrences k=0"),("occ1","occurrences k=2"),
                ("reads1","reads returned"),("nearby","GIAB vars within 60bp")):
    print(f"{label:<22}{avg(fail,k):>10.2f}{avg(ok,k):>12.2f}")

# how many failures returned only ONE occurrence -- i.e. the second contig was
# never found even with tolerance
solo=sum(1 for d in fail if d['occ1']<=1)
print()
print(f"failures where mm=2 still found <=1 occurrence: {solo} of {len(fail)}")
print(f"failures with ZERO alt reads despite >1 occurrence: "
      f"{sum(1 for d in fail if d['na']==0 and d['occ1']>1)}")
