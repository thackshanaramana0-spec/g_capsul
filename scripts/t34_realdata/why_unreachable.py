#!/usr/bin/env python3
"""
WHY is the alternate contig unreachable at the nine failing sites?

Tolerance is ruled out (probe length swept to 120 bp, k to 9, no change). The
ALT context is present in the pseudogenome for 6 of 9. So the alt contig exists
and our probe still does not land on it.

The probe is 40 bp of reference ending 5 bp BEFORE the variant -- entirely
upstream. Two things could then go wrong, and they need different fixes:

  DIVERGENCE      the upstream window exists on the alt contig but differs from
                  the reference probe in more than k places -> a tolerance
                  problem after all, just a bigger one than tested.

  BOUNDARY        the upstream window is not on the alt contig at all, because
                  the contig starts at or after the variant -> no mismatch
                  tolerance can ever help, and the fix has to be to anchor the
                  search somewhere other than upstream.

This measures which, by locating the ALT context in the pg and then comparing
the 40 bp immediately upstream OF THAT against the reference probe.

Run from ~/t34real.
"""
COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))

REF0=2999001
FAIL=[3036066,3053969,3059523,3163883,3172647,3172774,3325265,3345618,3349190]

ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')

V={}
for line in open("giab_win.vcf"):
    if line.startswith("#"): continue
    f=line.split("\t"); V[int(f[1])]=(f[3],f[4])

probes={}
for line in open("probes_real.tsv"):
    pos,r,a,p,v = line.rstrip("\n").split("\t")
    probes[int(pos)]=p

def mism(a,b):
    return sum(1 for x,y in zip(a,b) if x!=y) if len(a)==len(b) else None

def analyse(pos):
    r,a = V[pos]
    i = pos-REF0
    probe = probes[pos]
    HALF=20
    # ALT-haplotype context: reference around the variant, with ALT substituted
    ctx = ref[i-HALF:i] + a + ref[i+1:i+HALF+1]
    rows=[]
    for orient,needle in (('+',ctx),('-',rc(ctx))):
        start=0
        while True:
            at = pg.find(needle, start)
            if at<0: break
            start = at+1
            # where does the 40bp upstream window sit, on this orientation?
            if orient=='+':
                up = pg[at-HALF-45 : at-HALF-5]     # aligns with the ref probe
                cmpprobe = probe
            else:
                # on the minus strand the upstream window is downstream in pg
                up = pg[at+HALF+1+5 : at+HALF+1+45]
                cmpprobe = rc(probe)
            if len(up)!=40:
                rows.append((orient,at,'OFF THE END OF THE PG',None)); continue
            rows.append((orient,at,None,mism(up,cmpprobe)))
    return r,a,rows

print("For each failing site: every place the ALT context occurs in the pg, and")
print("how far the 40bp upstream window there is from the reference probe.")
print()
print(f"{'pos':>10}{'ref':>4}{'alt':>4}{'strand':>8}{'pg pos':>10}{'upstream mismatches vs probe':>30}")
summary={'boundary':0,'divergent':0,'close':0,'absent':0}
for pos in FAIL:
    r,a,rows = analyse(pos)
    if not rows:
        print(f"{pos:>10}{r:>4}{a:>4}{'-':>8}{'-':>10}{'ALT CONTEXT NOT IN PG':>30}")
        summary['absent']+=1
        continue
    best=None
    for orient,at,note,m in rows:
        shown = note if note else f"{m}"
        print(f"{pos:>10}{r:>4}{a:>4}{orient:>8}{at:>10}{shown:>30}")
        if m is not None and (best is None or m<best): best=m
    if best is None:                     summary['boundary']+=1
    elif best<=2:                        summary['close']+=1
    else:                                summary['divergent']+=1
print()
print("verdict tally:", summary)
print()
print("close     = upstream window is within k=2, so the search SHOULD have found it")
print("divergent = upstream window exists but is too far -> tolerance problem")
print("boundary  = upstream window runs off the pseudogenome -> anchoring problem")
