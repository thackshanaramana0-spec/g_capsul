#!/usr/bin/env python3
"""
For the nine residual sites, separate STORAGE from RETRIEVAL.

The archive round-trips byte-identically, so the ALT-bearing reads are
necessarily inside it. This confirms that directly and then asks the next
question down: is the ALT context also present in the exported PSEUDOGENOME?

  ALT in decoded reads, ALT in pg      -> the pg has the alternate contig, so
                                          the search should have reached it
  ALT in decoded reads, NOT in pg      -> the alternate allele survives only as
                                          per-read deviations, so a probe
                                          matched against the pg consensus can
                                          never locate it by sequence, at any
                                          mismatch tolerance

The second case is a genuine structural limit of sequence addressing, not a
tuning failure, and it is a different mechanism from the one already fixed.

Run from ~/t34real.
"""
COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))

REF0=2999001
FAIL=[3036066,3053969,3059523,3163883,3172647,3172774,3325265,3345618,3349190]

ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')
reads=open("decr/reads.out").read().split("\n")

V={}
for line in open("giab_win.vcf"):
    if line.startswith("#"): continue
    f=line.split("\t"); V[int(f[1])]=(f[3],f[4])

succ=[int(l.split("\t")[0]) for l in open("probes_real.tsv")]
succ=[p for p in succ if p not in FAIL][:9]

def probe_ctx(pos, allele, half=15):
    i=pos-REF0
    return ref[i-half:i] + allele + ref[i+1:i+half+1]

def count_in_reads(ctx):
    r=rc(ctx); n=0
    for s in reads:
        if ctx in s or r in s: n+=1
    return n

def in_pg(ctx):
    return (pg.find(ctx)>=0) or (pg.find(rc(ctx))>=0)

def report(title, positions):
    print(f"=== {title} ===")
    print(f"{'pos':>10}{'ALT reads':>11}{'REF reads':>11}{'ALT in pg':>11}{'REF in pg':>11}")
    nalt_pg=0
    for pos in positions:
        r,a=V[pos]
        ca=probe_ctx(pos,a); cr=probe_ctx(pos,r)
        na=count_in_reads(ca); nr=count_in_reads(cr)
        pa=in_pg(ca); pr=in_pg(cr)
        if pa: nalt_pg+=1
        print(f"{pos:>10}{na:>11}{nr:>11}{str(pa):>11}{str(pr):>11}")
    print(f"-> ALT context present in pseudogenome: {nalt_pg} of {len(positions)}")
    print()

report("the 9 residual sites", FAIL)
report("control: 9 sites that succeeded", succ)
