#!/usr/bin/env python3
"""
Is the damage one-sided?

The upstream window at the alt occurrences diverges by 6-13 bases from the
reference probe -- far more than one haplotype difference. The pseudogenome is
built by greedy overlap chaining, so a contig's neighbourhood is an artefact of
which reads overlapped, not the genomic neighbourhood. A reference-derived probe
anchored upstream can then miss the alt contig at any tolerance.

If the chaining break is on ONE side, the other side is still reference-like,
and a probe placed there would find the same occurrence. This measures both
sides at every alt occurrence, for the failing sites AND for a control group of
successful ones, so "downstream is better" cannot be mistaken for a property of
all sites rather than of the failures.

Neither probe contains the variant, so neither can rig the result.

Run from ~/t34real.
"""
COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))

REF0=2999001
FAIL=[3036066,3053969,3059523,3163883,3172647,3172774,3325265,3345618,3349190]
HALF=20; PL=40; GAP=5

ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')
V={}
for line in open("giab_win.vcf"):
    if line.startswith("#"): continue
    f=line.split("\t"); V[int(f[1])]=(f[3],f[4])

def mism(a,b): return sum(1 for x,y in zip(a,b) if x!=y) if len(a)==len(b) else 999

def sides(pos):
    """-> list of (strand, pg_at, up_mismatches, down_mismatches)"""
    r,a = V[pos]; i=pos-REF0
    up_ref   = ref[i-GAP-PL : i-GAP]          # 40bp ending 5bp BEFORE variant
    down_ref = ref[i+1+GAP : i+1+GAP+PL]      # 40bp starting 5bp AFTER variant
    ctx = ref[i-HALF:i] + a + ref[i+1:i+HALF+1]
    out=[]
    for orient, needle in (('+',ctx), ('-',rc(ctx))):
        s=0
        while True:
            at = pg.find(needle, s)
            if at<0: break
            s = at+1
            if orient=='+':
                up   = pg[at-HALF-GAP-PL+0 : at-HALF-GAP]
                # variant sits at at+HALF, so downstream window starts at +HALF+1+GAP
                down = pg[at+HALF+1+GAP : at+HALF+1+GAP+PL]
                mu, md = mism(up,up_ref), mism(down,down_ref)
            else:
                # minus strand: forward-upstream maps AFTER the block, and
                # forward-downstream maps BEFORE it, both revcomped
                up   = pg[at+HALF+1+GAP : at+HALF+1+GAP+PL]
                down = pg[at-HALF-GAP-PL : at-HALF-GAP]
                mu, md = mism(up,rc(up_ref)), mism(down,rc(down_ref))
            out.append((orient,at,mu,md))
    return out

def summarise(title, positions):
    print(f"=== {title} ===")
    print(f"{'pos':>10}{'occurrences':>13}{'best up':>9}{'best down':>11}  reachable at k=2 from")
    tally={'up':0,'down':0,'both':0,'neither':0,'no occ':0}
    for pos in positions:
        occ=sides(pos)
        if not occ:
            print(f"{pos:>10}{0:>13}{'-':>9}{'-':>11}  (ALT context absent from pg)")
            tally['no occ']+=1; continue
        bu=min(o[2] for o in occ); bd=min(o[3] for o in occ)
        u_ok, d_ok = bu<=2, bd<=2
        who = 'both' if u_ok and d_ok else 'down only' if d_ok else 'up only' if u_ok else 'NEITHER'
        if u_ok and d_ok: tally['both']+=1
        elif d_ok: tally['down']+=1
        elif u_ok: tally['up']+=1
        else: tally['neither']+=1
        print(f"{pos:>10}{len(occ):>13}{bu:>9}{bd:>11}  {who}")
    print("  tally:", tally)
    print()

summarise("the 9 FAILING sites", FAIL)
succ=[int(l.split('\t')[0]) for l in open("probes_real.tsv")]
succ=[p for p in succ if p not in FAIL][:15]
summarise("control: 15 SUCCEEDING sites", succ)
