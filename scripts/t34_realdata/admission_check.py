#!/usr/bin/env python3
"""
IS THE DENOMINATOR FAIR?

Sites are admitted to the experiment if the UPSTREAM probe is found exactly in
the pseudogenome. That rule was inherited from the published run and kept
deliberately, so that bilateral and upstream-only are compared over the same
400 sites.

But it invites an objection: if admission depends on the upstream probe, the
site set might be biased toward sites where upstream anchoring works, and
bilateral's gain could be an artefact of that.

This measures the objection directly. For every het SNV in the window it asks
which probes anchor at all:

   up only     admitted today
   down only   EXCLUDED today -- these are the sites the objection is about
   both        admitted today
   neither     excluded under any rule

If "down only" is large, the current denominator hides sites where bilateral
would help most, and the reported gain is an UNDERSTATEMENT rather than an
inflation. If it is ~0, the denominator is not doing any work and the comparison
is clean either way.

Run from ~/t34real.
"""
REF0=2999001; PL=40; GAP=5
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg_real.fa") if l[0]!='>')

n=up_only=down_only=both=neither=0
for line in open("sites.tsv"):
    pos,r,a=line.split(); pos=int(pos); i=pos-REF0
    u=ref[i-GAP-PL:i-GAP]; d=ref[i+1+GAP:i+1+GAP+PL]
    if len(u)!=PL or len(d)!=PL: continue
    if set(u)-set("ACGT") or set(d)-set("ACGT"): continue
    n+=1
    hu = pg.find(u)>=0
    hd = pg.find(d)>=0
    if hu and hd: both+=1
    elif hu: up_only+=1
    elif hd: down_only+=1
    else: neither+=1

print(f"het SNVs with well-formed probes : {n}")
print(f"  both probes anchor exactly      : {both}")
print(f"  upstream only  (admitted today) : {up_only}")
print(f"  DOWNSTREAM ONLY (excluded today): {down_only}   <- the objection")
print(f"  neither anchors                 : {neither}")
print()
admitted = both + up_only
either   = both + up_only + down_only
print(f"denominator under the current rule (upstream anchors): {admitted}")
print(f"denominator under an either-anchor rule             : {either}")
print(f"sites the current rule excludes that bilateral could serve: {down_only}")
