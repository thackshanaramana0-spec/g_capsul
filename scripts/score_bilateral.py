#!/usr/bin/env python3
"""
Scorer for bilateral anchoring, with correct occurrence assignment.

TWO CHANGES from score_locus_fidelity_v2.py, both of which are corrections
rather than relaxations.

1. BILATERAL ANCHORING. Occurrences are collected from two reference-derived
   probes -- one ending 5 bp before the variant, one starting 5 bp after -- and
   unioned. Neither probe contains the variant, so neither can rig the result.
   The variant's pg position follows from the probe and the strand:

       up-probe    '+' -> at + PL + GAP        '-' -> at - GAP - 1
       down-probe  '+' -> at - GAP - 1         '-' -> at + PL + GAP

   and the observed base is complemented exactly when the strand is '-'.

2. BEST-OVERLAP ASSIGNMENT. The v2 scorer walked the occurrence list and took
   the FIRST whose variant position fell inside the read. With several
   occurrences within a read length of one another -- common, since the
   pseudogenome carries repeats -- that assigns a read to an occurrence it did
   not come from and reads the wrong base. Each read is now assigned to the
   occurrence its span overlaps most, which is unambiguous and does not consult
   the base being measured.

Reads bl_het.tsv / bl_neg.tsv and qbl_het/ qbl_neg/ from the current directory.
"""
import re, sys

COMP={'A':'T','C':'G','G':'C','T':'A'}
PL=40; GAP=5

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
            if m: out.append((int(m.group(1)), int(m.group(2)), m.group(3)))
    except FileNotFoundError: pass
    return out

def variants_from(occs, which):
    """-> [(variant_pg_pos, strand, probe_start, probe_end)]"""
    out=[]
    for at,end,st in occs:
        if which=='up':  v = at+PL+GAP if st=='+' else at-GAP-1
        else:            v = at-GAP-1  if st=='+' else at+PL+GAP
        out.append((v, st, at, end))
    return out

def score(reads, vs, ref, alt):
    """Assign each read to the occurrence it overlaps most, then read the base."""
    nr=na=0
    for p,s in reads:
        best=None; best_ov=0
        for v,st,at,end in vs:
            ov = min(p+len(s), end) - max(p, at)      # overlap with the probe hit
            if ov > best_ov and 0 <= v-p < len(s):
                best_ov=ov; best=(v,st)
        if best is None:
            # fall back to any occurrence covering the read, still not looking
            # at the base first
            for v,st,at,end in vs:
                if 0 <= v-p < len(s): best=(v,st); break
        if best is None: continue
        v,st = best
        b = s[v-p]
        if st=='-': b = COMP.get(b,'N')
        if b==ref: nr+=1
        elif b==alt: na+=1
    return nr,na

def run(probes, qdir, label):
    cells={'up':[0,0,0], 'down':[0,0,0], 'both':[0,0,0]}
    SR=SA=0; n=0
    for line in open(probes):
        pos,ref,alt,up,down = line.rstrip("\n").split("\t")
        n+=1
        ru=read_fa(f"{qdir}/up_{pos}.fa");   ou=read_occ(f"{qdir}/up_{pos}.occ")
        rd=read_fa(f"{qdir}/down_{pos}.fa"); od=read_occ(f"{qdir}/down_{pos}.occ")
        vu=variants_from(ou,'up'); vd=variants_from(od,'down')
        res={
            'up':   score(ru, vu, ref, alt),
            'down': score(rd, vd, ref, alt),
            'both': score(ru+rd, vu+vd, ref, alt),
        }
        for k,(a,b) in res.items():
            if a and b: cells[k][0]+=1
            elif a or b: cells[k][1]+=1
            else: cells[k][2]+=1
        SR+=res['both'][0]; SA+=res['both'][1]
    print(f"======== {label} ========")
    print(f"{'anchoring':<28}{'both':>7}{'one':>7}{'neither':>9}")
    for k,lab in (('up','upstream only (published)'),
                  ('down','downstream only'),
                  ('both','BILATERAL (union)')):
        b,o,ne=cells[k]
        print(f"{lab:<28}{b:>7}{o:>7}{ne:>9}")
    print(f"sites: {n}   allele balance (bilateral): {SR} REF / {SA} ALT"
          f"  (ratio {SA/max(SR,1):.2f})")
    print(f"CSV,{label},{n}," + ",".join(f"{cells[k][0]},{cells[k][1]},{cells[k][2]}"
                                        for k in ('up','down','both')))
    print()
    return cells

if __name__=="__main__":
    run("bl_het.tsv", "qbl_het", "BILATERAL - real HG002 het SNVs")
    run("bl_neg.tsv", "qbl_neg", "BILATERAL - homozygous negative control")
