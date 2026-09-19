#!/usr/bin/env python3
"""
T3.4-v2 scorer. Separated from run_locus_fidelity_v2.sh so it can be unit
tested against synthetic inputs of known ground truth -- the scoring rule is
the thing most able to silently flatter the result, so it is the thing that
most needs to be testable on its own.

Reads, from the current directory:
  probes.tsv            pos \t ref \t alt \t probe \t pgvariantpos
  q/s0_<pos>.fa .occ    sequence query, exact search
  q/s1_<pos>.fa .occ    sequence query, mismatch-tolerant search
  q/c_<pos>.fa          coordinate query

Emits the 2x2 decomposition plus the coordinate control.
"""
import sys, re

COMP = {'A':'T','C':'G','G':'C','T':'A'}

def rc(s):
    return "".join(COMP.get(c,'N') for c in reversed(s))

def read_fa(path):
    """-> [(placement_pos, sequence)] from '>rN pos=P len=L' records."""
    out=[]; p=None
    try:
        for l in open(path):
            if l.startswith('>'):
                m=re.search(r'pos=(\d+)', l); p=int(m.group(1)) if m else None
            elif p is not None:
                s=l.strip()
                if s: out.append((p,s))
    except FileNotFoundError:
        pass
    return out

def read_occ(path):
    """-> [(start,end,strand)] from '[query] occ A B +/-' lines on stderr."""
    out=[]
    try:
        for l in open(path):
            m=re.match(r'\[query\] occ (\d+) (\d+) ([+-])', l)
            if m: out.append((int(m.group(1)), int(m.group(2)), m.group(3)))
    except FileNotFoundError:
        pass
    return out

def score_containment(reads, probe, ref, alt):
    """Original rule: the read must CONTAIN the reference probe."""
    rcp=rc(probe); nr=na=0; off=len(probe)+5
    for _p,s in reads:
        i=s.find(probe)
        if i>=0 and i+off<len(s):
            b=s[i+off]
        else:
            j=s.find(rcp)
            if j>=0 and j-6>=0: b=COMP.get(s[j-6],'N')
            else: continue
        if b==ref: nr+=1
        elif b==alt: na+=1
    return nr,na

def score_offset(reads, occs, ref, alt, plen=40):
    """Placement-offset rule, identical in form to the coordinate arm.
    The probe ends 5 bp before the variant, so given a probe of length `plen`
    the variant's pg position follows from each resolved occurrence:
      + strand: probe occupies [at, at+plen), variant sits at at+plen+5
      - strand: the probe is revcomp there, so the variant sits at at-6 --
                independent of plen, since the gap is on the other side -- and
                the base must be complemented.
    plen is passed rather than assumed, because probe length is a parameter of
    the experiment: the mismatch tolerance a probe can support is bounded by
    its length (seeds must stay >= MINSEED), so a longer probe is the only way
    to search with a larger k."""
    vs=[(at+plen+5,'+') if st=='+' else (at-6,'-') for at,_e,st in occs]
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

def tally(sites, mm_label="2"):
    cells={k:[0,0,0] for k in ("coord","s0c","s1c","s0o","s1o")}
    SR=SA=0; n=0
    for pos,ref,alt,probe,v in sites:
        n+=1
        r0=read_fa(f"q/s0_{pos}.fa"); o0=read_occ(f"q/s0_{pos}.occ")
        r1=read_fa(f"q/s1_{pos}.fa"); o1=read_occ(f"q/s1_{pos}.occ")
        rcd=read_fa(f"q/c_{pos}.fa")
        P=len(probe)
        res={
            "s0c": score_containment(r0, probe, ref, alt),
            "s0o": score_offset(r0, o0, ref, alt, P),
            "s1c": score_containment(r1, probe, ref, alt),
            "s1o": score_offset(r1, o1, ref, alt, P),
            "coord": (sum(1 for p,s in rcd if 0<=v-p<len(s) and s[v-p]==ref),
                      sum(1 for p,s in rcd if 0<=v-p<len(s) and s[v-p]==alt)),
        }
        for k,(nr,na) in res.items():
            if nr and na: cells[k][0]+=1
            elif nr or na: cells[k][1]+=1
            else: cells[k][2]+=1
        SR+=res["s1o"][0]; SA+=res["s1o"][1]
    return cells, n, SR, SA

def load_probes(path="probes.tsv"):
    sites=[]
    for line in open(path):
        pos,ref,alt,probe,v = line.rstrip("\n").split("\t")
        sites.append((pos,ref,alt,probe,int(v)))
    return sites

def main():
    IND = sys.argv[1] if len(sys.argv)>1 else "HG002"
    MM  = sys.argv[2] if len(sys.argv)>2 else "2"
    cells,n,SR,SA = tally(load_probes(), MM)
    print()
    print(f"======== T3.4-v2 — {IND}: where the residue is lost ========")
    print(f"{'configuration':<34}{'both':>7}{'one':>7}{'neither':>9}")
    for key,label in (("coord","coordinate (control)"),
                      ("s0c",  "content: exact   + containment"),
                      ("s1c",  f"content: mm={MM}    + containment"),
                      ("s0o",  "content: exact   + offset"),
                      ("s1o",  f"content: mm={MM}    + offset")):
        b,o,ne=cells[key]
        print(f"{label:<34}{b:>7}{o:>7}{ne:>9}")
    print(f"sites tested: {n}")
    print(f"allele balance (mm+offset arm): {SR} REF / {SA} ALT  (ratio {SA/max(SR,1):.2f})")
    print("=============================================================")
    print(f"T34V2,{IND},{n},"
          + ",".join(f"{cells[k][0]},{cells[k][1]},{cells[k][2]}"
                     for k in ("coord","s0c","s1c","s0o","s1o")))

if __name__ == "__main__":
    main()
