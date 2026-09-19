#!/usr/bin/env python3
"""
CLOSING THE EXACT-MATCH GAP WITHOUT PAYING 41.6 ms/site.

The naive closure scans all 16,481 deviation-carrying reads per query. That is
linear in the subset for every site and is the wrong shape: the subset does not
change between queries, so the work belongs at BUILD time, not query time.

Index: a k-mer -> read-id map over the subset only.

Why one k-mer suffices. A read that CONTAINS a 40 bp probe contains every one of
the probe's k-mers. So a single k-mer lookup yields a candidate set guaranteed
to include every true hit, and verification then removes false candidates. No
recall is given up -- the lookup is a superset filter, not an approximation.

Both strands are handled by querying the probe and its reverse complement.

Reports build cost, memory, per-site query cost and recall, against the same
brute-force ground truth over ALL reads.
"""
import re, sys, time, os

K = int(os.environ.get("IDXK", "20"))
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 50

COMP = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return "".join(COMP.get(c,'N') for c in reversed(s))
def canon(s):
    r = rc(s); return s if s <= r else r

pg = "".join(l.strip() for l in open("pg_real.fa") if l[0] != '>')
reads = [l.strip() for l in open("decr/reads.out") if l.strip()]

# ---- the subset: reads whose literal sequence is absent from the consensus ----
sub_cache = "subset_absent.txt"
if os.path.exists(sub_cache):
    dev = [l.strip() for l in open(sub_cache) if l.strip()]
    print(f"subset loaded from cache: {len(dev)} reads")
else:
    t0 = time.time()
    dev = [s for s in reads if pg.find(s) < 0 and pg.find(rc(s)) < 0]
    open(sub_cache, "w").write("\n".join(dev) + "\n")
    print(f"subset derived: {len(dev)} reads in {time.time()-t0:.1f}s (cached)")

print(f"subset is {100.0*len(dev)/len(reads):.1f}% of {len(reads)} reads")

# ---- BUILD: k-mer -> read ids, over the subset only ----
t0 = time.time()
idx = {}
for i, s in enumerate(dev):
    for j in range(0, len(s) - K + 1):
        idx.setdefault(s[j:j+K], []).append(i)
t_build = time.time() - t0
npost = sum(len(v) for v in idx.values())
print(f"index built: {len(idx)} distinct {K}-mers, {npost} postings, {t_build:.1f}s")
try:
    import resource
    print(f"peak RSS after build: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f} MB")
except Exception:
    pass

def lookup(probe):
    """every subset read containing `probe`, via one k-mer then verify"""
    out = set()
    for p in (probe, rc(probe)):
        if len(p) < K: continue
        cands = idx.get(p[:K])
        if not cands: continue
        for i in cands:
            if p in dev[i]: out.add(canon(dev[i]))
    return out

def read_fa(path):
    out = set(); pos = None
    try:
        for l in open(path):
            if l.startswith('>'):
                m = re.search(r'pos=(\d+)', l); pos = int(m.group(1)) if m else None
            elif pos is not None:
                x = l.strip()
                if x: out.add(canon(x))
    except FileNotFoundError: pass
    return out

sites = []
for line in open("bl_het.tsv"):
    p, ref, alt, up, down = line.rstrip("\n").split("\t")
    sites.append((p, up, down))
    if len(sites) >= NS: break

tt = fq = faug = 0
t_q = 0.0
for pos, up, down in sites:
    pats = [up, rc(up), down, rc(down)]
    truth = set()
    for s in reads:
        for q in pats:
            if q in s: truth.add(canon(s)); break
    got = read_fa(f"qbl_het/up_{pos}.fa") | read_fa(f"qbl_het/down_{pos}.fa")
    t0 = time.time()
    extra = lookup(up) | lookup(down)
    t_q += time.time() - t0
    tt += len(truth); fq += len(truth & got); faug += len(truth & (got | extra))

n = len(sites)
print()
print(f"reads containing a probe (brute force, ALL reads): {tt}")
print(f"  query alone            : {fq}   recall {fq/max(tt,1):.4f}")
print(f"  query + INDEXED lookup : {faug}   recall {faug/max(tt,1):.4f}")
print(f"  indexed lookup cost    : {t_q/n*1000:.3f} ms/site   (naive scan was 41.6 ms/site)")
print(f"  speedup on the closure : {41.6/max(t_q/n*1000, 1e-9):.0f}x")
