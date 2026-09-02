import collections, math
L=151
# ---- ground truth ----
meta={}
for line in open('ecoli_sim.meta'):
    i,s,st=map(int,line.split()); meta[i]=(s,st)
truthE=set()
for line in open('ecoli_sim.truth'):
    a=line.split()
    if a[2]=='E': truthE.add((int(a[0]),int(a[1])))
varloci=set(int(line.split()[0]) for line in open('ecoli_sim.varloci'))
def gpos(rid,pos):
    s,st=meta[rid]; return s+(L-1-pos if st else pos)

# ---- Pass 1: per-column base counts; store per-read (base,e) for candidates ----
cnt=collections.defaultdict(lambda:[0,0,0,0])
qsum=collections.defaultdict(lambda:[0.0,0.0,0.0,0.0])  # sum of (1-e) weight per base
reads_at=collections.defaultdict(list)   # col -> list of (base, e)  (all reads)
mm_rows=[]                                # (rid,pos,col,rb)
depth_at={}
for line in open('sim_pileup.tsv'):
    if line[0]=='r': continue
    v=line.split()
    rid=int(v[0]); q=int(v[1]); mm=int(v[8]); pos=int(v[9]); col=int(v[10]); rb=int(v[12])
    if col<0 or rb>3: continue
    e=10**(-q/10.0)
    cnt[col][rb]+=1
    reads_at[col].append((rb,e))
    if mm: mm_rows.append((rid,pos,col,rb))

# ---- Genotype-likelihood variant caller per column ----
def loglik(reads, alleles):
    # alleles: list of allele bases; genotype = uniform mixture (hom=1 allele, het=2)
    m=len(alleles); ll=0.0
    for b,e in reads:
        p=0.0
        for a in alleles:
            p += (1-e) if b==a else (e/3.0)
        ll += math.log(max(p/m,1e-300))
    return ll

def call_variants(QUAL_thr):
    var_alleles={}   # col -> set(alleles) if variant
    for col,rl in reads_at.items():
        c=cnt[col]; d=sum(c)
        if d<6: continue
        order=sorted(range(4), key=lambda b:-c[b])
        M,mn=order[0],order[1]
        if c[mn]<1: continue
        ll_homM = loglik(rl,[M])
        ll_het  = loglik(rl,[M,mn])
        ll_homm = loglik(rl,[mn])
        best=max(ll_het,ll_homm)
        qual = 10.0*(best-ll_homM)/math.log(10)
        if qual>=QUAL_thr:
            var_alleles[col]= {M,mn} if ll_het>=ll_homm else {mn}
    return var_alleles

def near(g,S): return any((g+d) in S for d in (-1,0,1))
print("QUAL | VarP  VarR | ErrP  ErrR | #varcols")
for thr in [10,20,30,50,80]:
    va=call_variants(thr)
    # classify mismatches
    callE=set()
    for rid,pos,col,rb in mm_rows:
        if col in va and rb in va[col]: pass   # variant, not an error
        else: callE.add((rid,pos))
    tp=len(callE&truthE); fp=len(callE-truthE)
    pe=tp/max(tp+fp,1); re=tp/max(len(truthE),1)
    # variant loci (genome)
    cv=set()
    for col in va:
        gs=collections.Counter(gpos(r,p) for (r,p,c,b) in mm_rows if c==col)
        if gs: cv.add(gs.most_common(1)[0][0])
    tpv=sum(1 for g in cv if near(g,varloci))
    pv=tpv/max(len(cv),1); rv=sum(1 for g in varloci if near(g,cv))/len(varloci)
    print(f" {thr:3d} | {pv:.3f} {rv:.3f} | {pe:.3f} {re:.3f} | {len(cv)}")

# callable ceiling: fraction of truth at depth>=6
callableV=sum(1 for g in varloci)  # all in region; check depth via any col mapping ~ skip
print(f"\ntrue variant loci={len(varloci)}  injected errors={len(truthE):,}")

# ---- depth-stratified error recall: separate method-limit from coverage-limit ----
# a truth error is 'callable' if that read-position sits at depth>=3 in the pileup
pos_depth={}   # (rid,pos) -> depth  (from pileup rows; need rid,pos,depth)
for line in open('sim_pileup.tsv'):
    if line[0]=='r': continue
    v=line.split()
    rid=int(v[0]); depth=int(v[7]); pos=int(v[9])
    pos_depth[(rid,pos)]=depth
va=call_variants(30)
callE=set()
for rid,pos,col,rb in mm_rows:
    if not (col in va and rb in va[col]): callE.add((rid,pos))
callable_err=[e for e in truthE if pos_depth.get(e,1)>=3]
d1_err=[e for e in truthE if pos_depth.get(e,1)<3]
rec_callable=len(callE & set(callable_err))/max(len(callable_err),1)
print(f"\n[QUAL>=30] error recall breakdown:")
print(f"  total injected errors: {len(truthE):,}")
print(f"  at callable depth>=3:  {len(callable_err):,}  -> recall {rec_callable:.3f}  (METHOD ceiling)")
print(f"  at depth<3 (invisible): {len(d1_err):,}  ({100*len(d1_err)/len(truthE):.1f}% of errors, un-detectable by any pileup)")
