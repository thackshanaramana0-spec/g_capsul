import re,sys
# Lift an ARCS contig-coordinate call VCF to genome coordinates (allele-aware) so a
# standard tool (rtg vcfeval / hap.py) can score it. Uses the contig->reference
# alignment (C2R, eval-only) to place each call, fetches the genome REF base from the
# reference FASTA, and complements alleles for reverse-strand contigs. Emits a valid
# genome-coordinate VCF (GT=0/1, or 1/2 when neither observed allele is the ref).
#   usage: lift_vcf.py <calls.contig.vcf> <c2r.sam> <ref.fa> <chrom> <out.vcf> [contigs.fa]
# The optional contigs.fa (the caller's ARCS_DUMP_CONTIGS sequences, ">contig_N") enables
# the exact haplotype-flank indel lift: the genome equals exactly ONE of the two bubble
# haplotypes, so matching each allele-with-flanks against the genome resolves both the
# insertion/deletion polarity and the strand uniformly, repeat-robustly.
CALLS,C2R,REF,CHROM,OUT=sys.argv[1:6]
CONTIGS=sys.argv[6] if len(sys.argv)>6 else None
comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}

# reference sequence for CHROM (1-based indexing via seq[pos-1])
seq=[]
cur=None
for line in open(REF):
    if line[0]=='>':
        cur=line[1:].split()[0]
        continue
    if cur==CHROM: seq.append(line.strip())
refseq=''.join(seq)
def gref(pos1): return refseq[pos1-1].upper() if 1<=pos1<=len(refseq) else 'N'

# contig -> (genome_start0, cigar, reverse, contig_len)
cmap={}
for line in open(C2R):
    if line[0]=='@':continue
    f=line.split('\t');flag=int(f[1]);mapq=int(f[4])
    if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20:continue
    if f[0] in cmap:continue
    cmap[f[0].split("|")[0]]=(int(f[3])-1,f[5],bool(flag&0x10),len(f[9]))
def c2g(rn,cp0):
    if rn not in cmap:return None
    rs,cig,rev,clen=cmap[rn];target=(clen-1-cp0) if rev else cp0
    refidx=rs;readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            if readidx<=target<readidx+n: return (refidx+(target-readidx)+1, rev)  # 1-based genome pos
            refidx+=n;readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
    return None

def cig_op_at(rn,cp0):
    # Return the CIGAR operation ('M'/'I'/'D'/'S') that covers original-contig position
    # cp0 in the reference contig's alignment to the genome. This tells us whether a body
    # base is ALIGNED to the genome (M -> the genome carries it -> deletion in the sample)
    # or INSERTED relative to the genome (I -> genome lacks it -> insertion) — robustly,
    # from the actual alignment, not a repeat-fooled sequence search.
    if rn not in cmap:return None
    rs,cig,rev,clen=cmap[rn];target=(clen-1-cp0) if rev else cp0
    readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            if readidx<=target<readidx+n:return 'M'
            readidx+=n
        elif op in 'IS':
            if readidx<=target<readidx+n:return op
            readidx+=n
        # D/N consume reference only, not the query index
    return None

def rc(s):
    return ''.join(comp.get(c,'N') for c in reversed(s))

# contig sequences (for the exact haplotype-flank indel lift)
contigseq={}
if CONTIGS:
    cur=None
    for line in open(CONTIGS):
        line=line.rstrip('\n')
        if not line: continue
        if line[0]=='>': cur=line[1:].split()[0]; contigseq[cur]=[]
        elif cur is not None: contigseq[cur].append(line)
    contigseq={k:''.join(v) for k,v in contigseq.items()}

def normalize_indel(pos,ref,alt):
    # left-trim shared suffix then shared prefix (keep 1 anchor base) so representation is
    # canonical; rtg still haplotype-matches, but this keeps our own records tidy/valid.
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]:
        ref=ref[:-1]; alt=alt[:-1]
    return pos,ref,alt

def hapflank_lift(rn,cp1,cref,calt):
    # Exact lift using contig flanks: build each bubble haplotype (leftflank + allele +
    # rightflank) from the reference contig, find which one the genome carries (in either
    # orientation), and read the genome-frame REF/ALT + position straight off that match.
    # Returns (pos,ref,alt) or None. cref is always the reference-contig's allele.
    if rn not in contigseq: return None
    R=contigseq[rn]; a0=cp1-1; K=25
    if a0<K or a0+len(cref)+K>len(R): return None
    if R[a0:a0+1]!=cref[0]: return None                  # sanity: contig carries cref here
    lf=R[a0-K:a0]; rf=R[a0+len(cref):a0+len(cref)+K]
    hapR=lf+cref+rf; hapO=lf+calt+rf                     # the two full local haplotypes
    gp0=c2g(rn,cp1-1)
    if gp0 is None: return None
    gpos=gp0[0]
    lo=max(0,gpos-300); hi=min(len(refseq),gpos+300); gw=refseq[lo:hi].upper()
    # match whichever haplotype the genome carries; that allele becomes REF, the other ALT
    for hap,mine,other in ((hapR,cref,calt),(hapO,calt,cref)):
        for ori in (hap, rc(hap)):
            k=gw.find(ori)
            if k<0: continue
            gstart=lo+k                                  # 0-based genome start of matched haplotype
            if ori==hap:                                 # forward strand: allele sits after lf
                apos=gstart+K                            # 0-based anchor position
                ref_al=mine; alt_al=other                # genome frame == contig frame
                p,ref_al,alt_al=normalize_indel(apos+1,ref_al,alt_al)
            else:                                        # reverse strand: window holds rc(hap)
                # rc(hap)=rc(rf)+rc(mine)+rc(lf); the rc(mine) block starts at gstart+K
                anchor0=gstart+K-1                        # genome base just left of the block
                if anchor0<0: return None
                ganch=refseq[anchor0].upper()
                # rc(mine)/rc(other) are right-anchored; re-anchor on ganch (drop shared tail)
                ref_al=ganch+rc(mine)[:-1]; alt_al=ganch+rc(other)[:-1]
                p,ref_al,alt_al=normalize_indel(anchor0+1,ref_al,alt_al)
            if not ref_al or not alt_al: return None
            if refseq[p-1:p-1+len(ref_al)].upper()!=ref_al.upper(): return None
            return (p,ref_al,alt_al)
    return None

rows=[]
indel_rows=[]
indel_rows_rev=0
for line in open(CALLS):
    if line[0]=='#':continue
    f=line.rstrip('\n').split('\t')
    rn=f[0];cp1=int(f[1]);cref=f[3].upper();calt=f[4].upper()
    # multi-allelic SNV (polyploid): ALT="C,G" — REF single base, every ALT single base.
    # Lift each allele independently (same contig->genome mapping) and emit one record with
    # the comma-joined ALTs in genome frame (complemented for reverse-strand contigs).
    if len(cref)==1 and ',' in calt and all(len(a)==1 for a in calt.split(',')):
        r=c2g(rn,cp1-1)
        if r is None: continue
        gpos,rev=r; gR=gref(gpos)
        if gR=='N': continue
        # observed alleles = the contig CONSENSUS base (cref) PLUS the ALT alleles — cref is
        # itself one of the real haplotype alleles, and at a multi-allelic site the consensus
        # may be a non-reference allele, so it must be counted or that allele is lost.
        obsc=[cref]+calt.split(',')
        obs=[(comp.get(a,'N') if rev else a) for a in obsc]
        alts=sorted({a for a in obs if a!=gR and a!='N'})
        if alts:
            gt="/".join(str(i+1) for i in range(len(alts))) if len(alts)>1 else "0/1"
            rows.append((gpos,gR,",".join(alts),gt))
        continue
    if len(cref)!=1 or len(calt)!=1:
        # ── indel: the caller labels the LONGER haplotype contig "reference", so its
        # ref/alt polarity is arbitrary; the contig may also map reverse-strand. The GENOME
        # is authoritative. `body` = the differing bases (longer allele minus its anchor),
        # oriented to the forward genome (reverse-complement for reverse-strand contigs).
        # c2g already returns the exact genome coordinate of ANY contig position on either
        # strand, so we test the body against the genome at that EXACT adjacent position
        # (no fuzzy window — short bodies match spuriously in a window and misplace calls):
        #   forward: body sits just AFTER the anchor;  reverse: just BEFORE it (rc'd).
        # If the genome carries body there the alt DELETED it, else the alt INSERTED it.
        # Primary: exact haplotype-flank lift (needs contigs.fa) — resolves polarity+strand.
        hf=hapflank_lift(rn,cp1,cref,calt)
        if hf is not None:
            indel_rows.append((hf[0],hf[1],hf[2],"0/1")); continue

        gp0=c2g(rn,cp1-1)                                 # genome pos of the contig anchor
        if gp0 is None: continue
        gpos,rev=gp0
        del_type = len(cref)>len(calt)                    # caller expressed a DELETION (body on THIS contig)
        L = cref if len(cref)>=len(calt) else calt        # anchor + body (contig frame)
        body = L[1:]
        if body=='' or 'N' in body: continue
        d=len(body)
        fwd = rc(body) if rev else body                   # forward-genome orientation

        # ── Determine polarity vs the genome. When the caller emitted a DELETION, the
        # body bases live on THIS (reference) contig, so the bwa CIGAR tells us directly
        # whether the genome carries them (M -> deletion) or they are inserted (I ->
        # insertion) — repeat-robust. When the caller emitted an INSERTION, the body is on
        # the OTHER contig; fall back to a genome-adjacency probe. ──
        is_del=None
        if del_type:
            body0 = (cp1-1)+1                             # 0-based position of first body base on contig
            op = cig_op_at(rn, body0)
            if op=='M':   is_del=True                     # body aligns to genome
            elif op in ('I','S'): is_del=False            # body inserted vs genome
        if is_del is None:                                # insertion-type, or CIGAR indeterminate
            if not rev: is_del = (refseq[gpos:gpos+d].upper()==fwd)
            else:       is_del = (refseq[gpos-d:gpos].upper()==fwd)

        if is_del:
            # deleted bases sit AFTER anchor (fwd strand) or BEFORE it (rev strand)
            if not rev:
                apos=gpos
            else:
                apos=gpos-d-1                              # left anchor before the deleted run
            if apos<1 or apos>len(refseq): continue
            gA=refseq[apos-1].upper()
            if gA=='N': continue
            # genome-consistency guard: REF body must actually be in the genome here
            if refseq[apos:apos+d].upper()!=fwd:
                # try a small left-align search so repeats still classify correctly
                hit=None
                for p in range(max(1,apos-d-2),apos+d+3):
                    if refseq[p:p+d].upper()==fwd and refseq[p-1].upper()!='N': hit=p;break
                if hit is None: continue
                apos=hit; gA=refseq[apos-1].upper()
            indel_rows.append((apos, gA+fwd, gA, "0/1"))
        else:
            apos = gpos if not rev else gpos-1            # left anchor for the insertion
            if apos<1 or apos>len(refseq): continue
            gA=refseq[apos-1].upper()
            if gA=='N': continue
            indel_rows.append((apos, gA, gA+fwd, "0/1"))
        continue
    r=c2g(rn,cp1-1)
    if r is None: continue
    gpos,rev=r
    gR=gref(gpos)
    if gR=='N': continue
    # observed alleles in genome frame (major=cref, minor=calt), complement if reverse
    a1 = comp.get(cref,'N') if rev else cref
    a2 = comp.get(calt,'N') if rev else calt
    if a1=='N' or a2=='N': continue
    # standard VCF representation: REF=gR, ALT=non-ref observed allele(s)
    if   a1==gR and a2!=gR: alt=a2; gt="0/1"
    elif a2==gR and a1!=gR: alt=a1; gt="0/1"
    elif a1!=gR and a2!=gR:
        # neither observed allele is the reference -> both are ALT
        if a1==a2: alt=a1; gt="1/1"
        else:      alt=a1+","+a2; gt="1/2"
    else:
        continue                                        # both == ref, not a variant
    rows.append((gpos,gR,alt,gt))

allrows=rows+indel_rows
allrows.sort(key=lambda r:(r[0],len(r[1]),len(r[2])))
with open(OUT,'w') as o:
    o.write("##fileformat=VCFv4.2\n")
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write("##contig=<ID=%s>\n"%CHROM)
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    seen=set()
    for gpos,gR,alt,gt in allrows:
        key=(gpos,len(gR),len(alt))                     # allow a SNV and an indel to coexist
        if key in seen: continue
        seen.add(key)
        o.write("%s\t%d\t.\t%s\t%s\t30\tPASS\t.\tGT\t%s\n"%(CHROM,gpos,gR,alt,gt))
print("lifted %d calls (%d SNV, %d indel; both strands) -> %s"
      %(len(seen),len(rows),len(indel_rows),OUT))
