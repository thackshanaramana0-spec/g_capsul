import re,sys,os
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
    # Match whichever haplotype the genome carries; that allele becomes REF, the
    # other ALT.
    #
    # POLARITY FIX (2026-09-03). This loop used to take the FIRST haplotype that
    # matched anywhere in the +/-300bp window, which silently made hapR win by
    # loop order. Inside a tandem repeat BOTH haplotypes match -- shifted copies
    # of the repeat unit exist on either side -- so the winner was decided by
    # iteration order rather than by evidence, and at repeats it lands on the
    # wrong one. Measured on the tetraploid benchmark: truth TTTTA->T (deletion)
    # came out as T->Tttta (insertion) at 20:3332481, TTTTATTTA->T as
    # T->TTTTATTTA at 20:3346020, GCA->G as G->GCA at 20:3097933 -- the same two
    # haplotypes each time, with REF and ALT swapped.
    #
    # The alignment already knows where this locus belongs: c2g() returns the
    # genome coordinate of the contig anchor, so a TRUE haplotype match starts
    # at about (gpos-1-K) while a repeat-shifted spurious match starts a whole
    # number of repeat units away. Collect every candidate match and keep the
    # one closest to that expected start. Ties keep the previous hapR-first
    # order, so any locus where only one haplotype matches is unchanged.
    _cands=[]
    for hap,mine,other in ((hapR,cref,calt),(hapO,calt,cref)):
        for ori in (hap, rc(hap)):
            st=0
            while True:
                k=gw.find(ori,st)
                if k<0: break
                _cands.append((abs((lo+k)-(gpos-1-K)), lo+k, hap, ori, mine, other))
                st=k+1
    _cands.sort(key=lambda c:c[0])
    # AMBIGUITY GUARD (2026-09-03). Inside a tandem repeat BOTH haplotypes match
    # the genome window, so a sequence search cannot decide which one the genome
    # actually carries -- and picking either way is a coin flip on insertion vs
    # deletion polarity. When both hapR and hapO match, refuse to decide here and
    # fall through to the CIGAR-based path below, which has real alignment
    # evidence (a D operation in the reference contig's own alignment means the
    # genome carries those bases, i.e. the event is a deletion). Measured on the
    # tetraploid benchmark: contig_97 aligns 7M1I433M8D294M -- an explicit 8bp
    # deletion -- at the (TTTA)n locus 20:3346020 where truth is TTTTATTTA->T and
    # the sequence-search path emitted the inverse T->TTTTATTTA.
    if len({id(c[2]) for c in _cands}) > 1:
        _hapset = {('R' if c[2] is hapR else 'O') for c in _cands}
        if len(_hapset) > 1:
            return None
    for _d,_gs,hap,ori,mine,other in _cands:
        if True:
            gstart=_gs                                   # 0-based genome start of matched haplotype
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

_feat={}                      # eval-only: genome pos -> caller INFO
rows=[]
indel_rows=[]
indel_rows_rev=0
# EVAL-ONLY, OPT-IN (LIFT_KEEP_INFO=1). The lift hardcodes INFO to "." and
# deduplicates output keys, which is correct for scoring a caller's VCF but
# destroys two things a CANDIDATE channel needs: the per-record evidence fields
# (DP/AF/MLEN/MMCNT) and the RECURRENCE -- how many source records mapped to the
# same genome locus. For mem_extmm that recurrence IS the evidence: one record
# is a sequencing error, many independent pg regions disagreeing at one locus is
# a real variant. Counted here, emitted as RC=. Default path is byte-unchanged.
_KEEPINFO = bool(os.environ.get('LIFT_KEEP_INFO'))
_rc = {}
_inf = {}
def _note(gpos, gR, alt, info):
    if not _KEEPINFO: return
    k = (gpos, len(gR), len(alt))
    _rc[k] = _rc.get(k, 0) + 1
    if k not in _inf: _inf[k] = info if info else '.'
for line in open(CALLS):
    if line[0]=='#':continue
    f=line.rstrip('\n').split('\t')
    rn=f[0];cp1=int(f[1]);cref=f[3].upper();calt=f[4].upper()
    # Support evidence the caller already computed (SVTYPE/AF/DP/ANCHORS). Used
    # ONLY to arbitrate between competing indel hypotheses at the same genome
    # locus (see the commit step at the end); never to filter a call out.
    _info=f[7] if len(f)>7 else ''
    def _iv(k, d=0.0):
        m=re.search(k+r'=([0-9.]+)', _info)
        return float(m.group(1)) if m else d
    _support=(_iv('ANCHORS'), _iv('DP'), _iv('AF'))
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
            _note(gpos,gR,",".join(alts),_info)
            _feat[gpos]=_info
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
            indel_rows.append((hf[0],hf[1],hf[2],"0/1",_support)); continue

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
        else:
            # INSERTION-TYPE, CIGAR EVIDENCE (2026-09-03). Previously the CIGAR
            # was consulted ONLY when the caller had already guessed deletion,
            # so an insertion-labelled call never got the benefit of the
            # alignment -- even though the caller's ref/alt polarity is
            # arbitrary (it labels the longer contig "reference"). If the
            # reference contig's own alignment carries a D/N operation of the
            # same length within a repeat-unit's reach of the anchor, the GENOME
            # holds those bases and the event is really a deletion.
            if rn in cmap:
                _rs,_cig,_rev,_clen = cmap[rn]
                _t = (_clen-1-((cp1-1))) if _rev else (cp1-1)
                _ri, _qi = 0, 0
                for _n,_op in re.findall(r'(\d+)([MIDNSHP=X])', _cig):
                    _n=int(_n)
                    if _op in 'M=X': _ri+=_n; _qi+=_n
                    elif _op in 'IS': _qi+=_n
                    elif _op in 'DN':
                        if _n==d and abs(_qi-_t) <= max(2*d, 12):
                            is_del=True; break
                        _ri+=_n
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
            indel_rows.append((apos, gA+fwd, gA, "0/1", _support))
        else:
            apos = gpos if not rev else gpos-1            # left anchor for the insertion
            if apos<1 or apos>len(refseq): continue
            gA=refseq[apos-1].upper()
            if gA=='N': continue
            indel_rows.append((apos, gA, gA+fwd, "0/1", _support))
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
    _note(gpos,gR,alt,_info)
    _feat[gpos]=_info          # eval-only: carry caller features for TP/FP analysis

# COMMIT TO ONE INDEL CALL PER LOCUS.
# At an ambiguous repeat the extractor can propose several mutually exclusive
# indels at the same position (measured on the tetraploid benchmark: 86 indel
# records over only 81 distinct positions, including three contradictory calls
# at one (CA)n locus -- an insertion of CA, an insertion of CACACACA, and a
# deletion of CACACA). A caller has to commit: at most one of those can be
# true, and every extra one is a false positive by construction. DiscoSNP++
# emits exactly one record per position (60 calls / 60 positions) and that is
# part of why its indel precision is higher.
# Arbitration uses the evidence the caller ALREADY produced -- junction anchor
# support first, then depth, then allele fraction -- and is applied blind to
# the truth set. Nothing is filtered out: a locus with one hypothesis keeps it
# unchanged, so this can only remove self-contradictions.
def _leftalign(pos, ref, alt):
    # Standard VCF left-alignment against the reference: while the alleles end
    # in the same base and both are >1, trim that base and shift left. Two
    # different representations of the same repeat-locus event converge to the
    # same (pos,ref,alt) after this, which is what makes per-locus arbitration
    # below actually collapse them. Without it bcftools norm does the shift
    # LATER, at scoring time, and re-creates the duplicates this step removes
    # (observed at the (CA)n locus 20:3097933).
    guard=0
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1] and guard<200:
        ref=ref[:-1]; alt=alt[:-1]; guard+=1
    while ref and alt and ref[-1]==alt[-1] and (len(ref)==1 or len(alt)==1) and pos>1 and guard<200:
        b=refseq[pos-2].upper()
        if not b or b=='N': break
        ref=b+ref[:-1]; alt=b+alt[:-1]; pos-=1; guard+=1
    return pos, ref, alt

# REFUTED BY MEASUREMENT 2026-09-03 -- OPT-IN (CAPS_ONE_INDEL_PER_LOCUS=1),
# OFF by default. The reasoning is sound (a caller should not emit mutually
# exclusive indels at one locus; DiscoSNP++ emits exactly 1 record per
# position, 60/60, and CAPSULE emitted 86 over 81 positions) and it does
# improve PRECISION (tetraploid FP 13 -> 10/11). But it costs more true
# positives than it removes false ones: tetraploid indel F1 0.555 -> 0.548.
# The multi-hypothesis emission was net-positive because one of the competing
# hypotheses was often right. Kept on the record, not shipped.
if os.environ.get('CAPS_ONE_INDEL_PER_LOCUS'):
    indel_rows=[( *_leftalign(r[0],r[1],r[2]), r[3], r[4] if len(r)>4 else (0.0,0.0,0.0)) for r in indel_rows]
    _best={}
    for r in indel_rows:
        pos=r[0]; sup=r[4] if len(r)>4 else (0.0,0.0,0.0)
        if pos not in _best or sup > _best[pos][1]:
            _best[pos]=(r,sup)
    indel_rows=[v[0][:4] for v in _best.values()]
else:
    indel_rows=[r[:4] for r in indel_rows]

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
        if _KEEPINFO:
            _ii = _inf.get(key,'.')
            o.write("%s\t%d\t.\t%s\t%s\t30\tPASS\t%s;RC=%d\tGT\t%s\n"
                    %(CHROM,gpos,gR,alt,_ii,_rc.get(key,1),gt))
        else:
            o.write("%s\t%d\t.\t%s\t%s\t30\tPASS\t.\tGT\t%s\n"%(CHROM,gpos,gR,alt,gt))
import os as _os
with open(_os.path.splitext(OUT)[0]+".feat.tsv","w") as _ff:
    for _p,_i in sorted(_feat.items()): _ff.write("%d\t%s\n"%(_p,_i))
print("lifted %d calls (%d SNV, %d indel; both strands) -> %s"
      %(len(seen),len(rows),len(indel_rows),OUT))
