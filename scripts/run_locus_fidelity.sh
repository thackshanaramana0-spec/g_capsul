#!/usr/bin/env bash
# T3.4 — LOCUS RETRIEVAL FIDELITY: does a locus query return the evidence?
#
#   usage: run_locus_fidelity.sh <capsule_decode> <archive> <ref.fa> [IND] [outdir] [N]
#
# T3.1/T3.2/T3.3 measure what retrieval COSTS. This measures whether what comes
# back is USABLE, which turned out to be the more important question.
#
# WHY THIS TABLE EXISTS. A locus query was returning reads that agreed with each
# other perfectly -- 121 overlapping pairs, zero mismatching bases -- which is
# impossible for real reads at a heterozygous site. The cause is structural: in
# a pseudogenome built to minimise bits, a het locus is not one place but N
# parallel places, because two internally-consistent haplotype contigs compress
# better than one contig plus a column of disagreements, and each is stored
# again reverse-complemented. Measured: median 4 parallel loci per site, with
# the two haplotypes of GIAB 20:3001343 sitting 18.6 Mb apart in pg space.
#
# So a COORDINATE names one of those places and returns one haplotype, with the
# variation gone. Content addressing resolves every parallel representative at
# once and returns their union, which is a real pileup. This script measures
# both on the same archive and the same sites, so the comparison is internal
# and cannot be attributed to anything else.
#
# The phenomenon itself is known in the DE NOVO ASSEMBLY literature (Purge
# Haplotigs 2018, Redundans 2016). What is measured here is its consequence in
# a compression archive, where it is an addressability defect rather than a
# quality defect -- and where the assembly field's remedy (purge the redundant
# haplotigs) is inadmissible, because deleting a haplotig deletes an allele,
# a variant, and losslessness.
set -u
DEC="${1:?usage: run_locus_fidelity.sh <capsule_decode> <archive> <ref.fa> [IND] [outdir] [N]}"
ARC="${2:?}"; REF="${3:?}"
IND="${4:-HG002}"; OUT="${5:-/tmp/t34_$IND}"; N="${6:-100}"
TRUTH="$HOME/giab_truth/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
CHROM=20; LO=3000000; HI=3600000          # a window with dense GIAB truth
PAR="${T34_PAR:-6}"                        # correctness only: parallel is safe

# Resolve to absolute paths BEFORE the cd below. The script cd's into its
# output directory, so a relative archive/ref/decoder path would silently stop
# resolving there -- which failed as "export failed" and looked like a decoder
# problem rather than a path one.
abspath(){ case "$1" in /*) printf '%s' "$1";; *) printf '%s/%s' "$(pwd)" "$1";; esac; }
DEC="$(command -v "$DEC" 2>/dev/null || abspath "$DEC")"
ARC="$(abspath "$ARC")"; REF="$(abspath "$REF")"
case "$OUT" in /*) ;; *) OUT="$(abspath "$OUT")";; esac

for t in tabix samtools python3; do command -v $t >/dev/null || { echo "MISSING: $t" >&2; exit 1; }; done
[ -x "$DEC" ] || { echo "MISSING/not executable: $DEC" >&2; exit 1; }
[ -s "$REF" ] || { echo "MISSING ref: $REF" >&2; exit 1; }
[ -s "$ARC" ]   || { echo "MISSING archive: $ARC" >&2; exit 1; }
[ -s "$TRUTH" ] || { echo "MISSING truth: $TRUTH" >&2; exit 1; }
mkdir -p "$OUT/q"; cd "$OUT"

echo "=== T3.4 locus retrieval fidelity — $IND chr$CHROM:$LO-$HI ==="
[ -s pg.fa ] || "$DEC" export "$ARC" pg.fa >/dev/null 2>&1
[ -s pg.fa ] || { echo "export failed: $DEC export $ARC" >&2; exit 1; }
[ -s ref.txt ] || samtools faidx "$REF" $CHROM:$((LO-1000))-$((HI+1000)) 2>/dev/null \
    | tail -n +2 | tr -d '\n' | tr 'acgt' 'ACGT' > ref.txt
tabix "$TRUTH" $CHROM:$LO-$HI 2>/dev/null \
  | awk -F'\t' '$10~/^0[|\/]1|^1[|\/]0/ && length($4)==1 && length($5)==1 {print $2"\t"$4"\t"$5}' > hets.tsv
echo "  GIAB het SNVs in window: $(wc -l < hets.tsv)"

# Build the probe list. A probe is 40 bp of REFERENCE ending 6 bp BEFORE the
# variant, so it never contains the variant itself -- otherwise it could only
# ever match the haplotype it was taken from, which would rig the result.
python3 - "$N" "$((LO-1000))" > probes.tsv <<'PY'
import sys
N=int(sys.argv[1]); REF0=int(sys.argv[2])
ref=open("ref.txt").read().strip()
pg="".join(l.strip() for l in open("pg.fa") if l[0]!='>')
n=0
for line in open("hets.tsv"):
    pos,r,a=line.split(); pos=int(pos)
    lo=pos-45-REF0; hi=pos-5-REF0
    if lo<0 or hi>len(ref): continue
    p=ref[lo:hi]
    if len(p)!=40 or set(p)-set("ACGT"): continue
    at=pg.find(p)
    if at<0: continue                      # need a coordinate control as well
    print(f"{pos}\t{r}\t{a}\t{p}\t{at+45}")
    n+=1
    if n>=N: break
PY
echo "  probes usable: $(wc -l < probes.tsv)"

# Both query modes, same archive, same sites. Parallel is legitimate here: this
# measures WHAT comes back, not how fast, so no timing is recorded.
cat probes.tsv | xargs -P "$PAR" -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r POS R A P V <<< "{}"
  '"$DEC"' query '"$ARC"' q/s_$POS.fa "$P"                   >/dev/null 2>&1
  '"$DEC"' query '"$ARC"' q/c_$POS.fa $((V-100))-$((V+100))  >/dev/null 2>&1'

python3 - "$IND" <<'PY'
import sys
IND=sys.argv[1]
comp={'A':'T','C':'G','G':'C','T':'A'}
def rc(s): return "".join(comp.get(c,'N') for c in reversed(s))
sb=so=sn=cb=co=cn=0; SR=SA=0; n=0
for line in open("probes.tsv"):
    pos,r,a,probe,v=line.rstrip("\n").split("\t"); v=int(v); n+=1
    seqs=[]
    try:
        for l in open(f"q/s_{pos}.fa"):
            if l[0]!='>': seqs.append(l.strip())
    except FileNotFoundError: pass
    rcp=rc(probe); nr=na=0
    for s in seqs:
        i=s.find(probe)
        if i>=0 and i+45<len(s): b=s[i+45]
        else:
            j=s.find(rcp)
            if j>=0 and j-6>=0: b=comp.get(s[j-6],'N')
            else: continue
        if b==r: nr+=1
        elif b==a: na+=1
    SR+=nr; SA+=na
    sb+= 1 if (nr and na) else 0; so+= 1 if ((nr or na) and not (nr and na)) else 0
    sn+= 1 if not (nr or na) else 0
    cr=ca=0; p=None
    try:
        for l in open(f"q/c_{pos}.fa"):
            if l[0]=='>': p=int(l.split("pos=")[1].split()[0])
            else:
                s=l.strip(); off=v-p
                if 0<=off<len(s):
                    b=s[off]
                    if b==r: cr+=1
                    elif b==a: ca+=1
    except FileNotFoundError: pass
    cb+= 1 if (cr and ca) else 0; co+= 1 if ((cr or ca) and not (cr and ca)) else 0
    cn+= 1 if not (cr or ca) else 0
print()
print(f"======== T3.4 LOCUS RETRIEVAL FIDELITY — {IND} ========")
print(f"{'addressing mode':24}{'both alleles':>14}{'one allele':>12}{'neither':>9}")
print(f"{'coordinate':24}{cb:>14}{co:>12}{cn:>9}")
print(f"{'content (sequence)':24}{sb:>14}{so:>12}{sn:>9}")
print(f"sites tested: {n}   content recovers both at {100*sb/max(n,1):.1f}%, coordinate at {100*cb/max(n,1):.1f}%")
print(f"allele balance across all returned evidence: {SR} REF / {SA} ALT  (ratio {SA/max(SR,1):.2f})")
print("=======================================================")
print(f"T34,{IND},{n},{cb},{co},{cn},{sb},{so},{sn},{SR},{SA}")
PY
