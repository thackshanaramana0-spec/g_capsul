#!/usr/bin/env python3
# Simulate a diploid genome with known het SNVs + het indels, emit ~cov x reads.
# Used to smoke-test arcs indel calling: the truth (het del + het ins) is printed so
# we can check `arcs call` recovers it. Deterministic.
import random, sys
random.seed(7)
L = 8000
cov = 80
rl = 150
bases = "ACGT"
hapA = [random.choice(bases) for _ in range(L)]

# hapB = hapA with a few het SNVs + one het deletion + one het insertion
hapB = list(hapA)
snv_sites = [1200, 2600, 4100, 5500, 6900]
for p in snv_sites:
    o = [b for b in bases if b != hapB[p]]
    hapB[p] = random.choice(o)

# het deletion of 3 bp at 3000 (hapB loses bases 3000..3002)
del_pos, del_len = 3000, 3
del_ref = ''.join(hapA[del_pos-1:del_pos+del_len])   # anchor + deleted (VCF-style)
hapB_str = ''.join(hapB[:del_pos]) + ''.join(hapB[del_pos+del_len:])

# het insertion of "GAT" at 5000 (in the already-deleted hapB coordinate frame)
ins_pos, ins_seq = 5000, "GAT"
hapB_str = hapB_str[:ins_pos] + ins_seq + hapB_str[ins_pos:]

hapA_str = ''.join(hapA)

def emit(fh, hap, n, tag):
    for i in range(n):
        s = random.randint(0, len(hap) - rl)
        r = hap[s:s+rl]
        if random.random() < 0.5:
            comp = {'A':'T','C':'G','G':'C','T':'A'}
            r = ''.join(comp[c] for c in reversed(r))
        fh.write(f"@{tag}_{i}\n{r}\n+\n{'I'*rl}\n")

with open(sys.argv[1], "w") as fh:
    nreads = (L * cov) // rl
    emit(fh, hapA_str, nreads // 2, "A")
    emit(fh, hapB_str, nreads // 2, "B")

print(f"TRUTH het SNVs (hapA coord, 0-based): {snv_sites}")
print(f"TRUTH het DEL: pos0={del_pos} len={del_len} ref={del_ref} alt={hapA[del_pos-1]}")
print(f"TRUTH het INS: near hapA coord ~{ins_pos+del_len} seq={ins_seq}")

# reference = hapA; write FASTA + truth VCF (genome coords, 1-based) for the lift test
if len(sys.argv) > 2:
    with open(sys.argv[2], "w") as fa:
        fa.write(">chrSIM\n")
        for i in range(0, len(hapA_str), 70): fa.write(hapA_str[i:i+70] + "\n")
    with open(sys.argv[3], "w") as tv:
        tv.write("##fileformat=VCFv4.2\n##contig=<ID=chrSIM>\n")
        tv.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        tv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        recs=[]
        for p in snv_sites:
            recs.append((p+1, hapA[p], hapB[p], "0/1"))
        # deletion: hapB lost bases del_pos..del_pos+del_len-1 (0-based) -> anchor at del_pos-1
        recs.append((del_pos, del_ref, hapA[del_pos-1], "0/1"))
        # insertion: hapB gained ins_seq before hapA coord del_pos... after del, at ins_pos in
        # hapB frame == genome(hapA) coord ins_pos+del_len. anchor = hapA[that-1].
        gp = ins_pos + del_len
        recs.append((gp, hapA_str[gp-1], hapA_str[gp-1] + ins_seq, "0/1"))
        for pos, r, a, gt in sorted(recs):
            tv.write(f"chrSIM\t{pos}\t.\t{r}\t{a}\t30\tPASS\t.\tGT\t{gt}\n")
    print(f"wrote reference {sys.argv[2]} + truth {sys.argv[3]}")
