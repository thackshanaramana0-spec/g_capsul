#!/usr/bin/env python3
# Real, full-pipeline decoder for the locked sequence+order scope.
# Reconstructs every original read's sequence from all layers and writes
# them out in original file order, for direct diff against the true
# original FASTQ's sequence column. This is the actual test this project's
# "round trip: VERIFIED" claims never did until this was built.
#
# STAGE 100 REWRITE: matches the simplified encoder that dropped the whole
# rank/bucket/uidOrder scheme in favor of PgRC2's real design (confirmed
# from source, SeparatedPseudoGenomePersistence.cpp:446, singleFileMode):
# position/strand/length are a flat array, indexed DIRECTLY by original
# read index -- no rank layer, no bucket-of-duplicates, no separate order
# stream at all. orig2uid.bin (original -> unique-id) survives only to
# correlate a read with its mismatch data, since mismatches are computed
# once per unique sequence and shared by any duplicate original reads.
import sys, struct

WORK = sys.argv[1]           # directory containing the dumped files
PG_LEN = int(sys.argv[2])
MAIN_PG_END = int(sys.argv[3])
OUT = sys.argv[4]

COMP = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
COMP_B = {ord(k):ord(v) for k,v in COMP.items()}
def revcomp(s):
    return ''.join(COMP[c] for c in reversed(s))

def read_varints(path):
    with open(path,'rb') as f:
        data = f.read()
    vals = []
    i = 0
    while i < len(data):
        v = 0; shift = 0
        while True:
            b = data[i]; i += 1
            v |= (b & 0x7f) << shift
            if not (b & 0x80): break
            shift += 7
        vals.append(v)
    return vals

def read_u32(path):
    with open(path,'rb') as f: d = f.read()
    return list(struct.unpack(f'<{len(d)//4}I', d))

def read_u16(path):
    with open(path,'rb') as f: d = f.read()
    return list(struct.unpack(f'<{len(d)//2}H', d))

def read_bits(path, n):
    with open(path,'rb') as f: d = f.read()
    bits = []
    for byte in d:
        for k in range(7,-1,-1):
            bits.append((byte>>k)&1)
    return bits[:n]

# ---- layer 1: literal.txt ----
with open(f'{WORK}/literal.txt','rb') as f:
    literal = f.read()

# ---- layer 3: mem_triples.bin ----
# STAGE 100 REWRITE: the encoder now does a real trim-on-overlap pass
# (mirroring SimplePgMatcher.cpp:99-131) before ever writing literal.txt or
# these triples, so the two are GUARANTEED consistent by construction --
# no overlaps, no duplicates, no reconstruction-from-a-lossy-dump needed.
# Also now carries a real is_rc flag (4th field) per match -- previously
# absent, which was the root cause of ~22% of matches failing a direct-copy
# consistency check (they were RC-discovered and need reverse-complementing,
# not a plain copy).
with open(f'{WORK}/mem_triples.bin','rb') as f: raw = f.read()
n_triples = len(raw)//13   # dst(4)+src(4)+len(4)+is_rc(1), no padding, read as raw bytes
triples = []
off = 0
for _ in range(n_triples):
    dst, src, ln = struct.unpack_from('<III', raw, off); off += 12
    is_rc = raw[off]; off += 1
    triples.append((dst, src, ln, is_rc))

pg = bytearray(PG_LEN)
lit_i = 0
pos = 0
for dst, src, ln, is_rc in triples:
    # gap before this match: literal
    pg[pos:dst] = literal[lit_i:lit_i+(dst-pos)]
    lit_i += dst - pos
    if is_rc:
        seg = bytes(COMP_B[b] for b in reversed(pg[src:src+ln]))
        pg[dst:dst+ln] = seg
    else:
        pg[dst:dst+ln] = pg[src:src+ln]
    pos = dst + ln
pg[pos:PG_LEN] = literal[lit_i:lit_i+(PG_LEN-pos)]
lit_i += PG_LEN - pos
assert lit_i == len(literal), f"literal not fully consumed: {lit_i} vs {len(literal)}"
pg = pg.decode('ascii')
print(f"[decode] pg reconstructed: {len(pg)} bytes, {n_triples} clean refs, literal fully consumed ({lit_i})", file=sys.stderr)

# ---- position/strand/length: flat, indexed by UNIQUE id (0..n_unique-1) ----
# Round-2 (position-sorted + delta-coded + rank permutation) was tried and
# reverted: shrank pos_abs.bin back to the old scheme's size, but needed a
# new permutation array costing ~log2(n!) (2.27MB on E. coli) to stay
# invertible -- net effect was flat, not an improvement. Reverted to this
# simpler direct-indexing version.
# STAGE 100 SIZE FIX ROUND 3: pos_abs.bin is now fixed-width uint32 (was
# varint) -- xz'd fixed-width consistently beat xz'd varint on real data.
positions = read_u32(f'{WORK}/pos_abs.bin')
lengths = read_u16(f'{WORK}/read_lengths.bin')
n_unique_pos = len(positions)
assert len(lengths) == n_unique_pos
strands = read_bits(f'{WORK}/pos_strand.bin', n_unique_pos)

# ---- orig2uid.bin: original -> unique-id, delta-coded ----
# delta=0 means "first occurrence" (id = running counter, counter++);
# nonzero delta means "duplicate, back-reference = counter - delta".
# Self-describing, no ambiguity -- mirrors the encoder exactly.
with open(f'{WORK}/orig2uid.bin','rb') as f: _d = f.read()
_deltas = struct.unpack(f'<{len(_d)//4}i', _d)
orig2uid = []
_exp = 0
for _dv in _deltas:
    if _dv == 0:
        orig2uid.append(_exp)
        _exp += 1
    else:
        orig2uid.append(_exp - _dv)
n_orig = len(orig2uid)

# ---- mismatch streams: mm_count_per_read.bin is indexed by UNIQUE id
# (0..n_unique-1), mm_ref/obs/pos are sequential in that same ascending order ----
mmcount = read_u16(f'{WORK}/mm_count_per_read.bin')
n_unique = len(mmcount)
with open(f'{WORK}/mm_ref.bin','rb') as f: mm_ref = f.read()
with open(f'{WORK}/mm_obs.bin','rb') as f: mm_obs = f.read()
with open(f'{WORK}/mm_pos.bin','rb') as f: mm_pos = f.read()
mm_offset = [0]*(n_unique+1)
for i in range(n_unique):
    mm_offset[i+1] = mm_offset[i] + mmcount[i]
assert mm_offset[n_unique] == len(mm_ref) == len(mm_obs)

# ---- reconstruct each ORIGINAL read's sequence via its unique-id lookup ----
out_seqs = [None]*n_orig
for o in range(n_orig):
    u = orig2uid[o]
    p, L, s = positions[u], lengths[u], strands[u]
    slice_ = pg[p:p+L]
    seq = list(revcomp(slice_)) if s else list(slice_)
    cnt = mmcount[u] if u < n_unique else 0
    if cnt:
        off = mm_offset[u]
        for m in range(cnt):
            obsc = chr(mm_obs[off+m]); j = mm_pos[off+m]
            if j == 255:
                continue  # capped position, disclosed limitation for reads >255bp
            seq[j] = obsc if not s else COMP[obsc]
    out_seqs[o] = ''.join(seq)

# ---- layer 6: N-reads ----
with open(f'{WORK}/n_reads.txt') as f:
    n_reads = [l.rstrip('\n') for l in f]
n_indices = read_u32(f'{WORK}/n_indices.bin')
assert len(n_reads) == len(n_indices)
n_by_index = dict(zip(n_indices, n_reads))

# ---- final: walk original file order ----
# out_seqs[o] is indexed by NON-N-READ ENCOUNTER ORDER (o=0..n_orig-1), not
# the true original file index -- N-reads were skipped entirely during that
# numbering. The true original index for the o-th non-N read is recoverable
# as the sorted complement of n_indices (N-read positions) within the full
# range.
n_total = n_orig + len(n_indices)
n_index_set = set(n_indices)
non_n_true_indices = [i for i in range(n_total) if i not in n_index_set]
assert len(non_n_true_indices) == n_orig, f"{len(non_n_true_indices)} vs {n_orig}"

out = [None]*n_total
for o in range(n_orig):
    out[non_n_true_indices[o]] = out_seqs[o]
for idx, seq in n_by_index.items():
    out[idx] = seq

missing = sum(1 for x in out if x is None)
print(f"[decode] total reads={n_total} missing={missing}", file=sys.stderr)

with open(OUT,'w') as f:
    for x in out:
        f.write((x or '') + '\n')
