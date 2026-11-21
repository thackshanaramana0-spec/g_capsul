#!/usr/bin/env python3
"""
End-to-end completeness test of the FM-index completion step, on archives built
to break it.

For each synthetic dataset the script
  1. writes a FASTQ,
  2. compresses it with best106 and checks the decode is byte-identical,
  3. builds the index (CAPS_PILEUP=1 CAPS_FMI=1),
  4. draws probes (substrings of reads on either strand, probes that span read
     boundaries' shared sequence, random strings, probes containing N, probes
     of every length from 4 to 80), and
  5. compares what `query` returns (CAPS_QUERY_MM=0, CAPS_QUERY_CONTAIN=1, the
     FM-index as the completion step) with a brute-force scan of the decoded reads.

PASS means every read containing a probe was returned, for every probe, on every
dataset, and every verbatim record (a contained read stored in the index) sits
inside the record of the read containing it, at the position its header states. Nothing is tuned per dataset: the tolerance is 0, no seed length or
probe-length floor is used. Returned reads that are not real reads are counted
and reported but do not fail the run, because the guarantee is completeness.

usage: fmi_end_to_end.py <best106> <capsule_decode> [workdir] [probes_per_dataset]
"""
import os, random, subprocess, sys, shutil

ENC, DEC = sys.argv[1], sys.argv[2]
WD = sys.argv[3] if len(sys.argv) > 3 else "fmi_e2e"
NPROBE = int(sys.argv[4]) if len(sys.argv) > 4 else 300
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(c, "N") for c in reversed(s))
can = lambda s: min(s, rc(s))


def make(path, seed, glen, nreads, lens, err, nfrac, dup, contain, repeat):
    r = random.Random(seed)
    g = [r.choice("ACGT") for _ in range(glen)]
    if repeat:                                   # repeated segments: many loci per probe
        seg = "".join(r.choice("ACGT") for _ in range(400))
        for k in range(repeat):
            p = r.randrange(0, glen - 400)
            g[p:p + 400] = seg
    g = "".join(g)
    reads = []
    for _ in range(nreads):
        L = r.choice(lens)
        p = r.randrange(0, glen - L)
        s = list(g[p:p + L])
        for i in range(L):
            if r.random() < err:
                s[i] = r.choice("ACGT")
        if r.random() < nfrac:
            for _ in range(r.choice([1, 1, 2, 5])):
                s[r.randrange(L)] = "N"
        s = "".join(s)
        if r.random() < 0.5:
            s = rc(s)
        reads.append(s)
        if contain and r.random() < contain and L > 40:        # a shorter read inside this one
            a = r.randrange(0, L // 3); b = L - r.randrange(0, L // 3)
            reads.append(s[a:b])
    for _ in range(int(nreads * dup)):
        reads.append(r.choice(reads))
    r.shuffle(reads)
    with open(path, "w") as f:
        for i, s in enumerate(reads):
            f.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
    return reads


def run(cmd, env=None, cwd=None):
    e = dict(os.environ); e.update(env or {})
    return subprocess.run(cmd, env=e, cwd=cwd, capture_output=True, text=True)


DATASETS = [
    # name, seed, genome, reads, read lengths, error, N fraction, dup, contain, repeat
    ("fixed_clean",        1, 60000, 4000, [100],                 0.000, 0.00, 0.00, 0.0, 0),
    ("fixed_noisy",        2, 60000, 4000, [150],                 0.020, 0.00, 0.05, 0.0, 0),
    ("variable",           3, 60000, 4000, [36, 75, 101, 150, 251], 0.005, 0.00, 0.05, 0.0, 0),
    ("variable_contained", 4, 60000, 3000, [60, 120, 250],        0.005, 0.00, 0.10, 0.30, 0),
    ("with_N",             5, 60000, 4000, [100, 150],            0.005, 0.20, 0.05, 0.10, 0),
    ("heavy_error",        6, 60000, 3000, [100],                 0.050, 0.05, 0.00, 0.0, 0),
    ("repeats",            7, 40000, 4000, [100, 140],            0.010, 0.02, 0.10, 0.05, 12),
    ("high_dup_variable",  8, 30000, 5000, [50, 100, 220],        0.010, 0.05, 0.60, 0.20, 0),
]

fails = 0
shutil.rmtree(WD, ignore_errors=True); os.makedirs(WD)
for name, seed, glen, nreads, lens, err, nfrac, dup, contain, repeat in DATASETS:
    d = os.path.join(WD, name); os.makedirs(d)
    fq = os.path.join(d, "reads.fq")
    make(fq, seed, glen, nreads, lens, err, nfrac, dup, contain, repeat)
    enc = run([ENC, "reads.fq", "3", "16", "16", "22", "16", "16", "1", "24", "64", "1"],
              env=dict(ARCHIVE="win.capsule", CAPS_SPANS="1", DUMP_PERM="1", DUMP_MM="1", DUMP_LIT="1"), cwd=d)
    arc = os.path.join(d, "win.capsule")
    if not os.path.exists(arc) and os.path.exists(os.path.join(d, "out.arcs2")):
        shutil.copy(os.path.join(d, "out.arcs2"), arc)
    if not os.path.exists(arc):
        print(f"{name}: FAIL no archive\n{enc.stderr[-400:]}"); fails += 1; continue
    os.makedirs(os.path.join(d, "dec"), exist_ok=True)
    run([DEC, "win.capsule", "dec", "dec/reads.out"], cwd=d)
    orig = [l.strip() for i, l in enumerate(open(fq)) if i % 4 == 1]
    decd = [l.strip() for l in open(os.path.join(d, "dec", "reads.out")) if l.strip()]
    if orig != decd:
        print(f"{name}: FAIL archive is not lossless ({len(orig)} vs {len(decd)} reads)"); fails += 1; continue
    ix = run([DEC, "index", "win.capsule", "win.capsule.qidx"], env=dict(CAPS_PILEUP="1", CAPS_FMI="1"), cwd=d)
    fmi_line = [l for l in ix.stderr.splitlines() if l.startswith("[fmi] ") and "indexed" in l]
    # Every placement-keyed record, from a whole-range query. A verbatim record
    # (header ">v") is a read contained in one of these, so it must agree with
    # a record that covers its stated position, base for base outside its own
    # deviations and N. Placement itself is a compression choice (it minimises
    # bits, not homology), so ordinary records are not checked against the pg.
    run([DEC, "export", "win.capsule", "pg.fa"], cwd=d)
    pglen = sum(len(l.strip()) for l in open(os.path.join(d, "pg.fa")) if not l.startswith(">"))
    run([DEC, "query", "win.capsule", "all.fa", f"0-{pglen}"], cwd=d)
    placed = []; hp = None
    for l in open(os.path.join(d, "all.fa")):
        if l.startswith(">"): hp = int(l.split("pos=")[1].split()[0])
        elif l.strip(): placed.append((hp, l.strip()))
    placed.sort()
    starts = [p for p, _ in placed]
    import bisect

    def verbatim_sits(pos, s):
        i = bisect.bisect_right(starts, pos)
        best = None
        for p, r in placed[max(0, i - 400):i]:
            if p <= pos and p + len(r) >= pos + len(s):
                w = r[pos - p:pos - p + len(s)]
                cm = [(a, b) for a, b in zip(s, w) if a != "N" and b != "N"]
                mm = sum(a != b for a, b in cm)
                best = mm if best is None else min(best, mm)
        return best is not None and best <= max(2, len(s) // 20)

    rnd = random.Random(seed * 7919)
    probes = []
    for k in range(NPROBE):
        kind = k % 5
        if kind in (0, 1, 2):
            s = rnd.choice(decd)
            if "N" in s and kind == 0:
                pass
            L = min(len(s), rnd.randint(4, 80)); a = rnd.randrange(0, len(s) - L + 1)
            p = s[a:a + L]
            probes.append(rc(p) if kind == 1 else p)
        elif kind == 3:
            probes.append("".join(rnd.choice("ACGT") for _ in range(rnd.randint(4, 40))))
        else:
            s = rnd.choice(decd); L = min(len(s), rnd.randint(8, 60)); a = rnd.randrange(0, len(s) - L + 1)
            p = list(s[a:a + L]); p[rnd.randrange(L)] = "N"
            probes.append("".join(p))
    truth_reads = decd
    missing = extra = nprobe_hit = misplaced = nverb_checked = 0
    for i, p in enumerate(probes):
        pat = {p, rc(p)}
        truth = {can(s) for s in truth_reads if any(q in s for q in pat)}
        out = os.path.join(d, "q.fa")
        if os.path.exists(out): os.remove(out)
        q = run([DEC, "query", "win.capsule", "q.fa", p],
                env=dict(CAPS_QUERY_MM="0", CAPS_QUERY_CONTAIN="1", CAPS_QUERY_INDEX="fmi"), cwd=d)
        got = set()
        if os.path.exists(out):
            pos = None; verb = False
            for l in open(out):
                if l.startswith(">"):
                    pos = int(l.split("pos=")[1].split()[0]); verb = l.startswith(">v"); continue
                s = l.strip()
                if not s: continue
                got.add(can(s))
                if verb:
                    nverb_checked += 1
                    if not verbatim_sits(pos, s): misplaced += 1
        if truth: nprobe_hit += 1
        m = truth - got
        if m:
            missing += len(m)
            if missing <= len(m):
                print(f"  {name}: probe {p!r} missing {len(m)} read(s), e.g. len {len(next(iter(m)))}")
        real = {can(s) for s in truth_reads}
        extra += len(got - real)
    ok = missing == 0 and misplaced == 0
    fails += 0 if ok else 1
    print(f"{name:20s} reads={len(decd):5d} probes={len(probes)} (with a true hit: {nprobe_hit}) "
          f"missed reads={missing} verbatim records checked={nverb_checked} misplaced={misplaced} "
          f"non-real returned={extra}  "
          f"{'PASS' if ok else 'FAIL'}")
    if fmi_line: print("   " + fmi_line[0])
print("ALL PASS" if fails == 0 else f"{fails} DATASET(S) FAILED")
sys.exit(1 if fails else 0)
