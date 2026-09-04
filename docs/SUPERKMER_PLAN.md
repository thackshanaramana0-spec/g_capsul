# Superkmer spill — implementation plan

Written 2026-09-04 after two full-scale runs failed by exhausting a 233 GB disk
(spill reached 25 GB, twice). This is the plan to close the 19x disk gap, not a
sketch: every step names what is built, what it depends on, and how it is
verified.

---

## 1. The gap, stated exactly

| | bytes written per k-mer |
|---|---|
| GATB superkmer | `(L + k + 3)/4 / L` with L≈11, k=31 → **~0.95** |
| our spill today | 8-byte k-mer + 4-byte count → **12** |
| our spill after dropping the count | 8-byte k-mer per occurrence → **8** |

Full chr20 has ~1.4 billion k-mer occurrences:

| format | spill volume |
|---|---|
| (kmer,count) — what failed twice | ~25 GB |
| key-only (built, unverified) | ~11 GB |
| **superkmer** | **~1.3 GB** |

## 2. Why the count-drop alone is not enough

It takes 25 GB to ~11 GB, which fits today's 35 GB free — but it is fragile: a
larger input, or a fuller disk, fails again. Superkmers are the only version
that is bounded in the way GATB's is.

## 3. What superkmers require — the dependency chain

Each step is a precondition for the next. **This is why my key-range version
gained nothing: it broke the chain at step 1.**

    1. minimizer per k-mer          (consecutive k-mers usually share one)
    2. group consecutive k-mers by shared minimizer  -> superkmer
    3. partition BY MINIMIZER       (so a whole superkmer lands in one file)
    4. store as 2-bit packed bases  (L + k - 1 bases, not L k-mers)
    5. read back: expand -> canonicalise -> sort -> count

Partitioning by anything other than the minimizer scatters neighbouring k-mers
and makes step 2 impossible, which is exactly what top-bits partitioning did.

## 4. Implementation, step by step

### 4.1 Rolling minimizer (the one genuinely new algorithm)

For each k-mer, the minimizer is the smallest m-mer inside it (GATB uses m=10
for k=31). Computing it naively is O(k-m) per k-mer; a **monotonic deque**
makes it O(1) amortised:

    for each position i:
        push m-mer(i) onto the back, popping any back element >= it
        pop the front if it has slid out of the window
        front = current minimum

Ordering: GATB ranks minimizers by **measured frequency** (5% sampling pass,
`RepartitionAlgorithm.cpp:345-380`) so rare m-mers win, making superkmers longer
and partitions balanced. **v1 uses lexicographic order** — simpler, and the cost
is shorter superkmers and some partition skew, both measurable. Frequency
ranking is a later refinement, not a precondition.

### 4.2 The strand problem — the part that is genuinely hard

`kc` stores CANONICAL k-mers. Two consecutive k-mers can canonicalise to
opposite strands, so their packed bases are not contiguous in one orientation
and cannot be stored as one overlapping sequence.

GATB handles this with `which()` and the kxmer machinery: it tracks runs of
k-mers on the same strand and treats a strand flip as a run boundary
(`SortingCountAlgorithm.cpp:830-845`).

**Plan:** store superkmers in FORWARD read orientation (no canonicalisation on
write) and canonicalise only on read-back, when each k-mer is expanded. The
minimizer is then computed on the forward sequence too. This sidesteps the
strand problem entirely at the cost of storing both strands' k-mers
separately — which we already do, since every read is scanned once forward.

### 4.3 On-disk format

    [uint8 L]  [ceil((L + k - 1) / 4) bytes of 2-bit packed bases]

`L` <= 255 caps a superkmer at 255 k-mers; longer runs split, costing one extra
header. Bases pack 4/byte, matching `Model.hpp:1401`.

### 4.4 Read-back and the index

Partitions are keyed by minimizer, so they are **not** key-ranges and
concatenating them does not give a sorted `kc`. Two options:

* **(a)** sort `kc` globally at the end — one 140M-element sort, expensive;
* **(b)** keep one sorted array per partition and make `kc_find(kmer)` compute
  the query's minimizer, select that partition, and binary search inside it.

**Choose (b).** It is O(log n) exactly as now, needs no global sort, and the
minimizer of a query k-mer is computed in O(k-m) once per lookup. It does mean
`kc_find` changes, and every caller of it must keep working.

## 5. Verification — three independent checks

This is the part that must not be skipped, because the whole accuracy result
rides on `kc` being correct.

1. **k-mer set identity.** `kc nodes=` must equal 1,063,607 on the r2 window and
   140,719,632 at full chr20, exactly as today.
2. **Count identity.** Dump `(kmer,count)` with and without superkmers and
   `cmp` the sorted files. Counts must match exactly, not approximately —
   `MINC` and the coherence filter both read them.
3. **End-to-end F1.** The r2 window must still score **TP=331 FP=18 F1=0.881**.
   Any deviation means the k-mer set changed and the accuracy claim is at risk.

Step 3 is the one that protects the claim; steps 1-2 localise a failure if it
occurs.

## 6. Risk to the accuracy win — and how it is contained

The accuracy result (full chr20 F1 0.876 vs DiscoSNP++ 0.847) depends only on
`kc`'s contents. Superkmers change **how `kc` is built**, never what it
contains, so a correct implementation is invisible to accuracy — and check 2
proves correctness directly rather than by inference.

**Containment:** superkmer spill is a third mode behind
`CAPS_KC_SUPERKMER=1`, leaving both the in-RAM path (default) and the current
key-only spill untouched. If it fails verification it is disabled, not debugged
in place, and the accuracy claim is never exposed.

## 7. Order of work

| # | step | verify |
|---|---|---|
| 1 | rolling minimizer, standalone | brute-force min over the window on 10K reads |
| 2 | superkmer grouping + packing/unpacking | round-trip: pack then unpack must return the same k-mers |
| 3 | wire into spill behind the flag | `kc nodes=` identical on the window |
| 4 | per-partition `kc_find` | end-to-end F1 = 0.881 on r2 |
| 5 | full chr20 | disk volume, peak RSS, wall time, F1 |

Nothing proceeds to the next row until the current one verifies.

## 8. Honest expectation

Disk ~11 GB -> ~1.3 GB. **RAM and wall time are NOT expected to improve much
from this alone** — the spill already delivered its RAM win (12.35 GB vs
30.04 GB without it). Superkmers make the spill *sustainable*, which is what
turns it from something that fills a disk into something that can be claimed.

The separate speed lever is radix pre-binning (item 1 of the earlier list),
which attacks the 70 s sort directly and does not depend on any of this.

---

# EXECUTION RESULT (2026-09-04)

All five items built. Every gate in §5 passed, and accuracy never moved.

## Verification, in the order the plan specified

| step | check | result |
|---|---|---|
| 1 | rolling minimizer vs brute force, 1,180,000 positions | **0 mismatches** |
| 2 | superkmer pack/unpack, k-mer multiset identity | **identical** (2,360,000 both) |
| 3 | `kc nodes` in-RAM vs superkmer spill | **1,063,607 both** |
| 4 | end-to-end F1 on r2 | **TP=331 FP=18 F1=0.881, identical** |
| 5 | full chr20 | spill **2.4 GB** vs 25 GB (running) |

## Measured mechanism numbers

| | before | after |
|---|---|---|
| spill bytes/k-mer | 12 | **1.161** (10.3x) |
| full-chr20 spill volume | 25 GB (filled the disk twice) | **~2.4 GB** |
| mean superkmer length | -- | 9.64 k-mers |
| read storage | ~192 B/read | ~80 B/read (2-bit packed) |
| quality storage | 8 bits/base | **1 bit/base** |
| batch sort | one large n log n | 256 cache-resident radix bins |

Our 1.161 B/k-mer against GATB's ~0.95 is explained: they rank minimizers by
measured frequency so superkmers run longer (L≈11 vs our 9.64 with lexicographic
order). That refinement is not built and is the remaining gap on this axis.

## The bug the identity gate caught

The first superkmer run produced **1,063,514** k-mers against 1,063,607 -- 93
short, 0.009%, and invisible in F1. Cause: after an `N` the minimizer deque
refills 10 bases later and was assigning minimizers to k-mers whose 31-base
window still contained that `N`; those superkmers were then rejected wholesale.
Fixed by tracking the last invalid position (`lastN`) and skipping any k-mer
whose window has not cleared it.

**This is the argument for multiset identity as the gate rather than "F1 looks
fine".** A 0.009% silent divergence between code paths would have shipped.

## What is NOT done

- **Frequency-ranked minimizers** (GATB's 5% sampling pass). Would lengthen
  superkmers toward their ~0.95 B/k-mer and balance partitions.
- **Multi-pass** on `minimizer % nbPass`. Without it spill volume is bounded
  only by input size; at chr20 scale 2.4 GB is fine, at WGS scale it would need
  revisiting.
- **Flat-buffer reads.** Packing cut ~192 -> ~80 B/read, but the 32-byte
  `std::string` header still dominates a 37-byte payload; a flat buffer with
  offsets would remove it.
