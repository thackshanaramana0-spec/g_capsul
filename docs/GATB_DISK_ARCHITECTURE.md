# How GATB/DiscoSNP++ actually counts k-mers on disk

Read from source 2026-09-04, after my own disk spill wrote 21 GB and filled the
disk. This documents what they do, step by step, and what my version got wrong.
Every claim cites the file and line it came from.

---

## The full chain

### Step 0 — declare a memory ceiling
`ConfigurationAlgorithm.cpp:334-343`

    if (_max_memory == 0) _max_memory = System::info().getMemoryProject();
    if (_max_memory == 0) _max_memory = 5000;                  // MB
    if (_max_memory > (system_mem*2)/3) _max_memory = (system_mem*2)/3;

Their low RAM is a **budget they choose**, not a property they achieve. The rest
of the pipeline is derived from it.

### Step 1 — sample the reads and rank minimizers BY FREQUENCY
`RepartitionAlgorithm.cpp:321, 345-380`

    nbseq_sample = min(estimateSeqNb * 0.05, 50000000);   // 5% of reads
    ... count every m-mer ...
    sort(_counts.begin(), _counts.end());                 // by frequency
    for (i) _freq_order[_counts[i].second] = i;
    repartitor.setMinimizerFrequencies(_freq_order);

A preliminary pass measures m-mer frequencies, then minimizers are chosen by
**frequency rank, not lexicographic order**
(`ComparatorMinimizerFrequencyOrLex`). **I missed this entirely.** It matters
twice over:

* a **rare** minimizer changes less often along a read, so superkmers are
  LONGER, so fewer bytes per k-mer;
* a **frequent** m-mer (`AAAA...`) would otherwise own a huge partition and
  break the memory budget -- the exact skew my top-bits partitioning hit
  (partition 0 held ~24x its uniform share).

### Step 2 — derive the partition count from the budget
`ConfigurationAlgorithm.cpp:407`

    _nb_partitions = ((volume_per_pass * _nb_partitions_in_parallel) / _max_memory) + 1;

Partitions are sized so one fits in the ceiling. If that exceeds
`max_open_files`, they halve parallelism, then add passes (`:413-415`).

### Step 3 — route each superkmer to its partition
`SortingCountAlgorithm.cpp:806-812`

    if ((superKmer.minimizer % _nbPass) == _pass && superKmer.isValid()) {
        size_t p = _repartition (superKmer.minimizer);
        superKmer.save (_superkmerFiles, p);
    }

Two things at once: `% _nbPass` selects this pass's minimizer subset
(**multi-pass** when one pass will not fit), and `_repartition` maps minimizer
to partition through the balanced table from step 1.

### Step 4 — store superkmers 2-bit packed
`Model.hpp:1401`

    int required_bytes = (superKmerLen + kmerSize + 3) / 4;

**This is the whole disk story.** A superkmer holding `L` k-mers spans
`L + k - 1` bases and is written as packed bases, 4 per byte.

| | bytes per k-mer |
|---|---|
| GATB superkmer (L≈11, k=31) | `(11+31+3)/4 / 11` ≈ **0.95** |
| my spill (`uint64 kmer` + `uint32 count`) | **12** |

**~12.6x.** My 21 GB would be ~1.3 GB in their format.

### Step 5 — process one partition at a time
Each partition owns a disjoint minimizer set, so it is counted independently and
never merged against the others. Peak = one partition, not the sum.

---

## What my implementation got wrong, and why

I partitioned on the **top bits of the canonical k-mer** and thought this was a
simplification of their scheme -- it preserves sort order, so concatenating
sorted partitions yields a sorted `kc` with no merge at all. That reasoning was
correct in isolation and wrong in consequence:

1. **It makes superkmers structurally impossible.** Consecutive k-mers in a read
   have unrelated top bits, so they scatter across partitions. Superkmers
   require consecutive k-mers to land *together*, which is precisely what a
   shared minimizer guarantees. **Minimizers are not a partitioning
   convenience -- they are the precondition for superkmer packing.** I traded a
   12.6x disk win for a sort-order convenience worth almost nothing.
2. **Canonical k-mers are skewed toward low values** (`min(v, rc(v))`), so
   partition 0 got ~24x its share. Measured: extrapolating total distinct
   k-mers from partition 0 overestimated by 27x (1.9B vs 69M actual).
3. **No frequency balancing**, so even with more partitions the skew persists.
4. **No multi-pass**, so disk volume is unbounded in input size. It hit 21 GB
   and filled the disk.

## Cost of the honest comparison

| | GATB | mine |
|---|---|---|
| partition key | frequency-ranked minimizer | top bits of canonical k-mer |
| consecutive k-mers | same partition | scattered |
| unit stored | superkmer (L k-mers as L+k-1 bases) | one (kmer,count) per k-mer |
| encoding | 2-bit packed | 8-byte kmer + 4-byte count |
| bytes/k-mer | ~0.95 | 12 |
| balancing | measured m-mer frequencies | none |
| overflow | multi-pass on `minimizer % nbPass` | none |
| disk at full chr20 | ~1.3 GB | **21 GB (filled the disk)** |

## What this means for us

**Porting the disk path is a replacement, not a patch.** It needs the whole
chain -- sampling pass, frequency ranking, minimizer selection, superkmer
grouping, 2-bit packing, balanced repartition, multi-pass. Each step depends on
the one before it; adopting any subset gains nothing, which is exactly what my
partial version demonstrated.

**The transferable insight is 2-bit packing, and it applies in RAM without any
of the rest.** `call_seqs` currently holds 12.6M `std::string`s (~2.42 GB) to
store 4-symbol data at 8 bits per base. Packed 4-per-byte that is ~466 MB, and
it needs no partitioning, no minimizers, and no disk. Same idea as their
`required_bytes` line, applied where we actually hold the memory.

**Already verified from the same insight:** quality only ever answers "is this
base above threshold?", so it is now stored as 1 bit per base instead of an
8-bit phred character -- measured entry RSS 6,821 -> 5,282 MB at full chr20,
with byte-identical calls.

## Status of my disk spill — CORRECTED

I first recorded this as "refuted, a net negative". **That was wrong**, and the
measurement that settled it is:

| | with spill | without spill |
|---|---|---|
| `kc_H_build` | **69.7 s** | 186.4 s |
| peak RSS | **12.35 GB** | 30.04 GB |

**The spill is a 2.4x RAM and 2.7x speed win.** Without it every run stays
resident and the merge runs over many large in-RAM runs under memory pressure.

The disk exhaustion that made me condemn it was **housekeeping, not the
mechanism**: each run created a fresh spill directory and none of the previous
ones were deleted (`kcs3`, `kcs4`, `kcspill`, `kcspill2`, ~21 GB each).
Removing them recovered 35 GB.

So the accurate verdict: the spill WORKS and is worth keeping. What it lacks is
GATB's volume efficiency -- 12 bytes per k-mer against their ~0.95 -- so it
needs ~21 GB of scratch where they need ~1.3 GB. That is a reason to adopt
superkmer packing, not a reason to remove the spill.

**Operational requirement:** the spill directory must be cleaned before each
run. That is a real failure mode -- it silently filled a 233 GB disk in four
runs.
