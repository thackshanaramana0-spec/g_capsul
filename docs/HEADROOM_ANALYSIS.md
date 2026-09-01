# Genuine headroom: where PgRC2's RAM and speed actually come from

File-and-line, with the arithmetic. Nothing here is an impression.

## 1. RAM — the dominant item is one line

**Ours,** `stages/106_inprocess.cpp:279`:

    std::unordered_map<uint64_t, std::vector<uint32_t>> seen; seen.reserve(1u<<21);

One hash node per distinct read (key 8 + vector 24 + hash 8 + next 8 = 48 B)
**plus a separately heap-allocated inner vector** (~32 B of malloc header alone).
That is **80 B per distinct read**, and it exists only to detect duplicates.

**Theirs,** `pseudogenome/generator/ParallelGreedySwipingPackedOverlapPseudoGenomeGenerator.cpp:160-194`:

    if (srIt != blockEnd && compareReads(*(srIt - 1), *srIt) == 0) { ... duplicatesCount++; }

They sort read *indices* into blocks keyed by leading symbol, then compare
adjacent entries. The whole structure is `sortedReadsIdxs`, one flat
`vector<uint32_t>` -- **4 B per read**. No map, no per-key allocation.

And the same sorted-blocked array is then the index the greedy overlap search
walks. **One structure, two jobs.** We build a hash map for dedup AND a separate
prefix index for overlap.

| dataset | reads | our map | theirs | headroom | measured load peak |
|---|---|---|---|---|---|
| E. coli | 1,235,646 | 94 MB | 4.7 MB | **90 MB** | 228 MB |
| L. major | 4,739,289 | 362 MB | 18 MB | **344 MB** | 459 MB (map = 79%) |
| yeast | ~5.9M | 450 MB | 23 MB | **428 MB** | 830 MB peak RSS |
| C. elegans | ~30.7M | 2,342 MB | 117 MB | **2.2 GB** | running |
| T. cacao | ~61.3M | 4,677 MB | 234 MB | **4.4 GB** | running |

Named structures at load total 119 MB against a 228 MB measured peak on
E. coli; the map accounts for essentially the whole 109 MB difference.

### Second RAM item: A4's per-thread slot arrays
`4 B x reads x 12 threads` -- 0.21 GB on L. major, projected 2.74 GB on
T. cacao. Introduced deliberately to kill a 1152 MB candidate-list blow-up, so
it was a good trade at 1.7 GB, but it scales with reads AND threads.

### Third RAM item, only visible at scale: A3's fork
Observed live on C. elegans: parent 3,615 MB and child 2,748 MB resident
simultaneously. Copy-on-write duplicates every page the child writes, and the
child writes round 2, mapping and coding. Invisible below ~2 GB inputs.

## 2. Speed — we are already ahead where it counts

**Their hottest symbol, 14.87%**,
`coders/SymbolsPackingFacility.cpp compareSuffixWithPrefix`:

    while (true) {
        int cmp = (int) reverse[*sufSeq][sufIdx] - reverse[*preSeq][preIdx];
        ...
    }

When the suffix is not byte-aligned they compare **one base per iteration**
through a lookup table. Only the aligned case reaches `compareSequences`.

**Ours,** `stages/106_inprocess.cpp` `rcmp`:

    while(L>=32){ if(w32(pa)!=w32(pb)) return false; pa+=32; pb+=32; L-=32; }

**32 bases per 64-bit comparison at arbitrary offset.** This is why our assembly
is 1.64 s against their 2.75 s despite their better parallel utilisation.

**Conclusion: there is no speed headroom in comparison or assembly. We are
1.7x ahead there and should not touch it.**

The speed headroom is entirely in two places:

| | ours | PgRC2 | headroom |
|---|---|---|---|
| compression as share of profile | ~16.7% | ~2.5% | ~14 pp = ~5.4 CPU-s of 38.5 on E. coli |
| cores utilised | 3.76 of 12 | 5.79 of 12 | 54% more throughput at their efficiency |

Both trace to the same root: the 7-way probe, which does 1.86x their total CPU
work and generates the page-fault traffic.

## 3. Where this can realistically end

- **RAM:** replacing the dedup map with a sorted-blocked index removes 90 MB
  (E. coli) to 4.4 GB (T. cacao) and makes the overlap index free. That alone
  moves us from ~2.5x their peak toward parity.
- **Speed:** if coding falls from 38% of wall to their ~10%, E. coli goes from
  8.35 s to roughly 6 s against their 3.58 s. Combined with utilisation, a
  realistic floor is **1.5-1.7x their wall**, not parity -- because their
  remaining advantage is breadth of parallelism (31 OpenMP regions to our 1),
  not any single hot loop.
- **Assembly: leave alone.** We win it.

## 4. What is NOT headroom, measured
- `pos_abs` coding: already at 0.92-0.99x its set-plus-permutation bound.
- Symbol packing: both sides are exactly 2 bits/base, 56 MB on E. coli.
- Comparison inner loop: ours is 32x theirs per operation.
