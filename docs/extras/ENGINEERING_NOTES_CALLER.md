# Engineering notes: the variant caller (include/caps_caller.h)

Long explanatory comments from `include/caps_caller.h` that were condensed in the source during the documentation revision. They are kept here in their original wording (with a few phrases edited for tone) because they record measurements, rejected approaches and the reasoning behind design choices. The source keeps a short statement of each conclusion. Measurement labels such as `HG002 r2` refer to the chr20 tuning windows described in `docs/claim2/`. Working notes cited by file name (for example `INDEL_LOSS_SKELETAL.md`) are not part of this repository; the numbers quoted beside each citation are the evidence.


## Tolerant flank re-convergence (flank_match_tol)

Source: `include/caps_caller.h`, lines 237-248 before the documentation revision.

Original comment:

```
// REFUTED, kept for the record (not wired). Tolerant re-convergence test.
// Tried because:  on HG002 r2, 29 of the 32
// MISSED indels already have a contig PAIR covering them -- the substrate is
// fine, the extractor is what fails. It demanded a byte-identical 15 bp flank,
// so a single heterozygous SNV anywhere in that flank destroyed the bubble.
// DiscoSNP++ explicitly models this case (its own header calls a bubble
// "an indel+n SNPs"), which is part of why its indel recall beats ours.
// Tolerance reuses the project's already-frozen HDMAX=2 rather than
// introducing a new constant.
// RESULT: net NEGATIVE -- r2 indel F1 0.602->0.595, r5 0.600->0.582,
// HG005 0.542->0.508. Recall rose slightly but precision fell more: the
// bubble GEOMETRY is not what limits us, and loosening it only admits FPs.
```

What the code now says:

```
// Not called. Re-convergence test that tolerates up to `tol` mismatches in the
// flank, so that a heterozygous SNV next to an indel does not break the bubble.
// It lowered indel F1 on every chr20 test window (precision fell more than recall
// rose), so extract_bubble keeps the exact flank_match. Measurements: docs/extras/ENGINEERING_NOTES_CALLER.md.
```

## Homopolymer run-length branch

Source: `include/caps_caller.h`, lines 275-287 before the documentation revision.

Original comment:

```
    // MEASURED MOTIVATION (working note INDEL_LOSS_SKELETAL.md): every one of the 9
    // truth indels DiscoSNP++ finds on HG002 r2 that we miss is a 1 bp indel
    // inside a homopolymer run of 7-8 identical bases. The generic path below
    // cannot represent them: a 1 bp change in an 8-mer run shifts its LENGTH
    // (8 -> 7 or 9), the two haplotypes stay identical through the whole run,
    // and the byte-identical `flank_match` never succeeds at any offset g
    // because the flanks are themselves shifted by the length difference.
    // DiscoSNP++ resolves this by extending one graph path base-by-base and
    // taking the shortest closing extension, plus an explicit ambiguity
    // allowance (max_ambigous_indel, default 20).
    // Here the equivalent is direct: when the divergence sits at a
    // homopolymer, measure the RUN LENGTH on each side and emit the difference
    // as the indel, then verify re-convergence AFTER both runs.
```

What the code now says:

```
    // A 1 bp indel inside a homopolymer changes the run LENGTH, so the two
    // haplotypes stay identical through the run and their flanks are offset by the
    // length difference. No shift g in the generic search below then gives an exact
    // flank match. Here the run length is measured on each side, the difference is
    // emitted as the indel, and re-convergence is checked after both runs. (DiscoSNP++
    // reaches the same events by extending a graph path base by base.) These events
    // were the whole of the 1 bp indels missed on the HG002 tuning window; see
    // docs/extras/ENGINEERING_NOTES_CALLER.md.
```

## Short tandem repeat run-length branch

Source: `include/caps_caller.h`, lines 323-360 before the documentation revision.

Original comment:

```
    // Generalizes the homopolymer branch above from unit length 1 to unit
    // lengths 2..6. MEASURED MOTIVATION (working note HET_INDEL_FRESH_SCAN.md, 2026-09-03):
    // a fresh scan of the raw FP/FN records found a real disagreement at
    // chr20:3332481-3332512, a (TTTA)n tetranucleotide repeat that GIAB's own
    // truth annotates `difficultregion=AllTandemRepeats_lt51bp_slop5`. The
    // homopolymer branch above cannot reach it: its guard requires a SINGLE
    // repeating nucleotide (`rc0 == B[ib-1] && b2i(rc0) >= 0`), so a 4 bp
    // repeat UNIT never qualifies. The generic flank_match loop below cannot
    // reach it either, for exactly the reason the homopolymer comment already
    // gives: when a whole repeat unit is gained or lost, the flanks on the two
    // sides are themselves shifted by the length difference, so a
    // byte-identical flank match fails at every offset g.
    //
    // The fix is the same run-length measurement the homopolymer branch already
    // uses and that is already validated -- count how far each side continues
    // the periodic unit and emit the difference -- just tiled by a k-mer unit
    // instead of a single character. diff is computed from the run lengths, so
    // its SIGN (insertion vs deletion) is correct by construction, unlike the
    // generic loop, which accepts whichever g happens to satisfy flank_match
    // first and has no length/direction guarantee inside a repeat.
    //
    // Deliberately starts at U=2: U=1 is the homopolymer branch's territory and
    // is left completely untouched so this cannot regress that already-tuned
    // path. Requires >=2 unit copies before the divergence so this only fires
    // at genuine repeats, not at coincidental k-mer equality.
    //
    // ── REFUTED BY MEASUREMENT, 2026-09-03. OPT-IN (CAPS_STRBUBBLE=1), OFF by
    // default. The reasoning above is sound and the branch DOES engage --
    // instrumented at 682 firings on HG002 r2 -- but the caller's output is
    // BYTE-IDENTICAL with it on and off (`cmp` on calls.vcf), and the scored
    // result is unchanged to the digit (INDEL P=0.704 R=0.576 F1=0.633 both
    // ways). So the STR loci it reaches were already being resolved
    // equivalently by the generic loop, or its candidates die in the same
    // downstream filters every other candidate dies in -- consistent with the
    // standing structural finding (working note HOW_DISCOSNP_WINS.md sec 4) that the
    // indel gap is not in bubble GEOMETRY but in what the substrate can offer
    // a filter to reject. Kept behind a flag, not deleted, per this project's
    // rule that refuted ideas stay on the record with their measurement.
```

What the code now says:

```
    // Generalizes the homopolymer branch above from unit length 1 to repeat units of
    // 2 to 6 bases: count how far each side continues the unit and emit the length
    // difference, so the sign (insertion or deletion) is correct by construction.
    // U=1 is left to the homopolymer branch. At least two copies of the unit must
    // precede the divergence, so the branch fires only inside a repeat.
    //
    // OPT-IN (CAPS_STRBUBBLE=1), OFF by default: the branch fires, but the caller's
    // VCF was byte-identical with it on and off on the HG002 tuning window, so the
    // loci it reaches are already resolved by the generic search. See
    // docs/extras/ENGINEERING_NOTES_CALLER.md.
```

## Whole-pair gapped alignment scan (scan_pair)

Source: `include/caps_caller.h`, lines 407-416 before the documentation revision.

Original comment:

```
// MEASURED MOTIVATION (working note INDEL_LOSS_SKELETAL.md §8): both allele-contigs
// DO exist (84% of placed reads match their contig perfectly), but the two are
// never paired at the event, because near a homopolymer every 25-mer
// overlapping the run differs between haplotypes, so the shared anchors that
// survive sit 213-135,453 bp away. Walking outward from such an anchor dies on
// ~335 bp contigs.
// This removes both constraints at once: take the diagonal implied by ANY
// shared anchor and scan the ENTIRE co-linear region, resolving each divergence
// with a gap and continuing past it on the updated diagonal. Events anywhere
// along the pair are then reachable from one distant anchor.
```

What the code now says:

```
// Both allele contigs usually exist, but near a homopolymer every 25-mer that
// overlaps the run differs between haplotypes, so the nearest shared anchor can be
// hundreds of bases away and an outward walk from it ends on short contigs. This
// takes the diagonal implied by ANY shared anchor and scans the whole co-linear
// region, resolving each divergence with a gap and continuing on the updated
// diagonal, so an event anywhere along the pair is reachable from one anchor.
```

## Rolling canonical 25-mer

Source: `include/caps_caller.h`, lines 499-503 before the documentation revision.

Original comment:

```
// Every hot k-mer walk in this file called pack25() and rc25() at each
// position: 25 loop iterations each, 50 per position, to produce a value that
// differs from its predecessor by ONE base. The read scan in the indel pass
// alone walks 4M reads x ~85 positions -- ~17 billion iterations to compute
// 340M k-mers.
```

What the code now says:

```
// Packing each 25-mer from scratch costs 25 steps per strand per position, although
// consecutive windows differ by one base. Roll25 keeps both strands and updates them
// in O(1) per base.
```

## Calling substrate rebuild: motivation

Source: `include/caps_caller.h`, lines 548-565 before the documentation revision.

Original comment:

```
// MEASURED MOTIVATION (2026-09-02, HG002 chr20 r2, 422 truth het-SNV sites):
//   * the reads contain the variants: 99.5% of truth sites carry BOTH alleles
//     with >=3 reads each at depth >=6 -- the ceiling is not the data;
//   * per-contig depth is NOT the blocker: 96.2% of sites already have >=6
//     reads stacked on one contig;
//   * the blocker is READ->CONTIG ASSIGNMENT. CAPSULE assigns each read to the
//     chain it EXACTLY overlaps, so ref-allele and alt-allele reads land on
//     different chains and never meet in a pileup. Only 124 of 75,115 reads
//     ever reached the mismatch-tolerant pigeonhole mapper -- chaining claimed
//     the rest. Re-assigning reads mismatch-tolerantly makes 70.1% of truth
//     sites immediately show a both-allele stack.
//   * reads covering one site are additionally spread over a MEDIAN of 5
//     contigs, so even a correct assignment splits the evidence -- hence the
//     collapse pass below, which keeps one representative contig per locus.
//
// Both passes are gated and default ON for calling only; they never touch the
// compression path. CAPS_NO_REMAP=1 / CAPS_NO_COLLAPSE=1 restore the previous
// behaviour exactly, so every measurement has a revert.
```

What the code now says:

```
// Why the substrate is rebuilt. The encoder assigns each read to the chain it
// overlaps EXACTLY, so reads carrying the reference and the alternate allele land on
// different chains and never meet in one pileup, and reads covering one site are
// spread over several redundant contigs. The rebuild therefore (1) collapses
// redundant contigs to one representative per locus and (2) re-places every read
// with mismatches allowed. On the HG002 tuning window the reads themselves carried
// both alleles at almost every truth site; the assignment was the limit. Figures in
// docs/extras/ENGINEERING_NOTES_CALLER.md.
//
// Both passes run for calling only and never touch the compression path.
// CAPS_NO_REMAP=1 / CAPS_NO_COLLAPSE=1 disable them.
```

## Flat open-addressing k-mer set

Source: `include/caps_caller.h`, lines 587-604 before the documentation revision.

Original comment:

```
// `collapse_contigs` was 51.8 s of a 207.6 s archive-path run (25%), and a
// benchmark of its exact operation mix showed it is essentially 100% hash
// INSERTION at std::unordered_set's floor: 140M inserts + 28M queries measured
// 69.6 s there.
//
// The same operation mix on a flat open-addressing table measures **4.4 s** --
// 15.8x. unordered_set allocates a node per element and chases a pointer per
// probe; this is one contiguous array, one cache line per probe, no allocation.
//
// SEMANTICS ARE IDENTICAL: incremental insert and query, which is what the
// collapse walk needs (each contig queries the set as it stands, then adds to
// it). This is a container change, not an algorithm change -- the accept/reject
// decision is bit-for-bit the same.
//
// Two earlier attempts at this stage failed BECAUSE they targeted the queries
// (4% of the cost): caching the k-mers gained 5 s for +1.7 GB, and a Bloom
// pre-filter cost 140M adds to guard 28M lookups and made the run WORSE. The
// container itself was the cost all along.
```

What the code now says:

```
// The collapse walk is almost entirely hash insertion. A flat open-addressing table
// (one contiguous array, one cache line per probe, no per-element allocation) runs
// the same insert/query mix about 16x faster than std::unordered_set. Semantics are
// identical: incremental insert and query, so the accept/reject decision for each
// contig is unchanged. Measurements: docs/extras/ENGINEERING_NOTES_CALLER.md.
```

## Huge-page hint ordering

Source: `include/caps_caller.h`, lines 607-611 before the documentation revision.

Original comment:

```
// ORDER MATTERS, and this project has already shipped an inert version of this
// hint by calling madvise AFTER the array was filled -- pages already faulted
// in are not remapped. Every call site must sit immediately after the
// allocation and before the first write, and the log line says which arrays
// actually got it.
```

What the code now says:

```
//
// ORDER MATTERS: madvise only affects pages not yet faulted in, so each call must sit
// immediately after the allocation and before the first write. The log line names
// the arrays that received the hint.
```

## FlatKmerSet sizing

Source: `include/caps_caller.h`, lines 648-658 before the documentation revision.

Original comment:

```
    // `init` used to round total k-mer POSITIONS over ALL contigs up to the
    // next power of two and double it. On the 4M-read subset that is a 2^29
    // slot table -- 4.29 GB -- holding only the DISTINCT k-mers of the contigs
    // that actually survive the collapse. Measured load factor: ~8%. So every
    // one of ~125M probes was a cold miss over a 4 GB region that is 92% empty,
    // and 4 GB of RSS was reserved for it, twice.
    //
    // Start from a modest table and DOUBLE when it passes 50% load. A hash
    // set's contents do not depend on its capacity, so the accept/reject
    // decisions -- and therefore which contigs survive -- are bit-for-bit
    // unchanged; only the footprint the probes range over shrinks.
```

What the code now says:

```
    //
    // The table starts small and doubles past 50% load. A hash set's contents do not
    // depend on its capacity, so which contigs survive the collapse is unchanged;
    // only the memory the probes range over shrinks. (A table sized from total k-mer
    // positions ran at about 8% load over several GB.)
```

## FlatKmerSet initial size

Source: `include/caps_caller.h`, lines 663-678 before the documentation revision.

Original comment:

```
        // START WHERE IT ENDS.
        //
        // Capping the initial size at 2^26 meant the table always grew exactly
        // once: 512 MB allocated, then 1024 MB allocated while the 512 MB is
        // still live to be rehashed out of -- a 1.5 GB transient to reach a
        // 1.0 GB table. Measured, the collapse step was the largest single
        // jump in the whole run at +2076 MB.
        //
        // `cap` is the total k-mer POSITIONS over the input contigs and the
        // distinct count cannot exceed it, so the largest power of two not
        // exceeding cap is a slot count the set always fits in, and in
        // practice lands near 45% load (k-mers drawn from overlapping windows
        // repeat). It is also exactly where the doubling ended up anyway --
        // 2^27 both before and after -- so this changes only where the table
        // STARTS, never its final size or its contents. grow() remains as the
        // safety net.
```

What the code now says:

```
        // Start at the largest power of two not exceeding `cap` (total k-mer positions,
        // an upper bound on distinct k-mers). This lands near 45% load and avoids a
        // grow-and-rehash step with both tables live. Only the starting size changes,
        // never the contents; grow() remains as the safety net.
```

## Collapse sort tie-break

Source: `include/caps_caller.h`, lines 716-723 before the documentation revision.

Original comment:

```
    // NOTE: adding an index tie-break here (to remove a dependence on
    // std::sort's unspecified order for equal-length contigs) was MEASURED and
    // REVERTED. It is a real latent nondeterminism, but breaking the tie by
    // index CHANGES WHICH CONTIGS SURVIVE -- 140,561 -> 140,659 on the 4M
    // subset -- and every downstream stage then does more work, eating ~130 s
    // of the 152 s this function's other fixes save. Correctness-neutral it is
    // not, so it is not bundled with a speed change. If it is wanted, it must
    // be measured on its own against F1.
```

What the code now says:

```
    // KNOWN NONDETERMINISM, left as is: equal-length contigs are ordered by
    // std::sort's unspecified order. Adding an index tie-break changes which contigs
    // survive the collapse and therefore the downstream work and calls, so it is a
    // behaviour change that would need its own F1 measurement, not a cleanup.
```

## Collapse: one k-mer walk per contig

Source: `include/caps_caller.h`, lines 747-765 before the documentation revision.

Original comment:

```
        // ── COMPUTE THE CANONICAL k-MERS ONCE PER CONTIG ────────────────
        // This used to walk the contig TWICE -- once on a stride-5 sample to
        // query, then again over every position to claim -- recomputing
        // pack25 + rc25 both times. The claim pass is a superset of the query
        // pass, so one walk produces both: collect the canonical k-mers in
        // order, and the query is every 5th entry of what was collected.
        //
        // Structural, not a tuning change: the SET of k-mers examined and the
        // SET claimed are unchanged, so `hit`, `tot` and the accept/reject
        // decision are bit-identical.
        // ONE walk: collect canonical k-mers AND flag which came from a
        // stride-5 position. The old code walked twice (stride-5 query, then
        // full claim), recomputing pack25+rc25 both times.
        //
        // Exactness detail: the query loop stepped the CONTIG index by 5 and
        // skipped positions where pack25 failed, so "every 5th surviving
        // k-mer" is NOT the same set. The source position is therefore
        // recorded per entry, keeping `hit`/`tot` and the accept/reject
        // decision bit-identical.
```

What the code now says:

```
        // ONE walk per contig collects every canonical k-mer (all are claimed) and
        // flags the ones at stride-5 contig positions (those are queried). The source
        // position is recorded per entry because positions where pack25 fails are
        // skipped, so "every 5th collected k-mer" would not be the same set. `hit`,
        // `tot` and the accept/reject decision are those of a separate query walk
        // followed by a claim walk.
```

## gapped_indel_scan: dominant starts

Source: `include/caps_caller.h`, lines 839-842 before the documentation revision.

Original comment:

```
        // Take the two DOMINANT implied starts rather than demanding that every
        // seed agree. A single spurious repeat hit inside the same contig is
        // enough to break a strict all-agree test -- which is why the first
        // version of this scan found exactly 1 event in a whole window.
```

What the code now says:

```
        // Take the two DOMINANT implied starts rather than requiring every seed to
        // agree: one spurious repeat hit inside the same contig would otherwise
        // break the test.
```

## build_substrate timers

Source: `include/caps_caller.h`, lines 903-905 before the documentation revision.

Original comment:

```
    // [BS-SPLIT] temporary internal timers: build_substrate is 55 s of the
    // 128 s archive-path run across its two calls, and "collapse vs index vs
    // place" is the split that decides which lever is real.
```

What the code now says:

```
    // [BS-SPLIT] stderr timers for the collapse, index and placement steps. They
    // only print and do not affect results.
```

## Seed index over surviving contigs

Source: `include/caps_caller.h`, lines 949-979 before the documentation revision.

Original comment:

```
    // PACKED SORTED ARRAY, not a hash map, for the same reason the encoder's
    // MEM/mapping indices were converted earlier: a hash map here means
    // pointer-chasing + poor locality on the loop the caller's own comment
    // calls "the single largest phase" (738 s of the first full-chr20 run).
    // Sorted-by-key + binary search replaces `cidx.find(can)` with a single
    // std::lower_bound and a bounded linear scan over the matching run.
    //
    // ORDER-PRESERVING BY CONSTRUCTION, which matters: the original capped
    // each key's hits at 64 with `if (vec.size() < 64) push_back(...)`, so
    // which 64 survive for a repetitive k-mer depends on INSERTION ORDER
    // (contig, then position, ascending). A flat array in that same iteration
    // order, stable_sort'ed by key, reproduces the identical per-key ORDER
    // (stable_sort preserves relative order of equal keys) -- so truncating
    // each run to its first 64 entries afterward selects EXACTLY the same 64
    // the hash map would have, not merely a same-SIZE but different set.
    // ── SEED INDEX BUILT IN PLACE ───────────────────────────────────────────
    //
    // This was the run's high-water mark. The old shape was: fill `flat`
    // (16 B x ~106M = 1.70 GB), hand it to __gnu_parallel::stable_sort -- which
    // allocates a FULL second copy, measured +2.43 GB -- and then copy the
    // survivors out into three more arrays.
    //
    // Now there is ONE array. Entries are bucketed on the way in, each bucket
    // is sorted in place, and the 64-cap compacts in place. No sort temporary,
    // no second copy.
    //
    // The stable sort's guarantee was only ever "equal keys keep contig-then-
    // position order", because that is the order the fill produces. Saying that
    // outright in the comparator -- (mix, key, ci, pos) -- yields the identical
    // sequence without needing stability, which is what lets each bucket use
    // in-place std::sort. The 64-cap therefore still keeps exactly the same 64.
```

What the code now says:

```
    //
    // PACKED SORTED ARRAY, not a hash map: one binary search and a bounded scan per
    // lookup instead of pointer chasing, on the caller's largest phase.
    //
    // INVARIANT: each key keeps at most 64 hits, and they must be the FIRST 64 in
    // (contig, position) order, the order in which the hits are produced. The
    // comparator (mix, key, ci, pos) states that order explicitly, so an unstable
    // in-place sort inside each bucket gives the same sequence a stable sort would,
    // and the in-place 64-cap keeps exactly those 64. One array is used throughout,
    // with no sort temporary or second copy (this step set the run's peak memory).
```

## Gapped-scan index

Source: `include/caps_caller.h`, lines 1123-1126 before the documentation revision.

Original comment:

```
    // idx (the hash map) is now built ONLY for gapped_indel_scan, which is
    // opt-in (CAPS_GAPSCAN) and not part of any measured pipeline -- building
    // it unconditionally, every call, for a consumer that is off by default
    // was pure waste on every production run.
```

What the code now says:

```
    // idx (a hash map) serves only gapped_indel_scan, which is opt-in (CAPS_GAPSCAN),
    // so it is built only then.
```

## Read-level gapped indel channel (CAPS_GAPSCAN)

Source: `include/caps_caller.h`, lines 1194-1205 before the documentation revision.

Original comment:

```
    // REFUTED, and kept only behind an explicit opt-in (CAPS_GAPSCAN=1).
    // The idea: a read spanning an indel should show two different implied
    // contig starts. Measured: 1 event per window with a strict all-seeds-
    // agree test, 0 with a robust two-dominant-groups test.
    // The reason is architectural, not a bug. This substrate keeps the two
    // haplotypes as SEPARATE contigs -- which is precisely what makes the
    // bubble passes work -- so an indel-carrying read is placed on its OWN
    // haplotype's contig, where it matches with no gap at all. No read ever
    // needs a gapped placement, and the indel evidence lives in the contig
    // PAIR rather than in any read-vs-contig discrepancy.
    // A read-level gapped channel is therefore only meaningful for a caller
    // built on a single merged consensus; it cannot add recall here.
```

What the code now says:

```
    // OPT-IN (CAPS_GAPSCAN=1), not used by any reported configuration. This substrate
    // keeps the two haplotypes as separate contigs, so a read carrying an indel is
    // placed on its own haplotype's contig with no gap; the indel evidence is in the
    // contig PAIR, not in a read-versus-contig gap. A read-level gapped channel can only
    // add recall to a caller built on one merged consensus. It found 0-1 events per
    // test window here.
```

## Per-phase timing

Source: `include/caps_caller.h`, lines 1218-1221 before the documentation revision.

Original comment:

```
    // TEMPORARY phase timing (2026-09-03), added to find where wall time
    // actually goes before optimizing blind -- OpenMP over contigs measured
    // SLOWER at both 400kb and 5Mb scale (+11%, +13%), so the assumption
    // that the per-contig loop dominates needs checking, not more tuning.
```

What the code now says:

```
    // Per-phase wall time and RSS, printed to stderr. Diagnostic only; results do not
    // depend on it.
```

## Per-phase RSS

Source: `include/caps_caller.h`, lines 1226-1229 before the documentation revision.

Original comment:

```
    // Real per-phase RSS, not another structure-by-structure guess. Reads
    // VmRSS straight from /proc/self/status (Linux-only, this project's own
    // server) -- no third-party dependency, no sampling thread, negligible
    // cost next to the phases themselves.
```

What the code now says:

```
    // Reads VmRSS from /proc/self/status (Linux only; elsewhere it reports -1).
```

## Read placements carried into Method B

Source: `include/caps_caller.h`, lines 1263-1277 before the documentation revision.

Original comment:

```
        // ── READ PLACEMENTS, CARRIED THROUGH INSTEAD OF DISCARDED ───────────
        // cd_in already holds, for every original read, which contig it sits
        // on and at what offset -- the encoder derives it from ppos[u], the
        // pseudogenome position it computed in order to COMPRESS. Method B
        // used to copy only `contigs` and drop the rest, so the caller threw
        // away the one piece of information k-mer counts cannot reconstruct:
        // where each molecule actually came from.
        //
        // This is the same information a linked de Bruijn graph stores.
        // McCortex builds it in a dedicated threading pass (20 GiB of links on
        // top of a 50 GiB graph, Turner et al. 2018); LueVari stores read
        // provenance as colors and spends succinct data structures on it. The
        // idea of using read provenance to resolve repeats is THEIRS and is
        // published -- what is different here is only that compression already
        // paid for it, so carrying it costs three vector copies.
```

What the code now says:

```
        // ── READ PLACEMENTS, CARRIED THROUGH ─────────────────────────────────
        // cd_in holds, for every original read, its contig and offset, derived by the
        // encoder from the pseudogenome position it computed in order to compress.
        // k-mer counts cannot reconstruct this provenance. Linked de Bruijn graphs
        // (McCortex, Turner et al. 2018) and colored graphs (LueVari) use the same
        // idea of read provenance to resolve repeats; here compression has already
        // computed it, so carrying it costs three vector copies.
```

## Read seed index (eBWT2SNP-style read-substring test)

Source: `include/caps_caller.h`, lines 1302-1336 before the documentation revision.

Original comment:

```
    // ── READ SEED INDEX -- DISABLED BY DEFAULT (2026-09-03), see below ──────
    // EBWT2SNP's precision mechanism, from its paper: every emitted fragment
    // of length 2k+1 must be an actual SUBSTRING of at least C real READS
    // (within Hamming distance 2). It reports 99.13% precision against
    // DiscoSNP++'s 77.80% -- BUT on SIMULATED chr22 data (29x). CORRECTED
    // 2026-09-04 after reading the primary source (docs/EBWT_LITERATURE_
    // CORRECTION.md): on the paper's own REAL chr1 data (43-47x) eBWT2SNP is
    // 66.62% precision and DiscoSNP++ is actually MORE precise at 74.57%. The
    // guarantee below was real and worth testing, but its justification here
    // was a simulated-data number, not a real one, and the real-data result
    // for the source method is BELOW what this project already measures.
    // Our previous test only asked whether individual 31-mers appear in the
    // read k-mer TABLE -- far weaker, because a CHIMERIC junction (two
    // paralogs joined at a shared anchor) can have every one of its 31-mers
    // present, each contributed by a different read, while no single read
    // contains the whole fragment. That is precisely our false-positive class.
    //
    // THE FLIP: this index was already recorded as near-vacuous for us
    // (working note INDEL_PRECISION_ROOT_CAUSE.md, "Group B -- read-support tests
    // are near-vacuous for us": every candidate we emit is read-supported BY
    // CONSTRUCTION, since it was built FROM reads in the first place -- the
    // test EBWT2SNP needs to reject chimeras built from a GRAPH has nothing
    // to reject here). That was documented, but the index was still built
    // every run because a config flag (CAPS_NO_READSUB) existed to skip it,
    // not because the default did. Verified fresh, not just trusted from the
    // doc: ran all 5 standard windows (r2,r3,r4,r5,na) with CAPS_NO_READSUB=1
    // and diffed every TP/FP/FN against the enabled run -- BYTE-IDENTICAL on
    // all 5, SNV and indel both. So this is not a shrink like ridx's sorted-
    // array rewrite below (kept for the record, now dead code under the
    // opt-in flag) -- it is a full elimination of the largest remaining
    // single structure (~900 MB at 971K reads, projected multi-GB at full
    // chr20), for a change already measured to move nothing.
    // CAPS_FORCE_READSUB=1 restores the old default for re-measurement if
    // this ever needs revisiting (e.g. a future substrate change makes
    // candidates less strictly read-derived than they currently are).
```

What the code now says:

```
    // ── READ SEED INDEX -- DISABLED BY DEFAULT ─────────────────────────────
    // Read-substring support test after eBWT2SNP: an emitted fragment must occur
    // as a substring of at least C reads (Hamming distance <= 2), which rejects
    // chimeric junctions whose k-mers each occur in some read while no single read
    // contains the whole fragment. Here every candidate is already built from reads,
    // so the test has nothing to reject: on the five chr20 tuning windows the calls
    // were byte-identical with and without it, and the index is the largest single
    // structure. It is therefore off by default. CAPS_FORCE_READSUB=1 builds it.
    // Literature figures and the measurement: docs/extras/ENGINEERING_NOTES_CALLER.md.
```

## Chunked k-mer counting

Source: `include/caps_caller.h`, lines 1441-1456 before the documentation revision.

Original comment:

```
    // CHUNKED COUNTING (2026-09-03) -- the properly-done version of the
    // sorted-array idea that was tried and REVERTED once already (see git
    // history / prior comment here): a single-shot "collect everything, then
    // sort" spikes to ~11.9 GB transient at full chr20 scale, which can be
    // WORSE at peak than the hash map it replaces. This version bounds that
    // transient to one batch: collect and RLE-compress BATCH_KMERS raw
    // k-mers at a time into a small sorted run, discard the raw batch, then
    // k-way merge the runs (summing counts for a key that appears in more
    // than one run). Peak extra memory is O(batch size + distinct k-mers),
    // never O(total k-mer occurrences) -- the number that made the earlier
    // attempt fail. Output (the (kmer,count) set) is identical either way.
    // Batch size is DIVIDED BY THREAD COUNT so aggregate batch memory is the
    // same as the serial version. Each thread reserving the full 32M entries
    // (268 MB) would cost 268 MB x nthreads -- 3.2 GB on 12 cores -- turning
    // the parallel speedup into a RAM regression. Smaller batches simply mean
    // more, shorter runs, which the k-way merge handles natively.
```

What the code now says:

```
    // CHUNKED COUNTING. Collecting every k-mer occurrence and sorting once would
    // need memory proportional to total occurrences (about 12 GB at full chr20).
    // Instead BATCH_KMERS occurrences at a time are sorted and run-length counted
    // into a run, the raw batch is discarded, and the runs are k-way merged with
    // counts summed per key. Peak extra memory is O(batch + distinct k-mers) and the
    // (kmer, count) set is the same. The batch size is divided by the thread count so
    // total batch memory does not grow with threads.
```

## kc_H_build timers

Source: `include/caps_caller.h`, lines 1485-1487 before the documentation revision.

Original comment:

```
    // [KC-SPLIT] temporary internal timers: kc_H_build is the largest single
    // item left (11.9 s of 81.1 s) and "spill vs partition sort vs k-way
    // merge" is the split that decides which part of it to attack.
```

What the code now says:

```
    // [KC-SPLIT] stderr timers for the spill, partition sort and merge steps of
    // kc_H_build. Diagnostic only.
```

## Huge pages for kc

Source: `include/caps_caller.h`, lines 1494-1507 before the documentation revision.

Original comment:

```
    // ── HUGE PAGES FOR kc ───────────────────────────────────────────────────
    // kc reaches ~2.1 GB at full chr20 (140,719,632 x 16 B) and the bubble
    // traversal BINARY-SEARCHES it, so that walk is dTLB-bound. This box runs
    // /sys/kernel/mm/transparent_hugepage/enabled = [madvise]: a region gets
    // huge pages only if it asks, and nothing here ever asked.
    //
    // ORDER MATTERS. This project has already shipped an inert huge-page hint
    // by calling madvise AFTER the array was filled, which does nothing for
    // pages already faulted in. Every call site below sits immediately after a
    // reserve() and before the fill, so the hint lands on untouched memory.
    //
    // THERE ARE THREE reserve() SITES (in-RAM, spill/superkmer, and the k-way
    // merge) and the spill path is the one that runs at full scale. Hinting
    // only one of them is how this silently did nothing the first time.
```

What the code now says:

```
    // ── HUGE PAGES FOR kc ───────────────────────────────────────
    // kc is a multi-GB array that the bubble traversal binary-searches, so that
    // walk is TLB-bound. With transparent_hugepage=[madvise] a region gets huge
    // pages only if it asks. The hint must follow reserve() and precede the fill
    // (madvise does not affect pages already faulted in), and all three reserve()
    // sites (in-RAM, spill/superkmer, k-way merge) carry it.
```

## Singleton k-mer drop

Source: `include/caps_caller.h`, lines 1508-1515 before the documentation revision.

Original comment:

```
    // NOTE: dropping cnt==1 nodes from kc was implemented and MEASURED here,
    // and it is a NO-OP at full scale -- kc stayed at exactly 140,719,632 nodes.
    // The "60.4% of nodes are singletons" figure in the MINC comment below is a
    // WINDOW-scale measurement; at full chr20 (30x) the superkmer counting path
    // has already summed duplicate keys before any filter could see them, so
    // there are essentially no singletons left to drop. Reverted rather than
    // kept as dead weight. If this is ever revisited, measure the cnt==1
    // fraction AT THE SCALE YOU RUN, not from the window figure.
```

What the code now says:

```
    // Dropping cnt==1 nodes from kc does not shrink it at full chr20 scale: the
    // superkmer counting path has already summed duplicate keys, so almost no
    // singletons remain (the singleton share quoted in the MINC comment is a
    // window-scale figure). Measure the cnt==1 fraction at the scale being run
    // before revisiting this.
```

## Declared memory ceiling

Source: `include/caps_caller.h`, lines 1556-1573 before the documentation revision.

Original comment:

```
    // ── LAYER 3: DECLARED MEMORY CEILING ────────────────────────────────────
    // Standard systems practice (GATB/KMC both take a `-max-memory`), and the
    // reason it is here is a concrete failure, not tidiness: the spill was
    // opt-in through an environment variable, so a full-chr20 run launched
    // without it counted entirely in RAM and peaked at 28.8 GB instead of the
    // configured ~12 GB, at 190.9 s instead of 25.4 s for kc_H_build. The
    // k-mer set was identical and every accuracy number was unaffected -- the
    // run was simply 4x heavier than the tool is capable of, because a
    // performance-critical default lived in the caller's environment.
    //
    // A budget the tool DECLARES cannot be forgotten. The ceiling defaults to
    // 60% of MemAvailable (measured at run time, not assumed), the projected
    // in-RAM cost of counting is estimated from the input, and the spill turns
    // itself on when the estimate exceeds the ceiling.
    //
    // This changes no output. Spilling alters WHERE k-mers are counted, never
    // which ones exist -- the merge sums counts per key across runs either way
    // -- and that identity is what the k-mer-count gate verifies.
```

What the code now says:

```
    // ── LAYER 3: DECLARED MEMORY CEILING ──────────────────────────────────
    // As in GATB and KMC (`-max-memory`), the counter declares a memory budget,
    // estimates the in-RAM counting cost from the input, and turns the disk spill on
    // by itself when the estimate exceeds the budget. Without this a full-chr20 run
    // started without the spill variable counted entirely in RAM at about 4x the
    // memory. Spilling changes where k-mers are counted, never which ones exist:
    // the merge sums counts per key either way, which the k-mer-count gate checks.
```

## Absolute default memory budget

Source: `include/caps_caller.h`, lines 1585-1612 before the documentation revision.

Original comment:

```
    // ── AN ABSOLUTE DEFAULT BUDGET, NOT A FRACTION OF WHAT IS FREE ──────────
    // This defaulted to 60% of MemAvailable, which silently DISABLES the spill
    // on a large machine: on a 76 GB box the ceiling came out at 45,548 MB, the
    // 22,695 MB projection fit under it, counting ran in RAM, and the caller
    // peaked at 28.3 GB taking 191.7 s in kc_H_build -- against 25.4 s and ~7 GB
    // for the spill path this file already implements. The bigger the machine,
    // the worse the tool behaved.
    //
    // GATB does the opposite, and it is why DiscoSNP++ holds 3.29 GB on this
    // same box and dataset (measured 2026-09-06, 75.72 s / 3.29 GB, matching
    // its published 76.5 s / 3.45 GB). ConfigurationAlgorithm.cpp:334-343:
    //
    //     if (_max_memory == 0) _max_memory = 5000;                  // MB
    //     if (_max_memory > (system_mem*2)/3) _max_memory = (system_mem*2)/3;
    //
    // An absolute default FIRST, with system memory used only as a CAP. The
    // budget is a choice the tool makes, not a resource it discovers -- which
    // is the whole point of the "declared ceiling" this file already argues
    // for a few lines above, applied to its own default.
    //
    // 5000 MB is GATB's value, adopted deliberately rather than invented: it is
    // the number their partitioning has been tuned against for years, and our
    // spill derives partition count from the ceiling the same way. Explicit
    // CAPS_MAXRAM_MB still wins over it.
    //
    // THIS CHANGES NO OUTPUT. Spilling alters WHERE k-mers are counted, never
    // which ones exist -- the merge sums counts per key either way, and the
    // k-mer-count gate verifies exactly that identity.
```

What the code now says:

```
    // ── AN ABSOLUTE DEFAULT BUDGET, NOT A FRACTION OF WHAT IS FREE ──────────
    // A budget of 60% of MemAvailable disabled the spill on large machines (the
    // projection fit, so counting ran in RAM at several times the memory). The
    // default is therefore an absolute 5000 MB, capped by system memory, which is
    // GATB's rule (ConfigurationAlgorithm.cpp:334-343):
    //
    //     if (_max_memory == 0) _max_memory = 5000;                  // MB
    //     if (_max_memory > (system_mem*2)/3) _max_memory = (system_mem*2)/3;
    //
    // An explicit CAPS_MAXRAM_MB overrides it. This changes no output, for the
    // reason given above.
```

## Superkmer spill is the default

Source: `include/caps_caller.h`, lines 1669-1685 before the documentation revision.

Original comment:

```
    // SUPERKMER IS THE DEFAULT WHENEVER WE SPILL, AND THAT IS A CORRECTNESS
    // REQUIREMENT, NOT A PERFORMANCE ONE.
    //
    // It used to be opt-in on top of the spill, which left the key-only spill
    // reachable as the default spilling path. Measured on HG002 r2 with the
    // ceiling forced low enough to trigger a spill:
    //
    //     in-RAM             kc nodes = 1,063,607   SNV F1 = 0.886
    //     spill, key-only    kc nodes =   741,350   SNV F1 = 0.000
    //     spill, superkmer   kc nodes = 1,063,607   SNV F1 = 0.886
    //
    // The key-only path silently drops 30% of the k-mer set. docs/ called it
    // "built, unverified"; it is in fact broken, and it was one environment
    // variable away from being used. Superkmer reproduces the in-RAM multiset
    // exactly, so it is what a spill turns on.
    //
    // CAPS_KC_NOSUPERKMER=1 restores the old path for re-measurement only.
```

What the code now says:

```
    // SUPERKMER IS THE DEFAULT WHENEVER THE COUNTER SPILLS, FOR CORRECTNESS.
    // The key-only spill path loses part of the k-mer set (HG002 tuning window,
    // spill forced: 741,350 kc nodes against 1,063,607 in RAM, SNV F1 0.000 against
    // 0.886). The superkmer path reproduces the in-RAM multiset exactly.
    // CAPS_KC_NOSUPERKMER=1 selects the key-only path and exists only to reproduce
    // that measurement; do not use it for results.
```

## Frequency-ranked minimizers

Source: `include/caps_caller.h`, lines 1689-1728 before the documentation revision.

Original comment:

```
    // Standard practice in KMC2/GATB, adopted here with attribution; the
    // mechanism is theirs, the derivation from our own data is ours.
    //
    // The minimizer was the LEXICOGRAPHICALLY smallest 10-mer in each window.
    // That is the worst possible choice for low-complexity DNA: poly-A packs to
    // the value 0, so it wins every window it touches, and the most frequent
    // m-mers in a genome end up selected constantly. Frequent minimizers make
    // SHORT runs -- the minimizer changes as soon as the window slides past the
    // repeat -- and each run pays a length byte plus k-1 bases of overlap. Long
    // runs are what make a superkmer cheap.
    //
    // Ranking by measured frequency inverts it: a RARE m-mer stays the window
    // minimum for longer, so runs lengthen and bytes/k-mer falls. The ordering
    // key is (frequency, m-mer), so it is derived from the input rather than
    // fixed, and the m-mer breaks ties to keep the result deterministic.
    //
    // The table is 4^10 counts = 4 MB, filled from a strided SAMPLE of reads.
    // A sample is sufficient because only the ORDER matters, not the counts,
    // and the stride is deterministic so runs are reproducible.
    //
    // THIS CANNOT CHANGE THE K-MER SET. A superkmer is a container: whichever
    // minimizer is chosen, the span still yields exactly its L k-mers and the
    // merge still sums counts per key. Grouping changes, membership does not --
    // which is precisely what the k-mer-count gate checks.
    // MEASURED AND REJECTED AS A VOLUME LEVER -- OPT-IN (CAPS_KC_FREQMIN=1).
    // On HG002 r2 with an identical k-mer set (1,063,607 both ways) and
    // identical SNV F1 (0.886 both ways):
    //     lexicographic     10,832,693 spill bytes   10.185 B/k-mer
    //     frequency-ranked  12,398,016 spill bytes   11.657 B/k-mer  (+14%)
    //
    // The reasoning that motivated it was backwards. Lexicographic ordering
    // makes poly-A win, and because it KEEPS winning across a repetitive
    // stretch the minimizer run is LONG -- which is exactly what makes a
    // superkmer cheap. Frequency ranking elects whichever rare m-mer is in the
    // window, and rare m-mers are transient, so the minimizer changes more
    // often and runs get SHORTER. KMC2/GATB adopt frequency ranking for
    // PARTITION BALANCE, not for volume; the "~1.16 -> ~0.95 B/k-mer"
    // expectation was carried over without checking which property it buys.
    // Kept, off, because partition balance may still matter if the merge is
    // ever made memory-bound per partition.
```

What the code now says:

```
    // Frequency-ranked minimizers, as in KMC2 and GATB: order m-mers by sampled
    // frequency (ties broken by the m-mer, so the order is deterministic) instead of
    // lexicographically.
    //
    // THIS CANNOT CHANGE THE K-MER SET: a superkmer is a container, and the merge
    // sums counts per key whatever the grouping.
    //
    // OPT-IN (CAPS_KC_FREQMIN=1), OFF by default. It does not reduce spill volume:
    // on the HG002 tuning window it wrote 14% more bytes per k-mer, because rare
    // minimizers change more often and give shorter runs. KMC2 and GATB use it for
    // partition balance, which would matter only if the merge became memory-bound
    // per partition. Measurements: docs/extras/ENGINEERING_NOTES_CALLER.md.
```

## Spill write failure handling

Source: `include/caps_caller.h`, lines 1760-1766 before the documentation revision.

Original comment:

```
        // A full spill disk previously failed silently: every `if (f) {...}`
        // below has no `else`, so a fopen("ab") failure (disk full, quota,
        // permissions) just dropped that batch's k-mer occurrences with no
        // diagnostic -- for the one input-shape (spill engaged) this project's
        // own exactness claims most need to hold. Recorded here instead of
        // aborting mid-parallel-region (which would need per-thread unwind
        // logic for no benefit); checked once, right after the region joins.
```

What the code now says:

```
        // A failed spill write (disk full, quota, permissions) is recorded in a flag
        // rather than aborting inside the parallel region, and checked once after the
        // region joins, so a short spill is an error instead of silently dropped
        // k-mer occurrences.
```

## Spill partitioned by key range

Source: `include/caps_caller.h`, lines 1785-1802 before the documentation revision.

Original comment:

```
                // The first version of this spill wrote each batch as one
                // sorted RUN. That was the wrong shape: every run spans the
                // whole key space, so the k-way merge has to hold all runs
                // open simultaneously and peak stays proportional to the total.
                //
                // GATB does not do that. SortingCountAlgorithm.cpp:804 routes
                // each superkmer to a partition chosen by its MINIMIZER --
                //     size_t p = _repartition(superKmer.minimizer);
                // so every partition owns a DISJOINT slice of key space and
                // partitions are processed one at a time, with no merge across
                // them at all. Their peak is one partition; ours was the sum.
                //
                // We can do it more simply than they do: partition on the TOP
                // BITS of the canonical k-mer. Then partition order IS key
                // order, so concatenating the sorted partitions yields a
                // globally sorted kc with no merge step whatsoever. (GATB needs
                // minimizers because superkmers must keep adjacent k-mers
                // together; we store plain k-mers, so we do not.)
```

What the code now says:

```
                // Partition on the TOP BITS of the canonical k-mer, so every partition owns
                // a disjoint key range and partitions are processed one at a time (GATB
                // partitions by minimizer for the same reason, SortingCountAlgorithm.cpp:804).
                // Partition order is key order, so concatenating the sorted partitions gives a
                // globally sorted kc with no merge. (GATB needs minimizers because superkmers
                // must keep adjacent k-mers together; kc stores plain k-mers.)
```

## Spill only the key

Source: `include/caps_caller.h`, lines 1804-1812 before the documentation revision.

Original comment:

```
                    // SPILL ONLY THE KEY, NOT THE COUNT.
                    // A (kmer,count) record is 12 bytes; the count is
                    // recomputed when the partition is read back and sorted,
                    // so writing it is pure waste. Measured: the count made
                    // the full-chr20 spill 25 GB, which exhausted a 233 GB
                    // disk. Dropping it is a 33% cut for free. (GATB gets far
                    // more -- ~0.95 B/k-mer -- by packing SUPERKMERS, which
                    // needs the whole minimizer chain; this is the part
                    // available without it.)
```

What the code now says:

```
                    // SPILL ONLY THE KEY: the count is recovered when the partition is read back
                    // and sorted (one entry per occurrence), so writing it would waste a third
                    // of the spill.
```

## Rolling 31-mer window (counting)

Source: `include/caps_caller.h`, lines 1837-1842 before the documentation revision.

Original comment:

```
            // The previous loop called pack31() and canon31() at every
            // position, and BOTH are O(k): pack31 loops 31 times to build the
            // k-mer, rc31 loops 31 times to reverse-complement it. At full
            // chr20 that is ~1.4 BILLION k-mer positions x ~62 operations
            // = ~87 billion operations, and it is the reason this phase costs
            // 72.8 s -- almost none of that work is necessary.
```

What the code now says:

```
            // Packing and reverse-complementing each k-mer from scratch is O(k) per
            // position, at about 1.4 billion positions for full chr20.
```

## Superkmer spill format

Source: `include/caps_caller.h`, lines 1856-1859 before the documentation revision.

Original comment:

```
            // before wiring: 1.161 bytes/k-mer against 12 for (kmer,count)
            // records, with the k-mer multiset proven identical.
            // This is what makes the spill sustainable -- the record format hit
            // 25 GB at full chr20 and exhausted the disk twice.
```

What the code now says:

```
            // The format was verified at 1.161 bytes per k-mer, against 12 for
            // (kmer, count) records, with an identical k-mer multiset.
```

## Minimizer window after N

Source: `include/caps_caller.h`, lines 1898-1900 before the documentation revision.

Original comment:

```
                        // which silently dropped 93 k-mers on the r2 window
                        // (1,063,514 vs 1,063,607) because those superkmers were
                        // then rejected wholesale.
```

What the code now says:

```
                        // which would drop k-mers next to an N (their superkmers are
                        // rejected whole).
```

## Heap merge, two passes

Source: `include/caps_caller.h`, lines 1995-2015 before the documentation revision.

Original comment:

```
    // HEAP MERGE, TWO PASSES (2026-09-03). The comment above ("a heap would be
    // faster but is not needed at this run count") was correct about the run
    // count and wrong about the cost, because the cost is driven by the OUTPUT
    // length, not the run count. The previous version did two linear scans over
    // all runs per output element: at full chr20 that is roughly
    //     ~400M distinct 31-mers x 47 runs x 2 scans = ~37 BILLION comparisons,
    // which is essentially all of this phase's 166.7 s. A min-heap over run
    // heads makes it O(output x log runs) instead.
    //
    // Pass 1 counts distinct keys so pass 2 can allocate `kc` exactly once.
    // The tempting one-pass alternative -- reserve the summed run lengths --
    // is an exact upper bound but a useless one: every run carries its own copy
    // of each common key, so that sum is ~47 x 32M entries (~18 GB reserved) to
    // hold ~400M (~4.8 GB) of real output. Two heap passes cost ~2x the merge
    // but the merge is now ~8-17x cheaper, so this is still several times
    // faster than the code it replaces AND removes the ~14.4 GB transient the
    // un-reserved doubling caused at the final reallocation.
    //
    // Byte-identical by construction: a k-way merge is a pure function of its
    // runs, and both passes visit keys in the same ascending order and sum the
    // same per-run counts.
```

What the code now says:

```
    // HEAP MERGE, TWO PASSES. The cost is driven by the output length, not the run
    // count, so a linear scan of all run heads per output key is too slow at full
    // chr20; a min-heap over run heads makes the merge O(output x log runs).
    // Pass 1 counts distinct keys so that pass 2 allocates kc exactly once (reserving
    // the summed run lengths would over-reserve several-fold, because every run holds
    // its own copy of each common key). Both passes visit keys in the same ascending
    // order and sum the same per-run counts, so the output is unchanged.
```

## Parallel superkmer read-back

Source: `include/caps_caller.h`, lines 2033-2038 before the documentation revision.

Original comment:

```
        // This was once reverted after the parallel version produced 2,273
        // extra k-mers -- but that was a MISDIAGNOSIS. The extra k-mers came
        // from 2-bit read packing mapping non-ACGT bases to 'A', which made
        // k-mers spanning an N look valid. The parallelism was never at fault;
        // reverting it removed a correct optimisation. Restored once the N bug
        // was fixed and kc identity (1,063,607) re-verified.
```

What the code now says:

```
        //
        // The parallel read-back reproduces the serial k-mer set exactly (kc node
        // count checked on the HG002 tuning window). An earlier discrepancy was caused
        // by 2-bit packing mapping non-ACGT bases to 'A', not by the parallelism.
```

## Parallel merge by key range: motivation

Source: `include/caps_caller.h`, lines 2041-2043 before the documentation revision.

Original comment:

```
        // The k-way merge that used to follow this loop was 8.21 s of an
        // 11.88 s kc_H_build: one serial priority_queue over 256 files with a
        // 12-byte fread per record, ~100M times.
```

What the code now says:

```
        // A single serial priority-queue merge over all partition files was most of
        // kc_H_build's time.
```

## Two-pass exact-allocation merge

Source: `include/caps_caller.h`, lines 2149-2158 before the documentation revision.

Original comment:

```
        // The first version of this merge deduped each bin into bout[b] and
        // then concatenated all of bout into kc, so the merged result existed
        // TWICE. That is invisible on the indel path, where the peak is
        // elsewhere, but on the graph-only path (the published Claim 2 config)
        // kc IS the peak, and it cost +1.25 GB there -- a regression on one
        // path bought with a win on another.
        //
        // Allocating kc at the pre-dedup bound instead was also tried and is
        // worse: resize() down does not release capacity, so the waste is held
        // for the rest of the run.
```

What the code now says:

```
        // On the graph-only path (the published Claim 2 configuration) kc sets the peak
        // memory, so the merged result must never be resident twice, and allocating kc
        // at the pre-dedup bound would hold the excess for the rest of the run (resize()
        // down does not release capacity).
```

## Reserve from first partition

Source: `include/caps_caller.h`, lines 2265-2269 before the documentation revision.

Original comment:

```
            // kc.push_back was growing by doubling, and at full chr20 the
            // final array is ~1.7 GB, so the last reallocation held old+new
            // ~3.4 GB -- a pure transient with no phase attached to it
            // (measured at 2M reads: RSS flat at 2772 MB, PEAK 4127 MB).
            // The two-pass merge this replaced used an exact reserve; the
```

What the code now says:

```
            // Growing kc by doubling would hold old and new arrays together at the
            // last reallocation (several GB at full chr20).
```

## H estimator binning

Source: `include/caps_caller.h`, lines 2189-2199 before the documentation revision.

Original comment:

```
        // H ESTIMATOR BUG (found 2026-09-02): cnt_max was the MAXIMUM k-mer
        // count, which is set by repeats and pinned at the 5000 cap, so
        // bw = 5000/200 = 25. The true single-copy peak (~11 at 27x read depth:
        // 27 * (148-31+1)/148 / 2) then lands in bucket 0, which the valley
        // logic excludes, forcing the peak search to start at bucket 3 and
        // yielding H ~= 87 -- about 8x too high. With H that wrong, the depth
        // guard (d > DHI*H*1.25) and the k-mer sanity gate (> KHI*H) can never
        // fire, so two of the frozen filters were silently inert.
        // Fix: bin at FULL RESOLUTION (bw = 1) over a bounded range. Repeat
        // k-mers above the range are irrelevant to locating the single-copy
        // peak, and 2000 covers any realistic per-haplotype depth.
```

What the code now says:

```
        // H ESTIMATION bins the count histogram at FULL RESOLUTION (bw = 1) over a
        // bounded range. Bin width derived from the maximum count (set by repeats) would
        // put the single-copy peak in bucket 0, which the valley logic excludes, and
        // overestimate H several-fold; the depth guard (d > DHI*H*1.25) and the k-mer
        // sanity gate (> KHI*H) would then never fire. Counts above the range are repeats
        // and do not affect the single-copy peak; 2000 covers any realistic
        // per-haplotype depth.
```

## Rejected histogram-bimodality ploidy test

Source: `include/caps_caller.h`, lines 2218-2225 before the documentation revision.

Original comment:

```
            // (A histogram-bimodality ploidy test was tried here and REMOVED:
            // it did not discriminate. Measured ratios: HG002 0.125 (diploid),
            // SARS-CoV-2 0.331 (haploid), E. coli 0.088 (haploid) -- the
            // haploids fell both above AND below the diploid, so no threshold
            // separates them. The reasoning was wrong: at 30x with k=31 the
            // located peak sits near the HOMOZYGOUS depth, so "mass at 2H"
            // measures repeats, not ploidy. Ploidy is detected instead by the
            // pair test below, the mechanism already validated 5/5 by HETSCAN.)
```

What the code now says:

```
            // (Histogram bimodality does not identify ploidy here: at 30x with k=31 the
            // located peak sits near the homozygous depth, so mass at 2H measures repeats.
            // The measured ratios for haploid samples fell both above and below HG002.
            // Ploidy is detected by the k-mer pair test below.)
```

## kcount without a temporary string

Source: `include/caps_caller.h`, lines 2234-2239 before the documentation revision.

Original comment:

```
    // Same thing without materialising the 31-mer. The junction scans in the
    // indel emit path call this once per offset as
    // `kcount(hap.substr(q, 31))`, so every probe allocated a fresh
    // std::string purely to be read 31 times and thrown away -- inside the
    // section that measured 1,018 s of an 1,093 s stage. `pack31` reads a raw
    // pointer, so the copy was never needed. Identical result by construction.
```

What the code now says:

```
    // Same as kcount(std::string) without building a 31-character string per probe;
    // the junction scans in the indel emit path call it once per offset.
```

## Contig-batched pileup

Source: `include/caps_caller.h`, lines 2247-2261 before the documentation revision.

Original comment:

```
    // CONTIG-BATCHED (2026-09-03). Previously col/C/FLmaj/FLmin/majoff/minoff
    // were built ONCE over every read in the whole chromosome and held
    // resident simultaneously -- measured as the dominant RAM sink on a real
    // full chr20 run (tens of GB). Nothing in this stretch crosses a contig
    // boundary: colkey always carries cid, and every lookup into these five
    // structures is keyed by a (cid,pos) that only ever matches the contig
    // being scanned. So the identical computation is re-scoped to run one
    // contig at a time, with all five structures declared fresh inside the
    // loop and freed before the next contig starts -- peak RAM becomes
    // O(largest single contig's pileup), not O(whole chromosome).
    // Output is byte-identical to the old global-pass version: the original
    // did one `std::sort(kept)` over (cid,pos) pairs from every contig before
    // emitting; emitting per-contig in increasing cid order and concatenating
    // produces the exact same total order, since cid is the primary sort key
    // and each contig's own `kept` is still locally sorted by pos.
```

What the code now says:

```
    // CONTIG-BATCHED. Nothing in this stretch crosses a contig boundary: every key
    // into the pileup structures carries the contig id. So the pileup is built one
    // contig at a time and freed before the next, and peak memory is that of the
    // largest single contig's pileup rather than the whole chromosome's (tens of GB
    // at full chr20). Output order is unchanged: contigs are emitted in increasing
    // cid, cid is the primary sort key, and each contig's `kept` is sorted by pos.
```

## H override hoisted

Source: `include/caps_caller.h`, lines 2299-2305 before the documentation revision.

Original comment:

```
    // Hoisted out of the per-contig loop (2026-09-03): this override does not
    // depend on cid0 at all -- it re-read the same env var on every contig
    // for no reason, and once the loop below runs on multiple threads,
    // writing the shared `H` from inside it would be a data race. Reading it
    // once here, before any thread starts, is both a correctness fix for the
    // parallel version and a no-op for the sequential one (same value either
    // way).
```

What the code now says:

```
    // Read once, before any thread starts: the override does not depend on the
    // contig, and writing the shared `H` inside the parallel loop would be a data
    // race.
```

## What is and is not batched

Source: `include/caps_caller.h`, lines 2308-2327 before the documentation revision.

Original comment:

```
    // WHAT IS AND IS NOT BATCHED, corrected 2026-09-03 after a review pass
    // caught an over-reach in the first version of this refactor.
    // Only `col` (the pileup: every read base at every covered position) and
    // `recs` are large -- measured shapes put col in the tens of GB at full
    // chr20 scale, recs at a few GB. `C` (candidate columns) and the flank
    // maps are SMALL: 585 candidates on a 400 kb window scales to roughly
    // 1e5 genome-wide, a few hundred MB with the flank contexts included.
    // The first version of this refactor scoped C and the flank maps per
    // contig too, which bought nothing and silently broke `med_cand_depth`
    // -- the depth guard's median, which the original computes over ALL
    // candidates chromosome-wide, became a per-contig median. That guard is
    // inert on the default path (it is gated on !collapse_ran) so no
    // recorded result moved, but it is exactly the kind of silent
    // methodology change that must not happen unmeasured.
    // So: build col/recs per contig and free them there, accumulate C and
    // the flank maps into thread-local buffers, then merge and run the
    // filter/emit stages EXACTLY as the original wrote them, over the
    // complete global C. Keys are (cid,pos) and each contig is handled by
    // exactly one thread, so no two threads ever produce the same key and
    // merging is a pure move with no count reconciliation.
```

What the code now says:

```
    // WHAT IS AND IS NOT BATCHED. Only `col` (every read base at every covered
    // position) and `recs` are large, so only they are built and freed per contig.
    // The candidate columns `C` and the flank maps are small and are accumulated into
    // thread-local buffers, then merged, and the filter and emit stages run over the
    // complete global C exactly as in the unbatched form. This matters for
    // `med_cand_depth`, the depth guard's median, which is defined over ALL candidates
    // chromosome-wide and must not become a per-contig median. Keys are (cid,pos) and
    // each contig is handled by one thread, so the merge is a move with no count
    // reconciliation.
```

## DBG bubble channel on kc

Source: `include/caps_caller.h`, lines 2382-2410 before the documentation revision.

Original comment:

```
    // THE REALIZATION. A de Bruijn graph is exactly two things: a set of
    // k-mers, and the implicit edge X->Y whenever Y = (X shifted left one base
    // + one new base) is also in the set. `kc` above IS that set -- every
    // canonical 31-mer of every read, with its count -- and `kc_find` IS the
    // membership test. We spend 166.7 s building a full de Bruijn graph and
    // then use it to compute ONE SCALAR (H, the coverage threshold below).
    //
    // So DiscoSNP++'s bubble finder does not need a new index, a new pass over
    // the reads, or a second assembly. Every primitive it uses maps onto a
    // kc_find:
    //
    //   their graph.successors(node)        -> 4 kc_find probes, keep cnt >= minc
    //   their branching node (>= 2 succ)    -> same test
    //   their graph.successors(node1,node2) -> for b in ACGT, BOTH extensions
    //                                          present (their getNodesCouple in
    //                                          GATB Graph.cpp:1602-1640 requires
    //                                          the SAME nt on both paths)
    //   their closure (succ.first==.second) -> node1 == node2
    //
    // This is a faithful port of Bubble.cpp:300-347 (start: branching node,
    // every successor PAIR) and Bubble.cpp:509-527 (expand: lockstep walk,
    // close on node equality), not an approximation of it. Two earlier attempts
    // failed precisely because they lacked this closure test:
    // chaining ties are exact-match and so blind to heterozygosity
    // (the GRAPH HARVEST block in docs/extras/ENGINEERING_NOTES_ENCODER.md), and near-miss pairs have no graph to
    // walk against so they fall back on a string-flank test that was already
    // measured net-negative (flank_match_tol).
    //
    // Opt-in while it is being measured: CAPS_DBG=1.
```

What the code now says:

```
    // A de Bruijn graph is a set of k-mers plus the implicit edge X->Y whenever Y
    // (X shifted left one base plus one new base) is also in the set. `kc` is that
    // set (every canonical 31-mer of the reads, with its count) and `kc_find` is the
    // membership test, so DiscoSNP++'s bubble finder needs no new index, no new pass
    // over the reads and no second assembly. Each primitive maps onto kc_find:
    //
    //   their graph.successors(node)        -> 4 kc_find probes, keep cnt >= minc
    //   their branching node (>= 2 succ)    -> same test
    //   their graph.successors(node1,node2) -> for b in ACGT, BOTH extensions
    //                                          present (their getNodesCouple in
    //                                          GATB Graph.cpp:1602-1640 requires
    //                                          the SAME nt on both paths)
    //   their closure (succ.first==.second) -> node1 == node2
    //
    // The traversal reimplements DiscoSNP++ Bubble.cpp:300-347 (start: a branching
    // node, every successor PAIR) and Bubble.cpp:509-527 (expand: lockstep walk,
    // close on node equality). The closure test is the essential part: approaches
    // without it (chaining ties, which are exact-match; string-flank tests on
    // near-miss pairs) cannot see heterozygosity or lowered F1.
    //
    // Enabled by CAPS_DBG=1, which `capsule_decode call` sets.
```

## Per-path support from kc

Source: `include/caps_caller.h`, lines 2413-2417 before the documentation revision.

Original comment:

```
    // every read back onto every bubble. We do not need that pass: a canonical
    // k-mer's count in kc IS the number of reads containing it, so walking the
    // path and taking the minimum gives per-path read support directly. The
    // minimum rather than the mean because a bubble is only as well supported
    // as its weakest link -- one unsupported k-mer means no read spans it.
```

What the code now says:

```
    // every read back onto every bubble. Here that pass is not needed: a canonical
    // k-mer's count in kc is the number of reads containing it, so the minimum along
    // the path is per-path read support. The minimum rather than the mean, because
    // one unsupported k-mer means no read spans that part of the path.
```

## MINC floor

Source: `include/caps_caller.h`, lines 2427-2453 before the documentation revision.

Original comment:

```
        // This was previously DERIVED from the k-mer histogram valley (capped
        // at H/3), which sounds principled and is optimising the wrong thing.
        // The valley is the right threshold only when the floor is your ONLY
        // filter. It is not: every bubble is validated downstream against real
        // reads (the coherence pass). In a two-stage design -- permissive
        // generation, then evidence-based filtering -- ALL precision must come
        // from the evidence stage, because precision bought upstream costs
        // recall that no downstream stage can recover.
        //
        // MEASURED, on the window's k-mer count distribution:
        //     cnt=1      60.4% of nodes   (singleton errors)
        //     cnt=2-9     3.9%
        //     cnt>=10    35.6%            (real signal)
        // Nodes surviving:  >=2: 39.6%   >=3: 38.2%   >=4: 37.8%   >=5: 37.6%
        // So raising the floor from 2 to 5 removes 5% of nodes -- no
        // meaningful compute saving -- while deleting exactly the 2-4 count
        // band where LOW-COVERAGE HETEROZYGOUS ALLELES live. The floor's only
        // real job is discarding the 60% singletons.
        //
        // This also explains a failure of the derived version: it moved UP
        // (4 -> 5) at full chr20, precisely when more data makes a lower floor
        // SAFER, because absolute error counts grow and drag the valley with
        // them. That is the wrong direction for the wrong reason.
        //
        // DiscoSNP++ runs -c 3 for the same reason: a permissive graph, with
        // kissreads2 supplying precision. 2 is the smallest defensible value
        // ("seen more than once"), so the floor now does only what it should.
```

What the code now says:

```
        // The floor is not the precision filter: every bubble is validated downstream
        // against the reads (the coherence pass), and in a two-stage design (permissive
        // generation, then evidence-based filtering) precision bought upstream costs recall
        // that no later stage can recover. On the tuning window, 60.4% of graph nodes have
        // count 1 (singleton errors), 3.9% have 2-9 and 35.6% have >= 10. Raising the
        // floor from 2 to 5 removes only 5% of nodes while deleting the 2-4 band where
        // low-coverage heterozygous alleles are. A floor derived from the histogram valley
        // also rises with depth, the wrong direction. MINC = 2 ("seen more than once")
        // discards the singletons and nothing else. DiscoSNP++ uses a permissive -c 3 for
        // the same reason, with kissreads2 supplying precision.
```

## Interleaved successor search: justification

Source: `include/caps_caller.h`, lines 2476-2483 before the documentation revision.

Original comment:

```
        // This is worth doing only because the loop is LATENCY-bound rather
        // than bandwidth-bound -- a distinction this project has been burned by
        // before (a filter that won in a latency-bound loop lost in a
        // bandwidth-bound one). Verified two ways before writing it: a model of
        // the access pattern predicts ~113 s of pure latency against an
        // observed 53.5 s traversal (same order, i.e. latency-dominated), and a
        // standalone benchmark on a table of the same shape and size measured
        // 2.123 s sequential vs 1.393 s interleaved -- 1.52x.
```

What the code now says:

```
        // This helps only because the loop is latency-bound rather than bandwidth-bound:
        // a model of the access pattern and the observed traversal time agree on that,
        // and a standalone benchmark on a table of the same shape ran 1.52x faster
        // interleaved.
```

## Successor search seeded from the directory

Source: `include/caps_caller.h`, lines 2499-2504 before the documentation revision.

Original comment:

```
            // This 4-way interleaved binary search was seeded with the FULL
            // array on every call -- 140,719,632 entries at full chr20, so ~27
            // steps x 4 keys, and it runs twice for every solid node in the
            // traversal. It is the traversal's real successor lookup; kc_find
            // is not on this path, which is why giving kc_find a directory
            // moved the traversal only 11.6%.
```

What the code now says:

```
            // This 4-way interleaved search is the traversal's successor lookup (kc_find is
            // not on this path).
```

## Lockstep walk support accumulation

Source: `include/caps_caller.h`, lines 2544-2546 before the documentation revision.

Original comment:

```
        // Also accumulates the weakest k-mer count seen on each path (see the
        // DbgBubble comment): this is per-path read support, obtained during
        // the walk we are already doing rather than from a second pass.
```

What the code now says:

```
        // Also accumulates the weakest k-mer count seen on each path (see the
        // DbgBubble comment), so per-path read support comes from this walk rather than
        // from a second pass.
```

## Bounded recursion over ambiguous continuations

Source: `include/caps_caller.h`, lines 2549-2555 before the documentation revision.

Original comment:

```
        // The previous version bailed the instant the two paths had more than
        // one common successor ("if (++nc > 1) return false"). DiscoSNP++ does
        // NOT do this: expand() (Bubble.cpp:598) recurses over every successor
        // pair through expand_heart. So we were MORE permissive than them at
        // the bubble entrance (no checkBranching) and STRICTLY LESS permissive
        // during extension -- and the abort counters showed the cost: of 2,628
        // entrances only 462 closed, 1,856 died mid-walk.
```

What the code now says:

```
        // When the two paths have more than one common successor, the search recurses over
        // every successor pair, as DiscoSNP++ expand() does (Bubble.cpp:598, through
        // expand_heart). Stopping at the first ambiguity lost most bubble entrances
        // mid-walk.
```

## MAXPOLY

Source: `include/caps_caller.h`, lines 2565-2598 before the documentation revision.

Original comment:

```
        // MEASURED, NOT JUST ARGUED (2026-09-05). This was justified on
        // structural grounds alone, which is weaker than it needed to be --
        // DiscoSNP++ ships P=3, so "1 is right" is a claim that has to be
        // tested against their default. Swept over the tuning window and four
        // held-out windows, SNV F1:
        //     window      MAXPOLY=1   MAXPOLY=2   MAXPOLY=3 (their P=3)
        //     r2 (tune)     0.886       0.881       0.881
        //     r3            0.892       0.886       0.879
        //     na            0.958       0.958       0.956
        //     r4            0.894       0.895       0.897
        //     r5            0.937       0.923       0.916
        //     held-out mean 0.9203      0.9155      0.9120
        // Monotone, and 1 wins. Relaxing it does buy recall (+6 TP on r2) but
        // costs roughly double in precision (+12 FP), because a second
        // polymorphism inside one walk is far more often two error branches
        // than two real het sites. Our narrower setting beats their default on
        // this pipeline; it is now evidence rather than an argument.
        //
        // WHY 1 IS STRUCTURAL, not a fitted value. MAXPOLY caps how many EXTRA
        // differences a single bubble may absorb. At k=31 and human
        // heterozygosity ~1/1000, the chance a second het site falls inside the
        // same 31 bp window is small and a third is rare, so 1 already captures
        // essentially all genuine clustering; larger values mostly admit paths
        // through repeats. The sweep agrees rather than defines it -- recall
        // saturates by P=2 (0.808 at both 2 and 3) while false positives keep
        // rising (F1 0.861 / 0.850 / 0.850 with coherence, 0.828 / 0.794 /
        // 0.765 without). A genome with much higher heterozygosity would
        // justify a higher cap, which is the property that makes this a
        // reasoned bound rather than a benchmark constant.
        // DiscoSNP++'s own -P default is 3. Measured here, 1 is better on this
        // data both WITHOUT coherence (F1 0.828 / 0.794 / 0.765 for P=1/2/3)
        // and WITH it (0.861 / 0.850 / 0.850): recall saturates by P=2 while
        // false positives keep accumulating. Defaulting to 3 would have made a
        // bare run silently worse than every number reported for this channel.
```

What the code now says:

```
        // MAXPOLY = 1, MEASURED. SNV F1 on the tuning window (r2) and four held-out chr20
        // windows:
        //     window      MAXPOLY=1   MAXPOLY=2   MAXPOLY=3 (DiscoSNP++ -P default)
        //     r2 (tune)     0.886       0.881       0.881
        //     r3            0.892       0.886       0.879
        //     na            0.958       0.958       0.956
        //     r4            0.894       0.895       0.897
        //     r5            0.937       0.923       0.916
        //     held-out mean 0.9203      0.9155      0.9120
        // A higher cap adds recall (+6 TP on r2) but about twice as many false positives,
        // because a second polymorphism inside one walk is more often two error branches
        // than two real heterozygous sites. The structural reason: at k=31 and human
        // heterozygosity of about 1/1000, a second het site inside one 31 bp window is
        // rare, so larger values mostly admit paths through repeats. A genome with much
        // higher heterozygosity would justify a higher cap.
```

## Multi-polymorphism extension

Source: `include/caps_caller.h`, lines 2623-2631 before the documentation revision.

Original comment:

```
            // first. Requiring identical bases at every step kills the bubble
            // outright, and measurement says that is our single largest source
            // of missed variants: of 115 missed truth SNVs on HG002 r2, 35
            // (30%) die exactly here as "dead end".
            // DiscoSNP++ does not stop: expand_heart is called with
            // nb_polymorphism+1 over successor pairs carrying DIFFERENT bases
            // (Bubble.cpp:636), up to max_polymorphism. We benchmark them at
            // -P 3, so this is a capability we were being measured against
            // without having ported it.
```

What the code now says:

```
            // first. Requiring identical bases at every step kills the bubble
            // outright; on the HG002 tuning window 30% of missed truth SNVs ended here.
            // DiscoSNP++ continues: expand_heart is called with nb_polymorphism+1 over
            // successor pairs carrying DIFFERENT bases (Bubble.cpp:636), up to
            // max_polymorphism.
```

## Unitig extension: motivation

Source: `include/caps_caller.h`, lines 2679-2686 before the documentation revision.

Original comment:

```
        // MEASURED MOTIVATION, not a guess: without this a bubble sequence is
        // flank(31) + path, mean 63 bp on HG002 r2. bwa cannot place a 63-mer
        // uniquely in a 63 Mb chromosome, so found bubbles were being lost at
        // the LIFT rather than at the calling -- 397 bubbles produced only 243
        // lifted calls. DiscoSNP++ does not have this problem because it
        // extends each bubble into its surrounding unitig before output
        // (expand_one_simple_path / its -t traversal); this is a layer where we
        // were strictly worse than the tool we are competing with.
```

What the code now says:

```
        // Without extension a bubble sequence is its 31 bp flank plus the path (mean
        // 63 bp on the tuning window), which bwa cannot place uniquely in chr20, so
        // bubbles were lost at the LIFT rather than at calling. DiscoSNP++ likewise
        // extends each bubble into its surrounding unitig before output
        // (expand_one_simple_path, its -t traversal).
```

## Why indel recall is bounded

Source: `include/caps_caller.h`, lines 2671-2690 before the documentation revision.

Original comment:

```
        // Traced every missed indel on HG002 r2 back into kc. At the missed
        // sites the branching node has GOOD coverage (12-30) but exactly ONE
        // successor -- and that holds with NO coverage floor at all: raw
        // successor counts come back as [0,0,21,0]. The alternate haplotype's
        // k-mers are not below MINC, they are ABSENT.
        //
        // The reason: 66% of the missed indels are homopolymer-length changes
        // and most of the rest are tandem-repeat expansions. Inserting a base
        // into a run of the same base does not create a distinct k-mer path --
        // it makes the run longer. Probing the reads directly confirms it: for
        // a G->GT call inside a T-run, the REF and ALT probes are the SAME
        // STRING (TCTGGGTTTTTGTTTTTCGGGTTTTTTTTTTT), each found in 6 reads.
        //
        // There is no second path, so no bubble exists -- for us OR for any
        // de Bruijn caller, DiscoSNP++ included. This is a representational
        // limit of the graph, not a defect in the traversal, and it is why
        // nine successive attempts at filters, gates and orderings all failed
        // to move indel recall. Resolving these events needs length-aware
        // evidence (read pileup depth over the run), which is a different
        // mechanism from bubble finding.
```

What the code now says:

```
        // At the missed indel sites on the tuning window the branching node has good
        // coverage (12-30) but exactly ONE successor even with no coverage floor: the
        // alternate haplotype's k-mers are absent, not below MINC. Most of these events
        // are homopolymer-length changes and the rest mostly tandem-repeat expansions.
        // Inserting a base into a run of that base makes the run longer without creating
        // a distinct k-mer path, so no bubble exists for any de Bruijn caller, DiscoSNP++
        // included. This is a representational limit of the graph, not a traversal
        // defect, and filters, gates or orderings cannot recover these events. They need
        // length-aware evidence such as read depth over the run.
```

## Onodera superbubble detection

Source: `include/caps_caller.h`, lines 2761-2777 before the documentation revision.

Original comment:

```
        // 2013. The earlier version here walked each branch as a unitig and
        // took the earliest node all branches reached -- that is NOT this
        // algorithm. It stops at any junction, so it can only ever find
        // parallel simple paths and never explores a structure with internal
        // branching (which is why multi-allelic sites came out as 0).
        //
        // The load-bearing condition is "a node enters the frontier only once
        // ALL of its parents have been visited". That is what allows nodes
        // reachable by several internal paths to be handled correctly.
        //   S      = frontier: discovered AND all parents visited
        //   seen   = discovered but not yet ready
        //   sink   = |S| == 1 and |seen| == 0  (only the exit remains)
        // Aborts on a tip (no successors) or on returning to the entrance.
        // abort-reason counters: 2,628 branching nodes yield only 282 bubbles,
        // so ~89% fail somewhere. Which reason dominates decides where the
        // missing recall is (ours 0.468 vs DiscoSNP++ 0.763 from the same
        // k-mer set), rather than guessing.
```

What the code now says:

```
        // 2013. A unitig walk that takes the earliest node all branches reach can only
        // find parallel simple paths; this is the full algorithm.
        //
        // The load-bearing condition is "a node enters the frontier only once
        // ALL of its parents have been visited". That is what allows nodes
        // reachable by several internal paths to be handled correctly.
        //   S      = frontier: discovered AND all parents visited
        //   seen   = discovered but not yet ready
        //   sink   = |S| == 1 and |seen| == 0  (only the exit remains)
        // Aborts on a tip (no successors) or on returning to the entrance.
        // The ab_* abort-reason counters are diagnostics.
```

## Independent indel BFS (CAPS_DBG_IBFS)

Source: `include/caps_caller.h`, lines 2780-2783 before the documentation revision.

Original comment:

```
        // OPT-IN (CAPS_DBG_IBFS=1) and never on by default. The indel claim is
        // withdrawn -- see working note CLAIM2_FINAL.md -- so this path is kept for the
        // record and for future work, but it is unmeasured at scale and must
        // not be able to affect a shipped number. IBFS_MAXINDEL mirrors D=100.
```

What the code now says:

```
        // OPT-IN (CAPS_DBG_IBFS=1), off by default and not used by any reported result;
        // unmeasured at scale. IBFS_MAXINDEL mirrors DiscoSNP++'s D=100.
```

## Both orientations

Source: `include/caps_caller.h`, lines 2824-2826 before the documentation revision.

Original comment:

```
        // This code had no notion of strand at all, which is the single
        // largest reason graph-only recall (0.512) trailed DiscoSNP++'s
        // (0.763) on the very same k-mer set.
```

What the code now says:

```
        // Extending both sides is required for recall: without it graph-only recall was
        // 0.512 against DiscoSNP++'s 0.763 on the same k-mer set.
```

## Parallel superbubble traversal

Source: `include/caps_caller.h`, lines 2836-2854 before the documentation revision.

Original comment:

```
        // ── LAYER 15: PARALLELISE THE SUPERBUBBLE TRAVERSAL ─────────────
        // This region used to be `if(!WANT_SB)`, i.e. the ENTIRE traversal ran
        // single-threaded whenever superbubble mode was on -- which is every
        // Method B run. Traversal is the dominant cost (55 s of 98.5 s at full
        // chr20), so the caller was giving up ~11 of 12 cores for its longest
        // stage. The stated reason was that find_sb's abort counters are
        // shared; that is a two-line fix (atomics), not a reason to serialise.
        //
        // What is shared, and how each is handled:
        //   * dbg_bubbles  -- already per-thread (`local_bubbles`), merged once
        //   * sb_sites     -- now per-thread (`local_sb`), merged once
        //   * ab_* counters-- diagnostics only, made atomic
        // Nothing else in the loop writes shared state; kc, MINC and the walk
        // lambdas are read-only.
        //
        // Bubble and site ORDER changes with thread scheduling. That is already
        // documented as harmless -- every record is self-contained and the VCF
        // is sorted downstream -- and the gate for this change is that the
        // bubble COUNT and the resulting F1 are unchanged.
```

What the code now says:

```
        // ── PARALLEL SUPERBUBBLE TRAVERSAL ─────────────────────────────────
        // Traversal is the caller's longest stage, so it runs in parallel in superbubble
        // mode too. Shared state:
        //   * dbg_bubbles  -- per-thread (`local_bubbles`), merged once
        //   * sb_sites     -- per-thread (`local_sb`), merged once
        //   * ab_* counters-- diagnostics only, atomic
        // Nothing else in the loop writes shared state; kc, MINC and the walk lambdas are
        // read-only. Bubble and site ORDER depend on thread scheduling, which is harmless:
        // every record is self-contained and the VCF is sorted downstream. The bubble
        // COUNT and F1 do not change.
```

## Superbubble exit

Source: `include/caps_caller.h`, lines 2885-2891 before the documentation revision.

Original comment:

```
                // Locate the exit with the REAL superbubble algorithm
                // (find_sb, the Onodera/BubbleGun port), then read each
                // allele's string off the branch trace that reaches it.
                // The previous version picked "the earliest node every branch
                // reaches", which is an approximation: it ignores the
                // all-parents-visited condition and therefore never explores a
                // structure with internal branching.
```

What the code now says:

```
                // The exit is located with the superbubble algorithm (find_sb, the
                // Onodera/BubbleGun port), which respects the all-parents-visited condition,
                // and each allele's string is read off the branch trace that reaches it.
```

## Variant groups (SKA lo)

Source: `include/caps_caller.h`, lines 2893-2929 before the documentation revision.

Original comment:

```
                // MEASURED PROBLEM: of 755 branching nodes with 3-4 successors,
                // find_sb fails on 734 (97%) and none of the 21 successes has
                // more than 2 branches converging. Multi-allelic sites die
                // inside find_sb, whose exit test needs the frontier to drain
                // to exactly one node with nothing pending -- a condition that
                // 3-4 branches through repeat-rich sequence essentially never
                // satisfy. That is structural to the Onodera formulation.
                //
                // SKA lo solves the same problem by relaxing WHAT COUNTS AS A
                // SITE: "a variant group is identified when AT LEAST TWO paths
                // start from the same entry node and end at the same exit
                // node". Not all paths -- at least two. The traces already hold
                // every node each branch visits, so the group can be read
                // straight out of them: pick the earliest node reached by the
                // most branches (>=2). A 4-way node where 3 branches meet is
                // then a 3-allele site instead of a discarded one.
                //
                // Tried first and only when the strict superbubble finds
                // nothing, so simple two-path bubbles keep the tighter
                // definition and SNV calling is unaffected.
                //
                // OUTCOME, MEASURED: this DOES form variant groups (sites
                // 818 -> 827) but ZERO of them have more than two allele paths.
                // The group node is selected to maximise the number of
                // converging branches, so a 3-branch meeting point would win if
                // one existed. None does.
                //
                // The reason is biological, not algorithmic: at 30x DIPLOID
                // coverage there is no third haplotype for a third branch to
                // converge on. The 3-4 successor nodes are error or repeat
                // branches, which die out -- which is exactly what MINC and the
                // tip filter exist to make them do. T5.2 is therefore not
                // reachable on diploid human data by any traversal change; it
                // needs a sample that genuinely carries >2 haplotypes (the
                // tetraploid construction of T5.3, or a pooled/polyploid
                // sample). Kept because it is strictly more general than the
                // strict superbubble and costs nothing when it finds nothing.
```

What the code now says:

```
                // find_sb's exit test needs the frontier to drain to exactly one node with
                // nothing pending, which 3-4 branches through repeat-rich sequence almost never
                // satisfy, so multi-allelic sites fail there. SKA lo relaxes what counts as a
                // site: "a variant group is identified when AT LEAST TWO paths start from the
                // same entry node and end at the same exit node". The traces already hold every
                // node each branch visits, so the group is the earliest node reached by the most
                // branches (>= 2). This is tried only when the strict superbubble finds nothing,
                // so two-path bubbles and SNV calling are unaffected.
                //
                // On 30x diploid human data no group has more than two allele paths, because no
                // third haplotype exists to converge on; the extra successors are error or repeat
                // branches. Multi-allelic recovery therefore needs a sample with more than two
                // haplotypes (the tetraploid construction, or a pooled sample).
```

## Shortest closure ordering

Source: `include/caps_caller.h`, lines 2945-2974 before the documentation revision.

Original comment:

```
                    // This previously maximised BRANCH COUNT first and used
                    // depth only as a tiebreak, so a distant node reached by 3
                    // branches beat a near one reached by 2. DiscoSNP++ orders
                    // it the other way round and never proposes the distant
                    // exit at all: their indel BFS stops at the smallest
                    // closing event (Bubble.cpp:238-243, "if (insert_size ==
                    // found_del_size-1) break" and "if (insert_size >
                    // found_del_size) continue"), which is why long spurious
                    // indels are absent from their output rather than filtered
                    // out of it.
                    //
                    // MEASURED AFTERWARDS: this ordering change is INERT in
                    // practice. Variant groups only fire when ns > 2, which is
                    // rare, so nearly every indel comes from find_sb on a
                    // 2-successor node -- and find_sb returns the FIRST node
                    // satisfying its drain condition, which is already the
                    // earliest in traversal order. Closure distances measured
                    // 709 of 827 in the 32-39 bp bin, i.e. about one k-mer past
                    // the divergence, which is the correct shape for real
                    // bubbles rather than evidence of distant closure. The
                    // long-indel false positives therefore do NOT come from
                    // choosing a far exit, and this ordering is kept only
                    // because it matches the reference implementation.
                    //
                    // Adopting their ORDERING (a generative bound) rather than
                    // adding another post-hoc length filter -- four of those
                    // were tried and every one cost more true positives than
                    // false ones, twice catastrophically (F1 0.537 -> 0.035).
                    // Preferring the nearest reconvergence keeps short events
                    // and simply does not generate the long ones.
```

What the code now says:

```
                    // Among candidate exits the nearest reconvergence is preferred, as in DiscoSNP++,
                    // whose indel BFS stops at the smallest closing event (Bubble.cpp:238-243), so
                    // long spurious events are not generated rather than filtered afterwards. In
                    // practice the ordering rarely matters: groups fire only when ns > 2, and
                    // find_sb already returns the earliest exit.
```

## Partial convergence

Source: `include/caps_caller.h`, lines 2998-3024 before the documentation revision.

Original comment:

```
                    // PARTIAL CONVERGENCE (kept, but it is NOT what blocks T5.2).
                    // MEASURED: of 755 branching nodes with 3-4 successors,
                    // find_sb FAILS on 734 (97%), and of the 21 that succeed
                    // none has more than 2 branches converging. Multi-allelic
                    // sites are lost inside find_sb itself -- its exit
                    // condition needs the frontier to drain to exactly one node
                    // with nothing pending, which 3-4 branches through
                    // repeat-rich sequence essentially never satisfy. That is a
                    // structural property of the Onodera formulation, not a
                    // threshold, so T5.2 is not reachable by relaxing this
                    // check. Kept anyway because it is strictly more
                    // information and cannot affect SNVs.
                    //
                    // This used to require EVERY branch to reach the exit
                    // (`if (!found) ok = false`) and discarded the whole site
                    // otherwise. Measured consequence: 755 branching nodes have
                    // 3-4 successors, yet not one produced a 3+-allele site --
                    // every multi-allelic candidate was thrown away because one
                    // of its branches died in a tip or a repeat.
                    //
                    // A 4-successor node where 3 branches converge and 1 dies
                    // is a 3-ALLELE SITE, not a failure. Keeping the branches
                    // that DO converge is strictly more information than
                    // discarding all of them, and it cannot affect SNV calling
                    // because SNVs come from the pairwise lockstep walk, which
                    // is untouched by this. At least 2 convergent branches are
                    // required -- one path is not a bubble.
```

What the code now says:

```
                    // PARTIAL CONVERGENCE. A 4-successor node where 3 branches converge and 1 dies
                    // in a tip or repeat is a 3-allele site, so the branches that do converge are
                    // kept instead of discarding the site. At least 2 convergent branches are
                    // required (one path is not a bubble). SNV calling is unaffected, because SNVs
                    // come from the pairwise lockstep walk. Multi-allelic loss on diploid data is
                    // inside find_sb (see above), not here.
```

## Pairwise walk is not skipped

Source: `include/caps_caller.h`, lines 3059-3076 before the documentation revision.

Original comment:

```
                // This used to `continue`, so enabling the superbubble path
                // REPLACED the pairwise walk. Measured cost of that
                // substitution on HG002 r2: SNV F1 0.881 -> 0.791 (-29 TP,
                // +42 FP) -- the superbubble finds indels the lockstep walk
                // structurally cannot, but it is WORSE at plain SNVs, because
                // its per-branch unitig traces are a coarser view of a simple
                // two-path bubble than walking the two paths in lockstep.
                //
                // The two are complementary, not alternatives:
                //   pairwise walk -> SNVs   (lockstep, same base both paths)
                //   superbubble   -> INDELS (paths of different length)
                // so both run at every branching node and each contributes the
                // variant class it is actually good at. Superbubble records are
                // tagged so the emitter can keep only its different-length
                // ones; see SB_ONLY_INDEL below.
                //
                // (`CAPS_DBG_SB_REPLACE=1` restores the old substituting
                // behaviour for re-measurement.)
```

What the code now says:

```
                //
                // The pairwise walk and the superbubble path are complementary:
                //   pairwise walk -> SNVs   (lockstep, same base both paths)
                //   superbubble   -> INDELS (paths of different length)
                // Both run at every branching node. Replacing the pairwise walk with the
                // superbubble path lowered SNV F1 (0.881 -> 0.791 on the tuning window), because
                // per-branch unitig traces are a coarser view of a two-path bubble. Superbubble
                // records are tagged so the emitter keeps only the different-length ones (see
                // SB_ONLY_INDEL). CAPS_DBG_SB_REPLACE=1 reproduces the substituting variant for
                // re-measurement only.
```

## Pre-walk coverage reject (rejected)

Source: `include/caps_caller.h`, lines 3082-3091 before the documentation revision.

Original comment:

```
                    // REFUTED: a cheap coverage-ratio reject before the walk.
                    // Rationale was sound -- 1,226,535 branching nodes yield only
                    // 111,499 bubbles, so 91% of walks are wasted and traversal is
                    // 55.5 s of 98.5 s. But measured at full chr20 it changed
                    // nothing (98.7 s vs 98.5 s, traversal 55.25 s vs 55.48 s) and
                    // cost 300 bubbles. Reason: the test was gated on ns > 2, and
                    // almost every branching node has exactly 2 successors, so it
                    // essentially never fired. The waste is not in walks that a
                    // coverage test can predict -- it is in walks that fail to
                    // CLOSE, which is only knowable by walking.
```

What the code now says:

```
                    // No cheap pre-walk reject is applied: most failed walks fail to CLOSE, which
                    // only walking reveals, and a coverage-ratio pre-test changed neither run time
                    // nor results at full chr20.
```

## Indel discovery (start_indel_prediction)

Source: `include/caps_caller.h`, lines 3111-3146 before the documentation revision.

Original comment:

```
            // Ported from kissnp2 Bubble.cpp:200-285, reimplemented from the
            // published algorithm. THE recall gap, found by tracing their
            // control flow rather than their filters.
            //
            // THEIR STRUCTURE (Bubble.cpp:310-346): at every branching node,
            // for every successor PAIR, they call start_snp_prediction() and
            // then start_indel_prediction() -- INDEPENDENTLY. The indel search
            // does not require the SNP search to succeed, and it does not
            // require any superbubble to exist.
            //
            // OURS, BEFORE THIS: SNVs came from the pairwise lockstep walk
            // (runs at every branching node) but indels came ONLY from
            // find_sb, the Onodera superbubble detector. Measured at full
            // chr20, find_sb SUCCEEDS ON 15.9% OF NODES:
            //
            //     ok=195,133  exhausted=960,541  budget=51,942  tip=17,409
            //
            // so at 1,030,061 branching nodes we never attempted an indel at
            // all. That is the same asymmetry as the results: SNV recall 0.827
            // from the path that runs everywhere, indel recall 0.414 from the
            // path that runs on a sixth of the graph. It is not the
            // homopolymer limit -- that explains the hardest subset, this
            // explains the population.
            //
            // THE ALGORITHM. Extend ONE path while the other waits, until the
            // extended path's last base equals the base that would let the two
            // reconverge, then close in LOCKSTEP with the existing walk() --
            // which is exactly their expand(). Note that with their default
            // `authorised_branching=0` the breadth-first queue never holds more
            // than one live node (they abandon the search at the first
            // junction), so this is a unitig extension, not a search.
            //
            // The shape test is applied inline here for the same reason it is
            // applied to superbubble records: lockstep closure guarantees the
            // tail agrees, but the inserted block must still be ONE contiguous
            // event.
```

What the code now says:

```
            // Reimplemented from the published DiscoSNP++ algorithm (kissnp2
            // Bubble.cpp:200-285). At every branching node and every successor PAIR,
            // DiscoSNP++ runs its SNP search and its indel search INDEPENDENTLY
            // (Bubble.cpp:310-346): the indel search needs neither the SNP search nor a
            // superbubble to succeed. Indels taken only from find_sb would be attempted at
            // the small share of nodes where find_sb succeeds (15.9% at full chr20).
            //
            // THE ALGORITHM. Extend ONE path while the other waits, until the extended path's
            // last base equals the base that would let the two reconverge, then close in
            // LOCKSTEP with walk(), which is DiscoSNP++'s expand(). With their default
            // `authorised_branching=0` the queue never holds more than one live node (the
            // search stops at the first junction), so this is a unitig extension, not a search.
            //
            // The shape test is applied inline, as for superbubble records: lockstep closure
            // guarantees the tail agrees, but the inserted block must be ONE contiguous event.
```

## Read coherence

Source: `include/caps_caller.h`, lines 3082-3100 before the documentation revision.

Original comment:

```
        // THE FP MECHANISM, named in this repo's own docs: a CHIMERIC JUNCTION.
        // Every k-mer along a bubble path can be present in kc while the FULL
        // path exists in no single read -- each k-mer contributed by a
        // different read, joined at a shared anchor. k-mer counts are blind to
        // this by construction, which is why our precision (0.865) trails
        // DiscoSNP++'s (0.951) even though our recall is now higher.
        //
        // eBWT2SNP's precision guarantee is exactly this test: every emitted
        // fragment of length 2k+1 must be an actual SUBSTRING of at least C
        // real reads. DiscoSNP++ gets the equivalent from kissreads2, a whole
        // separate tool that maps every read back onto every bubble.
        //
        // WE DO NOT NEED THAT TOOL. The reads are already in memory (`seqs`),
        // so one pass with a small probe index gives a stronger check than
        // theirs at a fraction of the cost: index each bubble by the 31-mer
        // SPANNING its variant, sweep the reads once, and verify the full
        // fragment really occurs. This is a case where having the reads
        // resident -- a consequence of being a compressor -- buys something a
        // standalone caller has to pay a separate pass for.
```

What the code now says:

```
        // The false-positive mechanism this targets is a CHIMERIC JUNCTION: every k-mer
        // along a bubble path can be present in kc while the FULL path occurs in no single
        // read (each k-mer contributed by a different read, joined at a shared anchor).
        // k-mer counts cannot see this. The test, as in eBWT2SNP, is that the fragment
        // must be a substring of real reads. DiscoSNP++ obtains the equivalent from
        // kissreads2, a separate pass mapping every read onto every bubble. Here the
        // reads are already in memory (`seqs`), so bubbles are indexed by the 31-mer
        // spanning the variant, the reads are swept once, and the full fragment is
        // verified.
```

## COHC derived from coverage

Source: `include/caps_caller.h`, lines 3105-3123 before the documentation revision.

Original comment:

```
            // This was a constant 2, swept on HG002 r2 at 30x. That does not
            // generalise: "how many reads must contain this fragment" scales
            // with depth. At 10x, 2 reads is a large fraction of a haplotype's
            // coverage and the filter becomes punitive; at 100x, 2 reads is
            // noise and the filter stops filtering.
            //
            // H is the measured haploid k-mer depth (the single-copy histogram
            // peak, computed above). A heterozygous allele is covered by ~H
            // reads, so requiring a small FRACTION of H tracks the data.
            // H/10 with a floor of 2 reproduces the swept value exactly at
            // this benchmark's depth (H=20 -> 2) while adapting elsewhere,
            // which is what makes it a formula rather than a fitted constant
            // that happens to be right once.
            //
            // NOTE the contrast with MINC, where a per-dataset derivation
            // (the k-mer histogram valley) generalised WORSE than a structural
            // constant: it drifted upward at full scale and cost 0.057 recall.
            // Derivation is right when the quantity genuinely scales with the
            // data -- read support does, the singleton-error floor does not.
```

What the code now says:

```
            // The number of reads that must contain a fragment scales with depth, so COHC is a
            // fraction of H, the measured haploid k-mer depth (a heterozygous allele is
            // covered by about H reads): H/10 with a floor of 2. At the benchmark depth (H=20)
            // this equals 2, the value swept on the tuning window. Unlike the error floor
            // MINC, which does not scale with depth, read support does.
```

## Fragment half-width (HALFW): derivation

Source: `include/caps_caller.h`, lines 3131-3161 before the documentation revision.

Original comment:

```
            // Fragment half-width. The containment test requires a read to
            // contain the WHOLE fragment (2*HALFW+1 = 63 bp at HALFW=31), which
            // for indels is the binding constraint -- a 148 bp read must cover
            // 63 bases centred on the event, so genuine indels near read ends
            // are lost. SNVs tolerate it because they are far more numerous.
            // ── FRAGMENT HALF-WIDTH, DERIVED FROM MEASURED READ LENGTH ──
            // The containment test asks whether a READ contains the whole
            // 2*HALFW+1 fragment, so HALFW controls how much unique context
            // each allele must be anchored in. That is the read-length signal
            // a linked/read-coloured graph exploits: a longer anchor spans a
            // repeat that a 63 bp fragment sits inside, so repeat/paralog
            // collapses -- our dominant false-positive mechanism, carrying 12x
            // the allele depth of true sites -- stop finding support.
            //
            // It was a constant 31 (63 bp) against ~148 bp reads, which wastes
            // more than half the available context. Measured on HG002 r2:
            //     HALFW=31   63 bp   TP=336  FP=20  P=0.944  F1=0.886
            //     HALFW=45   91 bp   TP=336  FP=16  P=0.955  F1=0.891
            //     HALFW=60  121 bp   TP=327  FP=14  P=0.959  F1=0.880
            // Precision rises monotonically; recall is untouched until the
            // fragment stops fitting comfortably in a read, at which point
            // genuine variants near read ends are lost (9 TPs at HALFW=60).
            //
            // So the ceiling is a property of READ LENGTH, not of this dataset,
            // and the parameter is a fraction of it rather than a fitted
            // number: 0.30 x median read length puts the fragment at ~60% of a
            // read, keeping ~40% as slack for the variant to sit anywhere in
            // it. On 148 bp reads that gives 44; on 100 bp reads 31 (the old
            // default, which was right for the read length it was tuned on);
            // on 250 bp reads 75. Floor of 31 because the probe must straddle
            // a 31-mer junction.
```

What the code now says:

```
            // FRAGMENT HALF-WIDTH. A supporting read must contain the WHOLE fragment
            // (2*HALFW+1 bases), so HALFW sets how much unique context anchors each allele.
            // More context raises precision (repeat and paralog collapses stop finding
            // support) but loses true variants near read ends, so the useful ceiling depends
            // on read length. halfw_auto below is 0.30 x the median read length (floor 31,
            // because the probe must straddle a 31-mer junction): 44 on 148 bp reads.
```

## Fragment half-width (HALFW): held-out validation

Source: `include/caps_caller.h`, lines 3176-3194 before the documentation revision.

Original comment:

```
            // HELD-OUT VALIDATION REJECTED IT. On the tuning window the
            // derived value looked like a free win -- 4 false positives removed
            // and not one true positive lost. It did not transfer:
            //     window   HALFW=31 -> derived 44
            //     r2 (tune)  0.886  ->  0.891   +0.005
            //     r3         0.892  ->  0.892    0.000
            //     na         0.958  ->  0.954   -0.004
            //     r4         0.894  ->  0.891   -0.003
            //     r5         0.937  ->  0.935   -0.002
            // Mean over the four HELD-OUT windows: -0.002. On r2 the longer
            // anchor removed false positives for free; elsewhere it removes
            // almost none and costs true positives near read ends. The
            // mechanism is real -- precision does rise monotonically with
            // context -- but the operating point that pays is a property of
            // that window, not of read length, so the formula does not rescue
            // it. Default stays 31, the value every validated measurement used.
            // The derivation is kept because it is the right SHAPE if a future
            // dataset has materially different read lengths; opt in with
            // CAPS_DBG_HALFW.
```

What the code now says:

```
            // The derived value is NOT the default: it raised F1 on the tuning window but not
            // on the four held-out chr20 windows (mean change -0.002), so HALFW stays 31
            // (63 bp), the value of every reported measurement. CAPS_DBG_HALFW opts in.
            // Measurements: docs/extras/ENGINEERING_NOTES_CALLER.md.
```

## Probe straddles the junction

Source: `include/caps_caller.h`, lines 3220-3241 before the documentation revision.

Original comment:

```
                // koff was vpos (=31), i.e. exactly the start of the path, so
                // the probe k-mer covered path+rext and NEVER included the
                // flank. For an SNV that is harmless -- the variant base is at
                // position 0 of the k-mer, so a read carrying the allele still
                // matches. For an INDEL it removes the whole discriminating
                // signal: what separates a real indel from a repeat's distant
                // closure is whether reads actually cross the FLANK->PATH
                // junction, and that junction was never in the probe.
                //
                // Centring the k-mer on the junction (half in the flank, half
                // in the path) means a read only matches if it genuinely spans
                // the event. This is the same principle Lancet and Scalpel
                // apply -- "evaluate the support of each alternate haplotype
                // with respect to the RAW READ DATA" -- applied at the place
                // that actually distinguishes indels.
                // TWO PROBES: one at the path start (as before) and one
                // centred on the junction. A read matches if EITHER is present,
                // so the junction evidence is available without demanding that
                // every supporting read span 15 bases of flank as well as the
                // event -- which is what a single centred probe required, and
                // it cost more recall than the precision it bought
                // (P 0.733 -> 0.791 but F1 0.611 -> 0.535, measured).
```

What the code now says:

```
                //
                // A probe that starts at the path never includes the flank. For an SNV that is
                // harmless (the variant base is position 0 of the k-mer), but for an INDEL the
                // discriminating signal is whether reads cross the FLANK->PATH junction. As in
                // Lancet and Scalpel, support is evaluated against the raw reads at the place that
                // distinguishes the event.
                // TWO PROBES: one at the path start and one centred on the junction. A read
                // matches if EITHER is present. A single centred probe alone required every
                // supporting read to span 15 flank bases as well as the event and cost more
                // recall than the precision it bought (measured F1 0.611 -> 0.535).
```

## Placement dispersion: not a filter

Source: `include/caps_caller.h`, lines 3273-3276 before the documentation revision.

Original comment:

```
            // NOT FILTERED YET. Filters justified by mechanism alone have
            // repeatedly measured negative in this caller, so this
            // pass only COUNTS, and the separation is checked against vcfeval
            // before any bubble is rejected.
```

What the code now says:

```
            // This pass only COUNTS; it rejects nothing unless PLACEMODE (a diagnostic) is
            // set.
```

## Exclusive read assignment

Source: `include/caps_caller.h`, lines 3310-3318 before the documentation revision.

Original comment:

```
            // For a LENGTH-CHANGING bubble it is not harmless. The two
            // fragments share a flank and a right context and differ by the
            // indel, so in repetitive or homopolymeric sequence one read can
            // satisfy both `find`s -- and then a site with no allelic evidence
            // at all still clears COHC on both paths. That is a plausible
            // source of the precision asymmetry we actually observe
            // (SNV 0.950 vs INDEL 0.658, against DiscoSNP++'s indel 0.936):
            // kissreads2 assigns each read to its BEST path rather than to
            // every path it contains.
```

What the code now says:

```
            // For a LENGTH-CHANGING bubble it is not harmless. The two
            // fragments share a flank and a right context and differ by the
            // indel, so in repetitive or homopolymeric sequence one read can
            // satisfy both `find`s, and a site with no allelic evidence then
            // clears COHC on both paths. kissreads2 assigns each read to its BEST
            // path instead of every path it contains. This is a candidate cause of the
            // lower indel than SNV precision.
```

## Quality-aware support

Source: `include/caps_caller.h`, lines 3353-3365 before the documentation revision.

Original comment:

```
                                // kissreads2 scores each bubble path by read
                                // coverage AND average phred; we were counting
                                // containment only, and ignoring `quals`
                                // entirely even though it is passed in. A
                                // sequencing error is a low-quality base by
                                // definition, so a read whose variant base is
                                // low quality is not evidence for that allele.
                                // This is the discriminator that separates a
                                // real low-coverage het allele (few reads, high
                                // quality) from an error (few reads, low
                                // quality) -- something a coverage floor cannot
                                // do at all, which is why precision has to live
                                // here rather than upstream.
```

What the code now says:

```
                                // kissreads2 scores each bubble path by read
                                // coverage AND average phred. A sequencing
                                // error is a low-quality base by definition,
                                // so a read whose variant base is low quality
                                // is not evidence for that allele. This
                                // separates a real low-coverage het allele
                                // (few reads, high quality) from an error (few
                                // reads, low quality), which a coverage floor
                                // cannot do.
```

## Per-position coverage test (rejected)

Source: `include/caps_caller.h`, lines 3425-3437 before the documentation revision.

Original comment:

```
            // REFUTED: a per-position coverage test, ported from kissreads2
            // (fragment.cpp:86-91, "for(i..stop) if(local_coverage[i]<min)
            // return false"). It measured EXACTLY no change at full chr20
            // (kept=111,499 killed=44,248, identical), and the reason is
            // structural, not a bug: our counter only increments after
            // r.find(f) succeeds -- i.e. after a read already contains the
            // WHOLE fragment -- so every position necessarily has the same
            // count as the containment tally and the test cannot ever fire.
            // kissreads2 differs in that reads OVERLAPPING a fragment
            // partially still contribute coverage to the positions they span,
            // which is what exposes a thin interior position. Reproducing that
            // needs partial-overlap alignment, not a stricter test over the
            // counts we already have.
```

What the code now says:

```
            // No per-position coverage test (kissreads2 fragment.cpp:86-91) is applied: the
            // counter here increments only after a read contains the WHOLE fragment, so every
            // position has the containment count and such a test cannot fire. kissreads2
            // also counts partially overlapping reads; reproducing that needs partial-overlap
            // alignment.
```

## Bubble chaining

Source: `include/caps_caller.h`, lines 3487-3507 before the documentation revision.

Original comment:

```
        // WHY, and this is a measured motivation rather than a design
        // preference. Emitting each bubble as its own short sequence produced
        // a mean 63 bp dcontig; extending each one into its unitig only
        // reached 238 bp and moved F1 not at all (0.914 -> 0.913). The reason
        // is in the assembly literature rather than in our code: unitigs
        // "break at heterozygotes", so in a heterozygous genome the very next
        // junction that stops the extension IS the next variant. Unitigs are
        // the wrong extension unit for het data.
        //
        // A BUBBLE CHAIN is the right one: BubbleGun defines it as a linear
        // stretch of bubbles where the sink of one is the source of the next.
        // Extending THROUGH consecutive bubbles instead of stopping at them
        // gives a genuinely long, uniquely-mappable sequence carrying several
        // variants -- which fixes the mapping problem and the closure
        // granularity in one change. It is also what makes phasing possible at
        // all (PHASM builds haplotypes exactly this way), since consecutive
        // bubbles spanned by the same reads are by definition linked.
        //
        // Implementation: index sites by entry node; a site whose right unitig
        // extension lands on another site's entry node is its successor. Walk
        // each maximal chain once, marking members consumed.
```

What the code now says:

```
        // Unitigs "break at heterozygotes": in a heterozygous genome the junction that
        // stops a unitig extension is the next variant, so extending each bubble into its
        // unitig barely lengthens it (mean 63 -> 238 bp on the tuning window, F1
        // unchanged). A BUBBLE CHAIN (BubbleGun: a linear stretch of bubbles where the
        // sink of one is the source of the next) extends THROUGH consecutive bubbles and
        // gives a long, uniquely mappable sequence carrying several variants. It is also
        // the basis of phasing (PHASM), since consecutive bubbles spanned by the same
        // reads are linked.
        //
        // Implementation: index sites by entry node; a site whose right unitig extension
        // lands on another site's entry node is its successor. Walk each maximal chain
        // once, marking members consumed.
```

## checkRepeatSize provenance

Source: `include/caps_caller.h`, lines 3532-3535 before the documentation revision.

Original comment:

```
                    // The structural gate their closure has and ours did not.
                    // Ported from kissnp2 Bubble.cpp:952 (`checkRepeatSize`),
                    // AGPL — reimplemented from the published algorithm, not
                    // copied.
```

What the code now says:

```
                    // Reimplemented from the published DiscoSNP++ algorithm (kissnp2
                    // Bubble.cpp:952, `checkRepeatSize`, AGPL); no code is copied.
```

## checkRepeatSize: why a graph test

Source: `include/caps_caller.h`, lines 3551-3560 before the documentation revision.

Original comment:

```
                    // THIS IS WHY THE EARLIER FILTERS FAILED. They asked
                    // whether the REFERENCE context is homopolymeric, and it
                    // does not separate -- measured over all 5,001 scored
                    // indels, true positives are MORE homopolymeric than false
                    // ones (frac run>=4: 0.531 vs 0.427). This asks a different
                    // question, about the GRAPH: not "is this sequence
                    // repetitive" but "is this indel's position determined".
                    // A real indel in a homopolymer still has a determined
                    // position when the flanks pin it down; a spurious one does
                    // not.
```

What the code now says:

```
                    // This is a question about the GRAPH ("is this indel's position determined"),
                    // not about the sequence ("is the reference context repetitive"). The latter does
                    // not separate: over all 5,001 scored indels, true positives are MORE
                    // homopolymeric than false ones (fraction with a run >= 4: 0.531 against 0.427).
                    // A real indel in a homopolymer still has a determined position when the flanks
                    // pin it down; a spurious one does not.
```

## Clean-indel shape

Source: `include/caps_caller.h`, lines 3570-3599 before the documentation revision.

Original comment:

```
                    // The real structural difference between their indel
                    // closure and ours, found by putting the two side by side.
                    //
                    // THEIRS: `start_indel_prediction` extends ONE path until
                    // it can close, then `expand` walks BOTH paths in LOCKSTEP
                    // and closes only when `nextNode1 == nextNode2`. So their
                    // two alleles are identical everywhere except one
                    // contiguous inserted block -- the geometry is enforced by
                    // the traversal and a non-indel shape can never be emitted.
                    //
                    // OURS: two unitig traces that happen to reach a common
                    // node via find_sb. Nothing requires them to agree
                    // anywhere in between, so paths differing at SEVERAL
                    // positions are emitted as one multi-base REF/ALT -- which
                    // vcfeval scores as an INDEL. This file already recorded
                    // the symptom ("two same-length paths can still differ at
                    // several positions ... those are emitted as a multi-base
                    // REF/ALT") without connecting it to the cause.
                    //
                    // The shape test is exact and parameter-free: s_long is
                    // s_short with ONE contiguous block inserted iff the common
                    // prefix and common suffix together already cover the whole
                    // of s_short.
                    //
                    //     LCP + LCS >= |s_short|
                    //
                    // Anything failing it differs by more than a single indel
                    // and is a traversal artefact, not a variant. This asks
                    // nothing about the sequence being repetitive -- the
                    // question that was already measured not to separate.
```

What the code now says:

```
                    // DiscoSNP++ extends ONE path until it can close and then walks BOTH paths in
                    // LOCKSTEP, closing only when the next nodes are equal, so its two alleles differ
                    // in exactly one contiguous block by construction. Two unitig traces that reach a
                    // common node via find_sb carry no such guarantee: paths differing at SEVERAL
                    // positions would be emitted as one multi-base REF/ALT, which vcfeval scores as an
                    // INDEL. The shape test is exact and parameter-free: s_long is s_short with ONE
                    // contiguous block inserted iff the common prefix and common suffix together
                    // cover the whole of s_short:
                    //
                    //     LCP + LCS >= |s_short|
                    //
                    // Anything failing it differs by more than a single indel and is a traversal
                    // artefact, not a variant.
```

## Shortest-closure rule

Source: `include/caps_caller.h`, lines 3650-3671 before the documentation revision.

Original comment:

```
                    // Their indel search is a BFS that STOPS at the smallest
                    // closing indel (Bubble.cpp:238-243):
                    //     if (insert_size == found_del_size-1) { clear(); break; }
                    //     if (insert_size >  found_del_size)   continue;
                    // so long indels are never GENERATED. Ours takes whatever
                    // exit the superbubble traversal reached, which through a
                    // repeat is often a distant one -- precisely the 21-32 bp
                    // false-positive cluster (448 of them at full chr20, with
                    // essentially zero true positives in that range).
                    //
                    // The structural difference matters: filtering long events
                    // after the fact removes them but also removes the real
                    // ones caught in the same net, which is why a blanket
                    // length rejection was already measured NEGATIVE on the
                    // other channel (precision 0.738 -> 0.813 but recall
                    // 0.562 -> 0.407, F1 0.637 -> 0.542). Preferring the
                    // SHORTEST closure keeps short events and simply never
                    // proposes the long ones.
                    //
                    // Applied here as the closest available analogue: when a
                    // site yields several indel-shaped allele pairs, keep the
                    // one with the smallest length difference.
```

What the code now says:

```
                    // DiscoSNP++'s indel search is a BFS that STOPS at the smallest closing indel
                    // (Bubble.cpp:238-243):
                    //     if (insert_size == found_del_size-1) { clear(); break; }
                    //     if (insert_size >  found_del_size)   continue;
                    // so long indels are never generated. A superbubble exit reached through a repeat
                    // is often a distant one, which produces the 21-32 bp false-positive mode. A
                    // blanket length rejection after the fact also removes real events (measured on
                    // the other channel: F1 0.637 -> 0.542), whereas preferring the shortest closure
                    // does not propose the long ones. The closest analogue here: when a site yields
                    // several indel-shaped allele pairs, keep the one with the smallest length
                    // difference.
```

## Homopolymer/STR artefact filter

Source: `include/caps_caller.h`, lines 3672-3700 before the documentation revision.

Original comment:

```
                    // ── HOMOPOLYMER / STR ARTEFACT FILTER ────────────────
                        // After the size bound the remaining indel false
                        // positives are dominated by 1 bp events, and they skew
                        // 2:1 toward DELETIONS (184 vs 91) while true 1 bp
                        // indels are balanced (986 del / 942 ins). Their
                        // context is homopolymeric -- CC, GG, TT are the most
                        // common. That is the classic homopolymer-length
                        // artefact, the known hard case for every
                        // reference-free caller.
                        //
                        // is_str_event() already encodes exactly this test and
                        // is used by the older bubble channel; it was simply
                        // never wired into the superbubble path. It asks
                        // whether the inserted/deleted sequence extends a run
                        // or tandem unit already present in the flank -- a
                        // measured property of the sequence, not a rule about
                        // this dataset.
                        // REFUTED AND DISABLED BY DEFAULT (opt in with
                        // CAPS_DBG_STR=1). It was added on a window-scale
                        // observation -- 1 bp false positives skewing 2:1
                        // toward deletions in homopolymeric context -- and it
                        // measured WORSE: it drops 9 true indels to remove 7
                        // false ones (F1 0.601 -> 0.535).
                        //
                        // The older channel had already measured exactly this
                        // and recorded it three hundred lines away in the same
                        // file: "Homopolymer context, by contrast, does NOT
                        // separate: 37.9% vs 34.4%". The filter was therefore
                        // a refuted idea re-derived from a smaller sample.
```

What the code now says:

```
                    // ── HOMOPOLYMER / STR ARTEFACT FILTER: OPT-IN (CAPS_DBG_STR=1) ──
                    // Rejects an indel whose inserted or deleted sequence extends a run or tandem
                    // unit already present in the flank (is_str_event). Off by default: it removed
                    // more true indels than false ones (F1 0.601 -> 0.535), consistent with the
                    // earlier finding that homopolymer context does not separate true from false
                    // indels here.
```

## Size bound: window-scale note

Source: `include/caps_caller.h`, lines 3712-3725 before the documentation revision.

Original comment:

```
                    // EARLIER NOTE, kept because it explains why this looked
                    // useless before: bounding indel length measured as pure
                    // loss at WINDOW scale (F1 0.582/0.601/0.611 at 8/12/20 bp)
                    // because that window contained only 16 indel FPs and none
                    // of them were long. Full scale has 448 long FPs. A filter
                    // must be evaluated where the failure mode actually
                    // occurs.
                    // The scored false positives look long (19-38 bp) but the
                    // RAW emitted indels max out at 20 bp -- the long shapes
                    // are the lift's left-alignment expanding short calls.
                    // Measured: capping at 8 / 12 / 20 bp left FP flat at 16
                    // and only removed true positives
                    // (F1 0.582 / 0.601 / 0.611). The 16 false positives are
                    // SHORT indels, so length carries no signal here.
```

What the code now says:

```
                    // At window scale a length cap only removed true positives (the window's few
                    // indel false positives were short); the long false-positive mode appears only at
                    // full chr20, so the cap is evaluated there.
```

## Graph-channel indel emission

Source: `include/caps_caller.h`, lines 3726-3746 before the documentation revision.

Original comment:

```
                    // INDEL EMISSION IS OFF BY DEFAULT -- MEASURED, NOT ASSUMED.
                    // The superbubble genuinely FINDS real indels: held-out
                    // indel TP rose on every window (r3 33->41, na 35->37,
                    // r4 18->25, r5 30->34). But it adds far more false
                    // positives than it recovers, and held-out indel F1 falls
                    // on 3 of 4 windows:
                    //     r3 0.518->0.471   na 0.654->0.514
                    //     r4 0.522->0.575   r5 0.659->0.586   mean -0.051
                    // An earlier tuning-window reading suggested +0.012, but
                    // that depended on a fitted MINC=3 and did not survive
                    // held-out testing -- two independent reasons to distrust
                    // it, both measured. The machinery stays (it is what lifts
                    // SNV precision 0.925 -> 0.957); only the indel OUTPUT is
                    // gated, behind CAPS_DBG_INDEL=1 for further work.
                    // NOTE the test that matters is the shape of the EMITTED
                    // RECORD, not the path lengths. Two same-length paths can
                    // still differ at several positions (e.g. ACGT vs TGCA);
                    // those are emitted as a multi-base REF/ALT, which vcfeval
                    // scores as an INDEL. Gating on s0.size()!=s1.size() alone
                    // let 23 such records through per window and was why indel
                    // F1 fell on all five windows even with "indels disabled".
```

What the code now says:

```
                    // GRAPH-CHANNEL INDEL EMISSION IS OFF BY DEFAULT (CAPS_DBG_INDEL=1 enables it).
                    // The superbubble path finds real indels (held-out indel TP rose on every window)
                    // but adds more false positives than it recovers: held-out indel F1 fell on 3 of
                    // 4 windows (mean -0.051). The machinery stays because it raises SNV precision
                    // (0.925 -> 0.957); only its indel OUTPUT is gated. The reported indel results
                    // come from the full caller (indel_pass), selected in `capsule_decode call` by
                    // CAPS_CALL_INDELS=1.
                    // The gate tests the shape of the EMITTED RECORD, not the path lengths: two
                    // same-length paths can differ at several positions (e.g. ACGT vs TGCA) and would
                    // be emitted as a multi-base REF/ALT, which vcfeval scores as an INDEL.
```

## One record per differing position

Source: `include/caps_caller.h`, lines 3921-3926 before the documentation revision.

Original comment:

```
                // Previously any bubble with more than one difference was
                // emitted as one multi-base REF/ALT and then discarded by the
                // clean-SNV gate, which is why raising -P changed nothing at
                // all: the walk found those bubbles and the emitter threw them
                // away. DiscoSNP++ writes one VCF record per SNP inside a
                // bubble; this does the same.
```

What the code now says:

```
                // DiscoSNP++ writes one VCF record per SNP inside a bubble, and so does this;
                // emitting one multi-base REF/ALT would be discarded by the clean-SNV gate.
```

## Coverage ceiling: wrong tail

Source: `include/caps_caller.h`, lines 3941-3942 before the documentation revision.

Original comment:

```
                // H. Every threshold raised before this was filtering the
                // wrong tail, which is why they all cost recall for nothing.
```

What the code now says:

```
                // H. Low-coverage thresholds therefore filter the wrong tail.
```

## Size bound on the emitted record

Source: `include/caps_caller.h`, lines 4000-4020 before the documentation revision.

Original comment:

```
                    // The bound applied during the superbubble walk measures
                    // the PATH length difference, and the caller's own indels
                    // already max out around 20 bp -- so it fires on almost
                    // nothing (dropped=3 of 117). The false positives are
                    // 20-37 bp only AFTER the lift's left-alignment expands
                    // them, so the discriminating quantity is the size of the
                    // record we actually emit.
                    //
                    // Measured on the tetraploid window: capping at 15 bp
                    // removes 12 of 16 indel false positives and ZERO true
                    // positives, taking precision 0.746 -> 0.922 and F1
                    // 0.537 -> 0.577, past DiscoSNP++'s 0.553. On full chr20
                    // the same split holds: true indels are 1-8 bp while the
                    // false ones form a second mode at 21-32 bp with
                    // essentially no true positives in it.
                    // Event size = the length DIFFERENCE between the two
                    // paths. (Bounding max(path1,path2) instead was tried and
                    // was catastrophic -- F1 0.537 -> 0.035 -- because the
                    // paths carry the whole bubble, flank and closure
                    // included, so a 15 bp cap on their total length rejects
                    // essentially every real event.)
```

What the code now says:

```
                    // The bound during the superbubble walk measures the PATH length difference and
                    // rarely fires; the false positives become 20-37 bp only after the lift's
                    // left-alignment, so the bound is applied to the emitted record. On the
                    // tetraploid window a 15 bp cap removed 12 of 16 indel false positives and no true
                    // positives (F1 0.537 -> 0.577); at full chr20 true indels are 1-8 bp and the false
                    // ones form a second mode at 21-32 bp.
                    // Event size is the length DIFFERENCE between the two paths. Bounding the path
                    // lengths themselves rejects almost every real event (F1 0.537 -> 0.035), because
                    // the paths carry the whole bubble, flank and closure included.
```

## Long-allele rejection (rejected)

Source: `include/caps_caller.h`, lines 4027-4036 before the documentation revision.

Original comment:

```
                    // REFUTED: rejecting records where BOTH alleles are long.
                    // The long-FP records look like 15->17, 17->18, 18->20 --
                    // paired-long alleles with a 1-2 bp difference -- so
                    // bounding the SHORTER allele looked like the discriminator.
                    // Measured: F1 0.537 -> 0.035. Real indels have both alleles
                    // long too, because path1/path2 carry the entire bubble
                    // (flank + event + closure), not just the event. Length of
                    // the allele strings therefore carries no signal here at
                    // all -- the same reason bounding max(path1,path2) failed
                    // identically.
```

What the code now says:

```
                    // Allele-string length is not used as a filter: real indels also have long
                    // alleles (path1/path2 carry flank, event and closure), and bounding the shorter
                    // allele gave the same collapse in F1 as bounding the path lengths.
```

## REF is the assembled base

Source: `include/caps_caller.h`, lines 4051-4055 before the documentation revision.

Original comment:

```
                        // REF must be the ASSEMBLED base; ALT is whichever path
                        // disagrees with it. Emitting path1 as REF regardless was
                        // wrong whenever the assembly carried path2, and turned
                        // every such call into a false positive (measured: FP
                        // 17 -> 182 before this fix).
```

What the code now says:

```
                        // REF must be the ASSEMBLED base; ALT is whichever path disagrees with it.
                        // Taking path1 as REF regardless is wrong whenever the assembly carries path2.
```

## Depth guard

Source: `include/caps_caller.h`, lines 4242-4252 before the documentation revision.

Original comment:

```
    // DEPTH GUARD, restored to its original semantic. ARCS's own description
    // is "reject columns deeper than 2.5x the MEDIAN CANDIDATE DEPTH
    // (collapsed repeats)", but this port had written it as DHI*H*1.25,
    // substituting the haploid-depth estimate for the median. That silently
    // couples the guard to H, and after the collapse pass a legitimate pileup
    // is DEEPER than raw coverage (reads from several duplicate contigs now
    // stack on one), so an H-derived ceiling rejects real sites. Measured:
    // forcing H to 10/20/40/80 moved average SNV F1 to 0.346/0.858/0.879/0.880
    // -- i.e. the "better" H values were simply the ones that disabled this
    // gate. Computing the median candidate depth restores the intended
    // behaviour and self-calibrates to any coverage, with no H dependence.
```

What the code now says:

```
    // DEPTH GUARD: reject columns deeper than 2.5x the MEDIAN CANDIDATE DEPTH
    // (collapsed repeats), as in the ARCS caller this pileup follows. It uses the
    // median of the candidates, not a ceiling derived from H, because after the
    // collapse pass a legitimate pileup is deeper than raw coverage (reads from
    // several duplicate contigs stack on one); an H-derived ceiling rejects real
    // sites. The median self-calibrates to any coverage.
```

## k-mer sanity threshold

Source: `include/caps_caller.h`, lines 4300-4311 before the documentation revision.

Original comment:

```
        // K-MER SANITY, re-derived after the H estimator was corrected.
        // The gate exists to reject REPEATS, whose flanks occur at multiples
        // of the homozygous depth (2H, 3H...). A legitimate homozygous-context
        // flank occurs at ~1.0H, so a ceiling of KHI=1.1H leaves only 10%
        // margin and rejects real sites on noise alone. That constant was
        // calibrated against the OLD, upward-biased H estimator (which
        // returned ~4x the true peak and so supplied the margin by accident);
        // correcting the estimator without re-deriving its companion threshold
        // would be the actual error. The repeat boundary is at 2H, so the
        // decision point is the midpoint 1.5H -- "closer to two copies than to
        // one" -- which is a derivation from the gate's stated purpose, not a
        // fitted value. KHI itself is left untouched for the uncollapsed path.
```

What the code now says:

```
        // K-MER SANITY. The gate rejects REPEATS, whose flanks occur at multiples of the
        // homozygous depth (2H, 3H...). A legitimate flank occurs at about 1.0H, so a
        // ceiling of 1.1H would reject real sites on noise. The repeat boundary is 2H, so
        // the decision point is the midpoint 1.5H ("closer to two copies than to one"),
        // derived from the gate's purpose. KHI is kept for the uncollapsed path.
```

## Two substrates

Source: `include/caps_caller.h`, lines 4366-4393 before the documentation revision.

Original comment:

```
    // DUAL VIEW. The collapse pass is what makes the SNV pileup work (it puts
    // both alleles' reads in one frame), but it is actively HARMFUL to the
    // bubble passes below, which need the two haplotype copies to still exist
    // as SEPARATE contigs -- collapsing deletes exactly the pair a bubble is
    // made of. Measured on HG002 held-out windows: indel F1 fell 0.357 -> 0.293
    // (r3) / 0.159 (na) once collapse was enabled.
    // So: SNVs are called on the collapsed+remapped substrate (above), and the
    // bubble passes run on the ORIGINAL uncollapsed contigs (below). One
    // assembly, two views, each used where it is correct.
    // A/B TESTED, both directions, same reads (HG002 r2): bubbles on the
    // COLLAPSED substrate score INDEL F1 0.422; on the ORIGINAL uncollapsed
    // contigs 0.349. My dual-view hypothesis -- that bubbles need the
    // duplicate pair collapse deletes -- is REFUTED: collapse helps the bubble
    // passes too, because a bubble needs the two HAPLOTYPES to differ, not the
    // same haplotype duplicated. Default is therefore the collapsed substrate;
    // CAPS_BUBBLE_UNCOLLAPSED=1 restores the other view for re-measurement.
    // MEASURED TENSION (HG002 r2, same reads, dup_frac swept):
    //     dup=0.45 -> SNV F1 0.882, INDEL F1 0.217
    //     dup=0.80 -> SNV F1 0.821, INDEL F1 0.422
    //     no collapse -> SNV F1 0.856, INDEL F1 0.349
    // The two variant classes want DIFFERENT amounts of collapse, and neither
    // setting is good for both. The SNV pileup wants both haplotypes' reads in
    // ONE frame (aggressive collapse); a bubble needs the two haplotypes to
    // still exist as TWO comparable contigs (mild collapse -- enough to remove
    // same-haplotype duplicates, not enough to merge the haplotypes).
    // So build TWO substrates from the one assembly, each at its own optimum.
    // This is the dual view done correctly; an earlier attempt paired the
    // pileup view with the UNCOLLAPSED contigs and was refuted (0.349).
```

What the code now says:

```
    //
    // TWO SUBSTRATES FROM ONE ASSEMBLY. SNVs and indels need different amounts of
    // contig collapse. The SNV pileup wants both haplotypes' reads in ONE frame
    // (aggressive collapse, dup_frac 0.45). A bubble needs the two haplotypes to
    // remain TWO comparable contigs, with only same-haplotype duplicates removed (mild
    // collapse, bub_dup 0.92). On the tuning window (HG002 r2):
    //     dup=0.45    -> SNV F1 0.882, INDEL F1 0.217
    //     dup=0.80    -> SNV F1 0.821, INDEL F1 0.422
    //     no collapse -> SNV F1 0.856, INDEL F1 0.349
    // so the SNV pileup and the bubble passes each run on their own substrate.
    // CAPS_BUBBLE_UNCOLLAPSED=1 runs the bubble passes on the uncollapsed contigs
    // instead (re-measurement only).
```

## indel_pass timing start

Source: `include/caps_caller.h`, lines 4394-4397 before the documentation revision.

Original comment:

```
    // Start the indel_pass accounting HERE, before the second build_substrate.
    // It was set after that call, so a ~52 s substrate build sat inside the
    // stage total but outside every lap -- and a timer labelled
    // "2nd build_substrate + setup" was measuring only the setup after it.
```

What the code now says:

```
    // indel_pass timing starts here, so the second build_substrate is included.
```

## Local read realignment (not used)

Source: `include/caps_caller.h`, lines 4421-4430 before the documentation revision.

Original comment:

```
    // A per-candidate "local read realignment" check was drafted here and
    // reverted before being wired in (2026-09-03): tracing through what it
    // would actually test showed it collapses into one of two things already
    // on record -- either DROP-COV (per-contig read coverage, already
    // implemented) if scoped to each contig separately, or the extended-
    // context matching already MEASURED AND REFUTED (INDEL_PRECISION_ROOT_CAUSE.md
    // Group A: our contigs are too short to supply 60+ bp of matching context)
    // if scoped to verifying rc_/ac are truly the same locus. Recorded here,
    // not silently dropped -- this is a real negative
    // result, not an abandoned draft.
```

What the code now says:

```
    // No per-candidate local read realignment is applied: scoped per contig it
    // reduces to DROP-COV (per-contig read coverage, implemented), and scoped to
    // cross-locus verification it needs 60+ bp of matching context that these
    // contigs are too short to supply.
```

## kidx layout

Source: `include/caps_caller.h`, lines 4468-4485 before the documentation revision.

Original comment:

```
        // KIDX AS SORTED FLAT ARRAY (2026-09-03), same overhead removal as
        // ridx/kc above -- the last remaining hashmap-of-vectors, and the
        // measured dominant contributor to the indel pass's RAM (this whole
        // pass was 62% of RSS growth at real mid-scale measurement, and
        // kidx is its only structure sized by TOTAL CONTIG BASES with a
        // per-key heap vector). Unlike ridx, no stride subsampling and no
        // per-key cap are needed: kidx indexes contig bases, not
        // reads x coverage, so total size is bounded by total contig
        // length, and any one key's true occurrence count is exactly what
        // the MAXOCC filter below already needs -- capping would risk
        // silently changing which keys pass that filter, so this reproduces
        // the ORIGINAL uncapped semantics exactly, just laid out flat.
        // 16 BYTES, NOT 24 -- same packing already applied to the pcluster
        // anchor record. {uint64 kmer; uint32 ci, pos; uint8 orient} pads to
        // 24 B, and kidx holds one entry per contig k-mer position (~106M
        // here), so 2.55 GB plus the parallel stable_sort's full copy. That
        // pair is what takes the run's high-water mark from 8.3 GB to 10.5 GB.
        // A canonical 25-mer is 50 bits, so the orientation bit rides in bit 63.
```

What the code now says:

```
        // KIDX AS A SORTED FLAT ARRAY (one record per contig k-mer position), not a
        // hash map of vectors; it is the indel pass's largest structure. No stride
        // subsampling or per-key cap: each key's true occurrence count is what the
        // MAXOCC filter below needs, so capping could change which keys pass.
        // 16-BYTE RECORDS: {uint64 kmer; uint32 ci, pos; uint8 orient} would pad to
        // 24 B. A canonical 25-mer is 50 bits, so the orientation bit rides in bit 63.
```

## kidx built in place

Source: `include/caps_caller.h`, lines 4493-4507 before the documentation revision.

Original comment:

```
        // BUILT IN PLACE, like the seed index above and for the same reason:
        // kidx is 16 B x ~106M = 1.70 GB and handing it to
        // __gnu_parallel::stable_sort allocates a full second copy. With the
        // seed index and pv fixed, THIS was the run's remaining high-water
        // mark -- measured +2394 MB in one step at `indel setup`.
        //
        // Entries are bucketed on the way in by the top bits of the RAW key
        // (not of mix(key), because this index is binary-searched by key and
        // must stay in key order), then each bucket is sorted in place with an
        // explicit (key, ci, pos) comparator. That comparator states what the
        // stable sort was actually relying on -- the fill is contig-major,
        // position-minor, so equal keys were already in (ci,pos) order -- so
        // the sequence is identical and stability is no longer needed.
        // Unlike the seed index there is no per-key cap here, so buckets stay
        // contiguous and no compaction pass is required.
```

What the code now says:

```
        // BUILT IN PLACE, like the seed index, so no sort temporary doubles it. Entries
        // are bucketed by the top bits of the RAW key (not of mix(key): this index is
        // binary-searched by key and must stay in key order), and each bucket is sorted
        // in place with an explicit (key, ci, pos) comparator. The fill is contig-major,
        // position-minor, so this is the order a stable sort would give. With no per-key
        // cap, buckets stay contiguous and no compaction is needed.
```

## Anchor multiplicity

Source: `include/caps_caller.h`, lines 4546-4553 before the documentation revision.

Original comment:

```
        // ANCHOR MULTIPLICITY. The rule used to be `occ.size() != 2 -> skip`,
        // i.e. an anchor was only usable if its 25-mer occurred EXACTLY twice
        // in the whole contig set. A genuine hap1/hap2 pair whose k-mer also
        // appears in any third fragment was therefore discarded silently.
        // Now every cross-contig PAIR among a small number of occurrences is
        // tried; the bubble extractor itself rejects pairs that do not
        // diverge-and-reconverge, so admitting more candidates costs
        // specificity only where the geometry genuinely looks like a bubble.
```

What the code now says:

```
        // ANCHOR MULTIPLICITY. Every cross-contig PAIR among a small number of
        // occurrences is tried, not only anchors that occur exactly twice: a real
        // hap1/hap2 pair whose k-mer also occurs in a third fragment must not be lost. The
        // bubble extractor rejects pairs that do not diverge and reconverge, so extra
        // candidates cost specificity only where the geometry looks like a bubble.
```

## indel_pass lap timing

Source: `include/caps_caller.h`, lines 4591-4594 before the documentation revision.

Original comment:

```
        // the second build_substrate -- which runs before the first lap --
        // landed in the stage total but in no timer, and 128 s of a 323 s
        // stage looked "unaccounted" through four wrong guesses at where it
        // was (im emit loop 0.02 s, lvotes dead code, 6b2 25.6 s, XSNV 14.3 s).
```

What the code now says:

```
        // the second build_substrate, which runs before the first lap, is counted.
```

## Keep only runs that do work

Source: `include/caps_caller.h`, lines 4606-4612 before the documentation revision.

Original comment:

```
        // MEASURED FAILURE, first attempt: threading over ALL runs gave 2.9%
        // (1153.8 -> 1120.3 s) at 482% CPU of a possible 1200%. There are
        // 86,487,848 runs and the body bails immediately on
        // `occ.size() < 2 || > MAXOCC`, so the overwhelming majority are
        // no-ops -- ~1.35M dynamic scheduling handoffs for nothing, which
        // swamps the real work. Parallelising the outer loop was the wrong
        // GRANULARITY, not the wrong idea.
```

What the code now says:

```
        // Most runs are no-ops (the body returns at once when occ.size() < 2 or
        // > MAXOCC), so threading over ALL runs spends its time on scheduling.
```

## Order-preserving merge

Source: `include/caps_caller.h`, lines 4636-4639 before the documentation revision.

Original comment:

```
            // replays those semantics explicitly -- which also fixes the
            // measured content difference from the first attempt (5 records),
            // where the merge relied on thread-visit order and
            // schedule(dynamic) does not match serial order.
```

What the code now says:

```
            // replays those semantics explicitly; relying on thread-visit order would not
            // match the serial result, because schedule(dynamic) does not preserve order.
```

## Multi-start extraction

Source: `include/caps_caller.h`, lines 4697-4706 before the documentation revision.

Original comment:

```
            // nearest difference. Measured on the 9 indels DiscoSNP++ finds and
            // we miss (working note INDEL_LOSS_SKELETAL.md §6): usable shared anchors
            // are ABUNDANT for these loci (49-1133 per contig pair), but the
            // nearest one lies 213-135,453 bp from the event -- so the walk
            // always terminates on some other difference first and the event is
            // never examined. Anchors are neither scarce nor over-repeated;
            // they are simply in the wrong PLACE.
            // Fix: after taking the first bubble, continue walking past it and
            // extract again, so one anchor pair can yield several events along
            // its length instead of only the closest.
```

What the code now says:

```
            // nearest difference. For the missed indels on the tuning window, usable shared
            // anchors were abundant but far from the event (213-135,453 bp), so the walk ended
            // on another difference first. After taking a bubble the walk continues past it
            // and extracts again, so one anchor pair can yield several events along its
            // length.
```

## Closing anchor

Source: `include/caps_caller.h`, lines 4734-4745 before the documentation revision.

Original comment:

```
            // THE structural difference from DiscoSNP++. Its bubble leaves a
            // shared graph node and must RE-CONVERGE on another shared node,
            // so both ends of both paths are k-mers that exist in the reads.
            // Ours opened on one 25-mer anchor and closed on nothing stronger
            // than 15 bases of contig-vs-contig identity -- a pair of paralogs
            // sharing one anchor can satisfy that by luck, which is why our
            // indel precision sat at 0.729 against their 0.910 while our
            // recall was already comparable (0.497 vs 0.522).
            // Require the re-convergence point to be a real closing anchor: a
            // full BK-mer identical in both contigs AND present in the reads.
            // This is the both-ends-anchored property obtained on our own
            // contig set, with no graph.
```

What the code now says:

```
            // DiscoSNP++'s bubble leaves a shared graph node and must RE-CONVERGE on another
            // shared node, so both ends of both paths are k-mers present in the reads. An
            // opening 25-mer anchor plus 15 bases of contig identity at the close can be
            // satisfied by chance by two paralogs, which lowers indel precision. So the
            // re-convergence point must be a real closing anchor: a full BK-mer identical in
            // both contigs AND present in the reads. This is the both-ends-anchored property,
            // obtained on the contig set without a graph.
```

## Homopolymer-aware closing anchor

Source: `include/caps_caller.h`, lines 4759-4763 before the documentation revision.

Original comment:

```
                // exactly that amount. Every one of the 9 truth indels
                // DiscoSNP++ finds and we miss is such an event
                // (working note INDEL_LOSS_SKELETAL.md). Slide the alt side by up to
                // the event length and re-test; a genuine run-length change
                // closes exactly there.
```

What the code now says:

```
                // exactly that amount (all of the 1 bp indels missed on the tuning window are
                // such events). The alt side is slid by up to the event length and re-tested; a
                // real run-length change closes exactly there.
```

## SNV-tolerant closing anchor

Source: `include/caps_caller.h`, lines 4776-4795 before the documentation revision.

Original comment:

```
                // MEASURED MOTIVATION: instrumenting every drop reason on the
                // tetraploid benchmark shows this test is not merely dominant,
                // it is nearly the whole story -- 2,512 of 2,865 candidate
                // drops (88%) are DROP-CLOSE, against 323 DROP-COV, 21
                // DROP-ANCH and 9 DROP-JUNC. A byte-identical BK=25-mer is
                // required on BOTH contigs, so ONE heterozygous SNV anywhere
                // in those 25 bases destroys an otherwise-valid bubble, and at
                // ~1 het/1.3kb in real human data that is common.
                //
                // Why tolerance here is principled rather than just looser:
                // the very NEXT test requires the closing anchor to be
                // observed in the READS (kcount >= MC). DiscoSNP++'s design is
                // exactly this -- propose permissively from the graph, then
                // let kissreads2 reject against reads -- and the standing
                // finding of this project (working note HOW_DISCOSNP_WINS.md sec 4) is
                // that OUR read-validation has nothing to reject because our
                // proposals are already read-derived and strict. Loosening the
                // proposal is what gives the read test something to do.
                // Tolerance reuses the frozen HDMAX (=2) rather than
                // introducing a new tunable constant.
```

What the code now says:

```
                // On the tetraploid benchmark 88% of candidate drops were at the closing-anchor
                // test: a byte-identical BK=25-mer is required on BOTH contigs, so one
                // heterozygous SNV in those 25 bases (common at about 1 het per 1.3 kb in human
                // data) destroys an otherwise valid bubble. Tolerance is justified because the next
                // test requires the closing anchor to be observed in the READS (kcount >= MC):
                // propose permissively, then validate against reads, as DiscoSNP++ does with
                // kissreads2. Tolerance reuses the frozen HDMAX (=2) instead of a new constant.
```

## Unique closing anchor

Source: `include/caps_caller.h`, lines 4813-4818 before the documentation revision.

Original comment:

```
                // Measured NEGATIVE and therefore opt-in: requiring the
                // closing anchor to be unique as well costs more than it buys
                // (five-window mean indel F1 0.5976 -> 0.5932, HG003
                // 0.672 -> 0.640). A repeated closing k-mer still closes a real
                // bubble; demanding uniqueness at BOTH ends over-constrains a
                // fragmented contig set. Kept for re-measurement elsewhere.
```

What the code now says:

```
                // OPT-IN: requiring a unique closing anchor lowered indel F1 (five-window mean
                // 0.5976 -> 0.5932); a repeated closing k-mer still closes a real bubble.
```

## Extended agreement

Source: `include/caps_caller.h`, lines 4828-4837 before the documentation revision.

Original comment:

```
            // EXTENDED AGREEMENT. Read-support tests cannot separate our false
            // bubbles from true ones: every candidate is built FROM contigs,
            // which are built from reads, so all of it is read-supported by
            // construction (measured -- requiring full junction coherence
            // rather than any-point support moved indel F1 by +0.001).
            // What actually separates a real haplotype pair from a paralog
            // pair that happens to share one 25-mer is HOW FAR the two contigs
            // agree either side of the event: true haplotypes agree over
            // hundreds of bases, paralogs diverge quickly. Verify a long span
            // rather than the 15 bp flanks the bubble extractor checks.
```

What the code now says:

```
            // EXTENDED AGREEMENT. Read-support tests cannot separate false bubbles here,
            // because every candidate is built from contigs built from reads. What separates
            // a real haplotype pair from a paralog pair sharing one 25-mer is how far the two
            // contigs agree on either side of the event: haplotypes agree over hundreds of
            // bases, paralogs diverge quickly. This verifies a long span rather than the
            // 15 bp flanks the bubble extractor checks.
```

## Extended agreement: measurement

Source: `include/caps_caller.h`, lines 4839-4845 before the documentation revision.

Original comment:

```
                // REFUTED and therefore OFF by default (EXT=0). Swept across
                // five windows: mean indel F1 0.5930 (off), 0.5864 (60bp/tol3),
                // 0.5644 (120/4), lower still at 200/6. The idea is sound in
                // principle but our contigs are too SHORT to supply the context
                // -- demanding 60+60 bp either side rejects true bubbles that
                // simply run out of contig. Kept behind CAPS_BUB_EXT for
                // re-measurement on a longer-contig substrate.
```

What the code now says:

```
                // OFF by default (EXT=0, CAPS_BUB_EXT to enable): the contigs are too short to
                // supply 60+ bp either side, so it rejected true bubbles (five-window mean indel
                // F1 0.5930 -> 0.5864 at 60 bp).
```

## Targeted paralog filter

Source: `include/caps_caller.h`, lines 4846-4855 before the documentation revision.

Original comment:

```
                // TARGETED PARALOG FILTER. Extended agreement was measured
                // NET NEGATIVE when applied to ALL events (five-window mean
                // 0.5930 -> 0.5864 at 60bp) -- because most events are 1-2bp
                // and our contigs are too short to supply the context, so it
                // rejected true short indels.
                // But the FP burden is concentrated in LONG events: pooled over
                // all 8 evaluations, len>=3 has TP=32 FP=26 (precision 0.55)
                // against len==1's TP=97 FP=20 (0.83). A long indel also
                // implies long contigs, so the context IS available exactly
                // where it is needed. Apply the check only there.
```

What the code now says:

```
                // TARGETED PARALOG FILTER: extended agreement applied only to long events. The
                // false-positive burden is concentrated at len >= 3 (pooled over 8 evaluations:
                // precision 0.55 against 0.83 at len 1), and a long indel implies longer contigs,
                // so the context is available where it is needed.
```

## Targeted paralog filter: measurement

Source: `include/caps_caller.h`, lines 4857-4860 before the documentation revision.

Original comment:

```
                // MEASURED NEUTRAL at 40bp on r2/r4, so OFF by default and
                // opt-in via CAPS_LONG_EXT. Paralogs agree over 40bp easily, so
                // the span needed to discriminate exceeds what a ~335bp contig
                // can supply even for long events.
```

What the code now says:

```
                // OFF by default (CAPS_LONG_EXT): neutral at 40 bp, because paralogs agree over
                // 40 bp easily and longer spans exceed the contig lengths.
```

## k-mer frequency band

Source: `include/caps_caller.h`, lines 5142-5153 before the documentation revision.

Original comment:

```
                // This is neither a context test, a read-support test, nor a
                // matching constraint -- all of which failed -- but a frequency
                // test, and it reuses the frozen KHI against the corrected H.
                // REFUTED, so opt-in only. Swept across five windows at
                // 1.0xH: 0.624/0.692/0.613/0.566 against a baseline of
                // 0.655/0.699/0.635/0.569 -- consistently WORSE, i.e. the cap
                // removes true indels rather than repeat-derived ones.
                // Diagnosis: our alt-junction k-mer counts already sit near the
                // heterozygous level because the candidate came from a contig
                // built out of those very reads, so a frequency band has no
                // repeat-specific signal left to cut on. Same underlying reason
                // the read-support tests were vacuous here.
```

What the code now says:

```
                // It reuses the frozen KHI against H.
                // OPT-IN only: across five windows at 1.0xH it lowered indel F1 on every window
                // (it removed true indels), because alt-junction k-mer counts already sit near the
                // heterozygous level when the candidate comes from a contig built from those
                // reads.
```

## Positional-clustering indel channel

Source: `include/caps_caller.h`, lines 5205-5223 before the documentation revision.

Original comment:

```
        // The eBWT2SNP principle, ported. Their 99.13% SNP precision comes from
        // never CONSTRUCTING a candidate: reads covering one genome position
        // share a right context, suffix-sorting groups them into a cluster, and
        // the differing allele sits immediately before it -- so both haplotype
        // fragments are lifted straight out of the reads. Their paper leaves
        // the INDEL case explicitly unimplemented ("extract the left-context
        // and perform a local alignment ... future work"); this is that.
        //
        // Why it should succeed where every filter failed: our false bubbles
        // come from SPLICING one contig's flank onto another's allele, and no
        // read-support test, context test or matching can detect a chimera
        // whose pieces are all read-derived (see working note INDEL_PRECISION_ROOT_CAUSE.md).
        // Here the evidence is a real read carrying a real gap against the
        // contig; the contig only supplies coordinates and the REF allele.
        // Anchoring on the RIGHT context is also what makes this different from
        // the earlier read-gap scan that found nothing: that one asked a read
        // where it placed best, and an alt-haplotype read places perfectly on
        // its OWN contig. Here every read covering an anchor is compared to the
        // SAME contig, so alt-haplotype reads must reveal their gap.
```

What the code now says:

```
        // After eBWT2SNP: reads covering one genome position share a right context, so
        // grouping reads by a right-context anchor puts the differing allele immediately
        // before it, and both haplotype fragments come from the reads rather than from a
        // constructed candidate. The eBWT2SNP paper leaves the indel case unimplemented
        // ("extract the left-context and perform a local alignment ... future work"); this
        // channel implements it. The evidence is a real read carrying a real gap against
        // the contig; the contig supplies only coordinates and the REF allele. Bubble
        // false positives come from splicing one contig's flank onto another's allele,
        // which read-support tests cannot detect when every piece is read-derived. Every
        // read covering an anchor is compared with the SAME contig (not placed where it
        // fits best), so alt-haplotype reads reveal their gap.
```

## Positional-clustering substrate

Source: `include/caps_caller.h`, lines 5227-5239 before the documentation revision.

Original comment:

```
            // SUBSTRATE CHOICE. The bubble passes need the two haplotypes to
            // survive as separate contigs, which is why they run on the mildly
            // collapsed substrate (dup=0.92). This channel has no such need --
            // it compares READS against ONE contig -- so it can use the
            // aggressively collapsed pileup substrate (dup=0.45), where reads
            // stack far better. Selectable so the two can be compared.
            // MEASURED: the mildly-collapsed BUBBLE substrate is the right one
            // here, despite this channel not needing haplotype separation the
            // way the bubble passes do. On the aggressively collapsed pileup
            // substrate (dup=0.45) indel F1 falls 0.655 -> 0.526 on r2, because
            // that collapse merges the two haplotypes into ONE contig: a read
            // carrying the indel then often matches a contig that already
            // contains it, and no gap appears at all.
```

What the code now says:

```
            // SUBSTRATE CHOICE: the mildly collapsed bubble substrate (dup=0.92). On the
            // aggressively collapsed pileup substrate (dup=0.45, CAPS_PCLUSTER_COLSUB) the two
            // haplotypes merge into one contig, a read carrying the indel often matches a
            // contig that already contains it, and no gap appears (indel F1 0.655 -> 0.526 on
            // the tuning window).
```

## pkidx as a flat sorted map

Source: `include/caps_caller.h`, lines 5249-5259 before the documentation revision.

Original comment:

```
            // FLAT SORTED MAP, NOT A HASH MAP. `pv` is already sorted by key
            // and the selection walk below emits entries in that order, so the
            // result is inherently sorted -- it was being poured into an
            // unordered_map (21.6M nodes, one allocation each) and then
            // ITERATED, which is pointer-chasing through 21.6M nodes.
            //
            // Same reasoning as the collapse_contigs container swap that took
            // the run 207.6 -> 170.9 s: the semantics here need ordered append
            // plus membership, both of which a sorted array does better.
            // Lookups are binary search, but they are guarded by the Bloom
            // bitset so only a small fraction ever reach it.
```

What the code now says:

```
            // FLAT SORTED MAP, NOT A HASH MAP. `pv` is sorted by key and the selection walk
            // emits entries in that order, so the result is inherently sorted; ordered append
            // plus membership is what the semantics need.
```

## One-shot hash directory

Source: `include/caps_caller.h`, lines 5276-5294 before the documentation revision.

Original comment:

```
                // The indel-pass read scan probes this map ~340M times (4M
                // reads x ~85 anchors). It did so via a two-hash Bloom filter
                // -- 2 random touches of a 43 MB bitset -- and then, for every
                // survivor, index_of, a binary search over 21.6M keys: ~25 more
                // dependent probes. That pair was 6.66 s of a 73.4 s run.
                //
                // This resolves the same query in TWO touches whether it hits
                // or misses, and needs no filter in front of it: one probe into
                // a 2^23-entry offset table, then one cache line of a bucket
                // array holding (key, index) pairs -- ~2.6 entries per bucket,
                // 32 B, so the scan almost never leaves that line. It strictly
                // dominates the old path: same cost on a miss, and on a hit it
                // also returns the index instead of starting a binary search.
                //
                // Buckets are taken from mix(key), not the key itself, because
                // raw canonical k-mers are strongly non-uniform and a
                // prefix-keyed directory would concentrate into a few huge
                // buckets. `k` is NOT reordered -- the anchor loop's iteration
                // order over pkidx is load-bearing and stays exactly as it was.
```

What the code now says:

```
                // The indel-pass read scan probes this map about 85 times per read. This resolves
                // a query in TWO memory touches whether it hits or misses: one probe into a
                // 2^23-entry offset table, then one cache line of a bucket array of (key, index)
                // pairs (about 2.6 entries per bucket). Buckets come from mix(key), because raw
                // canonical k-mers are strongly non-uniform. `k` is NOT reordered: the anchor
                // loop's iteration order over pkidx is load-bearing.
```

## Anchor selection by sorting

Source: `include/caps_caller.h`, lines 5325-5349 before the documentation revision.

Original comment:

```
                // MEASURED: this block is ~92% of indel_pass (480 s of 522 s on
                // a 4M-read subset), which is itself 84% of the archive-path
                // caller -- and it yields ~10% of the indels. The cost is
                // structural, not incidental: three unordered_maps taking ~65M
                // insertions each (reserved at 2M, so they rehash repeatedly),
                // followed by two full traversals doing erase() DURING
                // iteration with extra hash lookups per element.
                //
                // Identical semantics, computed by sorting instead:
                //   * count occurrences of the CANONICAL 25-mer across BOTH
                //     orientations (counting forward-canonical only measured a
                //     real precision loss, r2 indel P 0.760 -> 0.586);
                //   * the FIRST occurrence supplies (contig, pos, orient) --
                //     the original's `if (++pcount[cn] == 1)` -- and contigs
                //     are visited in increasing ci, so "first in sort order"
                //     IS "first seen";
                //   * survivors are the keys with count 1 that are
                //     forward-oriented. (The original ran two erase passes; the
                //     second, `pcount != 1`, subsumes the first at the default
                //     PUNIQ=1. Both are applied here so any PUNIQ>1 behaves as
                //     before.)
                //
                // Build is the count -> prefix-sum -> parallel-fill-by-index
                // pattern already used for kidx above, so entry order matches
                // the serial push order exactly.
```

What the code now says:

```
                // Computed by sorting rather than with hash maps (this block dominated indel_pass
                // time). Semantics:
                //   * count occurrences of the CANONICAL 25-mer across BOTH orientations
                //     (counting forward-canonical only lowered indel precision);
                //   * the FIRST occurrence supplies (contig, pos, orient); contigs are visited
                //     in increasing ci, so "first in sort order" IS "first seen";
                //   * survivors are the keys with count 1 that are forward-oriented. Both erase
                //     conditions are applied, so PUNIQ > 1 keeps its meaning (at the default
                //     PUNIQ=1 the second subsumes the first).
                // The build is count -> prefix-sum -> parallel fill by index, as for kidx, so entry
                // order matches the serial push order.
```

## pv record layout

Source: `include/caps_caller.h`, lines 5352-5362 before the documentation revision.

Original comment:

```
                // {uint64 k; uint32 ci, pos; uint8 fwd} pads to 24 B, and this
                // array holds one entry per contig k-mer position -- 106.4M on
                // the 4M-read subset, so 2.55 GB, which the parallel stable_sort
                // then copies. That pair is the run's memory PEAK: 13.54 GB,
                // measured, with the sort alone adding 3.3 GB.
                //
                // A canonical 25-mer occupies 50 bits, so the orientation flag
                // rides in bit 63 of the key and the record is exactly two
                // 8-byte words: 1.70 GB, and 1.70 GB for the sort's temporary.
                // The sort masks the flag off, so the ordering -- and the
                // stability the run-grouping below depends on -- is unchanged.
```

What the code now says:

```
                // {uint64 k; uint32 ci, pos; uint8 fwd} would pad to 24 B on an array with one
                // entry per contig k-mer position. A canonical 25-mer occupies 50 bits, so the
                // orientation flag rides in bit 63 and the record is two 8-byte words. The sort
                // masks the flag off, so the ordering and the stability the run grouping depends
                // on are unchanged.
```

## pv derived from kidx

Source: `include/caps_caller.h`, lines 5370-5390 before the documentation revision.

Original comment:

```
                // pv WAS kidx, REBUILT -- and it is the run's memory PEAK.
                //
                // Read the two builds side by side: both walk the same contigs
                // (pc_cd IS cdb unless CAPS_PCLUSTER_COLSUB substitutes the
                // collapsed SNV substrate), both take 25-mers with the same
                // Roll25 and the same ok() predicate, both pack
                // canon | (f<=r ? 0 : bit63), both store the same (ci,pos), and
                // both are stable-sorted by the masked key. They are the same
                // array, computed twice.
                //
                // kidx is built at the top of the indel pass and is still live
                // here (its last read is the cross-contig SNV loop, far below),
                // so the second copy is pure duplication: 106M x 16 B = 1.70 GB
                // for the array plus 1.70 GB for the parallel sort's temporary,
                // and that pair is exactly where VmHWM peaks.
                //
                // Deriving the anchors from kidx is EXACT, not approximate: the
                // fill order is the same (contig 0,1,2..., ascending position),
                // the sort is stable on the same key, so each run's first-seen
                // entry -- the one the occ==1 test keeps -- is the same entry.
                // CAPS_PV_REBUILD=1 forces the old path for A/B checking.
```

What the code now says:

```
                // pv is derived from kidx instead of being rebuilt. Both walk the same contigs
                // (pc_cd IS cdb unless CAPS_PCLUSTER_COLSUB substitutes the collapsed SNV
                // substrate), take 25-mers with the same Roll25 and ok() predicate, pack
                // canon | (f<=r ? 0 : bit63), store the same (ci,pos) and are stable-sorted by the
                // masked key, so they are the same array. kidx is still live here, so a rebuild
                // would duplicate it at the run's memory peak. The derivation is exact: same fill
                // order, same stable sort, so each run's first-seen entry (the one the occ==1 test
                // keeps) is the same entry. CAPS_PV_REBUILD=1 rebuilds it for A/B checking.
```

## pcov aliases cov

Source: `include/caps_caller.h`, lines 5496-5506 before the documentation revision.

Original comment:

```
            // pcov WAS cov, RECOMPUTED.
            //
            // This loop was character-for-character the coverage build at the
            // top of the indel pass -- same contigs, same read_cid/read_pos/
            // read_clip, same 60000 clamp -- because pc_cd IS cdb unless
            // CAPS_PCLUSTER_COLSUB puts the collapsed SNV substrate in its
            // place. So on every default run it built a second ~300 MB copy of
            // an array already in hand, with another ~590M serial increments.
            //
            // Alias when the two are literally the same object; keep the build
            // for the opt-in case, where they genuinely differ.
```

What the code now says:

```
            // pcov aliases cov when pc_cd IS cdb (the default): the coverage build would be
            // identical (same contigs, read_cid/read_pos/read_clip and 60000 clamp). It is
            // rebuilt only for the opt-in CAPS_PCLUSTER_COLSUB substrate.
```

## Only pkidx anchors stored

Source: `include/caps_caller.h`, lines 5528-5544 before the documentation revision.

Original comment:

```
            // ONLY ANCHORS PRESENT IN pkidx ARE EVER LOOKED UP (2026-09-03).
            // The sole consumer of this map is the `for (auto& kv : pkidx)`
            // loop below, whose first act is `rc_reads.find(kv.first)`. Any
            // key absent from pkidx is therefore hashed, stored, and never
            // read -- and that was the overwhelming majority of it, because
            // pkidx has just been filtered down to contig-UNIQUE anchors
            // while this loop indexes every 25-mer of every read, including
            // the error k-mers that dominate the distinct-key count (one
            // sequencing error spawns up to 25 novel k-mers). At full chr20
            // that is 12,604,917 reads x 84 offsets = 1.06 BILLION insertions
            // to serve a lookup set orders of magnitude smaller.
            //
            // Skipping non-pkidx keys is output-identical BY CONSTRUCTION,
            // not by measurement: a loop that only ever reads keys in pkidx
            // cannot observe the absence of keys that are not in pkidx. pkidx
            // is fully built and filtered above before this point, so the
            // membership test is well-defined here.
```

What the code now says:

```
            //
            // Only anchors present in pkidx are stored. The sole consumer is the loop over
            // pkidx below, which looks up rc_reads.find(kv.first), so a key absent from pkidx
            // is never read; most read 25-mers (including error k-mers) are such keys. This
            // is output-identical by construction, and pkidx is complete before this point.
```

## Bitset pre-filter: motivation

Source: `include/caps_caller.h`, lines 5546-5551 before the documentation revision.

Original comment:

```
            // MEASURED: this loop is `pcluster: scan + emit`, 251 s of a 600 s
            // stage (42%). It walks ~85 k-mers of each of 4M reads -- ~340M
            // iterations -- and every one does `pkidx.count(cn)`, a hash probe
            // into a 21.6M-entry map. The insertions after it are already
            // filtered ("provably dead otherwise"), so the PROBES are the cost,
            // not the work they admit.
```

What the code now says:

```
            // Each read probes about 85 k-mers against pkidx, and almost all probes miss.
```

## Parallel read scan fused with anchor work

Source: `include/caps_caller.h`, lines 5568-5600 before the documentation revision.

Original comment:

```
            // ── PARALLEL READ SCAN ──────────────────────────────────────────
            // 4M reads x ~85 k-mers, serial, and the bitset above already
            // removed ~99% of the hash traffic -- so what remains is a
            // compute-bound walk over independent reads. Measured 151.4 s of a
            // 339.6 s stage.
            //
            // Reads are independent, so each thread collects (key, read, pos)
            // into its own buffer and they are merged afterwards IN THREAD
            // ORDER, then by read index. That reproduces the serial insertion
            // order exactly, which matters because `vec.size() < 200` caps each
            // key -- a different arrival order would keep a different 200.
            // ── READ SCAN AND ANCHOR WORK, FUSED ────────────────────────────
            //
            // This used to be two passes with a CSR between them: scan every
            // read, buffer every (anchor, read, pos) hit per thread, merge the
            // buffers into rc_off/rc_val grouped by anchor, then walk the
            // anchors and do the indel comparison. That intermediate cost
            // 81.3M x 12 B of per-thread buffers plus a 650 MB value array plus
            // an 86 MB offset array -- measured as the run's largest single
            // allocation step at +1710 MB.
            //
            // The grouping bought exactly one thing: a 200-hits-per-anchor cap,
            // which is the only part of this that depends on the order hits
            // arrive. MEASURED, it never binds -- the longest row over
            // 21,605,670 anchors is 70 hits and at_cap is 0, at 30x coverage on
            // contig-UNIQUE 25-mers, which is what an anchor is. With the cap
            // inert, nothing downstream can observe the order: pvotes and panch
            // are set inserts and ploc keeps the smallest (contig,pos), so all
            // three are order-free and merge by union/min.
            //
            // So the comparison is done inline as each hit is found. Same work,
            // same per-thread accumulators, same merge -- with no intermediate
            // at all, and one less pass over 81.3M hits.
```

What the code now says:

```
            // ── PARALLEL READ SCAN, FUSED WITH THE ANCHOR WORK ──────────────────
            // Reads are independent, so each thread scans its own reads and does the indel
            // comparison inline as each anchor hit is found, with per-thread accumulators
            // merged afterwards. No per-anchor grouping is built. That grouping existed only
            // for a 200-hits-per-anchor cap, which never binds on contig-unique 25-mers at 30x
            // (longest row measured: 70 hits). With the cap inert nothing downstream depends
            // on hit order: pvotes and panch are set inserts and ploc keeps the smallest
            // (contig,pos), so all three merge by union or minimum.
```

## Parallel read placement

Source: `include/caps_caller.h`, lines 1015-1036 before the documentation revision.

Original comment:

```
    // PARALLEL (2026-09-03). This loop was the single largest phase in the
    // caller -- 738 s of the first full-chr20 run, and it runs TWICE (here for
    // the pileup substrate, again at the bubble substrate inside indel_pass).
    // It was serial only because it had never been profiled at scale: the
    // enclosing phase marker is called "ridx_build", but ridx is disabled by
    // default, so all of that time was actually being spent right here under
    // a name that pointed at dead code.
    //
    // Safe to parallelise by inspection, not by hope: iteration `o` writes
    // ONLY S.read_cid[o] / read_pos[o] / read_rc[o] / read_clip[o], which are
    // disjoint across iterations and pre-sized above; `idx` and `S.contigs`
    // are read-only once step 2 has finished (bound through a const ref below
    // so a stray insert cannot compile); every best_* accumulator is declared
    // inside the body and is therefore private. `placed`/`improved` are the
    // only shared scalars and are a reduction. Output is order-independent,
    // so the result is byte-identical to the serial version, not merely
    // equivalent -- which is the gate this change is held to.
    //
    // NOTE the prior negative result this does NOT contradict: OpenMP over
    // CONTIGS in the pileup measured slower (+11%, +13%, comment at the top of
    // run_variant_call). That loop has few, wildly unequal items. This one has
    // 12.6M uniform items; dynamic scheduling covers the tail.
```

What the code now says:

```
    // PARALLEL. This placement loop is the caller's largest phase and runs twice
    // (here for the pileup substrate, again for the bubble substrate inside
    // indel_pass). It is safe to parallelise: iteration `o` writes ONLY
    // S.read_cid[o] / read_pos[o] / read_rc[o] / read_clip[o], which are disjoint
    // across iterations and pre-sized above; `idx` and `S.contigs` are read-only
    // once step 2 has finished (bound through a const ref below so a stray insert
    // cannot compile); every best_* accumulator is declared inside the body and is
    // private. `placed`/`improved` are the only shared scalars and are a reduction.
    // The output does not depend on order, so it is byte-identical to a serial
    // run. (OpenMP over CONTIGS in the pileup is slower, because that loop has few,
    // very unequal items; this one has millions of uniform items and dynamic
    // scheduling covers the tail.)
```

## PG anchoring: idea

Source: `include/caps_caller.h`, lines 3959-3964 before the documentation revision.

Original comment:

```
    // We do not have to rely on that walk. The encoder has already assembled a
    // 63 Mb pseudogenome, and cd_in.contigs are long, real, well-anchored
    // sequence. Locating the bubble's flank k-mer inside one of them lets the
    // call be emitted in THAT contig's coordinates instead. DiscoSNP++ cannot
    // do this -- it has no assembly -- so this is our own structure paying for
    // itself again.
```

What the code now says:

```
    // The alternative is to locate the bubble's flank k-mer in one of the encoder's
    // long assembled contigs (cd_in.contigs) and emit the call in that contig's
    // coordinates, which DiscoSNP++ cannot do without an assembly.
```

## Maximum-weight matching: motivation

Source: `include/caps_caller.h`, lines 4955-4962 before the documentation revision.

Original comment:

```
        // That constraint is exactly what we lacked. A true heterozygous site
        // has exactly ONE alternate haplotype, but our anchor search lets one
        // locus pair with several alt contigs (and one alt contig serve many
        // loci), which is how paralogs and repeats enter as false bubbles.
        // Every precision idea tried before this was LOCAL -- more flanking
        // context, longer fragments, read substrings -- and all of them failed
        // for the same reason: our contigs are too short to supply context.
        // Matching is GLOBAL and needs no context at all.
```

What the code now says:

```
        // A true heterozygous site has exactly ONE alternate haplotype, but the anchor
        // search lets one locus pair with several alt contigs (and one alt contig serve
        // many loci), which is how paralogs and repeats can enter as false bubbles.
        // Matching is GLOBAL and needs no flanking context, which these short contigs
        // cannot supply.
```

## Maximum-weight matching: outcome

Source: `include/caps_caller.h`, lines 4975-4979 before the documentation revision.

Original comment:

```
        // Informative failure: the matching constraint assumes competing
        // pairings are the problem. They are not -- our false bubbles are
        // mostly UNCONTESTED, i.e. a locus with exactly one (wrong) partner,
        // which a matching cannot reject because there is nothing to compete
        // against it.
```

What the code now says:

```
        // The matching constraint assumes competing pairings are the problem, but false
        // bubbles here are mostly UNCONTESTED (a locus with exactly one wrong partner),
        // which a matching cannot reject.
```

## Read-level junction support: motivation

Source: `include/caps_caller.h`, lines 5056-5060 before the documentation revision.

Original comment:

```
            // Until now an indel was accepted on CONTIG-coverage proxies
            // (covwin / medcov) -- never on evidence that any READ actually
            // carries the alt allele. That is precisely the job DiscoSNP++
            // gives kissreads2, and precisely why its indel precision is
            // 0.88-0.97 against our 0.32-0.63.
```

What the code now says:

```
            // Contig-coverage proxies (covwin / medcov) do not show that any READ carries the
            // alt allele. DiscoSNP++ gives that job to kissreads2.
```

## Junction fragment length

Source: `include/caps_caller.h`, lines 5077-5082 before the documentation revision.

Original comment:

```
                // fragment genuinely rejects TRUE indels. Reason: our alt
                // fragment splices the REF contig's flanks onto the ALT allele,
                // so any nearby heterozygous SNV makes it a sequence no read
                // ever contains. EBWT2SNP does not have this problem because it
                // reads both haplotype fragments straight out of the eBWT
                // cluster instead of constructing one.
```

What the code now says:

```
                // fragment rejects TRUE indels: the alt fragment splices the REF contig's flanks
                // onto the ALT allele, so any nearby heterozygous SNV makes it a sequence no read
                // contains. eBWT2SNP avoids this by reading both haplotype fragments out of the
                // eBWT cluster instead of constructing one.
```

## Adaptive junction context

Source: `include/caps_caller.h`, lines 5089-5092 before the documentation revision.

Original comment:

```
                // not bad candidates being caught, it was good candidates near
                // contig ends being discarded unexamined. Our contigs are short,
                // so take as much context as EXISTS, down to a floor, and judge
                // the candidate on that.
```

What the code now says:

```
                // good candidates near contig ends being discarded unexamined. So the context is
                // as much as exists, down to a floor, and the candidate is judged on that.
```
