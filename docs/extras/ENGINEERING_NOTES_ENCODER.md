# Engineering notes: the encoder (src/encoder.cpp)

Long explanatory comments from `src/encoder.cpp` that were condensed in the source during the documentation revision. They are kept here word for word because they record measurements, rejected approaches and the reasoning behind design choices. The source keeps a short statement of each conclusion. `STAGE n` labels in these notes are the project's own numbering of the experiment that introduced a block. They are not stages of the program.


The quoted text is unedited, so it can name working documents (for example under `docs/_removable/`, or a `docs/*.md` that is not in this repository) and sibling code (ARCS, PgRC2) that are not shipped here. Those names are provenance only. Nothing in the code or the paper depends on them.

## Sweep seed width and the opt-in LOWCOV mode

Source: `src/encoder.cpp`, lines 636-668 before the documentation revision.

Original comment, verbatim:

```
    // Sweep seed width. 32 bases fills a uint64 exactly, which is why it was
    // fixed there -- but it also silently floors the shortest overlap the sweep
    // can ever see, and the measured trend (MINOV 50 -> 40 -> 32 giving 29.88M
    // -> 28.91M -> 28.28M) runs straight into that floor. A narrower seed costs
    // specificity (more candidates to verify) but lets round 2 pair reads that
    // overlap by less than 32, which after the division are still reads known
    // to tile well.
    // LOWCOV: an explicitly-requested low-coverage mode. It changes NO default
    // and is never auto-enabled, because the signal that would gate it
    // (leftover_frac) is only known AFTER round 1 -- the very stage this
    // controls -- so auto-detection would mean re-running chaining. Cheap
    // sample-based coverage estimators were tried and REJECTED: a 25-mer
    // frequency histogram gives mode=2 for every dataset regardless of true
    // depth, and unique-kmer fraction orders P. aeruginosa at 26x (0.929)
    // ABOVE Drosophila at 2.6x (0.919), because a fixed-size sample conflates
    // genome size with depth.
    //
    // Measured on Drosophila SRR40104920 (2.6x), sweeping SEEDW/MINOV together
    // with everything else fixed -- monotonic all the way to the floor:
    //     16/16  36,851,932     12/12  36,783,282
    //     10/10  36,558,158      8/8   36,103,845  (-2.38% vs tuned baseline)
    // 8 is the hard floor ("below this a seed is noise"). leftover_frac falls
    // 0.812 -> 0.785 and the second region shrinks by 6.3 M bases, i.e. more
    // reads chain instead of being appended raw -- attacking the cost at its
    // source rather than re-encoding it.
    //
    // NOT made a default, and deliberately NOT fitted to the locked suite: on
    // the three datasets checked it wins on H. salinarum (-1.33%) and
    // Drosophila (-2.38%) but LOSES on E. coli (+0.86%), and no validated
    // property separates those cases. Choosing a threshold that happens to
    // suit the 7 locked files would be exactly the overfit this project
    // forbids -- an eighth dataset could land the wrong side of it. So it stays
    // opt-in, evidenced, and off.
```

What the code now says:

```
    // Sweep seed width SW. 32 bases fills a uint64 exactly, which is why it is the
    // default, but it also floors the shortest overlap the sweep can see, and the
    // measured trend (MINOV 50 -> 40 -> 32 giving 29.88M -> 28.91M -> 28.28M) runs
    // into that floor. A narrower seed costs specificity (more candidates to verify)
    // but lets round 2 pair reads that overlap by less than 32.
    // LOWCOV is an explicitly requested low-coverage mode (SW = 8, the floor, below
    // which a seed is noise). It changes no default and is never auto-enabled: the
    // signal that would gate it (leftover_frac) is known only after round 1, the stage
    // it controls, and cheap coverage estimators were tried and rejected (a 25-mer
    // histogram gives mode 2 for every dataset). Measured, it gives -2.38% on
    // Drosophila SRR40104920 (2.6x) and -1.33% on H. salinarum but +0.86% on E. coli, and
    // no validated property separates those cases, so it stays opt-in instead of
    // being fitted to the locked set. The sweep and the rejected estimators are in
    // docs/extras/ENGINEERING_NOTES_ENCODER.md.
```

## Graph harvest: the idea, the measurement, and why it was refuted

Source: `src/encoder.cpp`, lines 1219-1256 before the documentation revision.

Original comment, verbatim:

```
    // ── GRAPH HARVEST -- REFUTED BY ITS OWN MEASUREMENT, kept instrumented ──
    //
    // THE IDEA: the overlap graph is already computed inside the sweep below
    // and thrown away. For each open tail `a` at level L, cand[i*CCAP +
    // 0..ccnt[i]-1] is a's OUT-EDGE SET -- every read whose first L bases are
    // verified equal (rcmp) to a's last L bases -- and the commit keeps the
    // first untaken entry and discards the rest. Recording the discarded
    // alternative looked free, and structurally matched DiscoSNP++'s own
    // seeding step (Bubble.cpp:309: successors.size()>=2, then every pair).
    //
    // THE MEASUREMENT (HG002 r2, 75,115 reads, 2026-09-04):
    //     commits=57,553   ccnt==1: 54,966 (95.5%)   ccnt==2: 2,453 (4.3%)
    //     surviving branches: 85     of which 73 (86%) at L = 100% of Lmax
    // Against 424 truth het-SNVs and 95 truth het-indels in the same window.
    //
    // WHY IT CANNOT WORK, and this is structural, not tuning:
    //   1. `rcmp` is an EXACT match. Heterozygosity means the two haplotypes
    //      DIFFER. Two haplotype reads can therefore never both exactly match
    //      the same suffix unless the variant lies outside the overlap window.
    //      A het site does not produce a tie here -- it produces exactly one
    //      surviving candidate, which is why ccnt==1 is 95.5% of commits.
    //   2. Of the ties that do occur, 86% sit at L = Lmax, and per the comment
    //      at the head of the level loop below, "At L = rlen[a] the offset is 0
    //      ... rcmp compares the two reads in full -- exactly the duplicate
    //      test." They are duplicate reads, not branches.
    //   3. 2,368 of the 2,453 two-candidate commits had the alternative
    //      already claimed (prv[other]!=NONE): two chains meeting, not a branch.
    //
    // So this harvests ~12 non-duplicate ties per 519 truth variants. It
    // cannot serve as a calling channel (accuracy), and it cannot replace
    // rc_reads/pkidx as a candidate source (RAM/time) because it finds ~0.15%
    // as many loci. Both halves of docs/_removable/GRAPH_HARVEST_EXACT_PLAN.md are dead.
    //
    // KEPT (not deleted) because it is CAPS_CALL-gated, verified bit-identical
    // on the Claim 1 path (ERR5181310 ARCHIVE_TOTAL=352997 and every stream
    // byte-for-byte), costs one predictable branch when off, and the printed
    // histogram is the evidence for the paragraph above. Full write-up:
    // docs/GRAPH_HARVEST_REFUTED.md.
```

What the code now says:

```
    // ── GRAPH HARVEST: instrumentation of a refuted idea, active only with CAPS_CALL ──
    // The idea: the sweep's overlap graph is computed and then discarded. For each open
    // tail `a` at level L, cand[i*CCAP + 0..ccnt[i]-1] is a's out-edge set (every read
    // whose first L bases equal a's last L bases under the exact comparison rcmp), and the
    // commit keeps the first untaken entry. Recording the discarded alternative looked free
    // and matched DiscoSNP++'s own seeding step (Bubble.cpp:309).
    // Measured (HG002 r2, 75,115 reads): of 57,553 commits 95.5% had exactly one
    // candidate, and 85 branches survived, 73 of them (86%) at L = Lmax, which are duplicate
    // reads, against 424 truth het-SNVs and 95 het-indels in the window. It cannot work,
    // structurally: rcmp is an exact match, so two haplotype reads cannot both match one
    // suffix unless the variant lies outside the overlap, and a het site yields one candidate.
    // It finds about 0.15% as many loci as rc_reads/pkidx, so it can serve neither as a
    // calling channel nor as a candidate source.
    // Kept, not deleted, because it is gated by CAPS_CALL, was verified bit-identical on
    // the archive path (ERR5181310 ARCHIVE_TOTAL=352997 and every stream byte for byte), and
    // costs one predictable branch when off. The printed histogram is the evidence above.
```

## Mapping stage: seed geometry, acceptance, and why MAXMAP is Lmax/13

Source: `src/encoder.cpp`, lines 2431-2516 before the documentation revision.

Original comment, verbatim:

```
        // Part width must SCALE with the tolerance, not stay pinned at 32.
        // Pigeonhole needs MAXMM+1 parts to fit inside one read, so a fixed
        // 32-base part caps tolerance at 150/32-1 = 3 mismatches -- ask for more
        // and the parts no longer fit, no seed is indexed, and mapping silently
        // does nothing. PgRC2 makes the part length the parameter and derives
        // the tolerance from it (targetMismatches = readLength/partLength - 1).
        // Same here: parts are as wide as they can be, capped at 32 because that
        // is what a uint64 holds exactly at 2 bits per base.
        // STAGE 17. The above ties acceptance to the seed geometry: NPARTS =
        // MAXMM+1 parts means a read is rejected past 3 mismatches, because that
        // is all the pigeonhole lemma can certify with parts this wide.
        //
        // PgRC2 does not tie them together. ReadsMatchers.cpp:700 reads
        //     uint8_t maxMismatches = readLength / minCharsPerMismatch;
        // with minCharsPerMismatch = 3 at CODER_LEVEL_NORMAL, so a 150-base read
        // is accepted at up to FIFTY mismatches. The seed (their
        // readsExactMatchingChars = 38) only proposes candidates; acceptance is a
        // separate, far looser test. Their log shows what that buys:
        //     290,521 LQ reads -> 11,298 survivors   (96.1% placed)
        // against ours at 74.8%. The survivor pseudogenome that holds 94% of the
        // remaining gap is made almost entirely of reads this limit rejected.
        //
        // Their seed is 38 bases; ours stays 32 because that is what a uint64
        // holds exactly at 2 bits per base. A shorter seed is more sensitive, not
        // less, so this is not a handicap -- it proposes more candidates, and
        // acceptance is what decides.
        //
        // Not copied: they re-scan and keep the BEST match per read (the
        // "better-matches" column in their log, ReadsMatchers.cpp:315-328). That
        // only shrinks their mismatch stream, which this progression does not
        // measure, so first-acceptable wins here. Noted rather than silently
        // skipped.
        // STAGE 20. Seeds were indexed at fixed multiples of SEEDW (0,32,64,96),
        // which has two defects at a 150-base read: a read whose only clean
        // 32-mer starts at, say, offset 17 is invisible, and bases 128-149 are
        // in no seed at all. PgRC2 does not do this -- at CODER_LEVEL_NORMAL the
        // matcher is mode 'c', CopMEMReadsApproxMatcher, which finds maximal
        // exact matches at ANY offset (mapReadsIntoPg, ReadsMatchers.cpp:715-735).
        // Their fixed 38 is a minimum MEM length, not a position.
        //
        // Note their pigeonhole guarantee is only
        //     targetMismatches = readLength/exactMatchingChars - 1 = 150/38-1 = 2
        // while maxMismatches is 50, so past 2 mismatches their seed is a
        // heuristic filter too. The difference is purely that their filter looks
        // everywhere and ours looked in four places.
        //
        // Sliding by SEEDSTRIDE recovers most of that for a proportional index
        // cost and no extra scan cost -- the pg is still swept once, one hash
        // probe per position, only the candidate lists get longer.
        // MAXMAP: mismatch ceiling for accepting a read placement.
        // Lmax/3 copies ReadsMatchers.cpp:700, but that ceiling was chosen when
        // the mismatch POSITION and COUNT streams were not being counted in our
        // total (see the Stage 27 note below -- it predicted exactly this).
        // With those streams counted, the real economics are: a mismatch costs
        // ~0.88 B (position+count+symbol, measured) while a base stored as pg
        // literal costs ~0.229 B (measured), so for a length-L read the
        // break-even is m ~= 0.26*L -- BELOW Lmax/3. Accepting past that point
        // buys a shorter pg with a mismatch stream that costs more than the
        // bases it saved. Env-overridable so the real optimum can be swept
        // against the CORRECTED total rather than assumed.
        // Re-swept on REAL full files after the mm_pos/mm_cnt transforms made
        // mismatches cheaper, which shifts this optimum upward (cheaper
        // mismatches -> accept more of them -> shorter pg). Real optima now:
        // E. coli 20, P. aeruginosa 27, but validated across all five real
        // files a fixed 20 LOSES 1.2% on P. falciparum, making the aggregate a
        // wash (+0.04% worse than a fixed 12). Minimax regret picks 12: worst
        // case 0.54% versus 20's 1.24%. So 12 stays, now justified on real
        // full files rather than the subsamples it was first derived from.
        // The property was identified: READ LENGTH. Every dataset the constant
        // was tuned on is 100-151 bp, a range over which a constant and a ratio
        // are indistinguishable. PgRC2 has always used the ratio form --
        // ReadsMatchers.cpp:700, `maxMismatches = readLength / minCharsPerMismatch`
        // -- i.e. one mismatch permitted per N bases, not a fixed count.
        //
        // A fixed ceiling is wrong in both directions: at 251 bp it is 4.8%
        // tolerance where 151 bp gets 8%, so long reads are under-mapped and
        // spill into the pg as literal; at short read lengths it would be too
        // permissive. Expressed as a ratio it is scale-free.
        //
        // The divisor comes from our own measured economics, not from theirs:
        // a mismatch costs ~0.88 B (position+count+symbol) against ~0.229 B for
        // a pg-literal base, so the break-even is ~0.26*L. The tuned optimum of
        // 12 at L=151 is L/12.6; L/13 reproduces it (11.6 at 151, 7.7 at 100,
        // 19.3 at 251) while scaling correctly, and stays well under the
        // 0.26*L cost ceiling. MIN_MAXMAP floors it so very short reads still
        // get a usable tolerance.
```

What the code now says:

```
        // MAPPING SEEDS AND TOLERANCE.
        // Pigeonhole needs MAXMAP+1 parts to fit inside a read, so a part width fixed at
        // 32 would cap the tolerance at 150/32-1 = 3 mismatches, and asking for more would
        // index no seed at all. So the seed width is a parameter, at most 32 (what a uint64
        // holds at 2 bits per base), and acceptance is decoupled from it. A seed only
        // PROPOSES a candidate placement. Acceptance is a separate and looser test, as in
        // PgRC2, where the seed is a heuristic filter beyond what the lemma certifies
        // (ReadsMatchers.cpp:700 accepts up to readLength/3 mismatches).
        // Seeds are taken from each read every SEEDSTRIDE bases, not only at fixed multiples of
        // SEEDW, so a read whose only clean stretch falls between them is still found, and the
        // pseudogenome is swept once with one hash probe per position. A read is found if it
        // shares an exact stretch of at least SEEDW + SEEDSTRIDE*K2 - 1 bases with the pseudogenome.
        // MAXMAP, the mismatch ceiling for accepting a placement, is Lmax / 13, at least 6.
        // It is a ratio, not a constant, because a read of L bases is worth mapping only up to
        // the break-even where its stored mismatches cost as much as the bases they save. Measured:
        // a mismatch costs about 0.88 B (position, count and symbol) against about 0.229 B for a
        // pseudogenome literal base, so the break-even is about 0.26 * L, and accepting past it buys
        // a shorter pseudogenome with a mismatch stream that costs more than it saved. The best fixed
        // value differs by dataset (E. coli 20, P. aeruginosa 27, and a fixed 20 loses 1.2% on
        // P. falciparum). Minimax regret over the five real files chose 12, and the divisor 13 reproduces
        // that at 151 bp (11.6) while scaling with read length (7.7 at 100 bp, 19.3 at 251 bp). MAXMAP_DIV
        // overrides the divisor, for sweeps. The derivation and the PgRC2 comparison are in
        // docs/extras/ENGINEERING_NOTES_ENCODER.md.
```

## Order permutation analysis and the switch to direct positions

Source: `src/encoder.cpp`, lines 3441-3471 before the documentation revision.

Original comment, verbatim:

```
    // ── STAGE 22: the order permutation (their stage 6) ─────────────────────
    // To restore the original file order the archive must carry, for every
    // original read, where it sits in emission order. PgRC2 stores the inverse
    // permutation as a raw uint32 array and LZMAs it
    // (SeparatedPseudoGenomePersistence.cpp:226-233) -- no permutation-specific
    // coding at all in single-end mode.
    //
    // Measured, that costs them 7,063,459 - 4,898,620 = 2,164,839 B on this
    // input. The floor for a RANDOM permutation of 999,340 elements is
    // log2(n!) = 2,309,000 B -- so their output is 144 KB BELOW the random
    // floor, which means the permutation is not random and LZMA is finding real
    // structure in it. Any coder that assumes uniformity (Lehmer / factorial
    // base) would therefore do WORSE than what they already have, not better.
    // That is why this stage dumps the permutation for analysis instead of
    // going straight to a clever coder.
    //
    // Note positions from the RC pass are in reverse-complement coordinates.
    // That is fine here: only the ORDER matters for this measurement, and a
    // real encoder would store the strand flag alongside (ARCS already does).
    // STAGE 100 REWRITE -- scrapped the whole rank/bucket/uidOrder scheme
    // after reading PgRC2's real source line by line
    // (SeparatedPseudoGenomePersistence.cpp:446, compressReadsPgPositions,
    // singleFileMode branch -- our exact case, single-end). Their design for
    // SE mode is a single flat array: position[original_read_index] =
    // absolute pg position, written directly in original-read order. No
    // rank layer, no duplicate-count bucket, no separate order/permutation
    // stream at all -- position alone restores both WHERE a read's sequence
    // is AND, by construction (reads emitted in original index order),
    // its place in the file. This replaces layers 2 (perm.u32), 4, 7, and
    // 8 (orig2uid_ranks.bin/uidorder.bin) all at once with one simpler,
    // directly-indexed design that has no indirection left to get wrong.
```

What the code now says:

```
    // ── Read order and placement ───────────────────────────────────────────
    // To restore the original file order the archive must say, for every original read,
    // where its sequence is. PgRC2's design for single-end input is one flat array,
    // position[original read] = absolute pg position, written in original order
    // (SeparatedPseudoGenomePersistence.cpp:446), so position alone restores both a
    // read's sequence and its place in the file, with no separate permutation. Here the
    // positions are stored once per UNIQUE read, and orig2uid maps each original read to its
    // unique read. An earlier scheme with a stored order permutation was scrapped: a random
    // permutation of n reads costs about log2(n!) (2,309,000 B for 999,340 reads), and
    // PgRC2's own permutation coded 144 KB BELOW that floor, which shows the order is not
    // random and that a coder assuming uniformity would do worse. DUMP_PERM gates this
    // block and is required for a usable archive.
```

## Position storage: per unique read, and fixed width instead of varints

Source: `src/encoder.cpp`, lines 3474-3514 before the documentation revision.

Original comment, verbatim:

```
        // STAGE 100 SIZE FIX -- position/strand/length were being stored once
        // PER ORIGINAL READ (999,308 for E. coli), redundantly repeating the
        // exact same value for every duplicate of a unique sequence (150,199
        // of them). `orig2uid.bin` (original -> unique-id) is needed anyway
        // for mismatch correlation and already gives decode everything it
        // needs to look up ANY original read's data -- so position/strand/
        // length only need to exist ONCE PER UNIQUE READ (849,109 here, real
        // ~15% fewer entries, before even counting that xz now sees no
        // redundant repeated values to (imperfectly) find on its own).
        // Confirmed this doesn't reintroduce the old rank/bucket bug: unlike
        // that scheme, orig2uid.bin stores the FINAL ANSWER directly (which
        // unique-id), no further indirection or missing bucket-size info
        // needed -- this is exactly the design already proven byte-identical
        // on E. coli, just no longer wastefully duplicated per original read.
        // STAGE 100 SIZE FIX ROUND 2 -- TRIED AND REVERTED. Sorting by
        // position + delta-coding did shrink pos_abs.bin back to the old
        // scheme's size (526,888 B, confirmed) -- but required a new
        // rank_of_uid[] permutation array to stay invertible, and that
        // permutation cost 2,275,236 B (close to log2(n!) raw entropy,
        // since assembly order has little correlation with position-sort
        // order). Net effect on real E. coli data: 5,784,467 -> 5,790,447,
        // essentially flat (+0.1%) -- the cost was relocated, not removed.
        // This is the same fundamental cost the old broken perm.u32 scheme
        // was paying, just moved to a different array. Reverted to the
        // simpler direct-indexing version below (no permutation array
        // needed at all), which gives the same real total with less
        // complexity and no permutation-cost risk on other datasets.
        // STAGE 100 SIZE FIX ROUND 3 -- real, established, "not flashy"
        // technique (confirmed this session: DEFLATE itself picks per-block
        // between stored/static-Huffman/dynamic-Huffman, keeping whichever
        // is smaller -- the exact same "try a plain encoding, keep it if
        // smaller" philosophy, decades-proven, not invented here). Tested
        // directly on real data: fixed-width uint32 positions, xz'd,
        // consistently beat the varint encoding this stream used before --
        // E. coli 2,549,216 vs 2,766,428 (-7.9%), P. aeruginosa 232,668 vs
        // 261,680 (-11.1%). Varint's per-value length-prefix bits apparently
        // break up the byte-alignment patterns xz's LZ77 stage would
        // otherwise find across positions -- fixed-width keeps that
        // structure intact. Switched from variable-length varint to a
        // fixed uint32 array; decode_locked_seqorder.py updated to match
        // (read_u32 instead of read_varints for this stream).
```

What the code now says:

```
        // Position, strand and length are stored once per UNIQUE read, not once per original
        // read: orig2uid already maps any original read to its unique read, so repeating the
        // values for every duplicate was pure redundancy (150,199 of 999,308 entries on E. coli).
        // Two layouts were measured. Sorting by position and delta-coding shrank pos_abs back to
        // the old scheme's size (526,888 B) but needs a rank permutation costing 2,275,236 B,
        // close to log2(n!), so the total was flat (+0.1%) and it was reverted. Fixed-width uint32
        // positions under xz beat varints by 7.9% on E. coli (2,549,216 vs 2,766,428) and 11.1% on
        // P. aeruginosa, because a varint's length bits break the byte alignment that xz's LZ77
        // finds across positions, so the stream is a fixed uint32 array.
```

## MEM match merging, the memory fix, and the extension-tolerance rule

Source: `src/encoder.cpp`, lines 4094-4150 before the documentation revision.

Original comment, verbatim:

```
        // STAGE 100 REWRITE -- the two separate per-region c[]/cr[] coverage
        // bitmaps below (used only for nm_main/nm_second reporting now) were
        // ALSO independently driving literal.txt's writing -- but a decoder
        // trying to reconstruct literal.txt's coverage from the separately-
        // dumped `allrefs` triples could NOT reliably reproduce the same
        // positions (confirmed directly: byte COUNTS matched by coincidence,
        // but actual covered POSITIONS did not -- every reconstructed read
        // came out wrong, 0/15424 on a real decode attempt, even for reads
        // with zero mismatches). Root cause: `allrefs` contains redundant,
        // overlapping candidate matches (this session found 17,343 exact
        // duplicate destinations on real E. coli data) with no way to
        // recover, from the dump alone, which ones the ORIGINAL c[]/cr[]
        // merge actually kept.
        //
        // Fix, matching PgRC2's real, proven design exactly (SimplePgMatcher
        // .cpp:85-160, read line by line this session): sort ALL matches
        // (both regions, forward+RC, now unified in one `allrefs` with a
        // real is_rc flag -- previously missing entirely) by destination,
        // then walk them TRIMMING overlap with already-consumed territory
        // instead of just OR-ing a bitmap. This guarantees, by construction,
        // that literal.txt's mark positions and the emitted (dst,src,len,
        // is_rc) triples describe the EXACT same non-overlapping partition
        // of `pg` -- nothing to reconstruct or get out of sync later.
        // RAM FIX (round 2, C. elegans scale) -- found by bracketing the
        // MEM-matching phase with fine-grained RSS checkpoints: internal
        // checkpoints never exceeded 926MB, yet the kernel-tracked peak
        // (/usr/bin/time -v) was consistently 1.27-1.30GB -- a transient
        // spike DURING this block, invisible to before/after sampling.
        // Concrete, measured overlap found: `c`/`cr` (main-pg coverage
        // bitmaps, main_pg_end bytes each -- 135MB combined at C. elegans
        // scale) were never freed after their last real use (the rem_main
        // tally just below) -- they stayed allocated for the ENTIRE rest
        // of the function, fully overlapping the second pass's own `c2`/
        // `cr2`/`Q`/`R` allocations (qlen bytes/chars each -- another
        // ~360MB at this scale). Fixed by scoping `c`/`cr` so they free
        // before the second pass begins, same principle as the earlier
        // allrefs/cleanRefs fix: don't hold data alive past its last use.
        // Extension mismatch tolerance, derived from a measured property, not a
        // per-dataset switch: MEM_MAXMM=0 whenever leftover_frac is small (every
        // one of the 7 locked datasets), rising toward REF_MAXMM as more of the
        // file fails to chain at all in round 1 -- exactly the population that
        // needs a candidate overlap to survive one sequencing error to be
        // captured at all. Applied only to the main self-match and the
        // second-region CROSS match; the SECOND_SELF pass is left at 0
        // (unaffected) since it was already measured to lose economically at
        // MAXMM=0 and more matches there would only make that worse.
        // Measured on the 7 locked datasets: leftover_frac ranges 0.222-0.518.
        // It is NOT near zero at normal coverage as first assumed -- duplicate
        // reads and reads whose only overlap falls at a tail edge also end up
        // leftover, not just true low-coverage gaps. So a bare proportional
        // formula fired everywhere and changed every locked archive. Rebuilt
        // with a floor above the measured locked maximum (0.518) plus real
        // margin, and a ceiling near Drosophila's 2.6x-coverage value (0.812):
        // T0=0.60 keeps MEM_MAXMM=0 on every locked dataset with room to
        // spare; T1=0.85 is where the tolerance reaches its cap. This is still
        // a function of the measured property, not a per-dataset switch --
        // any future dataset lands wherever ITS leftover_frac places it.
```

What the code now says:

```
        // Matches are sorted by destination (both regions, forward and reverse complement, in
        // one list with a real is_rc flag) and walked while TRIMMING overlap with territory
        // already consumed, as PgRC2's SimplePgMatcher.cpp:85-160 does. That guarantees, by
        // construction, that literal.txt's mark positions and the emitted (dst, src, len, is_rc)
        // references describe the same non-overlapping partition of `pg`. An earlier version merged
        // coverage bitmaps instead and could not be reconstructed by a decoder: the redundant
        // overlapping candidates (17,343 exact duplicate destinations on E. coli) left no way to
        // tell which ones the bitmap merge had kept.
        // Memory: the main-region coverage bitmaps c and cr (about 135 MB at C. elegans scale) are
        // scoped so they are freed before the second pass allocates its own, because the
        // kernel-tracked peak (1.27-1.30 GB) was above every internal checkpoint (926 MB).
        // Extension mismatch tolerance MEM_MAXMM is derived from leftover_frac, a measured property,
        // not set per dataset. It is 0 up to T0 = 0.60 and reaches its cap REF_MAXMM at T1 = 0.85
        // (Drosophila at 2.6x is 0.812). The locked datasets span 0.222-0.518, so they all stay at 0.
        // A bare proportional formula fired everywhere, because duplicate reads and reads whose only
        // overlap falls at a tail edge also end up leftover, not just true low-coverage gaps. It
        // applies to the main self-match and the second-region CROSS match only. The SECOND_SELF
        // pass stays at 0: it already loses economically at 0 and more matches would worsen that.
```
