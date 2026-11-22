# Engineering notes: stream coders

Long explanatory comments from `include/coders_inproc.h`, `seqpar_core.h`, `names_coder.h` and `quality_coder.h` that were condensed in the source during the documentation revision. They are kept here in their original wording, with references to individual working sessions removed, because they record measurements and the reasoning behind design choices. The source keeps a short statement of each conclusion.


The quoted text is unedited, so it can name working documents (for example under `docs/_removable/`, or a `docs/*.md` that is not in this repository) and sibling code (ARCS, PgRC2) that are not shipped here. Those names are provenance only. Nothing in the code or the paper depends on them.

## MEM source coder: per-component profile against PgRC2 and the derivation of the bound

Source: `include/coders_inproc.h`, lines 293-333 before the documentation revision.

Original comment:

```
// Reference-stream coder: sources bounded by their own destination.
//
// Profiled per component at MINMEM 24, RC-only self-matching:
//
//   component     ours    bits/match      PgRC2   bits/match    delta
//   sources    186,852         25.38    128,671        23.83  +58,181
//   lengths     63,988          8.69     48,509         8.99  +15,479
//   dest gaps   43,932          5.97          0            -  +43,932
//
// Two things that profile settles.
//
// The gap stream is NOT a loss. Their destination positions are in-band: a '%'
// MATCH_MARK sits in the pg literal and rides along inside the VarLenDNA
// codebook, which carries phrases like "T%", "AT%", "GG%". So the fair unit is
// literal plus positions -- ours 2,985,570 + 43,932 = 3,029,502 against their
// 3,056,474, and we are 26,972 AHEAD. An earlier note in this file called the
// gap stream pure overhead they avoid; that was double-counting.
//
// Sources ARE the loss, and the reason is a varint. log2(23,233,953) = 24.47
// bits is the naive floor and we spend 25.38, because a varint for values near
// 23M needs four bytes. Theirs sits at 23.83 -- BELOW their own naive floor --
// which means LZMA is finding structure in it.
//
// The structure is not subtle once looked for: in a self-match the source is
// always earlier than the destination, and matches are stored in destination
// order. So src is uniform in [0, dst), not in [0, pg_len), and coding it that
// way costs log2(dst), which averages about 1.4 bits below log2(pg_len) for
// destinations spread across the pg. That is ~23.0 bits/match, under what they
// manage, and it needs no new information -- the decoder already knows dst when
// it reads src.
//
// Cross-matches (survivor pg into main pg) are the one exception: there the
// source lives in the main pg while the destination is past main_pg_end, so the
// bound is main_pg_end rather than dst. The decoder can tell the two apart from
// dst alone, so nothing extra is stored to distinguish them.
//
// Range coder is the same LZMA-style one used for the permutation in stage 23.
// Round trip is verified before any size is reported.
//
//   g++ -O3 -march=native -o refc 37_ref_coder.cpp
//   ./refc mem_triples.bin <pg_len> <main_pg_end>
```

What the code now says:

```
// Coder for the SOURCE position of each MEM reference (the mem_triples stream).
// A reference is (dst, src, len, rc). dst, len and rc travel in their own streams
// (mem_dstgap, mem_len, mem_rc), so this stream holds only src. References are
// coded in destination order, each src uniform within a bound derived from its
// own dst, which the decoder already knows when it reads src:
//   self-match inside the main pg:           src is earlier than dst, so src in [0, dst)
//   cross-match (source in the main pg, destination past main_pg_end):
//                                            src in [0, main_pg_end)
//   second-region self-match:                src coded relative to main_pg_end,
//                                            flagged in the mem_self stream
// Coding src in [0, dst) instead of [0, pg_len) costs log2(dst) bits, on average
// about 1.4 bits fewer per match than log2(pg_len) when destinations are spread
// across the pg. The measured per-component profile against PgRC2 (sources,
// lengths, destination gaps, and why the gap stream is not a loss) is in
// docs/extras/ENGINEERING_NOTES_CODERS.md.
```

## mm_pos buckets coded in parallel: why and how much

Source: `include/coders_inproc.h`, lines 818-832 before the documentation revision.

Original comment:

```
    // ── BUCKETS CODED IN PARALLEL ───────────────────────────────────────────
    //
    // MEASURED: this function is 1.58 s of a 2.13 s mm_pos job, and mm_pos is
    // the coding pool's LONGEST job -- the pool runs at its Amdahl floor, so
    // this loop alone sets the encoder's floor. It is a nested search: per
    // bucket it builds two layouts and probes FOUR coders on each, all serial.
    //
    // Buckets are independent -- each reads only its own vector and produces
    // its own record -- so they are coded concurrently and the RESULTS ARE
    // EMITTED IN THE ORIGINAL std::map ORDER afterwards. Byte-identical output:
    // same candidates, same smaller-wins rule, same emission order.
    //
    // Concurrency is bounded because each in-flight bucket holds an LZMA
    // encoder state; buckets are few (one per distinct mismatch count) and the
    // pool already has this thread, so the bound is deliberately modest.
```

What the code now says:

```
    // Buckets are independent, so they are coded concurrently (this loop was the
    // longest job in the encoder's coding pool) and emitted afterwards in the
    // original std::map order, which keeps the output byte-identical to serial
    // coding. Concurrency is deliberately modest because each in-flight bucket
    // holds an LZMA encoder state. Measurements: docs/extras/ENGINEERING_NOTES_CODERS.md.
```

## seqpar_core.h: header history (the stage 75 and stage 81 notes)

Source: `include/seqpar_core.h`, lines 2-72 before the documentation revision.

Original comment:

```
// STAGE 81 -- encode-only mode on top of stage 75, unchanged otherwise.
//
// User's bigger-scope question after stage 80: does the earlier "beat
// SPRING on speed" combined-pipeline comparison fairly compare compress
// time to compress time? Diagnostic done FIRST across every remaining
// coder (per user direction: think and understand before executing),
// measuring encode-only vs full (encode+decode+self-verify) for each:
//
//   tool        encode-only   full     verify overhead
//   seqpar      0.313s        0.729s   0.416s (57%)
//   namebq      0.203s        1.238s   1.035s (84%)
//   qualbq      0.219s        1.622s   1.403s (87%)
//   permcoder   0.028s        0.176s   0.148s (84%)
//   refcoder    ~0            0.004s   negligible
//
// The percentages look uniformly dramatic, but only ONE of these sits on
// the real critical path. names/quality (namebq/qualbq) already run
// CONCURRENT with assembly (stage 76's overlap, ~3.4-3.6s) -- their full
// 1.2-1.6s runtime finishes comfortably inside that window regardless, so
// their verify overhead is already fully hidden and cutting it would not
// change the pipeline's total wall clock at all. permcoder's full time
// (0.176s) is dwarfed by seqpar's (0.729s) in the post-assembly phase, so
// fixing permcoder alone would not move the total either -- seqpar would
// still be the bottleneck. seqpar is the one tool whose 0.416s of
// self-verification sits on the critical path (it runs AFTER
// assembly, not absorbed by the overlap), so it is the one this stage
// fixes -- not all five, because the diagnostic showed the other four
// would not move the real number.
//
// None of these tools ever wrote a compressed archive to disk (they are
// measurement instruments throughout this project, confirmed by reading
// their own main()), so the earlier combined-pipeline comparison against
// real SPRING/Genozip was, on this specific point, charging itself for
// internal self-verification work no real archiver performs during normal
// compression -- SPRING/Genozip don't decode what they just encoded
// before writing it out.
//
// ENCODE_ONLY=1 env var skips the round-trip decode+compare block;
// default (unset) behavior is UNCHANGED from stage 75, verified: same
// coded byte count (3,007,051 B) in both modes, encode output is
// completely unaffected since verify runs strictly after and never
// touches the encoded bytes.
//
// Real, measured, repeated result (3 runs each, real data):
//   default (with verify): 0.743s / 0.749s / 0.750s
//   ENCODE_ONLY=1:          0.399s / 0.393s / 0.384s
// Roughly halves the tool's own time, matching the diagnostic's 0.416s
// prediction almost exactly.
//
//   g++ -O3 -march=native -pthread -o seqpar81 81_seqpar_encode_only.cpp
// ============================================================================
//
// STAGE 75 -- threaded, block-chunked sequence/literal coder. Real,
// disclosed extension of the same proven pattern already used for names
// (stage 69) and quality (stage 71): split into fixed blocks, fresh model
// per block, dispatch across threads. Applied here to stage 74's exact
// context-mixing model (constexpr orders, fused gather+mixer loop,
// TBITS=16) -- no algorithm change, only parallel framing, same discipline
// as every other threaded stage.
//
// Real, KNOWN, disclosed risk before testing: this model's match component
// specifically exploits LONG-RANGE repeats (up to the full literal
// stream's length via a 2^24-slot hash table) -- unlike names/quality
// where cross-block correlation was already shown to not exist across
// blocks (X coordinate has no cross-read structure; quality reflects
// per-read sequencing noise, not cross-read pattern), sequence DOES have
// real, exploited long-range structure (repeats), so chunking is expected
// to cost MORE here than it did for names/quality -- tested for real
// below, not assumed away.
//
//   g++ -O3 -march=native -pthread -o seqpar 75_dna_mix_par.cpp
```

What the code now says:

```
// DNA context-mixing coder. Codes a sequence of 2-bit base codes (A=0, C=1,
// G=2, T=3) and decodes it again. The encoder uses seq_encode_mem and the
// decoder seq_decode_mem, both defined here so the model cannot drift between
// the two. It codes the `literal` stream: the pseudogenome bases that no MEM
// reference covers.
//
// Model, rebuilt from scratch for each chunk (see run_chunk):
//   * Each base is two binary decisions: the high bit (node 0), then the low bit
//     given the high bit (node 1 or 2). Probabilities are 12-bit.
//   * Ten context models, one per order in ORD (the previous 1, 2, 3, 4, 6, 8,
//     11, 14, 18 and 22 bases). Each hashes its context and the node into a table
//     of 2^TBITS 16-bit probabilities. The learning rate starts at 1/2 and settles
//     at 1/32 as a slot's count grows. There are no collision checks.
//   * A match model: the last MINLEN = 24 bases are hashed to the position that
//     followed their previous occurrence, and the base found there is the
//     prediction. Its confidence is learned per match length, capped at 31. It
//     reaches back only to the start of its own chunk.
//   * A logistic mixer combines the eleven stretched predictions. The weight set
//     is chosen by the node and the previous two bases.
//   * A final APM (secondary estimation) stage refines the mixed probability,
//     using the node and the previous four bases. The coded probability is one
//     part mixer output to three parts APM output.
//   * A 32-bit carry-less binary arithmetic coder (BinEnc, BinDec) codes each
//     decision.
//
// Stream layout from seq_encode_mem: [nchunks:u32][total_symbols:u64], then per
// chunk [symbols:u32][coded_bytes:u32], then the chunk payloads back to back.
// Chunks are independent, so they are coded and decoded in parallel. The chunk
// count (4 by default, env SEQT in the encoder) is a COMPRESSION parameter, not
// only a threading one: this model's match component exploits long-range repeats,
// so more chunks mean less context per chunk. On the E. coli literal, 12 chunks
// gave 1,809,073 B, 4 chunks 1,804,223 B and 1 chunk 1,801,710 B.
//
// The development history of this coder, including the measurement that decided
// where its verification overhead sits, is in docs/extras/ENGINEERING_NOTES_CODERS.md.
```

## seqpar_core.h: why the match-model table is sized from the chunk

Source: `include/seqpar_core.h`, lines 130-148 before the documentation revision.

Original comment:

```
// STAGE 101 -- real root cause found via direct code inspection + layer-
// by-layer profiling (LAYER_BY_LAYER_ANALYSIS.md): mmtab was a FIXED
// 2^24-slot (16.7M entries x 4B = 64MB) hash table, allocated fresh per
// CHUNK regardless of that chunk's actual data size. With 12 parallel
// chunks (this pipeline's default), that is ~768MB, confirmed to match
// the pipeline's observed ~820MB peak RAM almost exactly -- and confirmed
// FIXED, not input-scaled, by measuring it unchanged (819,612 KB vs
// 821,536 KB) across a 5MB and a 7MB literal.txt.
//
// Verified safe to resize before touching anything: mmtab is a pure
// match-finder heuristic (used only to PREDICT a repeat for the context
// model, lines ~224-227 below) -- on a miss or collision it just falls
// back to mlen=0, no correctness dependency. Confirmed against a real,
// established precedent, not invented: zstd's own ZSTD_adjustCParams
// scales hashLog/chainLog down from a srcSizeHint for exactly this
// reason (small input doesn't need a hash table sized for a huge one).
// Same principle as the earlier alphabet-adaptive quality-
// coder fix -- a formula based on a MEASURED property of the actual
// chunk, not a fixed constant or a per-dataset guess.
```

What the code now says:

```
// mmtab, the match model's hash table, has 2^b slots with b = ceil(log2(chunk
// symbols)), at least 16 (64K slots, 256 KB, so a tiny chunk does not get a
// degenerate table) and at most MMBITS_MAX. It used to be a fixed 2^24 slots (64 MB)
// per chunk whatever the chunk's size. With 12 parallel chunks that was about
// 768 MB, which matched the pipeline's observed ~820 MB peak. The table is only a
// match-finder heuristic: a miss or collision gives no prediction (mlen = 0), so
// sizing it to the chunk has no effect on correctness. zstd scales its hash tables
// from a source-size hint for the same reason (ZSTD_adjustCParams).
```

## names_coder.h: header history (what the stage 86 coder lacked)

Source: `include/names_coder.h`, lines 2-41 before the documentation revision.

Original comment:

```
// Names / read-ID coder, shared by the encoder (106_inprocess.cpp) and the
// archive decoder (capsule_decode.cpp) so the two cannot diverge -- the same
// reason seqpar_core.h exists.
//
// Algorithm: stage 86 (= stage 70's global frozen dictionary + stage 85's
// bounded-memory streaming queue, merged). Namespaced `nmc` following this
// project's existing per-coder convention (mmc, refc, pgc), each of which
// carries its own RangeEnc/RangeDec.
//
// FOUR THINGS THIS FIXES relative to stage 86, all of which a mechanical
// copy of that file would have carried into the archive:
//
// 1. Stage 86 streams from a names.txt that someone extracted by hand. Here
//    pass 1 and pass 2 both stream the ID column straight out of the FASTQ,
//    so nothing accumulates. Holding the headers instead (the obvious
//    MemStream route, matching every other stream in this pipeline) would be
//    ~5 GB on T. cacao and would defeat the entire point of merging stage 85
//    in -- the bounded queue would sit behind a fully-materialized input.
//
// 2. Stage 86's "dict_header~=25356 B" is `entries*4`, arithmetic, never
//    serialized. Here the dictionary is really serialized, and returned RAW
//    for the caller to run through best_encode() like every other stream --
//    so its cost is real bytes in the archive, and it gets the coder
//    selector instead of raw uint32s.
//
// 3. Stage 86 has NO block index: its round trip re-reads the original file
//    to learn how many names are in each block. A decoder holding only the
//    archive cannot do that. Block byte-lengths and name-counts are a real
//    stream here and are really counted. This is the exact bug class that
//    has already cost this project twice (mm_pos/mm_count omitted from
//    totals; mem_triples' destination/length/rc dropped as "diagnostics
//    only" leaving an archive that could not decode) -- a cost landing in a
//    stream nobody measures.
//
// 4. Decode streams block-by-block to a FILE* rather than returning a
//    vector<string> of every name, so the decoder keeps the same bounded
//    memory property as the encoder.
//
// The '@' is not stored: it is structurally constant on every FASTQ header
// line and is re-added on output. SPRING and Genozip both do the same.
```

What the code now says:

```
// Names (read-ID) coder, shared by the encoder (src/encoder.cpp) and the archive
// decoder (src/decoder.cpp) so the two cannot diverge, for the same reason
// seqpar_core.h is shared. Namespace `nmc`.
//
// Approach (SPRING's): split each FASTQ header into tokens and code each token
// against the token at the SAME position in the previous header. A token is a run
// of letters, a run of digits, a run of zeros, or one other character. The '@' is
// not stored, since it is constant on every header line, and is re-added on output.
// Token types (TokType), each with its own adaptive model per token position:
//   ID_MATCH   identical to the previous header's token
//   ID_DELTA   a number 1..255 above the previous one
//   ID_ZDELTA  a signed difference within +-32,768, offered only for a field where,
//              after at least 20 observations, 30% or more of the earlier
//              differences were within +-2,048
//   ID_DIGIT   any other number: looked up in that position's value dictionary,
//              or coded as four bytes if the value is not in it
//   ID_ZEROS, ID_ALPHA, ID_CHAR   a run of zeros, a run of letters, one character
//   ID_SEQLEN  the number is this read's own length, which the archive already
//              holds in read_lengths, so nothing is stored
//   ID_END     end of header
// A digit run is split once its value reaches 2^28, so every number fits 32 bits.
//
// Two streaming passes over the FASTQ ID column, so the headers are never all
// resident:
//   pass 1 (dict_pass) measures, for each token position, whether a value
//     dictionary is cheaper than coding the number's bytes (priced from the
//     pass's own histograms), which positions hold identical text in every read,
//     which always equal the read length, and whether line 3 of each record is
//     '+' alone or '+' followed by the header.
//   pass 2 codes blocks of BLOCK names (250,000 by default) on a thread pool, each
//     block with fresh adaptive models, behind a bounded queue.
//
// Result (struct Encoded), which becomes three archive streams:
//   body   names_body   the blocks' range-coded payloads, concatenated
//   dict   names_dict   the token layout, the constant tokens and the value
//                       dictionaries, serialized by serialize_dict
//   index  names_index  per block, the payload byte length and the name count, as
//                       varints, so a decoder holding only the archive can find
//                       every block
// dict and index are returned RAW so the caller codes them through best_encode
// like every other stream. decode_to_file streams one block at a time, so the
// decoder has the encoder's bounded-memory property.
```

## names_coder.h: the two defects fixed in the dictionary frequency model

Source: `include/names_coder.h`, lines 143-166 before the documentation revision.

Original comment:

```
// TWO REAL DEFECTS FIXED HERE, both only visible on high-cardinality fields
// (the X/Y coordinates of raw Illumina headers), which is exactly where
// Genozip beats us:
//
// 1. The rescale ceiling was a hardcoded 60000 while `total` STARTS at N (all
//    frequencies initialised to 1). So any dictionary with more than ~60k
//    symbols was already over the limit before coding its first symbol, and
//    halved every frequency on every symbol thereafter -- the model could
//    never accumulate a statistic, and the dictionary degenerated to a flat
//    log2(N) code (~16.6 bits/value at N=98k) while still paying dictionary
//    overhead. Between ~3.7k and 60k symbols it was impaired rather than
//    dead. The range coder's ACTUAL constraint is tot <= 2^24 (it computes
//    range/=tot after normalising range >= 2^24), so 60000 was needlessly
//    conservative by two orders of magnitude. The ceiling now scales with the
//    alphabet, giving every dictionary a real adaptation window.
//
// 2. encSym/decSym computed the cumulative frequency by LINEAR SCAN, O(N) per
//    symbol. At N=98k over 500k names that is ~10^11 operations -- this is
//    what looked like a hang on raw-Illumina headers (it was quadratic
//    blow-up, not a deadlock), and in decSym the unbounded
//    `while(lo+freq[sym]<=target)` walk could also run past the array.
//    Replaced with a Fenwick tree (binary indexed tree): O(log N) prefix sum
//    to encode, O(log N) binary-lifting search to decode. Standard,
//    decades-proven for large-alphabet adaptive range coding.
```

What the code now says:

```
// Adaptive model over one token position's value dictionary, with cumulative
// frequencies kept in a Fenwick tree so each symbol costs O(log N). The rescale
// ceiling is max(60000, 8N), capped at 2^22, and scales with the alphabet because
// the range coder only needs the total to stay at or below 2^24. Two defects of an
// earlier version are why: a fixed 60000 ceiling made any dictionary of more than
// ~60k symbols halve its frequencies on every symbol, so it never learned and
// degraded to a flat log2(N) code, and a linear-scan cumulative frequency made raw
// Illumina headers (high-cardinality X/Y coordinates) quadratic.
```

## quality_coder.h: the trial grid, how it runs, and the one-strategy experiment that was reverted

Source: `include/quality_coder.h`, lines 189-232 before the documentation revision.

Original comment:

```
    // ── THE 8 TRIALS RUN CONCURRENTLY ────────────────────────────────────
    // This tries every (offset, strategy) pair and keeps the smallest, which
    // is 2 x 4 = 8 FULL compressions of the same block. That is why this coder
    // measured 16.6 MB/s where fqzcomp alone runs an order of magnitude faster
    // -- the throughput was never the coder, it was doing the coder's work
    // eight times.
    //
    // The eight trials are completely independent: each reads a shifted copy
    // of the block and writes its own buffer. Running them concurrently is
    // therefore FREE SPEED WITH NO SIZE COST -- the same eight candidates are
    // evaluated and the same one wins.
    //
    // DETERMINISM. The serial loop kept the first STRICTLY smaller result, so
    // ties resolved to the lowest (ci, strat). Results are collected into a
    // fixed-index array and scanned in that same order, so the selection is
    // identical regardless of completion order -- verified byte-identical.
    // CAPS_QSTRAT / CAPS_QOFF limit the trial grid, for measuring what the 8
    // trials actually buy. Unset, both are the full grid and the output is
    // bit-identical to before. Strategy order matters: the scan below keeps the
    // first strictly-smaller result, so restricting to the first k strategies
    // evaluates a PREFIX of the same ordered candidate list.
    // ── ONE STRATEGY, NOT FOUR. MEASURED 2026-09-05 ─────────────────────
    // The four fqzcomp strategies were being tried in full on every block and
    // buy NOTHING. Measured by varying the trial grid and comparing archives:
    //
    //   E. coli   2 offsets x 4 strategies  68,677,977   2 x 2  68,677,977  (same)
    //             1 offset  x 4 strategies  68,771,735   1 x 2  68,771,735  (same)
    //             1 offset  x 1 strategy    68,771,735               (same)
    //   HG002     2 x 4  574,014,786        2 x 2  574,014,786       (same)
    //
    // So strategy contributes 0.000% on both, while the OFFSET dimension is
    // worth 0.137% and is kept. That turns 8 full compressions of every 256 MB
    // block into 2 -- this coder's own comment already identified the 8x as why
    // it ran at 16.6 MB/s when fqzcomp alone is an order of magnitude faster.
    //
    // CAPS_QSTRAT=4 restores the full strategy search.
    // REVERTED TO 4 AFTER MEASUREMENT. Dropping to one strategy looked free on
    // E. coli and HG002 (both 0.000%) and is NOT: on S. acidocaldarius it costs
    // +476 B (+0.003%), isolated exactly --
    //     reference            15,159,145
    //     new code, QSTRAT=4   15,159,145  IDENTICAL  (mext + pipelining clean)
    //     new code, QSTRAT=1   15,159,621  DIFFERS
    // Two datasets are not evidence for a default that touches 82% of the
    // archive. CAPS_QSTRAT=1 remains available for anyone accepting that cost.
```

What the code now says:

```
    // Trial grid: every (offset, fqzcomp strategy) pair is compressed and the
    // smallest kept, up to 2 offsets x 4 strategies = 8 full compressions of the
    // block. The trials are independent and run concurrently, and the results are
    // scanned in a fixed order (ties go to the lowest offset index, then the lowest
    // strategy), so the winner does not depend on completion order. The strategy
    // dimension is worth 0.000% on E. coli and HG002 but 0.003% (476 B) on
    // S. acidocaldarius, so all four stay by default. The offset dimension is worth
    // about 0.14% on E. coli. CAPS_QSTRAT=1..4 and CAPS_QOFF=1..2 shrink the grid,
    // for measurement. The measurements are in docs/extras/ENGINEERING_NOTES_CODERS.md.
```
