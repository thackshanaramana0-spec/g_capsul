# GPT2026 structural COMPACT investigation

## Objective and isolation

Target: a structural, reproducible win in archive bytes, compression wall time,
and peak physical memory against SPRING, Genozip, and PgRC2 across the exact
locked set. Preserve the assembly, placements, and archive-only Claims 2/3.
The caller is out of scope. Small engineering wins are building blocks, not
completion of the publication objective.

Baseline: `dff89010a74d43c7ab79b7aeaddecce2397f7e57`, branch `gpt2026` created
from `c_star_pg_advance`. Untouched encoder built with `scripts/build106.sh`
at `/tmp/gpt2026_baseline106`. Existing benchmark status/projection files
were subsequently changed externally; preserve those changes.

All new files and flags use GPT2026/gpt2026. Substantial replacements remain
opt-in. No existing path, stream, document, or Git history is removed.

## Investigation sequence

1. Read the encoder end to end and derive archive consumers from the decoder.
   Skim the supplied Claude transcript for measurements and invocation traps.
   Its conclusions are hypotheses. Separate measured facts from code-derived
   facts and expected effects.
2. Establish the full-FASTQ baseline on E. coli and HG002, with the supported
   adaptive wrapper, `CAPS_SPANS=1 CAPS_NAMES=1 CAPS_QUAL=1`, same thread count,
   one timed job at a time, >=40 GB free disk. Record archive bytes, wall/user/
   system time, faults, kernel peak RSS and stage logs. Audit process-tree
   memory: a forked compressor's single-process RSS is not total physical use.
3. Decompose load, all three greedy sweeps, mapping index/scan/reduction, MEM,
   and coding. Use compile-time counters for hot loops; instrumentation is a
   separate build and its timing cannot establish the final speed claim.
   `perf` is currently unavailable (`perf_event_paranoid=4`, task-clock denied);
   do not describe counter measurements as hardware profiles.
4. Explain the dominant repeated work and representations before changing
   them. Read relevant original papers AND implementation for that operation.
   Implement the smallest experiment that discriminates between mechanisms,
   retaining the unmodified path. If its gate fails, record and disable it.
5. For output preservation, require nonempty archives and `cmp`, first on
   E. coli and HG002, then the locked set. Decode the archive; compare names,
   sequence, line 3 and qualities in original order. Preserve 15/15 tests.
   For output changes, run every locked dataset, accept no per-file size
   regressions and require aggregate improvement. A LOSSY result halts work
   on optimizations until its cause is established.
6. Exercise archive-only export, coverage and query; stream/pseudogenome/
   placement changes additionally require the complete HG002 Claim 2/3 gate.
   A byte-identical archive gives exactly the same consumer input but does
   not replace the requested operational checks.
7. Once a promising change is verified, repeat uncontaminated A/B timings in
   alternating order and compare competitors under equivalent preservation
   scopes. PgRC2 sequence/order results cannot be compared to full-FASTQ
   archive totals. Report unsupported datasets explicitly. Then re-profile
   and attack the next limiting structure; the three-axis goal remains open.

## Code-derived map and leads (not measured findings yet)

* Load: a complete FNV/hash-set pass estimates duplication, followed by a
  second parse/pack/dedup pass. Names and quality independently reread the
  input in a concurrent worker. Variable reads add a lexical containment
  sort. Packed reads and offsets are copied again for huge-page placement.
* Greedy sweep: sorted prefix seeds plus an open-addressed lookup; at every
  overlap length, compact open tails, verify up to eight candidates per tail
  in parallel, then commit in original order with serial fallback. Round 1
  produces the admission mask. Round 2 discards and rebuilds the links.
* Adaptive candidates fork after round 1, grouped by MINOV; members of each
  group share round 2 and mapping. The default wrapper has two MINOV groups.
  Each group currently repeats the common high-overlap part of round 2.
* Mapping: sparse read seeds are sorted and truncated by key frequency. A
  hash table and bit filter route a dense forward/RC pseudogenome scan to
  read candidates. Every logical text chunk allocates a slot for every read.
  It retains one minimum per read in hit/hmm arrays, then reduces by chunk
  order after joining. `readMM` is unchanged during each strand scan.
* Leading exact-output experiment: feed each chunk's current best mismatch
  count back into verification. Formula: candidate must beat both the
  pre-scan bound and the already observed local minimum. A zero local minimum
  terminates that read's local search. Preserve first-on-equal handling,
  chunk identities, forward-before-RC order, and existing candidate coverage.
  No extra candidate may be discarded merely because a seed is frequent.
* Broader mapping representation experiment: compact the state to indexed
  leftovers, then investigate read ownership instead of per-text-chunk
  ownership. Aim to remove replicated state and repeated verification
  together. Memory and ordering consequences need measurements first.
  **Refined after the weak atomic-state result:** transpose the candidate
  join itself. Index only pseudogenome positions whose keys occur in the
  already capped read-seed index, then assign reads to workers. Each read's
  words and best bound can stay local while testing its candidate positions;
  this removes repeated random access to the much larger packed-read store.
  Preserve exactly the original chunk seed coverage, retained read-seed
  edges and `(mismatches, seed coordinate, part)` ordering. A filtered
  position index must be measured; indexing every pseudogenome position
  blindly could outweigh the removed per-chunk state. This has not yet been
  implemented and is not an established speed/RAM result.
* MEM: sparse main-pseudogenome index, forward/RC parse with optional mismatch
  extension, candidate cap, lazy and backwards extension, then destination
  overlap trimming. Literal text and reference metadata remain addressable.
* Coding: raw streams in open_memstream buffers are copied to vectors, then
  independently transformed/coded. Names and quality bodies are copied from
  inherited globals into candidate-local structures and again into results.
  Peak memory must account for those lifetimes and copy-on-write sharing.

## Authoritative archive dependencies read in capsule_decode.cpp

* Assembly: header PGLEN/MAINEND/MINMEM; `literal`, `mem_triples`,
  `mem_dstgap`, `mem_len`, `mem_rc`, optional `mem_self`, and optional
  `mem_extmm_cnt`, `mem_extmm_pos`, `mem_extmm_obs`. Extension mismatches must
  be applied before a later reference reads the corrected source.
* Placement: `pos_abs`, optional `pos_sec`/`pos_region`, `pos_strand`,
  `read_lengths`, `orig2uid_flags`, `orig2uid_vals`. Lengths index ORIGINAL
  reads, positions index UNIQUE reads; containment requires longest unique
  lengths when restoring a reverse-strand contained original.
* Read differences: `mm_sym`, `mm_pos`, `mm_cnt` OR `mm_cnt_flags`/
  `mm_cnt_vals`; `n_indices`, `n_cnt`, `n_pos` restore ambiguous bases.
* Full FASTQ: `names_body`, `names_dict`, `names_index`, `qual_body`,
  `qual_index`; line-3 mode lives inside the names dictionary layout.
* Claim 2 full path additionally needs `contig_spans` and qualities; archive
  export also uses spans when individual contigs are requested. Coverage has
  an early exit using placements/lengths/aliases and header scalars.

## Verification issues to resolve explicitly

The existing sanity runner checks sorted sequences, not the complete FASTQ
in order. Add a separate gate; do not overwrite historical results. Read
length >255 also exceeds the existing one-byte N-position/count fields;
check real inputs for reachability before declaring any full dataset lossless.
The mapping index packs the seed-part number in three bits; NPARTS can exceed
eight for long reads. Determine the runtime domain and archive effect before
changing it. These are code observations, not confirmed dataset failures.

## Literature consulted

* [Kiełbasa et al., adaptive seeds (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3044862/):
  rarity-based seeding addresses repetitive candidate growth. It does not
  prove that replacing CAPSULE's capped seeds preserves its chosen alignments
  or archive size; sensitivity must be measured here.
* [Kowalski and Grabowski, PgRC2 (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11908645/):
  use the paper's stage decomposition to compare operation-level mechanisms.
  Local source is available at `/root/arcs-clean/method_c`; inspect it without
  modifying it. In `matching/ReadsMatchers.cpp`, best mismatch counts are
  updated during candidate processing. Check the parallel variant's scope
  before attributing any performance difference to that observation.

### Quality research follow-up (2026-09-09)

* [CRAM 3.1 (Bonfield, 2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8896640/)
  distinguishes the higher CPU cost of adaptive arithmetic coding from rANS
  and describes FQZ's configurable 16-bit context. The relevant local question
  is how much time is parameter selection, context construction and symbol
  coding. The new standalone quality probe times parameter selection and
  coding separately, using the unmodified vendored implementation. It must
  reproduce the archive's quality body before its profile is trusted.
* [Fqzcomp5 upstream](https://github.com/jkbonfield/fqzcomp5) documents sequence
  context gains on long-read technologies, periodic codec trials, and faster
  bit-packing/rANS modes. These results do not establish gains on our locked
  short-read data. They motivate measuring codec/model costs and testing
  conditional information on the actual inputs, without inheriting either
  positive or negative transfer assumptions.
* [ACO (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9175485/) investigates
  row-mean context, base context and serpentine quality traversal. Changing
  traversal is a representation experiment; any retained ordering information
  must be counted in archive size and available to the archive-only decoder.
  Published ratio gains are not evidence of a simultaneous speed/RAM win here.
