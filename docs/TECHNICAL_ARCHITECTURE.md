# CAPSULE — complete technical architecture

**Compact, Addressable, Pseudogenome-Structured, Unified Lossless Encoder.**
Written 2026-09-02 as a from-scratch, layer-by-layer account of every
mechanism currently in the archive: what it does, why it exists, and where its
code lives. This is the reference to read before touching any stage file.

Scope covered here: sequence, read order, names, quality, line 3 — the full
FASTQ (§1-8) — **plus, as of 2026-09-03, Claim 2 (variant calling, §9) and
Claim 3 (addressability, §10), both of which now live in this sandbox**
(`include/caps_caller.h` and `stages/capsule_decode.cpp`'s export/coverage/
query modes respectively), not only in the outer `/root/arcs-clean` project.
This scope line was wrong for several weeks after Claims 2 and 3 were built
here — corrected now, not carried forward.

---

## 0. The one-sentence description

A FASTQ's reads are greedily chained by suffix-prefix overlap into one
long **pseudogenome**; reads that cannot be chained are placed onto it
cheaply (pigeonhole mapping) or, failing that, assembled into a **second
region**; the two regions are then self-matched to remove residual
redundancy; and everything left over — the pseudogenome's literal bases,
every read's placement, every mismatch against its placement, and the
read order needed to reproduce the original file exactly — is split into
independent streams and each is entropy-coded by whichever of several
real coders is smallest for that stream's actual statistics.

Names and quality are separate, additive columns riding in the same
container, each with its own coder, gated on their own environment flags
so the sequence-only archive is untouched when they are off.

---

## 1. Container format

File: `stages/106_inprocess.cpp`, struct `Archive`. On-disk layout:

```
magic     "CAPSULE\0"           8 B
version   uint16                2 B
pg_len, main_pg_end  uint64 x2 16 B
minmem    uint32                4 B
n_streams uint16                2 B
per stream: namelen u8, name, len u64, payload
```

Streams are self-identifying by name, not by fixed order — the decoder
(`stages/capsule_decode.cpp`) reads them into a `std::map<string, vector<uint8_t>>`
and looks each one up by name (`S["literal"]`, `S["pos_abs"]`, etc.). Adding a
stream is additive and cannot shift any other stream's offset — the exact
property that a positional format lacked and that broke silently once (see
`docs/FAILURES_AND_REFUTED_IDEAS.md` §"positional archive format").

Header cost is ~300 B against archives of 2.6–150 MB (0.001–0.01%), and it is
counted in every reported `ARCHIVE_TOTAL`.

---

## 2. Assembly: building the pseudogenome

### 2.1 Round 1 — greedy exact suffix-prefix chaining

`stages/106_inprocess.cpp`, the sweep/chain-emission block (ported from
stages 01/05/08/09/45).

- Every read is packed 2 bits/base into a `uint64_t`-backed array.
- A seed index (default width `SEEDW`, dataset-swept) maps every `SEEDW`-mer
  to the list of reads containing it at that offset.
- For each read, the algorithm looks for another read whose prefix exactly
  matches this read's suffix for at least `MINOV` bases (also swept), and
  chains them: `read_i ... read_j` become one contiguous span in the
  pseudogenome instead of two separate literal copies.
- This is single-pass and exact — no approximate matching, no minimizers in
  the shipped path (a minimizer/tolerant variant was tried at stage 02 and
  found to trade specificity for candidates with no net win at this stage).
- **The single most important fix in the project's history** (commit
  `3e06957`): the sweep used to start at `Lmax-1`, one below the read length,
  which makes an EXACT DUPLICATE (two identical reads) structurally invisible
  to chaining — two identical 251-base reads do not overlap at L=250, only at
  L=251. Starting the sweep at `Lmax` instead let duplicates chain for free
  instead of needing a separate pre-assembly dedup pass with its own
  `orig2uid`-shaped cost. On S. acidocaldarius (9.4% duplication) this cut the
  pseudogenome 9,134,100 → 6,157,270 bytes and flipped the project's one
  remaining loss into a win.

### 2.2 Round 2 — division and second-pass chaining

Reads that chained in round 1 are divided from the leftover pool; a second,
looser sweep (`MAXMAP`, `MINOV` re-applied) tries to attach remaining reads to
the ends of already-formed chains. This is the "two-round division" of stage
09, refined by stage 14's relaxed-division variant.

### 2.3 Pigeonhole mapping

Reads that still have not chained are mapped onto the pseudogenome built so
far using a cheap, read-relative placement (stage 03): find an exact k-mer
seed hit inside the existing pseudogenome, verify against a bounded
mismatch tolerance (`MAXMM`, default 3), and record `(dst, src, len, rc)` — a
reference — instead of storing the read's bases again.

### 2.4 The MAXMAP coverage ramp

`106_inprocess.cpp`, just after the base `MAXMAP` computation. `MAXMAP` (the
per-read candidate-mapping ceiling) is widened automatically as a function of
`leftover_frac` — the fraction of reads that failed round-1 chaining, known
only after round 1 completes:

```
leftover_frac > 0.60  ->  ramp MAXMAP toward Lmax/5 (from the default Lmax/13)
```

This is a genuinely algorithmic, non-dataset-specific lever: it is keyed on a
measured property of THIS input, computed fresh every run, with a floor
(T0=0.60) set safely above every locked dataset's measured leftover_frac so
normal-coverage behaviour is provably unchanged (byte-identical, verified).
It only engages on inputs that measurably need it — e.g. low-coverage data —
and was the mechanism that closed part (not all) of the low-coverage gap to
SPRING (see `docs/COVERAGE_AWARE_MAXMAP.md`).

### 2.5 The second region

Reads that fail even pigeonhole mapping are not stored raw. They are
themselves run through a second, independent instance of the round-1/round-2
chaining process, producing a **second pseudogenome region**, appended after
the main region's end (`main_pg_end` in the container header marks the
boundary). This is what turns "unmatched leftovers" into more compressible
sequence instead of literal padding.

### 2.6 MEM self-match (cross-region and self-region deduplication)

Once both regions exist, a maximal-exact-match (MEM) pass — the real
workhorse of the size result — matches the whole pseudogenome (both regions)
against itself using an exact-match index (`copMEM`-style: hash `k1`/`k2`-mers,
verify to a full MEM), looking for repeated stretches at ANY position, not
just chain-adjacent ones. Every accepted MEM becomes a `(dst, src, len, rc)`
reference exactly like a pigeonhole placement, and the underlying pg bytes it
covers are removed from the literal stream.

**Cost-aware acceptance (`COST_GATE`, `docs/COST_AWARE_ACCEPTANCE.md`):** a
candidate MEM is only accepted if its actual encoded cost (vint-coded
src/len/gap, computed exactly — not a fixed length threshold like the older
`MINMEM`) is smaller than the literal cost ceiling (2.0 bits/base — the
information-theoretic max for 4 symbols) it would replace. This is a real,
if small, algorithmic improvement over the blunt fixed-`MINMEM` rule, off by
default where it was measured not to help and on where it was measured to.

**Second-region self-match against itself specifically was tried and
REJECTED** (`SECOND_SELF`, off by default) — see the failures document; it is
correct and lossless but removes literal bytes worth less than the references
it adds, because LZMA already finds that redundancy implicitly at a much
lower marginal cost than an explicit reference.

### 2.7 Extension mismatch tolerance

`extendTol` lambda in `parse_range`: extends a MEM match past its first
mismatch, up to `REF_MAXMM` (4) tolerated substitutions, trimming a trailing
wasted mismatch. Each accepted mismatch inside a reference is recorded as
`(position-within-match, observed-base)` rather than terminating the match —
this trades a few extra mismatch-stream bytes for a much longer, and
therefore cheaper-per-base, reference. **Measured and shipped OFF by
default** (`MEM_MAXMM = 0`): on the dataset it was tested on it cost 490,763 B
more than it saved (36,982,418 → 37,473,181 B). The code stays, correctly
gated behind `MEM_MAXMM_OVERRIDE`, purely to reproduce that measurement — see
the failures document.

---

## 3. The archive streams — what each one holds and how it is coded

`stages/106_inprocess.cpp`'s `jobs` vector is the authoritative list; each
entry is `{name, encode_fn}`, executed in parallel across a thread pool
(`NT = min(hardware_concurrency, jobs.size())`), and the resulting bytes are
what actually lands in the container. This section walks the same list.

| stream | what it holds | encoder |
|---|---|---|
| `literal` | pseudogenome bytes not covered by any reference, as 2-bit codes | `seq_encode_mem` — a DNA-specific coder (adaptive order-k context model on ACGT, shared with `seqpar_core.h` so the standalone and in-process paths cannot diverge) |
| `mem_triples` | one entry per accepted reference: `(dst, src, len, rc)` | `refc::encode` — src is delta-coded against a running "last end" cursor; the alternative reference-encoding schemes tried are in `docs/HOW_PGRC2_CODES_REFERENCES.md` |
| `mem_self` | present only if the second-region self pass fired | `best_encode` (see §4); absent stream costs nothing, existing archives without it are unaffected |
| `mem_extmm_cnt/pos/obs` | extension-tolerance mismatches (§2.7), present only if `MEM_MAXMM>0` | `best_encode` / `mmc::encode` |
| `mem_dstgap`, `mem_len`, `mem_rc` | per-reference destination gap / length-minus-MINMEM / reverse-complement flag, split into separate streams because PgRC2's own log shows the identical split ("Mismatches counts (zero flags)" / "(non-zero values)") and because each has different, much more regular statistics alone than mixed | `best_encode` each |
| `pos_abs` | per-UNIQUE-read absolute placement in the pseudogenome (uint32) | `best_encode_chunked`, `u32shaped=true` — fixed-width uint32 + byte-plane split + xz, chosen over varint after measurement (E. coli −7.9%, P. aeruginosa −11.1% vs varint; varint's length-prefix bits break the byte alignment xz's LZ77 stage exploits) |
| `pos_strand` | 1 bit per unique read: forward/reverse-complement placement | `best_encode` |
| `read_lengths` | per-ORIGINAL-read length, uint16 | `const_or_encode` — detects a file-constant length exactly and stores it ONCE (E. coli: 3,106,518 identical values → a handful of bytes, entropy 0) |
| `orig2uid_flags` / `orig2uid_vals` | maps each ORIGINAL read to its UNIQUE read id, split into a 1-bit-per-read "is this a duplicate" flag plus the sparse non-zero alias payload | `best_encode` each — the split matters because 79.55% of E. coli's deltas are exactly zero and no single model fits both the all-zero majority and the sparse 273,172-distinct-value minority well |
| `mm_ref` / `mm_obs` | for every mismatch between a placed read and its reference: the pg's byte and the read's actual byte | `mmc::encode` — an ADAPTIVE model keyed on `ref` (4 possible reference bases → 3-symbol alphabet for "not-ref"), NOT an independent per-mismatch code; this is where 40% of PgRC2's mismatch-symbol cost was beaten (124,280 vs 208,234 B on the same data) |
| `mm_pos` | position within each read where each mismatch occurs | `mmpos_encode_buckets` (bucket by per-read mismatch count, each bucket its own period-N range coder — the direct analogue of PgRC2's `compressRlMisRevOffDest` per-period streams) OR flat `best_encode`, whichever measures smaller; one leading byte says which. **Delta-coded (MMDELTA) when `Lmax<=256`**, absolute-varint-of-delta above that (fixed 2026-09-02, see the failures document — this used to be a fixed byte that silently clamped at 255) |
| `mm_cnt` | per-unique-read mismatch count | either split into `mm_cnt_flags`+`mm_cnt_vals` (same zero/nonzero split rationale as orig2uid) or flat, whichever is smaller |
| `n_pos`, `n_indices`, `n_cnt` | positions of `N` bases within N-containing reads, which original reads contain any N, and how many N's each has | `best_encode` each — N's are substituted with `A` before every other stage runs, then corrected back on the way out, so N-handling never has to touch the assembly, mapping or mismatch logic at all |
| `names_body` / `names_dict` / `names_index` | the read-ID column (§5) | `nmc::encode_from_fastq` + `best_encode` on dict/index |
| `qual_body` / `qual_index` | the quality column (§6) | `qlc::encode_from_fastq` (fqzcomp) + `best_encode` on the index |

---

## 4. The general-purpose stream selector — `best_encode`

`include/coders_inproc.h`. Given a raw byte buffer, tries several real coding
methods and keeps the smallest, prefixing one method byte so the decoder
knows which inverse to apply:

```
0/1  xz (LZMA)                          -- general purpose baseline
2    PPMd7 (LZMA SDK, public domain)    -- context-mixing, good on text-like data
3    FSE (Yann Collet, BSD)             -- tANS entropy coder, good on skewed small alphabets
4    the project's own adaptive range coder, period=1
5/6  u32 byte-plane split + xz / lzma   -- for u32-shaped streams (positions, mismatch deltas)
7    chunked: split into K spans, each independently best_encode'd and concatenated
CONST_MARKER  a stream that is one value repeated N times, stored as (count, value) -- ~9+width bytes total regardless of N
```

`const_or_encode` (used for `read_lengths`) checks the constant case FIRST,
before falling through to `best_encode`, because detecting it exactly is
cheaper and smaller than letting a general coder discover it.

**A real, general bug in this selector was found and fixed 2026-09-02**: FSE's
1-byte RLE output is not the repeated symbol, and the decoder's existing
workaround for that case only happened to be correct when the repeated
symbol was zero. See `docs/FAILURES_AND_REFUTED_IDEAS.md` for the full account
— this is now fixed by rejecting FSE/HUF's RLE result at encode time.

---

## 5. Names / read-ID column — `include/names_coder.h`, namespace `nmc`

Gated on `CAPS_NAMES=1`; unset, the archive is byte-identical to one built
without the feature.

### 5.1 Tokenizer

Ported from SPRING's real algorithm (`id_compression.cpp`, read line by line,
not inferred): each header is split left-to-right into typed tokens by
character class — alphabetic runs, digit runs, runs of the literal character
`'0'`, and single non-alphanumeric characters — and each token is compared
positionally against the same token-index slot in the PREVIOUS read's header.

### 5.2 Token types

`enum TokType { ID_ALPHA, ID_DIGIT, ID_CHAR, ID_MATCH, ID_ZEROS, ID_DELTA,
ID_END, ID_ZDELTA, ID_SEQLEN }`

- `ID_MATCH` — this token is byte-identical to the same slot in the previous
  read; costs one type symbol, no payload.
- `ID_ALPHA` / `ID_CHAR` — raw bytes, for tokens with no better representation.
- `ID_ZEROS` — a run of `'0'` characters, coded as a length.
- `ID_DELTA` — a numeric token whose value is `prev+delta` for `0<delta<256`
  (SPRING's original scheme).
- `ID_ZDELTA` — this project's own addition: a WIDER signed delta
  (`-32768..32767`, zigzag-coded into two bytes), gated behind a **self-
  learning per-token-index track record** (`hit`/`seen` counters — fires only
  once ≥20 observations show ≥30% of deltas would fall in the wide-but-not-
  ID_DELTA range for that specific token index). Without this gate a first
  naive attempt got stuck in a chicken-and-egg trap: whichever candidate
  happened to be tried first looked artificially cheap.
- **`ID_SEQLEN`** — this project's second real addition, added 2026-09-02
  after reading Genozip's real source and finding it splits `length=NNN`
  into a context separate from the rest of the header. The value of such a
  token IS the read's own sequence length, which the archive already stores
  in `read_lengths` — coding it again inside the name pays twice for one
  fact. Carries NO payload, just the type symbol; an index where this holds
  for EVERY read in the file is hoisted into the header for zero per-read
  cost at all (same admission rule as file-constant-token elision, §5.4).
  Measured: −96.5% on a dataset where this token dominated the whole names
  stream (39,416 → 1,363 B), −20.1% on another.

### 5.3 The value dictionary (`GlobalDict` / `LocalDictFreq`)

A global, per-token-index dictionary mapping observed numeric values to
compact symbol ids, built in a streaming first pass (`dict_pass`) and
serialized as part of `names_dict`. Whether a token index uses the dictionary
at all is a MEASURED COST COMPARISON in that same pass — price
`N*H(values) + header` (dictionary) against `N * sum of the four byte models'
order-0 entropies` (raw fallback) and keep the cheaper, per index. A losing
index is never registered, so an empty `GlobalDict` costs the range coder
**zero bits** (a one-symbol alphabet) — no flag, no format change needed.

**Two real defects in `LocalDictFreq` were found and fixed 2026-09-02:**
the adaptive frequency model rescaled at a hardcoded ceiling of 60,000 while
`total` STARTS at the alphabet size N — so any dictionary above ~60k symbols
was already over the limit before its first symbol and halved every
frequency on every symbol thereafter, degenerating to a flat `log2(N)` code
while still paying dictionary overhead; and cumulative-frequency lookup was
an O(N) linear scan, which is what looked like a "hang" on large-alphabet
inputs (quadratic blow-up, not a deadlock). Fixed with a Fenwick tree
(O(log N) encode/decode) and a ceiling that scales with the alphabet, bounded
by the range coder's real constraint (`tot <= 2^24`). Measures zero bytes'
difference on the currently-locked datasets because the dictionary GATE
already switches the dictionary off where it was previously broken — the fix
matters for future, larger-alphabet inputs and for the standalone stage
binaries, which predate the gate.

### 5.4 File-constant token elision

A token index whose text is IDENTICAL in every single read of the file (the
accession, instrument, flowcell, `length=` prefix text, etc.) is detected in
the same pass-1 scan and hoisted into the header once, costing zero bits per
read thereafter — not even a token-type symbol. Only engaged when every read
has the SAME token count, which is what makes the token walk terminable
without needing an explicit `ID_END` per read; files with a variable token
count (mixed instrument runs) keep the per-read path exactly as before.

### 5.5 Line 3

Genozip's real source (`qname.c`) treats line 3 (`+...`) as its own thing;
this project detects, in the same pass-1 scan, whether line 3 is (a) the bare
character `+` in every read, (b) `+` followed by an exact repeat of that
read's own header in every read, or (c) neither — and stores which as ONE
BYTE for the whole file (`Layout::line3_mode`). Case (c) is not yet supported
and the archive must not claim losslessness if it occurs (checked, not
assumed — a deliberately mixed test file was verified to correctly report
mode 2 rather than falsely claiming (a) or (b)).

### 5.6 Block structure and streaming

Both dictionary-build and encode passes stream the FASTQ header column
directly off disk (`fgets` in a loop) — never materializing the full column,
which would be ~5 GB on the largest locked dataset. Encoding is blocked
(default 250,000 names/block) and threaded via a bounded blocking queue
(reader thread pushes chunks, blocks once the queue is full; worker threads
pop and compress), bounding peak memory by queue capacity rather than file
size — the same property Genozip's real `dispatcher.c` has (read directly,
not inferred: a bounded streaming pipeline, not static up-front partitioning
like SPRING).

`names_index` carries, per block, `(compressed byte length, name count)` as
varints — a decoder holding only the archive needs both to find block
boundaries and to know how many names to decode from each block; omitting
this is the exact bug class (see failures doc) that has cost this project
twice before in other streams.

---

## 6. Quality column — `include/quality_coder.h`, namespace `qlc`

Gated on `CAPS_QUAL=1`. Wraps VENDORED `thirdparty/htscodecs/fqzcomp_qual.c`
(BSD 3-clause, James Bonfield / Genome Research Ltd), not a reimplementation
— see `docs/REIMPL_NOTES.md`'s quality section for the full measurement that
justified vendoring over reimplementing.

### 6.1 Why vendored, not reimplemented

This project's OWN quality coder (stages 67→92, an fqzcomp-style adaptive
context model: `(value history, position bucket, volatility bucket)` →
adaptive frequency table → range code) beats SPRING (7/8 datasets, −2.88%)
and Genozip (8/8, −5.24%) but loses to the REAL fqzcomp on all 8 tested
datasets by 1.1–4.5%, while also being ~2× slower single-threaded-vs-12-
threaded at the same RAM. `htscodecs` is BSD, unlike PgRC2 (GPL-3, why THAT
had to be reimplemented from scratch) — so there is no license reason not to
vendor. The vendored closure is small: `fqzcomp_qual.c` needs only
`htscodecs_tls_alloc/free` from `utils.c`, no rANS, no arith_dynamic, no
pack/rle — 2 `.c` files, 6 headers.

Every cross-column signal this project holds that fqzcomp structurally
cannot see was tested by HELD-OUT entropy (train/score split, not in-sample,
which was caught overfitting once — see failures doc) and refuted: tile
(+12.20% worse), the base call itself (+1.04% worse), an is-N flag (0.00%,
nothing). fqzcomp's own paper (CRAM 3.1, Bonfield, *Bioinformatics* 38(6),
2022) concedes a 16-bit context ceiling and small blocks, both concessions to
CRAM's random access that this project does not need — but neither costs
anything at these data volumes: the codec's own 300 MB block cap already
covers the largest locked quality column (145 MB) in a single block, and this
project's own stage-92 sweep found a BIGGER context is 5.5–12.5% WORSE at
real data volumes (under one observation per context).

### 6.2 What is fed to fqzcomp, and what is deliberately withheld

- **Record lengths** (`fqz_slice.len[]`) come from the archive's OWN quality
  scan per block, not from `read_lengths` directly, but on the decode side
  they are supplied from `read_lengths` (already decoded before the quality
  column is reached) rather than trusting anything fqzcomp itself stored —
  fqzcomp's own length storage auto-detects constant length and stores ONE
  value per block in that case, so this already costs ~nothing on
  fixed-length files; suppressing it for variable-length files was measured
  not worth patching vendored code to chase.
- **`flags[]` is deliberately left at 0**, NOT fed from `pos_strand`.
  `FQZ_FREVERSE` exists because in BAM/CRAM a reverse-strand alignment's
  quality is ALREADY stored reversed and fqzcomp un-reverses it to restore
  the true 5'→3' quality-degradation correlation. In FASTQ the quality is
  ALREADY in original read orientation, and `pos_strand` describes
  PSEUDOGENOME placement, not read orientation — feeding it would reverse
  strings that were never reversed and destroy the correlation fqzcomp
  relies on. `FQZ_FREAD2` likewise applies to paired files, which this
  single-file path does not produce.
- **The value alphabet is offset, not raw ASCII.** fqzcomp sizes its models
  on `max_sym`; raw ASCII (35–74 on real data) makes it carry nearly double
  the symbol space of the phred range it actually needs. Two candidate
  offsets — the block's own observed minimum, and the FASTQ-standard 33 —
  are tried per block and the smaller output kept (fqzcomp's own
  qmap/qshift selection is not monotone in `max_sym`, so "tighter" is not
  reliably "smaller"). The winning offset rides in the index, costing
  nothing extra in format.
- **All 4 of fqzcomp's internal strategies are tried per block**, smallest
  kept; its own default alone left 1.5–3.0% on the table against its own
  best-of-4 on 3 of 8 tested datasets.

### 6.3 Block/index shape

Mirrors names exactly: `qual_body` is fqzcomp's compressed block payloads
concatenated as-is (each block is independently, self-describingly
compressed — `fqz_decompress` needs no external strategy flag), `qual_index`
carries per block `(byte length, read count, offset byte)` as varints through
`best_encode`, so its cost is real archive bytes, never an estimate.

---

## 7. Decoder — `stages/capsule_decode.cpp`

Reads the container, inverts every general-purpose stream via
`capsule_decode_stream` (dispatches on the method byte from §4), then:

1. Rebuilds the pseudogenome: replays every `mem_triples` reference onto a
   `pg` buffer of size `pg_len`, applying any extension mismatches
   IMMEDIATELY after each reference's own copy and before the NEXT
   reference's copy can read a not-yet-corrected byte — this ordering
   subtlety (`mmc::StreamDecoder`, one symbol at a time through the same
   adaptive state the batch coder used) was a real bug found earlier in the
   project: a batch decode-then-apply approach read stale bytes whenever a
   later match's source fell inside an earlier match's destination, which is
   routine in tandem repeats.
2. Reconstructs every UNIQUE read's sequence as a slice of `pg` (forward) or
   its reverse complement (`pos_strand`), applies its recorded mismatches,
   and expands every ORIGINAL read via `orig2uid`.
3. Restores `N` characters from `n_pos`/`n_indices`/`n_cnt`.
4. If present, decodes `names_body`/`names_dict`/`names_index` — feeding the
   already-decoded `read_lengths` in as the `ID_SEQLEN` reference — to a
   `.names` side file.
5. If present, decodes `qual_body`/`qual_index` — again feeding
   `read_lengths` for record boundaries — to a `.qual` side file.

A full 4-line FASTQ is then assembled from the four decoder outputs
(`sequence`, `.names`, `.qual`, and the recovered `line3_mode`) with no
reference to the original file, and has been verified BYTE-IDENTICAL (same
MD5) against the original input.

---

## 8. Build

```bash
scripts/build106.sh   /tmp/capsule_enc     # encoder (must link -fopenmp)
scripts/build_decode.sh /tmp/capsule_dec   # decoder
```

Both scripts compile the vendored C sources (`thirdparty/htscodecs/*.c`) with
`gcc`, not `g++` — the vendored code relies on implicit `void*` conversions
that C++ rejects — then link the resulting objects into the C++ binary.

```bash
CAPS_NAMES=1 CAPS_QUAL=1 DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 \
  INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/capsule_enc \
  bash scripts/encode_adaptive.sh
/tmp/capsule_dec out.capsule <outdir> <outdir>/reads.txt
```

`DUMP_PERM=1 DUMP_MM=1` are REQUIRED, not optional debug flags — see the
failures document — they gate the per-read streams that make the archive
decodable at all.

---

## 9. Variant calling — `include/caps_caller.h` (Claim 2)

Added 2026-09-03 to this document; the code itself and its numbers predate
this addition — see `docs/CLAIM2_FINAL_VERDICT.md` for the full evaluation
history and `docs/CLAIM2_TABLES_AND_INDEL_SCAN.md` for the locked table
structure. This section covers architecture only.

### 9.0 The one-sentence description

Reference-free heterozygous SNV and indel calling from the SAME
pseudogenome and read-placement data the encoder already computed for
compression (`CallData`, populated pre-MEM in `stages/106_inprocess.cpp`,
gated on `CAPS_CALL=1`) — no separate assembler, no alignment to a
reference genome, activated as an in-process side effect of compressing.

### 9.1 `CallData` — what the encoder hands the caller

Captured **before** the MEM self-match step removes redundancy (§2.6),
because calling needs the pre-collapse contig structure, not the
maximally-deduplicated archive form: `contigs[]` (the pre-MEM chained
sequences), `read_cid[]`/`read_pos[]`/`read_rc[]`/`read_clip[]` — each
original read's contig id, offset, strand, and clip amount. A second FASTQ
pass reloads sequence+quality once `CAPS_CALL=1` is set, since the encoder's
own in-memory reads may already be partially consumed by that point.

### 9.2 Dual-substrate design — the central architectural decision

A single collapse level cannot serve both SNV and indel calling well, so
the caller builds the SAME contig set at **two different collapse
aggressiveness levels**, controlled by `dup_frac` in `collapse_contigs()`:

* **Aggressive collapse (`dup_frac≈0.45`)** for SNV pileup — merges more
  contigs together, which concentrates read depth and makes minor-allele
  fraction (`MAF=0.20`, ploidy-scaled via `MAF_K = MAF*2/PLOIDY`) estimation
  reliable on a diploid sample.
* **Mild collapse (`dup_frac≈0.92`)** for bubble/indel calling — keeps
  contigs closer to their pre-collapse form, which preserves the two-path
  bubble structure (`extract_bubble`, `extract_snv_bubble`) an indel call
  needs; over-collapsing here merges the very two haplotype paths a bubble
  call depends on distinguishing.

`build_substrate()` is called twice per input, once per collapse level, and
`Substrate` (defined at line 353) is the shared result type both paths
consume downstream.

### 9.3 The three candidate-generation channels

1. **SNV pileup** (aggressive substrate) — per-position allele counting
   against the ploidy-scaled MAF threshold.
2. **Bubble/indel extraction** (mild substrate) — `extract_bubble` walks
   right from a shared anchor comparing two contigs, `scan_pair` does a
   whole-pair gapped alignment scan for indels the anchor-walk alone
   misses, and `extract_snv_bubble` is the cross-contig SNV analogue used
   when a variant sits between contigs rather than within pileup depth on
   one.
3. **Positional-clustering indel channel** (`pcluster`, line ~1891) —
   eBWT2SNP-inspired: clusters reads by a right-context anchor that is
   unique across the contig set, catching indels neither pileup nor bubble
   extraction reaches. This channel, plus read-level junction support and
   closing-anchor (both-ends-anchored bubble) detection, are the three
   mechanisms that moved indel F1 from ~0.36 to 0.637 over this project's
   history (`docs/HOW_DISCOSNP_WINS.md` §3).

### 9.4 Multi-allelic emission

`CAPS_PLOIDY=k` admits up to `k` co-occurring alleles at a site and emits a
native multi-allelic VCF record (`ALT=C,G`), rather than DiscoSNP++'s
behavior of emitting separate biallelic records plus invalid alt-vs-alt
rows for the same site — a real capability difference, not a scoring
artifact (verified real-vs-synthetic distinction in
`docs/POLYPLOID_BENCHMARK.md`).

### 9.5 Why filtering alone cannot close the remaining indel gap

The structural finding from `docs/HOW_DISCOSNP_WINS.md` §4, load-bearing
enough to repeat here: CAPSULE's indel candidates are built FROM contigs,
which are built FROM reads — every candidate is read-supported by
construction, so a read-validation filter (the mechanism that gives
DiscoSNP++'s `kissreads2` its precision) has nothing to reject. DiscoSNP++'s
candidates are graph-traversal *hypotheses* that can be unsupported by any
single read; ours cannot be. Four filter attempts confirmed this
empirically (net negative or neutral, every time) — closing the gap
requires a different candidate *generator* (both-ends k-mer anchoring, the
eBWT/dBG route in `docs/CALLER_ARCHITECTURE_PLAN.md`), not another filter
on top of the current one.

### 9.6 Testing

`scripts/test_claim2.sh` (added 2026-09-03) — a synthetic diploid genome
with known het SNVs and clean (non-homopolymer) het indels, run through the
real caller → bwa-alignment → `lift_vcf.py` pipeline, checked for recall.
This is the first automated test this 2132-line file has ever had; see
`docs/INDUSTRIAL_CHECKLIST_CLAIM2.md` for the bug this test itself had (and
caught in itself) before being trusted.

---

## 10. Addressability — `stages/capsule_decode.cpp` export/coverage/query (Claim 3)

Added 2026-09-03; code and numbers predate this addition — see
`docs/CLAIM3_LOCKED.md` for the full evaluation, bug history, and prior-art
position. This section covers architecture only.

### 10.0 The one-sentence description

Three operations a conventional pipeline computes from scratch — assemble a
genome, align reads to compute per-base depth, index reads for coordinate
lookup — are instead served by decoding archive streams the compressor
already wrote for its own purposes, because the pseudogenome (§2) already
**is** the assembly, and `pos_abs`/`read_lengths` (§3) already **are** the
placement index.

### 10.1 Three early-exit modes, one function

`capsule_decode_all()` (`stages/capsule_decode.cpp:164`) takes a `mode`
argument and returns as soon as the streams that mode needs are decoded —
none of the three pays for full read reconstruction:

* **`export`** (line 318) — decodes `literal` + `mem_triples` (+ the
  mismatch-override streams) to rebuild the pseudogenome byte array, emits
  it as FASTA, stops. No per-read stream touched.
* **`coverage`** (line 195) — hoisted **above** the pseudogenome rebuild
  entirely: needs only `PGLEN` (header) and `pos_abs`/`read_lengths`, not
  one byte of pg content. A difference-array depth computation,
  O(reads + PGLEN), not O(reads × length).
* **`query`** (line 390) — the one mode that DOES need the rebuilt pg (it
  must return actual sequence), then does a linear overlap scan of
  `pos_abs` against the requested range, one record per UNIQUE read
  (duplicates share a placement and would return byte-identical records).

### 10.2 The `orig2uid` indexing invariant — and the bug it caused

`pos_abs` is indexed by UNIQUE read; `read_lengths` is indexed by ORIGINAL
read, always (an invariant `106_inprocess.cpp` states explicitly). Walking
both with one shared loop counter is only correct when there are zero
duplicate reads — the first version of `coverage` did exactly this and
silently dropped every duplicate read's contribution to depth, a 20%
undercount on E. coli, found and fixed this session
(`stages/capsule_decode.cpp:208-239`, expanding through `orig2uid_flags`/
`orig2uid_vals` exactly as the main read-reconstruction path already did).
Kept here as a permanent architectural note, not just a changelog entry,
because this exact asymmetry is a standing trap for any future code that
touches both streams together.

### 10.3 Limitation: `query` is O(archive), not O(range)

No range index exists over `pos_abs` (no interval tree, no sorted-offset
table) — reads are stored as slices of the pg, so `query` must rebuild the
whole pg before answering any range. Measured, corrected 2026-09-03 (an
earlier figure compared against the wrong baseline — see `CLAIM3_LOCKED.md`
§6.4): query and full reconstruction cost about the same in time (1.62×
saving), and the real advantage is **output selectivity** — 132× fewer
reads, 112× fewer bytes returned for a narrow range on E. coli — not
asymptotic speed.

### 10.4 Testing

`scripts/test_claim3.sh` (added 2026-09-03) — invariant checks (export is
pure ACGT, coverage's total covered-base-units equals the true sum of read
lengths, query's full-range record count equals the encoder's own unique
count) rather than hand-computed expected pg bytes. Verified to actually
catch the §10.2 bug: running the pre-fix binary against this test's
synthetic input reproduces an 18.8% coverage undercount, the same failure
class measured on real data.
