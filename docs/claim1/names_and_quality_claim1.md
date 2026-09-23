---
Date: 2026-09-20
Title: Names, Quality, and Line 3 — the Three Non-Sequence Columns of Claim 1
Purpose: `mechanism_insight_claim1.md` covers Claim 1's DNA-sequence stream
  and its novel-insight comparison to PgRC2 in depth, but treats names and
  quality as a single passing table row ("Full columns — SPRING-derived
  tokenizer; vendored fqzcomp"). This file is the missing depth: how the
  names tokenizer, the quality coder, and line 3 actually work, condensed
  and verified against `docs/TECHNICAL_ARCHITECTURE.md` §5-6 (the project's
  own from-scratch reference), not rewritten from memory.
When to refer to this file: Writing the Materials and Methods description
  of the full 4-column FASTQ reconstruction (sequence + names + quality +
  line 3); explaining why quality is vendored rather than reimplemented;
  explaining the two real named/token-coding contributions
  (`ID_ZDELTA`, `ID_SEQLEN`) that are this project's own, not SPRING's.
Keywords: names, read IDs, tokenizer, ID_SEQLEN, ID_ZDELTA, quality,
  fqzcomp, htscodecs, vendored, line 3, CAPS_NAMES, CAPS_QUAL
---

# Names, quality, and line 3

Claim 1's scope, precisely: **sequence + read order + names + line 3** are
this project's own implementation; **quality** is the one column vendored
rather than reimplemented (`include/quality_coder.h`, wrapping
`thirdparty/htscodecs/fqzcomp_qual.c`, BSD 3-clause). Both are gated
(`CAPS_NAMES=1`, `CAPS_QUAL=1`) and the archive is byte-identical to the
sequence-only path with them unset — nothing here can move a Claim 1
sequence-only number.

## Names — `include/names_coder.h`, namespace `nmc`

**Tokenizer**: ported from SPRING's real algorithm (`id_compression.cpp`,
read line by line, not inferred). Each FASTQ header is split left-to-right
into typed tokens by character class, and each token is compared
positionally against the same slot in the *previous* read's header — a
`ID_MATCH` token (byte-identical to the previous read, one symbol, no
payload) is the cheap common case for the run-length structure real
instrument headers have.

**Two additions that are this project's own, not SPRING's:**

- **`ID_ZDELTA`** — a wider signed delta (zigzag-coded into two bytes,
  range −32768..32767) than SPRING's original `ID_DELTA` (`prev+delta` for
  `0<delta<256` only). Gated behind a self-learning per-token-index track
  record (`hit`/`seen` counters, fires only after ≥20 observations show
  ≥30% of deltas would fall in the wide-but-not-`ID_DELTA` range for that
  specific index) — a naive first attempt without this gate got stuck in a
  chicken-and-egg trap where whichever candidate was tried first looked
  artificially cheap.
- **`ID_SEQLEN`** — added 2026-09-02 after reading Genozip's real source and
  finding it splits a `length=NNN` field into its own context. The value of
  such a token IS the read's own sequence length, which the archive already
  stores in `read_lengths` elsewhere — coding it again inside the name pays
  twice for one fact. Carries **no payload**, just the type symbol; an index
  where this holds for every read is hoisted into the header for zero
  per-read cost. **Measured −96.5% on the dataset where this token
  dominated the names stream** (39,416 → 1,363 B), −20.1% on another.

**Value dictionary** (`GlobalDict`/`LocalDictFreq`): a per-token-index
dictionary of observed numeric values, built in a streaming pass. Whether an
index uses the dictionary at all is a measured cost comparison — dictionary
cost vs. raw order-0 entropy fallback — per index, so a losing index costs
the range coder zero bits. Two real defects found and fixed 2026-09-02: the
frequency model's rescale ceiling (hardcoded 60,000) was below the starting
`total` for large alphabets, degenerating the model silently; and
cumulative-frequency lookup was O(N) linear, which looked like a hang on
large-alphabet inputs (quadratic blow-up, not a deadlock) — fixed with a
Fenwick tree. Zero byte impact on the currently-locked datasets (the gate
already switched the dictionary off wherever it was broken); matters for
future larger-alphabet inputs.

**File-constant token elision**: a token index identical in every read
(accession, instrument, flowcell text) is hoisted into the header once,
costing zero bits per read thereafter — only engaged when every read has the
same token count.

**Streaming, not in-memory**: both the dictionary-build and encode passes
stream the header column directly off disk, never materializing the full
column (~5 GB on the largest locked dataset). Blocked (250,000 names/block),
threaded via a bounded queue — the same streaming-pipeline property Genozip's
real `dispatcher.c` has, verified by reading it, not inferred.

## Quality — `include/quality_coder.h`, namespace `qlc`

**Why vendored, not reimplemented**: this project's own quality coder
(stages 67→92, an fqzcomp-style adaptive context model) beats SPRING (7/8
datasets, −2.88%) and Genozip (8/8, −5.24%) but **loses to the real fqzcomp
on all 8 tested datasets by 1.1–4.5%**, while running roughly 2× slower
single-threaded against fqzcomp's 12-threaded path. `htscodecs` is BSD —
unlike PgRC2 (GPL-3, the reason PgRC2's assembler specifically had to be
reimplemented from scratch, not vendored) — so there is no license barrier to
vendoring here. The vendored closure is small: `fqzcomp_qual.c` needs only
`htscodecs_tls_alloc/free` from `utils.c` — 2 `.c` files, 6 headers, no rANS,
no arith_dynamic, no pack/rle.

**Cross-column signals were tested and refuted**, by held-out entropy
(train/score split — an in-sample version of this test overfit once and was
caught, see the project's refuted-ideas record): tile number (+12.20%
worse), the base call itself (+1.04% worse), an is-N flag (0.00%, no
effect). fqzcomp's own published concessions (CRAM 3.1, Bonfield,
*Bioinformatics* 38(6), 2022 — a 16-bit context ceiling and small blocks,
both trade-offs for CRAM's random-access requirement) cost nothing here: the
codec's 300 MB block cap already covers the largest locked quality column
(145 MB) in one block, and this project's own sweep found a *bigger* context
is 5.5–12.5% worse at these data volumes (under one observation per
context) — the opposite direction CRAM's constraint would suggest mattered.

**What is fed to fqzcomp, and what is deliberately withheld** — three
real, verified decisions:
1. Record lengths are supplied from the archive's own `read_lengths` on
   decode, not trusted from fqzcomp's own auto-detected block storage.
2. **`flags[]` is deliberately left at 0**, not fed from `pos_strand`.
   `FQZ_FREVERSE` exists in BAM/CRAM because a reverse-strand alignment's
   quality is *already* stored reversed there. In FASTQ, quality is already
   in original read orientation, and `pos_strand` describes *pseudogenome*
   placement, not read orientation — feeding it would incorrectly reverse
   strings that were never reversed and destroy the 5'→3' degradation
   correlation fqzcomp relies on.
3. The value alphabet is offset (two candidate offsets tried per block, the
   smaller output kept), since fqzcomp's own model sizing is not monotone in
   the symbol-space size — "tighter offset" is not reliably "smaller
   output." All 4 of fqzcomp's internal strategies are tried per block, best
   kept; its own default alone left 1.5–3.0% on the table against its own
   best-of-4 on 3 of 8 tested datasets.

## Line 3 — one byte for the whole file

Genozip's real source (`qname.c`) treats line 3 (`+...`) as its own thing.
This project detects, in the same pass-1 header scan, whether line 3 is (a)
the bare `+` in every read, (b) `+` followed by an exact repeat of that
read's own header, or (c) neither — and stores which as **one byte for the
whole file**. Case (c) is not yet supported, and the archive is checked
(not assumed) to correctly report mode 2 rather than falsely claim (a) or
(b) on a deliberately mixed test file.

## A real, fixed bug found by testing all three claims' features together

Verified against `docs/FINAL_ALGORITHMIC_SCAN.md` (2026-09-03). Every
per-claim test built its own archive with only its own claim's flags set;
none had exercised `CAPS_NAMES=1 CAPS_QUAL=1` through the *recommended*
compression path (`encode_adaptive.sh`'s candidate sweep) together with a
relative input path — the normal, documented invocation
(`INPUT=reads.fq`). Doing so produced an archive that silently decoded to
**zero names, zero quality bytes** — no error, no warning, normal exit
code. Root cause: the candidate-sweep fork `chdir()`s into a per-candidate
scratch directory before the names/quality encoders reopen the input file
to extract that data; two other environment-derived paths (`CALL_VCF`,
`CAPS_DUMP_CONTIGS`) had already been fixed for this exact hazard once
before, but the input-path variable the names/quality path depends on had
not received the same treatment — the second occurrence of the identical
bug *pattern*, not a new class of bug. Fixed by resolving the input path
to an absolute path once, immediately after reading `argv[1]`, before any
fork can occur — the one change that covers both the `CANDIDATES` and
`GSEARCH` sweep code paths at once. **Verified, not assumed**: byte-identical
output pre-fix vs post-fix for every path that does not exercise the
sweep-plus-names/quality combination — meaning **no currently-published
number in this project was ever affected**, since no published result uses
names/quality through a candidate sweep — and `scripts/verify_lossless.sh`
on a real locked dataset with the fixed binary reproduces the exact
previously-measured archive size byte for byte.

## Decode order, and the one real ordering bug it had to avoid

`src/decoder.cpp` rebuilds the pseudogenome first (replaying every
`mem_triples` reference, applying extension mismatches *immediately* after
each reference's own copy — a batch decode-then-apply approach was found to
read stale bytes whenever a later match's source fell inside an earlier
match's not-yet-corrected destination, which is routine in tandem repeats),
then reconstructs sequence, restores `N`s, and only then decodes names
(feeding it the already-decoded `read_lengths` for `ID_SEQLEN`) and quality
(same, for record boundaries). A full 4-line FASTQ assembled from all four
decoder outputs, with no reference to the original file, has been verified
**byte-identical (same MD5)** against the original input.
