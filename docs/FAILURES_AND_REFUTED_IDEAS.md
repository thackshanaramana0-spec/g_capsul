# CAPSULE — failures, refuted ideas, and real bugs found

Written 2026-09-02. Every entry below was IMPLEMENTED AND MEASURED, not
argued away — that is the standing discipline for this repo (`CLAUDE.md`
rule: "a change that fails its gate is reverted and recorded, not tuned
until it passes"). This document exists so none of them is retried.

Two categories, kept separate: **real correctness bugs** (silent data loss,
now fixed) and **algorithmic ideas that were correct but did not pay**
(measured, reverted, kept off).

---

## PART A — Real correctness bugs (silent data loss), found and fixed 2026-09-02

All four were found the same day, while wiring the quality column, when the
user pushed back on one dataset's failure ("are you sure that's a bug") and
that led to a full round-trip audit of all 14 locked datasets — the first
time anything had ever DECODED an archive and compared it to the original
FASTQ, dataset by dataset. Before this, Phase 1 and Phase 2b had measured
archive SIZE only.

### A.1 — Mismatch position byte clamp (commit `23be207`)

**Where:** `stages/106_inprocess.cpp`, the extension-mismatch emission loop
(`fputc((uint8_t)(j>255?255:j), fp)`), and its decoder counterpart in
`stages/capsule_decode.cpp`.

**What:** mismatch positions within a read were stored in ONE byte. For any
read longer than 256 bases, positions past base 255 were silently CLAMPED to
255 — the archive recorded a wrong position for every such mismatch, with no
error anywhere. The code's own prior comment called this "the pre-existing
>255bp cap limitation" — known, but its consequence (which locked datasets
actually exceed 256 bp) was not checked.

**Trigger:** any read >256 bp. Two locked datasets: ERR552797 (M.
tuberculosis, up to 301 bp) and SRR40271341 (H. pylori, 300 bp fixed).

**Effect measured:** both datasets were reported as compression WINS from
archives that could not reproduce their own input.

**Fix:** the same per-read delta is now written as a VARINT instead of a
fixed byte, so any position is representable. The bucketed `mm_pos` form
(which indexes one byte per mismatch) is forced off above Lmax=256, since
that indexing does not survive varints; the flat form is used instead.

**Verification:** the `Lmax<=256` code path is untouched by construction
(same byte-per-mismatch encoding as before) — verified BYTE-IDENTICAL
archives on E. coli and SARS-CoV-2 pre/post fix. H. pylori: FAIL → LOSSLESS.

---

### A.2 — Mismatches emitted for orphaned unique reads (commit `121fea9`, bug 1)

**Where:** `stages/106_inprocess.cpp`, the mismatch-emission loop — it
iterated over every UNIQUE read index without checking whether any ORIGINAL
read actually maps to it.

**What:** under containment (a shorter read is fully contained within, and
therefore represented by, a longer read at the same pseudogenome position),
a "unique" id can end up with NO original read pointing to it — the shorter
read got absorbed differently, or its own placement was chosen instead. Such
an orphaned unique is never emitted by the decoder — but the ENCODER still
walked it, computed its mismatches with an assumed length of 0 (no original
to take a length from), and wrote those wrong `(ref, obs, pos)` triples into
the shared `mm_ref`/`mm_obs`/`mm_pos` streams anyway.

**Why this corrupted OTHER reads:** the mismatch symbol coder (`mmc::encode`)
is ADAPTIVE — its frequency tables update after every symbol, in stream
order. A handful of wrong `ref` bytes (computed from a bogus zero-length
placement) desynchronised the model's state for every symbol emitted after
them. On a 100,000-read C. jejuni sample: 901 of 81,648 unique reads were
orphaned, producing just 8 bad `ref` bytes out of 18,774 total — but those 8
bad refs cascaded into 6,892 wrong `obs` bytes and ultimately 3,712 wrong
reads (3.7% of the file), most differing from their true sequence in only
1–3 bases (61% of the wrong reads had exactly one base wrong — the
signature of a coder that had drifted, not a placement that was grossly
wrong).

**Trigger:** containment, which can only occur with variable-length reads
(two identical-length reads cannot contain one another). This is exactly
why every fixed-length dataset round-tripped from the start and only
variable-length ones did not.

**Fix:** build a `uid_referenced` bitmap from `orig2uid` before the mismatch
loop runs, and skip any unique index with no original pointing to it.
Strictly better on two axes at once: it removes a correctness break AND
removes pure waste (coding mismatches for a read that is never output).

**Verification:** C. jejuni sample 3,712 wrong reads → 95 (see A.3 for the
remaining 95). 901 orphaned uniques correctly skipped.

---

### A.3 — Reverse-strand contained reads decoded from the wrong end (commit `121fea9`, bug 2)

**Where:** `stages/capsule_decode.cpp`, the parallel read-reconstruction loop
— it indexed the pseudogenome slice with the ORIGINAL read's own length `L`
for both forward and reverse-complement reads.

**What:** `pg[pp .. pp+RLu)` holds the reverse complement of the UNIQUE
read (length `RLu`), not of any one original read that maps to it. A
contained original read is a PREFIX of the unique read in READ space — but
because reverse-complementing reverses the string, that same contained
region is a SUFFIX of the slice in PG space. Indexing with the original's
own (shorter) length `L` instead of the unique's length `RLu` therefore read
bases from the WRONG END of the reverse-complemented slice for every
contained, reverse-strand read.

**Signature, found by direct inspection (not guessed):** of the 95 reads
still wrong after fix A.2, ALL 95 were reverse-complement AND ALL 95 were
strictly shorter than their unique read's length — a 100% clean correlation
that made the mechanism unambiguous before writing the fix.

**Same bug exists in `scripts/decode_105.py`** (the Python reference
decoder, used by `verify_lossless.sh`) — confirmed by observing it fail on
exactly the same 95 reads, not a coincidence.

**Fix:** compute `RLu = max(rlenU[u], L)` and index the reverse-complement
slice with `RLu`, not `L`. Forward reads are prefixes in BOTH read space and
pg space and were never affected — the fix is a no-op for them.

**Verification:** C. jejuni sample: 95 → 0. Full sample now LOSSLESS.

---

### A.4 — FSE/HUF RLE round-trip corruption (commit `f3ab0c2`)

**The most general of the four — not tied to read length, strand, or any
dataset-specific property.**

**Where:** `include/coders_pgrc.h`, `fse_encode`/`huf_encode` (the accept
condition `r==0 || r>=n` let `r==1` through) and the decoder's existing
special case for it.

**What:** `FSE_compress` (and `HUF_compress`) return exactly 1 byte when the
input is a single repeated symbol (RLE). That ONE byte is a compressed-form
artifact, NOT the repeated symbol itself. An existing decoder workaround
reconstructed such a stream as "`rawlen` copies of `payload[0]`" — a fix that
had been written for, and only tested against, `orig2uid_flags`, which
happened to be 57,563 ZERO bytes; `payload[0]` there was coincidentally `0`,
the correct symbol, so the workaround looked correct. For ANY OTHER constant
value, this silently decoded the entire stream as ZEROS.

**Trigger, found in the wild:** M. tuberculosis (ERR552797) has 3 reads that
are ENTIRELY `N`. Their per-read N-count stream (`n_cnt`) was therefore 3
bytes of the constant value `0x23` (35 — the read length, since every base
is N). `FSE_compress` returned 1 byte, `0x00`. The archive decoded that
column as three ZEROS, and all three all-N reads came back as all-`A` — a
complete, silent loss of every N in those reads, with the archive reporting
no error whatsoever.

**Why this is more dangerous than A.1–A.3:** it required no specific read
length, no specific strand, no containment — only ANY stream, anywhere in
the archive, that happens to be constant and non-zero. It could in principle
have hit any of the other `best_encode`-coded streams under different data.

**Fix:** reject FSE/HUF's own RLE result at encode time (`r<=1` instead of
`r==0`), so the stream selector falls through to a coder that round-trips
correctly. This costs essentially nothing: a genuinely constant stream is
ALREADY handled correctly and far more cheaply by `const_or_encode`
(9 bytes regardless of length) wherever that path is used, and for streams
not routed through `const_or_encode`, any other coding method handles a
1–3-byte input just as compactly. The decoder's existing "fill with
`payload[0]`" special case is left in place (harmless, and needed to decode
any archive already written through the old, zeros-only-safe path).

**Verification:** M. tuberculosis 20k-read sample: FAIL → LOSSLESS.
E. coli, which also had a zeros-constant stream taking the old RLE path,
changed by exactly +33 bytes (+0.0009%) and was never actually corrupted
(the old path happened to be correct there) — confirming the fix's cost is
negligible and its scope is exactly the streams that were actually broken.

---

### Why none of these four were caught earlier — the process gap

Two structural gaps, both now understood precisely:

1. **Phase 1 and Phase 2b benchmarks measured archive SIZE only.** A wrong
   archive and a correct archive of similar size look identical on a size
   table. Nothing in either benchmark ever decoded an archive and compared
   it to the source FASTQ.
2. **The one lossless check that existed, `verify_lossless.sh`, does not
   test the archive.** It decodes the encoder's DUMPED intermediate streams
   (`decode_105.py` reading `.bin` files written alongside the archive), not
   the `.capsule` container itself — so the entire entropy-coding layer,
   where bug A.4 lived, was never exercised by it at all. This is explicitly
   flagged in `CLAUDE.md` §6.1 as a known gap, dated BEFORE these four bugs
   were found — the warning existed; the check that would have caught it did
   not.

That check had also only ever been run on 7 of the then-locked datasets
(E. coli, SARS-CoV-2, P. aeruginosa, H. salinarum, L. major, S. aureus, S.
acidocaldarius), none of which are among the four that turned out broken.

**The fix to the process, not just the code:** a full-audit script now
decodes every dataset and diffs its sequence column against the original
FASTQ before any size number is allowed to be reported. All 14 locked
datasets confirmed LOSSLESS after the four fixes above.

---

## PART B — Algorithmic ideas: implemented, measured, and REJECTED

Each of these is a real, working implementation that was measured to make
the archive larger, slower, or otherwise worse, and is kept in the codebase
(usually behind an environment-variable gate defaulting OFF) purely as the
evidence for why it is not the shipped default.

### B.1 — Second-region self-match

`run(Q, qlen, CROSS, ...)` matches the second region against the MAIN
pseudogenome only, never against itself. Implemented twice (two different
mechanisms), both structurally correct and lossless. Measured: removes
literal bytes worth only ~212 KB on one dataset while the references needed
to remove them cost ~369 KB — net loss, because LZMA already captures that
same redundancy implicitly at ~0.24 bits/base, far cheaper than an explicit
21.5-bit-average reference. On SARS-CoV-2, 74.3% of the second region's
32-mers are repeat occurrences, yet the self-match still nets negative once
reference cost is counted honestly. Behind `SECOND_SELF`, off by default.

### B.2 — PgRC2's both-side-overlap admission rule

PgRC2 only admits a read to the pseudogenome if it overlaps well-tiling
neighbours on BOTH sides. Implemented and measured on 3 datasets: worse on
all 3 (+2,845, +45,746, +78,597 bytes). Behind `BOTHSIDE`, off.

### B.3 — Extension mismatch tolerance in MEM matching (MEM_MAXMM)

Extending a MEM match past its first mismatch, tolerating up to
`REF_MAXMM=4` substitutions per match. Measured: 36,982,418 → 37,473,181
bytes, a 490,763-byte LOSS on the tested dataset — the extra mismatch-stream
cost exceeds the reference-length savings. Code kept, hardcoded OFF
(`MEM_MAXMM=0`), reproducible only via `MEM_MAXMM_OVERRIDE`.

### B.4 — Read-scale (SPRING-style) conversion of region-scale MEM matches

Investigated as the leading hypothesis for closing the low-coverage gap to
SPRING: could region-scale MEM references (21–38 bits each) be replaced with
cheap SPRING-style read-relative mismatch coding (~2 bits/mismatch)?
**Refuted structurally, not just numerically**: measured average MEM match
length is 44.9 bases against a ~150-base read — most matches are shorter
than a read, so there is no single "read-relative" frame to convert to. The
mechanisms operate at genuinely different scales and cannot be unified by
this route.

### B.5 — BSC-on-raw as a literal-coder sanity check

Compared this project's literal coder against BSC (independent, BWT-based,
general-purpose) applied directly to the raw pseudogenome bytes, to check
whether the DNA-specific coder was adding real value or just doing what a
strong general compressor would do anyway. Result: we already beat raw BSC
by 1,575,156 bytes on the tested file — the DNA-specific model earns its
keep. (This was a validation that PASSED, listed here because it was one of
the "is there a cheap universal alternative" checks run during the
low-coverage investigation, alongside B.4 which failed.)

### B.6 — Per-region MINMEM (separate MEM thresholds for main vs second region)

Monotonically worse at every tested value: +30,652 bytes at threshold 32,
+2,010,039 bytes at threshold 128. The 431K second-region references each
pay for themselves individually; splitting the threshold does not help.

### B.7 — Splitting positions by region

Checked whether coding `pos_abs` separately for main-region vs second-region
reads would beat the unified stream. It does not — the unified stream
already codes BELOW the region-split bound (3,629,158 vs 3,699,725 bytes on
the tested file).

### B.8 — MINOV as a search-optimizable lever

`MINOV` (minimum overlap for round-1/round-2 chaining) was checked for
whether a search procedure could find a better value than the current
default. Already optimal on both files tested, AND the response curve is
NOT unimodal — no local or gradient-style search applies to it at all; only
an exhaustive sweep can be trusted.

### B.9 — Minimum-degree-first matching order

Hypothesis: processing chain candidates in minimum-degree-first order (a
classic graph-matching heuristic) might resolve ambiguous chaining better.
Measured: byte-identical output to the existing order. 96.8% of chain tails
have ZERO competing candidates at any given overlap level, and under 1% have
more than one — there is essentially no contention for this heuristic to
resolve on real data.

### B.10 — Fixing the tail-eligibility bug

A REAL bug was found: reads shorter than the sweep's current starting length
`L` are permanently ineligible to ever be considered as chain tails, even
later in the sweep when `L` has shrunk to their size. Confirmed as a genuine
defect. But FIXING it measured WORSE on both datasets tested (C. jejuni
+26,109 bytes, SARS-CoV-2 +4,178 bytes) — such short reads are cheaper
handled by mapping onto an already-formed chain than by being assembled
into it. Left unfixed, deliberately, with this measurement as the reason.

### B.11 — Full-length links as containment

A guard intended to treat a full-length overlap as containment (skip
further chaining) was checked: the guard never actually fires on real data.
Archives are identical with and without it.

### B.12 — Parallel coder probes (running all `best_encode` candidate methods concurrently per stream)

6–11% faster wall-clock, for +54–86% peak RAM. Not shipped as the default
given the RAM cost was judged disproportionate to the speed gain (this
project's benchmarks are not primarily RAM-constrained at current scale, but
the trade was not judged worth taking blindly).

### B.13 — Architecture rewrite (a full internal restructuring, earlier session)

A significant internal architecture change was implemented and measured:
0% size change, 0% speed change, +54 MB RAM. The REAL fix that actually
produced a size improvement in the same investigation window was unrelated
and much smaller: aligning LZMA's `lp` parameter to integer stride
boundaries, which alone was worth +1.46% vs PgRC2 (3/5 wins) — a reminder
that a structural rewrite is not a substitute for finding the actual
bottleneck.

### B.14 — X/Y coordinate fixed-width columns for the names coder (2026-09-02)

Proposed after reading Genozip's real `dyn_int.c` (numeric fields stored as
a width-minimised binary column into a general codec) and because this
project had already proven the same PRINCIPLE on `pos_abs` (§ Technical
Architecture, fixed-width uint32 + xz beating varint by 7.9–11.1%).
**Measured before building anything**, using a standalone column-cost
harness: best case −0.4% against the current inline coder, WORSE on both Y
columns tested. Entropy analysis explains why: the current inline coder is
already BELOW order-0 entropy on 3 of 4 tested X/Y columns and within
0.3–4.6% of DELTA entropy on 3 of 4 — there is no room left for a
re-representation trick to find, because X+Y are 99.8% of the entire names
stream on Illumina headers and are already coded near the theoretical bound.
By the same argument, Genozip cannot be beating this project through a
better coordinate representation either — there is no representational
headroom for either side to claim.

### B.15 — Magnitude-bucket unification of ID_DELTA/ID_ZDELTA

To remove the per-token TYPE-SYMBOL cost (measured by direct bit
attribution to be 7.8–10.5% of the names stream, 105–125 KB on tested
files) — a real, quantified cost that Genozip pays none of, since each of
its fields is its own dedicated context. Tried unifying the two delta token
types into one type with a magnitude-bucket-coded payload. Measured
**+28,826 bytes WORSE**. The type bit is NOT redundant: the specialised
per-type models it lets the coder select between are worth more than the
bit costs. (Removing `ID_ZDELTA` outright, as a simpler test of the same
idea, was far worse still: 1,354,682 → 1,698,147 bytes.)

### B.16 — Raising `Model`'s rescale ceiling (the OTHER adaptive model, not `LocalDictFreq`)

Tried scaling `Model`'s (used for token types, alpha bytes, delta bytes,
etc. — separate from the dictionary's `LocalDictFreq`) rescale ceiling from
60,000 upward, on the theory that a slower-decaying model might capture
more structure. Result: −4,790 bytes on one dataset, **+1,878 bytes on
another** — fails the project's own no-regression gate. Correctly classified
and rejected as a PARAMETER, not an algorithmic fix, precisely because its
sign was not consistent across datasets.

### B.17 — Quality: reimplementing fqzcomp instead of vendoring it

Considered explicitly, and the paper's own conceded limitations (16-bit
context ceiling, small blocks — both concessions to CRAM's random-access
requirement) were checked as potential attack surfaces for a from-scratch
reimplementation that could beat it. Neither pays at this project's data
volumes: the codec's own 300 MB block cap already covers every locked
dataset's quality column (largest measured: 145 MB) in ONE block, and this
project's own stage-92 sweep independently found that a BIGGER context is
5.5–12.5% WORSE at real data volumes (models become too sparse). Combined
with htscodecs' permissive BSD license (making vendoring legally free,
unlike PgRC2's GPL-3), reimplementation was rejected in favour of vendoring.

### B.18 — Cross-column signals for quality: tile, base call, is-N flag

Tested by HELD-OUT entropy (train on even-indexed reads, score odd-indexed
reads, Laplace-smoothed) whether any signal outside the quality column
itself — specifically ones fqzcomp structurally cannot see, since it only
ever receives the quality string — could improve on it:

    conditioning on TILE (from the name column)   +12.20% WORSE
    conditioning on the base call (ACGTN)          +1.04% WORSE
    conditioning on an is-N flag                    0.00%  nothing

**A methodology trap recorded here because it nearly produced a false
positive:** the FIRST measurement of tile conditioning, done IN-SAMPLE (same
data used to build and score the model), showed −3.02% — apparently better.
This is pure overfitting: adding any context strictly lowers in-sample
entropy, always, regardless of whether the context carries real signal. A
naive first held-out attempt (simple first-half/second-half split) then
showed +227% — ALSO wrong, because tile values are sequentially ordered in
the file, so a first-half/second-half split tests the model on tile values
it has never seen at all. Only an INTERLEAVED held-out split (alternating
reads into train/test) is a fair test, and it is what produced the
+12.20% figure reported above and used to reject the idea.

---

## Index: which fixed/rejected item lives in which file

| item | commit / stage file |
|---|---|
| A.1 clamp fix | `23be207`, `stages/106_inprocess.cpp` + `stages/capsule_decode.cpp` |
| A.2 orphaned uniques | `121fea9`, `stages/106_inprocess.cpp` |
| A.3 RC contained-read indexing | `121fea9`, `stages/capsule_decode.cpp` |
| A.4 FSE/HUF RLE | `f3ab0c2`, `include/coders_pgrc.h` |
| B.1 second-region self-match | `docs/SECOND_REGION_SELF_MATCH.md`, `SECOND_SELF` flag |
| B.3 MEM_MAXMM | `stages/106_inprocess.cpp`, `MEM_MAXMM_OVERRIDE` |
| B.14–B.16 names levers | `docs/REIMPL_NOTES.md`, section appended 2026-09-02 |
| B.17–B.18 quality levers | `docs/REIMPL_NOTES.md`, quality section |
