---
Date: 2026-09-20
Title: Technical Architecture — Container Format and the Assembly Pipeline
  (Rounds 1/2, Pigeonhole Mapping, MAXMAP Ramp, Second Region, MEM Self-Match)
Purpose: The one file that explains HOW the pseudogenome and archive are
  actually built, layer by layer, independent of any single claim's
  results. Condensed and verified against `docs/TECHNICAL_ARCHITECTURE.md`
  §1-2 (the project's own from-scratch architecture reference, written
  2026-09-02), not rewritten from memory.
When to refer to this file: Explaining the encoder's mechanism to a reader
  who has not seen the code; writing the Materials and Methods section;
  understanding WHY the allele-splitting mechanism in
  `central_insight_and_mechanism.md` happens (it falls directly out of
  round-1 chaining + MEM self-match, described here).
Keywords: architecture, container format, assembly, chaining, pseudogenome,
  MEM, maximal exact match, pigeonhole mapping, MAXMAP, second region,
  copMEM, cost-aware acceptance
---

# Technical architecture — container format and assembly

## Container format

File: `stages/106_inprocess.cpp`, struct `Archive`. On-disk layout: an 8-byte
magic, a version field, pseudogenome length + main-region boundary, and then
a sequence of self-identifying streams (name, length, payload). The decoder
(`stages/capsule_decode.cpp`) reads every stream into a name-keyed map and
looks each one up by name rather than by position — adding a stream is
additive and cannot shift any other stream's offset, a property an earlier
positional format lacked (it broke silently once; see
`docs/FAILURES_AND_REFUTED_IDEAS.md`). Header cost is ~300 B against archives
of 2.6–150 MB (0.001–0.01% of archive size).

## Assembly: building the pseudogenome, in order

**Round 1 — greedy exact suffix-prefix chaining.** Every read is packed 2
bits/base. A seed index maps every fixed-width k-mer to the reads containing
it at that offset. For each read, the algorithm looks for another read whose
prefix exactly matches this read's suffix for at least `MINOV` bases, and
chains them into one contiguous pseudogenome span instead of two separate
literal copies. This is single-pass and **exact** — no approximate matching
in the shipped path. **The single most important fix in the project's
history** (commit `3e06957`): the sweep used to start one base below the read
length, which makes an exact duplicate read structurally invisible to
chaining (two identical 251-base reads only overlap at L=251, not L=250).
Starting the sweep at the true max length let duplicates chain for free; on
S. acidocaldarius (9.4% duplication) this cut the pseudogenome from 9,134,100
to 6,157,270 bytes and flipped the project's one remaining compression loss
against PgRC2 into a win.

**Round 2 — division and second-pass chaining.** Reads that chained in round
1 are set aside; a second, looser sweep tries to attach the remaining reads
to the ends of already-formed chains.

**Pigeonhole mapping.** Reads that still have not chained are mapped onto the
pseudogenome built so far via a cheap read-relative placement: find an exact
k-mer seed hit inside the existing pseudogenome, verify against a bounded
mismatch tolerance (default 3), and record `(dst, src, len, rc)` — a
reference — instead of storing the read's bases again.

**The MAXMAP coverage ramp.** `MAXMAP` (the per-read candidate-mapping
ceiling) widens automatically as a function of `leftover_frac` — the fraction
of reads that failed round-1 chaining, known only after round 1 completes.
This is a measured-input-property lever, not a fitted constant: it is keyed
fresh on every run, with a floor set safely above every locked dataset's
measured leftover fraction, so normal-coverage behavior is provably
unchanged (byte-identical, verified). It engages only on inputs that
measurably need it — e.g. low-coverage data.

**The second region.** Reads that fail even pigeonhole mapping are not
stored raw. They are run through their own independent instance of the
round-1/round-2 chaining process, producing a **second pseudogenome region**,
appended after the main region's boundary. This turns unmatched leftovers
into more compressible sequence instead of literal padding.

**MEM self-match — the real workhorse of the size result.** Once both
regions exist, a maximal-exact-match pass (`copMEM`-style: hash k-mers,
verify to a full MEM) matches the whole pseudogenome against itself, looking
for repeated stretches at *any* position, not just chain-adjacent ones. Every
accepted MEM becomes a `(dst, src, len, rc)` reference exactly like a
pigeonhole placement, and the pg bytes it covers are removed from the literal
stream. Acceptance is cost-aware: a candidate MEM is accepted only if its
actual encoded cost is smaller than the literal cost ceiling (2.0 bits/base)
it would replace — a real algorithmic improvement over a blunt fixed-length
threshold, gated on where it was measured to help.

**Extension mismatch tolerance exists in the code but ships OFF.** A MEM
match can in principle extend past its first mismatch, trading a few extra
mismatch-stream bytes for a longer, cheaper-per-base reference. Measured and
rejected: it cost 490,763 B more than it saved on the dataset it was tested
on. Kept gated behind an override flag purely to reproduce that measurement,
per this project's standing rule that refuted ideas are recorded, not
deleted.

**Why the shipped encoder searches a 4-point candidate grid, not a wider
one — measured, not assumed.** `scripts/encode_adaptive.sh` tries several
`(MAXMAP, MINOV)` candidates per encode and keeps whichever produces the
smallest archive. Verified against the raw measurement data
(`docs/GRID_COST_MEASURED.txt`, `docs/GRID2_COST_MEASURED.txt`, 14 real
locked datasets, sizes 8.7 MB to 574 MB): widening the search from 4 points
to 8 costs **+0.000% archive size on 10 of 14 datasets**, with the 4 real
deviations tiny (SRR2584863 +0.012%, SRR40271341 +0.032%, SRR29296997
+0.028–0.472% depending on which grid variant, SRR39257532 the one real
outlier at +0.357%) — mean **+0.031%** across all 14. A single fixed
candidate point, by contrast, costs a mean of **+0.41%**, swinging as high
as +0.782% on individual datasets — an order of magnitude worse and highly
dataset-dependent, exactly as predicted by this project's own finding that
no fixed ratio of read length can express the optimal mapping ceiling
(§"MAXMAP coverage ramp" above). The 4-point grid is therefore the measured
sweet spot: it captures nearly all of the 8-point grid's benefit (mean
+0.031% vs +0.000%, i.e. leaves 0.031% on the table) while searching half
the space, and searching only one point would cost more than 13× as much
archive size for the time saved.

## Why this section is the mechanism, not just the implementation

`central_insight_and_mechanism.md` states the paper's central finding as a
consequence of "a size-minimizing, reference-free compressor" without
re-deriving the specific code path. This file is that code path: round-1
exact chaining is what forces two divergent alleles of a het site onto
separate contigs the moment their difference exceeds one seed's worth of
exactness, and MEM self-match is what then independently re-discovers and
collapses any allele that *does* still share enough sequence to be
worth referencing — the two mechanisms operating together are why the
paper's ablation (F1 0.431 → 0.888) and Claim 3's bilateral anchoring both
exist as corrections to the *same* underlying assembly decision, not two
unrelated fixes.
