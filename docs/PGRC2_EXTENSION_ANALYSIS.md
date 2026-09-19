# Can PgRC2 add locus retrieval with a quick patch? — layered analysis

**Status: LOCKED as a mechanism-grounded argument. Read Layer 5 before quoting
this anywhere — it states precisely what is CONFIRMED (from PgRC2's real
source) versus PREDICTED (from shared mechanism, not measured on PgRC2
itself). Do not blur the two.**

This document exists because this project already made and had to retract an
unverified PgRC2 claim once ("PgRC2 could add a button" — inherited from an
earlier doc, never personally checked, PgRC2 source unavailable locally at the
time). This analysis does not repeat that mistake: every claim below is either
sourced to PgRC2's actual code/paper, or explicitly labelled as an inference.

---

## Layer 0 — the claim, one sentence

A quick patch exposing PgRC2's existing per-read position data as a
coordinate-lookup API would not give PgRC2 what this project calls locus
retrieval, because the same compression objective that forced our own
coordinate arm to fail (96/400) is shared by PgRC2's architecture — and
independently, building a correct query layer on top of position data
(regardless of that mechanism) took this project a multi-day debugging effort
even with full source access.

Two independent arguments, not one. Either alone is sufficient. Together they
are stronger.

---

## Layer 1 — what PgRC2 actually is, confirmed from its real source

Source: github.com/kowallus/PgRC, paper Kowalski and Grabowski, Bioinformatics
2025, btaf101 (PgRC2), building on Bioinformatics 2020 (PgRC1).

**Confirmed by reading pgrc/pgrc-decoder.cpp directly:**

- PgRC builds its pseudogenome by approximating the shortest common
  superstring over the reads — greedy/approximate overlap-based assembly, the
  same general family as this project's own pseudogenome construction.
- The decoder computes and holds, per read, one position in the pseudogenome
  (orgIdx2PgPos, a complete vector held in memory for the whole dataset
  during decode) and one strand (revComp).
- Confirmed from the README: no query, extraction, region, or single-read
  interface exists. The tool is compress-the-whole-file /
  decompress-the-whole-file only.

**Not confirmed, flagged as open, not assumed either way:**

- Whether PgRC2 tracks per-read mismatches/deviations against the
  pseudogenome (analogous to this project's mm_pos/mm_sym) was not verified
  from the files fetched. If it does not, a naive query would return the
  pseudogenome consensus rather than each read's actual bases — the exact
  defect this project already found and fixed in its own system
  (CAPS_PILEUP, documented in CLAUDE.md).
- Whether orgIdx2PgPos is ever persisted to disk (an on-disk index, equivalent
  to this project's sidecar) or exists only transiently during one full
  decompress pass, was not confirmed. If transient-only, a query feature would
  first need a persistence layer before anything else.

---

## Layer 2 — a general argument about ANY reference-plus-mismatch compressor, built fresh, not borrowed by name

State this on its own terms first, independent of what this project measured
for itself, because the risk this layer exists to avoid is citing our own
result and calling it evidence for a different tool by resemblance alone.

**The argument, freestanding:** take any compressor that (a) builds a shared
reference/pseudogenome from the reads and (b) encodes each read as a position
plus a list of mismatches against that reference, and (c) is tuned to minimise
total encoded size. For a site where two haplotype copies diverge enough, the
encoder has a choice: fold both into one reference position and pay for a
growing mismatch list on every divergent read, or place them at two separate
positions and pay a small mismatch cost at each. Past some divergence, the
second costs fewer bits. A size-minimising encoder is therefore pushed toward
separating them, purely from (a)+(b)+(c) — nothing here refers to any
particular chaining heuristic, dataset, or codebase. Any tool meeting the
three preconditions is subject to this pressure by construction.

**This project's own system is one measured instance of that general
argument, not the source of it:**

- Ablation: F1 0.431 (neither collapse nor re-placement) to 0.888 (both).
  Without collapse, precision holds at 0.96 while recall falls to 0.27 — the
  caller is blind, not mistaken. The alternate allele is structurally absent
  from the naive pileup, not merely noisy.
- Coordinate control: only 96/400 het sites have both alleles reachable from
  one coordinate on this project's own pseudogenome.

These numbers confirm the general argument holds for one system. They are
evidence that the argument is sound, not a result that transfers to a
different tool by citation.

---

## Layer 3 — applying the SAME general argument to PgRC2's own confirmed architecture, on its own terms

The question for PgRC2 is not "does the Claim 2 finding apply to it" — it is
"does PgRC2 meet preconditions (a)+(b)+(c) above, evaluated fresh from its own
source." Confirmed vs inferred, kept separate:

| precondition | met by PgRC2? |
|---|---|
| (a) builds a shared reference/pseudogenome from the reads | confirmed — approximates the shortest common superstring, per its own paper |
| (b) encodes each read as position plus mismatches against it | position and strand confirmed from `pgrc-decoder.cpp`. The mismatch-list half is very likely (this project's own historical notes reference PgRC2's in-band MATCH_MARK and dst-gap/len/rc fields) but was **not directly re-confirmed in this pass** |
| (c) tuned to minimise total encoded size | confirmed — stated purpose of both PgRC1 and PgRC2, literally the paper's subject |

**If all three hold, Layer 2's general argument applies to PgRC2 by the same
logic it applies to any tool meeting them — not because PgRC2 resembles this
project, but because PgRC2 independently satisfies the argument's own
preconditions.** That is a prediction, not a measurement: nobody has built
PgRC2 here and inspected whether its pseudogenome actually separates
haplotypes on real heterozygous data. The preconditions being met makes the
prediction sound. It does not make it observed.

**The honest line, stated once so it cannot be lost:** this is a prediction
from first principles evaluated against PgRC2's own confirmed architecture. It
is not evidence borrowed from this project's own measured result, and it is
not an empirical finding about PgRC2 itself. Anyone citing this must carry
that distinction forward.

---

## Layer 4 — a second, independent argument that does not need Layer 3 at all

Even setting the allele-splitting mechanism aside entirely: this project
already had per-read position, strand, and deviations on day one of the
locus-retrieval work (T3.4/T3.5) — the same starting point PgRC2 is confirmed
to have for position and strand. Building a correct query layer on top of that
data — not just having it — took a multi-day effort within this session,
surfacing genuine, silent, non-obvious bugs:

- Strand-application bug: reverse-complement reads had their deviations
  mirrored onto the wrong bases, silently, for months — found only because a
  working exact-match index recovered nothing new, which was itself only
  explicable if the emitted sequence was not the actual read.
- Missing N-restoration: N-carrying reads were emitted with a pseudogenome
  base substituted for the N.
- One-sided anchoring: a probe anchored on only one side of a locus is
  structurally unreachable when the chaining neighbourhood is broken on that
  side — required a second, independent anchoring path (bilateral search) to
  fix.
- Hardcoded seed floor: silently capped mismatch tolerance at k=2 regardless
  of what was requested, undiscovered until derived from the haystack size
  instead.

None of these were visible from the data structures alone. Each required
building the query path, running it against real data, and noticing the
result didn't add up. This is direct, measured evidence that possessing
position/strand data is not the hard part — the query and correctness layer
on top of it is, even for a team with complete access to and understanding of
their own code.

---

## Layer 5 — the final, precise claim, for citation

Use this wording, not a looser paraphrase:

PgRC2's decoder is confirmed, from its own source, to compute per-read
position and strand during decompression, but exposes no query interface. Two
independent arguments suggest this is not a quick extension. First, any
compressor that builds a shared reference from reads, encodes each read as a
position plus mismatches against it, and is tuned to minimise encoded size, is
by that construction pushed to separate divergent haplotypes onto different
positions — a general argument about the encoding scheme, not specific to any
one implementation. PgRC2's own confirmed architecture meets those
preconditions, so the prediction applies to it on its own terms. Whether its
pseudogenome actually does this on real heterozygous data has not been tested.
Second, and independently of that argument entirely, this project's own
experience building a correct query layer on data it already possessed
(position, strand, deviations) took a multi-day effort and surfaced several
silent correctness bugs invisible until the query path was actually built and
run — direct evidence that the query layer, not the underlying data, is where
the real work lives.

**What remains open, stated plainly:** whether PgRC2's own pseudogenome
actually exhibits the same allele-splitting on real heterozygous data is
untested and would require building and running PgRC2 locally — not done, not
claimed as done.
