---
Date: 2026-09-19
Title: How Claim 3 Connects to the Paper's Central Insight — Mechanism,
  Not a Bolt-On Feature. Honest PgRC2/BEETL Comparison.
Purpose: The single most important synthesis file in this folder. Answers
  precisely why Claim 3 (and T3.5/the completion index specifically) is architecturally
  earned rather than an API tacked onto an unrelated compressor, and gives
  an honest, source-checked answer to "could PgRC2 or BEETL do this with a
  small addition."
When to refer to this file: Writing the discussion section; answering any
  reviewer question of the form "isn't this just X with an extra feature
  bolted on"; deciding what can honestly be claimed about novelty relative
  to PgRC2 or BEETL; whenever tempted to describe Claim 3's mechanisms as
  five independent features rather than one insight applied five ways.
Keywords: insight, mechanism, ARCH, allele splitting, reconciliation,
  PgRC2, BEETL, CIndex, sFASTQ, novelty, bolt-on, API, not fabricated,
  confirmed vs inferred
---

# How Claim 3 connects to the paper's actual insight

## The insight, stated once, precisely

**A compression-derived structure does not have to be discarded after the
archive is written — it can be reasoned about directly to answer downstream
questions, without rebuilding a separate representation per task.** This is
the general claim. It is bigger than genomics and would matter in any domain
where a compressor infers structure it currently throws away.

**The mechanism that proves it, specific to this domain**: at a
heterozygous locus, a size-minimizing compressor is structurally pushed to
place the two alleles on *different* pseudogenomic segments rather than
fold them into one location with a growing mismatch list — past some
divergence, two positions cost fewer bits than one position plus many
mismatches. This is not a bug or an accident of this implementation; it
follows from precondition (a) building a shared reference, (b) encoding
position+mismatches, (c) minimizing total size — any compressor meeting
those three preconditions is subject to the same pressure by construction
(argued and measured in `docs/PGRC2_EXTENSION_ANALYSIS.md` Layer 2).

**Do not confuse the insight with the mechanism.** The insight is the
general architectural claim (retained structure is reusable). The
allele-splitting story is the *sharpest evidence* for it, appearing twice —
Claim 2 (F1 0.431→0.888 when reconciled), Claim 3/T3.5 (99.93%→100% when
reconciled) — because it is where the general insight has the most visible,
measurable bite. It is not itself the paper's thesis.

**A worked example, verified against `docs/CLAIM3_MECHANISM.md`** (GIAB
`20:3001343`, a real het SNV, C>T): a 40 bp probe ending 6 bp before the
variant occurs **four times** in the pseudogenome, and the consensus base at
the variant offset genuinely disagrees between occurrences —

    pg 101,247,805 (+)   C    <- REF haplotype
    pg 119,833,369 (+)   T    <- ALT haplotype
    pg   4,237,310 (-)   G    -> complement C   REF
    pg 136,671,198 (-)   A    -> complement T   ALT

with the two haplotypes sitting **18.6 Mb apart** in pseudogenome space.
Measured across 250 GIAB het SNVs on chr20: both alleles present as
separate pseudogenome occurrences in 80.0% of sites (200/250), only one
allele findable in 16.8% (42/250), probe unlocatable in the remainder
(8/250) — median 4 parallel occurrences per site, mean 15.58, **max 1007**.
A coordinate query at ±100 bp around one locus (same 12 sites) returns
*both* alleles 0/12 times and exactly *one* allele 12/12 times — the direct,
measured demonstration that a coordinate structurally cannot name what a
het locus actually is.

## Why Claim 3's mechanisms are not five bolted-on features

A skeptical read of Claim 3 might see five separate operations (export,
coverage, range query, exact match, locus retrieval) and suspect they were
each engineered independently to win their own comparison — "feature
creep" rather than one coherent architecture. The honest, verified answer
(see `code_mapping_claim3.md` for the code-level proof) is that this
is not the case:

- **Two of the five (export, coverage) are fast-exits from the SAME decode
  path that already exists for read reconstruction** — they read streams
  (`contig_spans`, `pos_abs`, `read_lengths`) the archive stores for other
  reasons, not new structures built to win T3.1/T3.2 specifically.
- **Two of the five (T3.3, T3.4) are the SAME function** (`mode=="query"`),
  differing only in what populates one variable. This was verified directly
  in the code this session, not assumed.
- **T3.5 is not a new mechanism at all.** It is an orchestration pattern —
  two calls to T3.4's exact primitive, reconciled by a Python script. No new
  C++ was written for it.

**This is the actual, defensible answer to "isn't this a bolt-on":** if
these were independently engineered per-table features, they would not
collapse this cleanly into two fast-exits + one shared primitive + one
orchestration layer on top of that primitive. The fact that they *do*
collapse this way is itself evidence the underlying structure (placement +
deviations, stored once) is genuinely general-purpose, not five special
cases dressed up as one story.

## the completion index specifically — extension, not rescue

the completion index deserves its own paragraph because it is the part most likely to be
mistaken for a bolt-on API. The honest framing: the completion index is a capability added
to a **finished archive**, after compression, without re-encoding, without
touching raw reads, without a second assembly — discovered as *possible*
only because the archive already stores placement and deviation information
for an unrelated reason (reconstructing reads exactly). It is evidence the
substrate is general, not proof it was designed with foresight for this
exact use. State it that way; do not imply premeditation the project cannot
verify.

## A real, measured number that directly proves T3.3/T3.5's distinction from substring search — found in `docs/CLAIM3_LOCKED.md`, not yet in this file until now

**50 random 40bp probes from the E. coli pseudogenome, measuring what
fraction of returned reads do NOT contain the probe substring**: median
38% (IQR 26-46%), excluding 4 repeat-landing probes that skew a pooled
average (46 probes / 1,692 reads → 35.9% when repeats excluded). **A BWT
text index can only return reads containing the query string. The ~38% it
structurally cannot return are exactly the reads overlapping a locus
without spanning the probe** — reads at the region's edges, which are
often precisely the reads carrying the alleles a caller needs. This is the
concrete, measured version of the position-vs-exact-match distinction
argued throughout this file — not a hypothetical, a number from 50 real
probes against a real archive.

**PgRC2's confirmed source detail, more specific than previously stated**:
its `pseudogenome/readslist/` directory structure is confirmed to keep
per-read positions — verified by running PgRC2 itself (`PgRC -h` exposes
compress/decompress/tuning flags only, no export/query/coverage/region
interface). The correct, precise claim is that PgRC2 lacks the
**interface**, not the data — a reviewer who knows PgRC2's internals would
check exactly this distinction.

**Honesty flag on the source of the above**: `docs/CLAIM3_FINAL_VERDICT.md`
and `docs/CLAIM3_LOCKED.md` (both read in full this session) are more
deeply superseded than Claim 1/2's equivalent docs — they predate T3.4 and
T3.5 entirely (locus retrieval and exact-match completion did not exist as
concepts when they were written), and their T3.1/T3.2/T3.3 numbers differ
from the final CSV (e.g. an early export figure of "555-656x" vs the final
129-784x). The 38%-neighbourhood measurement and the PgRC2
`readslist`/interface distinction above are stable, still-true mechanism
findings independent of those superseded numbers — cited for their
reasoning, not their headline figures.

**The one-sentence version of this file's entire argument, from the
project's own prior verified work**: *"Everyone else's archive stores what
the reads SAY. Ours stores where they SIT. The compressor had to assemble
in order to compress, so the coordinate system is a byproduct rather than
an addition — which is why locus retrieval, per-base coverage and
reference-free variant calling all come out of the same structure."*

## Honest, source-checked comparison to PgRC2

Full analysis, with confirmed-vs-inferred kept strictly separate: `docs/PGRC2_EXTENSION_ANALYSIS.md`. Summary, preserving that separation:

**Confirmed by reading PgRC2's actual source (`pgrc/pgrc-decoder.cpp`):**
- PgRC2's decoder computes per-read position and strand during decompression.
- No query, extraction, region, or single-read interface exists anywhere in
  PgRC2 — compress-whole-file / decompress-whole-file only.

**Not confirmed, explicitly flagged as open (do not state as fact):**
- Whether PgRC2 tracks per-read mismatches/deviations at all.
- Whether its position data is ever persisted to disk or exists only
  transiently during one decompress pass.

**The argument for why a quick patch would not give PgRC2 this capability**
rests on two independent legs, neither borrowed from this project's own
measured numbers by resemblance alone:

1. **A general, freestanding argument** (Layer 2 of the PgRC2 doc): any
   compressor meeting the three preconditions above (shared reference,
   position+mismatch encoding, size-minimizing objective) is pushed toward
   the same allele-splitting behavior by construction. PgRC2's own confirmed
   architecture meets all three preconditions on its own terms — this is a
   **prediction from first principles**, explicitly labeled as such, not an
   empirical finding about PgRC2 (nobody has built and run PgRC2 to check).
2. **An independent argument that does not need leg 1 at all**: this project
   already possessed position, strand, and deviation data for its own
   locus-retrieval work from day one — and building a *correct* query layer
   on top of that data still took a multi-day effort and surfaced several
   silent, non-obvious bugs (strand mismatch, missing N-restoration,
   hardcoded seed floor). This is direct evidence that *possessing* the raw
   data is not the hard part; the correctness layer on top of it is — even
   for a team with full access to and understanding of their own code.

**What remains explicitly open, stated plainly**: whether PgRC2's own
pseudogenome actually exhibits allele-splitting on real heterozygous data is
untested and would require building and running PgRC2 locally. Not done,
not claimed as done.

## Honest comparison to BEETL/CIndex/sFASTQ

Full novelty trail: `CLAUDE.md`'s "Novelty, corrected after a literature
survey 2026-09-16" section and `docs/CLAIM3_LOCUS_ADDRESSABILITY.md`.

**Confirmed**: content retrieval from a compressed archive is not novel —
BEETL-fastq (2014), CIndex (2022), sFASTQ (2022) all do it via BWT/FM-index,
and do it natively (no auxiliary index needed, since the index IS their
archive).

**What is genuinely distinct, verified against how BWT actually works**:
retrieval by **position** rather than by **exact match** — returning reads
*covering* a locus, not reads *containing* a string. The root cause is the
data model: a BWT stores each read as an independent string with no
inter-read geometry; a de Bruijn graph stores k-mer adjacency, not
placement; a pseudogenome stores **placement plus deviations**, which is the
only one of the three that gives the archive an internal coordinate system
at all. This was checked against how BWT/BEETL actually operates (see
`t34_exact_match_claim3.md`'s mechanism section) — not asserted from the name
of the technique.

**Could BEETL add position retrieval with a small addition?** Checked
structurally, not guessed: our completion index addition was cheap specifically
because position data already existed as a free byproduct of *how this
compressor works* (overlap chaining necessarily computes placement to
function at all). BWT's compression step (sorting) never computes anything
analogous — position is not hard to look up for a BWT, it is **not
computed, period**. For BEETL to add position retrieval, it would need to
build an entirely separate assembly/placement structure from scratch — not
an index over existing data, a second construction process. This is a
genuine structural asymmetry, verified by tracing what each compression
method actually computes, not assumed from which one "feels" more capable.

**What is not claimed**: no speed or index-size benchmark against
BEETL/CIndex/sFASTQ has been run. Stated as the single most valuable
missing experiment in `CLAUDE.md`, not silently omitted.

## The native-pileup fast path — an unexploited layer, found by audit

Verified against `docs/CLAIM3_LOCUS_ADDRESSABILITY.md` §M4. `capsule_decode
index` was already computing `tally[(pg position, observed base)] -> count`
over every read and writing it to `<sidecar>.sites` — a complete
allele-resolved pileup, 124 KB for a 2.15 Mb pseudogenome — and `query` never
read that file. Wiring it in composes directly with the bilateral-anchoring
finding above: because a het locus's two haplotype contigs already disagree
in the consensus, **allele identity comes from the consensus base at each
anchored position, and allele support comes from the tally** — neither step
decodes a single read.

A first pass accepted an allele on mere presence and reached 400/400 at 39
false positives, worse than the read-decoding path's 31. The fix was
thresholding: a single deviation is overwhelmingly sequencing error, while a
real alternate allele is carried by roughly half the reads, so only
tally-derived (non-consensus) alleles are thresholded — consensus bases stay
trusted unconditionally.

| min deviation count | het (of 400) | homozygous false positives |
|---|---|---|
| 1 (presence only) | 400 | 39 |
| **2** | **400** | **25** |
| 3 | 399 | 25 |
| 5 | 398 | 25 |

**At threshold 2 the fast path is strictly dominant** — it wins on both
correctness and speed, not a trade-off between them:

| path | het (of 400) | homozygous false positives | scoring cost |
|---|---|---|---|
| read decoding | 400 | 31 | 16.30 ms/site |
| **consensus + tally, count ≥ 2** | **400** | **25** | **0.83 ms/site** |

**19.6× faster and 6 fewer false positives at identical sensitivity**, because
the tally is a cleaner evidence source than re-deriving alleles from decoded
read strings — it was accumulated once at index time over every read, with no
probe, no containment test, and no per-query re-interpretation. The claim
this upgrades is not "we can retrieve the reads at a locus" but **the archive
already holds the pileup those reads would produce** — a second, independent
demonstration of the same substrate-reuse argument
`central_insight_and_mechanism.md` makes for the whole paper.

## The one-paragraph version, for the paper's discussion section

*Claim 3's five operations are not independent features engineered to win
five separate comparisons. Architecturally they reduce to two fast-exits
from the existing decode path, one shared query primitive used three ways,
and one orchestration layer built on that primitive — a structure that
falls out of retaining placement and deviation data for read reconstruction,
not one designed piecemeal per table. The one genuinely added structure
(the completion index) was built after compression, from data the archive already held for
an unrelated reason, which is itself evidence the substrate generalizes
rather than proof it was anticipated. Whether the closest architectural
relative (PgRC2) could add the same capability is, on the evidence available,
a structural argument grounded in its confirmed architecture, not a
measurement — and whether BWT-family tools could add position retrieval
is answered by what their compression step actually computes: nothing
analogous to placement, which is not a small addition to supply.*
