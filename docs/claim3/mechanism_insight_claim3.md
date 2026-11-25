# Claim 3 mechanism: why the operations follow from the archive structure

This file explains why the operations of Claim 3 (export, coverage, range query, exact-match and
locus retrieval) follow from the structure that the archive already keeps, where the FM-index fits,
and whether PgRC2 or the BWT-family tools could reproduce the result with a small addition. It
relates to the overview in [`overview_claim3.md`](overview_claim3.md). "The development notes" are
the author's internal working record, which is not distributed with the repository.

## The insight

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
(argued and measured in the development notes, in the PgRC2 extension analysis).

The insight and the mechanism are different things. The insight is the general architectural
claim that retained structure is reusable. Allele splitting is the clearest evidence for it, and it
appears twice: in Claim 2 (F1 0.431 to 0.888 when reconciled) and in Claim 3, T3.5 (99.93% to 100%
when reconciled). It is where the insight has the most measurable effect, and it is not itself the
thesis of the paper.

**A worked example, verified against the Claim 3 mechanism note** (GIAB
`20:3001343`, a real het SNV, C>T): a 40 bp probe ending 6 bp before the
variant occurs **four times** in the pseudogenome, and the consensus base at
the variant offset disagrees between occurrences —

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
creep" rather than one coherent architecture. The verified answer
(see `code_mapping_claim3.md` for the code-level proof) is that this
is not the case:

- **Two of the five (export, coverage) are fast-exits from the SAME decode
  path that already exists for read reconstruction** — they read streams
  (`contig_spans`, `pos_abs`, `read_lengths`) the archive stores for other
  reasons, not new structures built to win T3.1/T3.2 specifically.
- **Two of the five (T3.3, T3.4) are the SAME function** (`mode=="query"`),
  differing only in what populates one variable. This was verified directly
  in the code.
- **T3.5 is not a new mechanism at all.** It is an orchestration pattern —
  two calls to T3.4's exact primitive, reconciled by a Python script. No new
  C++ was written for it.

**This answers the objection that the operations are bolt-ons.** If
these were independently engineered per-table features, they would not
collapse this cleanly into two fast-exits + one shared primitive + one
orchestration layer on top of that primitive. The fact that they *do*
collapse this way is itself evidence the underlying structure (placement +
deviations, stored once) is general-purpose, not five special
cases dressed up as one story.

## The FM-index specifically — extension, not rescue

The FM-index is the part most likely to be mistaken for a bolt-on API. It is a capability added
to a **finished archive**, after compression, without re-encoding, without
touching raw reads, without a second assembly — possible
only because the archive already stores placement and deviation information
for an unrelated reason (reconstructing reads exactly). It also indexes only
the reads that deviate from the pseudogenome, because the pseudogenome
already holds the others. It is evidence that the
substrate is general. It is not evidence that the substrate was designed with this use in
mind.

## Reads that a substring search cannot return

Take 50 random 40 bp probes from the E. coli pseudogenome and measure what fraction of the returned
reads do not contain the probe as a substring. The median is 38% (IQR 26 to 46%), excluding 4
probes that land in repeats and skew a pooled average (46 probes and 1,692 reads give 35.9% with
the repeats excluded). A BWT text index can return only reads that contain the query string, so
it cannot return these reads. They are the reads that overlap a locus without spanning the probe,
at the edges of the region, and they are often the reads that carry the alleles a caller needs.
This is a measured version of the distinction between position and exact match.

**PgRC2 lacks the interface and not the data.** Its `pseudogenome/readslist/` directory keeps
per-read positions. Running PgRC2 itself showed that `PgRC -h` exposes only compress, decompress and
tuning flags, with no export, query, coverage or region interface.

The two earlier Claim 3 notes in the development notes (the final-verdict note and the first
results note) predate T3.4 and T3.5, and their T3.1 to T3.3 numbers differ from the final CSV (for
example an early export figure of 555 to 656 times against the final 129 to 784 times). The 38%
measurement and the PgRC2 `readslist` observation are independent of those numbers and remain
valid.

The argument of this file in one sentence: other archives store what the reads say, and this one
stores where they sit. The compressor has to assemble in order to compress, so the coordinate
system is a byproduct, and locus retrieval, per-base coverage and reference-free variant calling
all come from the same structure.

## Source-checked comparison to PgRC2

The full analysis, with confirmed and inferred statements kept separate, is in the PgRC2 extension analysis of the development notes. In summary:

**Confirmed by reading PgRC2's actual decoder source:**
- PgRC2's decoder computes per-read position and strand during decompression.
- No query, extraction, region, or single-read interface exists anywhere in
  PgRC2 — compress-whole-file / decompress-whole-file only.

**Not confirmed (open):**
- Whether PgRC2 tracks per-read mismatches/deviations at all.
- Whether its position data is ever persisted to disk or exists only
  transiently during one decompress pass.

**Why a quick patch would not give PgRC2 this capability** rests on two independent arguments:

1. **A general, freestanding argument** (Layer 2 of the PgRC2 doc): any
   compressor meeting the three preconditions above (shared reference,
   position+mismatch encoding, size-minimizing objective) is pushed toward
   the same allele-splitting behavior by construction. PgRC2's own confirmed
   architecture meets all three preconditions on its own terms — this is a
   **prediction from first principles**, which is a prediction from first principles and not an
   empirical finding about PgRC2, because the allele-splitting behavior was not checked by running PgRC2.
2. **An independent argument that does not need leg 1 at all**: this project
   already possessed position, strand, and deviation data for its own
   locus-retrieval work from day one — and building a *correct* query layer
   on top of that data still took a multi-day effort and surfaced several
   silent, non-obvious bugs (strand mismatch, missing N-restoration,
   hardcoded seed floor). This is direct evidence that *possessing* the raw
   data is not the hard part; the correctness layer on top of it is — even
   for a team with full access to and understanding of their own code.

**Open:** whether PgRC2's own
pseudogenome actually exhibits allele-splitting on real heterozygous data is
untested, and testing it would require building and running PgRC2 locally.

## Comparison to BEETL/CIndex/sFASTQ

The literature survey behind this comparison is in the development notes (the section on novelty corrected after the survey of September 2026) and in the Claim 3 locus-addressability write-up. [`../novelty_and_prior_art.md`](../novelty_and_prior_art.md) summarizes the related work.

**Confirmed**: content retrieval from a compressed archive is not novel —
BEETL-fastq (2014), CIndex (2022), sFASTQ (2022) all do it via BWT/FM-index,
and do it natively (no auxiliary index needed, since the index IS their
archive).

**What is distinct:**
retrieval by **position** rather than by **exact match** — returning reads
*covering* a locus, not reads *containing* a string. The root cause is the
data model: a BWT stores each read as an independent string with no
inter-read geometry; a de Bruijn graph stores k-mer adjacency, not
placement; a pseudogenome stores **placement plus deviations**, which is the
only one of the three that gives the archive an internal coordinate system
at all. This was checked against how BWT and BEETL operate (see the mechanism section of
[`t34_exact_match_claim3.md`](t34_exact_match_claim3.md)).

**Could BEETL add position retrieval with a small addition?** Not easily, for a structural reason. G_CAPSUL's FM-index addition was cheap specifically
because position data already existed as a free byproduct of *how this
compressor works* (overlap chaining necessarily computes placement to
function at all). BWT's compression step (sorting) never computes anything
analogous — position is not hard to look up for a BWT, it is **not
computed, period**. For BEETL to add position retrieval, it would need to
build an entirely separate assembly/placement structure from scratch — not
an index over existing data, a second construction process. This is a
structural asymmetry that follows from what each compression method computes.

**What is not claimed**: no speed or index-size benchmark against
BEETL/CIndex/sFASTQ has been run. It is recorded as the most valuable
missing experiment in the development notes.

## The native pileup fast path

An audit found an unused layer. `capsule_decode
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

**At threshold 2 the fast path is better on both correctness and speed**, with no trade-off between them:

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

## Summary

*Claim 3's five operations are not independent features engineered to win
five separate comparisons. Architecturally they reduce to two fast-exits
from the existing decode path, one shared query primitive used three ways,
and one orchestration layer built on that primitive — a structure that
falls out of retaining placement and deviation data for read reconstruction,
not one designed piecemeal per table. The one added structure
(the FM-index) was built after compression, from data the archive already held for
an unrelated reason, which is itself evidence the substrate generalizes
rather than proof it was anticipated. Whether the closest architectural
relative (PgRC2) could add the same capability is, on the evidence available,
a structural argument grounded in its confirmed architecture, not a
measurement — and whether BWT-family tools could add position retrieval
is answered by what their compression step actually computes: nothing
analogous to placement, which is not a small addition to supply.*
