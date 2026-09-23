---
Date: 2026-09-20
Title: Architecture vs DiscoSNP++, Layer by Layer — Where We Match the
  Textbook, Where We Differ, and One Self-Correction
Purpose: `mechanism_insight_claim2.md` covers the allele-splitting mechanism
  and T2.1-T2.5 results; this file is the missing layer-by-layer
  architecture comparison against DiscoSNP++ specifically — an 18-row
  correspondence table plus the four real differentiators, condensed and
  verified against `docs/ARCHITECTURE_VS_DISCOSNP.md` (frozen 2026-09-10),
  not rewritten from memory.
When to refer to this file: Writing the Methods/Discussion comparison to
  DiscoSNP++ at implementation-level detail (not just result numbers);
  explaining precisely what is shared textbook algorithm vs. genuine
  contribution; citing the "assembly is the byproduct" self-correction.
Keywords: DiscoSNP++, kissnp2, GATB, architecture comparison, self-correction,
  retained assembly, kissreads2, COVCAP, HETSCAN, PLACEMENTS_AS_LINKS_REFUTED
---

# Architecture vs DiscoSNP++, layer by layer

## The shape of the comparison

DiscoSNP++'s caller (`kissnp2`, built on the GATB library) and this project's
Method B caller share a genuine amount of textbook bubble-calling machinery:
k-mer counting shape, graph membership testing, bubble start/closure
conditions, and branching policy are the *same algorithm*, arrived at
independently in several cases and verified by reading DiscoSNP++'s own
source. Of an 18-layer correspondence table, roughly half are marked
**textbook** (shared, not a contribution either way), a handful are
**theirs** (DiscoSNP++'s design, sometimes adopted, sometimes measured and
rejected), and the rest are **ours** — genuine differences, each with a
measured reason, not asserted superiority.

**Layers marked textbook (shared, same algorithm independently arrived at):**
k-mer counting shape, memory-budget declaration, graph membership oracle,
bubble start condition (branching node, both strands), SNV closure (lockstep
walk), branching policy (stop at first junction), low-complexity filtering
(DUST, off by default in both).

**Layers marked theirs (DiscoSNP++'s specific design choice):**
- Minimizer ordering (frequency-ranked) — tried on our side, **measured
  +14% volume, rejected** in favor of lexicographic.
- `-P 3` multi-polymorphism tolerance vs our `MAXPOLY=1` — narrower on
  purpose (structural: k=31 vs ~1/1000 het rate), and we still win on SNV
  F1 despite the narrower tolerance.
- Indel proposal: DiscoSNP++ attempts a BFS call at **every** branching node
  (`start_indel_prediction`, unconditional); ours fires only at superbubble
  nodes, **15.9% of nodes**. This is a real, honestly-acknowledged
  architectural deficit — not a place we claim to win, and it is the
  concrete reason the indel claim (T2.3) is reported as a real trade-off,
  not argued away.
- Indel ambiguity rejection (`checkRepeatSize`) — ported from their
  threshold, then measured to be a monotone loss on our data, so kept at
  their default rather than re-tuned (a case of adopting their parameter
  because re-deriving it did not pay).

**Layers marked ours (genuine differences, each measured):**
- **Read coherence.** DiscoSNP++ uses `kissreads2`, a *separate binary* that
  re-maps every read onto every candidate bubble. This project keeps reads
  in memory via one indexed sweep instead — a deliberate second FASTQ pass
  for this purpose specifically, **not** reads left resident by compression
  (the no-`CAPS_CALL` memory footprint is unaffected). This is what makes a
  *per-base* quality test affordable (233 MB bitmap) where DiscoSNP++ uses a
  per-path mean phred (2.27 GB of phred strings).
- **Depth-derived parameters, not fixed constants.** `COHC = max(2, H/10)`
  and `COVCAP = 2×PLOIDY×H` are functions of the sample's own measured
  haploid depth and ploidy, which is why precision holds 0.911–0.950 across
  a 10–30× coverage sweep and the tetraploid arm (T2.5) works without
  retuning. DiscoSNP++ exposes `-b`, `-P`, `-D`, `-max_ambigous_indel` as
  fixed numbers the user sets once.
- **`HETSCAN`** — a pair-fraction gate that declines to call on haploid
  input (measured 0.016 on E. coli), something DiscoSNP++ has no equivalent
  of.
- **Output**: a lossless archive from the same pass, not just VCF.

## The self-correction — "the ASSEMBLY is the byproduct, not the graph"

**An earlier version of this project's own documentation claimed the k-mer
table used by the caller was "already built by the compressor."** This was
checked against the actual code and found **false**: `kc_H_build` runs
inside `run_variant_call` (`include/caps_caller.h`), and
`src/encoder.cpp` (the encoder) builds no k-mer table at all — the
k-mer graph is paid for at call time, same as DiscoSNP++.

**What IS genuinely free is the assembly.** The encoder builds contigs in
order to compress, and the caller reuses them directly rather than
re-placing all reads from scratch. Measured: the alternative
(`build_substrate`, re-placing all 12.6M reads) takes **738 s serial at
full chr20**; skipping it via reuse is the difference between a ~57 s and a
~150 s caller run. This is the precise, measured version of the paper's
"retained assembly" thesis — real, but narrower than the withdrawn claim it
replaced.

## One tested-and-refuted headroom candidate, kept for the record

**"Free read threading" does not exist.** The hypothesis: the archive's
`ppos`/`read_cid` placements might substitute for the dedicated
read-to-path threading pass McCortex and LueVari pay for explicitly. Built
and measured twice, withdrawn — precision is *identical* across every split
tested (alleles sharing a contig vs. disjoint, short vs. long pseudogenome
span; 0.916–0.952 across all splits, no signal). **Root cause, structural**:
a pseudogenome is built by greedy overlap chaining, which *deliberately
merges* near-identical sequence — repeat copies are the most compressible
thing in a genome, so they collapse first. `ppos` records where a read was
**stored**, not where it **came from**. The repeat information was not left
unexploited on the table; it was **spent** to achieve Claim 1's compression
ratio, and recovering it would mean not collapsing repeats — directly
opposed to what makes Claim 1 win. (Note: the underlying idea of
read-coloured-graph calling was never claimed as novel — LueVari already
does reference-free SNP calling this way.)

## Why a de Bruijn graph cannot simply be built over our contigs

A tempting shortcut was tried and refuted (`docs/PG_AS_GRAPH_REFUTED.md`,
2026-09-04, on the full-chr20 HG002 run): transform the *existing*
pseudogenome into a de Bruijn graph rather than re-reading raw reads, since
DiscoSNP++'s bubble caller could then apply directly (a node with exactly 2
outgoing edges is a het site; 3+ is a repeat/paralog to reject). A first
probe seemed to refute this emphatically (97.6% false-positive rate at
degree ≥3) but was itself wrong — it measured *contig multiplicity*
(how many overlapping contigs tile the same locus), not true de Bruijn
*branch degree* (distinct successor bases, at most 4). Kept in the record
specifically because the wrong number was persuasive and pointed at the
same eventual conclusion — the exact circumstance in which a wrong
measurement gets believed.

**The correct measurement, on real GIAB-scored calls**: true positives sit
at a branching node only 3.7% of the time; false positives, 10.8% — a real,
2.9× enrichment signal, larger than most mechanisms tried against the
het-indel gap. But filtering on it makes F1 *worse* (0.592 → 0.580): the
true-positive pool is 7.6× larger than the false-positive pool, so a 2.9×
enrichment still removes far more true calls than false ones in absolute
terms.

**The finding that actually matters**: 97.7% of true-positive anchors and
92.0% of false-positive anchors have out-degree 1. **The contig set barely
branches at all** — it is a collection of near-disjoint linear paths, not a
branching graph, and this is not a tuning artifact. Greedy overlap chaining
*linearizes the data on purpose*, because that is what makes the
compression ratio good — and it is precisely what destroys the branch
structure a variant bubble needs to live in. DiscoSNP++'s graph branches
because it is built from raw reads, where both haplotypes remain separate
evidence at every position; this project's pseudogenome is built by
*choosing one path* through that evidence, and the redundancy removal that
gives Claim 1 its compression ratio is the same operation that removes the
branches a bubble caller needs. **This is a third, independent instance of
the paper's central tension** (the representation that compresses a signal
best is the one that conceals it — see `central_insight_and_mechanism.md`)
— here appearing as "graph branch structure" rather than "coordinate
addressability" or "heterozygous allele visibility," refuted by direct
measurement rather than assumed.

## A general representational limit, not specific to this project's caller

Verified against `docs/INDEL_BOUND.md`. Tracing every indel this project's
caller misses into its own k-mer graph found the branching node at each
missed site has real coverage (12–30 reads) but exactly **one** successor
base, with no coverage floor at play — raw successor counts return literal
`[0,0,21,0]`. The alternate haplotype's k-mers are genuinely *absent* from
the graph, not filtered out by a threshold. **66% of missed indels are
homopolymer-length changes**; most of the rest are tandem-repeat
expansions. Inserting one base into a run of that same base does not
create a second graph path at all — probing reads directly for a `G→GT`
call inside a poly-T run found the REF and ALT probe strings are *literally
identical* (`TCTGGGTTTTTGTTTTTCGGGTTTTTTTTTTT`), each present in 6 reads.
**No bubble exists for this event in a de Bruijn representation — for this
project's caller or for DiscoSNP++, or for any other de Bruijn-graph-based
caller**, since the limitation is in what a de Bruijn graph can represent
at a homopolymer, not in either implementation. This bounds what any
filter-based precision fix could achieve: at this project's measured
TP=3,290 on full chr20, beating DiscoSNP++'s F1 would require FP<366 (a
79% false-positive reduction with zero true-positive loss); even *perfect*
precision (FP=0, TP unchanged) yields only F1 0.595 — still short — because
the missing recall is a representational ceiling, not a threshold that
can be tuned away.

## Eight precision mechanisms tried against the indel gap — two helped, six did not, in a consistent pattern

Verified against `docs/INDEL_PRECISION_ROOT_CAUSE.md` (2026-09-02), written
after reading DiscoSNP++'s, eBWT2SNP's, and Kmer2SNP's precision mechanisms
directly from source/paper. Eight distinct precision-improving mechanisms
were implemented and swept, not merely proposed: (1) read-level junction
support via k-mer counts (kissreads2's principle) — **+0.18 precision, the
only large win**; (2) closing-anchor bubbles at both ends (DiscoSNP++) —
+0.03 precision, the only other mechanism that helped; (3) unique closing
anchor; (4) tolerant re-convergence flank; (5) full junction coherence;
(6) read-substring guarantee at 31/61/91bp (eBWT2SNP's principle); (7)
maximum-weight matching (Kmer2SNP's principle); (8) extended contig
agreement — all six were neutral or negative.

**The six failures group into three structural causes, not six unrelated
misses**: context-based tests (4, 6, 8) fail because this project's contigs
are short relative to fixed-context requirements (mean ~335 bp against
148 bp reads — a test needing 60–200 bp of flanking context runs out of
contig or rejects true events near contig ends); read-support tests (1
partly, 5, 6) are near-vacuous because every candidate here is built from
contigs which are built from reads, so it is read-supported *by
construction* — unlike DiscoSNP++, where a graph-traversal path can be
proposed that no read actually supports, giving `kissreads2` something real
to reject; the matching constraint (7) fails for a structural reason tied to
how this project's false bubbles specifically arise (repeat collapse, not
ambiguous pairing). **The honest conclusion, stated plainly in the source
document**: closing this gap is a substrate change of known shape (the
generator must read candidates out of the data rather than construct them),
not a tuning exercise — consistent with the layer-4/5 architectural
divergence already documented above, not a contradiction of it.

**A follow-up self-correction, 2026-09-02 (later the same day), kept on the
record**: the document's own §4 initially pointed at eBWT positional
clustering as the strongest fix candidate, citing its 99.13% precision
figure — the same figure `docs/EBWT_LITERATURE_CORRECTION.md` (2026-09-04)
later found to be simulated-data-only (real-data figure 66.62%, see
`novelty_and_prior_art.md`). Running each indel channel in isolation
(`CAPS_INDEL_ANCH=9999` suppresses the bubble channel, `CAPS_NO_PCLUSTER=1`
suppresses the positional-clustering channel) found the positional-
clustering channel alone is *weak* (F1 0.083–0.208 on two windows,
substantially below the bubble channel's 0.623–0.636 alone) — **partly
refuting the document's own earlier framing**. The channel was kept in the
shipped caller regardless, because combined with the bubble channel it adds
real, measured recall the bubble channel alone misses (0.655→0.699 combined
vs bubble-only, on the windows tested) — the honest reason it survives in
the final architecture is a measured combined-channel gain, not the
(subsequently corrected) precision figure that originally motivated trying
it.

## Net honest read

On the layers that are shared textbook algorithm (5–8 in the source table),
this project beats DiscoSNP++ on SNV F1 (0.879 vs 0.847 in the isolated
comparison) — the fair, apples-to-apples layer. On the one layer where
DiscoSNP++'s design is genuinely richer (unconditional per-node indel
proposal vs. our superbubble-gated 15.9%), that is stated as a real
architectural deficit, which is the honest reason the indel claim is framed
as a disclosed trade-off rather than an unqualified win.
