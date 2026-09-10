# A fresh caller architecture — the layers, designed from scratch

Written 2026-09-02, after reading DiscoSNP++'s source layer by layer
(`docs/DISCOSNP_INTERNALS.md`) and measuring where our current attempt loses
(`docs/CLAIM2_RESULTS.md`: ours 0.419 vs DiscoSNP++ 0.840 SNV F1 on identical
reads).

**This is a design document, not results.** Every layer is marked with what
is EVIDENCE-BACKED (measured this session), what is INFERRED (from reading
DiscoSNP++'s working implementation), and what is UNPROVEN (my design guess).
Do not cite anything here as a finding.

The framing: stop treating the caller as something bolted onto the
compression assembly. Design the calling substrate *for calling*, and share
with compression only what is genuinely shareable.

---

## What the measurements already forced us to conclude

Three hypotheses were formed and killed this session, and they constrain the
design:

1. Contigs too short → **refuted.** Ours are *longer* than DiscoSNP-era
   ARCS's (mean 335 vs 249) and still score half.
2. Anchor supply starved → **refuted.** We have 2.75× *more* exactly-2-
   occurrence 25-mers than ARCS (253,945 vs 92,163).
3. What survives: our substrate carries **3.8× the window in contig bases**
   (ARCS 2.7×, diploid ideal ~2×). The excess is *same-allele duplication* —
   one haplotype smeared across several contigs — which both fails to form
   bubbles and *consumes the 2-occurrence slot* the true hap1/hap2 pair needs.

So the defect is **substrate composition**, not contig length, not anchor
count, not filter thresholds (precision is 0.94 — what we call is right).

---

## The layer stack

### L0 — Ingest and k-mer index *(exists; shared with compression)*

Read loading, 2-bit packing, canonical k-mer counting. Already in the
encoder and in `caps_caller.h` (`pack31`/`canon31`). **This is the only layer
that should be shared with the compression path** — everything downstream
diverges.

### L1 — Coverage model *(exists, but under-used)*

Count histogram → error/heterozygous/homozygous coverage thresholds. We
compute `H` (haploid depth) already, by histogram valley-then-peak. Currently
it is used only in two filter comparisons.

**Fresh design: `H` should drive L2 and L3, not just filter at the end.**

### L2 — Error-aware k-mer filtering *(MISSING ENTIRELY — build this)*

DiscoSNP++ takes `-c 3`: k-mers below coverage 3 never enter the graph.
*(INFERRED from its interface and defaults, and consistent with its 0.971
precision.)*

We do **no** error filtering before building the calling substrate. Every
sequencing error that breaks an overlap creates a chain break and therefore
an extra contig fragment — a direct, mechanical contributor to the measured
3.8× redundancy. *(EVIDENCE-BACKED that the redundancy exists; UNPROVEN how
much of it is error-driven vs het-driven — measure this first, it is cheap:
histogram contig multiplicity against k-mer coverage.)*

### L3 — The calling substrate *(THE layer that decides the outcome)*

This is where 0.419 vs 0.840 is won or lost. Two candidate designs:

**(a) Haplotype-collapsed contig set** — keep our assembly, then collapse
contigs that are identical or near-identical (≤1–2 mismatches), keeping
*differing* pairs apart. Targets the measured 3.8×→~2× directly. Cheap to
build; reuses everything. *(UNPROVEN, but it is the minimal intervention
aimed squarely at the one surviving explanation.)*

**(b) A purpose-built compacted de Bruijn graph** — what DiscoSNP++ does
(GATB). Structurally the right substrate for bubbles, but it is a second
full data structure and effectively makes us a DiscoSNP++ reimplementation,
which is neither novel nor a fight we would obviously win.

**Recommendation: (a) first.** It is one experiment, ~90 s per window, and it
tests the surviving hypothesis directly. Fall back to (b) only if (a)
plateaus well below 0.840.

**The design constraint that keeps this OURS:** the collapsed substrate should
be *derivable from the archive*, so calling does not require re-assembly. That
is the actual novelty available to us — see "What makes this not a clone".

### L4 — Bubble enumeration *(exists, but structurally incomplete)*

Current state: `extract_snv_bubble` finds **isolated single-SNV bubbles
only**, and only where a canonical 25-mer occurs *exactly twice*.

Two known gaps, both visible in DiscoSNP++'s own output:

- **Multiple/close variants per bubble.** DiscoSNP++ emits IDs like `9_1`,
  `9_2`, `102_1` — several variants sharing one bubble. We cannot represent
  that at all; any het site within ~15 bp of another is invisible to us.
  *(EVIDENCE-BACKED: those IDs are in the VCF we scored.)* In a human sample
  at ~1 het/1000 bp this is not the dominant class, but it is not negligible.
- **`occ.size() != 2 → skip`.** A true haplotype pair whose k-mer also
  appears in a duplicate fragment has occ 3–4 and is silently discarded.
  Relaxing to 2–4 and testing all cross-contig pairs is a few lines.
  *(UNPROVEN gain, but it is free to test and directly implicated by the
  redundancy finding.)*

### L5 — Read-coherence and genotyping *(exists, but weak — this is the
precision layer)*

DiscoSNP++ dedicates an entire binary to this (`kissreads2`): it maps reads
back onto both bubble paths, keeps only "coherent" predictions, and derives
per-sample genotype/coverage. *(INFERRED to be the main reason its precision
is 0.971.)*

Ours is much thinner: window-mean coverage (`covwin`) and per-contig median
(`medcov`), plus an AF band. Notably we are already at precision 0.939 — so
this layer is **not currently our bottleneck** and should not be worked on
before L2/L3. It becomes the bottleneck the moment L3 raises recall, because
a richer substrate produces more candidate bubbles to reject.

**Fresh design opportunity:** we hold per-read *quality strings* (the archive
stores them), and mismatch positions are already computed during compression.
A quality-weighted coherence test is available to us essentially for free and
is not something DiscoSNP++'s coverage-based model does. *(UNPROVEN.)*

### L6 — Filters *(exists, frozen — do not touch)*

`HDMAX=2, MAF=0.20, DHI=2.5, KHI=1.1, MC=3, TRI=0.12, HALF=15`, ported
verbatim and validated held-out. Precision 0.939/0.948 across tuning and
held-out windows says these are working. **Any temptation to loosen these in
response to low recall is the wrong lever** — the variants are not reaching
the filters at all.

### L7 — Emission and lift *(exists, works)*

Contig-coordinate VCF, plus BWA-based lift for scoring only. Working and
verified. One lesson already banked: DiscoSNP++'s own `zero2one.py` shows how
easily an off-by-one hides here — ours is validated by the fact that lifted
calls score 0.94 precision, which a coordinate error would destroy.

---

## What makes this a new tool rather than a DiscoSNP++ clone

If we build L2–L5 as described, we have a competent bubble caller — but
DiscoSNP++ already is one, and it is fast, dependency-light and well
engineered. Matching it is not a contribution. The two things we can do that
it structurally cannot:

1. **The calling substrate is part of a lossless archive.** DiscoSNP++
   rebuilds its graph from reads every run and keeps nothing. If our
   collapsed substrate is stored in the archive, then re-calling, re-
   genotyping with different parameters, or calling a new variant class costs
   no re-assembly — and the archive is *still* the smallest lossless
   representation of the reads (Claim 1, 14/14). That is "one artifact, two
   products" in a form that survives the measurements, unlike the original
   "zero-cost byproduct" framing, which our own numbers have now weakened.
2. **Multi-sample calling across archives.** Cohort/joint calling is the
   documented route to the one class every single-sample reference-free
   caller misses — homozygous variants. Archives are addressable, so
   intersecting substrates across samples is tractable in a way that re-
   running a dBG per sample is not.

Neither is built. Both are honest, checkable claims rather than ratio
improvements, which is what this project's Claim 2 needs — because on pure
het-SNV F1, DiscoSNP++ at 0.840 is a genuinely strong incumbent and we should
not promise to beat it before the substrate work is done.

---

## Build order (each step measurable in ~90 s on one window)

1. **Measure** how much of the 3.8× redundancy is error-driven vs
   het-driven. Decides whether L2 alone is worth anything.
2. **L4 cheap fix:** relax `occ==2` to 2–4, test all pairs. Free, directly
   implicated.
3. **L3(a):** collapse near-identical contigs pre-calling. The main event.
4. Re-measure on r2 (tuning) and r3/na/r4/r5 (held out). If SNV F1 does not
   move materially past ~0.6, L3(b) — a real graph — is the honest next step
   and should be taken rather than tuning around.
5. **L5** quality-weighted coherence, only once recall has moved.
6. Only then: the archive-persisted substrate (novelty item 1).

Stop rule: if steps 2–4 leave us far below 0.840 and step L3(b) means
reimplementing a dBG caller, the right call may be to reframe Claim 2 around
novelty items 1–2 (persisted substrate, cohort calling) rather than claiming
a het-SNV ratio win. Decide that on the numbers, not in advance.
