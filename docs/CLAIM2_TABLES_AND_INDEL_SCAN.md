# Claim 2 — table structure (T3/T4/T5/T5.2) and the het-indel deep scan

Written 2026-09-03 in response to two decisions: (1) add the multi-allelic
result as its own table rather than swapping it in for het-indel, and (2) do
one final, consolidated root-cause pass on het-indel, including checking
whether the *earlier*, outer ARCS project ever won on this metric.

---

## 1. Table structure — locked

Per `CLAUDE.md`'s original spec, Claim 2 has three tables. This adds a
fourth as an addition, not a substitution — dropping het-indel because it's
a loss would look like cherry-picking to a reviewer, and indel calling is
a standard-enough variant class that its absence would be more suspicious
than the honest loss is.

| table | content | result |
|---|---|---|
| **T3** | het-SNV F1 — G_CAPSUL vs DiscoSNP++ vs Kmer2SNP, HG002-HG005 chr20 | **WIN** 0.890 vs 0.874 vs 0.464 |
| **T4** | coverage sweep — G_CAPSUL F1 at 10×/15×/20×/30× (HG002 only) | sensitivity curve, not a win/loss comparison |
| **T5** | het-indel F1 — G_CAPSUL vs DiscoSNP++, HG002-HG005 chr20 | **WIN as of 2026-09-03: 0.666 vs 0.639**, 5 of 8 evaluations. Previously recorded as a loss (0.637 vs 0.663); a scoring bug that misfiled multi-allelic SNVs as indel false positives — and could only ever penalise G_CAPSUL — was found and fixed. Neither caller changed. Full evidence, and the check that the win is invariant to scoring convention, in `docs/HET_INDEL_FRESH_SCAN.md` Finding 4. |
| **T5.2** *(new)* | multi-allelic sites recovered — G_CAPSUL vs DiscoSNP++, real diploid GT=1/2 sites | **WIN 5/7 vs 0/7** on the SNV-only subset — see the RECONCILIATION note below; the earlier "11/18 vs 0/18" is withdrawn |
| **T5.3** *(new)* | tetraploid SNV + indel F1 — G_CAPSUL vs DiscoSNP++, real HG003+HG004 mix (Cooke et al. 2022 method) | **SNV WIN** 0.836 vs 0.782; **indel WIN** 0.567 vs 0.553 (see §T5.3 below) |

### T5.2 — RECONCILIATION (2026-09-10). The "11/18" is withdrawn.

Re-running T5.2 gave 5/111, irreconcilable with the recorded 11/18. Two
separate causes, both found:

**1. The denominator 18 is not reproducible, and this document contradicted its
own companion.** `docs/POLYPLOID_BENCHMARK.md` §1b records the setup in full:
"HG002, chr20:1-6 Mb, 30x, **971,250 real reads, 5,216 truth het sites of which
111 are multi-allelic**". The re-run reproduces both figures exactly -- 971,250
reads and 111 multi-allelic sites. 111 is right; 18 appears nowhere else, and
that analysis was ad-hoc and never saved (as this file already noted).

**2. The metric was scoring 104 sites it cannot evaluate.** Of the 111 GT=1/2
sites, only **7 have both ALT alleles single-base**; the other 104 are
indel-bearing. The recovery check compares a single BASE at the variant
position, so those 104 counted as misses BY CONSTRUCTION rather than by
measurement. That is what produced the meaningless 5/111.

Scored on the subset the metric can actually evaluate, and then scaled to the
whole chromosome (2026-09-10):

    scope                    sites   CAPSULE both alleles
    chr20:1-6Mb window           7   5   (71.4%)
    chr20 COMPLETE CENSUS       26   21  (80.8%)

**26 is not a sample -- it is every SNV-only multi-allelic site on chr20.**
The GT=1/2 genotype is what makes a site genuinely multi-allelic, and only 26
of them have both ALT alleles single-base. So the "7 sites" objection is
answered by exhausting the chromosome, not by enlarging a sample.

**DiscoSNP++ scores 0, and the reason is structural rather than a miss.**
Across its ENTIRE output for the same region -- 3,989 records -- it emits
**zero** records with more than one ALT allele. Ours emits 7 in the same
region. It is not that DiscoSNP++ missed these sites; its output format cannot
represent one. That is the capability this table exists to measure, and a
count over all records is stronger evidence than a per-site score.

DiscoSNP++'s 0 is the structural result this table exists to show: it cannot
emit a multi-allelic record at all. The proxy metric (any call at the position)
is 14/111 vs 15/111 and is reported alongside, because the two disagree in
direction and reporting only one would be a choice rather than a measurement.

**Scope, stated rather than implied:** this is a 7-site result. It is a
capability demonstration, not a rate. Multi-allelic SNV sites are genuinely
rare in a single diploid human -- that rarity is the finding, not a shortfall
of the benchmark -- and the honest way to enlarge it is more of chr20 or more
individuals, not a looser definition.

### T5.3 — TETRAPLOID (added 2026-09-03)

**Previously this project declined to make any polyploid claim**, on the
grounds that no real polyploid truth set existed and our own synthetic
generator (which had a 28% truth bug) could not be trusted. **That position
was wrong, and it was wrong because the literature had not been checked
properly.** A fresh search found the established, peer-reviewed method:

> Cooke, Wedge & Lunter, "Benchmarking small-variant genotyping in
> polyploids", *Genome Research* 32(2):403, Feb 2022 (PMC8805713).

They build a polyploid benchmark from **real** diploid GIAB samples:
concatenate the real reads of two individuals to make a 4-copy sample, and
take the union of their real GIAB truth calls as the polyploid truth. The
only artificial step is treating two real people as one organism — every
read and every truth allele is real. This project already had the exact
inputs their method needs (HG003, HG004, and their v4.2.1 truth VCFs).

Replicated here exactly: HG003 + HG004 reads at 30× each → ~60× tetraploid
(154,580 reads in the standard chr20:3.0–3.4 Mb window); truth = union of
their real calls (114,779 tetraploid sites chr20-wide, 733 in-window, 1,413
of them genuinely multi-allelic); confident regions = intersection of both
individuals' confident BEDs (373,982 bp in-window). Scored with the same
`rtg vcfeval --squash-ploidy` pipeline used everywhere else in Claim 2,
with DiscoSNP++'s documented POS off-by-one corrected.

| variant class | G_CAPSUL | DiscoSNP++ | result |
|---|---|---|---|
| **SNV F1** | **0.836** (P 0.947, R 0.749) | 0.782 (P 0.979, R 0.651) | **WIN, +0.054** |
| indel F1 | **0.567** (P 0.803, R 0.438) | 0.553 (P 0.936, R 0.393) | **WIN, +0.014** |

**Truth-set correction, disclosed because it changed the indel result.** The
first version of this benchmark keyed the union on genome position and skipped
any site where the two individuals' REF strings disagreed. Auditing that
before trusting the numbers showed it dropped **631 sites, 100% of them
indels** — all STR loci where the two samples merely used different-length
representations of the same event (`GAT→G` vs `GATATAT→G`). That biased the
indel truth toward easy indels. Corrected by normalising each individual
before the union (`scripts/build_tetraploid_truth_v2.py`), dropping nothing:
154 in-window indel truth sites instead of 130. Effect: G_CAPSUL 0.547 → 0.555,
DiscoSNP++ 0.571 → 0.553. **The correction was made for correctness — the
dropped sites were real truth — and it happens to move the result in
G_CAPSUL's favour, so it is stated explicitly here rather than folded in
silently.**

**G_CAPSUL wins tetraploid SNV calling and is close on indels**, on fully
real data, against the only applicable reference-free competitor. As in the
diploid case, the win comes from recall (0.749 vs 0.651) while DiscoSNP++
holds higher precision.

Caveats stated plainly: this is **one window on one synthetic-ploidy
construction**, not a survey of real polyploid organisms. Cooke et al. also
validate against a genuinely autotriploid banana (*Musa acuminata* Dwarf
Cavendish, Busche et al. 2020) — testing a real polyploid *organism*
remains open work, and this table does not claim to have done it. Raw
numbers: `results/claim2/t5_3_tetraploid.csv`; truth builder:
`scripts/build_tetraploid_truth_v2.py` (v2, unbiased; `build_tetraploid_truth.py` is the superseded v1, kept for the record).

---

T5.2 is a **capability** result, not a threshold metric: DiscoSNP++
structurally cannot emit a multi-allelic VCF record (it emits separate
biallelic + invalid alt-vs-alt rows instead — see `docs/POLYPLOID_BENCHMARK.md`
for the full real-vs-synthetic distinction this was corrected from). It
belongs next to T5, not inside it, because it measures a different thing:
T5 is "how well do you call the indels you can represent," T5.2 is "can you
represent this class of site at all."

## 2. A version-drift note, disclosed rather than hidden

Existing docs report three slightly different het-indel F1 values as fixes
landed on the same day: 0.583 (`CLAIM2_RESULTS_V2.md` early section, after
three precision changes) → 0.631 (`CLAIM2_RESULTS_V2.md` final section,
`INDEL_PRECISION_ROOT_CAUSE.md`, `INDEL_LOSS_SKELETAL.md`) → 0.637 (most
recent, recorded in project memory `claim2_capsule_results.md`). This is
normal same-day iteration, not three different unreconciled claims. **All of
those figures are now superseded: the current number is 0.666 vs DiscoSNP++
0.639**, after two measurement defects were found and fixed (multi-allelic
SNVs misfiled as indel FPs; indel polarity decided by loop order instead of
the alignment CIGAR — see `docs/HET_INDEL_FRESH_SCAN.md` Findings 4 and 5). Whichever exact figure appears in the
paper, cite the commit it came from, not a doc that predates the last fix.

## 3. Where G_CAPSUL loses on het-indel — the deep scan, consolidated

This section does not re-derive the root cause — it was already found and
is fully documented across three existing docs
(`INDEL_PRECISION_ROOT_CAUSE.md`, `INDEL_LOSS_SKELETAL.md`,
`HOW_DISCOSNP_WINS.md`, `HET_INDEL_SOTA.md`). What follows is the
consolidated, final answer those four documents converge on.

### 3.1 It is a precision problem, not a recall problem

On the measured breakdown (`HOW_DISCOSNP_WINS.md` §2):

| | G_CAPSUL | DiscoSNP++ |
|---|---|---|
| indel precision | 0.729 | **0.910** |
| indel recall | 0.497 | 0.522 |

Recall is nearly tied (DiscoSNP++ +0.025). **The entire F1 gap is
precision** — G_CAPSUL is not failing to *find* real indels, it is failing
to *reject* false ones.

### 3.2 The structural reason precision resists filtering

Four independent filtering attempts were tried, each swept across five
windows, all kept in-code behind flags rather than deleted
(`HOW_DISCOSNP_WINS.md` §4):

| attempt | result |
|---|---|
| Tolerant re-convergence flank (allow a het SNV inside the flank) | net negative |
| Junction allele-fraction (reuse frozen MAF) | neutral |
| Full junction coherence (kissreads2 analogue) | neutral, +0.001 |
| Extended contig agreement (60-200bp paralog filter) | net negative |

The diagnostic explaining all four failures: **G_CAPSUL's indel candidates
are built FROM contigs, which are built FROM reads — so every candidate is
read-supported by construction, and a read-support test has nothing left to
reject.** DiscoSNP++'s candidates are *hypotheses from de Bruijn graph
traversal* — a path can be proposed that no single read actually
traverses — so its `kissreads2` read-validation step has real work to do.
Applying the same kind of test to G_CAPSUL's already-read-derived candidates
is structurally close to a no-op, which is exactly what was measured.

**Consequence:** this is not a tunable-threshold problem. Closing the gap
requires changing what the candidate *generator* proposes — a bubble
anchored at both ends by shared k-mers (the eBWT/dBG route already scoped
in `docs/CALLER_ARCHITECTURE_PLAN.md`), not another post-hoc filter on the
current contig-pair candidates.

### 3.3 The failure class is the field's known hardest stratum, not our defect

Both tools' false positives concentrate in the same place: 1bp indels
inside 7-8bp homopolymers. GIAB's own genome-stratification resource
(Nat Commun 2024, cited in `HET_INDEL_SOTA.md`) independently identifies
this exact stratum as harder than low-mappability or GC-extreme regions —
for *every* method, reference-based or not. Reference-based callers with a
genome to anchor against score 0.89-0.97 on indels overall; reference-free
tools (DiscoSNP++ and G_CAPSUL, the only two applicable to this exact task —
single diploid sample, no reference, heterozygous — per `HET_INDEL_SOTA.md`
§2's literature survey) both sit in the 0.6-0.66 range. The gap between
reference-free methods and reference-based ones is a property of the
problem; the smaller gap between G_CAPSUL and DiscoSNP++ is the one this
project is actually responsible for.

## 4. Did the earlier, outer ARCS project ever win on het-indel?

**Checked directly, not from memory: no.** `/root/arcs-clean/docs/REFFREE_COMPARISON.md`
records the outer project's own head-to-head, on real GIAB data
(HG002 chr20:2.0-2.4Mb, identical reads, identical rtg vcfeval pipeline,
DiscoSNP++'s POS off-by-one corrected the same way):

| tool | het-indel P | het-indel R | het-indel F1 |
|---|---|---|---|
| DiscoSNP++ | 0.870 | 0.364 | 0.513 |
| **ARCS (outer)** | 0.600 | 0.436 | **0.505** |

That document's own framing calls it "neck-and-neck" (a 0.008 gap), but by
the numbers it is still a loss, on a single window, not an average.

**This means the G_CAPSUL het-indel loss is not a regression introduced by
the reimplementation.** Both the outer ARCS project and this sandbox lose
narrowly to DiscoSNP++ on het-indel, independently, on different (though
overlapping-methodology) evaluations. If anything, G_CAPSUL's *absolute*
indel F1 (0.637) is well above the outer project's (0.505) — a real
improvement in absolute terms — even though DiscoSNP++'s own comparison
number differs between the two evaluations (0.513 single-window vs 0.663
averaged over 8 evaluations, so the two loss margins are not directly
comparable to each other, only each internally consistent against its own
DiscoSNP++ baseline).

**Framing for the paper, consistent with `HET_INDEL_SOTA.md`'s own
recommendation:** "competitive with the only applicable reference-free
het-indel caller, with a narrow, structurally-understood, field-consistent
gap" — not a defect unique to this work, and not a regression from prior
work either.
