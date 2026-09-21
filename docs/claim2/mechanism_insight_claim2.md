---
Date: 2026-09-19
Title: How Claim 2 Connects to the Paper's Central Insight — Mechanism,
  Not a Separate Achievement. Honest DiscoSNP++/Kmer2SNP Capability Check.
Purpose: The synthesis file for Claim 2. Explains why archive-native
  calling working at all is not a coincidence but the same mechanism
  Claim 3/T3.5 demonstrates for retrieval, and gives an honest,
  capability-checked (not guessed) account of why DiscoSNP++ and Kmer2SNP
  cannot do what T2.4/T2.5 measure.
When to refer to this file: Writing the discussion section connecting
  Claim 2 to Claim 3 and the paper's thesis; answering "why does calling
  from a compression-optimized archive work at all"; checking whether a
  DiscoSNP++/Kmer2SNP limitation claim is verified or assumed.
Keywords: insight, mechanism, reconciliation, collapse, allele splitting,
  build_substrate, DiscoSNP++, Kmer2SNP, bubble, multi-allelic, structural
---

# How Claim 2 connects to the paper's actual insight

## The insight, restated from the Claim 3 mechanism file, applied here

The paper's central mechanism finding: **a size-minimizing, reference-free
compressor is structurally pushed to place the two alleles of a
heterozygous site on different pseudogenomic segments, once divergence
passes some threshold** — folding both into one location costs a growing
mismatch list; splitting them costs two small, cheap placements. This is
not specific to retrieval (Claim 3) or to calling (Claim 2) — it is a
property of the *archive*, and both claims are separately measured
consequences of the same underlying fact.

**Why this matters for Claim 2 specifically**: a variant caller reading the
archive naively — treating one pseudogenomic position as one biological
site — would see only one allele at a split het locus, the same way a
naive single-coordinate query (T3.3's mechanism) only reaches one allele in
Claim 3. **The manuscript's own stated ablation makes this precise**: F1
0.431 with neither collapse nor re-placement, rising to 0.888 with both.
Without collapse specifically, precision holds at 0.96 while recall falls
to 0.27 — **the caller is blind to the second allele, not merely noisy
about it.** This is the same "the caller cannot see what it isn't told is
one site" pattern as T3.3-vs-T3.5 in Claim 3.

## The category-error correction behind why archive-native calling was worth trying at all

Verified against `docs/COMPRESSION_DERIVED_CALLING.md`. Before building a
caller from the archive, this project checked the literature on
assembly-derived variant calling and found a real concession worth citing:
Xu et al. 2017 (*Scientific Reports* 7:10963), simulating >3M SNVs
genome-wide, found "the assembly-based approach had a much lower recall
rate and precision compared to the alignment-based approach" at 30×
coverage — alignment-based recovering ~99% of imputed SNVs. Taken at face
value this reads as evidence *against* the whole project. **It is a
category error to apply it here, and finding that error is what justified
proceeding**: Xu et al.'s comparison is alignment-to-a-reference-genome
versus assemble-then-call — their 99% figure belongs to the arm that has a
reference. This project's caller has no reference at all; in their own
taxonomy it is assembly-based on *both* sides of what would be their
comparison, so their finding about reference alignment's advantage over
assembly does not transfer. (One further correction on the record: an
earlier draft of the source document quoted a "heterozygous SNVs have no
more than 50% chance to be called correctly" figure attributed to two URLs;
re-searching could not reproduce that sentence, and it is explicitly
withdrawn, not cited anywhere in this project's own claims.)

## A free candidate-generation channel already inside the archive, at ~1 byte per candidate

Verified against the same source. Independent of any caller logic, the
encoder's own MEM-extension mismatch streams (`mem_extmm_cnt/pos/obs`) are
`(position, reference base, observed base)` triples recorded so decode can
reproduce a read that was stored as a near-copy of an earlier pseudogenome
region — **a variant call set that already exists in the archive for an
unrelated reason**, the same "the archive already holds it" pattern as
Claim 3's native-pileup fast path. Measured on full HG002 chr20: 2,526,555
archive bytes yield 2,518,219 SNV candidate records — **about one byte per
candidate, 0.46% of the 547 MB archive** — and this is not an addable cost:
removing these streams makes the archive no longer lossless, so the calling
signal rides on bytes the format is already obligated to carry regardless
of whether a caller ever reads them. A second, related free signal exists
in the chaining loop itself: a read that shares a chaining seed but differs
by exactly one base over a long verified overlap is normally discarded
(`if(!rcmp(...)) continue;`) — and that discarded read is, structurally,
the *other* haplotype, already resident in registers at the point of
rejection.

## Why archive-native calling working at all is not surprising once this is understood

Before this project's own investigation, one might expect calling from a
compression-optimized structure to underperform calling from a
correctness-optimized alignment, since compression's objective (minimize
bits) and calling's objective (find true variants) seem unrelated. **The
actual finding inverts this**: because compression's own bit-minimization
pressure is what causes allele-splitting in the first place, *undoing* that
splitting (via `build_substrate`/`collapse_contigs` — see
`code_mapping_claim2.md`) recovers exactly the signal a reference-based
aligner would have had by construction, without needing a reference at all.
The mechanism that makes compression efficient is the same mechanism that,
once reversed, makes reference-free calling possible. This is the sentence
that belongs in the paper's discussion connecting Claim 1's efficiency to
Claim 2's accuracy — they are not independent achievements, one enables the
other's fix to be necessary and sufficient.

## T2.4 (multi-allelic) is the sharpest instance of this, not a bonus table

A multi-allelic site is the same phenomenon at higher order: not two
versions of a locus to reconcile, but three or more. The manuscript's own
text confirms this framing: "Multi-allelic sites test a different property
of the retained representation... whether several supported alleles can be
recovered as alternative forms of one locus rather than appearing only as
independent pairwise differences." **This is literally the same
reconciliation claim as T2.1's het-SNV ablation, generalized from 2 to N
alleles.** T2.5 (tetraploid) generalizes it again, along the ploidy axis
rather than the allele-count axis — same underlying mechanism, third
independent confirmation.

## Historical scale gap, and its resolution — verified against `docs/CLAIM2_FINAL_VERDICT.md` in full, not just its superseded-numbers banner

That document's own "Final verdict" section (written before the 2026-09-10
sweep) states plainly: *"Claim 2 should be described as validated on a
representative sample, not validated at the stated scale... the full 30x
chr20 per individual run has never been executed."* At that point every
Claim 2 number was measured on a ~75K-read chr20 *window*, not a full
individual. **This is resolved, not still open**: the same document's own
superseded-numbers table (top of file) labels the final figures explicitly
as *"full chr20, called from the archive"*, and this session's own T2.1-T2.5
CSV cross-check (`t21_snv_claim2.md` etc.) confirms full-chromosome-scale
counts (e.g. 37,011 TP on HG002, not a ~75K-window-sized count). State this
resolution explicitly if the paper's methods section needs to address
"was this validated at claimed scale" — the honest answer is yes, now,
where it was not at an earlier point in the project's history.

## A real, undissolved limit on how deeply this mechanism has been code-audited

The same document states: *"no algorithmic audit at the depth of Claim 3's
was performed on `caps_caller.h`"* — at 2,132 lines when written. **This
session found the file at 6,908 lines** — real growth since that audit
gap was recorded, meaning even the limited prior scrutiny does not cover
everything currently in the file. This directly reinforces, rather than
contradicts, this file's own earlier caveat: `build_substrate` and
`collapse_contigs` being the ablation's mechanism is a strong,
naming-and-commit-history-based inference, not a confirmed fact — and the
honest reason a deeper confirmation hasn't happened is that no audit at
that depth has ever been performed on this file, by anyone, at any point in
the project's history. This is worth stating plainly rather than
implying more rigor was applied than actually was.

## Two indel failure-mode findings, verified from a from-scratch re-scan

`docs/HET_INDEL_FRESH_SCAN.md` (2026-09-03) deliberately re-ran the het-indel
benchmark and inspected raw FP/FN VCF records by hand, without starting from
any prior analysis's conclusions — specifically to catch confirmation bias.

**Finding 1 — repeat context drives false positives, replicated via an
independent signal.** Prior docs identified homopolymer/short-tandem-repeat
context as the dominant false-positive driver using GIAB's own
`difficultregion=` truth annotations. This scan found the *same*
underlying signal through an unrelated path: the reference FASTA is
soft-masked (repeats in lowercase), and both `bcftools norm` and `rtg
vcfeval`'s reference lookup inherit that case when normalizing an indel —
so a call's case in the FP/FN report is a free, independent repeat marker.
Measured: false positives are lowercase (repeat-masked) 62.5% of the time
against true positives' 31.6% — **~2× enrichment**, replicating the
existing conclusion a second, unrelated way (small sample, n=16 FP/n=38 TP
on one window — suggestive, not definitive on its own, but consistent with
every prior measurement of the same phenomenon).

**Finding 2 — a distinct, previously undocumented failure mode: right
site, wrong polarity.** Cross-referencing FP and FN positions directly (not
done in prior analyses) found one exact positional match:
`chr20:3332481`, where GIAB truth calls a 4bp deletion (`TTTTA→T`) and this
project's caller called a 4bp *insertion* (`T→Tttta`) at the *identical*
position — not a missed candidate, not a locus mismatch, but a
length-change-direction error at what is almost certainly a poly-T
homopolymer, where two locally-plausible interpretations both fit the read
data and the wrong one was resolved. Occurred once in this window (1 of 16
FPs), reported as a genuinely distinct failure class and a concrete,
explicitly unresolved follow-up question — whether `extract_bubble`/
`scan_pair` (`caps_caller.h`) has a systematic insertion-vs-deletion bias
specifically at tandem-repeat loci — not as a conclusion this scan had the
sample size to reach.

## Verification discipline, evidenced concretely — worth citing for credibility, not just the results

The same document records two real measurement bugs caught and fixed
*before* being reported as caller results: a benchmark script's synthetic-
truth generator silently shifted indel positions (root-caused, fixed, not
a caller defect), and an apparent catastrophic runtime signal from
small-sample testing turned out to be a sampling artifact once the trend
across sample sizes was checked rather than trusting one data point. Both
are good, concrete evidence of the project's actual practice — checking an
invariant before trusting a number — rather than a claim made once and
left unexamined. **A third instance, verified against
`docs/PREFLIGHT_CHECKLIST.md`**: a pre-full-scale-run audit of every tunable
for accidentally-fitted constants found that the dBG bubble channel's
multi-polymorphism tolerance (`CAPS_DBG_MAXPOLY`) still defaulted to 3
(DiscoSNP++'s own `-P` default) in the *compiled binary*, while every
number this project had actually reported for that channel was measured
with `MAXPOLY=1` passed explicitly on the command line — a bare,
undocumented invocation would silently have used a configuration that
scores measurably worse (F1 0.861 at P=1 with coherence checking vs 0.850
at P=3, on the same window). Caught and fixed *before* the full-scale run,
not after a number was published from it — the default was changed to
match every validated configuration, and a bare run was re-verified to
reproduce the exact previously-reported TP/FP/FN afterward.

## Honest capability check: could DiscoSNP++ or Kmer2SNP be extended to do this?

**Kmer2SNP**: no indel model, no multi-allelic model, by construction — it
is a k-mer-graph SNV caller. This is stated as a structural limitation in
`docs/HET_INDEL_SOTA.md`'s survey and confirmed by its own absence from
T2.3/T2.4/T2.5 (`NOT_APPLICABLE`, not a measured 0). No claim is made here
about whether it *could* be extended — that would require reading its
source, which has not been done for this project (unlike the PgRC2 check
in Claim 3's mechanism file). **Do not claim more than this** — Kmer2SNP's
limitation is documented as observed capability, not verified architectural
impossibility, and the paper should reflect that distinction.

**DiscoSNP++**: its 0/26 on T2.4 is explained by output schema (bubble
records are inherently biallelic — one REF, one ALT per bubble), which is
a real, checkable architectural property of a bubble-based caller, not
guessed. **This has NOT been verified this session by reading DiscoSNP++'s
own source** — it is inferred from its documented output format and its
observed behavior (0 multi-allelic records across 3,989 total emitted, on
data where the truth set demonstrably contains 26 multi-allelic sites).
This is a reasonable, evidence-based inference, but it is not the same
rigor as the PgRC2 source-code check in Claim 3's mechanism file, and the
paper should not claim it as verified from DiscoSNP++'s own code unless
that check is actually done.

## The honest, precise paragraph for the paper's discussion section

*Claim 2's calling results are not independent of Claim 1's compression
mechanism — they are its direct consequence, reversed. The same
bit-minimization pressure that causes a compressor to place divergent
alleles at separate pseudogenomic addresses is what a naive caller would
read as two unrelated single-allele sites; reconciling those addresses
before calling (build_substrate/collapse_contigs) is what recovers the
signal, raising F1 from 0.431 to 0.888 in ablation. T2.4's multi-allelic
result and T2.5's tetraploid result are the same mechanism observed at
higher allele count and higher ploidy respectively, not separate
achievements. DiscoSNP++'s inability to represent multi-allelic sites
follows from its bubble-based, inherently biallelic output format; this is
inferred from its documented design and observed output, not verified by
reading its source, and is stated with that caveat.*
