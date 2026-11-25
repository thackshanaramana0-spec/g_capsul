# How Claim 2 connects to the paper's actual insight

> Mechanism of Claim 2: how archive-native calling relates to the mechanism behind Claim 3 retrieval, with the capabilities of comparator tools checked against their documentation.

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

Before building a
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

Independent of any caller logic, the
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

## T2.4 (multi-allelic) is the clearest instance of this, not a bonus table

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

## Evaluation scale

An earlier version of the Claim 2 record described the results as validated on a representative
sample and not at the stated scale, because every number then came from a window of about 75,000
reads of chromosome 20 and not from a full individual. That is resolved. The final figures are from
full chromosome 20, called from the archive, and the T2.1 to T2.5 CSV files give counts at that
scale (for example 37,011 true positives on HG002, which a 75,000-read window could not produce).

## How closely the mechanism has been code-audited

An earlier review noted that `caps_caller.h`, then 2,132 lines, had not had an algorithmic audit
at the depth applied to the Claim 3 code. The file now has 5,970 lines, so even that limited
scrutiny does not cover all of it. The statement that `build_substrate` and `collapse_contigs`
are the mechanism behind the ablation therefore rests on their names, their commit history and the
ablation itself, and not on a line-by-line audit.

## Two indel failure modes found by a fresh re-scan

The heterozygous-indel benchmark was re-run and the raw false-positive and false-negative VCF
records were inspected by hand, without starting from earlier conclusions, to avoid confirmation
bias.

**Finding 1: repeat context drives false positives.** Earlier analysis had identified
homopolymer and short-tandem-repeat context as the dominant false-positive driver, using GIAB's
`difficultregion=` annotations. The re-scan found the same signal by an unrelated route. The
reference FASTA is soft-masked (repeats in lowercase), and both `bcftools norm` and the reference
lookup of `rtg vcfeval` keep that case when normalizing an indel, so the case of a call in the
report marks a repeat. False positives were lowercase 62.5% of the time against 31.6% for true
positives, about a twofold enrichment. The sample is small (16 false positives and 38 true
positives in one window), so the result is suggestive and consistent with every earlier measurement
of the same effect.

**Finding 2: right site, wrong polarity.** Comparing false-positive and false-negative positions
directly gave one exact match, `chr20:3332481`. GIAB truth has a 4 bp deletion (`TTTTA` to `T`)
and the caller called a 4 bp insertion (`T` to `Tttta`) at the same position. This is a
length-change-direction error at what is almost certainly a poly-T homopolymer, where two local
interpretations fit the reads and the wrong one was chosen. It occurred once in this window (1 of
16 false positives). It is reported as a distinct failure class. Whether `extract_bubble` and
`scan_pair` in `caps_caller.h` have a systematic insertion-versus-deletion bias at tandem repeats
is an open question, and the sample is too small to answer it.

## Checks made before results were reported

Two measurement bugs were caught and fixed before any caller result was reported. A synthetic-truth
generator in a benchmark script silently shifted indel positions, which was a script defect and not
a caller defect. An apparent runtime blow-up in small-sample tests proved to be a sampling artifact
once the trend across sample sizes was examined.

A third check covered the tunable parameters. An audit before the full-scale run looked for
accidentally fitted constants. It found that the multi-polymorphism tolerance of the de Bruijn
bubble channel (`CAPS_DBG_MAXPOLY`) still defaulted to 3 in the compiled binary (DiscoSNP++'s own
`-P` default), while every reported number for that channel had been measured with `MAXPOLY=1`
passed on the command line. A bare invocation would have used a configuration that scores worse
(F1 0.861 at P=1 with coherence checking against 0.850 at P=3, on the same window). The default was
changed to 1 (the current value in `include/caps_caller.h`), and a bare run was confirmed to
reproduce the previously reported TP, FP and FN.

## Could DiscoSNP++ or Kmer2SNP do this?

**Kmer2SNP** has no indel model and no multi-allelic model, because it is a k-mer-graph SNV caller.
This follows from its design as described in the heterozygous-indel literature survey, and is
consistent with its absence from T2.3 to T2.5 (`NOT_APPLICABLE`, not a measured 0). Whether it
could be extended was not examined, because that would require reading its source. The limitation
is therefore an observed capability and not a verified architectural impossibility.

**DiscoSNP++** produced 0 of 26 on T2.4. This is explained by its output schema. A bubble record
has one REF and one ALT, so a bubble-based caller cannot express more than one alternate allele at
a locus. The explanation is inferred from the documented output format and from the observed
behavior (no multi-allelic record among 3,989 emitted, on data whose truth set contains 26
multi-allelic sites). DiscoSNP++'s source was not read, so the explanation is not verified at the
code level, which is a weaker check than the PgRC2 source check in
[`../claim3/mechanism_insight_claim3.md`](../claim3/mechanism_insight_claim3.md).

## Summary

The calling results of Claim 2 follow from the compression mechanism of Claim 1. The pressure to
minimize bits that makes a compressor place divergent alleles at separate pseudogenome addresses is
what a naive caller reads as unrelated single-allele sites. Reconciling those addresses before
calling (`build_substrate` and `collapse_contigs`) recovers the signal and raises the F1 from 0.431
to 0.888 in the ablation. T2.4 (multi-allelic) and T2.5 (tetraploid) show the same mechanism at a
higher allele count and a higher ploidy. DiscoSNP++'s inability to represent multi-allelic sites
follows from its biallelic bubble output format, as inferred above.
