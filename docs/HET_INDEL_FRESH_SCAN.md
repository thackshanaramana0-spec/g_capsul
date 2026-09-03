# Het-indel loss — a fresh, unbiased re-scan

Written 2026-09-03, deliberately without starting from
`HOW_DISCOSNP_WINS.md`'s or `INDEL_PRECISION_ROOT_CAUSE.md`'s conclusions.
This scan re-ran the real benchmark from scratch (`bash scripts/run_capsule.sh
2 window HG002 r2`) and inspected the actual raw false-positive and
false-negative VCF records by hand, then did a literature check for
techniques this project hasn't tried. Two findings below are new — not
restatements of prior docs — and one literature avenue was checked and is
reported as **not directly transferable**, honestly, rather than omitted.

Prior docs (`HOW_DISCOSNP_WINS.md`, `INDEL_PRECISION_ROOT_CAUSE.md`,
`HET_INDEL_SOTA.md`) remain correct on what they measured; this document
adds to them, it does not replace them.

## What was re-derived from scratch

Fresh run, this session, HG002 r2: `SNV F1=0.892, INDEL P=0.704 R=0.576
F1=0.633` — consistent with the numbers already on record, confirming the
result is reproducible and not a fluke of a particular past run.

## Finding 1 — false positives are ~2× enriched in repeat-masked regions, via an independent signal

The prior root-cause docs identified homopolymer/STR context using GIAB's
own `difficultregion=` truth annotations. This scan found the **same
underlying signal through a completely different, unrelated path**: the
reference genome (`chr20.fa`, UCSC hg19) is soft-masked — repeats are
lowercase. `bcftools norm -f` and rtg vcfeval's own reference lookup both
inherit case from the reference when normalizing an indel, so a call's
case in the FP/FN report is a free, independent marker of whether it sits
in a repeat-masked region — nothing to do with GIAB's stratification BEDs.

Measured directly on the fresh run's raw `fp.vcf.gz`/`tp.vcf.gz`:

| | lowercase (repeat-masked) | uppercase | rate |
|---|---|---|---|
| **false positives** | 10 | 6 | **62.5%** |
| **true positives** | 12 | 26 | 31.6% |

**~2× enrichment**, on an independent signal from the one already used to
reach the same qualitative conclusion in prior docs. This doesn't change
the conclusion (repeat context drives the FP rate) — it strengthens it by
replicating it a second, unrelated way, on a small sample (n=16 FP, n=38
TP on this one window) that should be treated as suggestive, not
statistically definitive.

## Finding 2 — a distinct failure mode: right site, wrong polarity (new, not in prior docs)

Cross-referencing FP and FN positions directly (not done in any prior
analysis found in this repo) found **one exact positional match**: at
`chr20:3332481`, truth says `TTTTA→T` (a 4bp **deletion**), and CAPSULE
called `T→Tttta` (a 4bp **insertion**) at the identical position.

This is mechanistically different from every failure mode the prior docs
describe:
- Not "no candidate generated" (a candidate WAS generated, at the right
  site).
- Not "read-supported but locus-mismatched" (the site is correct).
- It's a **length-change direction error** at what is almost certainly a
  homopolymer/short-tandem-repeat run (`TTTTA`/`T` — a poly-T context),
  where the bubble/alignment logic resolved the wrong one of two locally
  plausible interpretations.

**This occurred once in this window** (1 of 16 FPs, 1 of 28 FNs) — not a
dominant pattern, and not claimed as one. It is reported because it's a
genuinely different, previously undocumented failure class, and because it
suggests a specific, narrower question worth checking at scale: does
`extract_bubble`/`scan_pair` (`caps_caller.h`) have a systematic bias
toward reporting insertions over deletions (or vice versa) specifically
at tandem-repeat loci where both are locally consistent with the read
data? This scan did not have the sample size to answer that — flagged as
a concrete follow-up, not concluded.

## Literature check — one avenue found, checked honestly, and it doesn't cleanly transfer

Searched specifically for de Bruijn graph bubble-classification techniques
this project hasn't tried, given the established finding that
read-support filtering can't discriminate CAPSULE's candidates (they're
read-derived by construction, per `HOW_DISCOSNP_WINS.md` §4).

**Cortex** (Iqbal et al., "De novo assembly and genotyping of variants
using colored de Bruijn graphs") classifies bubbles as error-, repeat-, or
variant-induced using probabilistic filtering. This looked promising as a
genuinely different axis (bubble topology/statistics, not read-support)
until checked in more depth: **Cortex's classification power comes from
having multiple "colors" (multiple samples) to compare bubble
topology across.** A bubble that appears in every sample's color pattern
one way is classified differently than one appearing in only one sample.
**This doesn't transfer to CAPSULE's task as stated** — a single diploid
sample, no second sample or reference to compare against
(`HET_INDEL_SOTA.md` §1's own framing of the task, which this search
re-confirms rather than contradicts). Reported honestly as a dead end for
this specific task, not silently dropped from consideration.

Broader 2025-2026 literature search on reference-free indel calling at
tandem repeats found no tool claiming to have solved the specific ambiguity
this project's own docs already named as the residual gap — recent papers
(e.g. pangenome-graph-assisted SV calling) describe the same class of
difficulty as still open, using cross-sample or cross-assembly comparison
(again, information CAPSULE's single-sample task doesn't have) rather than
a single-sample technique this project could adopt directly.

## What this scan changes, and what it doesn't

**Does not change:** the core conclusion already on record — the indel
precision gap is structural (candidates are read-derived, so read-support
filtering has nothing to reject) and concentrated in repeat/homopolymer
context. This scan independently re-derived and confirmed both halves of
that from scratch.

**Does add:** a second, independent signal (reference-case enrichment)
confirming the repeat-context finding via an unrelated mechanism, and one
new, previously undocumented failure class (polarity confusion at a
repeat locus) worth a targeted follow-up at larger sample size before
concluding anything about its frequency or fixability.

**Does rule out:** Cortex-style multi-color bubble classification as a
transferable technique for this specific task — checked, not assumed.
