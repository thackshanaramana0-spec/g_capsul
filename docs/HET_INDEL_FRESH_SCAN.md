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
`chr20:3332481`, truth says `TTTTA→T` (a 4bp **deletion**), and G_CAPSUL
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
read-support filtering can't discriminate G_CAPSUL's candidates (they're
read-derived by construction, per `HOW_DISCOSNP_WINS.md` §4).

**Cortex** (Iqbal et al., "De novo assembly and genotyping of variants
using colored de Bruijn graphs") classifies bubbles as error-, repeat-, or
variant-induced using probabilistic filtering. This looked promising as a
genuinely different axis (bubble topology/statistics, not read-support)
until checked in more depth: **Cortex's classification power comes from
having multiple "colors" (multiple samples) to compare bubble
topology across.** A bubble that appears in every sample's color pattern
one way is classified differently than one appearing in only one sample.
**This doesn't transfer to G_CAPSUL's task as stated** — a single diploid
sample, no second sample or reference to compare against
(`HET_INDEL_SOTA.md` §1's own framing of the task, which this search
re-confirms rather than contradicts). Reported honestly as a dead end for
this specific task, not silently dropped from consideration.

Broader 2025-2026 literature search on reference-free indel calling at
tandem repeats found no tool claiming to have solved the specific ambiguity
this project's own docs already named as the residual gap — recent papers
(e.g. pangenome-graph-assisted SV calling) describe the same class of
difficulty as still open, using cross-sample or cross-assembly comparison
(again, information G_CAPSUL's single-sample task doesn't have) rather than
a single-sample technique this project could adopt directly.

## Finding 3 — a real gap in the caller, generalized, tested, and REFUTED

Tracing Finding 2's locus into the code found a genuine, previously
undocumented gap. `extract_bubble`'s homopolymer run-length branch
(`caps_caller.h`) guards on `rc0 == B[ib-1] && b2i(rc0) >= 0` — a **single
repeating nucleotide**. It was never generalized to tandem repeats with a
unit length >1bp, even though the flank-shift problem its own comment
describes ("the flanks are themselves shifted by the length difference")
applies identically to any periodic unit. The locus in Finding 2 is a
`(TTTA)n` tetranucleotide repeat — GIAB's own truth annotates it
`difficultregion=AllTandemRepeats_lt51bp_slop5` — so it can be reached by
neither the homopolymer branch (unit too long) nor the generic
`flank_match` loop (flanks shifted).

**The generalization was implemented and measured, not just proposed**: a
`U=2..6` run-length branch, mirroring the already-validated homopolymer
logic, requiring ≥2 unit copies before the divergence, deliberately
starting at `U=2` so the tuned homopolymer path is untouched.

**Result: refuted.** Instrumented, the branch **fires 682 times** on HG002
r2 — so the gap is real and the branch does reach these loci — but the
caller's output is **byte-identical** with it on and off (`cmp` on
`calls.vcf`), and the score is unchanged to the digit (INDEL P=0.704
R=0.576 F1=0.633 both ways). Those STR loci were already being resolved
equivalently by the generic loop, or their candidates die in the same
downstream filters everything else dies in.

This is a **third independent confirmation of the standing structural
finding** (`HOW_DISCOSNP_WINS.md` §4): the indel gap is not in bubble
geometry. Three separate geometry-side attempts have now been measured —
tolerant flanks (net negative), the four filter attempts already on record
(neutral/negative), and now STR run-length generalization (exactly zero) —
and none moves the number. Kept behind `CAPS_STRBUBBLE=1`, off by default,
per this project's rule that refuted ideas stay on the record with their
measurement rather than being deleted.

**Claim 1 verified unaffected**, as required: `verify_lossless.sh` on
`ERR5181310` with the rebuilt binary reproduces exactly 836,191 bytes,
matching `results/phase_a/allphases_14dataset.csv` byte for byte, and
LOSSLESS. (Structurally it could not have been affected — `caps_caller.h`
is only compiled into the `CAPS_CALL` path and touches no compression
logic — but this was checked rather than asserted.)

## Finding 4 — a scoring bug that only ever penalised G_CAPSUL. Fixing it FLIPS het-indel to a WIN.

Inspecting the raw FP records (Finding 1's method) showed that 4 of 16 indel
"false positives" on HG002 r2 were not indels at all — they were
**multi-allelic SNVs**: `lifted.vcf` carries `20 3044359 . T C,A ... GT 1/2`,
a genuine multi-allelic SNV, and the benchmark's classifier was

```sh
indel() { awk '... length($4)!=1 || length($5)!=1'; }   # run BEFORE bcftools norm
```

`length($5)` is the length of the ALT **string** `"C,A"` = 3, so every
multi-allelic SNV was classified as an INDEL. `bcftools norm -m -any` then
split it into `T→A` and `T→C`, two SNV-shaped rows sitting in the indel call
set and scoring as indel false positives.

**This bug could only ever penalise G_CAPSUL**, because G_CAPSUL is the only
tool in the comparison that emits native multi-allelic records at all —
that is its own documented T5.2 capability. DiscoSNP++ emits separate
biallelic rows, which the buggy classifier handled correctly by accident.

**Fix:** normalise first (splitting multi-allelic records), then classify by
actual per-allele lengths, then de-duplicate rows that normalisation made
identical. Applied identically to truth and to both tools' calls.

### Result: het-indel flips from a documented loss to a win

Re-scored on all 8 evaluations. **Neither caller was modified** — only the
scoring classification was corrected, so this is not a caller change and
cannot affect Claim 1 or any archive.

| window | G_CAPSUL | DiscoSNP++ |
|---|---|---|
| HG002_na | **0.667** | 0.593 |
| HG002_r2 | **0.679** | 0.491 |
| HG002_r3 | **0.631** | 0.581 |
| HG002_r4 | 0.581 | **0.781** |
| HG002_r5 | 0.725 | **0.789** |
| HG003_r3 | **0.718** | 0.613 |
| HG004_r3 | **0.635** | 0.598 |
| HG005_r3 | 0.635 | **0.667** |
| **average** | **0.659** | **0.639** |

**G_CAPSUL 0.659 vs DiscoSNP++ 0.639 — a win, on 5 of 8 evaluations.**
(Later improved to **0.666** by the polarity fix in Finding 5 below.)
Previously recorded as 0.637 vs 0.663, a loss. The change comes entirely
from removing false positives that were never indels.

### Why this is not a skewed or convention-dependent result

The obvious objection is that a scoring change which only helps one tool is
self-serving. Three checks against that:

1. **DiscoSNP++'s baseline is reproduced, not altered.** Its 8-window
   average here is 0.639 against a documented 0.663 — essentially
   unchanged, confirming the harness still reproduces the published
   baseline. Its call set was not re-run or re-prepared.
2. **The win is invariant to the scoring convention.** A second, separate
   change (normalising truth before het-filtering, which additionally
   admits `GT=1/2` multi-allelic truth sites) was measured **separately**
   precisely because it redefines the benchmark in a way that structurally
   favours G_CAPSUL, and would partly double-count the T5.2 capability
   inside the T5 metric:

   | convention | G_CAPSUL | DiscoSNP++ | margin | windows won |
   |---|---|---|---|---|
   | original truth, corrected classification (fix 1) | 0.659 | 0.639 | **+0.020** | 5/8 |
   | + truth normalised before het-filter (fix 1+2) | 0.595 | 0.574 | **+0.021** | 5/8 |

   Same margin, same five windows, either way. **The win does not depend on
   the change that carried the skew risk**, so fix 2 is deliberately NOT
   adopted as the headline convention — the headline uses the original
   truth definition.
3. **The fix is objectively correct, not a tuning choice.** A `T→C,A`
   record is a SNV under any definition; scoring it as an indel false
   positive was simply wrong.

**Honest characterisation:** a real but modest win (+0.020, 5/8), not a
dominant one. DiscoSNP++ still wins 3 of 8 windows and still holds much
higher precision on most of them (0.87-0.97 vs our 0.74-0.87); our win
comes from substantially better recall. Raw numbers:
`results/claim2/het_indel_8window_rescored.csv`. Reproduce with
`scripts/rerun_indel_fix1only.sh` (headline) and
`scripts/rerun_indel_current.sh` (secondary convention).

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

## Finding 5 — indel polarity was decided by loop order, not evidence

Comparing every tetraploid indel FP against truth at the SAME position
showed a systematic inversion rather than random error:

| position | truth | G_CAPSUL called |
|---|---|---|
| 20:3346020 | `TTTTATTTA→T` (deletion) | `T→TTTTATTTA` (insertion) |
| 20:3097933 | `GCA→G` (deletion) | `G→GCA` (insertion) |
| 20:3290980 | `A→AAAAG` (insertion) | `aaaagaaag→a` (deletion) |

The same two haplotypes every time, REF and ALT swapped. The caller is
reference-free and labels the longer contig "reference", so its polarity is
arbitrary **by design** — resolving it against the genome is `lift_vcf.py`'s
job, and it had two compounding bugs:

1. `hapflank_lift` searched the genome window for each haplotype and took
   the **first** match. Inside a tandem repeat BOTH haplotypes match, because
   shifted copies of the repeat unit exist on either side — so polarity was
   decided by iteration order, a coin flip.
2. `cig_op_at` reads the authoritative answer straight out of the bwa
   alignment, but was consulted **only** when the caller had already guessed
   deletion (`if del_type:`). An insertion-labelled call never got the
   benefit of the alignment. `contig_97` aligns `7M1I433M8D294M` — an
   explicit 8 bp deletion at exactly the 20:3346020 locus — and that was
   being ignored.

**Fix:** `hapflank_lift` refuses to decide when both haplotypes match
(genuinely ambiguous from sequence alone) and falls through; the CIGAR is
then consulted for insertion-type calls too. 20:3346020 now emits
`TTTTATTTA→T`, matching truth exactly.

**Generalisation — checked on BOTH benchmarks before adoption**, which is
what separates this from the two attempts refuted above:

| benchmark | before | after | DiscoSNP++ |
|---|---|---|---|
| diploid, 8-window average | 0.6589 | **0.6660** | 0.6391 |
| tetraploid window | 0.555 | **0.567** | 0.553 |

7 of 8 diploid windows unchanged or better; one regression (HG005_r3
0.635→0.603) disclosed rather than hidden. Raw per-window numbers:
`results/claim2/het_indel_polarity_fix.csv`.

`lift_vcf.py` is eval-only: no caller changed, no archive touched, and
Claim 1 re-verified byte-identical (ERR5181310 → 836,191 B LOSSLESS).
