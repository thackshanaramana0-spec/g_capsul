# LIMITATIONS — what this work does not establish

> **Status:** frozen 2026-09-10 at tag `v1.0-capsule` — code and results final.
> Authoritative numbers live in `../benchmark/results/`; verification status in
> [`../AUDIT.md`](../AUDIT.md). Where this file and a result file disagree, the result file wins.

Written from the measurements, not from caution. Each item states the limit, the
evidence for it, and whether a cause is known. Nothing here is hedging language
around a result that is actually fine; every entry is a real boundary.

---

## 0. Sequencing technology: short-read only. Long-read is refused, not attempted.

**This tool does not support long-read sequencing (Oxford Nanopore, PacBio) in
any form.** It was designed and validated exclusively against short-read
Illumina-shaped data (40–301 bases; see the 19 locked datasets). This is a
scope decision, not a partial or degraded capability.

**Why it cannot silently be tried anyway.** Every per-read decode path in the
encoder unpacks one read into a fixed-size stack buffer (`stages/106_inprocess.cpp`,
`MAX_READ_LEN = 1023`), a size chosen for short-read technology. A read past
that length cannot be represented by the current format at all — this is a
structural property of the archive, not a tunable parameter.

**What used to happen, and what happens now.** Before 2026-09-10, an oversize
read was silently dropped with no warning; a file made entirely of such reads
(i.e., any real long-read run) still produced a structurally valid, empty
archive at exit code 0 — a user would see "success" on data the tool did
nothing with. This is fixed: the encoder now refuses unconditionally, before
any read is touched, with the exact count and length of the offending reads
and the structural reason. Verified: constructing a synthetic 2000-base-read
file and running the real binary against it reproduces the refusal on demand,
5/5 repeated trials, with no crash. See `industry/README.md` for the full
account, including a second, independent defect (silent corruption of read
*names*, unrelated to read length, on pathologically long header lines) found
and fixed the same way.

**The decision, stated plainly: long-read data is out of scope. The fix was
not to add support — it was to make the tool say so instead of pretending to
succeed.**

---

## 1. Scope of the human evaluation

**Claims 2 and 3 are chr20 only, at 30×.** Four GIAB individuals (HG002–HG005),
the whole chromosome rather than windows, but one chromosome — because the
archives are chr20. Nothing here is evidence about whole-genome behaviour.

### Why chr20 is the scope, stated against what this field actually does

This is a real limit and is not argued away. But it is worth recording what
the comparable published work evaluates on, because the answer is not
"whole genome":

| work | evaluation data |
|---|---|
| **eBWT2SNP** (Prezza et al., *Algorithms Mol Biol* 2019) | **simulated** chr22 (29×) · **simulated** chr16 (22×) · real chr1 (43–47×) |
| **DiscoSNP++** (its own published indel figures) | **simulated** human chr1 |
| **SKA lo** (*MBE* 2025) | 55 *S. aureus* samples — bacteria, multi-sample |
| **PgRC / PgRC2 / SPRING** | compression only; no variant-calling evaluation at all |

Single-chromosome evaluation is the norm in reference-free variant calling,
and two of eBWT2SNP's three experiments are on **simulated** reads. Against
that baseline this work uses **real data throughout, four individuals rather
than one, and the complete chromosome rather than a window** — plus 19 real
datasets (55.0 GB) for Claim 1, which the calling papers do not attempt.

The honest position: whole-genome would strengthen Claims 2 and 3, and it is
the obvious next experiment. It is not a gap that puts this work below its
field's standard.

**What whole-genome would actually cost**, projected from measured HG002
chr20 figures (chr20 is 1/49 of the genome; the projection is linear and the
non-linearity caveat below is the important part):

    input FASTQ    ~190 GB per individual   (~0.8 TB for four)
    archive         ~26 GB per individual   (~105 GB for four)
    compress         ~2.9 h per individual
    call             ~1.3 h per individual  (~17 h compute for four)

**Compute is not the blocker; k-mer spill disk is.** The caller spills ~25 GB
at chr20 scale in the current format, which is ~1.2 TB at whole-genome scale
and exceeds the 233 GB disk these results were produced on. A superkmer spill
format that would cut this to roughly 64 GB was specified in detail
(`docs/_removable/SUPERKMER_PLAN.md`) and **never built**. Anyone attempting
whole-genome should expect to build that first; it is the single concrete
engineering prerequisite, and it is named rather than left to be discovered.

**No non-human diploid variant validation.** The other 15 datasets carry Claim 1
only. There is no truth set of comparable quality for them, which is a property
of the field rather than a shortcut taken here.

**Sample sizes, stated plainly.** T2.1/T2.3 rest on **four** individuals.
T2.4's census is **26** sites — the complete population for chr20, not a
sample, but 26 is a small number and a reviewer is entitled to say so. T3.4 is
**400** sites (4 × 100), of which only HG002's 100 were independently re-run
during the 2026-09-10 audit. None of these are underpowered *for what they
claim*, but none of them is large.

## 2. The published loss, and why it is published

**HG005 loses both T2.1 (0.834 vs 0.863) and T2.3 (0.593 vs 0.605).** The cause
is identified by controlled experiment, not guessed:

| | read length | distinct lengths | TP | FP | F1 |
|---|---|---|---|---|---|
| HG002/3/4 | 148 bp fixed | 1 | — | ~1,700 | 0.888–0.891 |
| HG005 as-is | 250 bp trimmed | **216** | 30,404 | **3,426** | 0.834 |
| HG005 truncated to 148 bp fixed | 148 bp fixed | 1 | 32,825 | **1,291** | **0.897** |

Same individual, same depth, same reads, same caller, same truth — only the
geometry changes. False positives fall 2.7× and HG005 becomes our *best*
individual. **The caller is tuned for fixed-length reads**, and HG005 is the only
variable-length set in the benchmark.

**The truncated figure is NOT reported as an HG005 result and must never be.**
Truncating discards 40% of every read and scores a different dataset. The table
keeps 0.834. The truncation exists to identify the cause.

*Actionable form:* if variable-length input matters for a deployment, the
caller's length-dependent parameters are the place to look. Scoped work with a
known target.

## 3. Indels: a win, with a hard ceiling underneath it

We win 3 of 4 individuals — and the *class* has a bound that no implementation
escapes.

**Recall is capped by representation, not by thresholds.** At the missed sites
the branching node has coverage 12–30 but exactly ONE successor, and that holds
with **no coverage floor at all** (raw successor counts `[0,0,21,0]`). The
alternate haplotype's k-mers are absent, not filtered. **66% of missed indels are
homopolymer-length changes**: for a `G→GT` inside a T-run the REF and ALT probes
are literally the same string. No bubble exists — for us or for **any** de Bruijn
caller, DiscoSNP++ included.

**Reference context cannot fix precision**, measured over all 5,001 scored
indels: true positives are **more** homopolymeric than false ones (frac run ≥ 4:
0.531 vs 0.427), so any homopolymer/STR filter removes TPs preferentially. That
is exactly what a window-scale test measured (dropped 9 TPs to remove 7 FPs).

**Thirteen mechanisms were implemented and measured; three helped.** Further
filter work on this class should be considered refuted in advance.

Both tools sit far below reference-based callers (indel F1 0.89–0.97 for
DeepVariant/DRAGEN/GATK). **That gap is the price of needing no reference**, and
it is a property of the problem, not of either implementation.

## 4. Compression is slower than SPRING

**2.99× on compress, 1.79× on decompress.** We are faster than SPRING on 3 of 19
datasets and lighter on 9 of 19. Assembly costs time; that is the trade, and it
is stated rather than buried.

The single largest remaining cost is a deliberate size/speed trade: the adaptive
candidate grid is **56% of runtime for 0.586% of archive size** (0.010% against
the second-best point). That is a decision, not a defect, and it is the user's
to make.

## 5. `query` is not a speed win

`genocat --head=100` extracts in **0.19 s** against our **0.46 s**. With the
optional sidecar index we reach 0.03 s, but Genozip is the honest comparator.

**98.6% of a query is spent rebuilding the pseudogenome** in order to hand back
1.6 MB of it. The answer itself is free (1.3% of the time).

**The sidecar is a cache, not a mechanism.** It sidesteps the reason query is
slow instead of removing it, and it is not part of the archive.

## 6. The addressability residue, and a real architectural limit

Content addressing returns both alleles at **345 of 400** sites. The 13.8%
residue is **explained, not unknown**: every one-allele site resolved MULTIPLE
loci that agreed, so it is not a missing haplotype. It is the mechanism's own
complement — a het variant is stored either split across haplotype contigs
(recovered by content addressing) or as **per-read deviations** on one contig.

**True-read retrieval requires the optional sidecar**, and this is a genuine
architectural limit rather than an unfinished feature: `mm_sym` is coded with an
adaptive model in **read order**, so read *k*'s deviations require decoding all
*k−1* before it. They cannot be reached from the archive on demand.

## 7. What the export ratio does and does not measure

`export` is 129–784× faster than SPAdes — and **our export emits 2 records (the
pseudogenome) where SPAdes emits tens of thousands of biological contigs.** The
ratio measures time-to-a-reference-free-coordinate-system from an archive that
had to exist anyway. It is **not** "the same output, faster", and the CSV carries
output bytes and row counts for both sides so a reader can see that. T3.2 and
T3.3 carry no such caveat.

## 8. Claims deliberately not made

- **Polyploid calling.** A synthetic-triploid win was measured, then retracted:
  under symmetric normalisation DiscoSNP++ scores 1.000 on that data. What is
  claimed instead is **multi-allelic** calling on real GIAB sites, where
  DiscoSNP++ emits zero such records in 3,989 structurally.
- **Sequence-compression state of the art.** PgRC2's own paper dismisses the
  streaming-regime tools; mstcom wins only 2 of 15 in its own benchmark. We beat
  PgRC2 by 1.88% on the DNA stream and say exactly that, not more.
- **Statistical significance.** Not reported, consistent with the field: no paper
  in FASTQ compression or reference-free variant calling reports significance for
  these comparison types.

## 9. Where the compression margin does *not* hold

At **2.6× coverage** (a dataset outside the locked suite) we lose to SPRING by
3.97% after all available fixes. The gap is architectural: SPRING matches
read-to-read at ~2 bits per mismatch; we match region-to-region at ~21 bits.
Both routes around it were built and refuted — the matches are the wrong *shape*
to convert (mean MEM length 44.9 bases against a 150 bp read), and removing the
reference metadata entirely costs more than it saves (BSC on the raw region is
1.58 MB **worse**).

Context: PgRC2 crashes outright on that dataset, and SPRING's own authors never
benchmarked the regime. Every locked dataset is normal coverage, where we win
19/19.

## 10. Known-unfixed, stated rather than discovered by a reader

- **Line 3 mode (c)** — a `+` line that is neither bare nor an exact repeat of
  the header is **not supported**, and the archive correctly refuses to claim
  losslessness if it occurs (checked with a deliberately mixed test file, not
  assumed).
- **`SECOND_SELF`, `MEM_MAXMM`, `LOCALITY`, `COST_GATE`** ship OFF. Each is
  correct and lossless; each was measured not to pay. They remain in the tree as
  the evidence for that.
- **The format has no payload checksum.** Truncation and corruption are caught
  structurally (every stream length must read completely), not cryptographically.
- **No random access by read id.** `query` is by pseudogenome coordinate or by
  content.
- **T2.4's `capsule_one` / `capsule_none` columns read `NOT_MEASURED`.** The
  scoring script's CSV schema had those columns from the start but nothing ever
  computed them — the values published there before 2026-09-10 (`4`, `1`) were
  never measured. The script now computes them; the CSV records
  `NOT_MEASURED` rather than back-filling a number no run produced. Only
  `capsule_both` (17) is a measurement.
- **T3.4 requires an optional sidecar built with `CAPS_PILEUP=1`.** Without it
  the coordinate arm scores 0 instead of 18, because `query` emits the
  consensus rather than the reads. The flag is now documented in
  `benchmark/documentation/REPRODUCE_EVERYTHING.md` and `RESULT_CODE.md`; it
  was documented nowhere before the 2026-09-10 audit, and cost three false
  "unreproducible" verdicts during that audit before the cause was found.
- **One published figure was withdrawn during the final audit.** T2.4 was
  published as 21/26 with no supporting raw log; three independent
  re-derivations gave 17/26 and the number was corrected everywhere. Recorded
  because the base rate of such defects in this work is demonstrably not zero:
  one in ten tables, found only because someone went looking. See `AUDIT.md`.
