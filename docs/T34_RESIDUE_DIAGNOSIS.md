# T3.4 — where the 13.8% residue actually goes

**Status: SUPERSEDED IN PART — the residue was closed. See
[`CLAIM3_LOCUS_ADDRESSABILITY.md`](CLAIM3_LOCUS_ADDRESSABILITY.md) for the
finished claim.**

This document is the diagnosis trail and is kept for the reasoning and the
refuted hypotheses, per the repo rule that retractions stay in place. It stops
at **391/400**. A third defect was found afterwards — the probe anchored
upstream only — and correcting it reaches **400/400** on the primary archive
and 398/400 on an independently assembled replicate. Sections 1-9 below remain
accurate as far as they go.

**Status at the time of writing: CONFIRMED ON REAL HG002 READS — at 391/400.**

The mechanism reproduces on real data. The headline number does not.

| input | exact+cont | mm+cont | exact+offset | **mm+offset** |
|---|---|---|---|---|
| simulated reads, real variants | 365 | 365 | 368 | **400** |
| **real HG002 reads** | 352 | 352 | 364 | **391** |

The synergy — neither correction doing anything alone, both together closing
most of the gap — holds on real reads. **400/400 was an artifact of simulated
reads and is withdrawn as the headline.** The real-data number is **391/400**,
and §8 records that the nine remaining sites are still unexplained after three
refuted hypotheses.

Controls in §6, real-reads detail in §8, limits in §7.

The published T3.4 figure — **345/400 by content vs 81/400 by coordinate**,
from `scripts/run_locus_fidelity.sh` — is **unchanged and unaffected**. That
script is untouched, and the modified decoder was verified byte-identical to
git HEAD on the default path across all 400 sites (§6).

---

## 1. The explanation on file was insufficient

T3.4 reports 345 of 400 heterozygous sites recovering **both** alleles by
content addressing. The residue of 55 sites has been explained as "variants
stored as per-read deviations on a single contig."

That explanation does not survive being taken seriously. If a site's variant
is a per-read deviation, then the sidecar — which exists precisely to apply
per-read deviations rather than return the consensus (`CAPS_PILEUP=1`, see
`PIPELINE.md`) — reconstructs the read carrying it. The site would score
**both**, not **one**. So the stated mechanism predicts the opposite of the
observed result, and something else is consuming those 55 sites.

This matters beyond bookkeeping. A residue that is a property of the data is a
ceiling. A residue that is a property of our own query and scoring code is a
defect, and defects are fixable.

---

## 2. Two defects, found independently

### Defect A — SEARCH: the probe is matched exactly against a consensus

`query` takes a sequence probe and locates it in the pseudogenome with
`std::string::find`, an **exact** match.

At a heterozygous locus the two haplotypes sit on **separate contigs** — that
is the Claim 2 finding, restated. A 40 bp reference probe therefore matches the
reference-agreeing contig and fails on the alternate contig wherever a *second*
variant happens to fall inside the same 40 bp window. Exact matching resolves
one haplotype and is structurally incapable of resolving the other.

The existing long-probe fallback does not rescue this. It is gated on probes
longer than 60 bp, and T3.4's probes are 40 bp, so it never fires.

### Defect B — SCORING: the two arms were scored by different rules

The v1 scorer requires each returned read to **contain** the reference probe
(`s.find(probe)`) before reading the base 45 positions along.

A read that carries the variant of interest may also differ from the probe
elsewhere — a sequencing error, or that same second variant. Containment throws
away exactly the reads that carry the evidence.

The **coordinate** arm never had this restriction. It reads the base at
`off = v - p` straight from the read's stored placement. So the two arms of the
comparison were being scored by different rules against the same archive, and
the stricter rule was applied to the arm we were claiming was better. The
asymmetry happened to cost us rather than flatter us, which is why it went
unnoticed, but it is not defensible in either direction and it had to go.

These two defects are **independent**. A is in the decoder, B is in the scorer.
Neither was introduced to explain the other.

---

## 3. The prediction, and why it is falsifiable

Both corrections are necessary and neither is sufficient:

- Fixing **search** alone returns the alternate-haplotype reads, but
  **containment** then discards them, because those reads differ from the probe
  in the window. Net effect: nothing.
- Fixing **scoring** alone removes the containment filter, but exact search
  never returned the alternate-haplotype reads in the first place, so there is
  nothing new to score. Net effect: nothing.

So the prediction is a **synergy** — a 2×2 in which three cells are flat and one
moves. This is the sharpest available form of the claim, because almost any
wrong mechanism predicts something additive instead. If either single
correction moves the number on real data, the mechanism stated here is wrong
and the experiment says so.

Note also what this rules out. Offset scoring alone changing nothing is direct
evidence that the new scorer is **not** merely permissive. A metric that had
been loosened into flattering us would move on its own.

## 4. Measured, 400 simulated sites per cell

`scratchpad/t34_sim2.py` — 148 bp reads, 30×, 40 bp probe, diploid, variable
rate of a *second* het variant falling inside the probe window.

| extra-het rate in window | exact + containment | mm=2 + containment | exact + offset | **mm=2 + offset** |
|---|---|---|---|---|
| 0 (clean) | 400/400 | 400/400 | 400/400 | **400/400** |
| 1/1200 (genome average) | 386/400 | 386/400 | 386/400 | **400/400** |
| 1/400 (clustered) | 363/400 | 363/400 | 363/400 | **400/400** |
| 1/200 (highly clustered) | 318/400 | 319/400 | 318/400 | **400/400** |

The predicted shape holds exactly. Three columns are flat across a 26-point
spread in difficulty and the fourth is saturated.

This is the same signature Claim 2's ablation already carries — collapse and
re-placement give 0.431 / 0.426 / 0.648 / 0.888, neither alone, both together.
The recurrence is not a coincidence. Both are cases of the same underlying
fact: **a heterozygous locus is not one place**, so any single-address,
single-representation assumption fails, and undoing one such assumption while
leaving the other in force changes nothing.

## 5. Implementation and verification

**Decoder** (`stages/capsule_decode.cpp`):

- `caps_query_mm()` reads `CAPS_QUERY_MM`, defaulting to **0**.
- `find_occurrences(hay, pat, k, out)` does pigeonhole seed-and-extend — split
  the probe into k+1 seeds, at least one of which must match exactly under k
  mismatches, then verify each candidate. This is the same matching strategy
  the **encoder** already uses to place reads onto the pseudogenome, so it is
  the archive's own primitive applied to reading rather than to writing. It is
  not a new algorithm bolted on, and it is not a brute-force scan.
- A seed-length floor of 12 bp caps the effective k downward for short probes,
  rather than degrading into meaningless 3 bp seeds.
- Both query paths (sidecar and main) call it, and each emits
  `[query] occ <start> <end> <strand>` on stderr, so a consumer can score by
  placement offset without re-deriving the occurrences.

**Verification, all re-run 2026-09-15 and passing:**

| check | result |
|---|---|
| `scratchpad/t34_verify.py` | 4/4 PASS, including the offset scorer finding no ALT on homozygous sites |
| `scratchpad/test_matcher.cpp` (helper copied verbatim from the shipped source) | 4/4 PASS over 8,000 randomised cases against brute force |
| `scratchpad/test_scorer_v2.py` end-to-end | PASS at 10/10/10/**20**, allele balance 1.00, no cell flattered |
| full-TU syntax check of the edited decoder | 0 errors, and mutation-tested at **both** call sites to prove the check is live there |

**The default path is provably unchanged.** With `CAPS_QUERY_MM` unset, k=0,
and matcher test [3] confirms k=0 is byte-identical to `std::string::find`.
`run_locus_fidelity.sh` discards `query`'s stderr entirely and reads only the
`.fa`, so the new stderr line is invisible to it. The published 345/400 cannot
move.

The decoder now also **builds and runs on Linux** with the change in place, so
the earlier Windows shim-header syntax check is superseded by a real build.

---

## 6. CONFIRMED ON A REAL ARCHIVE — 2026-09-15

Before the run, three prerequisites were **verified rather than assumed**, each
of which would have invalidated the result silently:

| prerequisite | check | result |
|---|---|---|
| reference slice and coordinate frame are correct | every GIAB REF base compared against the fetched sequence | **870/870 match** |
| the archive is a valid capsule, not a partial one | decode, then `cmp` the sequence column against the input FASTQ | **byte-identical, LOSSLESS** |
| the sidecar carries real placements | `capsule_decode index` output | **122,026 placements, 28,337 reads carrying 33,300 deviations** |

The third check caught a real mistake in the first attempt. The initial encode
produced `pos_abs` at 0 B and a sidecar reporting **0 placements**, which would
have made every arm score zero. The cause was a missing `DUMP_PERM`, the
harness variable `PIPELINE.md` already documents — the encoder's whole
position/strand/length block is behind `if(getenv("DUMP_PERM"))`. Had the
sidecar not been inspected, this would have been read as a failure of the
mechanism.

### The headline result, 400 real GIAB HG002 het SNVs

| configuration | both | one | neither |
|---|---|---|---|
| coordinate (control) | 54 | 346 | 0 |
| content: exact + containment | 365 | 35 | 0 |
| content: mm=2 + containment | 365 | 35 | 0 |
| content: exact + offset | 368 | 32 | 0 |
| **content: mm=2 + offset** | **400** | **0** | **0** |

The predicted shape holds on real data. Search alone moves nothing at all
(365 → 365). Scoring alone moves almost nothing (365 → 368). Together they
close the entire residue. The corrections are **synergistic, not additive**,
exactly as §3 predicted and as Claim 2's ablation behaves.

A worked instance, from the first site queried: exact search returned three
occurrences, and `mm=2` returned a fourth at pg position 1,134,029 — roughly
300 kb away from the other three. That is the alternate-haplotype contig, found
exactly where §2 says it must be.

### Negative control — 400 homozygous positions

A method that finds two alleles everywhere would also score 400/400, so the
same pipeline was run over 400 positions with no GIAB variant within ±200 bp
and identical bases on both simulated haplotypes.

| configuration | both (false positives) |
|---|---|
| exact + containment | 4 / 400 |
| exact + offset | 6 / 400 |
| mm=2 + offset | **18 / 400** |

Allele balance separates the two sets by two orders of magnitude: **1.17** at
the real het sites versus **0.01** here.

The 18 false positives are disclosed rather than tuned away. They are a
property of the inherited **"both" rule**, which counts a site as recovering
two alleles on the strength of a *single* supporting read — at 30× with a 0.2%
substitution rate, one miscalled base satisfies it.

### Does 400/400 survive a stricter rule?

This is the question that decides whether the result is signal or a loosened
metric. If the positive set only holds because the rule is permissive, the
negative set will hold with it.

| rule | het "both" | homozygous "both" |
|---|---|---|
| minor allele ≥ 1 read (the inherited rule) | 400 | 18 |
| **minor allele ≥ 2 reads** | **400** | **9** |
| minor allele ≥ 3 reads | 399 | 8 |
| minor allele ≥ 5 reads | 395 | 6 |
| minor-allele fraction ≥ 0.05 | 393 | 4 |
| minor-allele fraction ≥ 0.10 | 383 | 3 |

**400/400 survives at ≥ 2 reads while the false positives halve.** The
underlying distributions are effectively disjoint:

| minor-allele fraction | p5 | p25 | median | p75 | p95 |
|---|---|---|---|---|---|
| het sites | 0.121 | 0.333 | **0.409** | 0.458 | 0.500 |
| homozygous sites | 0.000 | 0.000 | **0.000** | 0.000 | 0.000 |

A true heterozygous site sits near 0.5 and the real sites cluster there. The
homozygous set has no minor allele at all through the 95th percentile. The
recovery is real signal.

### No regression on the default path

All 400 sequence queries were re-run with a decoder built from **unmodified git
HEAD** and compared byte-for-byte against the modified binary under default
environment:

```
identical=400  differing=0  missing=0
```

The published 345/400 cannot move.

### Reproducing it

`scripts/t34_realdata/` — see its README for the exact command sequence and the
checks that must pass first. About four minutes end to end, including the
downloads, on a 7 GB box.

On the benchmark box with the real HG002 archives, the equivalent experiment is
`scripts/run_locus_fidelity_v2.sh`.

## 7. Limits of this result — read before citing

**The reads are simulated.** The reference is real, the variants and their
phasing are real GIAB, and the encoder, archive, sidecar, query and scorer are
all the real shipped code — but HG002's actual FASTQ is ~200 GB and the box
used here has 7 GB of RAM. So this run confirms the *mechanism* on a real
archive. It is not a substitute for the full HG002 measurement, and it should
not be reported as one.

Also scoped out:

- **SNVs only.** Indels were not modelled in haplotype construction, so the 400
  sites are real GIAB het SNVs. T2.3's indel behaviour is untouched by this.
- **One 602 kb window** of chr20, the same window T3.4 publishes over.
- **Uniform read sampling**, no real coverage bias, no PCR duplicates, no
  platform-specific error profile.
- **18/400 false positives** on homozygous sites under the inherited rule, 9/400
  at ≥ 2 reads. This is a genuine cost of mismatch-tolerant search and belongs
  next to the 400/400 wherever it is quoted.

Per this repo's standing rules, this result is recorded as measured, including
the false-positive rate and the missing-`DUMP_PERM` misstep that preceded it.

**§7 as written above applied to the simulated-read run. The "uniform read
sampling" limitation was then removed — see §8, which supersedes it.**

---

## 8. REAL HG002 READS — 391/400, and 400/400 withdrawn

The simulated-read limitation turned out to be unnecessary. T3.4 queries only
chr20:3.0–3.6 Mb, so the reads for exactly that window can be range-fetched
from the GIAB 300× Illumina BAM over HTTP in about 30 seconds. No 200 GB
download is involved, and the earlier claim that one was is withdrawn — the
project's own published input (`HG002_pooled.fq`, `DATASET_LOCKED.md`) is chr20
at 30×, roughly 4 GB, not a whole genome.

**117,565 real HG002 reads**, 148 bp, 28.9× over the window, real error
profiles and real coverage bias. Archive verified **LOSSLESS byte-identical**.
Sidecar: 117,565 placements, 16,415 reads carrying 42,800 deviations.

| configuration | both | one | neither |
|---|---|---|---|
| coordinate (control) | 96 | 303 | 1 |
| content: exact + containment | 352 | 48 | 0 |
| content: mm=2 + containment | 352 | 48 | 0 |
| content: exact + offset | 364 | 36 | 0 |
| **content: mm=2 + offset** | **391** | **9** | **0** |

Allele balance 1.20 at het sites versus 0.01 at the homozygous control.
False positives on 400 homozygous positions: **21/400**.

The coordinate arm lands at 96/400, close to the published 81/400 — a useful
sign that this setup behaves like the real benchmark rather than like a
simulation.

**What holds:** the synergy. Search alone still moves nothing (352 → 352).
Scoring alone moves a little (352 → 364). Together they recover 39 sites.

**What does not hold:** 400/400. That was simulated reads. The real number is
**391/400** and the manuscript should say 391.

### The nine remaining sites — three hypotheses, all refuted

They are not noise, and they are not yet explained. What was ruled out, each by
measurement rather than argument:

1. **"They need more mismatch tolerance."** Their local GIAB variant density is
   3.5× that of successful sites (1.78 vs 0.50 within ±60 bp) and the search
   returns 3.2 occurrences against 14.5, which made this the obvious reading.
   **Refuted.** Raising `CAPS_QUERY_MM` to 3, 4, 6, 8 changed nothing at all —
   the seed floor silently caps a 40 bp probe at k=2, since pigeonhole needs
   k+1 seeds of ≥ 12 bp. Sweeping probe length to lift that cap does not
   recover them either:

   | probe | k_max | het both | homozygous false |
   |---|---|---|---|
   | 40 | 2 | 391 | 21 |
   | 60 | 4 | 391 | 15 |
   | 80 | 5 | 389 | 11 |
   | 100 | 7 | 383 | 12 |
   | 120 | 9 | 379 | 10 |

   Longer probes buy **specificity** (21 → 11 false positives) and never buy
   the nine sites. Tolerance is not the limit.

2. **"The ALT allele is absent from the 30× subsample."** **Refuted.**
   `samtools mpileup` on the input BAM shows every one of the nine carrying 6–15
   ALT reads at depth 15–30, indistinguishable from the control group.

3. **"The ALT context is missing from the pseudogenome."** **Does not
   discriminate.** Present for 6 of 9 failing sites and 6 of 9 succeeding ones.

What *is* established: the archive round-trips byte-identically, so the
ALT-bearing reads are provably inside it. These nine are a **retrieval**
failure, not a storage failure. Characterising them is open work, and the
`k_max = floor(P/12) − 1` cap is worth recording independently — it means
mismatch tolerance on a 40 bp probe cannot exceed 2 whatever the flag says.

### Reproduction

`scripts/t34_realdata/` for the pipeline. Real reads come from:

```bash
B=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam
samtools view -b -F 0x900 -s 42.10 "$B" 20:2999001-3601000 > win30x.bam
samtools fastq -n win30x.bam > real_reads.fq
```

Numbers in `results/T34_REALDATA_20260915/REAL_READS.csv`.

---

## 9. What locus addressability costs, and how to account for it

**The accounting rule, because it decides how every number here is quoted.**
The sidecar is built FROM the capsule, after compression has finished and the
archive is sealed. It is a downstream operation with its own costs, exactly as
Claim 2's variant calling is. It does not change the archive, it does not change
the compression ratio, and its bytes must never be added to a Claim 1 number.

This is the same accounting the repo already uses: `claim2_T2.1_snv.csv` carries
`wall_s` and `peak_ram_kb` for the calling step and no size column at all,
because the VCF is derived output rather than stored cost. If encoding peaks at
4 GB and a downstream step peaks at 5 GB, the compressor is not a 5 GB tool.
Two operations, two rows.

**The headline row**, on 117,565 real HG002 reads:

| | wall | peak RAM | needs a reference |
|---|---|---|---|
| **ours — `capsule_decode index`** | **0.31 s** | **17 MB** | **no** |
| bwa index + mem + sort + index | 3.03 s | 117 MB | yes |

**9.8× faster, 6.7× lighter, and no reference genome.** Index build time is a
selling point in this literature rather than an embarrassment — CRAM advertises
exactly this comparison against BAM.

The sidecar produced is 1,589,156 B. It is reported, not charged against
Claim 1, and it is rebuildable in 0.31 s.

**Per-query latency**, 50 loci:

| arm | ms/query | content-addressed |
|---|---|---|
| samtools view (coordinate) | 4.8 | no |
| ours (coordinate) | 11.2 | no |
| ours (sequence, exact) | 25.2 | yes |
| ours (sequence, mm=2) | 52.0 | yes |

Like-for-like we are **2.3× slower**, and that is the honest weakness. Both
figures are instantaneous in absolute terms, and 3.8 ms of ours is per-invocation
process startup, since every query relaunches the binary and reloads the index.
An earlier "11× slower" figure was wrong: it compared our SEQUENCE search
against their COORDINATE seek, which are different operations.

Content addressing has no baseline at all. `samtools` cannot answer "give me the
reads containing this sequence" without a reference genome to align against.

Numbers in `results/T34_REALDATA_20260915/SIDECAR_COST.csv`, scripts
`sidecar_cost2.sh` and `query_breakdown.sh`.

**Caveat.** This is one 602 kb window. bwa's 117 MB is largely fixed overhead
that will not grow much with input size, while our 17 MB will, so the memory
margin is the figure most likely to shrink at chr20 or WGS scale. Re-measure
there before putting it in the manuscript.
