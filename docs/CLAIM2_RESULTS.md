# CAPSULE Claim 2 — first real GIAB results (2026-09-02)

Real reads, real GIAB v4.2.1 truth, `rtg vcfeval` (the GA4GH engine `hap.py`
wraps), het-restricted, inside the GIAB confident regions — the same
methodology and the same chr20 windows the outer ARCS project used, so these
numbers are directly comparable to its published ones.

**Verdict up front: it works, it generalizes, and it currently LOSES to every
competitor. Precision is strong; recall is the problem, and the cause is
structural (assembly fragmentation), not parameters.**

## The numbers

HG002, chr20, 30× (streamed from the GIAB 300× BAM and downsampled),
~75,000 reads per 400 kb window, **2.4 s per call**.

| window | | TP | FP | FN | **P** | **R** | **F1** |
|---|---|---|---|---|---|---|---|
| r2 (3.0–3.4M) — *tuning window* | SNV | 108 | 7 | 292 | **0.939** | 0.270 | **0.419** |
| r2 | INDEL | 15 | 3 | 51 | 0.833 | 0.227 | 0.357 |
| r3 (4.0–4.4M) — *held out* | SNV | 127 | 7 | 383 | **0.948** | 0.249 | **0.394** |
| r3 | INDEL | 10 | 7 | 50 | 0.588 | 0.167 | 0.260 |

Truth in window: 422 het-SNV / 69 het-indel (r2), 533 / 70 (r3).

## Does it generalize? Yes.

Tuning-window F1 0.419 vs held-out 0.394 — a gap of 0.025. **No overfitting**:
the frozen parameters were ported verbatim from ARCS and never touched, and
the held-out window behaves like the tuning window. That was the specific
risk of developing on a small slice, and it did not materialise. (Precision
is actually *higher* on the held-out window.)

## Where it stands against ARCS — VERIFIED BY RERUN, not quoted

**Correction to an earlier version of this document.** It first compared
CAPSULE against ARCS's *published* 0.923/0.954, which was wrong twice over:
those numbers are **HG001**, while CAPSULE was run on **HG002**, and they
were quoted from `VARIANT_ANALYSIS_MASTER.md` rather than reproduced. ARCS
has now been rerun on the **identical `reads.fq`** through the **identical**
lift-and-score pipeline. The real picture:

| tool (same reads, same pipeline, HG002 r2) | SNV P | SNV R | **SNV F1** | INDEL F1 |
|---|---|---|---|---|
| ARCS, default config | — | — | **0.000** (0 SNVs called) | 0.46 (46 lifted) |
| ARCS, `ARCS_XSNV=1` | 0.981 | 0.635 | **0.771** | 0.710 |
| **CAPSULE (this repo)** | 0.939 | 0.270 | **0.419** | 0.357 |

Two findings that matter more than the ranking:

1. **ARCS's published 0.954 does not reproduce on this input.** Its default
   configuration produces `candidates=0, SNVs=0+0` here — zero SNV calls.
   The published numbers must come from a different read preparation
   (pairing, source, or depth), so **no ARCS figure should be quoted as a
   baseline for CAPSULE without rerunning it on the same input**, which is
   exactly the mistake this document originally made.
2. **CAPSULE beats ARCS's shipped default** (0.419 vs 0.000) because the
   cross-contig bubble pass is default-ON here and off there — the finding
   recorded in `docs/CLAIM2_BUILDER.md`. Against ARCS with that same pass
   enabled, CAPSULE still loses 0.419 vs 0.771. The honest gap is **0.771 vs
   0.419 on equal footing**, not 0.954 vs 0.419.

## Root cause: the contigs are far too short

Measured on r2's own output:

```
contigs = 4,537        for a 400,000 bp window
mean length = 335 bp   (~2.3 reads long)
>= 500 bp   = 11.4%
total contig bp = 1,520,128  (3.8x the window — haplotype fragments never merged)
```

CAPSULE's chaining requires an **exact** suffix-prefix overlap. Every
heterozygous site, every sequencing error at an overlap boundary, breaks the
chain. The result is an assembly that is excellent for *compression* (it is
what wins Claim 1 by 14/14) but is nearly read-length fragments for
*calling*:

- **Pileup calling is starved** — only 13 of 165 calls came from pileup
  columns; there is rarely enough depth at one position on one short contig.
- **Cross-contig bubble calling carries the load** (127 of 165) but needs a
  unique 25-mer anchor with clean flanks on *both* contigs, which a 335 bp
  fragment often cannot provide.
- **73% of true het sites produce no detectable bubble at all** → recall 0.25.

Precision 0.94 says the mechanism is *correct* — what it calls is right. It
simply cannot see most variants.

## Why ARCS does better with the same caller code

ARCS's Method B (`vodbg_pg`) does global greedy-overlap growth with
mismatch-tolerant read placement, producing long consensus contigs where
both haplotypes co-occupy one contig and het sites appear as pileup columns.
CAPSULE's exact-overlap chaining splits them instead. **Same caller, same
frozen parameters, different assembly — and the assembly is the whole
difference.** This confirms, from the opposite direction, the finding
recorded in `docs/CLAIM2_BUILDER.md`.

## The "contigs are too short" hypothesis — MEASURED AND REFUTED

An earlier version of this document blamed short contigs and proposed
mismatch-tolerant chain extension to lengthen them. **Measurement refutes
that.** Comparing the two assemblers' contigs on the same window:

| | contigs | mean len | ≥1 kb | total bp (window = 400 kb) |
|---|---|---|---|---|
| ARCS (`ARCS_XSNV=1`) — 0.771 F1 | 4,307 | **249** | 97 | 1,072,468 (2.7×) |
| CAPSULE — 0.419 F1 | 4,537 | **335** | 326 | 1,520,128 (3.8×) |

CAPSULE's contigs are **longer** on average, with **3× more** ≥1 kb contigs,
and it still gets half the recall. Lengthening contigs would not have fixed
anything — the proposed fix was aimed at the wrong quantity.

## The anchor-supply hypothesis — ALSO REFUTED

Next suspicion: the cross-contig pass needs a canonical 25-mer occurring
*exactly twice* (`occ.size() != 2 → skip`), so maybe fragmentation starves
the anchor pool. Measured 25-mer occurrence multiplicity across each
assembler's contigs:

| | distinct 25-mers | occ==1 | **occ==2** | occ==3 | occ≥4 |
|---|---|---|---|---|---|
| ARCS | 648,780 | 511,804 | **92,163 (14.2%)** | 17,876 | 26,937 |
| CAPSULE | 670,082 | 281,411 | **253,945 (37.9%)** | 61,100 | 73,626 |

CAPSULE has **2.75× more** exactly-2-occurrence anchors than ARCS. Anchor
supply is not the bottleneck either.

## What the evidence actually points to

CAPSULE carries 3.8× the window in contig bases against ARCS's 2.7×. For a
diploid the *ideal* is ~2× — one contig per haplotype. The excess is
**same-allele duplication**: because exact-overlap chaining breaks a chain at
every het site *and* every error, one haplotype's sequence ends up spread
across several separate contigs.

That predicts exactly the failure observed: most of CAPSULE's abundant
`occ==2` anchor pairs are **two copies of the same allele**, which walk
along matching and never diverge, so `extract_snv_bubble` finds no bubble —
and worse, a same-allele duplicate *consumes the 2-occurrence slot* that the
true hap1/hap2 pair needed, pushing the real pair to occ==3 or 4 where the
`occ.size()!=2` test discards it. That is consistent with every measurement
above: more anchors, longer contigs, more total sequence, fewer bubbles.

**Next experiments, in order (both cheap, neither yet run):**

1. **Collapse duplicate/contained contigs before calling.** Targets the
   excess redundancy directly. Note the irony: CAPSULE's MEM self-match
   already finds precisely this redundancy for compression — but the caller
   deliberately captures contigs *pre*-MEM to preserve haplotype separation,
   which preserves same-allele duplication along with it. A dedup pass that
   collapses *identical* fragments while keeping *differing* ones is exactly
   the needed middle ground.
2. **Relax `occ.size()==2` to a small range (2–4)** and test every
   cross-contig pair among the occurrences, keeping bubbles that actually
   diverge-and-reconverge. Directly recovers the true pairs currently
   discarded by the strict test.

**Do not** respond to these numbers by loosening the frozen filter
parameters. Precision 0.94 with recall 0.25 is a *visibility* problem, not a
threshold problem — the variants never reach the filters at all.

## Methodological note for the paper

Two hypotheses were formed and killed by measurement within one session
(short contigs; anchor starvation). Neither was argued away — each was
measured against the alternative assembler on identical input. This is the
same discipline recorded in `docs/FAILURES_AND_REFUTED_IDEAS.md` and it
belongs in the Methods section as evidence the tuning was not blind.

## Reproduce

```bash
bash scripts/run_window_bench_capsule.sh /tmp/caps_call_enc "$PWD/scripts" \
     ~/refs/chr20.fa HG002 r3
```

~90 s end to end per window (13 s stream, 2.4 s call, rest is BWA/rtg).
Windows: `r2` (tuning) and `r3`/`na`/`r4`/`r5` (held out); individuals
HG002/HG003/HG004/HG005. Requires `~/refs/chr20.fa` + its BWA index and
`~/giab_truth/` (both already present — see `docs/CLAIM2_DATA_AND_TOOLS.md`).

## Two bugs found and fixed while getting here

1. `CAPS_DUMP_CONTIGS` writes FASTA, but the pipeline copied ARCS's TSV→FASTA
   `awk` conversion, silently producing a garbage contigs file — symptom was
   `lifted 0 calls`. Fixed in `run_window_bench_capsule.sh`.
2. The pooled `/data/fastq/HG00*_pooled.fq` files are **tile-ordered, not
   coordinate-ordered** (read names run 1101→2212), so slicing them by line
   yields a thin genome-wide scatter, not a window. Windows must be streamed
   by region from the GIAB BAM (`samtools view <url> 20:LO-HI`), which is what
   the script does — 13 s for 75k reads.
