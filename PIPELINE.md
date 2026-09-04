# PIPELINE.md — G_CAPSUL is one pipeline

**This file is the entry point.** Claims 1, 2 and 3 are not three tools that
happen to live in one repository; they are three outputs of a single pass over
the reads. Any document that presents them as independent predates this file
and should be read against it.

---

## 1. The pipeline

```
                      reads (FASTQ)
                           |
        ┌──────────────────┴──────────────────┐
        │  TRUNK — built once, shared by all  │
        │  1. greedy overlap chaining         │
        │  2. pigeonhole mapping of leftovers │
        │  3. second-region assembly          │
        │     -> pseudogenome + contigs +     │
        │        per-read placements          │
        └──────────────────┬──────────────────┘
                           |
      ┌────────────────────┼────────────────────┐
      |                    |                    |
   CLAIM 1              CLAIM 2              CLAIM 3
   MEM self-match       k-mer graph (kc)     export / coverage / query
   + entropy coding     -> bubbles -> VCF    served from the archive
   -> .capsule archive
```

The trunk costs **3.6 min and 2.35 GB** at full chr20 and is paid once. Claim 2
consumes what the trunk already built; it never re-assembles.

**One caveat stated plainly:** Claim 2 currently calls at *compress time*
(`CAPS_CALL=1`), from in-memory state. Claim 3 reads a saved archive. So the
honest sentence is "compress once and get calls as a byproduct", not "open an
old archive and call variants". Calling from a stored archive is not wired.

## 2. Running it

```bash
scripts/build106.sh /tmp/capsule            # build (must include -fopenmp)

# Claim 1 only — archive, no calling. Bit-identical to the pre-caller encoder.
/tmp/capsule reads.fq 3 16 16 22 16 16 1 24 64 1

# Claim 1 + Claim 2 — archive AND variant calls, one pass
CAPS_CALL=1 CALL_VCF=out.vcf /tmp/capsule reads.fq 3 16 16 22 16 16 1 24 64 1

# Claim 2 graph caller alone (Method B), for head-to-head benchmarking
CAPS_CALL=1 CAPS_DBG=1 CAPS_DBG_ONLY=1 CALL_VCF=out.vcf \
  /tmp/capsule reads.fq 3 16 16 22 16 16 1 24 64 1
```

Benchmarks:

| script | what it runs |
|---|---|
| `scripts/run_pipeline_bench.sh` | **the whole pipeline** — compresses all 19, calls variants on the 4 human sets, then export/coverage/query on those same archives |
| `scripts/run_claim1_bench.sh` | compression + lossless + SPRING/Genozip, **all 19 datasets** |
| `scripts/run_fullchr20_bench_capsule.sh` | Claim 2 at full chr20 scale |
| `scripts/run_fullchr20_bench_disco.sh` | DiscoSNP++ arm, identical methodology |
| `scripts/run_window_bench_capsule.sh` | Claim 2 on one 400 kb window (fast iteration) |
| `scripts/run_claim3.sh` | Claim 3 |

## 3. The graph caller (Method B) — what it is

`kc`, the canonical 31-mer table, is built **by the caller** (`kc_H_build`) and
was previously used for a single scalar (the coverage threshold `H`). It **is**
a de Bruijn graph: a k-mer set plus `kc_find` as the membership oracle, so
bubble calling reuses a table that already had to exist for `H` — but the table
is a CALLER cost, not a compressor byproduct. Do not write "the graph is free
because the compressor built it"; it is not true and `106_inprocess.cpp`
contains no k-mer table.

What the compressor genuinely hands over is the **assembly**: Method B reuses
the encoder's contigs and skips `build_substrate`, which re-places all 12.6M
reads at a measured 738 s serial on full chr20.

Stages: `kc` build -> branching nodes, both strand orientations -> lockstep
bubble walk with bounded multi-polymorphism -> read-coherence filter (containment
+ base quality) -> coverage ceiling -> VCF.

### Parameters — all derived or structural, none fitted

| parameter | value | basis |
|---|---|---|
| `MINC` | 2 | structural: discard singletons only. The count histogram is bimodal (60.4% at 1, gap, 35.6% at >=10); raising it to 5 removes 5% of nodes for no compute saving while deleting the 2-4 band where low-coverage het alleles live. |
| `MINQ` | 20 | phred convention (1% error) — a property of the platform, not the sample |
| `COHC` | `max(2, H/10)` | **derived from measured haploid depth.** Required read support scales with coverage; a constant 2 is right at 30x and wrong elsewhere. |
| `COVCAP` | `4 x H` | **derived.** False positives carry 12x the allele depth of true positives (168.3 vs 13.9) — they are repeat/paralog collapses pooling coverage from every genomic copy. A real het allele would need 4x its expected depth to trip this. |
| `MAXPOLY` | 1 | structural: at k=31 and ~1/1000 heterozygosity, a second het site inside one k-mer window is uncommon and a third is rare |
| ploidy gate | 0.024 | HETSCAN pair-fraction threshold, validated 5/5 |

**A cautionary result kept deliberately:** `MINC` was once *derived per-dataset*
from the k-mer histogram valley. That generalised **worse** — it drifted upward
(4 -> 5) at full scale, because absolute error counts grow with data and drag
the valley with them, costing 0.057 recall. Derive what genuinely scales with
the data (read support does); keep structural what does not (the singleton
floor).

### The 19 datasets

15 non-human (`NEW_DATASET_LOCKED.md`) plus the 4 GIAB human chr20 sets
(HG002, HG003, HG004, HG005) = **19**. **The
human sets are not a separate category** — they are FASTQ, they compress
through the identical code path, and they belong in the Claim 1 table with
everything else. What distinguishes them is only that they are diploid human,
so Claim 2 also calls variants on them, from that same compression pass. An
earlier version of `run_claim1_bench.sh` excluded them purely because it looked
for `<name>_1.fq` and the GIAB files are `<name>_pooled.fq` — a filename
artefact, not a property of the data.

## 4. Results — full chr20, measured

| | HG002 | HG003 (held out) | DiscoSNP++ |
|---|---|---|---|
| **SNV F1** | **0.879** | **0.877** | 0.847 |
| precision | **0.951** | 0.945 | 0.951 |
| recall | **0.818** | **0.818** | 0.763 |
| caller time | 98.5 s | 98.9 s | 76.5 s |
| peak RAM | 6.99 GB | 7.12 GB | 3.45 GB |
| archive | **26 MB lossless** | 26 MB | none |

No Method B parameter was set on HG003. Recall transfers identically; F1 within
0.002.

**Where we lose, stated plainly:** RAM ~2x and wall time ~1.3x. About 3.6 GB of
the 6.99 GB is encoder state Method B never reads — it stays resident because
compression continues after the caller returns. Standalone the caller is ~3.4 GB,
at parity. Method B is also **SNV-only**; superbubble indel emission measured
-0.051 F1 held-out and is gated off.

**Claim 2 is a het-SNV claim** — see `docs/CLAIM2_FINAL.md` for the final
verdict, and `docs/ARCHITECTURE_VS_DISCOSNP.md` for the layer-by-layer
comparison and where the remaining headroom is.

**Indels are a closed loss, not an open task** (`docs/INDEL_BOUND.md`). At
TP=3,290 against 7,768 truth indels, beating DiscoSNP++'s 0.576 requires
FP < 366 — a 79% cut with zero TP loss — and *perfect* precision yields only
0.595. Recall is capped by representation (ALT k-mers absent, 66% homopolymer)
and reference context cannot fix precision (TPs are more homopolymeric than FPs,
0.531 vs 0.427 at n=5,001). Our indel *recall* already beats DiscoSNP++
(0.424 vs 0.393); the gap is precision alone.

## 5. Verification gates — what every change must pass

| gate | requirement |
|---|---|
| Claim 1 | archive **bit-identical** with `CAPS_CALL` unset — every stream, every assembly statistic |
| k-mer identity | `kc nodes` identical across code paths (1,063,607 on r2; 140,719,632 at full chr20) |
| accuracy | window F1 unchanged; held-out windows and a held-out individual for anything claimed |

**k-mer multiset identity is the gate, not F1.** Twice this session a change
passed an F1 check while silently altering `kc` — packing that mapped non-ACGT
to `A` added 2,273 k-mers and moved no calls. F1 is insensitive to divergences
that matter.

## 6. Refuted, with measurements — do not retry

| idea | result |
|---|---|
| graph over finished contigs | 97.7% of anchors have out-degree 1 — chaining linearises the data, which is what makes compression good and what removes the branches |
| chaining-tie harvest | 85 ties, 86% of them duplicate reads, ~12 real against 519 truth variants |
| near-miss (1-mismatch) pairs | +20 rescued, +91 FP: F1 worse |
| bubble chaining | 34-54% coverage, data-dependent |
| unitig extension | F1-neutral; unitigs break at heterozygotes |
| superbubble indels | -0.051 F1 on held-out windows |
| tangle recursion | monotonic loss (0.804 -> 0.760 as tangles allowed) |
| cheap pre-walk reject | exactly zero effect at full scale |
| per-position coverage | inert by construction — our counter only fires after full containment |
| higher `COHC` | 7.6 TPs lost per FP removed |
| histogram-bimodality ploidy test | did not discriminate (haploids fell both sides of the diploid) |

Full detail in `docs/` — each has its own file with the numbers.

## 7. Where the numbers live

| | |
|---|---|
| `results/claim2/t3_t5_full_chr20_HG002.csv` | first full-scale Claim 2 measurement |
| `results/claim2/t3_t5_full_chr20_headtohead.csv` | head-to-head vs DiscoSNP++ |
| `docs/METHOD_B.md` | the graph caller's scope and standing |
| `docs/GATB_DISK_ARCHITECTURE.md` | how DiscoSNP++/GATB counts k-mers, from source |
| `docs/SUPERKMER_PLAN.md` | superkmer spill: plan, execution, verification |
| `docs/PREFLIGHT_CHECKLIST.md` | pre-run audit and the defect it caught |
