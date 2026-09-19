# PIPELINE.md — G_CAPSUL is one pipeline

> **Status:** frozen 2026-09-10 at tag `v1.0.2-capsule` — code and results final.
> Authoritative numbers live in `../benchmark/results/`; verification status in
> [`../AUDIT.md`](../AUDIT.md). Where this file and a result file disagree, the result file wins.

> **START HERE INSTEAD, unless you specifically want the command reference and the one-pipeline argument.**
> `README.md` for what the project is and what it measured; `CLAUDE.md` for the
> rules and the refuted-ideas list; `benchmark/documentation/RESULT_CODE.md` for
> where any number came from; `benchmark/documentation/REPRODUCE_EVERYTHING.md`
> to re-run it. **This document predates the final 2026-09-09 sweep, so where it
> disagrees with a file in `benchmark/results/`, the result file wins.**

**This file is the command reference.** (It used to call itself the entry
point; `README.md` is.) Claims 1, 2 and 3 are not three tools that
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

## 1a. Every flag that changes a published number — the complete list

Checked exhaustively against the source (150+ `getenv()` flags exist; most
are internal `CAPS_DBG_*` research knobs, self-documented at their point of
use in `scripts/*.sh` and irrelevant to reproducing a published table).
These are the ones that are NOT internal plumbing — get one of these wrong
or omit it, and the tool runs successfully but produces a **different,
silently narrower or differently-configured result**, not an error:

| flag | required for | what happens if you omit it |
|---|---|---|
| `CAPS_QUAL=1` | full lossless round trip (Claim 1) | archive still valid, but cannot reproduce the original quality scores |
| `CAPS_NAMES=1` | full lossless round trip (Claim 1) | archive still valid, but cannot reproduce original read identifiers |
| `CAPS_SPANS=1` | Claim 2/3 prerequisite, and **T3.1's per-contig export specifically** | records contig-span metadata (`g_contig_spans`) that `export`/`coverage`/`query` need; cheap (two push_backs per contig), always safe to set |
| `CAPS_CALL=1` | Claim 2, calling inline during compression | `CAPS_SPANS` is implied by `CAPS_CALL` (`stages/106_inprocess.cpp:337`), but `CAPS_CALL` ALSO runs the full variant caller inline at compress time — measured ~20x heavier on HG002 (913s/13.2GB vs ~250s/~4GB). **If all you need is `contig_spans` for T3.1's export, set `CAPS_SPANS=1` alone — do not set `CAPS_CALL=1` for this, it pays for an unrelated, much heavier pass.** Without either, `capsule` compresses only — no `.vcf`, no per-contig export. |
| `CAPSULE_EXPORT_CONTIGS=1` | **T3.1 correctness (genome fraction, mismatch rate, all QUAST metrics).** Requires `contig_spans` in the archive (see `CAPS_SPANS`/`CAPS_CALL` above). | Set at `export` time, not encode time: `CAPSULE_EXPORT_CONTIGS=1 capsule_decode export <archive> <out.fa>`. Without it, `export` concatenates every contig into 2 records (the main/second pseudogenome regions) instead of writing each of the archive's real contigs separately. **This is not cosmetic — on E. coli it was 41,603 real contigs concatenated into 2, and scoring the concatenated form against a reference gives 80.8% genome fraction / 90% unaligned; scoring the true per-contig form gives 98.7% genome fraction, beating SPAdes.** Every T3.1 correctness number in `docs/T3.1_CORRECTNESS_FINAL_20260919.md` requires this flag. The two-record form is still used by the internal ploidy-gate path (`CAPS_DBG_ONLY`) deliberately — that path wants the concatenated form, not a bug. |
| `CAPS_CALL_INDELS=1` | **T2.3 (het-indel), T2.4 (multi-allelic), T2.5 (tetraploid)** | `capsule_decode call` silently takes the graph-only SNV path instead of erroring — a different, narrower configuration, not a failure. Every published indel/multi-allelic/tetraploid number depends on this being set. |
| `CAPS_PILEUP=1` | T3.4 (locus fidelity), the sidecar built by `capsule_decode index` | without it, `query`'s coordinate arm returns the pseudogenome **consensus** instead of each read's own deviations — every het site reads as 0/0 by construction, not by measurement |
| `CAPS_QUERY_MM=k` | **not required for any published number.** Claim 3 work only, see `docs/CLAIM3_LOCUS_ADDRESSABILITY.md` | defaults to 0, which is byte-identical to the exact `std::string::find` used for the published 345/400. Set k>0 and a sequence `query` tolerates k mismatches via pigeonhole seed-and-extend, resolving the alternate-haplotype contig at a het locus. **Capped by the seed floor at `k_max = floor(P/12) - 1`**, so a 40 bp probe cannot exceed k=2 whatever value is set. Do not set it when reproducing T3.4. |
| `CAPS_QUERY_MINSEED=n` | **not required for any published number.** Diagnostic only | defaults to 12, the pigeonhole seed-length floor. It is what bounds the tolerance: `k_max = floor(probe_len / MINSEED) - 1`, so a 40 bp probe caps at k=2 however large `CAPS_QUERY_MM` is set. Lowering it raises the cap at the cost of more candidates to verify. Values below 4 are ignored. **Measured to change nothing on exact-match recall** (0.979 at MINSEED 12/8/6/5), which is how that hypothesis was refuted — see `docs/CLAIM3_LOCUS_ADDRESSABILITY.md`. |
| `CAPS_XMI=1` | optional, at `index` time. Requires `CAPS_PILEUP=1` | builds `<sidecar>.xmi`, a k-mer index over the deviation-carrying reads, so `query` can also answer "which reads CONTAIN this probe" completely. Only deviation-carrying reads can be missed by the ordinary search, because a read with none is emitted as pure pseudogenome. `CAPS_XMI_STRIDE=n` (default 8) trades index size for the minimum probe length the guarantee holds at, `k_max + stride - 1`. Unset, no index is written and the sidecar is byte-identical. |
| `CAPS_QUERY_CONTAIN=1` | optional, at `query` time. Requires a `.xmi` | returns the UNION of "reads covering the locus" and "reads containing the probe". Refuses loudly, rather than degrading silently, if the probe is shorter than the index guarantee allows. Unset, `query` output is unchanged. **Not optional for T3.4/T3.5's published multi-individual numbers** — plain `query` alone reached 1.0000/100% on HG002 by chance (no test read exceeded the tolerance ceiling there); on HG003/HG004/HG005 plain `query` only reaches 0.86-0.96 / leaves one T3.5 site unresolved. `CAPS_QUERY_CONTAIN=1` is required to reach 1.0000 and 100% in general — see `docs/T34_T35_MULTI_INDIVIDUAL_20260919.md`. Raises false positives on homozygous controls by ~35-75% relative — a real, disclosed cost, not a free win. |
| `CAPS_PLOIDY=N` | Claim 2 on non-diploid input (T2.5 uses 4) | defaults to 2; wrong ploidy silently mis-scores heterozygous sites |
| `GSEARCH=1` | optional, smaller archives | selects golden-section search over MAXMAP instead of the default 4-point grid (~7-10 probes vs 4, -154,223 B over 7 files measured) |
| `MAXMAP=N` | optional, manual override | bypasses the automatic candidate sweep entirely; used only for targeted debugging, never for a published number. **Do not use to chase T3.1 assembly-quality metrics** — tested this session (`docs/T3.1_CORRECTNESS_FINAL_20260919.md`): raising it to 90 costs ~11.6% archive size for a genome-fraction result that is actually slightly *worse* than the default (`MAXMAP=30`, `MINOV=52`) configuration already gives. The archive-optimal setting is also the correctness-optimal one; do not retune for T3.1. |

`capsule_decode query` also accepts **comma-separated probes** in one
invocation (`query <archive> <out.fa> "SEQ1,SEQ2"`). The sidecar load, placement
decode and haystack build then happen once rather than per probe — measured
1.53x faster for two probes at identical peak memory, with occurrence sets
verified identical to separate invocations. Each `[query] occ A B ± i` line
carries the index of the probe that produced it. A single probe behaves exactly
as before.
| `ARCS_AUTOCHUNK_MB=N` | large inputs (>2 GB) on a memory-constrained box | suppresses/adjusts the automatic chunking threshold; the 19-dataset sweep did not need this (12 vCPU / ~90 GB box), but a smaller box attempting C. elegans or T. cacao-scale input will |

Every other environment variable referenced anywhere in `scripts/*.sh` is
either (a) harness plumbing that a fresh run of the actual benchmark script
inherits automatically — `DUMP_LIT`, `DUMP_PERM`, `DUMP_MM`, `GRID4`,
`GRID8`, `CAPS_SKIP_SCOPE_CHECK`, `CAPS_ENCODER_PATH` — or (b) a
`CAPS_DBG_*` research/ablation knob never used to produce a published
number. If you are hand-typing a command instead of running the provided
script, the table above is everything that can silently change your answer.

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

# Claim 2 calling FROM AN EXISTING ARCHIVE (no FASTQ read) -- what T2.1/T2.3/
# T2.4/T2.5 actually run. CAPS_CALL_INDELS=1 is REQUIRED: without it the
# decoder's `call` subcommand takes the graph-only (SNV-only) path silently --
# not an error, just a different, narrower configuration -- and every
# published het-indel, multi-allelic and tetraploid number depends on it.
CAPS_CALL_INDELS=1 /tmp/capsule_decode call archive.capsule out.vcf callwk/
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
| sequence layer | 26 MB | 26 MB | none |
| **full archive (lossless)** | **547 MB** | — | none |

No Method B parameter was set on HG003. Recall transfers identically; F1 within
0.002.

**Read the two size rows carefully — they measure different things.** An earlier
version of this table listed a single "archive 26 MB lossless" row, which was
wrong twice over: 26 MB is the SEQUENCE LAYER only (literal + mem_triples +
mm_*), and a 26 MB file cannot be lossless because it carries no names and no
quality. Measured on HG002 (4.28 GB input, full FASTQ, CAPS_NAMES + CAPS_QUAL):

| stream | size | share |
|---|---|---|
| qual_body | 470.6 MB | **82.0%** |
| pos_abs | 43.4 MB | 7.6% |
| names_body | 30.6 MB | 5.3% |
| literal | 14.5 MB | 2.5% |
| mm_pos, mem_triples, mm_sym, extmm | ~11 MB | 1.9% |
| **total** | **547 MB** | 13.41% of input |

Quality is 82% of a real archive, so any statement about "our archive size" that
omits it is off by a factor of twenty. The comparable, measured numbers on the
same file are SPRING 570 MB (13.98%) and Genozip 908 MB (22.24%) -- we win both,
by 4.0% and 39.7%, verified LOSSLESS.

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
| `docs/_removable/SUPERKMER_PLAN.md` | superkmer spill: plan, execution, verification |
| `docs/PREFLIGHT_CHECKLIST.md` | pre-run audit and the defect it caught |
