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

## Where it stands against the field — honestly, it loses

Same metric, same engine, same windows (outer ARCS's published numbers):

| tool | r2 SNV F1 | r3 SNV F1 |
|---|---|---|
| ARCS (outer project) | 0.923 | **0.954** |
| DiscoSNP++ | 0.826 | 0.886 |
| Kmer2SNP | 0.550 | 0.544 |
| **CAPSULE (this repo)** | **0.419** | **0.394** |

CAPSULE currently places **last**, behind even Kmer2SNP. This is not a
tuning miss and not a scoring artifact — see the cause below.

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

## What would actually fix it (not yet attempted)

In rough order of expected value per unit of work:

1. **Mismatch-tolerant chain extension in calling mode.** Allow a chain to
   extend across ≤1 mismatch when `CAPS_CALL` is set. This is the direct
   analogue of `MEM_MAXMM` (already implemented for MEM references, shipped
   off because it costs *compression* bytes) — but for calling, longer
   contigs are worth far more than the bytes cost, and calling mode does not
   need to preserve the archive. Would merge both haplotypes onto one
   contig, restoring pileup depth.
2. **Post-assembly contig merge for calling only** — join contigs whose ends
   overlap within ≤1–2 mismatches, before the caller runs. Cheaper to
   implement, does not touch the compression path at all.
3. **Accept the split and improve bubble detection** — shorter anchors than
   25-mers, or allow anchors occurring >2 times. Weakest option: it fights
   the fragmentation rather than fixing it, and short anchors cost precision.

**Do not** respond to these numbers by loosening the frozen filter
parameters. Precision 0.94 with recall 0.25 is a *visibility* problem, not a
threshold problem — the variants are not being presented to the filters at
all.

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
