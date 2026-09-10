# RESULTS — every table, as executed

> **Status:** frozen 2026-09-10 at tag `v1.0-capsule` — code and results final.
> Authoritative numbers live in `../benchmark/results/`; verification status in
> [`../AUDIT.md`](../AUDIT.md). Where this file and a result file disagree, the result file wins.

All numbers come from one run: `benchmark_1_run.sh` at commit `21ee619`,
2026-09-09 02:46:07 → 07:53:21 (5 h 07 m), 19 datasets, 0 failures, one timed
job at a time on an idle machine. Raw CSVs: `../benchmark/results/`.
Traceability, table by table: `../benchmark/documentation/RESULT_CODE.md`.

**Where this file and a CSV disagree, the CSV wins.**

---

## Claim 1 — COMPACT

**19/19 wins against SPRING. 19/19 against Genozip. 57/57 archives LOSSLESS.**

| | bytes | vs ours |
|---|---:|---:|
| **G_CAPSUL** | **6,011,370,147** | — |
| SPRING | 6,396,897,280 | **−6.03%** |
| Genozip | 10,594,522,604 | **−43.26%** |

19 datasets, 55.0 GB of FASTQ, spanning bacteria, archaea, viruses, fungi,
protists, plants, animals and human. Losslessness is verified by decoding the
**archive** — not the encoder's intermediate dumps, which bypass the entropy
layer where one whole class of bug lived.

Speed and memory, stated plainly: **2.99× SPRING on compress, 1.79× on
decompress.** We are faster than SPRING on 3 of 19 and lighter on 9 of 19.
Assembly costs time; that is the trade, and it is not hidden.

### Against PgRC2 — a different measurement, deliberately kept separate

**+1.88%** smaller on the DNA stream. It is not in the CSV above and must not be
mixed with it, because **PgRC2 is not a FASTQ compressor** — verified by running
it, not inferred: given a 440,048 B FASTQ it returns 151,000 B of bare DNA, one
line per read. No names, no `+` line, no quality. It cannot reproduce the input.

| | DNA | names | line 3 | quality | round-trips a FASTQ |
|---|---|---|---|---|---|
| PgRC2 | yes | no | no | no | **no** |
| SPRING | yes | yes | yes | yes | yes |
| Genozip | yes | yes | yes | yes | yes |
| **ours** | yes | yes | yes | yes | **yes** |

Three separate statements, all true: we beat the DNA specialist on the column it
specialises in; we are a complete lossless FASTQ compressor; and we beat both
complete compressors on size.

---

## Claim 2 — FAITHFUL

**Every number here is called FROM THE ARCHIVE** — `capsule_decode call
<archive> <out.vcf>`, no FASTQ, no reference, no separate assembly graph.
Competitors get identical reads, truth, confident regions, normalisation order
and scorer (`rtg vcfeval --squash-ploidy`, GIAB v4.2.1, het-restricted). Only
the caller differs.

### T2.1 — het-SNV F1

| individual | ours | DiscoSNP++ | Kmer2SNP |
|---|---:|---:|---:|
| HG002 | **0.888** | 0.847 | 0.464 |
| HG003 | **0.891** | 0.849 | 0.460 |
| HG004 | **0.891** | 0.851 | 0.469 |
| HG005 | 0.834 | **0.863** | 0.505 |
| **mean** | **0.876** | 0.853 | 0.475 |

**3 wins, 1 loss** — and the loss is published rather than averaged away. See
Limitations for why HG005 loses and what it costs to say so.

### T2.2 — coverage sweep (HG002)

| depth | 10× | 15× | 20× | 30× |
|---|---:|---:|---:|---:|
| F1 | 0.487 | 0.704 | 0.797 | **0.888** |

Monotone. Each depth is a full re-compress and re-call, so depth is the only
variable.

### T2.3 — het-indel F1

| individual | ours | DiscoSNP++ |
|---|---:|---:|
| HG002 | **0.620** | 0.576 |
| HG003 | **0.637** | 0.587 |
| HG004 | **0.632** | 0.595 |
| HG005 | 0.593 | **0.605** |
| **mean** | **0.621** | 0.591 |

**3 wins, 1 loss.** Kmer2SNP records `NOT_APPLICABLE_SNP_ONLY` — it has no indel
model, so the blank is a property of the tool, not a gap in the benchmark.

> **This result reverses an earlier conclusion of our own, and the reversal is
> instructive.** `docs/INDEL_BOUND.md` and `docs/CLAIM2_FINAL.md` withdrew the
> indel claim as a closed loss (0.516 vs 0.576) and proved, with the bound
> `F1 = 2·TP/(TP+FP+truth)`, that **no filter** could flip it at TP = 3,290.
> That proof was correct, and it named the way out: *a different candidate
> generator, not another filter.* The generator changed, and both halves of the
> arithmetic moved at once — **TP 3,290 → 3,871 and FP 1,711 → 824** — which is
> exactly what a filter cannot do. The bound itself still holds and still
> reproduces the CSV.

### T2.4 — multi-allelic sites

**17 / 26 (65.4%)** — and **26 is the complete census** for chr20, every
SNV-only multi-allelic site there is, not a sample.

> **Corrected 2026-09-10.** Published earlier as 21/26. Withdrawn during a
> reproducibility audit: no raw log on this server supported that figure, and
> three independent re-derivations (two encoder builds, one exact repeat)
> each gave 17/26 with identical intermediate counts. The denominator (26)
> is independently confirmed and unaffected.
> `benchmark/documentation/T2.4_reproduction_20260910/README.md`.

**DiscoSNP++ recovers 0 of 26**, because across its entire 3,989-record output
it emits no record with more than one ALT allele. That is a structural inability,
not a miss.

### T2.5 — tetraploid (real reads, not simulated)

HG003+HG004 concatenated, ploidy 4, truth = the union of their GIAB calls.

| class | **ours** | DiscoSNP++ |
|---|---:|---:|
| SNV | **0.897** | 0.782 |
| INDEL | **0.597** | 0.553 |

---

## Claim 3 — ADDRESSABLE

### T3.1 / T3.2 / T3.3 — export, coverage, query

| operation | vs | speedup |
|---|---|---:|
| `export` | SPAdes 4.0.0 | **129–784×** |
| `coverage` | bwa + samtools sort + index + mosdepth, timed as one pipeline | **16–54×** |
| `query` | — | **not a speed claim** (see below) |

**The export ratio carries a caveat, which is why the CSV records output bytes
and row counts for both sides.** Our export emits 2 records (the pseudogenome);
SPAdes emits tens of thousands of biological contigs. The ratio measures
time-to-a-reference-free-coordinate-system from an archive that had to exist
anyway — **not "the same output, faster."** T3.2 and T3.3 carry no such caveat.

**`query` is not faster than the alternative and we say so.** `genocat
--head=100` extracts in 0.19 s against our 0.46 s. With the optional sidecar
index we reach 0.03 s, but Genozip is the honest comparator and Claim 3 does not
lead with speed.

### T3.4 — locus retrieval fidelity, the mechanism table

400 GIAB het SNV sites, 4 individuals, **same archive, same sites, both
addressing modes** — an internal control, because no other tool can produce a row.

| addressing mode | both alleles returned |
|---|---:|
| by coordinate | **81 / 400** (20.2%) |
| by content (ours) | **345 / 400** (86.2%) |
| pooled allele balance | 1.00 |

**Probe design matters and is stated:** 40 bp of *reference* ending 6 bp
**before** the variant, so the probe never contains the variant. A probe taken
from one haplotype's own sequence could only match that haplotype and would rig
the result.

> **A correction against ourselves that must not be lost.** This table was
> published with coordinate at **0/400**. That was an artifact of our own query
> emitting the *consensus* rather than the reads, so no allele difference could
> ever appear and coordinate addressing looked absolutely incapable. The query
> now applies each read's deviations (commit `d58fc23`). The honest figure is a
> **4.3× gap, not an infinite one.** Any document still saying 0/400 is
> superseded.

**Cost of all addressability:** the `contig_spans` stream, 232,509 B — **0.041%**
of a 573 MB archive.

---

## Reproducibility

Verified 2026-09-10 by re-running, not assumed.

**Byte-identical regardless of machine configuration:**

| knob | tested |
|---|---|
| core count | 2, 3, 7, 12 → identical archive, identical VCF |
| candidate concurrency | `CAPS_CAND_SEQ=1`, `CAPS_CAND_PAR=4` → identical archive |
| caller memory budget | 5,000 MB default vs `CAPS_MAXRAM_MB=2000` → identical VCF |

**The one number to check on any machine:**

    CLAIMS=1 bash benchmark/scripts/sanity_archive_one.sh SRR2584863
    →  archive 68,429,027 B exactly · ratio 9.844% · lossless 3/3

If that differs on any machine of any size, it is a real regression, not an
environment difference.

**What legitimately varies:** wall times and peak RSS (machine speed and
allocator), and Genozip's archive size, which embeds run metadata and drifts
~37 B (0.005%) run to run — never gate on it.
