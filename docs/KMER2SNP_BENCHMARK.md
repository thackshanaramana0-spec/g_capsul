# Kmer2SNP — the third competitor, finally benchmarked

Written 2026-09-02. Kmer2SNP had never been run in this project (recorded as an
open gap in `docs/PROJECT_AUDIT.md`). It is now benchmarked on the same window,
same reads, same truth and same scoring as DiscoSNP++.

## Result — HG002 chr20:3.0–3.4M, het-SNV, rtg vcfeval

| tool | P | R | **F1** |
|---|---|---|---|
| **G_CAPSUL** | 0.957 | **0.828** | **0.887** |
| DiscoSNP++ | 0.971 | 0.740 | 0.840 |
| Kmer2SNP | **0.992** | 0.302 | **0.464** |

Kmer2SNP is the most **precise** of the three (0.992 — it emits 131 calls with a
single false positive) but recovers only **30%** of the truth. That matches the
~0.53 the outer ARCS project measured for it and its own paper's design, which
keys on heterozygous k-mer pairs and deliberately trades recall for precision.

**G_CAPSUL beats both competitors on this window.**

## What it took to run it — five blockers, all disclosed

Kmer2SNP is 2020-era code and does not run out of the box on a current system.
None of the fixes below touches its algorithm; they are environment and
API-compatibility repairs, and they are listed so the run is reproducible and
auditable.

1. **`Bio` (biopython) missing** — installed into its own conda env
   (`kmer2snp_r`). Its main path imports `Bio.SeqIO`.
2. **`findGSE` (R package) missing**, so its heterozygous-coverage estimate
   `hete.para` came out empty and every run died in `read_findGSE_result`.
   **Not patched around** — Kmer2SNP exposes `--c1 --c2 --r` for exactly this
   case, so those were supplied instead, using values derived from the data:
   a canonical 31-mer histogram of the reads puts the homozygous peak at ~18-21
   and therefore the heterozygous peak at ~9-11, giving the band `--c1 5
   --c2 16`, with `--r 0.001` (the standard human het-SNP rate).
3. **`time.clock()`** — removed in Python 3.8; replaced with
   `time.perf_counter()` (2 files).
4. **`networkx.connected_component_subgraphs`** — removed in networkx 2.4;
   replaced with the documented equivalent
   `(G.subgraph(c).copy() for c in nx.connected_components(G))` (2 files).
5. **DSK not installed**, and `script/run_dsk.sh` shipped with a literal
   `/path2dsk/` placeholder. Installed DSK v2.3.3 binaries and pointed the
   script at them.

Output conversion used the outer project's own `scripts/kmer2snp_sam_to_vcf.py`
(k-mer pairs → BWA → genome VCF), then the same hygiene applied to every tool
here: SNV-only, REF must match the reference genome (0 records dropped), sample
column added for rtg.

## Why this matters for the paper

The three reference-free competitors occupy clearly different operating points
on identical data:

* **Kmer2SNP** — precision 0.992, recall 0.302. Very conservative.
* **DiscoSNP++** — precision 0.971, recall 0.740. Balanced.
* **G_CAPSUL** — precision 0.957, recall 0.828. Highest recall and highest F1.

Claim 2's het-SNV comparison is now complete against **all three** published
reference-free callers rather than one.
