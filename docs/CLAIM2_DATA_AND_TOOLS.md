# Claim 2 — data, tools, and files inventory (verified 2026-09-02)

Everything needed to score G_CAPSUL's caller against real GIAB truth is
**already on this server** — nothing was missing, nothing downloaded today.
This file traces exactly what exists, where, and which new file in this repo
uses it. Not run yet (per instruction) — this is the inventory + wiring, not
a benchmark result.

## 1. Data — all present, verified by direct `ls`

| what | path | status |
|---|---|---|
| GRCh37 chr20 reference | `~/refs/chr20.fa` (64,286,035 B) | present |
| chr20 rtg SDF (pre-built) | `~/refs/chr20.sdf/` | present |
| GIAB v4.2.1 truth VCF, HG002 | `~/giab_truth/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` (+`.tbi`) | present |
| GIAB v4.2.1 truth VCF, HG003 | `~/giab_truth/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` (+`.tbi`) | present |
| GIAB v4.2.1 truth VCF, HG004 | `~/giab_truth/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` (+`.tbi`) | present |
| GIAB v4.2.1 truth VCF, HG005 | `~/giab_truth/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` (+`.tbi`) | present |
| GIAB confident-region BEDs, all 4 individuals | `~/giab_truth/HG00{2,3,4,5}_..._benchmark*.bed` | present |
| Standardized 30× chr20 FASTQ, HG002 | `/data/fastq/HG002_pooled.fq` | present |
| Standardized 30× chr20 FASTQ, HG003 | `/data/fastq/HG003_pooled.fq` | present |
| Standardized 30× chr20 FASTQ, HG004 | `/data/fastq/HG004_pooled.fq` | present |
| Standardized 30× chr20 FASTQ, HG005 | `/data/fastq/HG005_pooled.fq` | present |

**Verdict: Claim 2 data is 100% complete. No download needed.**

## 2. Tools — all present, verified by direct `command -v` / `ls`

| tool | status | used for |
|---|---|---|
| `bwa` | OK | contig→genome coordinate lift (eval-only) |
| `rtg` (rtg-tools) | OK | `rtg vcfeval` gold-standard scoring |
| `bcftools` | OK | left-align/normalize before scoring |
| `tabix` / `bgzip` | OK | indexed VCF handling, remote truth fetch |
| DiscoSNP++ | built at `~/DiscoSnp/build/` | competitor comparison |
| Kmer2SNP | present at `/root/Kmer2SNP/kmer2snp.py`, conda env `kmer2snp_r` at `/root/miniconda3/envs/kmer2snp_r` | competitor comparison |

**Verdict: every tool the outer project's Claim 2 pipeline needs is already
installed. Nothing to install.**

## 3. Files added to THIS repo (`c_star_pg_advance`) today, and what each does

Claim 2 previously existed only in the outer `/root/arcs-clean` project. This
repo now has its own complete, self-contained copy of the pipeline, wired to
G_CAPSUL's own caller instead of ARCS's:

| file | role | source |
|---|---|---|
| `include/caps_caller.h` | the caller itself (pileup/SNV/indel/polyploid, frozen params) | ported from `/root/arcs-clean/src/caller.cpp` |
| `stages/106_inprocess.cpp` (`CAPS_CALL=1` gate) | contig-span capture + `CallData` population + caller invocation | new, written for G_CAPSUL's own assembly (see `docs/_removable/CLAIM2_BUILDER.md` for why this differs structurally from ARCS's `vodbg_pg`) |
| `scripts/lift_vcf.py` | contig→genome coordinate lift (`hapflank_lift`), tool-agnostic | copied unchanged from outer project (no ARCS-specific coupling) |
| `scripts/eval_caller.py` | quick window-based P/R/F1 scorer (secondary to rtg) | copied unchanged |
| `scripts/extract_vcfeval_metrics.py` | parses `rtg vcfeval`'s `summary.txt` | copied unchanged |
| `scripts/sim_indel_bench.py`, `scripts/sim_indel.py` | synthetic diploid indel-truth generators | copied unchanged |
| `scripts/sim_polyploid.py` | synthetic k-ploid truth generator | copied unchanged |
| `scripts/run_indel_bench_capsule.sh` | synthetic indel pipeline: **G_CAPSUL call → BWA → lift → rtg** | adapted from outer `run_indel_bench.sh` — only the caller-invocation line changed (`CAPS_CALL=1 CALL_VCF=... CAPS_DUMP_CONTIGS=... "$CAPS" reads.fq <args>` instead of `"$ARCS" call reads.fq calls.vcf`) |
| `scripts/run_giab_indel_capsule.sh` | **real-GIAB** HG002 chr20:2.0–2.4M SNV+indel pipeline, DiscoSNP++-comparable | adapted the same way from outer `run_giab_indel.sh` |
| `scripts/run_polyploid_bench_capsule.sh` | polyploid (`CAPS_PLOIDY=k`) pipeline | adapted the same way from outer `run_polyploid_bench.sh` |

**Not copied (deliberately):** `scripts/run_giab_indel.sh`,
`scripts/run_indel_bench.sh`, `scripts/run_polyploid_bench.sh` themselves
(the originals) — no need to duplicate them un-adapted; the `_capsule.sh`
versions above are the ones this repo actually uses. `scripts/kmer2snp_sam_to_vcf.py`
was not copied either — it's a shared, tool-agnostic VCF converter already
usable in place from the outer project's `scripts/` if a Kmer2SNP comparison
is run from here.

## 4. One-line verdict — the 14-dataset Phase1/Phase2b/Phase3 benchmark

**Complete: 14/14 datasets finished, all round-trip-verified LOSSLESS, whole-file
Phase 3 wins 14/14 against SPRING and 14/14 against Genozip** — full numbers
in `docs/SOTA_COMPARISON.md` (12-row table there predates the final 2; the
authoritative source is `/tmp/allph/r.csv`, all 14 rows `p1ok/p2ok/p3ok=OK`).
The 15th dataset (Utricularia gibba, `SRR10676752`) is still queued for
download — not part of this run.

## 5. What running these scripts would need (for later, not now)

- `scripts/run_giab_indel_capsule.sh /tmp/caps_call_enc scripts /data/fastq/HG002_pooled.fq ~/refs/chr20.fa`
  — real-GIAB SNV+indel F1, directly comparable to the outer project's own
  0.936 (SNV) / 0.505 (indel) numbers.
- `scripts/run_indel_bench_capsule.sh <workdir> /tmp/caps_call_enc scripts`
  (after `python3 scripts/sim_indel_bench.py <workdir> 60000 50 12345`) —
  synthetic sanity check (outer project's own equivalent scores SNV F1=1.000).
- `scripts/run_polyploid_bench_capsule.sh <workdir> /tmp/caps_call_enc scripts 3`
  (after `python3 scripts/sim_polyploid.py <workdir> 3 40000 60 42`).

All three are held for a future "run it" instruction — this pass only
confirms the data/tools exist and wires the G_CAPSUL-specific scripts into
this repo, per the user's explicit "no need to run" instruction.
