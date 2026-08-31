# CLAUDE.md — ARCS Benchmark Project

This file governs all AI-assisted work in this repository. Read it fully before touching any script, dataset, or result.

---

## Project Identity

**ARCS** = Assemble → Retain → Compress → Serve  
**Paper title:** "ARCS: a retained assembly for unified lossless FASTQ compression and reference-free genomic analysis"  
**Target journal:** Nature Methods  
**Server:** Ubuntu 24.04, SDC3 Chennai, 12 vCPU, ~90 GB RAM, 250 GB SSD  
**Repo path on server:** `/root/arcs-clean`

---

## Absolute Rules — Never Break These

1. **Do not change the 10 benchmark accessions.** They are locked in `benchmark/DATASET_LOCKED.md`. Do not substitute, re-verify, or re-derive them. If a download fails, report which accession failed and stop.
2. **ARCS compression uses zero flags.** Always: `arcs compress INPUT.fq OUTPUT.arcs` — nothing else. No `--fast`, no `ARCS_PAR_SHARDS`, no `ARCS_CHUNK_THREADS`. Auto-chunk for >2 GB inputs is expected default behaviour, not a problem.
3. **Never run two timed benchmark jobs concurrently.** One tool at a time. Concurrent jobs contaminate timing and RAM measurements.
4. **Never silently alter methodology.** If something fails, stop, print the error, and wait for instructions.
5. **ARCS must be fully lossless before Claim 2 starts.** If any `TOOL=ARCS ... lossless=LOSSY` appears in Claim 1 output, halt immediately. SPRING `lossless=LOSSY` is expected (it strips `+` line IDs, which are spec-redundant); the lossless check normalizes the `+` line on both sides so SPRING comparison is fair on data recovery.
6. **Projected numbers are sanity checks only.** Fresh server measurements are authoritative. Do not reject a server result because it differs from a projection.
7. **ARCS must pass 7/7 ctests before any benchmark phase runs.**

---

## The 10 Locked Datasets

Source of truth: `benchmark/DATASET_LOCKED.md`. Reproduced here for fast reference.

**10 full SRA datasets (primary benchmark):**

| # | Accession | Organism | Kingdom | _1 size (uncompressed) | Peak RAM | Auto-chunk |
|---|-----------|----------|---------|------------------------|---------|-----------|
| 1 | SRR2584863 | E. coli B REL606 | Bacteria | ~576 MB | ~5 GB | No |
| 2 | ERR552797 | M. tuberculosis H37Rv | Bacteria | ~217 MB | ~3 GB | No |
| 3 | SRR554369 | P. aeruginosa PAO1 | Bacteria | ~334 MB | ~4 GB | No |
| 4 | ERR5181310 | SARS-CoV-2 | Virus | ~30 MB | ~1 GB | No |
| 5 | ERR17740259 | S. aureus WGS | Bacteria | ~970 MB | ~6 GB | No |
| 6 | SRR065390 | C. elegans N2 WGS | Animalia | ~11 GB | ~18 GB | Suppressed† |
| 7 | DRR976266 | S. cerevisiae | Fungi | ~1.67 GB | ~8 GB | No |
| 8 | SRR870667 | Theobroma cacao WGS | Plantae | ~15 GB | ~28 GB | Suppressed† |
| 9 | SRR36741279 | Leishmania major | Protista | ~1.70 GB | ~7 GB | No |
| 10 | SRR37283774 | P. falciparum | Protista | ~669 MB | ~5 GB | No |

**† `ARCS_AUTOCHUNK_MB=25000` is set in run_block1.sh — an env var, not a command-line flag. Ensures single-pass assembly on the 90 GB server.**

**Claim 2 only — HG002-HG005 chr20 at standardized 30× (NOT in Claim 1):**

| # | File | Individual | Standardized depth |
|---|------|-----------|-------------------|
| C2-1 | HG002_pooled.fq | NA24385 Ashkenazi son | 30× chr20 |
| C2-2 | HG003_pooled.fq | NA24149 Ashkenazi father | 30× chr20 |
| C2-3 | HG004_pooled.fq | NA24143 Ashkenazi mother | 30× chr20 |
| C2-4 | HG005_pooled.fq | NA24631 Han Chinese son | 30× chr20 |

Sourced from GIAB S3 WGS BAMs (chr20 stream via samtools), downsampled to 30× for fair T3 comparison. GIAB S3 source coverage varies per individual (60-300×); standardization removes this confound.

**Slots 6 and 8 exceed 2 GB; auto-chunk is suppressed via `ARCS_AUTOCHUNK_MB=25000` in run_block1.sh.**

**Banned accessions (never use):** SRR390728, SRR988075, SRR327342, SRR1663585, SRR1296601, ERR015526, SRR1294122, ERR174310, SRR16357346, SRR1945765

**Disk budget:** ~82 GB peak on 250 GB SSD — safe.

**Cross-validation:** Server SPRING bpb must match published gold values ±2%:
- SRR554369 = 0.2416 bpb (P. aeruginosa)
- SRR870667 = 1.2621 bpb (T. cacao — SPRING's worst; ARCS expected to win largest margin here)

If either deviates >2%, SPRING benchmark setup is wrong — stop and investigate.

---

## 3 Claims and 6 Tables

### Claim 1 — COMPACT (T1 + T2)
- **T1:** Archive size (bytes) — ARCS vs SPRING vs Genozip, all 10 datasets
- **T2:** Compress time (s) + decompress time (s) + peak RAM (VmHWM KB), all tools, all datasets
- ARCS expected to win ratio on all 10; speed slower than SPRING/Genozip (assembly cost is expected)

### Claim 2 — FAITHFUL (T3 + T4 + T5)
- **T3:** Het-SNV F1 — ARCS vs DiscoSNP++ vs Kmer2SNP on HG002/HG003/HG004/HG005 chr20
- **T4:** Coverage sweep — ARCS F1 at 10×/15×/20×/30× (HG002)
- **T5:** Het-indel F1 — ARCS vs DiscoSNP++ on HG002 chr20
- Expected: ARCS avg F1 ≈ 0.936 > DiscoSNP++ 0.918 > Kmer2SNP 0.532

### Claim 3 — ADDRESSABLE (T6)
- **T6:** `arcs export` vs SPAdes, `arcs coverage` vs BWA+mosdepth, `arcs query` (unique)
- Expected: export ~50-200× speedup vs SPAdes, coverage ~2-5× speedup vs BWA+mosdepth

---

## Phase Structure

Run phases strictly in order. Do not proceed to the next phase if the current one has failures.

### Phase 0 — Build and test (5-10 min)

```bash
cd /root/arcs-clean
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc)
ctest --test-dir build
```

**Pass condition:** `7/7 tests passed`  
**Stop if:** any test fails — do not proceed

```bash
# Tool check
for t in prefetch fasterq-dump aws spring genozip; do
    command -v "$t" && echo "$t OK" || echo "$t MISSING"
done
```

**Stop if:** any tool is MISSING — install before continuing

### Phase 1 — Download (2-4 hours)

```bash
bash benchmark/download.sh /data/fastq 2>&1 | tee /data/download.log
tail -30 /data/download.log
```

**Pass condition:** Last line contains `DOWNLOAD COMPLETE`  
**Stop if:** any line contains `FAIL` — paste the failing phase output for diagnosis

**Large file notes:**
- SRR065390 (C. elegans, ~22 GB SRA) and SRR870667 (T. cacao, ~15 GB SRA) are the slow downloads
- **AWS speed-up for SRA files:** download from S3 Open Data mirror instead of NCBI prefetch:
  ```bash
  aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR065390/SRR065390 /data/fastq/prefetch/SRR065390/SRR065390.sra
  aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR870667/SRR870667 /data/fastq/prefetch/SRR870667/SRR870667.sra
  ```
  Then re-run download.sh (it will skip prefetch and proceed to fasterq-dump)
- Phase 5 of download.sh also fetches Claim 2 prerequisites:
  - `~/refs/chr20.fa` — GRCh37 chr20 reference (~65 MB, chromosome named "20")
  - `~/giab_truth/` — GIAB truth VCFs v4.2.1 for HG002-HG005 (~2.8 GB total)

### Phase 2 — Claim 1: Compression benchmark (2-4 hours)

```bash
bash benchmark/benchmark.sh claim1 /data/fastq /root/arcs-clean/build/arcs ./results 2>&1 | tee ./results/claim1.log
```

After completion, run both checks:
```bash
grep -h "TOOL=ARCS" ./results/claim1/*.log
grep "TOOL=ARCS.*lossless=LOSSY" ./results/block1/*.log && echo "BUG — ARCS LOSSY" || echo "ARCS ALL LOSSLESS"
```

**Pass condition:** Every log line ends with `lossless=LOSSLESS`  
**Stop if:** any `LOSSY` — this is a bug, do not proceed to Claim 2  
**Cross-validate:** SPRING bpb for SRR554369 and SRR870667 must match published values ±2%

### Phase 3 — Claim 2: Variant calling (2-3 hours)

```bash
export CHR20_FA=~/refs/chr20.fa
export CHR20_SDF=~/refs/chr20.sdf
export GIAB_TRUTH=~/giab_truth
export DISCO_DIR=~/DiscoSnp
export CONDA_ENV=kmer2snp_r

bash benchmark/benchmark.sh claim2 /data/fastq /root/arcs-clean/build/arcs ./results 2>&1 | tee ./results/claim2.log
cat ./results/claim2/t3_snv_f1.csv
```

**Pass condition:** T3 CSV has F1 values for all 4 GIAB individuals  
**Expected range:** ARCS F1 0.90–0.95, DiscoSNP++ F1 0.88–0.94, Kmer2SNP F1 0.48–0.56  
**Stop if:** ARCS F1 < 0.85 — check VCF liftover and truth region

### Phase 4 — Claim 3: Archive analysis (1-2 hours)

```bash
bash benchmark/benchmark.sh claim3 /data/fastq /root/arcs-clean/build/arcs ./results 2>&1 | tee ./results/claim3.log
cat ./results/claim3/t6_results.csv
```

**Pass condition:** T6 CSV shows export speedup ≥ 40× for E. coli vs SPAdes  
**Stop if:** `arcs export` output is empty FASTA or `arcs coverage` TSV has zero rows

---

## Critical Command Syntax

```bash
# ARCS — zero flags always
arcs compress reads.fq out.arcs                          # compress
arcs decompress out.arcs decoded.fq                      # decompress
arcs compress --call reads.fq out.arcs out.vcf           # compress + call variants
arcs export out.arcs contigs.fa                          # export pseudogenome
arcs coverage out.arcs coverage.tsv                      # per-contig coverage
arcs query out.arcs 0-100000 region.fq                  # NOTE: hyphen, not colon

# SPRING — always needs -g for FASTQ output on decompress
spring -c -i in.fq -o out.spring -t $(nproc) -g         # compress
spring -d -i out.spring -o decoded.fq -t $(nproc) -g    # decompress — -g required

# Genozip
genozip --force -o out.genozip in.fq
genounzip --force -o out.fq in.genozip

# DiscoSNP++ — -G flag is mandatory
run_discoSnp++.sh ... -G chr20.fa                        # -G required, no exceptions

# Lossless check — order-free, CRLF-safe
paste - - - - < orig.fq | tr -d '\r' | sort > /tmp/a
paste - - - - < dec.fq  | tr -d '\r' | sort > /tmp/b
cmp -s /tmp/a /tmp/b && echo LOSSLESS || echo LOSSY
```

---

## What to Paste Back for Debugging

| Situation | Paste this |
|-----------|-----------|
| Phase 0 test failure | Full `ctest` output |
| Download failure | Lines around `FAIL` in download.log |
| LOSSY result in Claim 1 | The full `__ARCS.log` for the failing dataset |
| ARCS F1 < 0.85 in Claim 2 | Last 50 lines of claim2.log + t3_snv_f1.csv |
| SPAdes failure in Claim 3 | Last 30 lines of spades.log |
| Any unexpected output | The exact phase line that failed + 10 lines above it |

---

## Correct Binary Path

```bash
# CORRECT
/root/arcs-clean/build/arcs

# WRONG — old binary, do not use
/usr/local/bin/arcs
```

Always build from source on the server. Never use a system-installed arcs binary.
