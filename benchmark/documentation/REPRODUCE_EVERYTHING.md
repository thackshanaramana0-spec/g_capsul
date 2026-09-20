# REPRODUCE_EVERYTHING — a reviewer's complete path

If you want to verify every claim from nothing, this is the whole procedure.
Nothing here is aspirational: each step has been executed on the machine that
produced the results.

## 0. What you need

| | requirement | why this number |
|---|---|---|
| cores | 8+ (12 used) | fewer works; **output is byte-identical**, only time changes |
| RAM | 32 GB min, 64 GB comfortable | measured peaks: 10.9 GB compress, 19.9 GB caller |
| disk | **250 GB** | 55 GB data + 7 GB archives + ~36 GB peak transient |
| OS | Ubuntu 24.04 | the SPRING patch is GCC-13-specific |
| time | ~2 h setup, ~4 h download, **5 h benchmark** | measured |
| network | yes | SRA/S3, GIAB truth, and Genozip's activation |

External tools you must install: SPRING, Genozip **(interactive activation)**,
rtg-tools, DiscoSNP++, Kmer2SNP + conda, KMC, SPAdes, bwa, samtools, mosdepth,
seqtk, bcftools, tabix, SRA toolkit, AWS CLI. Exact versions and provenance:
`../../server/TOOLS.md`. Every install gotcha: `../../server/SETUP_FULL.md`.

## 1. Build the machine

    bash server/verify_environment.sh      # 27 checks, seconds, before you download 55 GB

Follow `server/SETUP_SHORT.md`. When something fails, the matching section of
`server/SETUP_FULL.md` has the exact error text and the fix — every failure
listed there has already happened here.

## 2. Get the data

    bash /root/arcs-clean/benchmark/download.sh /data/fastq 2>&1 | tee /data/download.log

19 datasets, 55.0 GB, plus `~/refs/chr20.{fa,sdf}` and `~/giab_truth/`.
Manifest with per-file byte counts: `../../server/DATASETS.md`. Pull the two
>10 GB SRA files from the S3 open-data mirror; NCBI prefetch is far slower.

## 3. Gate before you spend five hours

    bash benchmark/scripts/run_tests.sh              # 15/15
    bash benchmark/scripts/benchmark_0_preflight.sh  # must print VERDICT: GO

The preflight rebuilds both binaries from source, checks every tool, dataset,
reference and truth set, exercises the archive-calling path end to end on a
synthetic input, re-derives the validated k-mer count and F1, and prints the
projected runtime per phase. **A NO-GO names the missing item. Do not proceed
past it.**

## 4. One dataset end to end (~45 s) before the sweep

    CLAIMS=1 bash benchmark/scripts/sanity_archive_one.sh SRR2584863

**Checkpoint — the single most useful check in the project:**

    archive        68,429,027 B   exactly
    ratio          9.844%
    lossless       3/3 tools

If the archive size differs **on any machine of any size**, that is a real
regression, not an environment difference. It has been verified byte-identical
across 2/3/7/12 cores and every candidate-concurrency setting.

## 5. The full sweep (~5 h)

    bash benchmark/scripts/benchmark_1_run.sh

Runs one timed job at a time on an idle box. It prints per-dataset progress
with actual-vs-projected timing, writes durable checkpoints to `PROGRESS.txt`,
and **halts the entire run on any LOSSY archive** (a correctness gate, not a
tolerance). Produces all eight tables into
`results/benchmark_1_<timestamp>/`.

## 6. Compare against the recorded results

    diff <(sort benchmark/results/claim1/claim1_T1.1_T1.2.csv) \
         <(sort results/benchmark_1_*/claim1_T1.1_T1.2.csv)

**What must match exactly**

- every `archive_bytes` for CAPSULE
- every `lossless` verdict (57/57 LOSSLESS)
- every F1, TP, FP, FN in Claim 2
- T3.4's coordinate/content counts

**What will legitimately differ, and why**

| field | why |
|---|---|
| `compress_s`, `decompress_s`, `wall_s` | machine speed. Ratios vs SPRING/Genozip should hold, since all tools see the same machine |
| `peak_ram_kb` | machine and allocator |
| Genozip's `archive_bytes` | it embeds run metadata; drifts ~37 B (0.005%) run to run. **Never gate on it** |
| SPAdes / bwa baseline times | machine speed |

## 7. Verify individual claims without the full sweep

    # Claim 2, one individual, from the archive:
    bash benchmark/scripts/run_fullchr20_archive_capsule.sh \
         /tmp/capsule_bin/capsule_decode scripts ~/refs/chr20.fa <archive> HG002 /tmp/out

    # T3.4 (the Claim 3 mechanism), 100 sites.
    # THE SIDECAR MUST BE BUILT FIRST, AND WITH CAPS_PILEUP=1 -- without that
    # flag the sidecar carries pg + placements but NOT each read's own
    # deviations, `query` then emits the CONSENSUS rather than the reads, no
    # allele difference can appear, and the coordinate arm scores 0/100
    # instead of 18/100. Verified 2026-09-10: this is exactly what happens,
    # three times, if the flag is omitted. `query` picks the sidecar up
    # automatically from <archive>.qidx.
    CAPS_PILEUP=1 /tmp/capsule_bin/capsule_decode index <archive> <archive>.qidx
    bash benchmark/scripts/run_locus_fidelity.sh \
         /tmp/capsule_bin/capsule_decode <archive> ~/refs/chr20.fa HG002 /tmp/t34 100

    # T2.4 multi-allelic:
    PATH="$HOME/DiscoSnp:$PATH" bash benchmark/scripts/run_multiallelic_bench_capsule.sh \
         /tmp/capsule_bin/best106 scripts ~/refs/chr20.fa HG002 20:1000000-6000000 /tmp/t24

## 8. How configuration changes the answer — and where it does not

**Does not change the output at all** (verified byte-identical):

| knob | tested |
|---|---|
| core count | 2, 3, 7, 12 → identical archive; identical VCF |
| candidate concurrency | `CAPS_CAND_SEQ=1`, `CAPS_CAND_PAR=4` → identical archive |
| caller memory budget | default 5,000 MB vs `CAPS_MAXRAM_MB=2000` → identical VCF |

The caller's k-mer budget is an **absolute 5,000 MB default** (GATB's value)
with system RAM only as a *cap*, so any machine above ~7.5 GB free uses the
same budget the published runs used. Below that the cap binds and the spill
partitions differently — which changes *where* k-mers are counted, never
*which ones exist*.

**Does change the output — do not vary these when reproducing:**

| knob | effect |
|---|---|
| `CAPS_CALL=1` instead of `CAPS_SPANS=1` | runs the caller inline during compression, ~20× slower, and Claim 2 then calls again from the archive. It once put HG002 compression at 913 s against SPRING's 52 s |
| omitting `CAPS_SPANS=1` | archive lacks `contig_spans`; Claim 2 fails with `ARCHIVE LACKS contig_spans` |
| omitting `CAPS_NAMES=1` / `CAPS_QUAL=1` | a smaller archive that is **not** a lossless FASTQ archive — not comparable to SPRING or Genozip |
| substituting a dataset | invalidates comparison against the locked numbers; see the banned list in `../../server/DATASETS.md` |
| running two timed jobs at once | contaminates wall time and peak RAM. The harness enforces one at a time |

## 9. Full detail

`REPRODUCIBILITY.md` in this folder carries the raw evidence for §8 — the
actual measurements across core counts and memory budgets, and the three things
that legitimately vary.
