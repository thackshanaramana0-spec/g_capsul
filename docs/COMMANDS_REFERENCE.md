# CAPSULE — commands reference, by claim

Written 2026-09-03. Pure commands: what to run and what it does, nothing
else. For *why* any of this works, see `docs/TECHNICAL_ARCHITECTURE.md`
(Claim 1), `docs/CLAIM2_TABLES_AND_INDEL_SCAN.md` (Claim 2), and
`docs/CLAIM3_LOCKED.md` (Claim 3). Every command below is a real,
already-used invocation, not a hypothetical interface.

## The 3-line version — start here

```bash
bash scripts/run_capsule.sh 1     # Claim 1 (COMPACT):     verify losslessness on the locked E. coli dataset
bash scripts/run_capsule.sh 2     # Claim 2 (FAITHFUL):    fast synthetic regression test, no downloads needed
bash scripts/run_capsule.sh 3     # Claim 3 (ADDRESSABLE): fast synthetic regression test, no downloads needed
```

Each claim is **independent** — there's no requirement to run them in
order, and no reason a reviewer checking Claim 3 should have to run 1 or 2
first. Each command builds whatever binaries it needs on first use and
prints which dataset it used.

**No dataset path is hardcoded anywhere in this chain.** Every script reads
from `scripts/capsule_config.sh` — edit that one file if your datasets move,
or override for a single run without touching any file:
```bash
CAPSULE_DATA_DIR=/mnt/other/fastq bash scripts/run_capsule.sh 1
```

**Sub-phases**, once the default has run cleanly, per claim:
```bash
bash scripts/run_capsule.sh 1 compress   # keep the archive, not just verify
bash scripts/run_capsule.sh 1 sweep      # the real 14-dataset comparison (not yet a single script -- prints where to find the manual steps)
bash scripts/run_capsule.sh 2 giab       # real GIAB het-SNV+indel benchmark (needs chr20.fa + truth VCFs, see SERVER_SETUP_AND_DOWNLOADS.md)
bash scripts/run_capsule.sh 2 window HG002 r2   # single-individual, single-window benchmark
bash scripts/run_capsule.sh 3 full       # the real export/coverage/query benchmark vs SPAdes/bwa+mosdepth
```

Everything below this point is the direct, lower-level command for each
operation, for anyone who wants to skip the dispatcher and call the
underlying binaries/scripts themselves.

---

Build the two binaries once, shared by all three claims:

```bash
scripts/build106.sh /tmp/best106        # encoder
scripts/build_decode.sh /tmp/capsule_decode   # decoder + Claim 3 operations
```

---

## Claim 1 — COMPACT (compression)

**Compress** (the default sweep — tries multiple MAXMAP/MINOV candidates,
keeps the smallest):
```bash
INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 bash scripts/encode_adaptive.sh
```

**Compress with names/quality/N-handling turned on** (full FASTQ, not
sequence-only):
```bash
CAPS_NAMES=1 CAPS_QUAL=1 INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 \
    bash scripts/encode_adaptive.sh
```

**Decompress** (full 4-line FASTQ reconstruction):
```bash
/tmp/capsule_decode out.capsule <outdir> <outdir>/reads.fq
```

**Verify losslessness** (encode, decode, diff against the original file's
own sequence column — one command, no manual steps):
```bash
scripts/verify_lossless.sh reads.fq
```

**Full compress/decompress round trip via one file** (what a real user
would actually run):
```bash
scripts/capsule_roundtrip.sh reads.fq
```

---

## Claim 2 — FAITHFUL (variant calling)

**Call variants directly from an encode pass** (no separate decode step —
calling happens in-process during compression):
```bash
CAPS_CALL=1 CALL_VCF=calls.vcf CAPS_DUMP_CONTIGS=contigs.tsv \
    /tmp/best106 reads.fq 3 16 16 22 16 16 1 24 64 1
```

**Multi-allelic calling** (admit up to k co-occurring alleles per site,
emits `ALT=C,G`-style records):
```bash
CAPS_CALL=1 CAPS_PLOIDY=3 CALL_VCF=calls.vcf /tmp/best106 reads.fq 3 16 16 22 16 16 1 24 64 1
```

**Full real-GIAB het-SNV + het-indel benchmark** (caller → bwa-lift contigs
to genome coordinates → rtg vcfeval against GIAB truth, one command):
```bash
bash scripts/run_giab_indel_capsule.sh /tmp/best106 scripts reads.fq ~/refs/chr20.fa
```

**Window/individual benchmark** (T3-style: downloads/streams that
individual's chr20 window from GIAB S3 at the standardized 30× depth, calls,
lifts, scores against truth — `window` is `r2`/`r3`/`r4`/`r5`/`na` or an
explicit `LO:HI` range; `r2` is the ONLY window ever used for parameter
tuning, per this project's own held-out discipline):
```bash
bash scripts/run_window_bench_capsule.sh /tmp/best106 scripts ~/refs/chr20.fa HG002 r2
```

**Real multi-allelic validation** (synthetic-triploid or real-pooled
workdir must already contain `reads.fq`; `K` is the ploidy passed to
`CAPS_PLOIDY`):
```bash
bash scripts/run_polyploid_bench_capsule.sh /path/to/workdir /tmp/best106 scripts 3
```

**Synthetic regression test** (deterministic, no GIAB download needed —
run this first on any caller change, before spending time on real data):
```bash
bash scripts/test_claim2.sh
```

---

## Claim 3 — ADDRESSABLE (export / coverage / query)

**Export the assembled pseudogenome as FASTA** (no assembler run — the
archive already has it):
```bash
/tmp/capsule_decode export out.capsule contigs.fa
```

**Per-base coverage depth** (no alignment, no reference — reads its own
placement streams):
```bash
/tmp/capsule_decode coverage out.capsule coverage.tsv
```

**Reads overlapping a coordinate range** (reference-free coordinate query;
note this decodes the whole archive first — see the O(archive) limitation
in `CLAIM3_LOCKED.md` §5):
```bash
/tmp/capsule_decode query out.capsule region.fa 0-100000
```

**One-command reproducible benchmark** (builds binaries, installs SPAdes if
absent, compresses, runs all three operations plus their conventional
baselines — the industrial-grade runner):
```bash
bash scripts/run_claim3.sh /data/fastq ./results/claim3
```

**Synthetic regression test** (catches the exact bug class this claim's
real audit found — a coverage-total invariant, verified against the actual
pre-fix binary):
```bash
bash scripts/test_claim3.sh
```

---

## Cross-claim: the outer project's own benchmark harness

Run from `/root/arcs-clean` (the outer repo), against the **outer** `arcs`
binary, not this sandbox's `best106`/`capsule_decode` — see
`docs/CLAIM3_LOCKED.md` §7 for why these are not yet the same binary:

```bash
bash benchmark/benchmark.sh claim1 /data/fastq /root/arcs-clean/build/arcs ./results
bash benchmark/benchmark.sh claim2 /data/fastq /root/arcs-clean/build/arcs ./results
bash benchmark/benchmark.sh claim3 /data/fastq /root/arcs-clean/build/arcs ./results
bash benchmark/benchmark.sh all    /data/fastq /root/arcs-clean/build/arcs ./results
```
