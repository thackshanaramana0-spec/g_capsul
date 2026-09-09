# Reproducibility: does the same input give the same result on a different box?

**Tested 2026-09-10, not assumed.** The question that motivates this: every
number in the paper was produced on one 12-core, 82 GB machine. If a reviewer
re-runs on 4 cores or 16 GB, do the archives, the variant calls and the tables
come back the same?

**Answer: yes for output, no for wall time — which is the correct behaviour.**

## What was tested and what it showed

### 1. Archive vs candidate concurrency

The encoder derives how many adaptive candidates run at once from measured free
RAM (`/proc/meminfo` MemAvailable), so machine state does enter the code. It
enters the SCHEDULE only, never the SELECTION: `cands` is fixed and `K` is the
wave size, so every candidate is always evaluated.

    E. coli SRR2584863            archive           K
      baseline                    68,429,027 B      2
      CAPS_CAND_SEQ=1 (K=1)       68,429,027 B      1     IDENTICAL
      CAPS_CAND_PAR=4 (K=4)       68,429,027 B      4     IDENTICAL

### 2. Archive vs core count

    2 cores    68,429,027 B   36.17 s   IDENTICAL
    3 cores    68,429,027 B   29.39 s   IDENTICAL
    7 cores    68,429,027 B   16.71 s   IDENTICAL
    12 cores   68,429,027 B   12.38 s   (baseline)

Time changes 2.9x across the range; **not one byte of the archive changes**.

### 3. Decoder determinism and losslessness

Decoded reads md5 identical at 2, 7 and 12 cores, and the result is LOSSLESS
against the original FASTQ (order-free comparison of the sequence column).

### 4. The caller, from the archive

    12 cores, default RAM     93.90 s   11.0 GB   75,250 variants
     4 cores, default RAM    200.68 s   10.4 GB   75,250 variants   IDENTICAL VCF
    12 cores, CAPS_MAXRAM_MB=2000
                              92.72 s   11.1 GB   75,250 variants   IDENTICAL VCF

Byte-identical VCFs across a 2.1x change in wall time and a memory ceiling
restricted from 5,000 MB to 2,000 MB.

## Why memory does not change the answer, by design

The caller's k-mer budget is an **absolute default of 5,000 MB** (GATB's value,
adopted deliberately), with system memory used only as a CAP:

    mem_ceil_default = 5000;
    if (MEM_AVAIL_MB && mem_ceil_default > (MEM_AVAIL_MB*2)/3)
        mem_ceil_default = (MEM_AVAIL_MB*2)/3;

So on any machine with roughly 7.5 GB or more available, the budget is the same
number the paper's runs used, and behaviour is identical. Below that the cap
binds and the spill partitions differently — which changes WHERE k-mers are
counted, never which ones exist, because the merge sums counts per key either
way. The k-mer-count identity assertion in the harness checks exactly that.

This was a deliberate correction: the budget used to be 60% of MemAvailable,
which silently DISABLED the spill on a large box (28.3 GB and 191.7 s instead
of ~7 GB and 25.4 s). A budget discovered from the machine makes a tool behave
worse the bigger the machine is; a declared budget does not.

## What legitimately varies

- **Wall time and peak RSS.** They are properties of the machine, and T1.2/T2/T6
  report them as measurements, not as constants.
- **Genozip archive bytes**, by about 37 B run to run (0.005%) — it embeds run
  metadata. Never gate on byte-identity of a Genozip archive.
- **The graph-only caller path** (`CAPS_DBG_ONLY`) renumbers `dcontig_N`
  non-deterministically, so gate that path on F1 after lifting, not on VCF
  bytes. The SNV+indel path tested above IS byte-deterministic.

## How to re-verify on a new machine

    bash scripts/benchmark_0_preflight.sh          # must print GO
    bash scripts/run_tests.sh                      # 15/15
    CLAIMS=1 bash scripts/sanity_archive_one.sh SRR2584863

The third reproduces T1.1/T1.2 for one dataset in about a minute; the archive
should be 68,429,027 B and every tool should report LOSSLESS. If the archive
size differs on a machine of any size, that is a real regression, not an
environmental difference — that is the point of this document.
