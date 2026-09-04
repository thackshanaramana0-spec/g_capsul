# Invocation errors — bugs that were MINE, not the code's

Every entry here is a case where a tool or script was declared broken, or a
result was declared wrong, when the real fault was **how it was invoked**.
They are recorded because each one cost real time, several nearly reached a
paper, and the same class of mistake recurred four times in one session.

**The rule they all violate:** *read the tool's own help/usage before
concluding the tool is broken.* Every single one of these was answerable in
one `--help` or one `grep` of a script header.

---

## 1. SPRING `-g` on compress — would have voided every SPRING number

**Symptom:** `spring -c` failed on every input with `Program terminated
unexpectedly with error: No reads found.` Reproduced on a 2-read hand-written
FASTQ, on 1,000 reads and on 10,000 reads. `strace` showed SPRING opening the
input and creating its temp files, then reporting no reads. Two independent
builds (`/usr/local/bin/spring` and `/root/SPRING/build/spring`) failed
identically. **I concluded SPRING was broken on the machine.**

**Actual cause:** `-g` means *"the input is gzipped"*:

    -g [ --gzipped-fastq ]   enable if compression input is gzipped fastq
                             or to output gzipped fastq during decompression

Our inputs are plain `.fq`. SPRING looked for gzip data, found none, and
reported exactly what was true.

**Where it came from:** CLAUDE.md says *"SPRING — always needs `-g` for FASTQ
output on decompress"*. That is correct and applies to DECOMPRESS. It was
copied onto the compress call.

**Blast radius:** the same wrong flag was in `scripts/run_claim1_bench.sh`,
the script that produces the published T1/T2 tables. Every SPRING row from
that script was a FAILURE being recorded, not a measurement.

**Correct invocation:**

    spring -c -i reads.fq -o out.spring -t $(nproc)          # NO -g: plain input
    spring -d -i out.spring -o dec.fq  -t $(nproc)           # NO -g: plain output

Use `-g` only when the input really is `.fq.gz`, or when a gzipped output is
wanted. `-g` on decompress writes gzip bytes into a file named `.fq`, which a
plain-text comparison then reads as garbage and reports LOSSY — a second,
independent false failure from the same flag.

**Verified after the fix** (SRR29296997, 460,501 reads): compress 4.41 s,
17.5 MB, 868 MB RAM; decompress 4 s; round trip byte-identical once the `+`
line is normalised on both sides.

---

## 2. SPRING `+` line — a lossless tool reported LOSSY

**Symptom:** after fixing #1, the SPRING round trip still differed from the
original.

**Actual cause:** SPRING drops the optional ID after the `+` line. The FASTQ
spec makes that field redundant (it must repeat line 1 if present), so no data
is lost. CLAUDE.md already specifies normalising it on both sides so the
comparison is fair on *data recovery*.

**Correct comparison:**

    norm(){ paste - - - - < "$1" | awk -F'\t' '{print $1"\t"$2"\t+\t"$4}' | tr -d '\r' | sort; }
    cmp -s <(norm orig.fq) <(norm dec.fq) && echo LOSSLESS || echo LOSSY

---

## 3. Running the encoder bare — "Claim 1 and Claim 3 are broken"

**Symptom:** an archive built by calling the binary directly decoded to
**zero reads**; `coverage` reported 0 placements and `query` returned nothing.
**I diagnosed a data-loss bug in Claim 1, traced it to `read_lengths` being
empty, and edited `stages/106_inprocess.cpp` on that premise.**

**Actual cause:** the encode was invoked without the required environment:

    DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1    # required, NOT debug flags
    CANDIDATES=...                       # the adaptive candidate sweep

`docs/TECHNICAL_ARCHITECTURE.md:456` states plainly: *"`DUMP_PERM=1 DUMP_MM=1`
are REQUIRED, not optional debug flags."* Every real encode path
(`encode_adaptive.sh`, `verify_lossless.sh`, the pipeline scripts) sets them.
The `read_lengths` write sits inside `if(getenv("DUMP_PERM"))`, so a bare
invocation produces a deliberately incomplete archive.

**Blast radius:** none — the false-premise edit was reverted before it was
committed. But it consumed a long debugging detour and produced a confident,
wrong claim that two of the three paper claims were broken.

**Correct invocation:** never call the binary directly for anything that will
be measured. Use:

    INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/capsule_bin/best106 \
        bash scripts/encode_adaptive.sh

---

## 4. Quoting a number measured in a different configuration

**Symptom:** `PIPELINE.md` says the caller takes 98.5 s on full chr20; every
log I looked at said 1131 s. **I announced the doc was wrong.**

**Actual cause:** two different configurations.

| configuration | traversal | produces | time |
|---|---|---|---|
| default (no `CAPS_DBG_SB`) | parallel, 12 cores | SNVs | ~98.5 s |
| `CAPS_DBG_SB=1` | **serial** (`#pragma omp parallel if(!WANT_SB)`) | + indels | ~1131 s |

1131 / 12 ≈ 94 s. Same work, same machine, one core versus twelve. I had been
running everything with `CAPS_DBG_SB=1` while chasing indels, then compared it
against a figure measured in the SNV configuration. The doc was right.

(That serial gate is real and was worth fixing — it is what layer 15 removed —
but it was not evidence of a wrong number.)

---

## 5. Argument order taken from memory, not from the script

**Symptom:** caught before running. The DiscoSNP++ arm in `benchmark_1_run.sh`
was called as `<scripts_dir> <ref> <reads> <individual> <outdir>`.

**Actual signature**, from the script's own header:

    run_fullchr20_bench_disco.sh <ref.fa> <reads.fq> [individual] [outdir]

A wrong argument order here does not crash — it benchmarks the competitor on
the wrong input and silently produces a number that looks real.

---

## What to do differently

1. **Read the tool's `--help` before calling it broken.** Entries 1 and 5 were
   each one command away from being avoided.
2. **`grep` the script header for the invocation before writing a caller.**
   Entry 3's answer was a documented line in `docs/`.
3. **When a measurement disagrees with a document, suspect the configuration
   first.** Entry 4 was two flags apart, not a wrong number.
4. **Never invoke a measured binary directly.** Go through the wrapper that
   sets the required environment.
5. **A tool failing on a 2-read hand-written file is a strong signal — but of
   the INVOCATION, not the tool.** A genuinely broken binary usually fails
   differently (crash, link error), not with a clean, accurate diagnostic.

## How they were caught

All five were caught by a **smoke test on one small dataset** before the full
run. That is the argument for always running `SANITY_ONLY=1` first: it cost
one minute and caught a bug that would otherwise have silently voided every
SPRING row in the paper's headline table.
