# T2.4 correction, 2026-09-10 — "21/26" was never reproducible; the verified figure is 17/26

Found while performing a technical-freeze reproducibility audit: walking every
published table backward to a raw log that actually supports it. T2.4's
"21/26 (80.8%), complete chr20 census" had no raw log anywhere on this server.
This document is that raw log, for the number that replaces it, plus the full
investigation of why the old number existed.

**Corrected figure: CAPSULE both-allele sites: 17 / 26 (65.4%). DiscoSNP++: 0 / 26 (unchanged — still structural).**

---

## What was checked, in order

### 1. Independent re-derivation of the denominator (26)

Streamed the full chr20 region (`20:1-63025520`) from the real GIAB HG002
300x BAM, extracted real truth restricted to GIAB confident regions, and
counted multi-allelic (`GT=1/2`) sites whose two ALT alleles are both
single-base. Result: **952 truth multi-allelic sites total, 26 of them
SNV-only** — matching the published denominator exactly. **This part of the
original claim is correct and independently confirmed.**

### 2. First full run: `run_A_full_chr20.log`

Ran `scripts/run_multiallelic_bench_capsule.sh` end to end on that same
full-chromosome stream, current encoder (post the three defect fixes in
commit `b28dd40`). Result: **17 / 26.**

### 3. Determinism check: `run_B_same_reads_repeat.log`

This project has documented caller nondeterminism before (bubble-traversal
order affecting `dcontig` numbering, which can flip which duplicate record
`lift_vcf.py`'s dedup keeps at a given `(gpos, len(ref), len(alt))` — see
`scripts/lift_vcf.py` line 183). To test whether run A's 17/26 was itself an
artifact of that, re-ran the identical pipeline on the **identical, already
-streamed** `reads.fq` (not re-sampled from the BAM). Every intermediate
count matched run A exactly: 38,887 SNVs, 7,022 indels, 151,286 contigs,
48,635 lifted calls. Final score: **17 / 26, again.** The caller is stable on
this input; nondeterminism is not the explanation.

### 4. Pre-fix binary check: `run_C_prefix_binary_same_reads.log`

Commit `b28dd40` (three defect fixes for out-of-scope input) is the only
commit touching the encoder between the original T2.4 measurement (commit
`124dc74`) and now. `include/caps_caller.h` — the actual variant-calling
logic — is byte-identical across that range; only `stages/106_inprocess.cpp`
changed, and only in ways gated on conditions no real short-read data
triggers. To confirm rather than assume this, built the encoder at
`124dc74^` (`bc722c7`, via `git worktree`) and re-ran the identical pipeline
on the identical `reads.fq`. Every intermediate count matched runs A and B
exactly, including the pre-fix binary. Final score: **17 / 26, a third time.**
**The encoder fix is confirmed inert here — it is not the cause.**

### 5. Searching for the original measurement's raw evidence

Found `/root/.claude/jobs/66d96ba8/tmp/t24fix/` — a working directory from
earlier the same day. Its `regions.bed` spans **chr20:1,000,000–6,000,000**
(a 5 Mb window), its `reads.fq` holds 1,284,795 reads (not the ~11.6M a full
chromosome stream gives), and its own `truth_multiallelic.vcf` has exactly
**111** records of which **7** are SNV-only — matching this project's
**earlier, already-superseded** window-based figures (`5/111`, `7 SNV-only`),
not 26.

**No raw log or intermediate file anywhere on this server supports a caller
run against the full 63 Mb chromosome that recovered 21 of 26 sites.** The
"21/26, complete chr20 census" commit message (`124dc74`) states real,
independently-reproducible truth-side numbers (952 sites, 26 SNV-only) —
which do require scanning the whole chromosome and were almost certainly
computed correctly, since a truth-only scan needs no read streaming and is
cheap — but the **caller-recovery figure paired with them does not have a
corresponding execution**, and every attempt to reproduce it, three times,
independently, gives 17/26.

## Conclusion

**21/26 is withdrawn as unreproduced.** The verified, reproducible figure —
three independent runs, two encoder versions, identical intermediate counts
each time — is **17/26 (65.4%)**. Updated everywhere the old figure appeared;
see each document's own retraction note for exactly what changed.

**What is unaffected:** DiscoSNP++'s 0/26 (a structural property — it emits
zero multi-allelic records in 3,989 total, independent of any recovery
metric), the 952/26 truth-side census, and every OTHER table in this project
— T1.1 through T2.3, and T3.1 through T3.4, were independently re-derived
from the raw sweep log during this same audit and matched their published
CSVs exactly (see `benchmark/documentation/RESULT_CODE.md`'s cross-check
note). This is a defect in one number in one table, found and fixed, not a
reason to distrust the others.

## Reproduce this correction

```bash
PATH="$HOME/DiscoSnp:$PATH" bash scripts/run_multiallelic_bench_capsule.sh \
    <capsule_encoder> "$(pwd)/scripts" ~/refs/chr20.fa \
    HG002 20:1-63025520 <outdir>
```

Needs `~/refs/chr20.fa`, GIAB truth for HG002, and network access to the GIAB
S3 BAM. Takes roughly 17 minutes on a 12-core box (12 minutes streaming the
BAM, 3 minutes compressing, under a minute calling, under a minute for
DiscoSNP++ and scoring).
