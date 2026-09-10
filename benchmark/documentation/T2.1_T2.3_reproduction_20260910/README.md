# T2.1 + T2.3 HG002 re-executed, 2026-09-10 — 12/12 fields exact

Done after a fair challenge: the claim that closing the tier-2 verification gap
needed "~19 h of compute" was lazy. It needed **90 seconds**, because an HG002
archive from an earlier run had survived on disk and the caller can be pointed
straight at it.

**One run verifies two tables**, because the caller emits both classes in the
same pass.

    published (claim2_T2.1_snv.csv, HG002/CAPSULE)
      TP=37011  FP=1721  FN=7564   P=0.956  R=0.830  F1=0.888
    re-run 2026-09-10
      SNV TP=37011 FP=1721 FN=7564 P=0.956 R=0.830 F1=0.888   EXACT

    published (claim2_T2.3_indel.csv, HG002/CAPSULE)
      TP=3871   FP=824   FN=3913   P=0.825  R=0.497  F1=0.620
    re-run 2026-09-10
      INDEL TP=3871 FP=824 FN=3913 P=0.825 R=0.497 F1=0.620    EXACT

All 12 fields across both tables reproduce exactly. Raw log:
`rerun_HG002_from_archive.log`.

## Why this counts as independent

The published numbers came from `benchmark_1_run.sh` on 2026-09-09. This run
used a **different archive file** (built 2026-09-10 00:45 by a separate job),
a **freshly built decoder**, and a direct invocation of the runner rather than
the sweep harness. It shares no shell variable, no process and no intermediate
file with the original — which is exactly the property the log-vs-CSV
cross-check lacked.

It also independently confirms HG002's T1.1 `archive_bytes`: the surviving
archive is **573,767,964 B**, matching the published value exactly.

## Reproduce

```bash
bash scripts/run_fullchr20_archive_capsule.sh \
     <capsule_decode> "$(pwd)/scripts" ~/refs/chr20.fa \
     <HG002.capsule> HG002 <outdir>
```

Needs `~/refs/chr20.fa`, `~/refs/chr20.sdf`, GIAB truth for HG002, and `rtg`.
About 90 s against an existing archive.

## Scope

HG002 only — 1 of the 4 individuals in each table. HG003/HG004/HG005 remain
transcription-checked, and their archives did not survive. Re-running them
needs a re-compress (~4 min each) before the ~90 s call.
