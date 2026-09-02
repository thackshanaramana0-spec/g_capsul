# Final algorithmic scan — CAPSULE as one product, 2026-09-03

Every prior audit this project ran (Claim 1's own bug-finding history, and
this session's Claim 2/3 checklists) tested **one claim's features at a
time**. This scan's specific purpose was different: exercise all three
claims' features **together, on the same archive**, since that combination
had never actually been run before — and a whole-product claim needs that
combination to work, not just each third of it in isolation.

## What was tested

```bash
export CAPS_CALL=1 CAPS_NAMES=1 CAPS_QUAL=1 CALL_VCF=calls.vcf
INPUT=reads.fq ARCHIVE=combined.capsule BEST=/tmp/best106 \
    bash scripts/encode_adaptive.sh
```

then running Claim 3's `export`/`coverage`/`query` and a full decode against
the resulting `combined.capsule` — the exact scenario a real user
compressing one file with every feature turned on would hit, and the exact
scenario none of the per-claim checklists exercised, because each one
built its own test archive with only its own claim's flags set.

## Real bug found: names/quality silently dropped through the sweep

**Symptom:** with `CAPS_NAMES=1 CAPS_QUAL=1` set, compressing via
`encode_adaptive.sh` (the documented, recommended compression path —
`README.md`, `TECHNICAL_ARCHITECTURE.md` §8) produced an archive that
decoded to **zero names, zero quality bytes**, silently — no error, no
warning, `ARCHIVE_TOTAL` printed normally, exit code 0. Direct
single-candidate invocation of the same binary with the same flags worked
perfectly (verified: 3176 B names, 43480 B quality, both round-tripped
correctly).

**Root cause, isolated precisely:** `encode_adaptive.sh` triggers
`stages/106_inprocess.cpp`'s `CANDIDATES` fork mechanism (§9 of
`TECHNICAL_ARCHITECTURE.md`'s companion, or read the code at
`106_inprocess.cpp:996-1200`). Each forked child `chdir()`s into a
per-candidate scratch directory (line 1197) before continuing the encode.
`CALL_VCF` and `CAPS_DUMP_CONTIGS` were already resolved to absolute paths
**before** this chdir, with an explicit comment explaining exactly why
("so a relative path the caller passed still refers to the directory they
meant") — this exact bug class had already been found and fixed **once**,
for those two variables. `g_input_path` — which `nmc::encode_from_fastq`
and `qlc::encode_from_fastq` reopen to extract names and quality — never
got the same treatment. After the chdir, a relative input path (the normal
case: `INPUT=reads.fq`, not `INPUT=/full/path/reads.fq`) pointed at a file
that doesn't exist in the scratch directory, so both encoders silently
processed zero records.

**This also affected `GSEARCH=1`**, which has its own, separate chdir
(`106_inprocess.cpp:1069`) with the identical hazard — not just the plain
`CANDIDATES` sweep.

**Fix** (`stages/106_inprocess.cpp`, `main()`): resolve `g_input_path` to
an absolute path immediately after reading it from `argv[1]`, before any
fork can happen — the single canonical point that covers both the
`CANDIDATES` and `GSEARCH` code paths at once, rather than patching each
fork branch separately.

**Verification, not assumption:**
1. The exact original repro (relative path, through `encode_adaptive.sh`,
   `CAPS_NAMES=1 CAPS_QUAL=1`) now produces 1000/1000 names and quality
   records, sequence still lossless.
2. Byte-identical output, pre-fix vs post-fix binary, for (a) plain
   single-candidate encode and (b) the sweep with names/quality OFF —
   confirmed via `cmp`, meaning **every existing locked result in this
   project is untouched** (this only ever engages when names/quality are
   requested through a candidate sweep, which no currently-published number
   in this project uses).
3. `scripts/verify_lossless.sh` on a real locked dataset (`ERR5181310`)
   with the fixed binary reproduces the exact previously-measured archive
   size, 836,191 bytes — matching `results/phase_a/allphases_14dataset.csv`
   byte for byte.
4. The full three-claims-combined scenario (calling + names + quality, all
   at once) now works end to end: VCF produced, names and quality
   round-trip, sequence lossless, `export`/`coverage` both run cleanly on
   the same archive.

**Why this matters beyond the one bug:** this is the second time this
exact bug pattern (a relative path silently broken by a post-chdir sweep)
has been found in this codebase — the first time produced the fix already
applied to `CALL_VCF`/`CAPS_DUMP_CONTIGS`. That fix was correct but local;
it didn't generalize to the next variable with the same hazard. The lesson
this scan is recording, not just the bug: **when a fix addresses "a
relative path breaks after this chdir," check every other variable that
gets reopened after the same chdir, not only the one that happened to be
reported.** No other such variable was found in this pass (grep for
`getenv` calls whose value is later used to `fopen`/reopen a file,
cross-checked against every chdir site) — `CALL_VCF`, `CAPS_DUMP_CONTIGS`,
and `g_input_path` are now the complete set that needed this treatment, and
all three now have it.

## Other cross-cutting checks run in this pass

- **Does `capsule_decode` (Claim 3's operations) handle an archive built
  with Claim 2's `CAPS_CALL=1` correctly?** Yes — `CAPS_CALL` only affects
  what the *encoder* does as an in-process side effect (writing a VCF to a
  separate path); it adds no new archive streams and changes no existing
  one, so the decoder's behavior is identical with or without it. Confirmed
  directly in the combined test above (export/coverage ran cleanly on an
  archive built with `CAPS_CALL=1`).
- **Does Claim 1's `verify_lossless.sh` remain valid as the trust anchor
  after this fix?** Yes, and it was used as exactly that in this scan —
  the fix was validated against it, not the other way around.
- **Is there a second candidate-sweep-vs-relative-path hazard anywhere in
  the Claim 2 or Claim 3 scripts?** Checked: `run_giab_indel_capsule.sh`,
  `run_window_bench_capsule.sh`, and `scripts/run_claim3.sh` all invoke the
  encoder directly (no `CANDIDATES`/`GSEARCH`), so none of them can trigger
  the fork/chdir path this bug lived in. `test_claim2.sh`/`test_claim3.sh`
  likewise call `build106.sh`/`encode_adaptive.sh` from within their own
  workdir with absolute-resolved paths already (`$W/synth.fq`), so neither
  test could have caught or been affected by this bug — which is exactly
  why running all three claims together, deliberately, was necessary to
  find it at all.

## Verdict

One real, previously-undiscovered, now-fixed cross-cutting bug — found
specifically because this scan tested claims together rather than in
isolation, which is the whole point of treating CAPSULE as one product
rather than three separate deliverables. No other cross-cutting defect was
found in this pass. This does not replace the per-claim algorithmic audits
already done (Claim 3's coverage-undercount and region-boundary bugs,
Claim 1's four silent data-loss bugs) — it is the complement to them: those
found bugs *within* a claim's own code path, this one found a bug in the
*interaction* between claims.
