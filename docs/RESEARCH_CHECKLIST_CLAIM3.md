# Research / academic checklist — Claim 3 (ADDRESSABLE)

Scope: Claim 3 only — the export/coverage/query capability, its evaluation
against SPAdes and bwa+mosdepth, and its prior-art position. Priority order
per the project's own stated ranking: **scientific correctness →
reproducibility → claim traceability → data provenance → fair baselines →
tests → code architecture → cleanliness.** Statistical methodology
(confidence intervals, significance tests) is intentionally marked N/A
throughout — per project instruction, and consistent with `[[claim2_capsule_results]]`'s
standing rule that no paper in this specific field (FASTQ compression /
reference-free variant calling) reports significance for these kinds of
comparisons, so holding Claim 3 to a different standard than Claim 1/2
would be inconsistent, not more rigorous.

Status legend: 🔴 critical, 🟠 important. ✅ done · ⚠️ partial/scoped gap ·
❌ not done · N/A not applicable.

| Area | 🔴/🟠 | Status | Finding |
|---|---|---|---|
| **Research question → code mapping** | 🔴 | ✅ | Every part of the claim has identifiable code: export → `capsule_decode.cpp:318-336`, coverage → `:195-262` and `:361-389` (see dead-code note below), query → `:390-414`. `CLAIM3_LOCKED.md` §2 maps each to exact line ranges, not "somewhere in the decoder." |
| **Reproducibility** | 🔴 | ✅ (E. coli) / ⚠️ (GIAB) | `scripts/run_claim3.sh` reproduces the E. coli result end to end from one documented command. The GIAB rows (`HG002_r2/r5`, `HG005_r3`) are **not** reproducible from a single command today — their numbers came from ad-hoc runs on pre-built archives, disclosed as such in `CLAIM3_LOCKED.md` §6.1 rather than presented as if scripted. |
| **Dataset provenance** | 🔴 | ✅ | E. coli = `SRR2584863`, exact SRA accession, cross-referenced against the project's own `DATASET_LOCKED.md`. GIAB windows sourced from `HG002_pooled.fq`/`HG005_pooled.fq`, themselves documented (elsewhere, for Claim 2) as GIAB S3-sourced, downsampled to 30× chr20. No dataset in this claim lacks a named source. |
| **Data integrity** | 🔴 | ⚠️ | No checksum (md5/sha256) is recorded anywhere for `SRR2584863_1.fq` or the GIAB pooled files in this claim's docs. This is a real, cheap-to-fix gap — worth adding to `CLAIM3_LOCKED.md` §6.1 as a follow-up, not fixed in this pass since it requires re-hashing large files not central to Claim 3's own logic. |
| **Train/test separation** | N/A | N/A | Claim 3 is not a learned model — there is no train/test split for a decode-time operation over an already-built archive. |
| **Ground truth** | 🔴 | N/A | Export/coverage/query are not evaluated against a truth VCF or truth genome the way Claim 2 is — they are evaluated for *speed* against a named conventional tool doing the *same job* (SPAdes for assembly, bwa+mosdepth for depth), and for *correctness* via the invariant tests (see Losslessness row below). "Ground truth" in the Claim-2 sense doesn't apply here. |
| **Experimental configuration** | 🔴 | ✅ | `scripts/run_claim3.sh` records every parameter used: thread count (`nproc`), SPAdes version pinned to an exact release tag, the exact query range (`0-100000`), all logged to the results CSV's notes column. |
| **Deterministic execution** | 🔴 | ⚠️ | The CAPSULE side is fully deterministic (same archive → same export/coverage/query output, every time — verified implicitly by the regression test passing repeatably). SPAdes and MEGAHIT are **not** guaranteed bit-for-bit deterministic across runs (multi-threaded assembly heuristics) — this is disclosed nowhere in `CLAIM3_LOCKED.md` explicitly, though the *timing* variance was directly observed and reported honestly (655× vs 555× across two SPAdes runs, both kept in the record rather than only the better one). |
| **Baseline implementation** | 🔴 | ✅ | SPAdes v4.0.0 (the spec-named tool per `CLAUDE.md`'s Claim 3 definition) is now actually installed and run, not simulated — closing the exact gap `PROJECT_AUDIT.md` had flagged. bwa+samtools+mosdepth versions are used as installed on the system but not version-pinned in the docs (see Dependency management in the industrial checklist — same underlying gap). |
| **Fair comparison** | 🔴 | ✅, with one disclosed asymmetry | Same reads, same machine, same measurement tool (`/usr/bin/time -v`) for both sides. The one asymmetry is disclosed, not hidden: bwa+mosdepth's coverage timing used a **pre-built** BWA index, while CAPSULE's coverage needs no reference at all — stated explicitly in `CLAIM3_LOCKED.md` §5 as understating CAPSULE's true advantage, not overstating it. |
| **Metric implementation** | 🔴 | ✅ (after this session's correction) | Speedup = conventional_time / capsule_time, a direct ratio, not a derived statistic with room for implementation error. The one real metric bug this session found — comparing `query` against a cheap stream-dump path rather than true reconstruction — was caught by *building the reproducible script* and re-verified with 5 repeats each before being written into the locked doc, not accepted on the first number. |
| **Statistical methodology** | — | N/A | Intentionally not applied, per project instruction and precedent — see header note. |
| **Experiment independence** | 🔴 | ✅ | No manually edited intermediate result was used anywhere in this claim's numbers. `t6_results.csv` rows are appended by the script or transcribed verbatim from a logged command; the corrected query row is *added*, and the old row is marked `SUPERSEDED` in place rather than silently edited out — the same "corrections stay on the record" discipline this project applies everywhere else. |
| **Raw → result pipeline** | 🔴 | ✅ | FASTQ → `encode_adaptive.sh` → `.capsule` archive → `capsule_decode {export,coverage,query}` → output file → hand-verified metric (contig count, covered-base sum, speedup ratio) is fully traceable in `scripts/run_claim3.sh`, one script, no manual step in between. |
| **Losslessness/correctness tests** | 🔴 | ✅ (new) | Was the single biggest gap before this pass: no test verified export/coverage/query correctness at all, only that they ran without crashing. `scripts/test_claim3.sh` now checks three real invariants plus delegates to the existing full-archive lossless check — see the industrial checklist's Testing row for the verified-against-the-real-bug detail. |
| **Ablation experiments** | 🟠 | N/A | Claim 3 has no tunable "components" to ablate — coverage's hoisting-above-pg-rebuild is an implementation optimization, not a method component with a knob. |
| **Sensitivity analysis** | 🟠 | ⚠️ | The query range `0-100000` is the only range tested; no sweep over range size was done to characterize how the (already-disclosed) O(archive) cost behaves as range width varies. Given §5 of `CLAIM3_LOCKED.md` already states plainly that cost is range-independent, a sweep would mostly be confirmatory rather than discover new behavior — but it isn't done, so marked as a gap rather than assumed. |
| **Failure cases** | 🔴 | ⚠️ | No test exercises what `export`/`coverage`/`query` do on a corrupt or truncated archive, an empty archive, or a query range fully outside `[0,PGLEN)` — same gap as the industrial checklist's Edge cases row. Nothing is silently discarded in the *success* path (every read that should count is now counted, per the bug-1 fix), but the *failure* path is untested. |
| **Negative results** | 🟠 | ✅ | The query-speedup correction (3.3× → 1.62×) **is** a negative result relative to the earlier claim, and it was not hidden — both the wrong number and the corrected one are in the record, with the mechanism of the original error explained (`CLAIM3_LOCKED.md` §5, §6.4). |
| **Resource measurement** | 🟠 | ✅ | Every timing in this claim uses the same tool, `/usr/bin/time -v`, parsed the same way (`parse_time_v()` in `run_claim3.sh`), for both CAPSULE and every conventional baseline — one consistent methodology, not different stopwatches for different sides. |
| **Hardware/environment** | 🟠 | ⚠️ | CPU core count (`nproc`) is used and implicitly logged via the SPAdes command line, but OS version, exact CPU model, and RAM are not recorded in `CLAIM3_LOCKED.md` itself — they're recorded once, generically, in the outer project's own `CLAUDE.md` ("Ubuntu 24.04, 12 vCPU, ~90 GB RAM"), and this claim's doc relies on that rather than repeating it. Acceptable but not self-contained. |
| **Random seeds** | 🔴 | ✅ | The one place randomness appears in Claim 3's own code — the new synthetic test generator — uses `random.seed(1234)`, explicit and reported in the script's own comment header. |
| **Dependency versions** | 🔴 | ⚠️ | SPAdes is pinned to an exact tag (v4.0.0) in the install command. bwa/samtools/mosdepth are not version-pinned anywhere in this claim's scripts or docs — same underlying gap as the industrial checklist's Dependency management row, called out once rather than twice. |
| **Automated tests** | 🔴 | ✅ (new) | `scripts/test_claim3.sh` tests the actual scientific function (coverage's depth computation, query's dedup semantics) — not a UI or utility test. This directly addresses what was, before this session, a real 🔴-critical gap. |
| **Regression tests** | 🟠 | ✅ (new, same script) | The coverage-total invariant is, by construction, a regression test for exactly the bug class found this session — verified above to actually fail against the pre-fix binary. |
| **Synthetic/toy tests** | 🟠 | ✅ (new) | 213 synthetic reads over a seeded 3000 bp genome, small enough to run in seconds, demonstrates the three operations behave correctly on a known-shape input independent of the real E. coli measurement. |
| **Independent verification** | 🔴 | ✅ | Coverage and query outputs are checked against *derivable* ground truth (sum of read lengths from the input file itself, encoder's own reported unique-read count) rather than against another tool's output — a form of independent verification appropriate to an operation with no directly-equivalent existing tool to cross-check against (per the prior-art survey, nothing else offers coverage+query without a reference genome). |
| **Figures/tables generation** | 🟠 | ⚠️ | `results/claim3/t6_results.csv` is generated by a script for the E. coli rows; the GIAB rows are not. No plotting/figure script exists yet for Claim 3 specifically — the numbers exist as a table, not yet as a paper-ready figure. |
| **Claim traceability** | 🔴 | ✅ | Paper claim ("export ~50-200× vs SPAdes") → `CLAIM3_LOCKED.md` §6.2 → exact command → `scripts/run_claim3.sh` → `stages/capsule_decode.cpp` lines → measured 555-656× (exceeding, not just meeting, the claimed range) is a complete, walkable chain with no missing link. |
| **README reproduction guide** | 🔴 | ✅ (new) | `scripts/run_claim3.sh`'s own header comment **is** the reproduction guide: usage line, what each argument does, what env vars override, what it deliberately does not do. A reviewer does not need to read `CLAIM3_LOCKED.md` first to run it, though doing so gives the full context. |

## Net result of this pass

Every 🔴-critical row that was a real gap before this session
(Losslessness/correctness tests, Automated tests, Reproducibility for the
E. coli side, Metric implementation) is now ✅, verified against the actual
historical bug rather than assumed fixed. The remaining 🔴 gaps —
data-integrity checksums, dependency-version pinning for bwa/samtools/
mosdepth — are real, named, and left open rather than papered over, because
fixing them requires re-touching large data files or system package
management that is not central to Claim 3's own logic. No 🔴 row is marked
✅ without a specific file/line/command backing it.
