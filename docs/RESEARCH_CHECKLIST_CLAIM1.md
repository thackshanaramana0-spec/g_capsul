# Research / academic checklist — Claim 1 (COMPACT)

> **SUPERSEDED 2026-09-10 — numbers only.** This document predates the final
> 19-dataset sweep of 2026-09-09, and it tracks the 14-dataset run and an unrun 15th dataset. Every one of those figures has been
> replaced:
>
> | | superseded | **final** |
> |---|---|---|
> | Claim 1 | 14 datasets, various margins | **19/19 vs SPRING (−6.03%), 19/19 vs Genozip (−43.26%), 57/57 LOSSLESS** |
> | het-SNV F1 | 0.890 | **0.876** (mean of 4 individuals, full chr20, called from the archive) |
> | het-indel F1 | 0.637 / 0.666 | **0.621** |
> | multi-allelic | 11/18 | **17/26**, a complete chr20 census (corrected 2026-09-10, was 21/26); DiscoSNP++ 0 of 3,989 records |
> | tetraploid SNV | 0.836 | **0.897** |
> | export vs SPAdes | 555–656x | **129–784x** |
> | T3.4 coordinate | 0/400 | **81/400** (the 0 was our query emitting consensus, not reads) |
>
> Citable source: `benchmark/results/*.csv`, traced in
> `benchmark/documentation/RESULT_CODE.md`. **Where this document and a result
> file disagree, the result file wins.**
>
> Kept because its *reasoning* is still the record of how the conclusion was
> reached — this repo does not delete superseded work, it marks it.

Scope: the 14/15-dataset SPRING/Genozip comparison and the separate 7-dataset
PgRC2 comparison. Statistical methodology marked N/A per project instruction
and precedent (established in the Claim 2/3 checklists; compression-ratio
comparisons in this field don't report significance either).

Status legend: 🔴 critical, 🟠 important. ✅ done · ⚠️ partial/scoped gap ·
❌ not done · N/A.

**Headline finding: this is the most evidence-complete of the three claims,
with one real, named, and now partially-addressed exception — the 15th
dataset was never run, and this pass does not run it (per instruction: the
verdict should not stop on this, but it stays visible, not hidden).**

| Area | 🔴/🟠 | Status | Finding |
|---|---|---|---|
| **Research question → code mapping** | 🔴 | ✅ | The claim (competitive lossless compression via pseudogenome assembly) maps directly to `stages/106_inprocess.cpp` — no ambiguity about what code produces this result. |
| **Reproducibility** | 🔴 | ✅ (14/14) / ❌ (15th) | `scripts/verify_lossless.sh` reproduces the lossless property on any input in one command. The size *comparison* against SPRING/Genozip covers 14 of the intended 15 datasets — Utricularia gibba (`SRR10676752`) was never downloaded/run, per `PROJECT_AUDIT.md`'s own gap list. This is a known, named gap, not a newly discovered one. |
| **Dataset provenance** | 🔴 | ✅ | `NEW_DATASET_LOCKED.md` documents not just the accession list but the *reasoning* for a mid-project swap (Drosophila out, Utricularia gibba in) — verified against NCBI eutils rather than assumed, and the rejected alternatives (Arabidopsis, other Plantae/Animalia candidates) are explained, not silently skipped. |
| **Data integrity** | 🔴 | ⚠️ | No checksums recorded for the 14 (soon 15) input FASTQ files — same gap already named in Claims 2 and 3. |
| **Train/test separation** | N/A | N/A | Not a learned model; no split applies. |
| **Ground truth** | 🔴 | N/A | Compression ratio is measured against archive size directly, not a truth set — the relevant integrity check is losslessness, which is covered (see Losslessness row). |
| **Experimental configuration** | 🔴 | ✅ | `MAXMAP`/`MINOV` sweep candidates, entropy coder selection, and every frozen threshold are named with their exact values in `CLAUDE.md` and the encoder's own comments. |
| **Deterministic execution** | 🔴 | ✅ | Same archive byte-for-byte on repeat runs — the precondition for every size number in this project, already established before this session. |
| **Baseline implementation** | 🔴 | ✅ | SPRING invoked with its documented `-g` decode flag (a real, previously-easy-to-miss requirement per `CLAUDE.md`'s own command reference), Genozip with `--force`, PgRC2 run from its own real GPL-3 source (never vendored, cloned separately) — not simulated or approximated. |
| **Fair comparison** | 🔴 | ✅ | Same input files, same machine, both tools measured fresh (`claim1_locked_result.md`: "both tools measured fresh on an idle machine") — the idle-machine qualifier itself shows awareness of a real confound (concurrent load skewing timing) even though timing isn't the headline metric for size. |
| **Metric implementation** | 🔴 | ✅ | Archive size is measured directly via `stat`/file size, not a derived computation with room for a bug — the lowest-risk metric of any claim in this project. |
| **Statistical methodology** | — | N/A | Intentionally not applied, consistent with the other two claims' checklists. |
| **Experiment independence** | 🔴 | ✅ | No hand-edited archive or intermediate file found in the pipeline between compress and measure. |
| **Raw → result pipeline** | 🔴 | ✅ | FASTQ → `encode_adaptive.sh` → archive → `stat` size + `verify_lossless.sh` → recorded number is fully traceable, no manual step. |
| **Losslessness/correctness tests** | 🔴 | ✅ | This is Claim 1's strongest row across all three checklists: not only does `verify_lossless.sh` exist, but the project's own history shows it being *distrusted and improved* — the discovery that it checked dumped streams rather than the real archive (`CLAUDE.md` §6.1) led to fixing the actual gap (`a81f55c`) rather than trusting a green check. That is exactly the skepticism this checklist row is trying to enforce, already demonstrated in this project's own past. |
| **Ablation experiments** | 🟠 | ✅ | Multiple structural changes (second-region self-match, PgRC2's both-side-overlap rule) were implemented, measured, and explicitly reverted with the losing numbers kept on record (`claim1_locked_result.md`'s "Structural attempts refuted" section) rather than silently dropped. |
| **Sensitivity analysis** | 🟠 | ✅ | MAXMAP and MINOV are both documented as swept, with the "interior optimum" property explicitly verified on the losing dataset (S. acidocaldarius) rather than assumed. |
| **Failure cases** | 🔴 | ✅ | The four silent-data-loss bugs (`CLAUDE.md` §6.3) were found precisely because the project went looking for failure cases (>256bp reads, orphaned unique reads, constant non-zero streams) rather than trusting a passing size-only benchmark — a genuinely strong precedent, not a gap. |
| **Negative results** | 🟠 | ✅ | S. acidocaldarius's -0.83% loss is reported plainly in the locked result, not omitted from the 7-dataset aggregate. |
| **Resource measurement** | 🟠 | ✅ | Compress/decompress time and peak RAM are both recorded with the same tool for both sides (`claim1_locked_result.md`'s side-by-side numbers), consistent with the methodology already established for Claims 2/3. |
| **Hardware/environment** | 🟠 | ⚠️ | Same as the other two checklists — relies on the outer `CLAUDE.md`'s one-time hardware description rather than restating it per result. |
| **Random seeds** | 🔴 | N/A | No randomness in the encoder's decision path to seed. |
| **Dependency versions** | 🔴 | ⚠️ | PgRC2's exact commit/version isn't pinned in the locked-result doc (cloned "separately," not vendored, per `CLAUDE.md`) — a real, if minor, reproducibility gap: a future PgRC2 change could shift the comparison without this project noticing. |
| **Automated tests** | 🔴 | ✅ | `verify_lossless.sh` predates this session and is real. |
| **Regression tests** | 🟠 | ✅ | Every output-preserving change in this project's history is required to be `cmp`-identical across all locked datasets before being accepted (`CLAUDE.md` §9 standing rule 2) — a regression-test discipline stronger than a single script, built into the workflow itself. |
| **Synthetic/toy tests** | 🟠 | ⚠️ | No small synthetic dataset exists purely for fast sanity-checking the encoder the way `test_claim2.sh`/`test_claim3.sh` now do for the caller and decoder — `verify_lossless.sh` always runs against a real (if possibly small) FASTQ file, not a generated one. |
| **Independent verification** | 🔴 | ✅ | Cross-validated against published gold values: SPRING bpb for SRR554369 and SRR870667 must match published figures within ±2% (`CLAUDE.md`'s own cross-validation rule) — an external anchor, not just internal consistency. |
| **Figures/tables generation** | 🟠 | ❌ → ⚠️ (partially fixed) | The master CSV was in `/tmp`, not scripted output committed to the repo — fixed this pass by committing `results/phase_a/allphases_14dataset.csv`. Still open: no script regenerates this CSV from a single command the way `scripts/run_claim3.sh` does for Claim 3; it was produced by an ad-hoc run, now merely rescued from `/tmp` rather than made re-runnable end to end. |
| **Claim traceability** | 🔴 | ✅ | "14/14 wins vs SPRING and Genozip, +1.88% vs PgRC2" traces directly to `results/phase_a/allphases_14dataset.csv` (now committed) and `claim1_locked_result.md`'s per-stream breakdown — a complete, walkable chain. |
| **README reproduction guide** | 🔴 | ✅ | `README.md`'s own "Build and run" + "Verify losslessness" sections are a genuine one-command reproduction guide, and predate this session. |

## Net result of this pass

Claim 1 needed the least new work of the three claims, because it already
had what the other two lacked before this session: a real lossless-verification
test, a documented sensitivity/ablation history, and external cross-validation
against published values. The one concrete action this pass took —
committing the previously `/tmp`-only results CSV — closes a real
provenance risk. **The one gap that remains open and is not fixed here, per
instruction, is the 15th dataset (Utricularia gibba) never having been
run** — named plainly, not silently dropped from the count, and not a
blocker to stating everything else in this checklist accurately.
