# Industrial checklist — Claim 1 (COMPACT)

Scope: the encoder (`stages/106_inprocess.cpp`, the shipped compressor),
`include/coders_inproc.h`/`coders_pgrc.h`/`seqpar_core.h`, the build/encode/
verify scripts, and the 14/15-dataset SPRING/Genozip/PgRC2 comparison.

Status legend: ✅ done · ⚠️ partial / real gap, scoped · ❌ not done · N/A.

**Claim 1 starts from a stronger position than Claims 2 and 3 did**: it
already had a real, working lossless-verification script
(`scripts/verify_lossless.sh`) before this pass, unlike the caller and the
decoder, which had zero tests until this session. The gaps here are
narrower, and one concrete fix landed: the master results CSV the paper
draft itself cites was sitting only in `/tmp/allph/r.csv` — gone on
reboot — and is now committed at `results/phase_a/allphases_14dataset.csv`.

| Area | Status | Finding |
|---|---|---|
| **Repository structure** | ✅ | `stages/` (progression), `include/` (coders), `scripts/` (build/encode/verify), `docs/`, `results/`, `thirdparty/` — the structure `CLAUDE.md` itself documents and that the rest of this project's checklists already relied on. |
| **Architecture** | ⚠️ | `stages/106_inprocess.cpp` is the shipped encoder and is large (comparable in scale to `caps_caller.h`); it's also the one file whose correctness everything else in this repo depends on (Claims 2 and 3 both read its output streams). Not split into modules, though the stage-by-stage `stages/` progression *is* a form of architectural history that a from-scratch rewrite would lose. |
| **Entry points** | ✅ | `scripts/build106.sh` + direct binary invocation is a real, documented CLI, not "edit source and recompile." |
| **Configuration** | ✅ | `ENC_ARGS`, `MAXMAP`/`MINOV` are all environment/argument driven (`encode_adaptive.sh`), not hardcoded. The frozen entropy-coder thresholds are a deliberate, documented choice (ported and measured, not guessed), same standing rule as the caller's frozen constants. |
| **Reproducibility** | ✅ | `scripts/verify_lossless.sh reads.fq` is a genuine one-command reproduction: encode, decode, diff against the original file's own sequence column — this is the strongest reproducibility story of the three claims, and it predates this session. |
| **Dependency management** | ⚠️ | Same repo-wide gap already named in the Claim 2/3 checklists: no `requirements.txt`/lockfile. The encoder's own dependencies (LZMA SDK, FSE, htscodecs) are vendored in `thirdparty/`, which sidesteps most of this for the C++ core — the real gap is competitor tool versions (SPRING, Genozip, PgRC2) not pinned anywhere. |
| **Determinism** | ✅ | Same archive, byte-identical output run to run — this is asserted throughout this project's own history (`GSEARCH` sweep is described as reproducing the same winner deterministically) and is the precondition every other size claim depends on. |
| **Testing** | ✅ | `scripts/verify_lossless.sh` already existed and is real: it decodes the DUMPED streams and diffs against the original file, not a coder-level round-trip. Its own documented limitation (§6.1 of `CLAUDE.md`, historical) — that it doesn't decode the *archive* itself — was found and fixed in an earlier session (`a81f55c`), so this gap is closed, not open. |
| **Edge cases** | ⚠️ | The project's own history (`CLAUDE.md` §6.3) documents four real silent-data-loss bugs found and fixed by testing edge cases that existed in real data (reads >256bp, orphaned unique reads, RC-contained reads, constant non-zero streams) — a genuinely strong track record. What's *not* covered by an automated test today: an empty input file, a single-read file, or a file with zero duplicate reads (the S. acidocaldarius near-loss case suggests low-duplication inputs are a real, sensitive edge, not a hypothetical one). |
| **Validation** | ⚠️ | No explicit check that an input file is well-formed FASTQ before encoding begins — malformed input behavior is undocumented. |
| **Error handling** | ⚠️ | Not re-audited at Claim 3's depth in this pass — `106_inprocess.cpp` is large and a dedicated pass was out of scope here. |
| **Logging** | ✅ | Per-stage timing (`[stage] ...`) and per-stream size reporting (`[archive] ...`) are both structured and already used throughout this session's own debugging — not scattered prints. |
| **Code quality** | ⚠️ | Consistent with the project's documented style (every non-obvious constant commented with provenance). A full quality pass on `106_inprocess.cpp` at Claim 3's depth was not performed here. |
| **Documentation in code** | ✅ | The project's single strongest documentation habit — refuted ideas are commented in place with their measured cost, not deleted (`CLAUDE.md` §7's table is a direct index into these comments). |
| **Dead/debug code** | ✅ | Same convention as the caller: refuted features live behind env-var flags (`SECOND_SELF`, `BOTHSIDE`, `GSEARCH`), documented as intentionally kept, not abandoned. |
| **Security** | ✅ | No credentials in encoder code or scripts. |
| **File handling** | ✅ | `encode_adaptive.sh` and `build106.sh` take paths as arguments with sensible defaults (`/tmp/best106`), not hardcoded machine-specific locations. |
| **Resource handling** | ⚠️ | `encode_adaptive.sh` runs multiple MAXMAP/MINOV candidates and discards all but the winner — by design (documented as the accepted cost of avoiding a fitted formula), but it does mean 3 of 4 candidate archives are computed and thrown away rather than the search being incremental. |
| **Performance** | ✅ | Explicitly measured, not assumed: `[stage]` timers per phase, and the `-fopenmp` build-flag story (`build106.sh`'s own header comment) is a real, documented performance fix (-45% wall time) that was previously invisible because the build command only lived in shell history. |
| **Scalability** | ✅ | Measured on real datasets from ~30MB (SARS-CoV-2) to ~15GB (T. cacao, per `DATASET_LOCKED.md`'s suppressed-autochunk tier) — the widest scale range of any of the three claims' evidence. |
| **Interfaces** | ✅ | `capsule compress`/`capsule decompress`-equivalent CLI unchanged across the whole stage progression's history. |
| **Backward compatibility** | N/A | No external consumers of the archive format outside this project yet. |
| **Data provenance** | ✅ | `DATASET_LOCKED.md` and `NEW_DATASET_LOCKED.md` both name exact accessions, coverage, and — critically — document *why* a swap was made (Drosophila out, Utricularia gibba in, verified live via NCBI eutils, not guessed) rather than silently substituting. |
| **Result provenance** | ⚠️ → ✅ (fixed this pass) | The master results table the paper draft (`PAPER_DRAFT_CLAIM1.md`) cites lived only at `/tmp/allph/r.csv` — outside the repo, in a location that doesn't survive a reboot and isn't version-controlled. **Fixed in this pass**: copied to `results/phase_a/allphases_14dataset.csv`, now tracked. |
| **CI** | ❌ | Same repo-wide gap as Claims 2 and 3 — no `.github/workflows` anywhere. |
| **Versioning** | ✅ | The commit history for this encoder is exceptionally granular and self-documenting (e.g. `3e06957`'s own commit message explains the exact mathematical reason a duplicate can only be detected at `L=Lmax`) — a high bar, already met. |
| **Licensing** | ⚠️ | Same gap already fixed for Claim 3 (`thirdparty/fse/LICENSE` added) applies identically here since it's the same vendored FSE code used by the encoder — already resolved, no new action needed. The repo's own top-level license remains an open, project-wide decision, not fixed in any of the three claim passes. |
| **Maintainability** | ✅ | The stage-by-stage `stages/01...106` progression is unusual but serves exactly the maintainability goal: each numbered stage is evidence for one decision, so a developer can trace *why* the shipped encoder looks the way it does, not just *what* it does. |

## Net result of this pass

One concrete fix: the at-risk results CSV is now committed
(`results/phase_a/allphases_14dataset.csv`), closing a real
reproducibility/provenance gap before it caused data loss. Everything else
in this checklist was already in stronger shape than Claims 2/3 started
in — Claim 1 had a genuine round-trip test before this session even began.
The open items (CI, top-level license, a full line-by-line audit of
`106_inprocess.cpp` at Claim 3's depth) are named, not new discoveries that
change the claim's standing.
