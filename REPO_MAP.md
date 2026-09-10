# REPO_MAP — every folder, every subfolder, what's in it and why

Read this to find where something lives without guessing. For *what to read
first* to understand the project, that's `RECAP.md`, not this file — this
file is a map of the filesystem, not a narrative.

Root-level files are covered first, then each top-level directory,
subfolder by subfolder.

---

## Root files

| file | what it is |
|---|---|
| `README.md` | the general entry point — what the tool is, quickstart, top-level layout |
| `RECAP.md` | **the deep narrative** — full project history, every claim's mechanism, every bug found and fixed, standing rules, current frozen state. Read this for depth; README is deliberately general. |
| `REPO_MAP.md` | this file |
| `AUDIT.md` | the reproducibility audit's own record — what was verified by re-execution vs by transcription-check only, the one figure withdrawn, what the audit did NOT establish |
| `CLAUDE.md` | this repo's own technical ground-truth reference — current standing, what's done/not done, refuted ideas, commit map |
| `LICENSE` | MIT |
| `CITATION.cff` | machine-readable citation record (GitHub auto-parses this into APA/BibTeX) |
| `Makefile` | build entry point — `make` builds both `capsule_encode` and `capsule_decode` |
| `PIPELINE.md` | **command reference** — every subcommand, every flag that changes a published number (see its §1a), which script produces which table |
| `DATASET_LOCKED.md` | the OUTER ARCS project's 17-accession list — does **not** govern this repo, kept for reference only |
| `NEW_DATASET_LOCKED.md` | **this repo's actual locked set** — 15 non-human + 4 GIAB human = 19 datasets, reconciled against the published CSVs |
| `.gitignore` | excludes data files, compiled binaries, regenerable stream dumps |

---

## `paper/` — 13 files — the condensed, citable manuscript source

This is what a reader who wants the *result*, not the process, should read.
Several files here are symlinks into `docs/` (the canonical copy lives in
`docs/`, `paper/` is the reading-order view of it) — do not edit the symlink
target's content via the `paper/` path with a bulk tool like `sed -i`
without checking `git diff --raw` after (this broke twice during the
freeze; see `RECAP.md` §4 Phase 9).

| file | covers |
|---|---|
| `PAPER.md` | the whole submission in one place — abstract, the three claims, provenance |
| `METHODS.md` | the architecture as implemented, all three claims' mechanisms |
| `RESULTS.md` | every table, as executed |
| `LIMITATIONS.md` | what this work does not establish — read before citing anything as unconditional |
| `DISCUSSION.md` | the field, what's ours, what's not, six refuted hypotheses |
| `ARCHITECTURE.md` | the same system in diagrams |
| `NOVELTY_FINAL.md` (symlink) | the novelty argument against prior art read as primary sources |
| `CLAIMS_FINAL.md` (symlink) | the three claims stated to be quoted |
| `CLAIM3_MECHANISM.md` (symlink) | the locus-fragmentation measurement in full |
| `FORMAT.md` (symlink) | archive format |
| `PIPELINE.md` (symlink to root `PIPELINE.md`) | command reference |
| `ARCHITECTURE_VS_DISCOSNP.md` (symlink) | layer-by-layer against the closest caller |
| `COMPRESSION_DERIVED_CALLING.md` (symlink) | the channels that come free from compression, and why they still lose |

---

## `stages/` — 98 files — every experiment, one file per decision

Numbered `01` through `106`. This is the project's lab notebook made of
compilable code: each stage is the evidence for a decision, and several were
later contradicted by measurement (kept anyway, per this project's rule of
not deleting the record). **`106_inprocess.cpp` is the one that ships** —
it's the actual encoder behind `capsule_encode`. Everything below it is
history.

`capsule_decode.cpp` (also in this folder) is the decoder binary — it
implements `decompress`, `export`, `coverage`, `query`, `index`, and `call`
(Claim 2's entry point when calling from an existing archive rather than
inline during compression).

---

## `include/` — 7 files — the actual library

The code that both the standalone encoder and any future embedding would
share:

| file | what it is |
|---|---|
| `coders_inproc.h` | stream coders, the coder selector, transforms |
| `coders_pgrc.h` | PPMd7 / FSE / range coder wrappers and their decoders |
| `seqpar_core.h` | the DNA coder itself, shared so the standalone binary and the in-process path cannot diverge |
| `caps_caller.h` | **Claim 2** — the entire reference-free variant caller, 6,908 lines, zero external dependencies (verified) |
| `names_coder.h` | identifier/name stream coder |
| `quality_coder.h` | wraps vendored fqzcomp for the quality-score stream |
| `caps_pack.h` | archive container format read/write |

---

## `thirdparty/` — 43 files across 3 subfolders — vendored, licensed, disclosed

| subfolder | what's in it | license |
|---|---|---|
| `fse/` | FSE and Huf0 entropy coders (Yann Collet) | BSD |
| `ppmd/` | PPMd7 (from the LZMA SDK) | public domain |
| `htscodecs/` | fqzcomp quality-score codec (Genome Research Ltd, James Bonfield) | BSD 3-clause |

None of these are called as external tools at runtime — all compiled
directly into the shipped binaries. See `RECAP.md` §2 for why this doesn't
weaken the "no external tools" claim (same category as linking `libc`).

---

## `scripts/` — 56 files — build and benchmark automation

Includes `build106.sh` (the canonical build command — exists because the
`-fopenmp` flag once lived only in shell history and silently vanished,
costing 45% wall time), `encode_adaptive.sh` (the candidate-sweep encoder
wrapper), and the per-claim benchmark drivers
(`run_pipeline_bench.sh`, `run_claim1_bench.sh`, `run_fullchr20_bench_capsule.sh`,
`run_multiallelic_bench_capsule.sh`, `run_tetraploid_bench_capsule.sh`,
`run_claim3.sh`, etc.). `benchmark/documentation/SCRIPT_INVENTORY.md` marks
which of these 56 are the ~17 "live" ones still used versus historical.
`test_0_scope_and_capability.sh` is the self-contained proof of the tool's
stated scope — no locked dataset needed, run this first on any new machine.

---

## `benchmark/` — 68 files across 3 subfolders — the citable evidence

| subfolder | what's in it |
|---|---|
| `results/` | **THE citable numbers.** 8 CSVs, one run (2026-09-09), 19 datasets, 0 failures. If a document anywhere disagrees with a file here, this file wins. |
| `documentation/` | `RESULT_CODE.md` (every number traced to its script + output file), `REPRODUCE_EVERYTHING.md` (the full from-nothing reproduction procedure), `MANIFEST.sha256` (checksums of everything a published number depends on), `SCRIPT_INVENTORY.md`, `RESULTS_INVENTORY.md`, plus four dated `*_reproduction_20260910/` folders — the raw logs from the audit's independent re-executions (T2.1/T2.3, T2.4, T2.5, T3.4), and `T1.1_reverification_20260910/` — the E. coli re-verification that also surfaced the Genozip nondeterminism finding |
| `scripts/` | the top-level `benchmark.sh` dispatcher and its per-claim sub-scripts |

---

## `results/` — 41 files across 13 run-directories — working history, NOT citable

Has its own `results/README.md` explaining this explicitly. One exception:
`FULL_SWEEP_20260909/` is the original run that `benchmark/results/` is a
`cmp`-verified byte-identical copy of. Two directories are named to look
authoritative and are not — `claim2/` (development CSVs with superseded
scoring conventions) and `phase_a/` (decision evidence; its
`allphases_14dataset.csv` is VOID, pre-bug-fix numbers) — both explained in
`benchmark/documentation/RESULTS_INVENTORY.md`.

---

## `docs/` — 127 files across 3 subfolders — the full working record

The bulk of the repository by file count. ~29 of these files are directly
referenced from the citable path (`paper/`, `AUDIT.md`, `README.md`,
`CLAUDE.md`); the rest are the project's complete history — every
experiment, every plan, every refuted idea, every superseded draft. Kept
rather than deleted per this project's standing rule. Three subfolders:

| subfolder | what's in it |
|---|---|
| `_removable/` (17 files) | plan documents explicitly marked as not load-bearing for the current result — moved here rather than deleted, verified byte-identical to their pre-move blobs |
| `session_transcripts/` (3 files) | full exported chat transcripts from earlier sessions, kept as the working record of *how* conclusions were reached, not just what they concluded |
| `working_notes/` (7 files, moved here from the repo root during the 2026-09-11 structure cleanup) | `ARCHIVE_CALLER_SPEED_RAM.md`, `ARCHIVE_PATH_TO_72S.md` (superseded by the former), `CLAIM1_STREAM_ANALYSIS.md`, `FINAL_VERDICT.md`, `HANDOVER.md` (2026-09-02 project handover, historical), `RUNBOOK.md`, `SHIPPING_READINESS.md` — root-level working notes that cluttered the top-level listing without being part of the citable core |

The other ~100 files directly in `docs/` are the per-claim analysis history:
`CLAIM1_*`, `CLAIM2_*`, `CLAIM3_*` prefixed files trace each claim's
development; `docs/INDEX.md` is the original navigation aid for this
directory (itself partially superseded by `paper/` and this file).

---

## `server/` — 7 files — rebuild the machine from nothing

`verify_environment.sh` (27 checks, run this before downloading 55GB of
data), `SETUP_SHORT.md` / `SETUP_FULL.md` (the latter has the exact error
text and fix for every setup failure that has actually happened here),
`TOOLS.md` (exact versions and provenance of every external competitor tool
needed for the benchmark — SPRING, Genozip, DiscoSNP++, etc.), `DATASETS.md`
(per-file byte counts for the 19 locked datasets), `patch_spring.sh` (the
GCC-13-specific patch SPRING needs on Ubuntu 24.04).

---

## `industry/` — 2 files — the input-validation layer

`check_input_scope.py` — streams any FASTQ, extracts the encoder's own
structural bounds from source at run time, refuses out-of-scope input with
a stated reason (not a bare error). `README.md` — documents the two-layer
design (the encoder itself refuses unconditionally; this script reports the
same refusal in advance) and the real defects this design found.

---

## `.github/workflows/` — 1 file — CI

---

## The one-paragraph version, if you only remember this

**Code that ships:** `stages/106_inprocess.cpp` + `stages/capsule_decode.cpp`
+ `include/*.h` + `thirdparty/`. **Numbers that are citable:**
`benchmark/results/*.csv`, traced by `benchmark/documentation/RESULT_CODE.md`.
**The story:** `RECAP.md`. **The manuscript source:** `paper/*.md`.
**Everything else** (`docs/`, `results/`, `stages/01`-`105`, `scripts/`'s
non-live scripts) is real, honest working history — not required to trust
the result, kept because this project's standing rule is to preserve
evidence rather than delete it.
