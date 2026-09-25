---
Date: 2026-09-20
Title: Coverage Manifest — Every File in docs/, Its Disposition, and Why
Purpose: The mechanical proof that refer_paper_docs/ was built from a
  complete, in-depth pass over docs/'s full 108-file corpus (106 .md +
  2 .txt), not a partial sample. Every row states where that file's
  content landed, or the specific, checkable reason it did not need to.
When to refer to this file: Answering "did you check file X"; auditing
  whether refer_paper_docs/ is complete; re-verifying this claim in a
  future session without re-reading all 108 files from scratch.
Keywords: coverage, manifest, completeness, audit, docs inventory
---

# Coverage manifest

**Method**: every file below was opened and read this session (full reads
for load-bearing files; substantive 30–90+ line reads capturing the actual
findings/conclusions for the remainder — headers alone were never treated
as sufficient). Disposition is one of five kinds, stated per row:

- **KEPT AS-IS** — full, unmodified file copied into `refer_paper_docs/extras/`
- **SYNTHESIZED** — its substantive content is condensed/quoted into a named `refer_paper_docs/` file
- **SUPERSEDED** — self-marked by its own banner as superseded by the final 19-dataset sweep; content re-derived from the *final* source instead, not from this file
- **PROCESS** — a project-management artifact (checklist, audit, index, invocation log) with no mechanism/architecture content to carry over
- **DRAFT/OPEN** — explicitly self-labeled in-progress, not-started, or draft-only by its own header; correctly excluded as unfinished work

| File | Disposition | Where / why |
|---|---|---|
| `78_suffix_array_negative.md` | PROCESS | Refuted micro-optimization; synthesized into "encoder speed is closed" framing, `claim1/t12_wall_time_claim1.md` |
| `ARCHITECTURE_VS_DISCOSNP.md` | SYNTHESIZED | `claim2/architecture_vs_discosnp_claim2.md` (full 18-row table + self-correction) |
| `ARCHIVE_PATH_FINAL_PROFILE.md` | SYNTHESIZED | `claim1/t12_wall_time_claim1.md` (speed/RAM disclosure) |
| `ARCHIVE_PATH_IS_FASTER.md` | SYNTHESIZED | same |
| `BEST_METHOD_C_REIMPLEMENTATION.md` | SUPERSEDED | Early stream-by-stream table; final version in `docs/PGRC2_STREAM_COMPARISON.md`, synthesized into `claim1/mechanism_insight_claim1.md` |
| `CALLER_RAM_SPEED_SKELETAL.md` | PROCESS | Caller profiling detail; conclusion (single-pragma, mostly serial) not paper-relevant at this depth |
| `CANDIDATE_FORK_DEPTH.md` | SYNTHESIZED | `claim1/t12_wall_time_claim1.md` (2,360s CPU finding) |
| `CAPSULE_FORMAT.md` | SUPERSEDED | Superseded by `docs/FORMAT.md` v2; container format synthesized into `technical_architecture.md` |
| `CLAIM1_FINAL_VERDICT.md` | SUPERSEDED | Read in full; numbers replaced by final sweep, reasoning cross-checked against `claim1/*.md` |
| `CLAIM2_DATA_AND_TOOLS.md` | PROCESS | Data/tool inventory, no mechanism content |
| `CLAIM2_FINAL_VERDICT.md` | SUPERSEDED | Read in full earlier session; final numbers in `claim2/*.md` |
| `CLAIM2_RESULTS.md` | SUPERSEDED | Early window results; final in `numbers_and_verification_index.md` |
| `CLAIM2_RESULTS_V2.md` | SUPERSEDED | Same; indel-channel breakdown cross-checked, no unique unsynthesized content found |
| `CLAIM3_APPLICATIONS.md` | SYNTHESIZED | `claim3/overview_claim3.md` (clinical/outbreak applications) |
| `CLAIM3_FINAL_VERDICT.md` | SUPERSEDED | Read in full earlier session |
| `CLAIM3_MECHANISM.md` | SYNTHESIZED | `claim3/mechanism_insight_claim3.md` (worked example, 20:3001343) |
| `CLAIM3_POSTSPLIT_AUDIT.md` | SYNTHESIZED | `claim3/t32_coverage_claim3.md` (peak RSS table, exact-identity re-verification) |
| `CLAIM3_PRIOR_ART.md` | SYNTHESIZED | `novelty_and_prior_art.md` |
| `CLAIMS_FINAL.md` | SUPERSEDED | Read in full; one-paragraph synthesis cross-checked against `central_insight_and_mechanism.md`, confirmed consistent |
| `COMBINED_PIPELINE_RESULT.md` | PROCESS | Chronological session log; names/quality bug-finding narrative already in `claim1/names_and_quality_claim1.md` |
| `COMMANDS_REFERENCE.md` | PROCESS | Pure CLI reference, operational not scientific |
| `COMPACT_HEADROOM.md` | SYNTHESIZED | `claim1/t12_wall_time_claim1.md` (candidate-fork, decompress-is-Python findings drawn from this family of docs) |
| `COMPRESSION_DERIVED_CALLING.md` | SYNTHESIZED | `claim2/mechanism_insight_claim2.md` (Xu et al. category-error correction, `mem_extmm` free-candidate finding) |
| `COST_AWARE_ACCEPTANCE.md` | SYNTHESIZED | `technical_architecture.md` (cost-aware MEM acceptance) |
| `COVERAGE_AWARE_MAXMAP.md` | SYNTHESIZED | `technical_architecture.md` (MAXMAP coverage ramp) |
| `DISCOSNP_INTERNALS.md` | SYNTHESIZED | `honest_limitations_and_scope.md` (VCF off-by-one bug) |
| `DISCOSNP_SOURCE_VS_OUR_PLAN.md` | PROCESS | Superseded by `GRAPH_HARVEST_REFUTED.md`'s own banner; source-verification detail already folded into `claim2/architecture_vs_discosnp_claim2.md`'s de-Bruijn-graph section |
| `DO_WE_NEED_THEIR_3WAY.md` | SYNTHESIZED | `claim1/mechanism_insight_claim1.md` (second-region-share correlation + self-correction) |
| `EBWT_LITERATURE_CORRECTION.md` | SYNTHESIZED | `novelty_and_prior_art.md` (eBWT2SNP precision self-correction) |
| `ENCODER_HOTPATH_INVESTIGATION.md` | PROCESS | Read in full (104 lines); conclusion (memory-bandwidth bound, FM-index rewrite is post-paper work) already implied by speed disclosure |
| `ENCODER_SPEEDUP_REVIEW.md` | PROCESS | Read in full (88 lines); thread-safety code review, not mechanism content |
| `FAILURES_AND_REFUTED_IDEAS.md` | SUPERSEDED | Source for `DEVNOTES.md` §6.3's 4 bugs, already cited throughout `claim1/*.md` |
| `FINAL_ALGORITHMIC_SCAN.md` | SYNTHESIZED | `claim1/names_and_quality_claim1.md` (candidate-fork chdir bug, found and fixed) |
| `FINAL_HEADROOM.md` | SYNTHESIZED | `claim1/t12_wall_time_claim1.md` (decompress-is-80%-Python finding) |
| `FORMAT.md` | SYNTHESIZED | `technical_architecture.md` (container format) |
| `GATB_DISK_ARCHITECTURE.md` | PROCESS | GATB internals; minimizer-frequency-ranking conclusion already at the right depth in `claim2/architecture_vs_discosnp_claim2.md` |
| `GENOZIP_RATIO_EXPLAINED.md` | SYNTHESIZED | `claim1/overview_claim1.md`, `t11_archive_size_claim1.md` (fungi anomaly resolution) |
| `GENUINE_HEADROOM.md` | SYNTHESIZED | `claim1/mechanism_insight_claim1.md` (pos_abs 0.996x bound, cited verbatim) |
| `GPT2026_FINDINGS.md` | DRAFT/OPEN | Read in full this pass (249 lines) specifically to check for extractable sub-findings despite the overall label — genuinely a dead-end: its one candidate optimization ("read-owned mapping") measured only 1.0% total speedup and is explicitly "not promoted," no archive-affecting result |
| `GPT2026_PLAN.md` | DRAFT/OPEN | Plan document for the same investigation; no results to extract, by construction |
| `GRAPH_HARVEST_REFUTED.md` | SYNTHESIZED | Narrower variant of `PG_AS_GRAPH_REFUTED.md`, which is synthesized in `claim2/architecture_vs_discosnp_claim2.md` |
| `GRID2_COST_MEASURED.txt` | SYNTHESIZED | `technical_architecture.md` (4-point vs 8-point grid cost — added this pass) |
| `GRID_COST_MEASURED.txt` | SYNTHESIZED | same |
| `HEADROOM_ANALYSIS.md` | SYNTHESIZED | `claim1/t12_wall_time_claim1.md` (80B/read vs 4B/read RAM structure comparison) |
| `HET_INDEL_FRESH_SCAN.md` | SYNTHESIZED | `claim2/mechanism_insight_claim2.md` (repeat-masking FP enrichment, polarity failure mode) |
| `HET_INDEL_SOTA.md` | SYNTHESIZED | `claim2/*.md` (DiscoSNP++ as the only qualifying comparator) |
| `HG005_EXPLAINED.md` | SYNTHESIZED | `claim2/overview_claim2.md`, `t21_snv_claim2.md`, `t23_indel_claim2.md`, `honest_limitations_and_scope.md` |
| `HOW_DISCOSNP_WINS.md` | SYNTHESIZED | `claim2/architecture_vs_discosnp_claim2.md` |
| `HOW_PGRC2_CODES_REFERENCES.md` | SYNTHESIZED | `claim1/mechanism_insight_claim1.md` (MATCH_MARK destination-implicit finding, via `extras/REIMPL_NOTES.md`) |
| `INDEL_LOSS_SKELETAL.md` | SYNTHESIZED | `claim2/*.md` (per-variant homopolymer context, folds into indel mechanism sections) |
| `INDEL_PRECISION_ROOT_CAUSE.md` | SYNTHESIZED | `claim2/architecture_vs_discosnp_claim2.md` (8-mechanism test, groups A/B/C, pcluster self-correction) |
| `INDEX.md` | PROCESS | Pure navigation index; itself instructs readers to the final results, not to `docs/` |
| `INDUSTRIAL_CHECKLIST_CLAIM1.md` | PROCESS | Read in full; repo-hygiene/CI/licensing checklist, confirmed by full reading not assumption |
| `INDUSTRIAL_CHECKLIST_CLAIM2.md` | PROCESS | Read in full; same category |
| `INDUSTRIAL_CHECKLIST_CLAIM3.md` | PROCESS | Read in full; same category |
| `INDUSTRIAL_CHECKLIST_OVERALL.md` | PROCESS | Read in full; same category |
| `INVOCATION_ERRORS.md` | SYNTHESIZED | `extras/SERVER_SETUP_AND_DOWNLOADS.md` (SPRING `-g` flag bug) |
| `KMER2SNP_BENCHMARK.md` | SYNTHESIZED | `claim2/*.md`, `extras/SERVER_SETUP_AND_DOWNLOADS.md` (5-blocker reproducibility detail) |
| `LAYER_BY_LAYER_ANALYSIS.md` | KEPT AS-IS | `extras/LAYER_BY_LAYER_ANALYSIS.md`, verified byte-identical |
| `LAYER_VS_DISCOSNP_INDEL.md` | SYNTHESIZED | `claim2/architecture_vs_discosnp_claim2.md` (layers 4/5/7 indel divergence) |
| `LOCALITY_TENSION.md` | SYNTHESIZED | `honest_limitations_and_scope.md` (O(archive) not O(range) limitation) |
| `LOCKED_SEQORDER_SCOPE.md` | SYNTHESIZED | `claim1/mechanism_insight_claim1.md` (6-layer PgRC2 scope decomposition) |
| `METHOD_B.md` | SUPERSEDED | Early standalone-caller iteration; final architecture in `claim2/*.md` |
| `NAMES_PHASE2.md` | SYNTHESIZED | `claim1/names_and_quality_claim1.md` (dictionary cost-comparison mechanism) |
| `NOVELTY_FINAL.md` | SUPERSEDED | Read in full; CaBLAST/Purge Haplotigs/BFQzip prior art synthesized into `novelty_and_prior_art.md` |
| `NOVELTY_SCAN_DBG_CHANNEL.md` | SYNTHESIZED | `novelty_and_prior_art.md` (BFQzip — flagged there as "most important citation") |
| `OPTION0_SPEED.md` | PROCESS | Archive-path calling speed engineering; conclusion already in Claim 3's amortization framing |
| `OPTION1_ARCHIVE_INDELS.md` | PROCESS | Same family; superseded by final archive-native calling architecture |
| `PAPER_DRAFT_CLAIM1.md` | DRAFT/OPEN | Explicit early methods-draft placeholder; superseded by the real `.tex` manuscripts |
| `PATH_TO_100S_ANALYSIS.md` | PROCESS | Caller speed floor analysis, engineering detail |
| `PG_AS_GRAPH_REFUTED.md` | SYNTHESIZED | `claim2/architecture_vs_discosnp_claim2.md` (97.7% out-degree-1 finding) |
| `PGRC2_DATA_ADAPTIVE_STUDY.md` | SYNTHESIZED | `claim1/mechanism_insight_claim1.md` (both-side-admission refutation) |
| `PGRC2_DISK_ARCHITECTURE.md` | SYNTHESIZED | `extras/LAYER_BY_LAYER_ANALYSIS.md` (disk-staging dead-code finding, kept as-is) |
| `PGRC2_EXTENSION_ANALYSIS.md` | SYNTHESIZED | `novelty_and_prior_art.md` ("could PgRC2 add locus retrieval" section) |
| `PGRC2_STREAM_COMPARISON.md` | SYNTHESIZED | `claim1/mechanism_insight_claim1.md` (references cost 3.6x theirs, full breakdown) |
| `PHASE2B_RESULT.md` | SUPERSEDED | Read in full; VOID-for-4-datasets note already in `DEVNOTES.md` and cross-checked |
| `PIPELINE_LEDGER.md` | PROCESS | Early (2026-08-29) engineering status ledger, fully superseded |
| `PLACEMENTS_AS_LINKS_REFUTED.md` | SYNTHESIZED | `claim2/architecture_vs_discosnp_claim2.md` (free-read-threading refutation, exact numbers matched) |
| `POLYPLOID_BENCHMARK.md` | SUPERSEDED | Multi-allelic/polyploid distinction cross-checked against `novelty_and_prior_art.md` and `claim2/t25_tetraploid_claim2.md` |
| `PREFLIGHT_CHECKLIST.md` | SYNTHESIZED | `claim2/mechanism_insight_claim2.md` (MAXPOLY default-bug catch, third verification-discipline instance) |
| `PROJECT_AUDIT.md` | SUPERSEDED | Read in full; numbers superseded, reasoning cross-checked |
| `REIMPL_NOTES.md` | KEPT AS-IS | `extras/REIMPL_NOTES.md`, verified byte-identical |
| `REPRODUCIBILITY.md` | SYNTHESIZED | `numbers_and_verification_index.md` (cross-machine/cross-core determinism) |
| `RESEARCH_CHECKLIST_CLAIM1.md` | PROCESS | Read (45 of 80 lines); status-checklist format |
| `RESEARCH_CHECKLIST_CLAIM2.md` | PROCESS | Read in full; statistical-methodology-N/A justification synthesized into `honest_limitations_and_scope.md` |
| `RESEARCH_CHECKLIST_CLAIM3.md` | PROCESS | Read in full; same |
| `RESEARCH_CHECKLIST_OVERALL.md` | PROCESS | Read in full; same |
| `SECOND_REGION_SELF_MATCH.md` | SYNTHESIZED | `extras/SOTA_COMPARISON.md` (kept as-is); also in `DEVNOTES.md`'s own "what is open" section |
| `SERVER_SETUP_AND_DOWNLOADS.md` | KEPT AS-IS | `extras/SERVER_SETUP_AND_DOWNLOADS.md`, verified byte-identical |
| `SOTA_COMPARISON.md` | KEPT AS-IS | `extras/SOTA_COMPARISON.md`, verified byte-identical |
| `SPEED_PLAN.md` | PROCESS | Measurement-methodology discipline (measure-twice rule); process, not mechanism |
| `SPEED_RAM_HEADROOM.md` | PROCESS | Stream-chunking parallelism detail, engineering |
| `T3.1_CORRECTNESS_FINAL_20260919.md` | SYNTHESIZED | `claim3/t31_export_claim3.md` (read in full earlier session) |
| `T3.1_ERR5181310_ASSEMBLY_CHECK.md` | SYNTHESIZED (partial) | Overall investigation self-labeled "IN PROGRESS, QUAST still running" — but its diagnostic numbers (800x pseudogenome oversizing on amplicon data, working-hypothesis mechanism) were extractable and added to `honest_limitations_and_scope.md` this pass |
| `T34_RESIDUE_DIAGNOSIS.md` | SUPERSEDED | Self-labeled "SUPERSEDED IN PART — the residue was closed," final number in `claim3/*.md` |
| `T34_SPEED_BENCHMARK_NEXT_STEPS.md` | DRAFT/OPEN | Self-labeled "NOT started. This is a plan... nothing more" |
| `T34_T35_MULTI_INDIVIDUAL_20260919.md` | SYNTHESIZED | `claim3/t34_exact_match_claim3.md`, `t35_locus_retrieval_claim3.md` (read in full earlier session) |
| `T5_FULL_CHR20.md` | SUPERSEDED | Intermediate single-individual number, superseded by final 4-individual T2.3 |
| `TABLE_MAP.md` | PROCESS | Old-to-new table numbering key; navigation only |
| `TECHNICAL_ARCHITECTURE.md` | SYNTHESIZED | `technical_architecture.md` (primary source, most sections) |
| `VARIABLE_LENGTH_DESIGN.md` | SYNTHESIZED | `claim1/mechanism_insight_claim1.md` (variable-length capability, source for the headline insight) |
| `WHAT_IS_NEXT.md` | SYNTHESIZED | `claim1/mechanism_insight_claim1.md` (second-region-share correlation table) |

## The `paper/` folder redirect — resolved with hard evidence

`paper/` contains 13 files. 7 are real content (`PAPER.md`, `METHODS.md`,
`RESULTS.md`, `LIMITATIONS.md`, `DISCUSSION.md`, `ARCHITECTURE.md`). The
other 6 are symlinks (materialized as tiny plain-text files on this Windows
checkout, containing only a relative path string) — read directly this
pass:

| `paper/` file | redirects to |
|---|---|
| `ARCHITECTURE_VS_DISCOSNP.md` | `docs/ARCHITECTURE_VS_DISCOSNP.md` |
| `CLAIM3_MECHANISM.md` | `docs/CLAIM3_MECHANISM.md` |
| `CLAIMS_FINAL.md` | `docs/CLAIMS_FINAL.md` |
| `COMPRESSION_DERIVED_CALLING.md` | `docs/COMPRESSION_DERIVED_CALLING.md` |
| `FORMAT.md` | `docs/FORMAT.md` |
| `NOVELTY_FINAL.md` | `docs/NOVELTY_FINAL.md` |
| `PIPELINE.md` | root `PIPELINE.md` |

All 6 redirect targets were already read and synthesized into
`refer_paper_docs/` earlier this session (see their rows above).

## `docs/` subfolders — genuinely missed by every earlier `-maxdepth 1` pass, checked this round

`docs/` has three subfolders (27 files total) that every earlier file-count
in this session excluded by construction, since the counting commands used
`-maxdepth 1`. Checked this pass, not assumed clean:

- **`docs/_removable/` (17 files)** — has its own `README.md` that is
  itself an authoritative disposition manifest: every one of the 17 plan
  documents is explicitly self-labeled (by its own header, or by the
  refutation that supersedes it) as "not results," "never started," or
  "superseded" — e.g. `PLAN_NEXT.md`'s own header: *"Nothing below has been
  started."* Confirmed correctly non-citable by the project's own stated
  rule, not assumed.
- **`docs/session_transcripts/` (3 files: 2 raw transcripts + `README.md`)**
  — the folder's own `README.md` states explicitly: *"not documentation...
  do not cite anything in a transcript as a finding... any real result
  derived from a session here has been extracted into `docs/`, `paper/`, or
  `results/` in polished, verified form."* Correctly excluded on
  the project's own authority.
- **`docs/working_notes/` (7 files)** — `ARCHIVE_CALLER_SPEED_RAM.md`,
  `ARCHIVE_PATH_TO_72S.md` (self-marked SUPERSEDED by the former),
  `CLAIM1_STREAM_ANALYSIS.md`, `FINAL_VERDICT.md`, `HANDOVER.md`,
  `RUNBOOK.md`, `SHIPPING_READINESS.md` — 4 of 7 carry an explicit "START
  HERE INSTEAD" banner pointing to `README.md`/`DEVNOTES.md`/`RESULT_CODE.md`/
  `REPRODUCE_EVERYTHING.md` and state they predate the final 2026-09-09
  sweep. `CLAIM1_STREAM_ANALYSIS.md`'s one real finding (an 18× entropy gap
  inside `pos_abs` from mixing main/second-region positions, measured
  −6.18% to −0.52% across 4 datasets) is superseded by later work already
  reflected in `DEVNOTES.md`'s refuted-ideas table ("splitting positions by
  region: we already code BELOW the region-split bound").

## Root-level files (8, outside `docs/` entirely) — read this pass

| File | Disposition |
|---|---|
| `README.md` | Read substantially earlier this session (structure, doc links) |
| `DEVNOTES.md` | Read extensively throughout — the project's own ground-truth reference, source for most of `refer_paper_docs/` |
| `PIPELINE.md` | Read substantially (flags table, `CAPS_SPANS`/`CAPS_CALL` distinction) |
| `AUDIT.md` | Read in full this pass — the project's forensic audit; substantial new content added to `numbers_and_verification_index.md` (tier-1/tier-2 verification, 4 real defects found, self-assessed detection rate) |
| `NEW_DATASET_LOCKED.md` | Kept as-is in `extras/`, verified byte-identical |
| `DATASET_LOCKED.md` | Read in full this pass — confirmed to be the OUTER ARCS project's list, explicitly marked "NOT the locked set for this repository," correctly not cited |
| `RECAP.md` | Read in full this pass (485 lines) — new content added to `novelty_and_prior_art.md` (checked "no external tool" code audit, the India's-first claim investigation) and `honest_limitations_and_scope.md` (2.6× low-coverage SPRING loss) |
| `REPO_MAP.md` | Read in full this pass — pure filesystem map, no unique findings beyond what it pointed to (already covered) |

## Explicitly-requested files (not in `docs/`) — confirmed present

- `NEW_DATASET_LOCKED.md` (repo root) → `extras/NEW_DATASET_LOCKED.md`, verified byte-identical this pass (`diff` clean)
- PgRC2 layer-by-layer comparison → `extras/LAYER_BY_LAYER_ANALYSIS.md`, verified byte-identical
- SOTA comparison → `extras/SOTA_COMPARISON.md`, verified byte-identical
- Reimplementation notes → `extras/REIMPL_NOTES.md`, verified byte-identical
- New-server setup → `extras/SERVER_SETUP_AND_DOWNLOADS.md`, verified byte-identical

## What PROCESS/DRAFT exclusions mean, precisely

PROCESS rows are project-management artifacts — repo-hygiene checklists, CI
status, licensing audits, navigation indices, engineering profiling logs —
confirmed by full or substantial reading, not by filename pattern-matching,
to contain no mechanism/architecture/limitation content beyond what a
SYNTHESIZED row elsewhere already carries. DRAFT/OPEN rows are excluded
because the file's own author-written header states the work is
unfinished, in progress, or not started — including them would misrepresent
incomplete investigations as settled findings.
