# Documentation index

114 documents. Per this project's own standing rule, nothing here is deleted
when superseded — retractions and refuted ideas stay on record, marked in
place. This index exists so a new reader can navigate that history without
reading all 114 files.

> ## READ THIS BEFORE ANY DOCUMENT IN THIS FOLDER
>
> **`docs/` is the project's HISTORY, not its results.** Almost everything here
> was written before the final 19-dataset sweep of 2026-09-09, and the numbers
> in it have been superseded. The current results live outside this folder:
>
> | you want | go to |
> |---|---|
> | the results | `../benchmark/results/` — 8 CSVs, one run |
> | where a number came from | `../benchmark/documentation/RESULT_CODE.md` |
> | to re-run it | `../benchmark/documentation/REPRODUCE_EVERYTHING.md` |
> | what the project is | `../README.md` |
> | the rules + refuted ideas | `../CLAUDE.md` |
> | the methods, for the paper | `../paper/METHODS.md` |
>
> **Where a document here and a result file disagree, the result file wins.**
>
> The final headline, so no one has to reconstruct it from the files below:
> Claim 1 **19/19 vs SPRING (−6.03%), 19/19 vs Genozip (−43.26%), 57/57
> LOSSLESS**; Claim 2 het-SNV **0.876** vs 0.853 vs 0.475, multi-allelic
> **21/26** vs **0/26**, tetraploid SNV **0.897** vs 0.782, all called from the
> archive; Claim 3 export **129–784x**, coverage **16–54x**, locus retrieval
> **345/400** by content vs **81/400** by coordinate.
>
> Documents whose numbers are stale now carry their own SUPERSEDED banner. Two
> retractions in particular: **T3.4 coordinate is 81/400, not 0/400** (the 0 was
> our query emitting the consensus rather than the reads), and **T2.4's "5/111"
> and "11/18" are both withdrawn** in favour of 21/26.

## Start here

| doc | what it is |
|---|---|
| [`../README.md`](../README.md) | product overview, results, quick start |
| [`TECHNICAL_ARCHITECTURE.md`](TECHNICAL_ARCHITECTURE.md) | the complete architecture — assembly, streams, coders, caller, decoder |
| [`COMMANDS_REFERENCE.md`](COMMANDS_REFERENCE.md) | every command, organized by claim |
| [`SERVER_SETUP_AND_DOWNLOADS.md`](SERVER_SETUP_AND_DOWNLOADS.md) | bootstrap a fresh server: every dataset and tool, exact verified commands |
| [`FINAL_ALGORITHMIC_SCAN.md`](FINAL_ALGORITHMIC_SCAN.md) | the cross-claim audit that found the one bug no per-claim test could see |

## Per-claim verdicts — read the status column, several are superseded

| doc | claim | verdict |
|---|---|---|
| [`CLAIM1_FINAL_VERDICT.md`](CLAIM1_FINAL_VERDICT.md) | COMPACT | **SUPERSEDED** — 14 datasets and an unrun 15th. Now 19/19, and SRR10676752 is in the sweep |
| [`CLAIM2_FINAL_VERDICT.md`](CLAIM2_FINAL_VERDICT.md) | FAITHFUL | **SUPERSEDED** — chr20 windows. Now full chr20, 4 individuals, called from the archive |
| [`CLAIM3_LOCKED.md`](CLAIM3_LOCKED.md) | ADDRESSABLE | **SUPERSEDED** — a different experiment (outer ARCS repo, MEGAHIT/SPAdes) |
| [`CLAIM3_MECHANISM.md`](CLAIM3_MECHANISM.md) | ADDRESSABLE | **CURRENT** — the mechanism: a het locus is N parallel places, so no coordinate can name it |
| [`NOVELTY_FINAL.md`](NOVELTY_FINAL.md) | all three | **CURRENT** — what is novel, with the prior art actually read |
| [`CLAIMS_FINAL.md`](CLAIMS_FINAL.md) | all three | **CURRENT** — the three claims as they stand |
| [`HG005_EXPLAINED.md`](HG005_EXPLAINED.md) | FAITHFUL | **CURRENT** — why the one published loss happens (variable read length) |
| [`CLAIM2_TABLES_AND_INDEL_SCAN.md`](CLAIM2_TABLES_AND_INDEL_SCAN.md) | FAITHFUL | the T3/T4/T5/T5.2/T5.3 table structure and the full het-indel root-cause history |
| [`HET_INDEL_FRESH_SCAN.md`](HET_INDEL_FRESH_SCAN.md) | FAITHFUL | the 5 findings (2 confirmations, 1 refuted generalization, 2 fixes) that flipped het-indel and tetraploid-indel to wins |

## Checklists (industrial + research grade)

| doc |
|---|
| [`INDUSTRIAL_CHECKLIST_OVERALL.md`](INDUSTRIAL_CHECKLIST_OVERALL.md) / [`RESEARCH_CHECKLIST_OVERALL.md`](RESEARCH_CHECKLIST_OVERALL.md) — whole-product synthesis |
| [`INDUSTRIAL_CHECKLIST_CLAIM1.md`](INDUSTRIAL_CHECKLIST_CLAIM1.md) / [`RESEARCH_CHECKLIST_CLAIM1.md`](RESEARCH_CHECKLIST_CLAIM1.md) |
| [`INDUSTRIAL_CHECKLIST_CLAIM2.md`](INDUSTRIAL_CHECKLIST_CLAIM2.md) / [`RESEARCH_CHECKLIST_CLAIM2.md`](RESEARCH_CHECKLIST_CLAIM2.md) |
| [`INDUSTRIAL_CHECKLIST_CLAIM3.md`](INDUSTRIAL_CHECKLIST_CLAIM3.md) / [`RESEARCH_CHECKLIST_CLAIM3.md`](RESEARCH_CHECKLIST_CLAIM3.md) |

## Architecture and format detail

| doc |
|---|
| [`CAPSULE_FORMAT.md`](CAPSULE_FORMAT.md) — the archive container format |
| [`CALLER_ARCHITECTURE_PLAN.md`](CALLER_ARCHITECTURE_PLAN.md) — the Claim 2 caller's layer design |
| [`SOTA_COMPARISON.md`](SOTA_COMPARISON.md) — every competitor table, all three claims |
| [`CLAIM3_PRIOR_ART.md`](CLAIM3_PRIOR_ART.md) — literature survey for addressability novelty |
| [`HET_INDEL_SOTA.md`](HET_INDEL_SOTA.md) — literature survey for reference-free indel calling |
| [`POLYPLOID_BENCHMARK.md`](POLYPLOID_BENCHMARK.md) — multi-allelic (real) vs polyploid (synthetic) distinction |
| [`FAILURES_AND_REFUTED_IDEAS.md`](FAILURES_AND_REFUTED_IDEAS.md) — the project's own refuted-ideas ledger |
| [`PROJECT_AUDIT.md`](PROJECT_AUDIT.md) — a point-in-time audit of every claim vs spec |
| [`PAPER_DRAFT_CLAIM1.md`](PAPER_DRAFT_CLAIM1.md) — methods/architecture draft for a paper |

## Competitor internals (read once, referenced often)

| doc |
|---|
| [`DISCOSNP_INTERNALS.md`](DISCOSNP_INTERNALS.md) — DiscoSNP++ source, layer by layer |
| [`HOW_DISCOSNP_WINS.md`](HOW_DISCOSNP_WINS.md) — where its advantage actually shows up, measured |
| [`KMER2SNP_BENCHMARK.md`](KMER2SNP_BENCHMARK.md) — the compatibility fixes needed to run it, and its result |
| [`HOW_PGRC2_CODES_REFERENCES.md`](HOW_PGRC2_CODES_REFERENCES.md), [`PGRC2_DATA_ADAPTIVE_STUDY.md`](PGRC2_DATA_ADAPTIVE_STUDY.md), [`PGRC2_DISK_ARCHITECTURE.md`](PGRC2_DISK_ARCHITECTURE.md), [`PGRC2_STREAM_COMPARISON.md`](PGRC2_STREAM_COMPARISON.md) — PgRC2 source analysis |
| [`GENOZIP_RATIO_EXPLAINED.md`](GENOZIP_RATIO_EXPLAINED.md) — diagnosing an apparent Genozip anomaly (it wasn't one) |

## Historical: plans, headroom analyses, and process (kept per project policy, not current status)

These record how conclusions were reached — several are superseded by the
verdict docs above, and say so where relevant. Read only if you need the
history of a specific decision.

<details>
<summary>Plans</summary>

`CLAIM2_BUILDER.md`, `CLAIM2_FINAL_PLAN.md`, `NEXT_LEVEL_PLAN.md`,
`PHASE_A_PLAN.md`, `PLAN_DERIVE_NOT_SEARCH.md`, `PLAN_FLIP_THE_LOSS.md`,
`PLAN_NEXT.md`, `PLAN_PG_VOLUME.md`, `SCOPE_AND_PLAN.md`, `WHAT_IS_NEXT.md`
</details>

<details>
<summary>Headroom / performance analyses</summary>

`FINAL_HEADROOM.md`, `GENUINE_HEADROOM.md`, `HEADROOM_ANALYSIS.md`,
`SPEED_RAM_HEADROOM.md`, `LAYER_BY_LAYER_ANALYSIS.md`,
`COST_AWARE_ACCEPTANCE.md`, `COVERAGE_AWARE_MAXMAP.md`,
`SECOND_REGION_SELF_MATCH.md`, `DO_WE_NEED_THEIR_3WAY.md`,
`78_suffix_array_negative.md`
</details>

<details>
<summary>Superseded results snapshots</summary>

`CLAIM2_RESULTS.md`, `CLAIM2_RESULTS_V2.md` (het-indel figures superseded,
see banner; het-SNV figures still valid), `PHASE2B_RESULT.md`,
`COMBINED_PIPELINE_RESULT.md`, `INDEL_LOSS_SKELETAL.md`,
`INDEL_PRECISION_ROOT_CAUSE.md`
</details>

<details>
<summary>Reimplementation notes and scope history</summary>

`BEST_METHOD_C_REIMPLEMENTATION.md`, `REIMPL_NOTES.md`,
`LOCKED_SEQORDER_SCOPE.md`, `VARIABLE_LENGTH_DESIGN.md`, `NAMES_PHASE2.md`,
`PIPELINE_LEDGER.md`, `CLAIM2_DATA_AND_TOOLS.md`
</details>
