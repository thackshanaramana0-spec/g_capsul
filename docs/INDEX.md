# Documentation index

68 documents. Per this project's own standing rule, nothing here is deleted
when superseded — retractions and refuted ideas stay on record, marked in
place. This index exists so a new reader can navigate that history without
reading all 68 files. If you only read five documents, read the five in
**Start here**.

## Start here

| doc | what it is |
|---|---|
| [`../README.md`](../README.md) | product overview, results, quick start |
| [`TECHNICAL_ARCHITECTURE.md`](TECHNICAL_ARCHITECTURE.md) | the complete architecture — assembly, streams, coders, caller, decoder |
| [`COMMANDS_REFERENCE.md`](COMMANDS_REFERENCE.md) | every command, organized by claim |
| [`SERVER_SETUP_AND_DOWNLOADS.md`](SERVER_SETUP_AND_DOWNLOADS.md) | bootstrap a fresh server: every dataset and tool, exact verified commands |
| [`FINAL_ALGORITHMIC_SCAN.md`](FINAL_ALGORITHMIC_SCAN.md) | the cross-claim audit that found the one bug no per-claim test could see |

## Current status, per claim (authoritative)

| doc | claim | verdict |
|---|---|---|
| [`CLAIM1_FINAL_VERDICT.md`](CLAIM1_FINAL_VERDICT.md) | COMPACT | 14/14 wins vs SPRING/Genozip; one dataset (Utricularia gibba) not yet run |
| [`CLAIM2_FINAL_VERDICT.md`](CLAIM2_FINAL_VERDICT.md) | FAITHFUL | 5/5 comparisons win; validated on chr20 windows, not the full-scale individual |
| [`CLAIM3_LOCKED.md`](CLAIM3_LOCKED.md) | ADDRESSABLE | locked — export/coverage/query, all beat their spec targets |
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
