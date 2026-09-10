# `_removable/` — proposals that were never executed

**Nothing here is deleted. Nothing here is citable.**

Every document in this folder is a PLAN, a DESIGN, or a PROJECTION whose
central content was never built, or was built and then refuted by the document
that superseded it. They contain numbers that were *predicted*, not measured,
and a reader who quotes one will be quoting an intention.

They are kept because this repository's standing rule is that reverted and
refuted work stays visible so it is not silently re-attempted — several of
these were re-proposed more than once before being measured. The folder name
says what it is: material that can be removed without losing a result.

**Where the real content is:**

| you want | go to |
|---|---|
| the results | `../../benchmark/results/` |
| where a number came from | `../../benchmark/documentation/RESULT_CODE.md` |
| the paper | `../../paper/` |
| what was refuted, *with its measurement* | `../FAILURES_AND_REFUTED_IDEAS.md` — that file stays in `docs/`, because every entry in it was implemented and measured |

## What is in here, and why each one is here

| document | why it is not citable |
|---|---|
| `CALLER_ARCHITECTURE_PLAN.md` | its own header: *"This is a design document, not results. Do not cite anything here as a finding."* |
| `GRAPH_HARVEST_EXACT_PLAN.md` | **superseded and declared wrong** by `../GRAPH_HARVEST_REFUTED.md`, which built it and measured 85 loci where thousands were needed |
| `CHAINING_BRANCH_TAG_PLAN.md` | same — the refutation names this file explicitly |
| `QUERY_WINDOWED_PLAN.md` | *"measured and planned, not implemented."* Its closure measurements are real and are repeated in `../LOCALITY_TENSION.md`, which stays |
| `INDEL_PASS_SPEED_PLAN.md` | a speed plan for a path the final benchmark does not use |
| `NEXT_LEVEL_PLAN.md` | plan written 2026-08-29, overtaken by everything after it |
| `PLAN_NEXT.md` | its own header: *"Nothing below has been started."* |
| `PLAN_PG_VOLUME.md`, `PLAN_FLIP_THE_LOSS.md` | pg-volume plans; the loss they target was later flipped by an unrelated fix (`3e06957`) |
| `PLAN_DERIVE_NOT_SEARCH.md` | proposes replacing the parameter sweep with a derivation. Never built; the sweep shipped |
| `PHASE_A_PLAN.md` | plan for work that was done and is recorded in its own result files |
| `SUPERKMER_PLAN.md` | disk-spill plan; the shipped fix was a declared budget instead (`../REPRODUCIBILITY.md`) |
| `SCOPE_AND_PLAN.md` | scope definition from an earlier dataset set (10+4), superseded by the locked 19 |
| `CLAIM2_FINAL_PLAN.md` | engineering plan; the caller that shipped differs from it |
| `CLAIM2_BUILDER.md` | build status from 2026-09-02, superseded |
| `HANDOFF_PROMPT.md` | a prompt for handing the project to another model, not a record of it |

## Two files deliberately NOT moved here

- **`../LOCALITY_TENSION.md`** — its proposed remedy is refuted, but its
  *measurements* (reference reach-back distribution, transitive closure at
  1 kb / 10 kb / 100 kb / 1 MB) are real and are the evidence retiring an
  earlier "windowed rebuild is impossible" claim. Measurement outranks
  proposal.
- **`../GPT2026_PLAN.md`, `../GPT2026_FINDINGS.md`** — a concurrent agent's
  investigation, not ours to file.
