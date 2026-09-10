# `results/` — working output. NOT the citable numbers.

> **The citable numbers are in [`../benchmark/results/`](../benchmark/results/).**
> This directory is the working record of 14 benchmark run directories, kept because this
> project does not delete superseded measurements. Several of these runs are
> known-wrong and are retained as the evidence for a retraction.

Full per-directory account, including which runs are void and why:
[`../benchmark/documentation/RESULTS_INVENTORY.md`](../benchmark/documentation/RESULTS_INVENTORY.md).

The one exception: `FULL_SWEEP_20260909/` is the run that produced every
published table. All **8 CSVs** there are `cmp`-verified byte-identical to
`../benchmark/results/` (checked, not asserted); `FULL_SWEEP_20260909/`
additionally carries its own `README.md` and `TABLES.md`, which
`benchmark/results/` does not. The `benchmark/results/` copy — not this one —
is what the paper cites.

**Two directories here are named to look citable and are not:**
`claim2/` (development CSVs with superseded scoring conventions) and
`phase_a/` (decision evidence, and its `allphases_14dataset.csv` is VOID —
pre-fix numbers). Both are explained in RESULTS_INVENTORY.md.

| if you want | read |
|---|---|
| the published numbers | `../benchmark/results/*.csv` |
| what script made each number | `../benchmark/documentation/RESULT_CODE.md` |
| how to re-run from nothing | `../benchmark/documentation/REPRODUCE_EVERYTHING.md` |
| what was verified, and what wasn't | `../AUDIT.md` |
