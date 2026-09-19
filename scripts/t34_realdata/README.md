# T3.4 / Claim 3 — real-data reproduction

Full writeup: [`docs/CLAIM3_LOCUS_ADDRESSABILITY.md`](../../docs/CLAIM3_LOCUS_ADDRESSABILITY.md)
— read its TERMS section first. It defines *locus*, *probe*, *placement*, the
difference between retrieval by exact match and by position, and what the
number 400 counts.
Diagnosis trail: [`docs/T34_RESIDUE_DIAGNOSIS.md`](../../docs/T34_RESIDUE_DIAGNOSIS.md)
Numbers: [`results/T34_REALDATA_20260915/`](../../results/T34_REALDATA_20260915/)

Everything here runs on real data: reference from UCSC, variants and phasing
from the GIAB HG002 benchmark VCF, reads range-fetched from the GIAB 300x
Illumina BAM. No whole-genome download is involved — T3.4 queries a 600 kb
window, and that window is fetched in seconds.

## One command

```bash
bash run_window.sh 3000000 3600000 ~/w1 400     # the published window
bash run_window.sh 4000000 4600000 ~/w2 400     # a held-out window
```

It gates itself and refuses to continue unless

- the reference slice matches the GIAB REF column base-for-base, and
- the archive round-trips byte-identically to its input.

Both gates exist because both failure modes actually happened. An early run
produced a sidecar with **0 placements** from a missing `DUMP_PERM`, which would
have read as the mechanism failing rather than as a missing flag.

## Generalization

```bash
bash heldout_local.sh                    # held-out sites + a different assembly
SEED=13 FRAC=67 bash heldout_local.sh    # ~19x independent subsample
```

`heldout_local.sh` is the important one. It re-subsamples the same real reads
with a different seed and re-encodes, producing a **different pseudogenome**.
Since the failure mode is an artefact of overlap chaining, a different assembly
is the sharpest retest available without a new download.

## Diagnosis scripts, in the order they were used

| script | question it answers |
|---|---|
| `diag_residual.py` | what distinguishes the failing sites |
| `k_sweep.sh` | does more mismatch tolerance help (no) |
| `probelen_sweep.sh` | does a longer probe lift the `k_max = floor(P/12) - 1` cap (no) |
| `check_input_alleles.py` | is the ALT allele even in the input reads (yes) |
| `where_is_alt.py` | is the ALT context in the pseudogenome (does not discriminate) |
| `why_unreachable.py` | tolerance problem or anchoring problem (anchoring) |
| `bilateral_check.py` | is the damage one-sided (yes: upstream 0/6, downstream 5/6) |
| `bilateral_run.sh` | the fix, measured |
| `position_vs_exact_match.py` | how much evidence an exact-match archive cannot reach (86%) |
| `admission_check.py` | is the site-admission rule biasing the denominator (no, it is conservative) |
| `heldout_local.sh` | held-out sites and an independently re-assembled archive |
| `run_window.sh` / `finish_window.sh` | a whole independent genomic window, fetch and analyse |

## Gates and cost

```bash
bash final_gates.sh        # 5 regression gates, all must pass
bash verify_multiprobe.sh  # multi-probe is an efficiency change only
bash sidecar_cost2.sh      # index build cost vs bwa, per-query latency
bash query_breakdown.sh    # separates process startup from search
```

## Flags that matter

| flag | why |
|---|---|
| `CAPS_SPANS=1` | `query`/`export`/`coverage` need contig spans |
| `DUMP_PERM=1` | **without it the encoder emits no positions and the sidecar has 0 placements** |
| `CAPS_PILEUP=1` | sidecar carries per-read deviations rather than the consensus |
| `CAPS_QUERY_MM=k` | mismatch tolerance, default 0 = exact. Capped at `floor(P/12) - 1` |
