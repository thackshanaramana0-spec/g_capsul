# Pre-flight checklist — before the full-scale run

Completed 2026-09-04. Every line is a check that was RUN, not reasoned about.

## 1. Verification results

| # | check | result |
|---|---|---|
| 1 | **Claim 1 bit-identical** with `CAPS_CALL` unset | **PASS** — every archive stream and assembly stat identical on ERR5181310 |
| 2 | every tunable audited for fitted constants | **1 DEFECT FOUND AND FIXED** (see §2) |
| 3 | thread-safety of all 11 OpenMP regions | **PASS** — every shared write is a `critical`, a `reduction`, or thread-local |
| 4 | ploidy gate rejects haploid input | **PASS** — M. tuberculosis frac=0.0077 -> haploid, **0 records emitted** |
| 5 | full pipeline (substrate NOT skipped) still works | **PASS** — SNV F1 0.921 |
| 6 | disk-spill path produces identical `kc` | **PASS** — 1,063,607 nodes / 5,131 branching / 687 bubbles, spill on and off |

## 2. The defect the checklist caught

`CAPS_DBG_MAXPOLY` defaulted to **3** (DiscoSNP++'s own `-P` default), but every
number reported for this channel was measured with `MAXPOLY=1` passed
explicitly. **A bare run would have used a configuration that scores worse than
anything reported.** Measured both ways:

| P | without coherence | with coherence |
|---|---|---|
| **1** | **0.828** | **0.861** |
| 2 | 0.794 | 0.850 |
| 3 | 0.765 | 0.850 |

Recall saturates by P=2 while false positives keep accumulating. Default is now
1, and `CAPS_DBG_FLANK` was likewise moved to 1000 to match every validated run.

**Verified after the fix:** a bare `CAPS_DBG=1 CAPS_DBG_ONLY=1` run reproduces
TP=317 FP=17 FN=85, F1 0.861 — identical to the fully-flagged run.

## 3. Two validated operating points

| | Method B (graph only) | Full pipeline | DiscoSNP++ |
|---|---|---|---|
| SNV F1 (r2 tuning) | 0.861 | **0.921** | 0.847 (full chr20) |
| SNV F1 (held-out mean) | **0.911** | not re-run | — |
| precision | 0.949 | 0.923 | 0.951 |
| recall | 0.789 | **0.920** | 0.763 |
| caller time | **0.424 s** | ~6 s | — |
| indel F1 | n/a (SNV only) | 0.636 | 0.576 (full chr20) |

Both beat DiscoSNP++ on SNV F1. Method B is the cheap point; the full pipeline
is the accurate one. They are different products, not competing versions.

## 4. Parameter provenance — all Method B tunables

| parameter | default | derived or fitted? |
|---|---|---|
| `MINC` | `min(KVALLEY, H/3)`, floor 2 | **derived** from this dataset's k-mer histogram valley, capped below het-allele depth |
| `MAXPOLY` | 1 | swept, flat-region choice, verified with and without coherence |
| `COHC` | 2 | swept; C=1 and C=2 identical, C=3 costs recall |
| `BAL` | 0.0 (off) | measured harmful at every level once coherence exists |
| `MAXEXT` | 60 | bubble-length bound; measured no effect at 120 |
| `TIPLEN` | `3 * 31` | expressed in k, not an absolute |
| `FLANK` | 1000 | mapping-only; measured F1-neutral |
| `SPILL_BITS` | 8 | partition count; correctness-invariant |

No fitted per-dataset constants remain.

## 5. What is MEASURED vs PROJECTED going into the run

**Measured:** everything in §1 and §3, on HG002 r2 plus held-out r3/na/r4/r5.

**Projected, and the run exists to settle these:**

| | projection | basis |
|---|---|---|
| Method B full-chr20 time | ~57 s | window components scaled; substrate now genuinely 0 |
| Method B full-chr20 RAM | ~3.6 GB | `seqs` 1.86 GB + `kc` 1.44 GB + partition buffer |
| DiscoSNP++ | 76.5 s / 3.45 GB | **measured, not projected** |

**The RAM projection has never been validated at any scale** — at 75K reads `kc`
is 12 MB, so the spill mechanism has nothing to save and its benefit is unproven.

## 6. Known limitations carried into the run

1. **SARS-CoV-2 misclassifies as diploid** (frac 0.854). A deeply-sequenced virus
   is a genuine quasispecies, so the boundary is real, not a bug — but it means
   the ploidy gate is validated for bacteria vs human, not for viruses.
2. **Method B is SNV-only.** Indel emission is gated off: measured -0.051 F1 on
   held-out windows.
3. **Full-scale behaviour has surprised us before**, in both directions — both
   tools lost F1 going window -> full chr20 (ours 0.890->0.849, theirs
   0.874->0.847). Expect attenuation, not a clean transfer.
4. **Nothing in this session is committed.**
