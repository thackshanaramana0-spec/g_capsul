# Research / academic checklist — Claim 2 (FAITHFUL)

Scope: het-SNV (T3), coverage sweep (T4), het-indel (T5), and multi-allelic
(T5.2, per `docs/CLAIM2_TABLES_AND_INDEL_SCAN.md`). Statistical methodology
marked N/A throughout, per project instruction and precedent (no paper in
this field reports significance for this comparison type).

Status legend: 🔴 critical, 🟠 important. ✅ done · ⚠️ partial/scoped gap ·
❌ not done · N/A.

**Headline finding of this pass: Claim 2 is scientifically sound but not
scale-validated.** The mechanism, the comparisons, and the honesty of the
reporting are all in good shape. The evidence base (window-only, not
full-chromosome) is the real, named gap — see Scalability/Reproducibility
below, and do not let anything ✅ elsewhere read as "therefore ready."

| Area | 🔴/🟠 | Status | Finding |
|---|---|---|---|
| **Research question → code mapping** | 🔴 | ✅ | SNV pileup, bubble/indel extraction, and the multi-allelic path are all named, separable code in `caps_caller.h`, not a black box. |
| **Reproducibility** | 🔴 | ⚠️ | `run_giab_indel_capsule.sh` and `run_window_bench_capsule.sh` are real, runnable scripts — but they reproduce results on a **~75K-read window** (~0.6% of one full 30× individual), and the full-scale run has never been executed. Reproducing the *existing* number is fine; reproducing the *claim as stated in the spec* (full 30× chr20) is not yet possible from any script in this repo. |
| **Dataset provenance** | 🔴 | ✅ | GIAB accessions (HG002-HG005), NIST v4.2.1 truth VCF URLs, and confident-BED sources are named exactly, with the exact chr20 sub-range used for the current window (`2000000-2400000`). |
| **Data integrity** | 🔴 | ⚠️ | No checksum recorded for the downloaded GIAB truth VCFs or the pooled FASTQ files — same gap as Claim 3's checklist, not fixed here either. |
| **Train/test separation** | 🔴 | ✅ | `dup_frac` and other thresholds were swept **once, on a single designated tuning window**, then frozen and evaluated on 4 other windows plus 3 unseen individuals (`CLAIM2_RESULTS_V2.md` §"8 independent evaluations") — this is exactly the held-out discipline the checklist asks for, already in place before this pass. |
| **Ground truth** | 🔴 | ✅ | NIST GIAB v4.2.1 benchmark VCFs, restricted to heterozygous sites inside the confident BED — a standard, named truth source, not an ad-hoc one. |
| **Experimental configuration** | 🔴 | ✅ | Every frozen constant (`MAF=0.20`, `HDMAX=2`, etc.) is named with its value and provenance in `caps_caller.h`'s own comments and `CLAUDE.md`. |
| **Deterministic execution** | 🔴 | ✅ (caller) / ⚠️ (competitors) | The caller itself is deterministic (no RNG in its decision path). DiscoSNP++ and Kmer2SNP's own determinism across runs was not independently re-verified this session. |
| **Baseline implementation** | 🔴 | ⚠️ | DiscoSNP++ is run correctly with its documented `-G` flag and its known POS off-by-one bug corrected before scoring (`CLAUDE.md` §5 explicitly documents this fix). Kmer2SNP, however, has only ever been run as **one aggregate number** (P=0.992/R=0.302/F1=0.464), not per-individual like the other two tools in T3 — if T3 is reported per-individual, Kmer2SNP's column would need separate runs that don't exist yet. |
| **Fair comparison** | 🔴 | ✅ | Same reads, same reference, same rtg vcfeval pipeline, same confident-BED restriction, for every tool compared — no asymmetric treatment found. |
| **Metric implementation** | 🔴 | ✅ | Uses `rtg vcfeval --squash-ploidy`, an established third-party scorer, not a custom F1 implementation — removes an entire class of "did we compute F1 correctly" risk. |
| **Statistical methodology** | — | N/A | Intentionally not applied, per project instruction and field precedent. |
| **Experiment independence** | 🔴 | ✅ | No manually-edited intermediate VCF found in the pipeline between calling and scoring — `lift_vcf.py` is a real, inspectable transformation step, not a hand-patch. |
| **Raw → result pipeline** | 🔴 | ✅ | FASTQ → caller (contig-space VCF) → bwa alignment → `lift_vcf.py` (genome-space VCF) → rtg vcfeval → F1 is fully traceable through named scripts, no manual step. |
| **Losslessness/correctness tests** | 🔴 | ✅ (new) | Was the single biggest gap before this pass — the caller had **zero** automated correctness check, only real-data benchmarks. `scripts/test_claim2.sh` closes this: synthetic diploid genome with known truth, run through the real pipeline, checked for recall (10/10 SNV, 5/5 clean indel on first correct run). |
| **Ablation experiments** | 🟠 | ✅ | The project's own docs record multiple ablations already: the multi-substrate design (aggressive vs mild collapse), the positional-clustering indel channel, and four indel-precision filter attempts, each independently toggleable and each measured (`HOW_DISCOSNP_WINS.md` §4). |
| **Sensitivity analysis** | 🟠 | ✅ | T4 (the coverage sweep, 10×/15×/20×/30×) **is** a sensitivity analysis by definition, and `dup_frac` was independently swept and shown to sit on a flat plateau (0.35-0.50), not a fragile peak (`HOW_DISCOSNP_WINS.md` §5). |
| **Failure cases** | 🔴 | ⚠️ | Not verified: what happens on a region with zero true heterozygous sites, or coverage below the caller's implicit minimum. Nothing in the current pipeline is known to *silently* discard a failed run, but this hasn't been actively tested either. |
| **Negative results** | 🟠 | ✅ | The het-indel loss is on the record, not buried — kept as its own table (T5) rather than replaced by the winning T5.2 result, exactly per this session's own decision. |
| **Resource measurement** | 🟠 | N/A | Timing explicitly out of scope for Claim 2 this session — no methodology to critique because none is being reported. |
| **Hardware/environment** | 🟠 | ⚠️ | Same as Claim 3: relies on the outer `CLAUDE.md`'s one-time hardware description rather than repeating it per-claim. |
| **Random seeds** | 🔴 | ✅ (new) | The new synthetic test's seed (`4242`) is explicit and reported. No other randomness exists in the caller's own decision path. |
| **Dependency versions** | 🔴 | ⚠️ | DiscoSNP++ version is named (v2.6.2-12, per `HOW_DISCOSNP_WINS.md`). Kmer2SNP's dependency chain (DSK, biopython, networkx) required real compatibility fixes this project documented (`docs/KMER2SNP_BENCHMARK.md`) but versions aren't pinned in a lockfile. |
| **Automated tests** | 🔴 | ✅ (new) | Directly addresses what was a real 🔴-critical gap — see Losslessness/correctness tests above. |
| **Regression tests** | 🟠 | ✅ (new, same script) | `scripts/test_claim2.sh` would catch a caller regression that broke SNV/indel detection outright — it is, by construction, exactly this kind of test. |
| **Synthetic/toy tests** | 🟠 | ✅ (new) | The 15-variant synthetic diploid genome is exactly this — small, known, fast. |
| **Independent verification** | 🔴 | ✅ | Every real-data F1 number is computed via `rtg vcfeval`, an independent, widely-used third-party tool, not a project-internal scorer — this is the strongest form of independent verification available for this task. |
| **Figures/tables generation** | 🟠 | ❌ | No script generates `results/claim2/t3_snv_f1.csv`-style output automatically from a single command the way `scripts/run_claim3.sh` does for Claim 3 — T3/T5/T5.2 numbers exist in docs as hand-assembled tables, not machine-generated from a benchmark run's output. |
| **Claim traceability** | 🔴 | ✅ | "ARCS expected avg F1 ≈ 0.936" (the spec's own projection) → measured 0.890 → traced to `CLAIM2_RESULTS_V2.md`'s 8-evaluation table → traced to `run_window_bench_capsule.sh` → traced to `caps_caller.h` is a complete, walkable chain, even though the measured number differs from the spec's original projection (disclosed, not hidden — `CLAUDE.md`'s own rule 6 says projections are sanity checks, not authoritative). |
| **README reproduction guide** | 🔴 | ⚠️ | `run_giab_indel_capsule.sh`'s header comment documents usage, but there is no single top-level "how to reproduce Claim 2" guide the way `scripts/run_claim3.sh`'s header serves that role for Claim 3 — a reviewer would need to read three separate scripts to reconstruct the full T3/T4/T5/T5.2 picture. |

## Net result of this pass

The 🔴-critical gap that mattered most — **zero automated tests on 2132
lines of caller code** — is now ✅, verified against a real (if initially
buggy, now fixed) synthetic ground truth. The scientific methodology itself
(held-out windows, third-party scoring, honest negative-result reporting)
was already sound before this pass and remains so. **The one gap this
checklist cannot mark ✅ no matter how the code is audited is Reproducibility
at the claimed scale** — every accuracy number on record is a window
result, and the full 30× chr20 run the project's own spec commits to has
not been executed. That is the actual blocker to calling Claim 2 done, not
anything in the code.
