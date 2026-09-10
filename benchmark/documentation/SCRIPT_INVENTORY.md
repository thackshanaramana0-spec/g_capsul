# SCRIPT_INVENTORY — every script in `scripts/`, and whether it is live

Generated from the tree, not written by hand. A reviewer opening `scripts/`
finds 55 files; 17 of them produced the published results and the rest are the
record of decisions that were later reversed. This says which is which, so no
one has to guess.


## Live — symlinked into `benchmark/scripts/` and used by the published sweep

| script | purpose |
|---|---|
| `benchmark_0_preflight.sh` | ═══════════════════════════════════════════════════════════════════════════ |
| `benchmark_1_run.sh` | ═══════════════════════════════════════════════════════════════════════════ |
| `build106.sh` | Build the in-process encoder (stage 106). |
| `build_decode.sh` | Build the archive decoder. Mirrors build106.sh, including the gcc/g++ split |
| `build_tetraploid_truth.py` | — |
| `capsule_config.sh` | Single source of truth for dataset/tool locations across every CAPSULE |
| `encode_adaptive.sh` | Adaptive mapping-ceiling selection. |
| `kmer2snp_sam_to_vcf.py` | Convert Kmer2SNP SNP-kmer-pair output to a genome-coordinate VCF so it can be |
| `lift_vcf.py` | Lift an ARCS contig-coordinate call VCF to genome coordinates (allele-aware) so a |
| `run_fullchr20_archive_capsule.sh` | FULL-SCALE Claim 2 benchmark: whole chr20 at 30x, not a 400 kb window. |
| `run_fullchr20_bench_disco.sh` | FULL-SCALE Claim 2 competitor arm: DiscoSNP++ on the whole chr20 at 30x. |
| `run_kmer2snp.sh` | FULL-SCALE Claim 2 competitor arm: Kmer2SNP on the whole chr20 at 30x. |
| `run_locus_fidelity.sh` | T3.4 — LOCUS RETRIEVAL FIDELITY: does a locus query return the evidence? |
| `run_multiallelic_bench_capsule.sh` | T5.2 — real multi-allelic benchmark, one command. |
| `run_tests.sh` | ============================================================================ |
| `run_tetraploid_bench_capsule.sh` | T5.3 — real tetraploid benchmark, one command. |
| `sanity_archive_one.sh` | ═══════════════════════════════════════════════════════════════════════════ |

## Superseded / historical — kept deliberately, not deleted

These are not broken; they are earlier approaches. Several were **reverted after
measurement**, and this repo's standing rule is that reverted work stays visible
so it is not re-attempted. Do not cite a number produced by one of these.

| script | purpose |
|---|---|
| `76_combined_pipeline.sh` | STAGE 76 -- the combined end-to-end pipeline driver: sequence + order + |
| `82_combined_pipeline_fast.sh` | STAGE 82 -- stage 76's driver, updated to use stage 81's encode-only |
| `86_combined_pipeline_realbound.sh` | STAGE 86 -- stage 82's driver, updated to use stage 84/85's real |
| `88_combined_pipeline_safe.sh` | STAGE 88 -- stage 86's driver, updated to require stage 87's atomically- |
| `89_combined_sequential.sh` | STAGE 89 -- combined pipeline, back to SEQUENTIAL execution (revert of |
| `91_conditional_reorder.sh` | STAGE 91 -- resolves the reorder-vs-no-reorder loss found on SARS-CoV-2 |
| `bench.sh` | Per-stage benchmark harness for the reimpl prototypes. |
| `benchmark_final.sh` | ============================================================================ |
| `benchmark_sanity_one.sh` | ═══════════════════════════════════════════════════════════════════════════ |
| `build_tetraploid_truth_v2.py` | — |
| `capsule_roundtrip.sh` | compress to ONE file, then decompress using ONLY that file. Single command |
| `decode_105.py` | Real, full-pipeline decoder for the locked sequence+order scope. |
| `decode_locked_seqorder.py` | Real, full-pipeline decoder for the locked sequence+order scope. |
| `encode_and_call_once.sh` | Compress with the adaptive grid, then call variants ONCE -- not once per |
| `eval_caller.py` | ---- ground truth ---- |
| `extract_vcfeval_metrics.py` | — |
| `fetch_claim3_refs.sh` | ═══════════════════════════════════════════════════════════════════════════ |
| `rerun_indel_current.sh` | Run the CURRENT caller on all 8 cached windows and score under the corrected |
| `rerun_indel_fix1only.sh` | FIX-1-ONLY: correct multi-allelic CALL classification, but keep the ORIGINAL |
| `rescore_indel_conventions.sh` | Re-score all 8 het-indel evaluations under BOTH the old and the corrected |
| `run_capsule.sh` | CAPSULE — one entry point, three claims. Written for a reviewer who wants |
| `run_claim1_bench.sh` | T1 (archive size) + T2 (compress/decompress time + peak RAM), whole-file |
| `run_claim3.sh` | Claim 3 — ADDRESSABLE: T6 archive-derived analysis vs conventional pipeline. |
| `run_full_pipeline.sh` | One-command, reproducible run of the full reimpl pipeline PLUS the two |
| `run_fullchr20_bench_capsule.sh` | FULL-SCALE Claim 2 benchmark: whole chr20 at 30x, not a 400 kb window. |
| `run_giab_indel_capsule.sh` | CAPSULE version of the outer ARCS project's run_giab_indel.sh: REAL-GIAB |
| `run_indel_bench_capsule.sh` | CAPSULE version of the outer ARCS project's run_indel_bench.sh: synthetic |
| `run_locked_seqorder.sh` | Single-command runner for the locked sequence+order scope |
| `run_pipeline_bench.sh` | ONE PIPELINE, ALL 19 DATASETS. |
| `run_polyploid_bench_capsule.sh` | CAPSULE version of the outer ARCS project's run_polyploid_bench.sh. |
| `run_window_bench_capsule.sh` | CAPSULE reference-free het-SNV/indel benchmark on ONE chr20 window of ONE |
| `sim_indel.py` | Simulate a diploid genome with known het SNVs + het indels, emit ~cov x reads. |
| `sim_indel_bench.py` | Large synthetic diploid benchmark for reference-free indel calling, scored by the SAME |
| `sim_polyploid.py` | Synthetic POLYPLOID genome for validating multi-allelic reference-free SNV calling. |
| `test_claim2.sh` | Claim 2 synthetic regression test — deterministic, no GIAB download needed. |
| `test_claim3.sh` | Claim 3 synthetic regression test — deterministic, no external data needed. |
| `validate_caller_opt.sh` | Gate for caller RAM/speed changes: prove the optimised binary produces the |
| `verify_lossless.sh` | Encode a FASTQ, then decode the streams back and compare against the ORIGINAL |

## The ones worth knowing about specifically

- `verify_lossless.sh` — **decodes the encoder's DUMPED intermediate streams, not
  the archive.** That gap hid four silent data-loss bugs (`CLAUDE.md` §6.1/6.3).
  The published lossless verdicts come from `benchmark_1_run.sh`, which decodes
  the real archive. Keep this distinction; it is the single most expensive
  lesson in the project.
- `run_fullchr20_bench_capsule.sh` — the FASTQ calling path. Superseded by
  `run_fullchr20_archive_capsule.sh`: every published Claim 2 number is called
  **from the archive**, with no FASTQ and no reference.
- `decode_105.py` — carries the same reverse-strand containment bug fixed in the
  C++ decoder (`121fea9`). Reference only.
- `sim_indel.py`, `sim_polyploid.py` — simulated data. The published T2.5 uses
  **real** HG003+HG004 reads, not simulation.
- `gpt2026_*` (not listed above) — a concurrent agent's experiment. Untracked,
  on no branch, and no part of this result.

