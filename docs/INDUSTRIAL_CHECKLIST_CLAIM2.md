# Industrial checklist — Claim 2 (FAITHFUL)

Scope: `include/caps_caller.h` (the reference-free variant caller, 2132
lines), the Claim 2 benchmark scripts (`run_giab_indel_capsule.sh`,
`run_window_bench_capsule.sh`, `run_polyploid_bench_capsule.sh`,
`lift_vcf.py`, `eval_caller.py`, `extract_vcfeval_metrics.py`), and the new
`scripts/test_claim2.sh`. Does not re-audit the encoder or Claim 3's
decoder, which have their own checklists.

Status legend: ✅ done · ⚠️ partial / real gap, scoped · ❌ not done · N/A.

**Unlike the Claim 3 pass, this one cannot end in "ready."** Claim 2's
biggest gap is not a code defect — it's that every accuracy number on
record was measured on small chr20 windows (~75K reads, ~0.6% of one
individual's full 30× file), not the full dataset the project's own spec
calls for. That is named here, not hidden, and it is the single largest
item blocking a "locked" verdict, independent of anything below.

One concrete fix landed in this pass, not just a finding:
`scripts/test_claim2.sh` — a real, verified regression test for a caller
that had zero automated tests before this.

| Area | Status | Finding |
|---|---|---|
| **Repository structure** | ✅ | Caller lives in `include/`, benchmark drivers in `scripts/`, results and root-cause analysis in `docs/` — consistent with the rest of the repo. |
| **Architecture** | ⚠️ | `caps_caller.h` is 2132 lines in one header. It has internal structure (substrate building, SNV pileup, bubble extraction, positional clustering are separable named functions), but a 2000+-line single file is large by the checklist's own "no giant monolithic files" bar. Splitting it was not attempted in this pass — it's a real refactor, not a bug fix, and risks introducing regressions in code that has no pre-existing test coverage to catch them (until this session). |
| **Entry points** | ✅ | `CAPS_CALL=1` env var + the same encoder binary is the entry point — documented in the governing comment at `caps_caller.h`'s top and used identically by every benchmark script. No separate binary to build or maintain. |
| **Configuration** | ⚠️ | Thresholds (`HDMAX=2, MC=3, HALF=15, MAF=0.20, DHI=2.5, KHI=1.1, TRI=0.12`) are compile-time constants in the header, not runtime-configurable — deliberate, per the project's own standing rule ("frozen constants... ported verbatim from ARCS, never re-tuned"), so this is a documented choice, not an oversight. The real configuration gap is in the *scripts*: `run_giab_indel_capsule.sh` hardcodes `WD=~/giab_indel_capsule` and `export PATH=~/miniconda3/bin:$PATH` — machine-convention paths, not parameterized. |
| **Reproducibility** | ⚠️ | The GIAB benchmark scripts are real, runnable commands (not manual steps) — but every existing result is on a small window, and the window-selection itself (`CHROM=20; RLO=2000000; RHI=2400000`) is hardcoded per-script rather than a parameter, so reproducing a *different* window requires editing the script, not passing an argument. |
| **Dependency management** | ⚠️ | Same gap as Claim 3's checklist: no `requirements.txt`/lockfile. `run_giab_indel_capsule.sh` depends on `bwa`, `tabix`, `curl`, `rtg` (vcfeval) and a conda environment for Kmer2SNP, none version-pinned in the scripts themselves. |
| **Determinism** | ✅ (new) | `scripts/test_claim2.sh`'s synthetic generator uses `random.seed(4242)` — verified deterministic by running it twice and getting identical 10/10 SNV, 5/5 indel results both times. The caller itself has no randomness in its decision path (same property already established for the encoder). |
| **Testing** | ✅ (new) | Was **zero** before this pass — 2132 lines of caller logic with no automated test, only real-data benchmarks requiring GIAB downloads. `scripts/test_claim2.sh` now exists: builds a synthetic diploid genome with 10 known het SNVs and 5 known clean (non-homopolymer) het indels, runs the real pipeline (caller → bwa alignment → `lift_vcf.py`, the exact same chain `run_giab_indel_capsule.sh` uses), and checks recall. Found and fixed a real bug in its own generator during first use (index drift from mutating a shrinking list left-to-right — see script comments) before landing at 10/10 SNV, 5/5 indel recall, confirmed deterministic. |
| **Edge cases** | ❌ | No test exercises a haploid region, a region with zero heterozygous sites, a region below minimum coverage, or reads shorter than `HALF`/`MC`'s implicit minimums. The new synthetic test is a happy-path sanity check, not an edge-case suite. |
| **Validation** | ⚠️ | The caller assumes well-formed FASTQ input; there's no explicit check for e.g. reads shorter than the anchor length before the algorithm runs — it would presumably just find no candidates rather than crash, but this isn't verified either way. |
| **Error handling** | ⚠️ | Not audited line-by-line in this pass (2132 lines is a large surface) — no obvious blanket exception-swallowing found in a targeted read, but a full pass matching the rigor of Claim 3's audit was not done here, and should be named as not-yet-done rather than assumed clean. |
| **Logging** | ✅ | `[CAPS-CALL]` prefixed stderr lines report substrate/candidate/call counts at each stage (visible in every log excerpt used in this session's own debugging) — structured, not scattered prints. |
| **Code quality** | ⚠️ | Consistent with the rest of the project's style (frozen constants named and commented with their provenance). Not reviewed function-by-function for naming/focus in this pass — the 2132-line size makes a full quality pass a separate, larger effort than this checklist round covers. |
| **Documentation in code** | ✅ | Every non-obvious threshold and every refuted experimental feature is commented with *why*, consistent with the project's own standing rules (`CLAUDE.md` §7's refuted-idea table mirrors flags left in this file). |
| **Dead/debug code** | ✅ | Refuted experimental features are deliberately kept behind env-var flags (`CAPS_GAPSCAN`, `CAPS_NO_HPBUBBLE`, `CAPS_MULTIBUBBLE`, etc.) rather than deleted — a project convention, documented as intentional, not abandoned debris. |
| **Security** | ✅ | No credentials/tokens in caller code or scripts. `run_giab_indel_capsule.sh` fetches public NCBI/GIAB URLs only. |
| **File handling** | ⚠️ | `run_giab_indel_capsule.sh`'s `WD=~/giab_indel_capsule` is a fixed path, not parameterized or a tempdir — a second concurrent run on the same machine would collide. `scripts/test_claim2.sh` correctly uses a caller-supplied or `mktemp -d` workdir. |
| **Resource handling** | ⚠️ | Same as Claim 3's finding: `test_claim2.sh` doesn't clean up its `mktemp -d` workdir on a passing run. Minor, same fix (`trap ... EXIT`) would apply to both. |
| **Performance** | N/A | Timing is explicitly out of scope for Claim 2 per project decision this session ("time we will not report") — no timing table exists or is planned. |
| **Scalability** | ❌ | This is the headline gap: **every existing accuracy number is on a ~75K-read window, not the ~12.6M-read full individual.** A background scaling investigation this session found that small subsamples of the full file exhibit pathologically low read-overlap density (5.9% at 200K reads, climbing to 35.2% at 1M) purely as an artifact of insufficient local coverage at small N — not a defect, but it does mean **no existing result has been validated at the scale the paper's own spec commits to.** |
| **Interfaces** | ✅ | `CAPS_CALL=1` + `CALL_VCF`/`CAPS_DUMP_CONTIGS` env vars is a stable, unchanged interface across every benchmark script in this repo. |
| **Backward compatibility** | N/A | No prior version of this caller interface exists to be compatible with. |
| **Data provenance** | ✅ | GIAB accessions, truth VCF URLs (NIST v4.2.1), and confident-BED sources are named explicitly in `run_giab_indel_capsule.sh`'s own header comments. |
| **Maintainability** | ⚠️ | A developer touching the caller needs to understand the `orig2uid`/`pos_abs` indexing invariant (same one that caused Claim 3's coverage bug) since the caller reads the same streams — this coupling is documented but real, and not something a newcomer would discover without reading `106_inprocess.cpp`'s comments first. |

## Net result of this pass

One real fix landed: `scripts/test_claim2.sh`, a genuine regression test for
previously-untested caller code, verified against its own (found-and-fixed)
generator bug before being trusted. Everything else marked ⚠️ or ❌ is a
named, specific gap — the single most consequential one being **Scalability**:
no result in this claim has been measured at the scale the project commits
to. That is a data-collection gap, not a code-quality one, and it is not
something this checklist pass can close by itself.
