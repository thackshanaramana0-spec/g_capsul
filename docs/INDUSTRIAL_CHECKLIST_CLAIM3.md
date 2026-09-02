# Industrial checklist — Claim 3 (ADDRESSABLE)

Scope: this checklist covers only the code, scripts and docs that implement
and evaluate Claim 3 — `stages/capsule_decode.cpp`'s export/coverage/query
modes, `scripts/run_claim3.sh`, `scripts/test_claim3.sh`, and
`docs/CLAIM3_LOCKED.md`/`CLAIM3_PRIOR_ART.md`. It does not re-audit the
encoder (`106_inprocess.cpp`) or the caller (Claim 2), which have their own
history and are out of scope here.

Status legend: ✅ done · ⚠️ partial / real gap, scoped · ❌ not done, no
excuse · N/A not applicable to this claim.

Two real gaps were fixed while writing this checklist rather than just
reported — both explained inline: `scripts/test_claim3.sh` and
`thirdparty/fse/LICENSE`.

| Area | Status | Finding |
|---|---|---|
| **Repository structure** | ✅ | `stages/` (implementation), `scripts/` (build + run + test), `docs/` (analysis), `results/` (raw output, outer repo) — Claim 3's files sit exactly where the existing structure says they should. |
| **Architecture** | ✅ | The three operations are early-exit modes inside one function, not three copy-pasted implementations. `export`/`coverage` stop before costs they don't need (see `CLAIM3_LOCKED.md` §2). Clear boundary: decode-only, no dependency on the caller or names/quality coders. |
| **Entry points** | ✅ | `capsule_decode export\|coverage\|query <archive> <out> [range]` is a real CLI, documented in the binary's own usage message (`capsule_decode.cpp:662-664`) and in `scripts/run_claim3.sh`'s header. No source editing required to run it. |
| **Configuration** | ✅ | `scripts/run_claim3.sh` takes `DATA_DIR`, `OUT_DIR`, `SPADES`, `REF_ECOLI` as arguments/env vars, all with documented defaults. No machine-specific path is hardcoded in `capsule_decode.cpp` itself (checked directly: zero `/root`, `/tmp`, `/data` literals in that file). |
| **Reproducibility** | ✅ (E. coli) / ⚠️ (GIAB) | `scripts/run_claim3.sh` reproduces the full E. coli export/coverage/query result from one command in a clean environment (builds its own binaries, installs SPAdes if absent). GIAB rows still require manually pointing at pre-built archives — documented as an open item, not hidden. |
| **Dependency management** | ⚠️ | No `requirements.txt`/lockfile exists anywhere in this repo (it's C++/bash, not Python, for this claim's code — `test_claim3.sh`'s one Python block uses only the stdlib `random`, no pinning needed). The real dependency is SPAdes/bwa/mosdepth versions, which **are** now pinned/logged: SPAdes v4.0.0 exact download URL in the script, bwa/samtools/mosdepth version-checked at runtime but not pinned to an exact version. Gap: pin bwa/samtools/mosdepth versions too, or at least log `--version` output into the results dir. |
| **Determinism** | ✅ | `scripts/test_claim3.sh`'s synthetic generator uses `random.seed(1234)` — same input every run. The encoder itself has no randomness in its decision path (verified elsewhere in this repo: `verify_lossless.sh` produces byte-identical archives run to run). |
| **Testing** | ✅ (new) | Was ❌ before this pass — there was no test for export/coverage/query at all. `scripts/test_claim3.sh` now exists: 6 checks (export non-empty/pure-ACGT, coverage-total invariant, query dedup-count invariant, query position-bounds invariant, full lossless roundtrip). **Verified to actually catch the real bug it targets**: built the pre-fix binary from commit `2b5437a~1` and confirmed the coverage-total check fails on it (10,380 vs true 12,780, an 18.8% undercount — same failure class as the real 20% found on E. coli). This is a regression test that would have caught bug 1 automatically had it existed before. |
| **Edge cases** | ⚠️ | Tested: duplicate reads (via the coverage invariant), reads at the main/second-region boundary (verified in `CLAIM3_LOCKED.md` §3 bug 2). Not tested: empty archive, single-read archive, a query range entirely outside `[0,PGLEN)`, malformed archive (truncated file). These are plausible inputs for a "run this on anything" tool and are not covered. |
| **Validation** | ⚠️ | `query`'s `START-END` parsing (`capsule_decode.cpp:391-394`) does not validate that `START<END`, that both are non-negative integers, or that the string actually contains a `-`. A malformed argument like `query 100` (no dash) is caught (`fprintf` + `return 2`), but `query abc-def` would call `strtoull` on non-numeric text and silently produce `qa=0, qb=0` rather than erroring — an unvalidated-input gap. |
| **Error handling** | ✅ mostly | Every `fopen` in the export/coverage/query paths is null-checked with a clear message and non-zero return (`capsule_decode.cpp:241,323,372,396`). No blanket `catch(...)` (this is C, not C++ exceptions, so N/A in that specific form) or swallowed error codes found in the Claim 3 code paths. Gap noted above under Validation. |
| **Logging** | ✅ | All three operations log to stderr with a consistent `[operation] ...` prefix (`capsule_decode.cpp:259,334,387,412`) reporting size/count/timing-relevant facts. `scripts/run_claim3.sh` uses structured `phase()`/`pdone()`/`pskip()` helpers, not scattered ad-hoc `echo`. |
| **Code quality** | ✅ | Functions are focused (each mode is a self-contained block inside one function, not spread across files); comments explain *why* (the indexing-invariant comment at `capsule_decode.cpp:201-205` is the textbook example — it names the exact prior bug it prevents). No naming inconsistency found in the reviewed lines. |
| **Documentation in code** | ✅ | The governing comment block at `capsule_decode.cpp:154-163` explains what each mode reads and why it's cheaper than the conventional equivalent — exactly the "why, not what" standard. |
| **Dead/debug code** | ⚠️ (found, not yet removed) | Two real findings from this session's audit, both still present: an entire superseded `coverage` implementation is unreachable (`capsule_decode.cpp:361-389`, dead because the hoisted fixed version at line 200 always returns first) and its `uid_of` lambda (`capsule_decode.cpp:354-360`) is never called. Documented in `CLAIM3_LOCKED.md` §3 with an explicit recommendation to delete, not done here because removing working (if dead) code is its own gated change per this repo's standing rules — **tracked, not silently left**. |
| **Security** | ✅ | No credentials, API keys, tokens or private endpoints anywhere in the Claim 3 code or scripts (checked: `run_claim3.sh` only references public GitHub release URLs for SPAdes). No sensitive dataset committed — `SRR2584863` is public SRA data, already governed by this repo's `DATASET_LOCKED.md`. |
| **File handling** | ✅ | All paths in `run_claim3.sh` are parameterized (`$DATA_DIR`, `$OUT_DIR`, `$WD`) — no absolute machine-specific path baked in except the *documented default* `/data/fastq`, which is overridable via the first argument. |
| **Resource handling** | ⚠️ | `run_claim3.sh` and `test_claim3.sh` write built binaries and intermediate files into `$WD` but never clean them up — by design for `run_claim3.sh` (caching across re-runs is a deliberate feature, see its `pskip` logic), but `test_claim3.sh`'s default `mktemp -d` workdir is never removed after a passing run, which will accumulate temp directories under repeated CI-style invocation. Minor, easy follow-up: `trap 'rm -rf "$W"' EXIT` when no explicit workdir is given. |
| **Performance** | ✅ | `coverage` is deliberately hoisted above the pseudogenome rebuild specifically to avoid an O(archive) cost it doesn't need (`capsule_decode.cpp:195-199` comment) — a real, measured architectural decision, not an assumption. |
| **Scalability** | ✅ | Measured on a real, non-toy dataset: E. coli, 1,553,259 reads, 576 MB input, not a synthetic handful of reads (`CLAIM3_LOCKED.md` §6). The new synthetic test (§Testing above) deliberately stays small (213 reads) because its job is fast invariant-checking, not scale demonstration — the scale claim rests on the real E. coli numbers, not the test. |
| **Interfaces** | ✅ | The three-mode CLI (`export\|coverage\|query <archive> <out> [range]`) is stable and was not changed by any bug fix this session — only internal logic changed, the contract stayed put. |
| **Backward compatibility** | N/A | Claim 3 introduced these operations this session; there is no prior version of this interface to be compatible with. |
| **Data provenance** | ✅ | `SRR2584863` is identified by exact SRA accession, cross-referenced against `DATASET_LOCKED.md`'s locked list; `HG002_pooled.fq`/`HG005_pooled.fq` are identified as GIAB-sourced in `CLAIM3_LOCKED.md` §6.1. |
| **Result provenance** | ✅ | Every number in `results/claim3/t6_results.csv` traces to a command in `CLAIM3_LOCKED.md` §6/§6.4 or to `scripts/run_claim3.sh` directly — including the corrected query number, whose exact 5-repeat measurement commands are printed in §6.4 rather than just the final ratio. |
| **CI** | ❌ | No `.github/workflows` exists anywhere in this repo. `scripts/test_claim3.sh` exists now and returns a real exit code (0/1) suitable for wiring into CI, but nothing runs it automatically on push/PR. This is a repo-wide gap, not Claim-3-specific, and is out of scope to fix unilaterally here — flagged as the single largest remaining industrial gap. |
| **Versioning** | ✅ | Every fix this session is its own commit with a specific, non-generic message (`2b5437a`, `a7987f8`, `ba37c90`, etc. — see `CLAIM3_LOCKED.md` §8's commit map). No `final_v2` style files anywhere in the Claim 3 code path. |
| **Licensing** | ⚠️ (one real gap fixed) | `thirdparty/fse/*.c` headers reference "the LICENSE file in the root directory of this source tree" which was **missing** from the vendored copy — a real compliance gap, now fixed: `thirdparty/fse/LICENSE` added (BSD text, matching this project's own stated choice of the BSD option over FSE's dual GPLv2 license). `thirdparty/ppmd/` correctly states "public domain" in-source, needs no license file. `thirdparty/htscodecs/LICENSE.md` already exists. **Still open**: this repo (`c_star_pg_advance`'s own code) has no top-level `LICENSE` file at all — this is a project-wide decision, not something to invent unilaterally on Claim 3's behalf; flagged for the user to decide, not fixed here. |
| **Maintainability** | ✅ | A developer touching only `capsule_decode.cpp` does not need to understand `106_inprocess.cpp`'s encoder internals to modify export/coverage/query — the file reads its own archive format via the same `read_capsule`/`capsule_decode_stream` helpers the rest of the decoder already used, no new coupling introduced. |

## Net result of this pass

Two real, concrete fixes landed (not just noted): `scripts/test_claim3.sh`
(a genuine regression test, verified against the actual historical bug) and
`thirdparty/fse/LICENSE` (a genuine compliance gap). Everything else marked
⚠️ or ❌ above is disclosed as exactly that — an open, scoped, named item —
rather than glossed over. The only ❌ (CI) and the one open licensing
decision (top-level project license) are repo-wide, not Claim-3-specific,
and are named explicitly rather than silently left for later.
