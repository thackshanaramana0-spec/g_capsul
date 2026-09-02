# Industrial checklist — CAPSULE, the whole product

Written 2026-09-03. CAPSULE is marketed as three claims (COMPACT, FAITHFUL,
ADDRESSABLE) but ships as **one archive format and one pair of binaries**
— `stages/106_inprocess.cpp` (encoder, all three claims' write paths) and
`stages/capsule_decode.cpp` (decoder, all three claims' read paths). This
document synthesizes `INDUSTRIAL_CHECKLIST_CLAIM{1,2,3}.md` into a
whole-product view rather than repeating each row three times — read the
per-claim files for row-level detail; read this one for what's true of the
product as a whole, including the one defect that only existed in the
*interaction* between claims and could not have been found by any per-claim
checklist alone (`docs/FINAL_ALGORITHMIC_SCAN.md`).

Status legend: ✅ done · ⚠️ partial / real gap, scoped · ❌ not done · N/A.

| Area | Status | Synthesis across all three claims |
|---|---|---|
| **Repository structure** | ✅ | Consistent across all three: `stages/`, `include/`, `scripts/`, `docs/`, `results/`, `thirdparty/`. No claim invented its own layout. |
| **Architecture** | ⚠️ | Two large files carry the whole product: `106_inprocess.cpp` (encoder, Claims 1+2's write path) and `caps_caller.h` (2132 lines, Claim 2's logic). Neither is split into modules. This is a real, shared architectural debt across claims, not three separate small ones. |
| **Entry points** | ✅ | One consistent pattern for all three: env vars (`CAPS_CALL`, `CAPS_NAMES`, `CAPS_QUAL`) gate what a single encoder invocation does, and `capsule_decode <mode>` gates what the decoder does. A user does not need three different tools. |
| **Configuration** | ✅ (product) / ⚠️ (scripts) | The binaries themselves are fully env/argument-driven. The *benchmark scripts* (Claim 2's GIAB scripts especially) hardcode machine-convention paths (`~/miniconda3`, `~/giab_indel_capsule`) — a real, named, non-blocking gap. |
| **Reproducibility** | ⚠️ | Claim 1: strong (one-command `verify_lossless.sh`, predates this session). Claim 3: strong for E. coli (`scripts/run_claim3.sh`, built this session). Claim 2: the weak point — every accuracy number is a small-window result, full-scale run not yet done. **The product's reproducibility story is only as strong as its weakest claim**, which is Claim 2. |
| **Dependency management** | ⚠️ | No `requirements.txt`/lockfile anywhere in the repo — a genuine, shared, product-wide gap, not three separate ones. `docs/SERVER_SETUP_AND_DOWNLOADS.md` (new, this session) at least makes every tool's install method explicit and verified, which is the practical mitigation until a real lockfile exists. |
| **Determinism** | ✅ | True product-wide: the encoder is fully deterministic, and every new test added this session (`test_claim2.sh` seed `4242`, `test_claim3.sh` seed `1234`) is explicitly seeded. |
| **Testing** | ✅ (new, product-wide) | Before this week: Claim 1 had a real test, Claims 2 and 3 had **zero**. Now all three have one (`verify_lossless.sh`, `test_claim2.sh`, `test_claim3.sh`), plus this session added the one test that matters most for a "one product" claim: the three-claims-combined manual integration check in `docs/FINAL_ALGORITHMIC_SCAN.md`, which is what actually found the one real bug in this pass. **This integration check is not yet a script** — it was run by hand and should become one; see Recommendations. |
| **Edge cases** | ⚠️ | Genuinely strong on the specific edge cases each claim's own history already found (>256bp reads, orphaned unique reads, RC-contained reads — Claim 1; region boundaries, duplicate reads — Claim 3). Weak on generic edges no claim has tested: empty input, single-read input, zero-heterozygosity regions. |
| **Validation** | ⚠️ | No claim validates malformed input before processing — consistent gap across all three, not a surprise found in only one. |
| **Error handling** | ⚠️ | Consistently reasonable (fopen null-checks, no blanket exception swallowing found anywhere) but not audited at the same forensic depth in all three — Claim 3 got the deepest pass, Claim 1 and 2 were not re-audited to that same depth this session. |
| **Logging** | ✅ | Consistent product-wide convention: `[STAGE]`/`[CAPS-CALL]`/`[coverage]`/`[export]`/`[query]` prefixed stderr lines everywhere, no scattered ad-hoc prints in any of the three claims. |
| **Code quality** | ⚠️ | Consistent style (documented constants, refuted-idea comments) across all three, but no claim's code got a full line-by-line quality pass this session — the checklists' Architecture and Code quality rows are honest about not having done this at Claim-3-audit depth everywhere. |
| **Documentation in code** | ✅ | The strongest, most consistent habit in the whole product — every claim's code comments explain *why*, and refuted approaches are kept in place with their measured cost rather than deleted. This is a product-wide culture, not a per-claim accident. |
| **Dead/debug code** | ✅ (with one named exception) | Refuted features behind flags everywhere, by convention. The one exception: Claim 3's superseded `coverage` implementation and its `uid_of` lambda are genuinely dead (unreachable, not flag-gated) — named in `CLAIM3_LOCKED.md` §3 as a small, scoped cleanup item, not fixed there or here. |
| **Security** | ✅ | No credentials/tokens anywhere in any of the three claims' code, scripts, or the new setup doc (which only references public URLs). |
| **File handling** | ⚠️ → mostly ✅ (fixed this session) | The one real file-handling defect found this session — `g_input_path` breaking after a sweep's `chdir()` — was product-wide in effect (it could have silently broken names/quality on ANY claim's archive, not just a specific one) and is now fixed. Smaller, still-open gaps: `run_giab_indel_capsule.sh`'s fixed workdir would collide on concurrent runs. |
| **Resource handling** | ⚠️ | `test_claim2.sh`/`test_claim3.sh` don't clean up their `mktemp -d` workdirs on a passing run — same small gap, found independently in both per-claim checklists, worth one shared fix rather than two. |
| **Performance** | ✅ | Consistently measured, not assumed, across all three: Claim 1's `-fopenmp` story, Claim 3's hoisted-coverage optimization, Claim 2's timing intentionally out of scope by explicit product decision. |
| **Scalability** | ⚠️ | Claim 1: measured across the widest real range (30MB-15GB). Claim 3: measured on a real 1.55M-read dataset. **Claim 2 is the product's scalability gap** — no result above ~75K reads. This is the single most consequential unresolved item for the whole product, not just for Claim 2 in isolation, since a paper presenting "one product" would need all three claims credible at the same scale. |
| **Interfaces** | ✅ | Stable across every claim's history — no claim's CLI/env-var contract has changed since being introduced. |
| **Backward compatibility** | N/A | No external consumers of any claim's interface yet. |
| **Data provenance** | ✅ | Consistently strong: every dataset in every claim is named by exact accession with documented reasoning for any substitution (Claim 1's Utricularia gibba swap; Claim 2's GIAB accessions; Claim 3's reuse of Claim 1's E. coli accession). |
| **Result provenance** | ✅ (fixed this session) | Claim 1's at-risk `/tmp` CSV rescued; Claim 3's query-speedup number corrected after being found to use the wrong baseline. Both fixes are the same underlying discipline applied twice: don't let a number stand unverified just because it was already written down. |
| **CI** | ❌ | Zero `.github/workflows` anywhere. This is the single largest true product-wide gap and the only checklist row identical across all three per-claim files — worth fixing once, for the whole product, rather than three times. |
| **Versioning** | ✅ | Exceptionally strong and consistent: every fix in every claim's history is its own descriptive commit, no `final_v2` anti-pattern anywhere in 100+ commits across the whole project. |
| **Licensing** | ⚠️ | `thirdparty/fse/LICENSE` fixed this session (applies product-wide, since FSE is used by the shared encoder, not one claim). Still open, product-wide, not per-claim: this repository's own code has no top-level `LICENSE` file — a decision for the user, not something any per-claim pass should invent. |
| **Maintainability** | ⚠️ | Genuinely good within a claim (a developer can modify `capsule_decode.cpp`'s export mode without understanding the caller). Genuinely coupled across claims in one specific way, now documented rather than hidden: the `orig2uid`/`pos_abs` indexing invariant is read by both the decoder (Claim 3) and, implicitly, anything that reasons about read identity — a developer changing that invariant needs to know it's shared, and `TECHNICAL_ARCHITECTURE.md` §10.2 now says so explicitly. |

## What only appears when the product is tested as one thing

This is the section a per-claim checklist structurally cannot produce:
**one real, previously-undiscovered bug** was found this session by running
Claim 2's calling, Claim 1's names/quality, and Claim 3's export/coverage
all against the same archive — a combination no prior test, in any claim,
had exercised. Full writeup: `docs/FINAL_ALGORITHMIC_SCAN.md`. Fixed,
verified byte-identical against every existing locked result, and confirmed
working end-to-end afterward.

## Recommendations, in priority order

1. **Turn the manual three-claims-combined check into a script**
   (`scripts/test_integration.sh` or similar) so this class of bug is
   caught automatically on future changes, not only when someone happens to
   run all three features together by hand again.
2. **Wire up CI** — even a minimal workflow running the three per-claim
   test scripts plus the integration check on every push would have caught
   this session's bug automatically instead of requiring a deliberate scan.
3. **Decide and add a top-level project LICENSE.**
4. **Close Claim 2's scalability gap** — this is the one item that blocks
   describing the whole product as validated at the scale its own spec
   commits to.
