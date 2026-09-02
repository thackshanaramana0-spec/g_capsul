# Research / academic checklist — CAPSULE, the whole product

Written 2026-09-03. Synthesizes `RESEARCH_CHECKLIST_CLAIM{1,2,3}.md` into
one product-level view. The paper frames CAPSULE as three claims about one
archive format; this document checks whether the *product-level* claim —
"one pseudogenome, three properties" — holds up, not just whether each
property holds up in isolation. Statistical methodology remains N/A
throughout, per established project instruction and precedent.

Status legend: 🔴 critical, 🟠 important. ✅ done · ⚠️ partial/scoped gap ·
❌ not done · N/A.

| Area | 🔴/🟠 | Synthesis across all three claims |
|---|---|---|
| **Research question → code mapping** | 🔴 | ✅ product-wide. All three claims map to named code: `106_inprocess.cpp` (Claim 1's assembly+coding, Claim 2's `CallData` capture), `caps_caller.h` (Claim 2), `capsule_decode.cpp`'s three modes (Claim 3). No claim is a black box relative to the others. |
| **Reproducibility** | 🔴 | ⚠️ — **this is where the product-level claim is weakest.** Claim 1 and Claim 3 are each independently reproducible from one command. Claim 2 is not, at the scale its own spec requires. A reviewer evaluating "one product, three claims" would correctly note that only two of three are fully reproducible today. |
| **Dataset provenance** | 🔴 | ✅ product-wide, and notably self-consistent: Claim 3 deliberately *reuses* Claim 1's E. coli accession rather than introducing a new one, and Claim 2's GIAB accessions are the same across all its own sub-tables (T3/T5/T5.2). The product doesn't juggle inconsistent datasets across claims. |
| **Data integrity** | 🔴 | ⚠️ — identical gap in all three per-claim checklists (no checksums recorded anywhere). Genuinely one product-wide gap, not three coincidentally identical ones — worth one fix (a checksum manifest) rather than three. |
| **Ground truth** | 🔴 | ✅ (Claim 2) / N/A (Claims 1, 3) — only Claim 2 has a ground-truth-scored metric; the other two are measured directly (archive size, wall time), so this row is only meaningfully evaluable for Claim 2, where it's sound (GIAB v4.2.1, third-party `rtg vcfeval`). |
| **Experimental configuration** | 🔴 | ✅ product-wide. Every frozen constant in every claim is named with its value and, where relevant, its provenance (ported from ARCS, measured as an interior optimum, etc.) — a consistent discipline, not an accident of one claim. |
| **Deterministic execution** | 🔴 | ✅ (CAPSULE itself, all claims) / ⚠️ (competitor tools) — the archive and every operation on it are deterministic product-wide; SPAdes/DiscoSNP++/MEGAHIT's own determinism was not independently re-verified in any claim. |
| **Baseline implementation** | 🔴 | ✅ product-wide with one shared, named gap: Claim 1's PgRC2 version isn't pinned, Claim 2's Kmer2SNP is only run as one aggregate (not per-individual), Claim 3's SPAdes is now the spec-exact baseline (fixed this session, was MEGAHIT). Three different specific gaps, same underlying category: baseline precision varies claim to claim. |
| **Fair comparison** | 🔴 | ✅ product-wide — same reads, same machine, same measurement tool, in every claim, with every known asymmetry disclosed rather than hidden (Claim 3's pre-built BWA index favoring the baseline; Claim 1's idle-machine qualifier). |
| **Metric implementation** | 🔴 | ✅ product-wide, and this session directly demonstrated the discipline working: Claim 3's query-speedup metric was found to be measuring the wrong thing (a stream-dump path, not real reconstruction) and was corrected with a real re-measurement rather than left standing. |
| **Statistical methodology** | — | N/A throughout, consistent with project precedent. |
| **Experiment independence** | 🔴 | ✅ product-wide — no hand-edited intermediate result found in any claim's pipeline. |
| **Raw → result pipeline** | 🔴 | ✅ product-wide, though Claim 2's pipeline is the least automated end-to-end (no single script produces its full T3/T4/T5/T5.2 table the way `run_claim3.sh` does for Claim 3). |
| **Losslessness/correctness tests** | 🔴 | ✅ (new, product-wide as of this session) — this is the row most transformed by this week's work. Before: Claim 1 had a real test, Claims 2/3 had none. Now all three do, **and this session's own cross-claim integration check caught a real bug none of the three individual tests could have caught**, which is itself evidence the testing discipline is working, not just present. |
| **Ablation experiments** | 🟠 | ✅ (Claims 1, 2) / N/A (Claim 3, no tunable components to ablate) — genuinely present and measured, not merely claimed, in both claims where it applies. |
| **Sensitivity analysis** | 🟠 | ✅ product-wide: Claim 1's MAXMAP/MINOV sweep, Claim 2's `dup_frac` plateau and the T4 coverage sweep itself, Claim 3's per-dataset speedup range (254-656× across different inputs, not one cherry-picked number). |
| **Failure cases** | 🔴 | ⚠️ product-wide — no claim tests empty/degenerate input systematically; this is a consistent, not claim-specific, gap. |
| **Negative results** | 🟠 | ✅ product-wide, and this is arguably the product's strongest single credibility signal: Claim 1 discloses its one PgRC2 loss at exact margin, Claim 2 keeps the het-indel loss as its own table rather than replacing it with the multi-allelic win, Claim 3 corrected its own overstated speedup down from 3.3× to 1.62× rather than leaving the bigger number standing. Three independent instances of the same honesty discipline. |
| **Resource measurement** | 🟠 | ✅ (Claims 1, 3) / N/A by decision (Claim 2) — consistent methodology (`/usr/bin/time -v`, parsed the same way) everywhere it's used. |
| **Hardware/environment** | 🟠 | ⚠️ product-wide — all three claims rely on the outer `CLAUDE.md`'s one-time hardware description rather than restating it; a real but minor, consistent gap. |
| **Random seeds** | 🔴 | ✅ product-wide — every source of randomness in the product (both new synthetic tests) is explicitly seeded and reported; the core encoder/caller have none to seed. |
| **Dependency versions** | 🔴 | ⚠️ product-wide — `docs/SERVER_SETUP_AND_DOWNLOADS.md` (new) at least makes every version explicit and verified against this server, which is real progress, but nothing is pinned in an enforceable lockfile. |
| **Automated tests** | 🔴 | ✅ (new, product-wide) — same transformation as Losslessness/correctness tests above. |
| **Regression tests** | 🟠 | ✅ product-wide: Claim 1's `cmp`-identical gating rule (a workflow-level discipline, not just a script), and both new synthetic tests are regression tests by construction, each verified against the real historical bug it targets. |
| **Synthetic/toy tests** | 🟠 | ✅ (Claims 2, 3, new) / ⚠️ (Claim 1, none exists — always tests against real, if small, FASTQ files) — two of three now have this, one gap named. |
| **Independent verification** | 🔴 | ✅ product-wide: Claim 1 cross-validates against published SPRING bpb values, Claim 2 uses third-party `rtg vcfeval`, Claim 3's coverage/query invariants are checked against values derivable from the input independent of the implementation being tested. |
| **Figures/tables generation** | 🟠 | ✅ (Claim 3 only, `run_claim3.sh` → `t6_results.csv`) / ❌ (Claims 1, 2 — tables are hand-assembled from ad-hoc runs). This is the least mature area product-wide and the most mechanical to fix. |
| **Claim traceability** | 🔴 | ✅ product-wide — every headline number in every claim traces to a specific command, script, or doc section, including where the traced number turned out to need correction (Claim 3's query speedup) rather than the trace stopping at the first plausible-looking figure. |
| **README reproduction guide** | 🔴 | ✅ (Claims 1, 3) / ⚠️ (Claim 2, spread across three scripts with no single entry point) — `docs/COMMANDS_REFERENCE.md` (new, this session) now gives every claim a single place listing its actual commands, partially closing this for Claim 2 specifically. |

## The one finding this document can make that no per-claim checklist could

**A whole-product claim needs the product tested as a whole, not three
times separately.** This session's own experience proved the point
directly: three individually-passing test suites (`verify_lossless.sh`,
`test_claim2.sh`, `test_claim3.sh`) coexisted with a real, silent,
product-level data-loss bug (names/quality dropped under specific
conditions) that none of them could see, because none of them combined all
three claims' features on one archive. The bug is fixed
(`docs/FINAL_ALGORITHMIC_SCAN.md`), but the more durable lesson for the
research claim itself is methodological: **"one product, three claims" is
a testable proposition, and until this session it had never actually been
tested.**

## Net result

The product-level research claim is now on firmer ground than it was
before this session, for a specific, checkable reason: it survived the one
test that could have falsified it (features tested together, not just
separately), after that test found and fixed a real defect rather than
passing trivially. The **one 🔴-critical row that remains genuinely
unresolved at the product level is Reproducibility**, and specifically
Claim 2's scale gap — not a code defect, a data-collection one, named
plainly here as it was in Claim 2's own verdict document.
