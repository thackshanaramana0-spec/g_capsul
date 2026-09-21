---
Date: 2026-09-19
Title: Claim 3 — Code-to-Table Mapping and Result File Index
Purpose: Exact function/line references in stages/capsule_decode.cpp for
  every Claim 3 operation, and the exact file path for every number cited
  anywhere in this folder. Traced directly from the actual repo this
  session, not recalled from memory.
When to refer to this file: Verifying a claim in another file of this
  folder against the real code; onboarding to the codebase; before citing
  any number, to find its exact source file; when a diagram or paper
  section needs a precise line reference.
Keywords: code mapping, capsule_decode.cpp, line numbers, function trace,
  contig_spans, find_occurrences, completion index, benchmark results, file index
---

# Code-to-table mapping

All line numbers are from `stages/capsule_decode.cpp` as committed at
`e2be4e5` ("Claim 3 multi-individual generalization + T3.1 export-bug fix,
flags documented"), 2026-09-19. If the file has since changed, re-verify
with `grep -n` before trusting a specific number below — the surrounding
prose describing the mechanism should still hold even if a line shifts.

## T3.1 — Export

- **Per-contig export via `contig_spans`**: `capsule_decode.cpp:903-957`.
  Reads the `contig_spans` stream (LEB128-encoded gap/length pairs relative
  to the already-in-memory `pg` buffer), slices `pg.data()+offset` directly
  into contigs — no FASTA round-trip, no re-decode of any stream. Gated by
  `getenv("CAPSULE_EXPORT_CONTIGS")`.
- **Two-record legacy form** (what was being used before the fix, still
  used deliberately by the internal ploidy-gate path): same code block,
  the `else if(out_contigs && !getenv("CAPSULE_EXPORT_CONTIGS"))` branch —
  emits `pg[0:MAINEND]` and `pg[MAINEND:end]` as two records.
- **`contig_spans` written at encode time**: `stages/106_inprocess.cpp`,
  `g_contig_spans` populated at lines 2154, 2162, 3115 — gated by
  `CAPS_SPANS` (line 337: `CAPS_SPANS = getenv("CAPS_SPANS") || CAPS_CALL`).
  Comment at lines 320-336 explains the CAPS_SPANS-vs-CAPS_CALL cost
  asymmetry (~20x heavier with CAPS_CALL, since it also runs the full
  inline variant caller).

## T3.2 — Coverage

- **Fast-exit before pg rebuild**: `capsule_decode.cpp:417-490`
  (`if(mode=="coverage")` block). Reads `pos_abs`/`pos_sec`/`pos_region`
  (via `caps_join_positions`), `read_lengths`, `orig2uid_flags`/`orig2uid_vals`.
- **Sweep-line mechanism**: lines ~458-465, `diff[a]++, diff[b]--` over a
  `PGLEN`-sized array, one pass to emit runs.
- **Duplicate-undercount bug fix**: documented inline at lines 424-432 of
  the same block (`orig2uid` expansion required, else 20% undercount
  measured on E. coli).

## T3.3 / T3.4 — Query (shared function)

- **Sidecar load** (`.qidx` file): `capsule_decode.cpp:518-586` approx —
  2-bit packed pg, placements `P`/`L`, optional deviations (`mmix`/`mmoffs`/
  `mmsyms`), strand bitmap (`qstr`), N positions (`nix`/`noffs`).
- **T3.3's `rr` population** (literal coordinates): lines 640-642,
  `strtoull(modearg, "-")`.
- **T3.4's `rr` population** (pg-search): `find_occurrences()`, defined at
  `capsule_decode.cpp` ~line 227 onward (pigeonhole seed-and-extend,
  `MINSEED` derived from haystack size at the lambda beginning ~line 265).
- **the completion index completion**: lines 644-720 approx. Reads `<qidx>completion index`, builds
  candidate set from stride-aligned k-mers of the probe (both strands),
  verifies each candidate by reconstructing the full read (deviation +
  strand + N applied) and doing a literal substring test — the actual
  completeness-guarantee code, not just documentation of one.
- **the completion index build (encode-adjacent, at `index` time)**: gated by `CAPS_XMI`
  and requiring `CAPS_PILEUP=1` first (see `capsule_decode.cpp` around the
  `index` command's pileup/variant-site section, same function that writes
  `.sites`).

## T3.5 — Locus Retrieval

- **Not in `capsule_decode.cpp` at all.** Built at the orchestration layer:
  `scripts/t34_realdata/run_window.sh` (the `query()` shell function near
  the bottom of the script — calls `x_decode query` twice per site, once
  per probe direction) and `scripts/score_bilateral.py` (reads the two
  `.occ`/`.fa` output sets, unions them, assigns each read to its
  best-overlapping occurrence, tallies REF/ALT).
- The parameterized copy used for HG003/4/5/new-locus this session:
  `~/run_window_multi.sh` (WSL-local, built from `run_window.sh` with
  BAM/VCF URLs pulled out as arguments — not committed as a separate repo
  file; the committed `run_window.sh` remains the HG002-hardcoded original
  and is the one to generalize from if reproducing this).

## The `.sites` pileup (referenced in T3.1's refuted consensus-polish
## attempt and T3.5's mechanism discussion)

- Tally built during `index`, `capsule_decode.cpp` ~line 1474-1512
  (`tally[{a0+j, obs[off+m]}]++`). **Known strand bug, deliberately left
  unfixed** — documented as a comment at the tally site (see the commit
  message and `CLAUDE.md`'s 2026-09-19 section for why it was not adopted).

---

# Result file index — exact source for every number cited in this folder

| Number / claim | File |
|---|---|
| T3.1/T3.2/T3.3 speed table (all rows in files 01-03) | `benchmark/results/claim3_T3.1_T3.2_T3.3.csv` |
| T3.1 final correctness numbers (98.698% genome fraction etc.) | `docs/T3.1_CORRECTNESS_FINAL_20260919.md`, "SUPERSEDES" section at the bottom |
| T3.1's withdrawn 80.8%/90%-unaligned numbers, and why | same file, body sections above the SUPERSEDES marker |
| Why ERR5181310/SARS-CoV-2 is excluded from T3.1's fair-comparison set | `docs/T3.1_ERR5181310_ASSEMBLY_CHECK.md` |
| T3.4/T3.5 multi-individual results (all tables in files 04-05) | `docs/T34_T35_MULTI_INDIVIDUAL_20260919.md` |
| T3.4/T3.5's original single-individual (HG002-only) real-data trail | `docs/CLAIM3_LOCUS_ADDRESSABILITY.md`, `docs/T34_RESIDUE_DIAGNOSIS.md` |
| Raw CSVs behind the HG002 real-data trail | `results/T34_REALDATA_20260915/` (BILATERAL.csv, BOTH_FULL.csv, QUERY_ALONE_FULL.csv, GENERALIZATION_FIXED.csv, etc.) |
| PgRC2 comparison (confirmed vs inferred) | `docs/PGRC2_EXTENSION_ANALYSIS.md` |
| Flags reference (`CAPS_SPANS`, `CAPSULE_EXPORT_CONTIGS`, `CAPS_QUERY_CONTAIN`, etc.) | `PIPELINE.md`, flags table |
| Session-level findings summary, dated | `CLAUDE.md`, "2026-09-19 session" section |
| Architecture diagram source (Mermaid) | `docs/CLAIM3_ARCHITECTURE_DIAGRAM.md` |
| Future work: BEETL speed/size comparison, not yet run | `docs/T34_SPEED_BENCHMARK_NEXT_STEPS.md` |

## What was actually benchmarked this session vs what is documentation-only

**Actually run, with real output captured this session** (not just cited
from prior sessions): T3.1 at multiple configurations (archive-optimal,
tuned, dedup, consensus-polish, strand-fixed consensus-polish) via QUAST
against `NC_012967.1`; T3.4/T3.5 on HG003/HG004/HG005 and a second HG005
locus, via `run_window_multi.sh` + `score_bilateral.py` +
`exact_match_recall.py`; a completion index size/index-existence check across all
four individuals; a `.sites` strand-fix build-and-gate cycle (reverted).

**Not run this session, cited from earlier work**: the original T3.1/T3.2/T3.3
CSV (`benchmark/results/claim3_T3.1_T3.2_T3.3.csv`) and the original
single-individual T3.4/T3.5 real-data trail — both pre-date this session and
were read, not regenerated.

**Not run at all, anywhere, explicitly future work**: any benchmark against
BEETL, CIndex, or sFASTQ (speed or index size).
