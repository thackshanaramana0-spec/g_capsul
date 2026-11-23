# Code-to-table mapping

> Functions and files in src/decoder.cpp behind each Claim 3 operation, and the source file of every number cited in this folder.

All line numbers are from `src/decoder.cpp` as committed at
`e447ab4` ("Claim 3 multi-individual generalization + T3.1 export-bug fix,
flags documented"). If the file has since changed, re-verify
with `grep -n` before trusting a specific number below — the surrounding
prose describing the mechanism should still hold even if a line shifts.

## T3.1 — Export

- **Per-contig export via `contig_spans`**: `decoder.cpp:903-957`.
  Reads the `contig_spans` stream (LEB128-encoded gap/length pairs relative
  to the already-in-memory `pg` buffer), slices `pg.data()+offset` directly
  into contigs — no FASTA round-trip, no re-decode of any stream. Gated by
  `getenv("CAPSULE_EXPORT_CONTIGS")`.
- **Two-record legacy form** (what was being used before the fix, still
  used deliberately by the internal ploidy-gate path): same code block,
  the `else if(out_contigs && !getenv("CAPSULE_EXPORT_CONTIGS"))` branch —
  emits `pg[0:MAINEND]` and `pg[MAINEND:end]` as two records.
- **`contig_spans` written at encode time**: `src/encoder.cpp`,
  `g_contig_spans` populated at lines 2154, 2162, 3115 — gated by
  `CAPS_SPANS` (line 337: `CAPS_SPANS = getenv("CAPS_SPANS") || CAPS_CALL`).
  Comment at lines 320-336 explains the CAPS_SPANS-vs-CAPS_CALL cost
  asymmetry (~20x heavier with CAPS_CALL, since it also runs the full
  inline variant caller).

## T3.2 — Coverage

- **Fast-exit before pg rebuild**: `decoder.cpp:417-490`
  (`if(mode=="coverage")` block). Reads `pos_abs`/`pos_sec`/`pos_region`
  (via `caps_join_positions`), `read_lengths`, `orig2uid_flags`/`orig2uid_vals`.
- **Sweep-line mechanism**: lines ~458-465, `diff[a]++, diff[b]--` over a
  `PGLEN`-sized array, one pass to emit runs.
- **Duplicate-undercount bug fix**: documented inline at lines 424-432 of
  the same block (`orig2uid` expansion required, else 20% undercount
  measured on E. coli).

## T3.3 / T3.4 — Query (shared function)

- **Decoded sidecar** (`.qidx` file): loaded in `mode=="query"` at `decoder.cpp:558` onward —
  2-bit packed pg, placements `P`/`L`, optional deviations (`mmix`/`mmoffs`/
  `mmsyms`), strand bitmap (`qstr`), N positions (`nix`/`noffs`). `L` holds the
  length per unique read, written by `caps_unique_lengths()` (`decoder.cpp:368`),
  because the stored `read_lengths` stream is per original read.
- **T3.3's `rr` population** (literal coordinates): `strtoull(modearg, "-")`
  in the same function, after the sequence-probe branch.
- **T3.4's `rr` population** (pg-search): `find_occurrences()`, defined at
  `decoder.cpp:229` onward (pigeonhole seed-and-extend,
  `MINSEED` derived from haystack size inside it).
- **FM-index completion** (`CAPS_QUERY_CONTAIN=1`): `decoder.cpp:698` onward. Loads
  `<qidx>.fmi`, runs one backward search per probe and strand
  (`FMIndex::find`, `include/caps_fmi.h`), and adds the reads it returns to
  the pseudogenome hits. Reads stored verbatim in the `.fmi` are emitted as
  themselves after the main emit loop (`decoder.cpp:876`).
  `CAPS_QUERY_INDEX=xmi` selects the older k-mer completion index instead;
  `CAPS_QUERY_INDEX=fmi` makes a missing `.fmi` an error.
- **FM-index build (at `index` time)**: gated by `CAPS_FMI=1`, which implies the
  deviations (`WANT_PILEUP`, `decoder.cpp:1064`). The block at
  `decoder.cpp:1747` indexes every unique read whose emitted record is not a
  plain pseudogenome slice, plus a verbatim copy of every original read the
  emitted record does not reproduce. The index itself (SA-IS, 2-bit BWT, rank
  checkpoints, save/load) is `include/caps_fmi.h`. The older k-mer index
  build (`CAPS_XMI=1`) is still in the same function, `decoder.cpp:1648`.

## T3.5 — Locus Retrieval

- **Not in `decoder.cpp` at all.** Built at the orchestration layer:
  `scripts/claim3_scripts/run_locus_retrieval_capsule.sh` (the `query()` shell function
  — calls `capsule_decode query` twice per site, once per probe direction) and `scripts/claim3_scripts/score_locus_retrieval.py` (reads the two
  `.occ`/`.fa` output sets, unions them, assigns each read to its
  best-overlapping occurrence, tallies REF/ALT).
- The same script serves all four individuals: it takes the individual, the
  chr20 window and a work directory (`run_locus_retrieval_capsule.sh HG003 3000000
  3600000 <workdir>`), downloads the public BAM, VCF and reference slice, encodes,
  checks the archive is lossless, builds the index and scores both tables.
  `benchmark_1_run.sh` runs it for all four individuals as phase 4.

## The `.sites` pileup (referenced in T3.1's refuted consensus-polish
## attempt and T3.5's mechanism discussion)

- Tally built during `index`, `decoder.cpp` ~line 1474-1512
  (`tally[{a0+j, obs[off+m]}]++`). **Known strand bug, deliberately left
  unfixed** — documented as a comment at the tally site (see the commit
  message for why it was not adopted).

---

# Result file index — exact source for every number cited in this folder

| Number / claim | File |
|---|---|
| T3.1/T3.2/T3.3 speed table (all rows in files 01-03) | `results/claim3/claim3_T3.1_T3.2_T3.3.csv` |
| T3.1 final correctness numbers (98.698% genome fraction etc.) | the T3.1 correctness note, "SUPERSEDES" section at the bottom |
| T3.1's withdrawn 80.8%/90%-unaligned numbers, and why | same file, body sections above the SUPERSEDES marker |
| Why ERR5181310/SARS-CoV-2 is excluded from T3.1's fair-comparison set | the ERR5181310 assembly-check note |
| T3.4/T3.5 multi-individual results (all tables in files 04-05) | `docs/claim3/t34_t35_multi_individual_claim3.md` |
| T3.4/T3.5's original single-individual (HG002-only) real-data trail | the Claim 3 locus-addressability write-up, the T3.4 residue diagnosis |
| Raw CSVs behind the HG002 real-data trail | not shipped. The published T3.4 and T3.5 numbers are in `results/claim3/claim3_T3.4_exact_match_recall.csv`, `claim3_T3.5_locus_fidelity_100site.csv` and `claim3_T3.5_locus_retrieval.csv` |
| PgRC2 comparison (confirmed vs inferred) | the PgRC2 extension analysis |
| Flags reference (`CAPS_SPANS`, `CAPSULE_EXPORT_CONTIGS`, `CAPS_QUERY_CONTAIN`, etc.) | the pipeline notes, flags table |
| Findings summary | `docs/claim3/t34_t35_multi_individual_claim3.md` |
| Architecture diagram source (Mermaid) | the Claim 3 architecture-diagram draft |
| Future work: BEETL speed/size comparison, not yet run | the T3.4 speed-benchmark plan |

## What was benchmarked in the final verification vs what is documentation-only

**Run, with output captured** (not just cited
from earlier work): T3.1 at multiple configurations (archive-optimal,
tuned, dedup, consensus-polish, strand-fixed consensus-polish) via QUAST
against `NC_012967.1`; T3.4/T3.5 on HG003/HG004/HG005 and a second HG005
locus, via `scripts/claim3_scripts/run_locus_retrieval_capsule.sh` +
`scripts/claim3_scripts/score_locus_retrieval.py` +
`scripts/claim3_scripts/score_exact_match_recall.py`; the FM-index tests
(`tests/fmi_selftest.cpp`, `tests/fmi_end_to_end.py`) and its size check across all
four individuals; a `.sites` strand-fix build-and-gate cycle (reverted).

**Not run, cited from earlier work**: the original T3.1/T3.2/T3.3
CSV (`results/claim3/claim3_T3.1_T3.2_T3.3.csv`) and the original
single-individual T3.4/T3.5 real-data trail — both pre-date the final
verification and were read, not regenerated.

**Not run at all, anywhere, explicitly future work**: any benchmark against
BEETL, CIndex, or sFASTQ (speed or index size).
