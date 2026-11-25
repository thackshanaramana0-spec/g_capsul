# Code-to-table mapping (Claim 3)

This file maps the Claim 3 operations (T3.1 to T3.5: export, coverage, range query, exact-match
retrieval and locus retrieval) to the code that performs them and lists the source file of every
number cited in this folder. Line numbers are given for `src/decoder.cpp` and `src/encoder.cpp` as
of the current source and drift as the files change, so search for the identifier. The overall
design is in [`../technical_architecture.md`](../technical_architecture.md), sections 5 and 6.

## T3.1 Export

- **Per-contig export from `contig_spans`** is in the export block of `decoder.cpp` (the
  `has("contig_spans") && getenv("CAPSULE_EXPORT_CONTIGS")` branch, line 1216). It reads the
  `contig_spans` stream (LEB128 gap and length pairs relative to the pseudogenome buffer already
  in memory) and slices `pg.data()+offset` into contigs, with no FASTA round trip and no re-decode
  of another stream.
- **The two-record form** is the `else if(out_contigs && !getenv("CAPSULE_EXPORT_CONTIGS"))` branch
  (line 1068). It emits `pg[0:MAINEND]` and `pg[MAINEND:end]` as two records, and the internal
  ploidy-gate path still uses it.
- **`contig_spans` at encode time.** `src/encoder.cpp` declares `g_contig_spans` (line 440), gated by
  `CAPS_SPANS`, defined as `getenv("CAPS_SPANS") || CAPS_CALL` (line 431). The comment above it
  (from line 420) explains the cost difference: `CAPS_CALL` is about 20 times heavier, because it
  also runs the full variant caller inline.

## T3.2 Coverage

- The `if(mode=="coverage")` block (`decoder.cpp:528`) exits before the pseudogenome is rebuilt. It
  reads `pos_abs`, `pos_sec`, `pos_region` (through `caps_join_positions`), `read_lengths`, and
  `orig2uid_flags` and `orig2uid_vals`.
- Depth is computed by a sweep line, `diff[a]++` and `diff[b]--` over an array of the pseudogenome
  length, followed by one pass that emits runs.
- A duplicate-undercount bug is documented inline in the same block. The `orig2uid` expansion is
  required, and without it the depth was 20% too low on E. coli.

## T3.3 and T3.4 Query (shared function)

- **The decoded sidecar (`.qidx`)** is loaded in the `mode=="query"` block (`decoder.cpp:615`). It
  holds the 2-bit packed pseudogenome, the placements `P` and lengths `L`, the optional deviations
  (`mmix`, `mmoffs`, `mmsyms`), the strand bitmap `qstr` and the N positions (`nix`, `noffs`). `L`
  has one length per unique read and is written by `caps_unique_lengths()` (line 423), because the
  stored `read_lengths` stream is per original read.
- **T3.3** fills `rr` from literal coordinates, with `strtoull(modearg, "-")` in the same block
  after the sequence-probe branch.
- **T3.4** fills `rr` by searching the pseudogenome with `find_occurrences()` (line 276), a
  pigeonhole seed-and-extend search whose seed length `MINSEED` is derived from the haystack size.
- **FM-index completion** (`CAPS_QUERY_CONTAIN=1`, `want_contain` near line 768) loads
  `<qidx>.fmi` and runs one backward search per probe and strand (`FMIndex::find`,
  `include/caps_fmi.h`). The reads it returns are added to the pseudogenome hits. Reads stored
  verbatim in the `.fmi` are emitted as themselves after the main emit loop.
  `CAPS_QUERY_INDEX=xmi` selects the older k-mer completion index, and `CAPS_QUERY_INDEX=fmi`
  makes a missing `.fmi` an error.
- **Building the FM-index at `index` time** is gated by `CAPS_FMI=1`, which implies the deviations
  (`WANT_PILEUP`, line 1158). It indexes every unique read whose emitted record is not a plain
  pseudogenome slice, plus a verbatim copy of every original read that the emitted record does not
  reproduce. The index itself (SA-IS, 2-bit BWT, rank checkpoints, save and load) is
  `include/caps_fmi.h`. The older k-mer index build (`CAPS_XMI=1`) is still in the same function.

## T3.5 Locus retrieval

- It is not in `decoder.cpp`. It is built at the orchestration layer by
  `scripts/claim3_scripts/run_locus_retrieval_capsule.sh` (its `query()` function calls
  `capsule_decode query` twice per site, once per probe direction) and
  `scripts/claim3_scripts/score_locus_retrieval.py` (reads the two `.occ` and `.fa` output sets,
  unions them, assigns each read to its best-overlapping occurrence and tallies REF and ALT).
- The same script serves all four individuals. Its arguments are the individual, the chromosome 20
  window and a work directory (`run_locus_retrieval_capsule.sh HG003 3000000 3600000 <workdir>`). It
  downloads the public BAM, VCF and reference slice, encodes, checks that the archive is lossless,
  builds the index and scores both tables. `benchmark_1_run.sh` runs it for all four individuals as
  phase 4.

## The `.sites` pileup

The pileup is built during `index` (the `tally[{a0+j,obs[off+m]}]++` statement, `decoder.cpp:1720`).
It is referenced in the refuted consensus-polish attempt for T3.1 and in the T3.5 mechanism
discussion. It has a known strand bug that was left unfixed on purpose, as a comment at the tally
site records.

---

# Source file of every number

| Number | Source |
|---|---|
| T3.1/T3.2/T3.3 speed table (all rows) | `results/claim3/claim3_T3.1_T3.2_T3.3.csv` |
| T3.1 final correctness numbers (98.698% genome fraction etc.) | the T3.1 correctness note, "SUPERSEDES" section at the bottom |
| T3.1's withdrawn 80.8%/90%-unaligned numbers, and why | same file, body sections above the SUPERSEDES marker |
| Why ERR5181310/SARS-CoV-2 is excluded from T3.1's fair-comparison set | the ERR5181310 assembly-check note |
| T3.4/T3.5 multi-individual results (all tables) | `docs/claim3/t34_t35_multi_individual_claim3.md` |
| T3.4/T3.5's original single-individual (HG002-only) real-data trail | the Claim 3 locus-addressability write-up, the T3.4 residue diagnosis |
| Raw CSVs behind the HG002 real-data trail | not shipped. The published T3.4 and T3.5 numbers are in `results/claim3/claim3_T3.4_exact_match_recall.csv`, `claim3_T3.5_locus_fidelity_100site.csv` and `claim3_T3.5_locus_retrieval.csv` |
| PgRC2 comparison (confirmed vs inferred) | the PgRC2 extension analysis |
| Flags reference (`CAPS_SPANS`, `CAPSULE_EXPORT_CONTIGS`, `CAPS_QUERY_CONTAIN`, etc.) | the pipeline notes, flags table |
| Findings summary | `docs/claim3/t34_t35_multi_individual_claim3.md` |
| Architecture diagram source (Mermaid) | the Claim 3 architecture-diagram draft |
| Future work: BEETL speed/size comparison, not yet run | the T3.4 speed-benchmark plan |

## What was run and what is cited

**Run, with the output captured:** T3.1 in several configurations (archive-optimal, tuned, dedup,
consensus-polish and strand-fixed consensus-polish) through QUAST against `NC_012967.1`. T3.4 and
T3.5 on HG003, HG004 and HG005 and a second HG005 locus, with
`scripts/claim3_scripts/run_locus_retrieval_capsule.sh`, `score_locus_retrieval.py` and
`score_exact_match_recall.py`. The FM-index tests (`tests/fmi_selftest.cpp`,
`tests/fmi_end_to_end.py`) and its size check on all four individuals. A `.sites` strand-fix
build-and-gate cycle, which was reverted.

**Cited and not regenerated:** the original T3.1 to T3.3 CSV
(`results/claim3/claim3_T3.1_T3.2_T3.3.csv`) and the original single-individual T3.4 and T3.5
real-data trail, both of which pre-date the final verification.

**Not run:** any benchmark against BEETL, CIndex or sFASTQ, for speed or index size. This is future
work.
