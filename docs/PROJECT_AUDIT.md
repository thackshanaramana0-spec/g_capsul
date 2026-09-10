# Project audit — what is done, what is missing

> **SUPERSEDED 2026-09-10 — numbers only.** This document predates the final
> 19-dataset sweep of 2026-09-09, and it opens with "14/14 whole-file wins" and reports het-SNV 0.890, het-indel 0.637 and multi-allelic 11/18. Every one of those figures has been
> replaced:
>
> | | superseded | **final** |
> |---|---|---|
> | Claim 1 | 14 datasets, various margins | **19/19 vs SPRING (−6.03%), 19/19 vs Genozip (−43.26%), 57/57 LOSSLESS** |
> | het-SNV F1 | 0.890 | **0.876** (mean of 4 individuals, full chr20, called from the archive) |
> | het-indel F1 | 0.637 / 0.666 | **0.621** |
> | multi-allelic | 11/18 | **17/26**, a complete chr20 census (corrected 2026-09-10, was 21/26); DiscoSNP++ 0 of 3,989 records |
> | tetraploid SNV | 0.836 | **0.897** |
> | export vs SPAdes | 555–656x | **129–784x** |
> | T3.4 coordinate | 0/400 | **81/400** (the 0 was our query emitting consensus, not reads) |
>
> Citable source: `benchmark/results/*.csv`, traced in
> `benchmark/documentation/RESULT_CODE.md`. **Where this document and a result
> file disagree, the result file wins.**
>
> Kept because its *reasoning* is still the record of how the conclusion was
> reached — this repo does not delete superseded work, it marks it.

Written 2026-09-02 at the end of the session. Every claim checked against its
spec in `/root/arcs-clean/CLAUDE.md`, not against memory.

---

## Claim 1 — COMPACT ✅ DONE (14 of 15 datasets)

**Result: 14/14 whole-file wins vs SPRING and Genozip; +1.88% vs PgRC2 on
sequence.** All round-trip verified LOSSLESS before any size was recorded.

| phase | vs SPRING | vs Genozip |
|---|---|---|
| Phase 1 (sequence+order) | 11/12 | 12/12 |
| Phase 2b (+names+line3) | 12/12 | 12/12 |
| Phase 3 (whole file) | **14/14** | **14/14** |

**Gap: dataset 15 of 15 (Utricularia gibba, SRR10676752) has never been run.**
The `.sra` (7.3 GB) is on disk at `/tmp/newdl/SRR10676752/`; conversion to FASTQ
was started at the end of this session. Everything else in `NEW_DATASET_LOCKED.md`
is complete.

---

## Claim 2 — FAITHFUL ✅ TWO OF THREE CLASSES WON

All measured on real GIAB data with DiscoSNP++ rerun on identical reads and
identical scoring (its POS off-by-one corrected — without that it scores 0.004).

| class | G_CAPSUL | DiscoSNP++ | verdict |
|---|---|---|---|
| **het-SNV** | **0.890** | 0.874 | **WIN**, 5/8 evaluations, generalises to unseen individuals |
| **multi-allelic** | **11/18 sites** | **0/18** | **WIN** — capability the competitor lacks, replicated on two independent regions |
| het-indel | 0.637 | 0.663 | behind by 0.026, winning 4/8 |

Session movement: het-SNV **0.419 → 0.890**; het-indel **~0.36 → 0.637**.

**Gaps:**
1. **Kmer2SNP was never benchmarked** — it needs the R package `findGSE` and
   `Rscript` is not installed. It is the weakest of the three competitors
   (~0.53 in its own paper and in ARCS's runs), so this does not change any
   conclusion, but the comparison is incomplete.
2. **Synthetic polyploid numbers are not usable** until regenerated — the
   generator had a 28% truth bug (fixed this session, commit `d92a5b1`). This
   also invalidates the outer project's published triploid figures.
3. The outer project's `run_claim2.sh` has still never been run end to end.

---

## Claim 3 — ADDRESSABLE ✅ DONE THIS SESSION, DOMINANT

Was **unbuilt in both projects** — no `export`/`coverage`/`query` existed
anywhere, and no assembler was installed to compare against. All three are now
implemented as early-exit modes in the G_CAPSUL decoder, and MEGAHIT was
installed as the de-novo baseline.

| dataset | export | MEGAHIT | speedup | coverage | bwa+mosdepth | speedup | query |
|---|---|---|---|---|---|---|---|
| HG002 r2 | 0.032 s | 8.14 s | **254×** | 0.052 s | 1.33 s | **26×** | 0.039 s |
| HG002 r5 | 0.028 s | 8.78 s | **314×** | 0.049 s | 1.23 s | **25×** | 0.035 s |
| HG005 r3 | 0.045 s | 12.31 s | **274×** | 0.064 s | 1.46 s | **23×** | 0.051 s |
| E. coli | 0.465 s | SPAdes 258.15 s | **555×** | — | — | — | — |

SPAdes (v4.0.0, default full pipeline: BayesHammer correction + K21/33/55/77
iterative assembly + repeat resolution) is the spec-exact export baseline,
installed and run 2026-09-03 — earlier figures used MEGAHIT as a substitute
because SPAdes was not yet installed. 239 real contigs produced (largest
243,716 bp), peak RAM 5.20 GB. **Gap #5 (export vs the spec-named tool) is now
closed.**

Both beat the spec's targets (≥40× export, 2–5× coverage) by a wide margin, and
the coverage figure is **conservative** — the conventional side used a
pre-built BWA index, while G_CAPSUL needs no reference at all. `query` has no
competitor: it requires per-read coordinates that only an addressable archive
retains.

**Gap: these three operations exist in G_CAPSUL only.** The outer
`/root/arcs-clean/build/arcs` binary still has no such subcommands, so if the
paper describes them as ARCS features they must be ported.

---

## Cross-cutting gaps (all pre-existing, none introduced this session)

1. **Utricularia gibba** — 15th dataset, conversion in progress.
2. **Kmer2SNP** — never run (needs R + findGSE).
3. **Genozip's fungi anomaly** — 162 MB where SPRING gets 24 MB, flagged in
   `CLAUDE.md` §6.2 and still undiagnosed. Must be explained or explicitly
   caveated before publication; an unexplained 6–7× outlier in a competitor
   looks like a methodology error to a reviewer.
4. **Sandbox → outer project reconciliation** — names, quality, line-3, the four
   data-loss fixes, the caller and the Claim 3 operations all live in
   `c_star_pg_advance`. The shipped `arcs` binary has none of them.
5. **Speed/RAM vs SPRING and Genozip** — never measured. Only size was recorded
   in the 14-dataset run, and only PgRC2 has timing comparisons.

## What was corrected during this session (kept on the record)

Four results were retracted or corrected after measurement contradicted them:
the polyploid "win" (scoring artifact), the "assembler collapses haplotypes"
diagnosis (disproved by the link-scan measurement), the "read-derived candidates
are more precise" claim (refuted by channel isolation), and DiscoSNP++'s
apparent 0.004 indel score (its POS off-by-one). Two benchmark-harness bugs and
one generator bug were fixed. Every refuted mechanism is kept in-code behind a
flag rather than deleted.
