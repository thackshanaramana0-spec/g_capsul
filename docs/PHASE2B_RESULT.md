# Phase 2b — sequence + read order + names, measured on 14 datasets

> **SUPERSEDED 2026-09-10 — numbers only.** This document predates the final
> 19-dataset sweep of 2026-09-09, and CLAUDE.md records this file as VOID for 4 of its 14 datasets. Every one of those figures has been
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

> **VOID for 4 of 14 datasets, superseded 2026-09-02.** The numbers below were
> measured before 4 silent data-loss bugs were found and fixed (commits
> `23be207`, `121fea9`, `f3ab0c2`) — the archives for **ERR552797,
> SRR40271341, SRR40402583, SRR32429602** did not actually decode to their
> input, so their G_CAPSUL sizes here are not valid comparisons. The other 10
> datasets were unaffected (all fixed-length, and the bugs are only reachable
> with variable-length or all-N reads). A full re-run with the fixed binaries,
> covering Phase 1 / Phase 2b / Phase 3 separately plus combined, each level
> round-trip verified before its number is recorded, is in progress — see
> `docs/PHASE3_RESULT.md` (or the latest ALLPHASES result) once it lands, and
> treat this file as historical until then.

**2026-09-02, original run.** The first head-to-head of the G_CAPSUL archive
against real SPRING and Genozip on the locked dataset set. Everything below is
a real file on disk, not a stream total or a log sum.

Dataset set: `NEW_DATASET_LOCKED.md` (15 accessions). 14 ran — Utricularia
gibba (SRR10676752, the Plantae entry) is downloaded but its `fasterq-dump`
failed on `disk-limit exceeded` and it is not yet converted.

PgRC2 is excluded from this table: it has no names stage at all, so it cannot
be compared on this axis. Its sequence-only comparison is in
`PHASE1_RESULT` / section below.

---

## What each column contains, so the comparison is honest

| tool | what its number covers | how |
|---|---|---|
| G_CAPSUL | sequence + read order + names + line-3 mode | `CAPS_NAMES=1` through `scripts/encode_adaptive.sh` (the locked config), real archive size on disk |
| SPRING | sequence + read order + names | `--no-quality` only, so IDs are retained. Real archive size on disk. |
| Genozip | sequence + read order + names | no such switch exists, so summed from its own `genocat --STATS`: SEQ + every `Parent=QNAME` context + `length` + `LINE3` |

`LINE3` is counted on both sides because we now store line 3 too (as a
one-byte file mode — see the line-3 section of `REIMPL_NOTES.md`).

---

## Result

| dataset | organism | G_CAPSUL | SPRING | Genozip | vs SPRING | vs Genozip |
|---|---|---|---|---|---|---|
| SRR39257532 | A. fumigatus | 54,463,702 | 94,085,120 | 165,844,879 | **−42.1%** | −67.2% |
| SRR40402583 | C. jejuni | 5,265,395 | 5,888,000 | 27,160,660 | −10.6% | −80.6% |
| SRR2584863 | E. coli | 12,246,611 | 15,749,120 | 56,400,715 | −22.2% | −78.3% |
| SRR40271341 | H. pylori | 4,338,943 | 7,362,560 | 23,229,546 | −41.1% | −81.3% |
| SRR29296997 | H. salinarum | 3,748,306 | 5,386,240 | 14,806,057 | −30.4% | −74.7% |
| SRR32429602 | HCMV | 34,703,101 | 36,925,440 | 86,067,752 | −6.0% | −59.7% |
| SRR36741279 | L. major | 39,252,323 | 43,735,040 | 121,337,119 | −10.2% | −67.7% |
| ERR552797 | M. tuberculosis | 5,618,490 | 8,611,840 | 39,085,612 | −34.8% | **−85.6%** |
| SRR554369 | P. aeruginosa | 8,985,372 | 9,154,560 | 39,043,073 | −1.8% | −77.0% |
| SRR37283774 | P. falciparum | 24,135,957 | 27,207,680 | 47,351,622 | −11.3% | −49.0% |
| ERR12954017 | S. acidocaldarius | 4,136,613 | 5,509,120 | 28,016,377 | −24.9% | −85.2% |
| ERR17740259 | S. aureus | 20,046,973 | 22,947,840 | 91,428,235 | −12.6% | −78.1% |
| DRR976266 | S. cerevisiae | 30,781,796 | 34,836,480 | 174,166,510 | −11.6% | −82.3% |
| ERR5181310 | SARS-CoV-2 | 839,750 | 1,290,240 | 964,712 | −34.9% | −13.0% |
| **TOTAL** | | **248,563,332** | **318,689,280** | **914,902,869** | **−22.00%** | **−72.83%** |

**14/14 wins against both tools.** Tightest margin is P. aeruginosa at −1.8%
vs SPRING; widest is A. fumigatus at −42.1%.

---

## Phase 1 (sequence + order only), same 14 datasets

Included for completeness; PgRC2 can be compared here because this is its scope.

    vs Genozip   14/14 wins   -78.48%
    vs SPRING    12/14 wins   -23.41%   (losses: HCMV +3.0%, C. jejuni +4.4%)
    vs PgRC2      7/8  wins    -1.37%   (on the 8 it could process)

**PgRC2 failed outright on 6 of 14** — real tool failures, not setup errors:
`Unsupported variable length reads` (SARS-CoV-2, HCMV, C. jejuni) and hard
crashes (`*** stack smashing detected ***` on H. pylori, `free(): invalid next
size` on A. fumigatus). Worth one sentence in the paper: we process inputs it
cannot.

Note that adding names IMPROVES our standing against SPRING (12/14 → 14/14),
because the names column is where our margin over it is largest.

---

## Two invocation traps that invalidated a first run

Recorded because both produce plausible-looking wrong numbers rather than
errors.

1. **`capsule_enc reads.fq 3 16` is not the locked configuration.** It silently
   produces a much larger archive — E. coli 9,809,930 vs the correct 8,225,993,
   which turns a 7.1% win over PgRC2 into a 10.8% loss. Always encode through
   `scripts/encode_adaptive.sh`, which supplies `3 16 16 22 16 16 1 24 64 1`
   plus the 4-candidate MAXMAP/MINOV sweep.
2. **`DUMP_PERM=1 DUMP_MM=1` are required, not optional.** They are named like
   debug flags but gate the per-read streams (`pos_abs`, `pos_strand`,
   `read_lengths`, `mm_*`, `orig2uid_*`). Without them the encoder writes them
   empty and the decoder reconstructs zero reads. `encode_adaptive.sh` sets
   them. `VERIFY_DUMP=1` additionally writes the raw `.bin` dumps.

## One reading trap

`capsule_decode` writes **one sequence per line, not 4-line FASTQ**. Running
`awk 'NR%4==2'` over its output samples every 4th line and produces garbage —
that mistake produced a false "the decoder is broken on every dataset" alarm
in this session. The decoder is correct: verified byte-identical against the
original sequence column on E. coli (1,553,259 reads) and SARS-CoV-2 (912,571
reads), with all 11 streams round-tripping byte-identically.

---

## What is still missing

**Quality**, which is roughly 70% of a losslessly compressed FASTQ. It is the
only column not yet in the archive. See the quality section of
`REIMPL_NOTES.md` for the assessment: vendor fqzcomp (BSD), do not reimplement.
