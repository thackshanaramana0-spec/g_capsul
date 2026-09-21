---
Date: 2026-09-19
Title: T3.1 — Pseudogenome Export vs De Novo Assembly (SPAdes)
Purpose: Complete, verified record of what T3.1 measures, how, the real bug
  found and fixed this session, the final locked numbers, and every attempt
  made (and refuted) to close the remaining correctness gap.
When to refer to this file: Writing/checking T3.1's table or text in the
  manuscript; answering "is our export correct"; deciding what correctness
  language is defensible; before citing any T3.1 genome-fraction/mismatch
  number anywhere.
Keywords: T3.1, export, SPAdes, QUAST, genome fraction, contig_spans,
  CAPSULE_EXPORT_CONTIGS, assembly, MAXMAP, MINOV, indels, mismatch rate,
  N50, duplication ratio
---

# T3.1 — Export

## What it measures

Two independent things, kept separate throughout:
1. **Speed**: time to recover the pseudogenome as a standalone FASTA vs time
   for SPAdes 4.0.0 to assemble the same reads from scratch.
2. **Correctness**: how well that exported FASTA scores as a genome assembly
   against the true reference (QUAST), compared to SPAdes's own assembly.

These are independent axes. A tool can be honestly reported as faster
without implying it is also more correct — standard field practice, not a
project-specific excuse (verified against real precedent this session: no
paper claims speed implies correctness, and correctness is always reported
separately when claimed at all).

## Speed — locked, dominant, unconditional

Source: `benchmark/results/claim3_T3.1_T3.2_T3.3.csv`, table=T3.1 rows.

| Dataset | Ours (s) | SPAdes (s) | Speedup |
|---|---|---|---|
| ERR5181310 | 0.40 | 51.73 | 129.3x |
| SRR29296997 | 0.13 | 83.95 | 645.8x |
| SRR2584863 | 0.37 | 205.04 | 554.2x |
| SRR37283774 | 1.49 | 380.88 | 255.6x |
| DRR976266 | 0.70 | 548.82 | 784.0x |
| HG002 | 5.67 | 2,678.83 | 472.5x |

No caveats on this half of the claim. Untouched by anything below.

## Correctness — the real story, in order

### The bug (found and fixed 2026-09-19)

The default `export` path concatenates the entire pseudogenome into **2
FASTA records** (the main region + the second/appended region), regardless
of how many actual overlap-chained contigs exist inside it. On E. coli at
the archive-optimal encoder settings, the true contig count is **287,125**
(at the earlier, since-abandoned tuned settings it was 41,603) — every one
of those real contig boundaries was being erased, producing one giant fake
"contig" with tens of thousands of false junctions where unrelated
sequences were glued together.

Every correctness number measured before this fix was scoring that
concatenation artifact, not the archive's actual content. Specifically:
genome fraction 80.8%, unaligned length 47-90% depending on config — these
numbers are **withdrawn** and must not be cited.

### The fix

Not new code. The archive already stores contig boundaries in a stream
called `contig_spans`, and the decoder already had a per-contig export path
gated behind an environment flag that was simply never being used for
correctness testing:

```
CAPS_SPANS=1  <best106 encode>              # writes contig_spans (cheap)
CAPSULE_EXPORT_CONTIGS=1 <x_decode> export <archive> <out.fa>
```

Full mechanism and exact code locations: `code_mapping_claim3.md`.

**Verified free**: encoding with `CAPS_SPANS=1` produces a byte-identical
archive to encoding without it (18,282,397 B both ways, E. coli, tested
directly). `CAPS_CALL=1` also implies spans but additionally runs the full
inline variant caller — ~20x heavier (913s vs ~250s on HG002) — and should
NOT be used just to get `contig_spans` for export.

**Verified safe**: archive re-checked LOSSLESS through the export-fix path
(3,106,518 reads in, 3,106,518 out, byte-identical). The fix only changes
what `export` chooses to write; it never touches encode or read
reconstruction.

### Final, locked correctness numbers (E. coli, shipping/archive-optimal config)

Source: `docs/T3.1_CORRECTNESS_FINAL_20260919.md`, "SUPERSEDES" section (the
final, authoritative numbers in that file — read the whole file for the full
trail, including two withdrawn intermediate configurations).

| Metric | Ours | SPAdes | Winner |
|---|---|---|---|
| Genome fraction | **98.698%** | 98.346% | **ours** |
| Indels /100kbp | **0.25** | 0.51 | **ours** |
| Mismatches /100kbp | 19.40 | 2.75 | SPAdes |
| Duplication ratio | 2.021 | 1.000 | SPAdes |
| N50 | 2,054 | 139,621 | SPAdes |
| Unaligned length | 1,895,813 | 9,582 | SPAdes |
| Misassemblies | 7 | 1 | SPAdes |

**This is the number to cite. Not "2 of 5" as a headline — state each metric
plainly, this table format, in the paper.** SPAdes wins the majority of
metrics; ours wins genome fraction and indel rate. Both true simultaneously.
Archive cost of reaching this configuration: **zero** (`contig_spans` is
free; the export flag is post-archive, changes no bytes of the archive).

### Why the remaining gap is structural, not a tuning miss — four independent attempts, all refuted

This matters for the paper's discussion section: the gap was not left open
by omission. It was searched for from four different angles this session,
documented in full in `docs/T3.1_CORRECTNESS_FINAL_20260919.md`:

1. **Region-split test** (exclude the "second region" appended pool):
   genome fraction and mismatch rate identical with or without it — refuted,
   the gap lives in the main region, not the appended pool.
2. **Full parameter search** (`MAXMAP`, `MINOV`, `MAXMM` — all three exposed
   knobs, properly controlled after finding and fixing a bug in the first
   attempt at this): a real lever (`MAXMAP`) was found and pushed to its
   measured plateau, buying a real, kept improvement (mismatch rate halved
   at one earlier configuration) — but the archive-optimal shipping config
   turned out to already be as good or better, so no retuning survives in
   the final recommendation. `MAXMM` had zero effect across its whole range.
3. **Contig deduplication** (self-align real contigs, drop >=90%-covered
   near-duplicates at >=95% identity): duplication ratio improved
   (1.99→1.32), but genome fraction dropped just below SPAdes — relocates
   the gap, does not close it.
4. **Consensus polishing** from the archive's own pileup (`.sites`):
   refuted, and informatively so — majority-vote polishing made mismatches
   *worse* (11.86→407/100kbp at one config), because placement in this
   archive is chosen to minimize encoded bits, not to reflect homology, so
   reads from a different genomic copy can vote with their own bases at a
   repeat-collapsed position. This is a real, non-obvious structural finding
   worth stating in discussion: **this pileup is not a biological pileup.**

### Why full parity is out of reach without new algorithm work

Closing the remaining metrics (mismatch rate, duplication ratio, N50) would
require explicit repeat resolution and read error correction — the actual
algorithmic content of a dedicated assembler like SPAdes (multi-k de Bruijn
graph, bubble-popping, BayesHammer-style correction). The pseudogenome is
built by **greedy overlap-chaining optimized to minimize compressed size**,
not optimized for assembly correctness. Both are genuinely "assembly" in the
mechanical sense (confirmed against the project's own `paper/METHODS.md`,
which says "the compressor assembles reads into a pseudogenome by greedy
suffix-prefix overlap") — they simply optimize for different objective
functions, and that objective difference is the entire, sufficient
explanation for the win/loss pattern above.

## What NOT to write in the paper

- Do not cite 80.8% genome fraction, 90% unaligned, or any number from
  before the export-bug fix — all withdrawn.
- Do not imply the pseudogenome is being offered as a competing de novo
  assembler. It is not. See `mechanism_insight_claim3.md`.
- Do not claim "future work will close the remaining gap" as if a fix is
  expected — four real attempts found none available without new algorithm
  work. State it as a structural, disclosed limitation, not a pending TODO.
