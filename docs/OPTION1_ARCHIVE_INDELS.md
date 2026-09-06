# Option 1 — the FULL caller (SNV + indels) from a stored archive

## The gap this closes

Option 0 serves SNVs from a `.capsule` with no FASTQ (F1 0.8766, verified). But
indels were unreachable: `capsule_decode call` hard-codes `CAPS_DBG_ONLY`, which
skips `indel_pass` entirely. So "everything from the archive" held for T3/T4 and
Claim 3, and was **false for T5/T5.2/T5.3** — all three use the full
`CAPS_CALL` path.

It was not a one-flag change. The full caller needs each read's
(contig id, offset); with empty placement arrays `build_substrate` places NO
reads and the indel pass emits nothing, **silently**.

## Three pieces

1. **`contig_spans` stream** (encoder, `CAPS_CALL`-gated so plain archives stay
   byte-identical). Chain boundaries are the one thing the decoder cannot
   recompute — it rebuilds `pg` faithfully but never sees where one chain ended
   and the next began. Delta + LEB128.
2. **Per-contig export.** `export` emits the pseudogenome as 2 concatenated
   records — right for Claim 3, wrong for the caller, which collapses and
   re-places reads PER CONTIG.
3. **Placement reconstruction:** `pos_abs` + `pos_strand` + `contig_spans` →
   `read_cid`/`read_pos`/`read_rc` by binary search into the span table.

## MEASURED — 4M-read subset, same reads both ways

| | A: FASTQ | B: archive | |
|---|---|---|---|
| input | 1.26 GB FASTQ | **190 MB .capsule** | 6.6x smaller |
| wall | 929.8 s | **595.4 s** | **-36%** |
| peak RAM | 15.07 GB | 15.13 GB | same |
| VCF records | 36,638 | 36,648 | +10 |
| **indels** | **3,081** | **3,094** | +13 |
| placements | in-memory | **4,000,000 / 4,000,000** | perfect rebuild |

**The full path runs from the archive and emits indels.** It is 36% faster
because it skips compression, at the same RAM.

Mechanism verified on real data, not inferred:
* `[SPANS] 334,816 contig spans -> 1,004,451 B` = **3.00 B/span**, matching the
  pre-run arithmetic to the digit (~0.25% of a 547 MB archive).
* `[export] 334,816 contigs` — not the 2-record pseudogenome.
* `4,000,000/4,000,000 read placements rebuilt` — a wrong `orig2uid` delta or
  `pos_abs` width would have shown as a partial count or reads in gaps.

## NOT YET EXACT — 78 of 36,638 records differ (0.21%)

Do **not** describe this as identical. Option 0's SNV path is byte-identical;
this one is not yet.

### The `read_clip` theory is REFUTED — diagnosed, not guessed

`read_clip` does not exist in the encoder at all (zero matches in
`106_inprocess.cpp`). It is computed INSIDE `build_substrate`
(`S.read_clip[o] = best_clip`, caps_caller.h:827) during re-placement — an
OUTPUT, not an input. Both paths therefore generate it identically, and it
cannot be the cause.

### The real cause: placement CARRY-FORWARD, not correctness

    A (FASTQ):   contigs 334816 -> 140561   placed=1,749,455  re-placed=2,011,297
    B (archive): contigs 334816 -> 140561   placed=1,689,820  re-placed=2,070,932

* Both collapse **334,816 -> the identical 140,561 contigs**.
* `H=55` identical. **SNVs identical: 24,769 vs 24,769.** Candidates differ by 1.
* The 78 differing records are **entirely in the indel channel**.

The only difference is WHICH reads carry their placement forward versus get
re-placed: A has the encoder's live per-read state, while B reconstructs from
`pos_abs`, which is stored **per UNIQUE read** — duplicates share one entry.
The re-placement step recovers them (totals match at ~3.76M either way), but a
read recovered by mismatch-tolerant search can land a base or two differently
from one carried forward.

**So this is a fidelity limit of the format, not a bug.** Closing it would mean
storing placements per ORIGINAL read rather than per unique — which is exactly
the redundancy the compressor exists to remove. The honest statement is:

> the archive reproduces the call set to within 0.21% (SNVs identical, 78 of
> 36,638 records differ, all indels)

and NOT "identically".

## Pre-run audit — 10 checks, 2 real bugs caught before spending machine time

Both would have produced a plausible WRONG answer, not an error:

| # | check | result |
|---|---|---|
| 1 | span codec round-trip | PASS — 451,760 synthetic spans, 0 mismatches |
| 2 | `orig2uid` expansion | **FAIL → FIXED** — flag sense inverted, value is a DELTA not an absolute id |
| 3 | `pos_abs` element width | **FAIL → FIXED** — passed w=4; proven consumers use w=1 |
| 4 | `read_pos` semantics | PASS — contig-relative, matches `cov[cid][pos+j]` |
| 5 | strand mirroring (chains) | PASS |
| 6 | all 3 span emission sites | PASS — incl. the second-pass region |
| 7 | mapped reads / `rc_inplace` | PASS — settled at the STORAGE transform, which covers both paths |
| 8 | `read_rc` sense | PASS |
| 9 | placement coverage | PASS — `build_substrate` re-places the remainder |
| 10 | `read_clip` absence | PASS — guarded by a size check |

## Two predictions I made and got wrong, recorded

1. **"Reverse-strand reads will be mis-mapped."** Raised on the chain path,
   retracted, then found the retraction rested on insufficient evidence
   (chain path only). Settled properly at the storage transform: the encoder
   mirrors `q` for `prc[u]` reads BEFORE storing, so `pos_abs` is forward-pg for
   every read and the lookup is correct.
2. **"`parallel_loop` at 1.9 s vs 210 s means B will produce far fewer SNVs."**
   Wrong — B produced slightly MORE (36,648 vs 36,638). Inference from a stage
   timing, corrected by the measurement.
