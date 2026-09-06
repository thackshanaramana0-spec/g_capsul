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

**Prime suspect: `read_clip` has no archive source.** The encoder derives it
during placement rather than storing it, so the archive path passes zeros.
`build_substrate` guards it (`o < cd.read_clip.size() ? ... : 0`) so zero
degrades gracefully rather than corrupting — which is consistent with a small
difference rather than a large one.

**Fix:** store `read_clip` as a stream (the encoder has it at placement time)
and re-run. If the 78 vanish, the claim becomes exact.

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
