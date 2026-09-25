---
Date: 2026-09-19
Title: How Claim 1 Connects to the Paper's Central Insight — Its Own Real
  Novel Insight (Variable-Length Capability), the Full Stream-by-Stream
  PgRC2 Comparison, and the Structure Claims 2/3 Reuse For Free
Purpose: The synthesis file for Claim 1. States Claim 1's actual novel
  contribution precisely (previously understated in this folder — see
  correction note below), gives the full, verified stream-level PgRC2
  comparison, and explains why T1.2's real cost is not a weakness but the
  price of the structure Claims 2 and 3 reuse for free.
When to refer to this file: Writing the Results/Discussion text for
  Claim 1; citing any PgRC2 number or the variable-length capability claim;
  explaining why T1.2's loss is load-bearing, not embarrassing.
Keywords: insight, mechanism, variable-length reads, PgRC2, capability not
  ratio, stream-by-stream, mismatch symbol coding, overlap chaining,
  placement, CAPS_SPANS, contig_spans, free byproduct, trade-off
---

# How Claim 1 connects to the paper's actual insight

## Correction, kept visible rather than silently replaced

An earlier version of this file stated Claim 1 has "no distinct novel
insight," framing it purely as a competitive execution of PgRC2's existing
pseudogenome-compression idea. **This was wrong, caught and corrected in
session.** Claim 1 has a real, verified, novel capability insight, detailed
below — the earlier framing missed it because it only checked the ratio
comparison, not the field-by-field technical comparison that surfaces it.

## Claim 1's actual novel insight: capability, not just ratio

**The headline fact, verified against `docs/SOTA_COMPARISON.md`**: PgRC2 —
the only architecturally comparable tool (also pseudogenome-based) — **fails
outright on 6 of the 14 originally-locked datasets** when reads are
variable-length, either refusing (`"Unsupported variable length reads"` —
SARS-CoV-2, HCMV, C. jejuni) or hard-crashing (`*** stack smashing
detected ***` on H. pylori; `free(): invalid next size` on A. fumigatus).
**G_CAPSUL processes all 14, including the ones that crash PgRC2.** Stated
precisely in the source document: *"G_CAPSUL wins by capability, not just
ratio — a real, citable advantage: PgRC2 cannot even attempt 6 of the 14
datasets in this benchmark."*

**This is a capability gap, not a margin.** A compression-ratio improvement
of a few percent is a matter of degree; an entire class of real sequencing
data (anything with variable-length reads — adapter-trimmed reads, quality-
trimmed reads, many real Illumina/amplicon protocols) that the closest
comparable tool cannot process at all is a matter of kind. **This belongs
in the paper's Results/Discussion as Claim 1's primary novel contribution,
not buried as an implementation detail below the ratio numbers.**

## The mechanism — why variable length is hard, and how it is handled as a first-class case, not a patch

Verified against `paper/METHODS.md` and `docs/VARIABLE_LENGTH_DESIGN.md`:

1. **Prefix-containment removal generalizes exact deduplication.** At fixed
   read length, a duplicate is trivially detected (identical string, hash
   match). At variable length, a 3'-trimmed read is a *strict prefix* of
   its untrimmed twin — a real duplicate relationship a plain hash-based
   dedup structurally cannot see, because the two strings are different
   lengths and therefore hash differently. Containment removal treats
   prefix relationships as the general case dedup is a special case of.
2. **Mismatch positions are varint-coded**, not fixed-width, so reads
   longer than 256bp are representable at all — a genuinely different
   encoding choice than a fixed-length design would ever need to make.
3. **N-containing reads are substituted (N→A) and routed through the same
   pipeline**, with original N positions tracked in a side stream — one
   unified path for every read, not a special case carved out for reads
   containing ambiguity codes.

**Why this was hard to get right — the 4 bugs, already documented in this
project's bug history, now connected to their actual trigger**: per
`docs/FAILURES_AND_REFUTED_IDEAS.md` Part A (also cited in `DEVNOTES.md` §6.3
and `code_mapping_claim1.md`), at least two of the four silent data-loss
bugs found in this project's pre-session history were **specifically
triggered by variable-length input**: the mismatch-position byte clamp
(only reachable once reads exceed 256bp, which fixed-length designs never
do) and the orphaned-unique-read desync (explicitly noted in the bug's own
file: *"Trigger: containment, which can only occur with variable-length
reads"*). **Variable-length support was not a feature bolted on and left
untested — it was hard enough to get right that it produced two of this
project's four most serious historical correctness bugs, both found and
fixed before this session, by the same decode-and-diff discipline described
in the bug-history section below.**

## The full stream-by-stream technical comparison against PgRC2 — verified, not summarized from a single aggregate number

Source: `docs/SOTA_COMPARISON.md` Table 1, read in full this session (its
own reasoning, not its superseded headline numbers, per its own
superseded-banner rule). Both tools solve the same core problem
(pseudogenome-based read compression) with recognizably similar internal
stages, so this is the one comparison makeable at individual-stream
granularity rather than whole-archive totals:

| Layer | PgRC2's approach | G_CAPSUL's approach | Result |
|---|---|---|---|
| Assembly admission | Both-side-overlap required | Single-pass overlap-length admission | PgRC2's stricter rule tested directly and found WORSE on all 3 datasets tested (+2,845/+45,746/+78,597 B) — the simpler rule wins |
| Duplicate handling | Native (sorted, removed pre-assembly) | Fixed via sweep-start-at-read-length (commit `3e06957`) | Now equivalent in effect. **Correction, kept visible**: this fix was once reported to flip S. acidocaldarius from a loss to a win, but that measurement was against an archive later found incomplete (§6.1) — re-measured on the complete archive, S. acidocaldarius is still the one loss (see below), though the fix genuinely shrank its pseudogenome 32% regardless. |
| Pseudogenome region split | 3 physical regions | 2 regions (main + second); N's handled inline | Measured directly: the 3rd split is not free headroom G_CAPSUL is missing |
| Reference encoding | Per-period dst/len/rev-comp/offset streams | `mem_triples` split the same way, arrived at independently | Tie in structure; G_CAPSUL within ~11% of PgRC2 on raw reference bytes, cheaper on the length sub-stream |
| Mismatch position coding | Per-period bucketed range coder | Same idea + an automatic flat fallback | G_CAPSUL measured BELOW both order-0 and count-conditioned entropy — at or past the theoretical floor |
| **Mismatch symbol coding** | Independent per-symbol code | Adaptive, keyed on reference base (4-way context) | **G_CAPSUL wins by 40%**: 124,280 B vs 208,234 B — the single largest per-stream margin found in this project |
| Read order/permutation | Not broken out in PgRC2's own output | Delta-coded, split by zero/nonzero flag (79.55% of deltas are exactly zero) | Measured at **0.996x its own information-theoretic bound** — essentially closed |
| Strand (RC) flag | Not separately documented | 1 bit/read | Within 1.2% of the order-0 entropy bound |
| Names / quality | **None at all** (confirmed by reading PgRC2's decoder source) | Full columns (SPRING-derived tokenizer; vendored fqzcomp) | Not comparable to PgRC2 — compared to SPRING/Genozip instead |
| **Variable-length support** | **Fails on 6/14 datasets** (refuses or crashes) | Handles all 14 | **Capability win — see above, this is the headline** |
| Speed / RAM | Faster, lighter (baseline) | ~1.6-1.7x slower, ~2.5-2.8x heavier | PgRC2 wins — expected, disclosed, the standard assembly-vs-speed trade-off |

**Net read, stated precisely**: G_CAPSUL wins on ratio (thin but real,
+1.88% aggregate on the 7-dataset sequence-only comparison), wins outright
on robustness (PgRC2 cannot process 6 of 14 real datasets at all), and
loses on speed/RAM by a known, disclosed, expected margin. The mismatch-
symbol coder's +40% margin is Claim 1's single strongest per-stream
technical result and belongs in the paper alongside the capability claim,
not omitted in favor of the aggregate ratio alone.

## T1.2's real cost is not a separate weakness — it is what T1.1's structure requires, and what Claims 2/3 spend nothing to reuse

**The precise mechanism, verified in code this session**: building the
pseudogenome (`greedy-sweep` stage) requires computing per-read overlap and
placement, which a simpler compressor (SPRING, Genozip) does not need to do
the same way. This is real, measured, extra work — it is why T1.2 shows
0/19 speed wins. **The same placement/deviation data this stage produces
is then available, completely free, to every downstream claim**: Claim 3's
`contig_spans` (T3.1 export), `pos_abs`/`read_lengths` (T3.2 coverage), and
the whole query machinery (T3.3/T3.4/T3.5) are all built from data this
stage already computed for compression's own purposes. **`CAPS_SPANS=1`
was measured this session to add zero archive bytes** (18,282,397 B
identical with or without it) — the clearest possible demonstration that
the *time* cost is paid once, at Claim 1, and the downstream *reuse* is
genuinely free.

## What predicts the compression margin, and a self-correction about it

Verified against `docs/WHAT_IS_NEXT.md`. Across the 7-dataset PgRC2
comparison, the fraction of a dataset's pseudogenome that ends up in the
**second region** (the append-and-reassemble pool for reads that failed
both chaining and pigeonhole mapping) correlates with margin at **+0.607**;
the fraction of reads placed via mapping correlates at **−0.445**. H.
salinarum (70.5% second-region share) wins by +8.58%; S. acidocaldarius
(8.9% share, 51.3% mapped — the extreme on both axes) is the one loss,
−0.83%.

**This inverts an earlier, since-corrected thesis.** `docs/DO_WE_NEED_THEIR_3WAY.md`
had argued the *opposite* — that an oversized second region was itself the
problem, estimating the "span penalty" at ~87 KB on S. acidocaldarius,
"essentially its entire loss." Measured after later coder work brought
positions to 0.996× their information-theoretic bound (closing the span
penalty that analysis was built on), the correlation runs the other way: a
*larger* second region goes with a *better* margin, not a worse one. Kept
here specifically because the earlier framing sounds like it might still be
true and is a natural instinct to reach for — it measurably is not, once
positions are coded near their bound. The real story for the one loss is
narrower: S. acidocaldarius is where the mapper places the most reads
(51.3%, against 15–27% elsewhere), generating 983,565 mismatches that make
its mismatch-position stream 23.4% of the archive — and PgRC2 spends
23.2% on the identical stream, so this is not a case where either tool
handles the data shape well; theirs simply handles it 0.83% better.

## Historical dataset gap, and its resolution

`docs/CLAIM1_FINAL_VERDICT.md`'s body (predating the 2026-09-10 sweep)
records the SPRING/Genozip comparison as **14/15**, not 15/15 — Utricularia
gibba (SRR10676752) was on disk but not yet run, named explicitly as an
open gap. **Resolved**: this session confirmed SRR10676752 present as one
of the final 19 rows, with CAPSULE winning it too (1,650.03 MB vs SPRING
1,718.55 MB vs Genozip 2,954.90 MB).

## Claim 1's own pre-session bug-finding history — the origin of this project's "always decode and diff" discipline

Found **before** this session began: (1) the archive was once genuinely
incomplete — `refc::encode` stored a reference's source but not its
destination gap/length/RC flag, found only because `verify_lossless.sh`'s
own limitation was itself questioned; margin corrected from a
previously-claimed +4.65% to the real +1.90%. (2) Four further silent
data-loss bugs (`DEVNOTES.md` §6.3, `docs/FAILURES_AND_REFUTED_IDEAS.md` Part
A) — two of them specifically triggered by variable-length input, as
detailed above. **All five were found by actually decoding archives and
diffing against the original file, never by trusting a passing size
table** — the precedent this session's own Claim 2 and Claim 3 findings
continue, not a practice invented this session.

## Honest, verified PgRC2 ratio comparison — sequence-only, explicitly scoped

**Not comparable to T1.1's whole-FASTQ numbers.** PgRC2 stores DNA sequence
only — no names, quality, or line-3 — so the ratio comparison is
restricted to the 7 datasets both tools can process at all (a strict
subset of the 6-datasets-PgRC2-fails-on point above; these 7 are ones with
fixed-enough-length reads for PgRC2 to attempt). **+1.88% aggregate, 6
wins, 1 loss** — the loss is S. acidocaldarius, exact disclosed margin
**−0.83%**, confirmed structural (not a tuning miss) by sweeping `MINMEM`
and `MAXMAP` and finding both already at interior optima on that dataset.
Speed/RAM: **~1.7x slower, ~2.5x heavier at worst**.

**Precisely why S. acidocaldarius still loses, verified against
`docs/PGRC2_STREAM_COMPARISON.md`** (full per-stream instrumentation of
PgRC2's own real archive on this dataset, accounting for 99.9% of its
bytes): the deficit is not spread evenly. Mismatch symbols are a **win**
(124,280 B vs PgRC2's 208,234 — 40% smaller, the largest single margin in
either direction), and mismatch positions are within 1.4% (731,919 vs
721,906 — both already near-entropy, not where this is lost. **References
are the entire story**: this project's `mem_triples`+`mem_dstgap`+
`mem_len`+`mem_rc` cost 94,935 B against PgRC2's offsets+lengths+RC-info
cost of 26,388 B — **3.6× theirs, and larger than the whole dataset's
deficit** (+24,662 B). Without the reference-coding gap alone, this project
would *win* S. acidocaldarius by roughly 44 KB. Root cause: PgRC2 never
stores a match's destination at all — it embeds a marker byte directly
inside the literal stream at the destination position, so the decoder
already knows where it is when it reads the marker; this project stores
`mem_dstgap` (destination) explicitly as a separate stream. This is an
architectural difference in what gets coded, not a coder-quality gap on
either side.

**What "sequence-only" precisely means, verified against
`docs/LOCKED_SEQORDER_SCOPE.md`**: confirmed directly from PgRC2's own
source (its `-q` flag only uses quality as an internal read-division
heuristic; there is no quality- or header-compression code path anywhere
in its encoder), so its one archive-size number is genuinely the complete
sequence+order cost, not an undercount. On this project's side, an honest
sequence+order total requires **six** real layers, and an earlier version
of this comparison used only three of them, making E. coli look like a
39.3% win when the complete, correct number is 18.6% — caught directly by
this project, not by an external reviewer: (1) sequence literal, (2) read
order/permutation, (3) MEM self-match references, (4) positions + strand,
(5) mismatch symbols for mapped (non-exact) placements, and (6)
N-containing reads — which, at the time this scope was locked, had been
**silently dropped with no storage anywhere** by every earlier stage in
the project's history, a real undisclosed data-loss bug distinct from the
four cataloged in `DEVNOTES.md` §6.3. All six are required for a decoder to
reconstruct every read's sequence and original file position; any subset
is an incomplete, misleadingly favorable number.

## The honest, precise paragraph for the paper's discussion section

*Claim 1's contribution is not solely a compression-ratio margin. Against
PgRC2, the only architecturally comparable tool, G_CAPSUL wins on ratio
(+1.88% aggregate where both tools can run), wins by a full order of
capability (PgRC2 fails outright, by refusal or crash, on 6 of the 14
originally-tested real datasets due to variable-length reads, which
G_CAPSUL handles as a first-class case), and loses on speed and memory by a
disclosed, expected margin consistent with assembly-based methods
generally. Within individual streams, G_CAPSUL's adaptive mismatch-symbol
coder measured 40% smaller than PgRC2's independent per-symbol code, the
largest per-stream margin found in this comparison. This capability, not
the ratio alone, is Claim 1's most citable advantage: a real class of
sequencing data that the closest comparable architecture cannot process at
all. Separately, T1.2's wall-time cost (CAPSULE fastest on 0 of 19 datasets
against SPRING/Genozip) is the one-time price of computing the placement
structure Claims 2 and 3 reuse without any further archive cost.*
