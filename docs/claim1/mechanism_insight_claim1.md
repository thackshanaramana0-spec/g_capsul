# How Claim 1 connects to the paper's insight

> Mechanism of Claim 1: the contribution of the DNA-sequence stream, and a stream-by-stream comparison with PgRC2.

## Note on an earlier version of this file

An earlier version of this file stated that Claim 1 has no distinct insight and is a competitive
execution of PgRC2's pseudogenome-compression idea. That was incorrect. It compared only the
ratio and not the field-by-field technical comparison, which shows a capability difference
(detailed below).

## Claim 1's contribution: capability, not only ratio

**The headline fact, verified against `docs/extras/SOTA_COMPARISON.md`**: PgRC2 —
the only architecturally comparable tool (also pseudogenome-based) — **fails
outright on 6 of the 14 originally-locked datasets** when reads are
variable-length, either refusing (`"Unsupported variable length reads"` —
SARS-CoV-2, HCMV, C. jejuni) or hard-crashing (`*** stack smashing
detected ***` on H. pylori; `free(): invalid next size` on A. fumigatus).
**G_CAPSUL processes all 14, including the ones that crash PgRC2.** Stated
precisely in the source document: *"G_CAPSUL wins by capability, not just
ratio — a citable advantage: PgRC2 cannot even attempt 6 of the 14
datasets in this benchmark."*

**This is a capability gap, not a margin.** A compression-ratio improvement
of a few percent is a matter of degree; an entire class of real sequencing
data (anything with variable-length reads — adapter-trimmed reads, quality-
trimmed reads, many real Illumina/amplicon protocols) that the closest
comparable tool cannot process at all is a matter of kind. **This belongs
in the paper's Results/Discussion as Claim 1's primary novel contribution,
not buried as an implementation detail below the ratio numbers.**

## The mechanism — why variable length is hard, and how it is handled as a first-class case, not a patch

Verified against the methods draft and the variable-length design note:

1. **Prefix-containment removal generalizes exact deduplication.** At fixed
   read length, a duplicate is trivially detected (identical string, hash
   match). At variable length, a 3'-trimmed read is a *strict prefix* of
   its untrimmed twin — a real duplicate relationship a plain hash-based
   dedup structurally cannot see, because the two strings are different
   lengths and therefore hash differently. Containment removal treats
   prefix relationships as the general case dedup is a special case of.
2. **Mismatch positions are varint-coded**, not fixed-width, so reads
   longer than 256bp are representable at all — a different
   encoding choice than a fixed-length design would ever need to make.
3. **N-containing reads are substituted (N→A) and routed through the same
   pipeline**, with original N positions tracked in a side stream — one
   unified path for every read, not a special case carved out for reads
   containing ambiguity codes.

**Variable-length support and the bug history**: per
the failures and refuted-ideas record Part A (also cited in the development notes §6.3
and `code_mapping_claim1.md`), at least two of the four silent data-loss
bugs found earlier in this project were **triggered by variable-length input**: the mismatch-position byte clamp
(only reachable once reads exceed 256bp, which fixed-length designs never
do) and the orphaned-unique-read desync (explicitly noted in the bug's own
file: *"Trigger: containment, which can only occur with variable-length
reads"*). **Variable-length support produced two of the four most serious historical correctness bugs in
this project. Both were found and fixed by the decode-and-diff check described in the
bug-history section below.**

## Stream-by-stream comparison with PgRC2

Source: `docs/extras/SOTA_COMPARISON.md` Table 1 (its
own reasoning, not its superseded headline numbers, per its own
superseded-banner rule). Both tools solve the same core problem
(pseudogenome-based read compression) with recognizably similar internal
stages, so this is the one comparison makeable at individual-stream
granularity rather than whole-archive totals:

| Layer | PgRC2's approach | G_CAPSUL's approach | Result |
|---|---|---|---|
| Assembly admission | Both-side-overlap required | Single-pass overlap-length admission | PgRC2's stricter rule tested directly and found WORSE on all 3 datasets tested (+2,845/+45,746/+78,597 B) — the simpler rule wins |
| Duplicate handling | Native (sorted, removed pre-assembly) | Fixed via sweep-start-at-read-length (commit `6596554`) | Now equivalent in effect. **Correction, kept visible**: this fix was once reported to flip S. acidocaldarius from a loss to a win, but that measurement was against an archive later found incomplete (§6.1) — re-measured on the complete archive, S. acidocaldarius is still the one loss (see below), though the fix shrank its pseudogenome 32% regardless. |
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

## T1.2's cost is the price of T1.1's structure, and Claims 2 and 3 reuse it at no further cost

**The mechanism, as traced in the code**: building the
pseudogenome (`greedy-sweep` stage) requires computing per-read overlap and
placement, which a simpler compressor (SPRING, Genozip) does not need to do
the same way. This is measured extra work, and it is why T1.2 shows
0/19 speed wins. **The same placement/deviation data this stage produces
is then available, completely free, to every downstream claim**: Claim 3's
`contig_spans` (T3.1 export), `pos_abs`/`read_lengths` (T3.2 coverage), and
the whole query machinery (T3.3/T3.4/T3.5) are all built from data this
stage already computed for compression's own purposes. **`CAPS_SPANS=1`
was measured to add zero archive bytes** (18,282,397 B
identical with or without it). The time cost is paid once, at Claim 1, and the downstream
reuse adds no archive cost.

## What predicts the compression margin, and a self-correction about it

Verified against the size-work roadmap. Across the 7-dataset PgRC2
comparison, the fraction of a dataset's pseudogenome that ends up in the
**second region** (the append-and-reassemble pool for reads that failed
both chaining and pigeonhole mapping) correlates with margin at **+0.607**;
the fraction of reads placed via mapping correlates at **−0.445**. H.
salinarum (70.5% second-region share) wins by +8.58%; S. acidocaldarius
(8.9% share, 51.3% mapped — the extreme on both axes) is the one loss,
−0.83%.

**This inverts an earlier, since-corrected thesis.** The PgRC2 three-way-split note
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

The Claim 1 final-verdict note's body (predating the 2026-09-10 sweep)
records the SPRING/Genozip comparison as **14/15**, not 15/15 — Utricularia
gibba (SRR10676752) was on disk but not yet run, named explicitly as an
open gap. **Resolved**: SRR10676752 is one
of the final 19 rows, with CAPSULE winning it too (1,650.03 MB vs SPRING
1,718.55 MB vs Genozip 2,954.90 MB).

## Earlier bug history: the origin of the decode-and-diff check

Before the final verification: (1) the archive was once
incomplete — `refc::encode` stored a reference's source but not its
destination gap/length/RC flag, found only because `verify_lossless.sh`'s
own limitation was itself questioned; margin corrected from a
previously-claimed +4.65% to the real +1.90%. (2) Four further silent
data-loss bugs (the development notes §6.3, the failures and refuted-ideas record Part
A) — two of them specifically triggered by variable-length input, as
detailed above. **All five were found by decoding archives and
diffing against the original file and not by trusting a passing size
table.** The Claim 2 and Claim 3 findings follow the same practice.

## PgRC2 ratio comparison (sequence-only)

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
the PgRC2 stream-comparison note** (full per-stream instrumentation of
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
the locked-scope note**: confirmed directly from PgRC2's own
source (its `-q` flag only uses quality as an internal read-division
heuristic; there is no quality- or header-compression code path anywhere
in its encoder), so its one archive-size number is the complete
sequence+order cost, not an undercount. On this project's side, a
sequence+order total requires **six** real layers, and an earlier version
of this comparison used only three of them, making E. coli look like a
39.3% win when the complete number is 18.6%. The error was found internally: (1) sequence literal, (2) read
order/permutation, (3) MEM self-match references, (4) positions + strand,
(5) mismatch symbols for mapped (non-exact) placements, and (6)
N-containing reads — which, at the time this scope was locked, had been
**dropped with no storage anywhere** by every earlier stage in
the project's history. This was a data-loss bug distinct from the
four cataloged in the development notes §6.3. All six are required for a decoder to
reconstruct every read's sequence and original file position; any subset
is an incomplete number that overstates the margin.

## Summary paragraph for the discussion section

*Claim 1's contribution is not solely a compression-ratio margin. Against
PgRC2, the only architecturally comparable tool, G_CAPSUL wins on ratio
(+1.88% aggregate where both tools can run), wins on
capability (PgRC2 fails outright, by refusal or crash, on 6 of the 14
originally-tested real datasets due to variable-length reads, which
G_CAPSUL handles as a first-class case), and loses on speed and memory by a
disclosed, expected margin consistent with assembly-based methods
generally. Within individual streams, G_CAPSUL's adaptive mismatch-symbol
coder measured 40% smaller than PgRC2's independent per-symbol code, the
largest per-stream margin found in this comparison. This capability, in addition to the ratio, is a
main advantage of Claim 1: a class of
sequencing data that the closest comparable architecture cannot process at
all. Separately, T1.2's wall-time cost (CAPSULE fastest on 0 of 19 datasets
against SPRING/Genozip) is the one-time price of computing the placement
structure Claims 2 and 3 reuse without any further archive cost.*
