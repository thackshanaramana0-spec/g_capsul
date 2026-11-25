# Claim 1 mechanism: variable-length support and the comparison with PgRC2

This document explains why G_CAPSUL's lossless archives are smaller than those of SPRING and
Genozip (Claim 1), what the compression stage contributes beyond the size ratio, and how the
sequence stream compares with PgRC2, the one architecturally similar tool. The numbers come from
`docs/extras/SOTA_COMPARISON.md` and the result files in `results/claim1/`. Where a source is
cited as "the development notes", it is the author's internal working record, which is not part
of this repository.

## What Claim 1 contributes besides the ratio

PgRC2 is also pseudogenome-based. On 6 of the 14 datasets of the original benchmark it fails on
variable-length reads. It either refuses (`"Unsupported variable length reads"` on SARS-CoV-2,
HCMV and C. jejuni) or crashes (`*** stack smashing detected ***` on H. pylori and
`free(): invalid next size` on A. fumigatus). G_CAPSUL processes all 14. This is a difference in
which inputs can be handled, and a ratio improvement of a few percent does not capture it. Reads
of variable length are common in practice, for example after adapter or quality trimming and in
many amplicon protocols.

## How variable length is handled

1. **Prefix-containment removal generalizes exact deduplication.** For fixed-length reads a
   duplicate is an identical string and a hash finds it. With variable length, a 3'-trimmed read
   is a strict prefix of its untrimmed twin. The two strings differ in length and hash
   differently, so a hash-based deduplication cannot see the relationship. Containment removal
   treats prefix relationships as the general case, with exact deduplication as a special case.
2. **Mismatch positions are varint-coded.** A fixed-width code cannot represent reads longer than
   256 bp, so a variable-length design needs a different encoding.
3. **Reads containing N are substituted (N to A) and pass through the same pipeline.** The original
   N positions are kept in a side stream, so there is one path for every read.

Variable-length input also triggered two of the four silent data-loss bugs found during
development. The mismatch-position byte clamp is reachable only once reads exceed 256 bp. The
orphaned-unique-read desync can occur only through containment, which exists only for
variable-length reads. Both were found by decoding an archive and comparing it with the original
file, which is the check described under "Bug history" below.

## Stream-by-stream comparison with PgRC2

Both tools solve the same core problem with similar internal stages, so this is the one
comparison that can be made at the level of individual streams. The source is Table 1 of
`docs/extras/SOTA_COMPARISON.md`. Its reasoning is used here, and its superseded headline numbers
are not.

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

The "Read order/permutation" row describes the order information of the measurement that produced
this table. In the current archive format the original read order is carried by the streams
`orig2uid_flags` and `orig2uid_vals`, and there is no separate permutation stream (see
[`../technical_architecture.md`](../technical_architecture.md), section 3.6). The remaining rows
describe streams that exist under the same names or the same roles in the current format.

**Net result.** G_CAPSUL is ahead on ratio, by a thin but real margin (+1.88% aggregate on the 7
datasets that both tools can process, sequence only; internal test only, not published). It is ahead on robustness, since PgRC2
cannot process 6 of the 14 datasets. It is behind on speed and memory, by a known margin that is
expected for an assembly-based method. The mismatch-symbol coder, 40% smaller than PgRC2's
per-symbol code, is the largest per-stream margin found.

## The time cost of Claim 1 is reused by Claims 2 and 3

Building the pseudogenome requires per-read overlap and placement, which SPRING and Genozip do not
compute in the same way. This is extra work, and it is why G_CAPSUL is the fastest tool on none
of the 19 datasets in T1.2. The placement and deviation data that this stage produces is then
available to the later claims without further archive cost. `contig_spans` (T3.1 export),
`pos_abs` and `read_lengths` (T3.2 coverage), and the query machinery (T3.3 to T3.5) are all
built from data that compression already computed. Setting `CAPS_SPANS=1` was measured to add
zero archive bytes (18,282,397 bytes with and without it). The time is paid once, at compression.

## What predicts the compression margin

Across the 7-dataset PgRC2 comparison, the fraction of a pseudogenome that lies in the second
region (the pool of reads that failed both chaining and pigeonhole mapping) correlates with the
margin at +0.607. The fraction of reads placed by mapping correlates at -0.445. H. salinarum has
a 70.5% second-region share and wins by +8.58%. S. acidocaldarius has an 8.9% share and 51.3% of
its reads mapped, which is the extreme on both axes, and it is the one loss, at -0.83%.

An earlier analysis had argued the opposite, that an oversized second region was itself the
problem, with a span penalty of about 87 KB on S. acidocaldarius. After later coder work brought
position coding to 0.996 times its information-theoretic bound, that penalty no longer exists,
and the correlation runs the other way. The loss on S. acidocaldarius has a narrower cause. The
mapper places the most reads there (51.3%, against 15 to 27% elsewhere) and generates 983,565
mismatches, so the mismatch-position stream is 23.4% of the archive. PgRC2 spends 23.2% on the
same stream, so the difference is 0.83% in PgRC2's favor on a data shape that is hard for both.

## Dataset coverage of the comparison

An earlier record of the SPRING and Genozip comparison listed 14 of 15 datasets, because
Utricularia gibba (SRR10676752) had not yet been run. It is one of the final 19 datasets, and
G_CAPSUL is smallest on it (1,650.03 MB, against 1,718.55 MB for SPRING and 2,954.90 MB for
Genozip).

## Bug history

Before the final verification, the archive was once incomplete. `refc::encode` stored the source
of a reference but not its destination gap, length and reverse-complement flag. This was found
only because the limitation of `verify_lossless.sh` was itself questioned, and it changed the
margin from a claimed +4.65% to the actual +1.90%. Four further silent data-loss bugs were found
later (two of them triggered by variable-length input, as described above). All five were found by
decoding the archives and comparing them with the original file, and none by a passing size table.
The Claim 2 and Claim 3 results were checked the same way.

## PgRC2 ratio comparison (sequence only)

PgRC2 stores only the DNA sequence, with no names, qualities or line 3. Its ratio therefore cannot
be compared with the whole-FASTQ numbers of T1.1, and the comparison is restricted to the 7
datasets that PgRC2 can process (those with reads of near-fixed length). The result is +1.88%
aggregate, with 6 wins and 1 loss (internal test only, not published). The loss is S. acidocaldarius at -0.83%, and sweeping `MINMEM`
and `MAXMAP` showed both already at interior optima on that dataset, so the loss is structural and
not a tuning miss. G_CAPSUL is about 1.7 times slower and 2.5 times heavier in memory at worst.

The deficit on S. acidocaldarius comes from the reference streams alone. It was established by
instrumenting every stream of PgRC2's own archive on this dataset, which accounts for 99.9% of its
bytes.

| Stream | G_CAPSUL (bytes) | PgRC2 (bytes) | Comment |
|---|---|---|---|
| Mismatch symbols | 124,280 | 208,234 | G_CAPSUL 40% smaller, the largest single margin |
| Mismatch positions | 731,919 | 721,906 | Within 1.4%, both near the entropy |
| References (`mem_triples`, `mem_dstgap`, `mem_len`, `mem_rc` against offsets, lengths and RC info) | 94,935 | 26,388 | 3.6 times PgRC2, larger than the whole deficit of 24,662 bytes |

Without the reference-coding gap, G_CAPSUL would win this dataset by roughly 44 KB. PgRC2 never
stores the destination of a match. It places a marker byte inside the literal stream at the
destination position, so the decoder already knows where it is when it reads the marker.
G_CAPSUL stores the destination explicitly in `mem_dstgap`. The difference is in what is coded,
and it does not reflect the quality of either coder.

### What "sequence only" means for this comparison

PgRC2's `-q` option uses quality only as an internal heuristic for dividing reads, and its encoder
has no path that compresses quality or headers (confirmed from its source). Its archive size is
therefore the complete sequence and order cost. The matching G_CAPSUL total needs six layers:

1. the sequence literal
2. the read order
3. the MEM self-match references
4. positions and strand
5. mismatch symbols of mapped (non-exact) placements
6. reads containing N

All six are needed to reconstruct every read and its position in the file, and any subset
overstates the margin. An earlier version of this comparison used only three layers, which made
E. coli look like a 39.3% win when the complete figure is 18.6%. Reads containing N had also been
dropped with no storage anywhere at an earlier stage of development. That was a data-loss bug
separate from the four listed under "Bug history".

## Summary

Against PgRC2, the only architecturally comparable tool, G_CAPSUL has a smaller sequence-only
archive (+1.88% aggregate where both can run; internal test only, not published). It processes variable-length reads, on which PgRC2
refuses or crashes for 6 of the 14 original datasets. It is slower and heavier, as is usual for
assembly-based methods. Its adaptive mismatch-symbol coder is 40% smaller than PgRC2's independent
per-symbol code. The wall-time cost seen in T1.2 (G_CAPSUL fastest on 0 of 19 datasets against
SPRING and Genozip) is the one-time price of computing the placement structure that Claims 2 and
3 reuse without further archive cost.
