# Archive-path caller: speed and RAM

What `capsule_decode call` costs, what was changed to get there, how each change
was gated, and what was tried and thrown away.

Every number is one job at a time on the 12-vCPU box (CLAUDE.md rule 4), on the
4M-read HG002 chr20 subset archive unless a section says otherwise.

---

## 1. Result

|                | before   | after      |         |
|----------------|----------|------------|---------|
| wall           | 128.5 s  | **55.1 s** | 2.33x   |
| peak RSS       | 13.3 GB  | **7.06 GB**| -47%    |
| SNV F1         | 0.4539   | 0.4539     | identical |
| INDEL F1       | 0.2613   | 0.2613     | identical |

and on the FULL chr20 archive (12.6M reads, never tuned against — see section 9):

|                | before   | after       |         |
|----------------|----------|-------------|---------|
| wall           | 381.3 s  | **148.4 s** | 2.57x   |
| peak RSS       | 33.3 GB  | **16.2 GB** | -51%    |
| VCF            | 58,433   | 58,433      | byte-identical |

"before" is commit `0851692`.

The full chain, so the numbers are not conflated — only the last row is this
document's work:

| stage | wall | what changed |
|---|---|---|
| FASTQ-path caller | 671.0 s | starting point (`docs/ARCHIVE_PATH_IS_FASTER.md`) |
| archive-path caller | 438.6 s | **switched input**, not an optimisation: `parallel_loop` 214.3 s -> 2.1 s because compression already computed the read placements |
| ... | 128.5 s | earlier optimisation work, up to commit `0851692` |
| **now** | **55.1 s** | this document |

So 671 -> 55.1 s is **12.2x end to end**, of which the 438.6 -> 55.1 s part
(**8.0x**) is optimisation and the 671 -> 438.6 s part is a change of input.

Three runs at the final commit, back to back:

    run1 wall=54.89s peak=7.00GB  SNV 13858/2624/30717 F1=0.4539  INDEL 1217/318/6564 F1=0.2613
    run2 wall=55.11s peak=7.14GB  SNV 13858/2624/30717 F1=0.4539  INDEL 1217/318/6564 F1=0.2613
    run3 wall=55.32s peak=7.05GB  SNV 13858/2624/30717 F1=0.4539  INDEL 1217/318/6564 F1=0.2613

Not one TP, FP or FN moved at any point in this work.

---

## 2. What was verified, and how

`caps_caller.h` is included by the ENCODER as well as the caller, so
output-identity had to be checked on the compression side too, not assumed.

| check | method | result |
|---|---|---|
| archive unchanged | encode the same FASTQ with the encoder built at `0851692` and at final HEAD, `cmp` the archives | **byte-identical**, both configurations: 183,509,577 B bare, and 196,555,953 B with `DUMP_LIT/DUMP_PERM/DUMP_MM` (the configuration that actually writes `pos_abs`) |
| encoder VCF unchanged | same two runs, `cmp` the VCFs | **byte-identical** |
| round trip | `scripts/verify_lossless.sh` (decodes and compares against the ORIGINAL FASTQ's sequence column) | **LOSSLESS** |
| decoder streams | decode one archive with both binaries, `cmp` every output stream | identical (only the timing log differs) |
| calls | full lift + `rtg vcfeval` against GIAB truth, 3 independent runs | F1 identical, VCF 36,647 records |
| k-mer set | `[KC-SUPERK] distinct` asserted every run | 84,181,752 every run |

Per-change gating followed CLAUDE.md rule 2: every change here is
output-preserving and was held to `cmp`-identity on the VCF, not to "F1 didn't
move". Where a change touched a structure with an ordering rule, the rule was
restated in the code and the identity re-checked.

---

## 3. The changes

### Speed

| lever | stage effect | wall |
|---|---|---|
| hash-bucket directory over the seed index | 52 probes/lookup -> 2 | 128.5 -> 126.9 |
| exact branch-and-bound in the alignment scan | place 33.3 -> 11.3 s | 126.9 -> 105.1 |
| rolling canonical 25-mer at six k-mer walks | collapse 13.9 -> 7.3 s, read scan 12.8 -> 6.7 s | 105.1 -> 93.1 |
| parallel in-memory quality-as-text decode | quality 11.0 -> 5.2 s | 93.1 -> 86.5 |
| contigs handed back from step 1 | export 5.6 -> 0.0 s | 86.5 -> 81.1 |
| parallel bin-wise kc merge | merge 8.2 -> 0.7 s | 81.1 -> 73.4 |
| one-shot anchor directory (drops the Bloom filter) | read scan 6.6 -> 4.8 s | 75.1 -> 73.1 |
| word-wise skip in the alignment scan | place 11.3 -> 10.9 s | 71.8 -> 71.0 |
| `pcov` was `cov` recomputed; `cov` partitioned by contig | pcov 1.24 -> 0.09 s | 71.3 -> 70.3 |
| `pv` was `kidx` rebuilt | also -1.58 GB | 68.9 -> 66.9 |
| seed index built in place | also -2.43 GB spike | 67.2 -> 62.5 |
| `kidx` built in place | setup 3.51 -> 2.82 s | 62.2 -> 62.5 (neutral) |
| collapse claimed set sized to start where it ends | also -0.5 GB | 62.2 -> 61.1 |
| read scan fused with the anchor work | 9.17 s -> 3.88 s | 61.0 -> **55.0** |

### RAM

| lever | peak |
|---|---|
| collapse claimed set sized from what goes in | table 4.29 GB -> 1.02 GB |
| 16-byte pcluster anchor record | 13.30 -> 11.32 GB |
| 16-byte `kidx` record | 11.30 -> 10.51 GB |
| `pv` was `kidx` rebuilt | 10.25 -> 8.67 GB |
| `kidx` built in place | indel setup step +2394 -> +883 MB |
| collapse set starts where it ends | collapse step +2076 -> +1564 MB |
| read scan fused with the anchor work | 8.74 -> **7.06 GB** |

---

## 4. Three shapes, not fourteen unrelated tricks

**(a) The same array computed twice.** Three separate cases, found by reading
two builds side by side rather than by profiling:

- `pcov` is `cov` — same contigs, same placements, same 60000 clamp.
- `pv` is `kidx` — same k=25, same `Roll25`, same canon+flag packing, same
  `(ci,pos)`, both stable-sorted on the masked key.
- the contigs were obtained by re-opening the archive and rebuilding a 148 MB
  pseudogenome that step 1 already had in memory.

These were the cheapest wins available and the largest RAM ones. None of them
is visible in a flat profile: each looks like legitimate work in its own stage.

**(b) A structure whose cost is dependent memory probes.** Hash-map nodes, a
binary search over tens of millions of keys, a Bloom filter in front of a
binary search, a priority queue fed by 12-byte `fread`s. Replacing the
structure while keeping contents and ordering identical is bit-identical by
construction.

**(c) A parallel `stable_sort` that allocates a full second copy.** Three
1.70 GB arrays went through `__gnu_parallel::stable_sort`. In every case the
stability requirement was really "equal keys keep contig-then-position order",
which is what the fill already produced — so stating it in the comparator as
`(key, ci, pos)` gives the identical sequence, needs no stability, and allows a
bucketed in-place sort with no temporary at all.

The largest single speed win — branch-and-bound in the alignment scan — is none
of these three. It came from counting what the loop did (1.04 billion alignment
attempts, 93-94% of them abandoned after 19-30 bytes) rather than reasoning
about it. Two earlier guesses at that loop, both "it must be cache misses",
were wrong.

---

## 5. Tried, measured, thrown away

Recorded per CLAUDE.md rule 3 so they are not retried.

| idea | measurement | why it failed |
|---|---|---|
| huge-page hints on the collapse table | 72.01 s vs 72.49 s | inside noise; kept only because it is free and correct |
| hashing the per-thread event maps | anchor loop 4.29 -> 4.37 s | the maps are tiny — only ~300 events survive — so the loop was never map-bound |
| parallel split of the quality strings | quality 5.35 -> 5.40 s | `QBLOCK_BYTES` is 256 MB, so 590 MB of quality is **3 blocks**; the cost is `fqz_decompress`, not the split. Encoder-side parameter, not touched |
| merging kc bins directly into `kc` | peak 8.73 -> 9.42 GB | `kc` must then be allocated at the pre-dedup bound and `resize()` down does not release capacity — ~753 MB held for the rest of the run to save a ~250 MB transient |
| per-read candidate dedup in placement | never completed | reverted unmeasured rather than kept on a correctness-sensitive path |
| interleaving `scid`/`spos`, flat contig tables | place 10.54 -> 10.24 s | real but small; superseded by the in-place seed index which subsumes it |

---

## 6. Measurement traps in this pipeline

**A full disk corrupts this pipeline rather than failing it.** The k-mer
counter spills to disk and a truncated spill silently drops k-mers. One A/B
reported INDEL F1 0.1722 against 0.2613 and looked exactly like a bad code
change; the code was provably neutral (0 disagreements over 81,278,428 real
lookups) and the tell was `[KC-SUPERK] distinct=72485292` instead of
84,181,752. The harness now refuses to start below 40 GB free and asserts the
k-mer identity after every run. **Do not remove those checks.**

**Run-to-run noise is +/-1.5 s** at this scale — four consecutive runs of one
binary spanned 71.32-74.11 s. A single before/after pair cannot resolve
anything smaller. Every timing here is from interleaved A/B/A/B runs, and
sub-second effects are read off the stage timers, not the wall clock.

**An archive built without `DUMP_PERM=1` cannot drive `capsule_decode call`.**
`pos_abs`, `pos_strand` and `read_lengths` are written inside
`if(getenv("DUMP_PERM"))` in `stages/106_inprocess.cpp` (~line 3296). Without
it the caller stops with `ARCHIVE LACKS pos_abs`, and a plain decode silently
reports `reads written: 0` instead of erroring. The tell is the FILE SIZE: the
same 4M input gives 183,509,577 B bare against 196,555,953 B with the dump
flags. Use the recipe in `scripts/encode_and_call_once.sh`, not a bare
invocation. This cost one wasted 13-minute encode and a void A/B in this work,
and it is the second time the same flag has done it.

**The `[INDEL-PROF] === accounted Xs of stage ===` line must match the
`indel_pass` total.** It has caught wrongly-placed timers four times.

**Contig-space VCF identity is only valid within one input.** The caller is
nondeterministic in contig numbering across different substrates; comparing
VCFs is valid here because every gate compares two binaries on the SAME
archive.

---

## 7. Instrumentation left in place

`[BS-SPLIT]`, `[KC-SPLIT]`, `[INDEL-PROF]`, `[DEC-RSS]` — each prints elapsed
time plus `VmRSS`/`VmHWM` at every phase boundary. The RAM work was only
possible because the high-water mark is attributable to a phase; three of the
seven RAM changes were aimed at a step that a flat profile does not show at
all. The cost is a handful of `/proc/self/status` reads per run.

---

## 8. Where the remaining 7.06 GB sits, and what is closed

No single dominant step is left; the peak is now reached by accumulation:

    decode leftovers (seqs + quals)   1356 MB
    collapse                          +1564 MB   (1024 MB claimed-set transient)
    kc merge                          +2056 MB
    seed index                         +747 MB
    indel setup (cov + kidx)           +883 MB
    anchor directory                   +367 MB
                                      -------
                                       7228 MB

**Quality as a Q>=20 bitmap is CLOSED, not untried.** It would take `quals`
from 768 MB to 74 MB, and quality is read in exactly one place in this path
(`medq < 20`, the median of the minor allele's base qualities). But the median
of an EVEN-sized list is the average of the two middle values, so `(19, 21)`
passes at 20.0 while `(19, 20)` fails at 19.5 — and those two have identical
`>=20` bit patterns. Bits cannot decide the even case, so a bitmap changes
calls. The DBG path at line ~3675 already uses packed bitmaps precisely because
the only question asked there IS a threshold; that reasoning does not carry to
the pileup.

**Still open, with the arithmetic:**

- `seqs` and `quals` are `vector<std::string>` — 4M x (32 B object + ~176 B
  heap) = ~830 MB each to hold ~590 MB of data. A flat buffer plus offsets
  would save ~200 MB each, at the cost of touching every `seqs[i]` use site.
- the kc merge's +2056 MB could drop to ~1.45 GB with a two-pass scheme (count
  distinct per bin, allocate exactly, then dedup into place) for about +1.5 s.
  The one-pass in-place version was tried and is refuted above.
- the 1024 MB collapse claimed-set transient is near-minimal for an exact set:
  60.7M distinct 64-bit keys is 486 MB at 100% load, and open addressing needs
  a power-of-two table.

---

## 9. Does it generalise? Full chr20, 3.15x the tuning subset

Everything above was measured on a 4M-read subset. That is the set the work was
done against, so it cannot also be the evidence that the work generalises. The
same two binaries were therefore run on the FULL HG002 chr20 archive --
12,604,917 reads, 3.15x as many, 592 MB archive, encoded with the required
`DUMP_PERM=1`.

|            | pre-session (`0851692`) | final    |        |
|------------|-------------------------|----------|--------|
| wall       | 381.30 s                | **148.39 s** | 2.57x |
| peak RSS   | 33.28 GB                | **16.21 GB** | -51%  |
| VCF        | 58,433 records          | 58,433   | **byte-identical** |
| kc distinct| 140,719,632             | 140,719,632 | identical |

**The wins hold, and are slightly larger at scale than on the subset** (2.57x
vs 2.33x, -51% vs -47%). That is the expected direction: several of the changes
remove an allocation or a pass that scales with input size, so the bigger the
input the more they save. Nothing here was tuned to 12.6M reads -- this run is
the first time that archive was used.

The RAM result is the one that changes what is possible: 33.28 GB is most of a
36 GB machine and over a third of this 82 GB box, while 16.21 GB fits
comfortably alongside other work.

Stage peaks at full scale keep the same shape as on the subset, so the profile
is not an artefact of the small input:

    decode leftovers                 4134 MB
    collapse                        +3690 MB
    kc merge                        +2209 MB
    seed index                      +2853 MB
    indel setup (cov + kidx)        +1599 MB
    anchor directory                 +431 MB
                                   --------
                                    16594 MB

(The script printed `WARNING: pos_abs not visible in archive`; that is a false
alarm from a crude `strings | head` probe, not a real condition -- both
decoders read the archive and produced 58,433 records with the expected k-mer
count.)

---

## 10. The PUBLISHED path is a different path, and it needed separate work

Everything in sections 1-9 is `CAPS_CALL_INDELS=1`, the full SNV+indel caller.
**Claim 2's published configuration is not that.** It is `CAPS_DBG_ONLY`, which
is the DEFAULT when `CAPS_CALL_INDELS` is unset, and it takes a different route
through the caller -- Method B skips `build_substrate` entirely, so the seed
index, the placement scan and the collapse set (the three biggest levers above)
do not run at all there.

The competitive picture, from `docs/INDEL_PASS_SPEED_PLAN.md` and
`docs/COMPRESSION_DERIVED_CALLING.md`:

| | wall | RAM | SNV F1 |
|---|---|---|---|
| DiscoSNP++ | 75.7 s | 3.29 GB | 0.847 |
| ours, `CAPS_DBG_ONLY` | 106.0 s | 6.03 GB | **0.8766** |

Accuracy wins; speed and RAM do not. Measured on the full chr20 archive on this
box, the pre-session binary is 140.5 / 139.6 s at 6.34 / 6.41 GB.

### 10.1 VCF byte-identity is an INVALID gate on this path

The pre-session binary was run twice on the same archive and **differs from
itself**: 223,524 differing records out of 111,766, and 230,966 differing
`##contig` header lines. That is wholesale `dcontig_N` relabeling -- the
parallel bubble traversal enumerates bubbles in nondeterministic order, so
every record's contig name shifts even when the calls are the same.

So a `cmp` difference between two binaries here means nothing on its own, and
one was observed and initially misread as a regression. **Gate this path on F1
after lifting to chr20 coordinates, never on the VCF bytes.** (The indel path
IS stable and was correctly gated on `cmp` -- the two paths differ in this.)

### 10.2 A win on one path was a regression on the other

The parallel bin-wise kc merge (section 3) is an 8.2 -> 0.7 s win on the indel
path. On the graph-only path it is a **RAM regression of ~1.3 GB**, because the
serial k-way merge it replaced STREAMED from disk holding almost nothing, while
the bin version builds the whole deduped result in staging vectors and then
copies it into `kc` -- resident twice. That is invisible where the peak is
elsewhere and decisive where `kc` IS the peak.

Fixed by making the merge two-pass with exact allocation: pass A sorts each bin
and counts its distinct keys, retaining nothing, which gives the exact final
size so `kc` is allocated once at exactly that size; pass B redoes the read and
sort and writes deduped runs straight into place. Cost is one extra read+sort
of cache-resident bins; the merged array is never resident twice.

**The lesson is the gating, not the bug.** Every change in section 3 was gated
on the indel path only, so a change that helped there and hurt elsewhere passed
cleanly. A caller with two configurations needs both measured.

### 10.3 Result on the published path

Two fixes: the kc merge rebuilt as two-pass with exact allocation (10.2), and
the in-memory contigs handover extended to this path in its correct TWO-RECORD
form -- the graph path's ploidy gate samples the concatenated pseudogenome, not
the 451k individual contigs, so handing it the per-contig form would have been
a different input rather than a faster route to the same one.

Full chr20, `CAPS_DBG_ONLY`, two runs each:

| | wall | peak |
|---|---|---|
| pre-session (`0851692`) | 140.50 / 139.63 s | 6.34 / 6.41 GB |
| session before these fixes | 131.63 / 129.99 s | 7.96 / 7.89 GB |
| **with both fixes** | **116.09 / 116.26 s** | **6.67 / 6.69 GB** |

**-17% wall against the pre-session baseline, and the RAM regression closes
from +1.55 GB to +0.30 GB.** The residual +0.30 GB is the per-thread bin
buffers the parallel merge needs and the serial streaming merge did not; it
buys the merge going from 8.2 s to ~1.5 s and is the trade being made
knowingly.

    export pseudogenome   13.94 s -> 0.00 s
    kc k-way merge         1.08 s -> 1.52 s   (the deliberate extra pass)

The indel path is unaffected: **VCF byte-identical** to the verified HEAD run,
53.83 s / 7.01 GB.

### 10.4 Where this leaves the comparison

| full chr20 | wall | RAM | SNV F1 |
|---|---|---|---|
| DiscoSNP++ | 75.7 s | 3.29 GB | 0.847 |
| ours, published path, now | 116.2 s | 6.68 GB | **0.8766** |

**We lead on accuracy and do not lead on cost.** ~1.5x slower and ~2x heavier.
This session closed a third of the speed gap and none of the RAM gap.

The remaining speed is not in anything touched here: of the 116 s, decode is
26.7 s, `kc_H_build` ~12 s, and **the DBG traversal is ~74 s**. That is Method
B's core algorithm and it is untouched. Any serious attempt at DiscoSNP's wall
time has to go there, and it is a larger piece of work than the structural
cleanups in this document.

The RAM gap has the same shape: DiscoSNP holds 3.29 GB because GATB streams its
k-mer partitions and never materialises the whole counter set, while `kc` here
is 140,719,632 x 16 B = 2.25 GB resident plus ~4.4 GB of decoded reads and
quality. Closing that means streaming kc, not shaving allocations.

### 10.5 The 1-FP difference: chased to the end, and it is not a call change

The F1 gate showed one extra false positive (2655 -> 2656) with TP and FN
identical. The tempting reading was "nondeterminism, ignore it". That was
tested and REFUTED -- repeat runs are stable per binary (2655/2655 and
2656/2656), so the difference is real and attributable.

The chase, each step measured rather than argued:

| step | result |
|---|---|
| Is it the kc merge? | **No.** `dec_final` (has the merge, not the handover) gives 2655, same as pre-session |
| Did kc change? | **No.** kc nodes 140,719,632, branching 1,226,535, bubbles 115,842, anchored 591,248 -- identical |
| Did emission change? | **No.** emitted=101,984 with identical filter drops (covcap dropped=13,858) |
| Is the FASTA round trip lossy? | **No.** 2 records, 457,043,834 bytes = PG_LEN exactly, zero non-ACGT |
| Are the dcontig sequences the same? | **Yes as a SET** (sorted md5 identical, 115,842 both) but in a **different ORDER** (unsorted md5 differs) |
| Is the call set the same? | **YES, proven.** Keying every record by its contig's SEQUENCE instead of its label: 111,766 records each, exact match including INFO |

So the caller's output is unchanged. What differs is only the `dcontig_N`
numbering, which changes the order contigs reach `bwa`, which changes which
duplicate `lift_vcf.py` keeps at 2-3 positions -- visible as ALT disagreements
at the SAME coordinate (e.g. 20:26,260,686 C>A against C>T) and a lifted record
count of 47,332 against 47,334.

**Two things follow, and the second matters beyond this change.**

1. The change is output-preserving at the strongest level this path admits.
2. **The published F1 on this path carries a +/-1-FP sensitivity to bubble
   enumeration order that is inherent to the tool.** The baseline binary has it
   too -- it differs from ITSELF in 223,524 VCF records. Any future comparison
   here must key on contig sequence, not on `dcontig_N`, or it will chase this
   artefact. `docs/` already warned that VCF byte-identity is not a valid gate;
   this is the concrete mechanism.
