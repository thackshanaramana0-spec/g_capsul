# METHODS — the validated architecture, as implemented

**Consolidated 2026-09-10** from the final code, the final benchmark harness and
the executed results. Where an older document disagreed with the code, the code
won; where the code disagreed with an executed result, the result won.

Results and their provenance live in `../benchmark/`, not here. This file
describes *what the system does*; `../benchmark/documentation/RESULT_CODE.md`
proves *what it produced*.

**Siblings in this folder:** [PAPER.md](PAPER.md) (start here) ·
[ARCHITECTURE.md](ARCHITECTURE.md) (the same system in diagrams) ·
[RESULTS.md](RESULTS.md) · [DISCUSSION.md](DISCUSSION.md) ·
[LIMITATIONS.md](LIMITATIONS.md).

---

## 1. Compression — the pipeline

The compressor assembles reads into a **pseudogenome** by greedy suffix-prefix
overlap and maps the remainder onto it — the approach introduced by PgRC
[Grabowski & Kowalski 2020; PgRC2 2025], itself building on assembly-based read
compression from Quip [Jones et al. 2012]. The greedy shortest-common-superstring
and q-gram matching primitives are standard.

What is built on that approach here is a substantially different system.

**Scope.** PgRC emits DNA only and cannot reproduce a FASTQ. This is a complete
lossless archive: identifiers, line 3 and quality as well as sequence, carried
in **19 independently coded streams**, each assigned a coder by measurement
rather than by rule — the selector evaluates LZMA, PPMd7, FSE, an adaptive
range coder, byte-plane-split variants, a chunked mode and a constant-run
encoding, and keeps whichever is smallest for that stream's actual statistics.

**Read geometry.** PgRC2 requires constant-length reads at most 255 bp, and
refuses or crashes otherwise — it cannot attempt 6 of the 14 datasets in an
earlier iteration of this benchmark, three refused outright
(`Unsupported variable length reads`) and two by memory corruption. Variable
length is supported here as a first-class case: prefix-containment removal
generalises exact deduplication (a trimmed read is a strict prefix of its
untrimmed twin, which hash dedup cannot see), mismatch positions are varint-
coded so reads beyond 256 bp are representable, and reads containing N are
substituted N→A and routed through the same pipeline with their positions in a
side stream rather than into a separate pseudogenome.

**Construction.** PgRC2 maintains three pseudogenomes; we use two regions.
Overlaps are found by a parallel hash sweep rather than their sort-merge — 11×
more comparisons, 3.2× faster on 12 cores. Our round-1 division is the
structural equivalent of their high/low-quality split (theirs is inactive by
default in the released binary, so both are in practice topological) and leaves
a main region 81,737 bases smaller on *S. acidocaldarius*. The mapping
tolerance adapts to measured leftover fraction rather than a fixed ratio.

**Coders.** Read order is coded as a Lehmer permutation over a Fenwick tree at
18.49 bits/read against PgRC2's 22.82, within 0.022% of the information-
theoretic floor. Mismatch symbols are conditioned on the reference base rather
than coded independently, 40% below theirs. The DNA coder is multi-order
context mixing with SSE. Identifiers extend SPRING's positional tokenizer with
a self-gating wide delta and a token that resolves against the archive's own
stored read length instead of re-coding it.

No source is shared with PgRC2, which is GPL-3. And Claims 2 and 3 have no
counterpart in that lineage at all: they depend on the pseudogenome being
retained and served, which is precisely what it is built to discard.

Implemented in `stages/106_inprocess.cpp` (the encoder) and
`stages/capsule_decode.cpp` (the decoder), with stream coders in
`include/coders_inproc.h` and the shared DNA coder in `include/seqpar_core.h`.

In order:

1. **Load, filter, optional dedup.** Duplicate fraction is measured; dedup is
   enabled only when it pays.
2. **Greedy suffix-prefix overlap chaining** builds a pseudogenome from
   well-tiling reads. The sweep runs from `Lmax` downward and revisits every
   still-open end at every overlap length, longest first — a read whose best
   partner is already taken still gets its next-best one length down.
3. **Pigeonhole mapping** places the remaining reads onto that pseudogenome,
   accepting up to `MAXMAP` mismatches, keeping the *best* placement rather
   than the first acceptable one.
4. **Second region.** Reads that still do not place are appended and assembled
   separately.
5. **MEM self-match** removes redundancy within the pseudogenome.
6. **Stream emission and entropy coding.** Sequence, positions, strands, read
   lengths, mismatch positions and symbols, names, quality and line 3 are coded
   as separate streams, each with its own coder chosen by a selector.

The archive is a lossless FASTQ archive: **DNA, read names, line 3 and quality**.

### Two design decisions that are load-bearing

**The sweep starts at `Lmax`, not `Lmax−1`.** A suffix-prefix overlap of
exactly the read length *is* an exact duplicate, and that is the only length at
which one can appear. Starting one below made duplicates invisible to chaining,
forcing a pre-assembly dedup pass. Starting at `Lmax` cut one pseudogenome from
9,134,100 to 6,157,270 bytes and turned the project's last size loss into a
win.

**Placement keeps the best match, and the current best doubles as the
early-exit bound.** A candidate is abandoned the moment it is worse than what
the read already has, so better placements also mean less work.

## 2. Variant calling — from the archive

`capsule_decode call <archive> <out.vcf>` — no FASTQ, no reference, no separate
assembly graph. Implemented in `include/caps_caller.h`.

### What is reused, and what that is worth

Reference-free callers build their own substrate from raw reads: DiscoSNP++
constructs a GATB de Bruijn graph, then runs `kissreads2`, a **separate binary**
that re-maps every read onto every candidate bubble. Here the assembly already
exists — it is what the compressor built in order to compress — and the caller
consumes the encoder's contigs directly. The alternative was implemented and
measured: re-deriving the substrate re-places all 12.6 M reads and costs **738 s
serially at full chr20**, against reusing what the archive already holds.

Read coherence is likewise done in one indexed sweep over reads held in memory,
which makes a **per-base** quality test affordable where a separate mapping tool
uses a per-path mean: quality enters as a 1-bit-per-base bitmap at 233 MB rather
than 2.27 GB of phred strings.

Two passes are added, **to the caller only; the compression path is untouched**:

1. **`collapse_contigs()`** — a greedy longest-first non-redundant contig set,
   putting both haplotypes of a locus into one frame.
2. **Mismatch-tolerant re-placement** of every read onto the surviving contigs,
   both strands, fewest mismatches wins.

### Why both passes are necessary — the mechanism

A pseudogenome built to minimise bits gives each haplotype of a heterozygous
site **its own contig**, because two internally-consistent contigs compress
better than one contig plus a column of disagreements. The ref-allele and
alt-allele reads are then never placed at the same coordinate, and the variant
is **not present in the data structure at all**. No read-out layer can recover
it; the pseudogenome has to be re-framed.

Measured, full chr20, HG002, called from the archive:

| configuration | F1 | precision | recall |
|---|---|---|---|
| neither pass | 0.431 | 0.967 | 0.278 |
| re-placement only | **0.426** | 0.962 | 0.274 |
| collapse only | 0.648 | 0.963 | 0.488 |
| **both** | **0.888** | 0.956 | 0.830 |

**The passes are not additive, and that is the evidence they are mechanistic
rather than fitted.** Collapse alone gains +0.217; re-placement alone *loses*
0.005 — it is worse than doing nothing. Additivity predicts 0.212; the measured
joint effect is **+0.457**. Neither pass is a filter that happens to help:
collapse puts the two haplotypes into one frame, re-placement puts the reads
onto that frame, and either alone leaves the other half of the operation
undone.

The **error shape** independently confirms the cause. Without collapse,
precision holds at 0.96 while recall falls to 0.27 — the caller is not
mistaken, it is **blind**, which is exactly what "the alt-allele reads are
filed on a different contig" predicts. A tuning artefact or a bad threshold
would cost precision instead; this costs only sensitivity.

**The correction costs nothing in compression ratio**, because it exists only
in the caller's view of the archive, not in the archive.

### Parameters are functions of measured depth, not fitted constants

The coverage ceiling is `2 × ploidy × H` and the coherence floor `max(2, H/10)`,
where `H` is the sample's own haploid depth read off the k-mer histogram. This
is why precision holds between 0.911 and 0.950 across a 10–30× sweep, and why
the tetraploid arm works without retuning. The comparable tool exposes `-b`,
`-P`, `-D` and `-max_ambigous_indel` as fixed numbers.

### Multi-allelic sites

`CAPS_PLOIDY=k` admits up to *k* co-occurring alleles and emits a native
multi-allelic VCF record. Across DiscoSNP++'s entire output for the same
region — 3,989 records — it emits **zero** with more than one ALT allele. That
is a property of its output model, not a miss, and it is why this is reported
as a capability rather than as a rate.

## 3. Addressability — export, coverage, query

Served directly from the archive: the work a conventional pipeline does at
query time — assembling, or aligning and indexing — was already done at
compression time. Each operation decodes only the streams it needs and stops.

| operation | decodes | never touches | cost |
|---|---|---|---|
| `export` | `literal` + `mem_triples` (+ companions) → the pseudogenome | any per-read stream | 1 pg rebuild, 36 MB |
| `coverage` | `pos_abs` + `read_lengths` → per-base depth by difference array | **no pseudogenome content at all** — the exit is hoisted above the rebuild | **0 pg rebuilds, 4 MB** |
| `query` | the assembly layer + placements, then stops | quality, identifiers, full read reconstruction | 1 pg rebuild, 42 MB |

`coverage` running in 9× less memory than the other two is the amortisation
claim made concrete rather than asserted: it needs only where the reads sit,
which the compressor stored, so it never materialises what they say.

### A heterozygous locus is not one place — and that is why `query` takes a sequence

`query` accepts a pseudogenome coordinate range **or a DNA sequence**, and the
second form is not a convenience. It is required by the same fragmentation the
caller has to undo. In a pseudogenome built to minimise bits, a het locus is
not one place but **N parallel places** — median 4, two haplotypes × two
strands — measured **up to 18.6 Mb apart**. A coordinate therefore names one
of them and returns a single haplotype, with the variation gone.

**This is not a missing feature, and no API can patch it.** The fragmentation
is not in the read-out layer — it is in the representation the compressor
chose, and it was chosen because it is the one that costs fewest bits. A
coordinate is a single integer into a single frame; the locus has no single
frame to be an integer into. Exposing a richer coordinate call, adding a
locus-to-offset table, or indexing the pseudogenome harder all answer the same
malformed question faster. The only fix is to address the archive by *content*,
which is a different question, and 81/400 against 345/400 is the size of the
difference between the two.

Measured on 400 GIAB het SNV sites across four individuals, both addressing
modes against the **same archive and the same sites** — an internal control,
because no other tool can produce a row of this table:

| addressing mode | both alleles returned | one allele | neither |
|---|---|---|---|
| by coordinate | **81 / 400** (20.2%) | 316 | 3 |
| by content (sequence) | **345 / 400** (86.2%) | 55 | 0 |

pooled allele balance 1.00 — the returned evidence is not skewed toward either
haplotype.

The probe is 40 bp of **reference** ending 6 bp *before* the variant, so it
never contains the variant: a probe taken from one haplotype's own sequence
could only match that haplotype and would rig the result.

**Why the control is internal.** SPRING addresses by read index; Genozip's
`--regions` is refused on FASTQ for want of coordinates; PgRC2 and NanoSpring
expose no read-out at all; BEETL-fastq returns reads *containing* a string,
which is a different object — measured, a median **38%** of the reads returned
at a locus here do not contain the probe (n=50). CRAM can answer a locus, but
only after aligning to an external reference, which is a different experiment.

**What this does not claim.** Not speed: `genocat --head=100` extracts in
0.19 s against 0.46 s here. The claim is that the question can be asked at all
of a reference-free archive, and — per the mechanism above — that it cannot be
asked by coordinate.

### The sidecar, and the limit it works around

An optional index (`capsule_decode index`) caches the pseudogenome, the
placements, and each read's deviations. It is **not part of the archive** and
its cost is not counted in any compression number; every operation works
without it. With it, `query` returns the true reads rather than the consensus
beneath them, and runs about 15× faster.

The deviations cannot be reached from the archive on demand, and the reason is
structural: `mm_sym` is coded with an **adaptive model in read order**, so read
*k*'s deviations require decoding all *k−1* before it. That is a property of
the coding — decompression only ever walks reads in order — and the sidecar
breaks the ordering dependency by decoding once.

**Cost of addressability in the archive itself:** the `contig_spans` stream is
232,509 B — **0.041%** of a 573 MB archive.

## 4. What is novel, and what is not

Stated in full, with the prior art actually read, in `NOVELTY_FINAL.md`.
In brief:

**Cited, not claimed.** Compressive genomics (Loh, Baym & Berger, *Nature
Biotechnology* 2012) established computing directly on compressed data — but
for genome/protein *databases* to accelerate *search*, not for read archives
and not for variant evidence. Assembly-based read compression is a family
(PgRC, NanoSpring, Minicom). Searchable compressed archives exist (BEETL-fastq
2014, the population BWT 2017). Haplotype fragmentation is documented in the
*assembly* literature (Purge Haplotigs 2018, Redundans 2016).

**Ours.** That the compression objective **deforms** the structure it builds:
compressing a heterozygous site optimally separates its alleles, which both
hides the variant from any caller and fragments the locus beyond what any
coordinate can name. Every tool in this family discards the assembly —
NanoSpring's own paper says it is "strictly a compression intermediate ... not
preserved or made available for downstream genomic analysis" — so the
deformation has not previously been observed, let alone corrected.

And the assembly field's remedy for fragmentation is to **purge** the redundant
haplotigs. For a lossless archive that is inadmissible: deleting a haplotig
deletes an allele, a variant, and losslessness.

## 5. Formats and interfaces

Container and stream layout: `FORMAT.md`. Command reference: `PIPELINE.md`.
Architecture contrasted with a de Bruijn caller: `ARCHITECTURE_VS_DISCOSNP.md`.
The compression-derived calling argument: `COMPRESSION_DERIVED_CALLING.md`.
