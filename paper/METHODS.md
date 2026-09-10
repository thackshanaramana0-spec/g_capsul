# METHODS — the validated architecture, as implemented

**Consolidated 2026-09-10** from the final code, the final benchmark harness and
the executed results. Where an older document disagreed with the code, the code
won; where the code disagreed with an executed result, the result won.

Results and their provenance live in `../benchmark/`, not here. This file
describes *what the system does*; `../benchmark/documentation/RESULT_CODE.md`
proves *what it produced*.

---

## 1. Compression — the pipeline

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

The caller reads the archive's assembly layer and placements and adds two
passes, **to the caller only; the compression path is untouched**:

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

    neither             F1 0.431   P 0.967   R 0.278
    re-placement only      0.426     0.962     0.274   <- worse than neither
    collapse only          0.648     0.963     0.488
    both                   0.888     0.956     0.830

Synergy, not additivity: +0.217 and −0.005 alone, **+0.457 together**. And the
error *shape* confirms the cause rather than merely fitting it — without
collapse, precision holds at 0.96 while recall collapses to 0.27. The caller is
not mistaken, it is **blind**, which is what "the alt reads are on another
contig" predicts. Noise or a bad threshold would cost precision instead.

**The correction costs nothing in compression ratio.**

## 3. Addressability — export, coverage, query

Served directly from the archive; the work a conventional pipeline does at
query time (assembling, or aligning and indexing) was already done at
compression time.

| operation | what it decodes | what it does NOT |
|---|---|---|
| `export` | `literal` + `mem_triples` (+ companions) → the pseudogenome | no per-read stream |
| `coverage` | `pos_abs` + `read_lengths` → per-base depth via a difference array | **no pseudogenome content at all** — hoisted above the rebuild |
| `query` | the assembly layer + placements, then stops | never touches quality or names, never reconstructs all reads |

`query` accepts either a pseudogenome coordinate range or a **DNA sequence**.
The sequence form matters because of the same fragmentation the caller has to
undo: a heterozygous locus is not one place in the pseudogenome but N parallel
places (median 4 — two haplotypes × two strands), so a coordinate names one of
them and returns a single haplotype. Content addressing resolves every parallel
representative at once.

An optional sidecar index (`capsule_decode index`) caches the pseudogenome,
the placements and each read's deviations. It is **not part of the archive**;
without it every operation still works. With it, `query` returns the true reads
rather than the consensus, and runs about 15× faster.

**Cost of addressability:** the `contig_spans` stream is 232,509 B — **0.041%**
of a 573 MB archive.

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
