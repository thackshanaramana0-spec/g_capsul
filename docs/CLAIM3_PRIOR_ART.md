# Claim 3 — prior art survey: who else offers these operations?

Literature check 2026-09-03, done before locking Claim 3, to establish what is
genuinely novel and what is not. **Random access to compressed FASTQ is NOT
novel** and the paper must not claim it is.

## 1. Prior art that DOES exist

| tool / work | what it offers | access key |
|---|---|---|
| **bgzip / razip** | block-compressed gzip with random access | byte offset |
| **BEETL-fastq** (Bioinformatics 2014) | searchable compressed archive; indexes every 1024th read | **k-mer search / read id** |
| **CIndex** (Bioinformatics 2022) | compressed indexes for fast FASTQ retrieval | **record / id** |
| **sFASTQ / SFQ** (Electronics 2022) | succinct representation, fully random access to individual records, searchable on disk | **record id** |
| **GPU LZ77 work** (arXiv 2026) | position-invariant random access, any region in 0.4 ms, blocks self-contained | **byte/block range** |
| **CRAM / BAM** | genuine region queries | **genomic coordinate — but requires a REFERENCE and prior alignment** |

So: record-level and byte-level random access on FASTQ is well established, and
coordinate-level access exists for *aligned* formats.

**None of the above have the three specific operations.** Checked directly,
not inferred from "has some access mechanism":

| tool | export (assembly output) | coverage (per-base depth) | query (coordinate range) |
|---|---|---|---|
| BEETL-fastq | no | no | no — retrieves by read id / k-mer, not a coordinate range |
| CIndex | no | no | no — k-mer to containing-reads lookup only |
| sFASTQ | no | no | no — record-id random access, not coordinate range |
| GPU LZ77 (2026) | no | no | no — byte/block-range decode, not a genomic coordinate |
| CRAM/BAM + samtools/mosdepth | no (no assembly step — it aligns to a supplied reference) | **yes** | **yes** |
| PgRC / Minicom / NanoSpring | no — build an internal pseudogenome/contigs but never expose it | no | no |

The assembly-based compressors (PgRC, Minicom, NanoSpring) are the closest
architectural relatives — they build the same kind of pseudogenome/contig set
CAPSULE does — but per their papers and released tools, none expose it as a
user-facing operation. That is the strongest point for novelty, not a weak
one: the capability was structurally available in that whole family and
nobody surfaced it. CRAM/BAM is the only prior art with genuine
coverage+query, and it requires a reference genome and a prior alignment step
that CAPSULE does not have.

## 2. What is actually different here

Every FASTQ tool above is keyed on a **record identifier, byte offset, or k-mer
string** — because a FASTQ file has no coordinates to key on. CRAM/BAM do offer
coordinate access, but only after aligning to a reference genome.

CAPSULE's three operations are keyed on **coordinates in an assembled
pseudogenome that the compressor built at compress time, with no reference
genome anywhere in the pipeline**:

* `query <START-END>` — reads overlapping a **coordinate range**, reference-free.
* `export` — the assembly itself, already built, no assembler run
  (E. coli: 0.465 s vs MEGAHIT 135.5 s = **291×**; vs the spec-exact
  baseline SPAdes 258.15 s = **555×**, peak RAM 5.20 GB for SPAdes).
* `coverage` — per-base depth with **no alignment and no index**
  (25–33× vs bwa+samtools+mosdepth).

The defensible novelty is therefore **not** "random access to a compressed
FASTQ". It is:

> **reference-free coordinate addressability** — plus assembly and depth
> obtained as archive reads rather than as computations, because the compressor
> already had to compute them.

`export` and `coverage` in particular have no counterpart in any FASTQ
compressor found: BEETL-fastq, CIndex and sFASTQ return *records*, not an
assembly or a depth profile.

## 3. Honest wording for the paper

**Do not write:** "the first compressed FASTQ archive supporting random access."
That is false — BEETL-fastq did it in 2014.

**Do write:** something of the form "an archive addressable by assembled-genome
coordinate without any reference, from which the assembly and per-base depth are
recovered by decoding rather than recomputation", and cite BEETL-fastq, CIndex
and sFASTQ as the record-level prior art, with CRAM/BAM as the reference-based
coordinate prior art.

## 4. Limitation to state alongside it

`query` is **O(archive), not O(range)**: it must rebuild the pseudogenome
because reads are stored as slices of it, so a 1 kb query costs the same as a
1 Mb one. Measured on E. coli: query 0.530 s vs full decompression 1.74 s —
a 3.3× saving in time, but the output is ~1/280th the size (11,103 reads
instead of a 695 MB FASTQ). The advantage is **selectivity, not asymptotics**,
and the GPU work above genuinely beats us on that axis (0.4 ms, block-local
decode). State this rather than implying constant-time region access.
