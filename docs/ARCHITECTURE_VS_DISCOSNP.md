# Architecture, layer by layer: G_CAPSUL vs DiscoSNP++

> **Status:** frozen 2026-09-10 at tag `v1.0.1-capsule` — code and results final.
> Authoritative numbers live in `../benchmark/results/`; verification status in
> [`../AUDIT.md`](../AUDIT.md). Where this file and a result file disagree, the result file wins.

Third column answers one question only: **is this layer textbook (free to any
implementer), theirs (a DiscoSNP++/GATB design decision we adopted), or ours?**
Written from `~/DiscoSnp/tools/kissnp2/src/*.cpp` and our own source.

## The table

| # | layer | DiscoSNP++ (kissnp2 / GATB) | G_CAPSUL (Method B) | provenance |
|---|---|---|---|---|
| 1 | k-mer counting | GATB: minimizer-partitioned superkmers, disk-partitioned | same shape, built **by the caller** (`kc_H_build`), not by the compressor — see the correction below | **textbook** (KMC2/GATB) |
| 2 | minimizer ordering | frequency-ranked (partition balance) | lexicographic; frequency ranking implemented, **measured +14% volume, rejected** | **theirs**, and measured not to pay here |
| 3 | memory budget | `-max-memory`, declared | ceiling = 60% of measured MemAvailable, spill self-enables | **textbook** systems practice |
| 4 | graph representation | cascading-Bloom dBG, GATB | flat hash table (`kc`) + `kc_find` as membership oracle | **textbook** |
| 5 | bubble start | branching node, >=2 successors, both strands | same, both orientations | **textbook** |
| 6 | SNV closure | lockstep `expand`, close on `nextNode1==nextNode2` | lockstep pairwise walk | **textbook** — same algorithm |
| 7 | multi-polymorphism | `-P 3` | `MAXPOLY=1` (structural: k=31 vs ~1/1000 het) | **theirs**, ours narrower and we still win |
| 8 | branching policy | `checkBranching`, `b=0` | `branch_walk` stops at first junction — equivalent | **textbook**, arrived at independently |
| 9 | indel proposal | `start_indel_prediction`: BFS per successor pair, **unconditional** | superbubble (Onodera) only — fires at **15.9%** of nodes | **theirs**; our gap, ported opt-in |
| 10 | indel closure shape | guaranteed by lockstep traversal | `LCP+LCS >= |s_short|` test, added after tracing theirs | **theirs** (geometry), **ours** (the explicit test) |
| 11 | indel ambiguity | `checkRepeatSize`, reject `k-2-min(ext) > 20` | ported; **monotone loss on our data**, kept at their default | **theirs**, measured |
| 12 | low complexity | DUST; **off by default** (`l="-l"`) | absent | **theirs**, not a gap |
| 13 | read coherence | **kissreads2 — a separate tool that re-maps every read onto every bubble** | reads held in memory (a deliberate second FASTQ pass, cheap I/O — *not* left over from compression); one indexed sweep + 1-bit quality bitmap | **ours** |
| 14 | quality use | mean phred per path | per-base bitmap (`MINQ=20`), 233 MB vs 2.27 GB of phred strings | **ours** |
| 15 | parallelism | GATB thread pool throughout | traversal was serial under superbubble mode; now parallel, **2.9x**, counts identical | **textbook**, was our deficit |
| 16 | coverage ceiling | none equivalent | `COVCAP = 2*PLOIDY*H`, derived from measured depth | **ours** |
| 17 | ploidy handling | none | HETSCAN pair-fraction gate, declines on haploids (E. coli 0.016) | **ours** |
| 18 | output | FASTA bubbles -> VCF_creator | VCF direct, plus a lossless archive from the same pass | **ours** |

## Why we are different from them, not just similar

**1. The ASSEMBLY is the byproduct — not the graph. (Corrected 2026-09-05.)**
An earlier version of this file claimed the k-mer table was "already built by
the compressor". That is false: `kc_H_build` is inside `run_variant_call`
(`include/caps_caller.h:1812`) and `stages/106_inprocess.cpp` builds no k-mer
table at all. We pay for the graph like everyone else.

What IS genuinely free is the **assembly**. The encoder builds contigs in order
to compress, and Method B reuses them directly — `[DBG-ONLY] substrate skipped
(451,578 encoder contigs reused)`. The alternative, `build_substrate`, re-places
all 12.6M reads and measured **738 s serial at full chr20**; skipping it is the
difference between a ~57 s and a ~150 s Method B run. That is a measured,
retained-assembly saving, and it is the project's actual thesis (Assemble →
Retain → Compress → Serve), not a graph claim.

**2. Read coherence needs no second tool.** kissreads2 is a separate binary
that re-maps every read onto every bubble. We keep the reads in memory and do
one indexed sweep instead, which also makes a *per-base* quality test affordable
where they use a per-path mean (layers 13, 14). Stated precisely: the reads are
loaded by a deliberate second FASTQ pass, NOT left resident by compression — the
encoder re-reads them specifically so the no-CAPS_CALL memory footprint is
unaffected. The saving is the separate mapping tool, not the I/O.

**3. Our parameters are derived from measured depth; theirs are constants.**
`COHC = max(2, H/10)` and `COVCAP = 2*PLOIDY*H` are functions of the sample's
own haploid depth and ploidy. That is why precision holds 0.911-0.950 across a
10-30x sweep and why the tetraploid arm works without retuning. DiscoSNP++
exposes `-b`, `-P`, `-D`, `-max_ambigous_indel` as fixed numbers.

**4. We win where the shared algorithm is shared, and lose where theirs is
richer.** Layers 5-8 are the same textbook bubble calling, and on SNVs we beat
them (0.879 vs 0.847) — that is the fair comparison, and it isolates our
contribution to layers 13-17. On indels their layer 9 attempts a call at every
branching node while ours attempts one at 15.9% of them; that is a real
architectural deficit, honestly a place where their design is better, and it is
why the indel claim is withdrawn rather than argued.

## Architectural headroom — the leading candidate, tested and refuted

**"Free read threading" does not exist. Built, measured twice, withdrawn.**
Full detail in `PLACEMENTS_AS_LINKS_REFUTED.md`.

The hypothesis was that our `ppos`/`read_cid` placements are the read-to-path
threading McCortex pays a dedicated pass and ~20 GiB of links for
([Turner et al. 2018](https://academic.oup.com/bioinformatics/article/34/15/2556/4938484)),
and that LueVari pays succinct structures for as read colours
([Bioinformatics 2020](https://academic.oup.com/bioinformatics/article/36/22-23/5275/5734643)).
The placements are real and were being discarded — but they carry no signal:

| split | TP | FP | precision |
|---|---|---|---|
| alleles share a contig | 76 | 7 | 0.916 |
| alleles disjoint | 260 | 13 | 0.952 |
| pseudogenome span <= 10 kb | 73 | 4 | 0.948 |
| span > 10 kb | 267 | 16 | 0.944 |
| baseline | 336 | 20 | 0.944 |

Precision is identical across every split. The reason is structural: a
pseudogenome is built by greedy overlap chaining, whose coordinates are
optimised for compressibility, and chaining **deliberately merges** near
identical sequence — repeat copies are the most compressible thing in a genome,
so they collapse first. `ppos` records where a read was **stored**, not where it
came **from**. The repeat information was not left unexploited; it was **spent**
to compress. Recovering it would mean not collapsing repeat copies, which is
directly opposed to what makes Claim 1 win.

Note also that the underlying idea was never novel: LueVari already does
reference-free SNP calling on a read-coloured graph for exactly this purpose.

**What remains genuinely free, and is measured:** the assembly (point 1 above) —
encoder contigs reused, `build_substrate` skipped, 738 s saved at full chr20.
That is the retained-assembly thesis, and unlike the threading idea it survives
measurement.

**Prior art to be careful about.** "Searchable compressed archive of reads" is
not new — BEETL-fastq
([Cox et al., *Bioinformatics* 2014](https://academic.oup.com/bioinformatics/article/30/19/2796/2422232))
already does selective read extraction and variant-related queries from a
compressed archive. Claim 3 should be positioned against it explicitly rather
than claiming addressability as novel. The defensible framing is the
*combination*: lossless compression, a retained assembly, and reference-free
calling from one pass — not any one of those alone.

**Two smaller levers, both measured as open rather than assumed:**
- Calling from a *stored* archive (today it is compress-time only) — this is
  what would make "addressable" a demonstrated property rather than a design
  intent.
- Partition balance from layer 2's frequency ranking is still unmeasured; it was
  rejected on volume, which is not what it is for. If the merge is ever made
  memory-bound per partition, re-measure it there.

Sources:
- [Integrating long-range connectivity information into de Bruijn graphs (McCortex)](https://academic.oup.com/bioinformatics/article/34/15/2556/4938484)
- [BEETL-fastq: a searchable compressed archive for DNA reads](https://academic.oup.com/bioinformatics/article/30/19/2796/2422232)
