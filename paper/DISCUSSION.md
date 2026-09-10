# DISCUSSION — what this work claims, against the field that already exists

Condensed from the project's literature record: `docs/NOVELTY_FINAL.md`,
`docs/CLAIM3_PRIOR_ART.md`, `docs/HET_INDEL_SOTA.md`,
`docs/EBWT_LITERATURE_CORRECTION.md`, `docs/NOVELTY_SCAN_DBG_CHANNEL.md`,
`docs/ARCHITECTURE_VS_DISCOSNP.md`. Every source below was read as a primary
source, not summarised from an abstract.

---

## 1. The one-paragraph claim

> Compressing sequencing reads optimally requires assembling them, and the
> assembly encodes far more than the compressor stores: where every read sits,
> how each deviates from consensus, which sequences are alternative haplotypes.
> Every tool in this family discards that structure — NanoSpring's own paper
> says the assembly is *"strictly a compression intermediate ... not preserved
> or made available for downstream genomic analysis."* Worse, the compression
> objective actively **deforms** it: compressing a heterozygous site optimally
> separates its alleles, which hides the variant from any caller and fragments
> the locus beyond what any coordinate can name. We identify that deformation,
> measure it, and correct it — recovering reference-free variant calling
> (F1 0.431 → 0.888) and reference-free locus retrieval (81/400 → 345/400) from
> the same archive, at no cost in compression ratio and 0.041% of archive size,
> while remaining the smallest lossless FASTQ compressor measured.

---

## 2. What is NOT ours, stated first

A novelty claim is only worth as much as the list of things it concedes.

| piece | prior art | status |
|---|---|---|
| computing on compressed data | Loh, Baym & Berger, *Nat Biotech* 2012 (CaBLAST/CaBLAT) | **our paradigm.** Instantiated for genome/protein *databases* to accelerate *search* — not read archives, not variant evidence |
| **the pseudogenome architecture** | PgRC (Grabowski & Kowalski 2020), PgRC2 (2025) | **theirs.** We follow the architecture, not the construction — their sort-merge sweep was built and measured here: 11× fewer comparisons, 3.2× slower on 12 cores (`docs/REIMPL_NOTES.md` §39) |
| assembly-based read compression | Quip (Jones et al., *NAR* 2012) — the first; PgRC; NanoSpring; Minicom | **established.** We are in this lineage |

**Why SPRING, Genozip and PgRC2.** SPRING and Genozip are the only competitors
that, like this work, compress a complete FASTQ and reproduce the input file.
PgRC2 is included on the narrower DNA-stream basis, stated as such — it stores
no identifiers, no `+` line and no quality, verified by running it (440,048 B
FASTQ in, 151,000 B of bare DNA out).

Minicom is excluded as being in PgRC2's category: its README states it *"only
compresses DNA sequences ... does not support to compress the whole FASTQ
file"*, and PgRC2's own benchmark already reports 18–21% over SPRING/Minicom.
NanoSpring is excluded because it targets Oxford Nanopore long reads, which
this format refuses (`LIMITATIONS.md` §0); it is cited for its architectural
statement, not benchmarked.

## 3. The nearest prior art, and why it must be cited explicitly

> **BFQzip** — Guerrini, Louza & Rosone, *"Lossy Compressor preserving variant
> calling through Extended BWT"*, arXiv:2304.08534.

Reference-free FASTQ compression on EBWT + positional clustering, modifying
bases *and* quality scores together, explicitly to preserve variant-calling
information. **That is one shared structure serving compression and variant
calling, reference-free — the same *shape* as our central claim.**

The distinctions are real and in our favour, which is exactly why stating them
plainly costs nothing:

| | BFQzip | this work |
|---|---|---|
| compression | **LOSSY** — alters bases and quality | **LOSSLESS**, verified by decode-and-diff against the original FASTQ, 57/57 |
| relation to calling | *preserves* signal so an external caller does better | *performs* the calling, from the archive |
| output | a compressed file | a compressed file **+ a VCF + export/coverage/query** |
| structure | EBWT + positional clustering | pseudogenome + read placements + k-mer counts |

"Preserves variant-calling signal for someone else's caller" and "is itself the
caller" are different claims. But a reviewer who knows BFQzip and does not see
it cited will assume we did not know it existed.

---

## 4. What IS ours — three parts, in increasing order of strength

### 4.1 Capability (the weakest)

No tool retrieves the reads at a locus from a reference-free read archive.
SPRING addresses by read index; Genozip's `--regions` is *refused* on FASTQ
because the format carries no coordinates to index; BEETL returns reads
*containing* a string — a different object, and measurably so: **a median 38%
of the reads we return at a locus do not contain the probe** (n=50); CRAM
answers a locus but only after aligning to an external reference;
PgRC2/NanoSpring expose no read-out at all — `PgRC -h` offers compress,
decompress, and nothing else.

Being first at a capability is the weakest of the three arguments, and we say so.

### 4.2 Instantiation of an established paradigm where it was not attempted

Compressive genomics said *compute on the compressed form*. The read-compression
field built assemblies and **discarded them, by its own account**. We keep the
assembly and compute on it. That is a bridge between two literatures that had
not been connected — but it is an engineering-economics claim, not an
algorithmic one, and it should be written as such:

> Existing assembly-based compressors build k-mer structures and discard them
> after compression; existing reference-free callers build equivalent structures
> and discard them after calling. We show these are the same structure.

### 4.3 Mechanism — the load-bearing part

**The representation that compresses a heterozygous site best is the one that
conceals it.** Two internally-consistent contigs compress better than one contig
plus a column of disagreements, so optimal compression puts each allele on its
own contig; ref-allele and alt-allele reads are then never placed at the same
coordinate, and the variant is **not present in the data structure at all**.

We found no prior statement of this tension in the compression literature. Three
properties make it a mechanism rather than a fitted result:

1. **It predicts the error SHAPE, not just the error rate.** Without the
   correction, precision *holds* at 0.96 while recall falls to 0.27. The caller
   is not mistaken, it is **blind** — which is what "the alt reads are on
   another contig" predicts. Noise or a bad threshold would cost precision.
2. **The two corrective passes are synergistic, not additive.** +0.217 and
   −0.005 alone; **+0.457 together.** Re-placement alone is *worse than doing
   nothing*.
3. **It reappears one layer out.** A het locus is not one place in a
   compression-optimal pseudogenome — it is N parallel places (median 4, up to
   18.6 Mb apart), so the archive's own coordinate system cannot name a locus,
   and **adding a coordinate API would not fix it.** The 81/400 measures exactly
   that: the obvious engineering answer failing for the same structural reason.

**And the assembly field's remedy is inadmissible here.** Purge Haplotigs and
Redundans *purge* the redundant haplotigs to keep a pseudo-haploid reference.
For a lossless archive, deleting a haplotig deletes an allele, a variant, and
losslessness.

---

## 5. The negative results, and why they belong in the paper

Six substantial hypotheses were built, instrumented and refuted. They are not
padding: three of them are the ones a reviewer would otherwise propose.

| hypothesis | why it was compelling | what killed it |
|---|---|---|
| **read placements are free de Bruijn links** | McCortex pays ~20 GiB for exactly this; LueVari pays succinct colours for it. We had it already | precision **identical** across every split (0.916/0.952 by contig sharing; 0.948/0.944 by ≤10 kb span). `ppos` records where a read was **stored**, not where it came **from** — and chaining merges repeat copies *on purpose*, because that is what makes the archive small. The information was not left unexploited; it was **spent** |
| **turn the pseudogenome into a graph** | a node with exactly 2 successors is a het site; 3+ is a repeat. Textbook, and aimed at our measured weakness | **97.7% of anchors have out-degree 1.** Greedy chaining linearises the data on purpose; the redundancy removal that yields the archive is the same operation that removes the branches. FP *are* 2.9× enriched at branching nodes — and filtering on it still loses F1, because the TP pool is 7.6× larger. Enrichment ratios do not survive contact with base rates |
| **the placement pileup can replace the graph** | we have placements; DiscoSNP++ has to build a graph | measured on full chr20: pileup **F1 0.807 vs bubbles 0.877**, at 8× the caller time and 2.3× the RAM. Merging both channels was *worse than either* — recall reached the ceiling and precision collapsed 0.933 → 0.793 |
| **the compressor's own mismatch streams are a free variant channel** | `mem_extmm` holds 2.5 M (pos, ref, obs) triples the archive must carry anyway, and they rescue **72.3%** of the graph's misses | the evidence cannot isolate them. Best marginal precision over nine filters: **14.6% against a 44% bar.** Long clean anchors select for loci the graph already calls; the rescues live in short messy ones, indistinguishable from 1.09 M errors. Two of my own predictions — that recurrence is evidence, and that the rescues would be enriched for clustered variants — were both refuted by the same data |
| **locality-aware reference selection** | 44% of references point ~15 MB further than an identical nearer copy; the source stream is our biggest coding loss vs PgRC2 | implemented: max-length ties are **0.035%** of candidates. The tie-break fired on 100% of them and changed the archive by **+37 B**. The nearer copies are real but not in the candidate list — the lever is index sampling density, which is already the mapping bottleneck |
| **second-region self-match** | 74.3% of that region's 32-mers are repeat occurrences, and 2.2% were being removed | correct at 94.8% removal and the archive still **grew** 180 KB. LZMA was already exploiting the redundancy at 0.24 bits/base against 21.5 bits per explicit reference. *The measurement established that redundancy existed; it never established it was unexploited* |

The recurring lesson, stated once: **a free channel is not automatically worth
adding.** In the merge experiment 2,083 rescued truth sites arrived with ~7,500
extra false positives.

---

## 6. Two corrections against ourselves

Both are in the record because a result that has been corrected is more
trustworthy than one that never was.

**T3.4 was published at 0/400 and is 81/400.** The zero was an artifact of our
own query emitting the *consensus* rather than the reads, so no allele
difference could appear and coordinate addressing looked absolutely incapable.
The honest figure is a 4.3× gap, not an infinite one.

**Het-indel was withdrawn as a closed loss, and is now a 3-win.**
`docs/INDEL_BOUND.md` proved with the bound `F1 = 2·TP/(TP+FP+truth)` that no
filter could flip it at TP = 3,290 — and named the way out: *a different
candidate generator, not another filter.* The generator changed; TP rose to
3,871 and FP fell to 824, which is precisely what a filter cannot do. The bound
still holds and still reproduces the CSV.

---

## 7. Honest framing, claim by claim

- **Claim 1 is incremental and we say so.** A strong engineering result in a
  mature field, not a new idea. Its job is to earn the right to make Claims 2
  and 3 — without it, "our archive is also analysable" is a consolation prize.
- **Claim 2 claims competitiveness per class and unification overall**, not
  per-class dominance: 3 wins and 1 loss on SNV, 3 and 1 on indel.
- **Claim 3 does not lead with speed.** `genocat --head=100` beats our query
  (0.19 s vs 0.46 s). It leads with the mechanism.
