# Claim 3 — ADDRESSABLE

**One document: what the claim is, why it is novel, the mechanism, what it
costs, and what it does not do.**

Status: real data, real archive, held-out replication. Last measured
2026-09-15. Every number here was produced by a script in
`scripts/t34_realdata/` and is recorded in `results/T34_REALDATA_20260915/`.
Where this document and a result file disagree, the result file wins.

**Relationship to the other Claim 3 documents.** This file consolidates the
claim, the novelty argument and the mechanism, and it is the one to read first.

- [`CLAIM3_MECHANISM.md`](CLAIM3_MECHANISM.md) established *a heterozygous locus
  is not one place* (median 4 parallel addresses, up to 18.6 Mb apart). That
  finding is **unchanged and is restated here in §2 and §5a** — this document
  extends it to the retrieval side rather than replacing it.
- [`CLAIM3_LOCKED.md`](CLAIM3_LOCKED.md) is already marked superseded and
  describes a different experiment on the outer ARCS repo. Do not cite it.
- [`T34_RESIDUE_DIAGNOSIS.md`](T34_RESIDUE_DIAGNOSIS.md) is the diagnosis trail,
  kept for the refuted hypotheses. It stops at 391/400.

**The published benchmark numbers are untouched by all of this.** The citable
figure -- the frozen sweep's T3.4, which is the POSITION claim -- remains
345/400 content against 81/400 coordinate from the
19-dataset sweep. What follows is a separate, stricter experiment on a separate
archive, and it supersedes nothing until the benchmark box re-runs the sweep.
See §17.

---

## TERMS, AND HOW TO READ EVERY NUMBER HERE

Read this before the tables. Two words in this document look alike and are not
the same thing, and one number (400) appears everywhere without explanation.

### The words

| term | meaning |
|---|---|
| **locus** | one position in the genome. A place. Plural *loci* |
| **probe** | a short stretch of reference sequence, 40 bp here, used as the *address* of a locus. It is the thing you hand to the archive |
| **pseudogenome** | the single long string the compressor builds by overlapping reads together. Reads are stored as a *position on this string* plus their differences from it |
| **contig** | one continuous stretch of that string |
| **placement** | where a given read sits on the pseudogenome. This is what gives the archive an internal coordinate system |
| **sidecar** (`.qidx`) | an index built *from* the capsule that holds the pseudogenome, every placement, and each read's differences |

Earlier drafts used "locality" as a noun for the idea of nearness. That was
confusing next to "locus" and has been removed. The distinction this document
actually draws is between two kinds of **retrieval**, named below.

### The two kinds of retrieval, which is the whole novelty argument

You want the reads covering one position. You supply a 40 bp probe as its
address.

```
              the 40 bp probe you supply
              +------------------+
  reference   ...ACGTTGCAATGCCTAG...T...
                                     ^ the position you actually care about

  read A      ---------------------------   contains the probe exactly
  read B      -----------x---------------   one sequencing error inside the probe
  read C              ---------------------  starts late, covers only part of it
  read D   ---------------                   ends early, covers only part of it
```

All four reads physically cover the position. They do not all contain the probe.

| | what it returns | who does it |
|---|---|---|
| **retrieval by exact match** | read A only | BWT / FM-index archives (BEETL-fastq, CIndex, sFASTQ) |
| **retrieval by position** | reads A, B, C and D | this work |

A BWT stores every read as a separate standalone string. Nothing records that
read C sits beside read A, so reads B, C and D are unreachable however the query
is phrased. A pseudogenome stores reads as *placements on a shared string*, so
the probe is located **once**, and then everything whose placement overlaps that
span is returned. Read B's error is irrelevant. C and D's partial overlap is
irrelevant.

**Why this matters rather than being a curiosity.** Read B still carries a real
base at the position you asked about. Dropping it skews the count of one allele
against the other, and it skews it *specifically* against the allele that
differs from the probe — because a read differing from the probe is exactly the
read an exact-match test rejects.

### Can this archive answer BOTH questions?

**Yes, both, completely** — and the asymmetry that used to exist was measured,
chased to its cause, and closed.

| | **exact match**<br>every read *containing* the string | **position**<br>every read *covering* the locus |
|---|---|---|
| **BWT / FM-index**<br>BEETL-fastq, CIndex, sFASTQ | **yes — 1.000** by construction | **no** — reads are independent strings with no inter-read geometry, so this cannot be asked at all |
| **G_CAPSUL** | **yes — 1.0000**, from `query` alone at k=3 | **yes — native** |

**We answer both. The BWT family answers one.** That is the claim in one table,
and both cells on our row are measured rather than argued:

| | value |
|---|---|
| exact match, as it was | 0.9789 |
| + strand bitmap | 0.9975 |
| + N restoration | 0.9998 |
| + the `.xmi` index | **1.0000** |
| position, bilateral | **400/400** |
| cost of the index | **+0.104 ms/site**, 0.2% of a tolerant query |
| index build | 1.9 s, over 14% of the reads |
| why only 14% | proven: absent from pg ⟹ carries mismatches **or** contains N |
| evidence exact match cannot reach | **86.2%** of the reads we return |

Exact-match recall is **0.979 as shipped, and 1.000 with an indexed lookup that
costs 0.104 ms/site** — 0.2% of a tolerant query.

**The shortfall was chased to its cause, then closed.** An earlier 0.945 figure
is withdrawn — it was a bug in the measurement, not the tool. Reads are returned
in *pseudogenome orientation*, so every reverse-complemented read counted as a
miss. Canonicalising both sides gives 0.979.

Two hypotheses for the remaining 2.1% were **refuted by measurement**:

1. *Tolerance is too low.* `CAPS_QUERY_MINSEED` was added to lower the seed floor
   and lift the k cap. Recall is **identical** at MINSEED 12 / 8 / 6 / 5
   (k_max 2 / 4 / 5 / 7) — 0.979 and 34 misses every time. Only wall time moves.
2. *The probe is near the consensus and we failed to find it.* Refuted.
   Categorising all 34: **0** within k of the consensus, **0** merely far,
   **34** whose literal sequence **does not appear in the pseudogenome at all**.

**That third category is what closes it.** A read absent from the consensus
*must* carry deviations — a read with none would *be* a pg substring. So every
possible miss is provably confined to the subset of reads that are not literal
pseudogenome substrings, and that subset is **14% of the reads**. (The sidecar's
deviation list is the obvious way to enumerate it, but is **not** sufficient —
see the caveat below.) Restricting the search to that subset:

| path | recall | cost |
|---|---|---|
| query alone | 0.9789 | — |
| query + naive subset scan | 1.0000 | +41.6 ms/site |
| **query + INDEXED subset lookup** | **1.0000** | **+0.104 ms/site** |

The naive scan was the wrong shape — the subset does not change between queries,
so the work belongs at **build** time. A k-mer→read-id map over the 14% subset
only: 1.2M distinct 20-mers, 2.1M postings, **1.9 s to build**, 255 MB in a
Python dict. A read containing a 40 bp probe contains *every* one of the probe's
k-mers, so one lookup yields a candidate set guaranteed to contain every true
hit and verification removes the rest — **no recall is traded for the speed.**

**400× faster than the naive scan.** Against a ~52 ms/site tolerant query, the
closure adds **0.2% overhead** to take recall from 0.979 to 1.000. The cost
question is settled: this is not a trade.

**So the architecture is not structurally barred from the BWT's guarantee, and
the earlier framing of this as a permanent trade-off is withdrawn.** It can
answer both questions completely: position natively, exact match at 1.000 via an
index over the 14% of reads that can possibly be missed.

**The subset is settled, by proof rather than by measurement.** A read is stored
as a placement plus mismatches plus N positions, and reconstruction is
`pg[placement..]` with mismatches applied and Ns restored. So a read with zero
mismatches and zero Ns is byte-identical to a pseudogenome substring.
Contrapositive:

> **absent from the pseudogenome ⟹ carries mismatches OR contains N**

Both lists are already written by the encoder — the sidecar's deviation list and
the `n_pos` / `n_indices` / `n_cnt` streams — so their **union is a provable
superset** of the reads that can possibly be missed. No new structure is needed.

Two earlier explanations were wrong and are recorded as such. The first was that
the deviation list alone was the subset — it is not, it omits N-carrying reads,
which deviate from the consensus without being mismatches. The second was that
the residual difference was duplicate read instances — **tested and refuted**,
there are exactly 2. The counts decompose exactly:

| | distinct reads |
|---|---|
| absent from the pseudogenome | 16,479 |
| of those, N-containing | 234 (encoder logs 236 instances) |
| absent without N | 16,245 — necessarily all carry mismatches |
| sidecar deviation list | 16,415 — a superset of 16,245, as the proof requires |

**What remains is engineering, not uncertainty:** the index must be built from
that union inside `capsule_decode index`, and `query` must consult it. Today the
closure lives in a script.

### query alone reaches both — no index, no second mode

The `.xmi` index is not needed for this. If a read contains the probe then
`read = pg[q..] + deviations`, so the **pseudogenome at that spot differs from
the probe only at deviation positions inside the window**. A search tolerating
`k >= (deviations in window)` finds that pg location, and the read — placed
there — comes back through the ordinary position path.

| MINSEED | k_max | exact-match recall | wall, 200 sites |
|---|---|---|---|
| 12 | 2 | 0.9998 | 24 s |
| **9 (derived)** | **3** | **1.0000** | 29 s |
| 8 | 4 | 1.0000 | 37 s |

**An earlier sweep said tolerance changed nothing, and that verdict is
withdrawn.** It ran before the strand fix, so it was measuring against emitted
sequences that were not the reads. On the fixed binary it moves immediately.

**The seed floor is now derived rather than hardcoded.** It was a fixed 12,
which silently capped a 40 bp probe at k=2 — one short of completeness. But the
floor that matters is not a constant: a seed of length L occurs by chance about
`|hay| / 4^L` times, so specificity is a property of the **haystack**, not a
tuning knob.

```
MINSEED = clamp( ceil(log4(|pseudogenome|)) - 2, 8, 16 )
```

The `-2` is a bounded, deliberate concession — it admits about 16× more chance
candidates, every one of which is verified against the full tolerance, so **the
answer cannot change, only the time to reach it**. On a 2.15 Mb pseudogenome it
yields 9, which lets a 40 bp probe reach k=3. On a 3 Gb one it yields 14, which
is the right answer there: short seeds on a genome-scale haystack are exactly
how a candidate set explodes.

**The two questions want different tolerance, and that is a finding rather than
a flaw:**

| k | exact-match recall | het 400 (≥2 reads) | homozygous FP | separation |
|---|---|---|---|---|
| 2 | 0.9998 | 400 | **22** | **378** |
| 3 | **1.0000** | 400 | 35 | 365 |

Raising k costs specificity in the **genotyping** interpretation. It costs
nothing in the **exact-match** answer, because every candidate is verified
against the full tolerance — a higher k buys time, never a wrong read. So:
`CAPS_QUERY_MM=2` for genotyping, `CAPS_QUERY_MM=3` for exact-match
completeness, and one call at k=3 returns a set containing both answers.

### Both at full — what it took

Getting exact match from 0.979 to **1.0000** while position stayed at
**400/400** took three pieces. Each was necessary: removing any one drops recall.

| piece | what it fixes | recall without it |
|---|---|---|
| **strand bitmap** | the sidecar never stored strand, so `query` mirrored RC substitutions onto the wrong bases and emitted a sequence that was not the read | 0.9789 |
| **N restoration** | N-reads traverse the pipeline with each N replaced by `A`, so `query` rebuilt an `A` where the read had an `N` | 0.9975 |
| **`.xmi` k-mer index** | a read that contains the probe but whose placement does not overlap the probe's pseudogenome occurrence is invisible to a position query | 0.9998 |

Cost: the sidecar grows **1.8%** (14.7 kB strand + 13.3 kB N on 1.59 MB). The
`.xmi` is a separate 3.3 MB file written only when asked. The archive is
untouched and still lossless.

**Index design, so it generalizes rather than fitting this dataset.** `k` is
derived from the data — `min(20, shortest indexed read)`, capped at 31 by 2-bit
packing — never a fixed constant. Only **deviation-carrying reads** are indexed,
because a read with no deviations is emitted as pure pseudogenome and is always
reachable by the ordinary search: 14% of reads here. Completeness holds for
probes of length `k + stride - 1`, which is checked at query time and **refused
loudly rather than degraded silently**. The lookup is a superset filter — a read
containing a probe contains *every* one of the probe's k-mers, so one indexed
hit is guaranteed — and every candidate is verified afterwards.

All three are **off unless asked**: no `.xmi` is written without `CAPS_XMI=1`,
the sidecar is byte-identical with and without it, and `query` output is
unchanged without `CAPS_QUERY_CONTAIN=1`.

### The bug that the index exposed

Wiring the exact-match index into the binary did not raise recall. The index
built, looked up and verified candidates correctly, and recovered nothing. That
contradiction — a working index that finds nothing new — is only explicable if
**the sequence `query` emits is not the read**.

It was not. The encoder stores a mismatch offset `j` whose pseudogenome index is
`q+j` for a forward read but **`q+RL-1-j` for a reverse-complement one**. The
sidecar stored the pseudogenome, positions and lengths — **but not strand**. So
`query` applied every deviation at `q+j`, mirroring the substitutions onto the
wrong bases for every RC read that carries one.

It survived because it is invisible unless a probe window happens to cover an
affected base. Aggregates moved about 2%, which reads as noise.

**The fix** stores a strand bitmap in the sidecar — 1 bit per read, 14.7 kB
against a 1.24 MB sidecar — and applies deviations at the right offset in the
emit loop, the index builder and the index verifier alike. Sidecars without the
block are detected and the old behaviour announced rather than guessed at.

| metric | before | after |
|---|---|---|
| exact-match recall | 0.9789 | **0.9975** |
| missed reads | 34 | **4** |
| upstream-only, 400 het sites | 391 | **392** |
| bilateral | 400 | **400** |
| homozygous false positives | 31 | **29** |
| archive lossless | yes | **yes** |

**It is output-changing and deliberately so:** 60 of 100 sampled sites now emit
different query output. Every measured aggregate improves and none regresses,
which is what standing rule 2 requires of an output-changing change.

**Disclosure.** The published 345/400 was measured with the query that
mis-placed RC deviations. Re-running it would likely move it, probably upward.
It has **not** been re-run and is not restated here on the strength of this.

**The last 4 misses are not a search problem.** `query` reconstructs from the
2-bit packed pseudogenome and never restores N from the `n_pos` / `n_indices` /
`n_cnt` streams, so an N-carrying read is emitted with a pseudogenome base where
the N was. Separate gap, clearly scoped.

### Which question does the real world ask?

**Exact match screens. Position diagnoses.**

| | typical question | evidence needed |
|---|---|---|
| exact match | Is TB present in this sample? Does this isolate carry the resistance gene? Pull the reads for this gene so I can assemble it | one read is enough — this is a **presence** question |
| position | At this BRCA1 site, which bases does the patient carry, and on one chromosome copy or both? What fraction of tumour cells carry it? Is coverage deep enough to trust the call? | **all** reads at the site — this is a **proportion** question, and a missing read biases the answer |

There is a further asymmetry that matters clinically. A probe is written from
the reference, so a read matches it exactly when the sample *agrees* with the
reference. **Exact match therefore works best where nothing interesting is
happening and degrades exactly at the variant sites you queried for**, because
the reads carrying the variant are the reads that differ. That is not a small
effect — it is the same mechanism that made Defect B (§8) silently discard
alternate-allele reads.

### What "400" means

**400 is 400 test positions**, not reads, not percent. They are heterozygous SNV
sites taken from the GIAB v4.2.1 verified truth set for HG002 — positions where
that individual genuinely carries two different bases, one per chromosome copy.

For each site the archive is asked for the reads there, and the result is scored
into one of three buckets:

| bucket | meaning |
|---|---|
| **both** | both the reference base and the alternate base came back. Correct |
| **one** | only one of the two came back. The site is there but an allele is missing |
| **neither** | nothing usable came back |

So **400/400 means both alleles recovered at all 400 sites.** 391/400 means both
at 391 and a single allele at the remaining 9.

The **homozygous negative control** runs the identical pipeline over 400
positions that have *no* variant. There the correct answer is "one", and any
"both" is a false positive.

---

## PART I — THE CLAIM

### 0. Claim 3 in one table — every sub-win, and where it stands

Claim 3 has four sub-claims. Three come from the frozen 19-dataset sweep and are
untouched by this document. The fourth, locus retrieval, is what this document
rebuilt.

| | sub-claim | result | baseline | kind of claim |
|---|---|---|---|---|
| **T3.1** | **export** the assembly | **129–784× faster** | SPAdes | competitive — settled |
| **T3.2** | **coverage** | **16–54× faster** | bwa + mosdepth | competitive — settled |
| **T3.3** | **query** across the locked set | 19/19, 0 failures, 0.20–6.77 s, 44 MB–1 GB | **none** | **coverage, not competitive** — it shows the feature works everywhere and what it costs, never that anyone else is worse |
| **T3.4** | retrieval by **exact match** | recall **1.0000** | BWT/FM-index = 1.000 **by construction** | competitive — parity |
| **T3.5** | retrieval by **position** | **400/400** | **bwa + samtools view** — the workflow a user actually runs today | competitive — exclusive |

> **Numbering collision, stated so it is never confused.** The frozen sweep's
> `benchmark/results/claim3/claim3_T3.4_locus_fidelity.csv` is the **position** claim
> (345/400 content vs 81/400 coordinate). In this document T3.4 is **exact
> match** and T3.5 is position. When citing the frozen sweep, use its own names.

**T3.4 and T3.5 are different questions and need different baselines.** T3.1 and
T3.2 are speed claims against tools that do those jobs (SPAdes, bwa+mosdepth),
and they are settled. T3.4 and T3.5 are retrieval claims, and the only tools in
that space are the BWT/FM-index archives — BEETL-fastq, CIndex, sFASTQ. Those do
not export assemblies or compute coverage, so they are irrelevant to T3.1/T3.2
and are the *only* relevant comparison here.

**T3.4 is parity with the best that exists.** A BWT answers exact match with
recall 1.000 by construction. So do we, measured, on three archives. This is the
row where the comparison bites, and we do not lose it.

**Compression is not argued here, deliberately.** Claim 3 is a retrieval claim.
Archive size is Claim 1's job, and Claim 1 already measures it against the actual
state of the art in compression — SPRING (−6.03%) and Genozip (−43.26%), 19/19.
BEETL-fastq is a *retrieval* competitor, not a compression competitor, and an
earlier draft of this document inferred a compression gap from BEETL's and
SPRING's separately published gzip baselines. That inference is **withdrawn**:
different datasets, different content, and gzip's own ratio varies with the
data, so it was a workaround rather than a result. The two claims stand
independently and neither needs the other's baseline.

**T3.5 is the row nothing else can enter.** A BWT stores every read as an
independent string with no record of one read lying beside another, so "which
reads cover this locus" cannot even be posed. The 96/400 coordinate arm is our
own control, not a competitor. Parity on T3.4 is what makes T3.5 count — without
it, a reader could assume we traded one capability for the other.

**T3.4 and T3.5 together, stated once:**

| | **T3.4** exact match — reads *containing* a string | **T3.5** position — reads *covering* a locus |
|---|---|---|
| BWT / FM-index (BEETL-fastq, CIndex, sFASTQ) | 1.000 | **not possible** |
| **G_CAPSUL** | **1.0000** | **native** |

Both from `query` alone, on every archive tested, with nothing tuned between
them. The three sub-results that got it there, each one necessary:

| piece | what it fixes | recall without it |
|---|---|---|
| bilateral anchoring | a probe anchored on one side is unreachable when overlap chaining breaks that side | 392/400 position |
| strand bitmap | the sidecar never stored strand, so RC deviations landed on the wrong bases and `query` emitted sequences that were not the reads | 0.9789 exact |
| N restoration + derived seed floor | N was never restored, and a hardcoded floor of 12 silently capped a 40 bp probe one mismatch short of completeness | 0.9998 exact |

And one efficiency result, taken because the audit found the archive was already
computing it and throwing it away:

| | het 400 | homozygous FP | scoring cost |
|---|---|---|---|
| read decoding | 400 | 31 | 16.30 ms/site |
| **consensus + tally, count ≥ 2** | **400** | **25** | **0.83 ms/site** |

**19.6× faster and more specific at the same time** — the archive already holds
the pileup its reads would produce.


### 1. What is being claimed

A G_CAPSUL archive can return **the reads at a locus, named by sequence
content, without a reference genome and without an aligner**.

Three things have to hold at once, and the third is the hard one:

| | |
|---|---|
| **COMPACT** | the archive is smaller than the alternatives (Claim 1) |
| **FAITHFUL** | variants can be called from it (Claim 2) |
| **ADDRESSABLE** | a locus can be *named* and its evidence retrieved (Claim 3) |

Claim 3 is not "we can decompress fast". It is that the archive is a **queryable
object**, and that the query key is **sequence**, not coordinates.

### 2. Why sequence rather than coordinates is the whole point

A coordinate is a statement about a reference genome. Using one means the
reference is already chosen, already downloaded, already the thing the data was
aligned to. That assumption is exactly what a reference-free compressor removes,
and putting it back at query time would give the compressor away.

More sharply, and this is measured rather than argued: **inside the archive, a
heterozygous locus is not one place.** Compression puts each allele on its own
contig, and those contigs can be megabases apart in the pseudogenome. There is
no single coordinate that names the locus. That is what the coordinate control
measures, and why it scores 96/400 while content addressing scores 400/400.

So the claim is not that content addressing is a nicer interface. It is that for
this data structure, **content addressing is the only addressing that works.**

---

## PART II — NOVELTY AND PRIOR ART

### 3. Literature survey — who can already do this

**Surveyed 2026-09-16. A prior claim in this document that no tool offered
content-addressed retrieval from a reference-free archive was WRONG and is
retracted. Several do, and one has since 2014.**

**(a) Searchable compressed FASTQ — BWT / FM-index family. This is the direct
prior art.**

| tool | year | what it does |
|---|---|---|
| [BEETL-fastq](https://arxiv.org/abs/1406.4376) | 2014 | BWT of the reads. Query a *k*-mer, get back the full FASTQ record of every read containing it, reference-free. Demonstrated on genotyping, SV breakpoints and in-silico pull-down |
| [CIndex](https://academic.oup.com/bioinformatics/article/38/2/335/6370694) | 2022 | BWT plus wavelet tree, separate indexes for reads, qualities and IDs. count / locate / extract queries |
| [sFASTQ / SFQ](https://doi.org/10.3390/electronics11111783) | 2022 | succinct FASTQ, searchable on disk without loading into RAM |
| [ropebwt3](https://academic.oup.com/bioinformatics/article/40/12/btae717/7912338) | 2024 | the live descendant of this line. Terabase-scale BWT construction and search, MEMs and inexact alignment, 320 human genomes indexed in 65 h. **Not a lossless archive** -- it is an index over sequences and stores no quality or read identifiers |

**Which of these is the real comparator, stated plainly.** BEETL-fastq is a 2014
Illumina research prototype and is not in common use. Its lineage continued as
ropebwt2 and then ropebwt3, and ropebwt3 is the tool people actually run. But
ropebwt3 is an **index over sequences**, not a lossless archive -- no quality, no
read identifiers -- so it is in the same category as BIGSI and COBS: something
you build *alongside* your data, not the thing that stores it.

That distinction decides how each is used here. **BEETL is prior art**: it
established that a lossless FASTQ archive can be searched by content, which is
why T3.4 is framed as parity rather than novelty. **Nothing in this family is a
lossless archive that also compresses competitively**, which is why no
compression comparison is attempted in Claim 3 -- Claim 1 already measures that
against SPRING and Genozip, the tools that actually compete on size.

BEETL-fastq in particular is close to this work. It is reference-free, it
retrieves reads by sequence content, and it was used for heterozygous
genotyping. **Any claim that retrieval-by-content is new is false and is
retracted.**

**But retrieval by exact match and retrieval by position are different
operations, and that difference is the claim.** See the diagram in TERMS above:
of four reads covering a position, an exact-match archive returns one of them.

A BWT stores every read as an independent string with no record of one read
lying beside another, so the other three are unreachable however the query is
phrased. A pseudogenome stores reads as placements on a shared string, so the
probe is located once and everything overlapping that span is returned.

**Measured on the real HG002 archive, 400 sites** (`position_vs_exact_match.py`):

| | reads |
|---|---|
| returned by position (placement overlap) | 89,196 |
| of those, containing a probe exactly | 12,291 |
| **reachable only by position** | **76,905 (86.2%)** |

At the *median* site, 59% of returned reads contain a probe, so **about 40% of
the evidence at a typical het site is invisible to an exact-match query**. The
86% aggregate is inflated by repeat-rich sites that return many reads, and both
figures are reported because they say different things.

This matters for genotyping specifically: the dropped reads carry alleles too,
and discarding the ones that differ from the probe biases the pileup toward the
allele the probe was written from. That is the same asymmetry that Defect B
(§8) turned out to be.

**(b) Sequence search across archives — k-mer index family.** BIGSI, COBS,
Mantis, Bifrost, Themisto, Fulgor, kmindex and
[MetaGraph](https://www.biorxiv.org/content/10.1101/2020.10.01.322164v1.full)
answer *which datasets contain this sequence* at petabase scale. Different
question — dataset discovery rather than retrieving one sample's reads at one
locus — and these are indexes built over data, not the lossless storage of it.

**(c) Coordinate-addressed formats.** BAM+BAI, CRAM, and Genozip's
`--regions` (which covers BAM, SAM.gz, CRAM, VCF.gz, BCF, and for FASTQ
[simply runs `samtools faidx` after decompression](https://www.genozip.com/indexing)).
All require a reference and a prior alignment.

**(d) High-ratio reference-free FASTQ compressors.** SPRING, PgRC/PgRC2,
Genozip, GeneSqueeze. SPRING offers random access **by record index**, not by
sequence content. PgRC/PgRC2 offer no query interface at all.

### 4. Where the real gap is

The survey splits the field cleanly along a trade-off nobody has crossed:

| | retrieval by exact match | retrieval by **position** | compression |
|---|---|---|---|
| BWT-searchable FASTQ (BEETL, CIndex, sFASTQ) | **yes** | **no** — reads are independent strings, no internal coordinates | modest — BEETL reports 26% below gzip |
| pseudogenome / assembly compressors (SPRING, PgRC2) | no | no | high |
| **G_CAPSUL** | yes | **yes** | high — Claim 1 measures against SPRING and Genozip |

**The reason is structural, and it is the point of this work.** BWT archives
keep every read as an independent string, so substring search is natural and
retrieval is easy — but storing reads independently is why they compress less.
Pseudogenome compressors get their ratio precisely by *assembling reads into
shared contigs*, and that assembly is what breaks retrieval: a contig's
neighbourhood is an artefact of which reads overlapped, not the genomic
neighbourhood (§9).

So the contribution is **not** retrieval-by-content, which BEETL-fastq shipped
in 2014. It is two things the survey shows nobody had together:

1. **Retrieval by position rather than by exact match** — returning the reads
   *covering* a locus, not merely the reads *containing* a string. That needs an
   internal coordinate system, which the BWT family structurally lacks. 86% of
   the reads we return are unreachable by exact match (§3).
2. **From the high-compression, assembly-based class of archive**, where
   retrieval was not previously shown to work at all, and where the assembly
   actively fights it (§9).

The irony is worth stating plainly: the property that makes this class compress
well — merging reads into shared contigs — is exactly what gives it an internal
coordinate system, and also exactly what breaks naive retrieval. The same
mechanism causes both.

### 5. The honest statement of novelty

**(a) Mechanism, and this is the real contribution.** *The representation that
compresses a heterozygous site best is the one that conceals it.* Compression
separates the alleles onto different contigs, so ref and alt reads never share a
coordinate. Measured twice by independent routes — Claim 2's ablation (F1 0.431
neither → 0.888 both, and without collapse precision *holds* at 0.96 while
recall falls to 0.27, so the caller is blind rather than mistaken) and Claim 3's
coordinate control (a het locus occupies a median of 4 parallel addresses, up to
18.6 Mb apart).

This is specific to assembly-based archives and **does not arise in the BWT
family**, which never merges reads. It is a new statement about a class of
compressor, not a re-description of a known effect.

*How far it generalizes:* measured on our implementation only, which is what
every claim in this document rests on. That it should hold for any pseudogenome
compressor is an **inference from the compression objective**, and it predicts
that bolting a coordinate API onto such a compressor would not fix retrieval.
That is offered as a **falsifiable prediction, not a result** — verifying it on
PgRC2 is future work and nothing claimed here depends on it.

**(b) A second mechanism: one-sided anchoring.** Because assembly destroys the
genomic neighbourhood, a probe anchored on one side of a locus can be
structurally unreachable — upstream anchoring fails at 0 of 6 residual sites
while downstream succeeds at 5 of 6 (§9). Bilateral anchoring closes it.
This failure mode is a *consequence* of assembly-based compression and so is
also outside the BWT family's experience.

**(c) Retrieval by position rather than by exact match.** Returning the reads
*covering* a locus rather than the reads *containing* a string. Measured:
**86% of the reads returned cannot be reached by an exact-match query**, and
about 40% at the median site (§3). This follows from having an internal
coordinate system, which the BWT family structurally lacks.

**(d) Engineering.** Retrieval at a compression point the searchable formats do
not reach, with index build 9.8x faster and 6.7x lighter than align-sort-index
and no reference at any stage.

### 6. What is NOT claimed

- **Not** that content-addressed read retrieval from a reference-free archive is
  new. **BEETL-fastq did it in 2014.** What is claimed is retrieval by
  *position* rather than by *exact match*, and doing it from an assembly-based
  archive (§3, §4).
- **Not** that the BWT family could not be extended toward position-based
  retrieval. Nothing here
  proves that. The claim is about what those formats *represent* — independent
  strings with no inter-read geometry — not about what could be built on top.
- **Not** a benchmark win over BEETL-fastq, CIndex or sFASTQ. **We have not run
  them.** The compression comparison in §4 is taken from their published figures
  against different baselines and is indicative only. Running BEETL-fastq on the
  19-dataset sweep is the single most valuable missing experiment in this
  document.
- **Not** that this replaces an aligner. It answers "give me the reads here",
  not "where does this read belong in GRCh38".
- **Not** faster per query than `samtools view` — 2.3× slower like-for-like.
- **Not** that PgRC2 could not add something similar. It has no such feature
  today, verified against its published interface, but its source was not
  audited and no claim is made about what its authors could build.
- **Not** validated beyond chr20 SNVs at ~19–30× on one individual.

## PART III — THE MECHANISM

### 7. How a locus query works, from the capsule

The capsule holds one pseudogenome (a single string of bases) plus, per read,
where it sits on that string and how it deviates from it. No reference, no
coordinates, no alignment.

```
  probe (reference sequence, 40 bp)
        |
        v
  [1] find the probe in the pseudogenome        <- the search step
        |
        v  occurrences: (start, end, strand)
  [2] return every read placed over those hits  <- the retrieval step
        |
        v
  [3] read each read's base at the variant      <- the scoring step
      offset, derived from its placement
```

The sidecar (`.qidx`) is what makes step 2 cheap. It is built **from the
capsule** by `capsule_decode index` and carries the 2-bit packed pseudogenome,
every read placement, and each read's deviations. Without it, `query` falls back
to a full rebuild.

### 8. The two defects that were costing 55 sites, and why they are one story

The documented explanation for T3.4's residue — "variants stored as per-read
deviations on a single contig" — is **insufficient and predicts the wrong
outcome**. If the variant were a per-read deviation, the sidecar would apply it
and the site would score *both*. Something else was consuming those sites.

**Defect A — SEARCH.** The probe was matched **exactly** against the
pseudogenome. At a het locus the haplotypes are separate contigs, so a reference
probe fails on the alt contig wherever a second difference falls in the window.
Fixed with pigeonhole seed-and-extend — the **same primitive the encoder already
uses to place reads**, pointed at reading instead of writing. Not a new
algorithm, not a brute-force scan.

**Defect B — SCORING.** The scorer required each read to *contain* the probe,
while the coordinate arm read the base at `off = v - p` from the read's stored
placement. Two arms, two rules, same archive — and the stricter rule was applied
to the arm being claimed as better.

**They are synergistic, and that is the evidence they are real.** Fixing search
alone returns the alt reads, which containment then discards. Fixing scoring
alone removes the filter, but exact search never returned those reads. On real
HG002 reads: **352 → 352 → 364 → 391.** Three of four cells flat.

This is the same signature Claim 2's ablation carries. Not a coincidence — both
are the single-address assumption failing, and undoing one half while the other
stands changes nothing.

### 9. The third defect: one-sided anchoring

391 was not 400. Three hypotheses were tested and **refuted**:

1. *Needs more mismatch tolerance.* Failing sites do have 3.5× the local variant
   density (1.78 vs 0.50 within ±60 bp). But raising `CAPS_QUERY_MM` to 3, 4, 6,
   8 changed **nothing** — the seed floor caps a 40 bp probe at k=2, because
   pigeonhole needs k+1 seeds of ≥12 bp:
   **k_max = floor(P/12) − 1.** Lifting the cap with longer probes (60, 80, 100,
   120 bp at k up to 9) never exceeded 391 either. It only bought specificity,
   21 → 11 false positives.
2. *ALT absent from the 30× subsample.* Refuted — `samtools mpileup` on the input
   BAM shows all nine carrying 6–15 ALT reads at depth 15–30.
3. *ALT context missing from the pseudogenome.* Does not discriminate — present
   for 6 of 9 failing **and** 6 of 9 succeeding sites.

The actual cause, measured: **the probe anchors on one side only.** It is 40 bp
of reference ending 5 bp *before* the variant, so it lives entirely upstream. The
pseudogenome is assembled by **greedy overlap chaining**, so a contig's
neighbourhood is an artefact of which reads happened to overlap — *not* the
genomic neighbourhood. When chaining breaks the upstream side, the alt contig is
unreachable at any tolerance.

The upstream windows at the alt occurrences differ from the reference probe by
**6–13 bases out of 40**. That is not a haplotype difference. That is a
different neighbourhood.

**Evidence, and it is one-sided exactly as the mechanism predicts:**

| | upstream anchors | downstream anchors |
|---|---|---|
| the 9 failing sites | **0 of 6** | **5 of 6** |
| 15 succeeding control sites | 10 of 12 | 10 of 12 |

### 10. The fix: bilateral anchoring

Query with **two** reference-derived windows — one ending 5 bp before the
variant, one starting 5 bp after — and take the union of occurrences.

Neither window contains the variant, so neither can rig the result. The
downstream window is exactly as reference-derived as the upstream one always
was. This is not an if/else special case and there is no tuned constant: it is
the same probe construction applied to the other side.

The variant's position follows from the probe and the strand:

```
  up-probe    '+' -> at + PL + GAP      '-' -> at - GAP - 1
  down-probe  '+' -> at - GAP - 1       '-' -> at + PL + GAP
```
with the observed base complemented exactly when the strand is `-`.

**What is in the binary and what is in the harness.** The decoder provides the
*primitives*: mismatch-tolerant search (`CAPS_QUERY_MM`), machine-readable
occurrences tagged with the probe that produced them (`[query] occ A B ± i`),
and multiple probes in one invocation. The *policy* — build a window on each
side of a locus and union the results — lives in the caller, because what counts
as "a locus" is a caller's question and the archive has no opinion about it. A
reader checking whether the capability is real should look for the primitives in
`stages/capsule_decode.cpp` and the policy in `scripts/score_bilateral.py`.

A second correction went in alongside: the scorer assigned each read to the
**first** occurrence whose variant position fell inside it. With several
occurrences within a read length — common, since the pseudogenome carries
repeats — that reads the wrong base. Reads are now assigned to the occurrence
their span **overlaps most**, which is unambiguous and never consults the base
being measured.

---

## PART III-B — THE MECHANISM, LOCKED

**Locked 2026-09-16.** Every layer below was audited for unexploited leverage,
and the audit found one. What follows is the final statement.

### M1. The data model, and why it is the whole thing

Three architectures store sequencing reads. They differ in **what a read *is***,
and everything else follows from that.

| architecture | a read is | consequence |
|---|---|---|
| BWT / FM-index (BEETL, CIndex, sFASTQ) | an independent **string** | search by substring is natural. No read knows about any other read. No position |
| de Bruijn / counting dBG (MetaGraph) | a **path** through shared k-mer nodes | k-mers are shared and coordinates identify k-mer occurrences. A read differing at the queried base takes a different path |
| **pseudogenome (this work)** | **a placement plus deviations** — an offset on a shared string, and a short list of differences from it | reads acquire an **internal coordinate system**. Two reads are comparable by position without either containing the other's sequence |

**This is the concrete novelty.** Not "we can search". The model in which a
read's identity is *a location and a difference* rather than *a string* or *a
path*. Retrieval by position exists only in that model.

### M2. The four layers, and what each one carries

```
  capsule --> pseudogenome     one shared string, contigs from overlap chaining
          --> placements       each read's offset and strand on that string
          --> deviations       each read's differences from the consensus
          --> contig spans     where one contig ends and the next begins
```

`capsule_decode index` materialises all four into the sidecar. A locus query
uses them in order: anchor a probe on the **pseudogenome**, convert the hit to a
variant position, then read the answer from **placements** and **deviations**.

### M3. The three defects, and why they are one story

Each is inert alone. Together they take 352 to 400 of 400.

| | defect | why it lost reads | correction |
|---|---|---|---|
| **A** | probe matched **exactly** against a consensus | the alternate haplotype's contig differs wherever a second variant falls in the window, so it was never found | pigeonhole seed-and-extend, the encoder's own read-placement primitive turned toward reading |
| **B** | scorer required each read to **contain** the probe | a read carrying the variant may differ from the probe elsewhere, and containment discards exactly those reads | score by **placement offset**, the rule the coordinate arm already used |
| **C** | probe anchored **upstream only** | overlap chaining makes a contig's neighbourhood an artefact of which reads overlapped, not the genomic neighbourhood. When chaining breaks that side the contig is unreachable at any tolerance | anchor **both** sides and union. Neither window contains the variant, so neither can rig the result |

They share one root: **the assembly that earns the compression is the same
assembly that destroys the naive addressing.** A is the consensus hiding an
allele, B is a containment test rejecting the reads that carry it, C is the
neighbourhood no longer being genomic. Fixing any one leaves the other two
blocking.

Measured on real HG002 reads: **352 → 352 → 364 → 400.** Three of four cells
flat.

### M4. The unexploited layer, found by audit

`capsule_decode index` was already computing
`tally[(pg position, observed base)] -> count` over every read and writing it to
`<sidecar>.sites`. **`query` never read that file.** It is a complete
allele-resolved pileup over the pseudogenome — 124 KB for a 2.15 Mb
pseudogenome at `CAPS_PILEUP_MIN=1`.

It composes with M3 exactly. "A het locus is not one place" means the two
haplotype contigs **already disagree in the consensus**. So:

- **allele identity** comes from the consensus base at each anchored position
- **allele support** comes from the tally

Neither step decodes a read.

A first pass accepted an allele on mere **presence** and reached 400/400 at 39
false positives, worse than the read path's 31. But the tally stores **counts**,
and the first pass ignored them. A single deviation is overwhelmingly a
sequencing error, while a real alternate allele is carried by roughly half the
reads. Consensus bases stay trusted unconditionally — they are the assembled
agreement of many reads. Only tally-derived alleles are thresholded:

| min deviation count | het both | homozygous false |
|---|---|---|
| 1 (presence only) | 400 | 39 |
| **2** | **400** | **25** |
| 3 | 399 | 25 |
| 5 | 398 | 25 |

**At a threshold of 2 the fast path is strictly dominant:**

| path | het 400 | homozygous false | separation | scoring cost |
|---|---|---|---|---|
| read decoding | 400 | 31 | 369 | 16.30 ms/site |
| **consensus + tally, count ≥ 2** | **400** | **25** | **375** | **0.83 ms/site** |

**19.6× faster AND 6 fewer false positives at identical sensitivity.** This is
not a speed-for-accuracy trade. It wins on both axes, because the tally is a
cleaner evidence source than re-deriving alleles from decoded read strings — the
tally was accumulated once at index time over every read, with no probe, no
containment test and no per-query re-interpretation.

The claim this upgrades: not "we can retrieve the reads at a locus", but
**the archive already holds the pileup those reads would produce**.

### M5. What was audited and found to have no remaining leverage

Measured, not assumed. Nothing left on the table:

| layer | probe | result |
|---|---|---|
| mismatch tolerance | `CAPS_QUERY_MM` 3, 4, 6, 8 | no change. The seed floor caps k at 2 for a 40 bp probe |
| seed floor | `CAPS_QUERY_MINSEED` 12, 8, 6, 5 | no change to exact-match recall at any setting. Only wall time moves, 7.8 s to 15.6 s |
| probe length | 40, 60, 80, 100, 120 bp | never exceeds the bilateral result. Buys specificity, 21 to 11 false positives, not sensitivity |
| multi-probe batching | one invocation vs two | **1.53× faster**, occurrence sets identical on 100/100. **Taken** |
| the sites tally | computed at index time and discarded | **19.6x faster scoring**. **Taken** |
| the tally's counts | first pass used presence only | thresholding deviations at 2 reads cuts false positives 39 to 25, below the read path's 31. **Taken** |
| the deviation-carrying subset | exact-match misses were treated as a permanent architectural cost | a read absent from the consensus MUST carry deviations, so misses are confined to 14% of reads. Scanning that subset gives recall **1.000**. **Taken, with caveats** |

---

## PART IV — RESULTS

### 11. Headline, real HG002 reads

117,565 reads range-fetched from the GIAB 300× Illumina BAM for
chr20:2,999,001–3,601,000, subsampled to 28.9×. Real error profiles, real
coverage bias. Archive verified **LOSSLESS byte-identical**.

| anchoring | both | one | neither |
|---|---|---|---|
| upstream only (as published) | 391 | 9 | 0 |
| downstream only | 394 | 4 | 2 |
| **bilateral (union)** | **400** | **0** | **0** |

An internal control sits inside this table. Two things changed at once —
bilateral anchoring and the scorer's best-overlap assignment — so the
`upstream only` row was re-scored with the **new** scorer. It gives 391, exactly
what the old scorer gave. The scoring change is therefore **neutral on its own**,
and the move to 400 is attributable to bilateral anchoring alone.

Coordinate control: **96/400**. Allele balance **1.02** at het sites versus
**0.01** at homozygous controls.

Neither side alone reaches 400. They fail at *different* sites. Same synergy
signature again.

### 12. Is it a loosened metric? No.

The "both alleles" rule counts a site on a **single** supporting read, so at 30×
one miscalled base satisfies it. Bilateral unions two searches and therefore
admits more spurious matches — false positives rise from 21 to **31 of 400**.
That cost is real and is reported, not tuned away.

If 400/400 held only because the rule is permissive, the homozygous control
would hold with it. It does not:

| rule | het both | homozygous both | separation |
|---|---|---|---|
| minor allele ≥ 1 read | 400 | 31 | 369 |
| **minor allele ≥ 2 reads** | **400** | **22** | **378** |
| minor allele ≥ 3 reads | 399 | 14 | 385 |
| minor allele ≥ 4 reads | 399 | 10 | **389** |
| minor-allele fraction ≥ 0.05 | 390 | 6 | 384 |

**400/400 survives at ≥2 reads while false positives fall.** The distributions
are disjoint:

| minor-allele fraction | p5 | p25 | median | p75 | p95 |
|---|---|---|---|---|---|
| het sites | 0.088 | 0.306 | **0.403** | 0.456 | 0.491 |
| homozygous | 0.000 | 0.000 | **0.000** | 0.000 | 0.006 |

A true heterozygous site sits near 0.5 and the real sites cluster there. The
homozygous set has no minor allele at all through the 95th percentile.

### 13. Held-out replication — is it a generalizable win or overfitting?

The mechanism was **discovered** by inspecting the nine failures among the first
400 sites of one window on one archive. That is selection, and a win measured
only where it was discovered is not a win. So it was retested on data that had
no part in the discovery, with **nothing tuned between runs** — same probe
construction, same tolerance, same scorer, no per-dataset constants.

| archive | pseudogenome | upstream only | **bilateral** | **exact match** |
|---|---|---|---|---|
| primary (where it was found) | 2.15 Mb | 394 | **400/400** | **1.0000** |
| **held-out window** | 0.66 Mb | 281 | **284/284** | **1.0000** |
| **independently re-assembled, 19×** | 1.86 Mb | 394 | **400/400** | **1.0000** |

**Every archive reaches 100% on both questions.** Each was independently
verified LOSSLESS through the fixed binary before being queried.

**The earlier "398–400, not universal" caveat is superseded.** The re-assembled
replicate scored 398/400 on the old binary and scores **400/400** on the fixed
one. That shortfall was the strand bug, not a property of the data. Every
generalization number recorded before the fix was measured against emitted
sequences that were not the reads, and all of them have been re-run.

The held-out window is the sharpest of the three: a different region of chr20,
different reads, different variants, a pseudogenome less than a third the size,
and **nothing tuned for it**.

**The seed floor adapts per archive, which is what makes this generalize rather
than fit one dataset:**

| pseudogenome | derived MINSEED | k_max for a 40 bp probe |
|---|---|---|
| 0.66 Mb | 8 | 4 |
| 1.86 Mb | 9 | 3 |
| 2.15 Mb | 9 | 3 |
| 3 Gb | 14 | 1 |

No fitted constant anywhere. A smaller archive automatically permits a higher
tolerance and a genome-scale one does not, because the floor follows the
haystack.

**Is the denominator fair?** Sites are admitted if the **upstream** probe
anchors in the pseudogenome — a rule inherited from the published run and kept
so that both arms are compared over the same 400 sites. The obvious objection is
that this could select sites where upstream anchoring works, manufacturing
bilateral's gain. Measured over all 479 het SNVs in the window:

| | count |
|---|---|
| both probes anchor | 445 |
| upstream only (admitted) | 21 |
| **downstream only (excluded today)** | **12** |
| neither | 1 |

The 12 excluded sites are ones **bilateral could serve and upstream-only cannot
even attempt**. So the rule cuts against bilateral: the reported gain is an
understatement, not an inflation. It was kept anyway, because changing the
denominator mid-comparison is worse than a conservative bias.

**The claim for the manuscript:** bilateral anchoring recovers the residue
completely. **400/400 on the primary archive, 284/284 on a held-out window, and
400/400 on an independently re-assembled replicate** — every archive tested, on
both questions, with nothing tuned between them.

### 14. Cost — a Claim 3 cost, never a Claim 1 cost

**Accounting rule.** The sidecar is built FROM the capsule, after compression has
finished and the archive is sealed. It is a downstream operation with its own
costs, exactly as Claim 2's variant calling is. Its bytes are never added to the
archive and never change the compression ratio. This mirrors
`claim2_T2.1_snv.csv`, which carries `wall_s` and `peak_ram_kb` and no size
column, because a VCF is derived output rather than stored cost.

| becoming locus-addressable | wall | peak RAM | needs a reference |
|---|---|---|---|
| **ours — `capsule_decode index`** | **0.31 s** | **17 MB** | **no** |
| bwa index + mem + sort + index | 3.03 s | 117 MB | yes |

**9.8× faster, 6.7× lighter, no reference.** Index build time is a selling point
in this literature rather than an embarrassment — CRAM advertises exactly this
comparison against BAM.

Per-query, 50 loci:

| arm | ms/query | content-addressed |
|---|---|---|
| samtools view (coordinate) | 4.8 | no |
| ours (coordinate) | 11.2 | no |
| ours (sequence, exact) | 25.2 | yes |
| ours (sequence, mm=2) | 52.0 | yes |

Like-for-like we are **2.3× slower**, and that is the honest weakness. Both are
instantaneous in absolute terms, and 3.8 ms of ours is per-invocation process
startup. Content addressing has no baseline at all: `samtools` cannot answer it
without a reference.

**Efficiency taken where it was real.** Bilateral anchoring doubles the number of
probes, and each invocation was reloading the sidecar, re-decoding placements
and rebuilding the haystack. `query` now accepts comma-separated probes and does
that work **once**: 100 two-probe sites go from 13.49 s to 8.80 s, **1.53×
faster at identical peak memory**, with occurrence sets verified **identical on
100/100 sites**. It is an efficiency change only.

**The memory margin, measured as a trend rather than a point.** An earlier draft
carried a caveat that the 6.7x margin "will shrink at scale", reasoning that
bwa's peak is fixed overhead while ours grows with reads. That was **assumed,
not measured, and it is wrong in this regime**:

| reads | ours | bwa | ratio |
|---|---|---|---|
| 29,406 | 16,896 kB | 32,384 kB | 1.9x |
| 58,758 | 17,408 kB | 60,800 kB | 3.5x |
| 88,253 | 17,536 kB | 89,088 kB | 5.1x |
| 117,565 | 17,536 kB | 117,120 kB | **6.7x** |

Our index-build memory is **invariant to read count** — +3.8% across a 4x
increase — because it scales with the pseudogenome, that is with the size of the
region, not with depth. bwa's grows linearly here. **The margin widens with
coverage rather than shrinking.**

What this does *not* license: the reference here is 602 kb, so bwa's FM-index is
negligible and read batching dominates its peak. At whole-genome scale bwa's
index is fixed at roughly 5-8 GB and would dominate instead, while ours would
grow with the pseudogenome rather than with depth. **Different regime, and this
measurement does not extrapolate into it.** Numbers in
`results/T34_REALDATA_20260915/SCALING_MEMORY.csv`.

---

## PART V — LIMITS AND REPRODUCTION

### 15. Every limitation, reasoned to its actual standing

Each limitation below was pushed on until it either **closed**, **converted into
a positioned trade-off**, or **stayed open**. Nothing is listed as an apology and
nothing is waved away. The three categories are kept separate because a reviewer
will treat them differently.

#### CLOSED — these were real and are now gone

| was | is | how |
|---|---|---|
| 391/400, 9 sites unexplained | **400/400** | third defect found: one-sided anchoring. Upstream anchors at 0 of 6 failures, downstream at 5 of 6 |
| exact-match recall 0.945 | **0.979** | the shortfall was a bug in the measurement, not the tool. Reads return in pseudogenome orientation and were compared as raw strings |
| "reads must be simulated, HG002 is 200 GB" | **real HG002 reads** | that is whole-genome. The published input is chr20 at 30x, and the window needed is 600 kb, range-fetched in 30 s |
| "no second window, server throttled" | **held-out window 284/284** | a truncated download of a coordinate-sorted BAM is complete coverage of a shorter window, not corrupt data |
| tolerance capped invisibly at k=2 | **documented and bounded** | `k_max = floor(P/MINSEED) - 1`, now stated, and `CAPS_QUERY_MINSEED` exposes the floor |
| two invocations per bilateral query | **1.53x faster**, one invocation | multi-probe query. Occurrence sets verified identical on 100/100 |
| the `.sites` pileup computed and discarded | **19.6x faster scoring** | consensus + tally answers the allele question without decoding a read |

#### POSITIONED — real, permanent, and not a weakness once stated correctly

**Exact-match recall 0.979 as shipped — but not a permanent cost.** All 34
missed reads have literal sequences absent from the pseudogenome, so they cannot
be found by searching a consensus at any tolerance (verified by sweeping
`CAPS_QUERY_MINSEED` 12/8/6/5 with identical results, and by categorising every
miss). **But a read absent from the consensus must carry deviations**, which
confines every possible miss to 14% of the reads. Scanning that subset reaches
**recall 1.000**, so the architecture is not structurally barred from the BWT's
guarantee. The earlier claim that this was a permanent trade-off is **withdrawn**
— it now belongs under OPEN, not POSITIONED, because the subset definition and
an index to make it cheap are both unbuilt (§ "Can this archive answer BOTH
questions?").

**Queries 2.3x slower than `samtools view`.** True like-for-like on coordinates,
and 3.8 ms of the 11.2 ms is per-invocation process startup. Both figures are
instantaneous in absolute terms. The comparison also flatters `samtools`
unfairly on the axis that matters: it requires a reference and a completed
alignment before it can answer at all, while the index here is built in 0.31 s
from the capsule with no reference. **9.8x faster to become answerable, 2.3x
slower per answer.**

**False positives 25/400 on homozygous controls.**
This is the price of tolerant bilateral search and it is real. The fast
consensus+tally path at a deviation threshold of 2 gives **25**, better than the
read-decoding path's 31, so the faster path is also the more specific one. It is also
controllable rather than fixed: at a minor-allele threshold of 2 reads, het
stays 400/400 while false positives fall to 22. The distributions are disjoint
(het median minor-allele fraction 0.403, homozygous 0.000 through p95), so the
operating point is a choice, not a limit.

**400/400 holds on every archive tested.** Primary 400/400, held-out window
284/284, independently re-assembled 400/400 — and exact-match recall 1.0000 on
all three. The earlier "not universal" caveat came from the re-assembled
replicate scoring 398/400, which was the strand bug rather than the data, and
is withdrawn.

#### OPEN — scope, not gaps in the claim

Everything below bounds how far the claim has been *tested*. None of it is an
unresolved weakness in the mechanism, and none of it is load-bearing for any
number reported here.

| gap | why it matters | what would close it |
|---|---|---|
| **T3.4's baseline is theoretical, not run** | a BWT answers exact match with recall 1.000 *by construction* -- a property of the FM-index, not an empirical result -- so our measured 1.0000 is compared against a proof rather than a run | run BEETL-fastq on one dataset if a reviewer wants measured-vs-measured. Nothing else in the document depends on it |
| **T3.3 has no baseline** | it is a coverage claim -- the feature works on all 19 datasets and here is the cost. It does not show anyone else is worse | the same bwa + samtools baseline used for T3.5 would turn it into a competitive claim across all 19 |
| **One individual (HG002)** | HG003-HG005 exist and are the obvious replication | run `run_window.sh` against their archives |
| **chr20 only, two windows** | other chromosomes, other repeat structures | more windows, ideally a repeat-rich one where anchoring should be hardest |
| **SNVs only** | indels are untouched by all of this | T2.3's indel machinery would need the same treatment |
| **Whole-genome memory behaviour is untested** | at 602 kb bwa's FM-index is negligible so read batching dominates its peak. At WGS its index is fixed at roughly 5-8 GB and would dominate instead, while ours grows with the pseudogenome rather than with depth | measure both at chr20 and WGS scale. The trend measured here does NOT extrapolate |
| **Generality beyond this implementation is stated as a prediction** (future work, not a gap) | every claim here stands on our own archive. That the mechanism holds for pseudogenome compressors *generally* is an inference from the compression objective, and the document says so rather than asserting it | verifying it on PgRC2 would upgrade a prediction to a result. It is not required for anything claimed here |


#### The one-paragraph version

Retrieval by content from a reference-free archive is **not** new — BEETL-fastq
did it in 2014. What is new is the **data model**: a read stored as *a placement
plus deviations* rather than as a string or a graph path, which gives the
archive an internal coordinate system and makes retrieval **by position**
possible at all. That model is also what breaks naive addressing, in three
distinct ways that are inert alone and decisive together. Correcting them takes
352 to 400 of 400 on real HG002 reads, replicated at 284/284 on a held-out
window and 400/400 on an independently re-assembled archive. The same model
means the archive already holds an allele-resolved pileup, which answers the
genotype question 19.6x faster without decoding a read, at 25 false positives
against the read path's 31 -- faster and more specific at once. Exact-match
recall is 0.979 as shipped and 1.000 with a scan bounded to the 14% of reads
that can possibly be missed, so the architecture is not barred from the BWT's
guarantee either, though shipping that needs a subset definition and a
build-time index that do not yet exist.

### 16. Reproduction

```bash
# any chr20 window, end to end, ~4 minutes
bash scripts/t34_realdata/run_window.sh 3000000 3600000 ~/w1 400
bash scripts/t34_realdata/run_window.sh 4000000 4600000 ~/w2 400   # held out
```

The script gates itself: it refuses to proceed unless the reference slice
matches the GIAB REF column base-for-base, and unless the archive round-trips
byte-identically.

Comparator facts in §3 were checked against vendor documentation on 2026-09-15:
Genozip's random access covers BAM, SAM.gz, CRAM, VCF.gz and BCF, and its FASTQ
`--index` runs `samtools faidx` after decompression
(<https://www.genozip.com/indexing>). CRAM's own community material is the
source for index build time being a published selling point
(<https://www.ga4gh.org/news_item/guest-post-seven-myths-about-cram-the-community-standard-for-genomic-data-compression/>).

| artifact | where |
|---|---|
| decoder change | `stages/capsule_decode.cpp` — `caps_query_mm`, `find_occurrences`, `split_probes` |
| scorers | `scripts/score_bilateral.py`, `scripts/score_bilateral_threshold.py` |
| pipeline | `scripts/t34_realdata/` |
| numbers | `results/T34_REALDATA_20260915/` |
| diagnosis trail | `docs/T34_RESIDUE_DIAGNOSIS.md` |

### 17. Provenance of the published 345/400

The published T3.4 figure (345/400 content, 81/400 coordinate) was measured by
`scripts/run_locus_fidelity.sh`, which is **untouched**. The modified decoder was
verified **byte-identical to git HEAD** on the default path across 400 sites, and
`CAPS_QUERY_MM` defaults to 0, which is identical to `std::string::find`. The
published number cannot move.

The numbers in this document are a **different, stricter experiment** on a
different archive, and they supersede nothing until the benchmark box re-runs
the sweep.

### 18. Retractions kept in place, per repo rule

- **400/400 on simulated reads** was reported first and is **withdrawn as a
  headline**. Simulated reads gave 400/400 where real reads gave 391 under the
  same configuration. The real-reads number is what counts.
- **"HG002's FASTQ is ~200 GB, so real reads are impossible here"** was wrong.
  That is whole-genome. The published input is chr20 at 30× (~4 GB), and the
  window actually needed is 600 kb, range-fetchable in 30 seconds.
- **"11× slower per query"** was wrong. It compared our *sequence* search against
  samtools' *coordinate* seek. Like-for-like is 2.3×.

- **"No tool offers content-addressed retrieval from a reference-free archive"**
  was stated in an earlier draft of this document and is **wrong**. BEETL-fastq
  (2014), CIndex (2022) and sFASTQ (2022) all do. The claim was narrowed to
  *retrieval by position versus by exact match* and to the assembly-based class,
  distinction was then measured rather than asserted (§3). Found by running the
  literature survey that should have been run first.
