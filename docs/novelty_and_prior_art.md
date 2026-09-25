---
Date: 2026-09-19
Title: Novelty and Prior Art — Cross-Claim, Honest, Confirmed vs Inferred
Purpose: Consolidates the honest comparator analysis from all three claim
  folders into one reference, keeping confirmed facts and inferred
  arguments strictly separate throughout.
When to refer to this file: Writing the paper's related-work discussion;
  answering any reviewer question of the form "why is X not compared
  directly"; checking whether a novelty claim is source-verified or
  argued from first principles.
Keywords: novelty, prior art, PgRC2, BEETL, CIndex, sFASTQ, DiscoSNP++,
  Kmer2SNP, SPAdes, bwa, confirmed, inferred, comparator justification
---

# Novelty and prior art — the honest, consolidated picture

## Rule followed throughout this file, and throughout the three claim folders it draws from

Every claim below is labeled **confirmed** (verified by reading source code
or a primary document this session or a prior, cited session) or
**inferred** (a reasonable argument from documented behavior, not verified
against the comparator's own source). This distinction was enforced
strictly in each claim folder's mechanism file and is preserved here, not
flattened into one confident voice.

## Compression comparators — Claim 1

- **SPRING, Genozip**: general-purpose FASTQ compressors, the correct
  primary comparators — both handle full FASTQ (sequence + names + quality
  + line-3), matching G_CAPSUL's own scope. **Confirmed**: 19/19 wins vs
  both, independently re-derived from raw data this session.
- **PgRC2**: the closest architectural relative (also pseudogenome-based).
  **Confirmed not directly comparable** to the whole-FASTQ numbers — PgRC2
  stores sequence only, no names/quality/line-3. Correctly scoped
  sequence-only comparison: +1.88%, 6 wins/1 loss, ~1.7x slower/2.5x
  heavier (7 datasets both tools can process). See
  `claim1/mechanism_insight_claim1.md`.

## Variant-calling comparators — Claim 2

- **DiscoSNP++**: **confirmed** (via literature survey, `docs/HET_INDEL_SOTA.md`)
  to be the only other general-purpose reference-free caller handling
  heterozygous indels from a single diploid sample with zero reference —
  the right primary comparator, not a weak strawman.
- **Kmer2SNP**: reference-free, SNV-capable, included as a second baseline
  even though it loses badly on recall — honest inclusion, not
  cherry-picking the stronger comparator alone. **Confirmed** no indel or
  multi-allelic model, by construction (its own documented design).
- **DiscoSNP++'s 0/26 on multi-allelic sites (T2.4)**: explained by its
  bubble-based output format being inherently biallelic. **This is
  INFERRED from documented design and observed output (0 multi-allelic
  records across 3,989 total emitted), NOT verified by reading DiscoSNP++'s
  own source.** State this distinction if citing the explanation, not just
  the 0/26 result. See `claim2/mechanism_insight_claim2.md`.
- **eBWT2SNP's precision figure — self-correction, verified applied.**
  This project once cited eBWT2SNP's "99.13% precision" as motivation for
  an abandoned design direction (read-support clustering). **Confirmed**
  (`docs/EBWT_LITERATURE_CORRECTION.md`, reading the primary source —
  Prezza, Pisanti, Sciortino, Rosone, *Algorithms Mol Biol* 2019,
  PMC6364478): that figure is from a **simulated**-data experiment; on the
  paper's one real-data experiment (chr1, real reads) the ordering flips
  and **DiscoSNP++ is more precise than eBWT2SNP** (74.57% vs 66.62%). The
  positional-clustering channel eBWT2SNP inspired was kept (it measurably
  improved indel recall), but the precision-superiority argument was
  dropped. **Verified applied**: the current manuscript
  (`capsul_paper/copy.tex` line 911) already states this honestly —
  "eBWT2SNP's own published evaluation uses two simulated chromosomes and
  one real one" — so this is not an outstanding action item, only a fact
  worth keeping visible for provenance.

## Assembly comparator — Claim 3, T3.1

- **SPAdes**: **confirmed** correct comparator for 5 of 6 T3.1 datasets;
  the 6th (ERR5181310/SARS-CoV-2, amplicon sequencing, >1000x PCR-duplicated
  coverage) is explicitly outside SPAdes's own documented scope too (read
  directly from the SPAdes paper this session's earlier work) — flagged,
  not silently included as a fair comparison.
- **Result, honestly mixed**: speed dominant (129-784x). Correctness: wins
  genome fraction and indel rate, loses mismatch rate/duplication/N50 —
  explained by SPAdes optimizing for assembly correctness (dedicated
  repeat resolution + error correction) while this system optimizes for
  compressed size, both genuinely "assembly" in the mechanical sense
  (**confirmed** against `paper/METHODS.md`'s own wording), just different
  objective functions. See `claim3/t31_export_claim3.md` and
  `claim3/mechanism_insight_claim3.md`.

## Alignment comparator — Claim 3, T3.2

- **bwa+samtools+mosdepth**: a standard, uncontroversial conventional
  pipeline for the same output (per-base coverage) — clean comparison, no
  caveats needed, dominant win (16-54x, identical numbers not approximated).

## Two named prior-art lines, verified as already correctly cited in the manuscript

- **Compressive genomics (CaBLAST/CaBLAT — Loh, Baym, Berger, *Nature
  Biotechnology* 2012)**: **confirmed** cited (`gcapsul_manuscript.tex` line
  118) as the paradigm this project belongs to (computing directly on a
  compressed representation rather than decompressing first) — but
  correctly distinguished as solving a *different* problem: accelerating
  sequence *search* across many reference genome/protein *databases*, not
  compressing and serving variant evidence from *one sample's reads*.
- **Haplotig purging (Purge Haplotigs, Roach et al. 2018; Redundans,
  Pryszcz & Gabaldón 2016)**: **confirmed** cited (line 122) as the
  established assembly-literature precedent for "heterozygous regions
  assemble into separate contigs" — the phenomenon Claim 2's central
  mechanism depends on is not a novel *observation*, only the response to
  it is inverted. Both tools *remove* the duplicate haplotig to recover a
  single pseudo-haploid reference, because redundancy is noise for an
  assembly whose goal is one representative genome. For a **lossless
  archive** that remedy is inadmissible — deleting a haplotig deletes an
  allele, a variant, and losslessness itself. The correct response in a
  retained-archive setting is the opposite: keep every representative and
  resolve them by content, not purge them. This inversion, not the
  underlying observation, is what is novel.

## BFQzip — the single closest prior work by intent, verified as already the manuscript's own centerpiece comparison

**Confirmed** (`docs/NOVELTY_SCAN_DBG_CHANNEL.md`, flagged there as "the
single most important citation we were missing," and confirmed already
resolved — `guerrini2023bfqzip` is cited in the manuscript with a
dedicated comparison table, `gcapsul_manuscript.tex` line 336). BFQzip
(Guerrini, Louza, Rosone, arXiv:2304.08534) modifies bases *and* quality
scores of a FASTQ archive built on the extended BWT and positional
clustering, explicitly to preserve variant-calling information for a
downstream caller — **the same shape of claim this project makes**: one
structure serving compression and variant evidence together. The
distinctions are load-bearing, not cosmetic:

| | BFQzip | G_CAPSUL |
|---|---|---|
| compression | **lossy** (alters bases and quality scores) | **lossless**, verified by decode-and-diff |
| relationship to calling | *preserves* signal for an external caller | *performs* the call itself, in-process |
| output | a compressed file | archive + VCF + export/coverage/query |
| structure | EBWT + positional clustering | pseudogenome + read placements + k-mer counts |

"Preserves signal for someone else's caller" and "performs the calling
itself, losslessly" are different claims answering the same underlying
question. Since BFQzip is the one paper that poses this exact question in
this exact shape, a reviewer who knows it and does not see it cited would
reasonably assume it was missed — it is not; it is the manuscript's
own centerpiece related-work comparison.

## Content-retrieval / exact-match comparators — Claim 3, T3.4

- **BEETL-fastq (2014), CIndex (2022), sFASTQ (2022)**: **confirmed** (via
  literature survey, `BTR_NOTES.md`'s "Novelty, corrected after a literature
  survey 2026-09-16" section) to already do content retrieval from
  compressed reads via BWT/FM-index, natively, no auxiliary index needed —
  content retrieval itself is **not novel**. What is distinct: retrieval by
  **position** rather than exact match, which none of the three can do,
  because their data model (independent strings for BWT, k-mer adjacency
  for de Bruijn graphs) has no coordinate system. This root-cause
  explanation was checked against how BWT/BEETL actually operates this
  session, not asserted from the technique's name.
- **No speed or index-size benchmark against any of the three has been
  run.** Stated as the single most valuable missing experiment
  (`BTR_NOTES.md`, `claim3/t34_exact_match_claim3.md`).
- **Could BEETL add position retrieval with a small addition?** **Confirmed
  structurally, by tracing what each compression method actually
  computes**: this project's the completion index addition was cheap because placement
  data already existed as a free byproduct of overlap chaining. BWT's
  compression step (sorting) computes nothing analogous — position is not
  hard to look up for a BWT, it is simply never computed. Adding it would
  require a second, separate construction process, not a small index. See
  `claim3/mechanism_insight_claim3.md`.

## The PgRC2 "could they add locus retrieval with an evening of work" question — the deepest-checked comparison in the whole project

**Confirmed by reading PgRC2's actual source** (`pgrc/pgrc-decoder.cpp`):
its decoder computes per-read position and strand during decompression, but
exposes no query/extraction interface anywhere.

**Not confirmed, explicitly flagged open**: whether PgRC2 tracks per-read
mismatches, and whether its position data persists to disk or is transient.

**The argument that a quick patch would not suffice** rests on two
independent legs: (1) a freestanding, general argument that any compressor
meeting three preconditions (shared reference, position+mismatch encoding,
size-minimization) is pushed toward the same allele-splitting behavior —
applied to PgRC2's own confirmed architecture as a **prediction, explicitly
labeled as such, not a measurement** (nobody has built and run PgRC2 to
check); (2) an independent argument needing no assumption about PgRC2 at
all — this project already possessed position/strand/deviation data from
day one of its own retrieval work, and building a *correct* query layer on
top still took a multi-day effort and surfaced several silent bugs,
demonstrating that possessing the raw data is not the hard part even for a
team with full access to its own code. Full analysis:
`docs/PGRC2_EXTENSION_ANALYSIS.md`, synthesized in
`claim3/mechanism_insight_claim3.md`.

## "No tool combines all three" — verified as a checked code-audit claim, not an assumption

Verified against root-level `RECAP.md` §2 (checked this session, not
previously in `refer_paper_docs/`). The "one archive, three properties"
claim was tested directly against every real candidate architecture, not
merely asserted:

| tool | compresses losslessly | calls variants reference-free | queries a locus reference-free | one representation |
|---|---|---|---|---|
| SPRING / PgRC / PgRC2 / Genozip / NanoSpring | yes | no | no | n/a |
| DiscoSNP++ / Kmer2SNP | no (raw FASTQ in) | yes | no | n/a |
| CRAM | yes | needs external reference | needs external reference | no |
| BEETL-fastq | yes | **no** — hands off to external BWA + external samtools + an external reference | content match only, a different object (median 38% of what it returns does not even contain the query string) | **no** — three separate tools glued by a script |
| population BWT (*Genome Research* 2017) | **no** — discards quality scores by design, stated in its own text | only at *pre-specified* sites, needing an external truth VCF or genotyping array to know what a variant even is | content (k-mer) only | **no** — needs a separate 4.75 TB RocksDB store plus external Cortex graphs plus an external reference for validation |
| **G_CAPSUL** | yes, 57/57 verified | yes, from the archive alone | yes, from the archive alone | yes — one archive file |

**The "one representation" claim was independently code-audited, not
inferred from the architecture diagram**: a direct source scan (this
session's verification, of `src/encoder.cpp`, `include/caps_caller.h`,
and `src/decoder.cpp` — the encoder, caller, and all of
export/coverage/query) found **zero `system()`/`popen()`/`exec()`/network
calls anywhere** in any of the three. The only `fork()` calls (3, in the
encoder) each continue running the binary's own compiled code under a
bounded thread budget — internal parallelism, not shelling out to an
external tool. Zero external reference genome files are ever opened by the
core pipeline. What is openly borrowed and does not weaken this claim:
general-purpose entropy coders (LZMA, PPMd7, FSE/Huf0 — the same category
every real competitor in this field uses, compiled in, not shelled out to)
and vendored `fqzcomp` for quality-score compression specifically, which
feeds Claim 1 only — confirmed zero calls into it from either the caller
or the export/coverage/query code paths.

## A claim explicitly investigated and left unresolved, on purpose

Verified against `RECAP.md` §3. "India's first reference-free FASTQ
compressor" was checked, not casually asserted, and the check did **not**
resolve cleanly: FQC (Dutta, Haque, Bose, Reddy, Mande — TCS Innovation
Labs, Pune, 2015, *J. Bioinformatics and Computational Biology*) is a real,
published, India-built FASTQ compressor that one search-derived
classification described as reference-free — but its full text was
paywalled everywhere checked (PubMed, ResearchGate, Semantic Scholar,
WorldScientific, no free PDF anywhere on Google Scholar's cluster for it),
so **it was never actually read**, and a later search summary directly
contradicted the "reference-free" classification (apparently contaminated
by a description of an unrelated 2023 paper found in the same search).
**Status: genuinely unresolved, and left that way rather than asserted
either direction without reading the source.** What *is* checked and holds,
because it is a narrower compound claim (lossless compression + reference-
free variant calling + reference-free locus retrieval, together, as one
mechanism) rather than a crowded single-property claim: no tool found
anywhere combines all three. That is the defensible framing — stated
"to the best of our knowledge," not as an absolute.

## What is explicitly NOT claimed anywhere in this project

- No speed or size comparison against BEETL/CIndex/sFASTQ.
- No claim that DiscoSNP++'s multi-allelic limitation was verified from its
  own source code (inferred from documented design + observed behavior
  only).
- No claim that PgRC2's pseudogenome actually exhibits allele-splitting on
  real data (the argument is a structural prediction from its confirmed
  architecture, not an empirical finding — PgRC2 has not been built and
  run locally for this check).
