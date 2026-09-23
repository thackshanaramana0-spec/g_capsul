---
Date: 2026-09-19
Title: Paper Overview — What This Paper Actually Is
Purpose: The single top-level entry point for understanding the whole
  paper before diving into any claim-specific folder. States the real
  title/name structure, the actual manuscript section layout (verified
  against the live tex file), and how the three claims fit together.
When to refer to this file: First file to read when starting any work on
  this paper; onboarding; deciding where something belongs (which claim
  folder, or here); before writing an abstract, cover letter, or any
  paper-wide summary.
Keywords: overview, G_CAPSUL, ARCH, COMPACT, FAITHFUL, ADDRESSABLE, paper
  structure, PLOS, title, backronym
---

# What this paper actually is

## The name structure — two names, one is the headline, one is secondary

**G_CAPSUL** (Compact, Addressable, Pseudogenome-Structured, Unified
Lossless) is the tool's name and the paper's actual title/contribution. It
is what appears in the title and what the abstract introduces as the thing
being evaluated.

**ARCH** (Construct, Retain, Compress, Harness — exact wording still being
finalized, see `refer_paper_docs/claim3/mechanism_insight_claim3.md` for
the fullest current draft) is a proposed name for the underlying design
pattern G_CAPSUL implements. **As of this session, ARCH has NOT been added
to the live manuscript** (`capsul_paper/capsul_manuscript.tex`) — it exists
only in draft discussion and in this `refer_paper_docs/` folder set. Do not
treat ARCH as already-published terminology; it is a proposed addition,
pending the user's decision on final wording and placement.

**Decided placement, if adopted** (per session discussion): name ARCH once,
folded into the sentence introducing G_CAPSUL, not given a separate
"we present ARCH" sentence of its own — verified against real precedent
(Transformer, ResNet, MapReduce, PgRC2 abstracts) that winning papers use
one headline name, with any secondary mechanism named in passing inside
the sentence about the main contribution, not launched separately. See
`central_insight_and_mechanism.md` for the reasoning behind this pattern.

## The paper's actual structure, verified against the live manuscript

Confirmed directly from `capsul_paper/capsul_manuscript.tex` this session
(not assumed from a template):

```
Abstract
Author summary
Introduction
Results
  COMPACT       (Claim 1 — T1.1, T1.2)
  FAITHFUL      (Claim 2 — T2.1-T2.5)
  ADDRESSABLE   (Claim 3 — T3.1-T3.5)
Discussion
Materials and methods
  Datasets
  Method overview
  Pseudogenome construction and read placement
  Lossless FASTQ encoding and reconstruction
  Heterozygous structure reconciliation
  Reference-free variant calling
  Pseudogenome export, coverage, and positional retrieval
  Experimental setup and baselines
  Evaluation metrics
  Implementation and reproducibility
Conclusion
Supporting information
Author contributions
Acknowledgments
```

This is a PLOS-format paper (structured abstract conceptually divided into
background/methodology/conclusions per PLOS guidelines, separate
plain-language Author Summary — both confirmed against PLOS Computational
Biology's actual submission guidelines this session, not assumed).

## The three claims, restated at the top level

| Claim | Section name | What it tests | Folder |
|---|---|---|---|
| 1 | COMPACT | Is the archive smaller than general-purpose FASTQ compressors, losslessly? | `refer_paper_docs/claim1/` |
| 2 | FAITHFUL | Does the retained structure still support reference-free variant calling? | `refer_paper_docs/claim2/` |
| 3 | ADDRESSABLE | Can the same archive serve export, coverage, range query, exact match, and locus retrieval? | `refer_paper_docs/claim3/` |

**The three are not independent achievements stacked side by side.** Claim
1 builds the structure; Claims 2 and 3 are two different, separately-
measured demonstrations that the structure is reusable. This is stated
directly in the paper's own Introduction (line ~336: *"Three aspects of
that question were examined separately, corresponding to the three claims
introduced above: COMPACT, FAITHFUL, and ADDRESSABLE"*) and in its own
Conclusion (line ~1303: *"a compression-derived representation can remain
useful as a computational substrate rather than serving only as an
intermediate for reconstruction"*) — full argument in
`central_insight_and_mechanism.md`.

## What "new/main" vs "old/waste" means in this repo, for navigation purposes

`refer_paper_docs/` contains 111 files as of this session — the large majority are
historical working notes, many explicitly marked SUPERSEDED or WITHDRAWN
in their own text (per this project's own stated rule: retractions stay in
place, never deleted). **The authoritative, current sources are**:
`CLAUDE.md` (top-level, always the first thing to read), the three
`refer_paper_docs/CLAIM{1,2,3}_FINAL_VERDICT.md` files (each explicitly marks its own
superseded sections), `results/*.csv` (the actual numeric ground
truth), and this `refer_paper_docs/` folder tree (built and verified this
session specifically to be the current, navigable synthesis layer). When in
doubt about whether a `refer_paper_docs/` file is current, check whether `CLAUDE.md`
or a `refer_paper_docs/` file cites it as a live source — if not, treat it
as historical record, not something to draw new numbers from.

## File index — where to go next

- `coverage_manifest.md` — the mechanical proof: every one of `docs/`'s 108
  files, its disposition, and exactly where its content landed (or why it
  didn't need to)
- `technical_architecture.md` — the container format and the assembly
  pipeline itself (chaining, pigeonhole mapping, MAXMAP ramp, second
  region, MEM self-match) — HOW the archive is built, independent of any
  one claim's results
- `central_insight_and_mechanism.md` — the one mechanism tying all three
  claims together
- `novelty_and_prior_art.md` — honest, cross-claim comparison to PgRC2,
  BEETL/CIndex/sFASTQ, DiscoSNP++/Kmer2SNP, SPAdes/bwa
- `honest_limitations_and_scope.md` — every real trade-off, withdrawn
  number, and scope limit, in one place
- `numbers_and_verification_index.md` — master index of every headline
  number with source and verification status
- `claim1/`, `claim2/`, `claim3/` — per-claim detail folders
