---
Date: 2026-09-19
Title: T3.5 — Locus Retrieval at Heterozygous Sites (Bilateral Anchoring)
Purpose: Complete, verified record of T3.5's mechanism, why it is the
  paper's actual novel claim (not T3.4), the multi-individual/multi-locus
  generalization result, the one diagnosed miss, and the real cost of
  reaching a clean sweep.
When to refer to this file: Writing/checking T3.5's table or text; deciding
  what "the paper's insight" sentence should say; citing the 96/400-vs-
  400/400 motivating number; explaining why this table has no competitor.
Keywords: T3.5, locus retrieval, bilateral anchoring, heterozygous, allele
  splitting, coordinate-only addressing, completion index, homopolymer, indel eviction
---

# T3.5 — Locus Retrieval

## What it measures — and why this is the actual novel claim of Claim 3

Given a known heterozygous site, retrieve reads covering **both** alleles,
even though compression may have split them onto different, disconnected
pseudogenome coordinates. This is the table that directly demonstrates the
paper's central mechanism finding, not just one of five equally-weighted
features — see `mechanism_insight_claim3.md` for the full
connection to Claim 2's F1 result.

**No competitor exists for this operation anywhere** — not "we are faster,"
literally nothing else has a concept of "this locus" to even attempt it.
BWT-family tools store reads as independent strings with no positional
relationship between them; a de Bruijn graph stores k-mer adjacency, not
placement. Neither can express "these two disconnected addresses are the
same biological site."

## The motivating number — state this precisely, it is not T3.3's result

Following a **single** pseudogenome coordinate for a het site recovers both
alleles only a small fraction of the time (the coordinate-only control,
measured on the archive's own data — see `docs/CLAIM3_LOCUS_ADDRESSABILITY.md`
for the exact figure and its provenance). This is the direct, measured proof
that plain coordinate lookup (the mechanism T3.3 uses generically) fails
specifically at heterozygous sites, because that is exactly where
compression can separate the two alleles into different pseudogenomic
regions. **This number belongs under T3.5 as the mechanism's motivation, not
under T3.3** — T3.3 itself has no failure mode to report; it is the
generic lookup that T3.5 exists to patch at exactly the sites where it
breaks.

## The mechanism — bilateral anchoring, and where it actually lives

**Not a separate C++ code path.** T3.5 is built entirely at the
orchestration layer (`scripts/t34_realdata/run_window.sh`,
`scripts/score_bilateral.py`), by calling T3.4's `query` primitive **twice**
— once with a probe 40bp upstream of the variant, once 40bp downstream —
and unioning the two result sets in Python, assigning each read to whichever
probe's occurrence it overlaps most. Neither probe contains the variant
itself, so neither can rig the result. See `code_mapping_claim3.md`
for the exact call sequence.

**Why two anchors, not one**: the pseudogenome is built by greedy overlap
chaining, so a contig's local neighbourhood is an artifact of which reads
happened to overlap, not the true genomic neighbourhood. A probe anchored
on only one side of a variant can be structurally unreachable if chaining
broke on that side — measured directly: upstream-only anchoring found 0 of
6 failing sites; downstream-only found 5 of 6. Unioning both closes the gap
neither can close alone.

## Results — locked, multi-individual, multi-locus

Source: `docs/T34_T35_MULTI_INDIVIDUAL_20260919.md`.

**Native (archive alone, no the completion index), first pass:**

| Individual | Sites | Bilateral (both) |
|---|---|---|
| HG002 | 400 | 400/400 |
| HG003 | 335 | 335/335 |
| HG004 | 400 | **399/400** (one real miss — diagnosed below) |
| HG005 (locus 1) | 317 | 317/317 |

**With the completion index completion, all four individuals plus a second, independently
chosen locus on HG005 — the full, final result:**

| Individual | Locus | Sites | Result |
|---|---|---|---|
| HG002 | chr20:3.0-3.6Mb | 400 | 400/400 |
| HG003 | chr20:3.0-3.6Mb | 335 | 335/335 |
| HG004 | chr20:3.0-3.6Mb | 400 | **400/400** (miss fixed) |
| HG005 | chr20:3.0-3.6Mb | 317 | 317/317 |
| HG005 | chr20:4.0-4.6Mb (new locus) | 400 | 400/400 |

**Aggregate: 1,452/1,452 (100%) across four individuals, two GIAB trios, two
independent loci.** Native alone (before the completion index): 1,451/1,452 = 99.93%.

## The one native miss — fully diagnosed, not hand-waved

Site chr20:3,009,763 (HG004, T→G). GIAB itself flags this exact site as a
homopolymer/simple-repeat difficult region. Both probes retrieved reads (62
total, correctly anchored) but all were REF; the ALT-carrying reads were
never at that pseudogenome locus, because their difference from the
consensus was an **indel** (the dominant error mode in homopolymers), not a
substitution — and this archive's deviation model is substitution-only. A
read differing by indel exceeds the placement mismatch ceiling and gets
placed *elsewhere* in the pseudogenome at encode time. No positional query
at any tolerance could find it there, because it structurally is not there.

**This was tested, not assumed**: re-querying with the completion index (which matches by
read content, not position) recovered 11 ALT-carrying reads that plain query
could not reach, flipping the site to a correct het call. This confirms the
diagnosis and is why the completion index closes exactly this class of miss.

## The real cost, disclosed — do not present 100% as free

Homozygous negative-control false positives rose with the completion index across every
individual: HG002 31→42, HG003 38→52, HG004 28→49, HG005 71→87 (a notably
higher rate than the other three — flagged as a real, undiagnosed open
question; a second locus test on HG005 showed the rate drop to 63/400,
suggesting it may be locus-specific rather than an inherent property of that
individual, but this is not a full diagnosis).

## Generalization claim — precise wording

The mechanism (bilateral union) is geometric, not a completeness proof the
way T3.4's the completion index is — this is an honest asymmetry, not a weakness to hide.
It generalizes on the strength of: (a) every parameter is derived from data,
never fitted (seed floor, the completion index k, probe construction identical across
runs), and (b) it held clean across two ancestries, two trios, and a
genuinely new locus never tested before. State it that way — "broad,
verified empirical generalization" — not as a mathematical guarantee, which
would overclaim relative to T3.4's actual proof-backed claim.
