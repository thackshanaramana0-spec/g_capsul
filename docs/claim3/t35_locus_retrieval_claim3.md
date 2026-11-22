# T3.5 — Locus Retrieval

> Record of the T3.5 mechanism, the multi-individual and multi-locus result, the one diagnosed miss, and the cost of reaching a clean sweep.

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
measured on the archive's own data — see the Claim 3 locus-addressability write-up
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
orchestration layer (`scripts/claim3_scripts/run_locus_retrieval_capsule.sh`,
`scripts/claim3_scripts/score_locus_retrieval.py`), by calling T3.4's `query` primitive **twice**
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

## Results — four individuals, one locus

Source: `results/claim3/claim3_T3.5_locus_retrieval.csv`, produced by phase 4
of `scripts/benchmark_scripts/benchmark_1_run.sh`. Detail of the earlier
replication: `docs/claim3/t34_t35_multi_individual_claim3.md`.

A heterozygous SNV is tested only if both flanking 40-bp reference sequences
are unambiguous and the upstream flank occurs in the pseudogenome, the same
rule for every individual, so the number of sites differs between them.

| Individual | Locus | Sites | Position query (tolerance 2) | Archive + FM-index (tolerance 3) |
|---|---|---:|---:|---:|
| HG002 | chr20:3.0-3.6 Mb | 400 | 400/400 | 400/400 |
| HG003 | chr20:3.0-3.6 Mb | 335 | 335/335 | 335/335 |
| HG004 | chr20:3.0-3.6 Mb | 400 | **399/400** (one real miss, diagnosed below) | **400/400** |
| HG005 | chr20:3.0-3.6 Mb | 317 | 317/317 | 317/317 |
| **All** | chr20:3.0-3.6 Mb | **1,452** | 1,451/1,452 (99.93%) | **1,452/1,452 (100%)** |

"Position query" is the bilateral query on the decoded pseudogenome and read
placements. It is the "archive" column of the result CSV, together with its
homozygous-control hits (29, 35, 29 and 66 of 400). Phase 4 of
`benchmark_1_run.sh` re-runs both columns from the downloaded data.

A separate additional check on HG005 at chr20:4.0-4.6 Mb (400/400, measured
with the earlier k-mer completion index) is not part of the 1,452 aggregate
and is not in the manuscript's table.

## The one position-query miss — fully diagnosed, not hand-waved

Site chr20:3,009,763 (HG004, T→G). GIAB itself flags this exact site as a
homopolymer/simple-repeat difficult region. Both probes retrieved reads (62
total, correctly anchored) but all were REF; the ALT-carrying reads were
never at that pseudogenome locus, because their difference from the
consensus was an **indel** (the dominant error mode in homopolymers), not a
substitution — and this archive's deviation model is substitution-only. A
read differing by indel exceeds the placement mismatch ceiling and gets
placed *elsewhere* in the pseudogenome at encode time. No positional query
at any tolerance could find it there, because it structurally is not there.

**This was tested, not assumed**: querying by read content, which matches
what a read contains instead of where it was placed, recovered the
ALT-carrying reads that the position query could not reach (11 with the
earlier k-mer index, flipping the site to a correct het call), and the
FM-index recovers the site too (HG004 400/400). That confirms the diagnosis
and is why a content index closes exactly this class of miss.

## The real cost, disclosed — do not present 100% as free

Homozygous negative-control false positives (400 control sites per
individual) rise when the content index is added: HG002 29 to 48, HG003 35 to
52, HG004 29 to 49, HG005 66 to 87. HG005 is notably higher than the other
three. That is an open question, not a diagnosis: the earlier additional
check at a second HG005 locus gave 63 of 400, which suggests the rate is
locus-specific, not a property of the individual.

## Generalization claim — precise wording

The mechanism (bilateral union) is geometric, not a completeness argument the
way T3.4's FM-index is, which is an asymmetry, not a weakness to
hide. It generalizes on the strength of: (a) every parameter is derived from
data, never fitted (seed floor, probe construction identical across runs),
and (b) it held clean across two ancestries, two trios and a locus never
tested before. State it that way, "broad, verified empirical
generalization", not as a mathematical guarantee, which would overclaim
relative to T3.4's actual construction-backed claim.
