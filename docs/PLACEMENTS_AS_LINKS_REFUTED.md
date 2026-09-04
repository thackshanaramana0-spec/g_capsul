# The compressor's read placements are NOT de Bruijn graph links

Built, instrumented, measured twice, refuted. This closes the "free linked
de Bruijn graph" idea, which was the strongest remaining headroom hypothesis.

## The hypothesis

To compress, the encoder maps every read onto the pseudogenome and records
`ppos[u]` — a position, which is how a read gets stored as a coordinate instead
of a sequence. `CallData` already carries the derived `read_cid`, `read_pos`,
`read_rc`, and Method B was discarding all of it.

That looks exactly like the read-to-path threading other tools pay for:
McCortex builds a Linked de Bruijn Graph in a dedicated pass costing ~20 GiB of
links on top of a 50 GiB graph (Turner et al., *Bioinformatics* 2018), and
LueVari stores read provenance as colours in succinct structures (Bioinformatics
2020). If our placements were equivalent, we would get for free the structure
that is the field's standard answer to repeat resolution — and repeat collapse
is our measured false-positive mechanism (false sites carry 12x the allele depth
of true ones), and `find_sb` exhausts on 960,541 of 1,225,194 branching nodes
precisely because it cannot resolve repeats.

## What was built

`read_ppos` (global pseudogenome offset) was added to `CallData` and filled from
`ppos[u]` in the encoder; the placements are carried into Method B instead of
dropped; and the read-coherence sweep records, per bubble, where the compressor
placed every supporting read. Cost is three vector copies. All of it stays in.

## Measurement 1 — contig-level: do the two alleles share a locus?

A real heterozygous site is one locus with two haplotypes, so both alleles'
reads should be placed together; a collapse pools two copies. Splitting the
bubbles and scoring each half on HG002 r2:

| group | TP | FP | precision |
|---|---|---|---|
| alleles share >=1 contig | 76 | 7 | 0.916 |
| alleles disjoint | 260 | 13 | **0.952** |
| baseline (all) | 336 | 20 | 0.944 |

The split is complete (76+260 = 336, 7+13 = 20). "Disjoint" is *slightly
better*, the opposite of the hypothesis, and far too weak to act on — rejecting
the shared group would cost 76 true positives to remove 7 false ones.

Diagnosis: contig id is a poor proxy for locus. The substrate is fragmented
(451,578 contigs at full chr20), so one locus routinely spans several.

## Measurement 2 — global coordinates: are supporting reads near each other?

Distance needs one coordinate system, so this repeats the test on `read_ppos`.
Span distribution over 737 bubbles: `<=200:42 <=1k:7 <=10k:165 <=100k:68
<=1M:375 >1M:21 none:59`.

| group | TP | FP | precision |
|---|---|---|---|
| span <= 10 kb | 73 | 4 | 0.948 |
| span > 10 kb | 267 | 16 | 0.944 |
| baseline | 336 | 20 | 0.944 |

**Precision is identical in both groups.** The signal is not weak; it is absent.

## Why it fails — the reason is structural, not a bug

A pseudogenome is built by greedy overlap chaining, and its coordinate system is
optimised for **compressibility**, not for genomic fidelity. Chaining
deliberately **merges** near-identical sequence — that merging is precisely what
makes the archive small. Repeat copies are the most compressible thing in a
genome, so they are exactly what chaining collapses first.

So a read from repeat copy A and a read from copy B are *placed together on
purpose*. The placement records where a read was **stored**, not where it came
**from**. The repeat information was not left unexploited; it was **spent** to
achieve compression.

This is the sharp distinction that the analogy hid:

| | records | survives repeats? |
|---|---|---|
| McCortex links | which path a read took **through the graph** | yes — that is their purpose |
| LueVari colours | which **read** each k-mer came from | yes |
| our `ppos` | where a read was **stored in a structure that merges identical sequence** | **no — collapsed by design** |

The three are not interchangeable, and the difference is not a matter of
resolution or coordinates. It is that one is a record of provenance and the
other is a record of storage.

## Consequence

- The "free linked de Bruijn graph" headroom claim is **withdrawn**. It should
  not appear in the paper.
- Getting the equivalent signal would mean retaining read-to-path provenance
  *through* chaining — i.e. NOT collapsing repeat copies — which directly
  opposes what makes Claim 1 win. That is a real trade-off, not an oversight,
  and it is the honest reason the free lunch does not exist.
- What is kept: `read_ppos` plumbing and the placement instrumentation, both
  inert by default (`CAPS_DBG_PLACEMODE=0`), so the measurement is reproducible.
- Nothing shipped changes: kc nodes 1,063,607 and SNV F1 0.886 before and after.

## The general lesson, since this session repeated it

Both formulations were argued from a mechanism that sounded compelling and then
measured to carry no signal — as did nine indel filters before them. The
instrumentation-first order (count, split, score, and only then filter) is what
kept this from shipping. It cost two runs to refute an idea that would otherwise
have gone into a paper as a headline claim.
