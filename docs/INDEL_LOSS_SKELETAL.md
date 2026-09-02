# Where, why and how we lose indels — per-variant skeletal analysis

Written 2026-09-02. Not a general account: every truth indel in HG002 r2 is
partitioned by which tool found it, and the events we uniquely miss are listed
individually with their sequence context.

## 1. The partition — 66 truth het-indels, same reads, same scoring

| | count |
|---|---|
| found by BOTH | 17 |
| found by **CAPSULE only** | **21** |
| found by **DiscoSNP++ only** | **9** |
| found by neither | 19 |

| | TP | FP | FN | P | R | F1 |
|---|---|---|---|---|---|---|
| CAPSULE | **38** | 12 | 28 | 0.760 | **0.576** | **0.655** |
| DiscoSNP++ | 26 | 14 | 40 | 0.650 | 0.394 | 0.491 |

**On this window CAPSULE detects 46% more true indels than DiscoSNP++ (38 vs
26) and has FEWER false positives (12 vs 14).** The aggregate 8-window loss
(0.631 vs 0.663) is therefore NOT uniform — it is driven by other windows, and
DiscoSNP++'s advantage is narrower than the average suggests.

## 2. Exactly what DiscoSNP++ catches that we do not

All 9 events, with 16 bp of reference context:

| position | event | context | longest homopolymer |
|---|---|---|---|
| 3042034 | T>TA | `ACACTGAT` **`AAAAAAAA`** | 8 |
| 3086683 | TG>T | `GCCCCCATG` **`GGGGGGG`** | 8 |
| 3114402 | TG>T | **`TTTTTTTT`**`GGGGGGGG` | 8 |
| 3163221 | C>CT | `GTAGGTCC` **`TTTTTTTT`** | 8 |
| 3169541 | A>AT | `ATCATAGAT` **`TTTTTTC`** | 7 |
| 3314333 | G>GT | `ATTAAAGG` **`TTTTTTTT`** | 8 |
| 3349776 | C>CA | `TCTGTCTC` **`AAAAAAAA`** | 8 |
| 3364448 | C>CT | `CTGGATTC` **`TTTTTTTT`** | 8 |
| 3371018 | CT>C | `TCTCTCACT` **`TTTTTTT`** | 8 |

**Every one is a 1 bp indel inside a homopolymer run of 7–8 identical bases.**
6 insertions, 3 deletions, zero exceptions.

**Why they win there, mechanically.** A 1 bp change inside an 8-mer homopolymer
shifts the run from 8 to 7 or 9. In our contig-pair bubble the two haplotypes
differ only in run LENGTH, so walking left from the anchor the sequences stay
identical for the whole run and the divergence point is ambiguous by up to 8
positions — `extract_bubble` cannot localise the event, and the closing-anchor
test then fails because both paths re-converge at slightly different offsets.
DiscoSNP++'s BFS extends one graph path base by base and accepts the SHORTEST
closing extension, which resolves a length change directly; it also carries an
explicit ambiguity allowance (`max_ambigous_indel`, default 20).

## 3. What we catch that they do not

| event class | CAPSULE-only |
|---|---|
| 1 bp | 7 (3 DEL, 4 INS) |
| 2 bp | 3 INS |
| 3 bp | 1 INS |
| 4 bp | 4 INS |
| ≥5 bp | 6 (2 DEL, 4 INS) |

**14 of our 21 exclusive finds are ≥2 bp.** The two tools are complementary by
event length: they own 1 bp homopolymer events, we own everything longer.

## 4. Is this class benchmarked elsewhere? — yes, and it is a known hard class

GIAB's own stratification resource (Nat Commun 2024) defines "challenging
sequence" as the union of *all tandem repeats, all homopolymers >6 bp, all
imperfect homopolymers >10 bp, and difficult-to-map regions*, and reports that
**homopolymers score lower than low-mappability or GC-extreme regions** — i.e.
they are the hardest stratum in the standard benchmark. GIAB now excludes
homopolymers >10 bp from difficult-to-map benchmark regions because even
PCR-free short reads are unreliable there.

Our 9 missed events sit in 7–8 bp homopolymers: inside GIAB's "challenging"
definition, just under the length at which the benchmark itself stops trusting
short reads.

## 5. Consequences

1. **The loss is one specific, named class**: 1 bp indels in ≥7 bp
   homopolymers. Not "precision", not "chimeric bubbles", not a threshold.
2. **It is fixable in principle** by length-aware bubble extraction — comparing
   homopolymer RUN LENGTHS between haplotypes instead of walking for a
   divergence point — which is a change to `extract_bubble`, not a substrate
   change. This is the first indel lead in this project with a named mechanism
   behind it rather than a refuted guess.
3. **The complementarity is a publishable result in itself**: on this window we
   find 46% more true indels, they find a class we cannot, and the union of the
   two tools would exceed either.
