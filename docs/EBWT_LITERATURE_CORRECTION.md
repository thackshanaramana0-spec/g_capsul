# eBWT2SNP, ebwt2InDel, SKA lo — what the papers actually say

Read 2026-09-04, primary sources, after the full-chr20 run and the
`PG_AS_GRAPH_REFUTED.md` finding. Two of the three change something in this
project: one is a **correction to our own documentation**, one is **prior art
we had not identified and must cite**.

---

## 1. CORRECTION: our eBWT2SNP precision figure is a SIMULATED-data number

This project cites "99.13% precision" for eBWT2SNP in at least five places
(`docs/HET_INDEL_SOTA.md:19`, `docs/INDEL_PRECISION_ROOT_CAUSE.md:18,86,136`,
`include/caps_caller.h:805-806,2076`), and uses it to justify a design
direction — `caps_caller.h:805` reads *"It reports 99.13% precision against
DiscoSNP++'s 77.80% on chr22, and that guarantee is why."*

The chr22 attribution is right. **The data is simulated.** From the paper
(Prezza, Pisanti, Sciortino, Rosone, *Algorithms Mol Biol* 2019, PMC6364478):

| experiment | data | chr | cov | ebwt2snp sens/prec | DiscoSNP++ sens/prec |
|---|---|---|---|---|---|
| 1 | **SIMULATED** | 22 | 29x | 88.47% / **99.13%** | 80.69% / 77.80% |
| 2 | **SIMULATED** | 16 | 22x | 81.29% / 99.64% | 80.85% / 68.40% |
| 3 | **REAL** | 1 | 43-47x | 77.94% / **66.62%** | 71.94% / **74.57%** |

**On real data the ordering flips: DiscoSNP++ is MORE precise than eBWT2SNP
(74.57% vs 66.62%).** The 99.13%/77.80% pairing we quote is real in the sense
that both numbers come from the same table — but it is the simulated table, and
the conclusion it supports does not survive on real reads.

**Why this matters here, concretely:**

1. It removes the main stated reason for pursuing an eBWT-style read-support
   guarantee. `INDEL_PRECISION_ROOT_CAUSE.md:86` argues eBWT2SNP is attractive
   because it has "(99.13% precision), needs no assembly at all." At 66.62% on
   real data that argument is gone.
2. It retroactively **supports** two decisions already taken on measurement
   alone: disabling `ridx`/`read_support` by default (measured byte-identical,
   i.e. vacuous for us) and abandoning the read-clustering probes (best fair
   result 29.9% recall / 36.2% precision). We were not failing to reproduce a
   99% mechanism; the mechanism is not a 99% mechanism on real data.
3. Our own **measured** full-chr20 indel precision is **0.839**, and
   DiscoSNP++'s is 0.923 — both above eBWT2SNP's real-data 66.62% (different
   variant class and truth set, so not directly comparable, but it rules out
   "adopt eBWT clustering to fix precision" as an obviously-winning move).

Also correct two secondary claims in our docs: the paper reports clustering at
~3 MB RAM but **SNP calling at 12 GB** and ~3 hours total on real chr1 — not
the "415 MB" this project's notes attribute to it.

## 2. The one mechanism that is structurally compatible with what we have

Last night's `PG_AS_GRAPH_REFUTED.md` established that our pseudogenome is
linearised by construction — 97.7% of anchors have de Bruijn out-degree 1 — so
any bubble/branch-based method (DiscoSNP++, SKA lo) cannot be ported onto our
contigs.

**eBWT positional clustering is the exception, and it is worth stating why.**
It does not use branches at all. Its Theorem 3.3 locates clusters as local
minima of the LCP array: bases covering the same genome position appear
contiguously in the eBWT of the reads, and the cluster boundary is found from
the LCP profile rather than from a fixed k. The paper's own emphasis is that
this is **independent of k** — it adapts to how unique each genomic locus
actually is.

That is a real difference from what this project tried. Our read-clustering
probes used a **fixed k-mer right-context anchor**. The paper's entire
contribution is *not* fixing k. So our probe's failure (29.9%/36.2%) is not
evidence against the paper's mechanism — it is evidence against a k-fixed
approximation of it. That distinction should be recorded honestly rather than
letting the probe stand as a refutation of the published method.

**It is still not worth building, for reasons that are about cost and payoff,
not about the idea being wrong:**

- It requires the eBWT + LCP + GSA of the read set. We build none of these; our
  compressor builds a pseudogenome. Constructing them for 12.6M reads is a
  separate pipeline — precisely the "combining two projects" the project has
  already rejected as a framing.
- Its real-data precision is 66.62%, below what we already achieve.
- Indel support is explicitly future work in this paper ("indels behave exactly
  like SNPs in Theorems 3.2 and 3.3" but typing them "remains future work
  requiring local alignment"); that is what `ebwt2InDel` was written for.
- **We currently win het-indel at full scale anyway** (F1 0.592 vs 0.576).

## 3. PRIOR ART WE HAD NOT IDENTIFIED — and it is close to our thesis

> Guerrini, Louza, Rosone et al., *"Lossy Compressor preserving variant calling
> through Extended BWT"*, arXiv:2304.08534.

From the abstract: it lossily compresses FASTQ **reference-free**, modifying
bases *and* quality scores together, "based on the Extended Burrows-Wheeler
Transform (EBWT) and positional clustering," and shows it "is able to achieve
good compression while preserving information relating to variant calling more
than the competitors."

**That is one shared structure serving both compression and variant calling,
reference-free — the same shape as this project's central claim.** It is the
nearest prior art to our one-pipeline story and it is not cited anywhere in our
docs. It must be, and the paper must distinguish from it explicitly.

**The distinction is real and strongly in our favour, which is why stating it
plainly costs nothing:**

| | arXiv:2304.08534 | G_CAPSUL |
|---|---|---|
| compression | **LOSSY** — alters bases and quality scores | **LOSSLESS**, verified by decode-and-diff against the original FASTQ |
| what is preserved | "information relating to variant calling" | the entire file, bit for bit |
| shared structure | EBWT + positional clustering | pseudogenome + read placements |
| also serves | variant calling | variant calling **+ export/coverage/query** (Claim 3) |

A lossy compressor that preserves *enough* for variant calling is a genuinely
different claim from a lossless compressor that *also* calls variants from the
structure it built to compress. But a reviewer who knows this paper and does
not see it cited will assume we did not.

## 4. SKA lo — does not apply, recorded so it is not revisited

> *Mol Biol Evol* 42(4):msaf077, 2025.

Colored de Bruijn graph over split k-mers; its novelty is **variant groups**
(all variant combinations between an entry and an exit node) rather than
individual bubbles, which is why it handles dense/overlapping variants better.
Indels come from variant groups where one path has at most 2(k-1) nt. Reports
~96% sensitivity, zero false-positive indels, 1.7 GB / 14 s on 55 *S. aureus*
samples.

**Not a competitor for T3/T5 and not portable to us**, for two independent
reasons: it is a **multi-sample, within-strain pathogen** tool (colors =
samples), not a single-diploid-individual caller; and it is branch-based, so
§2's linearisation finding applies to it exactly as it does to DiscoSNP++.

---

## Actions taken

- `docs/HET_INDEL_SOTA.md` and `include/caps_caller.h` comments corrected to
  mark 99.13% as simulated and to add the real-data row.
- arXiv:2304.08534 added to the prior-art comparison.
- No code behaviour changed by any of this.

Sources: [Prezza et al., Algorithms Mol Biol 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6364478/) ·
[Lossy Compressor preserving variant calling through EBWT](https://arxiv.org/abs/2304.08534) ·
[SKA lo, MBE 2025](https://academic.oup.com/mbe/article/42/4/msaf077/8103706) ·
[ebwt2snp/ebwt2InDel](https://github.com/nicolaprezza/ebwt2snp)
