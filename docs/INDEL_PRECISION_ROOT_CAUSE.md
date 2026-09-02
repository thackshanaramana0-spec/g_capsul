# The indel precision gap — root cause, after reading all three competitors

Written 2026-09-02 after reading DiscoSNP++'s source, the eBWT2SNP paper
(Algorithms Mol Biol 2019), and Kmer2SNP's method, then implementing and
measuring **seven** distinct precision mechanisms against our own caller.

**Standing:** het-SNV **we win** (0.890 vs 0.874). het-indel we lose
(0.592 vs 0.663), entirely on precision (0.759 vs 0.910); our indel recall
(0.491) is close to theirs (0.522).

---

## 1. How each competitor actually achieves precision

| tool | precision mechanism | reported |
|---|---|---|
| **DiscoSNP++** | bubble leaves a shared graph node and must RE-CONVERGE on another shared node (both ends anchored in read k-mers); then `kissreads2` re-validates every path against reads | 0.910 here |
| **eBWT2SNP** | every emitted fragment of length 2k+1 (k=30) must be an actual **substring of ≥C real reads** within Hamming 2; clusters with no consensus are discarded | 99.13% (their chr22) |
| **Kmer2SNP** | models calling as a **maximum weight matching** — each heterozygous k-mer gets at most ONE partner | — |

## 2. What we implemented from each, and what it measured

All swept across five windows, all kept in-code behind flags.

| # | mechanism | source | result |
|---|---|---|---|
| 1 | Read-level junction support (k-mer counts) | kissreads2 principle | **+0.18 precision** (0.37→0.93 on one window); indel F1 ~0.40→0.58 |
| 2 | Closing anchor — bubble anchored at BOTH ends | DiscoSNP++ | **+0.03 precision**, F1 0.5930→0.5976 |
| 3 | Unique closing anchor | DiscoSNP++ | negative, 0.5976→0.5932 |
| 4 | Tolerant re-convergence flank ("indel + n SNPs") | DiscoSNP++ | negative, 0.602→0.595 |
| 5 | Full junction coherence (min over k-mers) | kissreads2 | neutral, +0.001 |
| 6 | Read-substring guarantee, 31/61/91 bp fragments | **eBWT2SNP** | neutral at 31bp, **negative** at 61/91bp |
| 7 | Maximum-weight matching | **Kmer2SNP** | neutral, 0.5976→0.5962 |
| 8 | Extended contig agreement (60–200 bp) | paralog filter | negative, 0.5930→0.5644 |

Two of eight helped. Six did not — and they failed in a *consistent pattern*.

## 3. The root cause

Group the failures by why they failed:

**Group A — context-based tests (4, 6, 8) fail because our contigs are short.**
Mean contig ~335 bp against 148 bp reads. Every test needing 60–200 bp of
flanking context either runs out of contig or rejects true events near contig
ends. Adaptive context was built specifically to rule out "truncation" as the
explanation, and it did not recover the loss.

**Group B — read-support tests (1 partly, 5, 6) are near-vacuous for us.**
Our candidates are built FROM contigs, and contigs are built FROM reads, so
every candidate is read-supported *by construction*. `kissreads2` works for
DiscoSNP++ because its paths are **hypotheses from graph traversal** — a path
can be proposed that no read traverses, and coherence throws it out. Ours are
**observed sequence**; there is nothing for the test to reject.

**Group C — the matching constraint (7) fails because our false bubbles are
uncontested.** A matching can only reject a pairing that loses to a better
one. Our false positives are typically a locus with exactly one — wrong —
partner. Nothing competes with them.

### The unifying statement

> **All three competitors read their variant candidates DIRECTLY out of the
> read data** — graph paths traversed through reads, eBWT cluster fragments
> lifted out of the read collection, k-mer pairs observed in reads.
> **We CONSTRUCT ours**, by splicing one contig's flanks onto another
> contig's allele.

A constructed candidate can be chimeric in a way that:

- read-support tests cannot detect (all its pieces are read-derived), and
- context tests cannot detect (our contigs are too short to give the span), and
- matching cannot detect (nothing competes with it).

**That is why the gap is 0.759 vs 0.910, and why it did not close under any of
the seven mechanisms.** It is a property of the candidate GENERATOR, not of
any filter — which is also why mechanism #2 (closing anchor, the one change
that made the generator more graph-like) was the only structural one that
helped.

## 4. What would actually close it

Not another filter. The generator must stop constructing candidates and start
reading them out of the data:

1. **eBWT positional clustering** — the strongest option in the literature
   (99.13% precision), needs no assembly at all, and the ARCS family already
   has FM-index/suffix-array machinery (`fm_apsp.cpp`, `sa_apsp.cpp`). Both
   haplotype fragments come out of the same cluster, so they are never
   chimeric.
2. **Longer contigs** would rescue Group A, but not Group B or C.

## 5. Honest bottom line

- **het-SNV: won and generalized** — 0.890 vs 0.874, 5 of 8 evaluations,
  largest win on a held-out window, wins on 2 of 3 unseen individuals.
- **het-indel: lost, and now understood** — 0.592 vs 0.663. The cause is
  stated structurally and backed by eight measurements rather than asserted.
  Closing it is a substrate change of known shape, not a tuning exercise.


---

## 7. CORRECTION to section 4, measured 2026-09-02 (later)

Section 4 argues the gap is closable only by generating candidates from the
read data rather than constructing them from contigs, and points at the
positional-clustering channel as the embodiment of that. **The channel-isolation
measurement partly refutes my own framing and belongs on the record.**

Running each indel channel alone (bubble channel suppressed with
`CAPS_INDEL_ANCH=9999`; pcluster suppressed with `CAPS_NO_PCLUSTER=1`):

| window | both | bubble-only | pcluster-only |
|---|---|---|---|
| HG002 r2 | 0.655 | 0.636 | **0.083** |
| HG002 r5 | **0.699** | 0.623 | **0.208** |

and on r2, where the per-channel P/R is reliable:

| channel | precision | recall |
|---|---|---|
| bubble (CONSTRUCTED candidates) | **0.795** | 0.530 |
| pcluster (candidates READ from data) | **0.500** | 0.045 |
| both | 0.760 | 0.576 |

**So the data-derived channel is LESS precise standalone, not more.** Its value
is that it is ADDITIVE — it finds indels the bubble channel misses (recall
0.530 -> 0.576 on r2, and a large 0.623 -> 0.699 on r5) — but it does not
deliver the precision that section 4 predicted for it.

What survives of the original analysis: filters could not close the gap
(ten refuted mechanisms), and adding an independent detection channel could
partly close it. What does NOT survive: the claim that reading candidates out
of the reads is intrinsically more precise *in this implementation*. Our
pcluster candidates are still anchored on, and positioned by, a contig, so they
inherit the substrate's weaknesses; eBWT2SNP's guarantee comes from both
haplotype fragments emerging from one cluster, which this does not reproduce.

The remaining gap (0.631 vs 0.663) should therefore be attributed to
**substrate quality**, not to a missing filter and not to a solved-in-principle
generator.
