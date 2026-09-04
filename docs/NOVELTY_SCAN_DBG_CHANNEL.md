# Novelty scan — is the dBG-from-compressor-index channel actually new?

2026-09-04, after the dBG bubble channel measured SNV F1 0.888 -> 0.910 on
HG002 r2. Primary sources checked before claiming anything.

**Short answer: the individual pieces all exist. The specific combination does
not appear to — but the claim has to be narrowed, and one competitor is closer
than any doc in this repo currently admits.**

---

## 1. What is definitively NOT novel

| piece | prior art | verdict |
|---|---|---|
| de Bruijn bubble calling | Cortex / cortex_var, Bubbleparse, McCortex, DiscoSNP++ | **long established.** We ported DiscoSNP++'s `Bubble.cpp` deliberately; claiming the algorithm would be false. |
| assembly-based compression | **Quip** (Jones et al., *NAR* 40(22):e171, 2012) — first assembly-based compressor, de novo assembly via a probabilistic de Bruijn structure, reads stored as positions in contigs | **established 2012.** Our pseudogenome approach is in this lineage. |
| pseudogenome read compression | **PgRC** (Grabowski & Kowalski, *Bioinformatics* 36(7):2082, 2020) | **established.** Our direct architectural rival, already the project's main comparison target. |
| colored dBG compression | ESS-color (*AMB* 2024) | established; compresses the graph itself, a different object. |

Nothing here is disputed. The channel we just built runs *DiscoSNP++'s
algorithm*, and that must be stated as such in any write-up.

## 2. The closest prior art, and it is closer than our docs say

> **BFQzip** — Guerrini, Louza, Rosone, *"Lossy Compressor preserving variant
> calling through Extended BWT"*, arXiv:2304.08534.
> Code: https://github.com/veronicaguerrini/BFQzip

Reference-free FASTQ compression using EBWT + positional clustering, modifying
bases *and* quality scores together, explicitly to preserve variant-calling
information — **one shared structure serving compression and variant calling.**
That is the same *shape* as our central claim, and it is the single most
important citation we were missing (already flagged in
`EBWT_LITERATURE_CORRECTION.md` §3; repeated here because it bears directly on
this channel's novelty, not just on the eBWT discussion).

**The distinctions are real and in our favour, but they must be stated
explicitly rather than assumed:**

| | BFQzip | G_CAPSUL |
|---|---|---|
| compression | **LOSSY** — alters bases and quality scores | **LOSSLESS**, verified by decode-and-diff against the original FASTQ |
| relationship to calling | *preserves* information so an external caller does better | *performs* the calling itself, in-process |
| output | a compressed file | a compressed file **+ a VCF + export/coverage/query** |
| structure | EBWT + positional clustering | pseudogenome + read placements + k-mer counts |

"Preserves variant-calling signal for someone else's caller" and "is itself the
caller" are different claims. But a reviewer who knows BFQzip and does not see
it cited will assume we did not know it existed.

## 3. What does appear to be new — narrowly stated

Neither Quip nor PgRC nor NanoSpring calls variants. Verified for PgRC from the
primary source: it is a compression-only tool; the paper mentions variant
calling only to justify *lossy quality* compression, and its stated future work
is threading and variable-length reads, not analysis. Quip's contigs exist to
store read positions; the paper does not call variants from them.

Conversely, Cortex/McCortex/DiscoSNP++ build a dBG to call variants and
**discard it** — they emit no archive.

So the narrow, defensible statement is:

> **The k-mer index a lossless assembly-based compressor already builds is
> sufficient, unmodified, to run de Bruijn bubble calling — so reference-free
> variant calling can be obtained as a zero-cost byproduct of compression
> rather than from a separately constructed graph.**

Measured support: `kc` (built in `kc_H_build`, 166.7 s, previously consumed to
produce a single coverage threshold `H`) yields 397 bubbles in **0.14 s with no
additional memory**, lifting SNV F1 from 0.888 to 0.910 on HG002 r2.

**What makes that a claim about ENGINEERING ECONOMICS, not about algorithms.**
The novelty is that the structure was already paid for. That is a weaker and
more honest claim than "new variant-calling method", and it is the one the
evidence supports.

## 4. Caveats that must travel with the claim

1. **One window.** r2 is the designated *tuning* window. r3/na/r4/r5 and full
   chr20 are not yet run for this channel. Until then 0.910 is a single point.
2. **The algorithm is DiscoSNP++'s**, ported from `Bubble.cpp:300-347` and
   `509-527`. Cite it as a port.
3. **Precision fell** 0.960 -> 0.925 while recall rose 0.826 -> 0.895. Net
   positive, but it is a recall trade, not a free win.
4. **Indels unchanged** (0.636). This channel is SNV-only so far.
5. A **thorough** novelty claim needs a proper literature search of
   assembly-based compressors beyond Quip/PgRC/NanoSpring/SPRING/Genozip; this
   scan covered the ones this project already benchmarks plus the dBG-calling
   family. It is not exhaustive and should not be described as such.

## 5. Recommended framing for the paper

Do **not** write "we introduce a new variant calling algorithm." Write:

> Existing assembly-based compressors (Quip, PgRC) build k-mer structures and
> discard them after compression; existing reference-free callers (Cortex,
> DiscoSNP++) build equivalent structures and discard them after calling. We
> show these are the same structure, and that a lossless pseudogenome
> compressor's internal k-mer index supports bubble-based variant calling
> unmodified, at negligible additional cost.

That sentence is supported by measurement, concedes the algorithm to
DiscoSNP++, and still states something no cited work states.

Sources: [PgRC](https://academic.oup.com/bioinformatics/article/36/7/2082/5670526) ·
[Quip](https://academic.oup.com/nar/article/40/22/e171/1137542) ·
[BFQzip](https://arxiv.org/abs/2304.08534) ·
[McCortex](https://github.com/mcveanlab/mccortex) ·
[ESS-color](https://link.springer.com/article/10.1186/s13015-024-00254-6)
