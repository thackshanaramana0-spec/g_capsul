# Calling from what compression already stored

## The question

Claim 2's published path (`CAPS_DBG=1 CAPS_DBG_ONLY=1`) builds a 140.7 M-node
de Bruijn graph and hunts bubbles. That is DiscoSNP++'s algorithm. DiscoSNP++
runs it because it has nothing else: raw reads, no reference, no alignment.

We are not in that position. Compression has already placed every read.
So why run their method at all?

## What the literature actually says (checked 2026-09-06)

**Verified.** Xu et al. 2017, *Scientific Reports* 7:10963
(`s41598-017-10826-9`), simulating >3M SNVs from the whole human genome:
"The assembly-based approach had a much lower recall rate and precision
compared to the alignment-based approach that would recover 99% of imputed
SNVs" at 30x coverage.

**Documented, weaker form:** assemblers "often collapse two heterozygous
alleles (typically represented as bulges in the de Bruijn graph) into a single
copy to increase the contiguity of the consensus assembly" (LJA, Bankevich et
al.).

**WITHDRAWN.** An earlier draft of this file quoted "heterozygous SNVs would
have no more than 50% chance to be called correctly" and attributed it to two
URLs. Re-searching does not reproduce that sentence. It is withdrawn; do not
cite it.

### The citation does not apply to us -- a category error worth recording

Xu et al. compare ALIGNMENT-BASED (align reads to a REFERENCE genome) against
ASSEMBLY-BASED (assemble, then call). Their 99% recall belongs to reference
alignment.

Config B below has NO reference: it is a pileup over reads placed on our own
de novo pseudogenome. In their taxonomy **B is assembly-based, and so is A**.
Mapping a finding about reference alignment onto a reference-free pileup was
the error that motivated this experiment. The thing that gives their
alignment-based arm 99% recall is precisely the reference -- which a
reference-free tool does not have and cannot borrow.

The experiment settled the question empirically anyway (B lost by measurement,
not by citation), and the one literature-derived idea that survived is the one
that was MEASURED on our own data rather than quoted: the k-mer proximity
failure, below.

## What this codebase already emits, and never scored

Two channels come out of the ENCODER, before any caller runs.

### mem_extmm — already in the archive

When the MEM extension stage stores a read as a copy of an earlier pseudogenome
region, the bases that differ are recorded so decode can reproduce the original.
Those mismatches are (position, reference base, observed base) triples: a
variant call set. On full HG002 chr20:

| stream | bytes |
|---|---|
| `mem_extmm_cnt` | 225,508 |
| `mem_extmm_pos` | 1,832,378 |
| `mem_extmm_obs` | 468,669 |
| **total** | **2,526,555 B = 2.53 MB** |

2,518,219 SNV records for 2.53 MB — about **one byte per candidate**, 0.46% of
the 547 MB archive. This is not an added cost: **remove these streams and the
archive is no longer lossless.** The calls are a re-reading of bytes the format
is already obliged to carry.

### near-miss — free inside the overlap loop

Chaining probes on a 16-base seed and verifies with an exact comparison. A read
that shares the seed but differs by EXACTLY ONE base over a long verified
overlap is discarded at `if(!rcmp(...)) continue;` — and that read is the other
haplotype. It is already in registers when it is rejected. Full HG002 chr20:

```
1-mismatch observations : 13,019,987   (12,328,180 placed)
distinct (pos,alt)      :    513,342
  depth >=2: 309,679   >=3: 233,891   >=5: 142,203   >=10: 49,846
truth het-SNV in confident regions:      45,885
```

## The measurement gap found 2026-09-06

`mem_extmm` records are emitted into the official VCF as `mcontig_N`, indexed by
`g_contig_spans`. **Nothing ever dumped that contig space.** `CAPS_DUMP_CONTIGS`
writes the caller's `dcontig`/`lcontig`/`bcontig` spaces only, so every
`mcontig_` record hit `lift_vcf.py:36` (`if rn not in cmap: return None`) and was
dropped before scoring.

On the A-config full chr20 run: 2,632,599 records emitted, of which **2,520,833
(96%) are `mcontig_` and were discarded at evaluation**; 47,238 lifted, all from
the 111,766 `dcontig` records.

**The published F1 is unaffected** — the records were dropped, never counted, so
they could neither inflate nor deflate the score. But an entire
compression-derived channel had never been evaluated in either direction.

Fixed by `CAPS_DUMP_MCONTIGS`, which dumps the referenced spans only (mirroring
how the DBG-ONLY path already uses `lc_used`).

## Novelty

Searched twice. Every reference-free caller found — kSNP, DiscoSNP++, SKA lo,
Kmer2SNP, EBWT2INDEL, LueVari, the 2025 distributed dBG pipeline — builds an
index or graph from raw reads and hunts bubbles. ROVER is read-pair-overlap
aware but alignment-based and needs a reference. None derive calls from a
compressor's own overlap-rejection set or from streams the archive stores for
losslessness.

## Measured: does the pileup beat the bubbles? NO.

Full HG002 chr20, same box, sequential, identical scoring
(`bwa mem` lift + `rtg vcfeval --squash-ploidy`, GIAB v4.2.1 confident
regions, 45,885 truth het-SNVs).

| config | what it is | caller | peak RAM | P | R | **F1** |
|---|---|---|---|---|---|---|
| **A** | bubbles only (`CAPS_DBG_ONLY`) — DiscoSNP's method | **106.0 s** | **6.03 GB** | 0.9328 | 0.8268 | **0.8766** |
| **B** | pileup on encoder placements (`CAPS_NO_INDELS`) | 851.8 s | 14.09 GB | 0.9627 | 0.6950 | **0.8072** |
| **C** | both channels merged (`CAPS_DBG` + pileup) | ~950 s | ~14 GB | 0.7933 | **0.8736** | 0.8315 |

**The hypothesis that motivated this experiment is refuted.** The pileup is
worse on every axis: −0.069 F1, 8x the caller time, 2.3x the RAM. The
prediction taken from the literature — that alignment-style pileup beats dBG
bubbles on heterozygous sites — does not transfer to this implementation.

The mechanism is in the precision/recall split. B is MORE precise (0.963 vs
0.933) and much less sensitive (0.695 vs 0.827). The collapse pass that puts
both alleles into one frame is also what destroys sites: 116,249 candidates
become 43,708 SNVs, of which 36,041 lift, against 47,238 for bubbles. The
pileup is conservative, not better.

Consequence: **the k-mer graph stays.** Removing `kc` from the pileup path
(its `kcount` flank-repeat filter at `caps_caller.h:4252`, plus the scalar `H`)
was abandoned as pointless once the path it serves was measured to lose.

An honest correction to the framing that opened this work: the saving was
never going to come from replacing the graph with the pileup, because the
pileup costs MORE. It can only come from channels that cost nothing, which is
what the two encoder channels above are.

## Why the graph misses what it misses (measured, full chr20)

For every truth site, distance to the nearest OTHER truth variant:

| set | median nearest neighbour | within 31 bp (= k) |
|---|---|---|
| A **missed** (FN, n=7,720) | 94 bp | **30.4%** |
| A **found** (TP, n=36,855) | 224 bp | 10.2% |

A site the graph misses is **3x more likely to have another variant inside one
k-mer**. This is the documented dBG failure mode -- a second SNV inside the
k-mer destroys the shared flank, so no bubble can form -- and it is a
STRUCTURAL ceiling on bubble recall, not a tuning problem. It is also why
raising coverage or sweeping thresholds cannot fix those 7,720 sites.

This predicts which channel can rescue them:

- **near-miss requires EXACTLY ONE mismatch** over a 74-148 bp verified overlap,
  so it shares the blind spot -- arguably a worse version of it, since its
  window is longer than k. (The 2-mismatch observations are counted --
  `nm2 = 5,940,724` -- but not emitted.)
- **mem_extmm carries `MEM_MAXMM = 2` mismatches per extension**
  (`106_inprocess.cpp:3779`, gated on the heterozygosity pre-scan; array cap
  `REF_MAXMM = 4`). A clustered PAIR is representable in one extension, which is
  exactly the case that breaks the k-mer flank.

So the mechanistic prediction is that `mem_extmm`'s rescues should be ENRICHED
for short nearest-neighbour distance relative to A's missed set. `marginal.sh`
tests that directly; if the enrichment is absent the story is wrong whatever
the F1 says.

## A naive union is NOT free: measured

Before C ran, its TP count was predicted from the overlap of A and B
(B rescues 2,083 of A's 7,720 missed sites, so TP should reach ~38,938):

| | predicted | actual |
|---|---|---|
| TP | 38,938 | **38,941** |
| FP | ~3,655 | **10,144** |

The complementarity model was right and the FP model was badly wrong. Recall
rose to 0.8736 -- essentially the computed ceiling -- but precision COLLAPSED
from 0.933 to 0.793 and F1 fell BELOW A.

Two causes, one flagged in advance and one general:

1. C's bubble channel is not A's (131,664 bubbles vs 115,842), so C is a union
   of a NOISIER bubble set with the pileup, not a clean A+B. This was recorded
   before the result was seen, not after.
2. The general lesson: **a free channel is not automatically worth adding.**
   2,083 rescued truth sites arrived with ~7,500 extra false positives. Union
   by concatenation buys recall and pays more than it buys.

Consequence for the encoder channels: they must be THRESHOLDED, never merged
raw. `mem_extmm` emits 2.5M candidates against 45,885 truth sites; unfiltered it
would destroy precision outright. `compchan.sh` therefore sweeps a recurrence
threshold rather than taking the union.

**A remains the best configuration at F1 0.8766.** The paradigm question as
posed -- can the placement pileup replace the graph -- is answered NO on all
three configurations.

## MEASURED: the archive channel (E), full chr20

`mem_extmm` records lifted and scored for the first time (the contig-space dump
and the `LIFT_KEEP_INFO=1` evidence passthrough both had to be written to make
this measurable at all).

**Unfiltered, the channel has MORE recall than the graph:**

| | TP | FP | P | R |
|---|---|---|---|---|
| E, all records | 40,159 | 1,087,454 | 0.036 | **0.9009** |
| A, bubbles | 36,855 | 2,655 | 0.933 | 0.8268 |

and it **rescues 5,578 of A's 7,720 missed sites (72.3%)**. The variants ARE in
the archive.

**But the evidence cannot isolate them.** For a union with A to raise F1:

    F1' = (73710 + 2dTP) / (84085 + dTP + dFP) > 0.8766  =>  dTP/dFP > 0.78

i.e. **more than 44% of added calls must be genuine rescues.** Measured:

| filter | n | TP | FP | P | rescues | marginal P |
|---|---|---|---|---|---|---|
| RC>=2 | 180,727 | 23,428 | 103,360 | 0.185 | 2,767 | 2.6% |
| RC>=3 | 46,024 | 4,979 | 20,081 | 0.199 | 739 | 3.5% |
| RC>=5 | 7,672 | 107 | 2,030 | 0.050 | 32 | 1.6% |
| MMCNT=1 | 144,434 | 15,351 | 103,087 | 0.130 | 1,468 | 1.4% |
| MMCNT=1,RC>=2 | 27,850 | 9,142 | 11,676 | 0.439 | 676 | 5.5% |
| MLEN>=200 | 63,972 | 25,195 | 32,882 | 0.434 | 1,725 | 5.0% |
| MLEN>=500 | 20,266 | 7,324 | 11,062 | 0.398 | 508 | 4.4% |
| MMCNT=1,MLEN>=200 | 16,178 | 8,446 | 6,479 | 0.566 | 476 | 6.8% |
| MMCNT=1,MLEN>=500,RC>=2 | 1,637 | 1,173 | 333 | **0.779** | 57 | **14.6%** |

**Best marginal precision 14.6% against a 44% bar.** Best-case union:
F1 = 0.8739 < A's 0.8766. **No filter setting improves on the published config.**

### Why, mechanistically

**Anchor quality selects for EASY loci.** A long, clean, single-mismatch MEM
anchor sits in unique well-covered sequence -- exactly where the k-mer graph
already succeeds. Filtering hard raises precision to 0.779 but keeps
overwhelmingly the sites A ALREADY calls (at MMCNT=1,RC>=2: 9,142 TP of which
only 676 are rescues -- 8,466 are duplicates). The 5,578 rescues live in short,
messy, multi-mismatch anchors, statistically indistinguishable from the 1.09M
errors around them. Filtering for confidence discards precisely the sites that
made the channel interesting.

### Two hypotheses of mine, both REFUTED by this data

1. **Recurrence is evidence.** It is ANTI-correlated with truth: base rate 2.96%
   TP, RC>=5 gives 1.4%. The pseudogenome stitches both haplotypes AND repeat
   copies, so "many regions disagree here" is the signature of a REPEAT -- the
   exact thing a caller must reject. Predicted before measuring, wrong.
2. **The clustered-variant mechanism.** Predicted E would rescue the sites whose
   k-mer flank is destroyed by a neighbouring variant. Measured: A's missed set
   is 30.4% clustered, E's rescues 27.4% -- no enrichment, slightly less. E
   rescues broadly, for a reason still unknown.

### What is missing, concretely

The encoder channels carry **no base quality at all**. A sequencing error and a
het allele are literally identical to them. The caller separates exactly this
way (`QMIN=20`, `caps_caller.h:3122`), and quality is the LARGEST thing in the
archive (qual_body 470.6 MB, 82%). That is the one untested axis that is both
principled and already paid for.

## MEASURED: the near-miss channel (D), full chr20

First valid scoring (the initial run was void -- the lift discarded the INFO
fields the filters needed).

| filter | n | TP | FP | P | rescues | marginal P |
|---|---|---|---|---|---|---|
| all | 104,121 | 34,503 | 51,686 | 0.400 | 3,923 | 7.1% |
| DP>=3 | 75,605 | 25,920 | 36,787 | 0.413 | 2,850 | 7.2% |
| DP>=5 | 42,310 | 12,661 | 21,730 | 0.368 | 1,331 | 5.8% |
| DP>=10 | 9,617 | 1,045 | 5,267 | 0.166 | 159 | 2.9% |
| **AF 0.25-0.75** | 29,927 | 17,300 | 7,936 | **0.686** | **1,873** | **19.1%** |
| DP>=5 & AF band | 6,289 | 3,254 | 1,877 | 0.634 | 289 | 13.3% |
| DP>=10 & AF band | 379 | 11 | 215 | 0.049 | 1 | 0.5% |

**ALLELE FRACTION IS THE DISCRIMINATING STATISTIC.** It raises precision
0.400 -> 0.686 while KEEPING MORE rescues than depth filtering (1,873 vs 1,331).
Depth alone is actively harmful (DP>=10 -> P 0.166). This is the difference
between evidence describing the ALLELE and evidence describing the ANCHOR, and
it is exactly why E cannot work: `MLEN`/`MMCNT`/`RC` say nothing about the 50/50
signature of heterozygosity. (The AF denominator was added this session; the
channel had no denominator at all before.)

Best union at AF 0.25-0.75: F1 = 0.8249 < A's 0.8766.

## The quality + true-VAF lever (built and measured 2026-09-06)

Chasing "the channels carry no base quality" exposed a deeper bug first: the
`AF` field was **mathematically impossible** -- AF>1 on 43.7% of sites.

**Why, and it is structural, not arithmetic.** `nmcov` counts reads PLACED at a
pg position; near-miss support counts reads REJECTED there for differing by one
base. **The pseudogenome never contains a pileup** -- reads matching the
assembly are chained in, reads carrying the other allele are stitched elsewhere
-- so the numerator and denominator were measured on two different populations.
(This is the same structural fact that refuted PG-ANCHOR.)

Reformulated as a genuine variant allele fraction for this data structure:

    REF support = nmcov[p]  (reads placed here, i.e. agreeing with the pg)
    ALT support = distinct reads rejected here for one mismatch
    VAF = ALT / (ALT + REF)

bounded in [0,1], with a clean **heterozygous mode at 0.5-0.6 holding 40% of
sites** -- a signal completely invisible under the old formula. Also fixed:
support now counts DISTINCT reads, not (a,b) pair observations (one read a
pairs with many partners b), which removed 4.5% pure inflation.

Base-quality mask: one bit per base, "is this base >= Q20", reusing `woff` (the
2-bit packing's own offsets) so base j of read u is bit `woff[u]*32+j` -- no
second index, 240 MB at full chr20, built during the parse that already reads
the quality line. Built ONLY when `CAPS_NM_VCF` is set. Measured: 92,468 of
513,342 sites (18%) have NO high-quality observation at all.

| filter | n | TP | FP | P | rescues | marginal P |
|---|---|---|---|---|---|---|
| all | 100,200 | 34,456 | 48,349 | 0.416 | 3,910 | 7.5% |
| VAF 0.30-0.70 | 48,521 | 25,378 | 15,355 | 0.623 | 2,852 | 15.7% |
| VAF 0.35-0.65 | 32,175 | 19,239 | 7,757 | 0.713 | 2,169 | 21.9% |
| HQ>=3 alone | 72,082 | 25,303 | 34,754 | 0.421 | 2,752 | 7.3% |
| HQ>=3 + VAF | 34,055 | 21,596 | 7,510 | 0.742 | 2,316 | 23.6% |
| **HQ>=3 + VAF_hq** | 26,474 | 17,076 | 5,587 | **0.754** | 1,864 | **25.0%** |

**Marginal precision 7.5% -> 25.0% (3.3x); precision 0.416 -> 0.754.** Best
union: F1 0.8460 < A's 0.8766. **A still stands.**

**Quality is NOT inert** (a prediction recorded before the run, and wrong):
adding HQ>=3 to the VAF band moved marginal precision 21.9% -> 25.0%. But
quality ALONE is useless (7.3%). The two axes are independent and only work
together -- VAF asks "is the support balanced", quality asks "is the support
trustworthy". Encoder cost of both: 287.17s/6.03GB -> 279.46s/6.36GB, i.e. no
time penalty and +0.33 GB.

**Next lever this exposed, unbuilt:** near-miss requires EXACTLY ONE mismatch,
so a read spanning two nearby het sites is rejected outright -- it is blind to
clustered variants for the same reason the k-mer graph is, and that is 30.4% of
A's misses. `nm2 = 5,940,724` two-mismatch observations are counted and
discarded. Emitting them as candidate PAIRS targets the one failure mode both
current channels share.

## Summary of every configuration measured

| config | what | F1 | vs A |
|---|---|---|---|
| **A** bubbles (published) | k-mer graph, DiscoSNP's method | **0.8766** | -- |
| B pileup on placements | collapse + pileup | 0.8072 | worse |
| C A+B merged | | 0.8315 | worse |
| D near-miss, best filter (VAF+quality) | encoder overlap rejects | union 0.8460 | worse |
| E mem_extmm, best filter | archive mismatch streams | union 0.8739 | worse |
| DiscoSNP++ (same box) | | 0.847 | -- |

Best marginal precision achieved by ANY encoder channel: **25.0%** (near-miss
with true VAF + base-quality mask), against a **44%** bar. **A stands.**

## VERDICT

**As published, configuration A is SPRING + DiscoSNP bolted together.**

| claim | verdict |
|---|---|
| archive data improves A's calls | **No** -- one boolean (the ploidy gate); `CAPS_DBG_PLACEMODE` is off by default so placements are diagnostic only |
| calling data is stored for free | **No** -- +218,555 B, 0.038% |
| the archive channel beats/augments A | **No** -- best union F1 0.8739 < 0.8766 |
| compression lets us skip assembly | **No** -- WITHDRAWN; bubble finding never needed a substrate, and DiscoSNP has no ploidy gate to serve |

What IS defensible: reference-free calling at **F1 0.8766 vs DiscoSNP++ 0.847**
on the same box, in one pass, sharing the compressor's process and 2-bit read
store. The F1 edge comes from `MINC=2` vs their `abundance_min=3`, a parameter
choice, not from the archive.

**"The archive is the caller" is NOT supported by any measurement in this file.**

## Status

A, B and C measured. E measured, negative for F1. Fusion cost measured: +0.038%.
D (near-miss) still unmeasured -- its first scoring run was invalidated by the
lift dropping INFO fields, now fixed but not re-run. C (union of both) running. The two encoder-derived channels
(D near-miss, E mem_extmm) and the A+E union are scored by `compchan.sh`.
No F1 is claimed for D/E/F until it lands.
