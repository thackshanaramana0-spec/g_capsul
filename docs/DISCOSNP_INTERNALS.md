# DiscoSNP++ — layer-by-layer source analysis (cloned v2.6.2-12-gdd16ac4)

Read directly from the clone at `~/DiscoSnp`, not from its paper or docs.
Written 2026-09-02, after DiscoSNP++ turned out to be the strongest competitor
on our own benchmark window (SNV F1 0.840 — see `docs/CLAIM2_RESULTS.md`).

Two things this analysis settles: **what DiscoSNP++ actually needs to run**
(the biopython question), and **why its VCF positions are off by one** — which
we first found empirically and can now point to in its source.

---

## 1. The dependency question — DiscoSNP++ does NOT need biopython

Across the whole tree there are exactly **3 files** importing `Bio`:

```
scripts/community_scripts/discosnp2fasta.py      <- community convenience script
scripts/simulations/targeted_mut_fasta_corrected.py   <- test-data generator
scripts/simulations/random_mut_fasta.py               <- test-data generator
```

**None of them is on the calling path.** Likewise `scipy` / `pandas` /
`numpy` / `networkx` / `matplotlib` appear only in `scripts/k3000/` (a
separate variant-graph post-processing tool) and the same community script.

Confirmed against what our actual run invoked (from its own log):

```
kissnp2            (C++)
kissreads2         (C++)
bwa mem            (external, only because we passed -G)
VCF_creator.py     (python3, stdlib only)
zero2one.py        (python3, stdlib only)
```

Dependency closure of the two Python scripts that really run:

| script | imports |
|---|---|
| `VCF_creator.py` | `os, sys, subprocess, re, time, getopt` + its own `functionObjectVCF_creator`, `ClassVCF_creator` |
| `functionObjectVCF_creator.py` | `os, sys, subprocess, re, time` + own module |
| `ClassVCF_creator.py` | `re, os, sys, subprocess, time` + own module |
| `zero2one.py` | `argparse, fileinput` |

**All standard library.** DiscoSNP++ needs only `python3` (≥3.0, checked
explicitly at the top of `run_VCF_creator.sh`), a C++ toolchain, and `bwa`
if `-G` is used.

**The biopython requirement was Kmer2SNP's, not DiscoSNP++'s.** Kmer2SNP
imports `Bio.SeqIO` in `libprism/local/generate_VCF_without_reference.py`
on its main path, and additionally needs the R package `findGSE` — which is
what is currently blocking it (`hete.para` comes out empty, `Rscript` is not
on PATH). Do not conflate the two tools' requirements.

---

## 2. The layers

### Layer 1 — `kissnp2` (C++, `tools/kissnp2/src/`)

Bubble detection on a GATB de Bruijn graph (`thirdparty/gatb-core`, a git
submodule). `Bubble.hpp`'s own header comment: *"class that tries to build a
bubble from a starting node ... and potentially extend it with right and left
unitigs/contigs ... instantiated N times, one per thread."*

Key knobs visible in `Bubble.hpp`: `max_polymorphism`, `max_depth` /
`max_breadth` (unitig/contig extension limits), `max_sym_branches`. The
bubble is the whole variant model — a local branch in the graph where two
paths diverge and reconverge. **This is structurally the same signal CAPSULE's
cross-contig pass looks for**, which is exactly why our own caller's
`extract_snv_bubble` / `extract_bubble` are the right shape; the difference
is that DiscoSNP++ finds bubbles in a *graph* it built for this purpose,
while we look for them between *contigs* our compressor happened to produce.

Output: `*.fa` of predicted variant paths (higher/lower path per bubble).

### Layer 2 — `kissreads2` (C++, `tools/kissreads2/src/`)

Maps the original reads back onto the predicted bubble paths to establish
**read coherence** and per-sample genotypes/coverage. Files
(`fragment_index.*`, `interface_xhash.*`) show it is a seed-index + hash
approach, not a full aligner. Its `-co` / `-unco` outputs split predictions
into coherent and incoherent; only the coherent ones become the final VCF.
This is DiscoSNP++'s precision mechanism, and it is why its precision on our
window is 0.971.

### Layer 3 — `bwa mem` (only when `-G <ref>` is passed)

Purely for **placing** the bubble paths onto a reference so the VCF can carry
genome coordinates. Same role as BWA in our own pipeline: evaluation/lift, not
calling. Without `-G` the VCF is emitted in bubble-local coordinates and
carries no mapping info (its own `--help` says exactly this).

### Layer 4 — `VCF_creator.py` + `ClassVCF_creator.py` (python3, stdlib)

Parses the SAM produced above and turns bubbles into VCF records. Relevant
mechanics:

- `ClassVCF_creator.py:452` — `self.mappingPosition = abs(int(self.listSam[3]))`.
  `listSam[3]` is the **SAM POS field, which is 1-based by the SAM spec**.
- `ClassVCF_creator.py:369` — `mappingPositionCouple = mappingPosition +
  int(correctedPos[0]) + int(mappingPositionCouple) - 1`, folding in the
  variant's offset within the path, with a `-1`.
- It also parses the `MD` tag (`GetTag`, line ~719) to locate mismatches, and
  handles reverse-strand paths and close/multiple SNPs per bubble (hence the
  `9_1`, `9_2`, `102_1` style IDs in the output).

### Layer 5 — `zero2one.py` (python3, stdlib) — **the off-by-one**

The last step of `run_VCF_creator.sh` (lines 368–373). Its entire job:

```python
new_POS = int(vcf_line[1]) + 1          # and the same +1 applied to every XA= position
```

It is applied to **both** the main VCF and the `_for_IGV` VCF — which is why
those two files have byte-identical POS columns in our run, despite the
program's own final message still claiming the IGV file is "0-based". That
message is stale.

---

## 3. The off-by-one, confirmed three ways

We found this empirically before reading the source; the source now explains
it, and the version history shows it is the shipped state.

1. **Reference-base agreement.** Of 331 mapped SNV records, the VCF's REF
   base matches the GRCh37 chr20 sequence at **`POS-1` for 326**, at `POS`
   for only 80.
2. **Scoring.** Through `rtg vcfeval` against GIAB truth: as emitted,
   **F1 = 0.004** (TP=2, FP=303). With `POS-1`, **F1 = 0.840** (TP=296,
   FP=9). A coordinate shift is the only thing that moves a score that way.
3. **Provenance.** `zero2one.py` was added 2021-11-10 in commit `7f77dcc`
   *"prepare one-based VCF"*. Thirty commits later, at
   `v2.6.2-12-gdd16ac4`, no coordinate correction follows it. The `+1` is
   applied on top of a position already derived from a 1-based SAM field,
   producing output one too high.

**Operational rule: any benchmark that runs this DiscoSNP++ version with
`-G` and does not subtract 1 from POS is measuring noise, not the tool.**
Uncorrected it looks catastrophically bad (0.004) when it is in fact the
strongest caller on this window (0.840). The outer ARCS project's docs note
the same quirk for indels; this analysis pins it to the exact script and
commit, and shows it applies to SNVs too.

---

## 4. What this means for us

- **DiscoSNP++ is not fragile or hard to run** — no biopython, no R, no
  scientific stack; it ran our 400 kb window in ~1 second. It is a serious,
  well-engineered competitor, and the 0.840 bar is real.
- **Its bubble model is the same abstraction our cross-contig pass uses.**
  The gap between our 0.419 and its 0.840 is therefore not a difference of
  idea, but of the substrate the bubbles are found in: a purpose-built dBG
  versus contigs our compressor produced as a side effect. That is the honest
  framing for the paper, and it is consistent with the same-caller-different-
  assembly result in `docs/CLAIM2_RESULTS.md`.
- **Kmer2SNP is the one with the heavy dependency chain** (biopython +
  R/findGSE), and it remains blocked on `findGSE`.
