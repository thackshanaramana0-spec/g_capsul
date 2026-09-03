# Claim 2 — the final plan

Written 2026-09-02, after: reading DiscoSNP++'s source layer by layer, running
it and ARCS on identical reads, killing three of my own hypotheses by
measurement, and validating against the literature that a de Bruijn graph is
**not** required to reach state of the art.

**Status of the inputs to this plan:** the layer stack is settled
(`CALLER_ARCHITECTURE_PLAN.md`), the route options are settled (dBG / eBWT
positional clustering / explicit collapse — all three reach SOTA in the
literature), and the science is not the constraint. What follows is an
engineering plan.

---

## The high-level idea, stated once

Every route to reference-free calling is a different way of achieving **one**
property:

> at a variant site, the two alleles must end up directly comparable.

- DiscoSNP++ gets it from a **dBG**: shared sequence collapses onto one node,
  so a het site is a bubble by construction.
- EBWT2SNP gets it from **eBWT positional clustering**: shared sequence
  collapses by BWT sorting, with no assembly at all.
- We should get it from the thing **G_CAPSUL already computes and already
  stores** — and it turns out we do compute it.

## The realisation this plan is built on

G_CAPSUL's compressor already performs, as part of ordinary compression, the
exact operation a caller needs:

1. **Pigeonhole mapping** places a read onto the pseudogenome **tolerating
   mismatches**, and records every one as `(mm_pos, mm_ref, mm_obs)`. A read
   from haplotype 2 mapped onto a haplotype-1 pseudogenome yields, literally,
   the alternate allele at the right position. **This is a pileup already.**
2. **MEM self-match** finds near-identical regions *anywhere* in the
   pseudogenome and — when `MEM_MAXMM > 0` — extends a match **through**
   mismatches, recording each as `ref_mmpos` / `ref_mmref` / `ref_mmobs`
   (`106_inprocess.cpp:2709-2726`). Two haplotype copies that got assembled
   into separate chains are exactly "near-identical regions", so MEM pairs
   them and **records precisely where they differ.**

That second mechanism is a bubble detector. It was built for compression, it
is already written, tested, and shipped — and it is **switched off**
(`MEM_MAXMM = 0`, line 2548) because it costs compression bytes (+490,763 B
measured). For *calling* those bytes are irrelevant.

**So the plan is not to build a dBG, nor an eBWT index. It is to turn on
machinery that already exists and read streams we already produce.**

This is also exactly the "divergent view of the same assembly" idea: one
assembly, two configurations — compression mode (`MEM_MAXMM=0`) and calling
mode (`MEM_MAXMM>0`), differing in one integer.

---

## Why my first attempt failed, in one line

I captured contigs **pre-MEM** and gave each its own coordinate space,
fragmenting a 400 kb window into 4,537 separate frames — which threw away
the single global pseudogenome coordinate system in which reads from both
haplotypes already stack, and threw away MEM's pairing of the duplicate
chains. The 3.8× redundancy I measured and called a defect is, in the global
frame, **the haplotype pairs themselves**.

---

## The plan, validated layer by layer

| layer | plan | status |
|---|---|---|
| **L0** ingest / k-mer index | unchanged, shared with compression | exists |
| **L1** coverage model | `H` from the count histogram; used to set support thresholds | exists |
| **L2** error filtering | a difference must be supported by ≥ MC reads, not by one MEM extension — this is what separates a het allele from a sequencing error | **rides on L5, no separate pass needed** |
| **L3** substrate | the **global pseudogenome** + MEM references with `MEM_MAXMM>0`. NOT per-chain contigs. One coordinate space. | code exists, gated off |
| **L4** variant enumeration | two independent sources: (a) `ref_mm*` = differences between paired haplotype regions (bubble-equivalent); (b) `mm_*` per-read mismatches = pileup columns in global pg coordinates | **new, small — both are reading streams we already write** |
| **L5** coherence / genotyping | count reads supporting ref vs alt at each pg position from `ppos` + `mm_*`; require both alleles supported; AF band | data exists, needs the counting pass |
| **L6** filters | frozen constants, unchanged | exists, validated |
| **L7** emission + lift | VCF in pg coordinates; existing BWA lift for scoring | exists, works |

The layers hold. The only structural change is **L3: stop fragmenting, use
the global pg frame**, which also makes L4(b) work for free.

---

## Build order

1. **Switch the caller's coordinate frame from per-chain contigs to the
   global pseudogenome.** Delete the span-splitting; `read_pos = ppos[uid]`
   directly. Expected effect: pileup depth per position rises from
   near-nothing to real coverage, so L4(b) starts producing candidates.
   *(This alone may move recall substantially — it is the direct fix for the
   measured failure.)*
2. **Run calling mode with `MEM_MAXMM_OVERRIDE=2`** and emit the `ref_mm*`
   differences as candidate het sites (L4(a)).
3. **Score both sources on r2 (tuning), then r3/na/r4/r5 (held out).**
4. If recall moves but precision drops, tighten L5 support counting — *not*
   the frozen L6 constants.
5. Only if this plateaus far below 0.840: fall back to **eBWT positional
   clustering** (not a dBG) — it is the better-performing mechanism in the
   literature and closer to the suffix/FM machinery this codebase family
   already has.

Each step is ~90 s on one 400 kb window (`scripts/run_window_bench_capsule.sh`).

---

## What is honest to claim, and when

- Right now: **0.419 SNV F1 vs DiscoSNP++'s 0.840** on identical reads. We are
  behind, and no plan changes that until it is measured.
- If steps 1–2 work, the claim becomes attractive and is genuinely ours:
  **the variant evidence is a by-product of compression that is already in the
  archive** — mismatch streams the compressor had to compute anyway. Not
  "we run a caller too", but "the archive already contains the calls".
- The two claims that need no ratio win, and that no competitor can make,
  remain: **calling without re-assembly from a persisted archive**, and
  **cohort calling across archives** (the only route to homozygous variants,
  which every single-sample reference-free caller structurally misses).

**Stop rule unchanged:** if the numbers do not move, reframe Claim 2 around
those two structural novelties rather than promising to beat DiscoSNP++.
Decide on measurements, not on this document.
