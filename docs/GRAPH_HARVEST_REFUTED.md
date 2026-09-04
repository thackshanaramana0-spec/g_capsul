# Graph harvest — built, measured, refuted in one session

2026-09-04. **This supersedes `GRAPH_HARVEST_EXACT_PLAN.md` and
`CHAINING_BRANCH_TAG_PLAN.md` entirely. Both are now wrong.** They are kept
per this repo's standing rule that retractions are marked in place, not
deleted, but nothing in them should be acted on.

---

## 1. What was built

The plan was implemented, not argued about: `stages/106_inprocess.cpp` now
carries a `CAPS_CALL`-gated harvest at the sweep's commit site (inside the
existing `#pragma omp single`, so no race), recording every discarded
alternative when a tail commits with exactly two candidates, plus a histogram
of `ccnt` at commit and the `L` distribution of surviving ties.

**Claim 1 gate: PASSED, bit-identical.** ERR5181310 (SARS-CoV-2),
`ARCHIVE_TOTAL=352997` on both binaries, and a diff of every `[archive]`
per-stream byte count plus every assembly statistic (`round1:`, `links=`,
`second pg`, `leftovers=`) is empty. Only per-job timing lines differ, which
is thread-scheduling jitter.

## 2. What it measured — HG002 r2, 75,115 reads

    commits = 57,553
    ccnt histogram at commit:
        1: 54,966   (95.5%)
        2:  2,453   ( 4.3%)
        3:     99      4: 16    5: 7    6: 5    7: 4    8: 3
    surviving branches (ccnt==2, alternative unclaimed): 85
    tie L as % of Lmax:  100%: 73   90%: 7   80%: 1   70%: 2   30%: 1   10%: 1
    memory: 1 KB

Truth in the same window: **424 het-SNVs, 95 het-indels.**

## 3. Why it cannot work — three independent structural reasons

### 3.1 Exact matching and heterozygosity are mutually exclusive

`rcmp` (`106_inprocess.cpp:565-570`) is a bit-exact comparison. A heterozygous
site is *by definition* a position where the two haplotypes differ. Two reads
from opposite haplotypes therefore **cannot both exactly match the same
suffix** unless the variant falls outside the compared window.

**A het site does not produce a tie in this sweep. It produces exactly one
surviving candidate** — the haplotype that matches. That is precisely what
`ccnt==1` at 95.5% of commits is showing.

This is the reason, and it is not fixable by tuning: it is what exact-match
chaining *is*.

### 3.2 86% of the ties that do occur are duplicate reads

73 of 85 surviving ties sit at `L = 100% of Lmax`. The code's own comment at
the head of the level loop (`106_inprocess.cpp:850-852`) states what that
means:

> "At L = rlen[a] the offset is 0, so the seed is the read's own prefix and
> rcmp compares the two reads in full -- **exactly the duplicate test**."

So the overwhelming majority of harvested "branches" are duplicate reads. The
genuine non-duplicate residue is **~12 ties against 519 truth variants**.

### 3.3 Most two-candidate commits are chains meeting, not branching

2,368 of the 2,453 `ccnt==2` commits (96.5%) had `prv[other] != NONE` — the
alternative was already claimed by another chain. That is two chains
converging on the same read, which is a *merge*, not a branch point.

## 4. Both halves of the plan are dead, including the one I expected to survive

`DISCOSNP_SOURCE_VS_OUR_PLAN.md` split the plan into two claims and predicted A
would survive. **Neither did:**

| claim | prediction | measured outcome |
|---|---|---|
| **A.** Harvest replaces `rc_reads`/`pkidx`/`build_substrate`#2 as a candidate source | "INTACT — this is the RAM/time win" | **REFUTED.** 85 loci cannot replace a structure that supplies thousands. It finds ~0.15% as many candidate sites. |
| **B.** Harvest supplies DiscoSNP++'s accuracy mechanism | already refuted by source reading | **REFUTED again, for a second independent reason** — not only is the closure test missing, the seeds themselves are duplicates. |

The projected table in `GRAPH_HARVEST_EXACT_PLAN.md` §4 (~3.5 GB, ~4 min) is
**void**. It rested on deleting `rc_reads` and the `pkidx` trio, which this
measurement shows cannot be deleted, because nothing replaces what they find.

## 5. Where heterozygous evidence actually lives in this pipeline

Measured on the identical window, same run:

    [MM] placed=12,883  total_mismatches=30,153  mean=2.34/read  zero-mm=1,721 (13.4%)
    [MMTOL] extension mismatches recorded: 17,661 across 11,691 refs

**30,153 mismatch observations against 85 chaining ties** — nearly three orders
of magnitude more signal, and `mm_pos`/`mm_sym` are already stored losslessly in
the archive for Claim 1.

That is the correct answer to "where is the het information": **in the
mismatch streams, not in the chain topology.** Chaining is exact and therefore
blind to heterozygosity by construction; the pigeonhole mapper is
mismatch-tolerant (MAXMAP=11 here) and therefore sees it.

Two important qualifications, so this does not become the next over-claim:

1. **`mm_pos` covers only mapped leftovers** — 12,883 of 75,115 reads (17%)
   here. Reads absorbed into chains carry no mismatch record. So the mismatch
   streams are not a complete het picture on their own.
2. **The `mem_extmm` half of this was already tested this session and lost.**
   Best grid-searched configuration reached F1 57.2 held-out; union with the
   existing channel reached 85.8 against the existing channel's 87.3 alone.
   Recorded so it is not re-attempted as if new.

**And the honest conclusion that follows:** the pileup path already extracts
this signal — by re-deriving it from the reads — which is exactly why SNV F1 is
0.849 and competitive with DiscoSNP++'s 0.847. The information is not being
missed. It is being obtained expensively.

## 6. What this leaves standing

Unchanged and still true, because none of it depended on the harvest:

- Full-scale head-to-head (`results/claim2/t3_t5_full_chr20_headtohead.csv`):
  SNV 0.849 vs 0.847 (parity), indel **0.592 vs 0.576 (win)**.
- The five RAM/speed fixes in `include/caps_caller.h`, compiled and **still
  unvalidated** — projected caller 58 min -> ~34 min, peak 34.94 -> ~20 GB.
  These are ordinary engineering on structures that are now confirmed
  necessary, and they are the remaining real work.
- `CALLER_RAM_SPEED_SKELETAL.md`'s findings (one `#pragma omp` in 2,535 lines;
  `ridx_build` mislabelled; `rc_reads` mostly dead entries) are untouched by
  this and remain the actionable path.

**Cost of finding this out: about 40 minutes and one bit-identical, gated,
reverted-in-spirit code change** — against a plan that projected days of
integration work on a foundation that does not exist. That is the process
working, not failing.
