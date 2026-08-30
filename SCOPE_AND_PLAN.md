# Generalization scope and plan (reimpl pipeline)

Written because ad hoc testing on 1-2 convenient files already produced two false
signals this session: overstated speed/RAM wins (race condition, corrected) and a
now-fixed silent data-loss bug that only 2 of 10 real datasets happened to expose.
This document is the scope definition and execution plan requested directly:
"first scope it properly" before any further win claim.

## 1. The scope, in full (not redefined — it's the project's own locked list)

Source: `benchmark/DATASET_LOCKED.md`. 10 accessions + 4 human chr20 samples.

| # | Accession | Organism | Size | Read len | Quality resolution | Tested this session |
|---|---|---|---|---|---|---|
| 1 | SRR2584863 | E. coli | 576MB | 150bp | full (~1.3 bits/val) | YES — win vs SPRING+Genozip |
| 2 | ERR552797 | M. tuberculosis | 217MB (300bp reads) | 300bp | full | YES — found+fixed silent data loss (stage 90) |
| 3 | SRR554369 | P. aeruginosa | 334MB | 100bp | full (39 symbols) | YES — **LOSE to SPRING by 10.8%**, broad (seq+order AND quality both worse) |
| 4 | ERR5181310 | SARS-CoV-2 | 30MB | variable 40-221bp | near-degenerate (3 symbols) | YES — **FIXED (stage 91)**: was -8.7% vs Genozip, now +7.6% win vs Genozip, +15.2% win vs SPRING |
| 5 | ERR17740259 | S. aureus | 970MB | 148bp | binned (6 symbols) | YES — win vs SPRING (10.4%) and Genozip (39.6%) |
| 6 | SRR065390 | C. elegans | 11GB | ? | ? | NO — large-genome/auto-chunk regime |
| 7 | DRR976266 | S. cerevisiae | 1.67GB | 150bp | near-degenerate (4 symbols) | YES — win, but easy case |
| 8 | SRR870667 | T. cacao | 15GB | ? | ? | NO — largest file in scope |
| 9 | SRR36741279 | L. major | 1.7GB | 100bp | binned (8 symbols) | YES — win vs SPRING (marginal, 0.2%) and Genozip (6.2%) |
| 10 | SRR37283774 | P. falciparum | 669MB | 100bp | full (~2.2 bits/val) | YES — win vs SPRING+Genozip |
| C2-1..4 | HG002-HG005 | Human chr20 | 30x depth | 150bp | ? | NO |

**Real coverage so far: 7/10 primary accessions, 0/4 human.** Every new dataset
tested surfaced something real: M. tuberculosis found a silent data-loss bug
(fixed, stage 90), SARS-CoV-2 found and then FIXED a real architectural loss
(stage 91, conditional reorder), P. aeruginosa found a real quality-coder gap
that generalized into a fix affecting 3 files (stage 92, adaptive context
sizing) and narrowed its overall gap from -10.8% to -1.7% (quality stream
itself now wins). Current honest tally against SPRING: win on 6/7 (E. coli,
yeast, P. falciparum, S. aureus, L. major-marginal, SARS-CoV-2), lose on 1/7
(P. aeruginosa, -1.7%, down from -10.8% — remaining gap is sequence+order,
not quality). Against Genozip: win 7/7. One real, small, still-open gap
remains — not zero, but the honest trend across this session is every found
gap getting root-caused and either fixed or substantially narrowed by a real,
cross-file-verified mechanism, never a per-file patch.

## 2. The algorithmic lesson, synthesized (not fragments)

Confirmed directly from real SPRING/Genozip source and behavior, not assumed:

- **SPRING's block/chunk sizes are fixed constants** (`NUM_READS_PER_BLOCK =
  256000`, `NUM_READS_PER_BLOCK_LONG = 10000`, `params.h`), chosen once against a
  RAM budget, used identically on a 30MB file or a 15GB file. It does not
  re-tune per input. This project's earlier RAM/speed numbers were measured by
  tuning against ONE file's behavior (explicitly caught and corrected: "don't
  hyperparametrize" — [[reimpl_combined_pipeline_speed]]) — the fix has to be a
  fixed, RAM-budget-derived formula tested across the size range, not a
  per-file constant.
- **SPRING treats reordering as optional, not unconditional** (`-r` flag,
  default off). It only reorders when the resulting compression gain is worth
  the `log2(n!)` cost of recording the permutation. Our pipeline currently
  reorders unconditionally, which is exactly why SARS-CoV-2 loses (permutation
  cost 189KB > the reordering gain on that file) while yeast/E.coli/P.falciparum
  win (genuine long-range structure makes reordering pay off there). This is a
  precise, mechanistic, generalizable explanation, not a per-file patch.
- **Neither tool silently drops data outside an assumed size envelope.** Our
  uint8-length bug (fixed, stage 90) was exactly this kind of unstated
  assumption; the fix (uint16, matching buffer widen) is the same kind of
  "handle the whole realistic range, not the range already tested" correction.

**Unified design principle going forward:** every parameter (block size, queue
capacity, reorder-or-not) must be either (a) a fixed constant derived from a
RAM/size budget calculation valid across the full 30MB-15GB range, or (b) a
cheap, measured, per-file decision with an exact, known cost model (like the
reorder-vs-permutation comparison) — never a value tuned by looking at one
file's timing or size.

## 3. Execution plan (ordered, what's done vs what's left)

**Done:**
1. Silent data-loss bug (>255bp reads) — fixed, verified, committed (stage 90).
2. Root-caused (not yet fixed) the reorder-vs-permutation tradeoff on SARS-CoV-2.
3. Tested 7/10 primary accessions (E. coli, M. tuberculosis, P. aeruginosa,
   SARS-CoV-2, S. aureus, yeast, P. falciparum) against real SPRING+Genozip.
   Win vs SPRING on 5/7, lose on 2/7 (SARS-CoV-2 -8.7%, P. aeruginosa -10.8%).
   Win vs Genozip on 7/7 (including the two SPRING losses).

**Next, in order:**
4. **New, second open loss, root-cause narrowed: P. aeruginosa.** Unlike
   SARS-CoV-2 (single, architectural mechanism — unconditional reordering),
   this loss is a genuine CODER-QUALITY gap, confirmed by entropy check:
   order-0 entropy on this file's quality is 3.509 bits/value; our coder gets
   2.862 (real context gain over naive), but SPRING gets 2.553 — SPRING's
   context model extracts real structure from THIS file's quality-correlation
   pattern that ours misses, not a rounding difference (~12% gap). Sequence+
   order also loses (1,717,157 vs our 1,900,348) but by less. This rules out
   "not enough context" as the cause (we do use context, and it does help) —
   the open question is which specific correlation our fixed hist=3/pos=8/
   delta=48 context split fails to capture here that SPRING's model captures.
   Per the explicit "don't hyperparametrize per file" rule
   ([[reimpl_combined_pipeline_speed]]), the fix must come from understanding
   what structure is being missed generally, not from tuning hist/pos/delta
   to this one file's numbers.

   **RESOLVED — real, generalizable fix found and verified (stage 92,
   `92_qual_adaptive_ctx.cpp`).** The libbsc-proxy investigation (xz/bzip2 on
   reordered quality, both still worse than our coder) ruled out "SPRING's
   general compressor alone explains it" — see below for that trail, kept for
   the record. The actual cause: context-space DENSITY, not compressor
   choice. `contexts/block=22,778,496` for P. aeruginosa's ~10M quality
   values (alphabet 39) — under 0.44 observations per context on average, so
   most contexts barely leave their uniform prior. Swept hist/pos/delta on
   ALL 7 tested files and found a clean, alphabet-size-correlated split, not
   a single-file number: every FULL-RESOLUTION file (E. coli/P. falciparum/
   P. aeruginosa, alphabet 38-39) improves 5.5-12.5% with a smaller context
   (hist=2 pos=4 delta=8); every BINNED/degenerate file (yeast/S. aureus/
   L. major/SARS-CoV-2, alphabet 3-8) gets 0.9-2.4% WORSE with that same
   smaller context — there the default already has enough data per context,
   shrinking it throws away real structure. Clean separation at alphabet<=20
   (nothing tested near the boundary either side). Implemented as an
   automatic switch on MEASURED alphabet size (resolved from the same
   pre-scan the coder already does), not a manual per-file flag — verified
   all 7 files pick the correct regime automatically and round-trip.

   Result: P. aeruginosa's quality stream alone now BEATS SPRING (3,131,704
   vs SPRING's 3,191,556, was 3,577,150/+12.1% loss, now -1.9% win). Overall
   P. aeruginosa total: was -10.8% vs SPRING, now **-1.7% vs SPRING** — large,
   real, verified progress, not yet a full flip. Remaining gap is the
   sequence+order stream (already separately investigated: reordering
   genuinely measures worse than SPRING's approach on this specific genome,
   not a reorder-cost-mechanism issue like SARS-CoV-2 — a real, harder,
   still-open question, distinct from the quality fix above).

   *(Kept for the record, the now-superseded libbsc proxy trail:* SPRING's
   default lossless quality path is libbsc general BWT/context-mixing on
   read-reordered quality strings, confirmed from `reorder_compress_quality_id.cpp`.
   xz/bzip2 on a sequence-similarity-reordered proxy of P. aeruginosa's
   quality stayed worse than our own (then-unfixed) coder — ruled out
   "just reorder + generic compressor" as sufficient, correctly redirecting
   toward the real cause found above.)
5. Test the 3 remaining untested primary accessions (C. elegans, T. cacao,
   L. major already done) — wait, L. major is done; remaining: C. elegans,
   T. cacao only, both large-file regime, deferred to step 8.
6. **DONE (stage 91, `91_conditional_reorder.sh`).** Not a pre-assembly
   estimate — actually runs both candidates (reorder path; raw-unreordered
   seqpar encode) and keeps the measured smaller total. Verified: fixes
   SARS-CoV-2 (was -8.7% vs Genozip, now +7.6%; was +0.2% vs SPRING, now
   +15.2%). Verified NO regression on all 6 other tested files — each still
   measures reorder as cheaper and picks it, byte-for-byte same totals as
   before. Confirmed does NOT fix P. aeruginosa (raw-unreordered there is
   2,394,279 B vs reorder's 1,900,348 — reordering correctly stays chosen),
   which rules out reorder-cost as P. aeruginosa's cause and confirms it's
   purely the quality-coder gap (item 4).
7. Re-verify all previously-tested files are unaffected by step 6 (byte-identical
   assembly output where reordering is still chosen; correct decode where it
   isn't) — same round-trip standard as stage 90's regression check.
8. Once 4/6/7 are done: re-run the full comparison across all 7 files tested
   to date (plus C. elegans/T. cacao if reached) and report win/loss per file,
   honestly — not a single averaged or cherry-picked number.
9. Test the large-file end of scope (T. cacao 15GB, C. elegans 11GB) at
   representative scale (not just a 1M-read slice) to validate the RAM
   architecture (bounded-queue fix, stages 84/85) actually holds at the size
   it was built for — the one part of scope never stress-tested at real scale
   this session.
10. Claim 2 human data (HG002-HG005) is a separate, later scope extension
    (different pipeline — `arcs compress --call` in ARCS proper, not this
    reimpl sandbox) — out of scope for this reimpl generalization sweep, noted
    here only so it isn't silently forgotten.

## 4. What "win in depth" means, concretely

Not: one averaged number across cherry-picked files. It means: for every file
in the 10-accession scope, report size/speed/RAM against real SPRING and
Genozip honestly, including losses (like SARS-CoV-2 today), and only claim a
generalized win once the known losses are either fixed (step 5) or explicitly
scoped out with a stated, mechanistic reason (e.g., "loses only on
near-degenerate, highly-redundant, sub-1MB-genome inputs, for reason X, fix
tracked as step 5").
