# Claim 1 (COMPACT) — final verdict

Same five-part structure as `CLAIM2_FINAL_VERDICT.md` and
`CLAIM3_LOCKED.md`. Per instruction: the missing 15th-dataset run does not
block this verdict — it is named plainly in every section it touches,
exactly as the full-scale gap was named for Claim 2 without blocking that
verdict either.

---

## 1. The research idea — what it is, and that it was actually done

**Idea:** build a pseudogenome from the reads at compress time (greedy
suffix-prefix overlap chaining → pigeonhole mapping of the remainder →
second-region assembly of leftovers → MEM self-match to remove residual
redundancy), then code every stream — literal bases, match references, read
order, mismatches — against its own measured entropy bound, not a generic
compressor.

**Built, not proposed:** `stages/106_inprocess.cpp` is the shipped,
in-process encoder; the `stages/01...106` progression is not decoration —
each numbered stage is a real, previously-tried design point, several later
contradicted by measurement and kept in the record rather than deleted
(`CLAUDE.md` §7).

**Is it real, or overlapping existing work?** PgRC2 is the direct
architectural relative (pseudogenome-based, same general family) and is the
project's primary comparison target, run from its own real source, never
approximated. SPRING and Genozip are the general-purpose FASTQ compressors
this claim is benchmarked against for the SPRING/Genozip 14/15-dataset
comparison. The claim is not "no one has built a pseudogenome compressor" —
PgRC2 already exists — it is "this independent implementation wins on
sequence content against PgRC2 and wins whole-file against SPRING/Genozip,"
which is a real, falsifiable, and measured comparison, not an overlap
concern.

## 2. The result — dominant, on the datasets actually run

| comparison | scope | result |
|---|---|---|
| vs SPRING/Genozip | 14 of 15 datasets (Utricularia gibba not yet run) | **14/14 whole-file wins vs both** |
| vs PgRC2 (sequence only) | 7 datasets | **+1.88% aggregate, 6 wins, 1 loss** (S. acidocaldarius, -0.83%) |

This is dominant in the sense that matters most for a compression claim:
wins are not narrow or contested — every stream in the losing dataset was
independently checked against its own information-theoretic bound
(`claim1_locked_result.md`: read order at 0.996× its bound, mismatch
symbols 40% ahead of PgRC2) rather than the aggregate number being taken on
faith. The one loss is disclosed with its exact margin, not folded into an
average that hides it.

**What is not yet true:** the SPRING/Genozip comparison is 14/15, not
15/15. Utricularia gibba's `.sra` file is on disk (`/tmp/newdl/SRR10676752/`,
7.3 GB) but has not been converted or run. This does not change the
direction of the claim — 14 consecutive wins is already strong evidence —
but "14/15" and "15/15" are different sentences in a paper, and only one of
them is currently true.

## 3. Bugs found, and what running the numbers established

This claim's history already contains the deepest bug-finding track record
of the three, from *before* this session:

1. **The archive was incomplete, not just imperfect** (`a81f55c`) — `refc::encode`
   stored a reference's source but not its destination gap, length, or RC
   flag, so the shipped archive could not place a single reference. Found
   only because `verify_lossless.sh`'s own limitation (checking dumped
   streams, not the real archive) was itself questioned and fixed. The
   aggregate margin fell from a previously-claimed +4.65% to the real
   +1.90% once this was corrected — a real, disclosed retraction.
2. **Four silent data-loss bugs** (`CLAUDE.md` §6.3): a 1-byte mismatch
   position field silently corrupting reads over 256bp, orphaned
   unique-read mismatches desynchronizing the adaptive coder, a
   reverse-complement contained-read indexing bug, and an FSE-RLE decode
   path that silently zeroed any constant non-zero stream. All four were
   found by actually decoding archives and diffing against the original
   file — dataset by dataset — not by trusting a passing size table.

**What this establishes:** identical to the standard set for Claims 2 and
3 — a claim with zero bugs found on its first deep audit would be the
suspicious outcome. Claim 1 has the strongest version of this evidence
because the bugs were found *before* this session's own checklist pass,
by the same kind of skepticism this pass is applying to Claims 2 and 3 now.

**This session's own contribution to Claim 1** was narrower and
lower-risk than the caller/decoder work: one results file
(`/tmp/allph/r.csv`, the exact CSV `PAPER_DRAFT_CLAIM1.md` cites) was found
living outside the repo, at risk of being lost on reboot, and was committed
to `results/phase_a/allphases_14dataset.csv`. No code bug was found in this
pass, and none was expected — the code had already been through the kind of
audit this session performed on Claims 2 and 3 for the first time.

## 4. Industrial-grade checklist

See **`docs/INDUSTRIAL_CHECKLIST_CLAIM1.md`**. One concrete fix (the
at-risk CSV, now committed). Claim 1 enters this checklist in the strongest
position of the three — it already had a real round-trip test before this
session — so most rows were already ✅ rather than needing new work.

## 5. Research/academic-grade checklist

See **`docs/RESEARCH_CHECKLIST_CLAIM1.md`**. Every 🔴-critical row is ✅
except three named, scoped gaps: the 15th dataset (Reproducibility), data
integrity checksums (shared gap with Claims 2/3), and PgRC2's version not
being pinned (Dependency versions). None of these are code defects.

---

## Final verdict

**Claim 1 is the most evidence-complete of the three claims and is ready to
lock on everything it has actually measured — with one explicitly named
exception that does not change its direction.**

- The idea is implemented with the deepest bug-finding history in the
  project, self-applied before any external checklist asked for it.
- The result dominates on real, non-toy data (14/14 vs SPRING/Genozip,
  +1.88% vs PgRC2), with the one loss disclosed at its exact margin.
- One real provenance risk (a results CSV living only in `/tmp`) was found
  and fixed this session.
- **The 15th dataset (Utricularia gibba) has not been run.** Per
  instruction, this does not block the verdict — but it means any paper
  claiming "15/15" rather than "14/15" would be stating something not yet
  measured. The honest, currently-true sentence is 14/14 wins on the
  datasets actually run, with the 15th named as outstanding, not hidden.

**Recommended next step, stated plainly:** convert the Utricularia gibba
`.sra` already on disk and run it through the same SPRING/Genozip
comparison as the other 14. Everything else in this claim is ready now.
