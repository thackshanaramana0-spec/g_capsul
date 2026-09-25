---
Date: 2026-09-19
Title: Claim 1 — Code-to-Table Mapping and Result File Index
Purpose: Exact stage/function references in src/encoder.cpp and
  the shared coder headers for T1.1/T1.2, and the exact source file for
  every number cited in this folder. Traced directly from the actual repo
  this session.
When to refer to this file: Verifying a claim in another file of this
  folder against the real code; before citing any number, to find its
  exact source file.
Keywords: code mapping, encoder.cpp, coders_inproc.h, seqpar_core.h,
  pipeline stages, selector, benchmark results, file index, silent data
  loss bugs
---

# Code-to-table mapping

Line/stage numbers from `src/encoder.cpp` (4,961 lines),
`include/coders_inproc.h` (931 lines), `include/coders_pgrc.h` (213 lines),
`include/seqpar_core.h` (357 lines), as present in the repo this session.
Traced by targeted grep/read, not a full line-by-line read of all 6,462
combined lines — sufficient to ground the mechanism claims in
`mechanism_insight_claim1.md`, not a claim to have read every line.

## Pipeline stages, `src/encoder.cpp` — confirmed via `phase(...)` markers

1. **`greedy-sweep`** (line 2331) — the pseudogenome construction stage:
   suffix-prefix overlap chaining. This is the stage whose cost is the
   direct explanation for T1.2's 0/19 speed result, and whose output
   (placement data) is what Claims 2/3 reuse for free.
2. **`pre-call` / `call`** (lines 3130, 3247) — the inline variant-calling
   path, only run when `CAPS_CALL=1` (see Claim 2's code mapping file for
   the caller itself, `include/caps_caller.h`). Confirmed this session
   (Claim 3's T3.1 investigation) that `CAPS_CALL` costs ~20x more than
   `CAPS_SPANS` alone, because it runs this full stage inline.
3. **`MEM+positions`** (line 3513) — reference/position stream construction;
   comment at this line explicitly compares this project's coded position
   size against PgRC2's ("PgRC2 pays 683,370 B coded for its reads-list
   offsets"), direct evidence the position-encoding comparison to PgRC2 is
   a real, code-level concern in this project, not an afterthought.
4. **`assemble+map`** (line 3602) — pigeonhole mapping of the remainder
   after chaining (`MINMEM`/`MEMSEED`/`STEP` logged at this line).
5. **`pre-coding` / `CODING`** (lines 4912, 4956) — the stream-coding stage,
   where the per-stream selector (below) runs.

## Stream coder selector — `include/coders_inproc.h:556`

Comment: *"universal per-stream selector over the REAL coder set"* — this
is the mechanism behind T1.1's size result: each stream (literal bases,
match references, read order, mismatches, etc.) is coded by multiple
candidate coders and the smallest result kept, rather than a fixed
generic-compressor pass. `include/coders_pgrc.h` provides PPMd7/FSE/range
coder implementations that are among the candidates this selector chooses
from. `include/seqpar_core.h` is the shared DNA coder used identically by
both the standalone binary and the in-process path (per its own module
comment, confirmed at conversation start of this session, not re-verified
line-by-line this pass).

## Losslessness verification discipline — the real history behind the "LOSSLESS" column

Confirmed against `BTR_NOTES.md` §6.1/6.3 (read this session): this project
previously shipped an incomplete verification path (`verify_lossless.sh`
decoded dumped intermediate streams, not the actual archive) that hid four
real silent data-loss bugs — a mismatch-position byte overflow above 256bp
reads, two "orphaned unique read" reconstruction bugs, and an FSE-RLE
constant-stream decode bug. All four were found only by actually decoding
the shipped archives and diffing against source, not by trusting the dump-
level check. **This is the concrete, real reason T1.1's methodology note
("every archive was independently decoded and verified byte-identical")
is stated as carefully as it is** — it is a corrected practice, not a
default one.

---

# Result file index — exact source for every number cited in this folder

| Number / claim | File |
|---|---|
| T1.1 (archive size, all rows, aggregate) | `results/claim1/claim1_T1.1_T1.2.csv` |
| T1.2 (wall time, all rows) | same CSV, `compress_s`/`decompress_s` columns |
| 19/19, −6.03%, −43.26%, 57/57 LOSSLESS headline | `BTR_NOTES.md` (top section); independently re-derived this session from the CSV above |
| PgRC2 sequence-only comparison (+1.88%, 6 wins/1 loss, ~1.7x/2.5x) | `BTR_NOTES.md` lines ~248-268, citing `refer_paper_docs/CLAIM1_FINAL_VERDICT.md` §42 and `refer_paper_docs/FINAL_HEADROOM.md` |
| The four corrected silent data-loss bugs | `BTR_NOTES.md` §6.3 |
| Full session-by-session Claim 1 development history | `refer_paper_docs/CLAIM1_FINAL_VERDICT.md`, `refer_paper_docs/COMPACT_HEADROOM.md`, `refer_paper_docs/PAPER_DRAFT_CLAIM1.md` |
| `CAPS_SPANS` free-byproduct verification (zero archive-byte cost) | This session's Claim 3 work — see `refer_paper_docs/claim3/t31_export_claim3.md` and `PIPELINE.md`'s flags table |

## What was verified this session vs cited from earlier work

**Verified this session, directly, by re-deriving from raw data**: the full
T1.1 aggregate (19/19, −6.03%, −43.26%, 57/57 LOSSLESS) was independently
recomputed from the CSV, not just read from `BTR_NOTES.md`. T1.2's 0/19-both-
axes result was computed the same way, discovered fresh this session (not
previously stated anywhere this explicitly as "0/19" — the manuscript's own
"bold marks fastest per row" phrasing is consistent with it but does not
state the aggregate count directly).

**Not independently re-verified this session** (cited from prior,
documented work): the PgRC2 sequence-only comparison numbers (+1.88%, 6/1,
~1.7x/2.5x) — read from `BTR_NOTES.md`'s own citation of
`refer_paper_docs/CLAIM1_FINAL_VERDICT.md`/`refer_paper_docs/FINAL_HEADROOM.md`, not re-run this
session. The four silent-data-loss bug descriptions — read from `BTR_NOTES.md`
§6.3, not independently re-confirmed by re-decoding archives this session.
