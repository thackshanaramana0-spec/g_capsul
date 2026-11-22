# Code-to-table mapping

> Functions and files in src/encoder.cpp and the shared coder headers behind T1.1 and T1.2, and the source file of every number cited in this folder.

Line and stage references are to `src/encoder.cpp` (about 4,950 lines),
`include/coders_inproc.h` (about 940 lines), `include/coders_entropy.h` (about 240 lines) and
`include/seqpar_core.h` (about 320 lines). They were traced by targeted search and reading, not by a
full line-by-line read, and they ground the mechanism statements in
`mechanism_insight_claim1.md`. Line numbers drift as the files change, so search for the stage
name.

## Pipeline stages, `src/encoder.cpp` — confirmed via `phase(...)` markers

1. **`greedy-sweep`** (line 2331) — the pseudogenome construction stage:
   suffix-prefix overlap chaining. This is the stage whose cost is the
   direct explanation for T1.2's 0/19 speed result, and whose output
   (placement data) is what Claims 2/3 reuse for free.
2. **`pre-call` / `call`** (lines 3130, 3247) — the inline variant-calling
   path, only run when `CAPS_CALL=1` (see Claim 2's code mapping file for
   the caller itself, `include/caps_caller.h`). Confirmed
   (Claim 3's T3.1 investigation) that `CAPS_CALL` costs ~20x more than
   `CAPS_SPANS` alone, because it runs this full stage inline.
3. **`MEM+positions`** (line 3513) — reference/position stream construction;
   comment at this line explicitly compares this project's coded position
   size against PgRC2's ("PgRC2 pays 683,370 B coded for its reads-list
   offsets"), direct evidence the position-encoding comparison to PgRC2 is
   a code-level concern in this project, not an afterthought.
4. **`assemble+map`** (line 3602) — pigeonhole mapping of the remainder
   after chaining (`MINMEM`/`MEMSEED`/`STEP` logged at this line).
5. **`pre-coding` / `CODING`** (lines 4912, 4956) — the stream-coding stage,
   where the per-stream selector (below) runs.

## Stream coder selector — `include/coders_inproc.h:556`

Comment: *"universal per-stream selector over the REAL coder set"* — this
is the mechanism behind T1.1's size result: each stream (literal bases,
match references, read order, mismatches, etc.) is coded by multiple
candidate coders and the smallest result kept, rather than a fixed
generic-compressor pass. `include/coders_entropy.h` provides PPMd7/FSE/range
coder implementations that are among the candidates this selector chooses
from. `include/seqpar_core.h` is the shared DNA coder used identically by
both the standalone binary and the in-process path (per its own module
comment, not re-verified line by line).

## Losslessness verification discipline — the real history behind the "LOSSLESS" column

Confirmed against the development notes §6.1/6.3 (read): this project
previously shipped an incomplete verification path (`verify_lossless.sh`
decoded dumped intermediate streams, not the actual archive) that hid four
real silent data-loss bugs — a mismatch-position byte overflow above 256bp
reads, two "orphaned unique read" reconstruction bugs, and an FSE-RLE
constant-stream decode bug. All four were found only by actually decoding
the shipped archives and diffing against source, not by trusting the dump-
level check. **This is the concrete, real reason T1.1's methodology note
("every archive was independently decoded and checked against its input")
is stated as carefully as it is** — it is a corrected practice, not a
default one.

---

# Result file index — exact source for every number cited in this folder

| Number / claim | File |
|---|---|
| T1.1 (archive size, all rows, aggregate) | `results/claim1/claim1_T1.1_T1.2.csv` |
| T1.2 (wall time, all rows) | same CSV, `compress_s`/`decompress_s` columns |
| 19/19, −6.03%, −43.26%, 57/57 LOSSLESS headline | the development notes (top section); independently re-derived from the CSV above |
| PgRC2 sequence-only comparison (+1.88%, 6 wins/1 loss, ~1.7x/2.5x) | the development notes lines ~248-268, citing the Claim 1 final-verdict note §42 and the final-headroom note |
| The four corrected silent data-loss bugs | the development notes §6.3 |
| Claim 1 development history | the Claim 1 final-verdict note, the compactness-headroom note, the Claim 1 paper draft |
| `CAPS_SPANS` free-byproduct verification (zero archive-byte cost) | The Claim 3 work — see `docs/claim3/t31_export_claim3.md` and the pipeline notes' flags table |

## What was verified vs cited from earlier work

**Verified, directly, by re-deriving from raw data**: the full
T1.1 aggregate (19/19, −6.03%, −43.26%, 57/57 LOSSLESS) was independently
recomputed from the CSV, not just read from the development notes. T1.2's 0/19-both-
axes result was computed the same way, newly derived (not
previously stated this explicitly as "0/19" — the manuscript's own
"bold marks fastest per row" phrasing is consistent with it but does not
state the aggregate count directly).

**Not independently re-verified** (cited from prior,
documented work): the PgRC2 sequence-only comparison numbers (+1.88%, 6/1,
~1.7x/2.5x) — read from the development notes' own citation of
the Claim 1 final-verdict note/the final-headroom note, not re-run. The four silent-data-loss bug descriptions — read from the development notes
§6.3, not independently re-confirmed by re-decoding archives.
