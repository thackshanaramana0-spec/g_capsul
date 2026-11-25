# Code-to-table mapping (Claim 1)

This file maps the Claim 1 results (T1.1 archive size, T1.2 wall time) to the code that produces
them and lists the source file of every number cited in this folder. Line references are to
`src/encoder.cpp` (4,989 lines), `include/coders_inproc.h` (938), `include/coders_entropy.h` (240)
and `include/seqpar_core.h` (324). They drift as the files change, so search for the stage name.
The mechanism statements in [`mechanism_insight_claim1.md`](mechanism_insight_claim1.md) rest on
them. For the full description of the encoder see
[`../technical_architecture.md`](../technical_architecture.md).

## Pipeline stages in `src/encoder.cpp` (the `phase(...)` markers)

1. **`greedy-sweep`** (line 2511) is the pseudogenome construction, by suffix-prefix overlap
   chaining. Its cost explains the T1.2 result, and its output (the placement data) is what Claims
   2 and 3 reuse.
2. **`pre-call` and `call`** (lines 3251 and 3368) form the inline variant-calling path, run only
   when `CAPS_CALL=1` (the caller is `include/caps_caller.h`, see the Claim 2 code mapping).
   The T3.1 investigation found that `CAPS_CALL` costs about 20 times more than `CAPS_SPANS` alone,
   because it runs this stage inline.
3. **`MEM+positions`** (line 3583) builds the reference and position streams. A log line here
   compares the coded position size with PgRC2's ("PgRC2 pays 683,370 B coded for its reads-list
   offsets").
4. **`assemble+map`** (line 3672) maps the remaining reads by pigeonhole after chaining (`MINMEM`,
   `MEMSEED` and `STEP` are logged here).
5. **`pre-coding` and `CODING`** (lines 4940 and 4984) are the stream-coding stage, where the
   per-stream selector below runs.

## Stream coder selector (`include/coders_inproc.h`, line 560)

The comment reads "universal per-stream selector over the REAL coder set". Each stream (literal
bases, match references, mismatches and so on) is coded by several candidate coders and the
smallest result is kept, in place of one fixed generic compressor. `include/coders_entropy.h`
provides the PPMd7, FSE and range-coder implementations that are among the candidates.
`include/seqpar_core.h` is the DNA coder shared by the standalone binary and the in-process path.

## Losslessness verification

An earlier `verify_lossless.sh` decoded dumped intermediate streams and not the shipped archive.
That check hid four silent data-loss bugs: a mismatch-position byte overflow for reads above
256 bp, two "orphaned unique read" reconstruction bugs, and an FSE-RLE constant-stream decode
bug. They were found only by decoding the shipped archives and comparing them with the source. This
is the reason the T1.1 method note states that every archive was independently decoded and checked
against its input.

---

# Source file of every number

| Number | Source |
|---|---|
| T1.1 archive sizes and aggregates | `results/claim1/claim1_T1.1_T1.2.csv` |
| T1.2 wall times | same file, columns `compress_s` and `decompress_s` |
| 19/19, -6.03%, -43.26%, 57/57 lossless | the development notes, and recomputed from the CSV above |
| PgRC2 sequence-only comparison (internal test only, not published) (+1.88%, 6 wins and 1 loss, about 1.7x and 2.5x) | the development notes, citing the Claim 1 final-verdict note and the headroom analysis |
| The four silent data-loss bugs | the development notes, section 6.3 |
| `CAPS_SPANS` adds zero archive bytes | [`../claim3/t31_export_claim3.md`](../claim3/t31_export_claim3.md) |

## What was re-derived and what is cited

The T1.1 aggregate (19/19, -6.03%, -43.26%, 57/57 lossless) was recomputed from the CSV. The T1.2
count of 0 of 19 on both axes was computed the same way. The manuscript's "bold marks fastest per
row" is consistent with it but does not state the count.

The PgRC2 sequence-only numbers (+1.88%, 6 and 1, about 1.7x and 2.5x(internal test only, not published)) and the descriptions of the
four silent data-loss bugs are cited from the development notes. They were not re-run or
re-confirmed by decoding archives again.
