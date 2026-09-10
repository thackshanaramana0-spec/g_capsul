# industry/ — reproducibility, scope, and safety tooling

This folder holds tooling whose job is not to produce a result but to make
sure every result that IS produced is trustworthy: that the input was in
scope, that a failure fails loudly instead of silently, and that the same
input behaves the same way regardless of machine.

Everything here was built from real, empirically-confirmed defects, not from
a checklist. Each one is documented below with how it was found.

---

## The two-layer design

**Layer 1 — the encoder itself refuses.** `stages/106_inprocess.cpp` contains
an unconditional scope check, near the top of `main()`, before anything else
runs. It streams the input once, computes real statistics (record count,
length distribution, header length), and calls `return 2` with a full
explanation if the input cannot be processed correctly — before any
per-read buffer is touched. **This cannot be bypassed** by skipping a wrapper
script or calling the binary directly; it is not a shell-level convenience,
it is the actual safety mechanism.

**Layer 2 — this folder reports the same verdict in advance, in more
detail.** `check_input_scope.py` is wired into `scripts/encode_adaptive.sh`
(the production entry point) and runs before the expensive multi-candidate
encode sweep starts, so out-of-scope input is refused in seconds with a full
explanation rather than after however long the sweep would otherwise take.
It extracts its bounds — `MAX_READ_LEN`, `nmc::MAXTOK` — **directly from the
encoder's own source at run time**, via regex, rather than retyping the
numbers. If the encoder's bound ever changes, this script's next run picks
it up automatically; if it cannot find the constant, it refuses to state a
verdict rather than check against a guessed or stale number.

Layer 2 existing does not weaken Layer 1. A user who calls the binary
directly, or edits the wrapper, or bypasses this folder entirely, still hits
the encoder's own unconditional gate.

---

## What is genuinely in scope, and why

**Supported:** short-read FASTQ (Illumina-shaped), any read length up to
`MAX_READ_LEN` bases (currently 1023 — see `check_input_scope.py`'s printed
header for the live value), fixed or variable length, forward and
reverse-complement strands, ACGTN alphabet, any read count from one upward,
with or without paired quality scores in standard Phred encoding.

**Not supported, and refused rather than silently mishandled:**

| input | what happens | why |
|---|---|---|
| any read over `MAX_READ_LEN` bases (long-read: Oxford Nanopore, PacBio) | encoder refuses, exit 2, with the exact count and longest length found | every per-read decode path unpacks into a fixed-size stack buffer sized for short-read technology (see `stages/106_inprocess.cpp`'s `MAX_READ_LEN` comment) |
| a header over `nmc::MAXTOK` characters, with `CAPS_NAMES=1` | encoder refuses, exit 2 | the names tokenizer indexes per-token model arrays by token count; past the bound it silently decodes to something other than what was encoded (verified, see below) — sequence and quality are unaffected if `CAPS_NAMES` is unset |
| empty or malformed (non-FASTQ, truncated) input | encoder refuses, exit 2 | previously produced a structurally valid, completely empty archive at exit 0 — success on a file it did nothing with |
| non-ACGTN alphabet (protein, RNA with U, corrupted binary) | **not gated today** — `check_input_scope.py` flags it as a caveat, but the encoder does not refuse | a real, named gap, not a hidden one; see the caveat text this script prints |
| variable-length input | supported, with a measured accuracy caveat | this project's own published results (`paper/LIMITATIONS.md` §2) show the caller's precision is measurably worse on the one variable-length individual in its benchmark than on fixed-length ones — not a refusal, a documented cost |

## What was found empirically, not by inspection

Both defects below were found by actually constructing the input and running
the real binary, because reading the code alone did not surface either one:

1. **Silent total data loss on oversize reads.** `if(b.size()>1023)
   continue;` dropped every oversize read with no message; on a file made
   entirely of such reads, every read vanished and the pipeline still wrote
   a valid, empty archive at exit 0.
2. **A `std::terminate()` crash, exposed while fixing (1).** The
   names/quality background thread (`q_thread`) was started *before* the
   scope check ran; an early `return` skipped joining it, and a joinable
   `std::thread`'s destructor calling `std::terminate()` is undefined
   behavior. This reproduced 100% of the time under plain execution and 0%
   of the time under a debugger — the textbook signature of a
   thread/static-teardown race — which is why it was checked by running the
   exact failing command six times in a row, not by reasoning about it once.
3. **Silent names corruption on pathological headers.** A header producing
   more than `MAXTOK` tokens overran the tokenizer's fixed-size per-token
   arrays' *logical* bound (guarded against overflow, but not against
   producing a different result). Verified by round-tripping a real
   2000-token header through the actual build: it came back byte-different
   on every one of 20 test records.

All three are fixed in `stages/106_inprocess.cpp`, verified byte-identical
against real E. coli data before and after, and covered by dedicated,
repeated (not single-shot) tests in `scripts/run_tests.sh`.

## check_input_scope.py

```bash
python3 industry/check_input_scope.py reads.fq [reads_2.fq ...]
python3 industry/check_input_scope.py --json reads.fq   # machine-readable
```

Exit codes: `0` in scope, `1` in scope with caveats worth reading, `2`
out of scope (the encoder will refuse this file), `3` could not read the
file or could not confirm the encoder's real bounds.

Never trust a hardcoded "this looks like Nanopore" string match — this
script does not do that anywhere. Every verdict is derived from streaming
the actual file and comparing what it finds against the actual encoder's
actual limits, extracted from source.
