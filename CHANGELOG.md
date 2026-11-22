# Changelog

## 1.2.0

Release prepared for the manuscript submitted to *Computational Biology and Chemistry*.

### Added
- `scripts/build_scripts/verify_roundtrip.sh`: ordered, byte-for-byte FASTQ round trip.
- FM-index over the reads that the pseudogenome cannot represent (`include/caps_fmi.h`), built
  by `capsule_decode index`, giving exact-match recall 1.0000 on the four GIAB individuals.
  Tests: `tests/fmi_selftest.cpp` and `tests/fmi_end_to_end.py`, run by `ctest`.
- Multi-individual locus-retrieval and exact-match benchmarks (phase 4 of
  `scripts/benchmark_scripts/benchmark_1_run.sh`; results in `results/claim3/`).
- `docs/configuration_reference.md`: every environment variable and positional argument.
- `docs/extras/ENGINEERING_NOTES_*.md`: the measurement history behind design decisions,
  moved out of the source comments.
- README sections: Try it on your own FASTQ, Input requirements, Reproduce the paper tables,
  Archive format, Known issues.

### Changed
- Scripts reorganized into `scripts/build_scripts`, `scripts/benchmark_scripts`,
  `scripts/claim2_scripts` and `scripts/claim3_scripts`; `include/coders_pgrc.h` renamed to
  `include/coders_entropy.h`.
- Source comments and documentation revised; executable code is unchanged by that revision.
- `THIRD_PARTY_NOTICES.md` lists the published methods whose techniques are reimplemented.

### Fixed
- The encoder now refuses, before writing any archive, input that it would otherwise encode
  wrongly: CRLF line endings, bases other than A, C, G, T, N, an N after base 256 or more than
  255 N in a read, truncated or malformed records, and mixed line-3 content (with names stored).
- `capsule_decode` exits with status 2 when the output directory does not exist.
- `scripts/benchmark_scripts/check_input_scope.py` located `src/encoder.cpp` one directory too
  high, so `encode_adaptive.sh` refused every input; it now finds the source and reports the
  same refusals as the encoder.
- The top-level `Makefile` is a wrapper over the CMake build.
- Signed left shift in the zigzag delta encoding of two position streams (output unchanged).

### Known issues
See the Known issues section of `README.md`.

## 1.1.1 and earlier

See the GitHub releases page: https://github.com/thackshanaramana0-spec/g_capsul/releases
