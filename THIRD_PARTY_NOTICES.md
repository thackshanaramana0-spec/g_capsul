# Third-party notices

This repository's own code is Apache License 2.0 licensed (see [`LICENSE`](LICENSE)). It
vendors three third-party components under `thirdparty/`. Their licences
govern those directories and are reproduced there.

| component | licence | where |
|---|---|---|
| `thirdparty/fse/` | Finite State Entropy, Yann Collet. BSD 2-Clause. | [`thirdparty/fse/LICENSE`](thirdparty/fse/LICENSE) |
| `thirdparty/htscodecs/` | htscodecs (fqzcomp quality coder), Genome Research Ltd. BSD 3-Clause. | [`thirdparty/htscodecs/LICENSE.md`](thirdparty/htscodecs/LICENSE.md) |
| `thirdparty/ppmd/` | Ppmd7 from the LZMA SDK, Igor Pavlov. Public domain. | No separate licence file ships with these extracted sources — the LZMA SDK is placed in the public domain by its author. |

`liblzma` is a runtime dependency, linked with `-llzma`, and is public domain.

**Not vendored**: PgRC2 (GPL-3) is used as an external comparison binary, cloned
separately. No PgRC2 source file is included in or linked into this project.

## Algorithmic provenance

The components below are not vendored. They are this project's own code, written to
implement techniques published by, or studied in the source code of, the projects named.
The source comments at each location give the specific technique and, where one was
consulted, the file and line of the original.

| technique in G_CAPSUL | origin | where in this repository |
|---|---|---|
| Pseudogenome assembly by sorted suffix-prefix overlap, application of exact matches after sorting by destination, mismatch positions coded as reverse offsets, mismatch-position bucketing and transpose, mismatch-count zero-flag split, LZMA position-coder parameters | PgRC2 (Kowalski and Grabowski) | `src/encoder.cpp`, `include/coders_inproc.h`. Background in `docs/extras/REIMPL_NOTES.md` |
| Sampled seed index for pigeonhole read mapping, using the copMEM sampling bound (seed length plus stride minus one) | copMEM (Grabowski and Bieniecki) | `src/encoder.cpp` (pigeonhole mapping). The MEM self-match stage does not use copMEM sampling; it uses an exact stride, as its comment explains |
| Read-identifier tokenization against the previous header | SPRING (Chandak et al.) | `include/names_coder.h` |
| Bubble and indel discovery on the de Bruijn graph, indel positional-ambiguity check (`checkRepeatSize`) | DiscoSNP++ / kissnp2 (AGPL), reimplemented from the published algorithm | `include/caps_caller.h` |
| k-mer counting spilled to disk by key range, memory budget | GATB / KMC2 | `include/caps_caller.h` |
| Positional clustering of deviating reads as an indel channel | eBWT2SNP (inspired by, not a reimplementation) | `include/caps_caller.h` |
| Quality coding | fqzcomp, vendored through htscodecs (see the table above) | `include/quality_coder.h` |
