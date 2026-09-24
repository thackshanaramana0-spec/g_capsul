# Third-party notices

This repository's own code is MIT licensed (see [`LICENSE`](LICENSE)). It
vendors three third-party components under `thirdparty/`. Their licences
govern those directories and are reproduced there.

| component | licence | where |
|---|---|---|
| `thirdparty/fse/` | Finite State Entropy, Yann Collet. BSD 2-Clause. | [`thirdparty/fse/LICENSE`](thirdparty/fse/LICENSE) |
| `thirdparty/htscodecs/` | htscodecs (fqzcomp quality coder), Genome Research Ltd. BSD 3-Clause. | [`thirdparty/htscodecs/LICENSE.md`](thirdparty/htscodecs/LICENSE.md) |
| `thirdparty/ppmd/` | Ppmd7 from the LZMA SDK, Igor Pavlov. Public domain. | No separate licence file ships with these extracted sources — the LZMA SDK is placed in the public domain by its author. |

`liblzma` is a runtime dependency, linked with `-llzma`, and is public domain.

**Not vendored, for the record**: PgRC2 (GPL-3) is used only as an external
comparison binary, cloned separately. No PgRC2 source is included in or
linked into this project, so its GPL-3 terms do not attach to this
repository.
