# Technical architecture

This document describes how the encoder builds the pseudogenome and how the archive is laid out.
It is independent of any single claim. The implementation is in `src/encoder.cpp` and
`src/decoder.cpp`, and every switch mentioned here is listed in
[`configuration_reference.md`](configuration_reference.md).

## Pipeline overview

| Stage | Purpose | Result |
|---|---|---|
| 1. Input check | Refuse unsupported input before any work is done. | A clean error, or a validated read set. |
| 2. Duplicate collapse | Merge exact duplicate reads when they are frequent. | A unique read set and a map back to the original order. |
| 3. Overlap chaining | Join reads whose suffix and prefix overlap exactly. | The main region of the pseudogenome. |
| 4. Pigeonhole mapping | Place the remaining reads on the existing pseudogenome. | Placements with a small number of mismatches. |
| 5. Second region | Chain the reads that could not be mapped. | A second pseudogenome region. |
| 6. MEM self-match | Replace repeated stretches by references. | A shorter literal stream. |
| 7. Stream coding | Entropy-code each stream with the coder that suits it. | The archive. |

## Container format

An archive is a fixed header followed by self-identifying streams. The header holds an 8-byte magic,
a version, the pseudogenome length and the boundary between the two regions. Each stream stores its
name, its length and its payload. The decoder reads every stream into a map keyed by name instead of
relying on position, so adding a stream cannot shift any other stream's offset. The header costs
about 300 bytes, which is 0.001 to 0.01% of the archives measured (2.6 to 150 MB). The exact layout
is shown in the README, and each stream is described in the table above `main()` in
`src/encoder.cpp`.

## Building the pseudogenome

**Exact overlap chaining.** Reads are packed at two bits per base. A seed index maps fixed-width
k-mers to the reads that contain them. For each read, the encoder looks for another read whose
prefix equals this read's suffix over at least `MINOV` bases, and chains the two into one
contiguous span instead of storing both. Matching is exact. The overlap sweep starts at the full
read length, so that exact duplicate reads, which overlap only at that length, chain without extra
handling. On S. acidocaldarius, where 9.4% of reads are duplicates, this reduced the pseudogenome
from 9,134,100 to 6,157,270 bytes.

**Second chaining pass.** Reads that did not chain in the first pass are set aside, and a second,
less strict sweep tries to attach them to the ends of existing chains.

**Pigeonhole mapping.** Reads that are still unplaced are mapped onto the pseudogenome built so
far. The encoder finds an exact seed hit, checks the read against a bounded mismatch tolerance
(3 by default), and records a reference `(destination, source, length, strand)` instead of storing
the read's bases again.

**Mapping ceiling.** The per-read candidate ceiling (`MAXMAP`) widens automatically with the
fraction of reads left over after the first chaining pass. That fraction is measured on every run,
so the ramp is derived from the input and is not a fitted constant. Its floor lies above the
leftover fraction of every benchmark dataset, so output on normal-coverage data is unchanged
(checked byte for byte). The ramp engages only on inputs, such as low-coverage data, that
measurably need it.

**Second region.** Reads that cannot be mapped are not stored raw. They are chained by their own
instance of the two-pass procedure, which produces a second pseudogenome region after the boundary
of the main one. This turns unmatched leftovers into sequence that compresses well.

**MEM self-match.** With both regions built, a maximal-exact-match pass hashes k-mers and verifies
full matches to find repeated stretches anywhere in the pseudogenome. Each accepted match becomes a
`(destination, source, length, strand)` reference like a mapping placement, and the bases it
covers leave the literal stream. A match is accepted only if its encoded cost is below the literal
cost it replaces (2.0 bits per base). This stage accounts for much of the size reduction.

**Extension past a mismatch.** A match can in principle extend beyond its first mismatch. This
option is implemented but off. On the dataset where it was tested it cost 490,763 bytes more than
it saved, and it remains behind an override flag only to reproduce that measurement.

## Parameter search

`scripts/build_scripts/encode_adaptive.sh` encodes with several `(MAXMAP, MINOV)` candidates and
keeps the smallest archive. Four candidates were chosen from measurements on 14 datasets of 8.7 to
574 MB:

| Search | Mean archive size relative to the 8-point grid | Worst dataset |
|---|---|---|
| 4-point grid (shipped) | +0.031% | +0.47% |
| 1 fixed point | +0.41% | +0.782% |

Widening from four to eight candidates changed the archive size by 0.000% on 10 of the 14
datasets. A single fixed candidate was about thirteen times worse on average and varied widely
between datasets, which is consistent with there being no fixed read-length ratio that gives the
best mapping ceiling.

## Why this matters for the other capabilities

Exact chaining is what places two divergent alleles of a heterozygous site on separate contigs once
their difference exceeds one seed's worth of exactness. The MEM self-match then collapses alleles
that still share enough sequence to be worth referencing. These two stages are the origin of the
effect that variant calling and locus retrieval reconcile, as described in
[`central_insight_and_mechanism.md`](central_insight_and_mechanism.md).
