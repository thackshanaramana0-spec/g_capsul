# CAPSULE — the archive format

**C**ompact, **A**ddressable, **P**seudogenome-**S**tructured, **U**nified
**L**ossless **E**ncoder.

## Why this exists

Losslessness was verified by decoding the RAW intermediate streams the encoder
dumps alongside the archive. That is a real round trip of the algorithm --
assembly, mapping, mismatches and read order are all reconstructed and compared
byte-for-byte against the original FASTQ -- but it never exercised the entropy
layer, and nothing could read the archive back at all. A compressor whose
archive no decoder reads is not a format.

## Layout

    magic      "CAPSULE\0"        8 B
    version    uint16             2 B
    pg_len     uint64             8 B
    main_pg_end uint64            8 B
    n_streams  uint16             2 B
    per stream:
        namelen  uint8
        name     bytes
        length   uint64
        payload  bytes

Streams are self-identifying, so the decoder matches by NAME and adding a stream
(names, quality) is additive rather than a format break. The pseudogenome
parameters live in the header instead of a side file.

Each `best_encode` payload is `[method:1][raw_len:8][data]`. The raw length is
required because FSE, PPMd and the range coder cannot recover the uncompressed
size from their own output; it is written uniformly so one decode path serves
every method. A constant stream is `[0xC0][count:8][value]`.

Total format overhead: **210 bytes** (H. salinarum 2,606,073 -> 2,606,283,
+0.008%). All previously reported sizes shift by this amount.

## Status

Implemented and verified against the encoder's own streams:

    pos_abs        1,842,004 B   IDENTICAL   (byte-planes + B1 raw-plane mask)
    read_lengths     921,002 B   IDENTICAL   (constant-stream encoding)
    pos_strand        57,563 B   IDENTICAL
    n_pos, n_indices, n_cnt       IDENTICAL

That covers LZMA, the byte-plane inverse, the raw-plane mask, the constant
encoding and the container walk: **6 identical, 0 differ**.

PPMd, FSE and the adaptive range coder now have inverses written
(`ppmd_decode`, `fse_decode`, `range_decode` in `include/coders_pgrc.h`), but
they are only exercised when a stream selects those methods -- H. salinarum
chose LZMA variants for every stream above.

## Not yet implemented

Four streams still need their coder inverse:

    literal       seq_encode_mem   (the DNA coder)
    mem_triples   refc::encode
    mm_sym        mmc::encode
    mm_pos        mmpos_encode_buckets when the bucketed form wins

Until those exist, `capsule d` cannot reconstruct reads from the archive alone;
the end-to-end check still routes through the dumped streams. The container, the
header and the general-purpose coders are done.

## Commands

    capsule c reads.fq out.capsule      (today: the stage-106 encoder)
    capsule d out.capsule outdir        (today: emits raw streams + pg_params)
