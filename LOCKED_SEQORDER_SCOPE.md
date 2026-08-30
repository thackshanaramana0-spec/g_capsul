# Locked scope: sequence + read-order only

This is the exact, verified scope PgRC2 measures — confirmed directly from
its real source (`/tmp/pgrc_study/PgRC`), not assumed. PgRC2's `-q` flag only
uses quality as an internal heuristic for read-division decisions; there is
no quality-compression or header-compression code path anywhere in its
encoder. Its one archive-size number IS the complete sequence+order cost.

**Binary: `100_locked_seqorder.cpp`.** Built by merging three previously
separate, diverged pieces of this codebase — none of them alone was correct
or complete:

- Base: `46_position_stream.cpp` (has the positions+strand stream, missing
  from the 87→98 lineage entirely)
- Ported in: the read-length fix from `90_readlen_u16.cpp` (46's own base
  still had the old `uint8_t rlen` bug — silently dropped reads >255 bases,
  same bug class as below, independently present on this branch)
- Ported in: the mismatch-symbol dump from `47_mismatch_coder.cpp`
- New: N-containing-read storage (previously silently dropped, see below)

## The 6 real layers

A decoder needs ALL SIX of these, losslessly, to reconstruct every read's
sequence and its original file position. Any subset is an incomplete,
misleadingly favorable number — this was tested and caught directly this
session (an earlier 3-layer-only total made E. coli look like a 39.3% win;
the real, complete number is 18.6%).

| # | Layer | File(s) | Coder | What it's for |
|---|-------|---------|-------|----------------|
| 1 | Sequence (literal) | `literal.txt` | `seqpar` | Survivor bases not covered by self-match or mapping |
| 2 | Read order | `perm.u32` | `permcoder` | Permutation restoring original file order |
| 3 | MEM references | `mem_triples.bin` | `refcoder` | Self-match (gap, source, length) triples from assembly |
| 4 | Positions + strand | `pos_delta.bin`, `pos_strand.bin` | `xz -9` (no dedicated coder built yet) | Where each read starts in the pg and which strand, for every read in the pg |
| 5 | Mismatch symbols | `mm_ref.bin`, `mm_obs.bin` | `mmcoder` | The actual substituted bases for reads placed via mapping, not exact MEM match |
| 6 | N-containing reads | `n_reads.txt`, `n_indices.bin` | `seqpar` (reads) + `xz -9` (indices) | Reads containing 'N' — cannot be 2-bit packed, need separate lossless storage |

## Layer 6 — the real bug this session found and fixed

Every stage from `87_atomic_dump.cpp` through `98_o1_removal.cpp` silently
`continue`'d past any read containing 'N', with **no storage anywhere** —
real, undisclosed data loss, not a missing-cost accounting gap like layers
4/5. It went undetected because every "round trip: VERIFIED" claim in this
project has only ever been coder-level (does the entropy coder correctly
decode its own dumped intermediate file back to itself) — never a genuine
original-FASTQ-to-decoded-FASTQ comparison. Two binaries that both drop the
same reads still match each other byte-for-byte; that was mistaken for
correctness against the source file, which had never actually been tested.

Fixed in `100_locked_seqorder.cpp`: N-reads are written, full raw sequence,
to `n_reads.txt` (one per line — can't 2-bit-pack a 5th symbol), with their
original 0-based read index to `n_indices.bin` (raw LE uint32), gated on
`DUMP_LIT`. **Directly verified against the true original file** (not just
coder round-trip) on E. coli: all 692 N-reads' recorded index + stored
sequence match the real original FASTQ exactly, byte for byte — the first
genuine original-to-decoded check run in this project, done specifically
because this bug was found.

## Running it

```bash
BEST=/tmp/best100   # rebuild: g++ -O3 -march=native -fopenmp -pthread -std=c++17 100_locked_seqorder.cpp -o /tmp/best100
INPUT=your.fq

rm -f literal.txt perm.u32 mem_triples.bin pos_delta.bin pos_strand.bin \
      mm_ref.bin mm_obs.bin mm_pos.bin mm_ctx3.bin n_reads.txt n_indices.bin

DUMP_LIT=1 DUMP_PERM=1 DUMP_MM=1 "$BEST" "$INPUT" 3 40 16 22 16 16 1 24 64 1
# note PG_LEN and "MEM main ... -> X" (that X is main_pg_end) from stderr

ENCODE_ONLY=1 /tmp/seqpar literal.txt 12 12            # layer 1
/tmp/permcoder perm.u32                                  # layer 2
/tmp/refcoder mem_triples.bin <PG_LEN> <main_pg_end>      # layer 3
xz -9 -k -c pos_delta.bin  | wc -c                        # layer 4a
xz -9 -k -c pos_strand.bin | wc -c                        # layer 4b
/tmp/mmcoder mm_ref.bin mm_obs.bin                        # layer 5
/tmp/seqpar n_reads.txt 1 1                               # layer 6a
xz -9 -k -c n_indices.bin  | wc -c                        # layer 6b
```

Sum all six `coded=`/byte-count numbers. That total is directly comparable
to PgRC2's real archive size (`/tmp/pgrc_study/PgRC/build/PgRC -o -t 12 -i
your.fq out.pgrc`, then `stat -c%s out.pgrc`) — same scope, both sides,
verified from real code on both sides, not assumed.

**Not yet built:** a single wrapper script summing all six automatically
(next real step, not done yet — do not claim it exists until it does).
No dedicated coder exists yet for positions+strand or N-read indices (using
`xz -9` as a real, correct, un-optimized placeholder) — a real, disclosed
opportunity for a smaller number later, not a correctness gap.

## Results so far (real, all 6 layers, verified)

| Dataset | Ours (6 layers) | PgRC2 (real archive) | Result |
|---|---|---|---|
| E. coli (1M reads) | 4,996,977 | 6,139,045 | **WIN 18.6%** |
| yeast (1M reads, cross-check only) | — (not yet run with all 6 layers) | 7,063,585 | matches historical 7,063,459 within 126 B, validates PgRC2-side methodology |

Only E. coli has the complete, real, all-6-layer number so far. The other 9
locked-scope datasets still need this exact same treatment before any
"we win N/10" claim on this scope is real.
