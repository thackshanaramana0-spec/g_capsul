# CAPSULE — project reference

**C**ompact, **A**ddressable, **P**seudogenome-**S**tructured, **U**nified **L**ossless **E**ncoder.

Read this file first. It records what is done, what is verified, what is NOT
done, and which ideas have already been tested and refuted so they are not
attempted again.

**For the 15-dataset SPRING/Genozip comparison (this repo, `c_star_pg_advance`,
not the outer ARCS binary), the locked set is `NEW_DATASET_LOCKED.md`, not the
17-accession list in `DATASET_LOCKED.md` below.** It records one swap made
2026-09-02 (Drosophila out, Utricularia gibba/Plantae in) and why. Read it
before running or citing that comparison.

**`docs/PHASE2B_RESULT.md` is VOID for 4 of 14 datasets — do not quote it.**
Its numbers (14/14 vs SPRING -22.00%, 14/14 vs Genozip -72.83%) were measured
before 4 silent data-loss bugs were found and fixed same-day (`23be207`,
`121fea9`, `f3ab0c2`): the archives for ERR552797, SRR40271341, SRR40402583 and
SRR32429602 did not actually decode to their input. All four are fixed and
verified LOSSLESS now (see section 6.3). A full re-run with Phase 1 / Phase 2b
/ Phase 3 measured separately, each level round-trip verified before its number
is recorded, is in progress — check for `docs/PHASE3_RESULT.md` or the latest
ALLPHASES result before citing any size figure from this repo.

**Quality is now wired** (`include/quality_coder.h`, vendored fqzcomp/htscodecs,
BSD 3-clause, gated on `CAPS_QUAL=1`) — see commit `908b769` and section 6.3.
The archive can now reproduce a complete 4-line FASTQ from the archive alone;
verified byte-identical (same MD5) on at least one dataset, full sweep pending.

Repo: `github.com/thackshanaramana0-spec/c_star_pg_advance`, branch
`c_star_pg_advance`, 127 commits, working tree clean.

---

## 1. What this is

An independent from-scratch implementation of pseudogenome-based read
compression, benchmarked head-to-head against PgRC2 (`/root/arcs-clean/method_c`,
cloned separately, GPL-3, never vendored).

**Scope as of 2026-09-02: sequence + read order + names + line 3.** Names are
now really in the archive (`include/names_coder.h`, streams `names_body` /
`names_dict` / `names_index`, behind `CAPS_NAMES=1`; Phase 1 stays
byte-identical with it unset) and line 3 is stored as a one-byte file mode.
**Quality is the one column still out** -- see the quality section of
`docs/REIMPL_NOTES.md`: assessed and measured, the decision is to vendor
fqzcomp (BSD 3-clause) rather than reimplement. Stage 92 is our own quality
coder and beats SPRING 7/8 and Genozip 8/8, but loses to real fqzcomp 0/8.

Method, in order: greedy overlap chaining builds a pseudogenome from
well-tiling reads -> remaining reads are pigeonhole-mapped onto it -> unmapped
reads are appended and assembled as a second region -> the pseudogenome is
self-matched to remove redundancy -> everything is emitted as separate streams
and entropy-coded.

---

## 2. Current standing — corrected 2026-09-02, read the warning first

**Every size figure published before commit `a81f55c` was measured against an
INCOMPLETE archive** that could not be decoded. See section 6.1. These are the
corrected numbers, taken from the real file size on disk:

| dataset | kingdom | ours | PgRC2 | margin | previously claimed |
|---|---|---|---|---|---|
| H. salinarum | Archaea | 2,788,399 | 3,050,477 | +8.59% | +14.57% |
| E. coli | Bacteria | 8,199,540 | 8,864,420 | +7.50% | +11.38% |
| L. major | Protista | 27,860,951 | 28,272,652 | +1.46% | +2.83% |
| P. aeruginosa | Bacteria | 8,961,459 | 9,043,181 | +0.90% | +1.84% |
| S. aureus | Bacteria | 13,506,629 | 13,595,003 | +0.65% | +3.77% |
| P. falciparum | Protista | 17,118,655 | 17,219,695 | +0.59% | +5.38% |
| **S. acidocaldarius** | Archaea | 3,143,897 | 3,114,782 | **-0.93%** | +0.33% |
| **aggregate** | | **81,579,530** | **83,160,210** | **+1.90%** | +4.65% |

**6 wins, 1 loss.** S. acidocaldarius is a loss again -- the duplicate-chaining
fix (`3e06957`) was real and shrank that pseudogenome 32%, but the flip to a win
was measured against an incomplete archive and did not happen.

The three datasets PgRC2 cannot process have NOT been re-measured since the
correction and their figures are stale:

    SARS-CoV-2, A. fumigatus, S. cerevisiae -- re-run before quoting

Speed/RAM against PgRC2 (E. coli, idle machine): ~1.7x slower, ~2.5x heavier at
worst. Wall 243.78 s -> 130.60 s (-46.4%) across 7 files, worst-case peak RSS
1,388 MB -> 930 MB (-33.0%). These predate the correction; the added streams
cost encode time that has not been re-measured.

## 3. Repository layout, and why

    stages/       96 .cpp -- the full 01..106 progression, one file per experiment.
                  Kept whole: each stage is the evidence for a decision, and
                  several were later contradicted by measurement. 106_inprocess.cpp
                  is the shipped encoder.
    include/      coders_inproc.h (stream coders, selector, transforms),
                  coders_pgrc.h (PPMd7 / FSE / range coder + their decoders),
                  seqpar_core.h (the DNA coder, shared so the standalone binary
                  and the in-process path cannot diverge)
    scripts/      build106.sh, encode_adaptive.sh, decode_105.py,
                  verify_lossless.sh
    thirdparty/   ppmd/ (LZMA SDK, public domain), fse/ (Yann Collet, BSD)
    docs/         analysis and plans; retractions are marked in place, never deleted
    results/      raw measurement output, including runs that were reverted
    DATASET_LOCKED.md, CLAUDE.md   copied in because the work depends on both

Data files, compiled binaries and regenerable stream dumps are gitignored. The
`.gitignore` rule was originally `mem_*.bin` while the dumps are `mm_*.bin` --
one letter apart -- which left 18 MB unignored until it was fixed.

---

## 4. Commands

```bash
scripts/build106.sh /tmp/best106                 # build (must include -fopenmp)
INPUT=reads.fq ARCHIVE=out.capsule BEST=/tmp/best106 \
    bash scripts/encode_adaptive.sh              # encode
scripts/verify_lossless.sh reads.fq              # encode -> decode -> diff vs original
```

`GSEARCH=1` selects golden-section search over MAXMAP instead of the 4-point
grid: smaller archives (-154,223 B over 7 files) at ~7-10 probes instead of 4.

**build106.sh exists because the build command previously lived only in shell
history, which is exactly how `-fopenmp` went missing.** Without it the single
`#pragma omp parallel for` is silently discarded and the largest stage runs
serial; linking it correctly was worth -45% wall time.

---

## 5. What is verified, and how

`scripts/verify_lossless.sh` encodes, decodes, and compares against the
**original FASTQ's sequence column** -- not a coder-level round trip. Verified
LOSSLESS on S. acidocaldarius, E. coli, SARS-CoV-2, H. salinarum,
P. aeruginosa, L. major, S. aureus.

**Read this carefully: that check routes through the RAW intermediate streams
the encoder dumps, not through the archive.** It exercises the whole algorithm
-- assembly, mapping, mismatches, read order -- but not the entropy layer.

---

## 6. What is NOT done — checked, not assumed

### 6.1 The verification path does not read the product

`scripts/verify_lossless.sh` decodes the DUMPED streams, not the archive. The
encoder writes both: an archive (2.6 MB on H. salinarum) and 13 uncompressed
dump files (11.8 MB). Only the archive ships and is measured; only the dumps are
verified.

That gap hid a real defect for the whole session. `refc::encode` stored only a
reference's source; its destination gap, length and RC flag were computed on the
next lines and dropped as "diagnostics only". Nothing else carried them -- our
literal is pure ACGT with no marker, unlike PgRC2's in-band MATCH_MARK -- so the
archive could not place a single reference. The dump's mem_triples.bin is
2,299,531 B with all four fields; the archive stream was 489,798 B with one.
Fixed in `a81f55c` by adding mem_dstgap, mem_len and mem_rc: 2,286,493 B across
seven files, and the aggregate margin fell from +4.65% to +1.90%.

ARCHIVE_TOTAL was never wrong -- it matches the file on disk byte for byte. The
file was incomplete, not the count.

**Until verify_lossless.sh decodes the ARCHIVE, no size figure in this repo is
final.** This project had already recorded the identical failure once (totals
omitting mm_pos and mm_count); the lesson did not generalise because the check
that enforces it does not exist.

### 6.2 The archive decoder is incomplete — SUPERSEDED, see 6.3

**This section is historical.** As of 2026-09-02 (commits `908b769` onward)
`capsule_decode` reconstructs sequence, names, quality AND line 3 from the
archive alone, and a complete 4-line FASTQ rebuilt from decoder output only has
been verified byte-identical (same MD5) to the original file. Read 6.3 for the
current state; this section is kept for the historical record per this file's
own rule that retractions are marked in place, not deleted.

`stages/capsule_decode.cpp` reads the
CAPSULE container and inverts the general-purpose coders. Verified identical
against the encoder's own in-memory streams on H. salinarum: `pos_abs`
(1,842,004 B), `read_lengths` (921,002 B), `pos_strand`, and the three N
streams -- 6 identical, 0 differ.

**Four streams have NO inverse and are skipped.** Confirmed by running the
decoder, not by reading the code:

    literal       seq_encode_mem          NO inverse
    mem_triples   refc::encode            NO inverse
    mm_sym        mmc::encode             NO inverse
    mm_pos        mmpos_encode_buckets    NO inverse (bucketed form only)

So `capsule d` cannot reconstruct reads from an archive alone today. PPMd, FSE
and the range coder DO have inverses written, but H. salinarum selected LZMA
variants for every stream, so those three paths are untested on real data.

**Other gaps:**

- **T. cacao (22 GB) never completed** -- cancelled at ~28 minutes. Plantae is
  the only uncovered kingdom. C. elegans (11 GB) did complete: 179,949,366 B,
  31:01, peak RSS 6,117 MB.
- **SARS-CoV-2 regressed +3,233 B (+0.38%)** from the L=Lmax change. Disclosed,
  not tuned away. C. jejuni, the other variable-length file, improved -12,917.
- **M. tuberculosis, H. pylori, C. jejuni, D. melanogaster** on disk, not run.
  HCMV not downloaded.
- **The 7 core datasets have no SPRING/Genozip numbers** -- only the three above.
- **Genozip's fungi results are anomalous** (162 MB where SPRING gets 24 MB) and
  unexplained. Do not put them in a table until diagnosed.
- **Human is out of scope by an earlier decision** (`ERR174310` is on the banned
  list, "too large, human excluded"). All three published human benchmarks use
  files we either banned or do not have. chr20 at 30x is feasible; full WGS is
  not at ~700M reads.

### 6.3 Four silent data-loss bugs found and fixed, 2026-09-02

Found by actually decoding archives (something Phase 1/2b had never done — see
6.1) and comparing to the original FASTQ, dataset by dataset, not by trusting
size tables. All four are fixed, committed, and each is keyed on a measured
property of the input, not a per-dataset special case (standing rule 1, below).

1. **`23be207`** — mismatch positions were stored in ONE byte; above Lmax=256
   the encoder clamped every position past 255 to 255, silently corrupting any
   read longer than 256 bp. Now a varint. Two locked datasets exceed 256 bp
   (ERR552797 301bp, SRR40271341 300bp) and both were being reported as
   compression wins from archives that could not reproduce their input.
2. **`121fea9`, bug 1** — mismatches were still emitted for "orphaned" unique
   reads that no original read maps to (a containment artifact, variable-length
   only). The decoder can't derive such a read's length, computed a wrong `ref`
   byte, and — because the mismatch coder is adaptive — that desynchronised
   every mismatch symbol after it. 8 bad bytes of 18,774 became 3,712 wrong
   reads on a C. jejuni sample.
3. **`121fea9`, bug 2** — a contained reverse-strand read (a prefix in read
   space is a SUFFIX in pg space) was indexed with the ORIGINAL's length
   instead of the UNIQUE's, so it was rebuilt from the wrong end. 95 of 100,000
   on the same sample, all RC, all shorter than their unique — unambiguous
   signature. The same bug exists in `scripts/decode_105.py`.
4. **`f3ab0c2`** — the most general one. `FSE_compress` returns a 1-byte RLE
   result for a constant input, and that byte is NOT the repeated symbol. The
   decoder's existing "1 byte -> fill with that byte" handling only happened to
   be correct for an all-ZEROS stream (orig2uid_flags, the case it was written
   for); any OTHER constant stream silently decoded as zeros. M. tuberculosis
   has three all-N reads, `n_cnt` was three bytes of `0x23`, FSE wrote `0x00`,
   and every N in those reads was lost. Not dataset-specific — reachable by any
   constant non-zero stream on any input. Fixed by rejecting FSE/HUF's RLE
   result at encode time so the selector falls through to a coder that
   round-trips.

**Why none of this was flagged earlier:** Phase 1 and Phase 2b measured archive
SIZE only; nothing ever decoded an archive. The one lossless check that existed
(`verify_lossless.sh`) tests the DUMPED intermediate streams, not the archive —
see 6.1 — so the entropy layer, where bug 4 lived, was never exercised at all.
It had also only ever run on 7 of the 17 locked datasets, none of which were
the four that turned out to be broken.

**Cost of correctness:** +33 B on E. coli (+0.0009%), the only fixed-length
dataset affected (it had a zeros-constant stream that took the FSE-RLE path;
it decoded correctly before the fix by coincidence, so only its bytes changed,
not its correctness). Verified byte-identical otherwise. Full audit: 11/11
datasets checked so far are LOSSLESS post-fix; a complete 14-dataset sweep with
Phase 1/2b/3 measured separately is in progress.

---

## 7. Ideas already tested and REFUTED — do not redo

Each was implemented and measured, not argued away.

| idea | result |
|---|---|
| A2: work-derived pool scheduling | +3.2% slower. Source bytes do not predict coder time: `literal` has the most bytes and takes 0.43 s, `orig2uid` fewer and 4.76 s. |
| Cost model `m* = L(1-r)*bpb/bpm` | 55.5% mean error, under-predicts every dataset |
| ...plus a mem_triples reference term | 45.3%, errors now mixed in sign |
| Per-read pricing (LZMA GetOptimum style) | 12-18% WORSE than the swept optimum. The decisions are COUPLED: mapping a read removes it from the append pool and destroys redundancy for every other leftover. Per-item pricing cannot express that. |
| Per-region MINMEM | monotonically worse (+30,652 at 32, +2,010,039 at 128). The 431K second-region references each pay for themselves. |
| Splitting positions by region | we already code BELOW the region-split bound (3,629,158 vs 3,699,725) |
| MINOV as a pg-span lever | already optimal on both files tested, and the curve is NOT unimodal, so no search applies |
| Minimum-degree-first matching | byte-identical output. 96.8% of tails have ZERO candidates at a level and under 1% have more than one -- there is no contention for order to resolve. |
| Fixing the tail-eligibility bug | Real bug (reads shorter than the starting L are permanently dropped), but fixing it measured WORSE: C. jejuni +26,109, SARS +4,178. Short reads should be mapped, not assembled. |
| Full-length links as containment | the guard never fires; identical archives with and without |
| Parallel coder probes | 6-11% faster for +54-86% RAM |
| Architecture rewrite (earlier session) | 0% size, 0% speed, +54 MB RAM |

**Withdrawn conclusions:** libgomp barrier spin is NOT our bottleneck -- PgRC2
spins MORE (24.21% vs 22.50%) and still finishes in a third the time. The real
differentiator was page faults (1,089,328 vs 169,421), caused by a 64 MB LZMA dictionary
requested regardless of stream size.

---

## 8. What is open, with evidence

**Second-region self-match** -- the single largest quantified opportunity.
`run(Q, qlen, CROSS, ...)` matches the second region against the MAIN pg only,
never against itself. On SARS-CoV-2 that removes 2.2% of a 14.7 MB region while
**74.3% of its 32-mers are repeat occurrences**. Implemented twice; both work
structurally and give SARS -94,279 to -108,432, E. coli +61,792 to +95,373,
discriminated by how much the cross-pass already removed (97.8% / 39.2% / 12.7%
unremoved). **Not shipped: the archive is LOSSY.** Reconstruction does not
honour references whose source is in the second region; the RC path with a
non-zero SRCBASE is the first place to look. See `docs/SECOND_REGION_SELF_MATCH.md`.

**Where the remaining size sits** (S. acidocaldarius, at its optimum):
`pos_abs` 45.8% of the archive, `mm_pos` 22.7%, `literal` 19.2%. Positions
dominate every dataset (45.8 / 46.2 / 48.5%) and are already coded at 0.92-0.99x
their set-plus-permutation bound, so there is little left there.

**RAM/speed headroom, measured:** A3's fork holds parent and child resident
simultaneously (849 MB on C. elegans, 1,695 MB projected on T. cacao) -- an
in-process loop needs only `admit` saved, 29 MB and 58 MB. A4's per-thread slot
arrays cost 4 B x reads x 12 and are an active regression on low-coverage data
(-41 MB on E. coli, +1,052 MB on S. aureus); they should engage on measured hit
density rather than always.

---

## 9. Standing rules

1. **Fixes must be formulas over a measured input property, never per-dataset
   special cases and never a fitted constant.** Three cost models were attempted
   this session and all three were discarded rather than patched with a
   multiplier.
2. **Every change is gated.** Output-preserving changes must be `cmp`-identical
   on all 7 files; output-changing changes must not regress any file and must
   improve the aggregate. Lossless is re-verified either way.
3. **A change that fails its gate is reverted and recorded, not tuned until it
   passes.** `results/phase_a/04_gate_A2_REVERTED.txt` exists for this reason.
4. **Never run two timed jobs concurrently.**
5. **Retractions stay in the docs.** Several conclusions here were wrong and are
   marked in place so they are not re-derived.

---

## 10. Commit map

    c91ba02  CAPSULE format: container + archive decoder (partial)
    d1812f5  Second-region self-match: quantified, reverted as LOSSY
    987843b  Variable-length: two hypotheses refuted and reverted
    3e06957  Flip the last loss: sweep from L=Lmax so duplicates chain
    94d27b2  Refuted: min-degree-first matching
    a2f7e5b  Plan: flip the loss via chain quality
    1ccaf47  Measured answer: we do not need PgRC2's 3-way split
    8a2864f  Golden-section search over MAXMAP, behind GSEARCH
    f8db7be  Third failed derivation: per-read pricing
    4623ae6  Dedup table + constant-stream elimination
    0234f4c  Plan: derive the mapping ceiling
    22092fa  Headroom analysis, file-and-line
    b6e91eb  Plan: what is open and in what order
    71ef868  B3: split orig2uid into flags + values
    0efbc9b  B2: bound the LZMA dictionary by input length
    aa80964  B1: store incompressible byte planes verbatim
    d73fea8  Fix: adaptive children overwrote each other's dumps
    b74a961  Restructure as a standalone repository

The single most important change of the session is `3e06957`. A suffix-prefix
overlap of exactly the read length IS an exact duplicate, and it is the only
length at which one can appear -- two identical 251-base reads do not overlap at
L=250. The sweep started at Lmax-1, so duplicates were invisible to chaining and
had to be handled by a pre-assembly dedup pass costing an `orig2uid` array over
every original read; a threshold then chose between those two bad options. On
S. acidocaldarius (9.4% duplication, below the 15% threshold) all 48,694
duplicates went into the pseudogenome, while PgRC2's log reports the identical
count removed for free. Starting the sweep at `Lmax` cut that pseudogenome
9,134,100 -> 6,157,270 and turned the project's last loss into a win.
