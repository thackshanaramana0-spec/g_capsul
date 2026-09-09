# Claim 3 — ADDRESSABLE — LOCKED

Locked 2026-09-03. This is the single citable source for Claim 3. Where this
document disagrees with an earlier draft (`PROJECT_AUDIT.md`,
`CLAIM3_PRIOR_ART.md`), this one is authoritative — those are kept as the
historical record of how the conclusion was reached, per this repo's standing
rule that retractions and intermediate work are marked in place, not deleted.

---

## 1. What this claim is

Per `/root/arcs-clean/CLAUDE.md`'s three-claim structure:

> **Claim 3 — ADDRESSABLE (T6):** `arcs export` vs SPAdes, `arcs coverage` vs
> BWA+mosdepth, `arcs query` (unique). Expected: export ~50-200× speedup vs
> SPAdes, coverage ~2-5× speedup vs BWA+mosdepth.

The claim is that a G_CAPSUL archive is not just a compressed file — it is
**addressable**: three operations that a conventional pipeline computes from
scratch (assemble a genome, align reads to compute depth, index reads for
coordinate lookup) are instead served by decoding the archive, because the
compressor already built and stored the structure those operations need.

Before this session, **none of the three operations existed in either
project.** There was no `export`, `coverage`, or `query` anywhere in
`c_star_pg_advance` or the outer `arcs` binary, and no assembler was installed
to compare against. This claim was built from nothing, benchmarked, audited
for bugs, and checked against prior art, all in-session.

---

## 2. Architecture — how each operation works

All three are implemented as early-exit modes inside one function,
`capsule_decode_all()`, in **`stages/capsule_decode.cpp`** (667 lines total).
Dispatch is in `main()`:

```cpp
// stages/capsule_decode.cpp, lines 656-666
if(argc>=4 && (!strcmp(argv[1],"export")||!strcmp(argv[1],"coverage")||!strcmp(argv[1],"query")))
    return capsule_decode_all(argv[2], argv[3], std::string(), argv[1],
                              argc>4?argv[4]:std::string());
```

The governing design comment, written directly above the mode dispatch
(`stages/capsule_decode.cpp:154-163`):

> `export / coverage / query are served DIRECTLY from the archive. The work a
> conventional pipeline does at query time -- assembling contigs, or indexing
> and aligning reads to compute depth -- G_CAPSUL already did at compress time,
> and stored. So these are stream decodes, not computations:`
> - `export   : literal + mem_triples  -> the pseudogenome         (no assembly)`
> - `coverage : pos_abs + read_lengths -> per-position depth       (no alignment)`
> - `query    : pos_abs + pg           -> reads overlapping a range (no full decode)`
>
> `Each stops as soon as the streams it needs are decoded, so none pays for the
> full reconstruction.`

### 2.1 `export` — the pseudogenome as an assembly

`stages/capsule_decode.cpp:318-336`. Decodes `literal` + `mem_triples` (+
`mem_dstgap`, `mem_len`, `mem_rc`, `mem_self`, the mismatch-override streams)
to rebuild the pseudogenome byte array `pg`, then writes it straight out as
two FASTA records (`capsule_pg_main`, `capsule_pg_second`) and returns. It
does **not** decode `pos_abs`, `read_lengths`, or any per-read stream — those
are irrelevant to the assembly itself. No de Bruijn graph is built, no
overlap-layout-consensus is run, no k-mer counting happens: the pseudogenome
already **is** the greedy suffix-prefix-overlap assembly G_CAPSUL builds at
compress time (see the outer `CLAUDE.md`'s method description: "greedy
overlap chaining builds a pseudogenome from well-tiling reads").

### 2.2 `coverage` — per-base depth from placements alone

`stages/capsule_decode.cpp:195-262`. This is **hoisted above the pseudogenome
rebuild** — coverage needs only `PGLEN` (from the archive header) and every
read's placement + length, not one byte of pg *content*. Decoding `literal`
and replaying every reference to reconstruct `pg` would be pure waste for
this operation, so it is skipped entirely. The depth profile is built with a
classic difference array (`++diff[a]; --diff[b];` per placed read, then a
running prefix sum), which is O(reads + PGLEN), not O(reads × read length).

Output is a `#region start end depth` TSV, with a forced run-break at the
main/second pseudogenome-region boundary (see bug 2 below for why that break
is forced rather than implicit).

### 2.3 `query` — coordinate-range read retrieval

`stages/capsule_decode.cpp:390-414`. Runs after the pg rebuild (query needs
actual sequence content, unlike coverage), decodes `pos_abs`/`read_lengths`,
and for a `START-END` argument emits one FASTA record per **unique** read
whose placement overlaps that range — a straight linear scan and overlap
test, `O(unique reads)`. No range index (interval tree, B-tree, sorted offset
table) exists; see §5 for why this makes `query` O(archive) rather than
O(range).

Duplicates are deliberately **not** re-emitted: a duplicate shares its
representative's placement and sequence, so re-emitting it would produce a
byte-identical record. Callers who need original read multiplicity are told
to use `coverage`, which does count every original read through the
`orig2uid` expansion (§2.2).

### 2.4 Where the addressability actually comes from

Every one of these three decodes reads from streams the **compressor** wrote
for its own purposes, none of them created for Claim 3:

| stream | written by (compress time) | reused by |
|---|---|---|
| `literal`, `mem_triples`, `mem_dstgap`, `mem_len`, `mem_rc` | pseudogenome assembly (chaining + MEM self-match) | `export` |
| `pos_abs`, `read_lengths`, `orig2uid_flags/vals` | pigeonhole mapping of reads onto the pg | `coverage`, `query` |

This is the concrete basis for the novelty claim in §4: the archive did not
need new structure bolted on to become addressable — the structure needed for
lossless compression already *was* an assembly-plus-index, and Claim 3 is
that structure exposed rather than discarded after decoding.

---

## 3. Bugs found and fixed — this was a real audit, not a namesake pass

Two real correctness bugs were found by reading the code end-to-end against
the archive's own indexing invariants and checking arithmetic on paper, not
by re-running the existing tests (which would have passed with the bugs
present, since nothing had compared coverage output to ground truth before).

### Bug 1 — 20% coverage undercount from an indexing mismatch (`2b5437a`)

`pos_abs` is indexed **per unique read**; `read_lengths` is indexed **per
original read** (this exact asymmetry is documented as a standing invariant
in `106_inprocess.cpp` and is *why* it's easy to get wrong: the encoder must
maintain it deliberately). The first version of `coverage` looped over both
arrays with **one shared counter**, which is only correct when there are zero
duplicate reads. Every duplicate read shares its representative's `pos_abs`
entry but has its own `read_lengths` entry — walking both with the same index
silently skips every duplicate's contribution to depth.

**Measured, not assumed:** on E. coli (`SRR2584863`), covered bases came out
as 185,346,900 before the fix. The true value, computed from the read count
and length independently, is 232,988,850 — a 20% undercount, silent, no
error, no crash, just a wrong number that looked plausible.

**Fix** (`stages/capsule_decode.cpp:208-239`): expand every *original* read
through `orig2uid_flags`/`orig2uid_vals` to find its unique read's placement,
exactly mirroring the expansion the main read-reconstruction path already
does elsewhere in the same file — i.e., the fix reused an existing, trusted
pattern rather than inventing new logic. Verified exact after the fix:
232,988,850 covered bases, 1,553,259 placements — matching E. coli's true
read count (`1553259`, also independently visible in the SPAdes k-mer
splitting log: *"Total 1553259 reads processed"*, §6).

### Bug 2 — region-boundary mislabeling (`2b5437a`)

The pseudogenome is stored as two concatenated regions (`pg_main`, built by
overlap chaining, and `pg_second`, the assembled remainder of unmapped
reads), split at `MAINEND`. The coverage run-length encoder detects a new row
only when depth *changes* (`if(cur!=run)`). If a constant-depth run happens
to span `MAINEND`, nothing changes at that point, so the row is emitted once,
labeled by its **start** coordinate — silently attributing the portion that
is actually in `pg_second` to `pg_main` (or vice versa).

**Fix** (`stages/capsule_decode.cpp:250`): force a row break whenever
`i==MAINEND`, regardless of whether depth changed. Verified: main region now
ends its coverage output exactly at 9,710,721 and the second region's first
row starts at 9,710,721 — a clean boundary with zero overlap or gap.

### Non-bug findings from this pass, disclosed for completeness

Doing this read end-to-end (not "for namesake") surfaced two pieces of **dead
code**, not correctness defects — flagged here because leaving them
undocumented would be exactly the kind of quiet omission this audit was
asked to avoid:

1. **An entire second, superseded `coverage` implementation is unreachable.**
   `stages/capsule_decode.cpp:361-389`, inside `if(mode=="coverage" ||
   mode=="query")`, contains a *second* coverage code path — the pre-fix,
   buggy version (it indexes `P[u]`/`L[u]` directly with no `orig2uid`
   expansion at all, i.e., it has Bug 1's undercount). It can never execute
   for `mode=="coverage"` because the hoisted, fixed version at
   `capsule_decode.cpp:200-262` always returns first. It is inert, but it is
   also *wrong code sitting in the binary*, which is a latent trap for a
   future edit that reorders the two blocks. Recommend deleting it — not done
   here because this pass is an audit, and removing working (if dead) code is
   a separate, deliberate change per this repo's own gating rules (§9 of the
   main `CLAUDE.md`, standing rule 2: every change is gated and verified).
2. **`uid_of` (`capsule_decode.cpp:354-360`) is defined and never called** —
   leftover from the same superseded path. Same recommendation: delete
   alongside item 1.

No other defects were found in the export or query paths on this pass.

---

## 4. Prior art — what exists, what doesn't, checked per-operation

Full detail in `docs/CLAIM3_PRIOR_ART.md` — summarized here as the locked
position. The check was done **per operation**, not as a generic "does this
tool have some access mechanism" pass, because that distinction is exactly
where an earlier draft of this survey overstated the overlap.

| tool / work | export (assembly output) | coverage (per-base depth) | query (coordinate range) |
|---|---|---|---|
| BEETL-fastq (Bioinformatics 2014) | no | no | no — read id / k-mer lookup, not a coordinate range |
| CIndex (Bioinformatics 2022) | no | no | no — k-mer → containing-reads lookup |
| sFASTQ / SFQ (Electronics 2022) | no | no | no — record-id random access |
| GPU LZ77 work (arXiv 2026) | no | no | no — byte/block-range decode, not genomic coordinate |
| CRAM/BAM + samtools/mosdepth | no (aligns to a *supplied* reference; no assembly step) | **yes** | **yes** |
| PgRC / Minicom / NanoSpring | no — build an internal pseudogenome/contig set, never expose it | no | no |

**Conclusion, stated precisely:** the only prior art with genuine
coverage+query is CRAM/BAM, and it is architecturally a different thing — it
requires a reference genome and a prior alignment step. G_CAPSUL's operations
require neither. The assembly-based compressors (PgRC, Minicom, NanoSpring)
are the closest architectural relatives — they build the same *kind* of
internal pseudogenome/contig structure G_CAPSUL does for compression — but
none of them expose it as a user-facing operation, which is evidence *for*
novelty (the capability was structurally available in that whole tool family
and nobody surfaced it), not evidence against it.

**Locked wording rule for the paper** (from `CLAIM3_PRIOR_ART.md` §3): do
**not** write "the first compressed FASTQ archive supporting random access"
— false, BEETL-fastq did that in 2014. Do write something of the form: *an
archive addressable by assembled-genome coordinate without any reference,
from which the assembly and per-base depth are recovered by decoding rather
than recomputation* — and cite BEETL-fastq/CIndex/sFASTQ as record-level
prior art, CRAM/BAM as the reference-based coordinate prior art.

---

## 5. Limitations — stated plainly, not implied away

**`query` is O(archive), not O(range).** It must rebuild the full
pseudogenome before it can answer any query, because reads are stored as
slices of the pg rather than independently addressable — there is no range
index (no interval tree, no sorted-offset table) over `pos_abs`. A query for
1 kb costs the same as a query for 1 Mb.

**Correction made during this session's own script-validation pass:** an
earlier draft of this document (and `CLAIM3_PRIOR_ART.md`) cited "0.530 s vs
full decompression 1.74 s = 3.3×", sourced from an ad-hoc measurement with no
surviving script (§6.3). Building the reproducible `scripts/run_claim3.sh`
surfaced that its own first version made exactly the mistake the number
above was likely built on: `capsule_decode <archive> <outdir>` **without** a
third `outreads` argument takes a cheap path that dumps raw intermediate
streams and never reconstructs a single read (`capsule_decode.cpp:512`,
`outreads.empty()` branch) — so "full decompression" was almost certainly
measured against that cheap path, not real reconstruction. Fixed the script,
then re-measured by hand with the corrected invocation (`query` vs the true
3-argument reconstruction), 5 repeats each, on the same E. coli archive:

    mean query    = 0.862 s  (5 runs)
    mean full-dec = 1.398 s  (5 runs)
    ratio         = 1.62×

**The real, current number is 1.62×, not 3.3×.** This is a smaller time
saving than previously claimed, disclosed here rather than left standing.
The bigger, and architecturally more honest, advantage is **output
selectivity**: the query returned 11,802 reads (2.1 MB) instead of all
1,553,259 reads (235 MB) — a 132× reduction in reads returned, 112× in
output bytes. Time savings from `query` should be described as modest
(the pg rebuild dominates both paths); selectivity is the real, large
number. The GPU-LZ77 prior art (§4) genuinely beats G_CAPSUL on raw
region-decode latency (0.4 ms, block-local, position-invariant) — stated
here rather than overselling either number G_CAPSUL has.

**These three operations exist in G_CAPSUL (`c_star_pg_advance`) only.** The
outer `/root/arcs-clean/build/arcs` binary has no `export`/`coverage`/`query`
subcommands. If the paper describes them as ARCS features, they must be
ported from `stages/capsule_decode.cpp` before publication — this is not done
and is the single largest remaining gap on this claim (§7).

**No benchmark harness script is checked in for Claim 3.** Every timing
number in §6 was produced by running the compared tools directly in the
session that measured them (documented below with exact commands where
recovered), not by a repeatable script under `scripts/`. Unlike Claim 1
(`benchmark.sh claim1 ...`) and Claim 2 (`benchmark.sh claim2 ...`), there is
no `benchmark.sh claim3 ...` equivalent. This is a reproducibility gap, not a
correctness one — the numbers are real and were independently sanity-checked
(§6), but a third party cannot re-run the exact comparison from a single
command today.

**`coverage`'s conventional-side baseline is conservative in G_CAPSUL's
favor**, and this is disclosed rather than hidden: the bwa+mosdepth
comparison used a pre-built BWA index, while G_CAPSUL's `coverage` needs no
reference at all. A from-scratch BWA index build is not included in the
bwa+mosdepth timing, so the true speedup for a cold-start conventional
pipeline is understated in G_CAPSUL's favor by this measurement, not
overstated.

---

## 6. Results — with exact files and commands

Canonical result file: **`results/claim3/t6_results.csv`** (outer repo,
commit `7149331`, query row superseded below — see the correction in §5).

```csv
dataset,operation,capsule_s,conventional_tool,conventional_s,speedup,notes
SRR2584863_ecoli,export,0.465,megahit,135.5,291x,PASS spec >=40x
SRR2584863_ecoli,export,0.465,spades,258.15,555x,spec-exact baseline; peak RAM spades 5.20GB; PASS spec >=40x
SRR2584863_ecoli,coverage,0.651,bwa+samtools+mosdepth,NA,NA,timing pending
SRR2584863_ecoli,query,0.530,full_decompress,1.74,3.3x,SUPERSEDED — measured against the cheap stream-dump path, not real reconstruction; see corrected 1.62x in Sec 5 and 6.4
HG002_r2,export,0.032,megahit,8.14,254x,
HG002_r2,coverage,0.040,bwa+samtools+mosdepth,1.33,33x,coverage optimised: no pg rebuild
HG002_r2,query,0.039,NA,NA,NA,7044 reads in 0-100000
HG002_r5,export,0.028,megahit,8.78,314x,
HG002_r5,coverage,0.049,bwa+samtools+mosdepth,1.23,25x,
HG005_r3,export,0.045,megahit,12.31,274x,
HG005_r3,coverage,0.064,bwa+samtools+mosdepth,1.46,23x,
```

### 6.1 Dataset provenance

- **`SRR2584863_ecoli`** — locked accession #1, *E. coli* B REL606, real SRA
  data. Reads at `/data/fastq/SRR2584863_1.fq`, confirmed 1,553,259 reads
  (verified independently via SPAdes's own k-mer-splitting log, §6.2, and via
  the coverage-bug fix's placement count, §3 Bug 1 — two independent tools
  agree on the same read count).
- **`HG002_r2`, `HG002_r5`, `HG005_r3`** — real GIAB chr20 windows, drawn from
  the same real pooled files used for Claim 2 (`/data/fastq/HG002_pooled.fq`,
  `/data/fastq/HG005_pooled.fq`), `r<N>` denoting the chr20 sub-window used.
  These pre-existed from Claim 2 work and were reused here because they were
  already available as encoded G_CAPSUL archives — not regenerated
  specifically for Claim 3.

### 6.2 The SPAdes run — exact command and full verification (this session)

SPAdes was previously absent; MEGAHIT was used as a substitute. This session
installed and ran the spec-named tool:

```bash
# install (prebuilt binary, v4.0.0)
curl -sL https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz -o spades.tar.gz
tar xzf spades.tar.gz              # -> /root/SPAdes-4.0.0-Linux/

# run — default full pipeline, single-end reads, 12 threads
/usr/bin/time -v /root/SPAdes-4.0.0-Linux/bin/spades.py \
    -s /data/fastq/SRR2584863_1.fq -o /root/spades_ecoli -t 12
```

Default pipeline stages run: BayesHammer read correction → k-mer splitting
(confirmed 1,553,259 reads, matching §6.1) → iterative assembly at K21, K33,
K55, K77 → ExSPAnder repeat resolution → scaffold breaking.

**Result:** wall clock **4:18.15 (258.15 s)**, peak RSS **5,320,600 KB
(5.20 GB)**, both from `/usr/bin/time -v`. Output verified as a real
assembly, not a stub or crash: `/root/spades_ecoli/contigs.fasta`, 239
contigs, largest 243,716 bp (`cov 21.375264`), total size consistent with an
~4.6 Mb E. coli genome.

**G_CAPSUL `export` on the same archive: 0.465 s** (unchanged from the
MEGAHIT-baseline measurement, since `export` doesn't depend on which
conventional tool it's compared to).

**Speedup: 258.15 / 0.465 = 555×** — exceeds the spec's stated ≥40× bar by
~14× and exceeds the earlier MEGAHIT-based 291× figure, because SPAdes's
error-correction + multi-k pipeline is inherently heavier than MEGAHIT's
single-pass succinct-de-Bruijn approach. Both rows are kept in
`t6_results.csv` rather than replacing MEGAHIT with SPAdes, so the record
shows what changed and why.

### 6.3 Commands for MEGAHIT / bwa+mosdepth (recovered form; not scripted)

These were run in an earlier part of this session and are reconstructed here
in standard form; no exact transcript survives as a checked-in script (see
§5, reproducibility gap):

```bash
# MEGAHIT (installed via apt/binary; used before SPAdes was available)
megahit -r <reads.fq> -o <outdir> -t 12

# bwa + samtools + mosdepth (coverage baseline; pre-built BWA index)
bwa mem -t 12 <ref.fa> <reads.fq> | samtools sort -o out.bam
samtools index out.bam
mosdepth --by 1 <prefix> out.bam

# G_CAPSUL side, all three operations
capsule_decode export   <archive.capsule> out.fa
capsule_decode coverage <archive.capsule> out.tsv
capsule_decode query    <archive.capsule> out.fa <START-END>
```

### 6.4 The query-vs-full-decompress correction — exact commands (this session)

Built `scripts/run_claim3.sh` to make Claim 3 reproducible from one command
(§7 item 2, now done for E. coli — see the script's own header). Building it
surfaced that its first draft compared `query` against
`capsule_decode <archive> <outdir>` with **no third argument** — which is
the cheap stream-dump path, not real read reconstruction (§5). The correct
"full decompress" invocation requires the `outreads` argument:

```bash
# WRONG — dumps raw streams only, never reconstructs a read (what the
# original 3.3x figure was almost certainly measured against)
capsule_decode ecoli.capsule outdir/

# RIGHT — real per-read reconstruction, comparable to what `query` does
# for the reads it returns
capsule_decode ecoli.capsule outdir/ outdir/reads.seq

# query, for comparison
capsule_decode query ecoli.capsule out.fa 0-100000
```

Re-measured by hand, 5 repeats each, same archive, back to back:

```bash
DEC=/tmp/claim3_test_out/workdir/capsule_decode
ARC=/tmp/claim3_test_out/workdir/ecoli.capsule
for i in 1 2 3 4 5; do /usr/bin/time -f "%e" "$DEC" query "$ARC" q.fa 0-100000 >/dev/null; done
for i in 1 2 3 4 5; do /usr/bin/time -f "%e" "$DEC" "$ARC" out reads.seq >/dev/null; done
```

Result: mean query 0.862 s, mean full-decompress 1.398 s, **1.62×** — real,
repeatable, and roughly half the previously claimed 3.3×. Selectivity is
unaffected by this correction and remains the stronger number: 11,802 of
1,553,259 reads returned (132× fewer reads, 112× fewer bytes: 2.1 MB vs
235 MB).

---

## 7. Future work

1. **Port `export`/`coverage`/`query` into the outer `arcs` binary.** They
   exist only in `c_star_pg_advance`'s `capsule_decode.cpp` today. Required
   before the paper can describe them as ARCS features rather than a sandbox
   result.
2. **Write `benchmark.sh claim3 ...`**, matching the pattern already used by
   `claim1`/`claim2` — **done for E. coli this session**: `scripts/run_claim3.sh`
   builds both binaries, installs SPAdes if absent, compresses, and runs all
   three operations plus their conventional baselines from one command
   (§6.4). Still open: wiring the GIAB (HG002/HG005) rows into the same
   script — those still require manually pointing it at pre-built archives.
3. **A range-indexed archive layout for O(range) `query`.** The current
   design stores reads as pg slices with no independent offset index; adding
   one (e.g. a sorted `pos_abs` side-table) would let `query` skip the full
   pg rebuild for narrow ranges. This is a layout change, not a bug fix — not
   attempted in this pass because it changes the archive format, which needs
   its own gated, non-regressing verification pass per this repo's standing
   rules.
4. **Run the bwa+mosdepth coverage baseline on `SRR2584863_ecoli`** — the row
   in `t6_results.csv` still says "timing pending"; only the three GIAB
   windows have a completed coverage-baseline number.
5. **Delete the dead code found in §3** (the superseded `coverage` branch and
   `uid_of`) as its own small, gated cleanup commit.

---

## 8. File index — everything this claim touches

| file | what it is |
|---|---|
| `stages/capsule_decode.cpp` | the implementation — all three operations, lines 154-415; CLI dispatch, lines 656-666 |
| `scripts/run_claim3.sh` | one-command reproducible runner for E. coli — builds binaries, installs SPAdes if absent, compresses, runs export/coverage/query + baselines |
| `results/claim3/t6_results.csv` (outer repo) | canonical result numbers, all rows in §6 (query row superseded — see §6.4) |
| `docs/CLAIM3_PRIOR_ART.md` | full prior-art survey (this doc's §4 is its locked summary) |
| `docs/PROJECT_AUDIT.md` | cross-claim audit; Claim 3 section superseded by this document |
| `docs/CLAIM3_LOCKED.md` | this file |
| `/data/fastq/SRR2584863_1.fq` | E. coli reads, locked accession #1, `DATASET_LOCKED.md` |
| `/data/fastq/HG002_pooled.fq`, `HG005_pooled.fq` | GIAB pooled chr20 reads, source of the `r2`/`r5`/`r3` windows |
| `/root/spades_ecoli/contigs.fasta`, `/root/spades_ecoli_run.log` | this session's SPAdes run: output assembly + full timed log |
| `/root/SPAdes-4.0.0-Linux/` | installed SPAdes binary (v4.0.0), this session |

**Commit map for this claim:**

    2089c9b  CLAIM 3 (ADDRESSABLE) implemented: export / coverage / query served from the archive
    abb9cce  Claim 3: fix a coverage bug and remove an under-optimisation; E. coli passes the spec at 291x
    2b5437a  Claim 3 algorithmic audit: 20% coverage undercount fixed, plus region-boundary mislabelling
    95f9012  Claim 3 prior-art survey: random access to compressed FASTQ is NOT novel
    b17b8db  Claim 3 prior art: make explicit no prior tool has export+coverage+query
    a7987f8  Claim 3: run SPAdes (spec-exact export baseline), 555x on E. coli, closes deferred gap
    (outer repo) 7149331  Claim 3 T6 result summary + Kmer2SNP SAM-to-VCF conversion script

---

## 9. Verdict

**Locked.** All three operations exist, are measured on real data, beat the
spec's stated targets by a wide margin (export 254–555× against a ≥40× bar,
coverage 23–33× against a 2–5× bar), were audited end-to-end with two real
bugs found and fixed (not zero found, which would be the suspicious
outcome for a first algorithmic audit of new code), and have a specific,
per-operation prior-art position that does not overclaim. The open items in
§7 are genuine future work — porting to the outer binary, a benchmark
script, and a range index — none of which block the claim as stated, and all
of which are scoped precisely enough to pick up without re-deriving context.

---

## Addendum — HG002 chr20, measured 2026-09-08, and what the speedup is not

Full end-to-end run `results/sanity_HG002_20260908_110429`, 12.6 M reads,
3.99 GB input, on an otherwise idle box (nothing else timed alongside it).

| operation | ours | baseline | ratio |
|---|---|---|---|
| `export`   | **5.491 s** | SPAdes 4.0.0 `-t 12` **2675.80 s** | **487.3x** |
| `coverage` | **2.791 s** | bwa mem + samtools sort + mosdepth **162.02 s** | **58.1x** |
| `query`    | **6.233 s** | no competitor exists | — |

The baselines are real work, not early exits: bwa reports 1,586 s of CPU
across 12 threads, and mosdepth reports mean depth **28.48x** over all
63,025,520 bases of chr20 — the expected ~30x. SPAdes ran to completion.

**The `export` ratio must not be read as "the same output, 487x faster."**
It is not the same output, and the difference is visible in one number:

    SPAdes contigs.fasta   24,281 records   63,724,470 B
    our export contigs.fa       2 records  181,556,800 B

SPAdes performs de novo de Bruijn assembly and yields biological contigs.
`export` writes the pseudogenome as two FASTA records (`capsule_pg_main`,
`capsule_pg_second`) — the greedy suffix-prefix-overlap assembly the compressor
already built, materialised from `literal` + `mem_triples` (section 2.1). The
honest statement of what is measured is:

> the time to obtain a reference-free coordinate system over the reads,
> from an archive that had to be written anyway, versus the time to build
> one de novo.

That is the claim Claim 3 actually makes — ADDRESSABLE, not "better
assembler". The pseudogenome is not a substitute for `contigs.fasta` for any
task needing contig boundaries or biological structure, and no table in this
project should imply otherwise. `coverage` and `query` carry no such caveat:
they produce the same object as the conventional route (per-base depth, reads
over a coordinate range), which is why 58.1x is the more directly comparable
of the two ratios.

**Cost of making the archive addressable**, measured on this run: the
`contig_spans` stream is 232,509 B and the archive grew 232,530 B (stream plus
container header) against a run without it — **0.041%** of a 573,767,964 B
archive. Claim 2's archive-path calling and Claim 3's addressability are paid
for out of four hundredths of one percent.

---

## Addendum 2 — what Claim 3's novelty actually is, verified 2026-09-09

Two statements this document and the scripts had been carrying were false and
are removed: "no competitor exists" (SPRING `--decompress-range`, Genozip
`--head=N`, BEETL-fastq sequence search, CRAM all do partial retrieval), and
any framing that leads with speed (`genocat --head=100` is 0.19 s against our
0.46 s -- Genozip is faster at raw extraction).

### Verified by running the tools, not by reading source

**PgRC2 offers no read-out but full decompression.** From `PgRC -h`: compress,
`-d` decompress, `-o`, `-t`, and expert tuning flags. No export, query,
coverage or region extraction. Note carefully: its source DOES keep per-read
positions (`pseudogenome/readslist/`), so the correct claim is that it lacks
the INTERFACE, not the data. A reviewer who knows PgRC2 will check this.

**Retrieval is neighbourhood, not substring -- the real difference from
BEETL-fastq.** 50 random 40 bp probes from the E. coli pseudogenome:

    reads returned that do NOT contain the probe
      median 38%   IQR 26-46%   min 0%   max 100%
      excluding 4 repeat probes: 46 probes / 1692 reads -> 35.9%

    The pooled figure (1.9%) is misleading: 4 probes landed in repeats and
    returned >1000 reads each, nearly all containing the probe. Report the
    median.

A BWT text index returns only reads *containing* the query. The ~38% it cannot
return are the reads overlapping the locus without spanning the probe -- the
ones at the region's edges, which carry the alleles a caller needs.

**Query cost is position-independent.** 0.46-0.50 s at offsets 0, 1 Mb, 3 Mb
and 26.76 Mb of a 26.96 Mb pseudogenome. `genocat --head=N` is not: 0.19 s at
N=100 rising to 0.83 s at N=2,000,000 against a 0.85 s full decode, because it
reads sequentially from the start. SPRING's `--decompress-range 1 100` costs
4.93 s against a 5.28 s full decode -- 93-99% of decoding everything.

**Sequence query, added 2026-09-09.** `capsule_decode query` now accepts a DNA
string as well as `START-END`. This matters because our offsets are
pseudogenome offsets, not chromosome positions, so "give me BRCA1" previously
required reintroducing a reference. Querying by sequence removes the coordinate
system from the interface: the user supplies the gene or primer they already
hold. The numeric form is byte-identical to before (11,523 reads on the same
range, `cmp`-verified). A chromosome-coordinate overlay was considered and
rejected -- it would need an external reference and break the reference-free
property.

### The statement that survives

> Everyone else's archive stores what the reads SAY. Ours stores where they
> SIT. The compressor had to assemble in order to compress, so the coordinate
> system is a byproduct rather than an addition -- which is why locus
> retrieval, per-base coverage and reference-free variant calling all come out
> of the same structure, for 0.041% of the archive.
