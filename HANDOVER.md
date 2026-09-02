# HANDOVER — read this first in any new session

Written 2026-09-02. This is the complete project handover: what the server
is, what this project is, what has been done, and exactly how to reproduce
every piece of it. If you are Claude picking this up in a fresh chat with no
memory of prior sessions, read this document top to bottom before touching
any code.

**Companion documents, each cross-linked, read in this order after this
one:**
1. `docs/TECHNICAL_ARCHITECTURE.md` — every layer, every stream, every coder
2. `docs/FAILURES_AND_REFUTED_IDEAS.md` — every bug found+fixed, every idea
   tried+rejected, with numbers
3. `docs/SOTA_COMPARISON.md` — layer-by-layer vs PgRC2, general vs
   SPRING/Genozip, phase-by-phase
4. `docs/PHASE3_RESULT.md` — final head-to-head numbers, once the in-progress
   run completes (see status below)

---

## 1. The server

- Ubuntu 24.04.4 LTS, kernel 6.8.0-134-generic, x86_64
- 12 vCPU, ~82 GB RAM, ~233 GB disk (`/dev/vda2`, mounted at `/`)
- No GPU. All work is CPU-bound C++/Python.
- Disk usage as of 2026-09-02: `/data/fastq` 58 GB (locked datasets),
  `/tmp/newdl` 21 GB (extended-set downloads). A large cleanup on this date
  freed 169.7 GB of regenerable scratch (test extracts, old binaries,
  duplicate downloads) — see the `four_silent_dataloss_bugs_fixed` and
  session log around 2026-09-02 for what was deleted and why it was safe.

---

## 2. The two projects that share this server — do not confuse them

### 2.1 `/root/arcs-clean` — the OUTER project, real ARCS

The actual, real ARCS/CAPSULE product per its own `CLAUDE.md`. Built binary
at `/root/arcs-clean/build/arcs`. This is what the eventual paper claims
apply to, structured as three claims:

- **Claim 1 (COMPACT):** T1 archive size + T2 speed/RAM, ARCS vs SPRING vs
  Genozip, all 10 locked datasets.
- **Claim 2 (FAITHFUL):** T3/T4/T5 — het-SNV/indel calling F1 vs DiscoSNP++
  and Kmer2SNP on GIAB HG002-HG005 chr20.
- **Claim 3 (ADDRESSABLE):** T6 — `arcs export`/`arcs coverage`/`arcs query`
  speed vs SPAdes/BWA+mosdepth.

The outer project's own locked dataset list is `benchmark/DATASET_LOCKED.md`
(the ORIGINAL 10 accessions — do not confuse with the sandbox's 15/17, §3).

**Status of the outer project's own Phase 2 (real ARCS binary) benchmark:**
stalled at 8/10 datasets as of the last check (`results/block1.log`, stale
since Aug 26) — SRR065390 (C. elegans) and SRR870667 (T. cacao) never
finished. All 8 completed ones were `lossless=LOSSLESS`. This has NOT been
re-run since; it is a SEPARATE, UNMERGED code path from the sandbox below.

### 2.2 `/root/arcs-clean/c_star_pg_advance` — the SANDBOX, CAPSULE

**This is where ALL the work described in this document happened.** A
SEPARATE git repository (own `.git`, own remote
`github.com/thackshanaramana0-spec/c_star_pg_advance`, branch
`c_star_pg_advance`, 170 commits as of this writing), gitignored from the
outer repo (added via commit "Ignore the extracted standalone repo
directory"). It is explicitly intended to BECOME the new ARCS — a
from-scratch reimplementation, benchmarked head-to-head against PgRC2 and,
more recently, directly against SPRING and Genozip.

**The sandbox's own `CLAUDE.md`** (inside `c_star_pg_advance/`) is the
authoritative reference for this project's current state, section by
section — read it, not just this handover, for anything not covered here.

**The two real ARCS/CAPSULE binaries do NOT currently share code.** The
sandbox's names/quality/bugfix work described in this handover has NOT been
ported into the outer project's `src/` (real `arcs` binary). This is a real,
open reconciliation gap — see §8.

---

## 3. The locked datasets — exact list, exact paths, exact download method

**Authoritative file:** `c_star_pg_advance/NEW_DATASET_LOCKED.md` — supersedes
the older `DATASET_LOCKED.md` in the same repo for the 15-dataset
SPRING/Genozip comparison. One swap was made 2026-09-02 (Drosophila removed,
Utricularia gibba added as the Plantae entry) — read that file for why.

### 3.1 The 15 (14 currently usable)

| # | accession | organism | kingdom | path on disk | status |
|---|---|---|---|---|---|
| 1 | SRR2584863 | E. coli B REL606 | Bacteria | `/data/fastq/SRR2584863_1.fq` | on disk |
| 2 | ERR552797 | M. tuberculosis H37Rv | Bacteria | `/data/fastq/ERR552797_1.fq` | on disk |
| 3 | SRR554369 | P. aeruginosa PAO1 | Bacteria | `/data/fastq/SRR554369_1.fq` | on disk |
| 4 | ERR5181310 | SARS-CoV-2 | Virus | `/data/fastq/ERR5181310_1.fq` | on disk |
| 5 | ERR17740259 | S. aureus | Bacteria | `/data/fastq/ERR17740259_1.fq` | on disk |
| 6 | DRR976266 | S. cerevisiae | Fungi | `/data/fastq/DRR976266_1.fq` | on disk |
| 7 | SRR36741279 | Leishmania major | Protista | `/data/fastq/SRR36741279_1.fq` | on disk |
| 8 | SRR37283774 | P. falciparum | Protista | `/data/fastq/SRR37283774_1.fq` | on disk |
| 9 | SRR32429602 | Human betaherpesvirus 5 (HCMV) | Virus | `/tmp/newdl/SRR32429602_1.fastq` | on disk |
| 10 | SRR39257532 | Aspergillus fumigatus | Fungi | `/tmp/newdl/SRR39257532_1.fastq` | on disk |
| 11 | SRR29296997 | Halobacterium salinarum | Archaea | `/tmp/newdl/SRR29296997_1.fastq` | on disk |
| 12 | ERR12954017 | Sulfolobus acidocaldarius | Archaea | `/tmp/newdl/ERR12954017_1.fastq` | on disk |
| 13 | SRR40271341 | Helicobacter pylori | Bacteria | `/tmp/newdl/SRR40271341_1.fastq` | on disk |
| 14 | SRR40402583 | Campylobacter jejuni | Bacteria | `/tmp/newdl/SRR40402583_1.fastq` | on disk |
| 15 | SRR10676752 | Utricularia gibba | Plantae | NOT converted | `.sra` downloaded, `fasterq-dump` failed on `disk-limit exceeded`; disk has ~92 GB free now, safe to retry |

**Excluded by design (the "2 big ones"):** SRR065390 (C. elegans, 11 GB,
Animalia) and SRR870667 (T. cacao, 15–22 GB, Plantae) — part of the OUTER
project's original 10, deliberately excluded from the sandbox's 15/17 for
tractability. T. cacao's sandbox run never completed (cancelled ~28 min in).

**Standalone low-coverage controls (NOT part of the 15, kept separately):**
- SRR40104920 (Drosophila, 2.6× coverage) — dropped from the study entirely
  2026-09-02 once Drosophila's normal-coverage entry was also removed
- ERR17716639 (Arabidopsis thaliana, 2.8× coverage) — KEPT as the sole
  low-coverage control, at `/tmp/kd/ERR17716639.fastq`

**Kingdom split, the 15:** Bacteria 6, Virus 2, Fungi 2, Protista 2,
Archaea 2, Plantae 1 (once #15 converts). Animalia has zero representation
by design (both Animalia candidates were either excluded-as-big or
dropped-as-low-coverage).

### 3.2 Exact download method used

Standard SRA toolkit, verified against NCBI eutils before use (not guessed):

```bash
prefetch <ACCESSION> -O /tmp/newdl
fasterq-dump /tmp/newdl/<ACCESSION>/<ACCESSION>.sra -O /tmp/newdl
```

Accession discovery/verification for the extended 7 (accessions 9–15 above)
used NCBI eutils directly:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=<ORGANISM>+AND+wgs[strategy]+AND+illumina[platform]"
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=<ID>&rettype=runinfo"
```

Utricularia gibba specifically was chosen after checking its REAL NCBI
assembly size (100,688,548 bp, assembly U_gibba_v2 — the smallest known
sequenced plant genome) and its REAL computed coverage (~65× from
spots × per-mate length / genome size), not assumed — the first two
Utricularia gibba SRA runs found were checked and one was rejected for being
only ~7.7× coverage (same low-coverage tier as the already-banned
Arabidopsis run).

### 3.3 GIAB (Claim 2 prerequisite, outer project only)

Already downloaded, not part of the sandbox's 15:

```
/data/fastq/HG002_pooled.fq   (NA24385, Ashkenazi son,   3.99 GB)
/data/fastq/HG003_pooled.fq   (NA24149, Ashkenazi father, 3.99 GB)
/data/fastq/HG004_pooled.fq   (NA24143, Ashkenazi mother, 3.99 GB)
/data/fastq/HG005_pooled.fq   (NA24631, Han Chinese son,  6.34 GB)
~/giab_truth/                 (v4.2.1 truth VCFs + confident-region BEDs)
~/refs/chr20.fa                (GRCh37 chr20, ~65 MB)
~/refs/chr20.sdf               (rtg-format index, built via `rtg format`)
```

Streamed from GIAB S3 WGS BAMs via `samtools`, downsampled to 30× chr20
with `seqtk sample --seed 42` for apples-to-apples T3 comparison (source
coverage varies 60–300× across individuals).

### 3.4 Claim 2 tooling (outer project, built 2026-09-02, NOT yet exercised end-to-end)

- **DiscoSNP++**: cloned+built at `~/DiscoSnp/build/bin/{kissnp2,kissreads2,...}`.
  Required a source patch: `thirdparty/gatb-core/gatb-core/thirdparty/kff-cpp-api/kff_io.hpp`
  was missing `#include <cstdint>`, a GCC-13 strictness break in old vendored
  code — patched, builds clean.
- **DSK** (k-mer counter, GATB-based, a Kmer2SNP dependency): built at
  `~/dsk/build/bin/{dsk,dsk2ascii}`, same `cstdint` patch applied
  preemptively.
- **Kmer2SNP**: cloned from the REAL upstream (`github.com/yanboANU/Kmer2SNP`,
  Yang et al. 2020) at `~/Kmer2SNP`. Needs DSK (above), `findGSE` (R package)
  and Python `networkx` — both already installed in conda env `kmer2snp_r`
  (`/root/miniconda3/envs/kmer2snp_r`, verified working: `library(findGSE)`
  and `import networkx` both succeed).
- **Real gap, NOT yet closed:** `benchmark/run_claim2.sh` (outer project)
  expects `scripts/kmer2snp_to_vcf.py` as a 4-arg END-TO-END wrapper
  (`FASTQ ref chrom out.vcf`, i.e. run DSK → findGSE → Kmer2SNP's graph
  matching → BWA → VCF, all in one script). What actually exists on disk is
  `scripts/kmer2snp_sam_to_vcf.py`, restored from a deleted commit
  (`f8edd34`), which is only the LAST step (5-arg: `pairs.snp kmer1.sam
  ref.fa chrom out.vcf`) — it assumes the SNP-pairs and BWA alignment were
  already produced elsewhere. The end-to-end wrapper needs to be written
  before Claim 2 can run via `run_claim2.sh` as currently written.

---

## 4. What "Claim 1" means here, and its status — READ CAREFULLY, there are TWO Claim 1s

This is the single most important disambiguation in this document.

**Outer project's Claim 1** (COMPACT, real ARCS binary, the 10 original
datasets, vs SPRING+Genozip): stalled at 8/10, last run 2026-08-26, not
re-run since the sandbox work described here.

**Sandbox's own "Claim 1"** (its historical name for the sequence+order
comparison, `c_star_pg_advance/CLAUDE.md` §2): **THIS is the one that is
"completely done" as of 2026-09-02** — meaning:

- Sequence + read order: locked, verified, +1.88% vs PgRC2 aggregate
  (6 wins / 1 loss on the 7-dataset core set) — see `claim1_locked_result`
  memory and Table 1 of `docs/SOTA_COMPARISON.md`.
- Names: wired, gated `CAPS_NAMES=1`, 7/8 vs SPRING (−18.37%), 6/8 vs
  Genozip (−7.35%) on the correct measurement method.
- Quality: wired 2026-09-02, gated `CAPS_QUAL=1`, vendored real fqzcomp.
- Line 3: wired 2026-09-02, one-byte file mode.
- **The full 4-line FASTQ round-trips byte-identical (verified by MD5)
  from the archive alone**, on every dataset checked so far.
- **4 silent data-loss bugs were found and fixed the same day** (see
  `docs/FAILURES_AND_REFUTED_IDEAS.md` Part A) — a full audit confirms
  14/14 locked datasets now decode LOSSLESS.

**So: note this down precisely — as of 2026-09-02, the SANDBOX's Claim 1
(full sequence+order+names+quality+line3 FASTQ compression, vs SPRING and
Genozip) is functionally complete and verified lossless on all 14 usable
locked datasets.** What remains is: (a) the 15th dataset (Utricularia gibba)
conversion, (b) the final Phase1/2b/3 benchmark run completing (in progress,
§6), and (c) reconciling this sandbox work back into the OUTER project's
real `arcs` binary, which has none of it yet (§8).

**Do not report "Claim 1 done" to mean the outer project's 3-claim paper
structure is at the Claim 2 stage.** They are different scopes with the same
name, by historical accident, in two different repos.

---

## 5. Competitor tools — exact versions, exact invocations

| tool | location | invocation used |
|---|---|---|
| PgRC2 | `/root/arcs-clean/method_c/build/PgRC` (GPL-3, cloned, never vendored) | `PgRC -o -t <N> -i in.fq out.pgrc` (`-o` = preserve read order, matching CAPSULE's scope) |
| SPRING | `/root/SPRING/build/spring` | Phase 1: `-c -i in.fq -o out --no-ids --no-quality -t N -w <dir>`. Phase 2b: same minus `--no-quality`. Phase 3: no flags. Patched (`src/main.cpp`) to respect `SPRING_KEEP_TEMP=1` for inspecting its own per-stream temp files directly — this is how SPRING's real `id_1.*` stream sizes were verified against `--no-ids` differencing (they agree to 0.45%) |
| Genozip | system `genozip`/`genocat`, v15.0.87 | `genozip --force -o out.genozip in.fq`; sizes read via `genocat --STATS out.genozip` (CAPITAL — gives exact per-context bytes; lowercase `--stats` rounds to 2-3 significant figures and cannot resolve gaps under ~5%) |
| fqzcomp (real, standalone) | built from `/tmp/htscodecs` (CRAM's htscodecs library) at `/tmp/fqzcomp_qual` | `fqzcomp_qual -r -s <0-3> in.q > out.comp` (raw quality-only mode; `-r` strips ASCII-33 internally) |

---

## 6. Current in-flight work (as of this document's writing)

A benchmark script (`/tmp/allphases.sh`, output `/tmp/allph.log` /
`/tmp/allph/r.csv`) is running, producing Phase 1 / Phase 2b / Phase 3 sizes
for CAPSULE, SPRING, and Genozip on all 14 usable datasets, WITH a
round-trip verification gate at every phase level (a phase's number is only
recorded if that phase's own decode reproduced its input — sequence-only for
P1, +names for P2b, full-FASTQ-MD5 for P3). 8 of 14 datasets complete as of
this writing, all passing every phase's round-trip check. See
`docs/SOTA_COMPARISON.md` Table 3 for the interim numbers and
`docs/PHASE3_RESULT.md` for the final table once it completes.

**If you are resuming this session:** check `/tmp/allph.log` first
(`grep -c '^\[ph\]' /tmp/allph.log` against 14 total) before doing anything
else — if it finished, promote the final numbers into
`docs/PHASE3_RESULT.md` and update `docs/SOTA_COMPARISON.md`'s "PENDING"
markers; if the temp files are gone (a machine restart, etc.), the script
itself is reproducible from `docs/SOTA_COMPARISON.md`'s description of what
it does, or re-derive it from `scripts/encode_adaptive.sh` + the competitor
invocations in §5 above.

---

## 7. Key gotchas — read before running anything, all independently verified this session

1. **`capsule_enc reads.fq 3 16` is NOT the locked configuration** and
   silently produces a much larger, still-valid-looking archive. ALWAYS use
   `scripts/encode_adaptive.sh` (supplies the 4-candidate MAXMAP/MINOV sweep
   and the full parameter list).
2. **`DUMP_PERM=1 DUMP_MM=1` are REQUIRED environment variables**, not
   optional debug flags, despite their names — they gate the per-read
   streams (`pos_abs`, `read_lengths`, `mm_*`, `orig2uid_*`) without which
   the archive silently encodes as empty and decodes to zero reads.
3. **`capsule_decode` writes ONE sequence per line, not 4-line FASTQ.**
   Running `awk 'NR%4==2'` on its output samples every 4th line and produces
   garbage — this produced a false "decoder is broken on every dataset"
   alarm earlier in this project's history. Compare its `r` output directly
   to `awk 'NR%4==2' original.fq`, not to the raw file.
4. **`verify_lossless.sh` does NOT test the archive** — it decodes the
   encoder's separately-dumped intermediate `.bin` files, not the
   `.capsule` container. A size or lossless claim based only on this script
   has NOT exercised the entropy-coding layer at all. Always additionally
   run `capsule_decode` on the real archive and diff against the source
   FASTQ before trusting a number.
5. **Never measure a competitor's per-column cost by differencing "with
   real names" vs "with constant names" against Genozip specifically** — it
   reuses `length=` for its own read-length storage, so removing names can
   make its archive BIGGER, giving a negative and meaningless "cost". Use
   `genocat --STATS`'s own per-context byte breakdown instead (see §5).
6. **Vendored C sources (`thirdparty/htscodecs/*.c`) must be compiled with
   `gcc`, not `g++`** — they rely on implicit `void*` conversions C++
   rejects. `scripts/build106.sh` and `scripts/build_decode.sh` already
   handle this correctly (compile to `.o` with gcc, link with g++) — do not
   simply add the `.c` files to a g++ command line.
7. **Any in-sample entropy comparison for a proposed context (tile, base
   call, etc.) is presumptively overfitting** — always use a HELD-OUT,
   INTERLEAVED train/test split (not first-half/second-half if the signal
   being tested might be sequentially ordered in the file, as tile values
   are). See `FAILURES_AND_REFUTED_IDEAS.md` B.18 for a worked example of
   both failure modes.
8. **Standing algorithmic-only rule** (this project's own, stated
   explicitly by the user multiple times across sessions): any fix must be a
   formula over a MEASURED property of the input, never a per-dataset
   special case or a fitted constant, and must not regress any other
   dataset. Every fix in this handover was built and verified this way.

---

## 8. What is genuinely still open

1. **Utricularia gibba (dataset #15)** not yet converted from `.sra` — disk
   space is no longer the blocker (freed 169.7 GB this session), just needs
   `fasterq-dump` re-run.
2. **The in-progress Phase1/2b/3 benchmark** (§6) needs to finish and its
   numbers need promoting into `docs/PHASE3_RESULT.md`.
3. **Sandbox work is NOT reconciled into the outer project's real `arcs`
   binary.** Names, quality, line-3, and all 4 bug fixes exist ONLY in
   `c_star_pg_advance`. The outer project's `src/name_num_codec.h` was last
   touched 2026-07-23, over a month before any of this session's names work
   existed. This is a real, acknowledged gap — porting is unstarted.
4. **The outer project's Claim 2** (`run_claim2.sh`) cannot run end-to-end
   yet — see §3.4, the missing Kmer2SNP wrapper script.
5. **Genozip's fungi results are flagged as anomalous** in the sandbox's own
   `CLAUDE.md` (162 MB where SPRING gets 24 MB) and unexplained — do not put
   them in a paper table until diagnosed.
6. **The low-coverage gap to SPRING** (Arabidopsis 2.8×, the standalone
   control) is disclosed but not fully closed — see `low_coverage_weakness`
   memory. Explicitly, per this project's own rule, do NOT substitute a
   small/low-coverage run for a locked Animalia/Plantae accession; either
   run the real large file or report that kingdom as pending.
