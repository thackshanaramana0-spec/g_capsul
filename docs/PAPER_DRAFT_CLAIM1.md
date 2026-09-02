# ARCS / CAPSULE — paper draft, Claim 1 (Methods + Architecture)

Written 2026-09-02. This is a **methods-and-architecture draft**, not a full
paper — Results (beyond the numbers Claim 1 already has) and Limitations for
Claims 2/3 are deliberately left as placeholders since those benchmarks are
not finished. Everything under "Claim 1" below is written as if going
straight into a Nature Methods-style submission: precise, falsifiable,
citable. Source material: `docs/TECHNICAL_ARCHITECTURE.md`,
`docs/SOTA_COMPARISON.md`, `docs/FAILURES_AND_REFUTED_IDEAS.md`, and the
project's own measured numbers (`/tmp/allph/r.csv`, `CLAUDE.md`).

**Naming note:** the outer project is branded ARCS in
`/root/arcs-clean/CLAUDE.md`; the sandbox implementation here is called
CAPSULE. For the paper, CAPSULE's sequence/names/quality engine IS what
"ARCS compress" should eventually wrap — see HANDOVER.md's reconciliation
item. Use "ARCS" as the tool name in the paper and "the pseudogenome
assembly stage" / "CAPSULE" only when this draft needs to be specific about
which codebase a claim was measured on.

---

## 0. Paper skeleton (for orientation — do not include this section itself)

1. Abstract
2. Introduction (the 3-claim framing)
3. Methods
   - 3.1 Overview / one-sentence description
   - 3.2 Pseudogenome assembly (Claim 1 core)
   - 3.3 Archive format and stream-specific coding (Claim 1 core)
   - 3.4 Variant calling from the pseudogenome (Claim 2 — placeholder)
   - 3.5 Addressable export/coverage/query (Claim 3 — placeholder)
   - 3.6 Benchmark datasets and competitor invocation
4. Results
   - 4.1 Claim 1 — compression ratio, speed, RAM (data below)
   - 4.2 Claim 2 — variant-calling accuracy (placeholder, pending run)
   - 4.3 Claim 3 — addressability speedups (placeholder, pending run)
5. Discussion / Limitations
6. Data and code availability

Sections 3.1–3.3, 3.6, and 4.1 are written in full below. 3.4/3.5/4.2/4.3 are
left as one-paragraph placeholders describing what WILL go there, since
writing results prose before the numbers exist would misrepresent the work.

---

## 1. Introduction framing (the 3 claims, one paragraph each)

> Lossless FASTQ compression tools are typically single-purpose: they shrink
> a file and can only give it back unchanged. We present ARCS, a single
> archive format built around an explicit **pseudogenome** — an assembled
> approximation of the sample's genome built directly from the reads at
> compress time — that is simultaneously (1) a competitive general-purpose
> FASTQ compressor, (2) a substrate from which heterozygous variants can be
> called directly, without a separate alignment step, and (3) an addressable
> store supporting sub-second contig export, per-region coverage, and
> coordinate-range read extraction without full decompression. We call these
> three properties COMPACT, FAITHFUL, and ADDRESSABLE, and evaluate each
> against the relevant specialist tools: SPRING and Genozip for compression
> ratio (Claim 1), DiscoSNP++ and Kmer2SNP for variant-calling accuracy
> (Claim 2), and SPAdes / BWA+mosdepth for addressability speed (Claim 3).

This is the frame the whole paper hangs on. Claim 1 is the only one with a
complete, verified dataset as of this writing (2026-09-02); say so plainly in
the introduction rather than implying all three are equally mature — Nature
Methods reviewers will check.

---

## 2. Methods — 3.2 Pseudogenome assembly

This is the part of the paper that must be exact enough for a competent
bioinformatician to reimplement without reading the source. Write it as a
numbered pipeline, each step named, each parameter given with its default
and how it is chosen.

### 3.2.1 Definitions

- A **read** is one sequence record from the input FASTQ (`_1` mate only for
  paired data in this benchmark — state this explicitly, it is a scope
  decision, not an oversight).
- The **pseudogenome (pg)** is a single contiguous byte buffer, built
  incrementally, that every read is ultimately expressed relative to: either
  as a literal (bytes newly appended to pg) or as a reference
  `(dst, src, len, rc)` meaning "this read equals pg[src : src+len],
  reverse-complemented if rc, placed so that decoding walks pg[dst : dst+len]".

### 3.2.2 Step 1 — exact suffix-prefix chaining

For every read, using a `k`-mer seed index (`SEEDW`-mers, position-aware), find
another read whose prefix exactly matches this read's suffix over at least
`MINOV` bases, greedily chaining `read_i, read_{i+1}, ...` into one
contiguous pg span. Sweep starts at `L = Lmax` (maximum read length in the
file), not `Lmax - 1` — **an exact duplicate read only overlaps itself at the
full read length**, so starting one base short makes duplicate detection
structurally invisible and forces every duplicate into a separate pre-pass.
Both `MINOV` and the seed width are dataset-swept (grid or golden-section
search, `GSEARCH=1`), not fixed constants, and the paper should report the
swept range and the selection criterion (smallest resulting archive) rather
than a single "best" value, since different datasets pick different optima.

### 3.2.3 Step 2 — pigeonhole mapping

Reads that fail to chain are mapped onto the pg built so far: an exact
`k`-mer seed hit inside pg is extended and accepted if it matches within a
bounded mismatch tolerance (`MAXMM`, default 3 per read). Accepted mappings
become references; rejected candidates leave the read in the leftover pool.
`MAXMAP`, the per-read candidate ceiling, widens automatically as a function
of `leftover_frac` (fraction of reads still unmapped after step 1) once that
fraction exceeds 0.60 — report this formula explicitly, it is the paper's
one adaptive/input-dependent parameter and reviewers will ask what makes it
"algorithmic" rather than tuned: it is keyed on a property of the CURRENT
input, computed fresh every run, not fit to any specific dataset.

### 3.2.4 Step 3 — second-region assembly

Reads that still fail step 2 are not stored as literal fallback. They are run
through their own instance of steps 1–2, producing a second contiguous pg
region, appended after the first (the boundary is `main_pg_end`, stored in
the archive header). This converts the hardest-to-place reads — the ones
most likely to represent low-coverage or repeat-adjacent regions — into
additional compressible sequence instead of raw bytes.

### 3.2.5 Step 4 — MEM self-match

Both regions, once built, are matched against themselves for
maximal exact matches (MEM) at any offset, not just chain-adjacency, using a
copMEM-style seed-and-extend index. Every accepted MEM removes its covered pg
bytes from the literal stream and adds one `(dst, src, len, rc)` reference,
identically shaped to a pigeonhole reference. Acceptance is **cost-aware**:
a MEM is kept only if its exact encoded cost (computed, not estimated) beats
the literal cost ceiling (2 bits/base) it would otherwise occupy. State the
one negative result plainly in Methods, since it is informative: matching
each region against **itself alone** (as opposed to the whole assembled pg)
was tested and rejected — general-purpose LZMA on the literal stream already
recovers that redundancy at lower marginal cost than an explicit reference
would (see Limitations).

### 3.2.6 Step 5 — mismatch-tolerant extension (reported as OFF by default)

A MEM/reference match may optionally extend past its first mismatch up to a
tolerance (`REF_MAXMM`, default 4 tolerated substitutions per match),
recording `(position, observed base)` per tolerated mismatch instead of
terminating the match. This is implemented and measured, and the paper should
report the negative result: on the dataset tested, this cost more in
mismatch-stream bytes than it saved in fewer/longer references
(+490,763 B), so it ships OFF. This is exactly the kind of ablation reviewers
expect and it is already done — do not omit it for looking "cleaner".

### 3.2.7 N-handling

`N` bases are substituted with `A` before assembly and corrected afterward
from three side streams (`n_pos`, `n_indices`, `n_cnt`) so that no assembly,
mapping, or mismatch-coding logic needs an alphabet larger than {A,C,G,T}.

---

## 3. Methods — 3.3 Archive format and stream-specific coding

### 3.3.1 Container

A self-identifying, named-stream container (magic, version, `pg_len`,
`main_pg_end`, `minmem`, `n_streams`, then `(namelen, name, len, payload)`
per stream). Streams are looked up by name at decode time, not by position —
state this as a design choice for extensibility (adding a stream, e.g.
quality, cannot perturb any other stream's offset).

### 3.3.2 Stream inventory (Table — include as a paper table, not prose)

Reproduce the "archive streams" table from `docs/TECHNICAL_ARCHITECTURE.md`
§3 (stream name / contents / coder), trimmed to the columns relevant to
Claim 1 (drop `names_*`/`qual_*` rows from THIS table; they belong in a
"whole-file" supplementary table since Claim 1 is sequence+order only).
Key rows to keep: `literal`, `mem_triples`, `mem_dstgap`/`mem_len`/`mem_rc`,
`pos_abs`, `pos_strand`, `read_lengths`, `orig2uid_flags`/`orig2uid_vals`,
`mm_ref`/`mm_obs`, `mm_pos`, `mm_cnt`, `n_pos`/`n_indices`/`n_cnt`.

### 3.3.3 The general-purpose selector

For every stream without a bespoke coder, `best_encode` tries LZMA, PPMd7,
FSE (tANS), an in-house adaptive range coder, a u32 byte-plane split, and a
chunked variant, keeping whichever is smallest, with one method byte
recording the choice for the decoder. State the design rationale: different
streams have wildly different statistics (skewed small alphabets vs
near-random 32-bit deltas vs mostly-constant), and no single general-purpose
coder dominates across all of them — this is measured, not assumed (cite the
per-stream ablations in Supplementary Table X, sourced from
`docs/FAILURES_AND_REFUTED_IDEAS.md`).

### 3.3.4 The one specialist coder worth naming in Methods: mismatch symbols

`mm_ref`/`mm_obs` use an ADAPTIVE model keyed on the reference base (a
4-symbol context, each conditioning a 3-symbol "which of the other 3 bases"
alphabet), not an independent code per mismatch. This is the single largest
per-stream margin found against PgRC2 (124,280 B vs 208,234 B on identical
input, a 40.3% reduction) and is worth its own paragraph and possibly its own
supplementary figure (context-conditioned entropy vs order-0 entropy vs
achieved bytes).

### 3.3.5 The two explicit zero/nonzero splits

`orig2uid_flags`/`orig2uid_vals` and `mm_cnt_flags`/`mm_cnt_vals` each split a
stream into a dense boolean flag plus a sparse value payload because a
majority of values are exactly zero (e.g., 79.55% of `orig2uid` deltas on
E. coli) and no single distribution fits both the all-zero majority and the
long-tailed minority well. Report this as a general technique
("zero-inflated stream splitting"), not a dataset-specific hack — it is
tested and applied identically wherever the same statistical shape appears.

### 3.3.6 What is deliberately NOT built for Claim 1

State explicitly, in Methods or Limitations: no minimizer-based approximate
seeding (tested at an earlier stage, no net win here — cite
`docs/FAILURES_AND_REFUTED_IDEAS.md` if it lists this); no PgRC2-style 3-way
physical pg split (measured unnecessary, `docs/DO_WE_NEED_THEIR_3WAY.md`); no
approximate/lossy mode. Every design choice in this section has a citable
ablation behind it — this is a genuine strength of the current internal
documentation and should be leaned on for the Methods section's rigor.

---

## 4. Methods — 3.6 Benchmark datasets and competitor invocation

### 3.6.1 Datasets

15 SRA accessions spanning Bacteria (6), Virus (2), Fungi (2), Protista (2),
Archaea (2), Plantae (1) — no Animalia at full coverage (disclose this as a
scope limitation, with the reasoning: no small-genome, high-coverage
Animalia SRA run was found; see `NEW_DATASET_LOCKED.md`). List the full table
(accession, organism, kingdom, `_1` coverage) as Supplementary Table 1 —
already written in `NEW_DATASET_LOCKED.md`. State that all 10 of the original
ARCS project's locked accessions plus the 5 extension accessions are drawn
from public SRA/ENA with no re-derivation after locking (the project's own
`DATASET_LOCKED.md` rule).

### 3.6.2 Competitors and exact invocation

- **SPRING**: `spring -c -i in.fq -o out.spring -t <nproc> -g`; for the
  sequence-only phase, `--no-ids --no-quality` isolates the comparable column.
- **Genozip**: `genozip --force -o out.genozip in.fq`; sequence-only cost is
  read from `genocat --STATS` (capital — exact byte accounting, not the
  rounded `--stats`), not by differencing archives with/without a flag
  (differencing is demonstrably wrong for Genozip specifically, since it
  reuses `length=` internally — this methodological point is worth its own
  sentence in Methods, since a naive reviewer might suggest differencing).
- **PgRC2**: default `CODER_LEVEL_NORMAL` settings, run on the 8 (of 15)
  datasets it can process at all — report its 6 hard failures
  (unsupported-variable-length rejections and crashes) as a capability
  result, not just excluded rows.

### 3.6.3 Verification protocol (this is the paper's most important
methodological safeguard, given the project's own history — see §7 below)

Every reported archive size is preceded by an actual decode-and-compare
against the original file (sequence column at minimum; full 4-line FASTQ,
byte-for-byte MD5, for the whole-file condition). State plainly that an
earlier internal round of measurements reported sizes from archives that had
NOT been round-trip verified at the entropy-coding layer, and that four
specific silent data-loss defects were found and fixed as a direct result of
instituting this check — this is a methods-rigor point that strengthens the
paper's credibility rather than weakening it; state it in Methods, not bury
it. See `docs/FAILURES_AND_REFUTED_IDEAS.md` Part A for the technical
detail, cited in Supplementary Methods.

---

## 5. Results — Claim 1 (numbers as of 2026-09-02, 12/14 datasets confirmed
at time of writing; treat the remaining rows as "pending" until the
in-progress run finishes)

### 5.1 Sequence + read order only (Phase 1, vs PgRC2, 7 comparable datasets)

Aggregate: ARCS 81,631,156 B vs PgRC2 83,192,412 B — **+1.88% smaller**,
6 wins / 1 loss (S. acidocaldarius, −0.83%, confirmed structural, not a
tuning miss — see Table 1 of `docs/SOTA_COMPARISON.md`). PgRC2 additionally
fails outright on 6 of 14 datasets attempted here (variable-length rejections
and crashes) — report this as a capability result alongside the ratio.

### 5.2 Sequence + read order (Phase 1, vs SPRING/Genozip, 12 datasets
confirmed, 2 more in progress)

| dataset | ARCS (P1) | SPRING (P1) | Genozip (P1) | ARCS vs SPRING | ARCS vs Genozip |
|---|---|---|---|---|---|
| ERR5181310 (SARS-CoV-2) | 836,191 | 931,840 | 909,869 | −10.3% | −8.1% |
| SRR554369 (P. aeruginosa) | 8,981,037 | 9,123,840 | 39,033,286 | −1.6% | −77.0% |
| ERR552797 (M. tuberculosis) | 4,716,998 | 6,983,680 | 37,540,454 | −32.5% | −87.4% |
| SRR2584863 (E. coli) | 8,225,993 | 10,506,240 | 52,359,417 | −21.7% | −84.3% |
| SRR29296997 (H. salinarum) | 2,756,247 | 4,352,000 | 13,736,390 | −36.7% | −79.9% |
| ERR12954017 (S. acidocaldarius) | 3,141,277 | 4,331,520 | 26,761,465 | −27.5% | −88.3% |
| SRR40402583 (C. jejuni) | 3,806,837 | 3,624,960 | 25,633,015 | +5.0% | −85.2% |
| SRR40271341 (H. pylori) | 3,853,956 | 6,481,920 | 22,349,357 | −40.5% | −82.8% |
| ERR17740259 (S. aureus) | 13,522,261 | 16,332,800 | 84,382,969 | −17.2% | −84.0% |
| SRR37283774 (P. falciparum) | 17,161,828 | 18,083,840 | 40,518,738 | −5.1% | −57.6% |
| DRR976266 (S. cerevisiae) | 21,822,808 | 24,555,520 | 164,037,369 | −11.1% | −86.7% |
| SRR36741279 (L. major) | 27,954,749 | 29,214,720 | 110,192,473 | −4.3% | −74.6% |

11/12 wins vs SPRING (1 loss, C. jejuni, +5.0% — a variable-length dataset;
worth its own sentence discussing why variable-length data is the harder
case for ARCS's fixed-width position coding). 12/12 wins vs Genozip. These
numbers are freshly re-measured post-bug-fix (2026-09-02) and are the ones
to cite; do not cite `docs/PHASE2B_RESULT.md`'s earlier figures, which are
marked void.

### 5.3 Whole-file (Phase 3: sequence + order + names + quality) vs
SPRING/Genozip, same 12 datasets

| dataset | ARCS (P3) | SPRING (P3) | Genozip (P3) | ARCS vs SPRING | ARCS vs Genozip |
|---|---|---|---|---|---|
| ERR5181310 | 8,686,774 | 9,390,080 | 9,218,949 | −7.5% | −5.8% |
| SRR554369 | 57,320,645 | 59,166,720 | 88,813,037 | −3.1% | −35.5% |
| ERR552797 | 46,964,185 | 52,101,120 | 82,799,524 | −9.9% | −43.3% |
| SRR2584863 | 68,677,977 | 74,045,440 | 119,618,198 | −7.2% | −42.6% |
| SRR29296997 | 15,654,489 | 17,500,160 | 27,233,618 | −10.5% | −42.5% |
| ERR12954017 | 15,159,145 | 16,967,680 | 39,545,150 | −10.7% | −61.7% |
| SRR40402583 | 9,040,464 | 9,543,680 | 30,863,017 | −5.3% | −70.7% |
| SRR40271341 | 40,620,675 | 45,742,080 | 62,259,448 | −11.2% | −34.7% |
| ERR17740259 | 83,421,422 | 92,303,360 | 164,126,889 | −9.6% | −49.2% |
| SRR37283774 | 67,382,606 | 71,168,000 | 94,989,854 | −5.3% | −29.1% |
| DRR976266 | 56,684,465 | 61,716,480 | 201,140,420 | −8.2% | −71.8% |
| SRR36741279 | 106,342,010 | 117,565,440 | 194,786,333 | −9.5% | −45.4% |

**12/12 whole-file wins against BOTH SPRING and Genozip.** This is the
headline Claim 1 result and should be the paper's primary table — every
single-column loss (C. jejuni vs SPRING, sequence-only) disappears once
names+quality are added, because ARCS's names/quality margins are large
enough to dominate the aggregate.

### 5.4 Speed and RAM (qualify honestly)

ARCS/CAPSULE is ~1.6–1.7× slower and ~2.5–2.8× heavier in peak RAM than
PgRC2 on the datasets measured so far (E. coli). No head-to-head wall-clock
numbers against SPRING/Genozip have been collected yet in this run — get
them before submission; the CSV columns exist (`/tmp/allph/r.csv` currently
tracks size only, not time — note this as an action item, see §7).

---

## 6. Placeholders for Claims 2 and 3 (do not write results prose yet)

**3.4 / 4.2 — Claim 2 (FAITHFUL).** One paragraph in Methods can already be
written: variant calling proceeds by treating the assembled pseudogenome as
a draft reference, aligning original reads back to it (or reusing the
mismatch streams already computed during assembly — state which, once
decided) and calling heterozygous SNVs/indels from the resulting pileup,
compared against GIAB v4.2.1 truth VCFs for HG002–HG005 chr20 at
standardized 30×. Do not write a Results subsection until `t3_snv_f1.csv`
exists — the outer project's Claim 2 pipeline is not yet run end-to-end
(missing Kmer2SNP wrapper, per `HANDOVER.md`).

**3.5 / 4.3 — Claim 3 (ADDRESSABLE).** One paragraph: because ARCS's archive
already stores per-read placement (`pos_abs`) and the assembled pg directly,
contig export (`arcs export`), per-contig coverage (`arcs coverage`), and
coordinate-range extraction (`arcs query`) can all be served from the
archive without full FASTQ reconstruction, and should be faster than
SPAdes/BWA+mosdepth by construction. No numbers exist yet
(`t6_results.csv` not produced) — do not write Results until Phase 4 of the
outer project's benchmark runs.

---

## 7. Discussion / Limitations (draft points)

- Speed and RAM are real, disclosed trade-offs, not hidden ones — assembly-
  based compression costs time the read-relative competitors do not pay.
- No Animalia representation at full coverage in the 15-dataset set; a
  genuine scope gap, not a cherry-pick (explain the search that failed to
  find a small-genome high-coverage Animalia SRA run).
- One structural loss vs PgRC2 (S. acidocaldarius) exists and is disclosed,
  confirmed not a tuning artifact.
- One structural loss vs SPRING exists at the sequence-only phase (C. jejuni,
  variable-length reads) that disappears once names+quality are added — worth
  a sentence on why fixed-width position coding is the likely cause and
  whether a variable-width fallback is future work.
- Quality coding: ARCS uses a vendored, unmodified fqzcomp/htscodecs core
  (BSD-licensed) rather than a novel context model; the paper should be
  honest that this component's algorithmic contribution is "correct
  integration into a whole-FASTQ addressable archive," not "a better
  quality codec" — the project's own internal ablations found no context
  fqzcomp doesn't already exploit that beats it (tile, base-identity, is-N
  flag all tested and refuted, see `docs/FAILURES_AND_REFUTED_IDEAS.md`
  §B.18 for the held-out-vs-in-sample methodology, which is itself worth a
  sentence — it is a real methodological point about entropy-estimation
  pitfalls that a Methods reviewer would appreciate seeing handled
  correctly).
- Names coding beats SPRING (6/8 to 7/8, correct-methodology measurement)
  and Genozip (6/8) on non-human datasets but this project deliberately used
  a differencing-vs-per-tool-section-accounting methodological correction
  mid-project — worth one sentence acknowledging the correction happened,
  since the corrected numbers are what's being published.

---

## 8. ~25 additional points the paper will need (checklist, not prose)

1. Abstract needs one number per claim — Claim 1's is ready (+1.88% vs
   PgRC2, 12/12 whole-file wins vs SPRING/Genozip); Claims 2/3 are not.
2. A results figure (bar chart, ARCS vs SPRING vs Genozip, log-scale bytes,
   12 datasets) — data is ready in §5.2/5.3 above.
3. A per-stream breakdown figure/table for at least one dataset (e.g. E.
   coli) showing where bytes go (`pos_abs` / `mm_pos` / `literal` dominate,
   per `CLAUDE.md` §8) — useful for the mechanism story.
4. Wall-clock and peak-RSS numbers vs SPRING/Genozip specifically (currently
   only measured vs PgRC2) — an open action item, not yet collected.
5. A formal statement of what "lossless" means here (byte-identical 4-line
   FASTQ reconstruction, MD5-verified) and how it's tested — methods
   reviewers will ask.
6. Justify the `_1`-mate-only scope decision for paired data explicitly (is
   this a real limitation for paired-end genomic use? say so if it is).
7. A supplementary ablation table listing every "tested and rejected" idea
   (already written in full in `docs/FAILURES_AND_REFUTED_IDEAS.md` — needs
   trimming/reformatting for a supplement, not rewriting).
8. Explicit license disclosure for vendored code (fqzcomp/htscodecs BSD
   3-clause; LZMA SDK / PPMd public domain; FSE BSD) — required by most
   journals' code-availability policy.
9. A statement on reproducibility: exact tool versions (SPRING, Genozip,
   PgRC2 commit/release), exact command lines (already drafted in §4 above),
   hardware spec (12 vCPU / 90 GB RAM / 250 GB SSD, Ubuntu 24.04).
10. Data availability: all 15 accessions are public SRA/ENA; state the
    accession table as Supplementary Table 1, no restricted data in Claim 1.
11. Related-work paragraph placement: SPRING, Genozip, PgRC2, mstcom (dismiss
    briefly, 2/15 wins per `competitor_landscape_mstcom` memory) — don't
    over-claim "sequence SOTA"; PgRC2 already beats SPRING more than ARCS
    beats PgRC2 in the published literature (§ Discussion point above).
12. State the "why pseudogenome, why not just read-relative like SPRING"
    architectural motivation up front in Methods — it's what unlocks Claims
    2/3, and a reviewer will ask why bother with assembly for compression
    alone if the ratio margin over SPRING/Genozip were the only payoff.
13. Consider whether a schematic figure (pipeline diagram: reads → chaining
    → pigeonhole → second region → MEM self-match → streams → archive) is
    needed — very likely yes for Methods.
14. Define "coverage" precisely as used in the dataset table (`_1`-file
    read-bases / genome-length) since it's a specific, stated convention,
    not the standard paired-depth definition.
15. Decide and state the human-data policy explicitly (excluded from Claim 1
    by the locked-accession ban list; included at chr20/30× only for Claim
    2) — a reviewer will ask about human applicability.
16. A short paragraph on the four bugs (§ Methods 3.6.3) doubles as a
    "software correctness" section some journals now expect for
    computational tools — consider a dedicated "Software Quality Assurance"
    subsection rather than folding it into Methods.
17. Statistical significance / variance: current results are single runs per
    dataset-tool pair; state whether repeat runs were done (they were not,
    per the log) and whether that's acceptable for a deterministic
    compressor (it likely is — compression size is deterministic given fixed
    input and flags — but wall-clock timing claims need repeats or explicit
    single-run caveats).
18. Genozip's anomalous fungi numbers (still unexplained per `CLAUDE.md`
    §6.2) must be resolved or explicitly caveated before publication — an
    unexplained 6-7x outlier in a competitor's own tool looks like a
    methodology error to a reviewer even if it isn't.
19. A "future work" paragraph: closing the C. jejuni variable-length loss,
    Claim 2/3 completion, human WGS scale.
20. Code availability statement + exact commit hash to cite
    (`c_star_pg_advance`, branch `c_star_pg_advance`, plus the outer
    `/root/arcs-clean` binary once reconciled — see HANDOVER.md item 1).
21. A methods paragraph on parameter selection philosophy (§3.2.2's
    "swept, not fixed" point) — ties into standing rule 1 in `CLAUDE.md` and
    pre-empts a "did you overfit hyperparameters to your test set" question.
22. Consider whether MINOV/SEEDW sweep ranges should be reported per-dataset
    in a supplementary table (transparency) rather than only the winning
    value.
23. A clear figure or table contrasting Phase 1 / Phase 2b / Phase 3 (already
    drafted as Table 3 in `docs/SOTA_COMPARISON.md`) — shows the reader
    exactly what's being compared at each stage, since competitor tools
    don't expose per-column archives the same way.
24. Ethics/consent statement for GIAB and any human-adjacent data used in
    Claim 2 (standard GIAB consent language, publicly available).
25. Software/hardware determinism caveat: LZMA/xz thread count affects
    byte-exact output in some configurations — state the exact thread count
    used (`nproc`) since this affects reproducibility of exact archive sizes.
26. Consider a runtime-complexity statement (Methods): assembly is
    approximately linear in total bases with an index-lookup constant factor;
    MEM self-match is the potentially superlinear stage — worth stating
    Big-O or at least empirical scaling behavior across the dataset-size
    range already measured (2.6 GB to ~2.2 GB — actually check the largest,
    DRR976266 raw=2,248,042,490 B, is the biggest in the current 12).

---

## 9. Who to cite (grouped)

**Direct competitors / baselines (must cite):**
- SPRING — Chandak, Tatwawadi, Ochoa, Hernaez, Weissman, *Bioinformatics*
  2020, "SPRING: a next-generation compressor for FASTQ data."
- Genozip — Lan, Lam, et al. (Divon Lan et al.), *GigaScience* / *Bioinformatics*
  (verify exact venue/year from the tool's own citation file — check
  `genozip --help` or its GitHub README for the canonical BibTeX before
  submission, do not guess the year).
- PgRC2 / PgRC — Kowalski, Grabowski, et al., the paper this whole project
  benchmarks against by name (verify exact title/venue — it is referenced
  throughout this project's docs as "PgRC2's own paper"; pull the exact
  citation from `docs/HOW_PGRC2_CODES_REFERENCES.md` or the PgRC2 repo
  before submission).

**Quality coding:**
- fqzcomp / htscodecs — Bonfield, "CRAM 3.1: extreme compression of genomic
  data," or the specific fqzcomp paper if cited separately; the project's own
  notes reference "CRAM 3.1, Bonfield, *Bioinformatics* 38(6), 2022" — verify
  this exact reference before submission (this draft carries it forward from
  `docs/TECHNICAL_ARCHITECTURE.md` §6.1 but it has not been independently
  re-verified against the actual published bibliographic record).

**Variant calling (Claim 2, when written):**
- GIAB — Zook et al., "Extensive sequencing of seven human genomes to
  characterize benchmark reference materials," *Scientific Data* 2016, and
  the v4.2.1 truth-set update paper.
- DiscoSNP++ — Uricaru et al. / Peterlongo et al. (verify exact author order
  and venue from the tool's own citation).
- Kmer2SNP — cite its own publication if one exists; otherwise cite the tool
  repository directly.

**General-purpose entropy coders used internally (cite as tools, brief):**
- LZMA SDK (public domain, Igor Pavlov) — cite the SDK/spec, not a paper.
- PPMd (Dmitry Shkarin) — cite the algorithm paper if available, else the
  SDK.
- FSE / tANS — Duda, "Asymmetric numeral systems," and Yann Collet's FSE
  implementation (cite both the ANS theory paper and the implementation).

**Underlying assembly concept:**
- Any standard reference on overlap-layout-consensus assembly and MEM-based
  matching (e.g., the copMEM paper, since this project's MEM index is
  explicitly copMEM-style) — pull the exact copMEM citation before
  submission (Kokot, Grabowski, or similar authorship — verify, do not
  guess).

**Do before submission, not now:** every citation above marked "verify" is
a placeholder pulled from this project's own internal docs' informal
references, not independently checked against a bibliographic database.
Confirm exact author lists, years, and venues (e.g., via a literature search
tool) before the citation list is finalized — do not submit with guessed
bibliographic details.

---

## 10. What NOT to do in this paper (guardrails carried over from project rules)

- Do not claim "sequence-compression SOTA" outright — PgRC2 already beats
  SPRING by more than ARCS beats PgRC2 in the published literature
  (`docs/SOTA_COMPARISON.md` Table 2's own honest framing). Claim
  "competitive ratio + a Claim 2/3 capability PgRC2 lacks entirely."
- Do not cite `docs/PHASE2B_RESULT.md`'s original numbers — void, superseded.
- Do not present Claim 2/3 numbers that do not exist yet.
- Do not hide the speed/RAM trade-off or the one structural loss per
  competitor — both are already correctly disclosed internally; keep that
  standard in the paper.
