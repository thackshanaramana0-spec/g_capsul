# RECAP — read this first if you are a new session picking this project up cold

This file exists for one purpose: a fresh Claude session (or a fresh human)
with zero memory of this project should read this file top to bottom and come
away understanding what was built, why, what was found, what was wrong and
got fixed, and exactly what state the repo is in right now. It is written
chronologically and thematically together — the story, then the current
frozen truth. Every number in this file is checked against
`benchmark/results/` and `AUDIT.md`, not typed from memory.

**If this file and a result file ever disagree, the result file wins.**
**If this file and `AUDIT.md` disagree on verification status, `AUDIT.md` wins.**
**Canonical tag: `v1.0.2-capsule`. If a later tag exists, that one wins instead.**

---

## 0. What this project is, in one paragraph

**G_CAPSUL** (also called **CAPSULE**) is a from-scratch, solo-built,
lossless FASTQ compressor that also does two things no other tool in the
world does from the same archive: it calls variants reference-free, and it
retrieves the reads at a genomic locus reference-free — all from one retained
representation, with no external reference genome, no external tool, and no
external API called at runtime anywhere in the pipeline. This was verified
directly from the source code, not assumed (see §7).

Built by one person — Thackshanaramana Balashanmugam, age 20, Tamil Nadu,
India — over roughly two weeks (2026-08-28 to 2026-09-10), no lab, no
co-authors, no institutional infrastructure. ~400+ hours, 486+ commits (check `git log --oneline | wc -l` for the current count),
~87,000 lines of C++, ~189,000 words of working documentation.

---

## 1. The three claims, and what each is worth

### Claim 1 — COMPACT (compression)
**19/19 datasets smaller than SPRING (−6.03% aggregate), 19/19 smaller than
Genozip (−43.26% aggregate), 57/57 archives verified LOSSLESS** by
decode-and-diff against the original FASTQ (not a coder-level round trip —
see §5 for why that distinction mattered). Separately, +1.88% smaller than
PgRC2 on the DNA-sequence stream alone (a different, narrower comparison,
because PgRC2 stores no names/no `+`/no quality and cannot reproduce a
FASTQ at all).

*Honest self-assessment:* the ratio number is modest and the paper says so.
But "modest ratio" should not be read as "modest engineering" — PgRC2 is
DNA-only and refuses or crashes on 6 of 14 real-world variable-length
datasets (`Unsupported variable length reads`, and two crashes:
`stack smashing detected`, `free(): invalid next size`). This tool handles
19 independently-coded streams (identifiers, `+` line, quality, sequence,
read order) end to end, including variable-length reads as a first-class
case (prefix-containment dedup, varint mismatch positions, N-reads routed
through the main pipeline). Claim 1's real job is infrastructural: it is
what Claims 2 and 3 stand on. A DNA-only archive with no read-out could not
have supported either.

### Claim 2 — FAITHFUL (reference-free variant calling)
**Called directly from the archive — no FASTQ, no reference genome, no
second assembly pass.** Het-SNV mean F1 **0.876** vs DiscoSNP++ 0.853 vs
Kmer2SNP 0.475, across four real GIAB individuals (HG002–HG005) at 30x
chr20. Het-indel 0.621 vs DiscoSNP++ 0.591 (3 wins, 1 loss — HG005, see §6).
Multi-allelic **17/26** vs DiscoSNP++ **0/26** (complete chr20 census, not a
sample). Synthetic tetraploid (real HG003+HG004 reads) SNV F1 0.897 vs 0.782.

**The mechanism — this is the paper's centre.** Optimal lossless compression
puts each haplotype of a heterozygous site on its own internally-consistent
contig, because that costs fewer bits than one contig plus a column of
disagreements. The ref-allele and alt-allele reads then never share a
coordinate, and the variant is not present in the data structure at all — no
read-out layer, however clever, can recover it from that representation.
This was identified, measured, and corrected (collapse + mismatch-tolerant
re-placement) at **zero cost to compression ratio**.

Ablation, full chr20, proving this is a mechanism and not a fitted filter:

    neither pass          F1 0.431   P 0.967   R 0.278
    re-placement only        0.426     0.962     0.274   <- WORSE than nothing
    collapse only             0.648     0.963     0.488
    both                      0.888     0.956     0.830

Non-additive: +0.217 and −0.005 alone predict +0.212 together; the measured
joint effect is **+0.457**. And independently, the error SHAPE confirms the
cause: without collapse, precision holds at 0.96 while recall collapses to
0.27 — the caller is not mistaken, it is **blind**, exactly what "the
alt-allele reads are on a different contig" predicts. A bad threshold would
have cost precision instead.

### Claim 3 — ADDRESSABLE (reference-free locus retrieval)
`export` 129–784x faster than SPAdes, `coverage` 16–54x faster than
bwa+mosdepth. And the capability nothing else has: retrieving the actual
reads at a genomic locus from a reference-free compressed archive.

**Same mechanism, one layer out.** A heterozygous locus is not one place in
a compression-optimal pseudogenome — it is **N parallel places** (median 4:
two haplotypes x two strands), measured **up to 18.6 Mb apart**. A
coordinate therefore names one of them and returns a single haplotype, with
the variation gone. **This is not a missing feature and no API can patch
it** — the fragmentation is in the representation the compressor chose, not
in the read-out layer, so a richer coordinate call or a better index would
answer the same malformed question faster, not correctly.

Measured on 400 GIAB het-SNV sites across four individuals, same archive,
same sites, both addressing modes:

    addressing mode          both alleles    one allele    neither
    by coordinate            81/400 (20.2%)      316            3
    by content (sequence)   345/400 (86.2%)       55            0

Probe design: 40bp of *reference* sequence ending 6bp before the variant, so
the probe cannot contain the variant and rig its own result.

*Not claimed:* speed (`genocat --head=100` is faster, 0.19s vs 0.46s) — the
claim is that the question can be asked at all, and that it cannot be asked
by coordinate.

---

## 2. Why "no other tool does all three" is a checked claim, not an assumption

This was tested directly against every real candidate, not assumed:

| tool | compresses losslessly | calls variants reference-free | queries a locus reference-free | one representation |
|---|---|---|---|---|
| SPRING / PgRC / PgRC2 / Genozip / NanoSpring | yes | no | no | n/a |
| DiscoSNP++ / Kmer2SNP | no (raw FASTQ in) | yes | no | n/a |
| CRAM | yes | needs external ref | needs external ref | no |
| BEETL-fastq | yes | **no — hands off to external BWA + external Samtools + external human reference** | content match only (different object: returns reads *containing* a string; median 38% of what it returns doesn't even contain the probe) | **no — 3 separate tools glued by a script** |
| population BWT (*Genome Research* 2017) | **no — discards quality scores by design, checked in its own text** | only at *pre-specified* sites (needs an external truth VCF or genotyping array to know what a variant even is) | content (k-mer) only | **no — needs a separate 4.75 TB RocksDB store + external Cortex graphs + external reference for validation** |
| **G_CAPSUL** | **yes, 57/57 verified** | **yes, from the archive alone** | **yes, from the archive alone** | **yes — one `.capsule` file** |

Verified directly in the source, not asserted (checked this session, file by
file): zero `system()`/`popen()`/`exec()`/network calls anywhere in
`stages/106_inprocess.cpp`, `include/caps_caller.h`, or
`stages/capsule_decode.cpp` (the encoder, the caller, and export/coverage/
query all live there). The only `fork()` calls (3, in the encoder) each
continue running the binary's *own* compiled code with a bounded thread
budget — internal parallelism, not shelling out. Zero external reference
genome files ever opened by the core pipeline. One system library is linked
(`liblzma`, standard XZ Utils) and it is used only to decompress raw
sequence bytes off disk — `caps_caller.h` has zero lzma calls, and the
export/coverage/query logic itself calls zero lzma functions of its own.

**What IS borrowed, openly, and does not change the above:** LZMA, PPMd7,
FSE/Huf0 (general-purpose entropy coders, compiled-in libraries — the same
category every real competitor uses; nobody in this field writes their own
arithmetic coder) and fqzcomp (BSD 3-clause, vendored) for one specific
stream: quality-score compression, which feeds Claim 1 only. A custom
quality coder was built and measured first (stage 92: beat SPRING 7/8, beat
Genozip 8/8) before conceding fqzcomp beats it on all three axes (size,
speed, RAM) by 1.1%–4.5% at default — vendored the better one rather than
pretend to have beaten it. Neither Claim 2's caller nor Claim 3's
export/coverage/query ever touches fqzcomp (checked: zero calls in that
code path).

---

## 3. The "India's first" / "world's first" question — handled precisely

Do NOT casually claim "India's first reference-free FASTQ compressor" — it
was checked and the claim does not survive cleanly. **FQC** (Dutta, Haque,
Bose, Reddy, Mande — TCS Innovation Labs, Pune, 2015, *J. Bioinformatics and
Computational Biology*) is a real, published, India-built FASTQ compressor
that one search-derived classification described as reference-free — but
its full text is paywalled everywhere checked (PubMed cookie-wall,
ResearchGate 403, Semantic Scholar empty, WorldScientific 403, no free PDF
on Google Scholar's version cluster), so this was **never actually read**,
and a later search summary directly contradicted the "reference-free"
classification (appears to have been contaminated by a description of an
unrelated 2023 paper). **Status: genuinely unresolved. Do not assert either
way without reading the actual PDF.**

What IS checked and holds, because it is a compound/negative claim rather
than a crowded-field one (only a short, checkable list of tools even
attempt two of the three pieces): **no tool anywhere combines lossless
compression + reference-free variant calling + reference-free locus
retrieval as one mechanism in one archive.** See the table in §2. This is
the safe, defensible "world's first" framing — write "to the best of our
knowledge" per standard academic practice, not as an absolute.

---

## 4. The technical journey — what actually happened, in order

This is not a straight line. Read this to understand *why* the repo looks
the way it does, and to avoid re-deriving conclusions that were already
reached and are now settled.

**Phase 1 — PgRC2 head-to-head reimplementation (Aug 28 – Sep 1).**
Built an independent, from-scratch reimplementation of pseudogenome-based
compression to compare against the real PgRC2 binary (GPL-3, never
vendored, cloned separately at `/root/arcs-clean/method_c`). Numerous
refuted hypotheses along the way (cost models, per-read pricing, parallel
coder probes) — kept in the docs as the evidence they don't work, per this
project's standing rule: **fixes must be formulas over a measured input
property, never a fitted constant or a per-dataset special case.**
S. acidocaldarius flipped from a loss to a win (`3e06957`) by starting the
overlap sweep at L=Lmax instead of Lmax-1, so exact duplicates (which can
only overlap at exactly the read length) stop being invisible to chaining.

**Phase 2 — the incomplete-archive bug, and why it matters (~Sep 2).**
A whole session's PgRC2 comparison numbers (+4.65% aggregate) turned out to
be measured against an archive that could not actually be decoded —
`refc::encode` stored a reference's source but dropped its destination gap,
length, and RC flag as "diagnostics only." The verification script existing
at the time (`verify_lossless.sh`) tested the encoder's *dumped intermediate
streams*, never the *archive* itself — a gap that hid this for the whole
session. Fixed in `a81f55c`; corrected aggregate fell to +1.90%, later
settled at **+1.88%, 6 wins 1 loss vs PgRC2** as the final citable number.
**Lesson, recorded because it recurred:** a verification path that doesn't
read the actual shipped product proves nothing about the actual shipped
product.

**Phase 3 — four silent data-loss bugs found and fixed same day (Sep 2).**
Found by actually decoding archives and diffing against the original FASTQ
— something no prior phase had done. All four keyed on a measured property
of the input, not a per-dataset patch: (1) mismatch positions stored in one
byte, silently clamped and corrupted above 256bp — now varint; (2) mismatches
emitted for "orphaned" reads with no derivable length, desynchronizing the
adaptive coder after them; (3) a contained reverse-strand read indexed with
the wrong length; (4) FSE's constant-input RLE result was silently assumed
to always mean "fill with zero," which was only true by coincidence for the
one stream it was written for — any other constant non-zero stream (e.g.
three all-N reads in an M. tuberculosis sample) silently decoded as zeros.
Cost of the fix: +33B on E. coli. **This is why "byte-identical" and
"lossless" are treated as two different, both-required checks in this
project — one can hold while the other silently fails.**

**Phase 4 — quality wired, full lossless round trip achieved (Sep 2).**
`include/quality_coder.h` wraps vendored fqzcomp/htscodecs (BSD 3-clause).
A complete 4-line FASTQ rebuilt from decoder output alone was verified
byte-identical (same MD5) to the original file for the first time.

**Phase 5 — the 19-dataset full sweep (Sep 9).** One run,
`benchmark_1_run.sh` at a clean git worktree, 5h07m, 19 datasets, 0
failures. This is the run that produced every table in
`benchmark/results/`. `results/FULL_SWEEP_20260909/` holds the original;
`benchmark/results/` is a `cmp`-verified byte-identical copy.

**Phase 6 — Claim 2 mechanism discovered and Claim 3 extended (~Sep 9-10).**
The compression-conceals-the-variant finding, the ablation, the
locus-fragmentation measurement (median 4 places, up to 18.6 Mb apart), and
the 81/400 vs 345/400 addressing-mode result — all as described in §1.

**Phase 7 — the reproducibility audit (Sep 10), the project's most
important self-correction.** A systematic pass tracing every claim
backward: claim -> metric -> raw result -> script -> config -> commit ->
dataset -> environment. Found and fixed, in order of severity:

- **T2.4's published 21/26 had no supporting raw log anywhere on the
  server.** Three independent re-derivations (post-fix binary, exact repeat
  on identical reads, pre-fix binary via `git worktree`) all gave **17/26**
  with identical intermediate counts. Corrected across 25 files. This is
  the one published figure this project withdrew.
- **T3.4's coordinate arm was published as 0/400.** Cause: the sidecar
  index requires an undocumented flag (`CAPS_PILEUP=1`); without it, `query`
  emits the pseudogenome consensus instead of each read's own deviations, so
  a het site can never show disagreement — 0 by construction, not by
  measurement. Fixed and documented; real figure is **81/400**.
- `.gitignore` was silently excluding the audit's own evidence (`*.log`,
  then `*.vcf` — fixed with targeted exceptions).
- The MANIFEST itself was wrong on a first attempt — built from a working
  directory that still had the concurrent agent's `gpt2026_*` files in it;
  regenerated from `git ls-files` only, verified clean via an actual fresh
  clone (this is now the standing method, see §9).
- T2.5's "verification" was initially circular (comparing log to CSV, which
  were written from the same shell variables in the same run — proves
  transcription fidelity, not measurement correctness). Documented as a
  two-tier distinction in `AUDIT.md`: re-executed vs transcription-checked
  only.

**Phase 8 — attribution and framing correction (Sep 10).** The paper had
under- and then over-attributed the pseudogenome idea to PgRC at different
points; settled correctly: Quip (2012) was first for assembly-based read
compression; greedy shortest-common-superstring and q-gram matching are
textbook; PgRC/PgRC2 (2020/2025) is the specific prior art whose approach
this project's engineering is built on top of, and it must be cited as
such — while also stating precisely what's substantially different
(19 independently-coded streams vs PgRC2's DNA-only; first-class
variable-length support vs PgRC2's refuse-or-crash; two regions vs PgRC2's
three pseudogenomes).

**Phase 9 — the freeze and structural audit (Sep 10-11), same session as
this recap.** Repository frozen at a series of tags because the freeze
process itself had defects that had to be found and fixed:

- Five tags in one day all claimed to be "final" — diffed in **commit-date
  order** (not alphabetical — alphabetical order gives the WRONG supersede
  chain, this was caught mid-session). The last one with any actual
  code/result change was `results-final-20260910-corrected`; four more tags
  after it were prose-only.
- `NEW_DATASET_LOCKED.md` was titled "15-dataset" but the published CSV has
  19 rows (15 non-human + 4 GIAB human) — the file never mentioned the
  four human datasets at all (0 hits for "HG00" or "GIAB"). Fixed:
  reconciliation stated up front, human sets given their own section with
  every role claim checked against the actual CSVs (a first draft wrongly
  credited T3.4 to HG002 only — it's all four — caught before commit).
- The same lock file's tail still claimed two datasets (U. gibba, C.
  elegans) were never benchmarked; both are in the final sweep at 3/3
  LOSSLESS each. Corrected in place (struck through, not deleted).
- `RESULTS_INVENTORY.md` said 12 run directories; there are 14.
  `results/claim2/` and `results/phase_a/` were unlisted, and both hold
  files named to look authoritative (`method_b_full_chr20.csv`,
  `t3_t5_full_chr20_HG002.csv`, `allphases_14dataset.csv` — the last is VOID,
  pre-bug-fix numbers). Fixed with explicit rows.
- `sed -i` used to bulk-stamp a freeze-date header across `paper/*.md` twice
  silently converted 7 symlinks into regular files (mode 120000 -> 100644),
  which would have let `paper/` and `docs/` drift apart with no signal.
  Caught by checking `git diff --raw` before committing; restored from the
  prior tag's blobs both times.
- The E. coli re-verification evidence (proving CAPSULE reproduces
  byte-exact: 68,429,027 B, twice) was sitting in an 834 MB gitignored
  working directory — the exact check `paper/PAPER.md` tells a reviewer to
  run first did not actually ship in a clone. Copying it in (40 KB of text)
  surfaced a real finding: **Genozip's default archive size is not
  reproducible across machine states** — same binary, same flags, three
  fresh runs clustered at 119.6 MB against a published 114.6 MB (+4.4%).
  Root cause checked directly, not guessed: Genozip's default VBlock size is
  selected from machine conditions, and setting it explicitly moves the same
  file across a 120.0–86.2 MB range. Confirmed neither Genozip's paper nor
  its documentation states this. **Published Claim 1 numbers were NOT
  changed** — the published Genozip figure happens to be the one that
  favours Genozip, so the −43.26% aggregate margin is understated, not
  inflated. Documented as one bullet in `paper/LIMITATIONS.md` §10 plus the
  raw evidence directory. One commit attempt during this fix broke on shell
  quoting and left a tag pointing at the wrong commit — caught, the tag was
  deleted and correctly re-cut, and the mistake is recorded in the
  corrected commit's own message rather than hidden.
- `CITATION.cff` added for journal-grade citability (Nature Portfolio and
  FAIR-software guidance both name this as required and GitHub does not
  provide it automatically).

---

## 5. Standing rules — do not violate these in future work

1. **Fixes must be formulas over a measured input property, never a fitted
   constant or a per-dataset special case.** (e.g. `COVCAP = 2 x ploidy x H`,
   `COHC = max(2, H/10)` from the sample's own k-mer histogram — not a
   number picked to make one dataset look good.)
2. **Never run two timed benchmark jobs concurrently** — contaminates
   timing and RAM measurements.
3. **A LOSSY result anywhere halts everything until fixed** — do not
   proceed past it, do not paper over it.
4. **Retractions stay in the docs, marked in place, never deleted.** Several
   conclusions in this project's history were wrong; they are struck through
   or annotated where they were made, not scrubbed.
5. **The verification path must read the actual shipped product**, not an
   intermediate dump — this was violated once (§4 Phase 2) and cost a whole
   session's numbers.
6. **Before recommending anything from a tag, doc, or memory: verify it
   against the current file/result, not against what a past session said.**
   This recap itself will go stale — check `benchmark/results/`, `AUDIT.md`,
   and `git tag -l` before trusting any specific number here if significant
   time has passed.
7. **Regenerate `MANIFEST.sha256` from `git ls-files` only**, never from a
   live working directory that might contain the concurrent agent's
   untracked `gpt2026_*` files (see §6).
8. **Diff tag history in commit-date order, never alphabetical** — this
   project already got this wrong once mid-session.
9. **`sed -i` across `paper/*.md` will silently destroy the symlinks to
   `docs/`** — check `git diff --raw` before committing any bulk edit
   there, and restore from the prior tag's blob if it happened.

---

## 6. Known, named weaknesses — do not hide these, do not re-litigate them as new findings

- **HG005 loses both T2.1 (0.834 vs 0.863) and T2.3 (0.593 vs 0.605).**
  Root cause identified by controlled experiment: HG005 is the only
  variable-length read set in the benchmark (250bp trimmed, 216 distinct
  lengths) against 148bp fixed for the other three; truncating it to 148bp
  makes it the *best* individual (F1 0.897) — proving the caller is tuned
  for fixed-length geometry, not that HG005 is inherently harder. The
  truncated figure must NEVER be reported as HG005's actual result.
- **Indel recall has a hard ceiling shared by every de Bruijn caller**,
  DiscoSNP++ included: at missed sites the alternate haplotype's k-mers are
  simply absent (one successor, zero coverage floor), and 66% of misses are
  homopolymer-length changes where REF and ALT probes are the literal same
  string — no bubble exists for any caller. 13 filter mechanisms were tried;
  3 helped. Further filter attempts on this class should be considered
  refuted in advance.
- **Single chromosome (chr20), four individuals.** Defensible against the
  field's actual norm (eBWT2SNP uses two *simulated* chromosomes and one
  real one) but a real, named limit — not argued away.
- **Whole-genome is not attempted, and the blocker is disk, not compute.**
  Projected: ~17h compute for all four individuals, but the caller's k-mer
  spill scales to ~1.2 TB at whole-genome scale, which exceeds any
  reasonably-sized server disk. A superkmer spill format that would cut
  this to ~64 GB is specified in `docs/_removable/SUPERKMER_PLAN.md` and
  **never built** — this is the single concrete engineering prerequisite
  for whole-genome, named rather than left to be discovered.
- **`query` is not a speed win** — `genocat --head=100` is faster
  (0.19s vs 0.46s). The claim is capability, not speed.
- **The FQC "India's first" question is unresolved** (§3) — do not assert
  it either way without reading the actual paywalled paper.
- **At 2.6x coverage** (outside the locked 19-dataset suite), this tool
  loses to SPRING by ~3.97% — an architectural gap (region-to-region
  matching costs ~21 bits per mismatch vs SPRING's read-to-read ~2 bits) that
  was investigated and is real, not fixable by a parameter.

---

## 7. Current frozen state — verify this section against reality before trusting it

As of this writing:

    canonical tag:        v1.0.2-capsule
    commit:                e6cad2044236e5d45a523b8ed7d1112d8ce5fdbc  (then 1a6b846 adds CITATION.cff)
    remotes in sync:       origin/gpt2026, origin/c_star_pg_advance
    manifest:              179/179 verified (benchmark/documentation/MANIFEST.sha256)
    fresh-clone tested:    yes — 2.5 MB clone, 16/16 scope tests pass with zero setup
    code/results changed
      since v1.0-capsule:  0 files (everything since has been documentation)

Run this to re-confirm before doing anything else in a new session:

    cd /root/arcs-clean/c_star_pg_advance   # or wherever this repo now lives
    git log --oneline -1
    git tag -l | sort
    sha256sum -c benchmark/documentation/MANIFEST.sha256
    bash scripts/test_0_scope_and_capability.sh

**What is genuinely NOT done, as of the freeze:**
- **The manuscript itself does not exist.** `paper/` is ~23,000 words of
  source material and dense internal cross-references, not a submittable
  paper. Nature Methods format is ~3,000 words + 3–4 figures. **Zero figures
  exist anywhere in this repository.** This is real, separate work — writing
  and designing figures, not more science. The two mechanism findings
  (Claim 2's compression-conceals-the-variant, Claim 3's locus-fragmentation)
  are the two things that most need a diagram; right now both exist only as
  prose and a table.
- **Zenodo DOI archival is not complete.** `CITATION.cff` is in place;
  Zenodo-GitHub linking requires the repo owner's own login (OAuth), which
  no agent can complete. Once the Zenodo toggle is flipped for this repo, a
  GitHub Release against `v1.0.2-capsule` (or later) will auto-archive and
  mint a DOI — needed because Nature Portfolio explicitly states GitHub
  links alone do not satisfy their code-availability requirement.
- **`docs/` (98 files) has 71 files not referenced from the citable path**
  (`paper/`, `AUDIT.md`, `README.md`, `CLAUDE.md`); 61 of those carry no
  marker indicating they are historical/superseded. None of this affects
  the validity of the published results — checked, no withdrawn figure
  appears unmarked in the citable path — but it means a reader wandering
  into `docs/` directly, rather than through `README.md -> AUDIT.md ->
  paper/`, could land on a stale draft with no signal that it's stale. Low
  priority; noted so it is not rediscovered as if new.

---

## 8. If you are a new session starting work here, do this first

1. Read `AUDIT.md` in full — it is the verification record and supersedes
   any confidence level implied elsewhere.
2. Check `git tag -l` and diff the newest tag against `v1.0.2-capsule` in
   commit-date order before assuming anything has or hasn't changed.
3. Read `paper/PAPER.md` for the condensed claims, then `paper/METHODS.md`,
   `paper/LIMITATIONS.md`, `paper/DISCUSSION.md` for depth.
4. Do NOT re-run the full 19-dataset benchmark to "double check" unless a
   specific, named failure mode is suspected — this project's own standing
   feedback (`feedback_verification_diminishing_returns` in memory) is that
   re-verifying past the point where a failure mode is ruled out wastes a
   session; five independent exact reproductions already happened this
   session for different tables.
5. If asked to write the manuscript: this recap plus `paper/*.md` is the
   complete source material. The gap is compression into ~3,000 words plus
   figures for the two mechanism findings, not new research.
6. If asked to attempt whole-genome scale: read §6's disk-blocker note
   first. The superkmer spill format must be built before attempting it on
   any server smaller than ~1.5 TB free disk.

---

*This file was built at the end of the freeze session, 2026-09-11, from the
repository's own `CLAUDE.md`, `AUDIT.md`, `paper/*.md`, and this session's
full working transcript. It is a recap, not a new source of truth — where it
conflicts with a result file, the result file wins.*
