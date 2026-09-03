# G_CAPSUL — 15-dataset locked set for the SPRING/Genozip comparison

**Supersedes `DATASET_LOCKED.md`'s 17-accession list for this specific
comparison.** `DATASET_LOCKED.md` is unchanged and remains the source of truth
for the original 10 + extended 7. This file records one deliberate swap made
2026-09-02, why, and the resulting exact 15.

Scope: `c_star_pg_advance` only (the G_CAPSUL sequence+order encoder). Not the
outer `/root/arcs-clean` ARCS binary, not names/quality.

---

## The swap

**Removed: Drosophila melanogaster, both accessions.**
- SRR40104876 (~19x, was slot 15 of the 17) — dropped from the main list.
- SRR40104920 (2.6x) — this was a deliberate LOW-COVERAGE control, not part of
  the main list, kept on disk at `/tmp/kd/SRR40104920_1.fastq`. Dropped
  entirely now that Drosophila is not in the study at any coverage.

**Added: Utricularia gibba (Plantae), replacing Drosophila's slot.**
Verified live via NCBI eutils 2026-09-02, not guessed:
- Accession **SRR10676752**
- Genome: **100,688,548 bp** (NCBI assembly U_gibba_v2) — the smallest known
  sequenced plant genome, 4.3x smaller than T. cacao's ~430 Mb, which is why a
  small-but-real Plantae run is possible here where it wasn't for Arabidopsis
  or Drosophila (see `low_coverage_weakness` memory: coverage x genome size
  sets a floor, and searching SRA for "smallest run" in Plantae/Animalia just
  selects for low coverage).
- Coverage (_1 only, same convention as the rest of this table): **~65x**
  (spots=26,658,699, per-mate length ~246bp -> 26,658,699 x 246 / 100,688,548)
- Illumina MiSeq, PAIRED, WGS, GENOMIC source, full unsubsampled SRA run
- BioProject PRJNA595351, study SRP237333, publicly consented
- `_1` file size: ~3.5 GB (half of 6,973 MB total for both mates)
- NOT on disk yet — queued for download

**Second swap, 2026-09-03: Removed Campylobacter jejuni, added C. elegans (Animalia).**
The 15-set as first locked had ZERO Animalia representation (Drosophila was
removed in the first swap above, and C. elegans/T. cacao were excluded from
this repo's set as "the big datasets" per `DATASET_LOCKED.md`). Auditing
kingdom coverage found Bacteria over-represented at 6 of 15 while Animalia
had none. C. jejuni (SRR40402583, an "extended set" addition, not one of
the original primary-10 accessions) was the least structurally necessary of
the 6 Bacteria entries to drop.

- **Removed:** SRR40402583, Campylobacter jejuni, Bacteria, ~162x.
- **Added:** SRR065390, C. elegans N2 WGS, Animalia — one of the ORIGINAL
  primary-10 accessions (`DATASET_LOCKED.md` #6), not a new discovery. Was
  already on disk at `/data/fastq/SRR065390_1.fq` (11,282,985,734 B,
  135,234,184 lines = 33,808,546 reads) from the outer project's own
  earlier work — verified present, not re-downloaded. 100 bp reads, ~34×
  coverage (_1 only).
  **No auto-chunk env var applies here** -- `ARCS_AUTOCHUNK_MB` is the
  OUTER ARCS binary's own mechanism; checked directly against this repo's
  encoder (`stages/106_inprocess.cpp`, `grep -rn AUTOCHUNK`) and confirmed
  it has no auto-chunking behavior at all, so there is nothing to suppress.
  The only real implication of this file's size for THIS encoder is that
  it is simply a large single-pass input like any other -- expect
  proportionally higher RAM/time than the smaller datasets, not a special
  flag. (`DATASET_LOCKED.md`'s ~18 GB peak RAM figure for this accession is
  measured on the OUTER `arcs` binary and is not assumed to transfer
  unchanged to this repo's own encoder -- if this dataset's peak RAM here
  needs stating, it should be measured fresh, not copied.)
- SRR40402583 is NOT deleted from disk or from history — it remains a
  valid, already-measured dataset (`results/phase_a/allphases_14dataset.csv`
  row `SRR40402583`), just no longer one of the current 15 headline
  accessions. Its number stays on record, marked superseded, not erased.

**Kept as-is: Arabidopsis thaliana, low-coverage control.**
- ERR17716639, 2.8x coverage, on disk at `/tmp/kd/ERR17716639.fastq`
- This is the ONE remaining deliberate low-coverage dataset. It stays outside
  the "15" count below — it's a control for the low-coverage weakness
  (see `low_coverage_weakness` memory: we lose to SPRING here, -2.74%,
  disclosed honestly, not hidden), not a member of the main comparison set.

---

## The 15 (main SPRING/Genozip comparison set)

| # | Accession | Organism | Kingdom | Coverage (_1) | Status |
|---|-----------|----------|---------|---------------|--------|
| 1 | SRR2584863 | E. coli B REL606 | Bacteria | ~59x | on disk |
| 2 | ERR552797 | M. tuberculosis H37Rv | Bacteria | ~49x | on disk |
| 3 | SRR554369 | P. aeruginosa PAO1 | Bacteria | ~24x | on disk |
| 4 | ERR5181310 | SARS-CoV-2 | Virus | ~500x | on disk |
| 5 | ERR17740259 | S. aureus | Bacteria | ~158x | on disk |
| 6 | DRR976266 | S. cerevisiae | Fungi | ~63x | on disk |
| 7 | SRR36741279 | Leishmania major | Protista | ~22x | on disk |
| 8 | SRR37283774 | P. falciparum | Protista | ~12.6x | on disk |
| 9 | SRR32429602 | Human betaherpesvirus 5 (HCMV) | Virus | n/a | **queued for download** |
| 10 | SRR39257532 | Aspergillus fumigatus | Fungi | ~20x | on disk |
| 11 | SRR29296997 | Halobacterium salinarum | Archaea | ~54x | on disk |
| 12 | ERR12954017 | Sulfolobus acidocaldarius | Archaea | ~116x | on disk |
| 13 | SRR40271341 | Helicobacter pylori | Bacteria | ~139x | on disk |
| 14 | SRR065390 | C. elegans N2 WGS | Animalia | ~34x | on disk |
| 15 | SRR10676752 | Utricularia gibba | Plantae | ~65x | on disk |

**Kingdom split:** Bacteria 5, Virus 2, Fungi 2, Protista 2, Archaea 2,
Plantae 1, Animalia 1 — **all 6 standard kingdoms now represented**, one
gap (Animalia) closed 2026-09-03 by the C. jejuni -> C. elegans swap above.

**Low-coverage supplementary control (outside the 15):**

| Accession | Organism | Coverage | Status |
|---|---|---|---|
| ERR17716639 | Arabidopsis thaliana | 2.8x | on disk, `/tmp/kd/ERR17716639.fastq` |

---

## Downloads needed

**None. All 15 confirmed on disk at `/data/fastq/<accession>_1.fq`,
2026-09-03** (verified by direct file check, not assumed): the two
originally queued (SRR32429602 HCMV, SRR10676752 Utricularia gibba) were
completed and relocated from `/tmp/newdl/` — a non-persistent location,
moved into place specifically so `scripts/run_claim1_bench.sh` (which reads
from `$CAPSULE_DATA_DIR`, i.e. `/data/fastq/`) can find them — and
SRR065390 (C. elegans, added in the second swap above) was already present
from earlier work.

**Outstanding, separate from "downloaded":** SRR10676752 (Utricularia
gibba) has had three real G_CAPSUL-side encode attempts (sequence-only,
+names, +names+quality — the last completing at 1,655,447,728 B) but no
saved archive, no lossless verification, and no SPRING/Genozip comparison
yet. It is on disk and ready to run, not yet benchmarked. SRR065390 (C.
elegans) has never been run through this repo's own comparison at all.
