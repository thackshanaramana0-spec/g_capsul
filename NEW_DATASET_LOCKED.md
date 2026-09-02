# CAPSULE — 15-dataset locked set for the SPRING/Genozip comparison

**Supersedes `DATASET_LOCKED.md`'s 17-accession list for this specific
comparison.** `DATASET_LOCKED.md` is unchanged and remains the source of truth
for the original 10 + extended 7. This file records one deliberate swap made
2026-09-02, why, and the resulting exact 15.

Scope: `c_star_pg_advance` only (the CAPSULE sequence+order encoder). Not the
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
| 14 | SRR40402583 | Campylobacter jejuni | Bacteria | ~162x | on disk |
| 15 | SRR10676752 | Utricularia gibba | Plantae | ~65x | **queued for download** |

**Kingdom split:** Bacteria 6, Virus 2, Fungi 2, Protista 2, Archaea 2,
Plantae 1. Animalia now has zero representation in the main 15 (Drosophila
removed, C. elegans and T. cacao excluded as the "big" datasets per
`DATASET_LOCKED.md`) — disclose this honestly, do not backfill with another
low-coverage Animalia run (see `low_coverage_weakness` memory for why that
would be confounded, not cheap).

**Low-coverage supplementary control (outside the 15):**

| Accession | Organism | Coverage | Status |
|---|---|---|---|
| ERR17716639 | Arabidopsis thaliana | 2.8x | on disk, `/tmp/kd/ERR17716639.fastq` |

---

## Downloads needed

Two accessions, both queued 2026-09-02:
- SRR32429602 (HCMV)
- SRR10676752 (Utricularia gibba)
