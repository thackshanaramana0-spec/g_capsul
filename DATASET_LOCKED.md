# ARCS Benchmark Dataset — FINAL LOCKED 2026-08-25

**DO NOT CHANGE THIS FILE.** All accessions verified via NCBI SRA/DDBJ 2026-08-25.

**2026-08-30 addendum: 7 datasets added, deliberately, per explicit user
decision (not a silent change) — see "Extended Set" section below. The
original 10 above and their banned-accession list are untouched.**

---

## The 10 Full Datasets (Claim 1 primary + Claim 3)

All are complete, unsubsampled SRA runs as deposited by sequencing labs.

| # | Accession | Organism | Kingdom | Read len | _1 size (uncompressed) | Peak RAM | Auto-chunk |
|---|-----------|----------|---------|---------|------------------------|---------|-----------|
| 1 | SRR2584863 | E. coli B REL606 WGS | Bacteria | ~108 bp | ~576 MB | ~5 GB | No |
| 2 | ERR552797 | M. tuberculosis H37Rv WGS | Bacteria | ~301 bp | ~217 MB | ~3 GB | No |
| 3 | SRR554369 | P. aeruginosa PAO1 WGS | Bacteria | 100 bp | ~334 MB | ~4 GB | No |
| 4 | ERR5181310 | SARS-CoV-2 amplicon WGS | Virus | ~150 bp | ~30 MB | ~1 GB | No |
| 5 | ERR17740259 | S. aureus WGS (Firmicutes, 33% GC) | Bacteria | ~148 bp | ~970 MB | ~6 GB | No |
| 6 | SRR065390 | C. elegans N2 WGS | Animalia | 100 bp | ~11 GB | ~18 GB | Suppressed* |
| 7 | DRR976266 | S. cerevisiae WGS | Fungi | ~150 bp | ~1.67 GB | ~8 GB | No |
| 8 | SRR870667 | Theobroma cacao WGS | Plantae | 69 bp SE | ~15 GB | ~28 GB | Suppressed* |
| 9 | SRR36741279 | Leishmania major WGS | Protista | ~75 bp | ~1.70 GB | ~7 GB | No |
| 10 | SRR37283774 | P. falciparum WGS | Protista | ~100 bp | ~669 MB | ~5 GB | No |

**\* Auto-chunk suppressed via `ARCS_AUTOCHUNK_MB=25000` env var in run_block1.sh. Server has 90 GB RAM; both datasets fit single-pass assembly comfortably.**

**Diversity:** 4 bacteria · 1 virus · 1 animal · 1 fungus · 1 plant · 2 protists — 6 kingdoms  
**Note:** S. aureus (slot 5) adds Firmicutes / Gram-positive / 33% GC — distinct from slots 1–3 (all Gram-negative, 51–67% GC).  
**Note:** Slots 6 and 8 were previously banned due to RAM constraints on 16 GB servers. On 90 GB server, peak RAM is ~18 GB and ~28 GB respectively — both safe.

---

## Claim 2 Datasets — HG002-HG005 chr20 at 30× (NOT part of Claim 1)

GIAB chr20 reads are used exclusively for Claim 2 variant calling (T3/T4/T5).
They are NOT included in Claim 1 compression — at ~0.7× chr20 coverage, assembly
is impossible and SPRING would win trivially. Dropping them keeps Claim 1 honest.

Coverage varies across GIAB S3 BAMs (HG002/HG005 ~300×, HG003/HG004 ~60×).
All four are downsampled to 30× chr20 so T3 F1 comparison is apples-to-apples.

| # | File | Individual | Ancestry | Source coverage | Standardized |
|---|------|-----------|---------|----------------|-------------|
| C2-1 | HG002_pooled.fq | NA24385 (HG002) | Ashkenazi son | ~300× | 30× chr20 |
| C2-2 | HG003_pooled.fq | NA24149 (HG003) | Ashkenazi father | ~60× | 30× chr20 |
| C2-3 | HG004_pooled.fq | NA24143 (HG004) | Ashkenazi mother | ~60× | 30× chr20 |
| C2-4 | HG005_pooled.fq | NA24631 (HG005) | Han Chinese son | ~300× | 30× chr20 |

**Method:** stream chr20 reads from GIAB S3 WGS BAMs via samtools, downsample to 30× with seqtk sample --seed 42. Reproducible and exact.

---

## Spot counts and coverage

| # | Accession | Spots | Read len | Genome | Coverage |
|---|-----------|-------|----------|--------|----------|
| 1 | SRR2584863 | ~2.5M | 108 bp | 4.6 Mb | ~59× |
| 2 | ERR552797 | ~720K | 301 bp | 4.4 Mb | ~49× |
| 3 | SRR554369 | ~1.5M | 100 bp | 6.3 Mb | ~24× |
| 4 | ERR5181310 | ~100K | 150 bp | 30 kb | ~500× |
| 5 | ERR17740259 | ~3.0M | ~148 bp | 2.8 Mb | ~158× |
| 6 | SRR065390 | ~33.8M PE | 100 bp | 100 Mb | ~34× (_1 only) |
| 7 | DRR976266 | 5.08M | ~150 bp | 12 Mb | ~63× |
| 8 | SRR870667 | ~74M SE | 69 bp | 430 Mb | ~12× |
| 9 | SRR36741279 | 9.49M | ~75 bp | 32 Mb | ~22× |
| 10 | SRR37283774 | 2.91M | ~100 bp | 23 Mb | ~12.6× |
| C2-1..4 | GIAB HG002-HG005 | ~12.6M each | 150 bp | chr20 63 Mb | 30× (Claim 2 only) |

## Download Command

```bash
bash benchmark/download.sh /data/fastq
```

## ARCS Compression Command

```bash
ARCS_AUTOCHUNK_MB=25000 arcs compress INPUT.fq OUTPUT.arcs
```

No command-line flags. `ARCS_AUTOCHUNK_MB=25000` is a server-side env var (set in run_block1.sh) that raises the auto-chunk threshold to 25 GB, ensuring single-pass assembly on all 10 datasets including SRR065390 (~11 GB) and SRR870667 (~15 GB).

## Disk space estimate (250 GB SSD)

| Files | Size |
|-------|------|
| SRA downloads (compressed) | ~30 GB |
| Uncompressed FASTQ _1 files | ~32 GB |
| ARCS archives + intermediates | ~8 GB |
| SPRING/Genozip archives | ~12 GB |
| **Total peak** | **~82 GB** |

250 GB SSD has ~168 GB headroom.

## Banned accessions — never use again

- SRR390728 — H. sapiens RNA-Seq (NOT WGS)
- SRR1296601 — Glycine max ncRNA-Seq (NOT M. tuberculosis)
- ERR015526 — H. sapiens WGS (NOT E. coli 36bp)
- SRR1294122 — H. sapiens RNA-Seq stem cells (NOT WGS)
- ERR174310 — H. sapiens WGS NA12877 26.6 GB (too large, human excluded from scope)
- SRR988075 — D. melanogaster 8.7 GB _1 (no published gold bpb, replaced)
- SRR16357346 — C. elegans CeNDR 6.4× coverage (too low coverage, replaced by SRR065390)
- SRR1945765 — Arabidopsis thaliana 6.5× coverage (too low coverage, replaced by SRR870667)
- SRR327342 — S. cerevisiae 3.4 GB _1 (replaced by DRR976266 higher coverage)
- SRR1663585 — D. melanogaster 2.33 GB _1 (no published gold bpb, replaced by SRR36741279)

## Extended Set — 7 Additional Datasets, added 2026-08-30

All verified live via NCBI eutils (`esearch`/`efetch` against `db=sra`) on
2026-08-30, not guessed. None overlap the banned-accession list above
(different runs from the same-named banned D. melanogaster/Arabidopsis
entries were considered and Arabidopsis was dropped for coverage reasons —
see rationale column).

Rationale for the group structure: real testing this session established
that win/loss correlates with **coverage → overlap density**, not kingdom
(see `benchmark/reimpl/LAYER_BY_LAYER_ANALYSIS.md` §3b and memory
`reimpl_coverage_regime_factors`). These 7 were chosen to (a) cover a 4th
"novel" domain — Archaea — untested by any FASTQ-compression paper found in
literature search, alongside virus/fungi/protista, (b) fill the
animal/plant imbalance in the overlap-win group, and (c) test whether
P. aeruginosa's real, measured −6.1% loss (isolated to the order/position
layer, cause still unresolved — NOT explained by coverage or kingdom)
replicates on other small, high-coverage bacterial genomes.

| # | Accession | Organism | Kingdom/Domain | Group | Coverage (est.) | Status |
|---|-----------|----------|-----------------|-------|-----------------|--------|
| 11 | SRR32429602 | Human betaherpesvirus 5 (HCMV) | Virus | Novel ×8 | genomic-source shotgun WGS (not amplicon like ERR5181310) | untested |
| 12 | SRR39257532 | Aspergillus fumigatus | Fungi | Novel ×8 | ~20× | untested |
| 13 | SRR29296997 | Halobacterium salinarum | Archaea | Novel ×8 | ~54×, GC~.68 | untested |
| 14 | ERR12954017 | Sulfolobus acidocaldarius | Archaea | Novel ×8 | ~116×, GC~.37 | untested |
| 15 | SRR40104876 | Drosophila melanogaster | Animalia | Overlap-win ×6 | ~19× | untested |
| 16 | SRR40271341 | Helicobacter pylori | Bacteria | Mild-loss ×3 (candidate) | ~139× | **unverified — coverage argues WIN, not loss; see rationale** |
| 17 | SRR40402583 | Campylobacter jejuni | Bacteria | Mild-loss ×3 (candidate) | ~162×, GC~.30 | **unverified — same caveat** |

**Honest flag on #16/#17:** picked by genome-size similarity to P. aeruginosa,
but their coverage (139×, 162×) is deeper than every current winner except
S. aureus (160×) — the established coverage trend predicts these will WIN,
not lose. They stay in "mild-loss (candidate)" only until actually measured
in the locked seq+order scope; do not report them as loss cases in any
paper table until real numbers exist.

**Not yet done:** none of these 7 are downloaded or tested. This section
records the deliberate selection only.

## Published SPRING bpb (PgRC2 2025, for cross-validation)

| Accession | SPRING bpb | Notes |
|-----------|-----------|-------|
| SRR554369 | 0.2416 | P. aeruginosa — use to verify SPRING is working correctly |
| SRR870667 | 1.2621 | T. cacao — SPRING's worst dataset; ARCS expected biggest win here |

Cross-validate: server SPRING bpb must match ±2%. If SRR870667 deviates, SPRING benchmark setup is wrong.
