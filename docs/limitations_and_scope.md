# Limitations and scope

> Consolidated list of limitations, scope boundaries and withdrawn numbers across the three claims.

## Withdrawn numbers (do not cite)

| Withdrawn number | Claim | Correct replacement | Why withdrawn |
|---|---|---|---|
| 80.8% genome fraction, 90% unaligned | T3.1 | 98.698% genome fraction | Measured against a broken export (contigs concatenated into 2 fake records) |
| 345/400 locus retrieval | T3.4/T3.5 | 1,452/1,452 (4 individuals, 1 locus) | Measured with a pre-strand-fix, single-individual query |
| 0/400 coordinate-only | T3.4/T3.5 | 81/400 (old), or the 96/400-style internal control cited in the current mechanism file | The 0 came from the query emitting the consensus, not the reads |
| 21/26 multi-allelic | T2.4 | 17/26 | No supporting raw log; three independent re-derivations all gave 17/26 |
| 5/111, 11/18 multi-allelic | T2.4 | 17/26 | 5/111 used a check that cannot evaluate indel-bearing sites; 11/18 unreproducible |
| 0.637 / 0.666 het-indel F1 | T2.3 | 0.621 | Measured on 8 chr20 windows, not the full chromosome |
| 0.890 het-SNV F1 | T2.1 | 0.876 | Pre-19-dataset-sweep figure |
| 0.836 tetraploid SNV F1 | T2.5 | 0.897 | Pre-final-sweep figure |
| 555-656x export speedup | T3.1 | 129-784x | Pre-final-sweep figure |
| "plain query alone reaches 1.0000, no completion index needed" | T3.4 | "an index over the reads that deviate from the pseudogenome is required: archive alone 0.6047-0.8536, with the FM-index 1.0000" | The earlier figure held only on HG002 |

## Disclosed trade-offs

- **T1.2 (wall time): G_CAPSUL is fastest on 0 of 19 datasets, for both compression and
  decompression, against SPRING and Genozip.** Confirmed by recomputing the raw CSV and by
  the manuscript's own bold marking. This is the cost of the structure that Claims 2 and 3
  reuse.
- **FM-index size: 1.9-3.5x the archive itself** (0.9-2.6 MB per individual), growing with the
  share of reads that deviate from the pseudogenome. It is optional, built after the archive, and
  adds no cost to Claim 1, but it does take disk space. The k-mer completion index it replaces was
  7-13x.
- **Locus-retrieval false-positive cost:** homozygous-control hits rise from 29/35/29/66 to
  48/52/49/87 of 400 when the FM-index is added to the position query. Completeness is bought with
  specificity.
- **T3.1 correctness: G_CAPSUL wins 2 of about 7 QUAST metrics (genome fraction, indels) and loses
  the rest (mismatch rate, duplication, N50, misassembly count) to SPAdes.** This is structural.
  Closing it needs repeat resolution and error correction, which is an algorithm project and not a
  configuration change. Four attempts (region split, full parameter sweep, contig dedup, consensus
  polish) were tried and rejected, each for a documented reason.
- **T2.1/T2.3: DiscoSNP++ wins on HG005**, the one exception across het-SNV and het-indel calling.
- **PgRC2 comparison: about 1.7x slower and up to 2.5x heavier**, on the narrower sequence-only
  scope where the comparison is valid.

## Scope limits

- **T3.4/T3.5 ground truth exists only for the 4 GIAB human individuals** (HG002-005). GIAB is
  human-only by the scope of the NIST program, and no validated heterozygous-site truth set exists
  elsewhere.
- **T3.4/T3.5 were tested on chr20 windows, not the whole genome.** Chr20-subset validation and
  restriction to GIAB confident regions are standard conventions.
- **T3.1 and T3.2 use 6 of the 19 locked datasets. T3.3 uses all 19.** T3.1 and T3.2 need an
  external tool (SPAdes, bwa) to also succeed on the same dataset. T3.3 needs nothing external.
  This is disclosed in commit `c73753f`.
- **T2.5 (tetraploid) is a 400 kb region, not a full chromosome**, as the manuscript's table
  caption states. Do not imply chromosome scale for this table.
- **Human whole-genome scale is out of scope.** Full-WGS datasets were excluded as too large. All
  human-scale results are chr20-region or 19-dataset-sweep scale, not full 3.1 Gb genomes.
- **No speed or index-size benchmark against BEETL, CIndex or sFASTQ exists.** This is the most
  valuable missing experiment.
- **At low coverage (2.6x, outside the locked 19-dataset suite) G_CAPSUL loses to SPRING by about
  3.97%.** The cause is architectural and not a tunable parameter. At low coverage most reads fail
  to chain and are handled by region-to-region MEM matching, which costs about 21 bits per
  mismatch against about 2 bits for SPRING's read-to-read matching. A parameter sweep cannot close
  an order-of-magnitude per-mismatch difference that comes from which matching regime activates.
  This case is not part of the published 19-dataset comparison. It is noted here for anyone who
  extends the suite to lower-coverage input.
- **`query` is O(archive), not O(range).** There is no windowed or local decode, and the reason is
  structural. The matcher of a self-referential compressor optimises match length and ignores match
  distance. On a measured E. coli archive, 44% of references point to a copy about 15 MB further
  away than an identical nearer copy that exists but was outside the matcher's sampled candidate
  set (index sampling step 15 positions, so a nearer copy is reachable about 1 time in 15). This
  distance makes the transitive closure of any local window large (a 100 kb window's closure is
  179,322 B, 0.67% of the pseudogenome). It is also why every member of this compressor family
  (PgRC, NanoSpring, SPRING) offers only full compression and decompression. A proposed fix, a
  nearest-among-equal-length-candidates tie-break, was implemented and rejected: ties are only
  0.035% of candidates, so the tie-break had nothing to act on. The remaining lever is index
  sampling density, which trades directly against the mapping stage's measured bottleneck
  (candidate-walk cost) and is an open question. The decoded-index cache (17x faster query)
  avoids this cost and does not remove it.
- **DiscoSNP++'s multi-allelic limitation is inferred from its documented behaviour, not verified in
  its source.** See `novelty_and_prior_art.md`.

## Statistical methodology

No confidence intervals or significance tests are reported. Comparable work in this area (FASTQ
compression ratio comparisons, and reference-free variant-calling F1 comparisons against
DiscoSNP++ and Kmer2SNP-class tools) does not report them for these comparison types. Where
sensitivity to the input matters, it is addressed with sweeps: the MAXMAP/MINOV sweep of Claim 1,
the T2.2 coverage sweep and the `dup_frac` plateau of Claim 2, and the per-dataset speedup range
of Claim 3. The evidence offered is that the result holds across a swept range of inputs.

## Why SARS-CoV-2 / ERR5181310 is excluded from the T3.1 comparison

The diagnostic numbers below come from an assembly check on ERR5181310 that was never finalised.
They are cited for those numbers only and not as a completed investigation. On this dataset
(amplicon sequencing, more than 1000x PCR-duplicated coverage of a 30 kb genome) the exported
pseudogenome is about 800x oversized relative to the true genome (24,072,995 bp exported against
a 29,903 bp reference), although exact-duplicate collapse already removes 92% of the reads
(1,825,142 to 152,085 unique).

**Working hypothesis, not confirmed by tracing the code path.** Greedy overlap chaining assumes
that overlap length and quality reliably signal that a read extends the genome. At this ratio of
depth to genome size, many unrelated near-duplicate reads overlap the growing chain's tail well
enough by chance that the chain never terminates. It strings nearly the whole pool of 152,085
reads into one path (152,085 x about 150 bp is about 24,000,000, matching the observed `PG_LEN`
almost exactly). For this reason the dataset is excluded from the T3.1 comparison and not scored
against SPAdes. **The T3.1 speed claim (129-784x) is unaffected and independently measured.** The
limitation concerns only whether the correctness of the export could be claimed equal to SPAdes on
this data shape. It is not claimed.

## Code-level items

- **A strand bug exists in the `.sites` pileup tally**, inconsistent with the strand correction
  in `query`. A fix was written, built and gated, then reverted, because adopting it would change
  the published native-pileup numbers (400/400, 25 FP, 19.6x faster) and that re-gate was not done.
  The fix is documented as a comment at the tally site.
- **Consensus-polishing an export from the archive's own pileup makes correctness worse.** The
  negative result is kept in the code as a comment so that it is not repeated.
- **A retuning of the encoder (raising `MAXMAP`) was abandoned.** It would have cost about 11.6% of
  Claim 1's compression ratio to improve a Claim 3 metric that the paper does not claim to
  optimise.
- **DiscoSNP++'s VCF output has an off-by-one**, found by reading its source. The final step of its
  pipeline (the zero-to-one experiment script) adds `+1` to every `POS` value, and its final
  console message still describes the output as "0-based". Scoring DiscoSNP++'s output without
  correcting for this reduces its measured accuracy to essentially zero (0.004 in an earlier
  version of the results). All Claim 2 comparisons against DiscoSNP++ are scored with the
  correction applied.
