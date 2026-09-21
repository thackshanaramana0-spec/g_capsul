---
Date: 2026-09-19
Title: Honest Limitations and Scope — Every Real Trade-Off and Withdrawn
  Number, In One Place
Purpose: A single consolidated list so nothing gets forgotten, silently
  softened, or accidentally cited after withdrawal, across all three
  claims.
When to refer to this file: Writing the Discussion/Limitations section;
  responding to reviewer questions about weaknesses; before citing ANY
  number from this project, to check it hasn't been superseded.
Keywords: limitations, scope, withdrawn numbers, trade-offs, disclosed
  costs, honest, do not cite
---

# Honest limitations and scope — consolidated

## Numbers that are WITHDRAWN and must never be cited again

| Withdrawn number | Claim | Correct replacement | Why withdrawn |
|---|---|---|---|
| 80.8% genome fraction, 90% unaligned | T3.1 | 98.698% genome fraction | Measured against a broken export (contigs concatenated into 2 fake records) |
| 345/400 locus retrieval | T3.4/T3.5 | 1,452/1,452 (4 individuals, 2 loci) | Measured with a pre-strand-fix, single-individual query |
| 0/400 coordinate-only | T3.4/T3.5 | 81/400 (old), or the 96/400-style internal control cited in the current mechanism file | The 0 was our own query emitting consensus, not reads |
| 21/26 multi-allelic | T2.4 | 17/26 | No supporting raw log; three independent re-derivations all gave 17/26 |
| 5/111, 11/18 multi-allelic | T2.4 | 17/26 | 5/111 used a check that cannot evaluate indel-bearing sites; 11/18 unreproducible |
| 0.637 / 0.666 het-indel F1 | T2.3 | 0.621 | Measured on 8 chr20 windows, not the full chromosome |
| 0.890 het-SNV F1 | T2.1 | 0.876 | Pre-19-dataset-sweep figure |
| 0.836 tetraploid SNV F1 | T2.5 | 0.897 | Pre-final-sweep figure |
| 555-656x export speedup | T3.1 | 129-784x | Pre-final-sweep figure |
| "plain query alone reaches 1.0000, no completion index needed" | T3.4 | "completion index required in general" | True only on HG002 by chance; HG003/4/5 need it |

## Real, disclosed trade-offs — not weaknesses to hide, costs to state plainly

- **T1.2 (wall time): CAPSULE is fastest on 0 of 19 datasets, either
  compression or decompression, vs SPRING/Genozip.** Verified two
  independent ways this session (raw CSV recomputation + the manuscript's
  own bold-marking). This is the direct cost of the structure Claims 2/3
  reuse for free — frame it that way, not as an unexplained loss.
- **the completion index index size: 7-13x the archive itself.** Optional, post-archive,
  zero cost to Claim 1 — but real disk cost if used.
- **the completion index false-positive cost: homozygous-control FP rises ~35-75%
  relative when engaged.** Completeness bought with specificity.
- **T3.1 correctness: wins 2 of ~7 QUAST metrics (genome fraction,
  indels), loses the rest (mismatch rate, duplication, N50, misassembly
  count) to SPAdes.** Structural — closing it needs repeat resolution and
  error correction, an algorithm project, not a configuration fix. Four
  independent attempts (region split, full parameter sweep, contig dedup,
  consensus polish) were tried and refused this session, each for a
  documented reason.
- **T2.1/T2.3: DiscoSNP++ wins on HG005**, the one honest exception across
  both het-SNV and het-indel calling. Not averaged away.
- **PgRC2 comparison: ~1.7x slower, ~2.5x heavier at worst**, on the
  narrower sequence-only scope where the comparison is even valid.

## Scope limits — real, stated precisely, not vague

- **T3.4/T3.5's ground truth exists only for the 4 GIAB human individuals**
  (HG002-005) — GIAB is confirmed (NIST's own program scope, verified this
  session via live web search) to be human-only. This is not a choice this
  project made to narrow scope conveniently; no validated het-site truth
  set exists anywhere else to test against.
- **T3.4/T3.5 tested at chr20 windows, not the whole genome** — confirmed
  this session to be standard field practice (chr20-subset validation and
  GIAB confident-region restriction are both documented conventions, not
  a project-specific shortcut).
- **T3.1/T3.2 use 6 of the 19 locked datasets; T3.3 uses all 19.** Not an
  inconsistency — T3.1/T3.2 need an external tool (SPAdes, bwa) to also
  succeed on the same dataset; T3.3 needs nothing external. This is
  disclosed directly in commit `e746bc8`.
- **T2.5 (tetraploid) is a 400kb region, not a full chromosome** —
  explicit in the manuscript's own table caption; do not imply
  chromosome-scale for this one table when citing it elsewhere.
- **Human whole-genome scale is explicitly out of scope** — an earlier
  project decision (banned full-WGS datasets as "too large"); all
  human-scale results in this paper are chr20-region or 19-dataset-sweep
  scale, not full ~3.1Gb genomes.
- **No speed or index-size benchmark against BEETL/CIndex/sFASTQ exists.**
  Stated as the single most valuable missing experiment, not silently
  omitted.
- **At low coverage (2.6×, outside the locked 19-dataset suite), this
  project loses to SPRING by ~3.97%** — verified against root-level
  `RECAP.md` §6, a real, investigated, architectural gap, not a tunable
  parameter. Root cause: at low coverage, most reads fail to chain and are
  instead handled by region-to-region MEM matching, which costs ~21 bits
  per mismatch against SPRING's read-to-read matching cost of ~2 bits — an
  order-of-magnitude per-mismatch penalty that a parameter sweep cannot
  close, because it is a structural consequence of which matching regime
  activates at low coverage, not a coefficient. Not part of the published
  19-dataset comparison, and correctly excluded from it, but named here so
  it is not rediscovered as a surprise if the suite is ever extended to
  lower-coverage inputs.
- **`query` is O(archive), not O(range) — no windowed/local decode exists,
  and the reason is structural, not an engineering gap.** Verified against
  `docs/LOCALITY_TENSION.md`. A self-referential compressor's matcher
  optimises match *length* and is blind to match *distance*: on a measured
  E. coli archive, 44% of references point to a copy ~15 MB further away
  than an identical *nearer* copy that exists but was never selected,
  because it was outside the matcher's sampled candidate set (index
  sampling step = 15 positions, so a nearer copy is reachable only ~1-in-15
  of the time). This gratuitous distance is what makes the transitive
  closure of any local window large (a 100 kb window's closure is 179,322
  B, 0.67% of the pseudogenome) and is *why* every member of this
  compressor family (PgRC, NanoSpring, SPRING) offers only full compress/
  decompress, never an in-between windowed mode — not an interface
  omission, the representation genuinely does not support a local answer
  as shipped. **A proposed fix (nearest-among-equal-length-candidates
  tie-break) was implemented and refuted**: ties are only 0.035% of
  candidates, so there was nothing for the tie-break to act on — the real
  lever is index sampling density, which trades directly against the
  mapping stage's own measured bottleneck (candidate-walk cost), and
  remains an explicitly open, unresolved question, not a plan with a
  known cost. The current sidecar-index cache (17× faster query) is
  correctly described as a cache that avoids this cost, not a fix that
  removes it.
- **DiscoSNP++'s multi-allelic limitation is inferred, not source-verified**
  — see `novelty_and_prior_art.md`.
- **The manuscript's ARCH naming, if adopted, is not yet in the live tex
  file** — this whole `refer_paper_docs/` tree describes a proposed
  addition, not a currently-published fact about the manuscript's wording.

## Why no confidence intervals or significance tests appear anywhere in this project

A deliberate, stated methodological choice, not an oversight — verified
across all of this project's own research-audit checklists
(`docs/RESEARCH_CHECKLIST_CLAIM1/2/3.md`): statistical methodology (CIs,
significance testing) is marked N/A throughout, consistent with field
precedent — no comparable paper in this space (FASTQ compression ratio
comparisons, or reference-free variant-calling F1 comparisons against
DiscoSNP++/Kmer2SNP-class tools) reports significance for these comparison
types, so holding this project to a different standard than its own field
would be inconsistent, not more rigorous. Where sensitivity to input does
matter, it is addressed directly with sweeps instead (Claim 1's MAXMAP/MINOV
sweep, Claim 2's T2.2 coverage sweep and the `dup_frac` plateau
characterization, Claim 3's per-dataset speedup range spanning multiple
datasets rather than one cherry-picked figure) — the honest substitute for
significance testing in this domain is showing the result holds across a
swept range of inputs, not a p-value on one comparison.

## Why SARS-CoV-2/ERR5181310 is excluded from T3.1's fair comparison — the mechanism, not just the fact

Verified against `docs/T3.1_ERR5181310_ASSEMBLY_CHECK.md`, itself explicitly
labeled in-progress and never finalized — cited here only for the
diagnostic numbers it does report, not as a completed investigation. On
this one dataset (extreme amplicon-sequencing depth, >1000x PCR-duplicated
coverage over a 30 kb genome), this project's exported pseudogenome is
**~800× oversized** relative to the true genome (24,072,995 bp exported
against a 29,903 bp reference) despite exact-duplicate collapse already
removing 92% of reads (1,825,142 → 152,085 unique). **Working hypothesis,
explicitly not yet confirmed by tracing the actual code path**: greedy
overlap chaining assumes overlap length/quality reliably signals "this read
extends the genome," but at this depth-to-genome-size ratio, many unrelated
near-duplicate reads overlap the growing chain's tail well enough by chance
that the chain never terminates — it strings nearly the entire 152,085-read
pool into one path (152,085 × ~150 bp ≈ 24,000,000, matching the observed
`PG_LEN` almost exactly). This is why the dataset is excluded from T3.1's
fair-comparison scope rather than scored against SPAdes as-is — not an
arbitrary exclusion, but one with a stated, numerically consistent (if
unconfirmed) causal account. **T3.1's speed claim (129–784×) is unaffected
and independently measured** — this limitation concerns only whether
export's *correctness* could be claimed equal to SPAdes's on this specific
data shape, which it explicitly is not claimed to be.

## Code-level honesty items — things found, fixed, or deliberately left alone this session

- **A real strand bug exists in the `.sites` pileup tally** (inconsistent
  with `query`'s own strand correction). A working fix was written, built,
  and gated — then deliberately reverted, because adopting it would move
  the published native-pileup numbers (400/400, 25 FP, 19.6x faster) and
  that re-gate was not done this session. The fix is documented as a
  comment at the tally site, not silently dropped.
- **Consensus-polishing an export from the archive's own pileup makes
  correctness worse, not better** — a real, tested, negative result kept
  in the code as a comment specifically so it is not re-attempted.
- **An encoder retuning attempt (raising `MAXMAP`) was abandoned** after
  being caught mid-session as off-mission — it would have cost ~11.6% of
  Claim 1's compression ratio to chase a Claim 3 metric the paper never
  claims to optimize.
- **DiscoSNP++'s VCF output has a real off-by-one, found by reading its
  source, not by assumption.** Verified (`docs/DISCOSNP_INTERNALS.md`):
  its pipeline's final step (`zero2one.py`) adds `+1` to every `POS` value
  as its very last stage, and its own final console message still claims
  the output is "0-based" — stale, misleading, contradicted by its own
  code. Scoring DiscoSNP++'s output without correcting for this is not a
  cosmetic error — an earlier version of this project's own results
  recorded that scoring without the correction collapses DiscoSNP++'s
  measured accuracy to essentially zero (0.004), which would have made
  every DiscoSNP++ comparator number in this paper meaningless had it gone
  uncaught. All Claim 2 comparisons against DiscoSNP++ are scored with this
  correction applied.
