# Replacing the parameter search with a derivation

## The problem, stated exactly

`encode_adaptive.sh` runs FOUR full encodes over (MAXMAP, MINOV) and keeps the
smallest. That is a brute-force stand-in for a model of the data. It costs 4x
runtime and, worse, it cannot extrapolate: on an unseen dataset it interpolates
between four fixed guesses rather than responding to what the data is.

## What the literature actually does

**GenomeScope 2.0** (Bioinformatics 2017; Nat Commun 2020) fits
`f(x) = SUM_i alpha_i * NB(x, mu = i*lambda, size)` to the k-mer spectrum --
2p negative binomials at lambda, 2lambda, ... Repeat content is the mass at
copy number >= 2; sequencing errors form a separate peak near x = 1; genome
size is G*(k-1)/p. Requires >= 15x per homologue. **Stated failure mode: high
heterozygosity and high repetitiveness are mutually confusable** (Fragaria
iinumae was inferred tetraploid when diploid).

**KmerGenie** (Bioinformatics 2014) hash-samples k-mers -- keep those whose
hash is 0, a fraction epsilon, scale by 1/epsilon -- for an order-of-magnitude
speedup. It models genomic k-mers as ~30 Gaussians over copy number and errors
as Pareto(alpha); 5 free parameters haploid, 8 diploid, fitted by BFGS. Its
selection criterion, *maximise predicted distinct genomic k-mers*, is derived
from the model rather than tuned.

## Why we take KmerGenie's METHOD and not GenomeScope's MODEL

Our governing quantity is `r_second`: the fraction of bases removed by MEM
self-matching in the SECOND region -- i.e. the redundancy among the reads that
failed to map. That is a property of a **biased subset**, not of the global read
set. A global k-mer spectrum estimates global redundancy, which is a proxy at
best; and GenomeScope's own failure mode (heterozygous AND repetitive) is
precisely the regime we would be reaching for with human data.

Fitting a global spectrum to a subset property would be a coincidence, not a
model. What transfers is the *method*: estimate a distribution from a bounded
sample instead of computing it exhaustively.

## The derivation we already have

From PLAN_PG_VOLUME.md, the mapping ceiling is not a parameter at all. It is a
break-even between two costs a read can incur:

    cost_map(m)  = m * bytes_per_mismatch
    cost_append  = readLen * (1 - r_second) * bytes_per_base
    accept iff cost_map(m) < cost_append
    =>  m* = readLen * (1 - r_second) * bytes_per_base / bytes_per_mismatch

Position and strand are paid on BOTH paths, so they cancel and must not appear.
The only unknown is r_second, and the circularity -- r_second is only known
after the stage the ceiling governs -- is what forced the brute-force search.

## Breaking the circularity, KmerGenie-style

1. **Probe.** Map a bounded sample of reads at a permissive ceiling. Append that
   sample's unmapped reads and run the existing second-region MEM pass on the
   sample ONLY. Measure r_second directly. Cost O(sample), not O(n). This is
   KmerGenie's trick: estimate the distribution from a sample and scale.
2. **Derive** m* from the formula, using the measured r_second and the byte
   rates this run has already computed (bytes/base from the main region,
   bytes/mismatch from the mm streams).
3. **Decide** per read against m*, then run the real mapping pass.
4. **Verify** by construction: the decoder is untouched, because the choice only
   moves a read between two representations that both already decode.

This is not an if/else and not a threshold. It is one formula over two measured
rates, and it produces a different ceiling for every dataset because the rates
differ for every dataset.

## Evidence that a real signal exists

Measured optima, from the sweeps already run:

| dataset | MEM-second removal | optimum MAXMAP | as divisor of L |
|---|---|---|---|
| H. salinarum | 93.1% | ~8 | L/19 |
| S. acidocaldarius | 50.3% | ~45 | L/5.6 |

Optima differ 5.6x while read length differs 1.66x, so no constant divisor can
express it -- but the ratio (1 - r_second) is 0.069 vs 0.497, a 7.2x ratio
against an observed 5.6x. **The mechanism tracks.** That is the correlation the
formula predicts, and it is why this is a model and not a fit.

## What must be tested BEFORE building

The formula predicts m* from r_second. Two things have to hold, and neither is
assumed:

1. **Does a bounded sample estimate r_second accurately?** Measure r_second on
   samples of 1%, 5%, 25% and compare against the full-run value on all
   datasets. If sample estimates do not converge, the probe is invalid and this
   plan stops here.
2. **Does the derived m* land near the measured optimum?** Compare m* computed
   from the true r_second against the optimum found by sweeping, per dataset.
   If it does not, the formula is wrong and no amount of sampling saves it.

Only if both hold does the search get replaced. If (1) fails but (2) holds, the
probe is the problem and can be made cheaper differently. If (2) fails, the cost
model is wrong and must be re-derived, not patched.

## MINOV

MINOV is left in the search for now, deliberately. It interacts with MAXMAP --
at MAXMAP=50 on S. acidocaldarius, raising MINOV 16 -> 113 cut mem_triples 35%
-- so its optimum may collapse once the ceiling is derived rather than swept.
Deriving it before that is measured would be guessing. Expected outcome: 4
encodes -> 2, and possibly 1.

## What this is NOT
- Not a per-dataset or per-kingdom switch.
- Not a fitted constant. Every quantity in m* is measured during the run.
- Not GenomeScope. We are not fitting a negative-binomial mixture to a global
  spectrum to infer a subset property.
