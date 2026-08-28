# Independent reimplementation of PgRC2's assembler

Written from the algorithm, not copied. PgRC2's source was read closely to
understand what each stage does — and reading it is what found several of the
results below — but every line here is our own. `method_c/` is the unmodified
upstream clone used as the measurement target; it is GPL-3 and is never
vendored (see `.gitignore` and `../PGRC2_COMPARISON.md`).

Purpose: settle whether PgRC2's speed, memory and assembly quality are
reachable from textbook methods alone, and extract what transfers to ARCS.
These are throwaway measurement instruments, kept for the numbers in them.

Target (yeast_sub.fq, 1M reads, 851,275 unique after N-filter and dedup):
PgRC2 = 22.60M pseudogenome, 4.0 s, 232 MB.

## Results

| # | file | what it adds | pg | notes |
|---|------|--------------|-----|-------|
| 01 | greedy_exact | textbook greedy SCS, cycle avoidance | 49.84M | 40 s; matches their filter+dedup to the read (851,275) |
| 02 | minimizer_tolerant | minimizer-seeded tolerant overlap | 38.88M | tolerance replaces the quality split |
| 03 | pigeonhole_mapping | map leftovers into the pg (q-gram lemma) | 36.92M | + in-place RC rescan |
| 04 | rc_in_greedy | reverse complement during growth | 38.54M | **worse** — refuted |
| 05 | single_pass_overlap | packed 32-base seeds, rolling, outward walk | 34.46M | **40 s → 6 s**; the O(n·L²) fix |
| 06 | multi_candidate | keep backup partners | 35.46M | **worse** — weak merges cost bytes |
| 08 | true_sweep | full descending-length sweep | 36.34M | **worse** — same reason |
| 09 | two_round_division | round 1 labels, round 2 builds | 34.16M | chains 29.90M → 25.13M |
| 10 | scaled_seed_width | part width = readLen/(mm+1) | 31.59M | fixed 32 capped tolerance at 3 |
| 11 | second_pseudogenome | assemble survivors; admit[] bugfix | 29.06M | main pg **20.15M < their 21.10M** |
| 12 | classifier_sweep | round 1 links at widest range | 28.91M | best config |
| 14 | relaxed_division | admit non-singletons (86.7%) | 34.56M | **worse** — strictness is the point |
| 15 | parameterised_seed | sweep seed width | **27.49M** | optimum at 16; 12 and 10 regress |

Best: **27.49M, 11 s, 421 MB** against their 22.60M, 4.0 s, 232 MB.

## What this settled

Their speed and memory ARE reachable from textbook parts — same class on both,
with none of their code. Assembly quality was not fully reached (21.6% short),
though the main pseudogenome alone beats theirs; the deficit is that ours holds
141k fewer reads, so more survivors fail to place.

Every gain above came from fixing a defect in this code, not from adding
something PgRC2 has and we lacked:
  - a fixed 32-base seed silently ceilings mismatch tolerance at 3
  - a missing admit[] check emitted excluded reads twice (5 MB)
  - the survivor pseudogenome stage was absent entirely
  - the round-1 classifier was needlessly conservative
  - O(n·L²) re-hashing, not the algorithm, was the 40 s
  - MINOV had been swept in a different architecture and the conclusion carried over

## What transferred to ARCS

  - packed 32-base seeds as integer keys (hash_apsp.cpp)
  - index one side, stream the other
  - outward offset walk: first verified hit is the longest overlap
  - exact-size CSR table, which also fixes silent k-mer loss on bucket overflow

## Findings worth keeping

  - PgRC2's quality-based HQ/LQ split is INACTIVE by default; the division that
    does the work is structural, computed from round-1 chain topology
  - "weak merges cost bytes" is CONDITIONAL: true before the division, inverted
    after it. Once reads pass a tiling filter a marginal merge converts a
    150-byte survivor into a ~30-byte chain member
  - their read-to-pg matching is the q-gram/pigeonhole lemma, textbook
  - reverse complement in the greedy is refuted twice, in two architectures

Build any of these standalone:
    g++ -O3 -march=native -o scs <file>.cpp
    ./scs reads.fq [maxmm] [minov] [seedwidth]
