# Stage 78 (negative result) — suffix-array-based MEM matching, tried and rejected

User asked to actually try the suffix-array alternative for `pg MEM matching`
(32% of assembly's time, the single biggest phase) after stage 77's component
analysis flagged it as the one real algorithmic alternative in the literature
(what MUMmer/E-MEM use) but too risky to attempt without explicit direction.
This documents the real attempt and why it was rejected — not theorized,
built and measured on the real pseudogenome.

## What was built

Cloned `IlyaGrebnov/libsais` (Apache 2.0, actively maintained, SIMD-optimized
SA-IS construction, includes OpenMP-parallel variants) rather than
hand-rolling SA-IS from scratch — the same "read/use a proven real
implementation" discipline applied to fqzcomp/SPRING/Genozip all session,
since suffix array construction is notoriously easy to get subtly wrong.

Real pseudogenome dumped from a real run (`DUMP_PG=1`): 23,233,953 bases.

Built SA + PLCP + LCP via `libsais_omp`/`libsais_plcp`/`libsais_lcp_omp`
(12 threads), then a real (not theoretical) MEM finder: a greedy
left-to-right self-match parse using a capped bidirectional SA-rank-neighbor
scan with running-min LCP — a deliberately SIMPLER, lower-correctness-risk
variant of the textbook linear-time LZ-factorization-via-suffix-array
algorithm (which needs PSV/NSV precomputation plus a union-find structure —
avoided here specifically because it is easy to get subtly wrong and this
was meant to be a fast go/no-go check, not the final engineering).

## Real, measured result — decisive negative on two independent axes

| | current implementation (seed-index) | SA-based (this attempt) |
|---|---|---|
| Coverage | forward + reverse-complement | forward only (RC never added — already losing before this cost is even included) |
| Time | 0.60-1.05s (both directions, already parallel) | 0.90s (forward only; SA+LCP construction 0.43-0.46s + parse 0.44-0.45s) |
| Bases removed, main pg | 10,209,946 (47.2%) | 2,368,876 (10.2%) |

**Speed**: SA+LCP construction alone, even parallelized across 12 threads,
already costs about as much as the current implementation's ENTIRE phase
(both directions). Adding reverse-complement coverage (needed for parity)
would roughly double the input processed, pushing total time well past the
current implementation's.

**Match quality**: swept the neighbor-scan cap K from 64 to 512 (an 8x
widening) and match count moved from 41,537 to 41,541 — essentially no
change, proving the cap was not the limiting factor. The simplified scan has
a real structural gap: repeat regions with many copies can have their
closest EARLIER occurrence sitting far outside any local SA-rank window, so
a capped local scan misses most of the real matches the current seed-index
approach finds. Closing that gap needs the full textbook algorithm
(PSV/NSV + union-find), which would add more computation, making the
already-losing speed number worse, not better.

## Decision: not integrated

Both axes are real, measured, and unfavorable — not "needs more tuning."
Fixing the match-quality gap requires strictly more engineering and
strictly more time than what already loses on speed. The current
seed-index-based matcher (`pg MEM matching` / `pigeonhole mapping` in
`45_sweep_loop.cpp`, confirmed already parallel via
`std::thread::hardware_concurrency()`) stays as the right choice for both
phases. No code in the real pipeline changed as a result of this
investigation — this file exists so the negative result is on record and
this specific approach is not silently retried later.

Test harness (not part of the pipeline, kept as a record only):
`/tmp/sa_test/sa_mem.cpp`, `/tmp/sa_test/satime_omp.cpp`, linked against a
locally built `libsais` (cloned from
`github.com/IlyaGrebnov/libsais`, not vendored into this repo).
