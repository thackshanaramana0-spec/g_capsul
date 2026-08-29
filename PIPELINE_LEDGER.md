# Pipeline component ledger — status of every layer

Built 2026-08-29 to work through systematically, not skip around. Each row:
current status, whether it has been diagnosed for headroom, whether a web
search has been done specifically for it.

| # | component | current | headroom checked? | web search done? |
|---|---|---|---|---|
| 1 | load/filter/dedup | 0.60-0.68s | yes (stage 43-45) | no — not needed, mechanical |
| 2 | division round1/round2 | optimal (R1MINOV swept) | yes | no |
| 3 | chain emission | 0.22s | not deeply profiled | no |
| 4 | pigeonhole mapping | 0.59-0.65s | PARTIAL — only compared to PgRC2's own bug | NOT YET |
| 5 | survivor/second pg | measured | not deeply profiled | no |
| 6 | MEM matching, single pass | matched to greedy optimum via cycle-cover proof | yes (matching optimality) | yes (SCS approx literature) |
| 7 | MEM matching, ITERATIVE | BUILT & VERIFIED: -42,444 B, converges 2 passes | DONE (bug found+fixed) | n/a, own diagnostic |
| 8 | order/permutation coder | at log2(n!)+0.022% floor | yes | NOT YET (paired/structured order) |
| 9 | position/strand stream | near Poisson floor | partial | NOT YET |
| 10 | reference stream (gaps/src/len) | at bounded floor (stage 37) | yes | NOT YET |
| 11 | mismatch stream | BUILT & VERIFIED: 241,365 B (was 242,209 est) | DONE | done (Illumina context-error lit) |
| 12 | sequence/literal coder | CM 1.906 bpb vs xz 2-bit, no middle ground | partial | yes (DNA-COMPACT/GeCo3) |
| 13 | FWD self-match after RC converges | untested post-iteration | NO | no |
| 14 | MAXCAND sweep post-iteration-fix | untested | NO | no |

Working order: 11 (mismatch coder, biggest unknown) -> 4 (mapping, diagnose
properly this time) -> 7 (make iteration a real stage) -> 13/14 (cheap re-checks
now that the matcher bug is fixed) -> 8/9/10 (confirm floors with fresh eyes) ->
12 (middle-ground coder).

## Results from this session's deep pass (2026-08-29, continued)

| finding | net effect | verified | cost |
|---|---|---|---|
| mismatch stream: real adaptive coder (ref-conditioned exclusive) | -844 B | round-trip VERIFIED | negligible |
| MEM iteration: RC-only to convergence (5 passes) | -42,444 B | round-trip VERIFIED | ~1.1s multi-threaded |
| MEM iteration: +1 FWD pass on converged residual | -35,436 B | round-trip VERIFIED | ~0.8s |
| MEM iteration: further alternation (3 more rounds) | ~+680 B (net negative, noise) | round-trip VERIFIED | not worth it |
| **combined (mismatch + RC-converge + 1 FWD pass)** | **-78,044 B** | | ~2s |

Speed-config (fast xz coder) total: 6,567,830 -> 6,489,786.
Lead over PgRC2 (7,035,682): 6.65% -> **7.76%** -- matches what previously
required the 19s CM coder, now reached with the fast coder plus ~2s of
mechanical MEM iteration.

MAXCAND sweep (64/128/256/512) on the second pass: saturates almost
immediately (301,420 -> 301,784 bases, +0.05%). Confirms the multi-pass gain is
NOT a candidate-cap artifact -- it is the greedy left-to-right parse leaving
real matches unreachable until an earlier accepted match is removed and the
position becomes visible again. Classic LZ optimal-parsing effect, not a bug.

CRITICAL CHECK: the coordinate bug found and fixed in this session's own
verification tooling (real_mem.cpp) was NOT present in the actual pipeline
(46/47_*.cpp) -- the real pipeline uses a separate RC bitmap and a whole-array
positional flip afterward, which is the correct conversion. Confirmed by
reading the exact call site. Every previously reported MEM number in this
project stands unchanged.
