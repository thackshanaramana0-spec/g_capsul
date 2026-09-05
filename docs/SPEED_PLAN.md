# Encoder speed: state, plan, and the measure-twice rule

## Where we are (2026-09-05)

HG002 **914.80 s -> 513.16 s (1.78x)**, HG005 **2523.56 -> 1256.20 s (2.01x)**,
archives byte-identical on every dataset swept. The win came from fork DEPTH
(round 2 and the mapping search each ran 8x with 2 distinct results), not from
cache/TLB tuning.

Remaining 513 s on HG002:

```
shared prefix 235 s : load 14 | round 1 35 | round 2 66 | emit 7 | mapping 113
per-cand tail 278 s : second-region sweep 8-158 | MEM run() 80-116 | index 1 | trim 2 | coding
```

## Levers, ranked by (payoff x confidence) / risk

| # | lever | mechanism | payoff | risk | size cost |
|---|---|---|---|---|---|
| **A** | `rpk`/`woff`/`rlen` huge pages | post-fill `madvise` was inert on the 445 MB array the mapping loop random-accesses 1.93e9 times | unknown; touches round 2 + mapping = 179 s of the shared prefix | very low: memory layout only | none |
| **B** | share the MEM main self-match | `pg[0,main_pg_end)` (152.7 MB) is group-invariant; it is 70-87% of `run()` | ~80-120 s wall | medium-high: ~400-line hoist | none |
| **C** | raise candidate concurrency K | `per_child = max(peak, 1.5 x input)`; actual peak is 0.68 x input, so K=4 where 8 fits | unknown | low: env-tunable, no code change | none |
| **D** | shrink the candidate grid | 8 points cost the whole tail | large | none technical | **0.010% vs 2nd best, 0.586% vs worst** |

**D is not mine to decide.** We beat SPRING by 4.0%, so 0.586% is ~15% of the
margin. Surfaced with numbers; not acted on.

**Not shareable, checked:** the second-region sweep (8 s -> 158 s across
candidates). Low-MAXMAP candidates append 451k-477k reads vs 158k-172k, and
chains over a superset are not chains over a subset.

## The measure-twice rule (why this is written down)

Three separate wrong numbers were produced today by sloppy measurement: a
perf-attached run reported "22% faster" for a 1.8% change; a THP test launched
while another job was running reported a 4x regression; and a `grep -E "round 2
(assembly)"` matched nothing because `(...)` is an ERE group, reporting n=0 for
the very stage under study. Also two `pkill`/`pgrep -f` patterns matched their
own command line and killed the shell instead of the job.

So, for every change from here:

1. **Correctness gate first, timing second.** Byte-identity on a fixed 4-dataset
   subset (H. salinarum, S. acidocaldarius, E. coli, P. falciparum) before any
   timing is believed.
2. **Two independent full-scale runs**, not one, on a machine verified idle by
   `ps -eo stat,%cpu,time` -- a process at 0.0% CPU can still have been alive 13
   hours.
3. **No perf attached** to a run whose wall time will be quoted.
4. **Fixed-string extraction** (`grep -F`) for stage names containing `()`.
5. One change measured at a time; `A_only` source kept aside so A and B never
   share a measurement.
6. Full 19-dataset sweep before anything is called generalized.
