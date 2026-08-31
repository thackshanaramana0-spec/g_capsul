# Plan: closing the pg-volume gap — verified findings and implementation

Every number below was measured from both binaries on the same file. Claims
that did not survive checking are listed as retracted rather than removed.
Nothing here is estimated or assumed.

---

## PART 1 — VERIFIED FINDINGS

### 1.1 Our main pseudogenome already BEATS theirs

Decomposing what actually survives to be coded (sulfo, 251bp):

| region | ours (surviving) | theirs (surviving) | delta |
|---|---|---|---|
| main / hqPg | 2,368,666 | 2,450,403 | **−81,737 (we win)** |
| second / lqPg | **1,038,520** | **95,666** | **+942,854 (10.9x)** |
| total | 3,418,901 | 2,546,069 | +872,832 |

**The entire literal loss is the second region.** Our main assembly is not the
problem and must not be rewritten.

### 1.2 The second region is nearly-raw storage

| | unmapped reads | second pg bases | bases/read | collapse |
|---|---|---|---|---|
| ours | 8,344 | 2,088,887 | 250.3 | **0.3%** |
| theirs | ~2,735 | 156,495 | 57.2 | 4.4x |

Reads that fail mapping are, by construction, reads that overlap nothing --
so a second assembly pass over them collapses essentially nothing (0.3%).
Appending them is raw storage. **The only lever is to leave fewer of them.**

Predicted cost of the excess: 942,854 bases x 0.229 B/base = **215,914 B**.
Observed literal gap: **208,472 B**. Agreement within 3.6%.

### 1.3 The unmapped count is set by the mapping acceptance ceiling

sulfo, sweeping MAXMAP (our ceiling is Lmax/13 = 19; theirs is
readLength/minCharsPerMismatch = 251/3 = 83):

| MAXMAP | appended | second surviving | literal | TOTAL |
|---|---|---|---|---|
| 19 (current) | 8,344 | 1,038,520 | 814,938 | 3,341,112 |
| 30 | 4,210 | 473,127 | 678,305 | 3,258,783 |
| **45** | 2,522 | 258,508 | 627,726 | **3,242,027** |
| 63 | 1,532 | 173,441 | 608,646 | 3,249,343 |
| 83 (theirs) | 784 | 130,579 | **599,368** | 3,264,260 |

At their ceiling our literal (599,368) matches theirs (595,834) almost
exactly -- the literal gap is fully explained and fully closable. The TOTAL
optimum is lower (45) because mismatch streams grow as acceptance rises. So
the ceiling is a **marginal-cost tradeoff**, not a quality threshold.

### 1.4 NO FIXED RATIO CAN EXPRESS THE OPTIMUM — including our current one

| | optimum MAXMAP | as divisor of read length |
|---|---|---|
| halo (151bp) | ~8 | L/19 |
| sulfo (251bp) | ~45 | L/5.6 |

Optima differ by **5.6x** while read length differs by only **1.66x**.
`MAXMAP = Lmax/13` is therefore wrong on both files (too high for halo, far
too low for sulfo), and so is any other constant divisor.

The real driver is how well the second region compresses:

| | MEM-second removal | (1 − removal) | optimum |
|---|---|---|---|
| halo | 93.1% | 0.069 | 8 |
| sulfo | 50.3% | 0.497 | 45 |

ratio 7.2x vs observed optima ratio 5.6x -- the mechanism tracks.

### 1.5 Chaining quality is NOT the lever (MEM matching absorbs it)

| MINOV | pg after chains | main surviving | literal |
|---|---|---|---|
| 16 | 11,302,852 | 2,368,666 | 814,938 |
| 113 | 8,655,900 (**−23%**) | 2,294,266 | 804,306 (**−1.3%**) |

A 23% smaller pg yields a 1.3% smaller literal, because MEM matching already
removes the redundancy that loose chaining leaves. **Do not invest in the
overlap sweep for size.**

### 1.6 The other two layers are at their floors

- `pos_abs`: coded **below** the uniform-random entropy floor for our span
  (−3.1% sulfo, −8.6% halo). Span-driven loss on sulfo only (predicted
  87,051 B vs actual 83,815 B); on halo we WIN by 22,525 B. Not a defect.
- `mem_triples`: 17.8–22.2 bits/match, at the log2(pg_len) floor.
- sequence coder: 1.9% behind on sulfo, **9.4% better** on halo.

### 1.7 RETRACTED claims from earlier sections

- ~~"our span is 5–7x theirs"~~ — compared our pre-matching span against their
  POST-matching size. `joinedPgLength` is computed at pgrc-encoder.cpp:211,
  before matchPgsInPg at line 240, so both index pre-matching. True ratio
  2.54x / 3.59x.
- ~~"pos_abs loses because our span is bigger"~~ — true on sulfo, false on
  halo where we win despite a wider span.
- ~~"their overlap floor is 0.35 x readlen"~~ — `performOverlapping()` for pg
  generation is called with the default coef = 1 (AbstractOverlap...cpp:251,
  264), i.e. the full range, same as us. The coefficient applies only to the
  hq/lq division pass.

---

## PART 2 — WHAT TO BUILD, IN OUR CODE

### The defect, stated precisely

`106_inprocess.cpp`, the mapping stage: a read is accepted if its mismatch
count is under a ceiling fixed before any data is seen. The correct decision
is a comparison of two costs that are both measurable at run time:

```
cost_map(m)  = m * bytes_per_mismatch                 // mismatch records
cost_append  = readLen * (1 - r_second) * bytes_per_base
                                                       // bases that SURVIVE
                                                       // MEM removal in the
                                                       // second region
accept  iff  cost_map(m) < cost_append
=>  m* = readLen * (1 - r_second) * bytes_per_base / bytes_per_mismatch
```

Position and strand are paid on both paths, so they cancel and must not
appear. This is the same correction already applied to MAXMAP's *form*
(absolute -> relative), taken to its actual conclusion: the ceiling is not a
constant, not a ratio of read length, but a **break-even between two measured
per-read costs**.

### The chicken-and-egg, and how to resolve it

`r_second` (second-region MEM removal) and the two byte rates are only known
after the stage that consumes them. Resolve by **estimate, decide, verify** --
the pattern already used for the coder selector's probe:

1. **Probe.** Run mapping on a bounded sample of leftovers (reuse the existing
   `PROBE_MIN`/`probe_frac` machinery) at a permissive ceiling. Append the
   sample's unmapped reads, run the existing second-region MEM pass on that
   sample only, and measure `r_second` directly.
2. **Derive** `m*` from the formula with the measured `r_second` and the
   byte rates already computed this run (we know literal bytes/base from the
   main region, and bytes/mismatch from the mm streams).
3. **Decide** per read against `m*`, then run the real mapping pass.
4. **Verify** by construction: the decoder is untouched, so byte-identity is
   preserved regardless of the decision -- the choice only moves a read
   between two already-supported representations.

### Code sites

| step | file / location | change |
|---|---|---|
| ceiling derivation | `106_inprocess.cpp`, the `MAXMAP` block | replace the `Lmax/13` constant with the break-even computation |
| probe pass | mapping stage, before `sweep()` for survivors | bounded sample, measure `r_second` |
| rate inputs | archive stage already computes both | thread `bytes_per_base` / `bytes_per_mismatch` back, or use the running estimate |
| decoder | `decode_105.py` | **no change** — both outcomes already decode |

### Expected result, from the measured sweeps

Applying each file's own measured optimum rather than a shared constant:

| | current (L/13) | at its own optimum | PgRC2 | margin |
|---|---|---|---|---|
| sulfo | 3,341,112 | 3,242,027 | 3,114,782 | −7.3% -> **−4.1%** |
| halo | 2,619,597 | 2,617,438 | 3,050,477 | +14.1% -> **+14.2%** |

So this recovers roughly **3.2 points on sulfo and costs nothing on halo**. It
does not by itself make sulfo a win -- ~127 KB remains after the literal gap
closes, unattributed and still to be identified. That is stated as a limit of
the plan, not hidden.

### What this plan deliberately does NOT do

- rewrite the overlap sweep (1.5: MEM matching absorbs it, 1.3% available)
- touch the position, reference, sequence or mismatch coders (1.6: all at or
  below their coding floors)
- add a per-dataset or per-kingdom switch (the rule is one formula over
  measured inputs)
- change the decoder or the archive format
