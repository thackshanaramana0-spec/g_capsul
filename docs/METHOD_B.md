# Method B — the graph caller, standalone

**Standing directive, set 2026-09-04. Read this before touching the dBG code.**

## The rule

Method B is a caller built from scratch on the trunk. It consists of exactly:

    trunk (chaining + pseudogenome)  ->  kc (k-mer graph)  ->  bubbles  ->  VCF

It does **NOT** include, and must never be compared against or fall back on:

- the pileup (`parallel_loop`)
- `indel_pass`, the second `build_substrate`, `rc_reads`, `pkidx`

Those are the OLD caller. They are a different method at a different order of
RAM and runtime, and mixing the two produces meaningless numbers — a Method B
result quoted next to a 35 GB / 58 min figure from the old path says nothing
about either. **Report Method B on its own terms only.**

Run it with `CAPS_DBG=1 CAPS_DBG_ONLY=1`.

## Where it stands (HG002 r2, measured)

| | Method B | DiscoSNP++ (full chr20) |
|---|---|---|
| SNV F1 | **0.804** | 0.847 |
| precision | **0.946** | 0.951 |
| recall | **0.699** | 0.763 |
| caller time | **1.78 s** | — |
| caller peak RSS | **~111 MB** | — |

Progression, all measured on the same window:

| change | recall | F1 |
|---|---|---|
| initial (one strand, strict superbubble) | 0.468 | 0.636 |
| permissive pairwise walk | 0.512 | 0.666 |
| **+ both strand orientations** | **0.699** | **0.804** |

The strand fix was the single largest win: `kc` stores CANONICAL k-mers, so the
graph is bidirected and every node has two sides. The code walked only one.
GATB carries an explicit strand per node (`Graph.cpp:1599`); we had none.

## Remaining gap: recall 0.699 vs 0.763

The identified limiter is in `walk()`:

    if (++nc > 1) return false;    // ambiguous, give up

**We abandon a bubble the moment the two paths hit a tangle mid-walk.**
DiscoSNP++ does not: `expand()` (Bubble.cpp:598) RECURSES over every successor
pair via `expand_heart` rather than bailing. So at the entrance we are more
permissive than they are (no `checkBranching`), but during extension we are
strictly less permissive. That asymmetry is where the remaining recall is most
likely sitting.

Ranked candidates, to be tested one at a time against this baseline:

1. **`nc > 1` bailout** — port their bounded recursion over successor pairs.
2. **`is_tip`** may delete a real het branch whose local coverage dips.
3. **`MINC = 4`** — derived from the histogram valley, but still a floor that
   removes low-coverage het alleles.
4. **`BAL = 0.35`** — swept on this window; verify it is not costing recall.

## Rules that still apply

- Constants must be formulas over a measured input property, never fitted
  (CLAUDE.md rule 1). `MINC` already follows this; `BAL` does not yet.
- r2 is the tuning window. r3/na/r4/r5 are held out and must not be tuned on.
- Every change measured individually, reverted and recorded if it fails.
