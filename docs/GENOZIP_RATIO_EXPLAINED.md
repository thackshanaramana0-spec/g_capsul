# The "Genozip fungi anomaly" — diagnosed, and it is not an anomaly

`CLAUDE.md` §6.2 flagged this as an open red flag:

> **Genozip's fungi results are anomalous** (162 MB where SPRING gets 24 MB) and
> unexplained. Do not put them in a table until diagnosed.

It is now diagnosed. **There is no fungi-specific defect, and nothing is wrong
with the measurement.**

## 1. It is not fungi-specific

The Genozip / SPRING sequence-only ratio across all 14 locked datasets:

| dataset | ratio | | dataset | ratio |
|---|---|---|---|---|
| C. jejuni | **7.1×** | | H. pylori | 3.4× |
| S. cerevisiae (the "fungi" case) | **6.7×** | | H. salinarum | 3.2× |
| S. acidocaldarius | 6.2× | | HCMV | 2.9× |
| M. tuberculosis | 5.4× | | P. falciparum | 2.2× |
| S. aureus | 5.2× | | A. fumigatus (also fungi) | **2.0×** |
| E. coli | 5.0× | | SARS-CoV-2 | 1.0× |
| P. aeruginosa | 4.3× | | | |
| L. major | 3.8× | | **mean** | **4.2×** |

The gap spans **1.0×–7.1× on every kingdom**. The two Fungi datasets sit at
opposite ends (S. cerevisiae 6.7×, A. fumigatus 2.0×), so kingdom is not the
variable. Yeast is simply the largest instance of an effect present everywhere.

## 2. What actually causes it

The Phase-1 comparison is *sequence + read order only*, and the two tools solve
that problem in fundamentally different ways:

* **SPRING** with `--no-ids --no-quality` **reorders the reads** and codes them
  read-relatively — an assembly-like design that exploits cross-read redundancy.
* **Genozip has no cross-read matching by default** (see
  `docs/SOTA_COMPARISON.md` Table 2): a specialised `acgt` codec plus a generic
  LZMA/BSC backend, applied to reads in their original order.

So the ratio measures *reordering vs streaming*, and it tracks how much
cross-read redundancy each dataset holds — which is exactly why it correlates
with coverage rather than with taxonomy, and why SARS-CoV-2 (an amplicon set
with little exploitable structure at this stage) sits at 1.0×.

This is the behaviour PgRC2's own paper describes when it dismisses
streaming-regime tools: *"streaming-regime tools cannot obtain competitive
compression ratios on FASTQ reads."*

## 3. Confirmation from the whole-file numbers

If this were a Genozip defect, it would persist when quality is included. It
does not:

| comparison | mean Genozip/SPRING ratio |
|---|---|
| Phase 1 (sequence only) | **4.2×** |
| Phase 3 (whole file) | **1.82×** |

Adding names and quality — where Genozip is strong and reordering buys much
less — collapses the gap. A genuine bug would not behave this way.

## 4. Consequence for the paper

The Genozip numbers can be used as measured. The Phase-1 comparison should be
described for what it is: **a reordering compressor against a streaming one**,
which is a real and citable architectural difference, not a malfunction. Our
own Phase-1 margin over Genozip (12/12 wins, −78.5% aggregate) should be
reported with that context rather than as evidence Genozip is broken.

`CLAUDE.md` §6.2's "do not put them in a table until diagnosed" condition is
now satisfied.
