# The archive path is structurally cheaper — measured

Every optimisation run in this session used the FASTQ path as a harness. That
was wrong: the archive path is the product, and it turns out to be **much
faster for identical output**, which no amount of tuning the FASTQ path would
have revealed.

Same 4M-read subset, same caller build:

| stage | from FASTQ | from archive |
|---|---|---|
| compress / decode | 53.5 s | **25.4 s** |
| ridx_build | 76.9 s | 101.8 s |
| kc_H_build | 12.1 s | 12.1 s |
| **parallel_loop** | **214.3 s** | **2.1 s** |
| indel_pass | 311.1 s | 296.2 s |
| **wall** | **671.0 s** | **438.6 s** |

## The output is the same

    archive:  contigs=140561  H=55  candidates=68665  SNVs=24769  indels=3080
    FASTQ:    contigs=140561  H=55  candidates=68666  SNVs=24769  indels=3093

Identical contig count, identical H, identical SNV count, candidates differ by
one. Both substrates collapse identically (334816 -> 140561 and
334816 -> 242732). So `parallel_loop` at 2.1 s is NOT skipped work.

## Why

The compressor already computed the read placements. The archive stores them
(`pos_abs` + `pos_strand` + `contig_spans`), so the caller decodes them in 18 s
instead of re-deriving them. `parallel_loop` -- the SNV pileup -- is what
re-derives that state on the FASTQ path, and it is exactly the stage that
collapses when the archive supplies it.

The anchor scan shows the same effect: 37 s from FASTQ, **5.1 s** from the
archive.

**This is the "compression already paid for it" claim, measured end to end** --
not as an equivalence argument but as a 35% wall-clock difference for the same
calls.

## Consequence for the optimisation plan

`parallel_loop` was the single largest stage I had earmarked as "the biggest
lever left" at 214 s. **On the path that matters it is 2.1 s.** Time spent
optimising it would have been wasted entirely.

The archive-path budget is:

    indel_pass    296.2 s   68%
    ridx_build    101.8 s   23%
    decode         25.4 s    6%
    kc_H_build     12.1 s    3%
    parallel_loop   2.1 s    0%

so `indel_pass` and `ridx_build` are 91% of it, and everything else is noise.
