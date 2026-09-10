# T2.5 independently re-executed, 2026-09-10 — every value confirmed

T2.5 was the one table the raw-log cross-check could not verify: the sweep
log's end-of-run section prints it with `column -s, -t "$CSV7"`, which
pretty-prints the very CSV it is meant to corroborate. That is circular. So
it was re-executed instead.

**Result: all 16 published values reproduce exactly.**

    published CSV                          fresh re-run 2026-09-10
    CAPSULE     SNV   472  7  101          CAPSULE_SNV      TP=472 FP=7  FN=101
                      P=0.985 R=0.824                       P=0.985 R=0.824
                      F1=0.897                              F1=0.897
    DiscoSNP++  SNV   373  8  200          DiscoSNP++_SNV   TP=373 FP=8  FN=200
                      P=0.979 R=0.651                       P=0.979 R=0.651
                      F1=0.782                              F1=0.782
    CAPSULE     INDEL  51  8   61          CAPSULE_INDEL    TP=51  FP=8  FN=61
                      P=0.864 R=0.455                       P=0.864 R=0.455
                      F1=0.597                              F1=0.597
    DiscoSNP++  INDEL  44  3   68          DiscoSNP++_INDEL TP=44  FP=3  FN=68
                      P=0.936 R=0.393                       P=0.936 R=0.393
                      F1=0.553                              F1=0.553

Every TP, FP, FN, precision, recall and F1 identical for both tools and both
variant classes. Raw log: `rerun_20260910.log`.

Intermediate figures also reproduced: 154,580 mixed reads, 373,982 confident
bp in the window, 759 truth sites, 13,173 contigs, 609 lifted calls.

## What this is

Cooke, Wedge & Lunter, *Genome Research* 2022 (PMC8805713): two real diploid
GIAB individuals' real reads (HG003 + HG004) concatenated into one 4-copy
sample, truth = the union of their real GIAB v4.2.1 calls. Nothing simulated.
Called **from the archive** (the log's `[call] compressing to an archive, then
calling FROM it (no FASTQ read)` line), scored by `rtg vcfeval`, with
DiscoSNP++ run on the identical mixed reads.

## Reproduce

```bash
PATH="$HOME/DiscoSnp:$PATH" bash scripts/run_tetraploid_bench_capsule.sh \
    <encoder> "$(pwd)/scripts" ~/refs/chr20.fa \
    HG003 HG004 4 20:3000000-3400000 <outdir>
```

About 51 seconds on a 12-core box; needs network access to the GIAB S3/HTTP
endpoints for the read windows and truth VCFs.
