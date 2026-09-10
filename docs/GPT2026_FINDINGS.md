# GPT2026 measured findings

Status: investigation in progress. No three-axis win has been established.
The caller and production defaults have not been changed.

## Reproduced baseline (2026-09-09 IST)

Encoder built from `dff89010a74d43c7ab79b7aeaddecce2397f7e57` with the
shipped build script, before edits. Subsequent external commits through
`660706e` change benchmark scripts/tables, not encoder/decoder/coder sources.
Those commits and the external benchmark output changes are preserved.

Invocation: the shipped adaptive wrapper, 12 OpenMP threads,
`CAPS_SPANS=1 CAPS_NAMES=1 CAPS_QUAL=1 CAPS_PHASE=1 CAPS_SWEEP_TIMING=1`.
No concurrent timed job. Free-disk and machine-idle checks are recorded in
each experiment's metadata. These are single observations, not repeat-run
confidence intervals. Diagnostic instruments are excluded from these builds.

| Input | Archive bytes | Wall s | GNU time max RSS KiB | Sampled tree peak PSS KiB | Full FASTQ verification |
|---|---:|---:|---:|---:|---|
| SRR2584863 | 68,429,027 | 12.00 | 2,014,900 | 1,936,841 | 1,553,259 records; 695,163,748 bytes identical |
| HG002 pooled | 573,767,964 | 207.10 | 4,950,664 | 6,789,263 | 12,604,917 records; 4,279,197,941 bytes identical |

Both archive sizes match the handoff's measurements. The 15 existing tests
passed with no failures. Full-FASTQ verification independently reconstructs
line 3 from the archive's names dictionary (E. coli mode 1; HG002 mode 0),
and compares all columns without sorting. Source FASTQ bytes are an oracle
only, never used to reconstruct output.

Raw results: `/tmp/gpt2026_work/{ecoli,hg002}_baseline/` contain the archive,
encoder log, GNU time, memory samples, invocation metadata, decoder log, and
full-FASTQ verification result. Baseline test log:
`/tmp/gpt2026_work/baseline_tests.log`.

## Memory measurement changes the interpretation

GNU time's maximum RSS is a single-process high-water value, not the maximum
physical memory of a live fork tree. On HG002, summed sampled PSS peaks at
6.47 GiB during coding (~191 s), while the single-process maximum is 4.72
GiB. Summed RSS at that sample is 7.65 GiB and double-counts shared pages.
PSS is sampled approximately every half-second plus collection overhead;
its maximum is a lower bound on the true peak. On E. coli, a short peak was
missed, so sampled PSS falls below the exact single-process high-water mark.
Keep both metrics; never substitute one silently for the other.

The coding phase explicitly copies the inherited 470,556,575-byte HG002
quality body into a candidate-local `QL`, then copies it again into a job
result. The names body (30,645,045 bytes) follows the same path. This is a
code-derived allocation mechanism, not yet an A/B memory saving. A borrowed
immutable stream view could remove these copies while preserving every
archive byte. It would need lifetime and serialization-order verification.

## Mapping: evidence to distinguish two causes

The fresh HG002 mapping stages cost 59.45 and 61.21 seconds. Round 1 costs
42.60 s; the two round-2 sweeps cost 21.79 and 26.60 s. The whole encoder's
quality body is 470.6 MB: historical sequence-only stream shares do not
describe the current full-FASTQ archive.

The mapping verifier reads `readMM[rid]`, which remains fixed during each
strand scan. Its local `hmm[t][slot[rid]]` minimum is updated during that
scan, but is not fed back into candidate rejection. The new compile-time
`GPT2026_MAP_PROFILE` counters measure candidate visits, actual verifications,
previous local minima, zero minima, repeated winning positions, non-improving
accepted candidates, and word comparisons after an already-known local bound.
No atomic counters are added to the hot loop. The instrument's own timing
does not establish an optimization benefit.

The completed diagnostic reproduces the baseline archive byte for byte.
For the forward scan with 3,663,240,292 valid verifications, 2,069,419,808
(56.5%) encounter an existing local minimum. Of 11,940,457,154 packed-word
comparisons, 2,957,306,660 (24.8%) occur after its implied bound is reached.
356,663,213 candidates pass the loose verifier but cannot improve the local
minimum. Only 6,075,480 are already locally exact: a zero-only shortcut would
miss most of this opportunity.

The reverse scan is a different workload: 3,508,727,741 verifications, but
only 43,191,791 encounter a local hit. It starts with forward-scan bounds.
The final mismatch histogram is dominated by one-mismatch reads (892,491
of 1,951,461 placed reads in one candidate). Therefore the previous claim
that an 8-base filter is always inapplicable because MAXMAP is 29 is not a
proof: the actual `lim` is often `cur - 1`, not MAXMAP. Quantify that domain
before choosing a discriminator. This is an additional lead, not yet a
validated filtering implementation.

The diagnostic mapping stage takes ~106 s against ~61 s uninstrumented;
its counters have substantial cost. Only event counts and output equality
are evidence from this build. Its timings must not enter the speed table.

Relevant comparison: PgRC2's local `matching/copmem/CopMEMMatcher.cpp`,
`processApproxMatchQueryTight`, sets `maxMismatches = res - 1` immediately
after an improvement and returns at the allowed minimum. Its parallel caller
in `matching/ReadsMatchers.cpp` assigns independent reads to threads. This
motivates testing read-owned best-match state in CAPSULE. Merely increasing
seed specificity does not address an already-known good match being checked
again.

Exact-output restructuring investigated after that diagnostic:
represent a read's best candidate as one ordered atomic integer containing
the mismatch count and its original enumeration rank. Workers may prune
against that value; candidates tied in mismatch count must still compete by
the original rank. This could replace per-text-chunk slot and hit arrays,
avoid repeated verification, and make memory O(reads) instead of
O(workers * reads), without changing seeds or placements. Prove the rank
against actual chunk and seed-part ordering, including RC, before measuring.

## Read-owned mapping experiment: insufficient benefit

`GPT2026_MAP_OWNED` is compiled out normally and requires a runtime value of
1 in the experimental build. One ordered atomic per read replaces the
per-worker slot arrays and per-chunk hit lists. The ordering key encodes
mismatch count, original seed coordinate and seed part. This preserves
first-on-equal selection independently of worker scheduling; an earlier
candidate must still be allowed to tie a later incumbent.

| Input | Baseline wall s | Owned wall s | Owned max RSS KiB | Owned sampled peak PSS KiB | Archive |
|---|---:|---:|---:|---:|---|
| SRR2584863 | 12.00 | 11.68 | 2,044,200 | 1,967,453 | byte-identical |
| HG002 pooled | 207.10 | 204.95 | 4,934,252 | 6,933,642 | byte-identical |

HG002 mapping stages decrease from 59.45/61.21 to 57.95/58.99 s in these
single observations. Mapping live RSS drops by about 170 MB per candidate,
but whole-run peak physical memory does not decrease: coding still sets
the peak. Total time improves only 1.0%, within the range requiring repeated
measurements. This does **not** meet the structural performance target and
is not promoted. The earlier software counter reductions are not a measured
end-to-end speedup. Keep this experiment isolated and disabled while testing
the independent buffer ownership mechanism.

The combined experimental binary with **all runtime flags off** also
reproduces the E. coli archive exactly (11.65 s; max RSS 2,018,404 KiB).
This control reinforces why a 0.32-second E. coli difference cannot be
claimed as a reliable optimization.

## Buffer ownership experiment: memory saving, timing not yet accepted

An independent memory-only run releases the waiting
parent’s encoder snapshot after the final group fork, and moves already
compressed names/quality bodies through the child job pipeline. It changes
neither coding choices nor stream order. The mapping experiment is disabled
for this measurement.

| Input | Wall s | Max RSS KiB | Sampled peak tree PSS KiB | Archive |
|---|---:|---:|---:|---|
| SRR2584863 | 11.71 | 1,992,188 | 1,888,095 | byte-identical |
| HG002 pooled | 219.05 | 4,877,828 | 4,568,916 | byte-identical |

On HG002 the sampled tree peak decreases 32.7%, from 6.47 to 4.36 GiB.
The peak moves from coding (~191 s) to the initial quality/assembly overlap
(~49 s); parent RSS after release is 5 MB. The single-process high-water
RSS changes little, as expected when that earlier phase remains intact.
Wall time increases 5.8% in this first run, so this is a memory result,
**not an accepted simultaneous time/RAM win**. Alternating repeated timings
are required before deciding whether it introduces a time regression.

The memory archive was decoded afresh using the immutable baseline decoder:
12,604,917 HG002 FASTQ records / 4,279,197,941 bytes match completely.
`combined-tests` passed all 15 existing regression tests with owned mapping,
parent release and buffer moves enabled. The suite includes archive-only
export/coverage/query; its synthetic single-candidate encodes do not exercise
the outer group fork. The real E. coli/HG002 adaptive runs exercise that fork.

## Locked input domain audit

`stages/gpt2026_fastq_domain.cpp` scans the exact 15 locked inputs and HG002,
without encoding or normalizing them. The complete scan found zero reads
with N positions beyond 255, more than 255 Ns, or lengths beyond 1023.
ERR552797 and SRR40271341 do contain reads longer than 255, but none of their
Ns exceed the current position domain. Thus the suspected N-stream overflow
is **not reached by these locked files**. This is a domain audit, not a
substitute for full archive losslessness. Raw record counts and maxima:
`/tmp/gpt2026_work/locked_domain.tsv`.

## Seed-extension experiment: no demonstrated speed benefit

`GPT2026_MAP_EXTENSION` snapshots each posting's immutable strand-scan
mismatch bound beside an eight-base extension. Extension Hamming distance
is an exact lower bound on full-alignment distance; candidates without a
complete extension fail open. The original capped posting list and tie
order remain unchanged. A 32-bit posting word replaces repeated bound
lookups and permits rejection before loading the packed read.

E. coli flags off/on and HG002 flag on all reproduce baseline archives
exactly. E. coli off/on takes 12.06/12.27 s. HG002 takes 223.65 s with max
RSS 4,986,888 KiB and sampled tree peak PSS 7,002,618 KiB. This is **not a
speed win** in the first observation and is not promoted. The index adds
about 59-61 MB per candidate while a scan is live. Later unchanged stages
also vary from the original baseline, so repeated runs need stronger
contamination evidence before assigning the entire difference to this code.

## Run provenance and concurrency limitation

The checkout is being changed externally during this investigation, including
the decoder and benchmark tables. Those edits are preserved. The measured
baseline/experimental encoder binaries are fixed files in `/tmp`; the
full-FASTQ gates use the copied baseline decoder. External decoder edits
must not be attributed to GPT2026 or used to imply our caller changed.

The initial timing series checked machine idle only before each run. It
does not prove isolation throughout a run. The measurement helper now logs
file hashes and compares whole-machine busy CPU seconds against GNU time's
command CPU plus the monitor's own CPU. A large unexplained residual marks
possible contamination. A small residual is only supporting evidence, not
a proof of exclusive machine use. Repeated timings with this evidence are
required for promotion; the early single-run tables remain observations.

## Quality-stage decomposition and model capacity

The standalone diagnostic calls the existing quality scheduler and times
FQZ parameter selection and coding with thread CPU clocks. Its 56 HG002
trials spend 45.02 CPU s in parameter selection and 157.96 CPU s in coding
(including model setup). Its 470,556,575-byte quality body matches the actual
baseline archive exactly. All four strategies remain candidates; no trial
has been removed based on historical winners.

Each original quality context allocates a 256-symbol model even when the
selected alphabet has maximum symbol 39 or 41. A model occupies
`4 * capacity + 16` bytes, and each trial initializes 65,536 contexts.
The namespaced compact experiment dispatches to the smallest power-of-two
capacity covering **the actual selected model alphabet**, with a minimum
of four. For HG002 this selects 64: 17 MiB per model bank instead of 65 MiB.
Length and selector models retain capacity 256. Every active symbol,
frequency, ordering rule, normalization and arithmetic update is preserved.
Unused trailing zero-frequency entries cannot affect valid encoding.

Both uninstrumented isolated quality builds reproduce the body and index:

| HG002 quality stage only | Wall s | Max RSS KiB | Sampled peak PSS KiB | CPU s |
|---|---:|---:|---:|---:|
| Baseline | 32.92 | 3,029,356 | 2,653,507 | 206.71 |
| Compact model capacity | 32.16 | 2,439,380 | 2,041,914 | 206.29 |

Peak RSS is 19.5% lower in this isolated stage. Runtime is similar; the
small difference does not establish a speed benefit. Unexplained CPU load
is 0.083/0.058 core equivalents respectively, below the helper's contamination
threshold. These are **not full-encoder measurements**. Full integration,
capacity-boundary tests and the regression suite remain to be completed.

## Limits and outstanding work

Output equality has passed on two inputs for the mapping experiment; its
performance is insufficient. The combined experimental regression suite is
15/15; full locked-set archive validation remains outstanding. No competitor
superiority is established. The
three-axis, all-locked-datasets target remains open. `perf` is denied by the
current sandbox/kernel setting; current profiles are timers and software
counters only. PgRC2's supported scope must be matched explicitly: its
[upstream description](https://github.com/kowallus/PgRC) describes DNA-stream
compression and an order-preserving mode, not full-FASTQ preservation.
