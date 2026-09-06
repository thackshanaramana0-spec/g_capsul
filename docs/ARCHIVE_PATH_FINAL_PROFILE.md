# Archive path — complete profile after the structural work

## Result

    671.0 s   FASTQ path, unoptimised (where this work started)
    432.0 s   archive path, unoptimised
    213.5 s   archive path, six structural fixes

**3.1x, F1 bit-identical at every step:**

    SNV    TP=13858  FP=2624  FN=30717  P=0.8408  R=0.3109  F1=0.4539
    INDEL  TP=1217   FP=318   FN=6564   P=0.7928  R=0.1564  F1=0.2613

## The fully accounted profile

For the first time this session the timers sum EXACTLY to the stage total
(159.167 s accounted = 159.167 s `indel_pass`), so nothing is hidden.

| stage | time | share |
|---|---|---|
| build_substrate #2 | 73.1 s | 33.6% |
| pcluster scan | 47.6 s | 21.9% |
| build_substrate #1 (`ridx_build`) | 44.2 s | 20.3% |
| decode reads+qual | 18.4 s | 8.4% |
| emit + XSNV | 12.7 s | 5.8% |
| kc_H_build | 11.8 s | 5.4% |
| pcluster sort+select | 9.8 s | 4.5% |
| anchor scan | 7.0 s | 3.2% |
| everything else | ~15 s | 6.9% |
| **build_substrate TOTAL** | **117.3 s** | **53.9%** |

## Why build_substrate's 117 s is NOT waste

It is called twice with different `dup_frac` (0.45 for the SNV pileup, 0.92 for
the bubble passes), and the two calls are independent BY CONSTRUCTION -- both
take `cd_in`, the encoder's original placements, not each other's output.

Cost tracks surviving contigs almost exactly:

    call 1: 334,816 -> 140,561 contigs   44.2 s
    call 2: 334,816 -> 242,732 contigs   73.1 s
    ratio      1.73x contigs              1.65x time

The dominant term is read RE-placement, which is required whenever a read's
contig was collapsed away: 2.07M reads in call 1, 3.03M in call 2 (5.1M of 8M
total). Call 2 re-places MORE despite keeping more contigs, because at
dup=0.92 a different SET survives -- not a superset -- so fewer reads find
their original contig intact.

That loop is already `#pragma omp parallel for` and is genuine per-read work
(4M reads x 2 strands x offsets). There is no redundancy left to remove.

## Honest conclusion

**The structural levers are exhausted at 3.1x.** What remains is real
computation on 4M reads and ~65M k-mers:

* build_substrate 117 s -- two independent collapses, already parallel
* pcluster scan 47.6 s -- already threaded and bitset-filtered (was 112 s)
* decode 18.4 s -- entropy decode of 4M reads
* kc_H_build 11.8 s -- k-mer counting

Further gains would require either changing what is computed (a different
algorithm, or dropping the second substrate -- both change results) or
parameter tuning, which is explicitly out of scope.

## What the session actually removed

| fix | what it removed |
|---|---|
| `collapse_contigs` single walk | recomputing the same k-mers twice per contig |
| `claimed.reserve()` sized from input | rehashing a set from a 1M seed to tens of millions, twice per run |
| seed index parallel fill | serial `push_back` + serial sort over ~65M positions |
| pcluster anchor loop parallel | serial iteration over 21.6M anchors |
| per-pair read copy | a full `std::string` per candidate pair |
| 6b2 XSNV: filter then thread | 65M mostly-empty scheduler tasks |
| buffered contig dump | 1.05M unbuffered `fprintf` calls (benchmark-only) |

None is a parameter change.

## Instrumentation lesson

Four separate targets in this session were identified by READING CODE and were
wrong -- the `im` emit loop (0.02 s, a fix was written and ready), the `lvotes`
channel (dead code behind a flag), `6b2` (25 s not 142 s), and XSNV (14 s). Each
was caught by the `accounted X of stage` check, which prints the sum of the laps
against the stage total and says explicitly that a mismatch means the timers are
wrong.

The last and largest miss was structural: `g_ip_start` was set AFTER the second
`build_substrate` call, so the single most expensive item in the whole run --
73.1 s, 33.6% -- sat inside the stage total and outside every lap, while a timer
labelled "2nd build_substrate + setup" reported 4.7 s for the setup after it.
