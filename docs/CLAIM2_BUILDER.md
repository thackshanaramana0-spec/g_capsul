# CAPSULE Claim 2 builder — status, 2026-09-02

**Data/tools inventory + the CAPSULE-adapted benchmark scripts live in
`docs/CLAIM2_DATA_AND_TOOLS.md`** — read that first if you're about to
actually run a GIAB/synthetic benchmark rather than extend the caller.

Reference-free variant calling ported from the outer ARCS project's caller
(`/root/arcs-clean/src/caller.cpp`) and wired directly onto CAPSULE's own
assembly. This is the first working end-to-end version — smoke-tested, not
yet GIAB-validated. Read this before extending it.

## What's built

- `stages/106_inprocess.cpp`, gated on `CAPS_CALL=1` (zero cost unset, same
  convention as `CAPS_NAMES`/`CAPS_QUAL`):
  - `g_contig_spans` captured directly inside the two existing chain-emission
    loops (main pg round-1/round-2, and the second-region sweep) — every
    chain is its own contiguous `[start,end)` span of `pg`, recorded BEFORE
    the MEM self-match stage runs. MEM only changes which pg bytes get
    literal-encoded, never `pg`'s content, so these spans stay valid archive-
    wide.
  - After assembly + pigeonhole mapping finish, every unique read's `ppos`
    is resolved to a contig id via binary search over span starts, then
    expanded to every ORIGINAL read through the existing `orig2uid` array
    (no new dedup/containment logic — reused as-is; containment never flips
    orientation, verified by reading that code path, so `read_rc` for an
    aliased original read is just `prc[orig2uid[orig]]` directly).
  - A second, lightweight FASTQ pass reloads original seq+qual in file order
    (mirrors the `>1023bp` skip the main load pass applies, so indices stay
    aligned with `orig2uid`) — kept separate from the main load pass so the
    no-`CAPS_CALL` path's memory footprint is completely unaffected.
- `include/caps_caller.h` — the caller itself, ported from `caller.cpp`
  near-verbatim: internal canonical-31-mer counting (no DSK), frozen SNV
  filters (`HDMAX=2, MAF=0.20, DHI=2.5, KHI=1.1, MC=3, TRI=0.12, HALF=15` —
  **not re-tuned, ported exactly**), contig-bubble indel extraction, and the
  cross-contig SNV bubble pass.

## The one real architectural finding from porting

ARCS's own assembler (`vodbg_pg`, Method B) evidently merges most
heterozygous sites into ONE shared consensus contig, with the minor allele's
reads placed onto it via mismatch-tolerant mapping — so ARCS's pileup-based
SNV pass is the dominant signal and its cross-contig bubble pass
(`ARCS_XSNV`) is experimental/opt-in, off by default.

**CAPSULE's assembler does not behave the same way.** Its round-1/round-2
chaining is EXACT suffix-prefix overlap (no mismatch tolerance in the chain
step itself — tolerance only enters later, in pigeonhole mapping). That
means a chain breaks at EVERY heterozygous site by construction: the two
haplotypes fragment into separate short contigs there, not into one merged
consensus with mismatches. Confirmed on the synthetic smoke test below —
0 pileup-column candidates, 62 contigs from ~24,000 reads over a 60 kb
diploid genome. **So for CAPSULE, the cross-contig SNV bubble pass is the
PRIMARY signal, not a fallback** — shipped default-ON here (`CAPS_NO_XSNV=1`
to disable), the reverse of ARCS's own default. This is a structural
consequence of CAPSULE's stricter, no-mismatch-tolerance chaining, not a
tuning choice.

## Smoke test (synthetic, not yet GIAB)

60 kb random diploid genome, 30 planted het SNVs, 100 bp reads, ~40× total
(20×/haplotype), reads randomly RC'd 50%:

```
[CAPS-CALL] contigs=62 H=27 candidates=4 SNVs=0 indels=0
[CAPS-CALL] xcontig SNVs=20
[CAPS-CALL] 20 records
```

**20/30 planted SNVs recovered (66.7% recall on a first pass), correct VCF
format, allele fractions 0.29–0.69 (consistent with heterozygosity), no
crash, no false-position artifacts observed.** This is a real signal, not
noise — but it is NOT a GIAB-comparable F1 number. No precision/recall
scoring was done (no coordinate lift yet — see below), and the synthetic
generator here is a fresh one-off script, not `sim_indel_bench.py`/
`sim_polyploid.py` from the outer project.

Confirmed zero-cost when `CAPS_CALL` unset: plain archive build produces the
same `ARCHIVE_TOTAL` as before this change.

## What's NOT done yet (next steps, in order)

1. **Coordinate lift + real scoring.** Port `scripts/lift_vcf.py`
   (`hapflank_lift`) to place CAPSULE's contig-coordinate calls onto a real
   reference, then score with `rtg vcfeval` against GIAB truth — same
   pipeline the outer project already validated. Nothing here has been
   scored against ground truth yet; the 20/30 above is an eyeball count on
   synthetic data with no lift step.
2. **Real GIAB run.** Reuse `~/refs/chr20.fa` / `~/giab_truth/` (already on
   disk from the outer project's Phase 1 download) on a CAPSULE-encoded
   chr20 region, same regions the outer project already used, for a
   directly comparable number.
3. **Indel path is ported but untested.** The indel-bubble code compiled and
   is wired identically to the SNV bubble pass, but the smoke test above
   planted no indels — needs its own synthetic or real-GIAB pass before
   trusting it.
4. **Polyploid path (`CAPS_PLOIDY=k`) is ported but untested** — same
   caveat as indels.
5. **Decide the streaming-scale question.** The second FASTQ pass loads all
   original seq+qual into memory (`std::vector<std::string>`) — fine at
   chr20/GIAB-region scale (the outer project's own Claim 2 scope), not
   suitable for the largest locked datasets (C. elegans/T. cacao scale) if
   Claim 2 is ever run there. Not needed yet since Claim 2's scope is chr20
   regions, not whole-genome.
6. **`ARCS_INDEL_DEBUG`/`ARCS_DUMP_CONTIGS` env-var equivalents** are ported
   as `CAPS_DUMP_CONTIGS` only; indel debug logging was left out of the port
   (low priority, add if debugging real-GIAB indel misses later).

## Files

- `include/caps_caller.h` — the caller (new)
- `stages/106_inprocess.cpp` — contig capture + CallData population + FASTQ
  re-read + caller invocation, all behind `CAPS_CALL=1`
- No build-script changes needed (header-only; `scripts/build106.sh` already
  includes `-I"$HERE/include"`)

## Reproduce the smoke test

```bash
CAPS_CALL=1 CALL_VCF=/path/calls.vcf INPUT=/path/sim.fq \
  ARCHIVE=/path/sim.capsule BEST=/tmp/caps_call_enc \
  bash scripts/encode_adaptive.sh
grep -v '^#' /path/calls.vcf | wc -l
```

Note: `encode_adaptive.sh` forks 4 MAXMAP/MINOV candidates (A3) that all run
`CAPS_CALL` and all write to the SAME `CALL_VCF` path — they overwrite each
other rather than racing byte-for-byte (each write is a full, complete
`fopen`+write+`fclose`), so the file left on disk is whichever candidate's
write lands last, not necessarily the smallest-archive candidate's calls.
**This is a real, currently-unaddressed interaction** — for anything beyond
a smoke test, invoke the binary directly with one fixed MAXMAP/MINOV (skip
the sweep) so there's no ambiguity about which candidate's VCF you're
looking at.
