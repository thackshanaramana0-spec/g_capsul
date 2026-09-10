# T3.4 independently reproduced, 2026-09-10 — and the undocumented flag that blocks it

**Outcome: T3.4's published HG002 row is CONFIRMED, exactly, field for field.**
This is the opposite result from the T2.4 audit in the sibling directory, and
the contrast is the point: the same audit method that withdrew one number
confirmed another.

    published CSV row:  HG002,100,18,82,0,85,15,0,2771,3023,DONE
    fresh re-run:       T34,HG002,100,18,82,0,85,15,0,2771,3023,0,15
                                      ^^ ^^ ^  ^^ ^^ ^  ^^^^ ^^^^
                        every field identical

Reproduced from a **freshly built archive** (573,767,964 B, matching the known
production figure) using a **freshly built decoder** from current source.

---

## The finding that matters for reproducibility

Getting there took four attempts, and the first three all produced
**coordinate = 0**, not 18 — the same "0" that this project already
identified once as a consensus-emitting artifact and corrected. The cause is
not a defect. It is an **undocumented required flag**:

    CAPS_PILEUP=1 capsule_decode index <archive> <archive>.qidx

Without `CAPS_PILEUP=1`, `capsule_decode index` writes a sidecar carrying the
pseudogenome and the placements but **not each read's own deviations** (see
`stages/capsule_decode.cpp`, the `WANT_PILEUP` branch — without the flag it
`return 0`s before the deviation-append step). `query` then reconstructs each
read as the **consensus** at its position, so every returned read agrees with
every other, no allele difference can appear at a heterozygous site, and the
coordinate arm scores 0 **by construction rather than by measurement**.

This was verified as cause, not guessed:

| attempt | sidecar | coordinate arm | log |
|---|---|---|---|
| 1 | none (stale pre-fix decoder) | 0 | — |
| 2 | none (current decoder) | 0 | `run_without_CAPS_PILEUP_gives_0.log` |
| 3 | built WITHOUT `CAPS_PILEUP=1` (120,274,641 B) | 0 | — |
| 4 | built WITH `CAPS_PILEUP=1` (172,045,852 B, "+1,804,605 reads carrying 7,466,871 deviations") | **18** | `run_with_CAPS_PILEUP_matches_published.log` |

The sidecar build log for attempt 4 is `sidecar_build_with_CAPS_PILEUP.log`.

**`CAPS_PILEUP` appeared in zero documentation files before this audit.** A
reviewer following the documented reproduction steps would have hit exactly
this wall and concluded the published 18/100 was unreproducible — the same
wrong conclusion I reached three times. Now documented in
`benchmark/documentation/REPRODUCE_EVERYTHING.md` §7 and in `RESULT_CODE.md`'s
T3.4 entry as a REQUIRED SETUP row.

## Why the sidecar is required at all, and why that is honest

The per-read deviations cannot be reached from the archive on demand:
`mm_sym` is coded with an **adaptive model in read order**, so read *k*'s
deviations require decoding all *k−1* before it. That is a property of the
coding, not an oversight — decompression only ever walks reads in order. The
sidecar breaks that ordering dependency by decoding once, up front.

This is already stated as a real architectural limit in
`paper/LIMITATIONS.md` §6. The sidecar is **not part of the archive** and its
cost is not counted in any compression number.

## Reproduce

```bash
# 1. Build the archive (or use the one Phase 1 kept)
CAPS_SPANS=1 CAPS_NAMES=1 CAPS_QUAL=1 \
  INPUT=/data/fastq/HG002_pooled.fq ARCHIVE=HG002.capsule BEST=<encoder> \
  bash scripts/encode_adaptive.sh          # -> 573,767,964 B

# 2. Build the sidecar WITH the flag. This is the step that is easy to miss.
CAPS_PILEUP=1 <capsule_decode> index HG002.capsule HG002.capsule.qidx

# 3. Run the fidelity test (query finds the .qidx automatically)
bash scripts/run_locus_fidelity.sh <capsule_decode> HG002.capsule \
     ~/refs/chr20.fa HG002 <outdir> 100
```

Roughly 6 minutes for the archive, 15 seconds for the sidecar, 5 minutes for
the 200 queries on a 12-core box.

## Scope of what this confirms

The HG002 row, which is 100 of T3.4's 400 sites. HG003/HG004/HG005 were not
re-run in this audit; they use the identical script and the identical
mechanism, and nothing about them was found suspect, but they carry the same
"confirmed by inheritance, not by re-execution" caveat that any unre-run row
does.
