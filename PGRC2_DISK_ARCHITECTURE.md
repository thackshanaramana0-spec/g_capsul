# PgRC2's disk/persistence architecture — real source analysis + empirical verification

Prerequisite analysis for the structural rewrite. Everything here is read from
the real cloned source (`/tmp/pgrc_study/PgRC`) and verified by running the
real release binary — no assumptions.

---

## 1. THE DECISIVE FINDING: their disk-staging is DEAD CODE in the release build

The disk-staging / memory-release-valve architecture we identified in an
earlier session as "the actual, deeper gap" **does not execute in the binary
we benchmark against.**

### Source evidence

`pgrc/pgrc-encoder.cpp` gates all 6 stage checkpoints on one condition:

```cpp
if (params->disableInMemoryMode || params->endAtStage == stageCount) {
    persistX();
    data.disposeChainData();
}
```

Both operands are unreachable in the release build:

| Operand | Default | How it can be set | Reachable in release? |
|---|---|---|---|
| `params->disableInMemoryMode` | `false` (`pgrc-params.h:72`) | **No CLI flag exists at all** — `grep` over the whole tree finds only the declaration and these 6 reads | **NO** |
| `params->endAtStage` | `UINT8_MAX` (`pgrc-params.h:69`) | Only `-E` → `setEndAtStage()` (`PgRC.cpp:170`), which sits inside `#ifdef DEVELOPER_BUILD` | **NO** |

And `CMakeLists.txt:228` defines `DEVELOPER_BUILD` **only** for the separate
`PgRC-dev` target:

```cmake
target_compile_definitions(PgRC-dev PUBLIC "-DDEVELOPER_BUILD")
```

The `PgRC` release target — the one built at `build/PgRC`, the one used for
every comparison in this project and in their paper — does not define it.
So `disableInMemoryMode || endAtStage == stageCount` is **always false**, and
`persistX()` + `disposeChainData()` at those 6 checkpoints **never run**.

### Empirical confirmation

Ran the real release binary in a clean directory and listed every file after:

```
in.fq        (input)
out.pgrc     (output)
```

**Zero intermediate files. Zero temp directories.** Confirmed on two separate
datasets (H. salinarum, S. acidocaldarius). The `mkdir(tmpDirectoryName)` at
`pgrc-encoder.cpp:71` is itself gated on `params->extraFilesForValidation`
(the `-V` flag), which is off by default.

**Conclusion: PgRC2 in its shipped configuration is a pure in-memory
compressor that touches disk exactly twice — read input, write archive.**

---

## 2. What their architecture ACTUALLY is: in-memory incremental compress-and-release

The real memory management lives in the `else` branches — the paths that DO
execute. The pattern is: as soon as a structure's final compressed bytes are
written to the output stream, the in-memory structure is destroyed.

From `executePgRCChain()` (`pgrc-encoder.cpp:108-252`), the paths that run:

| Line | Action | What it frees |
|---|---|---|
| 157 | `data.hqPg->disposeReadsList()` | HQ pg's reads list, right after its mapping stage |
| 169 → 183 | `compressLQPgReadsList()` then `data.lqPg->disposeReadsList()` | LQ reads list, immediately after compressing it |
| 189 → 203 | `compressNPgReadsList()` then `data.nPg->disposeReadsList()` | N reads list, immediately after compressing it |
| 205-206 | `delete(data.divReadsSets); = nullptr` | The entire divided reads set |
| 221 / 226 | `data.orgIdx2PgPos.clear()` / `data.rlIdxOrder.clear()` | The order/position array, right after compressing it |

Note the shape: **compress → write to `pgrcOut` → free**, repeatedly, so peak
RAM is bounded by the largest single live structure rather than the sum of all
of them. No disk involved.

The `SeparatedPseudoGenomeOutputBuilder` supports both targets via one switch
(`SeparatedPseudoGenomePersistence.cpp:688-696`):

```cpp
void initDest(ostream *&dest, const string &fileSuffix) {
    if (dest == nullptr) {
        if (onTheFlyMode())   // == !pseudoGenomePrefix.empty()
            dest = new ofstream(...);   // disk
        else
            dest = new ostringstream(); // memory
    }
}
```

In release runs the prefix is empty, so **every one of their 10 streams**
(`pgDest, pgPropDest, rlPosDest, rlOrgIdxDest, rlRevCompDest, rlMisCntDest,
rlMisSymDest, rlMisPosDest, rlOffDest, rlMisRevOffDest`) is an in-memory
`ostringstream`.

### The resume mechanism (also dead in release)

Every `prepareForXxx()` opens with `if (!data.divReadsSets) { ...reload from
disk... }` (e.g. `pgrc-encoder.cpp:271-284, 298-307`). That is what makes
`disposeChainData()` safe in dev mode: the next stage transparently reloads
what it needs. In release, the data is always still in memory, so these
reload branches never fire either.

---

## 3. The premise was inverted: WE are the disk-staging tool, not them

Listing what each pipeline writes for one dataset (H. salinarum, 460,501 reads):

| | PgRC2 (release) | Ours (locked seq+order) |
|---|---|---|
| Intermediate files written | **0** | **18** |
| Files | — | literal.txt, mem_triples.bin, mem_gaps.bin, mem_lens.bin, mem_srcs.bin, mem_srcdelta.bin, mm_count_per_read.bin, mm_ref.bin, mm_obs.bin, mm_pos.bin, mm_ctx3.bin, pos_abs.bin, pos_direct.bin, pos_strand.bin, pos_strand_direct.bin, orig2uid.bin, read_lengths.bin, n_reads.txt/n_indices.bin |
| Process model | single process, streams held in memory | assembly process + N separate coder processes, each re-reading from disk |

We already stage everything through disk — far more aggressively than they do.
**The "we have no disk-staging architecture and PgRC2 does" premise is exactly
backwards**, and any rewrite aimed at *adding* disk staging would be moving in
the wrong direction.

---

## 4. Real RAM measurements — we are at parity or better at these scales

Same file, same machine, `/usr/bin/time -v` peak RSS:

| Dataset | PgRC2 (whole run) | Ours: assembly | Ours: seqpar (L1) | Ours: xz (L4) | Our effective peak |
|---|---|---|---|---|---|
| H. salinarum (460K reads, 210 MB) | 185 MB | 112 MB | 41 MB | 84 MB | **112 MB** |
| S. acidocaldarius (517K reads, 333 MB) | 194 MB | 208 MB | 41 MB | 86 MB | **208 MB** |

Because our coders run as *separate sequential processes*, our effective peak
is the max of the stages, not the sum — 112 MB vs their 185 MB on the first
file (we win), 208 MB vs 194 MB on the second (7% behind, near parity).

The `seqpar` 820 MB fixed allocation documented in `LAYER_BY_LAYER_ANALYSIS.md`
§1b is **already fixed** — the adaptive `mmbitsFor(n)` change landed, and the
binary now measures 41 MB on these files.

Caveat, stated honestly: the earlier C. elegans-scale measurement (4M reads,
peak 1.30 GB, ~1.86x behind PgRC2) is a different regime and is **not**
re-measured here. The parity claim above is scoped to the ~0.5M-read files
tested in this pass.

---

## 5. What this means for the structural rewrite

The rewrite target must be re-aimed. Ruled out by this analysis:

- ~~Add disk staging / memory-release valve~~ — they don't have one in the
  release build, and we already stage more than they do.
- ~~Their position array is structurally cheaper~~ — A/B tested directly in
  `LAYER_BY_LAYER_ANALYSIS.md` §3e: our unique-indexed scheme already beats
  their direct-per-original scheme on 5 of 6 files.
- ~~Coverage regime / duplication rate explains the losses~~ — both disproven
  against real data (§3e).

What survives as the real structural difference, and where the remaining
size gap lives:

**Their 3-way pseudogenome split.** They build *three* separate pseudogenomes
— `hqPg` (high-quality reads assembled by exact overlap), `lqPg` (low-quality
reads mapped onto the hqPg), and `nPg` (N-containing reads) — each with its
own independent reads list and its own set of 10 streams, each compressed and
released separately. We build one pg plus an undifferentiated leftover pool.

This matters because of the measured mismatch-layer finding (§3e): our L5 cost
tracks the fraction of reads falling through to lenient leftover-matching
(E. coli 22.8% leftover → win; H. salinarum 37.6% leftover → −36.8% loss).
Their quality-based division routes reads into the *right* pg up front, so
the reads that would become our expensive fuzzy-matched leftovers are instead
handled by a dedicated pg with its own coder. That is a genuine structural
difference, it is on the path where our measured loss actually lives, and it
has not been ruled out by anything tested so far.

**The no-pre-dedup change is already built and verified** (`103_no_predup.cpp`,
byte-identical on P. aeruginosa / H. salinarum / S. acidocaldarius): H.
salinarum improved 2,588,038 → 2,479,901 (−4.2%, gap −36.8% → −31.1%),
P. aeruginosa flat (+0.2%). Real but partial — consistent with the conclusion
that the remaining gap is in the quality-division/mismatch-routing structure,
not in the dedup or position layers.
