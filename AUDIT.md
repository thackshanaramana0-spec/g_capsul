# AUDIT — the technical freeze and reproducibility audit, 2026-09-10

> **Status:** frozen 2026-09-10 at tag `v1.0.1-capsule` — code and results final.
> Authoritative numbers live in `benchmark/results/`; this file is that verification
> record. Where it and a result file disagree, the result file wins.

The final pass before this project's machine is terminated. Its standard was:
*six months later, with no memory of the project and only this repository plus
the documented datasets, could someone reconstruct exactly what was built, why,
how it was evaluated, and where every published number came from?*

This document is the record of that audit: what was checked, how, what it
found, and what it did **not** establish. It is deliberately blunt about the
last part.

**Audited state: `0f7a145`, tagged `results-final-20260910-audited`.**

---

## 1. The headline

Every published table was checked, and **one number changed as a result —
downward**. But "checked" means two different things depending on the table,
and the difference matters: see the tier note directly below the table. Only
four tables were re-executed against an independent measurement; the rest were
checked for transcription fidelity only.

| table | how it was verified | outcome |
|---|---|---|
| T1.1 / T1.2 archive size, time, RAM | every row re-derived from `benchmark/results/run.log.txt` by an independently written parser, then diffed against the CSV | **57/57 rows**, including every LOSSLESS verdict |
| T2.1 het-SNV | same method | **12/12** |
| T2.2 coverage sweep | same method | **3/3** (30× is carried from T2.1, not re-run) |
| T2.3 het-indel | same method | **8/8** |
| **T2.4 multi-allelic** | **no raw log existed — re-executed 3×** | **21/26 WITHDRAWN → 17/26** |
| T2.5 tetraploid | log evidence was circular — **re-executed** | **all 16 values confirmed exactly** |
| T3.1 export | same method | **6/6** (output bytes, contig count *and* speedup) |
| T3.2 coverage | same method | **6/6** |
| T3.3 query | same method | **19/19** |
| T3.4 locus fidelity | not in the sweep log — **re-executed** | **HG002 row confirmed field for field** |

Raw evidence for every re-execution is committed:
`benchmark/documentation/T2.4_reproduction_20260910/`,
`T2.5_reproduction_20260910/`, `T3.4_reproduction_20260910/`.

### The limit of the "re-derived from the raw log" checks — read this before trusting the table above

The rows marked *re-derived* were checked by parsing `run.log.txt` and diffing
against the CSV. **That is weaker than it sounds, and the weakness is
structural.** In `benchmark_1_run.sh` the log line and the CSV row are written
from the *same shell variables in the same execution*:

    ok     "compress  archive=$(mbs $ARCH) ..."      # -> the log
    printf "%s,CAPSULE,%s,%s,..." ... "$ARCH" ...    # -> the CSV

So a log/CSV match proves **transcription fidelity**, not measurement
correctness. If `stat -c%s` had been pointed at the wrong file, or
`parse_time_v` had misparsed, both artefacts would carry the identical wrong
number and the check would still report a clean match.

The verification therefore has two tiers, and they are not equivalent:

| tier | what it proves | tables |
|---|---|---|
| **Re-executed** — independent measurement | catches wrong *values* | T2.4 (3×), T2.5, T3.4 (HG002 row), **T2.1 + T2.3 (HG002 rows, 12/12 fields)**, T1.1's E. coli and HG002 `archive_bytes` |
| **Transcription-checked** — log vs CSV | catches copying/typing errors only | T2.1/T2.3's **HG003, HG004, HG005** rows; all of T2.2, T3.1, T3.2, T3.3; T3.4's HG003/HG004/HG005 rows; and T1.1/T1.2's timing and RAM columns (size and lossless columns are independently checked — see below) |

**The single real defect this audit found (T2.4) was found by tier 1, and could
not have been found by tier 2.** That is the honest measure of what the two
tiers are worth. Tier 2 is not worthless — it would catch a mis-transcribed
table — but it cannot detect the failure mode that actually occurred here.

**But re-executing everything is the wrong way to close this**, and saying so
was lazy. Most of the gap closes with checks costing seconds, because the
tier-2 columns are not all equally weak:

| T1.1/T1.2 column | what it actually rests on |
|---|---|
| `raw_bytes` | **INDEPENDENT.** Re-checked 2026-09-10 against the real FASTQ files still on disk: **19/19 exact**. Nothing to do with the run's variables. |
| `ratio_pct` | **DERIVED**, and the arithmetic was verified: `100 × archive/raw` reproduces the stored value on **57/57** rows. |
| `lossless` | **INDEPENDENT.** It is a decode-and-compare, not a copied variable. |
| `archive_bytes` | The genuinely weak one — but the file that was `stat`-ed is the same file that then decoded LOSSLESS, so a wrong `stat` needs a *second, coordinated* error to survive. E. coli independently reproduces at 68,429,027 B, and the encoder is proven byte-identical across cores and settings. |

So T1.1/T1.2 is in far better shape than "transcription-checked" implies.

The F1 columns of T2.1/T2.3 had no equivalent cross-check — so they were
**re-executed too, and it took 90 seconds, not 19 hours**: an HG002 archive
from an earlier job had survived on disk, and the caller can be pointed
straight at it. One run verifies both tables, because the caller emits SNV and
INDEL in the same pass.

    T2.1 HG002  TP=37011 FP=1721 FN=7564 P=0.956 R=0.830 F1=0.888   EXACT
    T2.3 HG002  TP=3871  FP=824  FN=3913 P=0.825 R=0.497 F1=0.620   EXACT

**12/12 fields reproduce exactly**, using a different archive file, a freshly
built decoder, and a direct runner invocation — sharing no variable, process or
intermediate with the original. Evidence:
`benchmark/documentation/T2.1_T2.3_reproduction_20260910/`.

**The lesson, since the original claim here was wrong:** "re-execute
everything" was the expensive answer to the wrong question. Most of the gap
closed with checks costing seconds — stat the inputs, verify the arithmetic,
point the caller at a surviving archive. What remains genuinely unverified is
now small and specific, not "seven tables".

---

## 2. What the audit caught

Four real defects, none of which would have been found by reading the code.

### 2.1 T2.4's 21/26 was unreproducible — corrected to 17/26

Walking the claim → raw-log chain found that T2.4's published figure had **no
supporting raw log anywhere on this server**. Rather than assume either number
correct, it was investigated:

1. **Re-derived the denominator independently** by streaming the real GIAB
   HG002 BAM: 952 multi-allelic truth sites, 26 SNV-only. Matches the
   published denominator exactly — *that half of the claim was right*.
2. **Ran the caller end to end** on the full chromosome: **17/26**.
3. **Suspected known caller nondeterminism** (bubble-traversal order can flip
   which duplicate `lift_vcf.py`'s dedup keeps). Re-ran on the *identical*
   already-streamed reads: **17/26 again**, every intermediate count matching
   (38,887 SNVs, 7,022 indels, 151,286 contigs, 48,635 lifted calls).
   Nondeterminism refuted.
4. **Suspected my own earlier encoder fix** (`b28dd40`, the only commit
   touching the encoder in between). Built the encoder at the prior commit via
   `git worktree` and re-ran: **17/26, a third time**. The fix is confirmed
   inert.
5. **Found the original working directory.** Its `regions.bed` spans
   chr20:1,000,000–6,000,000 — a **5 Mb window**, not the chromosome — and its
   own truth file has 111 sites, 7 SNV-only, matching this project's *earlier,
   already-superseded* window figures.

**Conclusion:** the truth-side census (952/26) is real and was computed
correctly; the caller-recovery figure paired with it has no corresponding
execution. 21/26 withdrawn. The verified figure is **17/26 (65.4%)**,
propagated to all 25 files that cited it, each with a dated retraction note.

The *qualitative* claim is unaffected: DiscoSNP++ recovered **0/26** in all
three reproductions, structurally — it emits zero multi-allelic records in
3,989 total.

### 2.2 T3.4 needed an undocumented flag — three false failures before finding it

Reproducing T3.4 gave `coordinate = 0` three times, matching the *old*,
already-acknowledged-wrong figure. The cause was not a defect but a
**required flag documented nowhere**:

    CAPS_PILEUP=1 capsule_decode index <archive> <archive>.qidx

Without it the sidecar carries the pseudogenome and placements but **not each
read's deviations**, `query` emits the *consensus* rather than the reads, no
allele difference can appear at a heterozygous site, and the coordinate arm
scores 0 **by construction rather than by measurement**.

| attempt | sidecar | coordinate arm |
|---|---|---|
| 1–3 | absent, or built without the flag (120,274,641 B) | **0** |
| 4 | built **with** `CAPS_PILEUP=1` (172,045,852 B, carrying 7,466,871 deviations) | **18** ✓ |

`CAPS_PILEUP` appeared in **zero** documentation files. A reviewer following
the documented steps would have hit the same wall and concluded the published
number was fabricated — the same wrong conclusion reached here three times.
Now documented in `REPRODUCE_EVERYTHING.md` and as a REQUIRED SETUP row in
`RESULT_CODE.md`.

### 2.3 `.gitignore` was silently excluding the audit evidence itself

Twice. The blanket `*.log` rule excluded the preserved reproduction logs —
including ones already committed, so a README cited raw evidence **that was
never actually in the repository**. Then the same for `*.vcf`, which excluded
the T2.4 correction's own truth file. Both fixed with targeted exceptions;
`MANIFEST.sha256` extended to cover `.log` and `.vcf` so they are checksummed.

### 2.4 The manifest itself was wrong — found by actually cloning

`MANIFEST.sha256` was first generated with a plain `find` over the working
directory, which silently pulled in a concurrent agent's **untracked** files.
A genuine fresh clone — the exact scenario the manifest exists to support —
failed to verify it with 8 `No such file` errors. Regenerated from
`git ls-files` only.

### 2.5 T2.5's "verification" was circular

The sweep log's end-of-run section prints T2.5 with `column -s, -t "$CSV7"` —
it pretty-prints the CSV it is supposed to corroborate. Caught by reading the
harness rather than trusting the log's appearance. Re-executed instead; all 16
values confirmed exactly.

---

## 3. What was checked and found clean

- **Fresh-clone reproduction.** Cloned from GitHub into a scratch directory,
  built the encoder from that clone alone, ran it under a restricted 2-core
  profile: **byte-identical** archive to this checkout's binary on 12 cores.
- **The canonical checkpoint** (E. coli, `sanity_archive_one.sh`) reproduces
  **68,429,027 B** exactly, 3/3 LOSSLESS — before and after every code change
  in this session.
- **Seeds.** Every script using randomness has a fixed seed; zero unseeded
  uses.
- **Secrets.** Full-pattern scan across all tracked files: clean.
- **Licenses.** `thirdparty/ppmd` had no LICENSE file while its two siblings
  did — added, quoting the public-domain declaration from its own source.
- **Large/stray files.** A 2.4 MB session transcript at the repository root
  moved to `docs/session_transcripts/` with a README stating it is process
  record, not documentation, and must not be cited.
- **Documentation matches implementation.** `paper/` carries no TODO,
  placeholder or "pending" markers; README's documented commands were executed
  and work as written.
- **Tests.** 17/17 self-contained tests; 16/16 scope-and-capability checks.

---

## 4. What this audit did NOT establish

Stated plainly, because an audit that only lists successes is not an audit.

- **No second physical machine.** No container runtime is available in this
  environment and no second host exists. The strongest available substitute
  was used — fresh clone, fresh build, restricted core count — and it passed.
  That is not the same as a genuinely independent machine.
- **No independent human reviewer.** Single session, single party. Every
  finding here was produced by the same process that produced the work being
  audited.
- **T3.4's HG003/HG004/HG005 rows were not re-run.** Only HG002 (100 of 400
  sites). They use the identical script and mechanism and nothing about them
  was found suspect, but they are confirmed *by inheritance*, not by
  execution.
- **The six-month reconstruction test has not been performed by an actual
  outsider** — only simulated by the same party that built the repository.
- **Most tables were transcription-checked, not re-executed** (see §1's tier
  note). The one defect this audit found was in the re-executed group and was
  undetectable by transcription checking. That is not proof the
  transcription-only tables are wrong — it is proof the method used on them
  could not have told us either way.

## 4a. Why "all known gaps are closed" is not "there are no gaps"

Every gap listed in this document has been documented or fixed. That is a
statement about the **known** set, and it should not be read as a statement
about the total.

The evidence that the two differ is in this audit itself:

- T2.4's missing evidence existed for a full day before anyone looked for it,
  and was found only because someone deliberately walked the claim → log chain.
- The `.gitignore` evidence-exclusion bug was found **twice** — the second
  instance (`*.vcf`) was missed on the first pass that fixed `*.log`.
- The manifest bug was found by accident, while doing something else.
- T3.4's missing flag took three failed attempts to identify, and the first two
  produced confident wrong diagnoses.

**Detection rate here is demonstrably below 100%.** One bad number in ten
tables was found; the honest inference is not "there was exactly one" but
"the process that produced one can produce another, and the audit that caught
it is not exhaustive." The highest-value next step is not more documentation —
it is re-executing the tier-2 tables, and a second person looking.

---

## 5. Trail

| commit | what |
|---|---|
| `b28dd40` | three silent out-of-scope data-loss bugs found and fixed |
| `478a241` | scope self-test + `industry/` safety layer |
| `7410a6a` | long-read scope decision documented in the paper |
| `b4c71e6` | licenses, stray files, first checksums, provenance note |
| `2c7941c` | manifest fixed after the fresh-clone test exposed it |
| `7056c4f` | **T2.4 corrected 21/26 → 17/26** |
| `e1850b6` | Claim 3 audited; T2.5 disclosed as unverified |
| `0f7a145` | **T2.5 re-executed and confirmed; every table now has evidence** |

### The canonical tag, and why there are five that say "final"

**The canonical tag is `v1.0.1-capsule`. Use that one.** Everything else below is
history, kept because git tags are immutable and this project does not rewrite
its record.

Five tags were cut on 2026-09-10, each named as if it were the last. They are
not five candidate results — four of them are byte-identical in code and in
results, and differ only in documentation:

| tag | commit | code/results changed vs previous |
|---|---|---|
| `results-final-20260910` | `b4c71e6` | — (first) |
| `results-final-20260910-corrected` | `7056c4f` | **2 files** — the T2.4 21/26 → 17/26 correction |
| `results-final-20260910-audited` | `0f7a145` | 0 |
| `results-final-20260910-audited-v2` | `9e6722d` | 0 |
| `results-final-20260910-locked` | `d85dd67` | 0 |
| `v1.0-capsule` | `2ad4e14` | 0 |
| **`v1.0.1-capsule`** | HEAD | 0 — documentation only |

`v1.0.1-capsule` fixes four documentation defects found by auditing the repo
*structure* after the freeze: the locked-dataset file did not reconcile its 15
against the published tables' 19 (the difference is the four GIAB human sets,
now written out); it still listed two datasets as never benchmarked when both
are in the final sweep at 3/3 LOSSLESS; `RESULTS_INVENTORY.md` said twelve run
directories where there are fourteen, leaving `results/claim2/` and
`results/phase_a/` unlisted — both containing files named to look citable; and
the outer ARCS `DATASET_LOCKED.md` sat at the repo root with the more
authoritative-sounding name and no warning that it does not govern this repo.

Measured with `git diff --name-only <a> <b> -- benchmark/results stages include
scripts`, in commit-date order. **The last substantive change to code or
results was `-corrected`**; everything after it is prose. A reader six months
from now should check out `v1.0.1-capsule` and ignore the rest — the naming was
mine and it was bad, so it is documented rather than quietly re-tagged.

### The deleted backup tag

`backup-pre-purge-20260908` (commit `36aa7c9`) was deleted at the freeze. It
held 1.77 GiB and contained exactly two things not present at HEAD:

- `out.vcf`, 160,737,970 B — a caller output. Every reference to that filename
  in the repo is a **command-line argument in an example**
  (`capsule_decode call <archive> out.vcf`), never a citation of that file's
  contents. Regenerable from the archive.
- a 2.4 MB session transcript, referenced only from `docs/_removable/`.

The 14 planning documents that tag also held were **not** lost: they were
*moved* to `docs/_removable/` and every one was verified byte-identical by blob
hash before the tag was removed. The stale `docs/<PLAN>.md` paths that the move
left behind in ten files were repaired at the same time.

## 6. How to check this yourself

```bash
sha256sum -c benchmark/documentation/MANIFEST.sha256   # file integrity
bash scripts/test_0_scope_and_capability.sh            # scope, proven live
bash scripts/run_tests.sh                              # 17 self-contained tests
```

Then `benchmark/documentation/RESULT_CODE.md` traces every published number to
the script that produced it and the file that holds it.
