# G_CAPSUL vs the field — layer-by-layer, then general, then phase-by-phase

Written 2026-09-02. Three tables, in the order requested: (1) a stream-by-
stream technical comparison against PgRC2 specifically, since it is the only
tool close enough in architecture to compare layer-for-layer; (2) a general
comparison against SPRING and Genozip, which have no equivalent internal
layers to line up against (SPRING is read-relative, Genozip has no assembly
step at all); (3) what specifically changes between Phase 1 / Phase 2b /
Phase 3 and why each tool's number means something different at each phase.

**Status note (updated 2026-09-02, run complete):** the full 14-dataset
Phase1/Phase2b/Phase3 benchmark finished — every number in Table 3 below is
final and round-trip-verified, no PENDING rows remain. Source of truth:
`/tmp/allph/r.csv`. The 15th locked dataset (Utricularia gibba) is still
queued for download and not included.

---

## Table 1 — layer-by-layer vs PgRC2

Both tools solve the same core problem (pseudogenome-based read
compression) with recognizably similar internal stages, so this is the one
comparison that can be made at the level of individual streams rather than
whole-archive totals. Source: `docs/HOW_PGRC2_CODES_REFERENCES.md`,
`docs/claim1_locked_result` (memory), and PgRC2's own source
(`/root/arcs-clean/method_c`, read directly, GPL-3, never vendored).

| layer | PgRC2's approach | G_CAPSUL's approach | who wins / how effective |
|---|---|---|---|
| **Assembly admission** | Both-side-overlap required (`getBothSidesOverlappedReads`) — a read only joins the pseudogenome if it overlaps well on BOTH sides | Single-pass exact suffix-prefix chaining, admission by overlap length only | Tested PgRC2's rule directly (B.2 in failures doc): WORSE on all 3 datasets tested (+2,845 / +45,746 / +78,597 B). G_CAPSUL's simpler rule wins here. |
| **Duplicate handling** | Reads sorted, exact duplicates removed for free before assembly | Originally needed a separate pre-assembly dedup array (`orig2uid`) because the sweep started below read length; fixed by starting the sweep AT read length (commit `3e06957`) so duplicates chain naturally | PgRC2's native approach and G_CAPSUL's fixed approach are now equivalent in effect; G_CAPSUL's fix flipped S. acidocaldarius from a loss to a win (pg 9,134,100 → 6,157,270 B) |
| **3-way pseudogenome split** | Splits into 3 physical regions (main / off-target / N-reads-ish groupings) at the disk-persistence layer | Single pseudogenome + one second region (2 total); N's substituted and corrected inline instead of a 3rd structural split | Measured directly (`docs/DO_WE_NEED_THEIR_3WAY.md`): G_CAPSUL's 2-region design does NOT need the 3rd split to match or beat PgRC2 — the extra split PgRC2 does is not free headroom G_CAPSUL is missing |
| **Reference encoding** | `compressRlMisRevOffDest`: per-period streams (period = mismatch count per read), destination/length/rev-comp/offset each their own stream | `mem_triples` (dst/src/len/rc) split the same way (`mem_dstgap`, `mem_len`, `mem_rc` separate streams) — same insight, arrived at independently and later confirmed to match PgRC2's own log output | Tie in structure; G_CAPSUL measured within ~11% of PgRC2 on raw reference byte count and CHEAPER specifically on the length sub-stream |
| **Mismatch position coding** | Per-period range coder, bucketed by mismatch count | `mmpos_encode_buckets` — same "bucket by count, period-N range coder" idea, PLUS a flat fallback chosen automatically when smaller | G_CAPSUL's positions measured BELOW both order-0 entropy AND count-conditioned entropy — i.e. at or past the theoretical floor for this representation |
| **Mismatch symbol coding** | Independent per-mismatch symbol code | `mmc::encode` — ADAPTIVE, keyed on the reference base (4-way context → 3-symbol "not equal to ref" alphabet) | **G_CAPSUL wins by 40%**: 124,280 B vs PgRC2's 208,234 B on the same measured data — the single largest per-stream margin found in this project |
| **Read-order / permutation** | Stores enough to reconstruct original order; exact mechanism not re-derived here | Delta-coded `orig2uid` split into a 1-bit zero/nonzero flag stream + sparse alias-value stream, because 79.55% of deltas are exactly zero (non-duplicate reads) and no single model fits both the all-zero majority and the sparse ~273K-distinct-value minority well | G_CAPSUL's read order measured at 0.996× its own information-theoretic bound — essentially closed; not compared stream-for-stream against PgRC2's exact scheme since PgRC2's own paper does not break this stream out separately |
| **Strand (RC) flag** | Not separately documented in their released log output at this granularity | 1 bit/read, `pos_strand`, `best_encode` | G_CAPSUL: 16,314 B against a computed order-0 bound of 16,116 B (orders 1–8 tested down to 16,095 B) — within 1.2% of the bound at order-0 already |
| **Position coding** | Not the focus of their public per-stream breakdown | Fixed-width uint32 + byte-plane split + xz, chosen over varint by direct measurement (varint's length-prefix bits break byte alignment xz's LZ77 needs) | E. coli −7.9%, P. aeruginosa −11.1% vs the varint alternative G_CAPSUL itself used to use |
| **Sequence quality/coding config** | Fixed per their `CODER_LEVEL_NORMAL` defaults (`MINMEM`=45 in their param file) | `MINMEM` (24) and `MAXMAP` (49) both independently swept and confirmed as INTERIOR OPTIMA on the one dataset G_CAPSUL loses on (S. acidocaldarius) — i.e. the loss is not a tuning miss | Both tools' defaults are locally optimal for their own designs; G_CAPSUL's loss margin there is genuinely structural (see Table 2), not a missed parameter |
| **Names / read-ID column** | **None at all** — confirmed directly by reading their decoder source: sequence-only output, no ID handling whatsoever | Full column: SPRING-derived tokenizer + adaptive dictionary + `ID_SEQLEN` cross-column reference + line-3 elision (§ Technical Architecture) | Not a comparable axis — PgRC2 has nothing here. G_CAPSUL's names are compared against SPRING and Genozip instead (Table 2) |
| **Quality column** | **None at all** — same confirmation as names | Vendored fqzcomp (BSD), full block+index archive integration | Same as names: not comparable to PgRC2, compared to SPRING/Genozip/real fqzcomp instead |
| **Variable-length read support** | **Fails outright** on 6 of the 14 locked datasets: `Unsupported variable length reads` (SARS-CoV-2, HCMV, C. jejuni) and hard crashes (`*** stack smashing detected ***` on H. pylori, `free(): invalid next size` on A. fumigatus) | Handles all 14, including the ones that crash PgRC2 — after fixing the 4 bugs in `FAILURES_AND_REFUTED_IDEAS.md` that were specifically triggered by variable-length input | **G_CAPSUL wins by capability, not just ratio** — a real, citable advantage: PgRC2 cannot even attempt 6 of the 14 datasets in this benchmark |
| **Whole-archive aggregate** (7 datasets both tools could run, sequence+order only) | 83,192,412 B total | 81,631,156 B total | **G_CAPSUL +1.88% smaller, 6 wins / 1 loss** (S. acidocaldarius −0.83%, the one structural loss) |
| **Speed** | Faster: baseline | ~1.6–1.7× slower than PgRC2 | PgRC2 wins on speed — expected and disclosed; assembly-based lossless compression trades time for ratio, and the paper states this plainly rather than hiding it |
| **Peak RAM** | Lighter: baseline | ~2.5–2.8× heavier at worst | PgRC2 wins on RAM for the same reason as speed |

**Net read on Table 1:** G_CAPSUL wins on ratio (thin but real, +1.88%
aggregate), wins outright on robustness (PgRC2 cannot process 6 of 14 real
datasets), and loses on speed/RAM by a known, disclosed, expected margin.
Within individual streams, the mismatch-symbol coder (+40%) is G_CAPSUL's
single strongest technical result; the one structural loss
(S. acidocaldarius) is confirmed NOT a tuning miss.

---

## Table 2 — general comparison vs SPRING and Genozip (and the field)

SPRING and Genozip have no equivalent internal "layers" to compare
stream-for-stream — SPRING reorders reads and codes them read-relatively
(closer to a de Bruijn/overlap-graph design without an explicit
pseudogenome), and Genozip has NO cross-read matching mechanism at all by
default. So this table compares by COLUMN (sequence, names, quality) and by
CAPABILITY rather than by internal stage.

| axis | SPRING | Genozip | G_CAPSUL |
|---|---|---|---|
| **Core mechanism** | Hash-table seed index (32bp two-substring), Hamming-distance read reordering + read-relative mismatch coding | NO cross-read matching by default — a specialised `acgt` codec + generic LZMA/BSC/bzip2 backend | Explicit pseudogenome assembly (greedy chain + pigeonhole + second region + MEM self-match) |
| **Sequence, Phase 1 (14 datasets)** | 12/14 wins for G_CAPSUL, **−23.41%** aggregate (2 losses: HCMV +3.0%, C. jejuni +4.4%, both PRE the bug-fix numbers — pending re-confirmation) | 14/14 wins for G_CAPSUL, **−78.48%** aggregate | — |
| **Names column** | SPRING HAS a names coder (`id_compression.cpp`, dynamic per-character-class tokenizer, no format knowledge) | Genozip's names coder (`qname.c`) matches a hand-curated library of ~50 known header FORMATS, gives each field its own permanent context, codes numerics as width-minimised binary — architecturally different, not a "smarter coder", just format-aware where SPRING is format-blind | G_CAPSUL: 7/8 wins vs SPRING (**−18.37%**), 6/8 vs Genozip (**−7.35%**) by the CORRECT measurement method (each tool's own name-section bytes — see `docs/REIMPL_NOTES.md` for why differencing gives a wrong, sometimes negative, answer for Genozip specifically) |
| **Quality column** | Uses BSC on assembly-reordered reads (confirmed by reading `reorder_compress_quality_id.cpp` directly — NOT a custom quality context model) | LZMA/BSC/bzip2 backend, no quality-specific model | G_CAPSUL's own coder (stage 92): 7/8 vs SPRING (−2.88%), 8/8 vs Genozip (−5.24%); vendored real fqzcomp inside the archive matches or beats standalone fqzcomp exactly (0/8 for our own coder vs real fqzcomp before vendoring, now closed by vendoring) |
| **Robustness (variable-length, edge cases)** | Runs on all 14 (with appropriate flags) | Runs on all 14 | Runs on all 14 after the 4-bug fix session; PRE-fix, 4 of 14 archives silently failed to round-trip while still reporting plausible sizes |
| **Whole-FASTQ (Phase 3), 7 tested so far** | G_CAPSUL wins all 7 tested (see Table 3) | G_CAPSUL wins all 7 tested (see Table 3) | — |
| **Speed** | Fast (SPRING is a production, speed-conscious tool) | Very fast (streaming architecture, `dispatcher.c`) | Slower than both — expected, same trade-off as vs PgRC2 |
| **Where G_CAPSUL is NOT claiming SOTA** | — | — | Quality column vs REAL fqzcomp specifically (fqzcomp is a specialist CRAM codec, not a whole-FASTQ tool — see `REIMPL_NOTES.md`); X/Y coordinate coding in names is already at its entropy bound on both sides, no further win available there for either party |
| **PgRC2, for reference** | PgRC2 beats SPRING by 15.7% (SE) / 14.0% (PE) per PgRC2's own published paper | PgRC2's paper dismisses Genozip outright ("streaming-regime tools cannot obtain competitive compression ratios on FASTQ reads") | G_CAPSUL beats PgRC2 by a much thinner +1.88% (Table 1) — the honest published-literature context for how large G_CAPSUL's SPRING/Genozip margins actually are |

**Net read on Table 2:** the paper's defensible claim is "competitive
compression against SPRING and Genozip, PLUS an addressable, unified
archive" — not "sequence-compression SOTA" (PgRC2 already holds a larger,
published margin over SPRING than G_CAPSUL does; mstcom is a further,
though fragile, contender behind PgRC2 — see `competitor_landscape_mstcom`
memory). G_CAPSUL's genuinely defensible SOTA-adjacent claims are: (a) beats
PgRC2 itself by a real, if thin, margin while ALSO adding names+quality+
addressability that PgRC2 has none of, and (b) processes real datasets that
crash or are refused by PgRC2.

---

## Table 3 — Phase 1 / Phase 2b / Phase 3, what specifically differs

Each phase adds exactly one FASTQ column to the comparison, and — critically
— changes what "the same file" means for each competitor tool, since SPRING
and Genozip don't expose per-column archives the way G_CAPSUL's container
does. This table exists so a reader knows exactly what each phase's number
does and does not include.

| | **Phase 1** | **Phase 2b** | **Phase 3** |
|---|---|---|---|
| **FASTQ columns covered** | sequence + read order | + names + line 3 | + quality = the WHOLE file |
| **G_CAPSUL invocation** | default (`CAPS_NAMES`/`CAPS_QUAL` unset) | `CAPS_NAMES=1` | `CAPS_NAMES=1 CAPS_QUAL=1` |
| **SPRING invocation** | `--no-ids --no-quality` (isolates seq+order cleanly — a real flag, not a neutralisation hack) | `--no-quality` only (IDs retained) | full, no flags (whole FASTQ) |
| **Genozip invocation / measurement** | its own `--stats` SEQ field only | SEQ + every `Parent=QNAME` context + `length` + `LINE3` (Genozip has no `--no-quality`-equivalent flag, so its number is assembled from `genocat --STATS`'s exact per-context bytes, not a single archive) | the whole `.genozip` archive size, unmodified |
| **PgRC2** | included (only phase it can be compared on — no names/quality stage exists) | excluded (no names stage) | excluded (no quality stage) |
| **What changes G_CAPSUL's OWN number between phases** | nothing extra beyond the base archive | + `names_body`/`names_dict`/`names_index` streams | + `qual_body`/`qual_index` streams |
| **Verification required before a number counts** | decode → compare sequence column to original | decode → compare sequence AND `.names` output to original | decode → rebuild full 4-line FASTQ from decoder output alone → compare to original byte-for-byte (MD5) |
| **Known result pattern so far** | 12/14 vs SPRING, 14/14 vs Genozip, 7/8 vs PgRC2 (PRE bug-fix numbers, re-confirming) | Adding names IMPROVES the SPRING standing (12/14 → 14/14 in the earlier run) because names is where G_CAPSUL's margin over SPRING is largest | 7/7 tested so far beat BOTH SPRING and Genozip on the full file — see interim table below |
| **Datasets affected by the 4-bug-fix session** | ERR552797, SRR40271341 (the >256bp clamp, bug A.1) affected Phase 1 numbers directly | same 2, PLUS SRR40402583/SRR32429602 whose Phase 2b names/line3 columns were fine but whose SEQUENCE column (shared across all 3 phases) was not | all 4 affected, since Phase 3 = Phase 1 + Phase 2b's columns plus quality |

### Phase 3 final results — 14/14 datasets complete, run finished 2026-09-02

All three levels round-trip verified per dataset before being recorded.
**14/14 whole-file (Phase 3) wins vs both SPRING and Genozip.** Source of
truth: `/tmp/allph/r.csv` (all 14 rows `p1ok=p2ok=p3ok=OK`).

| dataset | Phase 1 | Phase 2b | Phase 3 (whole FASTQ) | vs SPRING P3 | vs Genozip P3 |
|---|---|---|---|---|---|
| ERR5181310 (SARS-CoV-2) | 836,191 | 838,696 | 8,686,774 | −7.5% | −5.8% |
| SRR554369 (P. aeruginosa) | 8,981,037 | 8,985,459 | 57,320,645 | −3.1% | −35.5% |
| ERR552797 (M. tuberculosis) | 4,716,998 | 5,983,492 | 46,964,185 | −9.9% | −43.3% |
| SRR2584863 (E. coli) | 8,225,993 | 12,246,611 | 68,677,977 | −7.2% | −42.6% |
| SRR29296997 (H. salinarum) | 2,756,247 | 3,748,338 | 15,654,489 | −10.5% | −42.5% |
| ERR12954017 (S. acidocaldarius) | 3,141,277 | 4,136,647 | 15,159,145 | −10.7% | −61.7% |
| SRR40402583 (C. jejuni) | 3,806,837 | 5,265,313 | 9,040,464 | −5.3% | −70.7% |
| SRR40271341 (H. pylori) | 3,853,956 | 4,693,232 | 40,620,675 | −11.2% | −34.7% |
| ERR17740259 (S. aureus) | 13,522,261 | 20,046,973 | 83,421,422 | −9.6% | −49.2% |
| SRR37283774 (P. falciparum) | 17,161,828 | 24,135,957 | 67,382,606 | −5.3% | −29.1% |
| DRR976266 (S. cerevisiae) | 21,822,808 | 30,781,796 | 56,684,465 | −8.2% | −71.8% |
| SRR36741279 (L. major) | 27,954,749 | 39,252,551 | 106,342,010 | −9.5% | −45.4% |
| SRR32429602 (HCMV) | 27,444,199 | 34,701,936 | 55,694,667 | −2.9% | −47.9% |
| SRR39257532 (A. fumigatus) | 37,904,777 | 54,463,895 | 108,558,425 | −29.3% | −52.2% |

Utricularia gibba (`SRR10676752`, the 15th locked dataset) is still queued
for download, not included here or in the 14-dataset count above.

---

## How to regenerate these tables

```bash
bash /tmp/allphases.sh    # NOT committed to scripts/ -- a real, named gap,
                          # see docs/INDUSTRIAL_CHECKLIST_CLAIM1.md
cat results/phase_a/allphases_14dataset.csv   # rescued from /tmp/allph/r.csv 2026-09-03
```

Column order in that CSV: `dataset,raw,p1,p1ok,p2,p2ok,p3,p3ok,spr1,spr2,spr3,gz1,gz2,gz3`.
A phase's `pN` value is only meaningful where the paired `pNok` column reads
`OK` — this is the round-trip gate described in Table 3.

---

## Table 4 — Claim 2 competitors (added 2026-09-03)

Not a compression comparison — DiscoSNP++ and Kmer2SNP are the only two
tools identified as applicable to this exact task (reference-free
heterozygous variant calling from a single diploid sample; see
`docs/HET_INDEL_SOTA.md` §2 for the literature survey ruling out
DeepVariant/GATK/Clair3/etc. as reference-based, and `ska lo`/eBWT2SNP as
solving a different problem). Full numbers and per-window breakdown live in
`docs/CLAIM2_TABLES_AND_INDEL_SCAN.md` — this table is the summary anchor.

| tool | het-SNV F1 | het-indel F1 | multi-allelic sites found | notes |
|---|---|---|---|---|
| **G_CAPSUL** | **0.890** | 0.637 | **11/18** | wins SNV + multi-allelic, loses indel |
| DiscoSNP++ v2.6.2-12 | 0.874 | **0.663** | 0/18 | strongest applicable competitor; cannot emit true multi-allelic records |
| Kmer2SNP | 0.464 | not applicable (SNP-only by construction) | not applicable | weakest of the three; real compatibility fixes required to even run it, see `docs/KMER2SNP_BENCHMARK.md` |

## Table 5 — Claim 3 competitors (added 2026-09-03)

Addressability has no direct single-tool competitor — it is compared
per-operation against the conventional pipeline that would otherwise
compute the same thing. Full numbers in `docs/CLAIM3_LOCKED.md` §6;
per-operation prior-art position (why none of these tools offer the same
*combination* of operations) in `docs/CLAIM3_LOCKED.md` §4.

| operation | conventional tool | G_CAPSUL speedup | spec target |
|---|---|---|---|
| export | SPAdes v4.0.0 (spec-exact) | **555-656×** | ≥40× |
| export | MEGAHIT 1.2.9 (earlier substitute, kept on record) | 254-314× | ≥40× |
| coverage | bwa 0.7.17 + samtools + mosdepth 0.3.6 | 23-33× | 2-5× |
| query | full decompression (same binary) | 1.62× (corrected — an earlier 3.3x figure was measured against the wrong baseline, see `CLAIM3_LOCKED.md` §6.4) | not spec'd; real advantage is 132× output selectivity, not time |
