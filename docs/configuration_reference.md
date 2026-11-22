# Configuration reference: environment variables and positional arguments

This page lists every switch that the G_CAPSUL programs read, and says which ones a normal run needs.
It was written from the source. Part B is generated from the source text itself (file, line and the
expression that reads the variable), so it cannot drift from what the code does. Re-generate it if
the code changes.

## How switches work

* Almost every switch is **presence-based**: the code tests `getenv("NAME") != nullptr`. Setting the
  variable to any value, including `0` or an empty string, turns the feature on. To turn it off,
  `unset NAME`. A few switches read a number (`atoi`, `atof`) and are marked "value" below.
  `NODEDUP`, `VERIFY_DUMP` and `FWD_SELF` are the exceptions that read the value as a flag (see below).
* A switch read at program start (`static const bool ... = getenv(...)`) is fixed for the whole run.
* Switches change what is written or how fast it is produced. The archive format (version 2) is the
  same. A decoder never needs the encoder's switches, with the exceptions of the query-time switches
  in A.3.

## A. Switches a normal run uses

### A.1 Compression (`best106`, `src/encoder.cpp`)

The command line is `best106 <in.fq> [argv[2] .. argv[11]]`. The positional arguments are documented
in the table above `main()` in `src/encoder.cpp` and are repeated here (all optional, defaults from
the source):

| arg | name | default | meaning |
|---|---|---|---|
| argv[2] | MAXMM | 3 | parsed and **not used** by any code |
| argv[3] | MINOV | 16 | shortest suffix-prefix overlap the builder pass chains |
| argv[4] | SW | 32 | seed width of the chaining sweep, at most 32 |
| argv[5] | SEEDSTRIDE | 8 | spacing of the seeds taken from each read for the mapping index |
| argv[6] | SEEDW | 32 | mapping seed width, clamped to 8..32 |
| argv[7] | R1MINOV | SW | shortest overlap accepted in the first chaining round |
| argv[8] | K2 | 1 | stride when sweeping the pseudogenome for seeds |
| argv[9] | MINMEM | 45 | minimum length of a MEM reference |
| argv[10] | MAXCAND | 64 | candidate matches examined per query position |
| argv[11] | LAZY | 1 | nonzero turns on lazy matching in the MEM parse |

The values in the README command (`3 16 16 22 16 16 1 24 64 1`) set argv[2] to argv[11] in that order.

| switch | kind | effect |
|---|---|---|
| `DUMP_PERM`, `DUMP_MM` | presence | **Required for a usable archive.** They gate the code that fills the placement, read-length, order and mismatch streams. The names are historical. `encode_adaptive.sh` and the Dockerfile set them. |
| `VERIFY_DUMP` | flag value | `VERIFY_DUMP=1` also writes the raw stream files (`literal.txt`, `mem_triples.bin`, `pos_abs.bin` and so on) into the working directory, for `scripts/build_scripts/verify_lossless.sh`. |
| `DUMP_LIT` | presence | Set by the scripts, never read by the encoder. |
| `CAPS_NAMES` | presence | Also code the read names (`names_body`, `names_dict`, `names_index`). |
| `CAPS_QUAL` | presence | Also code the quality values (`qual_body`, `qual_index`). |
| `CAPS_SPANS` | presence | Record the contig spans (`contig_spans`). Needed to call variants from the archive later. |
| `CAPS_CALL` | presence | Call variants at compress time. Implies `CAPS_SPANS`. About 20 times the cost of recording spans only. |
| `CALL_VCF` | value | Output path of the compress-time VCF. Default `out.vcf`. |
| `ARCHIVE` | value | Output archive path. Default `out.arcs2`. |
| `SEQT` | value | Number of chunks of the `literal` coder. Default 4. This is a compression parameter (see `include/seqpar_core.h`). |
| `NODEDUP` | flag value | `NODEDUP=1` skips the pre-assembly duplicate collapse. `NODEDUP=0` keeps it on. Unset, the encoder measures the duplicate fraction and skips the collapse when it is below 0.15. |
| `MAXMAP`, `MAXMAP_DIV` | value | `MAXMAP` forces the mapping mismatch ceiling. Otherwise it is Lmax divided by `MAXMAP_DIV` (default 13), at least 6, and it is raised only when the leftover fraction is above 0.60 (a low-coverage ramp, see the comment at the `MAXMAP` definition). |
| `MINOV` | value | Overrides the sweep's minimum overlap, for sweeps. |
| `LOWCOV` | presence | Low-coverage mode: sweep seed width 8 and minimum overlap at most 8. Never selected automatically. |
| `MEM_MAXMM_OVERRIDE` | value | Overrides the MEM extension-mismatch tolerance, which is otherwise derived from the leftover fraction. |
| `SECOND_SELF` | value | 1 or 2: also self-match the second region (forward only, or forward and reverse complement). Off when unset. |
| `FWD_SELF` | flag value | Forward self-matching in the MEM stage. Default on. `FWD_SELF=0` turns it off. |
| `CAPS_NO_QOVERLAP` | presence | Do not overlap the names and quality coding with the main pipeline. |
| `CAPS_QTHREADS` | value | Threads for the overlapped names and quality coding. |
| `CAPS_QBLOCK_MB` | value | Quality block size in MB. The default is `QBLOCK_BYTES` in `quality_coder.h`. |
| `CAPS_QSTRAT`, `CAPS_QOFF` | value | Limit the quality trial grid to the first N strategies (1..4) or offsets (1..2). For measurement. |
| `PROBE_FRAC` | value | Fraction of a stream used to choose its coder. Default 0.40. |
| `CHUNK_PAR` | value | Parallel chunks in the chunked coder. Default 2. |
| `CANDIDATES` | value | A comma-separated list of `MAXMAP:MINOV` pairs (for example `12:16,12:20`). Each pair is compressed as a candidate and the smallest archive is kept. |
| `CAPS_PHASE` | presence | Print the wall-clock time of each phase. |

### A.2 Index building (`capsule_decode index`, `src/decoder.cpp`)

| switch | kind | effect |
|---|---|---|
| `CAPS_FMI` | presence | Also build the FM-index `<archive>.qidx.fmi` over the deviating reads. Implies `CAPS_PILEUP`. |
| `CAPS_FMI_SAMPLE` | value | Suffix-array sampling interval of the FM-index. Default 256. |
| `CAPS_PILEUP` | presence | Also write `<archive>.sites`, the per-position base counts. |
| `CAPS_PILEUP_MIN` | value | Keep a (position, base) pair only if at least this many reads show it. Default 3. |
| `CAPS_XMI` | presence | Build the older k-mer completion index `.xmi`. Kept for comparison. |
| `CAPS_XMI_STRIDE` | value | Stride of the `.xmi` k-mers. Default 8. |

### A.3 Query (`capsule_decode query`, `src/decoder.cpp`)

| switch | kind | effect |
|---|---|---|
| `CAPS_QIDX` | value | Path of the decoded sidecar. Default `<archive>.qidx`. A nonexistent path gives a query with no index. |
| `CAPS_QUERY_CONTAIN` | presence | Return every read that contains the probe. Uses `<sidecar>.fmi` when present. |
| `CAPS_QUERY_INDEX` | value | `xmi` selects the older completion index. `fmi` makes a missing `.fmi` an error instead of a fallback. |
| `CAPS_QUERY_MM` | value | Mismatches tolerated by the probe seeding. Default 0, which is exact. |
| `CAPS_QUERY_MINSEED` | value | Minimum seed length, at least 4. Default is derived from the sidecar size. |
| `CAPS_SKIP_QUAL` | presence | Do not decode quality values. |

### A.4 Variant calling (`capsule_decode call`, `include/caps_caller.h`)

| switch | kind | effect |
|---|---|---|
| `CAPS_CALL_INDELS` | presence | Also call indels (read in `src/decoder.cpp`). |
| `CAPS_CALL_NOQUAL` | presence | Skip quality decoding. Sets `CAPS_SKIP_QUAL`. |
| `CAPS_PLOIDY` | value | 2 to 4. Default 2. |
| `CAPS_MAXRAM_MB` | value | Memory ceiling of the k-mer counter. Default is 5000 MB, or two thirds of the available memory if that is lower. |
| `CAPS_HAPLOID_COV` | value | Overrides the caller's haploid coverage estimate H, a positive integer. |

## B. Complete inventory (generated)

Every name read by `getenv` in `src/` and `include/`. The last column is the source text of the
first reading site. Switches that are not listed in Part A are experiments, ablations (the `NO_`
names turn a stage off) and debugging aids. They are **not part of the supported interface**, and
their effect is defined by the code at the line shown. Where a default appears in the expression
it is the number after the `:`.

| name | first site | source expression |
|---|---|---|
| `ARCHIVE` | `src/encoder.cpp:1896` (+3) | `const std::string base = getenv("ARCHIVE") ? getenv("ARCHIVE") : "out.arcs2"; bool child=false; size_t mine=0;` |
| `BOTHSIDE` | `src/encoder.cpp:2302` | `const bool BOTHSIDE = getenv("BOTHSIDE") != nullptr; for(uint32_t i=0;i<n;++i){` |
| `CALL_VCF` | `src/encoder.cpp:2217` (+5) | `const char* cve = getenv("CALL_VCF"); const std::string cv = cve ? pabso(cve) : (std::string(pcwd) + "/out.vcf");` |
| `CANDIDATES` | `src/encoder.cpp:1846` | `if(const char* cl = getenv("CANDIDATES")){ std::vector<std::pair<uint32_t,uint32_t>> cands;   // (MAXMAP, MINOV)` |
| `CAPSULE_EXPORT_CONTIGS` | `src/decoder.cpp:1068` (+1) | `} else if(out_contigs && !getenv("CAPSULE_EXPORT_CONTIGS")){ // TWO-RECORD FORM, for the graph (CAPS_DBG_ONLY) path.` |
| `CAPS_ALIGNPAIR_CREDIT` | `include/caps_caller.h:4923` | `if (const char* e = std::getenv("CAPS_ALIGNPAIR_CREDIT")) CREDIT = atoi(e); a3.anchors += std::min(pv.second.n, CREDIT);` |
| `CAPS_ALIGNPAIR_MAXEV` | `include/caps_caller.h:4903` | `if (const char* e = std::getenv("CAPS_ALIGNPAIR_MAXEV")) MAXEV = atoi(e); size_t added = 0;` |
| `CAPS_BUBBLE_UNCOLLAPSED` | `include/caps_caller.h:4347` | `const bool BUB_UNCOL = std::getenv("CAPS_BUBBLE_UNCOLLAPSED") != nullptr; // Swept on the tuning window only: INDEL F1 0.315 (0.65) / 0.413 (0.80) /` |
| `CAPS_BUB_DUP_FRAC` | `include/caps_caller.h:4354` | `if (const char* e = std::getenv("CAPS_BUB_DUP_FRAC")) bub_dup = atof(e); CallData cd_bub;` |
| `CAPS_BUB_EXT` | `include/caps_caller.h:4788` | `if (const char* e = std::getenv("CAPS_BUB_EXT"))     EXT = atoi(e); if (const char* e = std::getenv("CAPS_BUB_EXT_TOL")) EXT_TOL = atoi(e);` |
| `CAPS_BUB_EXT_TOL` | `include/caps_caller.h:4789` | `if (const char* e = std::getenv("CAPS_BUB_EXT_TOL")) EXT_TOL = atoi(e); if (EXT > 0) {` |
| `CAPS_BUCKET_HIST` | `src/encoder.cpp:2759` | `if(getenv("CAPS_BUCKET_HIST")){ const uint32_t CAPS[]={8,16,32,64,128,256,1024,4096};` |
| `CAPS_CALL` | `src/encoder.cpp:418` | `static const bool CAPS_CALL  = getenv("CAPS_CALL")  != nullptr; // ── SPANS WITHOUT CALLING ───────────────────────────────────────────────────` |
| `CAPS_CALL_INDELS` | `src/decoder.cpp:2249` | `const bool WANT_INDELS = std::getenv("CAPS_CALL_INDELS") != nullptr; if (!WANT_INDELS) {` |
| `CAPS_CALL_NOQUAL` | `src/decoder.cpp:2286` | `if(getenv("CAPS_CALL_NOQUAL")) setenv("CAPS_SKIP_QUAL","1",1); // Set BEFORE step 1 now, because step 1 is what hands the contigs back and` |
| `CAPS_CAND_PAR` | `src/encoder.cpp:2077` | `if(const char* cc = getenv("CAPS_CAND_PAR")) { long q=atol(cc); if(q>0) K=(size_t)q; } fprintf(stderr,"  [a3] candidates=%zu concurrency=%zu "` |
| `CAPS_CAND_SEQ` | `src/encoder.cpp:2044` | `if(getenv("CAPS_CAND_SEQ")) K = 1; else {` |
| `CAPS_CLOSE_TOL` | `include/caps_caller.h:4740` | `if (!same && std::getenv("CAPS_CLOSE_TOL")) { int mm = 0;` |
| `CAPS_CLOSE_UNIQ` | `include/caps_caller.h:4759` | `if (std::getenv("CAPS_CLOSE_UNIQ")) { uint64_t cv;` |
| `CAPS_DBG` | `src/encoder.cpp:3323` (+4) | `const bool QUAL_BITMAP = getenv("CAPS_DBG_ONLY") && getenv("CAPS_DBG"); const bool SEQ_PACK    = QUAL_BITMAP;   // same gate: graph-only caller` |
| `CAPS_DBG_BAL` | `include/caps_caller.h:3799` (+1) | `const double BAL = std::getenv("CAPS_DBG_BAL") ? atof(std::getenv("CAPS_DBG_BAL")) : 0.0;` |
| `CAPS_DBG_COHC` | `include/caps_caller.h:3086` (+1) | `const int COHC = std::getenv("CAPS_DBG_COHC") ? atoi(std::getenv("CAPS_DBG_COHC"))` |
| `CAPS_DBG_COVCAP` | `include/caps_caller.h:3816` (+1) | `const double CAPMULT = std::getenv("CAPS_DBG_COVCAP") ? atof(std::getenv("CAPS_DBG_COVCAP"))` |
| `CAPS_DBG_DEEPFIRST` | `include/caps_caller.h:2880` | `const bool SHORTEST_FIRST = std::getenv("CAPS_DBG_DEEPFIRST") == nullptr; int bestN = 0; size_t bestDepth = SIZE_MAX; uint64_t bestNode = 0;` |
| `CAPS_DBG_DUMPKC` | `include/caps_caller.h:3678` | `if (const char* kd = std::getenv("CAPS_DBG_DUMPKC")) { if (FILE* kf = fopen(kd, "w")) {` |
| `CAPS_DBG_EMITMAXINDEL` | `include/caps_caller.h:3913` (+1) | `const int EMAXI = std::getenv("CAPS_DBG_EMITMAXINDEL") ? atoi(std::getenv("CAPS_DBG_EMITMAXINDEL")) : 15;` |
| `CAPS_DBG_EXCL` | `include/caps_caller.h:3229` | `const bool EXCL = std::getenv("CAPS_DBG_EXCL") != nullptr; std::vector<uint32_t> ambig(dbg_bubbles.size(), 0);` |
| `CAPS_DBG_FA` | `include/caps_caller.h:3697` | `if (const char* bf = std::getenv("CAPS_DBG_FA")) { if (FILE* f = fopen(bf, "w")) {` |
| `CAPS_DBG_FLANK` | `include/caps_caller.h:2579` (+1) | `const int MAXFLANK = std::getenv("CAPS_DBG_FLANK") ? atoi(std::getenv("CAPS_DBG_FLANK")) : 1000;` |
| `CAPS_DBG_FORCE` | `include/caps_caller.h:2364` | `if (std::getenv("CAPS_DBG") && (looks_diploid \|\| std::getenv("CAPS_DBG_FORCE"))) { auto t_dbg = clk::now();` |
| `CAPS_DBG_HALFW` | `include/caps_caller.h:3118` (+1) | `const int HALFW = std::getenv("CAPS_DBG_HALFW") ? atoi(std::getenv("CAPS_DBG_HALFW")) : 31;` |
| `CAPS_DBG_IBFS` | `include/caps_caller.h:2752` | `const bool WANT_IBFS = std::getenv("CAPS_DBG_IBFS") != nullptr; const int IBFS_MAXINDEL = std::getenv("CAPS_DBG_IBFSMAX")` |
| `CAPS_DBG_IBFSMAX` | `include/caps_caller.h:2753` (+1) | `const int IBFS_MAXINDEL = std::getenv("CAPS_DBG_IBFSMAX") ? atoi(std::getenv("CAPS_DBG_IBFSMAX")) : 100;` |
| `CAPS_DBG_INDEL` | `include/caps_caller.h:3656` (+1) | `if (!clean_snv && !std::getenv("CAPS_DBG_INDEL")) { ++n_ind_drop; continue; } (void)is_indel;` |
| `CAPS_DBG_MAXAMB` | `include/caps_caller.h:3573` (+1) | `const int MAXAMB = std::getenv("CAPS_DBG_MAXAMB") ? atoi(std::getenv("CAPS_DBG_MAXAMB")) : 20;` |
| `CAPS_DBG_MAXEXT` | `include/caps_caller.h:2381` (+1) | `const int MAXEXT = std::getenv("CAPS_DBG_MAXEXT") ? atoi(std::getenv("CAPS_DBG_MAXEXT")) : 60;` |
| `CAPS_DBG_MAXPOLY` | `include/caps_caller.h:2476` (+1) | `const int MAXPOLY = std::getenv("CAPS_DBG_MAXPOLY") ? atoi(std::getenv("CAPS_DBG_MAXPOLY")) : 1;` |
| `CAPS_DBG_MINC` | `include/caps_caller.h:2378` (+1) | `const uint32_t MINC = std::getenv("CAPS_DBG_MINC") ? (uint32_t)atoi(std::getenv("CAPS_DBG_MINC"))` |
| `CAPS_DBG_MINQ` | `src/encoder.cpp:3325` (+5) | `const int  QMIN = getenv("CAPS_DBG_MINQ") ? atoi(getenv("CAPS_DBG_MINQ")) : 20; {` |
| `CAPS_DBG_MINSUP` | `include/caps_caller.h:3797` (+1) | `const uint32_t MINSUP = std::getenv("CAPS_DBG_MINSUP") ? (uint32_t)atoi(std::getenv("CAPS_DBG_MINSUP")) : 0u;` |
| `CAPS_DBG_NOAMB` | `include/caps_caller.h:3572` | `if (is_indel && !std::getenv("CAPS_DBG_NOAMB")) { const int MAXAMB = std::getenv("CAPS_DBG_MAXAMB")` |
| `CAPS_DBG_NOCHAIN` | `include/caps_caller.h:3392` | `if (WANT_SB && !sb_sites.empty() && !std::getenv("CAPS_DBG_NOCHAIN")) { std::unordered_map<uint64_t, size_t> by_entry;` |
| `CAPS_DBG_NOCOH` | `include/caps_caller.h:3077` | `if (!std::getenv("CAPS_DBG_NOCOH") && !dbg_bubbles.empty()) { auto t_coh = clk::now();` |
| `CAPS_DBG_NOMULTI` | `include/caps_caller.h:3458` | `if (st.alleles.size() > 2 && !std::getenv("CAPS_DBG_NOMULTI")) { // distinct allele sequences, order preserved` |
| `CAPS_DBG_NOSHAPE` | `include/caps_caller.h:3562` | `if (is_indel && !std::getenv("CAPS_DBG_NOSHAPE")) { const std::string& sl = (s0.size() > s1.size()) ? s0 : s1;` |
| `CAPS_DBG_NOVGROUP` | `include/caps_caller.h:2862` | `if (ns > 2 && !std::getenv("CAPS_DBG_NOVGROUP")) { std::unordered_map<uint64_t, std::pair<int,size_t>> reach; // node -> (count, earliest depth sum)` |
| `CAPS_DBG_ONLY` | `src/encoder.cpp:3323` (+5) | `const bool QUAL_BITMAP = getenv("CAPS_DBG_ONLY") && getenv("CAPS_DBG"); const bool SEQ_PACK    = QUAL_BITMAP;   // same gate: graph-only caller` |
| `CAPS_DBG_PGANCHOR` | `include/caps_caller.h:3962` | `if (std::getenv("CAPS_DBG_ONLY") && std::getenv("CAPS_DBG_PGANCHOR")) { auto t_pa = clk::now();` |
| `CAPS_DBG_PLACEMODE` | `include/caps_caller.h:3336` (+1) | `const int PLACEMODE = std::getenv("CAPS_DBG_PLACEMODE") ? atoi(std::getenv("CAPS_DBG_PLACEMODE")) : 0;` |
| `CAPS_DBG_SB` | `include/caps_caller.h:2692` | `const bool WANT_SB = std::getenv("CAPS_DBG_SB") != nullptr; const int SB_MAXPATH = std::getenv("CAPS_DBG_SBLEN")` |
| `CAPS_DBG_SBBUDGET` | `include/caps_caller.h:2696` (+1) | `const int SB_BUDGET = std::getenv("CAPS_DBG_SBBUDGET") ? atoi(std::getenv("CAPS_DBG_SBBUDGET")) : 200;` |
| `CAPS_DBG_SBLEN` | `include/caps_caller.h:2693` (+1) | `const int SB_MAXPATH = std::getenv("CAPS_DBG_SBLEN") ? atoi(std::getenv("CAPS_DBG_SBLEN")) : 60;` |
| `CAPS_DBG_SBMAXINDEL` | `include/caps_caller.h:3607` (+1) | `const int SBMAXI = std::getenv("CAPS_DBG_SBMAXINDEL") ? atoi(std::getenv("CAPS_DBG_SBMAXINDEL")) : 15;` |
| `CAPS_DBG_SB_KEEPSNV` | `include/caps_caller.h:3587` | `if (!is_indel && !std::getenv("CAPS_DBG_SB_KEEPSNV")) continue; // ── INDEL SIZE BOUND, derived from the TP distribution ──` |
| `CAPS_DBG_SB_REPLACE` | `include/caps_caller.h:2952` | `if (std::getenv("CAPS_DBG_SB_REPLACE")) continue; }` |
| `CAPS_DBG_STR` | `include/caps_caller.h:3629` | `if (std::getenv("CAPS_DBG_STR")) { const std::string& lng = (s0.size() > s1.size()) ? s0 : s1;` |
| `CAPS_DBG_TIPLEN` | `include/caps_caller.h:2598` (+1) | `const int TIPLEN = std::getenv("CAPS_DBG_TIPLEN") ? atoi(std::getenv("CAPS_DBG_TIPLEN")) : 3 * 31;` |
| `CAPS_DBG_WALKBUDGET` | `include/caps_caller.h:2474` (+1) | `const int WALK_BUDGET = std::getenv("CAPS_DBG_WALKBUDGET") ? atoi(std::getenv("CAPS_DBG_WALKBUDGET")) : 0;` |
| `CAPS_DUMP_CONTIGS` | `src/encoder.cpp:2223` (+3) | `if(const char* dce = getenv("CAPS_DUMP_CONTIGS")){ const std::string dc = pabso(dce);` |
| `CAPS_DUMP_MCONTIGS` | `src/encoder.cpp:4532` | `if(const char* mcp=getenv("CAPS_DUMP_MCONTIGS")){ if(FILE* mf=fopen(mcp,"w")){` |
| `CAPS_DUMP_REFS` | `src/decoder.cpp:981` (+1) | `if(getenv("CAPS_DUMP_REFS")){ FILE* rf=fopen(getenv("CAPS_DUMP_REFS"),"wb");` |
| `CAPS_DUP_FRAC` | `include/caps_caller.h:818` | `if (const char* e = std::getenv("CAPS_DUP_FRAC")) dup = atof(e); if (dup_override > 0) dup = dup_override;` |
| `CAPS_FMI` | `src/decoder.cpp:1159` (+1) | `(getenv("CAPS_PILEUP")!=nullptr \|\| getenv("CAPS_FMI")!=nullptr); if(mode=="index"){` |
| `CAPS_FMI_SAMPLE` | `src/decoder.cpp:1988` (+1) | `const uint32_t D = getenv("CAPS_FMI_SAMPLE") ? (uint32_t)atoi(getenv("CAPS_FMI_SAMPLE")) : 256u; std::string err;` |
| `CAPS_FORCE_READSUB` | `include/caps_caller.h:1276` | `const bool WANT_RSUB = std::getenv("CAPS_FORCE_READSUB") != nullptr; if (WANT_RSUB) {` |
| `CAPS_GAPSCAN` | `include/caps_caller.h:998` (+1) | `if (std::getenv("CAPS_GAPSCAN")) { idx.reserve(1u << 21);` |
| `CAPS_HAPLOID_COV` | `include/caps_caller.h:2267` | `if (const char* he = std::getenv("CAPS_HAPLOID_COV")) { int v = atoi(he); if (v > 0) H = (uint32_t)v; } ` |
| `CAPS_HPDBG` | `include/caps_caller.h:289` (+2) | `if (std::getenv("CAPS_HPDBG") && back >= 3) fprintf(stderr, "[hp] run='%c' back=%zu ea=%zu eb=%zu diff=%ld\n", rc0, back, ea, eb, (long)ea-(long)eb);` |
| `CAPS_INDEL_ANCH` | `include/caps_caller.h:4997` | `if (const char* e = std::getenv("CAPS_INDEL_ANCH")) MIN_ANCH = atoi(e); // LENGTH-SCALED EVIDENCE. Pooled over all 8 evaluations (TP=253,` |
| `CAPS_INDEL_SERIAL` | `include/caps_caller.h:4613` | `const bool IND_PAR = std::getenv("CAPS_INDEL_SERIAL") == nullptr; int n_thr = 1;` |
| `CAPS_INDEL_SUP` | `include/caps_caller.h:5096` | `if (const char* e = std::getenv("CAPS_INDEL_SUP")) MINSUP = atoi(e); if ((int)best_sup < MINSUP) {` |
| `CAPS_JUNC_FL` | `include/caps_caller.h:5056` | `if (const char* e = std::getenv("CAPS_JUNC_FL")) FL = atoi(e); // ADAPTIVE CONTEXT. Fixing the half-length and REJECTING any` |
| `CAPS_JUNC_MAX` | `include/caps_caller.h:5083` | `const bool COH = std::getenv("CAPS_JUNC_MAX") == nullptr; uint32_t best_sup = COH ? UINT32_MAX : 0;` |
| `CAPS_KC_FREQMIN` | `include/caps_caller.h:1574` | `const bool FREQMIN = SUPERK && std::getenv("CAPS_KC_FREQMIN") != nullptr; std::vector<uint32_t> mm_freq;` |
| `CAPS_KC_FREQSTRIDE` | `include/caps_caller.h:1581` (+1) | `const size_t STRIDE = std::getenv("CAPS_KC_FREQSTRIDE") ? (size_t)atoll(std::getenv("CAPS_KC_FREQSTRIDE")) : 32;` |
| `CAPS_KC_NOSPILL` | `include/caps_caller.h:1527` | `bool spill_forced_off = std::getenv("CAPS_KC_NOSPILL") != nullptr; if (!SPILLDIR && !spill_forced_off && MEM_CEIL_MB && PROJ_MB > MEM_CEIL_MB) {` |
| `CAPS_KC_NOSUPERKMER` | `include/caps_caller.h:1559` | `const bool  SUPERK = SPILL && std::getenv("CAPS_KC_NOSUPERKMER") == nullptr; ` |
| `CAPS_KC_SPILLBITS` | `include/caps_caller.h:1547` | `if (const char* sb = std::getenv("CAPS_KC_SPILLBITS")) { int v = atoi(sb); if (v >= 1 && v <= 20) SPILL_BITS = v;` |
| `CAPS_KC_SPILLDIR` | `include/caps_caller.h:1522` | `const char* SPILLDIR_ENV = std::getenv("CAPS_KC_SPILLDIR"); // Auto-enable when the projection exceeds the ceiling and a target` |
| `CAPS_KFREQ` | `include/caps_caller.h:5117` | `if (std::getenv("CAPS_KFREQ")) { double kcap = 1.5 * (double)H;` |
| `CAPS_KFREQ_CAP` | `include/caps_caller.h:5119` | `if (const char* e = std::getenv("CAPS_KFREQ_CAP")) kcap = atof(e) * (double)H; if ((double)best_sup > kcap) continue;` |
| `CAPS_LINKSCAN` | `include/caps_caller.h:5663` | `if (std::getenv("CAPS_LINKSCAN")) { const int TAIL_MIN = 20;          // bases of tail that must re-match` |
| `CAPS_LINKSCAN_MIN` | `include/caps_caller.h:5741` | `if (const char* e = std::getenv("CAPS_LINKSCAN_MIN")) LMIN = atoi(e); size_t n_ls = 0;` |
| `CAPS_LONG_ANCH` | `include/caps_caller.h:5016` | `if (const char* e = std::getenv("CAPS_LONG_ANCH")) LONGMULT = atoi(e); const int need = MIN_ANCH + LONGMULT * std::max(0, len - 2);` |
| `CAPS_LONG_EXT` | `include/caps_caller.h:4786` | `if (const char* e = std::getenv("CAPS_LONG_EXT")) EXT = atoi(e); }` |
| `CAPS_MAPDBG` | `src/encoder.cpp:1430` | `static const bool MAPDBG = (getenv("CAPS_MAPDBG")!=nullptr); #endif` |
| `CAPS_MATCHING` | `include/caps_caller.h:4957` | `const bool MATCHING = std::getenv("CAPS_MATCHING") != nullptr; if (MATCHING) {` |
| `CAPS_MAXOCC` | `include/caps_caller.h:4539` | `if (const char* e = std::getenv("CAPS_MAXOCC")) { int v = atoi(e); if (v >= 1 && v <= 64) MAXOCC = v;` |
| `CAPS_MAXRAM_MB` | `include/caps_caller.h:1500` (+1) | `const size_t MEM_CEIL_MB  = std::getenv("CAPS_MAXRAM_MB") ? (size_t)atoll(std::getenv("CAPS_MAXRAM_MB"))` |
| `CAPS_MBDBG` | `include/caps_caller.h:4677` (+1) | `if (std::getenv("CAPS_MBDBG")) fprintf(stderr,"[mb] enter ok=%d apos=%u rp=%u Rsz=%zu Asz=%zu\n",(int)b2.ok,b2.apos,rp2,R.size(),Aalt.size());` |
| `CAPS_MM_VAF` | `src/encoder.cpp:4412` | `const bool MMVAF = getenv("CAPS_MM_VAF")!=nullptr; std::vector<int32_t> refcov;` |
| `CAPS_MULTIBUBBLE` | `include/caps_caller.h:4674` | `if (const char* e = std::getenv("CAPS_MULTIBUBBLE")) MAXB = atoi(e); uint32_t rp2 = rp, qB2 = qB;` |
| `CAPS_NAMES` | `src/encoder.cpp:394` | `static const bool CAPS_NAMES = getenv("CAPS_NAMES") != nullptr; // ---- quality column (CAPS_QUAL=1) -------------------------------------------` |
| `CAPS_NM_DEPTH` | `src/encoder.cpp:2392` (+1) | `const int NMD=getenv("CAPS_NM_DEPTH")?atoi(getenv("CAPS_NM_DEPTH")):5; FILE* nv=NMV?fopen(NMV,"w"):nullptr;` |
| `CAPS_NM_MINL` | `src/encoder.cpp:1396` | `if(const char* e=getenv("CAPS_NM_MINL")) nm_minL=(uint32_t)atoi(e); std::vector<uint32_t> nxt,prv,ovl,ch_h,ch_t;` |
| `CAPS_NM_MINQ` | `src/encoder.cpp:435` (+1) | `static const int  NM_MINQ    = getenv("CAPS_NM_MINQ") ? atoi(getenv("CAPS_NM_MINQ")) : 20; // Emit 2-mismatch (clustered) partners as candidate pairs. Opt-in: it is a` |
| `CAPS_NM_PAIRS` | `src/encoder.cpp:439` | `static const bool NM_PAIRS   = getenv("CAPS_NM_PAIRS") != nullptr; static std::vector<std::pair<uint64_t,uint64_t>> g_contig_spans;   // [start,end) in pg coords` |
| `CAPS_NM_VCF` | `src/encoder.cpp:434` (+1) | `static const bool NM_QUAL    = getenv("CAPS_NM_VCF") != nullptr; static const int  NM_MINQ    = getenv("CAPS_NM_MINQ") ? atoi(getenv("CAPS_NM_MINQ")) : 20;` |
| `CAPS_NO_ALIGNPAIR` | `include/caps_caller.h:4901` | `if (!std::getenv("CAPS_NO_ALIGNPAIR")) { int MAXEV = 12;` |
| `CAPS_NO_CLOSE_ANCHOR` | `include/caps_caller.h:4705` | `if (!std::getenv("CAPS_NO_CLOSE_ANCHOR")) { uint32_t ra = bub.apos + (bub.type == 0 ? (uint32_t)bub.len : 0u);` |
| `CAPS_NO_COLLAPSE` | `include/caps_caller.h:796` (+1) | `const bool NO_COLLAPSE = std::getenv("CAPS_NO_COLLAPSE") != nullptr; const int  K           = 25;` |
| `CAPS_NO_HPBUBBLE` | `include/caps_caller.h:276` | `if (!std::getenv("CAPS_NO_HPBUBBLE")) { size_t ia = pA + d, ib = qB + d;` |
| `CAPS_NO_HPCLOSE` | `include/caps_caller.h:4721` | `if (!same && !std::getenv("CAPS_NO_HPCLOSE")) { for (int shift = -(int)bub.len; shift <= (int)bub.len && !same; ++shift) {` |
| `CAPS_NO_INDELS` | `include/caps_caller.h:4373` | `if (!std::getenv("CAPS_NO_INDELS") && cdb.contigs.size() >= 2) { constexpr int BK = 25, FLANK = 15;` |
| `CAPS_NO_JUNCTION` | `include/caps_caller.h:5041` | `if (!std::getenv("CAPS_NO_JUNCTION")) { // Fragment half-length. EBWT2SNP emits 2k+1 with k=30 (61 bp)` |
| `CAPS_NO_PCLUSTER` | `include/caps_caller.h:5180` | `if (!std::getenv("CAPS_NO_PCLUSTER")) { const int AK = 25, LW = 40;      // anchor k-mer, left window` |
| `CAPS_NO_PCLUSTER_AF` | `include/caps_caller.h:5591` | `if (!std::getenv("CAPS_NO_PCLUSTER_AF")) { int dp = (ccid < pcov.size() && apos2 < pcov[ccid].size())` |
| `CAPS_NO_PCLUSTER_HI` | `include/caps_caller.h:5613` | `if (!std::getenv("CAPS_NO_PCLUSTER_HI") && af > 1.0 - MAF) { ++pc_drop_afhi; continue; } }` |
| `CAPS_NO_QOVERLAP` | `src/encoder.cpp:760` | `if((CAPS_NAMES\|\|CAPS_QUAL) && !g_input_path.empty() && !getenv("CAPS_NO_QOVERLAP")){ q_thread = std::thread([&]{` |
| `CAPS_NO_REMAP` | `include/caps_caller.h:795` | `const bool NO_REMAP    = std::getenv("CAPS_NO_REMAP")    != nullptr; const bool NO_COLLAPSE = std::getenv("CAPS_NO_COLLAPSE") != nullptr;` |
| `CAPS_NO_XSNV` | `include/caps_caller.h:5785` | `if (!std::getenv("CAPS_NO_XSNV")) { constexpr int XSNV_FLANK = 15, XSNV_MIN_ANCHORS = 15;` |
| `CAPS_PCDBG` | `include/caps_caller.h:4927` (+5) | `if (std::getenv("CAPS_PCDBG")) fprintf(stderr, "[ap] pairs=%zu new-events=%zu total-events=%zu\n",` |
| `CAPS_PCLUSTER_AF` | `include/caps_caller.h:5602` | `if (const char* e = std::getenv("CAPS_PCLUSTER_AF")) afmin = atof(e); if (af < afmin) { ++pc_drop_aflo; continue; }` |
| `CAPS_PCLUSTER_ANCH` | `include/caps_caller.h:5578` | `if (const char* e = std::getenv("CAPS_PCLUSTER_ANCH")) PANCH = atoi(e); if ((int)panch[kv.first].size() < PANCH) continue;` |
| `CAPS_PCLUSTER_COLSUB` | `include/caps_caller.h:5187` | `const CallData& pc_cd = std::getenv("CAPS_PCLUSTER_COLSUB") ? cd : cdb; // 'cov' above is indexed by the BUBBLE substrate's contigs, so it` |
| `CAPS_PCLUSTER_MIN` | `include/caps_caller.h:5565` | `if (const char* e = std::getenv("CAPS_PCLUSTER_MIN")) PMIN = atoi(e); size_t n_pc = 0;` |
| `CAPS_PCLUSTER_UNIQ` | `include/caps_caller.h:5316` | `if (const char* e = std::getenv("CAPS_PCLUSTER_UNIQ")) PUNIQ = atoi(e); _iplap("pcluster: build flat array");` |
| `CAPS_PCLUSTER_VERIFY` | `include/caps_caller.h:5359` | `if (std::getenv("CAPS_PCLUSTER_VERIFY")) { std::unordered_map<uint64_t,uint32_t> ocount;` |
| `CAPS_PHASE` | `src/encoder.cpp:381` | `static const bool CAPS_PHASE = getenv("CAPS_PHASE") != nullptr; static std::chrono::steady_clock::time_point g_t0, g_tp;` |
| `CAPS_PILEUP` | `src/decoder.cpp:1159` | `(getenv("CAPS_PILEUP")!=nullptr \|\| getenv("CAPS_FMI")!=nullptr); if(mode=="index"){` |
| `CAPS_PILEUP_MIN` | `src/decoder.cpp:1728` (+1) | `const uint32_t KMIN = getenv("CAPS_PILEUP_MIN")?(uint32_t)atoi(getenv("CAPS_PILEUP_MIN")):3; FILE* vf=fopen((outdir+".sites").c_str(),"wb");` |
| `CAPS_PLOIDY` | `include/caps_caller.h:1262` | `if (const char* pe = std::getenv("CAPS_PLOIDY")) { int v = atoi(pe); if (v >= 2 && v <= 4) PLOIDY = v; } ` |
| `CAPS_PV_REBUILD` | `include/caps_caller.h:5283` | `const bool PV_FROM_KIDX = (&pc_cd == &cdb) && !std::getenv("CAPS_PV_REBUILD"); std::vector<PE> pv;` |
| `CAPS_QBATCH` | `include/quality_coder.h:321` | `if(const char* e=getenv("CAPS_QBATCH")){ long v=atol(e); if(v>0) NB=(size_t)v; } {` |
| `CAPS_QBLOCK_MB` | `include/quality_coder.h:270` | `if(const char* e = getenv("CAPS_QBLOCK_MB")){ long v=atol(e); if(v>0) BLOCK_BYTES=(size_t)v*1024u*1024u; } Encoded E;` |
| `CAPS_QIDX` | `src/decoder.cpp:616` (+1) | `std::string qidx = getenv("CAPS_QIDX") ? getenv("CAPS_QIDX") : (std::string(arcpath) + ".qidx"); FILE* qf = fopen(qidx.c_str(), "rb");` |
| `CAPS_QOFF` | `include/quality_coder.h:131` (+1) | `if(const char* e=getenv("CAPS_QOFF")){ int v=atoi(e); if(v>=1&&v<=J.nc) J.nc=v; } J.ntrial = J.nc * J.nstrat;` |
| `CAPS_QSTRAT` | `include/quality_coder.h:130` (+2) | `if(const char* e=getenv("CAPS_QSTRAT")){ int v=atoi(e); if(v>=1&&v<=4) J.nstrat=v; } if(const char* e=getenv("CAPS_QOFF")){ int v=atoi(e); if(v>=1&&v<=J.nc) J.nc=v; }` |
| `CAPS_QTHREADS` | `src/encoder.cpp:762` | `if(const char* e=getenv("CAPS_QTHREADS")){ int v=atoi(e); if(v>0) omp_set_num_threads(v); } if(CAPS_NAMES && !g_NM_done){ g_NM = nmc::encode_from_fastq(g_input_path.c_str` |
| `CAPS_QTIME` | `src/decoder.cpp:595` | `const bool QTIME = getenv("CAPS_QTIME")!=nullptr; auto _qt0 = std::chrono::steady_clock::now();` |
| `CAPS_QUAL` | `src/encoder.cpp:400` | `static const bool CAPS_QUAL  = getenv("CAPS_QUAL")  != nullptr; // ── NAMES AND QUALITY ARE COMPUTED ONCE, BEFORE THE CANDIDATE FORK ─────────` |
| `CAPS_QUERY_CONTAIN` | `src/decoder.cpp:768` | `const bool want_contain = sq && getenv("CAPS_QUERY_CONTAIN"); bool contain_done = false;` |
| `CAPS_QUERY_INDEX` | `src/decoder.cpp:767` | `const char* csel = getenv("CAPS_QUERY_INDEX"); const bool want_contain = sq && getenv("CAPS_QUERY_CONTAIN");` |
| `CAPS_QUERY_MINSEED` | `src/decoder.cpp:321` | `if(const char* e = getenv("CAPS_QUERY_MINSEED")){ int v = atoi(e);` |
| `CAPS_QUERY_MM` | `src/decoder.cpp:251` | `const char* e = getenv("CAPS_QUERY_MM"); if(!e) return 0;` |
| `CAPS_SEEDCAP` | `src/encoder.cpp:2713` | `const char* forced = getenv("CAPS_SEEDCAP"); if(forced) CAP = (size_t)atoll(forced);` |
| `CAPS_SEEDKEEP` | `src/encoder.cpp:2711` | `if(const char* e=getenv("CAPS_SEEDKEEP")){ double v=atof(e); if(v>0&&v<=1) KEEPFRAC=v; } size_t CAP = 0;` |
| `CAPS_SKIP_QUAL` | `src/decoder.cpp:495` (+3) | `if (has("qual_body") && (out_qtext \|\| out_qbits) && !getenv("CAPS_SKIP_QUAL")) { auto lenb_e = dec("read_lengths", 2);` |
| `CAPS_SPANS` | `src/encoder.cpp:431` | `static const bool CAPS_SPANS = (getenv("CAPS_SPANS") != nullptr) \|\| CAPS_CALL; // Near-miss base-quality mask: built only when the near-miss VCF is requested,` |
| `CAPS_STRBUBBLE` | `include/caps_caller.h:321` | `if (std::getenv("CAPS_STRBUBBLE")) { const size_t ia = pA + d, ib = qB + d;` |
| `CAPS_SWEEP_TIMING` | `src/encoder.cpp:1483` | `const bool SWTIME = getenv("CAPS_SWEEP_TIMING")!=nullptr; size_t pr_total=0; size_t sd_hit=0; size_t cand_seen=0;` |
| `CAPS_TRACE` | `include/caps_caller.h:4747` (+7) | `if (std::getenv("CAPS_TRACE")) fprintf(stderr, "[trace] DROP-CLOSE cid=%u apos=%u\n", rc_, bub.apos);` |
| `CAPS_XMI` | `src/decoder.cpp:1775` | `if(getenv("CAPS_XMI")){ uint16_t minL = 0xFFFF;` |
| `CAPS_XMI_STRIDE` | `src/decoder.cpp:1788` (+1) | `const int S = getenv("CAPS_XMI_STRIDE") ? atoi(getenv("CAPS_XMI_STRIDE")) : 8; if(K < 8 \|\| S < 1){` |
| `CHUNK_PAR` | `include/coders_inproc.h:746` | `if(const char* e=getenv("CHUNK_PAR")){ unsigned v=atoi(e); if(v) CPAR=v; } unsigned T = (unsigned)std::min<size_t>(spans.size(),` |
| `COST_GATE` | `src/encoder.cpp:3776` | `const bool COSTGATE = getenv("COST_GATE") != nullptr; // Greedy left-to-right parse over ONE slice of the destination. Matches` |
| `DBG_DST` | `src/encoder.cpp:4552` (+1) | `const uint64_t want = getenv("DBG_DST")?strtoull(getenv("DBG_DST"),nullptr,10):UINT64_MAX; size_t shown=0, mmiE=0;` |
| `DBG_FAILMODE` | `src/encoder.cpp:3775` (+1) | `const bool DBGFAIL = getenv("DBG_FAILMODE") != nullptr; const bool COSTGATE = getenv("COST_GATE") != nullptr;` |
| `DBG_O2U` | `src/encoder.cpp:3568` | `if(getenv("DBG_O2U")){ fprintf(stderr,"[dbg-o2u] size=%zu  first 32: ",orig2uid.size());` |
| `DBG_OVL` | `src/encoder.cpp:2344` | `if(getenv("DBG_OVL")){ // Overlap-length histogram over committed links, as a FRACTION of read` |
| `DBG_POS` | `src/encoder.cpp:4571` | `if(const char* dbg=getenv("DBG_POS")){ const uint64_t target=strtoull(dbg,nullptr,10);` |
| `DBG_STREAMS` | `src/encoder.cpp:4549` | `if(getenv("DBG_STREAMS")){ fprintf(stderr,"[DBGSTREAM] ref_mmcnt.size=%zu ref_mmpos.size=%zu ref_mmref.size=%zu ref_mmobs.size=%zu\n",` |
| `DUMP_MM` | `src/encoder.cpp:3274` (+1) | `if(!getenv("DUMP_MM")){ std::vector<uint64_t>().swap(rpk); std::vector<uint64_t>().swap(woff);` |
| `DUMP_PERM` | `src/encoder.cpp:3496` | `if(getenv("DUMP_PERM")){ const uint64_t PL=main_pg_end;   // RC scan ran before survivors were appended` |
| `DUMP_PG` | `src/encoder.cpp:4594` (+1) | `if(getenv("DUMP_PG")){ FILE* f=fopen("pg_full.txt","wb"); fwrite(pg.data(),1,pg.size(),f); fclose(f);` |
| `FWD_SELF` | `src/encoder.cpp:3679` (+1) | `const bool   FWD_SELF=(getenv("FWD_SELF")?atoi(getenv("FWD_SELF"))!=0:true); ` |
| `GSEARCH` | `src/encoder.cpp:1915` | `if(getenv("GSEARCH")){ uint32_t glo = Lmax/25; if(glo<6) glo=6;` |
| `LOG_PICK` | `include/coders_inproc.h:676` | `if(getenv("LOG_PICK")) fprintf(stderr,"  [pick] n=%-10zu u32=%d -> method %d\n", n, (int)u32shaped, best);` |
| `LOG_STREAMS` | `src/encoder.cpp:4881` (+1) | `if(getenv("LOG_STREAMS")) fprintf(stderr,"  [mmcnt] split=%zu (flags %zu + vals %zu)  flat=%zu  -> %s\n",` |
| `LOWCOV` | `src/encoder.cpp:813` (+1) | `if(getenv("LOWCOV")){ SW = 8; } const uint64_t SWMASK = (SW>=32)?~0ULL:((1ULL<<(2*SW))-1);` |
| `MAXMAP` | `src/encoder.cpp:2592` (+1) | `if(getenv("MAXMAP")) MAXMAP = (uint32_t)atoi(getenv("MAXMAP"));   // override for sweeps fprintf(stderr,"  [maxmap] Lmax=%u -> MAXMAP=%u (L/%u)\n", Lmax, MAXMAP, MAXMAP_D` |
| `MAXMAP_DIV` | `src/encoder.cpp:2560` (+1) | `const uint32_t MAXMAP_DIV = (uint32_t)(getenv("MAXMAP_DIV")?atoi(getenv("MAXMAP_DIV")):13); const uint32_t MIN_MAXMAP = 6;` |
| `MEM_MAXMM_OVERRIDE` | `src/encoder.cpp:4175` (+1) | `if(getenv("MEM_MAXMM_OVERRIDE")) MEM_MAXMM = atoi(getenv("MEM_MAXMM_OVERRIDE")); fprintf(stderr,"[MMTOL] leftover_frac=%.3f het_pair_frac=%.4f -> MEM_MAXMM=%d\n", leftove` |
| `MINOV` | `src/encoder.cpp:2279` (+1) | `if(getenv("MINOV")) sweep_minov=(uint32_t)atoi(getenv("MINOV"));  // override for sweeps sweep();` |
| `MMGOOD` | `src/encoder.cpp:2605` (+1) | `const uint8_t MMGOOD=(uint8_t)(getenv("MMGOOD")?atoi(getenv("MMGOOD")):0); const uint32_t NPARTS=(Lmax>=SEEDW)?((Lmax-SEEDW)/SEEDSTRIDE+1):1;` |
| `NMC_NODICT` | `include/names_coder.h:131` | `static const bool NO_DICT = getenv("NMC_NODICT")!=nullptr;  // drop the value dictionary static const uint32_t MAXTOK=1024;` |
| `NMC_NOZD` | `include/names_coder.h:130` | `static const bool NO_ZD   = getenv("NMC_NOZD")!=nullptr;    // drop ID_ZDELTA static const bool NO_DICT = getenv("NMC_NODICT")!=nullptr;  // drop the value dictionary` |
| `NODEDUP` | `src/encoder.cpp:514` (+3) | `bool NODEDUP = getenv("NODEDUP") && atoi(getenv("NODEDUP")); double DUPFRAC = -1.0;` |
| `PROBE_FRAC` | `include/coders_inproc.h:624` | `const char* e = getenv("PROBE_FRAC"); // PgRC2 uses 0.20. Measured on real files: 0.20 costs E. coli +50,149 B` |
| `SECOND_SELF` | `src/encoder.cpp:4207` (+1) | `if(getenv("SECOND_SELF")){ const int SSMODE=atoi(getenv("SECOND_SELF"));  // 1=fwd only, 2=fwd+rc` |
| `SEQT` | `src/encoder.cpp:4715` (+1) | `const unsigned SEQT = getenv("SEQT") ? (unsigned)atoi(getenv("SEQT")) : 4; ` |
| `TMPDIR` | `include/caps_caller.h:1529` | `const char* base = std::getenv("TMPDIR"); if (!base) base = "/tmp";` |
| `VERIFY_DUMP` | `src/encoder.cpp:4622` (+1) | `const bool VDUMP = getenv("VERIFY_DUMP") && atoi(getenv("VERIFY_DUMP")); auto vdump=[&](const char* fn, const MemStream& m){` |
