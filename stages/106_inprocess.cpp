// Greedy SCS + pigeonhole mapping, with the overlap search done in one pass.
//
// The earlier prototypes rebuilt a hash index once per candidate overlap length
// and re-hashed every read's L bytes each time -- O(n*L^2), about 19 billion byte
// reads on 851k reads, which is the entire reason they took 40 s while PgRC2's
// sort-and-merge takes 1.5 s. Nothing about the ALGORITHM was slow; the indexing
// was.
//
// Here each read is scanned once. A 32-base seed packs exactly into a uint64
// (2 bits/base), so the seed comparison is integer equality with no false
// positives from hashing, and sliding the seed one base is a shift-or. For a
// read A we walk its suffix offsets outward from the longest overlap; the first
// offset whose seed hits a read B and whose full suffix-prefix span verifies IS
// A's longest overlap, so we stop there. That makes the search O(n * offsets
// actually tried) instead of O(n * L^2).
//
// Stage C then maps whatever never chained into the finished pseudogenome using
// the pigeonhole lemma (a read with <= m mismatches cut into m+1 parts must have
// one part matching exactly), forward and then against the reverse-complemented
// text -- one index, one text, rather than a two-strand index.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <atomic>
#include <omp.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <map>
#include <dirent.h>
#if defined(__GLIBC__)
#include <malloc.h>
#endif

static const uint32_t NONE    = UINT32_MAX;
static const uint32_t SEED    = 32;               // exactly one uint64 at 2 bits/base
static const uint64_t SEEDMSK = ~0ULL;            // 32 bases * 2 bits = 64

static inline int b2(char c){
    switch(c){ case 'A':return 0; case 'C':return 1; case 'G':return 2; case 'T':return 3; }
    return -1; }
static inline uint64_t fnv(const char* p,uint32_t n){
    uint64_t h=1469598103934665603ULL;
    for(uint32_t i=0;i<n;++i){ h^=(uint8_t)p[i]; h*=1099511628211ULL; } return h; }
static void rc_inplace(std::string& s){
    auto comp=[](char c){ switch(c){case 'A':return 'T';case 'C':return 'G';
                                    case 'G':return 'C';case 'T':return 'A';} return c; };
    if(s.empty()) return;
    size_t i=0,j=s.size()-1;
    while(i<j){ char a=comp(s[i]),b=comp(s[j]); s[i]=b; s[j]=a; ++i; --j; }
    if(i==j) s[i]=comp(s[i]);
}
// pack the 32 bases at p; false if any non-ACGT
static inline bool pack(const char* p,uint64_t& out){
    uint64_t k=0;
    for(uint32_t i=0;i<SEED;++i){ int v=b2(p[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
    out=k; return true;
}

// ============================ STAGE 106 ============================
// SINGLE-PROCESS, IN-MEMORY, COMPRESS-AND-RELEASE.
//
// PgRC2's release build writes ZERO intermediate files and holds every stream
// as an in-memory ostringstream, compressing each and disposing of it right
// away (SeparatedPseudoGenomePersistence.cpp:688-696 initDest; the
// disposeReadsList/clear calls in pgrc-encoder.cpp:157-226). We wrote 18
// files and re-read them from separate coder processes -- 361 MB of I/O for a
// 456 MB input, all of it pure overhead (see PGRC2_DISK_ARCHITECTURE.md).
//
// Every stream here is an open_memstream FILE*, so all the existing fwrite/
// fputc calls are untouched and provably write the same bytes -- they just
// land in RAM. Each stream is then coded IN-PROCESS (coders_inproc.h, which
// reproduces the exact algorithms; verified to produce byte-identical sizes
// to the standalone binaries) and its raw buffer is freed immediately, so
// peak RAM is bounded by the largest live stream rather than their sum.
// Output is one archive; no intermediates ever touch disk.
#include <unordered_set>
#include <fstream>
#include <functional>
#include "coders_inproc.h"
#include "seqpar_core.h"
#include "names_coder.h"
#include "quality_coder.h"
#include "caps_pack.h"
#include "caps_caller.h"


// ── HUGE-PAGE HINT FOR THE BIG ARRAYS ──────────────────────────────────────
// MEASURED: the encoder runs at a 52.25% dTLB load-miss rate -- 6.7 BILLION
// misses out of 12.9 billion loads, against 1.6 trillion total cycles. A
// healthy program is under 1%.
//
// The cause is page size, not algorithm. One candidate's working set is
// ~1142 MB on HG002 (rpk 445 MB, seed/pent/ptab 96 MB each, nxt/prv/ovl/ch_h/
// ch_t 48 MB each). In 4 KB pages that is ~292,000 pages competing for a
// ~1500-entry TLB, so essentially every random access pays a page-table walk
// on top of its cache miss. At K=8 candidates the pressure is eightfold, which
// is precisely the 4x penalty over perfect scaling that round 2 shows
// (503.8 s measured against 126 s predicted).
//
// This machine has transparent huge pages set to [madvise]: they are available
// but only to callers that ASK. std::vector never asks. One madvise per big
// array raises TLB reach from ~6 MB to ~3 GB.
//
// PURELY A MEMORY-LAYOUT HINT: it changes no value the program computes, so
// the archive is identical by construction. If the kernel declines the hint
// the code is unaffected.
// Allocate-then-hint-then-fill. Calling madvise AFTER a vector is filled does
// nothing: the pages are already faulted in at 4 KB and MADV_HUGEPAGE will not
// retroactively collapse them (khugepaged may, far too slowly to matter).
// Measured proof that the naive version was inert: AnonHugePages stayed at
// 0.0 MB with hints on 16 arrays. reserve() takes the mapping without touching
// it, the hint lands on untouched memory, and assign() then faults it in as
// 2 MB pages.
// ── RAW-THREAD BUDGET ───────────────────────────────────────────────────────
// omp_set_num_threads() bounds OpenMP regions ONLY. Three stages fan out with
// raw std::threads sized by hardware_concurrency(), which ignores it entirely,
// so every candidate spawned 12 threads no matter how many candidates were
// running: the pigeonhole scan at 2 concurrent group children was 24 threads on
// 12 cores, and MEM run() at 4 concurrent candidates was 48 on 12.
//
// Set once at each fork point to that process's share of the machine.
static unsigned g_thr_budget = 0;                  // 0 = not forked: take the box
static inline unsigned thr_budget(){
    if(g_thr_budget) return g_thr_budget;
    unsigned h=std::thread::hardware_concurrency(); return h?h:1u;
}
static inline void hugehint(void* p, size_t bytes);
template<class V, class T>
static inline void hugefill(V& v, size_t n, const T& val){
    v.clear(); v.shrink_to_fit(); v.reserve(n);
    hugehint((void*)v.data(), n*sizeof(T));
    v.assign(n, val);
}
// REHOME a vector that was built incrementally. rpk/woff/rlen cannot use
// hugefill: their final size is not known until the FASTQ has been read, so
// they are grown with push_back/insert. They were therefore given a bare
// hugehint() after the fill -- which is precisely the inert pattern the comment
// above warns about, and it went unnoticed on the single largest array in the
// program. rpk is 445 MB on HG002, 39% of a candidate's working set, and it is
// one of the three arrays the pigeonhole loop dereferences at a random index
// 1.93 billion times (see the prefetch note there); it had no huge pages at
// all, and shrink_to_fit() re-faulted every page at 4 KB immediately before the
// hint that was supposed to help.
//
// Copying into a freshly reserved, hinted buffer is the only way to fix this
// after the fact: reserve() takes the mapping untouched, the hint lands on
// untouched memory, and the copy faults it in as 2 MB pages. Costs one pass
// over the array and holds both copies briefly -- rpk peaks at 706 MB RSS here
// against a 4.2 GB run peak, so the transient is affordable. Below the
// threshold it degrades to the plain shrink_to_fit() it replaces.
//
// PURELY A MEMORY-LAYOUT CHANGE: same bytes, same order, same values.
template<class V>
static inline void hugerehome(V& v){
    const size_t bytes = v.size()*sizeof(v[0]);
    if(bytes < (8u<<20) || v.empty()){ v.shrink_to_fit(); return; }
    V nv; nv.reserve(v.size());
    hugehint((void*)nv.data(), bytes);
    nv.insert(nv.end(), v.begin(), v.end());
    v.swap(nv);
}
static inline void hugehint(void* p, size_t bytes){
#ifdef MADV_HUGEPAGE
    if(!p || bytes < (4u<<20)) return;                 // not worth it below 4 MB
    const uintptr_t HP = 2u<<20;
    uintptr_t a=(uintptr_t)p, st=(a+HP-1)&~(HP-1), en=(a+bytes)&~(HP-1);
    if(en>st) madvise((void*)st, (size_t)(en-st), MADV_HUGEPAGE);
#else
    (void)p; (void)bytes;
#endif
}

struct MemStream {
    FILE*  f    = nullptr;
    char*  buf  = nullptr;
    size_t len  = 0;
    const char* name = "";
    void open(const char* nm){ name=nm; f=open_memstream(&buf,&len); }
    void close(){ if(f){ fclose(f); f=nullptr; } }
    void release(){ close(); if(buf){ free(buf); buf=nullptr; } len=0; }
    std::vector<uint8_t> bytes() const {
        return std::vector<uint8_t>((const uint8_t*)buf,(const uint8_t*)buf+len);
    }
};

// Every stream the decoder needs. Names match the old filenames one-for-one
// so the layer accounting stays directly comparable.
struct Streams {
    MemStream literal, mem_triples, pos_abs, pos_strand, read_lengths,
              orig2uid, mm_ref, mm_obs, mm_pos, mm_count, n_pos, n_indices, n_cnt;
    void openAll(){
        literal.open("literal"); mem_triples.open("mem_triples");
        pos_abs.open("pos_abs"); pos_strand.open("pos_strand");
        read_lengths.open("read_lengths"); orig2uid.open("orig2uid");
        mm_ref.open("mm_ref"); mm_obs.open("mm_obs"); mm_pos.open("mm_pos");
        mm_count.open("mm_count"); n_pos.open("n_pos");
        n_indices.open("n_indices"); n_cnt.open("n_cnt");
    }
};
static Streams STR;

// Archive writer: append a coded stream, then free the source immediately.
// CAPSULE container.
//
// The archive was a bare sequence of (uint64 length, bytes) in whatever order
// the job list happened to run, with the pseudogenome parameters written to a
// separate pg_params.txt. Nothing identified the format, its version, or which
// stream was which, so a decoder could only work by replicating the encoder's
// job order exactly -- and any added stream (names, quality) would silently
// shift every offset.
//
// Layout:
//   magic   "CAPSULE" + '\0'            8 B
//   version uint16                       2 B
//   pg_len, main_pg_end  uint64 x2      16 B
//   n_streams uint16                     2 B
//   per stream: namelen u8, name, len u64, payload
//
// Streams are self-identifying, so the decoder matches by NAME and a new stream
// is additive rather than a format break. Header cost is ~300 B against
// archives of 2.6-27 MB (0.001-0.01%), and it is counted in ARCHIVE_TOTAL.
struct Archive {
    FILE* out;
    size_t total=0;
    std::vector<std::pair<const char*,size_t>> parts;
    Archive(const char* path, uint64_t pg_len, uint64_t main_end, uint32_t minmem){
        out=fopen(path,"wb");
        if(!out) return;
        fwrite("CAPSULE\0",1,8,out);
        uint16_t ver=2; fwrite(&ver,2,1,out);
        fwrite(&pg_len,8,1,out); fwrite(&main_end,8,1,out);
        // MINMEM: mem_len stores (match length - MINMEM), and destination gaps
        // are measured from the END of the previous match, so both the lengths
        // and the destinations are unrecoverable without it.
        fwrite(&minmem,4,1,out);
        nstream_pos=ftell(out);
        uint16_t ns=0; fwrite(&ns,2,1,out);       // patched in finish()
        total += 8+2+8+8+4+2;
    }
    long nstream_pos=0; uint16_t nstreams=0;
    void put(const char* name, const std::vector<uint8_t>& coded){
        // (timing is recorded by the caller via tput)
        uint64_t n=coded.size();
        uint8_t nl=(uint8_t)strlen(name);
        fwrite(&nl,1,1,out); fwrite(name,1,nl,out);
        fwrite(&n,8,1,out);
        if(n) fwrite(coded.data(),1,n,out);
        total += 1 + nl + 8 + n;
        ++nstreams;
        parts.emplace_back(name,n);
    }
    void finish(){
        if(out){ fseek(out,nstream_pos,SEEK_SET); fwrite(&nstreams,2,1,out); }
        if(out) fclose(out);
        out=nullptr;
        for(auto& p:parts) fprintf(stderr,"  [archive] %-14s %10zu B\n",p.first,p.second);
        fprintf(stderr,"ARCHIVE_TOTAL=%zu\n",total);
        printf("ARCHIVE_TOTAL=%zu\n",total);
    }
};

// ---- phase timing (CAPS_PHASE=1) --------------------------------------------
// Speed work needs to know which phase owns the wall clock. perf attributes 60%
// to coder threads and 9% to libgomp, but that does not say how much of the run
// is assembly versus coding, which is the number that decides where to look.
static const bool CAPS_PHASE = getenv("CAPS_PHASE") != nullptr;
static std::chrono::steady_clock::time_point g_t0, g_tp;
static void phase_init(){ g_t0 = g_tp = std::chrono::steady_clock::now(); }

// ---- names / read-ID column (CAPS_NAMES=1) ----------------------------------
// Phase 2. Off by default so every Phase-1 archive stays byte-identical and
// the locked size result is provably untouched (standing rule 2: an
// output-preserving change must be cmp-identical on all 7 files). The coder
// re-reads the ID column straight from the input FASTQ rather than being fed
// headers captured during the main parse: capturing them would hold every
// header resident (~5 GB on T. cacao) and defeat the bounded-memory design
// that names_coder.h exists to provide. Two extra streaming passes over the
// ID column cost ~1 s per 5 M reads and no RAM.
static const bool CAPS_NAMES = getenv("CAPS_NAMES") != nullptr;
// ---- quality column (CAPS_QUAL=1) -------------------------------------------
// Wraps the vendored fqzcomp_qual codec; see include/quality_coder.h for why
// it is vendored rather than reimplemented. Same gating discipline as names:
// with CAPS_QUAL unset the archive is byte-identical to one built without it,
// so the locked Phase 1 result is untouched.
static const bool CAPS_QUAL  = getenv("CAPS_QUAL")  != nullptr;
// ── NAMES AND QUALITY ARE COMPUTED ONCE, BEFORE THE CANDIDATE FORK ─────────
// Neither column depends on MAXMAP or MINOV: they are functions of the input
// file alone. They were being encoded inside each candidate, so an 8-point
// sweep did the identical work eight times -- and quality is by far the most
// expensive single step in the encoder (measured 7.09 s of a 10.49 s
// candidate, 68%). Hoisting them above the fork computes them once with the
// whole machine, and every child inherits the result through copy-on-write.
static nmc::Encoded g_NM;  static bool g_NM_done = false;
static qlc::Encoded g_QL;  static bool g_QL_done = false;
// ---- Claim 2: reference-free variant calling (CAPS_CALL=1) ------------------
// Captures each chain (round-1/round-2 main pg, and the second-region sweep)
// as its own discrete contig BEFORE the MEM self-match stage runs -- MEM only
// changes which pg bytes get literal-encoded, never the pg content itself, so
// contig boundaries captured here stay valid for the caller regardless of
// what MEM later does. See include/caps_caller.h for why this matters (the
// indel-bubble mechanism needs genuinely separate haplotype contigs, which a
// single merged pg coordinate space would erase). Zero cost when unset.
static const bool CAPS_CALL  = getenv("CAPS_CALL")  != nullptr;
static std::vector<std::pair<uint64_t,uint64_t>> g_contig_spans;   // [start,end) in pg coords
static std::string g_input_path;
static void phase(const char* name){
    if(!CAPS_PHASE) return;
    auto now = std::chrono::steady_clock::now();
    double d = std::chrono::duration<double>(now - g_tp).count();
    double tot = std::chrono::duration<double>(now - g_t0).count();
    long rss = 0;
    if(FILE* f=fopen("/proc/self/status","r")){ char b[256];
        while(fgets(b,sizeof b,f)) if(!strncmp(b,"VmHWM:",6)){ sscanf(b+6,"%ld",&rss); break; }
        fclose(f); }
    fprintf(stderr,"[PHASE] %-14s %7.2fs  cum %7.2fs  peakRSS %6ld MB\n", name, d, tot, rss/1024);
    g_tp = now;
}

int main(int argc,char** argv){
    phase_init();
    if(argc>1) g_input_path = argv[1];
    // Resolve to absolute IMMEDIATELY, before any candidate/GSEARCH fork can
    // chdir() into a per-candidate scratch directory. CALL_VCF and
    // CAPS_DUMP_CONTIGS already get this treatment for exactly this reason
    // (see the CAPS_CALL child-setup block below) -- g_input_path did not,
    // which silently broke CAPS_NAMES/CAPS_QUAL whenever a relative input
    // path was used through encode_adaptive.sh's sweep: nmc::encode_from_fastq
    // and qlc::encode_from_fastq reopen g_input_path to extract names/quality,
    // and after the chdir a relative path no longer points at the file the
    // caller meant, so both silently encoded zero records while the archive
    // still reported success. Found 2026-09-03 by actually running all three
    // claims' features together on one archive, not by isolated per-feature
    // testing. A no-op when the path is already absolute or CAPSULE is run
    // without CANDIDATES/GSEARCH (the overwhelmingly common case), so this
    // cannot change any existing single-candidate archive.
    if(!g_input_path.empty() && g_input_path[0]!='/'){
        char cwdbuf[4096];
        if(getcwd(cwdbuf,sizeof cwdbuf)) g_input_path = std::string(cwdbuf) + "/" + g_input_path;
    }
    // RAM FIX -- confirmed real driver of the C. elegans-scale RSS gap
    // after three application-level hypotheses were measured and ruled
    // out (allrefs/cleanRefs double-holding, c/cr scope overlap, res[]
    // per-thread capacity -- see LAYER_BY_LAYER_ANALYSIS.md section 3c
    // for the full trail). glibc's default allocator gives each thread
    // its own malloc arena (up to 8x core count by default); freed memory
    // in one arena isn't consolidated across others, which shows up as
    // real RSS bloat the application's own logical lifetime never
    // predicts. Confirmed directly: MALLOC_ARENA_MAX=1 dropped C. elegans
    // peak RSS from 1,300,968 KB to 1,187,404 KB (-8.7%), byte-identical
    // correctness unaffected (0-diff, same clean-ref count). This is a
    // well-established, safe, zero-correctness-risk allocator tuning
    // (not a novel guess) -- setting it here via mallopt() makes the fix
    // robust to however the binary gets invoked, not dependent on an
    // external environment variable at every call site.
    mallopt(M_ARENA_MAX,1);
    if(argc<2){ fprintf(stderr,"usage: scs5 <in.fq> [maxmm] [minov]\n"); return 1; }
    const int      MAXMM = argc>2?atoi(argv[2]):3;
    // MINOV default 40 -> 16. Re-swept on REAL full files with the current
    // coders: lowering it collapses more of the pseudogenome, and with the
    // stream coding now efficient the shorter pg is worth far more than the
    // extra chain links cost. P. aeruginosa 9,214,602 -> 9,083,594 as MINOV
    // goes 40 -> 16; it floors at 16 because SEEDW=16 means shorter overlaps
    // cannot be found at all. Validated across all five real files:
    // aggregate vs PgRC2 +1.55% -> +2.89%, wins 3/5 -> 4/5, every file
    // BYTE_IDENTICAL. (MINOV=24 was also tried: +2.55%, worse.) S. aureus is
    // the one dataset that prefers a higher MINOV (+1.9% at 40 vs +0.7% at
    // 16) but the aggregate strongly favours 16.
    const uint32_t MINOV = argc>3?(uint32_t)atoi(argv[3]):16;
    // Sweep seed width. 32 bases fills a uint64 exactly, which is why it was
    // fixed there -- but it also silently floors the shortest overlap the sweep
    // can ever see, and the measured trend (MINOV 50 -> 40 -> 32 giving 29.88M
    // -> 28.91M -> 28.28M) runs straight into that floor. A narrower seed costs
    // specificity (more candidates to verify) but lets round 2 pair reads that
    // overlap by less than 32, which after the division are still reads known
    // to tile well.
    // LOWCOV: an explicitly-requested low-coverage mode. It changes NO default
    // and is never auto-enabled, because the signal that would gate it
    // (leftover_frac) is only known AFTER round 1 -- the very stage this
    // controls -- so auto-detection would mean re-running chaining. Cheap
    // sample-based coverage estimators were tried and REJECTED: a 25-mer
    // frequency histogram gives mode=2 for every dataset regardless of true
    // depth, and unique-kmer fraction orders P. aeruginosa at 26x (0.929)
    // ABOVE Drosophila at 2.6x (0.919), because a fixed-size sample conflates
    // genome size with depth.
    //
    // Measured on Drosophila SRR40104920 (2.6x), sweeping SEEDW/MINOV together
    // with everything else fixed -- monotonic all the way to the floor:
    //     16/16  36,851,932     12/12  36,783,282
    //     10/10  36,558,158      8/8   36,103,845  (-2.38% vs tuned baseline)
    // 8 is the hard floor ("below this a seed is noise"). leftover_frac falls
    // 0.812 -> 0.785 and the second region shrinks by 6.3 M bases, i.e. more
    // reads chain instead of being appended raw -- attacking the cost at its
    // source rather than re-encoding it.
    //
    // NOT made a default, and deliberately NOT fitted to the locked suite: on
    // the three datasets checked it wins on H. salinarum (-1.33%) and
    // Drosophila (-2.38%) but LOSES on E. coli (+0.86%), and no validated
    // property separates those cases. Choosing a threshold that happens to
    // suit the 7 locked files would be exactly the overfit this project
    // forbids -- an eighth dataset could land the wrong side of it. So it stays
    // opt-in, evidenced, and off.
    uint32_t SW = argc>4?(uint32_t)atoi(argv[4]):32;
    if(getenv("LOWCOV")){ SW = 8; }
    const uint64_t SWMASK = (SW>=32)?~0ULL:((1ULL<<(2*SW))-1);
    auto packSW=[&](const char* p,uint64_t& out)->bool{
        uint64_t k=0;
        for(uint32_t i=0;i<SW;++i){ int v=b2(p[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
        out=k; return true; };   // >= SEED
    auto T0=std::chrono::steady_clock::now();
    auto rss_mb=[]()->size_t{ FILE* f=fopen("/proc/self/statm","r"); if(!f) return 0;
        long pg_=0,res=0; if(fscanf(f,"%ld %ld",&pg_,&res)!=2) res=0; fclose(f);
        return (size_t)((double)res*(double)sysconf(_SC_PAGESIZE)/1048576.0); };
    // statm gives CURRENT rss, which misses a transient peak inside a stage.
    // VmHWM is the kernel's high-water mark -- the number a RAM claim needs.
    auto hwm_mb=[]()->size_t{ FILE* f=fopen("/proc/self/status","r"); if(!f) return 0;
        char ln[256]; size_t kb=0;
        while(fgets(ln,sizeof ln,f)) if(!strncmp(ln,"VmHWM:",6)){ sscanf(ln+6,"%zu",&kb); break; }
        fclose(f); return kb/1024; };
    auto lap=[&](const char* w){ auto t=std::chrono::steady_clock::now();
        fprintf(stderr,"  [stage] %-26s %6.2f s   rss=%zu MB  peak=%zu MB\n",w,
                std::chrono::duration<double>(t-T0).count(),rss_mb(),hwm_mb()); T0=t; };

    // ── load, drop N-reads, collapse exact duplicates ────────────────────────
    // ── STAGE 21: reads stored 2 bits/base ───────────────────────────────────
    // Profiled, peak RSS hit 224 MB at load, before any index existed -- 96% of
    // PgRC2's entire 232 MB peak, spent purely on holding the reads. 851,275
    // std::strings pay a 32-byte header plus allocator rounding on every 150
    // bases, ~192 B for 150 B of sequence.
    //
    // PgRC2 does not do that: SymbolsPackingFacility keeps reads packed and
    // compares packed against unpacked in place (countSequenceMismatchesVsUnpacked,
    // compareSequenceWithUnpacked). At 4 bases/byte the same reads are ~34 MB.
    //
    // Safe here because N-containing reads are already filtered out above, so
    // the alphabet really is 2-bit. Each read starts on a word boundary: it
    // wastes under a word per read and makes a prefix load aligned.
    std::vector<uint64_t> rpk;            // 32 bases per uint64, base 0 in the top bits
    std::vector<uint64_t> woff;           // starting WORD of each read
    std::vector<uint16_t> rlen;           // length in bases (uint16: >255bp reads are real, uint8 silently dropped them)
    // STAGE 22 (their stage 6, OrderInfo): to restore the original file order we
    // need, for every ORIGINAL read, where its sequence sits in the pg. Dedup
    // collapses duplicates, so the original->unique map has to be kept too.
    std::vector<uint32_t> orig2uid;       // original read (N-filtered out) -> unique id
    std::vector<uint16_t> origlen;        // TRUE length of each ORIGINAL read (see containment note)
    size_t n_in=0,n_filt=0;
    // ADAPTIVE DEDUP, keyed on a MEASURED property of the input.
    // Pre-assembly dedup collapses duplicate reads so position/length/strand/
    // mismatch data is stored once, but forces an orig2uid correlation array
    // over every ORIGINAL read -- our single biggest remaining loss to PgRC2,
    // who pay nothing there (a duplicate is just a 100%-length overlap in
    // their chain mechanism). Which side wins depends on how duplicated the
    // input actually is. Measured on real full files:
    //     E. coli        20.4% dup -> dedup WINS by 3.8%
    //     L. major       10.0% dup -> dedup LOSES by 0.40%
    //     P. aeruginosa   1.2% dup -> dedup LOSES by 0.49%
    // So the break-even sits between 10% and 20.4%; the threshold is set at
    // 15%. This is a formula over a measured input property, not a per-dataset
    // switch, per the standing algorithmic-first rule. A cheap hash-only
    // pre-pass measures the rate before the real load decides.
    bool NODEDUP = getenv("NODEDUP") && atoi(getenv("NODEDUP"));
    double DUPFRAC = -1.0;
    if(!getenv("NODEDUP")){
        std::ifstream pf(argv[1]);
        if(pf){
            std::unordered_set<uint64_t> seenh;
            std::string a1,b1,c1,d1; size_t tot=0, dup=0;
            while(std::getline(pf,a1)&&std::getline(pf,b1)&&std::getline(pf,c1)&&std::getline(pf,d1)){
                if(b1.find('N')!=std::string::npos) continue;
                uint64_t h=1469598103934665603ULL;
                for(char ch:b1){ h^=(uint8_t)ch; h*=1099511628211ULL; }
                ++tot;
                if(!seenh.insert(h).second) ++dup;
            }
            if(tot){ DUPFRAC=(double)dup/(double)tot; NODEDUP = (DUPFRAC < 0.15); }
            fprintf(stderr,"[dedup] measured duplicate fraction %.1f%% -> dedup %s\n",
                    100.0*DUPFRAC, NODEDUP?"OFF":"ON");
        }
    }
    // STAGE 100 -- real, previously-undisclosed data-loss bug: N-containing
    // reads were `continue`'d here with NO storage anywhere, in every stage
    // from 87 through 98. Caught only because a direct question ("how can
    // this be byte-identical if reads are being dropped") forced a check of
    // what "round trip: VERIFIED" actually meant throughout this project --
    // it was always coder-level (does the entropy coder correctly decode ITS
    // OWN dumped intermediate file), never a genuine original-FASTQ-to-
    // decoded-FASTQ comparison. Two different builds both silently dropping
    // the same reads still match each other byte-for-byte; that was mistaken
    // for correctness against the source file, which was never actually
    // tested. Fixed here for real: N-containing reads are stored (full raw
    // sequence text, since 2-bit packing can't represent a 5th symbol) in
    // n_reads.txt, one per line, plus their ORIGINAL read index (needed to
    // splice them back into the right position on decode) in n_indices.bin
    // as raw LE uint32 -- gated on DUMP_LIT since this is a literal-sequence
    // sibling stream, following the project's existing DUMP_* convention.
    // n_pos.bin: byte position of each N within its read (concatenated)
    // n_indices.bin: original index of each N-containing read
    // n_cnt.bin: how many Ns in each N-containing read
    STR.openAll();
    FILE* nrf = STR.n_pos.f;
    FILE* nif = STR.n_indices.f;
    FILE* ncf = STR.n_cnt.f;
    {
        std::ifstream f(argv[1]); std::string a,b,c,d;
        // Dedup index. This used to be
        //     unordered_map<uint64_t, vector<uint32_t>>
        // which costs a 48 B hash node PLUS a separately heap-allocated inner
        // vector (~32 B of malloc header before a single id is stored) for every
        // distinct read -- 80 B/read, and 1.24M individual allocations on
        // E. coli. Measured: named structures at load total 119 MB against a
        // 228 MB peak, and this map is the entire 109 MB difference; on L. major
        // it projects to 362 MB of a 459 MB load peak.
        //
        // PgRC2 has no such structure at all: it sorts read INDICES into blocks
        // by leading symbol and compares neighbours, one flat vector<uint32_t>
        // (ParallelGreedySwiping...Generator.cpp:160-194). We cannot adopt that
        // shape directly -- unique ids here must follow FIRST-OCCURRENCE order,
        // which sorting destroys, and it would force packing duplicate reads too.
        //
        // What carries over is the principle: a flat array, no per-key
        // allocation. Open addressing with a 32-bit hash tag beside the id, so a
        // probe rejects a mismatch without unpacking the stored read. Semantics
        // are unchanged -- same equivalence classes, same representative, same
        // orig2uid -- so the archive must stay byte-identical.
        struct DedupTable {
            std::vector<uint32_t> uid;      // UINT32_MAX = empty
            std::vector<uint32_t> tag;      // low 32 bits of the FNV hash
            size_t mask = 0, used = 0;
            void init(size_t cap){ size_t c=1; while(c<cap) c<<=1;
                uid.assign(c,UINT32_MAX); tag.assign(c,0); mask=c-1; used=0; }
            void grow(){ std::vector<uint32_t> ou=uid, ot=tag; size_t oc=uid.size();
                uid.assign(oc*2,UINT32_MAX); tag.assign(oc*2,0); mask=oc*2-1;
                for(size_t i=0;i<oc;++i) if(ou[i]!=UINT32_MAX){
                    size_t j=ot[i]&mask;
                    while(uid[j]!=UINT32_MAX) j=(j+1)&mask;
                    uid[j]=ou[i]; tag[j]=ot[i]; } }
        } seen;
        seen.init(1u<<20);
        std::vector<uint64_t> tmpw; char ubuf[1024];
        auto packstr=[&](const char* q,uint32_t L,std::vector<uint64_t>& out){
            out.clear();
            for(uint32_t i=0;i<L;i+=32){
                uint64_t v=0;
                for(uint32_t j=0;j<32;++j){
                    const uint32_t idx=i+j;
                    const uint64_t c=(idx<L)?(uint64_t)b2(q[idx]):0ULL;
                    v=(v<<2)|c;
                }
                out.push_back(v);
            }
        };
        while(std::getline(f,a)&&std::getline(f,b)&&std::getline(f,c)&&std::getline(f,d)){
            ++n_in;
            // STAGE 105 -- STRUCTURAL: N-containing reads used to bypass
            // assembly entirely and get dumped as raw text. Measured on real
            // P. aeruginosa that cost 586,364 B (567,140 seqpar + 19,224
            // indices) for 24,350 reads -- 70% of our whole gap to PgRC2 on
            // that file -- and it is NOT a coder problem (seqpar 1 1 already
            // beats xz 567,140 vs 644,304). The reads themselves are fine:
            // mean 1.23 N characters each, max 8, none all-N. They were being
            // excluded from overlap assembly over one or two ambiguous bases.
            // PgRC2 gives them a dedicated nPg so they still get assembled
            // (pgrc-encoder.cpp:184-204, runNPgGeneration). Here they are
            // substituted N->A and sent through the SAME pipeline as every
            // other read, with the N positions kept in a small side stream --
            // same goal as their nPg, without a second pseudogenome, and the
            // substituted base is simply corrected on the way out.
            if(b.find('N')!=std::string::npos){
                ++n_filt;
                if(nif){ const uint32_t idx=(uint32_t)(n_in-1); fwrite(&idx,4,1,nif); }
                uint8_t ncnt=0;
                for(uint32_t j=0;j<b.size();++j){
                    if(b[j]=='N'){
                        if(nrf){ fputc((int)(j>255?255:j),nrf); }
                        b[j]='A';
                        if(ncnt<255) ++ncnt;
                    }
                }
                if(ncf) fputc((int)ncnt,ncf);
            }
            if(b.size()>1023) continue;                       // uint16 length field, matches ubuf/b_ buffer size
            // NODEDUP=1: skip pre-assembly dedup entirely. orig2uid is our
            // single biggest remaining loss to PgRC2 (+427,060 B on real
            // P. aeruginosa, where they pay nothing because a duplicate is
            // just a 100%-length overlap in their chain mechanism). Retested
            // here because the earlier no-dedup measurement predates stage 105
            // and the mm_pos/mm_cnt transforms, which changed the economics.
            bool dup=false;
            if(NODEDUP){
                orig2uid.push_back((uint32_t)rlen.size());
                packstr(b.data(),(uint32_t)b.size(),tmpw);
                woff.push_back(rpk.size());
                rpk.insert(rpk.end(),tmpw.begin(),tmpw.end());
                rlen.push_back((uint16_t)b.size());
                continue;
            }
            const uint64_t h64=fnv(b.data(),(uint32_t)b.size());
            const uint32_t htag=(uint32_t)h64;
            size_t slot=htag&seen.mask;
            for(;;){
                const uint32_t id=seen.uid[slot];
                if(id==UINT32_MAX) break;                 // empty: not present
                if(seen.tag[slot]!=htag){ slot=(slot+1)&seen.mask; continue; }
                if(rlen[id]!=(uint16_t)b.size()){ slot=(slot+1)&seen.mask; continue; }
                const uint64_t* w=&rpk[woff[id]];
                const uint32_t nw=((uint32_t)rlen[id]+31)/32;
                for(uint32_t t=0;t<nw;++t){
                    for(uint32_t j=0;j<32;++j){
                        const uint32_t idx=t*32+j;
                        ubuf[idx]="ACGT"[(w[t]>>(2*(31-j)))&3ULL];
                    }
                }
                if(memcmp(ubuf,b.data(),b.size())==0){ dup=true; orig2uid.push_back(id); break; }
                slot=(slot+1)&seen.mask;
            }
            if(!dup){
                orig2uid.push_back((uint32_t)rlen.size());
                seen.uid[slot]=(uint32_t)rlen.size(); seen.tag[slot]=htag;
                if(++seen.used*2 > seen.uid.size()) seen.grow();
                packstr(b.data(),(uint32_t)b.size(),tmpw);
                woff.push_back(rpk.size());
                rpk.insert(rpk.end(),tmpw.begin(),tmpw.end());
                rlen.push_back((uint16_t)b.size());
            }
        }
        // These three were the hole in the huge-page work: built by push_back,
        // shrink_to_fit()ed (which re-faults every page at 4 KB) and only THEN
        // hinted -- inert. rpk alone is 445 MB. See hugerehome().
        hugerehome(rpk); hugerehome(woff); hugerehome(rlen);
        // VERIFY THE HINT ACTUALLY TOOK. The first version of this work shipped
        // completely inert and was only caught by reading smaps -- a madvise()
        // call in the source proves nothing. Print what the kernel actually gave
        // us so a regression is visible in every log rather than inferred.
        if(FILE* sm=fopen("/proc/self/smaps_rollup","r")){
            char ln[256]; unsigned long hp=0, rss=0;
            while(fgets(ln,sizeof ln,sm)){
                unsigned long v;
                if(sscanf(ln,"AnonHugePages: %lu kB",&v)==1) hp=v;
                else if(sscanf(ln,"Rss: %lu kB",&v)==1) rss=v;
            }
            fclose(sm);
            fprintf(stderr,"[HP] AnonHugePages=%lu MB of Rss=%lu MB after load "
                           "(rpk=%zu MB woff=%zu MB rlen=%zu MB)\n",
                    hp/1024, rss/1024,
                    rpk.size()*sizeof(rpk[0])/1048576,
                    woff.size()*sizeof(woff[0])/1048576,
                    rlen.size()*sizeof(rlen[0])/1048576);
        }
        // INVARIANT: read_lengths is indexed by ORIGINAL read, always.
        // It used to be indexed by unique id, which silently breaks whenever a
        // read's uid is later rewritten -- containment aliasing folds a shorter
        // read into a longer container (rlen[a] < rlen[b]), so the contained
        // read would inherit the container's length and decode too long.
        // Populated here, at load, so the invariant holds on every path rather
        // than only when the containment block happens to run.
        origlen.resize(orig2uid.size());
        for(size_t i=0;i<orig2uid.size();++i) origlen[i]=rlen[orig2uid[i]];
    }
    STR.n_pos.close(); STR.n_indices.close(); STR.n_cnt.close();
    const uint32_t n=(uint32_t)rlen.size();

    uint32_t Lmax=0; for(uint32_t i=0;i<n;++i) Lmax=std::max(Lmax,(uint32_t)rlen[i]);

    // 32 bases starting at an arbitrary ABSOLUTE base index, base 0 in top bits.
    auto w32=[&](uint64_t ab)->uint64_t{
        const size_t wi=ab>>5; const uint32_t sh=(uint32_t)(ab&31)*2;
        const uint64_t hi=rpk[wi];
        if(!sh) return hi;
        const uint64_t lo=(wi+1<rpk.size())?rpk[wi+1]:0ULL;
        return (hi<<sh)|(lo>>(64-sh));
    };
    // read i's bases [o, o+SWb) as an integer, SWb <= 32
    auto rseed=[&](uint32_t i,uint32_t o,uint32_t SWb)->uint64_t{
        return w32(woff[i]*32ULL+o)>>(2*(32-SWb));
    };
    // compare read a's [off,off+L) against read b's [0,L)
    auto rcmp=[&](uint32_t a,uint32_t off,uint32_t b,uint32_t L)->bool{
        uint64_t pa=woff[a]*32ULL+off, pb=woff[b]*32ULL;
        while(L>=32){ if(w32(pa)!=w32(pb)) return false; pa+=32; pb+=32; L-=32; }
        if(L){ const uint64_t m=~0ULL<<(2*(32-L)); if((w32(pa)^w32(pb))&m) return false; }
        return true;
    };
    // Same comparison, but COUNTS mismatching bases instead of stopping at the
    // first, bailing once the count exceeds `cap`. Returns the count (>cap
    // means "more than cap", value not exact).
    //
    // WHY THIS EXISTS. rcmp above is exact, and heterozygosity is by definition
    // a DIFFERENCE, so the exact test can never see a het partner -- it drops
    // it. But the partner is not invisible: the hash probe keys on a 16-base
    // SEED while the verified overlap L runs up to ~148, so whenever the
    // variant sits past the seed (measured ~89% of the time at L=148) the
    // partner read SHARES the seed, is enumerated by the probe loop, fails
    // rcmp by exactly one base, and is discarded unexamined at
    // `if(!rcmp(...)) continue;`.
    //
    // That discarded set -- seed-sharing, exactly-1-mismatch -- is the het
    // signal available at layer 1, and it costs no extra probing to look at:
    // the read is already in hand at the moment rcmp rejects it. This is a
    // DIFFERENT set from the ccnt>=2 ties measured in docs/GRAPH_HARVEST_
    // REFUTED.md, which counted only reads that PASSED rcmp (85 found, 86% of
    // them duplicate reads). Counting the near-misses was never measured.
    //
    // 2-bit-per-base popcount: XOR, fold the two bits of each base together so
    // any differing base becomes a single set bit, mask to one bit per base.
    // Returns the mismatch count (capped at 2 -- larger values are reported as
    // 3 and are not exact) and, when the count is exactly 1, writes the
    // 0-based offset of the differing base within the overlap into `mmpos`.
    auto rcmp_mm1=[&](uint32_t a,uint32_t off,uint32_t b,uint32_t L,int& mmpos)->int{
        uint64_t pa=woff[a]*32ULL+off, pb=woff[b]*32ULL;
        int mm=0, base=0; mmpos=-1;
        const uint64_t LO=0x5555555555555555ULL;
        auto acc=[&](uint64_t x,int blk){
            if(!x) return;
            const uint64_t nz=((x>>1)|x)&LO;
            mm+=__builtin_popcountll(nz);
            if(mm==1){
                // highest set 2-bit group = earliest differing base, because
                // w32 packs the first base in the most significant bits
                mmpos=blk+(__builtin_clzll(nz)>>1);
            }
        };
        while(L>=32){ acc(w32(pa)^w32(pb),base); if(mm>2) return 3;
                      pa+=32; pb+=32; L-=32; base+=32; }
        if(L){ const uint64_t m=~0ULL<<(2*(32-L)); acc((w32(pa)^w32(pb))&m,base); }
        return mm>2?3:mm;
    };
    // unpack read i into buf (caller guarantees >= rlen[i] bytes)
    auto runp=[&](uint32_t i,char* buf){
        const uint32_t L=rlen[i]; const uint64_t base=woff[i]*32ULL;
        for(uint32_t j=0;j<L;++j) buf[j]="ACGT"[(w32(base+j)>>62)&3ULL];
    };
    auto rappend=[&](std::string& dst,uint32_t i,uint32_t from){
        char b_[1024]; runp(i,b_); dst.append(b_+from,rlen[i]-from);
    };

    // ── STAGE 42: prefix-containment removal, generalising dedup ─────────────
    // Real variable-length data comes from 3' quality trimming, so a trimmed
    // read is a strict PREFIX of its untrimmed twin. Exact-hash dedup cannot see
    // that, and measured on 30% of reads trimmed to 130-149 bases the cost is
    // severe: unique reads 851,275 -> 918,598, both-side-overlapped 69.8% ->
    // 40.0%, survivors 12,065 -> 57,230, literal +24%.
    //
    // Myers' string graph (2005) and SGA (2012) both drop contained reads before
    // building, because a contained read carries nothing its container does not.
    // The important property here is that containment SUBSUMES dedup: at fixed
    // length one read contains another only when they are equal, so this
    // degenerates to exactly the existing behaviour and the fixed-length results
    // are unchanged -- which is checked, not assumed.
    //
    // Sorting lexicographically puts a prefix immediately before its extensions,
    // so one linear scan finds every prefix containment. A contained read is then
    // handled exactly like a duplicate: its originals are redirected to the
    // container (offset 0, so no extra bookkeeping) and it leaves the read set.
    std::vector<uint8_t> contained(n,0);
    size_t n_contained=0;
    // At constant read length one read can only contain another by equalling it,
    // which dedup has already handled -- so the whole pass is guaranteed to find
    // nothing and is pure cost. Measured, it was 0.31 s inside the load stage
    // (0.56 -> 0.87 s), which is 8% of the entire run spent proving a negative.
    bool varlen=false;
    for(uint32_t i=1;i<n;++i) if(rlen[i]!=rlen[0]){ varlen=true; break; }
    if(varlen){
        std::vector<uint32_t> ord(n);
        for(uint32_t i=0;i<n;++i) ord[i]=i;
        auto pfxcmp=[&](uint32_t x,uint32_t y)->int{     // compare min(len) bases
            uint32_t L=(rlen[x]<rlen[y])?rlen[x]:rlen[y];
            uint64_t px=woff[x]*32ULL, py=woff[y]*32ULL;
            while(L>=32){ const uint64_t a=w32(px),b=w32(py);
                          if(a!=b) return a<b?-1:1; px+=32; py+=32; L-=32; }
            if(L){ const uint64_t m=~0ULL<<(2*(32-L));
                   const uint64_t a=w32(px)&m,b=w32(py)&m;
                   if(a!=b) return a<b?-1:1; }
            return 0;
        };
        std::sort(ord.begin(),ord.end(),[&](uint32_t x,uint32_t y){
            const int c=pfxcmp(x,y);
            return c ? c<0 : rlen[x]<rlen[y];        // shorter first among prefixes
        });
        // a read is contained when the next entry shares all of its bases
        std::vector<uint32_t> alias(n,UINT32_MAX);
        for(size_t i=0;i+1<ord.size();++i){
            const uint32_t a=ord[i], b=ord[i+1];
            if(rlen[a]<=rlen[b] && pfxcmp(a,b)==0 && rlen[a]<rlen[b]){
                contained[a]=1; alias[a]=b; ++n_contained;
            }
        }
        // a container may itself be contained, so follow the chain to the end
        for(uint32_t i=0;i<n;++i){
            if(!contained[i]) continue;
            uint32_t t=alias[i]; int guard=0;
            while(t!=UINT32_MAX && contained[t] && ++guard<64) t=alias[t];
            alias[i]=t;
        }
        // CORRECTNESS: containment folds a SHORTER read into a longer one
        // (the condition above requires rlen[a] < rlen[b]). read_lengths was
        // indexed by unique id, so after aliasing a contained read inherited
        // its CONTAINER's length and decoded one or more bases too long. Every
        // regression dataset is fixed-length, so this never fired there; it
        // corrupted 27.6% of reads on the first variable-length file tested
        // (SARS-CoV-2 amplicon). Capture each ORIGINAL read's true length here,
        // before aliasing, and store lengths per original read instead.
        // origlen was captured at load time, before any aliasing, so the
        // rewrite below cannot corrupt it.
        for(uint32_t& u:orig2uid) if(u<n && contained[u] && alias[u]!=UINT32_MAX) u=alias[u];
        if(n_contained)
            fprintf(stderr,"contained reads folded into their container: %zu\n",n_contained);
    }
    // ord and alias are n x 4 bytes each and go out of scope above, but glibc
    // keeps the arena -- measured, the containment pass cost 12 MB of peak
    // (193,736 -> 205,844 KB) purely by being retained. Same lesson as stage 41.
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
    fprintf(stderr,"reads in=%zu N-reads(stored separately)=%zu unique=%u maxlen=%u\n",n_in,n_filt,n,Lmax);
    lap("load+filter+dedup");

    // `admit` selects which reads take part; hoisted above the index because
    // buildPref filters on it when rebuilding for the survivor pass.
    std::vector<uint8_t> admit(n,1);
    for(uint32_t i=0;i<n;++i) if(contained[i]) admit[i]=0;   // stage 42

    // ── STAGE 43: the prefix index, flat ─────────────────────────────────────
    // This was still an unordered_map<uint64,vector<uint32>> -- the exact
    // structure stages 18 and 40 replaced elsewhere for large wins, left in the
    // one place it is hit hardest: rounds 1 and 2 probe it 20.9M times and cost
    // 2.21 s of a 4.7 s run. 851k distinct keys each carry a node header plus a
    // separately allocated vector, ~95 MB, and every probe is a pointer chase.
    //
    // Flat instead: (key32 << 32) | rid in one sorted array, with an
    // open-addressing table pointing at each key's first entry. Sorting the whole
    // word sorts by key because the key is in the high bits, and rid ascends
    // within a key -- so candidates are visited in the same order the map's
    // insertion-ordered vector gave, and the output is unchanged. ~15 MB.
    //
    // The 32-bit key is exact at the default SW 16 and partial above it, which is
    // safe: a collision only proposes a candidate, and rcmp verifies every
    // candidate against the full L bases anyway.
    std::vector<uint64_t> pent;
    std::vector<uint32_t> ptab;
    size_t pmask=0;
    auto pmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
    std::vector<uint16_t> pext;
    auto buildPref=[&](bool useAdmit){
        pent.clear(); pent.shrink_to_fit(); pent.reserve(n);
        hugehint((void*)pent.data(), n*sizeof(uint64_t));
        for(uint32_t i=0;i<n;++i){
            if(useAdmit && !admit[i]) continue;
            if(rlen[i]<SW) continue;
            pent.push_back(((uint64_t)(uint32_t)rseed(i,0,SW)<<32)|(uint64_t)i);
        }
        std::sort(pent.begin(),pent.end());
        // ── SEQUENTIAL CANDIDATE DISCRIMINATOR ──────────────────────────────
        // MEASURED BOTTLENECK. Instrumenting the sweep on 3M human reads:
        //     probes 94,268,287   seed_hit 16.1%   cand_examined 371,742,373
        // i.e. 24.5 candidates are walked per seed hit, and each one touches
        // admit[b], rlen[b] and rcmp() at a RANDOM read index -- a cache miss
        // apiece, at IPC 0.62. The failed lookups are not the cost; the
        // candidate walk is.
        //
        // pext[q] holds the 8 bases FOLLOWING each entry's seed, stored in
        // index order so the walk reads it sequentially -- the same cache line
        // that already carries pent[q]. A candidate whose bases [32,40) differ
        // from the query's cannot pass rcmp, so it is rejected without ever
        // touching the read store.
        //
        // OUTPUT-PRESERVING, and only where it is valid: the check applies
        // only when L >= SW+8, because below that the extra bases lie beyond
        // the overlap being tested and comparing them would reject real links.
        pext.assign(pent.size(), 0xFFFFu);
        for(size_t i=0;i<pent.size();++i){
            const uint32_t b=(uint32_t)(pent[i]&0xFFFFFFFFULL);
            if(rlen[b] >= SW+8) pext[i]=(uint16_t)(rseed(b,SW,8)&0xFFFFu);
        }
        size_t sz=1; while(sz<pent.size()*2+1) sz<<=1;
        pmask=sz-1; hugefill(ptab,sz,UINT32_MAX);
        for(size_t i=0;i<pent.size();++i){
            const uint32_t k=(uint32_t)(pent[i]>>32);
            if(i && (uint32_t)(pent[i-1]>>32)==k) continue;
            size_t h=pmix(k)&pmask;
            while(ptab[h]!=UINT32_MAX) h=(h+1)&pmask;
            ptab[h]=(uint32_t)i;
        }
        hugehint(pent.data(), pent.size()*sizeof(pent[0]));
        hugehint(ptab.data(), ptab.size()*sizeof(ptab[0]));
    };
    auto pfind=[&](uint64_t key)->uint32_t{
        const uint32_t k=(uint32_t)key;
        size_t h=pmix(k)&pmask;
        while(ptab[h]!=UINT32_MAX){
            if((uint32_t)(pent[ptab[h]]>>32)==k) return ptab[h];
            h=(h+1)&pmask;
        }
        return UINT32_MAX;
    };
    buildPref(false);
    lap("prefix seed index");

    // ── longest verified suffix-prefix partner for each read, one pass ───────
    // ── TRUE descending-length sweep ────────────────────────────────────────
    // Storing one best partner per read is what broke this: sorted by overlap
    // and committed, a read whose single candidate was already consumed links
    // to nothing, and every unlinked read starts a chain that costs a full read
    // length. Measured, that was 78,327 chains x 150 b = 11.7 MB of pure chain
    // heads, against PgRC2's ~5,100 chains on the same read count -- the whole
    // gap, and nothing to do with tolerance or strandedness.
    //
    // The sweep instead revisits EVERY still-open end at every overlap length,
    // longest first, so a read whose best partner is taken still gets its
    // next-best one length down. Naively that re-reads L bytes per read per
    // length (O(n*L^2), the 40 s version). Here each open tail carries a rolling
    // 32-base seed: its suffix offset grows by exactly one as L drops by one, so
    // the seed advances with a single shift-or, and the L-byte memcmp runs only
    // after an exact 64-bit seed hit.
    std::vector<uint64_t> ppos;           // pg position of each unique read (stage 22)
    std::vector<uint8_t>  prc;            // stage 46: strand of that position
    size_t mm_total=0, mm_reads=0; std::vector<size_t> mm_hist(256,0);
    // STAGE 28. Best match, not first acceptable.
    // Stage 27 measured what first-acceptable costs: 10.74 mismatches per placed
    // read against PgRC2's 1.39. Every one of those has to be stored, so at
    // their coded rate (265,900 B for 1,369,413 mismatches = 1.55 bits each)
    // first-acceptable was quietly spending ~260 KB more than they do -- over
    // half the order-coding win, in a stream this progression never counted.
    //
    // PgRC2 re-scans and keeps the best placement (ReadsMatchers.cpp:315-328).
    // The mechanism that makes it cheap is that the current best doubles as the
    // early-exit bound: a candidate is abandoned the moment it is worse than
    // what the read already has, so better placements also mean less work.
    std::vector<uint8_t> readMM;          // 255 = unplaced
    // ── GRAPH HARVEST -- REFUTED BY ITS OWN MEASUREMENT, kept instrumented ──
    //
    // THE IDEA: the overlap graph is already computed inside the sweep below
    // and thrown away. For each open tail `a` at level L, cand[i*CCAP +
    // 0..ccnt[i]-1] is a's OUT-EDGE SET -- every read whose first L bases are
    // verified equal (rcmp) to a's last L bases -- and the commit keeps the
    // first untaken entry and discards the rest. Recording the discarded
    // alternative looked free, and structurally matched DiscoSNP++'s own
    // seeding step (Bubble.cpp:309: successors.size()>=2, then every pair).
    //
    // THE MEASUREMENT (HG002 r2, 75,115 reads, 2026-09-04):
    //     commits=57,553   ccnt==1: 54,966 (95.5%)   ccnt==2: 2,453 (4.3%)
    //     surviving branches: 85     of which 73 (86%) at L = 100% of Lmax
    // Against 424 truth het-SNVs and 95 truth het-indels in the same window.
    //
    // WHY IT CANNOT WORK, and this is structural, not tuning:
    //   1. `rcmp` is an EXACT match. Heterozygosity means the two haplotypes
    //      DIFFER. Two haplotype reads can therefore never both exactly match
    //      the same suffix unless the variant lies outside the overlap window.
    //      A het site does not produce a tie here -- it produces exactly one
    //      surviving candidate, which is why ccnt==1 is 95.5% of commits.
    //   2. Of the ties that do occur, 86% sit at L = Lmax, and per the comment
    //      at the head of the level loop below, "At L = rlen[a] the offset is 0
    //      ... rcmp compares the two reads in full -- exactly the duplicate
    //      test." They are duplicate reads, not branches.
    //   3. 2,368 of the 2,453 two-candidate commits had the alternative
    //      already claimed (prv[other]!=NONE): two chains meeting, not a branch.
    //
    // So this harvests ~12 non-duplicate ties per 519 truth variants. It
    // cannot serve as a calling channel (accuracy), and it cannot replace
    // rc_reads/pkidx as a candidate source (RAM/time) because it finds ~0.15%
    // as many loci. Both halves of docs/GRAPH_HARVEST_EXACT_PLAN.md are dead.
    //
    // KEPT (not deleted) because it is CAPS_CALL-gated, verified bit-identical
    // on the Claim 1 path (ERR5181310 ARCHIVE_TOTAL=352997 and every stream
    // byte-for-byte), costs one predictable branch when off, and the printed
    // histogram is the evidence for the paragraph above. Full write-up:
    // docs/GRAPH_HARVEST_REFUTED.md.
    struct GEdge { uint32_t a, b_taken, b_alt, L; };
    std::vector<GEdge> g_branch;
    // instrumentation: histogram of ccnt at the moment of commit (CCAP=8, so 0..8)
    size_t g_cc_hist[9] = {0};
    size_t g_commits = 0;
    // L distribution for 2-way ties, bucketed as a percentage of Lmax: a tie at
    // a long overlap is strong evidence (exact rcmp over many bases); a tie
    // near the sweep floor could be an error-induced collision. Bucketing lets
    // the threshold in any later filter be derived from measurement rather
    // than guessed (standing rule 1: formulas over measured properties).
    size_t g_L_bucket[11] = {0};
    // ── NEAR-MISS HARVEST accumulators (see rcmp_mm above for the rationale) ─
    // Declared out here, beside g_branch, because they are OpenMP reduction
    // targets inside sweep() AND are reported after round 2 returns.
    size_t nm1=0, nm2=0;
    // alt = the base the OTHER haplotype's read carries at that position (2-bit)
    struct NMObs { uint32_t read_a, pos_in_a; uint8_t alt; };
    std::vector<NMObs> nm_obs;
    // Minimum verified overlap for a near-miss to count, as a FRACTION of Lmax
    // rather than a fixed constant (standing rule 1): at Lmax=148 this is 74
    // bases, and one mismatch across 74+ otherwise-exact bases is not a chance
    // collision. Overridable so the threshold can be swept, not guessed.
    uint32_t nm_minL=Lmax/2;
    if(const char* e=getenv("CAPS_NM_MINL")) nm_minL=(uint32_t)atoi(e);
    std::vector<uint32_t> nxt,prv,ovl,ch_h,ch_t;
    hugefill(nxt,n,(uint32_t)NONE); hugefill(prv,n,(uint32_t)NONE);
    hugefill(ovl,n,(uint32_t)0);    hugefill(ch_h,n,(uint32_t)0);
    hugefill(ch_t,n,(uint32_t)0);
        hugehint(nxt.data(), nxt.size()*sizeof(nxt[0]));
        hugehint(prv.data(), prv.size()*sizeof(prv[0]));
        hugehint(ovl.data(), ovl.size()*sizeof(ovl[0]));
        hugehint(ch_h.data(), ch_h.size()*sizeof(ch_h[0]));
        hugehint(ch_t.data(), ch_t.size()*sizeof(ch_t[0]));
    std::vector<uint64_t> seed; hugefill(seed,n,(uint64_t)0);
    std::vector<uint8_t>  ok(n,0);
    std::vector<uint32_t> tails;
    std::vector<uint32_t> open_tails;   // still-open tails, independent of level
    size_t links=0, probes=0;
    // `admit` selects which reads take part; everything else is reset per round.
    // Round 1 only LABELS reads (which ones sit interior to a chain), so its
    // links cost nothing -- a link there buys classification, not bytes. It
    // therefore runs at the widest overlap the seed allows, while round 2, which
    // actually emits sequence, keeps the caller's floor. Running both at the
    // same conservative floor was leaving reads unlabelled that would have
    // classified fine: interior fraction 66% against PgRC2's 83%, which starves
    // the main pseudogenome and leaves the mapping stage with less to hit.
    uint32_t sweep_minov=MINOV;
    static size_t g_sd_hit=0, g_cand=0, g_probe=0;
    auto sweep=[&](){
        links=0; probes=0; g_sd_hit=0; g_cand=0; g_probe=0;
        // Round 1's result is discarded wholesale (nxt/prv/ovl are refilled on
        // the next call, and round 1 exists only to compute `admit`), so its
        // harvest must be discarded too -- only round 2's survives, which is
        // correct: round 2's nxt/ovl are what build pg and g_contig_spans.
        g_branch.clear(); g_commits=0; nm_obs.clear(); nm1=0; nm2=0;
        for(size_t z=0;z<9;++z) g_cc_hist[z]=0;
        for(size_t z=0;z<11;++z) g_L_bucket[z]=0;
        std::fill(nxt.begin(),nxt.end(),NONE); std::fill(prv.begin(),prv.end(),NONE);
        std::fill(ovl.begin(),ovl.end(),0);    std::fill(ok.begin(),ok.end(),0);
        for(uint32_t i=0;i<n;++i){ ch_h[i]=i; ch_t[i]=i; }
        tails.clear(); tails.reserve(n);
        open_tails.clear(); open_tails.reserve(n);
        for(uint32_t i=0;i<n;++i) if(admit[i]&&rlen[i]>=sweep_minov) open_tails.push_back(i);
        tails = open_tails;
    // STAGE 19. The sweep splits cleanly into a search and a commit.
    //
    // Everything expensive here -- the seed roll, the hash probe, and the
    // L-byte memcmp against every candidate -- reads only state that is fixed
    // for the whole of one L level: reads, pref, admit, and each tail's own
    // seed/ok slot. The only conditions that change WITHIN a level are
    // `prv[b]!=NONE` (b already claimed) and `ch_h[a]==b` (would close a cycle),
    // and both are a couple of array loads.
    //
    // So the search runs in parallel and the commit stays serial in the original
    // tail order. Because the parallel phase records EVERY candidate that passes
    // the static tests, in bucket order, the serial phase can apply the dynamic
    // tests to that same list and take the same first survivor the fully serial
    // loop would have taken. This is not an approximation of the old behaviour;
    // it is the old behaviour with the memcmp moved off the critical path.
    const uint32_t CCAP=8;               // candidates retained per tail per level
    std::vector<uint32_t> cand; std::vector<uint8_t> ccnt;
    unsigned NT=std::thread::hardware_concurrency(); if(!NT) NT=1;
    // A1: one team for the whole sweep, not one per overlap level.
    // The level loop runs 134 times in round 1 and 110 in round 2, and sweep()
    // is called three times, so entering the parallel region inside the loop
    // paid ~600 fork/joins. perf attributed 21.3% of all cycles to libgomp and
    // 16.0% to the kernel (against PgRC2's 1.7% kernel) -- that is team setup
    // and the IPIs it triggers, not compression. Hoisting keeps every write
    // disjoint exactly as before: the search is per-tail independent and the
    // commit below still walks i in order.
    bool sweep_done=false;
    // SWEEP PHASE TIMING (CAPS_SWEEP_TIMING=1). The level loop is
    // serial-prep -> parallel-search -> serial-commit, and only the middle
    // third is parallel. With 110 levels in round 2 the two serial passes are
    // paid 110 times over millions of reads, so Amdahl's law -- not the TLB --
    // is the candidate explanation for round 2 refusing to scale with threads.
    // Diagnostic only: unset, this is one branch per level and the archive is
    // bit-identical.
    double tt_prep=0, tt_par=0, tt_com=0; size_t nlev=0, tt_w=0, tt_open=0;
    const bool SWTIME = getenv("CAPS_SWEEP_TIMING")!=nullptr;
    size_t pr_total=0; size_t sd_hit=0; size_t cand_seen=0;
    #pragma omp parallel
    {
    const bool me0 = (omp_get_thread_num()==0);
    double ta=0, tb=0, tc=0, td=0;
    // Start at Lmax, not Lmax-1. A suffix-prefix overlap of exactly the read
    // length IS an exact duplicate, and it is the only length at which one can
    // appear: two identical 251-base reads do NOT overlap at L=250, because
    // that compares A[1..251) against B[0..250), equal only for a periodic read.
    // Excluding L=Lmax therefore made duplicates invisible to chaining, which is
    // why they had to be handled by a separate pre-assembly dedup pass -- and
    // that pass costs an orig2uid array over every ORIGINAL read, which is why a
    // threshold was then needed to decide between two bad options.
    //
    // PgRC2 needs no such pass and no such threshold: a duplicate is simply a
    // 100%-length overlap in their chain mechanism, contributing zero bases.
    // Measured on S. acidocaldarius, they remove 48,694 duplicates this way
    // while our dedup threshold (9.4% < 15%) turned dedup OFF and carried all
    // 48,694 into the pseudogenome.
    //
    // At L = rlen[a] the offset is 0, so the seed is the read's own prefix and
    // rcmp compares the two reads in full -- exactly the duplicate test. The
    // existing ch_h[a]==b guard already prevents a 2-cycle.
    for(uint32_t L=Lmax; L>=sweep_minov && L>=SW; --L){
        if(SWTIME&&me0) ta=omp_get_wtime();
        #pragma omp single
        {
        // Two different filters were being applied to the same array, and only
        // one of them is permanent.
        //
        //   nxt[a]!=NONE   the tail was extended -- it is done for good
        //   rlen[a]<L      the read is shorter than THIS level's overlap
        //
        // The second is not permanent: L decreases every iteration, so a read
        // too short at L becomes eligible at L-1. Dropping it from `tails`
        // removed it from every remaining level too. On constant-length reads
        // the test never fires and nothing is lost, which is why all 7
        // benchmark datasets are byte-identical either way -- they are all
        // constant length. On variable-length reads it discards nearly the
        // whole file at the first level: SARS-CoV-2 has 103 distinct lengths
        // from 31 to 221, so at L=221 every read but the longest was dropped
        // permanently.
        //
        // `open_tails` now holds the permanent filter and the per-level length
        // test selects into `tails` without destroying anything.
        size_t ow=0;
        for(uint32_t a:open_tails){
            if(nxt[a]!=NONE) continue;                       // extended: done for good
            open_tails[ow++]=a;
        }
        open_tails.resize(ow);
        tails.clear();
        for(uint32_t a:open_tails){
            if(rlen[a]<L) continue;                          // not at THIS level; keep for later
            const uint32_t off=(uint32_t)rlen[a]-L;          // grows by 1 per length
            if(off+SW>rlen[a]) continue;
            tails.push_back(a);
        }
        const size_t w=tails.size(); (void)w;
        // `break` cannot leave an OpenMP structured block, so the exit
        // condition becomes a shared flag tested after the implicit barrier.
        if(tails.empty()) sweep_done=true;
        else {
            if(cand.size()<w*CCAP) cand.resize(w*CCAP);
            ccnt.assign(w,0);
            pr_total=0; sd_hit=0; cand_seen=0;
        }
        }   // end single (implicit barrier: every thread now sees tails/cand)
        if(sweep_done) break;
        if(SWTIME&&me0) tb=omp_get_wtime();
        const size_t w = tails.size();

        // ── parallel: roll each seed, probe, keep candidates passing the
        //    static tests, in bucket order ──────────────────────────────────
        //
        // STAGE 45. This block runs once per overlap level -- 134 levels in
        // round 1, 110 in round 2 -- so anything paid per level is paid 244
        // times. Three such costs were here, none of them the actual search:
        //
        //   cand.assign(w*CCAP,NONE) memset 27 MB at the widest level and
        //   several MB at most others, to write a sentinel that is never read:
        //   the serial phase reads cand[i*CCAP+c] only for c < ccnt[i]. Only
        //   ccnt has to be cleared, because a tail whose seed misses the index
        //   `continue`s before assigning it. So the buffer is grown once to the
        //   high-water mark and never cleared again.
        //
        //   12 std::threads created and joined per level -- 2,928 spawns for a
        //   run whose whole budget is under four seconds. OpenMP keeps one pool
        //   alive across all of them, and schedule(static) over [0,w) hands each
        //   thread the same contiguous span the manual split did, so every tail
        //   is searched by exactly the work it was before.
        //
        //   ++pr[t] on a vector<size_t> puts 12 counters in two cache lines and
        //   increments them 20.9M times from 12 cores: every increment steals
        //   the line back. A local counter added once at the end is the same
        //   number with none of the traffic.
        //
        // The output is unchanged by construction -- the search is per-tail
        // independent and the serial commit below still walks i in order.
        {
            // num_threads(T) with T=(w<4096?1:NT) existed only to dodge fork
            // cost on small levels; with the team already open there is no
            // fork to dodge, and the writes were disjoint either way.
            #pragma omp for schedule(static) reduction(+:pr_total,nm1,nm2,sd_hit,cand_seen)
            for(long long ii=0;ii<(long long)w;++ii){
                const size_t i=(size_t)ii;
                const uint32_t a=tails[i];
                const uint32_t off=(uint32_t)rlen[a]-L;
                // packed reads are always ACGT, so the seed always exists
                seed[a]=rseed(a,off,SW); ok[a]=1;
                ++pr_total;
                const uint32_t pix=pfind(seed[a]);
                if(pix==UINT32_MAX) continue;
                ++sd_hit;
                const uint32_t pk=(uint32_t)seed[a];
                const uint16_t aext = (off+SW+8 <= rlen[a])
                                    ? (uint16_t)(rseed(a,off+SW,8)&0xFFFFu) : 0xFFFFu;
                uint8_t c=0;
                for(uint32_t q=pix;q<pent.size()&&(uint32_t)(pent[q]>>32)==pk;++q){
                    ++cand_seen;
                    // Sequential reject before any random read access.
                    if(L >= SW+8 && pext[q]!=0xFFFFu && pext[q]!=aext) continue;
                    const uint32_t b=(uint32_t)(pent[q]&0xFFFFFFFFULL);
                    if(b==a) continue;
                    if(!admit[b]) continue;                  // excluded reads are leftovers,
                                                             // never chain members -- without
                                                             // this they get emitted twice
                    if(rlen[b]<L) continue;
                    if(!rcmp(a,off,b,L)){
                        // ── NEAR-MISS HARVEST (CAPS_CALL only) ──────────────
                        // The het partner dies on the line above. It shares
                        // a's seed (that is why it is in this bucket at all)
                        // but differs somewhere in the verified overlap. If it
                        // differs by EXACTLY ONE base, that base is a
                        // candidate heterozygous site and this read is the
                        // alternate haplotype. No extra probe, no extra index:
                        // `b` is already in hand at the instant it is rejected.
                        //
                        // Only counted at long overlaps: a 1-mismatch hit over
                        // 100+ verified bases is strong evidence; the same over
                        // 20 bases is a coincidence. NMIN is expressed as a
                        // fraction of Lmax, not a constant, per standing rule 1.
                        if(CAPS_CALL && L>=nm_minL){
                            int mmpos=-1;
                            const int mm=rcmp_mm1(a,off,b,L,mmpos);
                            if(mm==1){
                                ++nm1;
                                // Record (tail read, offset of the differing
                                // base within read a). ppos[a] is not known
                                // until chains are emitted, so translation to
                                // pg coordinates happens there -- the same
                                // deferral g_contig_spans already uses.
                                // ~44k pushes across the entire run, against
                                // 20.9M probes: a critical section here is far
                                // cheaper than per-thread buffers plus a merge,
                                // and cannot race.
                                // b's base at the differing position: b's
                                // overlap starts at its own offset 0, so the
                                // differing base is b[mmpos].
                                const uint8_t altb=(uint8_t)((w32(woff[b]*32ULL+(uint32_t)mmpos)>>62)&3ULL);
                                #pragma omp critical(nmobs)
                                nm_obs.push_back({a,(uint32_t)(off+mmpos),altb});
                            }
                            else if(mm==2) ++nm2;
                        }
                        continue;
                    }
                    cand[i*CCAP+c]=b;
                    if(++c==CCAP) break;
                }
                ccnt[i]=c;
            }
            #pragma omp single
            { probes+=pr_total; g_sd_hit+=sd_hit; g_cand+=cand_seen; g_probe+=pr_total; }
        }
        if(SWTIME&&me0) tc=omp_get_wtime();

        // ── serial: same order, same first survivor ───────────────────────
        // A retained list that filled to CCAP may have been truncated, so if
        // every entry in a full list turns out to be taken, that tail is
        // re-probed serially. Without this the cap silently drops the 9th-and-
        // later candidate and the run loses links -- measured at CCAP=8, 14 of
        // 673,334 links and 1,288 bytes of literal. Rare enough to cost nothing,
        // and it makes the result identical to the fully serial sweep rather
        // than merely close to it.
        #pragma omp single
        {
        for(size_t i=0;i<w;++i){
            const uint32_t a=tails[i];
            bool done=false;
            for(uint8_t c=0;c<ccnt[i];++c){
                const uint32_t b=cand[i*CCAP+c];
                if(prv[b]!=NONE||ch_h[a]==b) continue;       // taken, or would cycle
                nxt[a]=b; prv[b]=a; ovl[a]=L;
                uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h; ++links;
                // ── graph harvest: record the discarded alternative ─────────
                // Placed here, inside the existing `omp single`, so it is
                // already serial: no race, no critical section, no atomic.
                // Reads only values the commit itself just used; writes only
                // to g_* which nothing else touches. When CAPS_CALL is unset
                // (every Claim 1 run) this is a single predictable branch and
                // the emitted archive is bit-identical.
                if(CAPS_CALL){
                    ++g_commits;
                    g_cc_hist[ccnt[i]<9?ccnt[i]:8]++;
                    if(ccnt[i]==2){
                        // the other entry of the pair, whichever slot it is in
                        const uint32_t other=cand[i*CCAP+(c==0?1:0)];
                        // It must be a genuine unclaimed alternative: if it is
                        // already in a chain this is two chains meeting, not a
                        // branch; ch_h guard mirrors the commit's own cycle test.
                        if(other!=b && prv[other]==NONE && ch_h[a]!=other){
                            g_branch.push_back({a,b,other,L});
                            const uint32_t pct = Lmax? (uint32_t)((uint64_t)L*10ULL/Lmax) : 0;
                            g_L_bucket[pct<11?pct:10]++;
                        }
                    }
                }
                done=true; break;
            }
            if(done||ccnt[i]<CCAP||!ok[a]) continue;         // list was complete
            const uint32_t off=(uint32_t)rlen[a]-L;
            const uint32_t pix=pfind(seed[a]);
            if(pix==UINT32_MAX) continue;
            const uint32_t pk=(uint32_t)seed[a];
            uint8_t seen=0;
            for(uint32_t q=pix;q<pent.size()&&(uint32_t)(pent[q]>>32)==pk;++q){
                const uint32_t b=(uint32_t)(pent[q]&0xFFFFFFFFULL);
                if(b==a) continue;
                if(!admit[b]) continue;
                if(rlen[b]<L) continue;
                if(!rcmp(a,off,b,L)) continue;
                if(seen++<CCAP) continue;                    // already tried above
                if(prv[b]!=NONE||ch_h[a]==b) continue;
                nxt[a]=b; prv[b]=a; ovl[a]=L;
                uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h; ++links;
                break;
            }
        }
        }   // end single (commit)
        if(SWTIME&&me0){ td=omp_get_wtime();
            tt_prep+=tb-ta; tt_par+=tc-tb; tt_com+=td-tc;
            ++nlev; tt_w+=w; tt_open+=open_tails.size(); }
    }       // end level loop
    }       // end parallel (one team for the whole sweep)
    if(SWTIME){
        const double tot=tt_prep+tt_par+tt_com;
        fprintf(stderr,"[sweep] levels=%zu  prep=%.2fs (%.1f%%)  par=%.2fs (%.1f%%)  "
                       "commit=%.2fs (%.1f%%)  total=%.2fs  serial=%.1f%%  "
                       "mean_tails=%zu mean_open=%zu\n",
            nlev, tt_prep, tot?100*tt_prep/tot:0, tt_par, tot?100*tt_par/tot:0,
            tt_com, tot?100*tt_com/tot:0, tot,
            tot?100*(tt_prep+tt_com)/tot:0,
            nlev?tt_w/nlev:0, nlev?tt_open/nlev:0);
    }
    };

    // ── round 1: label, do not build ────────────────────────────────────────
    // The first sweep exists only to find which reads sit INTERIOR to a chain --
    // overlapped on both sides. Those are the reads that tile the genome well.
    // Reads that end up at a chain end, or unlinked, overlap poorly; keeping
    // them in the pseudogenome is what drags the mean overlap down (measured
    // 119.9 against PgRC2's 133.7, where ideal tiling at this coverage is ~136).
    // They cost far less mapped into the finished sequence than assembled into
    // it, so round 2 rebuilds without them.
    // SW was doing two unrelated jobs: the width of the rolling seed, and the
    // minimum overlap round 1 will classify on. They do not have to be equal --
    // the seed is only a prefix filter, and the L-base memcmp behind it is what
    // decides. Measured, raising SW improves the main pg monotonically and
    // worsens the survivor pg monotonically (32 -> main 11,384,786 / survivor
    // 1,567,625; 16 -> 11,632,777 / 1,294,148), and the best half of each beats
    // PgRC2's total. Splitting the parameter is what makes both reachable at
    // once. Same shape as stage 17's split of acceptance from seed geometry.
    const uint32_t R1MINOV=(argc>7)?(uint32_t)atoi(argv[7]):SW;
    sweep_minov=(R1MINOV>=SW)?R1MINOV:SW;   // classifier floor, seed stays SW
    sweep();
    size_t both_sides=0;
    {
        std::vector<uint8_t> keep(n,0);
        for(uint32_t i=0;i<n;++i)
            if(nxt[i]!=NONE && prv[i]!=NONE){ keep[i]=1; ++both_sides; }
        admit.swap(keep);
    }
    fprintf(stderr,"round1: links=%zu both-sides-overlapped=%zu (%.1f%%)\n",
            links,both_sides,100.0*(double)both_sides/(double)n);
    lap("round 1 (division)");

    // ── A3: run every adaptive candidate off ONE shared prefix ──────────────
    // encode_adaptive.sh encodes the whole file once per (MAXMAP, MINOV) pair,
    // four times. But MINOV first takes effect at the round-2 sweep below, and
    // MAXMAP not until the mapping stage, so everything up to here -- load,
    // seed index, round 1 -- is identical for all four. Measured, that prefix
    // is 33.23 s of 133.94 s (24.8%), so four independent encodes repeat 100 s
    // of work that has exactly one answer.
    //
    // The candidates are run by forking here rather than by refactoring the
    // suffix into a re-runnable function. fork() gives each candidate an exact
    // copy-on-write snapshot of the post-prefix state, so no variable can leak
    // between candidates and the code path is bit-for-bit the one that runs
    // today -- the output is identical by construction rather than by review.
    // Children run one at a time, so peak RSS is unchanged.
    // ── TWO-LEVEL CANDIDATE FORK: state carried from level 1 to level 2 ─────
    // Level 1 forks per distinct MINOV; level 2, after `emit chains`, forks per
    // MAXMAP within that group. These live at main scope because the two fork
    // sites are ~450 lines apart.
    std::vector<uint32_t> g_l2_maxmaps;    // MAXMAP of each member of MY group
    std::vector<size_t>   g_l2_members;    // their GLOBAL candidate indices
    std::string           g_l2_base, g_l2_cwd;
    unsigned              g_l2_threads = 1;
    size_t                g_l2_conc    = 1;
    bool                  g_l2_active  = false;
    if(const char* cl = getenv("CANDIDATES")){
        std::vector<std::pair<uint32_t,uint32_t>> cands;   // (MAXMAP, MINOV)
        { const char* q=cl;
          while(*q){ unsigned a=0,b=0; if(sscanf(q,"%u:%u",&a,&b)==2) cands.push_back({a,b});
                     while(*q && *q!=',') ++q; if(*q==',') ++q; } }
        if(!cands.empty()){
            // ── GROUP BY MINOV: ROUND 2 IS A FUNCTION OF MINOV ALONE ────────
            //
            // The 8-point grid is 4 MAXMAP x 2 MINOV. MINOV is consumed by the
            // round-2 sweep; MAXMAP is not read until the mapping stage, ~300
            // lines later. So the four candidates sharing a MINOV compute a
            // BIT-IDENTICAL round-2 assembly and three of them throw it away.
            //
            // Not inferred from the code -- read off the production HG002 log,
            // where `round2:` prints exactly two distinct results across the
            // eight candidates:
            //     probes=165361129 links=10254332   x4   (MINOV=16)
            //     probes=173265221 links=10369827   x4   (MINOV=51)
            // and the per-candidate stage times were 785,785,785,782 s for one
            // group and 125,125,115,114 s for the other -- i.e. ~2,360 s of
            // recomputation of an assembly that was already in hand.
            //
            // So the fork is at the right PLACE (after round 1) but the wrong
            // DEPTH. Level 1 forks one child per distinct MINOV; that child
            // runs round 2 and emit chains once, with a real share of the
            // cores, and only then forks one grandchild per MAXMAP.
            //
            // OUTPUT-PRESERVING BY CONSTRUCTION: a grandchild's state at the
            // level-2 fork is exactly the state the corresponding candidate had
            // at the same point, because round 2 is a deterministic function of
            // (shared prefix, MINOV). It is also thread-count independent --
            // already established, since K=8/4/2 (1, 3 and 6 threads per
            // candidate) produced byte-identical archives.
            std::vector<uint32_t> gmin;                 // distinct MINOV, first-seen order
            std::vector<std::vector<size_t>> gidx;      // candidate indices per group
            for(size_t ci=0; ci<cands.size(); ++ci){
                size_t g=0; for(; g<gmin.size(); ++g) if(gmin[g]==cands[ci].second) break;
                if(g==gmin.size()){ gmin.push_back(cands[ci].second); gidx.push_back({}); }
                gidx[g].push_back(ci);
            }
            // fork() in a process that has live OpenMP threads deadlocks the
            // child: libgomp's internal locks are copied in whatever state they
            // held, and the child's first parallel region then waits forever on
            // threads that do not exist in it. POSIX only guarantees fork() in a
            // threaded program when the child execs. Measured here as a child
            // hung at 0.0% CPU. Tearing the runtime down first makes the process
            // genuinely single-threaded, so the snapshot is safe; the team is
            // rebuilt on demand inside each child.
            omp_pause_resource_all(omp_pause_hard);
            const std::string base = getenv("ARCHIVE") ? getenv("ARCHIVE") : "out.arcs2";
            bool child=false; size_t mine=0;
            // GOLDEN-SECTION over MAXMAP, when GSEARCH is set.
            //
            // Licensed by a measured structural property, not by tuning: the fine
            // sweep found total_size(MAXMAP) strictly descending then ascending on
            // every dataset tested (halo min at 8, ecoli at 20, sulfo at 45, no
            // violations). A unimodal objective admits golden-section search,
            // which reaches the true optimum in ~6 probes.
            //
            // Why bother when the fixed grid already lands within 0.4%: MAXMAP
            // swings archive size by 15-20% (halo 2.62 MB -> 3.14 MB across the
            // sweep). The grid is close on the datasets it was chosen against and
            // carries NO guarantee on unseen data. The search does.
            //
            // Each probe costs only a suffix, because A3 forks after the shared
            // prefix. Bounds come from the observed optima (L/19 .. L/5.6) with
            // margin; they bound the SEARCH, they do not choose the answer.
            bool gs_child=false;
            if(getenv("GSEARCH")){
                uint32_t glo = Lmax/25; if(glo<6) glo=6;
                uint32_t ghi = Lmax/4;  if(ghi<glo+4) ghi=glo+4;
                std::map<uint32_t,size_t> seen;
                uint32_t probe=0; int guard=0;
                auto need=[&](uint32_t m)->bool{ return !seen.count(m); };
                while(ghi>glo+2 && guard++<10 && !child){
                    uint32_t a = glo + (uint32_t)((ghi-glo)*0.382);
                    uint32_t b = glo + (uint32_t)((ghi-glo)*0.618);
                    if(a==b) b=a+1;
                    for(uint32_t m : {a,b}){
                        if(!need(m)) continue;
                        pid_t pid=fork();
                        if(pid<0){ perror("fork"); return 1; }
                        if(pid==0){ child=true; probe=m; break; }
                        int st=0; waitpid(pid,&st,0);
                        struct stat sb; std::string f=base+".cand"+std::to_string(m);
                        seen[m] = (stat(f.c_str(),&sb)==0) ? (size_t)sb.st_size : SIZE_MAX;
                    }
                    if(child) break;
                    if(seen[a] <= seen[b]) ghi=b; else glo=a;
                }
                if(!child){
                    size_t bs=SIZE_MAX; uint32_t bm=0;
                    for(auto& kv : seen) if(kv.second<bs){ bs=kv.second; bm=kv.first; }
                    char pcwd[4096]; if(!getcwd(pcwd,sizeof pcwd)) pcwd[0]='\0';
                    for(auto& kv : seen){
                        std::string f=base+".cand"+std::to_string(kv.first);
                        if(kv.first==bm) rename(f.c_str(), base.c_str()); else remove(f.c_str());
                        std::string d=f+".d";
                        if(DIR* dp=opendir(d.c_str())){
                            while(struct dirent* de=readdir(dp)){
                                if(de->d_name[0]=='.') continue;
                                std::string src=d+"/"+de->d_name;
                                if(kv.first==bm) rename(src.c_str(),(std::string(pcwd)+"/"+de->d_name).c_str());
                                else remove(src.c_str());
                            }
                            closedir(dp); rmdir(d.c_str());
                        }
                    }
                    fprintf(stderr,"  [gsearch] %zu probes, chose MAXMAP=%u -> %zu B\n",
                            seen.size(), bm, bs);
                    fprintf(stderr,"ARCHIVE_TOTAL=%zu\n",bs);
                    printf("ARCHIVE_TOTAL=%zu\n",bs);
                    return 0;
                }
                char cwdbuf[4096]; if(!getcwd(cwdbuf,sizeof cwdbuf)) cwdbuf[0]='\0';
                const std::string ab = (!base.empty()&&base[0]=='/')?base:(std::string(cwdbuf)+"/"+base);
                const std::string dd = ab+".cand"+std::to_string(probe)+".d";
                mkdir(dd.c_str(),0755);
                setenv("MAXMAP", std::to_string(probe).c_str(), 1);
                setenv("ARCHIVE",(ab+".cand"+std::to_string(probe)).c_str(),1);
                if(chdir(dd.c_str())!=0) perror("chdir");
                gs_child=true;
            }
            // ── CANDIDATES RUN CONCURRENTLY, BOUNDED BY MEASURED RAM ──────
            // This loop forked ONE candidate and waitpid()ed it before forking
            // the next, so the adaptive sweep was strictly serial. The
            // candidates are fully independent -- each writes its own
            // `.cand<N>` file and the parent then keeps whichever is smallest,
            // an order-independent choice -- so the serialisation bought
            // nothing. Measured on SARS-CoV-2 (471 MB input):
            //     1 candidate    10.42 s   557 MB peak   8,698,846 B
            //     4 candidates   39.63 s   557 MB peak   8,686,774 B
            // i.e. 3.8x the wall time for 0.14% smaller output, at IDENTICAL
            // peak RAM -- because only one child was ever resident.
            //
            // Running them concurrently spends RAM to buy that time back. The
            // concurrency is not a fixed number: it is derived from what the
            // machine actually has, so a small box still runs (serially if it
            // must) and a large one uses what it has.
            //
            //   per-child peak  ~= max(observed peak RSS, 1.5 x input bytes)
            //   K               = clamp(1, ncands, 0.60 x MemAvailable / peak)
            //
            // The 1.5x input term is the floor because a child continues to
            // grow after the fork; the observed peak covers the shared prefix
            // already resident. 60% of MemAvailable leaves headroom for the
            // page cache and the parent.
            //
            // OUTPUT IS UNCHANGED. Every child writes a distinct file and the
            // parent's selection is by size, so scheduling cannot alter which
            // archive wins. CAPS_CAND_SEQ=1 restores the serial behaviour for
            // A/B measurement.
            // Parent-side, before any fork: pay for names and quality ONCE, with
            // every core available. Inside the sweep each child holds only
            // NT/K threads, so the parallel trials inside the quality coder
            // would have no cores to use -- measured as a near-zero gain when
            // the coder was parallelised but left inside the candidates.
            if(!gs_child && !child){
                if(CAPS_NAMES && !g_NM_done && !g_input_path.empty()){
                    g_NM = nmc::encode_from_fastq(g_input_path.c_str()); g_NM_done = true;
                }
                if(CAPS_QUAL && !g_QL_done && !g_input_path.empty()){
                    g_QL = qlc::encode_from_fastq(g_input_path.c_str()); g_QL_done = true;
                }
                // TEAR THE RUNTIME DOWN AGAIN BEFORE FORKING.
                //
                // The teardown above (omp_pause_resource_all, ~100 lines up)
                // happens BEFORE this hoist, and the quality coder's trial
                // search is itself an OpenMP region -- so encoding here
                // RE-CREATES the thread pool, and the fork below then hits
                // exactly the deadlock that call exists to prevent: libgomp's
                // locks are copied in whatever state they held and the child
                // hangs at 0.0% CPU on its first parallel region.
                //
                // This is not theoretical. It was measured: the first version
                // of this hoist ran fine twice and then hung a verification
                // sweep for 13 minutes on the first dataset. A race that
                // sometimes passes is worse than a consistent failure, because
                // it would have stalled a 10-hour benchmark at random.
                if(g_NM_done || g_QL_done) omp_pause_resource_all(omp_pause_hard);
            }
            if(!gs_child){
                size_t K = cands.size();
                if(getenv("CAPS_CAND_SEQ")) K = 1;
                else {
                    size_t avail_mb = 0;
                    if(FILE* mi = fopen("/proc/meminfo","r")){
                        char k[64]; unsigned long v; char u[16];
                        while(fscanf(mi,"%63s %lu %15s",k,&v,u)>=2)
                            if(!strcmp(k,"MemAvailable:")){ avail_mb=v/1024; break; }
                        fclose(mi);
                    }
                    size_t in_mb = 0;
                    { struct stat sb; if(stat(argv[1],&sb)==0) in_mb=(size_t)sb.st_size/(1024*1024); }
                    size_t peak_mb = hwm_mb();
                    size_t per_child = peak_mb > (in_mb*3)/2 ? peak_mb : (in_mb*3)/2;
                    if(per_child < 64) per_child = 64;
                    size_t fit = avail_mb ? (avail_mb*3/5)/per_child : 1;
                    if(fit < 1) fit = 1;
                    if(fit < K) K = fit;
                    // K MUST DIVIDE THE CANDIDATE COUNT, or the last wave runs
                    // short-handed on a mostly idle machine.
                    //
                    // With N candidates run K at a time there are ceil(N/K)
                    // waves, and every child is pinned to P/K threads. If K
                    // does not divide N the final wave has fewer children than
                    // K, so it occupies (N mod K)*(P/K) cores and leaves the
                    // rest idle -- while still paying a full P/K-thread
                    // candidate's latency. N=4, P=12, K=3 is the worst case:
                    //     K=3 -> waves of 3 then 1, both at 4 threads = 2*T(4)
                    //     K=2 -> waves of 2 then 2, both at 6 threads = 2*T(6)
                    // Same number of waves, but every child gets 6 threads
                    // instead of 4, so K=2 STRICTLY BEATS K=3 here. Stepping
                    // down to the largest divisor of N that still fits memory
                    // is therefore never worse and is often better.
                    while(K > 1 && (cands.size() % K) != 0) --K;
                    if(const char* cc = getenv("CAPS_CAND_PAR")) { long q=atol(cc); if(q>0) K=(size_t)q; }
                    fprintf(stderr,"  [a3] candidates=%zu concurrency=%zu "
                                   "(avail=%zuMB per-child~%zuMB)\n",
                            cands.size(), K, avail_mb, per_child);
                }
                // DIVIDE THE CORES, DO NOT OVERSUBSCRIBE THEM.
                // Each candidate is itself OpenMP-parallel, so K concurrent
                // children each taking every core means K-fold oversubscription
                // and most of the win is lost to context switching. Measured on
                // SARS-CoV-2, 12 cores, 4 candidates:
                //     sequential                        39.54 s   1.00x
                //     K=4, 12 threads each (oversub)    25.80 s   1.53x
                //     K=4,  3 threads each (nproc/K)    11.43 s   3.46x
                // All three produce a BYTE-IDENTICAL archive. Splitting the
                // cores is what converts the concurrency into actual speed.
                // THE MODEL THIS IMPLEMENTS, so the win is predictable rather
                // than a lucky measurement. Per-candidate work is S (serial) +
                // Q/t (parallel on t threads):
                //     serial:      N * (S + Q/P)
                //     concurrent:  S + N*Q/P        (K=N, P/N threads each)
                // The total PARALLEL work is identical in both; concurrency
                // pays the SERIAL part once instead of N times. Hence
                //     speedup = (N*S + a) / (S + a),  a = N*Q/P
                // which tends to N when the work is serial-heavy and to 1 when
                // it is perfectly parallel. Measured 3.50x at N=4 implies
                // S ~= 5a, i.e. this encoder's per-candidate work is ~83%
                // serial -- consistent with the job pool's own diagnostic
                // (speedup 2.68x against a 1.59 s Amdahl floor).
                //
                // So the win generalises wherever the serial fraction stays
                // high, and degrades gracefully (never inverts) if a future
                // stage becomes more parallel. Where RAM forces K < N the
                // speedup falls out of the same formula with K substituted --
                // on the largest inputs expect ~1.9-2.6x rather than 3.5x, and
                // that is a memory bound, not a modelling error.
                // SPLIT THE MEASURED MEMORY BUDGET ACROSS THE TWO LEVELS, so
                // peak RAM is what it was. K is the number of candidate-sized
                // resident images the box was measured to afford; with two
                // levels the concurrent images are (groups x members), so the
                // product is held at K rather than each level taking K.
                const size_t Kg = std::min(gmin.size(), K ? K : (size_t)1);
                g_l2_conc = std::max<size_t>(1, (K ? K : 1) / (Kg ? Kg : 1));
                // Level-1 children each get a real share of the machine for
                // round 2 -- that is the whole point. Grandchildren then get
                // the per-candidate share they had before.
                const unsigned per_child_threads =
                    (unsigned)std::max<size_t>(1, (size_t)omp_get_max_threads() / (Kg ? Kg : 1));
                g_l2_threads =
                    (unsigned)std::max<size_t>(1, (size_t)omp_get_max_threads() /
                                                  ((Kg?Kg:1) * (g_l2_conc?g_l2_conc:1)));
                fprintf(stderr,"  [a3] two-level fork: %zu MINOV group(s) x members "
                               "(round 2 shared) conc=%zu x %zu threads=%u/%u\n",
                        gmin.size(), Kg, g_l2_conc, per_child_threads, g_l2_threads);
                std::vector<pid_t> running;
                for(size_t gi=0; gi<gmin.size() && !child; ++gi){
                    pid_t pid = fork();
                    if(pid < 0){ perror("fork"); return 1; }
                    if(pid == 0){
                        child=true; mine=gi;
                        // Child's share of the machine. Set before any parallel
                        // region in the child runs.
                        omp_set_num_threads((int)per_child_threads);
                        g_thr_budget = per_child_threads;   // bounds raw threads too
                        { char b[16]; snprintf(b,sizeof b,"%u",per_child_threads);
                          setenv("OMP_NUM_THREADS", b, 1); }
                        break;
                    }
                    running.push_back(pid);
                    // Keep at most Kg group children alive at once.
                    while(running.size() >= Kg){
                        int st=0; pid_t done = wait(&st);
                        if(done > 0){
                            running.erase(std::remove(running.begin(),running.end(),done),
                                          running.end());
                            if(!WIFEXITED(st) || WEXITSTATUS(st)!=0)
                                fprintf(stderr,"  [a3] a candidate failed (status %d)\n",st);
                        } else break;
                    }
                }
                if(!child) while(!running.empty()){
                    int st=0; pid_t done = wait(&st);
                    if(done <= 0) break;
                    running.erase(std::remove(running.begin(),running.end(),done),running.end());
                    if(!WIFEXITED(st) || WEXITSTATUS(st)!=0)
                        fprintf(stderr,"  [a3] a candidate failed (status %d)\n",st);
                }
            }
            if(!child && !gs_child){
                // parent: every candidate is finished; keep the smallest.
                size_t bestsz=SIZE_MAX, bi=0;
                for(size_t ci=0; ci<cands.size(); ++ci){
                    std::string f = base + ".cand" + std::to_string(ci);
                    struct stat sb;
                    if(stat(f.c_str(),&sb)==0 && (size_t)sb.st_size < bestsz){
                        bestsz=(size_t)sb.st_size; bi=ci; }
                }
                if(bestsz==SIZE_MAX){ fprintf(stderr,"[a3] all candidates failed\n"); return 1; }
                // Report EVERY candidate's size, not just the winner. The
                // smaller grids are subsets of this one (GRID4 is candidates
                // 2,3,6,7; the compiled-in single point is candidate 2), so one
                // 8-point run answers "what would a smaller grid have cost?"
                // exactly -- no extra encodes.
                for(size_t ci=0; ci<cands.size(); ++ci){
                    std::string f = base + ".cand" + std::to_string(ci);
                    struct stat sb2; size_t sz2 = (stat(f.c_str(),&sb2)==0)?(size_t)sb2.st_size:0;
                    fprintf(stderr,"  [a3] cand%zu MAXMAP=%u MINOV=%u size=%zu%s\n",
                            ci, cands[ci].first, cands[ci].second, sz2, ci==bi?"  <-- WINNER":"");
                }
                // Promote the winner's dumps into this directory, then drop
                // every candidate's scratch dir, so what remains on disk is the
                // archive AND the streams that actually produced it.
                char pcwd[4096]; if(!getcwd(pcwd,sizeof pcwd)) pcwd[0]='\0';
                const std::string pabs =
                    (!base.empty() && base[0]=='/') ? base
                                                    : (std::string(pcwd) + "/" + base);
                for(size_t ci=0; ci<cands.size(); ++ci){
                    std::string f = pabs + ".cand" + std::to_string(ci);
                    if(ci==bi) rename(f.c_str(), pabs.c_str()); else remove(f.c_str());
                    std::string d = f + ".d";
                    if(DIR* dp = opendir(d.c_str())){
                        while(struct dirent* de = readdir(dp)){
                            if(de->d_name[0]=='.') continue;
                            std::string src = d + "/" + de->d_name;
                            if(ci==bi) rename(src.c_str(),
                                              (std::string(pcwd) + "/" + de->d_name).c_str());
                            else remove(src.c_str());
                        }
                        closedir(dp); rmdir(d.c_str());
                    }
                }
                // Promote the WINNER's caller outputs to the paths the caller
                // actually asked for, and drop the losers'. Without this the
                // surviving VCF is whichever candidate finished last (a
                // scheduling race, not reproducible) -- see the matching
                // comment in the child block below.
                if(CAPS_CALL){
                    auto pabso=[&](const char* p)->std::string{
                        std::string s(p);
                        return (!s.empty() && s[0]=='/') ? s : (std::string(pcwd) + "/" + s);
                    };
                    const char* cve = getenv("CALL_VCF");
                    const std::string cv = cve ? pabso(cve) : (std::string(pcwd) + "/out.vcf");
                    for(size_t ci=0; ci<cands.size(); ++ci){
                        std::string f = cv + ".cand" + std::to_string(ci);
                        if(ci==bi) rename(f.c_str(), cv.c_str()); else remove(f.c_str());
                    }
                    if(const char* dce = getenv("CAPS_DUMP_CONTIGS")){
                        const std::string dc = pabso(dce);
                        for(size_t ci=0; ci<cands.size(); ++ci){
                            std::string f = dc + ".cand" + std::to_string(ci);
                            if(ci==bi) rename(f.c_str(), dc.c_str()); else remove(f.c_str());
                        }
                    }
                    fprintf(stderr,"  [a3] promoted candidate %zu's calls -> %s\n", bi, cv.c_str());
                }
                fprintf(stderr,"  [a3] chose MAXMAP=%u MINOV=%u -> %zu B (one shared prefix)\n",
                        cands[bi].first, cands[bi].second, bestsz);
                fprintf(stderr,"ARCHIVE_TOTAL=%zu\n",bestsz);
                printf("ARCHIVE_TOTAL=%zu\n",bestsz);
                return 0;
            }
            // child: adopt this candidate and fall through to the suffix.
            if(!gs_child){
            //
            // Each child must also get its OWN directory for the VERIFY_DUMP
            // stream files. They are written to the process CWD, so with every
            // child sharing one directory the surviving dumps belonged to
            // whichever child ran LAST -- not to the winner. Measured: 12/12
            // dumps matched the losing candidate and only 5/12 the winner (the
            // 5 being the shared-prefix streams). The archive was always the
            // winner's, so sizes were never wrong, but decoding those dumps
            // reconstructed a different encoding than the one shipped.
            char cwdbuf[4096];
            if(!getcwd(cwdbuf,sizeof cwdbuf)) cwdbuf[0]='\0';
            const std::string absbase =
                (!base.empty() && base[0]=='/') ? base
                                                : (std::string(cwdbuf) + "/" + base);
            // LEVEL 1 SETS ONLY WHAT ROUND 2 NEEDS. MAXMAP, ARCHIVE, the
            // scratch dir and the caller's per-candidate output paths are all
            // per-CANDIDATE, so they are deferred to the level-2 fork below --
            // which is also why this child must not chdir yet. Verified there
            // is no file written between here and `emit chains`.
            g_l2_base    = absbase;
            g_l2_cwd     = cwdbuf;
            g_l2_members = gidx[mine];
            g_l2_maxmaps.clear();
            for(size_t ci : gidx[mine]) g_l2_maxmaps.push_back(cands[ci].first);
            g_l2_active  = true;
            setenv("MINOV",  std::to_string(gmin[mine]).c_str(), 1);
            // THE MAPPING SEARCH RUNS ONCE, AT THE GROUP'S CEILING. See the
            // level-2 block for why one search at max(MAXMAP) answers every
            // member of the group.
            { uint32_t cap=0; for(uint32_t m : g_l2_maxmaps) if(m>cap) cap=m;
              setenv("MAXMAP", std::to_string(cap).c_str(), 1); }
            sweep_minov = gmin[mine];
            }
        }
    }

    // ── round 2: build over the well-overlapped reads only ──────────────────
    sweep_minov=MINOV;                // builder: caller's floor
    if(getenv("LOWCOV") && sweep_minov>8) sweep_minov=8;   // see LOWCOV note at SW
    if(getenv("MINOV")) sweep_minov=(uint32_t)atoi(getenv("MINOV"));  // override for sweeps
    sweep();
    fprintf(stderr,"round2: probes=%zu links=%zu\n",probes,links);
    fprintf(stderr,"[FUNNEL] probes=%zu  seed_hit=%zu (%.1f%%)  cand_examined=%zu  links=%zu (%.2f%% of probes)\n",g_probe,g_sd_hit,100.0*g_sd_hit/(g_probe?g_probe:1),g_cand,links,100.0*links/(g_probe?g_probe:1));
    lap("round 2 (assembly)");
    if(CAPS_CALL){
        fprintf(stderr,"[HARVEST] commits=%zu branches(ccnt==2)=%zu  ccnt hist:",
                g_commits,g_branch.size());
        for(int z=0;z<9;++z) if(g_cc_hist[z]) fprintf(stderr," %d:%zu",z,g_cc_hist[z]);
        fprintf(stderr,"\n[HARVEST] tie L as %% of Lmax:");
        for(int z=0;z<11;++z) if(g_L_bucket[z]) fprintf(stderr," %d0%%:%zu",z,g_L_bucket[z]);
        fprintf(stderr,"\n[HARVEST] memory=%zu KB (%zu B/edge)\n",
                g_branch.size()*sizeof(GEdge)/1024,sizeof(GEdge));
        fprintf(stderr,"[NEARMISS] minL=%u  1-mismatch=%zu  2-mismatch=%zu\n",
                nm_minL,nm1,nm2);
    }

    // ── emit chains; singletons held back for mapping ────────────────────────
    std::string pg; pg.reserve((size_t)n*40);
    ppos.assign(n,UINT64_MAX);
    prc.assign(n,0);
    std::vector<uint32_t> leftovers;
    uint32_t multi=0;
    const bool BOTHSIDE = getenv("BOTHSIDE") != nullptr;
    for(uint32_t i=0;i<n;++i){
        if(!admit[i]){ leftovers.push_back(i); continue; }   // excluded in round 1
        if(contained[i]) continue;                 // stage 42: its container carries it
        if(prv[i]!=NONE) continue;
        if(nxt[i]==NONE){ leftovers.push_back(i); continue; }
        ++multi;
        // PgRC2 admits a read to the pseudogenome only if it is overlapped on
        // BOTH sides (getBothSidesOverlappedReads, AbstractOverlapPseudoGenome
        // Generator.cpp:64-85): a read with both a predecessor and a successor
        // stays; everything else moves to the LQ set and is mapped instead. We
        // emit whole chains, endpoints included, which is why our pseudogenome
        // is 2.65x theirs on S. acidocaldarius -- 6,757,250 bases against
        // 2,546,054 -- and every read position costs log2(span) more.
        //
        // The rule is uniform; its EFFECT is data-dependent. On long chains it
        // drops almost nothing; on short chains it drops most reads. A chain of
        // length 2 contributes nothing, both reads being endpoints.
        if(BOTHSIDE){
            uint32_t len=1; for(uint32_t c=i; nxt[c]!=NONE; c=nxt[c]) ++len;
            if(len<3){
                for(uint32_t c=i;;){ leftovers.push_back(c);
                                     if(nxt[c]==NONE) break; c=nxt[c]; }
                continue;
            }
            leftovers.push_back(i);                       // head: one side only
            const uint64_t _cspan0=pg.size();
            uint32_t cur=nxt[i];
            ppos[cur]=pg.size(); rappend(pg,cur,0);
            while(nxt[cur]!=NONE && nxt[nxt[cur]]!=NONE){
                uint32_t o=ovl[cur]; cur=nxt[cur];
                ppos[cur]=pg.size()-o; rappend(pg,cur,o); }
            if(CAPS_CALL) g_contig_spans.push_back({_cspan0,pg.size()});
            if(nxt[cur]!=NONE) leftovers.push_back(nxt[cur]);   // tail
            continue;
        }
        const uint64_t _cspan0=pg.size();
        uint32_t cur=i; ppos[cur]=pg.size(); rappend(pg,cur,0);
        while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                               ppos[cur]=pg.size()-o; rappend(pg,cur,o); }
        if(CAPS_CALL) g_contig_spans.push_back({_cspan0,pg.size()});
    }
    if(getenv("DBG_OVL")){
        // Overlap-length histogram over committed links, as a FRACTION of read
        // length -- the quantity that decides how many new bases each link adds.
        size_t nb[11]={0}; size_t tot=0; double sumfrac=0;
        for(uint32_t i=0;i<n;++i){
            if(nxt[i]==NONE) continue;
            const double f=(double)ovl[i]/(double)rlen[i];
            int b=(int)(f*10); if(b>10)b=10; if(b<0)b=0;
            ++nb[b]; ++tot; sumfrac+=f;
        }
        fprintf(stderr,"[ovl] links=%zu mean_overlap_frac=%.3f\n",tot,tot?sumfrac/tot:0.0);
        for(int b=0;b<10;++b)
            fprintf(stderr,"[ovl]   %2d-%3d%% of readlen : %8zu links (%5.2f%%)  each adds ~%.0f bases\n",
                    b*10,(b+1)*10,nb[b],tot?100.0*nb[b]/tot:0.0,(1.0-(b*10+5)/100.0)*Lmax);
    }
    fprintf(stderr,"links=%zu chains(multi)=%u leftovers=%zu pg after chains=%zu\n",
            links,multi,leftovers.size(),pg.size());
    lap("emit chains");

    // ── NEAR-MISS: translate to pg coordinates and cluster ──────────────────
    // Each observation is (tail read a, offset of the differing base within a).
    // ppos[a] is a's start in the pseudogenome, so ppos[a]+offset is the pg
    // coordinate of the candidate variant. Reads covering the SAME genomic
    // site land on the SAME pg coordinate, which is what turns ~44k raw
    // observations into a much smaller set of distinct sites.
    //
    // The depth filter is the whole point: a heterozygous site is seen by many
    // independent read pairs and therefore accumulates a high count at one
    // coordinate; a sequencing error is seen once and sits alone. This is the
    // same logic the pileup applies -- the difference is that these candidate
    // pairs cost nothing to obtain, being a byproduct of chaining.
    if(CAPS_CALL){
        // key = pg position << 2 | alt base, so the same position with two
        // different alternate bases stays two candidates (a real multi-allelic
        // site) instead of being merged into one.
        std::vector<uint64_t> sites; sites.reserve(nm_obs.size());
        size_t unplaced=0;
        for(const auto& o : nm_obs){
            if(ppos[o.read_a]==UINT64_MAX){ ++unplaced; continue; }
            sites.push_back(((ppos[o.read_a]+o.pos_in_a)<<2)|(uint64_t)o.alt);
        }
        std::sort(sites.begin(),sites.end());
        size_t distinct=0; size_t d2=0,d3=0,d5=0,d10=0;
        const char* NMV=getenv("CAPS_NM_VCF");
        const int NMD=getenv("CAPS_NM_DEPTH")?atoi(getenv("CAPS_NM_DEPTH")):5;
        FILE* nv=NMV?fopen(NMV,"w"):nullptr;
        if(nv) fprintf(nv,"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for(size_t i=0;i<sites.size();){
            size_t j=i; while(j<sites.size()&&sites[j]==sites[i]) ++j;
            const size_t depth=j-i; ++distinct;
            if(depth>=2)  ++d2;
            if(depth>=3)  ++d3;
            if(depth>=5)  ++d5;
            if(depth>=10) ++d10;
            if(nv && depth>=(size_t)NMD){
                const uint64_t p=sites[i]>>2; const uint8_t alt=(uint8_t)(sites[i]&3ULL);
                // locate the contig containing pg position p
                auto it=std::upper_bound(g_contig_spans.begin(),g_contig_spans.end(),p,
                        [](uint64_t v,const std::pair<uint64_t,uint64_t>& s){ return v<s.first; });
                if(it!=g_contig_spans.begin()){
                    --it;
                    if(p<it->second && p<pg.size()){
                        const size_t cid=(size_t)(it-g_contig_spans.begin());
                        fprintf(nv,"contig_%zu\t%llu\t.\t%c\t%c\t.\tPASS\tSVTYPE=SNV;SRC=nearmiss;DP=%zu\n",
                                cid,(unsigned long long)(p-it->first+1),pg[p],"ACGT"[alt],depth);
                    }
                }
            }
            i=j;
        }
        if(nv){ fclose(nv); fprintf(stderr,"[NEARMISS] VCF (depth>=%d) -> %s\n",NMD,NMV);
            // Emit the MATCHING contig FASTA. CAPS_DUMP_CONTIGS writes the
            // caller's COLLAPSED substrate (different count, different ids),
            // so lifting this VCF against that file silently matches the wrong
            // contigs -- these records are indexed by g_contig_spans.
            std::string fa(NMV); fa+=".contigs.fa";
            if(FILE* cf=fopen(fa.c_str(),"w")){
                for(size_t k=0;k<g_contig_spans.size();++k)
                    fprintf(cf,">contig_%zu\n%s\n",k,
                            pg.substr(g_contig_spans[k].first,
                                      g_contig_spans[k].second-g_contig_spans[k].first).c_str());
                fclose(cf);
                fprintf(stderr,"[NEARMISS] matching contigs -> %s\n",fa.c_str());
            }
        }
        fprintf(stderr,"[NEARMISS] obs=%zu placed=%zu unplaced=%zu -> distinct (pos,alt)=%zu\n",
                nm_obs.size(),sites.size(),unplaced,distinct);
        fprintf(stderr,"[NEARMISS] sites by depth: >=2:%zu  >=3:%zu  >=5:%zu  >=10:%zu\n",
                d2,d3,d5,d10);
    }

    // The prefix index is not consulted anywhere in the mapping stage -- that
    // stage builds and queries its own seed index -- but it was staying live
    // across it, so the two large indexes coexisted at exactly the moment peak
    // RSS is set. It is next needed by the survivor sweep, which runs over ~13k
    // reads rather than 851k, so it is cheaper to drop it here and rebuild a
    // small one there than to carry the full one through.
    phase("greedy-sweep");
    // STAGE 41: return the freed arena to the OS.
    //
    // Three targeted reductions in a row measured ZERO -- ppos uint64->uint32
    // (-3.4 MB of array), the packed mapping index (13 B -> 8 B/entry and a
    // 30 MB temp removed), and the parallel loader. Each shrank something real
    // and none moved peak RSS, which is the signal: the peak is not the live
    // set. Live structures at mapping total ~155 MB against a measured 232 MB.
    //
    // The gap is glibc holding freed arenas. ARCS hit exactly this and says so
    // in src/vodbg_pg.cpp: "malloc_trim is what makes the release visible: glibc
    // keeps freed arenas by default, so an earlier measurement recovered only
    // 52 MB of a 312 MB free until the arena was returned explicitly." The
    // sweep's per-level candidate array and the 851k-entry prefix map are both
    // freed before mapping and both stay resident without this.
    { std::vector<uint64_t>().swap(pent); std::vector<uint32_t>().swap(ptab); }
#if defined(__GLIBC__)
    malloc_trim(0);
#endif

    // ── Stage C: pigeonhole mapping into the pg, forward then RC ─────────────
    std::vector<uint8_t> matched(n,0);
    hugefill(readMM,n,(uint8_t)255);
    hugehint(readMM.data(), readMM.size());
    size_t n_matched=0;
    {
        // Part width must SCALE with the tolerance, not stay pinned at 32.
        // Pigeonhole needs MAXMM+1 parts to fit inside one read, so a fixed
        // 32-base part caps tolerance at 150/32-1 = 3 mismatches -- ask for more
        // and the parts no longer fit, no seed is indexed, and mapping silently
        // does nothing. PgRC2 makes the part length the parameter and derives
        // the tolerance from it (targetMismatches = readLength/partLength - 1).
        // Same here: parts are as wide as they can be, capped at 32 because that
        // is what a uint64 holds exactly at 2 bits per base.
        // STAGE 17. The above ties acceptance to the seed geometry: NPARTS =
        // MAXMM+1 parts means a read is rejected past 3 mismatches, because that
        // is all the pigeonhole lemma can certify with parts this wide.
        //
        // PgRC2 does not tie them together. ReadsMatchers.cpp:700 reads
        //     uint8_t maxMismatches = readLength / minCharsPerMismatch;
        // with minCharsPerMismatch = 3 at CODER_LEVEL_NORMAL, so a 150-base read
        // is accepted at up to FIFTY mismatches. The seed (their
        // readsExactMatchingChars = 38) only proposes candidates; acceptance is a
        // separate, far looser test. Their log shows what that buys:
        //     290,521 LQ reads -> 11,298 survivors   (96.1% placed)
        // against ours at 74.8%. The survivor pseudogenome that holds 94% of the
        // remaining gap is made almost entirely of reads this limit rejected.
        //
        // Their seed is 38 bases; ours stays 32 because that is what a uint64
        // holds exactly at 2 bits per base. A shorter seed is more sensitive, not
        // less, so this is not a handicap -- it proposes more candidates, and
        // acceptance is what decides.
        //
        // Not copied: they re-scan and keep the BEST match per read (the
        // "better-matches" column in their log, ReadsMatchers.cpp:315-328). That
        // only shrinks their mismatch stream, which this progression does not
        // measure, so first-acceptable wins here. Noted rather than silently
        // skipped.
        // STAGE 20. Seeds were indexed at fixed multiples of SEEDW (0,32,64,96),
        // which has two defects at a 150-base read: a read whose only clean
        // 32-mer starts at, say, offset 17 is invisible, and bases 128-149 are
        // in no seed at all. PgRC2 does not do this -- at CODER_LEVEL_NORMAL the
        // matcher is mode 'c', CopMEMReadsApproxMatcher, which finds maximal
        // exact matches at ANY offset (mapReadsIntoPg, ReadsMatchers.cpp:715-735).
        // Their fixed 38 is a minimum MEM length, not a position.
        //
        // Note their pigeonhole guarantee is only
        //     targetMismatches = readLength/exactMatchingChars - 1 = 150/38-1 = 2
        // while maxMismatches is 50, so past 2 mismatches their seed is a
        // heuristic filter too. The difference is purely that their filter looks
        // everywhere and ours looked in four places.
        //
        // Sliding by SEEDSTRIDE recovers most of that for a proportional index
        // cost and no extra scan cost -- the pg is still swept once, one hash
        // probe per position, only the candidate lists get longer.
        // MAXMAP: mismatch ceiling for accepting a read placement.
        // Lmax/3 copies ReadsMatchers.cpp:700, but that ceiling was chosen when
        // the mismatch POSITION and COUNT streams were not being counted in our
        // total (see the Stage 27 note below -- it predicted exactly this).
        // With those streams counted, the real economics are: a mismatch costs
        // ~0.88 B (position+count+symbol, measured) while a base stored as pg
        // literal costs ~0.229 B (measured), so for a length-L read the
        // break-even is m ~= 0.26*L -- BELOW Lmax/3. Accepting past that point
        // buys a shorter pg with a mismatch stream that costs more than the
        // bases it saved. Env-overridable so the real optimum can be swept
        // against the CORRECTED total rather than assumed.
        // Re-swept on REAL full files after the mm_pos/mm_cnt transforms made
        // mismatches cheaper, which shifts this optimum upward (cheaper
        // mismatches -> accept more of them -> shorter pg). Real optima now:
        // E. coli 20, P. aeruginosa 27, but validated across all five real
        // files a fixed 20 LOSES 1.2% on P. falciparum, making the aggregate a
        // wash (+0.04% worse than a fixed 12). Minimax regret picks 12: worst
        // case 0.54% versus 20's 1.24%. So 12 stays, now justified on real
        // full files rather than the subsamples it was first derived from.
        // The property was identified: READ LENGTH. Every dataset the constant
        // was tuned on is 100-151 bp, a range over which a constant and a ratio
        // are indistinguishable. PgRC2 has always used the ratio form --
        // ReadsMatchers.cpp:700, `maxMismatches = readLength / minCharsPerMismatch`
        // -- i.e. one mismatch permitted per N bases, not a fixed count.
        //
        // A fixed ceiling is wrong in both directions: at 251 bp it is 4.8%
        // tolerance where 151 bp gets 8%, so long reads are under-mapped and
        // spill into the pg as literal; at short read lengths it would be too
        // permissive. Expressed as a ratio it is scale-free.
        //
        // The divisor comes from our own measured economics, not from theirs:
        // a mismatch costs ~0.88 B (position+count+symbol) against ~0.229 B for
        // a pg-literal base, so the break-even is ~0.26*L. The tuned optimum of
        // 12 at L=151 is L/12.6; L/13 reproduces it (11.6 at 151, 7.7 at 100,
        // 19.3 at 251) while scaling correctly, and stays well under the
        // 0.26*L cost ceiling. MIN_MAXMAP floors it so very short reads still
        // get a usable tolerance.
        const uint32_t MAXMAP_DIV = (uint32_t)(getenv("MAXMAP_DIV")?atoi(getenv("MAXMAP_DIV")):13);
        const uint32_t MIN_MAXMAP = 6;
        uint32_t MAXMAP = Lmax / (MAXMAP_DIV ? MAXMAP_DIV : 13);
        if(MAXMAP < MIN_MAXMAP) MAXMAP = MIN_MAXMAP;
        // Coverage-derived: L/13 was minimax-swept on the 7 locked datasets
        // (all leftover_frac 0.222-0.518) and stays exactly there for any of
        // them -- this only fires above that measured range, same T0/T1
        // margin as MEM_MAXMM's gate. At low coverage the ALTERNATIVE to
        // accepting a mismatch-heavy placement isn't "slightly more literal",
        // it's falling through to the far more expensive region-scale MEM
        // path entirely, so the break-even point shifts. Measured directly on
        // Drosophila SRR40104920 (leftover_frac=0.812): sweeping MAXMAP with
        // everything else fixed gives a clean interior optimum at 30 (150/5),
        // not at our break-even ESTIMATE of 39 -- 36,982,418 -> 36,836,970 B
        // (-145,448 B, -0.393%), 250,841 reads placed against 158,351 before.
        // Worse on both sides (24: 36,856,348; 39: 36,873,742; 50: 37,050,398,
        // above the untouched baseline) -- a real optimum, not "more is
        // better". Calibrated on this one dataset; the ramp's SHAPE (linear
        // in leftover_frac) is the honest part, its endpoint is not yet
        // cross-validated on a second low-coverage file.
        {
            const double lf = n ? (double)leftovers.size()/n : 0.0;
            const double T0=0.60, T1=0.85;
            if(lf>T0){
                const double t = std::min(1.0, (lf-T0)/(T1-T0));
                const uint32_t divLow = 5;                     // -> Lmax/5 = 30 at Lmax=150
                const uint32_t divNow = MAXMAP_DIV - (uint32_t)std::lround((double)(MAXMAP_DIV-divLow)*t);
                uint32_t adj = Lmax / (divNow?divNow:1);
                if(adj < MIN_MAXMAP) adj = MIN_MAXMAP;
                if(adj > MAXMAP) MAXMAP = adj;
            }
        }
        if(getenv("MAXMAP")) MAXMAP = (uint32_t)atoi(getenv("MAXMAP"));   // override for sweeps
        fprintf(stderr,"  [maxmap] Lmax=%u -> MAXMAP=%u (L/%u)\n", Lmax, MAXMAP, MAXMAP_DIV);
        // Sensitivity is SEEDW + SEEDSTRIDE - 1: the shortest exact read/pg
        // stretch guaranteed to be found. copMEM's K is 28 with k1*k2 = 10, so
        // theirs is 28+10-1 = 37 against our 32+8-1 = 39. A SHORTER seed with a
        // WIDER stride therefore improves sensitivity and shrinks the index at
        // the same time -- unlike shortening the stride, which improves
        // sensitivity and inflates the index (measured: stride 4 mapped 262 more
        // reads for 150 MB).
        uint32_t SEEDW=(argc>6)?(uint32_t)atoi(argv[6]):32;   // <= 32 (uint64 at 2 bits/base)
        if(SEEDW>32) SEEDW=32; if(SEEDW<8) SEEDW=8;
        const uint32_t SEEDSTRIDE=(argc>5)?(uint32_t)atoi(argv[5]):8;
        const uint32_t K2=(argc>8)?(uint32_t)atoi(argv[8]):1;   // query-side stride
        const uint8_t MMGOOD=(uint8_t)(getenv("MMGOOD")?atoi(getenv("MMGOOD")):0);
        const uint32_t NPARTS=(Lmax>=SEEDW)?((Lmax-SEEDW)/SEEDSTRIDE+1):1;
        if(SEEDW<8) SEEDW=8;                     // below this a seed is noise
        const uint64_t SW_MASK=(SEEDW>=32)?~0ULL:((1ULL<<(2*SEEDW))-1);
        fprintf(stderr,"  sensitivity floor: %u bases (SEEDW %u + k1 %u * k2 %u - 1)\n",SEEDW+SEEDSTRIDE*K2-1,SEEDW,SEEDSTRIDE,K2);
        auto packW=[&](const char* p,uint64_t& out)->bool{
            uint64_t k=0;
            for(uint32_t i=0;i<SEEDW;++i){ int v=b2(p[i]); if(v<0) return false;
                                          k=(k<<2)|(uint64_t)v; }
            out=k; return true; };
        fprintf(stderr,"  mapping: %u seeds x %u bases, stride %u, accepting up to %u mismatches\n",
                NPARTS,SEEDW,SEEDSTRIDE,MAXMAP);
        // Flat CSR, not unordered_map<uint64,vector<...>>. At stride 8 this index
        // holds ~4M entries over ~3.5M distinct keys, and a map would put a node
        // header plus a separately-allocated vector behind every one of them --
        // measured, that alone took peak RSS from 389 MB to 673 MB. Sorting
        // (key,rid,part) once and pointing an open-addressing table at each key's
        // first occurrence answers the same queries from three flat arrays.
        // STAGE 40: one packed array instead of three plus a temp.
        //
        // The modelled orientation flip (index the pg, stream reads) is
        // unaffordable: at 2.5x fewer operations it needs 141 MB of index and
        // peak RSS is already 233 MB against PgRC2's 232 MB. At EQUAL index size
        // the flip is 80% worse. So their mapping advantage is not orientation
        // but a denser index -- copMEM at their parameters samples the pg every
        // 5, giving ~4.3M entries against our 1.85M inside the same footprint.
        // The way to afford density is a smaller entry, which is this.
        //
        // Was: mkey(8) + mrid(4) + mpart(1) = 13 B/entry, built through a
        // vector<E> of 16 B/entry that coexists with all three during the copy --
        // 29 B/entry at the build moment, ~54 MB here.
        //
        // Now: (key32 << 32) | (rid << 3) | part in ONE uint64, sorted in place.
        // Sorting the whole word sorts by key, because the key occupies the high
        // bits -- so the run-detection and first-occurrence logic are unchanged.
        // 8 B/entry, no temp, no copy: ~15 MB.
        //
        // The key is truncated to 32 bits, which is exact for SEEDW <= 16 (the
        // default) and a partial key above it. Partial is still correct: a
        // colliding key only produces a candidate, and every candidate is
        // verified against the text anyway. It costs a few wasted verifies, never
        // a wrong placement.
        std::vector<uint64_t> ment;
        {
            // Exactly the upper bound, so push_back never doubles and no
            // shrink_to_fit copy is needed. Reserving half instead made the
            // vector reallocate mid-fill and then copy again on shrink, which
            // cost more peak than the packing saved: mapping RSS went 232 -> 240
            // MB on the first attempt.
            ment.reserve((size_t)leftovers.size()*NPARTS);
            hugehint((void*)ment.data(), (size_t)leftovers.size()*NPARTS*sizeof(uint64_t));
            for(uint32_t rid:leftovers){
                if(rlen[rid]<SEEDW) continue;
                for(uint32_t p=0;p<NPARTS;++p){
                    const uint32_t off=p*SEEDSTRIDE;
                    if(off+SEEDW>rlen[rid]) break;
                    const uint64_t k32=(uint64_t)(uint32_t)rseed(rid,off,SEEDW);
                    ment.push_back((k32<<32)|((uint64_t)rid<<3)|(uint64_t)p);
                }
            }
            std::sort(ment.begin(),ment.end());
        }
        auto MKEY =[&](size_t i)->uint32_t{ return (uint32_t)(ment[i]>>32); };
        auto MRID =[&](size_t i)->uint32_t{ return (uint32_t)((ment[i]>>3)&0x1FFFFFFFULL); };
        auto MPART=[&](size_t i)->uint32_t{ return (uint32_t)(ment[i]&7ULL); };
        size_t msize=1; while(msize < ment.size()*2+1) msize<<=1;
        const uint64_t MMASK=msize-1;
        std::vector<uint32_t> mtab; hugefill(mtab,msize,UINT32_MAX);
        // ment and mtab are the pigeonhole index -- the two largest random-access
        // structures in the stage perf attributes 32.4% of cycles to.
        hugehint(ment.data(), ment.size()*8);
        hugehint(mtab.data(), mtab.size()*4);
        fprintf(stderr,"[MEMRPT] ment %zu entries = %zu MB   mtab %zu slots = %zu MB\n",
                ment.size(), ment.size()*8/1000000, msize, msize*4/1000000);
        auto mmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
        for(size_t i=0;i<ment.size();++i){
            if(i && MKEY(i)==MKEY(i-1)) continue;
            size_t h=mmix(MKEY(i))&MMASK;
            while(mtab[h]!=UINT32_MAX) h=(h+1)&MMASK;
            mtab[h]=(uint32_t)i;
        }
        // STAGE 26. The mapping stage is probe-bound, not verify-bound and not
        // probe-COUNT-bound. 43M probes (21.6M pg positions x 2 directions) hit
        // an index of ~24 MB of CSR arrays plus a 16 MB open-addressing table --
        // past L3, so essentially every probe is a memory access. That is also
        // why stage 24's coprime sampling lost: it halved the probes but doubled
        // the table, and the table was the problem.
        //
        // Almost every probe MISSES: most pg positions carry no read seed. A one
        // bit per slot presence filter is 2 MB at this size, stays in cache, and
        // rejects those misses without ever touching the 40 MB structure. Only a
        // survivor pays for the real lookup. False positives cost one wasted
        // lookup and cannot change the result.
        std::vector<uint64_t> filt((size_t)1<<18,0);      // 2^24 bits = 2 MB
        const uint64_t FMASK=((uint64_t)1<<24)-1;
        for(size_t i=0;i<ment.size();++i){
            if(i && MKEY(i)==MKEY(i-1)) continue;
            const uint64_t h=mmix(MKEY(i))&FMASK;
            filt[h>>6]|=1ULL<<(h&63);
        }
        auto mmaybe=[&](uint64_t k)->bool{
            const uint64_t h=mmix(k)&FMASK;
            return (filt[h>>6]>>(h&63))&1ULL;
        };
        auto mfind=[&](uint64_t k)->uint32_t{
            size_t h=mmix(k)&MMASK;
            while(mtab[h]!=UINT32_MAX){ if(MKEY(mtab[h])==(uint32_t)k) return mtab[h]; h=(h+1)&MMASK; }
            return UINT32_MAX;
        };
        fprintf(stderr,"  seed index: %zu entries (packed, %.1f MB)\n",
                ment.size(),ment.size()*8.0/1048576.0);
        // The pg sweep is read-only against `seeds` and `reads`; the only shared
        // write is the matched flag. Threads take disjoint slices and collect
        // hits locally, so the merge is a union -- no atomics, no false sharing,
        // and the result cannot depend on thread count. Reading `matched` inside
        // the loop stays an optimisation only: a stale read costs a redundant
        // verify, never a wrong answer.
        bool scanning_rc=false;      // stage 46: which strand this scan is on
        auto scan=[&](const std::string& text){
            if(text.size()<SEEDW) return;
            // STAGE 25. The verify loop was the whole of the 2.6x mapping
            // deficit. Stage 21 packed the reads but then, for every candidate,
            // unpacked all 150 bases into a stack buffer and compared them
            // byte-by-byte -- ~300 operations per candidate, run ~26M times.
            //
            // PgRC2 does not unpack (SymbolsPackingFacility.cpp:344-362): it
            // packs the pattern a word at a time and compares whole packed
            // words, only falling back to per-symbol work when a word differs.
            // Packing the text ONCE per direction goes further -- both sides are
            // then packed, so a 150-base comparison is 5 word loads, 5 XORs and
            // 5 popcounts instead of 300 byte operations.
            std::vector<uint64_t> tpk((text.size()+31)/32+2,0);
            for(size_t i=0;i<text.size();i+=32){
                uint64_t v=0;
                for(uint32_t j=0;j<32;++j){
                    const size_t idx=i+j;
                    v=(v<<2)|(uint64_t)((idx<text.size())?std::max(0,b2(text[idx])):0);
                }
                tpk[i>>5]=v;
            }
            auto tw32=[&](uint64_t ab)->uint64_t{
                const size_t wi=ab>>5; const uint32_t sh=(uint32_t)(ab&31)*2;
                const uint64_t hi=(wi<tpk.size())?tpk[wi]:0ULL;
                if(!sh) return hi;
                const uint64_t lo=(wi+1<tpk.size())?tpk[wi+1]:0ULL;
                return (hi<<sh)|(lo>>(64-sh));
            };
            unsigned T=std::thread::hardware_concurrency(); if(!T) T=1;
            if(text.size()<(1u<<20)) T=1;
            std::vector<std::vector<std::pair<uint32_t,uint32_t>>> hit(T);
            std::vector<std::vector<uint32_t>> hmm(T);
            std::vector<std::thread> th;
            const size_t chunk=(text.size()+T-1)/T;
            // SAME T CHUNKS, FEWER WORKERS. T and the chunk partition are left
            // exactly as they were, so hit[t]/hmm[t] hold the same hits and the
            // serial merge below applies them in the same (thread, insertion)
            // order -- which matters, because that merge takes only strict
            // improvements, so a different enumeration order could pick a
            // different equal-mismatch placement. Only the number of chunks in
            // flight is bounded. That also drops the per-chunk `slot` array
            // (one uint32 per read) from T resident copies to W.
            const unsigned Wk = std::min<unsigned>(T, thr_budget());
            std::atomic<unsigned> nextc{0};
            th.reserve(Wk);
            for(unsigned w=0; w<Wk; ++w){
                th.emplace_back([&]{
                  for(unsigned t=nextc.fetch_add(1); t<T; t=nextc.fetch_add(1)){
                    const size_t lo=(size_t)t*chunk;
                    size_t hi=std::min(text.size(),lo+chunk);
                    if(lo>=hi) continue;
                    // A4. Keep only this thread's BEST placement per read instead
                    // of every accepted one. Measured on S. aureus: 59,187,817
                    // accepted hits for 1,705,714 reads -- 34.7 per read -- held
                    // as 1152 MB of hit/hmm capacity, which is the entire +877 MB
                    // that stage adds to peak RSS. The seed index is only 39.7 MB.
                    //
                    // Output is unchanged. The serial merge below applies hits in
                    // (thread, insertion) order and takes only strict improvements,
                    // so a read's final value is the minimum mismatch count and its
                    // final position is the FIRST hit achieving that minimum.
                    // Collapsing to the first strict minimum within each thread
                    // yields that identical pair, because a later equal-or-worse
                    // hit could never have displaced it in the merge either.
                    //
                    // Cost is one uint32 slot per read per thread, freed when the
                    // thread exits -- and it replaces a structure that grew with
                    // coverage depth rather than with the number of reads.
                    std::vector<uint32_t> slot((size_t)n, UINT32_MAX);
                    // start SEEDW-1 early so seeds spanning the boundary are not lost
                    const size_t from=(lo>=(size_t)(SEEDW-1))?lo-(SEEDW-1):0;
                    uint64_t k=0; uint32_t good=0;
                    for(size_t p=from;p<hi;++p){
                        int v=b2(text[p]);
                        if(v<0){ good=0; k=0; continue; }
                        k=((k<<2)|(uint64_t)v)&SW_MASK; ++good;
                        if(good<SEEDW) continue;
                        const size_t seedStart=p-SEEDW+1;
                        if(seedStart<lo) continue;      // owned by the previous slice
                        // STAGE 24: copMEM's lemma applies to the QUERY side too.
                        // Stage 20 took it for sensitivity (sampling the reads)
                        // but still probed every one of the pg's 21.6M positions,
                        // i.e. k2 = 1. Their k1*k2 <= L-K+1 with BOTH sides
                        // sampled is what makes their matching 0.41 s against our
                        // 1.06 s. Sampling the query every K2 and the reads every
                        // SEEDSTRIDE keeps the same floor as long as the product
                        // is unchanged, so K2=2 with stride halved is free
                        // sensitivity-wise and halves the probes.
                        if(K2>1 && (seedStart%K2)) continue;
                        if(!mmaybe(k)) continue;              // 2 MB filter, stays cached
                        const uint32_t ix=mfind(k); if(ix==UINT32_MAX) continue;
                        for(uint32_t q=ix;q<ment.size()&&MKEY(q)==(uint32_t)k;++q){
                            // ── SOFTWARE PREFETCH ───────────────────────────
                            // perf puts 32.4% of all cycles in this worker and
                            // the sweep runs at IPC 0.62 with 3.9 billion cache
                            // misses: the loop is memory-LATENCY bound, not
                            // bandwidth or compute bound. It walks 1.93 billion
                            // candidates on 3M human reads (59.5 per seed hit),
                            // and every one dereferences readMM[rid], rlen[rid]
                            // and woff[rid] at a random index.
                            //
                            // Those addresses are KNOWN several iterations
                            // ahead, because ment[] is walked sequentially --
                            // so the misses can be overlapped instead of
                            // serialised. woff is fetched one step further out
                            // than the rest because rpk[woff[rid]] is a
                            // dependent load and needs woff resident first.
                            //
                            // A prefetch has no semantic effect whatsoever: it
                            // cannot change which candidates are examined or
                            // accepted, so the archive is identical by
                            // construction, not by measurement.
                            if(q+8<ment.size() && MKEY(q+8)==(uint32_t)k){
                                const uint32_t r8=MRID(q+8);
                                __builtin_prefetch(&woff[r8],0,1);
                            }
                            if(q+3<ment.size() && MKEY(q+3)==(uint32_t)k){
                                const uint32_t r3=MRID(q+3);
                                __builtin_prefetch(&readMM[r3],0,1);
                                __builtin_prefetch(&rlen[r3],0,1);
                                __builtin_prefetch((const char*)&rpk[woff[r3]],0,1);
                            }
                            const uint32_t rid=MRID(q);
                            const size_t off=(size_t)MPART(q)*SEEDSTRIDE;
                            if(seedStart<off) continue;
                            const size_t st=seedStart-off;
                            const uint8_t cur=readMM[rid];
                            // Rescanning every already-placed read is what the
                            // best-match pass costs. A read already placed at or
                            // below MMGOOD is left alone: the remaining gain per
                            // read is small and the work per read is not.
                            if(cur<=MMGOOD) continue;
                            const uint32_t lim=(cur==255)?MAXMAP:(uint32_t)cur-1;
                            const uint32_t RL=rlen[rid];
                            if(st+RL>text.size()) continue;
                            // packed vs packed: 2 bits per base, so a XOR leaves a
                            // non-zero 2-bit field exactly where the bases differ.
                            // Folding the two bits together and masking to the low
                            // bit of each pair turns that into a popcount.
                            uint32_t mm=0;
                            const uint64_t rbase=woff[rid]*32ULL;
                            for(uint32_t j=0;j<RL;j+=32){
                                const uint32_t take=(RL-j<32)?(RL-j):32;
                                uint64_t d=w32(rbase+j)^tw32((uint64_t)st+j);
                                if(take<32) d&=~0ULL<<(2*(32-take));
                                d=(d|(d>>1))&0x5555555555555555ULL;
                                mm+=(uint32_t)__builtin_popcountll(d);
                                if(mm>lim) break;
                            }
                            if(mm<=lim){
                                uint32_t& sl = slot[rid];
                                if(sl==UINT32_MAX){
                                    sl=(uint32_t)hit[t].size();
                                    hit[t].push_back({rid,(uint32_t)st}); hmm[t].push_back(mm);
                                } else if(mm < hmm[t][sl]){      // strict: keeps the first minimum
                                    hit[t][sl].second=(uint32_t)st; hmm[t][sl]=mm;
                                }
                            }
                        }
                    }
                  }
                });
            }
            for(auto& x:th) x.join();
            {   // A4 diagnostic: is the +877 MB the candidate vectors, or the index?
                size_t hb=0,hn=0;
                for(size_t t=0;t<hit.size();++t){
                    hn+=hit[t].size();
                    hb+=hit[t].capacity()*sizeof(std::pair<uint32_t,uint32_t>)
                       +hmm[t].capacity()*sizeof(uint32_t);
                }
                fprintf(stderr,"  [a4] accepted hits=%zu  hit/hmm capacity=%.1f MB"
                               "  seedindex=%.1f MB  reads=%zu\n",
                        hn, hb/1048576.0, ment.size()*8.0/1048576.0, (size_t)n);
            }
            // Keep the fewest-mismatch placement per read. Threads may each
            // have found a different one; the merge is serial and authoritative,
            // so the readMM they raced on was only ever a bound hint.
            for(size_t t=0;t<hit.size();++t)
                for(size_t q=0;q<hit[t].size();++q){
                    const auto& pr=hit[t][q]; const uint32_t m=hmm[t][q];
                    if(m<readMM[pr.first]){
                        if(readMM[pr.first]==255) ++n_matched;
                        readMM[pr.first]=(uint8_t)m; matched[pr.first]=1;
                        ppos[pr.first]=pr.second; prc[pr.first]=scanning_rc?1:0;
                    }
                }
        };
        scanning_rc=false; scan(pg);
        rc_inplace(pg); scanning_rc=true; scan(pg); rc_inplace(pg);
        scanning_rc=false;
    }
    lap("pigeonhole mapping");

    // ── LEVEL-2 CANDIDATE FORK: split the MINOV group by MAXMAP ─────────────
    //
    // Everything above this line -- round 2, chain emission AND the pigeonhole
    // mapping search -- is now computed ONCE per MINOV group instead of once
    // per member.
    //
    // Mapping looks MAXMAP-dependent and is not. The search keeps the
    // MINIMUM-mismatch placement per read: `lim = (cur==255) ? MAXMAP : cur-1`
    // with a strict-improvement update, so MAXMAP is only the initial cap that
    // decides whether a read is placed at all, never which placement wins.
    // Running at the group's ceiling therefore yields every member's answer by
    // thresholding readMM.
    //
    // Proven from the production HG002 log before a line was changed. The
    // mismatch histograms of the four candidates in one group are bit-identical
    // in every bin below the smaller ceiling:
    //     MAXMAP=7   0:207172 1:928067 2:239064 3:117518 4:76881 5:57665
    //                6:45931 7:38044                        8+:0
    //     MAXMAP=11  ...same bins 0-7...   8:32688 9:28422 10:25114 11:22563
    //     MAXMAP=18  ...same bins 0-11...                    >=12:109346
    //     MAXMAP=29  ...same bins 0-11...                    >=12:195906
    // and the cumulative sums reproduce each candidate's mapped count exactly:
    //     1,710,342 = sum(bins 0..7)                    = placed(MAXMAP=7)
    //     1,819,129 = 1,710,342 + 108,787 (bins 8..11)  = placed(MAXMAP=11)
    //     1,928,475 = 1,819,129 + 109,346               = placed(MAXMAP=18)
    //     2,015,035 = 1,819,129 + 195,906               = placed(MAXMAP=29)
    // A larger ceiling only ADDS reads in higher bins; it never moves a read
    // that was already placeable at a lower count.
    if(g_l2_active){
        // TEAR THE OPENMP RUNTIME DOWN BEFORE FORKING. Round 2 ran parallel
        // regions in this process, so forking without this is exactly the
        // libgomp deadlock the level-1 fork already documents -- the child's
        // first parallel region waits forever on threads it did not inherit.
        // This project has already paid for that bug once, with a 13-minute
        // hang that showed up only after two clean runs.
        omp_pause_resource_all(omp_pause_hard);
        bool gchild=false; size_t mine2=0;
        if(g_l2_members.size()==1){ gchild=true; mine2=0; }   // nothing to split
        else {
            std::vector<pid_t> run2;
            auto reap=[&](size_t keep){
                while(run2.size() >= keep && !run2.empty()){
                    int st=0; pid_t d=wait(&st);
                    if(d<=0) break;
                    run2.erase(std::remove(run2.begin(),run2.end(),d),run2.end());
                    if(!WIFEXITED(st)||WEXITSTATUS(st)!=0)
                        fprintf(stderr,"  [a3] a candidate failed (status %d)\n",st);
                }
            };
            for(size_t j=0;j<g_l2_members.size() && !gchild;++j){
                pid_t pid=fork();
                if(pid<0){ perror("fork"); return 1; }
                if(pid==0){
                    gchild=true; mine2=j;
                    omp_set_num_threads((int)g_l2_threads);
                    g_thr_budget = g_l2_threads;            // bounds raw threads too
                    { char bb[16]; snprintf(bb,sizeof bb,"%u",g_l2_threads);
                      setenv("OMP_NUM_THREADS", bb, 1); }
                    break;
                }
                run2.push_back(pid);
                reap(g_l2_conc);
            }
            if(!gchild){
                reap(1);
                // This group child produced no archive of its own -- its
                // members did. _exit, not exit: it shares the parent's stdio
                // buffers and must neither flush them nor run destructors.
                fflush(stderr); _exit(0);
            }
        }
        // grandchild: adopt one candidate. Everything per-candidate that the
        // level-1 fork used to set is set HERE instead, including the chdir --
        // which is why level 1 must not chdir, and why it was verified that
        // nothing is written to the cwd between the two fork points.
        const size_t ci = g_l2_members[mine2];
        const std::string ddir = g_l2_base + ".cand" + std::to_string(ci) + ".d";
        mkdir(ddir.c_str(), 0755);
        setenv("MAXMAP",  std::to_string(g_l2_maxmaps[mine2]).c_str(), 1);
        setenv("ARCHIVE", (g_l2_base + ".cand" + std::to_string(ci)).c_str(), 1);
        // CAPS_CALL outputs have the SAME hazard the dump files had: every
        // candidate writes the caller's VCF (and contig dump) to one shared
        // path, so what survives belongs to whichever ran LAST, not to the
        // winner -- and which that is depends on scheduling, so the result was
        // not even reproducible run to run. Each candidate therefore writes its
        // own file at an absolute path resolved from the cwd captured BEFORE
        // any chdir (so a relative path the caller passed still refers to the
        // directory they meant), and the parent promotes the winner's.
        if(CAPS_CALL){
            auto abso=[&](const char* q)->std::string{
                std::string t(q);
                return (!t.empty() && t[0]=='/') ? t : (g_l2_cwd + "/" + t);
            };
            if(const char* cv = getenv("CALL_VCF"))
                setenv("CALL_VCF", (abso(cv) + ".cand" + std::to_string(ci)).c_str(), 1);
            else
                setenv("CALL_VCF", (g_l2_cwd + "/out.vcf.cand" + std::to_string(ci)).c_str(), 1);
            if(const char* dc = getenv("CAPS_DUMP_CONTIGS"))
                setenv("CAPS_DUMP_CONTIGS", (abso(dc) + ".cand" + std::to_string(ci)).c_str(), 1);
        }
        if(chdir(ddir.c_str())!=0) perror("chdir");

        // APPLY THIS CANDIDATE'S CEILING to the shared search. The search ran at
        // the group's maximum, so a read placed with more mismatches than this
        // candidate allows must be un-placed -- restoring exactly the state it
        // would have had if the search had run at this ceiling: readMM 255 and
        // matched 0 (as initialised), prc 0 (prc.assign(n,0)) and ppos
        // UINT64_MAX (ppos.assign(n,UINT64_MAX)), which the second-region chain
        // emission then sets for every admitted read. MAXMAP itself is not
        // re-read: it is never used after the mapping stage.
        {
            const uint32_t mycap = g_l2_maxmaps[mine2];
            size_t dropped=0;
            for(uint32_t rid : leftovers)
                if(matched[rid] && readMM[rid] > mycap){
                    matched[rid]=0; readMM[rid]=255; prc[rid]=0;
                    ppos[rid]=UINT64_MAX; --n_matched; ++dropped;
                }
            fprintf(stderr,"  [a3] cand %zu: MAXMAP=%u un-placed %zu of the shared "
                           "search's hits -> mapped=%zu\n", ci, mycap, dropped, n_matched);
        }
    }

    // Survivors get their own pseudogenome rather than being appended whole.
    // They are the reads that neither tiled well enough to assemble in round 2
    // nor matched anything already built -- but they still overlap EACH OTHER,
    // so a second sweep over just this set recovers most of their length.
    // Appending them raw costs a full read apiece (measured: 43,075 x 150 =
    // 6.46 MB); PgRC2 assembles its 11,298 survivors into 1,504,035 bytes.
    size_t appended=0, second_pg=0;
    const size_t main_pg_end = pg.size();
    // A REFERENCE NEEDS dst, src AND len. refc::encode stores only src
    // (`enc.encode(r.src, bound(r.dst))`), and the destination gaps and
    // lengths were computed beside it and then explicitly dropped as
    // "diagnostics only". Nothing else carried them: the literal stream holds
    // ACGT and nothing else, so unlike PgRC2 -- whose destinations ride in-band
    // as a MATCH_MARK inside the literal -- ours cannot be implicit. The
    // archive was therefore missing the data needed to place a single
    // reference, and every reported size was understated by that amount.
    // The read-level lossless check never caught it because it decodes the
    // dumped mem_triples.bin, which has all four fields.
    std::vector<uint8_t> ref_gaps, ref_lens, ref_rc;
    // Extension-mismatch streams: per-ref count, varint offsets (in the
    // same dst order as mem_triples), and the observed base at each --
    // coded against a `ref` side that is never stored, exactly like the
    // per-read mismatch mechanism: at decode time ref = pg[dst+offset]
    // read immediately after the reference copy, before this override.
    std::vector<uint8_t> ref_mmcnt, ref_mmpos, ref_mmref, ref_mmobs;

    {
        std::vector<uint8_t> keep(n,0);
        for(uint32_t rid:leftovers) if(!matched[rid]){ keep[rid]=1; ++appended; }
        admit.swap(keep);
        // Rebuild the prefix index over the survivors only (see the free above).
        buildPref(true);
        sweep();
        const size_t before=pg.size();
        for(uint32_t i=0;i<n;++i){
            if(!admit[i]||prv[i]!=NONE) continue;            // not a chain head here
            const uint64_t _cspan0=pg.size();
            uint32_t cur=i; ppos[cur]=pg.size(); rappend(pg,cur,0);
            while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                                   ppos[cur]=pg.size()-o; rappend(pg,cur,o); }
            if(CAPS_CALL) g_contig_spans.push_back({_cspan0,pg.size()});
        }
        second_pg=pg.size()-before;
        fprintf(stderr,"second pg: %zu reads -> %zu B (raw would be %zu B)\n",
                appended,second_pg,appended*(size_t)Lmax);
    }
    fprintf(stderr,"leftovers=%zu mapped=%zu appended=%zu\n",
            leftovers.size(),n_matched,appended);

    // ── Claim 2: reference-free variant calling (CAPS_CALL=1) ────────────────
    // g_contig_spans was captured pre-MEM in the two chain-emission loops above
    // (main pg + second region); pigeonhole-mapped reads land inside an
    // already-emitted span by construction, so every unique read's ppos/prc
    // resolves to exactly one contig via a binary search over span starts.
    if(CAPS_CALL){
        phase("pre-call");
        capscall::CallData cd;
        cd.contigs.reserve(g_contig_spans.size());
        for(auto& sp:g_contig_spans) cd.contigs.push_back(pg.substr(sp.first, sp.second-sp.first));
        std::vector<uint64_t> span_starts(g_contig_spans.size());
        for(size_t k=0;k<g_contig_spans.size();++k) span_starts[k]=g_contig_spans[k].first;
        auto find_cid=[&](uint64_t p)->uint32_t{
            auto it=std::upper_bound(span_starts.begin(),span_starts.end(),p);
            if(it==span_starts.begin()) return UINT32_MAX;
            size_t ci=(size_t)(it-span_starts.begin())-1;
            if(p>=g_contig_spans[ci].second) return UINT32_MAX;   // gap: not covered (shouldn't happen)
            return (uint32_t)ci;
        };
        std::vector<uint32_t> uid_cid(n,UINT32_MAX), uid_pos(n,0);
        for(uint32_t u=0;u<n;++u){
            if(ppos[u]==UINT64_MAX) continue;
            uint32_t ci=find_cid(ppos[u]);
            if(ci==UINT32_MAX) continue;
            uid_cid[u]=ci; uid_pos[u]=(uint32_t)(ppos[u]-g_contig_spans[ci].first);
        }
        const size_t n_orig=orig2uid.size();
        cd.read_cid.resize(n_orig); cd.read_pos.resize(n_orig); cd.read_rc.resize(n_orig);
        cd.read_ppos.resize(n_orig);
        for(size_t o=0;o<n_orig;++o){
            uint32_t u=orig2uid[o];
            cd.read_cid[o]=(u<n)?uid_cid[u]:UINT32_MAX;
            cd.read_pos[o]=(u<n)?uid_pos[u]:0;
            cd.read_rc[o] =(u<n)?prc[u]:0;
            // ppos[u] is the read's offset in the pseudogenome -- the number the
            // compressor computed in order to store the read as a position
            // instead of a sequence. Carrying it costs one vector; it is the
            // only global coordinate that exists for a read.
            cd.read_ppos[o]=(u<n)?ppos[u]:UINT64_MAX;
        }
        cd.valid=true;

        // Second FASTQ pass: original reads in original order, seq+qual only,
        // applying the SAME >1023bp skip the initial load pass uses so indices
        // stay aligned with orig2uid (see the load loop's `if(b.size()>1023)
        // continue;`). Loaded fresh rather than cached from the first pass to
        // keep the no-CAPS_CALL path's memory footprint completely unaffected.
        std::vector<std::string> call_seqs, call_quals;
        call_seqs.reserve(n_orig); call_quals.reserve(n_orig);
        // The graph-only caller consumes quality as a pass/fail bitmap; the
        // pileup path needs real phred characters, so this is gated on the
        // same flag the caller checks.
        const bool QUAL_BITMAP = getenv("CAPS_DBG_ONLY") && getenv("CAPS_DBG");
        const bool SEQ_PACK    = QUAL_BITMAP;   // same gate: graph-only caller
        const int  QMIN = getenv("CAPS_DBG_MINQ") ? atoi(getenv("CAPS_DBG_MINQ")) : 20;
        {
            // std::getline, NOT fgets with a fixed buffer: the main load pass
            // uses getline, and any read line longer than the buffer would be
            // split across two fgets calls, desynchronising the 4-line rhythm
            // and silently misindexing every subsequent read against
            // orig2uid. Same failure class as the four silent data-loss bugs
            // this project already found the hard way -- match the load pass
            // exactly rather than approximately.
            std::ifstream fin(g_input_path);
            if(!fin){ fprintf(stderr,"caps_caller: cannot reopen %s\n",g_input_path.c_str()); }
            else{
                std::string a,b,c,d;
                while(std::getline(fin,a)&&std::getline(fin,b)&&std::getline(fin,c)&&std::getline(fin,d)){
                    while(!b.empty()&&(b.back()=='\n'||b.back()=='\r')) b.pop_back();
                    while(!d.empty()&&(d.back()=='\n'||d.back()=='\r')) d.pop_back();
                    if(b.size()>1023) continue;                // mirror the load-pass skip
                    if(SEQ_PACK){
                        // capspack::pack_seq is the SINGLE definition of this
                        // format. capsule_decode rebuilds the same bytes from a
                        // stored archive, and if the two ever disagree by a byte
                        // the caller silently sees different reads. This project
                        // has already lost data to exactly that class of bug --
                        // packing that mapped non-ACGT to A added 2,273 k-mers
                        // and no F1 check could see it. Shared, never copied.
                        call_seqs.push_back(capspack::pack_seq(b));
                    } else {
                        call_seqs.push_back(b);
                    }
                    if(QUAL_BITMAP){
                        call_quals.push_back(capspack::pack_qual(d, QMIN));
                    } else {
                        call_quals.push_back(d);
                    }
                }
            }
        }
        if(call_seqs.size()!=n_orig)
            fprintf(stderr,"caps_caller: WARNING seq count %zu != orig2uid %zu (index skew, calls may be wrong)\n",
                    call_seqs.size(),n_orig);
        std::string call_vcf = getenv("CALL_VCF") ? getenv("CALL_VCF") : "out.vcf";
        int n_calls = capscall::run_variant_call(call_seqs, call_quals, cd, call_vcf);
        fprintf(stderr,"[CAPS-CALL] %d records -> %s\n", n_calls, call_vcf.c_str());
        phase("call");
    }

    // STAGE 100 -- ported from 47_mismatch_coder.cpp (58-line diff from this
    // file's base, 46_position_stream.cpp): dumps the real (ref, obs, pos,
    // ctx3) mismatch symbol streams instead of the estimated bits/mismatch
    // number used everywhere else this session. Gated on DUMP_MM, same
    // convention as DUMP_LIT/DUMP_PERM.
    if(getenv("DUMP_MM")){
        auto comp=[](char c)->char{ return c=='A'?'T':c=='T'?'A':c=='C'?'G':c=='G'?'C':'N'; };
        auto rbase=[&](uint32_t i,uint32_t j)->char{
            static const char L[4]={'A','C','G','T'};
            return L[(int)rseed(i,j,1)];
        };
        FILE* fr=STR.mm_ref.f, *fo=STR.mm_obs.f, *fp=STR.mm_pos.f;
        FILE* fx=fopen("/dev/null","wb");   // mm_ctx3 is diagnostic only, never decoded
        // LAYER 9 -- real gap found while designing the decoder: readMM[i]
        // (used only for reporting until now) is not guaranteed to equal
        // the ACTUAL number of mismatches emitted below (the dump loop
        // guards and skips reads with an out-of-range q, which readMM's
        // own count doesn't know about) -- using readMM directly for decode
        // would risk silently misaligning the mm_ref/obs/pos/ctx3 streams
        // the moment such a skip happens. Track the real, actually-emitted
        // per-read count instead, in the SAME raw-unique-id order (0..n-1)
        // this loop already iterates in -- so decode can walk mm_ref/obs in
        // lockstep with this count array and always know exactly how many
        // entries belong to each read. Also: a read needing zero corrections
        // (never mapped, or mapped with zero real mismatches) decodes
        // identically either way (just use the pg slice as-is) -- so this
        // single count array, with 0 covering both cases, is sufficient;
        // readMM's separate "255 = never mapped" sentinel isn't needed for
        // reconstruction at all.
        std::vector<uint16_t> mmcount(n,0);
        size_t emitted=0;
        // MISMATCH POSITIONS as REVERSE OFFSETS (PgRC2's real technique --
        // SeparatedPseudoGenomePersistence.cpp:823-905 codes rlMisRevOffDest,
        // an offset from the previous mismatch, never an absolute position).
        // Within a read the positions are strictly increasing, so the gap is
        // always smaller than the absolute value and the alphabet collapses
        // toward zero -- measured -10.4% on this stream in a direct prototype.
        // Guarded: only when Lmax<=256 is an absolute j always <256, which
        // makes gaps unambiguous under the existing 255 cap. Above that the
        // old absolute encoding is kept so the pre-existing >255bp cap
        // limitation behaves exactly as before. The mode is derivable by the
        // decoder from read_lengths.bin (max length), so no extra stream and
        // no format flag is needed.
        const bool MMDELTA = (Lmax<=256);
        // A unique read that NO original read maps to is never emitted by the
        // decoder, but its mismatches still occupied slots in mm_ref/mm_obs/
        // mm_pos -- and the decoder, having no original to take a length from,
        // computed its reverse-complement index with RL=0 and so derived the
        // WRONG `ref` byte. Since the mismatch coder is adaptive, those few
        // wrong refs desynchronised the model and corrupted every symbol after
        // them: on C. jejuni 8 bad refs of 18,774 turned into 6,892 wrong obs
        // bytes and 3,712 wrong reads. Such uniques arise from containment (a
        // shorter read absorbed into a longer one), which only happens with
        // variable-length input -- which is exactly why every fixed-length
        // dataset round-tripped and three variable-length ones did not.
        //
        // Emitting mismatches for a read nobody outputs is pure waste as well
        // as a correctness break, so they are skipped here. Files with no
        // orphaned uniques produce byte-identical streams.
        std::vector<uint8_t> uid_referenced(n, 0);
        for(uint32_t o : orig2uid) if(o < n) uid_referenced[o] = 1;
        size_t g_skipped_orphan=0;
        for(uint32_t i=0;i<n;++i){
            if(!uid_referenced[i]){ ++g_skipped_orphan; continue; }
            if(readMM[i]==255) continue;
            uint32_t prevj=0;
            const uint32_t RL=rlen[i];
            int64_t q; bool rc=prc[i];
            if(rc) q=(int64_t)main_pg_end-(int64_t)ppos[i]-(int64_t)RL;
            else   q=(int64_t)ppos[i];
            if(q<0||q+RL>main_pg_end) continue;               // guard: skip if out of range
            for(uint32_t j=0;j<RL;++j){
                char refc, obsc;
                if(!rc){ refc=pg[q+j]; obsc=rbase(i,j); }
                else   { refc=pg[q+RL-1-j]; obsc=comp(rbase(i,j)); }
                if(refc==obsc) continue;
                const char prevc = (j>0 && !rc) ? pg[q+j-1] :
                                    (j>0 &&  rc) ? comp(pg[q+RL-j]) : 'A';
                fputc(refc,fr); fputc(obsc,fo);
                // Above 256 the old code wrote `j>255?255:j` -- it CLAMPED,
                // silently emitting a wrong position for every mismatch past
                // base 255 and making the archive lossy for any read longer
                // than 256. Two of the locked datasets are (ERR552797 at 301
                // bp, SRR40271341 at 300 bp) and both failed to round trip.
                // Now the same delta is written, as a varint instead of a
                // byte, so any position is representable. The <=256 path is
                // untouched, so every existing archive stays byte-identical.
                if(MMDELTA){ fputc((uint8_t)(j-prevj),fp); prevj=j; }
                else { uint32_t d=j-prevj; prevj=j;
                       while(d>=0x80){ fputc((uint8_t)(d|0x80),fp); d>>=7; }
                       fputc((uint8_t)d,fp); }
                fputc(prevc,fx);
                ++emitted; ++mmcount[i];
            }
        }
        STR.mm_ref.close(); STR.mm_obs.close(); STR.mm_pos.close(); fclose(fx);
        { fwrite(mmcount.data(),2,mmcount.size(),STR.mm_count.f); STR.mm_count.close(); }
        fprintf(stderr,"[MM-DUMP] emitted=%zu mismatch symbols (readMM total was %zu -- reads with no valid q excluded)\n",
                emitted,mm_total);
    }
    {
        mm_total=0; mm_reads=0; std::fill(mm_hist.begin(),mm_hist.end(),0);
        for(uint32_t i=0;i<n;++i) if(readMM[i]!=255){ mm_total+=readMM[i]; ++mm_reads; ++mm_hist[readMM[i]]; }
        size_t zero=mm_hist[0];
        fprintf(stderr,"[MM] placed=%zu total_mismatches=%zu mean=%.2f/read  zero-mm=%zu (%.1f%%)\n",
                mm_reads,mm_total,mm_reads?(double)mm_total/mm_reads:0.0,
                zero,mm_reads?100.0*zero/mm_reads:0.0);
        fprintf(stderr,"[MM] PgRC2 for reference: 1,369,413 mismatches over 988,702 reads = 1.39/read\n");
        fprintf(stderr,"[MM] hist:"); for(int i=0;i<12;++i) fprintf(stderr," %d:%zu",i,mm_hist[i]);
        size_t tail=0; for(int i=12;i<256;++i) tail+=mm_hist[i];
        fprintf(stderr,"  >=12:%zu\n",tail);
    }

    // ── STAGE 22: the order permutation (their stage 6) ─────────────────────
    // To restore the original file order the archive must carry, for every
    // original read, where it sits in emission order. PgRC2 stores the inverse
    // permutation as a raw uint32 array and LZMAs it
    // (SeparatedPseudoGenomePersistence.cpp:226-233) -- no permutation-specific
    // coding at all in single-end mode.
    //
    // Measured, that costs them 7,063,459 - 4,898,620 = 2,164,839 B on this
    // input. The floor for a RANDOM permutation of 999,340 elements is
    // log2(n!) = 2,309,000 B -- so their output is 144 KB BELOW the random
    // floor, which means the permutation is not random and LZMA is finding real
    // structure in it. Any coder that assumes uniformity (Lehmer / factorial
    // base) would therefore do WORSE than what they already have, not better.
    // That is why this stage dumps the permutation for analysis instead of
    // going straight to a clever coder.
    //
    // Note positions from the RC pass are in reverse-complement coordinates.
    // That is fine here: only the ORDER matters for this measurement, and a
    // real encoder would store the strand flag alongside (ARCS already does).
    // STAGE 100 REWRITE -- scrapped the whole rank/bucket/uidOrder scheme
    // after reading PgRC2's real source line by line
    // (SeparatedPseudoGenomePersistence.cpp:446, compressReadsPgPositions,
    // singleFileMode branch -- our exact case, single-end). Their design for
    // SE mode is a single flat array: position[original_read_index] =
    // absolute pg position, written directly in original-read order. No
    // rank layer, no duplicate-count bucket, no separate order/permutation
    // stream at all -- position alone restores both WHERE a read's sequence
    // is AND, by construction (reads emitted in original index order),
    // its place in the file. This replaces layers 2 (perm.u32), 4, 7, and
    // 8 (orig2uid_ranks.bin/uidorder.bin) all at once with one simpler,
    // directly-indexed design that has no indirection left to get wrong.
    if(getenv("DUMP_PERM")){
        const uint64_t PL=main_pg_end;   // RC scan ran before survivors were appended
        // STAGE 100 SIZE FIX -- position/strand/length were being stored once
        // PER ORIGINAL READ (999,308 for E. coli), redundantly repeating the
        // exact same value for every duplicate of a unique sequence (150,199
        // of them). `orig2uid.bin` (original -> unique-id) is needed anyway
        // for mismatch correlation and already gives decode everything it
        // needs to look up ANY original read's data -- so position/strand/
        // length only need to exist ONCE PER UNIQUE READ (849,109 here, real
        // ~15% fewer entries, before even counting that xz now sees no
        // redundant repeated values to (imperfectly) find on its own).
        // Confirmed this doesn't reintroduce the old rank/bucket bug: unlike
        // that scheme, orig2uid.bin stores the FINAL ANSWER directly (which
        // unique-id), no further indirection or missing bucket-size info
        // needed -- this is exactly the design already proven byte-identical
        // on E. coli, just no longer wastefully duplicated per original read.
        // STAGE 100 SIZE FIX ROUND 2 -- TRIED AND REVERTED. Sorting by
        // position + delta-coding did shrink pos_abs.bin back to the old
        // scheme's size (526,888 B, confirmed) -- but required a new
        // rank_of_uid[] permutation array to stay invertible, and that
        // permutation cost 2,275,236 B (close to log2(n!) raw entropy,
        // since assembly order has little correlation with position-sort
        // order). Net effect on real E. coli data: 5,784,467 -> 5,790,447,
        // essentially flat (+0.1%) -- the cost was relocated, not removed.
        // This is the same fundamental cost the old broken perm.u32 scheme
        // was paying, just moved to a different array. Reverted to the
        // simpler direct-indexing version below (no permutation array
        // needed at all), which gives the same real total with less
        // complexity and no permutation-cost risk on other datasets.
        // STAGE 100 SIZE FIX ROUND 3 -- real, established, "not flashy"
        // technique (confirmed this session: DEFLATE itself picks per-block
        // between stored/static-Huffman/dynamic-Huffman, keeping whichever
        // is smaller -- the exact same "try a plain encoding, keep it if
        // smaller" philosophy, decades-proven, not invented here). Tested
        // directly on real data: fixed-width uint32 positions, xz'd,
        // consistently beat the varint encoding this stream used before --
        // E. coli 2,549,216 vs 2,766,428 (-7.9%), P. aeruginosa 232,668 vs
        // 261,680 (-11.1%). Varint's per-value length-prefix bits apparently
        // break up the byte-alignment patterns xz's LZ77 stage would
        // otherwise find across positions -- fixed-width keeps that
        // structure intact. Switched from variable-length varint to a
        // fixed uint32 array; decode_locked_seqorder.py updated to match
        // (read_u32 instead of read_varints for this stream).
        std::vector<uint32_t> dl; std::vector<uint8_t> sb; size_t bigd=0;
        std::vector<uint16_t> lenarr; lenarr.reserve(n);
        size_t placed=0;
        uint8_t acc=0; int nb=0;
        for(uint32_t u=0;u<n;++u){
            uint64_t q=0; uint8_t s=0; uint16_t L=0;
            if(ppos[u]!=UINT64_MAX){
                q=ppos[u]; s=prc[u]; L=(uint16_t)rlen[u];
                if(s){ const uint64_t rl=rlen[u]; q=(PL>=q+rl)?(PL-q-rl):0; }
                ++placed;
            }
            if(q>255) ++bigd;
            dl.push_back((uint32_t)q);
            lenarr.push_back(L);
            acc=(uint8_t)((acc<<1)|(s&1));
            if(++nb==8){ sb.push_back(acc); acc=0; nb=0; }
        }
        if(nb) sb.push_back((uint8_t)(acc<<(8-nb)));
        { fwrite(dl.data(),4,dl.size(),STR.pos_abs.f); STR.pos_abs.close(); }
        { fwrite(sb.data(),1,sb.size(),STR.pos_strand.f); STR.pos_strand.close(); }
        // MEASUREMENT ONLY, not yet wired to decode: PgRC2's own design stores
        // position DIRECTLY per ORIGINAL read (SeparatedPseudoGenomePersistence.
        // cpp:446), no unique-dedup indirection for position at all -- only our
        // orig2uid correlation (kept regardless, mismatches still need it) adds
        // indirection. Real archaea/P.aeruginosa losses on the order layer
        // suggest this direct scheme may beat ours at low duplication rate.
        // Dumping both so run_locked_seqorder.sh can measure real xz'd sizes
        // and report which wins BEFORE committing to a decoder rewrite --
        // same "try it, keep smaller if real" discipline as the varint/fixed32
        // and orig2uid-delta fixes, not applied blind this time.
        {
            std::vector<uint32_t> dd(orig2uid.size());
            std::vector<uint8_t> sd; sd.reserve((orig2uid.size()+7)/8);
            uint8_t acc2=0; int nb2=0;
            for(size_t o=0;o<orig2uid.size();++o){
                uint32_t u2=orig2uid[o];
                dd[o]=dl[u2];
                uint8_t bit=(sb[u2/8]>>(7-(u2%8)))&1;
                acc2=(uint8_t)((acc2<<1)|bit);
                if(++nb2==8){ sd.push_back(acc2); acc2=0; nb2=0; }
            }
            if(nb2) sd.push_back((uint8_t)(acc2<<(8-nb2)));
            // (pos_direct/pos_strand_direct were A/B diagnostics only -- dropped)
        }
        { // Always per ORIGINAL read -- one path, so encoder and decoder cannot
          // disagree about the indexing. A fallback here was the bug: on files
          // where containment never fired, origlen stayed empty and the encoder
          // silently reverted to per-unique while the decoder still read
          // per-original.
          fwrite(origlen.data(),2,origlen.size(),STR.read_lengths.f);
          STR.read_lengths.close(); }
        // orig2uid delta-coding: for a first-occurrence read, orig2uid[i] equals
        // the running "next new id" counter exactly (delta=0); only a duplicate's
        // back-reference deviates. Measured 96.1% zero-deltas on real data (low
        // dup-rate files), still 69-99.5% zero on the regression set -- xz over
        // this is a strict, generalizing win (34.6%-90.8% smaller), never worse,
        // because delta=0 is self-describing (unambiguously "new id, exp++") vs
        // any nonzero delta ("back-reference exp-delta, exp unchanged"). Fully
        // reversible, no ambiguity, verified on 4-file regression set + 1 new.
        {
            std::vector<int32_t> delta(orig2uid.size());
            if(getenv("DBG_O2U")){
                fprintf(stderr,"[dbg-o2u] size=%zu  first 32: ",orig2uid.size());
                for(size_t i=0;i<32 && i<orig2uid.size();++i) fprintf(stderr,"%u ",orig2uid[i]);
                fprintf(stderr,"\n");
            }
            uint32_t exp=0;
            for(size_t i=0;i<orig2uid.size();++i){
                uint32_t v=orig2uid[i];
                if(v==exp){ delta[i]=0; ++exp; }
                else delta[i]=(int32_t)exp-(int32_t)v;
            }
            fwrite(delta.data(),4,delta.size(),STR.orig2uid.f); STR.orig2uid.close();
        }
        fprintf(stderr,"[ORDER+POS] unique_placed=%zu (of %u unique reads)  positions raw=%zu B  strand raw=%zu B  lengths raw=%zu B  pos>255=%zu\n",
                placed,n,dl.size(),sb.size(),lenarr.size()*2,bigd);
        fprintf(stderr,"[POS] PgRC2 pays 683,370 B coded for its reads-list offsets\n"); phase("MEM+positions");
    }

    // Nothing below this point reads the reads, the prefix index, or any of the
    // sweep's per-read arrays -- the MEM stage works purely on `pg`. Holding
    // them through it was costing more than the stage itself uses: 851,275
    // std::strings at ~150 bases plus a 32-byte header apiece is ~168 MB, and
    // the prefix index another ~80 MB. Released before the index is built, so
    // the two allocations never coexist.
    { std::vector<uint64_t>().swap(rpk); std::vector<uint64_t>().swap(woff);
      std::vector<uint16_t>().swap(rlen); }
    { std::vector<uint64_t>().swap(pent); std::vector<uint32_t>().swap(ptab); }
    { std::vector<uint32_t>().swap(nxt);  std::vector<uint32_t>().swap(prv);
      std::vector<uint32_t>().swap(ovl);  std::vector<uint32_t>().swap(ch_h);
      std::vector<uint32_t>().swap(ch_t); std::vector<uint32_t>().swap(tails);
      std::vector<uint64_t>().swap(seed); std::vector<uint8_t>().swap(ok);
      std::vector<uint8_t>().swap(admit); std::vector<uint8_t>().swap(matched); }

    // ── MEM matching between and within pseudogenomes ────────────────────────
    // PgRC2's stage 7 (SimplePgMatcher::matchPgsInPg) and the single largest
    // thing this progression was missing. They build ONE matcher over the HQ pg
    // and then, reverse-complement matching enabled throughout:
    //   - match the LQ pg into it        (SimplePgMatcher.cpp:187)
    //   - match the N pg into it         (:194)
    //   - match the HQ pg into ITSELF    (:202)
    // replacing each maximal exact match with a mark plus (offset,length).
    //
    // Measured from their own stdout on yeast_sub.fq:
    //     HQ 21,104,500 -> 11,463,320   (45.8% removed)
    //     LQ  1,504,035 ->  1,230,006   (18.5%)
    //     N       2,100 ->      1,849   (12.1%)
    // Their real final literal is 12,695,175. The 22,610,635 this progression
    // has been calling "their 22.60M" is their pg BEFORE this stage runs.
    //
    // Their parameters at CODER_LEVEL_NORMAL (pgrc-params.h:138-145): minimum
    // match length 45, reverse complement on.
    //
    // NOT copied: they use CopMEM with sparse sampling, which their own paper
    // concedes can miss matches. This reuses the packed-32-base seed already in
    // this file and samples the SOURCE every STEP = 45-32+1 = 14. That is exact,
    // not heuristic: a match of length >= 45 spans >= 14 consecutive candidate
    // seed starts, which must contain a multiple of 14, so no qualifying match
    // can be missed.
    // STAGE 27. Whole-architecture check, not another layer tweak.
    // Stage 17 cut survivors 66,775 -> 15,296 by accepting up to readLen/3 = 50
    // mismatches and taking the FIRST acceptable placement. PgRC2 allows the
    // same 50 but re-scans and keeps the BEST (ReadsMatchers.cpp:315-328, the
    // "better-matches" column). Their accepted matches average 1,369,413
    // mismatches over 988,702 reads = 1.39 per read, coded to 265,900 B.
    //
    // Every mismatch we accept has to be STORED. This progression measures
    // pg literal and order and never counted that stream, so stage 17's win
    // could be partly borrowed against a cost that was never on the books.
    // This measures it before anything else is built on top.
    const size_t MINMEM  = (argc>9)?(size_t)atoi(argv[9]):45;
    const size_t MEMSEED = (MINMEM>14)?((MINMEM-14>32)?32:MINMEM-14):8;
    // Extension mismatches: a match that would otherwise TRUNCATE at the first
    // differing base instead continues past it, up to REF_MAXMM times, storing
    // WHERE it differed. What is stored there (the observed base) is never
    // written separately -- it is already sitting in `pg` at dst+offset, the
    // same way per-read mismatches never store `ref` because it is derivable
    // from the reconstructed pg. Offsets are destination-relative (0..len-1)
    // in FINAL genome-coordinate order, so decode can read them directly
    // against the already-copied (pre-override) pg content.
    static const int REF_MAXMM = 4;
    struct Ref {
        uint32_t dst, src, len; bool is_rc;
        uint8_t mmcnt = 0;
        uint32_t mmpos[REF_MAXMM] = {0,0,0,0};
    };
    std::vector<Ref> allrefs;
    size_t rem_main=0, rem_second=0, nm_main=0, nm_second=0;
    {
        // MINMEM was 45 because that is PgRC2's CODER_LEVEL_NORMAL default
        // (pgrc-params.h:145), copied without ever being swept here. With the
        // reference cost finally measured, the trade is computable: a match
        // costs ~40.5 bits (24.5 source + 7.4 gap + 8.6 length = 5.07 B) and
        // saves L * 1.9174/8 bytes of literal, so it pays whenever L >= ~21.
        // 45 is more than twice that.
        //
        // The memory objection dissolves too. STEP = MINMEM - SEED + 1, so
        // shrinking the seed alongside the threshold holds STEP -- and therefore
        // the index size -- constant: MINMEM 30 with a 16-base seed gives the
        // same STEP 15 as MINMEM 45 with 32.
        const size_t STEP     = MINMEM - MEMSEED + 1;
        auto packM=[&](const char* q,uint64_t& o)->bool{
            uint64_t k=0;
            for(size_t i=0;i<MEMSEED;++i){ int v=b2(q[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
            o=k; return true; };
        fprintf(stderr,"[MEM] MINMEM=%zu seed=%zu step=%zu\n",MINMEM,MEMSEED,STEP); phase("assemble+map");
        // Candidates are scanned in index order and the LONGEST is kept, so a cap
        // does not merely bound work -- it can hide the best match entirely. Our
        // matches average 179 bases at MINMEM 45 against copMEM's 250, and a
        // truncated candidate list is the obvious suspect.
        const size_t MAXCAND = (argc>10)?(size_t)atoi(argv[10]):64;
        const bool   LAZY    = (argc>11)?(atoi(argv[11])!=0):true;
        const bool   FWD_SELF=(getenv("FWD_SELF")?atoi(getenv("FWD_SELF"))!=0:true);

        // Flat CSR index instead of unordered_map<uint64_t,vector<uint32_t>>.
        // That map allocated a separate heap block per distinct key -- ~1.5M of
        // them here -- and paid a node header on each. Sorting (key,pos) pairs
        // once and pointing an open-addressing table at the first occurrence of
        // each key gives the same answers from two flat arrays: all occurrences
        // of a key are adjacent after the sort, so a lookup is one probe plus a
        // forward scan. This is the same exact-size-CSR fix already carried into
        // ARCS from this progression.
        lap("  [diag] entering MEM matching, before seed index");
        // The seed index was built inline over pg[0, main_pg_end) and could not
        // be rebuilt, so no pass could ever use a different source text. Wrapped
        // as a lambda over (T, tlen); calling it with pg.data()/main_pg_end
        // reproduces the original index exactly.
        std::vector<uint64_t> skey; std::vector<uint32_t> spos;
        std::vector<uint32_t> htab;
        size_t tsize=1; uint64_t TMASK=0;
        auto hmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
        auto build_index=[&](const char* T,size_t tlen){
            std::vector<std::pair<uint64_t,uint32_t>> tmp;
            tmp.reserve(tlen/STEP+16);
            for(size_t p=0;p+MEMSEED<=tlen;p+=STEP){
                uint64_t k; if(packM(T+p,k)) tmp.push_back({k,(uint32_t)p});
            }
            std::sort(tmp.begin(),tmp.end());
            skey.clear(); skey.shrink_to_fit(); skey.reserve(tmp.size());
            hugehint((void*)skey.data(), tmp.size()*8); skey.resize(tmp.size());
            spos.clear(); spos.shrink_to_fit(); spos.reserve(tmp.size());
            hugehint((void*)spos.data(), tmp.size()*4); spos.resize(tmp.size());
            hugehint(skey.data(), skey.size()*8); hugehint(spos.data(), spos.size()*4);
            for(size_t i=0;i<tmp.size();++i){ skey[i]=tmp[i].first; spos[i]=tmp[i].second; }
            tsize=1; while(tsize < skey.size()*2+1) tsize<<=1;
            TMASK=tsize-1;
            hugefill(htab,tsize,UINT32_MAX);
            hugehint(htab.data(), htab.size()*4);
            for(size_t i=0;i<skey.size();++i){
                if(i && skey[i]==skey[i-1]) continue;   // first occurrence only
                size_t h=hmix(skey[i])&TMASK;
                while(htab[h]!=UINT32_MAX) h=(h+1)&TMASK;
                htab[h]=(uint32_t)i;
            }
        };
        build_index(pg.data(), main_pg_end);
        auto lookup=[&](uint64_t k)->uint32_t{
            size_t h=hmix(k)&TMASK;
            while(htab[h]!=UINT32_MAX){ if(skey[htab[h]]==k) return htab[h]; h=(h+1)&TMASK; }
            return UINT32_MAX;
        };
        lap("  [diag] seed index built");

        enum Mode { CROSS, SELF_FWD, SELF_RC };
        // DBG_FAILMODE: classify every position that does NOT end in an
        // accepted match, so the real bottleneck is measured instead of
        // guessed. Atomics because parse_range runs multi-threaded; only
        // touched when the env var is set, so zero cost otherwise.
        static std::atomic<size_t> dbgNoKmer{0}, dbgNoSeedHit{0}, dbgCandExceeded{0}, dbgSeedButShort{0};
        const bool DBGFAIL = getenv("DBG_FAILMODE") != nullptr;
        const bool COSTGATE = getenv("COST_GATE") != nullptr;
        // Greedy left-to-right parse over ONE slice of the destination. Matches
        // are collected rather than written straight into the shared consumed
        // bitmap, so threads never touch the same bytes and the merge below is
        // an ordinary union.
        // STAGE 32. Every accepted match is a REFERENCE the archive must store:
        // where in the destination it starts, where in the source it points, and
        // how long it is. PG_LITERAL counts only the bases that survive, so that
        // cost has been invisible in every number this progression has reported.
        // PgRC2 pays 177,180 B for exactly these streams (offsets 114,824 +
        // lengths 45,485 + 13,847 + 3,024), and their raw 155,048 B offset stream
        // implies ~38,762 matches against our 59,013 -- so ours is plausibly
        // larger, on the order of the whole remaining deficit.
        //
        // Same failure stage 27 caught with mismatches: an aggressive stage looks
        // free because its cost lands somewhere nobody measures.
        // The source text was hardcoded to pg/main_pg_end in three places, so a
        // pass could only ever match against the MAIN region. Making it explicit
        // is a prerequisite for self-matching the second region; passing
        // pg.data()/main_pg_end reproduces today's behaviour exactly.
        auto parse_range=[&](const char* Q,size_t qlen,Mode mode,size_t lo,size_t hi,
                             std::vector<Ref>& out,
                             const char* S=nullptr,size_t slen=0,int MAXMM=0){
            if(!S){ S=pg.data(); slen=main_pg_end; }
            // Tolerant extension: keep going past a mismatch, up to MAXMM times,
            // recording each offset (query-relative, i.e. relative to qp -- run()
            // re-expresses these in destination-relative terms once dst is
            // final). A trailing mismatch that nothing profitable follows is
            // trimmed back off, so a match never gains length purely by ending
            // on a wasted substitution. MAXMM=0 reproduces the old exact-only
            // behaviour exactly -- verified byte-identical before any non-zero
            // MAXMM was enabled anywhere.
            auto extendTol=[&](const char* Q,size_t qlen,const char* S,size_t slen,
                               size_t qp,size_t s,size_t capL,int maxmm,
                               uint32_t* mmout,uint8_t& mmcnt)->size_t{
                // ── WORD-AT-A-TIME EXTENSION ────────────────────────────
                // This compared ONE BYTE per iteration while the pigeonhole
                // mismatch loop next door already compares 32 bases at a time
                // (w32 + popcount). perf puts 36.2% of all cycles in this
                // worker, MINMEM is 45, and matches routinely run into the
                // hundreds of bases -- so the overwhelming majority of these
                // iterations were single-byte compares of bases that MATCH.
                //
                // Q and S are plain ASCII here, so eight bases fit in a uint64
                // and a matching stretch advances eight at a time. On the
                // first differing word, ctz gives the byte index directly
                // (x86 is little-endian) and the original per-byte path takes
                // over for the mismatch itself.
                //
                // OUTPUT-PRESERVING: the same L is returned and the same
                // mismatch positions are recorded in the same order; only the
                // scanning of EQUAL bytes is done in wider steps.
                size_t L=0; mmcnt=0;
                for(;;){
                    size_t lim=capL; 
                    if(qlen-qp<lim) lim=qlen-qp;
                    if(slen-s  <lim) lim=slen-s;
                    if(L>=lim) break;
                    while(L+8<=lim){
                        uint64_t a,b; memcpy(&a,Q+qp+L,8); memcpy(&b,S+s+L,8);
                        const uint64_t d=a^b;
                        if(d){ L += (size_t)(__builtin_ctzll(d)>>3); break; }
                        L+=8;
                    }
                    if(L>=lim) break;
                    if(Q[qp+L]==S[s+L]){ ++L; continue; }
                    if((int)mmcnt>=maxmm) break;
                    mmout[mmcnt++]=(uint32_t)L; ++L;
                }
                while(mmcnt>0 && mmout[mmcnt-1]==L-1){ --mmcnt; --L; }
                return L;
            };
            // STAGE 34: lazy matching.
            // Taking the longest match at every position, left to right, is
            // GREEDY parsing, and the classic LZ result is that greedy is not
            // optimal: accepting a shorter match now can block a much longer one
            // starting a base later, and the fragments cost a reference each.
            // References are our only remaining loss (66,831 against their
            // ~43,190) and each costs ~40 bits, so fewer-and-longer is exactly
            // the objective.
            //
            // Lazy matching is the cheap form of lookahead, as used in DEFLATE:
            // having found a match at qp, look at qp+1, and if that is strictly
            // longer, leave qp as a literal and take the better one. One extra
            // probe per matched position, no extra memory, and nothing from
            // PgRC2 -- their parse is greedy too.
            auto bestAt=[&](const char* Q,size_t qlen,Mode mode,size_t qp,
                            size_t& bsrc)->size_t{
                uint64_t k; bsrc=0;
                if(!packM(Q+qp,k)) return 0;
                const uint32_t idx=lookup(k);
                if(idx==UINT32_MAX) return 0;
                size_t best=0, tried=0;
                for(uint32_t i=idx;i<skey.size()&&skey[i]==k;++i){
                    if(++tried>MAXCAND) break;
                    const size_t s=spos[i];
                    size_t capL=SIZE_MAX;
                    if(mode==SELF_FWD){ if(s>=qp) continue; capL=qp-s; }
                    else if(mode==SELF_RC){ if(s>=qlen-qp) continue; capL=(qlen-qp-s)/2; }
                    if(capL<MINMEM) continue;
                    uint32_t scratch[REF_MAXMM]; uint8_t sc=0;
                    size_t L=extendTol(Q,qlen,S,slen,qp,s,capL,MAXMM,scratch,sc);
                    if(L>best){ best=L; bsrc=s; }
                }
                return (best>=MINMEM)?best:0;
            };
            size_t qp=lo, lastend=lo;
            while(qp<hi && qp+MINMEM<=qlen){
                uint64_t k;
                if(!packM(Q+qp,k)){ if(DBGFAIL) dbgNoKmer.fetch_add(1,std::memory_order_relaxed); ++qp; continue; }
                const uint32_t idx=lookup(k);
                if(idx==UINT32_MAX){ if(DBGFAIL) dbgNoSeedHit.fetch_add(1,std::memory_order_relaxed); ++qp; continue; }
                size_t best=0, bestsrc=0, tried=0;
                uint32_t bestmm[REF_MAXMM]; uint8_t bestmmc=0;
                bool exceededCand=false;
                for(uint32_t i=idx;i<skey.size()&&skey[i]==k;++i){
                    if(++tried>MAXCAND){ exceededCand=true; break; }
                    const size_t s=spos[i];
                    // Cap usable length BEFORE extending, never after -- in a
                    // self-match the seed at qp also occurs at qp itself and an
                    // uncapped extension runs to end of text.
                    size_t capL=SIZE_MAX;
                    if(mode==SELF_FWD){ if(s>=qp) continue; capL=qp-s; }
                    else if(mode==SELF_RC){ if(s>=qlen-qp) continue; capL=(qlen-qp-s)/2; }
                    if(capL<MINMEM) continue;
                    uint32_t cmm[REF_MAXMM]; uint8_t cmmc=0;
                    size_t L=extendTol(Q,qlen,S,slen,qp,s,capL,MAXMM,cmm,cmmc);
                    if(L>best){ best=L; bestsrc=s; bestmmc=cmmc;
                        for(uint8_t t=0;t<cmmc;++t) bestmm[t]=cmm[t]; }
                }
                if(DBGFAIL && best<MINMEM){
                    if(exceededCand) dbgCandExceeded.fetch_add(1,std::memory_order_relaxed);
                    else dbgSeedButShort.fetch_add(1,std::memory_order_relaxed);
                }
                if(best>=MINMEM && LAZY && qp+1<hi && qp+1+MINMEM<=qlen){
                    size_t nsrc=0;
                    const size_t nb=bestAt(Q,qlen,mode,qp+1,nsrc);
                    if(nb>best){ ++qp; continue; }   // qp becomes literal; take the longer one next
                }
                if(best>=MINMEM){
                    // STAGE 33: extend BACKWARD before accepting.
                    // The source is sampled every STEP, so a true maximal match
                    // at source ss is only seen at the first sampled position
                    // p >= ss inside it. Extending forward from there captures a
                    // SUFFIX of the real match and abandons the leading p-ss
                    // bases -- about STEP/2 each -- to literal. That is why our
                    // matches average 189 bases against CopMEM's 250: copMEM
                    // reports maximal exact matches, this reported partial ones.
                    //
                    // Recovering them costs nothing: the match count is
                    // unchanged, so the reference streams are unchanged, and the
                    // recovered bases come straight off the literal.
                    size_t cap=qp-lastend;                       // no overlap with the previous match
                    if(bestsrc<cap) cap=bestsrc;                 // and none before the source start
                    if(mode==SELF_FWD){
                        // source must still end at or before the destination
                        const size_t slack=qp-bestsrc-best;      // >=0, enforced by capL above
                        if(slack<cap) cap=slack;
                    }
                    // SELF_RC needs no extra cap: its condition is
                    // bestsrc+best <= qlen-qp-best, which backward extension
                    // leaves unchanged on both sides.
                    size_t b=0;
                    // Backward extension, word-at-a-time (same reasoning as
                    // extendTol). Little-endian: loading 8 bytes ending at P
                    // puts P-1 in the HIGH byte, so scanning backward means
                    // finding the HIGHEST differing byte -- clz, not ctz.
                    while(b+8<=cap && qp>=b+8 && bestsrc>=b+8){
                        uint64_t x,y;
                        memcpy(&x,Q+qp-b-8,8); memcpy(&y,S+bestsrc-b-8,8);
                        const uint64_t d=x^y;
                        if(d){ b += (size_t)(7 - ((63-__builtin_clzll(d))>>3)); break; }
                        b+=8;
                    }
                    while(b<cap && Q[qp-b-1]==S[bestsrc-b-1]) ++b;
                    Ref nr; nr.dst=(uint32_t)(qp-b); nr.src=(uint32_t)(bestsrc-b);
                    nr.len=(uint32_t)(best+b); nr.is_rc=false;
                    nr.mmcnt=bestmmc;
                    for(uint8_t t=0;t<bestmmc;++t) nr.mmpos[t]=bestmm[t]+(uint32_t)b;
                    // COST-AWARE ACCEPTANCE (COST_GATE=1). MINMEM is a length
                    // proxy for "worth a reference" -- it does not know that a
                    // reference costs ~1 varint (src) + ~1 varint (len) + ~1
                    // varint (gap) + 1 bit (rc) regardless of length, while
                    // literal costs a near-fixed ~2 bits/base. A 24-base match
                    // barely clears that; a mismatch-laden extension (measured
                    // this session: +490,763 B net loss on Drosophila 2.6x)
                    // never does. This computes the SAME byte-cost formulas the
                    // real coders use (vint length, log2(dst) as an upper bound
                    // on the src coder's own bound -- true bound is <= dst in
                    // every case, self or cross, so this can only OVER-count
                    // src cost, never wrongly accept a match that is actually
                    // unprofitable) against a literal cost of 2.0 bits/base --
                    // the information-theoretic MAXIMUM for 4 symbols, so a
                    // match this gate rejects loses even under the most
                    // generous possible assumption about literal's cost.
                    bool accept = true;
                    if(COSTGATE){
                        auto vintBits=[](uint64_t v)->double{
                            size_t bytes=1; v>>=7; while(v){ ++bytes; v>>=7; } return (double)(bytes*8); };
                        const uint64_t gap = (nr.dst>=lastend)?((uint64_t)nr.dst-lastend):0;
                        const double srcBits = nr.dst>1 ? std::log2((double)nr.dst) : 1.0;
                        const double refBits = srcBits + vintBits((uint64_t)(nr.len-MINMEM)) + vintBits(gap) + 1.0;
                        const double litBitsSaved = (double)nr.len * 2.0;
                        accept = refBits < litBitsSaved;
                    }
                    if(accept){
                        out.push_back(nr);
                        lastend=qp+best; qp+=best;
                    } else ++qp;
                }
                else ++qp;
            }
        };
        // Split the destination across threads. The greedy parse is serial by
        // nature (it skips ahead by each accepted match), so what is parallelised
        // is the search, exactly as ARCS's repeat_elim does: every query reads
        // only the shared, already-built index. Slicing can only differ from a
        // single serial parse at the T-1 slice boundaries, bounded by ~45 bytes
        // each -- under 0.005% of a 21 MB pseudogenome.
        // DSTBASE: the cross-match is handed pg.data()+main_pg_end, so its qp is
        // relative to the survivor pg, not to the pseudogenome. Stored raw, a
        // cross-match destination of 0 claims to be the very start of the pg,
        // which is why 12,003 of 58,908 rows appeared to violate src < dst.
        // DSTBASE existed but SRCBASE did not: destinations were lifted into
        // pseudogenome coordinates while sources were left in whatever frame the
        // source text used. Invisible while the source IS the main pg (base 0),
        // and wrong the moment it is not.
        auto run=[&](const char* Q,size_t qlen,Mode mode,std::vector<uint8_t>& consumed,
                     bool RCDEST=false,size_t DSTBASE=0,
                     const char* S=nullptr,size_t slen=0,size_t SRCBASE=0,int MAXMM=0)->size_t{
            unsigned T=std::thread::hardware_concurrency(); if(!T) T=1;
            if(qlen < (1u<<20)) T=1;
            std::vector<std::vector<Ref>> res(T);
            std::vector<std::thread> th;
            const size_t chunk=(qlen+T-1)/T;
            // SAME T CHUNKS, FEWER WORKERS. T is left exactly as it was, so the
            // partition, the res[t] buckets and the merge order below are
            // unchanged -- the output is identical by construction. Only how
            // many chunks are in flight at once is bounded, which is what stops
            // 4 concurrent candidates putting 48 threads on 12 cores.
            const unsigned W = std::min<unsigned>(T, thr_budget());
            std::atomic<unsigned> nextc{0};
            for(unsigned w=0; w<W; ++w){
                th.emplace_back([&]{
                    for(unsigned t=nextc.fetch_add(1); t<T; t=nextc.fetch_add(1)){
                        const size_t lo=(size_t)t*chunk, hi=std::min(qlen,lo+chunk);
                        if(lo>=hi) continue;
                        parse_range(Q,qlen,mode,lo,hi,res[t],S,slen,MAXMM);
                    }
                });
            }
            for(auto& x:th) x.join();
            size_t nm=0;
            size_t capBytes=0; for(auto& v:res) capBytes+=v.capacity()*sizeof(Ref);
            fprintf(stderr,"    [diag] run() T=%u qlen=%zu res[] total capacity=%zu B (%.1f MB)\n",T,qlen,capBytes,capBytes/1e6);
            for(auto& v:res){
                nm+=v.size();
                for(auto& m:v){
                    for(uint32_t j=0;j<m.len && (size_t)m.dst+j<qlen;++j) consumed[m.dst+j]=1;
                    // An RC pass reports its destination in reverse-complement
                    // coordinates while the source stays in forward ones, so the
                    // two are not comparable until the destination is converted
                    // back: a match at RC position qp of length L covers forward
                    // [qlen-qp-L, qlen-qp). PgRC2 does the same conversion in
                    // correctDestPositionDueToRevComplMatching
                    // (SimplePgMatcher.cpp:58-60); we were storing the raw RC
                    // position, which left src and dst in different coordinate
                    // systems.
                    Ref r=m;
                    r.is_rc=RCDEST;
                    if(RCDEST){
                        r.dst=(uint32_t)(qlen-m.dst-m.len);
                        // Query-relative offset k becomes destination-relative
                        // len-1-k under reversal; decode walks the genome
                        // forward, so offsets must be re-expressed (and
                        // re-sorted ascending) in that frame, not the query's.
                        for(uint8_t t=0;t<r.mmcnt;++t) r.mmpos[t]=m.len-1-m.mmpos[t];
                        for(uint8_t a=0,b=r.mmcnt;a+1<b;++a,--b)
                            std::swap(r.mmpos[a],r.mmpos[b-1]);
                    }
                    r.dst=(uint32_t)(r.dst+DSTBASE);
                    r.src=(uint32_t)(r.src+SRCBASE);
                    allrefs.push_back(r);
                }
            }
            return nm;
        };

        // STAGE 100 REWRITE -- the two separate per-region c[]/cr[] coverage
        // bitmaps below (used only for nm_main/nm_second reporting now) were
        // ALSO independently driving literal.txt's writing -- but a decoder
        // trying to reconstruct literal.txt's coverage from the separately-
        // dumped `allrefs` triples could NOT reliably reproduce the same
        // positions (confirmed directly: byte COUNTS matched by coincidence,
        // but actual covered POSITIONS did not -- every reconstructed read
        // came out wrong, 0/15424 on a real decode attempt, even for reads
        // with zero mismatches). Root cause: `allrefs` contains redundant,
        // overlapping candidate matches (this session found 17,343 exact
        // duplicate destinations on real E. coli data) with no way to
        // recover, from the dump alone, which ones the ORIGINAL c[]/cr[]
        // merge actually kept.
        //
        // Fix, matching PgRC2's real, proven design exactly (SimplePgMatcher
        // .cpp:85-160, read line by line this session): sort ALL matches
        // (both regions, forward+RC, now unified in one `allrefs` with a
        // real is_rc flag -- previously missing entirely) by destination,
        // then walk them TRIMMING overlap with already-consumed territory
        // instead of just OR-ing a bitmap. This guarantees, by construction,
        // that literal.txt's mark positions and the emitted (dst,src,len,
        // is_rc) triples describe the EXACT same non-overlapping partition
        // of `pg` -- nothing to reconstruct or get out of sync later.
        // RAM FIX (round 2, C. elegans scale) -- found by bracketing the
        // MEM-matching phase with fine-grained RSS checkpoints: internal
        // checkpoints never exceeded 926MB, yet the kernel-tracked peak
        // (/usr/bin/time -v) was consistently 1.27-1.30GB -- a transient
        // spike DURING this block, invisible to before/after sampling.
        // Concrete, measured overlap found: `c`/`cr` (main-pg coverage
        // bitmaps, main_pg_end bytes each -- 135MB combined at C. elegans
        // scale) were never freed after their last real use (the rem_main
        // tally just below) -- they stayed allocated for the ENTIRE rest
        // of the function, fully overlapping the second pass's own `c2`/
        // `cr2`/`Q`/`R` allocations (qlen bytes/chars each -- another
        // ~360MB at this scale). Fixed by scoping `c`/`cr` so they free
        // before the second pass begins, same principle as the earlier
        // allrefs/cleanRefs fix: don't hold data alive past its last use.
        // Extension mismatch tolerance, derived from a measured property, not a
        // per-dataset switch: MEM_MAXMM=0 whenever leftover_frac is small (every
        // one of the 7 locked datasets), rising toward REF_MAXMM as more of the
        // file fails to chain at all in round 1 -- exactly the population that
        // needs a candidate overlap to survive one sequencing error to be
        // captured at all. Applied only to the main self-match and the
        // second-region CROSS match; the SECOND_SELF pass is left at 0
        // (unaffected) since it was already measured to lose economically at
        // MAXMM=0 and more matches there would only make that worse.
        // Measured on the 7 locked datasets: leftover_frac ranges 0.222-0.518.
        // It is NOT near zero at normal coverage as first assumed -- duplicate
        // reads and reads whose only overlap falls at a tail edge also end up
        // leftover, not just true low-coverage gaps. So a bare proportional
        // formula fired everywhere and changed every locked archive. Rebuilt
        // with a floor above the measured locked maximum (0.518) plus real
        // margin, and a ceiling near Drosophila's 2.6x-coverage value (0.812):
        // T0=0.60 keeps MEM_MAXMM=0 on every locked dataset with room to
        // spare; T1=0.85 is where the tolerance reaches its cap. This is still
        // a function of the measured property, not a per-dataset switch --
        // any future dataset lands wherever ITS leftover_frac places it.
        const double leftover_frac = n ? (double)leftovers.size()/n : 0.0;
        const double T0=0.60, T1=0.85;
        // MEASURED AND REJECTED: extension mismatch tolerance was built,
        // debugged (it exposed a real decoder bug -- see the streaming-decode
        // fix) and verified lossless, then measured honestly on the dataset it
        // was designed for: Drosophila SRR40104920 at 2.6x coverage went
        // 36,982,418 -> 37,473,181 B, a 490,763 B LOSS. Removing 11 M literal
        // bases saved 212,041 B while the references cost 369,424 B, because a
        // region-scale match's mismatch offset costs ~21 bits where SPRING's
        // read-scale equivalent costs ~2. The code stays (correct, and the
        // decoder bug it found was real) but it is OFF: enable with
        // MEM_MAXMM_OVERRIDE only to reproduce the measurement.
        int MEM_MAXMM = 0;
        // HETEROZYGOSITY PRE-SCAN (2026-09-03). Separate from leftover_frac
        // above, which is about COVERAGE gaps. This is about PLOIDY: a
        // haploid genome's pg, after dedup, holds each true stretch of
        // sequence ONCE -- there is nothing left to link, so MEM_MAXMM>0
        // would only cost bytes for zero calling benefit (het-SNV/indel
        // calling is not even a coherent operation on a haploid sample).
        // A heterozygous diploid genome's pg still holds BOTH haplotype
        // copies as separate literal regions (MEM_MAXMM=0 hasn't linked
        // them yet), so a canonical k-mer from their SHARED flanking
        // sequence appears roughly TWICE in pg, not once. Sampling (not
        // scanning every position) keeps this cheap: this must be "early,
        // cheap, and accurate" -- a wrong call here silently mis-gates
        // every downstream decision, so it is measured, not guessed.
        bool het_detected = false;
        double het_pair_frac = 0.0;
        {
            const size_t STRIDE = 41;   // sample every 41st position, not exhaustive
            const int K = 31;
            std::unordered_map<uint64_t,uint32_t> kc_sample;
            size_t sampled = 0;
            for (size_t i = 0; i + (size_t)K <= main_pg_end; i += STRIDE) {
                uint64_t v = 0; bool ok = true;
                for (int t = 0; t < K; ++t) { int bb = b2(pg[i+t]); if (bb < 0) { ok = false; break; } v = (v<<2) | (uint64_t)bb; }
                if (!ok) continue;
                // canonicalize against reverse complement (2*K bits, K<=31 fits in 62 bits)
                uint64_t rv = 0, vv = v;
                for (int t = 0; t < K; ++t) { rv = (rv<<2) | (3 - (vv & 3)); vv >>= 2; }
                ++kc_sample[std::min(v,rv)];
                ++sampled;
            }
            size_t pairs = 0, singles = 0, other = 0;
            for (auto& kv : kc_sample) {
                if (kv.second == 1) ++singles;
                else if (kv.second == 2) ++pairs;
                else ++other;
            }
            size_t distinct = pairs + singles + other;
            het_pair_frac = distinct ? (double)pairs / (double)distinct : 0.0;
            // Threshold, CALIBRATED against real measurement -- TWICE now,
            // both times caught by actually testing, not by trusting the
            // arithmetic. Attempt 1 (0.05) false-negatived on HG002
            // (0.0264). Attempt 2 (0.020) was a real ARITHMETIC ERROR --
            // 0.020 is BELOW M. tuberculosis's own measured 0.0209, so it
            // false-positived on a HAPLOID dataset and was caught only by
            // re-running M. tuberculosis and finding the archive had
            // actually changed (1,653,912 -> 1,600,645 B), not by
            // inspecting the number. Real measured values now on record:
            // haploid -- E. coli 0.0172, M. tuberculosis 0.0209 (the real
            // max), P. aeruginosa 0.0103; heterozygous -- HG002 0.0264.
            // 0.024 sits strictly between the two groups with real margin
            // on both sides (0.0031 below HG002, 0.0031 above TB). Still
            // n=1 positive, n=3 negative -- re-validate against HG003/4/5
            // and more haploid datasets before trusting this further.
            const double HET_THRESH = 0.024;
            het_detected = het_pair_frac >= HET_THRESH;
            fprintf(stderr,"[HETSCAN] sampled=%zu distinct=%zu pairs=%zu singles=%zu other=%zu pair_frac=%.4f het=%d\n",
                    sampled, distinct, pairs, singles, other, het_pair_frac, (int)het_detected);
        }
        if (het_detected) MEM_MAXMM = std::max(MEM_MAXMM, 2);
        if(getenv("MEM_MAXMM_OVERRIDE")) MEM_MAXMM = atoi(getenv("MEM_MAXMM_OVERRIDE"));
        fprintf(stderr,"[MMTOL] leftover_frac=%.3f het_pair_frac=%.4f -> MEM_MAXMM=%d\n", leftover_frac, het_pair_frac, MEM_MAXMM);
        {
            std::vector<uint8_t> c(main_pg_end,0), cr(main_pg_end,0);
            if(FWD_SELF) nm_main =run(pg.data(),main_pg_end,SELF_FWD,c,false,0,nullptr,0,0,MEM_MAXMM);
            { std::string R(pg,0,main_pg_end); rc_inplace(R);
              nm_main+=run(R.data(),main_pg_end,SELF_RC,cr,true,0,nullptr,0,0,MEM_MAXMM); }
            for(size_t i=0;i<main_pg_end;++i) if(cr[i]) c[main_pg_end-1-i]=1;
            for(size_t i=0;i<main_pg_end;++i) if(c[i]) ++rem_main;
        }   // c, cr freed here -- before the second pass allocates its own bitmaps/strings
        {
            const size_t qlen=pg.size()-main_pg_end;
            if(qlen>=MINMEM){
                std::vector<uint8_t> c2(qlen,0), cr2(qlen,0);
                std::string Q(pg,main_pg_end,qlen);
                nm_second =run(Q.data(),qlen,CROSS,c2,false,main_pg_end,nullptr,0,0,MEM_MAXMM);
                std::string R=Q; rc_inplace(R);
                nm_second+=run(R.data(),qlen,CROSS,cr2,true,main_pg_end,nullptr,0,0,MEM_MAXMM);
                for(size_t i=0;i<qlen;++i) if(cr2[i]) c2[qlen-1-i]=1;
                // The second region is matched against the MAIN pg only, so it
                // can never reference itself. Where the main pg is large that is
                // invisible; where it is small nothing is removed. SARS-CoV-2:
                // 14.8 MB second region, 2.4% removed, against 74.3% repeat
                // 32-mers -- the redundancy is there and unreachable by
                // construction.
                //
                // Decodability: the decoder applies references in dst order, so
                // a source must already be reconstructed, i.e. src+len <= dst.
                // SELF_FWD enforces exactly that at match time (s<qp, capL=qp-s)
                // in region-local coordinates, and DSTBASE/SRCBASE shift source
                // and destination by the same main_pg_end, so the invariant
                // survives the lift into pseudogenome coordinates.
                if(getenv("SECOND_SELF")){
                    const int SSMODE=atoi(getenv("SECOND_SELF"));  // 1=fwd only, 2=fwd+rc
                    std::vector<uint8_t> c3(qlen,0), cr3(qlen,0);
                    build_index(Q.data(), qlen);              // index the second region
                    nm_second+=run(Q.data(),qlen,SELF_FWD,c3,false,main_pg_end,
                                   Q.data(),qlen,main_pg_end);
                    if(SSMODE>=2){ std::string R3=Q; rc_inplace(R3);
                      nm_second+=run(R3.data(),qlen,SELF_RC,cr3,true,main_pg_end,
                                     Q.data(),qlen,main_pg_end); }
                    for(size_t i=0;i<qlen;++i) if(cr3[i]) c3[qlen-1-i]=1;
                    for(size_t i=0;i<qlen;++i) if(c3[i]) c2[i]=1;
                    build_index(pg.data(), main_pg_end);      // restore
                }
                for(size_t i=0;i<qlen;++i) if(c2[i]) ++rem_second;
            }
        }
        lap("  [diag] both run() passes done, allrefs built");

        // Trim-on-overlap pass, exactly mirroring SimplePgMatcher.cpp:99-131.
        //
        // RAM FIX -- found at C. elegans scale (LAYER_BY_LAYER_ANALYSIS.md
        // section 3c): the original version built a SEPARATE `cleanRefs`
        // vector while the raw, pre-trim `allrefs` was still fully alive,
        // only freeing it via swap() at the very end -- a real, measured
        // transient double-holding (assembly's peak RSS, 1.27GB, was far
        // above its own post-phase checkpoint of 861MB, localizing the
        // spike to exactly this window). Fixed via in-place compaction:
        // the scan is strictly sequential and every kept element's write
        // position is always <= its read position (nothing is ever
        // reordered or duplicated), so it's safe to write the trimmed
        // result back into `allrefs` itself as we go, then resize() down
        // -- zero second allocation, not a bigger or smarter one.
        std::sort(allrefs.begin(),allrefs.end(),[](const Ref&a,const Ref&b){ return a.dst<b.dst; });
        const size_t rawRefCount=allrefs.size();
        FILE* litf = STR.literal.f;
        {
            uint64_t pos=0; size_t w=0;
            for(size_t r=0;r<allrefs.size();++r){
                Ref m=allrefs[r];   // copy: safe even though allrefs[w] may alias allrefs[r] when w==r
                if(m.dst<pos){
                    const uint64_t overflow=pos-m.dst;
                    if(overflow>=m.len) continue;             // fully swallowed by an earlier match
                    m.len-=(uint32_t)overflow;
                    m.dst+=(uint32_t)overflow;
                    if(!m.is_rc) m.src+=(uint32_t)overflow;   // RC: src stays fixed on front-trim (proven design)
                    // Mismatch offsets are destination-relative from the match's
                    // OLD start: one inside the trimmed-off prefix belongs to
                    // literal now, not to this reference, and every surviving
                    // offset shifts left by the same overflow.
                    uint8_t w2=0;
                    for(uint8_t t=0;t<m.mmcnt;++t){
                        if(m.mmpos[t]<overflow) continue;
                        m.mmpos[w2++]=m.mmpos[t]-(uint32_t)overflow;
                    }
                    m.mmcnt=w2;
                }
                if(m.len<MINMEM) continue;                    // too short to be worth keeping once trimmed
                if(litf) for(uint64_t p=pos;p<m.dst;++p) fputc(pg[p],litf);   // literal gap before this match
                allrefs[w++]=m;
                pos=(uint64_t)m.dst+m.len;
            }
            if(litf) for(uint64_t p=pos;p<pg.size();++p) fputc(pg[p],litf);  // trailing literal
            allrefs.resize(w);
            allrefs.shrink_to_fit();   // actually release the now-unused tail, not just logically shrink
        }
        if(litf){ STR.literal.close(); fprintf(stderr,"[LIT] literal written (%zu clean refs, from %zu raw)\n",allrefs.size(),rawRefCount); }
        lap("  [diag] trim-on-overlap + literal.txt done");
        if(getenv("DBG_FAILMODE"))
            fprintf(stderr,"[FAILMODE] no_kmer=%zu no_seed_hit=%zu maxcand_exceeded=%zu seed_hit_too_short=%zu\n",
                    dbgNoKmer.load(),dbgNoSeedHit.load(),dbgCandExceeded.load(),dbgSeedButShort.load());
    }
    lap("pg MEM matching");

    const size_t lit_main   = main_pg_end-rem_main;
    const size_t lit_second = (pg.size()-main_pg_end)-rem_second;
    fprintf(stderr,"MEM main   : %zu -> %zu (removed %zu, %.1f%%, %zu matches)\n",
            main_pg_end,lit_main,rem_main,rem_main*100.0/(main_pg_end?main_pg_end:1),nm_main);
    fprintf(stderr,"MEM second : %zu -> %zu (removed %zu, %.1f%%, %zu matches)\n",
            pg.size()-main_pg_end,lit_second,rem_second,
            rem_second*100.0/((pg.size()-main_pg_end)?(pg.size()-main_pg_end):1),nm_second);
    fprintf(stderr,"PgRC2 for reference: main 21104500 -> 11463320 (45.8%%), "
                   "lq+n 1506135 -> 1231855, final literal 12695175\n");
    // Reference streams, in the form an archive would actually store: matches in
    // destination order, the start as a gap from the previous match's end (so the
    // values are small and positive), the source as an absolute pg offset, and
    // the length as a delta above MINMEM.
    {
        std::sort(allrefs.begin(),allrefs.end(),
                  [](const Ref&a,const Ref&b){ return a.dst<b.dst; });
        std::vector<uint8_t> srcs;
        // gaps/lens/rc are hoisted (declared before the MEM block) because they
        // MUST reach the archive: refc::encode stores only src.
        ref_gaps.clear(); ref_lens.clear(); ref_rc.clear();
        std::vector<uint8_t>& gaps = ref_gaps;
        std::vector<uint8_t>& lens = ref_lens;
        auto vint=[](std::vector<uint8_t>& o,uint64_t v){
            while(true){ uint8_t b=v&0x7f; v>>=7; o.push_back(b|(v?0x80:0)); if(!v) break; } };
        // Source offsets are the dominant component (182 KB of 300 KB, 24.7 bits
        // each = log2(pg length), i.e. the naive floor with no locality used).
        // Consecutive matches in destination order often point at nearby sources,
        // so a zigzag delta is also emitted and the smaller kept -- the same
        // keep-the-better discipline used for the aux layouts in ARCS.
        std::vector<uint8_t> srcd;
        uint64_t prev=0; int64_t psrc=0;
        for(const Ref& r:allrefs){
            const uint64_t g=(r.dst>=prev)?(r.dst-prev):0;
            vint(gaps,g); vint(srcs,r.src); vint(lens,(uint64_t)(r.len-MINMEM));
            ref_rc.push_back((uint8_t)(r.is_rc?1:0));
            // Repeat DISTANCE, not absolute position. Every self-match points
            // backwards, so dst-src is the LZ77 match distance, and distances
            // are strongly skewed towards small values where absolute offsets
            // into a 23 Mbase pg are near-uniform. (An earlier test used the
            // delta between CONSECUTIVE sources, which is a different quantity
            // and was worse than absolute -- sources have no locality with each
            // other, but each source is close to its own destination.)
            const int64_t d=(int64_t)r.dst-(int64_t)r.src;
            vint(srcd,(uint64_t)((d<<1)^(d>>63)));
            psrc=(int64_t)r.src;
            prev=(uint64_t)r.dst+r.len;
        }
        // mem_srcdelta: diagnostic only, dropped
        // STAGE 100 REWRITE: added is_rc (uint8, 4th field) -- previously
        // absent entirely, meaning a decoder had no way to know an RC-
        // discovered match needed reverse-complementing before copying.
        // allrefs here is the CLEAN, non-overlapping, trim-on-overlap
        // result (see above) -- exactly the runs literal.txt's marks
        // correspond to, one triple per mark, in the same dst order.
        { FILE* t=STR.mem_triples.f;
          ref_mmcnt.clear(); ref_mmpos.clear(); ref_mmref.clear(); ref_mmobs.clear();
          auto vintmm=[](std::vector<uint8_t>& o,uint64_t v){
              while(true){ uint8_t b=v&0x7f; v>>=7; o.push_back(b|(v?0x80:0)); if(!v) break; } };
          size_t totalmm=0;
          auto compb=[](uint8_t c)->uint8_t{ return c=='A'?'T':c=='C'?'G':c=='G'?'C':'A'; };
          for(const Ref& r:allrefs){
              fwrite(&r.dst,4,1,t); fwrite(&r.src,4,1,t); fwrite(&r.len,4,1,t);
              const uint8_t rc=(uint8_t)r.is_rc; fwrite(&rc,1,1,t);
              ref_mmcnt.push_back(r.mmcnt);
              for(uint8_t k=0;k<r.mmcnt;++k){
                  vintmm(ref_mmpos, r.mmpos[k]);
                  // ref: what the plain copy places at dst+offset BEFORE any
                  // override -- the exact value decode will independently
                  // compute at the same point in its own pg-rebuild pass.
                  const uint8_t refb = r.is_rc ? compb((uint8_t)pg[(size_t)r.src+r.len-1-r.mmpos[k]])
                                                : (uint8_t)pg[(size_t)r.src+r.mmpos[k]];
                  ref_mmref.push_back(refb);
                  ref_mmobs.push_back((uint8_t)pg[(size_t)r.dst+r.mmpos[k]]);
              }
              totalmm += r.mmcnt;
          }
          if(totalmm) fprintf(stderr,"[MMTOL] extension mismatches recorded: %zu across %zu refs\n",
                               totalmm, allrefs.size());
          // REAL PIPELINE INTEGRATION (2026-09-03): feed Claim 1's own
          // mem_extmm data (pg-internal cross-haplotype SNV mismatches,
          // computed above, gated by the same heterozygosity detection
          // that gates MEM_MAXMM) directly into Claim 2's VCF, additively.
          // This is the actual "call from what compression already
          // computed" pipeline, not a probe: real archive-computed data,
          // real VCF output, appended after run_variant_call's own
          // contig-bubble calls so nothing already emitted is disturbed.
          if (CAPS_CALL && totalmm) {
              std::string call_vcf = getenv("CALL_VCF") ? getenv("CALL_VCF") : "out.vcf";
              FILE* av = fopen(call_vcf.c_str(), "a");
              if (av) {
                  auto find_cid_mm=[&](uint64_t p)->int64_t{
                      auto it=std::upper_bound(g_contig_spans.begin(),g_contig_spans.end(),
                          std::make_pair(p,(uint64_t)0),
                          [](const std::pair<uint64_t,uint64_t>& a,const std::pair<uint64_t,uint64_t>& b){return a.first<b.first;});
                      if(it==g_contig_spans.begin()) return -1;
                      --it;
                      if(p>=it->second) return -1;
                      return (int64_t)(it-g_contig_spans.begin());
                  };
                  size_t mmi=0, emitted=0;
                  for(const Ref& r : allrefs){
                      for(uint8_t k=0;k<r.mmcnt;++k){
                          uint64_t dstpos = (uint64_t)r.dst + r.mmpos[k];
                          int64_t cid = find_cid_mm(dstpos);
                          if (cid>=0){
                              uint64_t local = dstpos - g_contig_spans[(size_t)cid].first;
                              char refc_ = (char)ref_mmref[mmi+k], obsc_ = (char)ref_mmobs[mmi+k];
                              fprintf(av, "mcontig_%ld\t%zu\t.\t%c\t%c\t.\tPASS\tSVTYPE=SNV;SRC=mem_extmm;MLEN=%u;MMCNT=%u\n",
                                      (long)cid, (size_t)local+1, refc_, obsc_, r.len, r.mmcnt);
                              ++emitted;
                          }
                      }
                      mmi += r.mmcnt;
                  }
                  // DIFFERENTIAL-GAP INDEL CHANNEL (2026-09-03) -- the one
                  // mechanism in this codebase that gets its indel evidence
                  // the way DiscoSNP++ gets its precision: from TWO
                  // independently-real anchors, not from searching for one
                  // coincidental shared k-mer between two separately-built
                  // contigs (the bubble channel's known weakness -- its
                  // false and true positives are statistically identical,
                  // docs/INDEL_PRECISION_ROOT_CAUSE.md §8).
                  // MEM's matcher is substitution-only, so a real indel
                  // FORCES it to terminate a match and re-seed just past the
                  // event. Two adjacent matches on the same haplotype pair
                  // therefore straddle the indel, and the DIFFERENCE between
                  // their src-gap and dst-gap IS the indel length -- the same
                  // way minimap2-class aligners infer indels between
                  // collinear blocks. Both ends are anchored by real extended
                  // alignments, which is precisely the "re-convergence on a
                  // shared node" property our bubble channel can only
                  // approximate by luck.
                  {
                      size_t n_gap = 0;
                      const int MAXGAPINDEL = 12;
                      for (size_t i2 = 1; i2 < allrefs.size(); ++i2) {
                          const Ref& a2 = allrefs[i2-1];
                          const Ref& b2 = allrefs[i2];
                          if (a2.is_rc != b2.is_rc) continue;
                          uint64_t a_dend = (uint64_t)a2.dst + a2.len;
                          uint64_t a_send = (uint64_t)a2.src + a2.len;
                          if (b2.dst < a_dend) continue;              // overlapping, not adjacent
                          if (b2.src < a_send) continue;
                          uint64_t dgap = (uint64_t)b2.dst - a_dend;
                          uint64_t sgap = (uint64_t)b2.src - a_send;
                          if (dgap == sgap) continue;                  // collinear, no indel
                          long long delta = (long long)dgap - (long long)sgap;
                          if (delta > MAXGAPINDEL || delta < -MAXGAPINDEL) continue;
                          if (dgap > 64 || sgap > 64) continue;        // too far apart to trust
                          int64_t cid2 = find_cid_mm(a_dend);
                          if (cid2 < 0) continue;
                          uint64_t local2 = a_dend - g_contig_spans[(size_t)cid2].first;
                          if (local2 == 0) continue;
                          char anchb = pg[(size_t)a_dend - 1];
                          if (delta > 0) {   // insertion in dst relative to src
                              std::string ins2 = pg.substr((size_t)a_dend, (size_t)delta);
                              fprintf(av, "mcontig_%ld\t%zu\t.\t%c\t%c%s\t.\tPASS\tSVTYPE=INDEL;SRC=mem_gap;LEN=%lld;DG=%llu;SG=%llu\n",
                                      (long)cid2, (size_t)local2, anchb, anchb, ins2.c_str(),
                                      (long long)delta, (unsigned long long)dgap, (unsigned long long)sgap);
                          } else {           // deletion in dst relative to src
                              long long dl = -delta;
                              if ((size_t)a_send + (size_t)dl > pg.size()) continue;
                              std::string del2 = pg.substr((size_t)a_send, (size_t)dl);
                              fprintf(av, "mcontig_%ld\t%zu\t.\t%c%s\t%c\t.\tPASS\tSVTYPE=INDEL;SRC=mem_gap;LEN=%lld;DG=%llu;SG=%llu\n",
                                      (long)cid2, (size_t)local2, anchb, del2.c_str(), anchb,
                                      (long long)dl, (unsigned long long)dgap, (unsigned long long)sgap);
                          }
                          ++n_gap;
                      }
                      fprintf(stderr,"[CAPS-CALL] mem_gap INDEL channel: %zu records\n", n_gap);
                  }
                  fclose(av);
                  fprintf(stderr,"[CAPS-CALL] mem_extmm channel: %zu SNV records appended to %s\n", emitted, call_vcf.c_str());
              }
          }
          if(getenv("DBG_STREAMS")){
              fprintf(stderr,"[DBGSTREAM] ref_mmcnt.size=%zu ref_mmpos.size=%zu ref_mmref.size=%zu ref_mmobs.size=%zu\n",
                      ref_mmcnt.size(), ref_mmpos.size(), ref_mmref.size(), ref_mmobs.size());
              const uint64_t want = getenv("DBG_DST")?strtoull(getenv("DBG_DST"),nullptr,10):UINT64_MAX;
              size_t shown=0, mmiE=0;
              for(size_t i=0;i<ref_mmcnt.size();++i){
                  const uint8_t c=ref_mmcnt[i];
                  const uint32_t rdst = allrefs[i].dst;
                  if(!c){ continue; }
                  const bool near = (want!=UINT64_MAX) && rdst+200>=want && rdst<=want+50;
                  if(near || (want==UINT64_MAX && shown<8)){
                      fprintf(stderr,"[DBGSTREAM] i=%zu dst=%u src=%u len=%u cnt=%u  ",
                              i,allrefs[i].dst,allrefs[i].src,allrefs[i].len,c);
                      for(uint8_t k=0;k<c;++k){
                          fprintf(stderr,"pos=%u ref=%c obs=%c  ",allrefs[i].mmpos[k],ref_mmref[mmiE+k],ref_mmobs[mmiE+k]);
                      }
                      fprintf(stderr,"\n");
                      ++shown;
                  }
                  mmiE+=c;
              }
          }
          if(const char* dbg=getenv("DBG_POS")){
              const uint64_t target=strtoull(dbg,nullptr,10);
              for(const Ref& r:allrefs){
                  if(target>=r.dst && target<(uint64_t)r.dst+r.len){
                      fprintf(stderr,"[DBGPOS] byte %llu covered by ref dst=%u src=%u len=%u is_rc=%d mmcnt=%u",
                              (unsigned long long)target,r.dst,r.src,r.len,(int)r.is_rc,r.mmcnt);
                      for(uint8_t t=0;t<r.mmcnt;++t) fprintf(stderr," mm[%u]=%u",t,r.mmpos[t]);
                      fprintf(stderr,"\n");
                      fprintf(stderr,"[DBGPOS] pg[dst..dst+len) = %.*s\n", (int)std::min((uint32_t)60,r.len), &pg[r.dst]);
                      fprintf(stderr,"[DBGPOS] pg[src..src+len) = %.*s\n", (int)std::min((uint32_t)60,r.len), &pg[r.src]);
                  }
              }
          }
          STR.mem_triples.close(); }
        // gaps and lens are NO LONGER dropped -- see the archive jobs below.
        fprintf(stderr,"[REF] matches=%zu  raw gaps=%zu srcs=%zu lens=%zu  total_raw=%zu\n",
                allrefs.size(),gaps.size(),srcs.size(),lens.size(),
                gaps.size()+srcs.size()+lens.size());
        fprintf(stderr,"[REF] PgRC2 pays 177,180 B coded for the same thing\n");
    }
    // DUMP_PG: the pseudogenome BEFORE MEM removal. Needed to test whether a
    // context model with a match model can capture the long repeats implicitly,
    // which would remove the reference streams entirely rather than shrink them.
    if(getenv("DUMP_PG")){
        FILE* f=fopen("pg_full.txt","wb"); fwrite(pg.data(),1,pg.size(),f); fclose(f);
        fprintf(stderr,"[PG] pg_full.txt written: %zu bases\n",pg.size());
    }
    printf("PG_LEN %zu\n",pg.size());
    printf("PG_LITERAL %zu\n",lit_main+lit_second);

    // ================= COMPRESS-AND-RELEASE, SINGLE ARCHIVE =================
    // Each stream is coded in-process and its raw buffer freed IMMEDIATELY
    // afterwards, so peak RAM tracks the largest live stream instead of the
    // sum of all of them -- PgRC2's disposeReadsList/clear discipline
    // (pgrc-encoder.cpp:157-226), which our multi-process design achieved only
    // by paying 361 MB of disk I/O per 456 MB of input.
    {
        const char* apath = getenv("ARCHIVE") ? getenv("ARCHIVE") : "out.arcs2";
        Archive ar(apath, (uint64_t)pg.size(), (uint64_t)main_pg_end, (uint32_t)MINMEM);
        auto __t0 = std::chrono::steady_clock::now();
        auto CODED = [&](const char* nm){
            auto n = std::chrono::steady_clock::now();
            fprintf(stderr,"  [codetime] %-14s %6.2f s\n", nm,
                    std::chrono::duration<double>(n-__t0).count());
            __t0 = n;
        };
        const size_t rss_before = rss_mb();
        // VERIFY_DUMP: also write the RAW (pre-coding) streams to files purely
        // so decode_105.py can prove the in-memory pipeline produces exactly
        // the bytes the file-based one did. Off by default -- the shipping
        // path never touches disk.
        const bool VDUMP = getenv("VERIFY_DUMP") && atoi(getenv("VERIFY_DUMP"));
        auto vdump=[&](const char* fn, const MemStream& m){
            if(!VDUMP) return;
            // open_memstream only publishes buf/len on flush. Without this the
            // still-buffered streams (mm_*, orig2uid) dumped as 0 bytes and the
            // decoder silently verified against empty input.
            if(m.f) fflush(m.f);
            FILE* g=fopen(fn,"wb"); if(!g) return;
            if(m.len) fwrite(m.buf,1,m.len,g);
            fclose(g);
        };
        vdump("literal.txt",STR.literal);       vdump("mem_triples.bin",STR.mem_triples);
        vdump("pos_abs.bin",STR.pos_abs);       vdump("pos_strand.bin",STR.pos_strand);
        vdump("mm_ref.bin",STR.mm_ref);         vdump("mm_obs.bin",STR.mm_obs);
        vdump("mm_pos.bin",STR.mm_pos);         vdump("mm_count_per_read.bin",STR.mm_count);
        vdump("n_pos.bin",STR.n_pos);           vdump("n_indices.bin",STR.n_indices);
        vdump("n_cnt.bin",STR.n_cnt);           vdump("read_lengths.bin",STR.read_lengths);
        vdump("orig2uid.bin",STR.orig2uid);
        // The decoder needs pg length and the main/second split point, which
        // were previously only recoverable by hand-reading the MEM main/second
        // lines out of stderr. Writing them beside the dumps makes the lossless
        // check reproducible instead of a manual step.
        if(VDUMP){
            if(FILE* pf=fopen("pg_params.txt","w")){
                fprintf(pf,"%zu %zu\n", pg.size(), main_pg_end);
                fclose(pf);
            }
        }

        // ---------------- PARALLEL STREAM CODING ----------------
        // The streams are independent, so coding them concurrently is
        // embarrassingly parallel and changes no bytes. Measured before this:
        // coding was 4.88 s of a 13.9 s run at ~100% CPU (fully serial) while
        // PgRC2 ran the whole compressor at 515%. RSS is flat across the
        // coding phase (195 MB before and after), i.e. every raw buffer is
        // already live before coding begins, so running them together does not
        // raise peak memory.
        //
        // Results are collected into a fixed-order slot table and written to
        // the archive in that order, so the output is byte-for-byte identical
        // and independent of thread scheduling.
        struct Job { const char* name; std::function<std::vector<uint8_t>()> fn; };
        std::vector<Job> jobs;
        std::vector<std::vector<uint8_t>> results;

        // mm_count must be decoded first: its values drive the mm_pos
        // bucketing. Do that one inline, then queue everything else.
        std::vector<uint16_t> mmcounts;
        std::vector<uint8_t> mmcnt_flags, mmcnt_vals, mmcnt_flat;
        bool mmcnt_is_split = false;
        {
            auto v = STR.mm_count.bytes();
            mmcounts.resize(v.size()/2);
            if(!mmcounts.empty()) memcpy(mmcounts.data(), v.data(), mmcounts.size()*2);
            STR.mm_count.release();
            mmcnt_flat = std::move(v);
        }

        // literal: symbol-pack now (cheap, frees the big raw buffer), code later
        std::vector<uint8_t> litsym;
        {
            litsym.reserve(STR.literal.len);
            for(size_t i=0;i<STR.literal.len;++i){
                switch((uint8_t)STR.literal.buf[i]){
                    case 'A': litsym.push_back(0); break; case 'C': litsym.push_back(1); break;
                    case 'G': litsym.push_back(2); break; case 'T': litsym.push_back(3); break; }
            }
            STR.literal.release();
        }
        lap("  [t] literal->symbols");
        auto v_tri  = STR.mem_triples.bytes(); STR.mem_triples.release();
        auto v_pos  = STR.pos_abs.bytes();     STR.pos_abs.release();
        auto v_str  = STR.pos_strand.bytes();  STR.pos_strand.release();
        auto v_mr   = STR.mm_ref.bytes();      STR.mm_ref.release();
        auto v_mo   = STR.mm_obs.bytes();      STR.mm_obs.release();
        auto v_mp   = STR.mm_pos.bytes();      STR.mm_pos.release();
        auto v_np   = STR.n_pos.bytes();       STR.n_pos.release();
        auto v_ni   = STR.n_indices.bytes();   STR.n_indices.release();
        auto v_nc   = STR.n_cnt.bytes();       STR.n_cnt.release();
        auto v_rl   = STR.read_lengths.bytes();STR.read_lengths.release();
        auto v_o2u  = STR.orig2uid.bytes();    STR.orig2uid.release();

        const uint64_t PGLEN_ = pg.size(), MAINEND_ = main_pg_end;
        // seqpar is itself threaded; give it fewer threads so it does not
        // oversubscribe against the other jobs running alongside it.
        unsigned HW = std::thread::hardware_concurrency(); if(!HW) HW = 12;

        // seqpar splits the literal into NCHUNKS coded independently, so the chunk
        // count is a COMPRESSION parameter, not just a threading one -- fewer
        // chunks means more context per chunk. Measured on E. coli literal:
        // 12 chunks 1,809,073 B; 4 chunks 1,804,223; 1 chunk 1,801,710. Four is
        // the knee: within 2.5 KB of the single-chunk optimum while still
        // parallel, and it also gained 15,862 B on P. aeruginosa.
        const unsigned SEQT = getenv("SEQT") ? (unsigned)atoi(getenv("SEQT")) : 4;

        // Names (Phase 2). Computed once here, emitted as three streams: the
        // block-coded body (already entropy coded, stored as-is), the global
        // dictionary and the block index. The latter two are RAW out of the
        // coder specifically so they go through best_encode() like every
        // other stream instead of carrying an ad-hoc serialization -- and so
        // their cost is real archive bytes, not an estimate. The block index
        // is what lets the decoder find block boundaries from the archive
        // alone; without it a decoder would need the original file, which is
        // exactly the shortcut that left mem_triples undecodable.
        static nmc::Encoded NM;
        if(CAPS_NAMES && !g_input_path.empty()){
            // Reuse what the parent computed before the fork; only recompute
            // if the hoist did not run (single-candidate / GSEARCH paths).
            NM = g_NM_done ? g_NM : nmc::encode_from_fastq(g_input_path.c_str());
            fprintf(stderr,"  [names] %llu names in %llu blocks: body %zu B, dict %zu B raw, index %zu B raw\n",
                    (unsigned long long)NM.n_names,(unsigned long long)NM.n_blocks,
                    NM.body.size(), NM.dict.size(), NM.index.size());
            jobs.push_back({"names_body",  [&]{ return NM.body; }});
            jobs.push_back({"names_dict",  [&]{ return best_encode(NM.dict.data(),  NM.dict.size());  }});
            jobs.push_back({"names_index", [&]{ return best_encode(NM.index.data(), NM.index.size()); }});
        }

        // Quality (Phase 3). Two streams, same shape as names: the fqzcomp
        // block payloads concatenated (already entropy coded, stored as-is)
        // and the block index RAW so it goes through best_encode() like every
        // other stream and its cost lands in the archive rather than in an
        // estimate. The index carries per block (byte length, read count,
        // symbol offset) -- all three are needed to invert the column from the
        // archive alone.
        static qlc::Encoded QL;
        if(CAPS_QUAL && !g_input_path.empty()){
            QL = g_QL_done ? g_QL : qlc::encode_from_fastq(g_input_path.c_str());
            fprintf(stderr,"  [qual] %llu reads / %llu quality bytes in %llu blocks: body %zu B, index %zu B raw\n",
                    (unsigned long long)QL.n_reads,(unsigned long long)QL.n_qbytes,
                    (unsigned long long)QL.n_blocks, QL.body.size(), QL.index.size());
            jobs.push_back({"qual_body",  [&]{ return QL.body; }});
            jobs.push_back({"qual_index", [&]{ return best_encode(QL.index.data(), QL.index.size()); }});
        }

        lap("  [t] .bytes() copies + pre-jobs");
        jobs.push_back({"literal",     [&]{ return seq_encode_mem(litsym, SEQT, SEQT); }});
        jobs.push_back({"mem_triples", [&]{ return refc::encode(v_tri, PGLEN_, MAINEND_); }});
        // Emitted only when a second-region self pass actually produced
        // references; absent otherwise, so every existing archive is unchanged.
        // Pushed only when self references exist: an empty stream still costs
        // its name and length in the container (17 B), and every existing
        // archive must stay byte-identical.
        { auto sf = refc::encode_self(v_tri, MAINEND_);
          if(!sf.empty())
              jobs.push_back({"mem_self", [&,sf]{ return best_encode(sf.data(), sf.size()); }}); }
        // Emitted only when any extension mismatch was actually recorded, same
        // discipline as mem_self: every archive without one is unaffected, byte
        // for byte -- which is every one of the 7 locked datasets, since
        // MEM_MAXMM is 0 there by construction (leftover_frac ~ 0).
        if(!ref_mmobs.empty()){
            jobs.push_back({"mem_extmm_cnt", [&]{ return best_encode(ref_mmcnt.data(), ref_mmcnt.size()); }});
            jobs.push_back({"mem_extmm_pos", [&]{ return best_encode(ref_mmpos.data(), ref_mmpos.size()); }});
            jobs.push_back({"mem_extmm_obs", [&]{ return mmc::encode(ref_mmref, ref_mmobs); }});
        }
        jobs.push_back({"mem_dstgap",  [&]{ return best_encode(ref_gaps.data(), ref_gaps.size()); }});
        jobs.push_back({"mem_len",     [&]{ return best_encode(ref_lens.data(), ref_lens.size()); }});
        jobs.push_back({"mem_rc",      [&]{ return best_encode(ref_rc.data(),   ref_rc.size());   }});
        jobs.push_back({"pos_abs",     [&]{ return best_encode_chunked(v_pos.data(), v_pos.size(), true); }});
        jobs.push_back({"pos_strand",  [&]{ return best_encode(v_str.data(), v_str.size(), false); }});
        jobs.push_back({"mm_sym",      [&]{ return mmc::encode(v_mr, v_mo); }});
        // mm_pos: the encoder picks whichever of the flat and bucketed forms is
        // smaller, but emitted no flag saying which -- the decoder had no way to
        // know. One leading byte fixes it: 0 = flat, 1 = bucketed.
        jobs.push_back({"mm_pos",      [&]{
            auto flat = best_encode(v_mp.data(), v_mp.size());
            // The bucketed form indexes one BYTE per mismatch; with varint
            // positions (Lmax>256) that mapping no longer holds, so the flat
            // form is forced there. Lmax<=256 is unaffected.
            auto buck = (Lmax<=256) ? mmpos_encode_buckets(v_mp, mmcounts)
                                    : std::vector<uint8_t>();
            const bool useB = (!buck.empty() && buck.size()+1 < flat.size()+1);
            std::vector<uint8_t> o; o.reserve((useB?buck.size():flat.size())+1);
            o.push_back(useB?1:0);
            const auto& src = useB?buck:flat;
            o.insert(o.end(), src.begin(), src.end());
            return o; }});
        jobs.push_back({"mm_cnt",      [&]{
            std::vector<uint8_t> cf, cv;
            if(mmcnt_split(mmcounts, cf, cv) && mmcnt_join(cf,cv,mmcounts.size())==mmcounts){
                auto a = best_encode(cf.data(), cf.size());
                auto b = best_encode(cv.data(), cv.size());
                auto flat = best_encode(mmcnt_flat.data(), mmcnt_flat.size());
                if(getenv("LOG_STREAMS"))
                    fprintf(stderr,"  [mmcnt] split=%zu (flags %zu + vals %zu)  flat=%zu  -> %s\n",
                            a.size()+b.size(), a.size(), b.size(), flat.size(),
                            (a.size()+b.size()<flat.size())?"SPLIT":"flat");
                if(a.size()+b.size() < flat.size()){
                    mmcnt_is_split = true;
                    mmcnt_flags = std::move(a); mmcnt_vals = std::move(b);
                    return std::vector<uint8_t>{};
                }
                return flat;
            }
            return best_encode(mmcnt_flat.data(), mmcnt_flat.size()); }});
        jobs.push_back({"n_pos",       [&]{ return best_encode(v_np.data(), v_np.size()); }});
        jobs.push_back({"n_indices",   [&]{ return best_encode(v_ni.data(), v_ni.size(), true); }});
        jobs.push_back({"n_cnt",       [&]{ return best_encode(v_nc.data(), v_nc.size()); }});
        // read_lengths is uint16 per original read; on constant-length data that
        // is one value repeated millions of times (E. coli: 3,106,518 B of "150",
        // entropy 0). Detected exactly, so variable-length input is unaffected.
        jobs.push_back({"read_lengths",[&]{ return const_or_encode(v_rl.data(), v_rl.size(), 2); }});
        // B3: orig2uid mixes two different signals in one stream. Measured on
        // E. coli: 79.55% of its deltas are ZERO (the reads that are not
        // duplicates) and 20.45% are sparse alias payloads over 273,172
        // distinct values. Its order-0 floor is 855,909 B against the 6,213,036
        // we hand the coder, and no single model fits both halves -- which is
        // precisely why this stream needed a 7-way probe and became the coding
        // pool's floor at 4.53 s, longer than `literal` which is BIGGER but
        // understood (0.38 s through the DNA coder).
        //
        // PgRC2 splits exactly this shape rather than searching for a coder
        // that copes; their log shows "Mismatches counts (zero flags)" and
        // "(non-zero values)" as separate streams. We already do it for mm_cnt.
        // Splitting hands the coder the structure instead of making it discover
        // it, and drops the data through the probe from 6.21 MB to ~1.46 MB.
        {
            const size_t cnt = v_o2u.size()/4;
            std::vector<uint8_t> ofl((cnt+7)/8, 0), ova;
            ova.reserve(v_o2u.size()/4);
            for(size_t i=0;i<cnt;++i){
                uint32_t x; memcpy(&x, v_o2u.data()+i*4, 4);
                if(x){
                    ofl[i>>3] |= (uint8_t)(0x80u >> (i&7));
                    ova.insert(ova.end(), v_o2u.data()+i*4, v_o2u.data()+i*4+4);
                }
            }
            jobs.push_back({"orig2uid_flags",[&,ofl]{
                return best_encode(ofl.data(), ofl.size()); }});
            jobs.push_back({"orig2uid_vals", [&,ova]{
                return best_encode(ova.data(), ova.size(), true); }});
            fprintf(stderr,"  [b3] orig2uid %zu B -> flags %zu B + vals %zu B (%zu nonzero of %zu)\n",
                    v_o2u.size(), ofl.size(), ova.size(), ova.size()/4, cnt);
        }

        results.resize(jobs.size());
        {
            std::atomic<size_t> next{0};
            unsigned NT = HW; if(NT > jobs.size()) NT = (unsigned)jobs.size();
            std::vector<double> jsec(jobs.size(),0.0);
            lap("  [t] job list built");
            auto POOL0 = std::chrono::steady_clock::now();
            phase("pre-coding");
            std::vector<std::thread> pool;
            for(unsigned t=0;t<NT;++t) pool.emplace_back([&]{
                for(;;){ size_t i = next++; if(i >= jobs.size()) break;
                         auto a=std::chrono::steady_clock::now();
                         results[i] = jobs[i].fn();
                         jsec[i]=std::chrono::duration<double>(
                                    std::chrono::steady_clock::now()-a).count(); }
            });
            for(auto& th : pool) th.join();
            double pool_wall=std::chrono::duration<double>(
                                std::chrono::steady_clock::now()-POOL0).count();
            double sum=0,mx=0; const char* mxn="";
            for(size_t i=0;i<jobs.size();++i){ sum+=jsec[i];
                if(jsec[i]>mx){ mx=jsec[i]; mxn=jobs[i].name; } }
            std::vector<size_t> ord(jobs.size());
            for(size_t i=0;i<ord.size();++i) ord[i]=i;
            std::sort(ord.begin(),ord.end(),
                      [&](size_t a,size_t b){ return jsec[a]>jsec[b]; });
            for(size_t k=0;k<ord.size();++k)
                fprintf(stderr,"  [job] %-14s %6.2f s\n",
                        jobs[ord[k]].name, jsec[ord[k]]);
            // wall cannot fall below the longest single job; that ratio is the
            // ceiling on any further threading of this phase.
            fprintf(stderr,"  [pool] wall %.2f s  cpu-sum %.2f s  longest=%s %.2f s"
                           "  speedup %.2fx  amdahl-floor %.2f s\n",
                    pool_wall,sum,mxn,mx,sum/(pool_wall>0?pool_wall:1),mx);
        }
        // Write in fixed job order -- deterministic regardless of scheduling.
        for(size_t i=0;i<jobs.size();++i){
            if(std::string(jobs[i].name) == "mm_cnt" && mmcnt_is_split){
                ar.put("mm_cnt_flags", mmcnt_flags);
                ar.put("mm_cnt_vals",  mmcnt_vals);
            } else {
                if(getenv("LOG_STREAMS"))
                    fprintf(stderr,"  [stream] %-18s coded=%zu B\n",
                            jobs[i].name, results[i].size());
                ar.put(jobs[i].name, results[i]);
            }
        }

        lap("  [t] pool + archive put");
        ar.finish();
        lap("stream coding");
        phase("CODING");
        fprintf(stderr,"[archive] rss before coding %zu MB, after %zu MB, peak %zu MB\n",
                rss_before,rss_mb(),hwm_mb());
    }
    return 0;
}
