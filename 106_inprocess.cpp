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
#include <omp.h>
#include <unistd.h>
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
struct Archive {
    FILE* out;
    size_t total=0;
    std::vector<std::pair<const char*,size_t>> parts;
    explicit Archive(const char* path){ out=fopen(path,"wb"); }
    void put(const char* name, const std::vector<uint8_t>& coded){
        // (timing is recorded by the caller via tput)
        uint64_t n=coded.size();
        fwrite(&n,8,1,out);
        if(n) fwrite(coded.data(),1,n,out);
        total += 8 + n;
        parts.emplace_back(name,n);
    }
    void finish(){
        if(out) fclose(out);
        out=nullptr;
        for(auto& p:parts) fprintf(stderr,"  [archive] %-14s %10zu B\n",p.first,p.second);
        fprintf(stderr,"ARCHIVE_TOTAL=%zu\n",total);
        printf("ARCHIVE_TOTAL=%zu\n",total);
    }
};

int main(int argc,char** argv){
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
    const uint32_t SW = argc>4?(uint32_t)atoi(argv[4]):32;
    const uint64_t SWMASK = (SW>=32)?~0ULL:((1ULL<<(2*SW))-1);
    auto packSW=[&](const char* p,uint64_t& out)->bool{
        uint64_t k=0;
        for(uint32_t i=0;i<SW;++i){ int v=b2(p[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
        out=k; return true; };   // >= SEED
    auto T0=std::chrono::steady_clock::now();
    auto rss_mb=[]()->size_t{ FILE* f=fopen("/proc/self/statm","r"); if(!f) return 0;
        long pg_=0,res=0; if(fscanf(f,"%ld %ld",&pg_,&res)!=2) res=0; fclose(f);
        return (size_t)((double)res*(double)sysconf(_SC_PAGESIZE)/1048576.0); };
    auto lap=[&](const char* w){ auto t=std::chrono::steady_clock::now();
        fprintf(stderr,"  %-26s %6.2f s   rss=%zu MB\n",w,
                std::chrono::duration<double>(t-T0).count(),rss_mb()); T0=t; };

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
        std::unordered_map<uint64_t,std::vector<uint32_t>> seen; seen.reserve(1u<<21);
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
            auto& bk=seen[fnv(b.data(),(uint32_t)b.size())];
            for(uint32_t id:bk){
                if(rlen[id]!=(uint16_t)b.size()) continue;
                const uint64_t* w=&rpk[woff[id]];
                const uint32_t nw=((uint32_t)rlen[id]+31)/32;
                for(uint32_t t=0;t<nw;++t){
                    for(uint32_t j=0;j<32;++j){
                        const uint32_t idx=t*32+j;
                        ubuf[idx]="ACGT"[(w[t]>>(2*(31-j)))&3ULL];
                    }
                }
                if(memcmp(ubuf,b.data(),b.size())==0){ dup=true; orig2uid.push_back(id); break; }
            }
            if(!dup){
                orig2uid.push_back((uint32_t)rlen.size());
                bk.push_back((uint32_t)rlen.size());
                packstr(b.data(),(uint32_t)b.size(),tmpw);
                woff.push_back(rpk.size());
                rpk.insert(rpk.end(),tmpw.begin(),tmpw.end());
                rlen.push_back((uint16_t)b.size());
            }
        }
        rpk.shrink_to_fit(); woff.shrink_to_fit(); rlen.shrink_to_fit();
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
    auto buildPref=[&](bool useAdmit){
        pent.clear(); pent.reserve(n);
        for(uint32_t i=0;i<n;++i){
            if(useAdmit && !admit[i]) continue;
            if(rlen[i]<SW) continue;
            pent.push_back(((uint64_t)(uint32_t)rseed(i,0,SW)<<32)|(uint64_t)i);
        }
        std::sort(pent.begin(),pent.end());
        size_t sz=1; while(sz<pent.size()*2+1) sz<<=1;
        pmask=sz-1; ptab.assign(sz,UINT32_MAX);
        for(size_t i=0;i<pent.size();++i){
            const uint32_t k=(uint32_t)(pent[i]>>32);
            if(i && (uint32_t)(pent[i-1]>>32)==k) continue;
            size_t h=pmix(k)&pmask;
            while(ptab[h]!=UINT32_MAX) h=(h+1)&pmask;
            ptab[h]=(uint32_t)i;
        }
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
    std::vector<uint32_t> nxt(n,NONE),prv(n,NONE),ovl(n,0),ch_h(n),ch_t(n);
    std::vector<uint64_t> seed(n,0);
    std::vector<uint8_t>  ok(n,0);
    std::vector<uint32_t> tails;
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
    auto sweep=[&](){
        links=0; probes=0;
        std::fill(nxt.begin(),nxt.end(),NONE); std::fill(prv.begin(),prv.end(),NONE);
        std::fill(ovl.begin(),ovl.end(),0);    std::fill(ok.begin(),ok.end(),0);
        for(uint32_t i=0;i<n;++i){ ch_h[i]=i; ch_t[i]=i; }
        tails.clear(); tails.reserve(n);
        for(uint32_t i=0;i<n;++i) if(admit[i]&&rlen[i]>=sweep_minov) tails.push_back(i);
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
    for(uint32_t L=Lmax-1; L>=sweep_minov && L>=SW; --L){
        // compact to tails still open and still long enough at this L
        size_t w=0;
        for(uint32_t a:tails){
            if(nxt[a]!=NONE) continue;                       // already extended
            if(rlen[a]<L) continue;
            const uint32_t off=(uint32_t)rlen[a]-L;          // grows by 1 per length
            if(off+SW>rlen[a]) continue;
            tails[w++]=a;                                    // still open: keep
        }
        tails.resize(w);
        if(tails.empty()) break;

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
        if(cand.size()<w*CCAP) cand.resize(w*CCAP);
        ccnt.assign(w,0);
        {
            const unsigned T=(w<4096)?1u:NT;
            size_t pr_total=0;
            #pragma omp parallel for schedule(static) num_threads(T) reduction(+:pr_total)
            for(long long ii=0;ii<(long long)w;++ii){
                const size_t i=(size_t)ii;
                const uint32_t a=tails[i];
                const uint32_t off=(uint32_t)rlen[a]-L;
                // packed reads are always ACGT, so the seed always exists
                seed[a]=rseed(a,off,SW); ok[a]=1;
                ++pr_total;
                const uint32_t pix=pfind(seed[a]);
                if(pix==UINT32_MAX) continue;
                const uint32_t pk=(uint32_t)seed[a];
                uint8_t c=0;
                for(uint32_t q=pix;q<pent.size()&&(uint32_t)(pent[q]>>32)==pk;++q){
                    const uint32_t b=(uint32_t)(pent[q]&0xFFFFFFFFULL);
                    if(b==a) continue;
                    if(!admit[b]) continue;                  // excluded reads are leftovers,
                                                             // never chain members -- without
                                                             // this they get emitted twice
                    if(rlen[b]<L) continue;
                    if(!rcmp(a,off,b,L)) continue;
                    cand[i*CCAP+c]=b;
                    if(++c==CCAP) break;
                }
                ccnt[i]=c;
            }
            probes+=pr_total;
        }

        // ── serial: same order, same first survivor ───────────────────────
        // A retained list that filled to CCAP may have been truncated, so if
        // every entry in a full list turns out to be taken, that tail is
        // re-probed serially. Without this the cap silently drops the 9th-and-
        // later candidate and the run loses links -- measured at CCAP=8, 14 of
        // 673,334 links and 1,288 bytes of literal. Rare enough to cost nothing,
        // and it makes the result identical to the fully serial sweep rather
        // than merely close to it.
        for(size_t i=0;i<w;++i){
            const uint32_t a=tails[i];
            bool done=false;
            for(uint8_t c=0;c<ccnt[i];++c){
                const uint32_t b=cand[i*CCAP+c];
                if(prv[b]!=NONE||ch_h[a]==b) continue;       // taken, or would cycle
                nxt[a]=b; prv[b]=a; ovl[a]=L;
                uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h; ++links;
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

    // ── round 2: build over the well-overlapped reads only ──────────────────
    sweep_minov=MINOV;                // builder: caller's floor
    sweep();
    fprintf(stderr,"round2: probes=%zu links=%zu\n",probes,links);
    lap("round 2 (assembly)");

    // ── emit chains; singletons held back for mapping ────────────────────────
    std::string pg; pg.reserve((size_t)n*40);
    ppos.assign(n,UINT64_MAX);
    prc.assign(n,0);
    std::vector<uint32_t> leftovers;
    uint32_t multi=0;
    for(uint32_t i=0;i<n;++i){
        if(!admit[i]){ leftovers.push_back(i); continue; }   // excluded in round 1
        if(contained[i]) continue;                 // stage 42: its container carries it
        if(prv[i]!=NONE) continue;
        if(nxt[i]==NONE){ leftovers.push_back(i); continue; }
        ++multi;
        uint32_t cur=i; ppos[cur]=pg.size(); rappend(pg,cur,0);
        while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                               ppos[cur]=pg.size()-o; rappend(pg,cur,o); }
    }
    fprintf(stderr,"links=%zu chains(multi)=%u leftovers=%zu pg after chains=%zu\n",
            links,multi,leftovers.size(),pg.size());
    lap("emit chains");

    // The prefix index is not consulted anywhere in the mapping stage -- that
    // stage builds and queries its own seed index -- but it was staying live
    // across it, so the two large indexes coexisted at exactly the moment peak
    // RSS is set. It is next needed by the survivor sweep, which runs over ~13k
    // reads rather than 851k, so it is cheaper to drop it here and rebuild a
    // small one there than to carry the full one through.
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
    readMM.assign(n,255);
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
        std::vector<uint32_t> mtab(msize,UINT32_MAX);
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
            std::vector<std::thread> th; th.reserve(T);
            const size_t chunk=(text.size()+T-1)/T;
            for(unsigned t=0;t<T;++t){
                const size_t lo=(size_t)t*chunk;
                size_t hi=std::min(text.size(),lo+chunk);
                if(lo>=hi) break;
                th.emplace_back([&,t,lo,hi]{
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
                            if(mm<=lim){ hit[t].push_back({rid,(uint32_t)st}); hmm[t].push_back(mm); }
                        }
                    }
                });
            }
            for(auto& x:th) x.join();
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

    // Survivors get their own pseudogenome rather than being appended whole.
    // They are the reads that neither tiled well enough to assemble in round 2
    // nor matched anything already built -- but they still overlap EACH OTHER,
    // so a second sweep over just this set recovers most of their length.
    // Appending them raw costs a full read apiece (measured: 43,075 x 150 =
    // 6.46 MB); PgRC2 assembles its 11,298 survivors into 1,504,035 bytes.
    size_t appended=0, second_pg=0;
    const size_t main_pg_end = pg.size();
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
            uint32_t cur=i; ppos[cur]=pg.size(); rappend(pg,cur,0);
            while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                                   ppos[cur]=pg.size()-o; rappend(pg,cur,o); }
        }
        second_pg=pg.size()-before;
        fprintf(stderr,"second pg: %zu reads -> %zu B (raw would be %zu B)\n",
                appended,second_pg,appended*(size_t)Lmax);
    }
    fprintf(stderr,"leftovers=%zu mapped=%zu appended=%zu\n",
            leftovers.size(),n_matched,appended);
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
        for(uint32_t i=0;i<n;++i){
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
                if(MMDELTA){ fputc((uint8_t)(j-prevj),fp); prevj=j; }
                else        fputc((uint8_t)(j>255?255:j),fp);
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
        fprintf(stderr,"[POS] PgRC2 pays 683,370 B coded for its reads-list offsets\n");
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
    struct Ref { uint32_t dst, src, len; bool is_rc; };
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
        fprintf(stderr,"[MEM] MINMEM=%zu seed=%zu step=%zu\n",MINMEM,MEMSEED,STEP);
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
        std::vector<uint64_t> skey; std::vector<uint32_t> spos;
        {
            std::vector<std::pair<uint64_t,uint32_t>> tmp;
            tmp.reserve(main_pg_end/STEP+16);
            for(size_t p=0;p+MEMSEED<=main_pg_end;p+=STEP){
                uint64_t k; if(packM(pg.data()+p,k)) tmp.push_back({k,(uint32_t)p});
            }
            std::sort(tmp.begin(),tmp.end());
            skey.resize(tmp.size()); spos.resize(tmp.size());
            for(size_t i=0;i<tmp.size();++i){ skey[i]=tmp[i].first; spos[i]=tmp[i].second; }
        }
        size_t tsize=1; while(tsize < skey.size()*2+1) tsize<<=1;
        const uint64_t TMASK=tsize-1;
        std::vector<uint32_t> htab(tsize,UINT32_MAX);
        auto hmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
        for(size_t i=0;i<skey.size();++i){
            if(i && skey[i]==skey[i-1]) continue;       // first occurrence only
            size_t h=hmix(skey[i])&TMASK;
            while(htab[h]!=UINT32_MAX) h=(h+1)&TMASK;
            htab[h]=(uint32_t)i;
        }
        auto lookup=[&](uint64_t k)->uint32_t{
            size_t h=hmix(k)&TMASK;
            while(htab[h]!=UINT32_MAX){ if(skey[htab[h]]==k) return htab[h]; h=(h+1)&TMASK; }
            return UINT32_MAX;
        };
        lap("  [diag] seed index built");

        enum Mode { CROSS, SELF_FWD, SELF_RC };
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
        auto parse_range=[&](const char* Q,size_t qlen,Mode mode,size_t lo,size_t hi,
                             std::vector<Ref>& out){
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
                    size_t L=0;
                    while(L<capL && qp+L<qlen && s+L<main_pg_end && Q[qp+L]==pg[s+L]) ++L;
                    if(L>best){ best=L; bsrc=s; }
                }
                return (best>=MINMEM)?best:0;
            };
            size_t qp=lo, lastend=lo;
            while(qp<hi && qp+MINMEM<=qlen){
                uint64_t k;
                if(!packM(Q+qp,k)){ ++qp; continue; }
                const uint32_t idx=lookup(k);
                if(idx==UINT32_MAX){ ++qp; continue; }
                size_t best=0, bestsrc=0, tried=0;
                for(uint32_t i=idx;i<skey.size()&&skey[i]==k;++i){
                    if(++tried>MAXCAND) break;
                    const size_t s=spos[i];
                    // Cap usable length BEFORE extending, never after -- in a
                    // self-match the seed at qp also occurs at qp itself and an
                    // uncapped extension runs to end of text.
                    size_t capL=SIZE_MAX;
                    if(mode==SELF_FWD){ if(s>=qp) continue; capL=qp-s; }
                    else if(mode==SELF_RC){ if(s>=qlen-qp) continue; capL=(qlen-qp-s)/2; }
                    if(capL<MINMEM) continue;
                    size_t L=0;
                    while(L<capL && qp+L<qlen && s+L<main_pg_end && Q[qp+L]==pg[s+L]) ++L;
                    if(L>best){ best=L; bestsrc=s; }
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
                    while(b<cap && Q[qp-b-1]==pg[bestsrc-b-1]) ++b;
                    out.push_back({(uint32_t)(qp-b),(uint32_t)(bestsrc-b),(uint32_t)(best+b),false});
                    lastend=qp+best; qp+=best;
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
        auto run=[&](const char* Q,size_t qlen,Mode mode,std::vector<uint8_t>& consumed,
                     bool RCDEST=false,size_t DSTBASE=0)->size_t{
            unsigned T=std::thread::hardware_concurrency(); if(!T) T=1;
            if(qlen < (1u<<20)) T=1;
            std::vector<std::vector<Ref>> res(T);
            std::vector<std::thread> th;
            const size_t chunk=(qlen+T-1)/T;
            for(unsigned t=0;t<T;++t){
                const size_t lo=(size_t)t*chunk, hi=std::min(qlen,lo+chunk);
                if(lo>=hi) break;
                th.emplace_back([&,t,lo,hi]{ parse_range(Q,qlen,mode,lo,hi,res[t]); });
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
                    if(RCDEST) r.dst=(uint32_t)(qlen-m.dst-m.len);
                    r.dst=(uint32_t)(r.dst+DSTBASE);
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
        {
            std::vector<uint8_t> c(main_pg_end,0), cr(main_pg_end,0);
            if(FWD_SELF) nm_main =run(pg.data(),main_pg_end,SELF_FWD,c);
            { std::string R(pg,0,main_pg_end); rc_inplace(R);
              nm_main+=run(R.data(),main_pg_end,SELF_RC,cr,true); }
            for(size_t i=0;i<main_pg_end;++i) if(cr[i]) c[main_pg_end-1-i]=1;
            for(size_t i=0;i<main_pg_end;++i) if(c[i]) ++rem_main;
        }   // c, cr freed here -- before the second pass allocates its own bitmaps/strings
        {
            const size_t qlen=pg.size()-main_pg_end;
            if(qlen>=MINMEM){
                std::vector<uint8_t> c2(qlen,0), cr2(qlen,0);
                std::string Q(pg,main_pg_end,qlen);
                nm_second =run(Q.data(),qlen,CROSS,c2,false,main_pg_end);
                std::string R=Q; rc_inplace(R);
                nm_second+=run(R.data(),qlen,CROSS,cr2,true,main_pg_end);
                for(size_t i=0;i<qlen;++i) if(cr2[i]) c2[qlen-1-i]=1;
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
        std::vector<uint8_t> gaps,srcs,lens;
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
          for(const Ref& r:allrefs){
              fwrite(&r.dst,4,1,t); fwrite(&r.src,4,1,t); fwrite(&r.len,4,1,t);
              const uint8_t rc=(uint8_t)r.is_rc; fwrite(&rc,1,1,t);
          }
          STR.mem_triples.close(); }
        // mem_gaps/mem_srcs/mem_lens: diagnostics only, dropped
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
        Archive ar(apath);
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

        jobs.push_back({"literal",     [&]{ return seq_encode_mem(litsym, SEQT, SEQT); }});
        jobs.push_back({"mem_triples", [&]{ return refc::encode(v_tri, PGLEN_, MAINEND_); }});
        jobs.push_back({"pos_abs",     [&]{ return best_encode(v_pos.data(), v_pos.size(), true); }});
        jobs.push_back({"pos_strand",  [&]{ return best_encode(v_str.data(), v_str.size(), false); }});
        jobs.push_back({"mm_sym",      [&]{ return mmc::encode(v_mr, v_mo); }});
        jobs.push_back({"mm_pos",      [&]{
            auto flat = best_encode(v_mp.data(), v_mp.size());
            auto buck = mmpos_encode_buckets(v_mp, mmcounts);
            return (!buck.empty() && buck.size() < flat.size()) ? buck : flat; }});
        jobs.push_back({"mm_cnt",      [&]{
            std::vector<uint8_t> cf, cv;
            if(mmcnt_split(mmcounts, cf, cv) && mmcnt_join(cf,cv,mmcounts.size())==mmcounts){
                auto a = best_encode(cf.data(), cf.size());
                auto b = best_encode(cv.data(), cv.size());
                auto flat = best_encode(mmcnt_flat.data(), mmcnt_flat.size());
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
        jobs.push_back({"read_lengths",[&]{ return best_encode(v_rl.data(), v_rl.size()); }});
        jobs.push_back({"orig2uid",    [&]{ return best_encode(v_o2u.data(), v_o2u.size(), true); }});

        results.resize(jobs.size());
        {
            std::atomic<size_t> next{0};
            unsigned NT = HW; if(NT > jobs.size()) NT = (unsigned)jobs.size();
            std::vector<std::thread> pool;
            for(unsigned t=0;t<NT;++t) pool.emplace_back([&]{
                for(;;){ size_t i = next++; if(i >= jobs.size()) break;
                         results[i] = jobs[i].fn(); }
            });
            for(auto& th : pool) th.join();
        }
        // Write in fixed job order -- deterministic regardless of scheduling.
        for(size_t i=0;i<jobs.size();++i){
            if(std::string(jobs[i].name) == "mm_cnt" && mmcnt_is_split){
                ar.put("mm_cnt_flags", mmcnt_flags);
                ar.put("mm_cnt_vals",  mmcnt_vals);
            } else {
                ar.put(jobs[i].name, results[i]);
            }
        }

        ar.finish();
        fprintf(stderr,"[archive] rss before coding %zu MB, after %zu MB\n",rss_before,rss_mb());
    }
    return 0;
}
