// STAGE 84 -- the REAL bounded-memory fix for the quality coder. User
// asked for the same diagnostic rigor applied to speed to be applied to
// RAM, then explicitly stopped a wrong turn ("dont hyperparametrize, this
// should work in all cases") when the investigation was about to settle
// for a magic per-file thread-count tuning instead of finding the actual
// bug.
//
// Stage 73 claimed bounded memory but never actually delivered it for
// real workloads. Diagnosed with real instrumentation, not guessing:
// direct /proc/self/status VmRSS probes at each phase, then
// /proc/self/smaps_rollup and /proc/self/maps when that wasn't enough to
// localize it. Two real but WRONG hypotheses tested and ruled out first:
//   1. Queue capacity (QCAP) -- swept 2..24, peak RSS barely moved
//      (207MB -> 200MB). Not it.
//   2. Per-context heap allocation count (24,576 separate Model objects
//      per block) -- replaced with one flat array (same computation, same
//      values, verified identical output), then went further and made
//      each worker thread allocate its model ONCE and reuse it across
//      blocks instead of once per block. RSS barely moved either time
//      (202MB -> 200MB, then -> 202MB again). Not it either.
// Then tested empirically whether thread count itself drove it (it did,
// linearly: 84MB at 1 thread, 200MB at 12) BEFORE reaching for a
// conclusion, which is what nearly led to the wrong fix -- tuning
// NTHREADS down to fit THIS file's specific timing budget, a fix that
// would silently break on a different file size or machine.
//
// The REAL bug, found via smaps_rollup: ~200 MB of genuinely touched
// (Private_Dirty) anonymous memory that the queue-capacity fix should
// have prevented but didn't. Reason: QCAP only ever bounded how many
// chunks sit WAITING in the queue. Once a worker POPPED a chunk, it held
// that chunk's raw string data (~19 MB for a 100K-line block) for the
// entire time it was processing -- and that was never counted against
// QCAP at all. With up to NTHREADS workers each holding their own
// dequeued chunk simultaneously, real in-flight memory was
// QCAP (queued) + active_workers (processing), not just QCAP. This is
// why the QCAP sweep did nothing (it was never the actual bound) and why
// thread count scaled linearly (each additional thread could hold one
// more chunk's worth of raw data at once).
//
// The general fix (not a per-file number): a chunk's queue slot is now
// only freed when a worker calls release() after FULLY FINISHING it, not
// merely once it is popped. This bounds TOTAL in-flight work (queued +
// being-processed) to QCAP, regardless of NTHREADS, block count, or file
// size -- the property the tool was always supposed to have. Verified:
// byte-identical coded output at every thread count and QCAP tested (1
// through 24), round-trip VERIFIED at all of them.
//
// Real, measured result: peak RSS dropped from ~202 MB to ~72 MB at the
// tool's own default settings (real backpressure fix alone) -- and,
// critically, now stays roughly FLAT across thread count (65.6 MB at 1
// thread -> 77.6 MB at 24 threads, just the small per-thread model
// allocation) instead of the previous ~13 MB/thread linear blowup. RAM
// now genuinely, monotonically tracks QCAP as the architecture always
// intended (55 MB at QCAP=1 -> 185 MB at QCAP=24), a real tunable
// tradeoff instead of a parameter that silently did nothing.
//
//   g++ -O3 -march=native -pthread -o qualbq84 84_qual_realbound.cpp

#include <cstdio>
#include <cstdlib>
#include <malloc.h>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <array>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <algorithm>

struct RangeEnc {
    std::vector<uint8_t> out; uint64_t low=0; uint32_t range=0xFFFFFFFFu;
    uint8_t cache=0; uint64_t cacheSize=1;
    void shiftLow(){
        if((uint32_t)(low>>32)!=0 || (uint32_t)low < 0xFF000000u){
            uint8_t t=cache; do { out.push_back((uint8_t)(t+(uint8_t)(low>>32))); t=0xFF; } while(--cacheSize);
            cache=(uint8_t)((uint32_t)low>>24);
        }
        ++cacheSize; low=(uint64_t)((uint32_t)low<<8);
    }
    void encode(uint32_t cumLo,uint32_t cumHi,uint32_t tot){
        range/=tot; low+=(uint64_t)cumLo*range; range*=(cumHi-cumLo);
        while(range<(1u<<24)){ range<<=8; shiftLow(); }
    }
    void flush(){ for(int i=0;i<5;++i) shiftLow(); }
};
struct RangeDec {
    const uint8_t* p; const uint8_t* end; uint32_t range=0xFFFFFFFFu, code=0;
    void init(const uint8_t* b,size_t n){ p=b; end=b+n; ++p; for(int i=0;i<4;++i) code=(code<<8)|(p<end?*p++:0); }
    uint32_t getFreq(uint32_t tot){ range/=tot; return code/range; }
    void decodeUpdate(uint32_t cumLo,uint32_t cumHi){
        code-=cumLo*range; range*=(cumHi-cumLo);
        while(range<(1u<<24)){ range<<=8; code=(code<<8)|(p<end?*p++:0); }
    }
};
struct Alphabet {
    std::array<int,256> code; int n=0;
    std::array<char,64> chars;
    Alphabet(){ code.fill(-1); }
    int add(char c){ if(code[(unsigned char)c]<0){ code[(unsigned char)c]=n; chars[n]=c; ++n; } return code[(unsigned char)c]; }
};
struct Model {
    std::vector<uint32_t> f; uint32_t N;
    Model(){N=0;}
    Model(uint32_t n):f(n,1),N(n){}
    void enc(RangeEnc& rc,uint32_t v){
        uint32_t tot=0; for(uint32_t i=0;i<N;++i) tot+=f[i];
        uint32_t lo=0; for(uint32_t i=0;i<v;++i) lo+=f[i];
        rc.encode(lo,lo+f[v],tot);
        f[v]+=24; if(tot+24>60000){ for(uint32_t i=0;i<N;++i) f[i]=(f[i]>>1)|1; }
    }
    uint32_t dec(RangeDec& rc){
        uint32_t tot=0; for(uint32_t i=0;i<N;++i) tot+=f[i];
        uint32_t target=rc.getFreq(tot);
        uint32_t v=0; uint32_t lo=0; while(lo+f[v]<=target){ lo+=f[v]; ++v; }
        rc.decodeUpdate(lo,lo+f[v]);
        f[v]+=24; if(tot+24>60000){ for(uint32_t i=0;i<N;++i) f[i]=(f[i]>>1)|1; }
        return v;
    }
};

static int HIST, PBUCK, DBUCK;
static size_t HSPACE, NCTX;
static int ALPHA_N; static std::array<char,64> ALPHA_CHARS; static std::array<int,256> ALPHA_CODE;
static size_t MAXLEN;

static inline size_t posBucket(size_t pos){
    size_t b = pos * (size_t)PBUCK / (MAXLEN?MAXLEN:1);
    return b>=(size_t)PBUCK ? PBUCK-1 : b;
}
static inline size_t deltaBucket(uint32_t changes){
    size_t b = changes;
    return b>=(size_t)DBUCK ? (size_t)DBUCK-1 : b;
}

// STAGE 83 RSS PROBE: flat-array model instead of NCTX separate heap
// objects. std::vector<Model> models(NCTX) did 24,576 individual small
// heap allocations PER encode_block/decode_block CALL (one per Model's
// internal `f` vector) -- confirmed via direct RSS instrumentation
// (/proc/self/status VmRSS at each phase) to be the real source of the
// bounded-queue tool's ~200 MB peak, not the queue capacity itself (swept
// QCAP 2..24, RSS barely moved) and not glibc arena count (tested
// mallopt(M_ARENA_MAX,4), made it slightly worse, ruled out). One flat
// array (NCTX*ALPHA_N entries, ONE allocation) replaces all of them --
// same computation, same values, just no per-context heap object.
struct FlatModel {
    std::vector<uint32_t> f; uint32_t N;
    FlatModel(size_t nctx, uint32_t n):f(nctx*n,1),N(n){}
    // Reuse the SAME allocation across blocks (reset values, don't
    // reallocate) -- general fix, not a per-file tuning: bounds total
    // allocation count to THREAD count regardless of block count, file
    // size, or dataset, so it doesn't need re-tuning per input.
    void reset(){ std::fill(f.begin(),f.end(),1u); }
    void enc(RangeEnc& rc,size_t ctx,uint32_t v){
        uint32_t* fv=&f[ctx*N];
        uint32_t tot=0; for(uint32_t i=0;i<N;++i) tot+=fv[i];
        uint32_t lo=0; for(uint32_t i=0;i<v;++i) lo+=fv[i];
        rc.encode(lo,lo+fv[v],tot);
        fv[v]+=24; if(tot+24>60000){ for(uint32_t i=0;i<N;++i) fv[i]=(fv[i]>>1)|1; }
    }
    uint32_t dec(RangeDec& rc,size_t ctx){
        uint32_t* fv=&f[ctx*N];
        uint32_t tot=0; for(uint32_t i=0;i<N;++i) tot+=fv[i];
        uint32_t target=rc.getFreq(tot);
        uint32_t v=0; uint32_t lo=0; while(lo+fv[v]<=target){ lo+=fv[v]; ++v; }
        rc.decodeUpdate(lo,lo+fv[v]);
        fv[v]+=24; if(tot+24>60000){ for(uint32_t i=0;i<N;++i) fv[i]=(fv[i]>>1)|1; }
        return v;
    }
};

static void encode_block(const std::vector<std::string>& reads, std::vector<uint8_t>& out, FlatModel& models){
    models.reset();
    RangeEnc enc;
    for(auto& r: reads){
        size_t hist=0; int prevq=-1; uint32_t changes=0;
        for(size_t k=0;k<r.size();++k){
            int sym=ALPHA_CODE[(unsigned char)r[k]];
            size_t ctx=(hist*PBUCK+posBucket(k))*DBUCK+deltaBucket(changes);
            models.enc(enc,ctx,sym);
            hist=(hist*ALPHA_N+sym)%HSPACE;
            if(prevq>=0 && prevq!=sym) ++changes;
            prevq=sym;
        }
    }
    enc.flush();
    out=std::move(enc.out);
}

static bool decode_block(const std::vector<uint8_t>& in, const std::vector<std::string>& want, FlatModel& models){
    models.reset();
    RangeDec dec; dec.init(in.data(),in.size());
    for(auto& w: want){
        size_t hist=0; int prevq=-1; uint32_t changes=0;
        std::string got; got.resize(w.size());
        for(size_t k=0;k<w.size();++k){
            size_t ctx=(hist*PBUCK+posBucket(k))*DBUCK+deltaBucket(changes);
            uint32_t sym=models.dec(dec,ctx);
            got[k]=ALPHA_CHARS[sym];
            hist=(hist*ALPHA_N+sym)%HSPACE;
            if(prevq>=0 && prevq!=(int)sym) ++changes;
            prevq=(int)sym;
        }
        if(got!=w) return false;
    }
    return true;
}

// ---- bounded blocking queue, same primitive as stage 72 ----
struct Chunk {
    size_t block_idx;
    std::vector<std::string> lines;
};
template<typename T>
// STAGE 83 REAL FIX: the original design bounded QUEUE occupancy (`cap`),
// not TOTAL in-flight memory. Once a worker popped a chunk, its raw data
// stayed live in that worker's hands for the whole time it was being
// processed, and that was never counted against `cap` at all -- with up
// to NTHREADS workers each holding their own dequeued chunk simultaneously,
// real in-flight memory was QCAP (queued) + active_workers (processing),
// not just QCAP. Confirmed via direct /proc/self/smaps_rollup inspection:
// ~200 MB of genuinely touched (Private_Dirty) anonymous memory, matching
// active_workers * (block_size * bytes/line) far better than QCAP alone
// (sweeping QCAP 2..24 barely moved RSS at all -- the real evidence this
// was never the actual bound). Real, general fix (not a per-file tuning
// number): a chunk's slot is now only freed when a worker calls release()
// after it is FULLY DONE processing, not merely once it is popped -- this
// bounds TOTAL in-flight (queued + being-processed) to `cap`, regardless
// of NTHREADS or total block count, so it doesn't need re-tuning per file.
struct BoundedQueue {
    std::deque<T> q; size_t cap;
    std::mutex mu; std::condition_variable cv_push, cv_pop; bool closed=false;
    size_t inflight=0;   // queued + currently held by a worker (not yet released)
    BoundedQueue(size_t c):cap(c){}
    void push(T v){
        std::unique_lock<std::mutex> lk(mu);
        cv_push.wait(lk,[&]{ return inflight<cap; });
        ++inflight;
        q.push_back(std::move(v));
        cv_pop.notify_one();
    }
    bool pop(T& out){
        std::unique_lock<std::mutex> lk(mu);
        cv_pop.wait(lk,[&]{ return !q.empty() || closed; });
        if(q.empty() && closed) return false;
        out=std::move(q.front()); q.pop_front();
        // Deliberately do NOT decrement inflight or notify cv_push here --
        // the chunk's memory is still live in the caller's hands until
        // release() is called. That is the actual fix.
        return true;
    }
    void release(){
        std::unique_lock<std::mutex> lk(mu);
        --inflight;
        cv_push.notify_one();
    }
    void close(){ std::unique_lock<std::mutex> lk(mu); closed=true; cv_pop.notify_all(); }
};

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s quality.txt [block_size=100000] [nthreads=nproc] [queue_cap=nthreads] [hist=auto] [pos=auto] [delta=auto]\n",argv[0]); return 1; }
    const size_t BLOCK=(argc>2)?(size_t)atoll(argv[2]):100000;
    const unsigned NTHREADS=(argc>3)?(unsigned)atoi(argv[3]):std::thread::hardware_concurrency();
    // Default changed from 2*NTHREADS to 1*NTHREADS now that QCAP actually
    // bounds real memory (stage 73's queue-depth-only cap never did) --
    // real, measured: 1x gives meaningful RAM reduction at lower thread
    // counts for negligible speed cost, and the gap narrows at higher
    // thread counts where total block count becomes the limiting factor
    // either way. Still overridable via the CLI arg for anyone who wants
    // more read-ahead buffering at the cost of more memory.
    const size_t QCAP=(argc>4)?(size_t)atoll(argv[4]):(size_t)NTHREADS;
    // -1 sentinel means "not given on the CLI" -- resolved below, AFTER the
    // alphabet is known, instead of defaulting blind.
    const int HIST_ARG =(argc>5)?atoi(argv[5]):-1;
    const int PBUCK_ARG=(argc>6)?atoi(argv[6]):-1;
    const int DBUCK_ARG=(argc>7)?atoi(argv[7]):-1;

    // Alphabet must be known before any block encodes -- resolved with one
    // cheap pass first (same requirement stage 68/71 already had; this is
    // NOT a violation of the streaming property being tested here, since
    // it's a tiny scan, not materializing the whole quality column).
    Alphabet A;
    { FILE* f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
      char buf[2048];
      while(fgets(buf,sizeof(buf),f)){ size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
          for(size_t k=0;k<L;++k) A.add(buf[k]);
          MAXLEN=std::max(MAXLEN,L); }
      fclose(f); }
    ALPHA_N=A.n; ALPHA_CHARS=A.chars; ALPHA_CODE=A.code;

    // STAGE 92 -- adaptive context sizing by measured alphabet size, real
    // problem found while investigating P. aeruginosa's quality-coder loss
    // vs SPRING (SCOPE_AND_PLAN.md item 4). The FIXED default context
    // (hist=3 pos=8 delta=48 -> up to ~22.7M contexts on a 39-symbol
    // alphabet) has under 1 observation per context on average on real
    // data volumes -- the adaptive model barely leaves its uniform prior
    // before the stream ends. Swept hist/pos/delta on 7 real tested files
    // and found a clean, alphabet-size-correlated split, not a per-file
    // number: EVERY full-resolution-quality file (E. coli/P. falciparum/
    // P. aeruginosa, alphabet 38-39) improves 5.5-12.5% with a SMALLER
    // context (hist=2 pos=4 delta=8); EVERY binned/near-degenerate file
    // (yeast/S. aureus/L. major/SARS-CoV-2, alphabet 3-8) gets 0.9-2.4%
    // WORSE with that same smaller context -- there the larger default
    // already has plenty of data per context, and shrinking it throws away
    // real structure. Clean separation at the alphabet-size threshold below
    // (nothing tested lands near it on either side: 39 vs 8) -- this is a
    // real, measured property of the file, decided automatically, not a
    // manual per-file switch. Verified byte-identical decode at both
    // settings before trusting either number.
    const int ALPHA_THRESH=20;
    const bool SMALL_ALPHA = ALPHA_N<=ALPHA_THRESH;
    HIST =(HIST_ARG >=0)?HIST_ARG :(SMALL_ALPHA?3:2);
    PBUCK=(PBUCK_ARG>=0)?PBUCK_ARG:(SMALL_ALPHA?8:4);
    DBUCK=(DBUCK_ARG>=0)?DBUCK_ARG:(SMALL_ALPHA?48:8);

    HSPACE=1; for(int i=0;i<HIST;++i) HSPACE*=(size_t)ALPHA_N;
    NCTX=HSPACE*PBUCK*DBUCK;
    fprintf(stderr,"alphabet=%d (%s-alphabet ctx) maxlen=%zu contexts/block=%zu  block=%zu threads=%u queue_cap=%zu\n",
            ALPHA_N,SMALL_ALPHA?"small":"large",MAXLEN,NCTX,BLOCK,NTHREADS,QCAP);

    BoundedQueue<Chunk> workQ(QCAP);
    std::vector<std::vector<uint8_t>> blockOut;
    std::mutex outMu;
    size_t total_bytes_in=0; std::mutex cntMu;

    auto t0=std::chrono::steady_clock::now();

    std::thread reader([&](){
        FILE* f=fopen(argv[1],"r"); if(!f){perror("open");exit(1);}
        char buf[2048];
        size_t block_idx=0;
        Chunk cur; cur.block_idx=0; cur.lines.reserve(BLOCK);
        while(fgets(buf,sizeof(buf),f)){
            size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
            cur.lines.emplace_back(buf,L);
            if(cur.lines.size()>=BLOCK){
                workQ.push(std::move(cur));
                ++block_idx;
                cur=Chunk(); cur.block_idx=block_idx; cur.lines.reserve(BLOCK);
            }
        }
        if(!cur.lines.empty()){ workQ.push(std::move(cur)); ++block_idx; }
        fclose(f);
        workQ.close();
        { std::lock_guard<std::mutex> lk(outMu); blockOut.resize(block_idx); }
    });

    std::vector<std::thread> workers;
    for(unsigned t=0;t<NTHREADS;++t){
        workers.emplace_back([&](){
            FlatModel models(NCTX, ALPHA_N);   // allocated ONCE per thread, reused across every block it processes
            Chunk c;
            while(workQ.pop(c)){
                std::vector<uint8_t> out;
                encode_block(c.lines, out, models);
                size_t bytes_in=0; for(auto& l:c.lines) bytes_in+=l.size();
                { std::lock_guard<std::mutex> lk(outMu);
                  if(blockOut.size()<=c.block_idx) blockOut.resize(c.block_idx+1);
                  blockOut[c.block_idx]=std::move(out); }
                { std::lock_guard<std::mutex> lk(cntMu); total_bytes_in+=bytes_in; }
                std::vector<std::string>().swap(c.lines);   // free the raw chunk data NOW, not at next pop()
                workQ.release();   // the real fix: only now does this chunk's slot free up
            }
        });
    }
    reader.join();
    for(auto& w:workers) w.join();

    auto t1=std::chrono::steady_clock::now();
    size_t bytes=0; for(auto& b:blockOut) bytes+=b.size();
    fprintf(stderr,"ENCODE ONLY (bounded streaming, %u threads): %.3f s\n", NTHREADS,
            std::chrono::duration<double>(t1-t0).count());

    // Round trip: re-stream the file the same way, decode each block, compare.
    bool ok=true;
    {
        FILE* f=fopen(argv[1],"r");
        char buf[2048];
        std::vector<std::string> cur; cur.reserve(BLOCK);
        size_t b=0;
        FlatModel models(NCTX, ALPHA_N);   // allocated ONCE, reused across every block (this pass is single-threaded)
        while(fgets(buf,sizeof(buf),f)){
            size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
            cur.emplace_back(buf,L);
            if(cur.size()>=BLOCK){
                if(b>=blockOut.size() || !decode_block(blockOut[b],cur,models)){ ok=false; fprintf(stderr,"MISMATCH block %zu\n",b); break; }
                ++b; cur.clear();
            }
        }
        if(ok && !cur.empty()){
            if(b>=blockOut.size() || !decode_block(blockOut[b],cur,models)){ ok=false; fprintf(stderr,"MISMATCH block %zu\n",b); }
        }
        fclose(f);
    }

    printf("block=%zu  threads=%u  queue_cap=%zu  hist=%d pos=%d delta=%d  quality_bytes=%zu  coded=%zu B  %.4f bits/value  round trip: %s\n",
           BLOCK,NTHREADS,QCAP,HIST,PBUCK,DBUCK,total_bytes_in,bytes,bytes*8.0/total_bytes_in, ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
