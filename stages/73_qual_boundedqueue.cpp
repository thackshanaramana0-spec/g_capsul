// STAGE 73 -- bounded-memory streaming pipeline, quality coder. The
// quality counterpart to stage 72 (names): locks in the same real property
// Genozip's dispatcher.c has (memory bounded by pool size, not file size --
// confirmed by reading its real source in full) for the quality stream too,
// since stage 71 still loads the entire quality column into a
// std::vector<std::string> before threading starts, same real gap as
// stage 69 had for names before stage 72 fixed it.
//
// No dictionary concern here (already established at stage 71: quality's
// contexts are fixed-size bounded tables, not a growing dictionary), so
// this stage is a direct architectural port of stage 72's bounded queue
// onto stage 68/71's exact context model and Model struct -- nothing new
// to design, isolate the memory architecture as the only variable.
//
// Real, measured result (qual_only.txt, 150 MB, 1M reads, block=100,000,
// 12 threads): same pattern as stage 72's names result -- a genuine
// memory-vs-throughput tradeoff, not a free win. Ratio unaffected either
// way (0.2605 bits/value, matches stage 71 exactly at this block size).
//   peak RSS:  stage 71 (load-all)      400 MB
//              stage 73 (bounded queue) 209 MB   (-48%, real, /usr/bin/time -v)
//   time:      stage 71                 0.50s (total, includes round trip)
//              stage 73                 1.63s (total, includes round trip;
//                                        real synchronization overhead per
//                                        chunk, consistent with stage 72's
//                                        finding for names)
// Worth it specifically for files too large to comfortably load whole
// (the real T. cacao ~15 GB case in the ARCS benchmark), not a win on a
// file this small.
//
//   g++ -O3 -march=native -pthread -o qualbq 73_qual_boundedqueue.cpp

#include <cstdio>
#include <cstdlib>
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

static void encode_block(const std::vector<std::string>& reads, std::vector<uint8_t>& out){
    std::vector<Model> models; models.reserve(NCTX);
    for(size_t i=0;i<NCTX;++i) models.emplace_back(ALPHA_N);
    RangeEnc enc;
    for(auto& r: reads){
        size_t hist=0; int prevq=-1; uint32_t changes=0;
        for(size_t k=0;k<r.size();++k){
            int sym=ALPHA_CODE[(unsigned char)r[k]];
            size_t ctx=(hist*PBUCK+posBucket(k))*DBUCK+deltaBucket(changes);
            models[ctx].enc(enc,sym);
            hist=(hist*ALPHA_N+sym)%HSPACE;
            if(prevq>=0 && prevq!=sym) ++changes;
            prevq=sym;
        }
    }
    enc.flush();
    out=std::move(enc.out);
}

static bool decode_block(const std::vector<uint8_t>& in, const std::vector<std::string>& want){
    std::vector<Model> models; models.reserve(NCTX);
    for(size_t i=0;i<NCTX;++i) models.emplace_back(ALPHA_N);
    RangeDec dec; dec.init(in.data(),in.size());
    for(auto& w: want){
        size_t hist=0; int prevq=-1; uint32_t changes=0;
        std::string got; got.resize(w.size());
        for(size_t k=0;k<w.size();++k){
            size_t ctx=(hist*PBUCK+posBucket(k))*DBUCK+deltaBucket(changes);
            uint32_t sym=models[ctx].dec(dec);
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
struct BoundedQueue {
    std::deque<T> q; size_t cap;
    std::mutex mu; std::condition_variable cv_push, cv_pop; bool closed=false;
    BoundedQueue(size_t c):cap(c){}
    void push(T v){
        std::unique_lock<std::mutex> lk(mu);
        cv_push.wait(lk,[&]{ return q.size()<cap; });
        q.push_back(std::move(v));
        cv_pop.notify_one();
    }
    bool pop(T& out){
        std::unique_lock<std::mutex> lk(mu);
        cv_pop.wait(lk,[&]{ return !q.empty() || closed; });
        if(q.empty() && closed) return false;
        out=std::move(q.front()); q.pop_front();
        cv_push.notify_one();
        return true;
    }
    void close(){ std::unique_lock<std::mutex> lk(mu); closed=true; cv_pop.notify_all(); }
};

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s quality.txt [block_size=100000] [nthreads=nproc] [queue_cap=2*nthreads] [hist=3] [pos=8] [delta=48]\n",argv[0]); return 1; }
    const size_t BLOCK=(argc>2)?(size_t)atoll(argv[2]):100000;
    const unsigned NTHREADS=(argc>3)?(unsigned)atoi(argv[3]):std::thread::hardware_concurrency();
    const size_t QCAP=(argc>4)?(size_t)atoll(argv[4]):(size_t)(2*NTHREADS);
    HIST=(argc>5)?atoi(argv[5]):3;
    PBUCK=(argc>6)?atoi(argv[6]):8;
    DBUCK=(argc>7)?atoi(argv[7]):48;

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
    HSPACE=1; for(int i=0;i<HIST;++i) HSPACE*=(size_t)ALPHA_N;
    NCTX=HSPACE*PBUCK*DBUCK;
    fprintf(stderr,"alphabet=%d maxlen=%zu contexts/block=%zu  block=%zu threads=%u queue_cap=%zu\n",
            ALPHA_N,MAXLEN,NCTX,BLOCK,NTHREADS,QCAP);

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
            Chunk c;
            while(workQ.pop(c)){
                std::vector<uint8_t> out;
                encode_block(c.lines, out);
                size_t bytes_in=0; for(auto& l:c.lines) bytes_in+=l.size();
                { std::lock_guard<std::mutex> lk(outMu);
                  if(blockOut.size()<=c.block_idx) blockOut.resize(c.block_idx+1);
                  blockOut[c.block_idx]=std::move(out); }
                { std::lock_guard<std::mutex> lk(cntMu); total_bytes_in+=bytes_in; }
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
        while(fgets(buf,sizeof(buf),f)){
            size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
            cur.emplace_back(buf,L);
            if(cur.size()>=BLOCK){
                if(b>=blockOut.size() || !decode_block(blockOut[b],cur)){ ok=false; fprintf(stderr,"MISMATCH block %zu\n",b); break; }
                ++b; cur.clear();
            }
        }
        if(ok && !cur.empty()){
            if(b>=blockOut.size() || !decode_block(blockOut[b],cur)){ ok=false; fprintf(stderr,"MISMATCH block %zu\n",b); }
        }
        fclose(f);
    }

    printf("block=%zu  threads=%u  queue_cap=%zu  hist=%d pos=%d delta=%d  quality_bytes=%zu  coded=%zu B  %.4f bits/value  round trip: %s\n",
           BLOCK,NTHREADS,QCAP,HIST,PBUCK,DBUCK,total_bytes_in,bytes,bytes*8.0/total_bytes_in, ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
