// STAGE 71 -- threaded, block-chunked quality coder.
//
// Before building, checked (per user direction) whether quality needed the
// same Genozip-style global-dictionary treatment that names did (stage 70).
// Reread SPRING's real quality source in full
// (/tmp/spring_clone/src/reorder_compress_quality_id.cpp, reorder_compress()
// at lines 127-191): quality has NO dictionary concept at all -- each block
// calls bsc::BSC_str_array_compress() completely independently, zero shared
// state, nothing analogous to the ID dictionary. This matches our own
// coder's structure: stage 68's contexts (value history + position +
// volatility, 24,576 of them) are FIXED-SIZE bounded tables, not a growing
// dictionary like names' ValueDict was -- there is no X-coordinate-style
// problem here to design around. So the simple SPRING-style block-chunked
// design (independent fresh model per block) is the right and ONLY real
// pattern that applies -- no stage-70-style two-pass dictionary needed.
//
// Also cross-checked SPRING's real default block size while there:
// NUM_READS_PER_BLOCK = 256,000 (params.h) -- remarkably close to what we
// independently found best for names (100K-250K sweet spot).
//
// Uses stage 68's exact context formula and Model unchanged -- parallelism
// is a call-site/framing change, not a new algorithm.
//
// Real result, 12 threads, sweeping block size, all round-trip verified
// (baseline single-threaded: 1.003s @ 0.2562 bits/value, matches stage 68
// exactly):
//   block=500,000 (2 blocks):   0.501s  0.2568 bits/value  (+0.2% ratio cost)
//   block=256,000 (4 blocks):   0.259s  0.2578              (+0.6%, SPRING's
//                                                             own real default)
//   block=100,000 (10 blocks):  0.107s  0.2605              (+1.7%)
//   block=50,000  (20 blocks):  0.108s  0.2640              (+3.0%)
// block=100,000 balances thread utilization at 12 cores well: 9.4x speedup
// for a 1.7% ratio cost, still comfortably beating both SPRING (0.2658)
// and Genozip (0.2628).
//
//   g++ -O3 -march=native -pthread -o qualpar 71_qual_blockpar.cpp

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <array>
#include <chrono>
#include <thread>
#include <atomic>
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

static void encode_block(const std::vector<std::string>& reads, size_t s, size_t e, std::vector<uint8_t>& out){
    std::vector<Model> models; models.reserve(NCTX);
    for(size_t i=0;i<NCTX;++i) models.emplace_back(ALPHA_N);
    RangeEnc enc;
    for(size_t r=s;r<e;++r){
        size_t hist=0; int prevq=-1; uint32_t changes=0;
        for(size_t k=0;k<reads[r].size();++k){
            int sym=ALPHA_CODE[(unsigned char)reads[r][k]];
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

static void decode_block(const std::vector<uint8_t>& in, const std::vector<uint16_t>& lens, size_t s, size_t e, std::vector<std::string>& out){
    std::vector<Model> models; models.reserve(NCTX);
    for(size_t i=0;i<NCTX;++i) models.emplace_back(ALPHA_N);
    RangeDec dec; dec.init(in.data(),in.size());
    for(size_t r=s;r<e;++r){
        size_t hist=0; int prevq=-1; uint32_t changes=0;
        std::string got; got.resize(lens[r]);
        for(size_t k=0;k<lens[r];++k){
            size_t ctx=(hist*PBUCK+posBucket(k))*DBUCK+deltaBucket(changes);
            uint32_t sym=models[ctx].dec(dec);
            got[k]=ALPHA_CHARS[sym];
            hist=(hist*ALPHA_N+sym)%HSPACE;
            if(prevq>=0 && prevq!=(int)sym) ++changes;
            prevq=(int)sym;
        }
        out[r-s]=std::move(got);
    }
}

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s quality.txt [block_size=256000] [nthreads=nproc] [hist=3] [pos=8] [delta=48]\n",argv[0]); return 1; }
    const size_t BLOCK=(argc>2)?(size_t)atoll(argv[2]):256000;
    const unsigned NTHREADS=(argc>3)?(unsigned)atoi(argv[3]):std::thread::hardware_concurrency();
    HIST=(argc>4)?atoi(argv[4]):3;
    PBUCK=(argc>5)?atoi(argv[5]):8;
    DBUCK=(argc>6)?atoi(argv[6]):48;

    std::vector<std::string> reads;
    { FILE* f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
      char buf[2048];
      while(fgets(buf,sizeof(buf),f)){ size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
          reads.emplace_back(buf,L); }
      fclose(f); }
    const size_t n=reads.size();
    size_t totalq=0; for(auto& r:reads) totalq+=r.size();
    MAXLEN=0; for(auto& r:reads) MAXLEN=std::max(MAXLEN,r.size());

    Alphabet A;
    for(auto& r:reads) for(char c:r) A.add(c);
    ALPHA_N=A.n; ALPHA_CHARS=A.chars; ALPHA_CODE=A.code;
    HSPACE=1; for(int i=0;i<HIST;++i) HSPACE*=(size_t)ALPHA_N;
    NCTX=HSPACE*PBUCK*DBUCK;

    const size_t nblocks=(n+BLOCK-1)/BLOCK;
    fprintf(stderr,"reads=%zu  total quality bytes=%zu  block=%zu  blocks=%zu  threads=%u  contexts/block=%zu\n",
            n,totalq,BLOCK,nblocks,NTHREADS,NCTX);

    std::vector<std::vector<uint8_t>> blockOut(nblocks);
    auto t0=std::chrono::steady_clock::now();
    {
        std::atomic<size_t> next{0};
        auto worker=[&](){
            size_t b;
            while((b=next.fetch_add(1))<nblocks){
                size_t s=b*BLOCK, e=std::min(n,s+BLOCK);
                encode_block(reads,s,e,blockOut[b]);
            }
        };
        std::vector<std::thread> pool;
        for(unsigned t=0;t<NTHREADS;++t) pool.emplace_back(worker);
        for(auto& th:pool) th.join();
    }
    auto t1=std::chrono::steady_clock::now();
    size_t bytes=0; for(auto& b:blockOut) bytes+=b.size();
    fprintf(stderr,"ENCODE ONLY (parallel, %u threads): %.3f s\n", NTHREADS,
            std::chrono::duration<double>(t1-t0).count());

    std::vector<uint16_t> lens(n); for(size_t r=0;r<n;++r) lens[r]=(uint16_t)reads[r].size();
    std::vector<std::vector<std::string>> blockGot(nblocks);
    {
        std::atomic<size_t> next{0};
        auto worker=[&](){
            size_t b;
            while((b=next.fetch_add(1))<nblocks){
                size_t s=b*BLOCK, e=std::min(n,s+BLOCK);
                blockGot[b].resize(e-s);
                decode_block(blockOut[b],lens,s,e,blockGot[b]);
            }
        };
        std::vector<std::thread> pool;
        for(unsigned t=0;t<NTHREADS;++t) pool.emplace_back(worker);
        for(auto& th:pool) th.join();
    }
    bool ok=true;
    for(size_t b=0;b<nblocks && ok;++b){
        size_t s=b*BLOCK, e=std::min(n,s+BLOCK);
        for(size_t r=s;r<e;++r) if(blockGot[b][r-s]!=reads[r]){ ok=false; fprintf(stderr,"MISMATCH block %zu read %zu\n",b,r); break; }
    }

    printf("block=%zu  threads=%u  hist=%d pos=%d delta=%d  coded=%zu B  %.4f bits/value  round trip: %s\n",
           BLOCK,NTHREADS,HIST,PBUCK,DBUCK,bytes,bytes*8.0/totalq, ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
