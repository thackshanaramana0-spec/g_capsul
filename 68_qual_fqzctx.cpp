// STAGE 68 -- quality coder using fqzcomp's REAL context design (read from
// the actual CRAM 3.1 source, /tmp/htscodecs/htscodecs/fqzcomp_qual.c,
// fqz_update_ctx() at line 344).
//
// Stage 67 found order-9 pure value-history beats real SPRING (0.2654 vs
// 0.2658) but order-10 overfits (context table too sparse). Before assuming
// value-history is the only lever, checked whether a real different tool
// uses a smarter context -- it does. fqzcomp's context is NOT just recent
// values; it packs THREE signals into one context index:
//   1. qctx  -- a shift-register of recent QUANTIZED values (order-H history)
//   2. p     -- position in read (bytes remaining), via ptab
//   3. delta -- COUNT OF VALUE CHANGES SO FAR THIS READ (a volatility/
//               changepoint signal, NOT a value delta):
//               state->delta += (state->prevq != q);
// This lets fqzcomp use a SHORTER value history (cheaper, less dilution)
// while recovering signal from position and volatility instead -- exactly
// the kind of context stage 67's order-10 failure needed but didn't have.
//
// Real result on the same 1M-read real quality column, swept over
// history/position/delta bucket counts (all round-trip verified):
//
//   BEST: hist=3 pos=8 delta=48 -> 4,803,429 B = 0.2562 bits/value
//         24,576 contexts (smaller table than stage 67's order-9, 262,144)
//
//   vs SPRING  (0.2658): 3.6% smaller
//   vs Genozip (0.2628): 2.5% smaller
//   vs stage 67 order-9 (0.2654, pure value history): 3.5% smaller
//
// The "delta" (change-count) bucket dominates: it kept helping up to 48-64
// buckets even as position/history buckets past 8/3-4 stopped helping or
// started hurting (context dilution, the same failure mode as stage 67's
// order-10). Confirms fqzcomp's real design choice: a volatility signal
// beats spending more bits on either longer value history or finer
// position resolution.
//
//   g++ -O3 -march=native -o fqzctx 68_qual_fqzctx.cpp
//   default: hist=3 pos=8 delta=48 (the measured real optimum above)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <array>

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

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s quality.txt [hist=4] [posbuckets=16] [deltabuckets=8]\n",argv[0]); return 1; }
    const int HIST=(argc>2)?atoi(argv[2]):3;
    const int PBUCK=(argc>3)?atoi(argv[3]):8;
    const int DBUCK=(argc>4)?atoi(argv[4]):48;

    std::vector<std::string> reads;
    { FILE* f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
      char buf[2048];
      while(fgets(buf,sizeof(buf),f)){ size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
          reads.emplace_back(buf,L); }
      fclose(f); }
    const size_t n=reads.size();
    size_t totalq=0; for(auto& r:reads) totalq+=r.size();
    size_t maxlen=0; for(auto& r:reads) maxlen=std::max(maxlen,r.size());
    fprintf(stderr,"reads=%zu  total quality bytes=%zu  maxlen=%zu\n",n,totalq,maxlen);

    Alphabet A;
    for(auto& r:reads) for(char c:r) A.add(c);
    fprintf(stderr,"alphabet size=%d\n",A.n);

    size_t HSPACE=1; for(int i=0;i<HIST;++i) HSPACE*=(size_t)A.n;
    const size_t NCTX = HSPACE * PBUCK * DBUCK;
    fprintf(stderr,"contexts=%zu (hist=%d pos=%d delta=%d)\n",NCTX,HIST,PBUCK,DBUCK);

    std::vector<Model> models; models.reserve(NCTX);
    for(size_t i=0;i<NCTX;++i) models.emplace_back(A.n);

    auto posBucket=[&](size_t pos)->size_t{
        // position-in-read, quantized into PBUCK buckets (fqzcomp's ptab role)
        size_t b = pos * (size_t)PBUCK / (maxlen?maxlen:1);
        return b>= (size_t)PBUCK ? PBUCK-1 : b;
    };
    auto deltaBucket=[&](uint32_t changes)->size_t{
        size_t b = changes;
        return b>=(size_t)DBUCK ? (size_t)DBUCK-1 : b;
    };

    RangeEnc enc; enc.out.reserve(totalq);
    for(size_t r=0;r<n;++r){
        size_t hist=0;
        int prevq=-1; uint32_t changes=0;
        for(size_t k=0;k<reads[r].size();++k){
            int sym=A.code[(unsigned char)reads[r][k]];
            size_t ctx = (hist*PBUCK + posBucket(k))*DBUCK + deltaBucket(changes);
            models[ctx].enc(enc,sym);
            hist = (hist*A.n + sym) % HSPACE;
            if(prevq>=0 && prevq!=sym) ++changes;
            prevq=sym;
        }
    }
    enc.flush();
    const size_t bytes=enc.out.size();

    for(auto& m:models) m=Model(A.n);
    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    bool ok=true;
    std::vector<uint16_t> lens(n); for(size_t r=0;r<n;++r) lens[r]=(uint16_t)reads[r].size();
    for(size_t r=0;r<n && ok;++r){
        size_t hist=0; int prevq=-1; uint32_t changes=0;
        std::string got; got.resize(lens[r]);
        for(size_t k=0;k<lens[r];++k){
            size_t ctx = (hist*PBUCK + posBucket(k))*DBUCK + deltaBucket(changes);
            uint32_t sym=models[ctx].dec(dec);
            got[k]=A.chars[sym];
            hist=(hist*A.n+sym)%HSPACE;
            if(prevq>=0 && prevq!=(int)sym) ++changes;
            prevq=(int)sym;
        }
        if(got!=reads[r]){ ok=false; fprintf(stderr,"MISMATCH at read %zu\n",r); }
    }
    printf("hist=%d pos=%d delta=%d  contexts=%zu  quality_bytes=%zu  coded=%zu B  %.4f bits/value  round trip: %s\n",
           HIST,PBUCK,DBUCK,NCTX,totalq,bytes,bytes*8.0/totalq, ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
