// STAGE 67 -- real quality-score coder, order-N Markov context, adaptive
// range coding. Built to get an HONEST, MEASURED number (not a theoretical
// entropy estimate) for the file-order baseline, before testing whether
// SPRING's real technique (reorder reads by assembly/similarity order, THEN
// compress with a strong general compressor -- confirmed by reading
// reorder_compress_quality_id.cpp: lossless mode calls bsc::BSC_str_array_
// compress on ASSEMBLY-ORDERED reads, not a custom quality context model at
// all) is what actually explains the real gap.
//
// Real data here: only 4 distinct quality symbols (binned Illumina quality),
// dominated 93%+ by one value.
//
// REORDERING TESTED AND REJECTED (real negative result, not a guess): two
// independent proxies for assembly order -- full-sequence lexicographic sort
// and sequence-minimizer sort -- were tried on real data (1M reads) with
// BOTH bzip2 (same BWT family as SPRING's real BSC backend) and this coder.
// Neither moved the number at all (bzip2: 0.352 bits/value in file order,
// 0.352 in either sorted order; this coder: 0.2764 file order vs 0.2768/
// 0.2771 sorted -- noise-level). Makes physical sense: quality reflects
// base-calling confidence from the SEQUENCING PROCESS (cycle, optics), not
// the DNA sequence, so sorting reads by sequence content does not cluster
// similar quality profiles. Reordering is SPRING's real mechanism but it is
// NOT why it beats naive Markov models here -- that gap was closed instead
// by simply sweeping the within-read context order:
//
//   order  bits/value  (real data, round-trip VERIFIED)
//   2      0.2823
//   3      0.2764
//   4      0.2726
//   5      0.2699
//   6      0.2679
//   7      0.2664
//   8      0.2655   <- beats real SPRING (0.2658)
//   9      0.2654   <- best; real optimum, 0.24% off real Genozip (0.2628)
//   10     0.2667   <- regresses: contexts too sparse (4^10 contexts over
//                       150M symbols dilutes the per-context sample count,
//                       adaptive rescale kicks in before contexts converge)
//
//   g++ -O3 -march=native -o qualcoder 67_qual_order3.cpp
//   default order is 9 (the measured real optimum above)
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

// Alphabet resolved at runtime from the actual data (small: this real data
// has 4 symbols, but keep it general up to 64 distinct quality chars).
struct Alphabet {
    std::array<int,256> code; int n=0;
    std::array<char,64> chars;
    Alphabet(){ code.fill(-1); }
    int add(char c){ if(code[(unsigned char)c]<0){ code[(unsigned char)c]=n; chars[n]=c; ++n; } return code[(unsigned char)c]; }
};

struct Model {
    std::vector<uint32_t> f; uint32_t N;
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
    if(argc<2){ fprintf(stderr,"usage: %s quality.txt [order=3]\n",argv[0]); return 1; }
    const int ORDER=(argc>2)?atoi(argv[2]):9;
    std::vector<std::string> reads;
    { FILE* f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
      char buf[2048];
      while(fgets(buf,sizeof(buf),f)){ size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
          reads.emplace_back(buf,L); }
      fclose(f); }
    const size_t n=reads.size();
    size_t totalq=0; for(auto& r:reads) totalq+=r.size();
    fprintf(stderr,"reads=%zu  total quality bytes=%zu\n",n,totalq);

    // Resolve alphabet from the real data first (small, cheap pass; needed
    // to size the context table -- not a violation of streaming, this is
    // exactly what the archive header already needs to store).
    Alphabet A;
    for(auto& r:reads) for(char c:r) A.add(c);
    fprintf(stderr,"alphabet size=%d: ",A.n); for(int i=0;i<A.n;++i) fprintf(stderr,"%c ",A.chars[i]); fprintf(stderr,"\n");

    // Context = last ORDER symbols (base-N digits packed into one index).
    size_t NCTX=1; for(int i=0;i<ORDER;++i) NCTX*=(size_t)A.n;
    std::vector<Model> models; models.reserve(NCTX);
    for(size_t i=0;i<NCTX;++i) models.emplace_back(A.n);

    // Read lengths are already carried by the real archive's length stream
    // (this measurement stores them separately, uncounted against quality's
    // own cost, matching how SPRING/Genozip's real numbers are also
    // quality-payload-only).
    RangeEnc enc; enc.out.reserve(totalq);
    std::vector<uint16_t> lens(n);
    for(size_t r=0;r<n;++r) lens[r]=(uint16_t)reads[r].size();
    for(size_t r=0;r<n;++r){
        size_t ctx=0;
        for(char c:reads[r]){
            int sym=A.code[(unsigned char)c];
            models[ctx].enc(enc,sym);
            ctx = (ctx*A.n + sym) % NCTX;
        }
    }
    enc.flush();
    const size_t bytes=enc.out.size();

    // round trip
    for(auto& m:models) m=Model(A.n);
    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    bool ok=true;
    for(size_t r=0;r<n && ok;++r){
        size_t ctx=0;
        std::string got; got.resize(lens[r]);
        for(size_t k=0;k<lens[r];++k){
            uint32_t sym=models[ctx].dec(dec);
            got[k]=A.chars[sym];
            ctx=(ctx*A.n+sym)%NCTX;
        }
        if(got!=reads[r]){ ok=false; fprintf(stderr,"MISMATCH at read %zu\n",r); }
    }
    printf("order=%d  quality_bytes=%zu  coded=%zu B  %.4f bits/value  round trip: %s\n",
           ORDER,totalq,bytes,bytes*8.0/totalq, ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
