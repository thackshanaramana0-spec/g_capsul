// STAGE 53 -- widen stage 50's mismatch coder from 4 contexts (ref base
// alone) to a cross-stream context: ref base x position-in-read-bucket x
// preceding pg base. This is the "joint/cross-stream coding" claim made
// concrete: mm_pos.bin and mm_ctx3.bin carry values from streams (position,
// sequence) that are ALREADY DECODED by the time the mismatch stream is
// coded in this pipeline's own stream order -- the same "condition on
// already-decoded prior data" pattern read out of real CABAC/channel-wise-
// autoregressive code, just with our existing division-based range coder
// (unchanged) and a wider context array (the only thing that changes).
//
// SAME coder mechanism as stage 50 -- RangeEnc/RangeDec, exclusive coding,
// adaptive frequency tables with the same rescale rule. Nothing structural
// is new; only which context selects which table.
//
//   g++ -O3 -march=native -o mmcoder53 53_crossstream_mismatch.cpp

// Real adaptive arithmetic coder for the mismatch symbol stream: exclusive
// coding (3 remaining symbols, since decoder knows the reference base)
// conditioned on the reference base, adapted online. Round trip verified.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
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
int main(int argc,char**argv){
    if(argc<5){ fprintf(stderr,"usage: %s ref.bin obs.bin pos.bin ctx3.bin\n",argv[0]); return 1; }
    FILE* fr=fopen(argv[1],"rb"); FILE* fo=fopen(argv[2],"rb");
    FILE* fp=fopen(argv[3],"rb"); FILE* fx=fopen(argv[4],"rb");
    fseek(fr,0,SEEK_END); size_t n=ftell(fr); fseek(fr,0,SEEK_SET);
    std::vector<uint8_t> ref(n), obs(n), pos(n), ctx3(n);
    if(fread(ref.data(),1,n,fr)!=n) return 1; if(fread(obs.data(),1,n,fo)!=n) return 1;
    if(fread(pos.data(),1,n,fp)!=n) return 1; if(fread(ctx3.data(),1,n,fx)!=n) return 1;
    fclose(fr); fclose(fo); fclose(fp); fclose(fx);
    auto code2=[](uint8_t c)->int{ return c=='A'?0:c=='C'?1:c=='G'?2:3; };

    // Context = ref base (4) x position-in-read bucket (3: <50, <100, >=100,
    // matching the 3' quality-drop pattern the earlier entropy estimate used)
    // x preceding pg base (4) = 48 contexts. SAME 3-slot adaptive table per
    // context, same update rule, as stage 50 -- only the context COUNT and
    // SELECTION changed.
    const int NCTX = 4*3*4;
    std::vector<std::array<uint32_t,3>> freq(NCTX);
    for(auto& f:freq) f={1,1,1};

    auto excludeIdx=[&](int refc,int obsc)->int{
        int k=0; for(int s=0;s<4;++s){ if(s==refc) continue; if(s==obsc) return k; ++k; } return -1; };
    auto includeSym=[&](int refc,int k)->int{
        int c=0; for(int s=0;s<4;++s){ if(s==refc) continue; if(c==k) return s; ++c; } return -1; };
    auto posBucket=[&](uint8_t p)->int{ return p<50?0:p<100?1:2; };
    auto ctxIdx=[&](int r,int pb,int pv)->int{ return (r*3+pb)*4+pv; };

    RangeEnc enc; enc.out.reserve(n);
    for(size_t i=0;i<n;++i){
        const int r=code2(ref[i]), o=code2(obs[i]);
        const int pb=posBucket(pos[i]), pv=code2(ctx3[i]);
        const int cx=ctxIdx(r,pb,pv);
        const int k=excludeIdx(r,o);
        auto& f=freq[cx]; uint32_t tot=f[0]+f[1]+f[2];
        uint32_t lo=0; for(int j=0;j<k;++j) lo+=f[j];
        enc.encode(lo,lo+f[k],tot);
        f[k]+=8; if(tot+8>65536){ for(int j=0;j<3;++j) f[j]=(f[j]>>1)|1; }
    }
    enc.flush();
    const size_t bytes=enc.out.size();

    // round trip
    for(auto& f:freq) f={1,1,1};
    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    bool ok=true;
    for(size_t i=0;i<n;++i){
        const int r=code2(ref[i]);
        const int pb=posBucket(pos[i]), pv=code2(ctx3[i]);
        const int cx=ctxIdx(r,pb,pv);
        auto& f=freq[cx]; uint32_t tot=f[0]+f[1]+f[2];
        uint32_t target=dec.getFreq(tot);
        int k=0; uint32_t lo=0; while(lo+f[k]<=target){ lo+=f[k]; ++k; }
        dec.decodeUpdate(lo,lo+f[k]);
        const int s=includeSym(r,k);
        const int wanted=code2(obs[i]);
        if(s!=wanted){ ok=false; fprintf(stderr,"MISMATCH at %zu: got %d want %d\n",i,s,wanted); break; }
        f[k]+=8; if(tot+8>65536){ for(int j=0;j<3;++j) f[j]=(f[j]>>1)|1; }
    }
    printf("symbols=%zu  coded=%zu B  %.4f bits/symbol  round trip: %s\n",
           n,bytes,bytes*8.0/n, ok?"VERIFIED":"FAILED");
    printf("stage 50 (ref-only, 4 ctx): 241,365 B  1.5446 bits/symbol\n");
    return ok?0:1;
}
