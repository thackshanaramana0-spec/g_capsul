// STAGE 50 -- real mismatch coder, replacing the estimate every prior number
// in this progression used (PgRC2's own observed 1.55 bits/mismatch).
//
// Exclusive coding (the decoder already knows the reference base, so the
// symbol is one of 3, not 4 -- same shape PgRC2 uses,
// reorderingSymbolsExclusiveMismatchEncoding), conditioned on the reference
// base, adapted online with a range coder. Round-trip verified: 241,365 B
// against the 242,209 B estimate on yeast_sub.fq's 1,250,111 mismatches.
//
// Real Illumina substitution bias is present in this data (C<->G rarer than
// the other ten, matching the literature) but mild here -- almost all the
// achievable gain over naive 2-bit coding is from exclusive coding alone
// (1.585 bits/mismatch); longer context (trinucleotide, position-in-read)
// adds under 2,000 B more by entropy estimate and was not worth a second
// coder for that little.
//
// Usage: ./mmcoder mm_ref.bin mm_obs.bin  (both produced by stage 47's
// DUMP_MM=1 -- one byte per mismatch, the reference base and observed base)
//
//   g++ -O3 -march=native -o mmcoder 50_mismatch_coder_real.cpp

// Real adaptive arithmetic coder for the mismatch symbol stream: exclusive
// coding (3 remaining symbols, since decoder knows the reference base)
// conditioned on the reference base, adapted online. Round trip verified.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
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
int main(int argc,char**argv){
    if(argc<3){ fprintf(stderr,"usage: %s ref.bin obs.bin\n",argv[0]); return 1; }
    FILE* fr=fopen(argv[1],"rb"); FILE* fo=fopen(argv[2],"rb");
    fseek(fr,0,SEEK_END); size_t n=ftell(fr); fseek(fr,0,SEEK_SET);
    std::vector<uint8_t> ref(n), obs(n);
    if(fread(ref.data(),1,n,fr)!=n) return 1; if(fread(obs.data(),1,n,fo)!=n) return 1;
    fclose(fr); fclose(fo);
    auto code2=[](uint8_t c)->int{ return c=='A'?0:c=='C'?1:c=='G'?2:3; };

    // 4 contexts (by ref base), each an adaptive frequency table over the 3
    // remaining symbols (indices 0,1,2 after excluding ref itself).
    uint32_t freq[4][3]; for(int r=0;r<4;++r) for(int k=0;k<3;++k) freq[r][k]=1;

    auto excludeIdx=[&](int refc,int obsc)->int{
        int k=0; for(int s=0;s<4;++s){ if(s==refc) continue; if(s==obsc) return k; ++k; } return -1; };
    auto includeSym=[&](int refc,int k)->int{
        int c=0; for(int s=0;s<4;++s){ if(s==refc) continue; if(c==k) return s; ++c; } return -1; };

    RangeEnc enc; enc.out.reserve(n);
    for(size_t i=0;i<n;++i){
        const int r=code2(ref[i]), o=code2(obs[i]);
        const int k=excludeIdx(r,o);
        uint32_t* f=freq[r]; uint32_t tot=f[0]+f[1]+f[2];
        uint32_t lo=0; for(int j=0;j<k;++j) lo+=f[j];
        enc.encode(lo,lo+f[k],tot);
        f[k]+=8; if(tot+8>65536){ for(int j=0;j<3;++j) f[j]=(f[j]>>1)|1; }
    }
    enc.flush();
    const size_t bytes=enc.out.size();

    // round trip
    for(int r=0;r<4;++r) for(int k=0;k<3;++k) freq[r][k]=1;
    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    bool ok=true;
    for(size_t i=0;i<n;++i){
        const int r=code2(ref[i]);
        uint32_t* f=freq[r]; uint32_t tot=f[0]+f[1]+f[2];
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
    return ok?0:1;
}
