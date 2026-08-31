// Order-k context model + range coder for the pseudogenome literal.
//
// The sequence gap against PgRC2 is 101,610 B, and it splits almost evenly:
//   232,022 extra literal bases  -> 55,863 B at their rate
//   coder deficit (1.9544 vs 1.9261 bpb) -> 45,747 B
// This attacks the second half.
//
// What they do: VarLenDNACoder maps 1-4 bases to one byte from a 242-phrase
// dictionary, then LZMA. The dictionary is ordered by the phrase's LAST base
// (VarLenDNACoder.cpp:200-219 -- all phrases ending A, then C, then G, then T),
// so a code byte's high bits carry the terminal base and the byte stream keeps
// structure LZMA can match on. That is also why it beats naive 2-bit packing:
// packing four bases per byte at fixed boundaries destroys the alignment of any
// repeat that does not start at a multiple of four, and LZMA then cannot see it.
//
// Why a context model should win here instead: this literal is what SURVIVED
// MEM removal, which already stripped every exact repeat of 45 bases or more.
// The stream is therefore close to repeat-free, so LZMA spends bits on
// match/literal signalling it can rarely use, while an order-k model spends
// everything on the only thing left -- the next base given its context.
//
// Textbook: order-k adaptive frequency counts with periodic halving, coded with
// the same LZMA-style range coder used in stage 23. Nothing from PgRC2; their
// path has no context model.
//
//   g++ -O3 -march=native -o cm 29_dna_cm.cpp && ./cm literal.txt [k]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>

struct RangeEnc {
    std::vector<uint8_t> out; uint64_t low=0; uint32_t range=0xFFFFFFFFu;
    uint8_t cache=0; uint64_t cacheSize=1;
    void shiftLow(){
        if((uint32_t)(low>>32)!=0 || (uint32_t)low < 0xFF000000u){
            uint8_t t=cache;
            do { out.push_back((uint8_t)(t+(uint8_t)(low>>32))); t=0xFF; } while(--cacheSize);
            cache=(uint8_t)((uint32_t)low>>24);
        }
        ++cacheSize; low=(uint64_t)((uint32_t)low<<8);
    }
    // encode symbol with cumulative frequency [cum, cum+f) out of tot
    void encode(uint32_t cum,uint32_t f,uint32_t tot){
        range/=tot; low+=(uint64_t)cum*range; range*=f;
        while(range<(1u<<24)){ range<<=8; shiftLow(); }
    }
    void flush(){ for(int i=0;i<5;++i) shiftLow(); }
};

int main(int argc,char** argv){
    if(argc<2){ fprintf(stderr,"usage: %s <literal.txt> [k]\n",argv[0]); return 1; }
    const int K = argc>2 ? atoi(argv[2]) : 12;          // context length in bases
    std::vector<uint8_t> sym;
    {
        FILE* f=fopen(argv[1],"rb"); if(!f){ perror("open"); return 1; }
        fseek(f,0,SEEK_END); const size_t sz=ftell(f); fseek(f,0,SEEK_SET);
        std::vector<char> buf(sz);
        if(fread(buf.data(),1,sz,f)!=sz){ fprintf(stderr,"short read\n"); return 1; }
        fclose(f);
        sym.reserve(sz);
        for(size_t i=0;i<sz;++i){
            switch(buf[i]){ case 'A':sym.push_back(0);break; case 'C':sym.push_back(1);break;
                            case 'G':sym.push_back(2);break; case 'T':sym.push_back(3);break;
                            default: break; }
        }
    }
    const size_t n=sym.size();
    const size_t CTX=(size_t)1<<(2*K);
    // uint16 counts, 4 per context. At K=12 that is 16M contexts = 128 MB.
    std::vector<uint16_t> cnt;
    cnt.assign(CTX*4,1);                                  // start uniform
    const uint32_t LIMIT=1<<13;                           // halve before overflow

    RangeEnc enc; enc.out.reserve(n/4);
    uint64_t ctx=0; const uint64_t CMASK=CTX-1;
    for(size_t i=0;i<n;++i){
        uint16_t* c=&cnt[(size_t)ctx*4];
        const uint32_t tot=(uint32_t)c[0]+c[1]+c[2]+c[3];
        const uint8_t s=sym[i];
        uint32_t cum=0; for(uint8_t j=0;j<s;++j) cum+=c[j];
        enc.encode(cum,c[s],tot);
        c[s]+=32;
        if(tot+32>=LIMIT){ for(int j=0;j<4;++j) c[j]=(uint16_t)((c[j]>>1)|1); }
        ctx=((ctx<<2)|s)&CMASK;
    }
    enc.flush();
    const size_t bytes=enc.out.size();
    printf("k=%-3d bases=%zu  coded=%zu B  %.4f bits/base\n",K,n,bytes,bytes*8.0/n);
    printf("      PgRC2 on their literal: 3,056,474 B over 12,694,903 = 1.9261 bits/base\n");
    printf("      2-bit + xz -9e on ours: 3,158,084 B over 12,926,925 = 1.9544 bits/base\n");
    return 0;
}
