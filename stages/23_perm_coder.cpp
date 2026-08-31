// Order information (PgRC2 stage 6) coded at the permutation entropy floor.
//
// PgRC2 stores the inverse permutation as a raw uint32 array and hands it to
// LZMA (SeparatedPseudoGenomePersistence.cpp:226-233) -- no permutation-specific
// coding at all in single-end mode. From its own log on yeast_sub.fq:
//
//     lzma ... compressed 4000000 bytes to 2852758 bytes (ratio 0.713)
//
// 2,852,758 B for 1M reads = 22.82 bits/read. A general-purpose compressor
// cannot do better than that in principle, because it does not know the array
// is a permutation. That single fact is worth a lot: once you know every value
// appears exactly once, the i-th element only has to be identified among the
// n-i values still unused, so the total content is
//
//     sum_{i=1..n} log2(i)  =  log2(n!)  =  2,309,466 B  =  18.49 bits/read
//
// which is 19% below what they spend. This encodes exactly that, by mapping the
// permutation to its Lehmer code (each element replaced by its rank among the
// values not yet used) and coding digit i as a uniform symbol in [0, n-i).
// Ranks come from a Fenwick tree over the still-available values, so the whole
// pass is O(n log n).
//
// Textbook parts throughout: the Lehmer/factorial-number-system code is
// classical, the Fenwick tree is standard, and the range coder is the usual
// LZMA-style carry-cached encoder. Nothing here is taken from PgRC2 -- their
// single-end path has no permutation coder to take.
//
// Round-trips by construction and it is checked: the decoder rebuilds the
// permutation and the result is compared element by element before any size is
// reported. A number that has not survived a decode is not evidence.
//
//   g++ -O3 -march=native -o pc 23_perm_coder.cpp
//   ./pc perm.u32
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <cmath>

// ── Fenwick tree over "value still available" ────────────────────────────────
struct Fenwick {
    std::vector<uint32_t> t; uint32_t n, LOG;
    void init(uint32_t n_){
        n=n_; t.assign(n+1,0); LOG=0; while((1u<<(LOG+1))<=n) ++LOG;
        // every value starts available; build in O(n)
        for(uint32_t i=1;i<=n;++i){
            t[i]+=1;
            const uint32_t j=i+(i&(-(int32_t)i));
            if(j<=n) t[j]+=t[i];
        }
    }
    // number of available values strictly below v  (v is 0-based)
    uint32_t rank(uint32_t v) const {
        uint32_t s=0; for(uint32_t i=v;i;i-=i&(-(int32_t)i)) s+=t[i]; return s;
    }
    void remove(uint32_t v){
        for(uint32_t i=v+1;i<=n;i+=i&(-(int32_t)i)) --t[i];
    }
    // smallest v with rank(v+1) == k+1, i.e. the (k+1)-th available value
    uint32_t select(uint32_t k) const {
        uint32_t pos=0, rem=k+1;
        for(int32_t pw=(int32_t)1<<LOG; pw>0; pw>>=1){
            const uint32_t np=pos+(uint32_t)pw;
            if(np<=n && t[np]<rem){ pos=np; rem-=t[np]; }
        }
        return pos;   // 0-based value
    }
};

// ── LZMA-style range coder, uniform symbols ──────────────────────────────────
struct RangeEnc {
    std::vector<uint8_t> out; uint64_t low=0; uint32_t range=0xFFFFFFFFu;
    uint8_t cache=0; uint64_t cacheSize=1;
    void shiftLow(){
        if((uint32_t)(low>>32)!=0 || (uint32_t)low < 0xFF000000u){
            uint8_t temp=cache;
            do { out.push_back((uint8_t)(temp+(uint8_t)(low>>32))); temp=0xFF; } while(--cacheSize);
            cache=(uint8_t)((uint32_t)low>>24);
        }
        ++cacheSize; low=(uint64_t)((uint32_t)low<<8);
    }
    void encode(uint32_t v,uint32_t m){          // v in [0,m)
        range/=m; low+=(uint64_t)v*range;
        while(range<(1u<<24)){ range<<=8; shiftLow(); }
    }
    void flush(){ for(int i=0;i<5;++i) shiftLow(); }
};
struct RangeDec {
    const uint8_t* p; const uint8_t* end; uint32_t range=0xFFFFFFFFu, code=0;
    void init(const uint8_t* b,size_t n){
        p=b; end=b+n; ++p;                        // first byte is always 0
        for(int i=0;i<4;++i) code=(code<<8)|(p<end?*p++:0);
    }
    uint32_t decode(uint32_t m){
        range/=m; uint32_t v=code/range; if(v>=m) v=m-1;
        code-=v*range;
        while(range<(1u<<24)){ range<<=8; code=(code<<8)|(p<end?*p++:0); }
        return v;
    }
};

int main(int argc,char** argv){
    if(argc<2){ fprintf(stderr,"usage: %s <perm.u32>\n",argv[0]); return 1; }
    std::vector<uint32_t> P;
    {
        FILE* f=fopen(argv[1],"rb"); if(!f){ perror("open"); return 1; }
        fseek(f,0,SEEK_END); const size_t sz=ftell(f); fseek(f,0,SEEK_SET);
        P.resize(sz/4);
        if(fread(P.data(),4,P.size(),f)!=P.size()){ fprintf(stderr,"short read\n"); return 1; }
        fclose(f);
    }
    const uint32_t n=(uint32_t)P.size();
    // The input must be a genuine permutation of 0..n-1, or the Lehmer mapping
    // is meaningless. Checked rather than assumed.
    {
        std::vector<uint8_t> seen(n,0);
        for(uint32_t i=0;i<n;++i){
            if(P[i]>=n||seen[P[i]]){ fprintf(stderr,"not a permutation at %u\n",i); return 1; }
            seen[P[i]]=1;
        }
    }
    double floorBits=0; for(uint32_t i=1;i<=n;++i) floorBits+=std::log2((double)i);

    // ── encode ───────────────────────────────────────────────────────────────
    Fenwick fw; fw.init(n);
    RangeEnc enc; enc.out.reserve((size_t)(floorBits/8)+1024);
    for(uint32_t i=0;i<n;++i){
        const uint32_t m=n-i;                     // values still unused
        if(m>1) enc.encode(fw.rank(P[i]),m);      // last digit carries no information
        fw.remove(P[i]);
    }
    enc.flush();
    const size_t bytes=enc.out.size();

    // ── decode and verify ────────────────────────────────────────────────────
    Fenwick fw2; fw2.init(n);
    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    bool ok=true;
    for(uint32_t i=0;i<n;++i){
        const uint32_t m=n-i;
        const uint32_t d=(m>1)?dec.decode(m):0;
        const uint32_t v=fw2.select(d);
        if(v!=P[i]){ fprintf(stderr,"MISMATCH at %u: got %u want %u\n",i,v,P[i]); ok=false; break; }
        fw2.remove(v);
    }

    printf("n              = %u\n",n);
    printf("raw uint32     = %zu B\n",(size_t)n*4);
    printf("log2(n!) floor = %.0f B  (%.2f bits/read)\n",floorBits/8,floorBits/n);
    printf("this coder     = %zu B  (%.2f bits/read)  overhead vs floor %.3f%%\n",
           bytes,bytes*8.0/n,100.0*(bytes-floorBits/8)/(floorBits/8));
    printf("round trip     = %s\n",ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
