// Reference-stream coder: sources bounded by their own destination.
//
// Profiled per component at MINMEM 24, RC-only self-matching:
//
//   component     ours    bits/match      PgRC2   bits/match    delta
//   sources    186,852         25.38    128,671        23.83  +58,181
//   lengths     63,988          8.69     48,509         8.99  +15,479
//   dest gaps   43,932          5.97          0            -  +43,932
//
// Two things that profile settles.
//
// The gap stream is NOT a loss. Their destination positions are in-band: a '%'
// MATCH_MARK sits in the pg literal and rides along inside the VarLenDNA
// codebook, which carries phrases like "T%", "AT%", "GG%". So the fair unit is
// literal plus positions -- ours 2,985,570 + 43,932 = 3,029,502 against their
// 3,056,474, and we are 26,972 AHEAD. An earlier note in this file called the
// gap stream pure overhead they avoid; that was double-counting.
//
// Sources ARE the loss, and the reason is a varint. log2(23,233,953) = 24.47
// bits is the naive floor and we spend 25.38, because a varint for values near
// 23M needs four bytes. Theirs sits at 23.83 -- BELOW their own naive floor --
// which means LZMA is finding structure in it.
//
// The structure is not subtle once looked for: in a self-match the source is
// always earlier than the destination, and matches are stored in destination
// order. So src is uniform in [0, dst), not in [0, pg_len), and coding it that
// way costs log2(dst), which averages about 1.4 bits below log2(pg_len) for
// destinations spread across the pg. That is ~23.0 bits/match, under what they
// manage, and it needs no new information -- the decoder already knows dst when
// it reads src.
//
// Cross-matches (survivor pg into main pg) are the one exception: there the
// source lives in the main pg while the destination is past main_pg_end, so the
// bound is main_pg_end rather than dst. The decoder can tell the two apart from
// dst alone, so nothing extra is stored to distinguish them.
//
// Range coder is the same LZMA-style one used for the permutation in stage 23.
// Round trip is verified before any size is reported.
//
//   g++ -O3 -march=native -o refc 37_ref_coder.cpp
//   ./refc mem_triples.bin <pg_len> <main_pg_end>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

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
    // A 32-bit range coder renormalises at 2^24, so `range` can be as small as
    // 2^24 when this is called. Dividing by a bound near 2^24.5 then yields
    // ZERO, low stops advancing and the renormalise loop pushes bytes forever.
    // Values this wide must be split: a high part and a 12-bit low part, each
    // small enough that range/m stays non-zero.
    void encodeRaw(uint32_t v,uint32_t m){
        if(m<2) return;
        range/=m; low+=(uint64_t)v*range;
        while(range<(1u<<24)){ range<<=8; shiftLow(); }
    }
    void encode(uint32_t v,uint32_t m){          // v uniform in [0,m)
        if(m<2) return;
        const uint32_t K=1u<<12;
        if(m<=K){ encodeRaw(v,m); return; }
        const uint32_t hi=v/K, lo=v%K, hm=(m+K-1)/K;
        encodeRaw(hi,hm);
        encodeRaw(lo,K);
    }
    void flush(){ for(int i=0;i<5;++i) shiftLow(); }
};
struct RangeDec {
    const uint8_t* p; const uint8_t* end; uint32_t range=0xFFFFFFFFu, code=0;
    void init(const uint8_t* b,size_t n){ p=b; end=b+n; ++p;
        for(int i=0;i<4;++i) code=(code<<8)|(p<end?*p++:0); }
    uint32_t decodeRaw(uint32_t m){
        if(m<2) return 0;
        range/=m; uint32_t v=code/range; if(v>=m) v=m-1;
        code-=v*range;
        while(range<(1u<<24)){ range<<=8; code=(code<<8)|(p<end?*p++:0); }
        return v;
    }
    uint32_t decode(uint32_t m){
        if(m<2) return 0;
        const uint32_t K=1u<<12;
        if(m<=K) return decodeRaw(m);
        const uint32_t hm=(m+K-1)/K;
        const uint32_t hi=decodeRaw(hm);
        const uint32_t lo=decodeRaw(K);
        return hi*K+lo;
    }
};

int main(int argc,char** argv){
    if(argc<4){ fprintf(stderr,"usage: %s <mem_triples.bin> <pg_len> <main_pg_end>\n",argv[0]); return 1; }
    const uint64_t PGLEN=strtoull(argv[2],nullptr,10);
    const uint64_t MAINEND=strtoull(argv[3],nullptr,10);

    // STAGE 100 REWRITE: mem_triples.bin now has a 4th field (is_rc, 1 byte)
    // added by the encoder to fix a real reverse-complement bug found while
    // building the first genuine full-pipeline decoder -- 13 bytes/record
    // now, not 12. A raw struct-fread at the old stride silently misaligned
    // every record after the first. Parsed manually (not via a packed
    // struct) to avoid any compiler-padding surprises. is_rc itself isn't
    // part of what this tool compresses (that's a separate, tiny, disclosed
    // stream, 1 bit/match) -- only dst/src/len feed the range coder below,
    // unchanged from before.
    struct Ref { uint32_t dst,src,len; };
    std::vector<Ref> R;
    {
        FILE* f=fopen(argv[1],"rb"); if(!f){ perror("open"); return 1; }
        fseek(f,0,SEEK_END); const size_t sz=ftell(f); fseek(f,0,SEEK_SET);
        const size_t n=sz/13;
        R.resize(n);
        std::vector<uint8_t> buf(sz);
        if(fread(buf.data(),1,sz,f)!=sz){ fprintf(stderr,"short read\n"); fclose(f); return 1; }
        fclose(f);
        // is_rc (byte 12 of each 13-byte record) isn't used by this tool --
        // it doesn't feed the src range coder below -- just skipped over.
        for(size_t i=0;i<n;++i){
            const uint8_t* p=&buf[i*13];
            memcpy(&R[i].dst,p,4); memcpy(&R[i].src,p+4,4); memcpy(&R[i].len,p+8,4);
        }
    }
    std::sort(R.begin(),R.end(),[](const Ref&a,const Ref&b){ return a.dst<b.dst; });
    const size_t n=R.size();

    // The bound the decoder can reconstruct without being told: a destination
    // inside the main pg is a self-match, so the source precedes it; a
    // destination past main_pg_end is a survivor cross-match, whose source lies
    // anywhere in the main pg.
    auto bound=[&](uint32_t dst)->uint32_t{
        return (uint64_t)dst<MAINEND ? (dst?dst:1) : (uint32_t)MAINEND;
    };

    double floorBits=0;
    for(const Ref& r:R) floorBits+=std::log2((double)std::max<uint32_t>(bound(r.dst),2));

    RangeEnc enc; enc.out.reserve(n*4);
    for(const Ref& r:R) enc.encode(r.src,bound(r.dst));
    enc.flush();
    const size_t bytes=enc.out.size();

    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    bool ok=true;
    for(size_t i=0;i<n;++i){
        const uint32_t v=dec.decode(bound(R[i].dst));
        if(v!=R[i].src){ fprintf(stderr,"MISMATCH at %zu: got %u want %u\n",i,v,R[i].src); ok=false; break; }
    }

    printf("matches        %zu\n",n);
    printf("naive floor    %.0f B  (%.2f bits/match, log2(pg_len)=%.2f)\n",
           n*std::log2((double)PGLEN)/8, std::log2((double)PGLEN), std::log2((double)PGLEN));
    printf("bounded floor  %.0f B  (%.2f bits/match)\n",floorBits/8,floorBits/n);
    printf("this coder     %zu B  (%.2f bits/match)  [12-bit split, so slightly above the bounded floor]\n",bytes,bytes*8.0/n);
    printf("varint + xz    186,852 B  (25.38 bits/match)   <- what this replaces\n");
    printf("PgRC2 offsets  128,671 B  (23.83 bits/match)\n");
    printf("round trip     %s\n",ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
