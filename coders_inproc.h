#pragma once
// Stage 106: in-process coder wrappers. Same algorithms, same bytes, but
// operating on memory buffers instead of files so the pipeline can run as a
// single process with no intermediate files (PgRC2's architecture -- their
// release build writes zero intermediates; see PGRC2_DISK_ARCHITECTURE.md).
// The RangeEnc/RangeDec model code lives in the original coder sources; the
// pieces needed for encoding are reproduced here in isolated namespaces so
// the standalone binaries keep working unchanged.
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <lzma.h>
#include <map>

// ---------- xz (liblzma), equivalent to `xz -9 -c` ----------
static std::vector<uint8_t> xz_compress(const void* data, size_t n, uint32_t preset=9){
    if(!n) return {};
    size_t cap = lzma_stream_buffer_bound(n) + 128;
    std::vector<uint8_t> out(cap);
    size_t pos = 0;
    lzma_ret r = lzma_easy_buffer_encode(preset, LZMA_CHECK_CRC64, nullptr,
                                         (const uint8_t*)data, n,
                                         out.data(), &pos, cap);
    if(r != LZMA_OK){ fprintf(stderr,"[xz] encode failed rc=%d\n",(int)r); return {}; }
    out.resize(pos);
    return out;
}
template<class T>
static std::vector<uint8_t> xz_compress_vec(const std::vector<T>& v, uint32_t preset=9){
    return xz_compress(v.data(), v.size()*sizeof(T), preset);
}

// ---------- mismatch symbol coder (from 50_mismatch_coder_real.cpp) ----------
namespace mmc {
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

static std::vector<uint8_t> encode(const std::vector<uint8_t>& ref,
                                   const std::vector<uint8_t>& obs){
    const size_t n=ref.size();
    if(!n) return {};
    auto code2=[](uint8_t c)->int{ return c=='A'?0:c=='C'?1:c=='G'?2:3; };
    uint32_t freq[4][3]; for(int r=0;r<4;++r) for(int k=0;k<3;++k) freq[r][k]=1;
    auto excludeIdx=[&](int refc,int obsc)->int{
        int k=0; for(int s=0;s<4;++s){ if(s==refc) continue; if(s==obsc) return k; ++k; } return -1; };
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
    return std::move(enc.out);
}
} // namespace mmc

// ---------- MEM-reference src coder (from 37_ref_coder.cpp) ----------
namespace refc {
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

struct Ref { uint32_t dst,src,len; };
// triples: raw mem_triples.bin bytes, 13 B/record (dst,src,len,is_rc).
static std::vector<uint8_t> encode(const std::vector<uint8_t>& triples,
                                   uint64_t PGLEN, uint64_t MAINEND){
    const size_t n=triples.size()/13;
    if(!n) return {};
    std::vector<Ref> R(n);
    for(size_t i=0;i<n;++i){
        const uint8_t* p=&triples[i*13];
        memcpy(&R[i].dst,p,4); memcpy(&R[i].src,p+4,4); memcpy(&R[i].len,p+8,4);
    }
    std::sort(R.begin(),R.end(),[](const Ref&a,const Ref&b){ return a.dst<b.dst; });
    auto bound=[&](uint32_t dst)->uint32_t{
        return (uint64_t)dst<MAINEND ? (dst?dst:1) : (uint32_t)MAINEND;
    };
    RangeEnc enc; enc.out.reserve(n*4);
    for(const Ref& r:R) enc.encode(r.src,bound(r.dst));
    enc.flush();
    return std::move(enc.out);
}
} // namespace refc

// ---------- mismatch-position bucketing + transpose (PgRC2's technique) ----------
// compressRlMisRevOffDest (SeparatedPseudoGenomePersistence.cpp:823-905) does
// not code mismatch positions as one flat stream. It routes each read's
// positions into a destination chosen by that read's mismatch COUNT, then
// optionally TRANSPOSES each bucket so column k holds the k-th mismatch of
// every read in it. Both moves make the statistics inside a stream far more
// uniform. Measured on real P. aeruginosa: -13.5% at MAXMAP=12, -12.5% at
// MAXMAP=40, on top of the reverse-offset coding already applied.
// Fully invertible from the counts alone, which the archive already carries.
static std::vector<uint8_t> mmpos_bucket_transpose(const std::vector<uint8_t>& pos,
                                                   const std::vector<uint16_t>& counts){
    std::map<uint16_t,std::vector<uint8_t>> buck;
    size_t off=0;
    for(uint16_t c : counts){
        if(c){ auto& b=buck[c]; b.insert(b.end(), pos.begin()+off, pos.begin()+off+c); }
        off += c;
    }
    std::vector<uint8_t> out; out.reserve(pos.size());
    for(auto& kv : buck){
        const uint16_t c=kv.first; const auto& v=kv.second;
        const size_t n=v.size()/c;
        for(uint16_t col=0; col<c; ++col)
            for(size_t r=0;r<n;++r) out.push_back(v[r*c+col]);
    }
    return out;
}
// Inverse -- used by the encoder's own round-trip self-check and by any decoder.
static std::vector<uint8_t> mmpos_untranspose(const std::vector<uint8_t>& t,
                                              const std::vector<uint16_t>& counts){
    std::map<uint16_t,size_t> nreads;
    for(uint16_t c : counts) if(c) ++nreads[c];
    std::map<uint16_t,std::vector<uint8_t>> buck;
    size_t off=0;
    for(auto& kv : nreads){
        const uint16_t c=kv.first; const size_t n=kv.second;
        std::vector<uint8_t> v(n*(size_t)c);
        for(uint16_t col=0; col<c; ++col)
            for(size_t r=0;r<n;++r) v[r*c+col]=t[off++];
        buck[c]=std::move(v);
    }
    std::map<uint16_t,size_t> cur;
    std::vector<uint8_t> out; out.reserve(t.size());
    for(uint16_t c : counts){
        if(!c) continue;
        auto& v=buck[c]; size_t& k=cur[c];
        for(uint16_t j=0;j<c;++j) out.push_back(v[k++]);
    }
    return out;
}

// ---------- mismatch-count zero-flag split (PgRC2's technique) ----------
// They do not code the per-read mismatch counts as one stream: it is split
// into "Mismatches counts (zero flags)" and "Mismatches counts (non-zero
// values)" (visible as two separate streams in their own stderr). Most reads
// have zero mismatches, so a 1-bit flag per read plus a dense value array for
// the rest beats a sparse uint16 array. Measured -10.6% on real
// P. aeruginosa. Values are uint8; the caller must confirm no count exceeds
// 255 (MAXMAP bounds it well below that) and fall back otherwise.
static bool mmcnt_split(const std::vector<uint16_t>& counts,
                        std::vector<uint8_t>& flags, std::vector<uint8_t>& vals){
    flags.clear(); vals.clear();
    flags.reserve((counts.size()+7)/8); vals.reserve(counts.size()/4+1);
    uint8_t acc=0; int nb=0;
    for(uint16_t c : counts){
        if(c>255) return false;
        acc=(uint8_t)((acc<<1)|(c?1:0)); ++nb;
        if(nb==8){ flags.push_back(acc); acc=0; nb=0; }
        if(c) vals.push_back((uint8_t)c);
    }
    if(nb) flags.push_back((uint8_t)(acc<<(8-nb)));
    return true;
}
static std::vector<uint16_t> mmcnt_join(const std::vector<uint8_t>& flags,
                                        const std::vector<uint8_t>& vals, size_t n){
    std::vector<uint16_t> out(n,0);
    size_t k=0;
    for(size_t i=0;i<n;++i){
        const bool nz = (flags[i>>3] >> (7-(i&7))) & 1;
        if(nz) out[i]= vals[k++];
    }
    return out;
}
