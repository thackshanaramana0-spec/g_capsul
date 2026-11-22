#pragma once
// In-process stream coders and the per-stream framing of the archive.
//
// Every coder here works on memory buffers, so the encoder runs as one process
// with no intermediate files (PgRC2's release build also writes none). Contents:
//   xz_compress, xz_compress_lzma     liblzma wrappers (preset 9, explicit lc/lp/pb)
//   mmc::                             mismatch symbol coder: the observed base given the base it replaced
//   refc::                            MEM source coder (the mem_triples stream) and its self-match flags
//   mmpos_encode_buckets, _decode_    mismatch positions, bucketed by the read's mismatch count
//   mmcnt_split, mmcnt_join           per-read mismatch counts as a zero-flag bitmap plus values
//   best_encode, best_encode_chunked, const_or_encode    per-stream method selection and framing
// The framing (method byte, constant marker, chunk table) is documented at
// capsule_decode_stream() in src/decoder.cpp, which is its inverse. The PPMd7,
// FSE and range coders the selector tries are in coders_entropy.h.
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <lzma.h>
#include <map>
#include <thread>
#include <atomic>
#include "coders_entropy.h"

// ---------- xz (liblzma), equivalent to `xz -9 -c` ----------
// Multi-threaded LZMA. PgRC2 passes a thread count to every LzmaCoderProps
// (PropsLibrary.cpp: noOfThreads = numberOfThreads>1 ? 2 : 1). The
// single-threaded lzma_easy_buffer_encode measured 145% CPU against PgRC2's
// 515% on a 12-core box. Falls back to the
// single-threaded path for small inputs, where MT block-splitting costs bytes
// and buys no time.
static std::vector<uint8_t> xz_compress(const void* data, size_t n, uint32_t preset=9){
    if(!n) return {};
    size_t cap = lzma_stream_buffer_bound(n) + 128;
    std::vector<uint8_t> out(cap);
    size_t pos = 0;
    const size_t MT_MIN = 1u<<21;                 // 2 MB: below this MT is a loss
    if(n >= MT_MIN){
        unsigned hw = std::thread::hardware_concurrency(); if(!hw) hw = 1;
        lzma_mt mt{};
        mt.threads      = hw > 8 ? 8 : hw;
        mt.block_size   = 0;                      // let liblzma choose
        mt.check        = LZMA_CHECK_CRC64;
        mt.preset       = preset;
        lzma_stream strm = LZMA_STREAM_INIT;
        if(lzma_stream_encoder_mt(&strm, &mt) == LZMA_OK){
            strm.next_in = (const uint8_t*)data; strm.avail_in = n;
            strm.next_out = out.data();          strm.avail_out = cap;
            lzma_ret r = lzma_code(&strm, LZMA_FINISH);
            pos = cap - strm.avail_out;
            lzma_end(&strm);
            if(r == LZMA_STREAM_END){ out.resize(pos); return out; }
        }
        pos = 0;                                  // MT failed; fall through
    }
    lzma_ret r = lzma_easy_buffer_encode(preset, LZMA_CHECK_CRC64, nullptr,
                                         (const uint8_t*)data, n,
                                         out.data(), &pos, cap);
    if(r != LZMA_OK){ fprintf(stderr,"[xz] encode failed rc=%d\n",(int)r); return {}; }
    out.resize(pos);
    return out;
}
// xz of a vector's raw bytes. Currently unused.
template<class T>
static std::vector<uint8_t> xz_compress_vec(const std::vector<T>& v, uint32_t preset=9){
    return xz_compress(v.data(), v.size()*sizeof(T), preset);
}

// LZMA with an explicit literal-context/literal-position/position-bit split.
// PgRC2 codes its position stream with lc=8,lp=2,pb=2 rather than the xz
// default lc=3,lp=0,pb=2 (SeparatedPseudoGenomePersistence.cpp:451-460,
// getReadsPositionsCoderProps). lp=2 aligns the context model to a 4-byte
// stride, which is exactly the shape of a uint32 position array. Measured on
// real P. aeruginosa pos_abs: default 5,082,152 -> lc=2,lp=2,pb=2 5,021,990
// (-1.2%).
static std::vector<uint8_t> xz_compress_lzma(const void* data, size_t n,
                                             uint32_t lc, uint32_t lp, uint32_t pb,
                                             uint32_t preset=9){
    if(!n) return {};
    lzma_options_lzma opt;
    if(lzma_lzma_preset(&opt, preset)) return {};
    opt.lc=lc; opt.lp=lp; opt.pb=pb;
    // B2: a dictionary larger than the input is unusable. LZMA matches can only
    // reference data already seen within this buffer, so any dict_size above n
    // reserves memory that no match can ever address. preset 9 asks for 64 MB
    // regardless of stream size, and with seven probes over six large streams
    // that is ~2.7 GB of allocate-touch-free per file -- measured as 1,089,328
    // minor page faults against PgRC2's 169,421 (6.4x), which is essentially
    // all of the encoder's 3.18 s system time.
    //
    // PgRC2 sizes dictionaries per stream by hand (1 MB, 8 MB, 16 MB, 64 MB).
    // The bound here is derived from the data instead: round n up to a power of
    // two, never above what the preset asked for, never below LZMA's minimum.
    // The .xz container records dict_size in the block header, so a decoder
    // adapts on its own, and the match set is unchanged by construction.
    if(opt.dict_size > n){
        uint32_t d = 1u<<12;                        // LZMA_DICT_SIZE_MIN
        while(d < n && d < opt.dict_size) d <<= 1;
        opt.dict_size = d;
    }
    lzma_filter filters[2];
    filters[0].id=LZMA_FILTER_LZMA2; filters[0].options=&opt;
    filters[1].id=LZMA_VLI_UNKNOWN;  filters[1].options=nullptr;
    size_t cap=lzma_stream_buffer_bound(n)+128;
    std::vector<uint8_t> out(cap); size_t pos=0;
    if(lzma_stream_buffer_encode(filters, LZMA_CHECK_CRC64, nullptr,
                                 (const uint8_t*)data, n, out.data(), &pos, cap) != LZMA_OK)
        return {};
    out.resize(pos);
    return out;
}
// Byte-plane (structure-of-arrays) split for a uint32 stream: emit every
// value's byte 0, then every byte 1, and so on. Each plane then has uniform
// statistics -- the high plane is highly skewed and compresses hard, the low
// planes are near-random. This beats stride-aligned lp tuning on real data
// (P. aeruginosa pos_abs: 4,983,604 vs 5,021,990 vs 5,082,152 for xz default).
// ── B1: do not spend coding effort on provably incompressible byte planes ──
//
// PgRC2 does not probe its position stream at all: getReadsPositionsCoderProps
// returns one LZMA coder whose context stride is DERIVED from the width the
// values need (SimplePgMatcher.cpp:233 picks DATAPERIOD 32 vs 64 from whether
// the pg length is standard). The coder is a consequence of what the data is.
//
// Our positions are indices into a multi-megabyte pseudogenome, so their low
// bytes are uniform noise by construction. Measured on pos_abs: planes 0 and 1
// compress to ratio exactly 1.000 on every file -- they come out LARGER than
// raw -- and they are half the stream. Shannon says a plane whose order-0
// entropy is 8.00 bits/byte cannot be coded below its raw size, so the decision
// needs no trial: compute the histogram once, in O(n), and store such planes
// verbatim.
//
// This is conditional on the data, not a fixed rule: orig2uid's low planes have
// entropy 2.36, and storing THOSE raw costs +2,141,527 B. The criterion selects
// correctly in both cases.
//
// The guard: a plane can have H0 = 8 and still be compressible if it has
// higher-order structure (a counter mod 256 is the textbook case). So a plane
// that looks incompressible by entropy is confirmed on a bounded 64 KB prefix
// -- PgRC2's own DEFAULT_MIN_PROBE_SIZE -- before it is skipped. That is O(64 KB),
// not O(n).
static double plane_entropy(const uint8_t* p, size_t n){
    if(!n) return 0.0;
    size_t c[256]={0};
    for(size_t i=0;i<n;++i) ++c[p[i]];
    double h=0.0, inv=1.0/(double)n;
    for(int k=0;k<256;++k) if(c[k]){ double q=c[k]*inv; h-=q*std::log2(q); }
    return h;
}

// Transposes a uint32 array into byte planes: out[k*cnt+i] is byte k of value i,
// with cnt = n/4. Any trailing n%4 bytes are dropped, so callers use it only when
// n is a multiple of 4 (best_encode checks this before offering methods 5 and 6).
static std::vector<uint8_t> u32_byteplanes(const uint8_t* d, size_t n){
    const size_t cnt=n/4;
    std::vector<uint8_t> out(cnt*4);
    for(size_t i=0;i<cnt;++i)
        for(size_t k=0;k<4;++k) out[k*cnt+i]=d[i*4+k];
    return out;
}
// Try every candidate encoding and keep the smallest, recording the choice in
// a 1-byte header so the decoder can invert it. "Measure, do not assume" --
// the same discipline used for every other transform here. PgRC2 does the
// equivalent with its selector coders.
//   0 = xz default, 1 = lzma lc2/lp2/pb2, 2 = byte-planes + xz,
//   3 = byte-planes + lzma lc0/lp0/pb0
// Currently unused. best_encode below covers the same candidates and is what the
// archive uses. These method ids are NOT the ones capsule_decode_stream() reads,
// so a stream coded here could not be decoded by it.
static std::vector<uint8_t> xz_best_for_u32(const void* data, size_t n){
    std::vector<std::vector<uint8_t>> cand(4);
    cand[0] = xz_compress(data,n);
    cand[1] = xz_compress_lzma(data,n,2,2,2);
    if(n>=4 && n%4==0){
        auto bp = u32_byteplanes((const uint8_t*)data,n);
        cand[2] = xz_compress(bp.data(),bp.size());
        cand[3] = xz_compress_lzma(bp.data(),bp.size(),0,0,0);
    }
    int best=0; size_t bs=SIZE_MAX;
    for(int i=0;i<4;++i) if(!cand[i].empty() && cand[i].size()<bs){ bs=cand[i].size(); best=i; }
    std::vector<uint8_t> out; out.reserve(bs+1);
    out.push_back((uint8_t)best);
    out.insert(out.end(), cand[best].begin(), cand[best].end());
    return out;
}

// ---------- mismatch symbol coder ----------
// Codes the base OBSERVED in a read at each deviation, given the pseudogenome base
// it replaced (`ref`). Because observed != ref, only three symbols are possible, so
// each symbol is coded among those three with an adaptive frequency table kept
// separately for each of the four ref bases (every frequency starts at 1, each
// coded symbol adds 8, and the table is halved, keeping at least 1, when the total
// would pass 65536). The decoder needs the identical ref sequence in the identical
// order for its tables to track. ref is never stored: it is pg[read position +
// offset], recomputed from the rebuilt pseudogenome. Coded form: raw range-coder
// bytes, no header; the first byte is the initial carry (zero) and the decoder
// skips it.
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

// Inverse. obs was coded against ref under an exclusion model (obs != ref, so
// three symbols), and the model adapts per ref base -- so the decoder must see
// the SAME ref sequence in the SAME order for its frequency tables to track.
// ref is never stored: it is pg[read_pos + mismatch_offset], which the decoder
// recomputes from the reconstructed pseudogenome.
static std::vector<uint8_t> decode(const uint8_t* d, size_t n,
                                   const std::vector<uint8_t>& ref){
    std::vector<uint8_t> obs(ref.size(), 0);
    if(!n || ref.empty()) return obs;
    auto code2=[](uint8_t c)->int{ return c=='A'?0:c=='C'?1:c=='G'?2:3; };
    const char SYM[4]={'A','C','G','T'};
    uint32_t freq[4][3]; for(int r=0;r<4;++r) for(int k=0;k<3;++k) freq[r][k]=1;
    RangeDec dec; dec.init(d,n);
    for(size_t i=0;i<ref.size();++i){
        const int r=code2(ref[i]);
        uint32_t* f=freq[r]; const uint32_t tot=f[0]+f[1]+f[2];
        const uint32_t target=dec.getFreq(tot);
        uint32_t lo=0; int k=0;
        for(; k<3; ++k){ if(lo+f[k]>target) break; lo+=f[k]; }
        if(k==3) k=2;
        dec.decodeUpdate(lo, lo+f[k]);
        int seen=0, sym=0;
        for(int t=0;t<4;++t){ if(t==r) continue; if(seen==k){ sym=t; break; } ++seen; }
        obs[i]=(uint8_t)SYM[sym];
        f[k]+=8; if(tot+8>65536){ for(int j=0;j<3;++j) f[j]=(f[j]>>1)|1; }
    }
    return obs;
}

// Streaming form of the same decoder, one symbol per call. Needed whenever
// `ref` for entry i+1 can depend on the OVERRIDE applied for entry i -- e.g.
// MEM matches chained through a tandem repeat, where a later match's SOURCE
// falls inside an earlier match's DESTINATION. Applying overrides only after
// a full batch decode (as `decode` above requires, since it needs the whole
// `ref` array up front) reads the earlier match's un-corrected byte in that
// case. This keeps the identical adaptive state as `decode` -- same freq
// tables, same update rule -- just exposed one call at a time so the caller
// can apply each override before computing the next entry's ref.
struct StreamDecoder {
    RangeDec dec; uint32_t freq[4][3];
    void init(const uint8_t* d, size_t n){
        dec.init(d,n);
        for(int r=0;r<4;++r) for(int k=0;k<3;++k) freq[r][k]=1;
    }
    uint8_t next(uint8_t refc){
        auto code2=[](uint8_t c)->int{ return c=='A'?0:c=='C'?1:c=='G'?2:3; };
        const char SYM[4]={'A','C','G','T'};
        const int r=code2(refc);
        uint32_t* f=freq[r]; const uint32_t tot=f[0]+f[1]+f[2];
        const uint32_t target=dec.getFreq(tot);
        uint32_t lo=0; int k=0;
        for(; k<3; ++k){ if(lo+f[k]>target) break; lo+=f[k]; }
        if(k==3) k=2;
        dec.decodeUpdate(lo, lo+f[k]);
        int seen=0, sym=0;
        for(int t=0;t<4;++t){ if(t==r) continue; if(seen==k){ sym=t; break; } ++seen; }
        f[k]+=8; if(tot+8>65536){ for(int j=0;j<3;++j) f[j]=(f[j]>>1)|1; }
        return (uint8_t)SYM[sym];
    }
};
} // namespace mmc

// ---------- MEM-reference src coder ----------
namespace refc {
// Coder for the SOURCE position of each MEM reference (the mem_triples stream).
// A reference is (dst, src, len, rc). dst, len and rc travel in their own streams
// (mem_dstgap, mem_len, mem_rc), so this stream holds only src. References are
// coded in destination order, each src uniform within a bound derived from its
// own dst, which the decoder already knows when it reads src:
//   self-match inside the main pg:           src is earlier than dst, so src in [0, dst)
//   cross-match (source in the main pg, destination past main_pg_end):
//                                            src in [0, main_pg_end)
//   second-region self-match:                src coded relative to main_pg_end,
//                                            flagged in the mem_self stream
// Coding src in [0, dst) instead of [0, pg_len) costs log2(dst) bits, on average
// about 1.4 bits fewer per match than log2(pg_len) when destinations are spread
// across the pg. The measured per-component profile against PgRC2 (sources,
// lengths, destination gaps, and why the gap stream is not a loss) is in
// docs/extras/ENGINEERING_NOTES_CODERS.md.

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
// A reference whose source sits past main_pg_end is a second-region self match;
// its source is coded relative to main_pg_end so it stays inside the model.
static inline uint32_t bound_sel(uint32_t dst, bool self, uint64_t MAINEND){
    if(self){ uint32_t b=(uint32_t)((uint64_t)dst-MAINEND); return b?b:1; }
    return (uint64_t)dst<MAINEND ? (dst?dst:1) : (uint32_t)MAINEND;
}
// bound_sel() with the self flag derived from (dst, src). Currently unused:
// encode() and decode() call bound_sel() with the flag they already hold.
static inline uint32_t bound_for(uint32_t dst, uint32_t src, uint64_t MAINEND){
    return bound_sel(dst, (uint64_t)src>=MAINEND && (uint64_t)dst>=MAINEND, MAINEND);
}
// One flag per reference, in the coder's own dst-sorted order. All zero whenever
// no self pass ran, so it costs nothing on every existing dataset.
static std::vector<uint8_t> encode_self(const std::vector<uint8_t>& triples, uint64_t MAINEND){
    const size_t n=triples.size()/13;
    std::vector<Ref> R(n);
    for(size_t i=0;i<n;++i){ const uint8_t* p=&triples[i*13];
        memcpy(&R[i].dst,p,4); memcpy(&R[i].src,p+4,4); memcpy(&R[i].len,p+8,4); }
    std::sort(R.begin(),R.end(),[](const Ref&a,const Ref&b){
        if(a.dst!=b.dst) return a.dst<b.dst;
        if(a.src!=b.src) return a.src<b.src;
        return a.len<b.len; });   // total order: std::sort is not stable, and the
                                  // flag pass and the source pass must agree exactly
    std::vector<uint8_t> f(n,0); bool any=false;
    for(size_t i=0;i<n;++i)
        if((uint64_t)R[i].src>=MAINEND && (uint64_t)R[i].dst>=MAINEND){ f[i]=1; any=true; }
    if(!any) return {};            // absent stream == no self references
    return f;
}
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
    std::sort(R.begin(),R.end(),[](const Ref&a,const Ref&b){
        if(a.dst!=b.dst) return a.dst<b.dst;
        if(a.src!=b.src) return a.src<b.src;
        return a.len<b.len; });   // total order: std::sort is not stable, and the
                                  // flag pass and the source pass must agree exactly
    // The rule above -- "a destination past main_pg_end implies a source inside
    // the main pg, so nothing extra need be stored" -- is an assumption about
    // WHICH passes exist, baked into the coder. A second-region self-match
    // breaks it: such a reference has dst >= main_pg_end AND src >= main_pg_end,
    // so its source is outside the model's range entirely and cannot be
    // represented. The reference then decodes to a legal-looking but wrong
    // position, which is why the destinations, the ordering invariant and the
    // literal accounting can all be exact while the recovered bases are wrong.
    //
    // A self reference is coded relative to main_pg_end instead, so its source
    // lands in [0, dst-main_pg_end). Whether a reference is one cannot be
    // derived from dst alone any more, so it is stored -- see encode_self().
    RangeEnc enc; enc.out.reserve(n*4);
    for(const Ref& r:R){
        const bool self=((uint64_t)r.src>=MAINEND && (uint64_t)r.dst>=MAINEND);
        const uint32_t v = self ? (uint32_t)((uint64_t)r.src-MAINEND) : r.src;
        enc.encode(v, bound_sel(r.dst,self,MAINEND));
    }
    enc.flush();
    return std::move(enc.out);
}

// Inverse. Each source was coded with a modulus derived from its OWN
// destination, so the destinations must be known first -- they come from the
// mem_dstgap stream, which is why that stream had to exist before this could.
// RangeDec above was already written and simply never called.
static std::vector<uint32_t> decode(const uint8_t* d, size_t n,
                                    const std::vector<uint32_t>& dst,
                                    uint64_t MAINEND,
                                    const std::vector<uint8_t>& selfflags={}){
    std::vector<uint32_t> src(dst.size(), 0);
    if(!n || dst.empty()) return src;
    RangeDec dec; dec.init(d,n);
    for(size_t i=0;i<dst.size();++i){
        const bool self = (i<selfflags.size() && selfflags[i]);
        const uint32_t v = dec.decode(bound_sel(dst[i],self,MAINEND));
        src[i] = self ? (uint32_t)(v+MAINEND) : v;
    }
    return src;
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
// mmpos_bucket_transpose and mmpos_untranspose, which flatten all buckets into one
// transposed stream, are currently unused. mmpos_encode_buckets (below) does the
// bucketing itself, decides per bucket whether to transpose, and frames the result.
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
// Inverse of mmpos_bucket_transpose. Currently unused (see above).
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
// Inverse of mmcnt_split: rebuilds n counts from the flag bitmap and the dense
// values (a flagged read takes the next value, an unflagged read is 0).
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

// ---------- universal per-stream selector over the REAL coder set ----------
// Tries the real coder set PgRC2's selector uses (PPMd, FSE/Huff0, range coders
// with per-position models, PropsLibrary.cpp getSelectorCoderProps) plus xz, and
// keeps the smallest, with a 1-byte method id so the result stays decodable.
//   0 xz default | 1 xz lc2lp2pb2 | 2 ppmd5 | 3 fse | 4 range(period=1)
//   5 byteplanes+xz | 6 byteplanes+lzma(0,0,0)
// Method 7 (chunked) is added by best_encode_chunked, and the 0xC0 marker for a
// constant stream by const_or_encode.
static std::vector<uint8_t> encode_method(const uint8_t* d, size_t n, int m){
    switch(m){
        case 0: return xz_compress(d,n);
        case 1: return xz_compress_lzma(d,n,2,2,2);
        case 2: return pgc::ppmd_encode(d,n,5,32);
        case 3: return pgc::fse_encode(d,n);
        case 4: return pgc::range_encode(d,n,1);
        case 5: { auto bp=u32_byteplanes(d,n); return xz_compress(bp.data(),bp.size()); }
        case 6: {
            // Byte planes coded together (one shared LZMA model), except that
            // any plane the data proves incompressible is stored verbatim.
            const size_t cnt=n/4; if(!cnt) return {};
            auto bp=u32_byteplanes(d,n);
            uint8_t rawmask=0;
            for(int k=0;k<4;++k){
                const uint8_t* pl=bp.data()+(size_t)k*cnt;
                if(plane_entropy(pl,cnt) < 7.999) continue;      // codeable
                // entropy says incompressible; confirm on a bounded prefix so a
                // high-entropy-but-structured plane is never thrown away
                const size_t PREF_CAP = 1u<<16;   // PgRC2 DEFAULT_MIN_PROBE_SIZE
                const size_t pref = cnt < PREF_CAP ? cnt : PREF_CAP;
                auto t = xz_compress_lzma(pl,pref,0,0,0);
                if(!t.empty() && t.size() < pref) continue;      // it IS codeable
                rawmask |= (uint8_t)(1u<<k);
            }
            if(!rawmask) {                                        // nothing to skip
                auto c=xz_compress_lzma(bp.data(),bp.size(),0,0,0);
                std::vector<uint8_t> o; o.reserve(c.size()+1);
                o.push_back(0); o.insert(o.end(),c.begin(),c.end());
                return o;
            }
            std::vector<uint8_t> keep; keep.reserve(bp.size());
            for(int k=0;k<4;++k) if(!(rawmask&(1u<<k)))
                keep.insert(keep.end(), bp.begin()+(size_t)k*cnt, bp.begin()+(size_t)(k+1)*cnt);
            auto c = keep.empty() ? std::vector<uint8_t>() 
                                  : xz_compress_lzma(keep.data(),keep.size(),0,0,0);
            std::vector<uint8_t> o; o.reserve(bp.size()+16);
            o.push_back(rawmask);
            uint32_t cl=(uint32_t)c.size();
            o.insert(o.end(),(uint8_t*)&cl,(uint8_t*)&cl+4);
            o.insert(o.end(),c.begin(),c.end());
            for(int k=0;k<4;++k) if(rawmask&(1u<<k))
                o.insert(o.end(), bp.begin()+(size_t)k*cnt, bp.begin()+(size_t)(k+1)*cnt);
            return o;
        }
    }
    return {};
}
// PgRC2 does NOT run every candidate over the whole stream -- SelectorCoderProps
// takes a probeFraction with a 64 KB floor (CodersLib.h:216-236,
// getSelectorCoderProps(..., 0.2, ...)), picks a winner on the probe, then
// encodes the full stream ONCE. Running all candidates at full length was
// costing 10.0 s of an 18.7 s run, 6.69 s of it on pos_abs alone.
static const size_t PROBE_MIN   = 1u<<16;   // 64 KB, their DEFAULT_MIN_PROBE_SIZE
// Fraction of a stream used to choose its method (override: PROBE_FRAC).
static double probe_frac(){
    const char* e = getenv("PROBE_FRAC");
    // PgRC2 uses 0.20. Measured on real files: 0.20 costs E. coli +50,149 B
    // (+0.63%) because a coder that wins on a 20% sample is not always the
    // winner on the full stream, while 0.40 recovers all but 60 B of full
    // selection at no measurable time cost, and makes no difference at all on
    // P. aeruginosa. Raised to 0.40 on that evidence.
    return e ? atof(e) : 0.40;
}

// Codes d[0..n) with the smallest candidate method and returns the framed result
// [method:1][raw_len:8][payload] (empty when n is 0). The method is chosen on a
// prefix of probe_frac() of the stream, at least PROBE_MIN bytes, and the whole
// stream is then coded once with the winner. A stream no longer than the probe is
// coded by every method in full and the smallest kept. u32shaped declares the data
// to be a uint32 array, which also offers the byte-plane methods 5 and 6 when n is
// a multiple of 4. If the chosen coder returns nothing, the stream falls back to xz.
static std::vector<uint8_t> best_encode(const uint8_t* d, size_t n, bool u32shaped=false){
    if(!n) return {};
    std::vector<int> methods = {0,1,2,3,4};
    if(u32shaped && n>=4 && n%4==0){ methods.push_back(5); methods.push_back(6); }

    size_t probe = (size_t)(n*probe_frac());
    if(probe < PROBE_MIN) probe = PROBE_MIN;
    if(probe > n) probe = n;
    if(u32shaped) probe &= ~(size_t)3;          // keep the 4-byte stride intact
    if(probe < 4) probe = n;

    int best = 0;
    if(probe < n){
        size_t bs = SIZE_MAX;
        for(int m : methods){
            auto c = encode_method(d, probe, m);
            if(!c.empty() && c.size() < bs){ bs = c.size(); best = m; }
        }
    } else {
        // stream is small enough that probing saves nothing -- decide exactly
        size_t bs = SIZE_MAX;
        std::vector<uint8_t> keep; int km = 0;
        for(int m : methods){
            auto c = encode_method(d, n, m);
            if(!c.empty() && c.size() < bs){ bs = c.size(); keep = std::move(c); km = m; }
        }
        // [method:1][raw_len:8][payload]. FSE, PPMd and the range coder cannot
        // recover the uncompressed length from their own output, so a decoder
        // has no way back without it. LZMA's xz container carries it, but the
        // field is written uniformly so one decode path serves every method.
        std::vector<uint8_t> out; out.reserve(bs+9);
        out.push_back((uint8_t)km);
        uint64_t rl=n; out.insert(out.end(),(uint8_t*)&rl,(uint8_t*)&rl+8);
        out.insert(out.end(), keep.begin(), keep.end());
        return out;
    }
    if(getenv("LOG_PICK"))
        fprintf(stderr,"  [pick] n=%-10zu u32=%d -> method %d\n", n, (int)u32shaped, best);
    auto coded = encode_method(d, n, best);
    if(coded.empty()) coded = xz_compress(d,n), best = 0;
    std::vector<uint8_t> out; out.reserve(coded.size()+9);
    out.push_back((uint8_t)best);
    uint64_t rl=n; out.insert(out.end(),(uint8_t*)&rl,(uint8_t*)&rl+8);
    out.insert(out.end(), coded.begin(), coded.end());
    return out;
}

// Code a stream as independent chunks, in parallel.
//
// Measured on all 7 datasets: the coder pool's wall clock equals its single
// longest job to within 0.06 s, so every other stream hides behind the biggest
// one and adding cores cannot help. The fix is to split WITHIN the largest
// stream, which is what `literal` already does.
//
// Method 7 nests: each chunk is itself a complete best_encode output, so the
// decoder needs no new per-method logic, only the chunk table.
//
//   [7][total_raw:8][nchunks:u32][coded_len:u32 x nchunks][chunk 0 payload][chunk 1 payload]...
//
// The split factor is DERIVED from the stream size, never fixed: chunking is a
// loss on small streams (H. salinarum pos_abs is slower chunked, 0.67 -> 0.82 s
// sequential) and the context reset costs real bytes. Below CHUNK_MIN the whole
// stream is coded as one piece and this path is not used at all.
static const size_t CHUNK_MIN   = 4u<<20;   // below this, splitting loses
static const size_t CHUNK_TARGET= 2u<<20;   // aim for chunks about this size

// Number of chunks for an n-byte stream: 1 below CHUNK_MIN, otherwise about
// n / CHUNK_TARGET, at least 2 and at most the hardware thread count.
static size_t chunk_count(size_t n){
    if(n < CHUNK_MIN) return 1;
    unsigned hw = std::thread::hardware_concurrency(); if(!hw) hw = 4;
    size_t k = n / CHUNK_TARGET;
    if(k > hw) k = hw;
    if(k < 2)  k = 2;
    return k;
}

static std::vector<uint8_t> best_encode_chunked(const uint8_t* d, size_t n,
                                                bool u32shaped=false){
    const size_t k = chunk_count(n);
    if(k < 2) return best_encode(d, n, u32shaped);
    size_t per = n / k;
    if(u32shaped) per &= ~(size_t)3;            // keep the 4-byte stride intact
    if(!per) return best_encode(d, n, u32shaped);

    std::vector<std::pair<size_t,size_t>> spans;
    for(size_t c=0;c<k;++c){
        size_t off = c*per, len = (c==k-1) ? (n-off) : per;
        if(off>=n) break;
        spans.push_back({off,len});
    }
    std::vector<std::vector<uint8_t>> parts(spans.size());
    {
        std::atomic<size_t> next{0};
        // Concurrency is bounded separately from the split. Every concurrent
        // chunk holds its own LZMA encoder state (~11.5x its dictionary), so
        // running all chunks at once is what costs RAM -- and the speed gain
        // stops as soon as this job is no longer the pool's bottleneck, while
        // the memory cost keeps growing. Splitting also helps on its own
        // (LZMA is superlinear in block size) at no memory cost at all.
        // Measured on L. major: 1 -> 41.69s/941MB, 2 -> 38.31s/989MB,
        // 4 -> 35.56s/1118MB, 12 -> 36.30s/1441MB. Twelve is both SLOWER and
        // 49% heavier than four, so the gain saturates well before the core
        // count while memory keeps climbing. Two takes most of the speed for
        // 2.5% memory; four costs 16% for the last few percent.
        unsigned CPAR = 2;
        if(const char* e=getenv("CHUNK_PAR")){ unsigned v=atoi(e); if(v) CPAR=v; }
        unsigned T = (unsigned)std::min<size_t>(spans.size(),
                        std::min<size_t>(CPAR, std::max(1u, std::thread::hardware_concurrency())));
        std::vector<std::thread> th;
        for(unsigned t=0;t<T;++t) th.emplace_back([&]{
            for(;;){ size_t i=next.fetch_add(1); if(i>=spans.size()) break;
                parts[i] = best_encode(d+spans[i].first, spans[i].second, u32shaped); }
        });
        for(auto& x:th) x.join();
    }
    size_t tot=0; for(auto& p:parts){ if(p.empty()) return best_encode(d,n,u32shaped); tot+=p.size(); }

    std::vector<uint8_t> out; out.reserve(tot + 13 + 4*parts.size());
    out.push_back((uint8_t)7);
    uint64_t rl=n; out.insert(out.end(),(uint8_t*)&rl,(uint8_t*)&rl+8);
    uint32_t nc=(uint32_t)parts.size(); out.insert(out.end(),(uint8_t*)&nc,(uint8_t*)&nc+4);
    for(auto& p:parts){ uint32_t L=(uint32_t)p.size();
        out.insert(out.end(),(uint8_t*)&L,(uint8_t*)&L+4); }
    for(auto& p:parts) out.insert(out.end(), p.begin(), p.end());

    // Splitting costs bytes (the coder context resets at each chunk) and this
    // function does not compare against coding the stream whole. The trade is set
    // once, by CHUNK_MIN: streams below it are never split, streams above it always are.
    return out;
}

// A stream with exactly ONE distinct element carries no information beyond that
// element and how many times it repeats. Measured: read_lengths on E. coli is
// 3,106,518 B holding the value 150 for all 1,553,259 reads -- one distinct
// value, order-0 entropy exactly 0. Its exact encoding is (count, value).
//
// This is detected, never assumed: variable-length input (SARS-CoV-2, <=221)
// takes the normal path unchanged. The marker 0xC0 cannot collide with
// best_encode's leading byte, which is a method id in 0..7, so the ordinary
// path pays nothing.
static const uint8_t CONST_MARKER = 0xC0;
static std::vector<uint8_t> const_or_encode(const uint8_t* d, size_t n, size_t w,
                                            bool u32shaped=false){
    if(n && w && n % w == 0 && n > w){
        bool same = true;
        for(size_t i=w;i<n;i+=w) if(memcmp(d, d+i, w)){ same=false; break; }
        if(same){
            std::vector<uint8_t> o; o.reserve(1+8+w);
            o.push_back(CONST_MARKER);
            uint64_t cnt = n / w;
            o.insert(o.end(), (const uint8_t*)&cnt, (const uint8_t*)&cnt + 8);
            o.insert(o.end(), d, d + w);
            return o;
        }
    }
    return best_encode(d, n, u32shaped);
}

// mm_pos: bucket by mismatch count, then code each bucket with the best real
// coder -- including the range coder at period = that bucket's count, which is
// exactly PgRC2's compressRlMisRevOffDest scheme (their per-bucket streams are
// the "period = N" lines in its stderr).
static std::vector<uint8_t> mmpos_encode_buckets(const std::vector<uint8_t>& pos,
                                                 const std::vector<uint16_t>& counts){
    if(pos.empty()) return {};
    std::map<uint16_t,std::vector<uint8_t>> b;
    size_t off=0;
    for(uint16_t c : counts){
        if(c){ auto& v=b[c]; v.insert(v.end(), pos.begin()+off, pos.begin()+off+c); }
        off += c;
    }
    // Frame: [nbuckets:u32], then per bucket
    //   [key:u16 (the mismatch count)][raw:u32][method:u8][transposed:u8][coded:u32][payload]
    // 12 B of header per bucket. The method ids are those of encode_method,
    // restricted to 0, 2, 3 and 4, and for method 4 the range coder's period is
    // the key. `transposed` says the bucket is stored column-major (see below).
    // The key and raw length could be derived from the per-read mismatch counts
    // but are stored so the stream can be decoded on its own. An earlier form wrote
    // the buckets back to back with no count, key or length, so nothing marked where
    // one ended and the stream could not be decoded at all.
    std::vector<uint8_t> out;
    auto put32=[&](uint32_t v){ out.insert(out.end(),(uint8_t*)&v,(uint8_t*)&v+4); };
    auto put16=[&](uint16_t v){ out.insert(out.end(),(uint8_t*)&v,(uint8_t*)&v+2); };
    put32((uint32_t)b.size());
    // Within a bucket every read has the SAME mismatch count, so the k-th
    // mismatch of each read forms a column with far more uniform statistics
    // than the row-major order: early columns hold small offsets, late columns
    // large ones. PgRC2 transposes for exactly this reason
    // (compressRlMisRevOffDest). Measured on S. acidocaldarius, the only
    // losing dataset, where mm_pos is 23.4% of the archive: 751,726 B in read
    // order against 710,820 B transposed, -5.4% by entropy.
    //
    // Chosen per bucket by coding both and keeping the smaller, so a bucket
    // that does not benefit is unaffected -- the flag costs one byte.
    // Buckets are independent, so they are coded concurrently (this loop was the
    // longest job in the encoder's coding pool) and emitted afterwards in the
    // original std::map order, which keeps the output byte-identical to serial
    // coding. Concurrency is deliberately modest because each in-flight bucket
    // holds an LZMA encoder state. Measurements: docs/extras/ENGINEERING_NOTES_CODERS.md.
    struct BRec { uint16_t c; uint32_t vsz; uint8_t bi, bt; std::vector<uint8_t> keep; };
    std::vector<BRec> recs(b.size());
    {
        std::vector<std::pair<const uint16_t, std::vector<uint8_t>>*> bv;
        bv.reserve(b.size());
        for(auto& kv : b) bv.push_back(&kv);
        unsigned hw = std::thread::hardware_concurrency(); if(!hw) hw = 4;
        unsigned nt = (unsigned)std::min<size_t>(bv.size(), std::max(1u, hw/2));
        std::atomic<size_t> next{0};
        auto work=[&]{
            for(;;){
                const size_t i = next.fetch_add(1);
                if(i >= bv.size()) break;
                const uint16_t c = bv[i]->first;
                std::vector<uint8_t>& v = bv[i]->second;
                std::vector<uint8_t> vt; vt.reserve(v.size());
                { const size_t rows=v.size()/c;
                  for(uint16_t col=0; col<c; ++col)
                      for(size_t r=0;r<rows;++r) vt.push_back(v[r*(size_t)c+col]); }
                int bi=0, bt=0; size_t bs=SIZE_MAX; std::vector<uint8_t> keep;
                for(int t=0;t<2;++t){
                    const std::vector<uint8_t>& src = t ? vt : v;
                    if(src.size()!=v.size()) continue;
                    std::vector<std::pair<int,std::vector<uint8_t>>> cand;
                    cand.emplace_back(0, xz_compress(src.data(),src.size()));
                    cand.emplace_back(2, pgc::ppmd_encode(src.data(),src.size(),5,32));
                    cand.emplace_back(3, pgc::fse_encode(src.data(),src.size()));
                    cand.emplace_back(4, pgc::range_encode(src.data(),src.size(),c));
                    for(auto& p:cand) if(!p.second.empty() && p.second.size()<bs){
                        bs=p.second.size(); bi=p.first; bt=t; keep=std::move(p.second); }
                }
                recs[i] = BRec{ c, (uint32_t)v.size(), (uint8_t)bi, (uint8_t)bt, std::move(keep) };
            }
        };
        if(nt <= 1){ work(); }
        else {
            std::vector<std::thread> th;
            for(unsigned t=0;t<nt;++t) th.emplace_back(work);
            for(auto& x:th) x.join();
        }
    }
    for(const BRec& r : recs){
        put16(r.c); put32(r.vsz);
        out.push_back(r.bi); out.push_back(r.bt);
        put32((uint32_t)r.keep.size());
        out.insert(out.end(), r.keep.begin(), r.keep.end());
    }
    return out;
}

// Inverse: rebuild the flat position stream from the framed buckets. The
// buckets hold positions grouped by mismatch count; `counts` says which read
// contributed how many, so the flat order is restored by walking the reads and
// drawing from the matching bucket in turn.
static std::vector<uint8_t> mmpos_decode_buckets(const uint8_t* d, size_t n,
                                                 const std::vector<uint16_t>& counts){
    std::vector<uint8_t> flat;
    if(n<4) return flat;
    uint32_t nb; memcpy(&nb,d,4); size_t off=4;
    std::map<uint16_t,std::vector<uint8_t>> b;
    for(uint32_t i=0;i<nb;++i){
        if(off+12>n) return {};
        uint16_t key; memcpy(&key,d+off,2); off+=2;
        uint32_t raw; memcpy(&raw,d+off,4); off+=4;
        const uint8_t meth=d[off]; off+=1;
        const uint8_t tflag=d[off]; off+=1;
        uint32_t cl; memcpy(&cl,d+off,4); off+=4;
        if(off+cl>n) return {};
        std::vector<uint8_t> v;
        switch(meth){
            case 0: { std::vector<uint8_t> o(raw); size_t ip=0,op=0; uint64_t ml=UINT64_MAX;
                      if(lzma_stream_buffer_decode(&ml,0,nullptr,d+off,&ip,cl,o.data(),&op,o.size())!=LZMA_OK) return {};
                      o.resize(op); v=std::move(o); } break;
            case 2: v = pgc::ppmd_decode(d+off,cl,raw); break;
            case 3: v = pgc::fse_decode(d+off,cl,raw); break;
            case 4: v = pgc::range_decode(d+off,cl,raw,key); break;
            default: return {};
        }
        if(v.size()!=raw) return {};
        if(tflag){                                  // column-major -> row-major
            std::vector<uint8_t> u(v.size());
            const size_t rows=v.size()/key; size_t k=0;
            for(uint16_t col=0; col<key; ++col)
                for(size_t r=0;r<rows;++r) u[r*(size_t)key+col]=v[k++];
            v.swap(u);
        }
        b[key]=std::move(v); off+=cl;
    }
    std::map<uint16_t,size_t> cur;
    for(uint16_t c : counts){
        if(!c) continue;
        auto it=b.find(c); if(it==b.end()) return {};
        size_t& k=cur[c];
        if(k+c>it->second.size()) return {};
        flat.insert(flat.end(), it->second.begin()+k, it->second.begin()+k+c);
        k+=c;
    }
    return flat;
}
