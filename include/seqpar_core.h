#pragma once
// DNA context-mixing coder. Codes a sequence of 2-bit base codes (A=0, C=1,
// G=2, T=3) and decodes it again. The encoder uses seq_encode_mem and the
// decoder seq_decode_mem, both defined here so the model cannot drift between
// the two. It codes the `literal` stream: the pseudogenome bases that no MEM
// reference covers.
//
// Model, rebuilt from scratch for each chunk (see run_chunk):
//   * Each base is two binary decisions: the high bit (node 0), then the low bit
//     given the high bit (node 1 or 2). Probabilities are 12-bit.
//   * Ten context models, one per order in ORD (the previous 1, 2, 3, 4, 6, 8,
//     11, 14, 18 and 22 bases). Each hashes its context and the node into a table
//     of 2^TBITS 16-bit probabilities. The learning rate starts at 1/2 and settles
//     at 1/32 as a slot's count grows. There are no collision checks.
//   * A match model: the last MINLEN = 24 bases are hashed to the position that
//     followed their previous occurrence, and the base found there is the
//     prediction. Its confidence is learned per match length, capped at 31. It
//     reaches back only to the start of its own chunk.
//   * A logistic mixer combines the eleven stretched predictions. The weight set
//     is chosen by the node and the previous two bases.
//   * A final APM (secondary estimation) stage refines the mixed probability,
//     using the node and the previous four bases. The coded probability is one
//     part mixer output to three parts APM output.
//   * A 32-bit carry-less binary arithmetic coder (BinEnc, BinDec) codes each
//     decision.
//
// Stream layout from seq_encode_mem: [nchunks:u32][total_symbols:u64], then per
// chunk [symbols:u32][coded_bytes:u32], then the chunk payloads back to back.
// Chunks are independent, so they are coded and decoded in parallel. The chunk
// count (4 by default, env SEQT in the encoder) is a COMPRESSION parameter, not
// only a threading one: this model's match component exploits long-range repeats,
// so more chunks mean less context per chunk. On the E. coli literal, 12 chunks
// gave 1,809,073 B, 4 chunks 1,804,223 B and 1 chunk 1,801,710 B.
//
// The development history of this coder, including the measurement that decided
// where its verification overhead sits, is in docs/extras/ENGINEERING_NOTES_CODERS.md.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <chrono>
#include <vector>
#include <thread>
#include <atomic>
#include <algorithm>

// Logistic function in fixed point: maps a stretched probability d, clamped to
// [-2047, 2047] (ln(p/(1-p)) scaled by 256), to a 12-bit probability, by
// interpolating a 33-point table. STRETCH, filled by initStretch(), is its inverse.
static short STRETCH[4096];
static inline int squash(int d){
    if(d>2047) d=2047; if(d<-2047) d=-2047;
    static const int t[33]={1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
        1546,2047,2549,2994,3348,3607,3785,3901,3975,4024,4050,4068,4079,4085,
        4089,4092,4093,4094};
    const int w=d&127; d=(d>>7)+16;
    return (t[d]*(128-w)+t[d+1]*w+64)>>7;
}
// Builds STRETCH as the inverse of squash. seq_encode_mem and seq_decode_mem call
// it, so the table is ready before any chunk is coded.
static void initStretch(){
    int pi=0;
    for(int x=-2047;x<=2047;++x){ const int v=squash(x);
        for(int p=pi;p<=v;++p) STRETCH[p]=(short)x; pi=v+1; }
    for(int p=pi;p<4096;++p) STRETCH[p]=2047;
}

// Carry-less binary arithmetic coder. The interval [x1, x2] is 32 bits wide and p
// is the 12-bit probability that the bit is 1, clamped to [1, 4094]. A byte is
// emitted whenever the top bytes of x1 and x2 agree. flush() writes four bytes.
// BinDec below is the exact inverse and must see the same probabilities.
struct BinEnc {
    std::vector<uint8_t> out; uint32_t x1=0,x2=0xFFFFFFFFu;
    void encode(int bit,int p){
        if(p<1) p=1; if(p>4094) p=4094;
        const uint32_t xmid=x1+(uint32_t)(((uint64_t)(x2-x1)*p)>>12);
        if(bit) x2=xmid; else x1=xmid+1;
        while(((x1^x2)&0xFF000000u)==0){ out.push_back((uint8_t)(x2>>24)); x1<<=8; x2=(x2<<8)|255; }
    }
    void flush(){ for(int i=0;i<4;++i){ out.push_back((uint8_t)(x2>>24)); x2<<=8; } }
};
struct BinDec {
    const uint8_t* p; const uint8_t* end; uint32_t x1=0,x2=0xFFFFFFFFu,x=0;
    void init(const uint8_t* b,size_t n){ p=b; end=b+n;
        for(int i=0;i<4;++i) x=(x<<8)|(p<end?*p++:0); }
    int decode(int pr){
        if(pr<1) pr=1; if(pr>4094) pr=4094;
        const uint32_t xmid=x1+(uint32_t)(((uint64_t)(x2-x1)*pr)>>12);
        const int bit=(x<=xmid);
        if(bit) x2=xmid; else x1=xmid+1;
        while(((x1^x2)&0xFF000000u)==0){ x1<<=8; x2=(x2<<8)|255; x=(x<<8)|(p<end?*p++:0); }
        return bit;
    }
};

// Context orders in bases of history (the 22-base order uses 44 of the 64 bits of
// `hist`). NM mixer inputs are the NORD orders plus the match model. MCTX is the
// number of mixer weight sets (3 nodes x the previous two bases), APMN the number
// of APM contexts (3 nodes x the previous four bases), MINLEN the match model's
// hashed context length in bases.
static constexpr int ORD[]={1,2,3,4,6,8,11,14,18,22};
static constexpr int NORD=(int)(sizeof(ORD)/sizeof(ORD[0]));
static constexpr int NM=NORD+1;
static constexpr int TBITS=16; static constexpr size_t TSZ=(size_t)1<<TBITS; static constexpr uint64_t TMASK=TSZ-1;
static constexpr int MMBITS_MAX=24;   // never worse than before
// mmtab, the match model's hash table, has 2^b slots with b = ceil(log2(chunk
// symbols)), at least 16 (64K slots, 256 KB, so a tiny chunk does not get a
// degenerate table) and at most MMBITS_MAX. It used to be a fixed 2^24 slots (64 MB)
// per chunk whatever the chunk's size. With 12 parallel chunks that was about
// 768 MB, which matched the pipeline's observed ~820 MB peak. The table is only a
// match-finder heuristic: a miss or collision gives no prediction (mlen = 0), so
// sizing it to the chunk has no effect on correctness. zstd scales its hash tables
// from a source-size hint for the same reason (ZSTD_adjustCParams).
static inline int mmbitsFor(size_t n){
    int b=1; while(((size_t)1<<b) < n && b<MMBITS_MAX) ++b;
    return std::max(16,b);   // floor of 16 (64K slots, 256KB) avoids a degenerately small table on tiny chunks
}
static constexpr int MINLEN=24;
static constexpr int MCTX=3*16;
static constexpr int APMN=3*256;

// Code one chunk with a FRESH model.
//   decode == false: ENCODES the n symbols sym0[0..n) into `out`.
//   decode == true:  DECODES n symbols from in_bytes[0..in_size) into `got`.
//                    sym0 is not read and may be null, and `out` is untouched.
static void run_chunk(const uint8_t* sym0, size_t n, bool decode,
                       const uint8_t* in_bytes, size_t in_size,
                       std::vector<uint8_t>& out, std::vector<uint8_t>& got){
    std::vector<uint16_t> tbl((size_t)NORD*TSZ, 2048);
    std::vector<uint8_t>  cnt((size_t)NORD*TSZ, 0);
    const int mmbits=mmbitsFor(n);
    const size_t mmsz=(size_t)1<<mmbits;
    std::vector<uint32_t> mmtab(mmsz,0);
    size_t mptr=0; uint32_t mlen=0;
    std::vector<uint16_t> mmconf(64*2, 2048);
    std::vector<int> wt((size_t)MCTX*NM, 65536/NM);
    std::vector<uint16_t> apm((size_t)APMN*33);
    for(int c=0;c<APMN;++c) for(int j=0;j<33;++j) apm[(size_t)c*33+j]=(uint16_t)(squash((j-16)*128)*16);

    BinEnc enc; if(!decode) enc.out.reserve(n/4);
    BinDec dec; if(decode) dec.init(in_bytes,in_size);
    if(decode) got.assign(n,0);
    uint64_t hist=0;
    int st[16]; size_t idx[16];

    for(size_t i=0;i<n;++i){
        const uint8_t s=decode?0:sym0[i];
        int dsym=0;
        for(int bp=0;bp<2;++bp){
            int bit=decode?0:((bp==0)?(s>>1):(s&1));
            const int node=(bp==0)?0:(1+(decode?dsym:(s>>1)));
            const int mc=(node*16)+(int)(hist&15);
            const int* w=&wt[(size_t)mc*NM];
            int dot=0;
            for(int m=0;m<NORD;++m){
                const int k=ORD[m];
                const uint64_t c=hist&((k>=32)?~0ULL:(((uint64_t)1<<(2*k))-1));
                uint64_t h=c*0x9E3779B97F4A7C15ULL;
                h^=(uint64_t)(k*0x2545F4914F6CDD1DULL);
                h^=(uint64_t)node*0xD1B54A32D192ED03ULL;
                h^=h>>29; h*=0xBF58476D1CE4E5B9ULL; h^=h>>32;
                idx[m]=(size_t)m*TSZ+(h&TMASK);
                st[m]=STRETCH[tbl[idx[m]]>>0];
                dot+=w[m]*st[m];
            }
            int mci=-1;
            if(mlen && mptr<i){
                const uint8_t pred=decode?got[mptr]:sym0[mptr];
                const int pbit=(bp==0)?(pred>>1):(pred&1);
                const bool usable=(bp==0)||((pred>>1)==(decode?dsym:(s>>1)));
                if(usable){
                    const int lc=(mlen>31)?31:(int)mlen;
                    mci=lc*2+pbit;
                    st[NORD]=STRETCH[mmconf[mci]];
                    if(!pbit) st[NORD]=-st[NORD];
                } else st[NORD]=0;
            } else st[NORD]=0;
            dot+=w[NORD]*st[NORD];
            dot>>=16;
            int p=squash(dot);
            const int ac=node*256+(int)(hist&255);
            const int sidx=(dot+2048)>>7, sw=(dot+2048)&127;
            const size_t a0=(size_t)ac*33+(size_t)((sidx<0)?0:(sidx>31?31:sidx));
            const int pa=(apm[a0]*(128-sw)+apm[a0+1]*sw)>>11;
            p=(p+3*pa)>>2;
            if(p<1) p=1; if(p>4094) p=4094;

            if(decode){ bit=dec.decode(p);
                        if(bp==0) dsym=bit; else got[i]=(uint8_t)((dsym<<1)|bit); }
            else        enc.encode(bit,p);

            const int err=((bit<<12)-p)*7;
            if(mci>=0){
                const int pbit=mci&1;
                const int tgt=(bit==pbit)?4095:0;
                mmconf[mci]=(uint16_t)((int)mmconf[mci]+((tgt-(int)mmconf[mci])>>5));
            }
            for(int m=0;m<NM;++m) wt[(size_t)mc*NM+m]+=(err*st[m])>>14;
            for(int m=0;m<NORD;++m){
                const int c=cnt[idx[m]];
                const int rate=(c<30)?(c+2):32;
                if(c<255) ++cnt[idx[m]];
                tbl[idx[m]]=(uint16_t)((int)tbl[idx[m]]+(((bit<<12)-(int)tbl[idx[m]])/rate));
            }
            const int g=(bit<<16)+(bit<<4)-bit-bit;
            apm[a0]  =(uint16_t)(apm[a0]  +((g-apm[a0])  >>6));
            apm[a0+1]=(uint16_t)(apm[a0+1]+((g-apm[a0+1])>>6));
        }
        const uint8_t cur=decode?got[i]:s;
        hist=(hist<<2)|cur;
        const uint8_t* refbuf = decode?got.data():sym0;
        if(mlen && mptr<i && refbuf[mptr]==cur){ ++mptr; if(mlen<65535) ++mlen; }
        else { mlen=0; mptr=0; }
        if(i+1>=(size_t)MINLEN){
            uint64_t hh=(hist&(((uint64_t)1<<(2*MINLEN))-1))*0x9E3779B97F4A7C15ULL;
            hh^=hh>>31; const size_t slot=(size_t)(hh&(mmsz-1));
            if(!mlen){ const uint32_t cand=mmtab[slot];
                       if(cand){ mptr=cand; mlen=1; } }
            mmtab[slot]=(uint32_t)(i+1);
        }
    }
    if(!decode){ enc.flush(); out=std::move(enc.out); }
}

// ---- in-memory encoder: codes `sym` (values 0..3) as NCHUNKS independent chunks of
// ceil(n/NCHUNKS) symbols on NTHREADS threads, and returns the chunk table followed
// by the payloads (layout in the header comment). No file I/O. This function does
// not decode its own output: round-trip checking belongs to the archive-level
// verification, not to this stream.
static std::vector<uint8_t> seq_encode_mem(const std::vector<uint8_t>& sym,
                                           unsigned NCHUNKS, unsigned NTHREADS){
    initStretch();
    const size_t n=sym.size();
    if(!n) return {};
    if(!NCHUNKS) NCHUNKS=1;
    if(!NTHREADS) NTHREADS=1;
    const size_t chunkSz=(n+NCHUNKS-1)/NCHUNKS;
    std::vector<std::vector<uint8_t>> chunkOut(NCHUNKS);
    {
        std::atomic<size_t> next{0};
        auto worker=[&](){
            size_t c;
            while((c=next.fetch_add(1))<NCHUNKS){
                size_t s=c*chunkSz, e=std::min(n,s+chunkSz);
                if(s>=e) continue;
                std::vector<uint8_t> dummy_got;
                run_chunk(sym.data()+s, e-s, false, nullptr, 0, chunkOut[c], dummy_got);
            }
        };
        std::vector<std::thread> pool;
        for(unsigned t=0;t<NTHREADS;++t) pool.emplace_back(worker);
        for(auto& th:pool) th.join();
    }
    // CHUNK TABLE. Without boundaries a decoder cannot tell where one chunk ends
    // and the next begins, so the stream would be undecodable whatever the coder.
    // Emit
    //   [nchunks:u32][total_syms:u64] then per chunk [raw_syms:u32][coded:u32]
    // before the payloads: 12 + 8*nchunks bytes, 44 B at the default 4 chunks,
    // against a literal of 0.6-7.8 MB.
    std::vector<uint8_t> all;
    size_t tot=0; for(auto& b:chunkOut) tot+=b.size();
    all.reserve(tot + 16 + 8*NCHUNKS);
    auto put32=[&](uint32_t v){ all.insert(all.end(),(uint8_t*)&v,(uint8_t*)&v+4); };
    auto put64=[&](uint64_t v){ all.insert(all.end(),(uint8_t*)&v,(uint8_t*)&v+8); };
    put32((uint32_t)NCHUNKS); put64((uint64_t)n);
    for(unsigned c=0;c<NCHUNKS;++c){
        const size_t st=(size_t)c*chunkSz, en=std::min(n,st+chunkSz);
        put32((uint32_t)(en>st?en-st:0));
        put32((uint32_t)chunkOut[c].size());
    }
    for(auto& b:chunkOut){ all.insert(all.end(), b.begin(), b.end()); b.clear(); b.shrink_to_fit(); }
    return all;
}

// Inverse of seq_encode_mem. Reads the chunk table, then replays run_chunk in
// decode mode over each chunk. run_chunk already had a `decode` parameter and a
// BinDec; it had simply never been called from an archive path.
static std::vector<uint8_t> seq_decode_mem(const uint8_t* d, size_t n){
    initStretch();
    if(n < 12) return {};
    uint32_t nch; memcpy(&nch,d,4);
    uint64_t total; memcpy(&total,d+4,8);
    if(!nch || nch > (1u<<20)) return {};
    const size_t tabOff = 12;
    if(n < tabOff + (size_t)nch*8) return {};
    std::vector<uint32_t> rawN(nch), codN(nch);
    for(uint32_t c=0;c<nch;++c){
        memcpy(&rawN[c], d+tabOff+(size_t)c*8,   4);
        memcpy(&codN[c], d+tabOff+(size_t)c*8+4, 4);
    }
    size_t off = tabOff + (size_t)nch*8;
    std::vector<std::vector<uint8_t>> got(nch);
    // Chunks are independent -- each has its own coded range and its own output,
    // and the chunk table gives every offset up front. Decoding them one after
    // another made seq_decode_mem 93% of the C++ decoder's wall clock on
    // L. major, while the table needed to parallelise it already existed: it was
    // added so the stream could be decoded at all.
    std::vector<size_t> offs(nch);
    { size_t o = off;
      for(uint32_t c=0;c<nch;++c){ offs[c]=o; o += codN[c]; }
      if(o > n) return {};
      off = o; }
    {
        std::atomic<uint32_t> next{0};
        unsigned T = std::min<unsigned>(nch, std::max(1u, std::thread::hardware_concurrency()));
        std::vector<std::thread> th;
        for(unsigned t=0;t<T;++t) th.emplace_back([&]{
            for(;;){
                uint32_t c = next.fetch_add(1);
                if(c >= nch) break;
                if(!rawN[c]) continue;
                std::vector<uint8_t> dummy_out;
                run_chunk(nullptr, rawN[c], true, d+offs[c], codN[c], dummy_out, got[c]);
            }
        });
        for(auto& x : th) x.join();
    }
    std::vector<uint8_t> all; all.reserve(total);
    for(auto& g : got) all.insert(all.end(), g.begin(), g.end());
    if(all.size()!=total) return {};
    return all;
}
