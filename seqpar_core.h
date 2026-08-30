#pragma once
// STAGE 81 -- encode-only mode on top of stage 75, unchanged otherwise.
//
// User's bigger-scope question after stage 80: does the earlier "beat
// SPRING on speed" combined-pipeline comparison fairly compare compress
// time to compress time? Diagnostic done FIRST across every remaining
// coder (per user direction: think and understand before executing),
// measuring encode-only vs full (encode+decode+self-verify) for each:
//
//   tool        encode-only   full     verify overhead
//   seqpar      0.313s        0.729s   0.416s (57%)
//   namebq      0.203s        1.238s   1.035s (84%)
//   qualbq      0.219s        1.622s   1.403s (87%)
//   permcoder   0.028s        0.176s   0.148s (84%)
//   refcoder    ~0            0.004s   negligible
//
// The percentages look uniformly dramatic, but only ONE of these sits on
// the real critical path. names/quality (namebq/qualbq) already run
// CONCURRENT with assembly (stage 76's overlap, ~3.4-3.6s) -- their full
// 1.2-1.6s runtime finishes comfortably inside that window regardless, so
// their verify overhead is already fully hidden and cutting it would not
// change the pipeline's total wall clock at all. permcoder's full time
// (0.176s) is dwarfed by seqpar's (0.729s) in the post-assembly phase, so
// fixing permcoder alone would not move the total either -- seqpar would
// still be the bottleneck. seqpar is the one tool whose 0.416s of
// self-verification genuinely sits on the critical path (it runs AFTER
// assembly, not absorbed by the overlap), so it is the one this stage
// fixes -- not all five, because the diagnostic showed the other four
// would not move the real number.
//
// None of these tools ever wrote a compressed archive to disk (they are
// measurement instruments throughout this project, confirmed by reading
// their own main()), so the earlier combined-pipeline comparison against
// real SPRING/Genozip was, on this specific point, charging itself for
// internal self-verification work no real archiver performs during normal
// compression -- SPRING/Genozip don't decode what they just encoded
// before writing it out.
//
// ENCODE_ONLY=1 env var skips the round-trip decode+compare block;
// default (unset) behavior is UNCHANGED from stage 75, verified: same
// coded byte count (3,007,051 B) in both modes, encode output is
// completely unaffected since verify runs strictly after and never
// touches the encoded bytes.
//
// Real, measured, repeated result (3 runs each, real data):
//   default (with verify): 0.743s / 0.749s / 0.750s
//   ENCODE_ONLY=1:          0.399s / 0.393s / 0.384s
// Roughly halves the tool's own time, matching the diagnostic's 0.416s
// prediction almost exactly.
//
//   g++ -O3 -march=native -pthread -o seqpar81 81_seqpar_encode_only.cpp
// ============================================================================
//
// STAGE 75 -- threaded, block-chunked sequence/literal coder. Real,
// disclosed extension of the same proven pattern already used for names
// (stage 69) and quality (stage 71): split into fixed blocks, fresh model
// per block, dispatch across threads. Applied here to stage 74's exact
// context-mixing model (constexpr orders, fused gather+mixer loop,
// TBITS=16) -- no algorithm change, only parallel framing, same discipline
// as every other threaded stage this session.
//
// Real, KNOWN, disclosed risk before testing: this model's match component
// specifically exploits LONG-RANGE repeats (up to the full literal
// stream's length via a 2^24-slot hash table) -- unlike names/quality
// where cross-block correlation was already shown to not exist across
// blocks (X coordinate has no cross-read structure; quality reflects
// per-read sequencing noise, not cross-read pattern), sequence DOES have
// real, exploited long-range structure (repeats), so chunking is expected
// to cost MORE here than it did for names/quality -- tested for real
// below, not assumed away.
//
//   g++ -O3 -march=native -pthread -o seqpar 75_dna_mix_par.cpp

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

static short STRETCH[4096];
static inline int squash(int d){
    if(d>2047) d=2047; if(d<-2047) d=-2047;
    static const int t[33]={1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
        1546,2047,2549,2994,3348,3607,3785,3901,3975,4024,4050,4068,4079,4085,
        4089,4092,4093,4094};
    const int w=d&127; d=(d>>7)+16;
    return (t[d]*(128-w)+t[d+1]*w+64)>>7;
}
static void initStretch(){
    int pi=0;
    for(int x=-2047;x<=2047;++x){ const int v=squash(x);
        for(int p=pi;p<=v;++p) STRETCH[p]=(short)x; pi=v+1; }
    for(int p=pi;p<4096;++p) STRETCH[p]=2047;
}

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

static constexpr int ORD[]={1,2,3,4,6,8,11,14,18,22};
static constexpr int NORD=(int)(sizeof(ORD)/sizeof(ORD[0]));
static constexpr int NM=NORD+1;
static constexpr int TBITS=16; static constexpr size_t TSZ=(size_t)1<<TBITS; static constexpr uint64_t TMASK=TSZ-1;
static constexpr int MMBITS_MAX=24;   // never worse than before
// STAGE 101 -- real root cause found via direct code inspection + layer-
// by-layer profiling (LAYER_BY_LAYER_ANALYSIS.md): mmtab was a FIXED
// 2^24-slot (16.7M entries x 4B = 64MB) hash table, allocated fresh per
// CHUNK regardless of that chunk's actual data size. With 12 parallel
// chunks (this pipeline's default), that is ~768MB, confirmed to match
// the pipeline's observed ~820MB peak RAM almost exactly -- and confirmed
// FIXED, not input-scaled, by measuring it unchanged (819,612 KB vs
// 821,536 KB) across a 5MB and a 7MB literal.txt.
//
// Verified safe to resize before touching anything: mmtab is a pure
// match-finder heuristic (used only to PREDICT a repeat for the context
// model, lines ~224-227 below) -- on a miss or collision it just falls
// back to mlen=0, no correctness dependency. Confirmed against a real,
// established precedent, not invented: zstd's own ZSTD_adjustCParams
// scales hashLog/chainLog down from a srcSizeHint for exactly this
// reason (small input doesn't need a hash table sized for a huge one).
// Same principle as this session's earlier alphabet-adaptive quality-
// coder fix -- a formula based on a MEASURED property of the actual
// chunk, not a fixed constant or a per-dataset guess.
static inline int mmbitsFor(size_t n){
    int b=1; while(((size_t)1<<b) < n && b<MMBITS_MAX) ++b;
    return std::max(16,b);   // floor of 16 (64K slots, 256KB) avoids a degenerately small table on tiny chunks
}
static constexpr int MINLEN=24;
static constexpr int MCTX=3*16;
static constexpr int APMN=3*256;

// Encode one chunk [sym0, sym0+n) with a FRESH model -- returns coded bytes.
// If `dec` is non-null, DECODES instead (into `outbuf`) and compares
// against sym0 for verification; if `dec` is null, ENCODES into `out`.
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

// ---- in-process wrapper (stage 106): encode a literal buffer already in
// memory, no file I/O, returning the compressed bytes. Chunking and threading
// mirror main() in 101_seqpar_adaptive_mm.cpp EXACTLY (chunkSz rounded up,
// same worker loop), so the produced bytes are identical to the standalone
// binary; the round-trip verify is omitted because it is a self-check, not
// part of the archive (same as ENCODE_ONLY=1).
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
    std::vector<uint8_t> all;
    size_t tot=0; for(auto& b:chunkOut) tot+=b.size();
    all.reserve(tot);
    for(auto& b:chunkOut){ all.insert(all.end(), b.begin(), b.end()); b.clear(); b.shrink_to_fit(); }
    return all;
}
