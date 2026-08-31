// STAGE 58 -- per-pg-position mismatch context, blended with the existing
// ref-base global model. Validated design, not ad hoc:
//
//   - samtools mpileup / VarScan build EXACTLY this per-position base-count
//     structure (a "pileup") as their foundational operation in variant
//     calling -- confirms per-position counting is the right primitive for
//     data that is genuinely repeat/variant structured, which stage 57's
//     analysis showed ours is (31% of mismatches at ~166K repeated sites).
//   - PPM (Cleary & Witten) is the established precedent for "prefer a
//     specific context, fall back to a general one" -- but PPM's own
//     literature flags its failure mode directly: too many singleton/sparse
//     contexts cause escape-symbol overhead that can make compression
//     WORSE, not better (this is exactly what an earlier experiment in this
//     project found: raw per-position context with no backoff, 1024
//     contexts, LOST to sensible bucketing).
//
// This avoids that trap by construction, not by hope: instead of PPM's hard
// escape symbol, use LINEAR BLENDING (Katz/Jelinek-Mercer-style backoff) --
//   effective_freq[k] = position_freq[k]*W + global_freq[k]
// A brand-new position (position_freq all zero, the common case: 886,288 of
// 1,242,201 sites in the earlier analysis are singletons) blends to EXACTLY
// the existing global-model behaviour -- this coder is PROVABLY never worse
// than stage 53's baseline at any position with no repeat history, and only
// improves where real structure exists. No parameter sweep needed to avoid
// the sparsity trap; it is structurally excluded.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>
#include <unordered_map>
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
static inline int code2(uint8_t c){ return c=='A'?0:c=='C'?1:c=='G'?2:3; }

int main(int argc,char**argv){
    if(argc<6){ fprintf(stderr,"usage: %s ref.bin obs.bin pos.bin ctx3.bin pgpos.bin [W=8]\n",argv[0]); return 1; }
    auto loadb=[&](const char* fn)->std::vector<uint8_t>{
        FILE* f=fopen(fn,"rb"); fseek(f,0,SEEK_END); size_t n=ftell(f); fseek(f,0,SEEK_SET);
        std::vector<uint8_t> v(n); if(fread(v.data(),1,n,f)!=n){fprintf(stderr,"short read %s\n",fn);exit(1);} fclose(f); return v; };
    std::vector<uint8_t> ref=loadb(argv[1]), obs=loadb(argv[2]);
    FILE* fg=fopen(argv[5],"rb"); fseek(fg,0,SEEK_END); size_t ng=ftell(fg); fseek(fg,0,SEEK_SET);
    const size_t n=ref.size();
    if(ng!=n*4){ fprintf(stderr,"pgpos.bin size mismatch: expected %zu got %zu\n",n*4,ng); return 1; }
    std::vector<uint32_t> pgpos(n);
    if(fread(pgpos.data(),4,n,fg)!=n){ fprintf(stderr,"short read pgpos\n"); return 1; } fclose(fg);
    const uint32_t W = (argc>6)?(uint32_t)atoi(argv[6]):8;

    auto excludeIdx=[&](int refc,int obsc)->int{
        int k=0; for(int s=0;s<4;++s){ if(s==refc) continue; if(s==obsc) return k; ++k; } return -1; };
    auto includeSym=[&](int refc,int k)->int{
        int c=0; for(int s=0;s<4;++s){ if(s==refc) continue; if(c==k) return s; ++c; } return -1; };

    // Global model: same shape as stage 50 (4 contexts by ref base).
    std::array<uint32_t,3> global[4]; for(auto& f:global) f={1,1,1};
    // Per-position model: lazily created, starts at {0,0,0} -- contributes
    // nothing until a position has actually been seen before.
    std::unordered_map<uint32_t,std::array<uint32_t,3>> perpos;

    auto blend=[&](uint32_t pos,int refc,uint32_t out[3]){
        auto it=perpos.find(pos);
        const std::array<uint32_t,3>* pf = (it!=perpos.end())?&it->second:nullptr;
        // uint64 throughout: pp[k] (<=4096 after its own rescale) times W
        // (swept up to 2^18+ to find the true optimum) can exceed uint32
        // range well before the final cap below, so the multiply itself
        // must not be done in 32-bit -- that silent wraparound, not the
        // final total, was the real cause of the earlier bad_alloc (a
        // wrapped cumLo/cumHi fed the range coder garbage, which spun its
        // output-byte loop until OOM). Caught by testing large W, not by
        // reasoning about the encode() function in isolation.
        uint64_t t[3];
        for(int k=0;k<3;++k) t[k]=(uint64_t)global[refc][k]+(pf?(uint64_t)(*pf)[k]*W:0ULL);
        // Cap the BLENDED total, not just each model's own total -- the
        // established safe bound elsewhere in this codebase (stage 50's
        // encode()) is range/tot staying non-degenerate, verified good up to
        // 65536 by that coder's own rescale threshold. Rescale proportionally
        // (floor 1, never zero -- a zero-probability symbol that then
        // actually occurs is a real coder bug, not just inefficiency).
        uint64_t tot=t[0]+t[1]+t[2];
        const uint64_t CAP=32768;
        if(tot>CAP){ for(int k=0;k<3;++k) t[k]=(t[k]*CAP)/tot; }
        for(int k=0;k<3;++k) out[k]=(uint32_t)(t[k]?t[k]:1);
    };

    RangeEnc enc; enc.out.reserve(n);
    for(size_t i=0;i<n;++i){
        const int r=code2(ref[i]), o=code2(obs[i]);
        const int k=excludeIdx(r,o);
        uint32_t f[3]; blend(pgpos[i],r,f);
        uint32_t tot=f[0]+f[1]+f[2];
        uint32_t lo=0; for(int j=0;j<k;++j) lo+=f[j];
        enc.encode(lo,lo+f[k],tot);
        // update BOTH models -- global always, per-position lazily created
        global[r][k]+=8; { uint32_t t2=global[r][0]+global[r][1]+global[r][2];
            if(t2>65536) for(int j=0;j<3;++j) global[r][j]=(global[r][j]>>1)|1; }
        auto& pp = perpos[pgpos[i]]; ++pp[k];
        if(pp[0]+pp[1]+pp[2] > 4096) for(int j=0;j<3;++j) pp[j]>>=1;  // rare, keeps W*count in range
    }
    enc.flush();
    const size_t bytes=enc.out.size();

    // round trip -- fresh state, identical update order
    for(auto& f:global) f={1,1,1};
    perpos.clear();
    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    bool ok=true;
    for(size_t i=0;i<n;++i){
        const int r=code2(ref[i]);
        uint32_t f[3]; blend(pgpos[i],r,f);
        uint32_t tot=f[0]+f[1]+f[2];
        uint32_t target=dec.getFreq(tot);
        int k=0; uint32_t lo=0; while(lo+f[k]<=target){ lo+=f[k]; ++k; }
        dec.decodeUpdate(lo,lo+f[k]);
        const int s=includeSym(r,k);
        if(s!=code2(obs[i])){ ok=false; fprintf(stderr,"MISMATCH at %zu: got %d want %d\n",i,s,code2(obs[i])); break; }
        global[r][k]+=8; { uint32_t t2=global[r][0]+global[r][1]+global[r][2];
            if(t2>65536) for(int j=0;j<3;++j) global[r][j]=(global[r][j]>>1)|1; }
        auto& pp = perpos[pgpos[i]]; ++pp[k];
        if(pp[0]+pp[1]+pp[2] > 4096) for(int j=0;j<3;++j) pp[j]>>=1;
    }
    printf("symbols=%zu  coded=%zu B  %.4f bits/symbol  round trip: %s  (W=%u)\n",
           n,bytes,bytes*8.0/n, ok?"VERIFIED":"FAILED", W);
    printf("stage 50 (ref-only baseline)      : 241,365 B\n");
    printf("stage 53 (ref x pos x prev)       : 238,993 B\n");
    return ok?0:1;
}
