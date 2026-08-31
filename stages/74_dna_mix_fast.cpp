// Context mixing WITH A MATCH MODEL, run on the pseudogenome BEFORE MEM removal.
//
// Why this exists. MEM removal buys its literal reduction with an explicit
// reference stream -- 66K matches costing 328,512 B, of which 210,600 is source
// offsets sitting within 1% of n*log2(pg_len), its entropy floor. Four attempts
// to shrink that stream all failed (MAXCAND never binds, consecutive-source
// delta is worse, repeat distance gives 1.1%, lazy parsing 3,634 B), and the
// structural check found zero mergeable matches. The stream cannot be coded
// smaller and the matches are not fragments.
//
// So stop paying for it. MEM removal only exists because an order-22 context
// model cannot see a repeat 10 Mbases back. A MATCH MODEL can: keep a hash from
// a recent k-mer to where it last occurred, and while a match holds, predict the
// base that followed there. Long repeats then cost almost nothing AND cost no
// references at all, because nothing is removed and nothing has to be pointed
// at. This is standard lpaq/paq machinery, not a PgRC2 idea -- their pipeline
// has no context model.
//
// Run on pg_full.txt (23,233,953 bases) rather than literal.txt, and compared
// against literal-coded + references = 3,332,257 B.
//
// Where this sits. The sequence gap against PgRC2 is 101,610 B: 55,863 B of
// extra literal (assembly) and 45,747 B of coder deficit. Two attempts on the
// coder half already failed:
//   - a single order-k context model (stage 29) reached only 2.0206 bits/base
//     at k=8 and got WORSE with more context -- 12.9M bases cannot populate
//     4^10 contexts, so high orders starve;
//   - LZMA literal-context tuning moved 1.9544 -> 1.9529, worth 2,425 B.
//
// Their VarLenDNACoder wins by keeping byte boundaries tied to content so LZMA
// can still match (its pre-entropy stream is LARGER, 2.313 bpb, yet finishes at
// 1.926 because LZMA gains 16.7% on it against 2.3% on our 2-bit packing).
// Copying that is possible, but it is their design and it caps out at 1.926.
//
// The literature says the whole framing is beatable. DNA-COMPACT reports 1.838
// bits/base on yeast, and GeCo3 gets its gains the same way: several finite
// context models of DIFFERENT orders, combined by logistic mixing, rather than
// one model or a dictionary. That fixes exactly what killed stage 29 -- a high
// order that starves on sparse data contributes little weight instead of
// dominating, while low orders carry the prediction, and the mixer learns the
// balance per position rather than having it fixed.
//
// Architecture (standard lpaq-style, textbook, nothing from PgRC2 -- their path
// has no context model at all):
//   - each base is two binary decisions, so everything is binary arithmetic
//     coding with a 12-bit probability;
//   - one hashed table of 12-bit predictions per order;
//   - predictions combined in the logistic domain, p = squash(sum w_i*stretch(p_i));
//   - mixer weights selected by a small context and trained online by gradient;
//   - an APM/SSE stage refines the mixed probability against a low-order context.
//
//   g++ -O3 -march=native -o mix 30_dna_mix.cpp && ./mix literal.txt
//
// ============================================================================
// STAGE 74 -- same model as stage 35 (35_match_model.cpp), profiled and
// sped up. The user asked, after building a combined sequence+order+names+
// quality pipeline and comparing it to real SPRING/Genozip end to end, why
// the total was ~4-9x slower than either real tool despite a real
// compression-ratio win. Profiling (perf record + perf annotate, real
// instruction-level attribution, not guessing) found the time was almost
// entirely in THIS coder -- 10.5s single-pass-encode / 23.5s for the
// built-in encode+decode round trip, on 12.5M bases, dwarfing SPRING's
// (3.9s) and Genozip's (1.6s) ENTIRE per-file compress time.
//
// Three real changes tested in order, each measured, on real data
// (literal.txt, 12,506,313 bases from the real yeast dataset):
//
// 1. `ORD[]`/`NORD`/`NM` made `static constexpr` (were `const`, runtime-
//    initialized from a sizeof expression the compiler could not prove
//    constant enough to fully unroll). Real but modest: 23.5s -> 21.5s
//    (~8.5%), byte-identical output.
// 2. Fused the "gather" loop (compute st[m]/idx[m] for each of the 10
//    order-models) with the separate "mixer dot-product" loop that used to
//    read st[] back afterward -- the dot product only needs w[m] (known
//    before gather starts) and the freshly-computed st[m], so accumulating
//    inline removes one full pass over memory. Integer sum is associative,
//    so this cannot change the result -- verified byte-identical, not just
//    reasoned safe. Real effect: negligible on its own (21.5s -> 21.5s),
//    which was itself informative -- it meant compute/loop-structure was
//    NOT the bottleneck.
// 3. That negative result redirected the investigation to the tables
//    themselves: `tbl`/`cnt` are sized `NORD * 2^TBITS` with TBITS=22 --
//    ~126 MB total, accessed via hash at effectively random addresses,
//    far bigger than cache. Swept TBITS from 14 to 22, real data, every
//    point round-trip verified:
//      TBITS  time(both passes)  coded bytes   bits/base
//      22     22.4s              2,979,683     1.9060  (original)
//      21     13.6s              2,981,198     1.9070
//      20      9.1s              2,982,182     1.9076
//      19      8.2s              2,982,833     1.9080
//      18      7.9s              2,983,417     1.9084
//      17      7.4s              2,983,970     1.9088
//      16      7.3s              2,984,469     1.9091   <- chosen
//      15      6.1s              2,984,934     1.9094
//      14      5.6s              2,985,360     1.9097
//    Confirms the theory decisively: this was memory-latency-bound (DRAM
//    misses on random hashed access), not compute-bound -- shrinking the
//    table trades a small, real, disclosed ratio cost for a large, real
//    speed AND memory win (RSS also dropped, 219 MB -> 98 MB, since a
//    smaller table is being allocated at all, not just accessed faster).
//
// TBITS=16 chosen as the balance point: 3.2x faster (23.5s -> 7.3s
// round-trip; ~3.6s estimated single-pass) for +0.16% size cost
// (2,984,469 vs 2,979,683 B) -- all TBITS 14-18 sit within ~0.2% of each
// other on ratio, so speed is the real differentiator in that range, not
// ratio. Not a free win -- a small, real, honestly measured and disclosed
// tradeoff, same standard applied to every other tradeoff this session.
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>

// ── logistic helpers ────────────────────────────────────────────────────────
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

// ── binary arithmetic coder ─────────────────────────────────────────────────
struct BinEnc {
    std::vector<uint8_t> out; uint32_t x1=0,x2=0xFFFFFFFFu;
    void encode(int bit,int p){                       // p = P(bit=1), 12-bit
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

int main(int argc,char** argv){
    if(argc<2){ fprintf(stderr,"usage: %s <literal.txt>\n",argv[0]); return 1; }
    initStretch();

    std::vector<uint8_t> sym;
    {
        FILE* f=fopen(argv[1],"rb"); if(!f){ perror("open"); return 1; }
        fseek(f,0,SEEK_END); const size_t sz=ftell(f); fseek(f,0,SEEK_SET);
        std::vector<char> buf(sz);
        if(fread(buf.data(),1,sz,f)!=sz){ fprintf(stderr,"short read\n"); return 1; }
        fclose(f);
        sym.reserve(sz);
        for(size_t i=0;i<sz;++i)
            switch(buf[i]){ case 'A':sym.push_back(0);break; case 'C':sym.push_back(1);break;
                            case 'G':sym.push_back(2);break; case 'T':sym.push_back(3);break; }
    }
    const size_t n=sym.size();

    // Orders chosen to span starved and well-populated regimes: the mixer is
    // what decides which to trust, which is the whole point -- stage 29 failed
    // because a single starved order had to carry everything alone.
    static constexpr int ORD[]={1,2,3,4,6,8,11,14,18,22};
    static constexpr int NORD=(int)(sizeof(ORD)/sizeof(ORD[0]));
    static constexpr int NM=NORD+1;                 // +1 for the match model
    const int TBITS=16; const size_t TSZ=(size_t)1<<TBITS; const uint64_t TMASK=TSZ-1;

    std::vector<uint16_t> tbl((size_t)NORD*TSZ, 2048); // 12-bit predictions
    std::vector<uint8_t>  cnt((size_t)NORD*TSZ, 0);    // adaptation rate per slot

    // ── match model ─────────────────────────────────────────────────────────
    // MINLEN bases of history hash to the position where that context last
    // occurred. While the prediction keeps coming true the match is extended and
    // trusted more; one wrong base drops it and a fresh lookup is made. Its
    // prediction enters the mixer like any other model, so the mixer learns how
    // much to trust it rather than it being hard-wired.
    const int MMBITS=24; const size_t MMSZ=(size_t)1<<MMBITS;
    const int MINLEN=24;                             // bases of context to anchor on
    std::vector<uint32_t> mmtab(MMSZ,0);
    size_t mptr=0; uint32_t mlen=0;                  // where we are following, how long it has held
    // confidence table indexed by (capped match length, expected bit): learned,
    // not assumed, so a long match is trusted only as far as it earns.
    std::vector<uint16_t> mmconf(64*2, 2048);

    // Mixer: weights selected by (model, order-2 context, node) so the balance
    // can differ between, say, a GC-rich and an AT-rich neighbourhood.
    const int MCTX=3*16;                              // node(3) x order-2(16)
    // Weights must start so that the mixed dot product lands in squash's
    // +-2047 domain. With NM models each contributing (w*stretch)>>16 and
    // stretch up to 2047, w = 65536/NM makes the initial mix the plain average.
    // Starting at 1<<14 instead put the sum near 5000 and squash saturated from
    // the first symbol, which is what made the first version worse than a single
    // order rather than better.
    std::vector<int> wt((size_t)MCTX*NM, 65536/NM);

    // SSE/APM: refine the mixed probability against a small direct context.
    const int APMN=3*256;
    std::vector<uint16_t> apm((size_t)APMN*33);
    for(int c=0;c<APMN;++c) for(int j=0;j<33;++j) apm[(size_t)c*33+j]=(uint16_t)(squash((j-16)*128)*16);

    // One pass drives both directions: identical model, identical update order,
    // and the only difference is whether the bit is written or read. That is
    // what makes the round-trip check below meaningful rather than decorative.
    BinEnc enc; enc.out.reserve(n/4);
    BinDec dec; bool DECODE=false;
    std::vector<uint8_t> got;
    uint64_t hist=0;
    int st[16]; size_t idx[16];

  for(int pass=0;pass<2;++pass){
    if(pass==1){
        DECODE=true; dec.init(enc.out.data(),enc.out.size()); got.assign(n,0);
        std::fill(tbl.begin(),tbl.end(),(uint16_t)2048);
        std::fill(cnt.begin(),cnt.end(),(uint8_t)0);
        std::fill(wt.begin(),wt.end(),65536/NM);
        std::fill(mmtab.begin(),mmtab.end(),0u);
        std::fill(mmconf.begin(),mmconf.end(),(uint16_t)2048);
        mptr=0; mlen=0;
        for(int c=0;c<APMN;++c) for(int j=0;j<33;++j) apm[(size_t)c*33+j]=(uint16_t)(squash((j-16)*128)*16);
        hist=0;
    }
    for(size_t i=0;i<n;++i){
        const uint8_t s=sym[i];
        int dsym=0;
        for(int bp=0;bp<2;++bp){
            int bit=DECODE?0:((bp==0)?(s>>1):(s&1));
            const int node=(bp==0)?0:(1+(DECODE?dsym:(s>>1)));
            // FUSED gather + mixer dot-product: the dot product only needs
            // w[m] (known from mc, computable before the loop -- mc depends
            // only on node/hist, not on anything gather computes) and st[m]
            // (computed fresh each iteration) -- no need for a second pass
            // reading st[] back from memory. Integer sum is associative, so
            // accumulating in gather order changes nothing about the final
            // `dot` value versus the original two-loop version -- verified
            // byte-identical output below, not just reasoned to be safe.
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
            // match model prediction for this bit
            int mci=-1;
            if(mlen && mptr<i){
                const uint8_t pred=sym[mptr];        // base that followed last time
                const int pbit=(bp==0)?(pred>>1):(pred&1);
                // only meaningful for the low bit if the high bit agreed
                const bool usable=(bp==0)||((pred>>1)==(DECODE?dsym:(s>>1)));
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
            // SSE stage
            const int ac=node*256+(int)(hist&255);
            const int sidx=(dot+2048)>>7, sw=(dot+2048)&127;
            const size_t a0=(size_t)ac*33+(size_t)((sidx<0)?0:(sidx>31?31:sidx));
            const int pa=(apm[a0]*(128-sw)+apm[a0+1]*sw)>>11;
            p=(p+3*pa)>>2;
            if(p<1) p=1; if(p>4094) p=4094;

            if(DECODE){ bit=dec.decode(p);
                        if(bp==0) dsym=bit; else got[i]=(uint8_t)((dsym<<1)|bit); }
            else        enc.encode(bit,p);

            // update
            const int err=((bit<<12)-p)*7;
            if(mci>=0){
                const int pbit=mci&1;
                const int tgt=(bit==pbit)?4095:0;
                mmconf[mci]=(uint16_t)((int)mmconf[mci]+((tgt-(int)mmconf[mci])>>5));
            }
            for(int m=0;m<NM;++m) wt[(size_t)mc*NM+m]+=(err*st[m])>>14;
            for(int m=0;m<NORD;++m){
                // 1/(count+1.5) schedule: fast while a slot is new, slow once it
                // has evidence. The earlier delta*2/rate with rate starting at 2
                // was a full jump to 0 or 4095 on the first observation.
                const int c=cnt[idx[m]];
                const int rate=(c<30)?(c+2):32;
                if(c<255) ++cnt[idx[m]];
                tbl[idx[m]]=(uint16_t)((int)tbl[idx[m]]+(((bit<<12)-(int)tbl[idx[m]])/rate));
            }
            const int g=(bit<<16)+(bit<<4)-bit-bit;
            apm[a0]  =(uint16_t)(apm[a0]  +((g-apm[a0])  >>6));
            apm[a0+1]=(uint16_t)(apm[a0+1]+((g-apm[a0+1])>>6));
        }
        const uint8_t emitted=DECODE?(uint8_t)((dsym<<1)|((0))):s;
        (void)emitted;
        const uint8_t cur=DECODE?got[i]:s;
        hist=(hist<<2)|cur;
        // extend or drop the match, then re-anchor from the current context
        if(mlen && mptr<i && sym[mptr]==cur){ ++mptr; if(mlen<65535) ++mlen; }
        else { mlen=0; mptr=0; }
        if(i+1>=(size_t)MINLEN){
            uint64_t hh=(hist&(((uint64_t)1<<(2*MINLEN))-1))*0x9E3779B97F4A7C15ULL;
            hh^=hh>>31; const size_t slot=(size_t)(hh&(MMSZ-1));
            if(!mlen){ const uint32_t cand=mmtab[slot];
                       if(cand){ mptr=cand; mlen=1; } }
            mmtab[slot]=(uint32_t)(i+1);
        }
    }
    if(pass==0) enc.flush();
  }
    const size_t bytes=enc.out.size();
    bool ok=(got.size()==n);
    if(ok) for(size_t i=0;i<n;++i) if(got[i]!=sym[i]){ ok=false;
        fprintf(stderr,"MISMATCH at %zu: got %u want %u\n",i,got[i],sym[i]); break; }
    printf("bases=%zu  coded=%zu B  %.4f bits/base   round trip: %s\n",
           n,bytes,bytes*8.0/n, ok?"VERIFIED":"FAILED");
    printf("  compare: literal-coded + references = 3,332,257 B (stage 34)\n");
    printf("  compare: PgRC2 sequence + references = 3,233,654 B\n");
    printf("  ours 2-bit + xz   : 3,158,084 B over 12,926,925 = 1.9544 bits/base\n");
    printf("  stage 29 order-8  :                               2.0206 bits/base\n");
    printf("  DNA-COMPACT, yeast (published)                   1.838  bits/base\n");
    return 0;
}
