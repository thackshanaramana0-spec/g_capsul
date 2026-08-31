// STAGE 48 -- iterate the pipeline's own MEM-removal matcher on its own
// residual literal until the gain drops below a threshold.
//
// The matcher in stage 47 is a single greedy left-to-right parse with a
// candidate cap (MAXCAND). Both are known sources of missed matches in the LZ
// literature: greedy parsing is provably non-optimal (a match taken now can
// block a longer one just past it), and a cap can hide the best candidate.
// Re-running the identical matcher on its own residual recovers exactly that
// class of miss. Verified: MAXCAND 64->512 barely moves the needle
// (+0.05%), so this is NOT a candidate-cap artifact -- it is the parse order.
//
// Usage: ./iter_mem literal.txt [MINMEM=24] [MAXCAND=64] [--dump]
// --dump writes literal2.txt (smaller residual) plus r2_gaps.bin/r2_srcs.bin/
// r2_lens.bin (the new reference stream, same vint format as stage 47's).
//
// A run bug was found and fixed while building this: an early version marked
// removed bytes at the raw reverse-complement-space position instead of the
// converted forward position. Byte COUNTS came out right by coincidence (a
// bijection preserves cardinality) so aggregate percentages looked plausible;
// only round-trip verification (content, not just length) caught it. The
// ACTUAL pipeline (46/47_*.cpp) does not have this bug -- confirmed by reading
// the exact call site, which uses a separate RC bitmap and a whole-array
// positional flip afterward, the correct conversion. This file has the fix.
//
//   g++ -O3 -march=native -pthread -o iter_mem 48_iterate_mem.cpp
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <string>
static void rc_inplace(std::string& s){
    for(char& c:s) c = c=='A'?'T':c=='C'?'G':c=='G'?'C':c=='T'?'A':'N';
    std::reverse(s.begin(),s.end());
}
static inline int b2(char c){ return c=='A'?0:c=='C'?1:c=='G'?2:c=='T'?3:-1; }

struct Ref { uint32_t dst, src, len; };

size_t self_rc_match(const std::string& pg, size_t MINMEM, size_t MEMSEED, size_t MAXCAND,
                      bool LAZY, std::vector<Ref>& allrefs, std::vector<uint8_t>& consumed){
    const size_t qlen=pg.size();
    const size_t STEP=MINMEM-MEMSEED+1;
    auto packM=[&](const char* q,uint64_t& o)->bool{
        uint64_t k=0; for(size_t i=0;i<MEMSEED;++i){ int v=b2(q[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
        o=k; return true; };
    std::vector<std::pair<uint64_t,uint32_t>> tmp;
    tmp.reserve(qlen/STEP+16);
    for(size_t p=0;p+MEMSEED<=qlen;p+=STEP){ uint64_t k; if(packM(pg.data()+p,k)) tmp.push_back({k,(uint32_t)p}); }
    std::sort(tmp.begin(),tmp.end());
    std::vector<uint64_t> skey(tmp.size()); std::vector<uint32_t> spos(tmp.size());
    for(size_t i=0;i<tmp.size();++i){ skey[i]=tmp[i].first; spos[i]=tmp[i].second; }
    size_t tsize=1; while(tsize<skey.size()*2+1) tsize<<=1; const uint64_t TMASK=tsize-1;
    std::vector<uint32_t> htab(tsize,UINT32_MAX);
    auto hmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
    for(size_t i=0;i<skey.size();++i){ if(i&&skey[i]==skey[i-1]) continue;
        size_t h=hmix(skey[i])&TMASK; while(htab[h]!=UINT32_MAX) h=(h+1)&TMASK; htab[h]=(uint32_t)i; }
    auto lookup=[&](uint64_t k)->uint32_t{ size_t h=hmix(k)&TMASK;
        while(htab[h]!=UINT32_MAX){ if(skey[htab[h]]==k) return htab[h]; h=(h+1)&TMASK; } return UINT32_MAX; };

    std::string R=pg; rc_inplace(R);
    auto bestAt=[&](size_t qp,size_t& bsrc)->size_t{
        uint64_t k; bsrc=0; if(!packM(R.data()+qp,k)) return 0;
        const uint32_t idx=lookup(k); if(idx==UINT32_MAX) return 0;
        size_t best=0,tried=0;
        for(uint32_t i=idx;i<skey.size()&&skey[i]==k;++i){
            if(++tried>MAXCAND) break;
            const size_t s=spos[i]; if(s>=qlen-qp) continue;
            const size_t capL=(qlen-qp-s)/2; if(capL<MINMEM) continue;
            size_t L=0; while(L<capL && qp+L<qlen && s+L<qlen && R[qp+L]==pg[s+L]) ++L;
            if(L>best){ best=L; bsrc=s; }
        }
        return (best>=MINMEM)?best:0;
    };
    auto parse_range=[&](size_t lo,size_t hi,std::vector<Ref>& out){
        size_t qp=lo,lastend=lo;
        while(qp<hi && qp+MINMEM<=qlen){
            uint64_t k; if(!packM(R.data()+qp,k)){ ++qp; continue; }
            const uint32_t idx=lookup(k); if(idx==UINT32_MAX){ ++qp; continue; }
            size_t best=0,bestsrc=0,tried=0;
            for(uint32_t i=idx;i<skey.size()&&skey[i]==k;++i){
                if(++tried>MAXCAND) break;
                const size_t s=spos[i]; if(s>=qlen-qp) continue;
                const size_t capL=(qlen-qp-s)/2; if(capL<MINMEM) continue;
                size_t L=0; while(L<capL && qp+L<qlen && s+L<qlen && R[qp+L]==pg[s+L]) ++L;
                if(L>best){ best=L; bestsrc=s; }
            }
            if(best>=MINMEM && LAZY && qp+1<hi && qp+1+MINMEM<=qlen){
                size_t nsrc=0; const size_t nb=bestAt(qp+1,nsrc);
                if(nb>best){ ++qp; continue; }
            }
            if(best>=MINMEM){
                size_t cap=qp-lastend; if(bestsrc<cap) cap=bestsrc;
                size_t b=0; while(b<cap && R[qp-b-1]==pg[bestsrc-b-1]) ++b;
                out.push_back({(uint32_t)(qp-b),(uint32_t)(bestsrc-b),(uint32_t)(best+b)});
                lastend=qp+best; qp+=best;
            } else ++qp;
        }
    };
    unsigned T=std::thread::hardware_concurrency(); if(!T) T=1; if(qlen<(1u<<20)) T=1;
    std::vector<std::vector<Ref>> res(T); std::vector<std::thread> th;
    const size_t chunk=(qlen+T-1)/T;
    for(unsigned t=0;t<T;++t){ const size_t lo=(size_t)t*chunk,hi=std::min(qlen,lo+chunk);
        if(lo>=hi) break; th.emplace_back([&,t,lo,hi]{ parse_range(lo,hi,res[t]); }); }
    for(auto& x:th) x.join();
    size_t nm=0;
    for(auto& v:res){ nm+=v.size();
        for(auto& m:v){
            Ref r=m; r.dst=(uint32_t)(qlen-m.dst-m.len);      // R-space qp -> forward position
            for(uint32_t j=0;j<r.len&&(size_t)r.dst+j<qlen;++j) consumed[r.dst+j]=1;   // mark FORWARD position
            allrefs.push_back(r); } }
    return nm;
}

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s <text> [MINMEM=24] [MAXCAND=64]\n",argv[0]); return 1; }
    FILE* f=fopen(argv[1],"rb"); fseek(f,0,SEEK_END); size_t n=ftell(f); fseek(f,0,SEEK_SET);
    std::string T(n,0); if(fread(&T[0],1,n,f)!=n) return 1; fclose(f);
    const size_t MINMEM=(argc>2)?(size_t)atoi(argv[2]):24;
    const size_t MEMSEED=(MINMEM>14)?((MINMEM-14>32)?32:MINMEM-14):8;
    const size_t MAXCAND=(argc>3)?(size_t)atoi(argv[3]):64;
    std::vector<Ref> refs; std::vector<uint8_t> consumed(n,0);
    size_t nm=self_rc_match(T,MINMEM,MEMSEED,MAXCAND,true,refs,consumed);
    size_t rem=0; for(auto c:consumed) rem+=c;
    printf("input=%zu  removed=%zu (%.3f%%)  matches=%zu  literal-out=%zu\n",
           n,rem,100.0*rem/n,nm,n-rem);
    if(argc>4 && std::string(argv[4])=="--dump"){
        FILE* lf=fopen("literal2.txt","wb");
        for(size_t i=0;i<n;++i) if(!consumed[i]) fputc(T[i],lf);
        fclose(lf);
        auto vint=[](std::vector<uint8_t>& o,uint64_t v){
            while(true){ uint8_t b=v&0x7f; v>>=7; o.push_back(b|(v?0x80:0)); if(!v) break; } };
        std::sort(refs.begin(),refs.end(),[](const Ref&a,const Ref&b){ return a.dst<b.dst; });
        std::vector<uint8_t> gaps,srcs,lens; uint64_t prev=0;
        for(auto& r:refs){ const uint64_t g=(r.dst>=prev)?(r.dst-prev):0;
            vint(gaps,g); vint(srcs,r.src); vint(lens,(uint64_t)(r.len-MINMEM)); prev=(uint64_t)r.dst+r.len; }
        FILE* g=fopen("r2_gaps.bin","wb"); fwrite(gaps.data(),1,gaps.size(),g); fclose(g);
        FILE* sf=fopen("r2_srcs.bin","wb"); fwrite(srcs.data(),1,srcs.size(),sf); fclose(sf);
        FILE* lnf=fopen("r2_lens.bin","wb"); fwrite(lens.data(),1,lens.size(),lnf); fclose(lnf);
        fprintf(stderr,"[DUMP] literal2.txt (%zu B), refs gaps=%zu srcs=%zu lens=%zu\n",
                n-rem,gaps.size(),srcs.size(),lens.size());
    }
    return 0;
}
