// STAGE 49 -- forward self-match, run AFTER stage 48's RC-only convergence.
//
// Stage 36 found forward self-matching a net LOSS on a single pass over the
// raw pg (extra literal saved does not cover the reference cost). It is a net
// WIN here, on the RC-CONVERGED residual: the residual's literal/reference
// economics genuinely differ from the raw pg's, because everything the RC
// pass could remove is already gone -- both conclusions are correct, scoped
// to different inputs, and stage 36's is not overturned by this.
//
// Usage: ./fwd_mem residual.txt [MINMEM=24] [MAXCAND=64]
// Always writes fwd_literal.txt, fwd_gaps.bin, fwd_srcs.bin, fwd_lens.bin.
//
//   g++ -O3 -march=native -pthread -o fwd_mem 49_fwd_after_converge.cpp

// Does forward self-matching pay off on the RC-converged residual? Stage 36
// found FWD self-match a net loss in one pass (extra literal saved does not
// cover the reference cost). Retest after RC has already been iterated to
// convergence -- the residual's economics may differ.
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
static inline int b2(char c){ return c=='A'?0:c=='C'?1:c=='G'?2:c=='T'?3:-1; }
struct Ref { uint32_t dst, src, len; };
int main(int argc,char**argv){
    FILE* f=fopen(argv[1],"rb"); fseek(f,0,SEEK_END); size_t n=ftell(f); fseek(f,0,SEEK_SET);
    std::string T(n,0); if(fread(&T[0],1,n,f)!=n) return 1; fclose(f);
    const size_t MINMEM=24, MEMSEED=10, STEP=15, MAXCAND=64;
    auto packM=[&](const char* q,uint64_t& o)->bool{
        uint64_t k=0; for(size_t i=0;i<MEMSEED;++i){ int v=b2(q[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
        o=k; return true; };
    std::vector<std::pair<uint64_t,uint32_t>> tmp;
    for(size_t p=0;p+MEMSEED<=n;p+=STEP){ uint64_t k; if(packM(T.data()+p,k)) tmp.push_back({k,(uint32_t)p}); }
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
    std::vector<uint8_t> consumed(n,0); std::vector<Ref> refs; size_t lastend=0;
    for(size_t qp=0;qp+MINMEM<=n; ){
        uint64_t k; if(!packM(T.data()+qp,k)){ ++qp; continue; }
        const uint32_t idx=lookup(k); if(idx==UINT32_MAX){ ++qp; continue; }
        size_t best=0,bestsrc=0,tried=0;
        for(uint32_t i=idx;i<skey.size()&&skey[i]==k;++i){
            if(++tried>MAXCAND) break;
            const size_t s=spos[i]; if(s>=qp) continue;      // FORWARD self-match: source precedes dest
            const size_t capL=qp-s; if(capL<MINMEM) continue;
            size_t L=0; while(L<capL && qp+L<n && T[qp+L]==T[s+L]) ++L;
            if(L>best){ best=L; bestsrc=s; }
        }
        if(best>=MINMEM){
            size_t cap=qp-lastend; if(bestsrc<cap) cap=bestsrc;
            const size_t slack=qp-bestsrc-best; if(slack<cap) cap=slack;
            size_t b=0; while(b<cap && T[qp-b-1]==T[bestsrc-b-1]) ++b;
            for(size_t j=0;j<best+b;++j) consumed[qp-b+j]=1;
            refs.push_back({(uint32_t)(qp-b),(uint32_t)(bestsrc-b),(uint32_t)(best+b)});
            lastend=qp+best; qp+=best;
        } else ++qp;
    }
    size_t rem=0; for(auto c:consumed) rem+=c;
    FILE* lf=fopen("fwd_literal.txt","wb");
    for(size_t i=0;i<n;++i) if(!consumed[i]) fputc(T[i],lf);
    fclose(lf);
    std::sort(refs.begin(),refs.end(),[](const Ref&a,const Ref&b){return a.dst<b.dst;});
    auto vint=[](std::vector<uint8_t>& o,uint64_t v){ while(true){ uint8_t b=v&0x7f; v>>=7; o.push_back(b|(v?0x80:0)); if(!v) break; } };
    std::vector<uint8_t> gaps,srcs,lens; uint64_t prev=0;
    for(auto& r:refs){ vint(gaps,(r.dst>=prev)?(r.dst-prev):0); vint(srcs,r.src); vint(lens,(uint64_t)(r.len-MINMEM)); prev=(uint64_t)r.dst+r.len; }
    FILE* g=fopen("fwd_gaps.bin","wb"); fwrite(gaps.data(),1,gaps.size(),g); fclose(g);
    FILE* sf=fopen("fwd_srcs.bin","wb"); fwrite(srcs.data(),1,srcs.size(),sf); fclose(sf);
    FILE* lnf=fopen("fwd_lens.bin","wb"); fwrite(lens.data(),1,lens.size(),lnf); fclose(lnf);
    printf("input=%zu removed=%zu (%.3f%%) matches=%zu literal-out=%zu\n",n,rem,100.0*rem/n,refs.size(),n-rem);
    return 0;
}
