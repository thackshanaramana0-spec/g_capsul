// Greedy SCS + pigeonhole mapping, with the overlap search done in one pass.
//
// The earlier prototypes rebuilt a hash index once per candidate overlap length
// and re-hashed every read's L bytes each time -- O(n*L^2), about 19 billion byte
// reads on 851k reads, which is the entire reason they took 40 s while PgRC2's
// sort-and-merge takes 1.5 s. Nothing about the ALGORITHM was slow; the indexing
// was.
//
// Here each read is scanned once. A 32-base seed packs exactly into a uint64
// (2 bits/base), so the seed comparison is integer equality with no false
// positives from hashing, and sliding the seed one base is a shift-or. For a
// read A we walk its suffix offsets outward from the longest overlap; the first
// offset whose seed hits a read B and whose full suffix-prefix span verifies IS
// A's longest overlap, so we stop there. That makes the search O(n * offsets
// actually tried) instead of O(n * L^2).
//
// Each read keeps up to MAXCAND partners, not just its single best. Greedy takes
// the globally longest overlap first, so by the time a read is reached its top
// choice is often already consumed -- with one candidate stored that read links
// to nothing at all (measured: 700,160 reads found a partner but only 628,143
// links committed, so 72,017 reads were dropped outright). Backups let those
// reads fall through to their next-longest overlap, which is what the classic
// per-length sweep does implicitly by revisiting every open end at every length.
//
// Stage C then maps whatever never chained into the finished pseudogenome using
// the pigeonhole lemma (a read with <= m mismatches cut into m+1 parts must have
// one part matching exactly), forward and then against the reverse-complemented
// text -- one index, one text, rather than a two-strand index.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <unordered_map>
#include <algorithm>

static const uint32_t NONE    = UINT32_MAX;
static const uint32_t SEED    = 32;               // exactly one uint64 at 2 bits/base
static const uint64_t SEEDMSK = ~0ULL;            // 32 bases * 2 bits = 64

static inline int b2(char c){
    switch(c){ case 'A':return 0; case 'C':return 1; case 'G':return 2; case 'T':return 3; }
    return -1; }
static inline uint64_t fnv(const char* p,uint32_t n){
    uint64_t h=1469598103934665603ULL;
    for(uint32_t i=0;i<n;++i){ h^=(uint8_t)p[i]; h*=1099511628211ULL; } return h; }
static void rc_inplace(std::string& s){
    auto comp=[](char c){ switch(c){case 'A':return 'T';case 'C':return 'G';
                                    case 'G':return 'C';case 'T':return 'A';} return c; };
    if(s.empty()) return;
    size_t i=0,j=s.size()-1;
    while(i<j){ char a=comp(s[i]),b=comp(s[j]); s[i]=b; s[j]=a; ++i; --j; }
    if(i==j) s[i]=comp(s[i]);
}
// pack the 32 bases at p; false if any non-ACGT
static inline bool pack(const char* p,uint64_t& out){
    uint64_t k=0;
    for(uint32_t i=0;i<SEED;++i){ int v=b2(p[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
    out=k; return true;
}

int main(int argc,char** argv){
    if(argc<2){ fprintf(stderr,"usage: scs5 <in.fq> [maxmm] [minov]\n"); return 1; }
    const int      MAXMM = argc>2?atoi(argv[2]):3;
    const uint32_t MINOV = argc>3?(uint32_t)atoi(argv[3]):40;   // >= SEED
    auto T0=std::chrono::steady_clock::now();
    auto lap=[&](const char* w){ auto t=std::chrono::steady_clock::now();
        fprintf(stderr,"  %-26s %6.2f s\n",w,std::chrono::duration<double>(t-T0).count()); T0=t; };

    // ── load, drop N-reads, collapse exact duplicates ────────────────────────
    std::vector<std::string> reads;
    size_t n_in=0,n_filt=0;
    {
        std::ifstream f(argv[1]); std::string a,b,c,d;
        std::unordered_map<uint64_t,std::vector<uint32_t>> seen; seen.reserve(1u<<21);
        while(std::getline(f,a)&&std::getline(f,b)&&std::getline(f,c)&&std::getline(f,d)){
            ++n_in;
            if(b.find('N')!=std::string::npos){ ++n_filt; continue; }
            auto& bk=seen[fnv(b.data(),(uint32_t)b.size())]; bool dup=false;
            for(uint32_t id:bk) if(reads[id]==b){ dup=true; break; }
            if(!dup){ bk.push_back((uint32_t)reads.size()); reads.push_back(b); }
        }
    }
    const uint32_t n=(uint32_t)reads.size();
    uint32_t Lmax=0; for(auto& r:reads) Lmax=std::max(Lmax,(uint32_t)r.size());
    fprintf(stderr,"reads in=%zu N-filtered=%zu unique=%u maxlen=%u\n",n_in,n_filt,n,Lmax);
    lap("load+filter+dedup");

    // ── index every read by the 32 bases at its start ────────────────────────
    std::unordered_map<uint64_t,std::vector<uint32_t>> pref;
    pref.reserve(n*2);
    for(uint32_t i=0;i<n;++i){
        if(reads[i].size()<SEED) continue;
        uint64_t k; if(pack(reads[i].data(),k)) pref[k].push_back(i);
    }
    lap("prefix seed index");

    // ── longest verified suffix-prefix partner for each read, one pass ───────
    struct Cand { uint32_t a,b,ov; };
    std::vector<Cand> cands; cands.reserve((size_t)n*2);
    const uint32_t MAXCAND=4;        // backups per read
    const uint32_t EXTRA=24;         // offsets to keep walking after the first hit
    size_t tried=0;
    for(uint32_t a=0;a<n;++a){
        const std::string& A=reads[a];
        if(A.size()<MINOV) continue;
        const uint32_t maxoff=(uint32_t)A.size()-MINOV;
        uint64_t k=0; bool have=false; uint32_t got=0, stop=maxoff;
        for(uint32_t off=1;off<=stop;++off){
            if(off+SEED>A.size()) break;
            if(!have){ if(!pack(A.data()+off,k)){ have=false; continue; } have=true; }
            else{
                int v=b2(A[off+SEED-1]);
                if(v<0){ have=false; continue; }
                k=((k<<2)|(uint64_t)v)&SEEDMSK;
            }
            ++tried;
            auto it=pref.find(k); if(it==pref.end()) continue;
            const uint32_t ov=(uint32_t)A.size()-off;
            for(uint32_t b:it->second){
                if(b==a) continue;
                const std::string& B=reads[b];
                if(B.size()<ov) continue;
                if(memcmp(A.data()+off,B.data(),ov)!=0) continue;
                cands.push_back({a,b,ov});
                if(++got==1) stop=std::min(maxoff,off+EXTRA);   // bound the extra walk
                if(got>=MAXCAND) break;
            }
            if(got>=MAXCAND) break;
        }
    }
    fprintf(stderr,"seed probes=%zu  candidate links=%zu\n",tried,cands.size());
    lap("overlap search (multi-cand)");

    // ── greedy: every candidate link, globally longest overlap first ─────────
    std::sort(cands.begin(),cands.end(),
              [](const Cand& x,const Cand& y){ return x.ov>y.ov; });
    std::vector<uint32_t> nxt(n,NONE),prv(n,NONE),ovl(n,0),ch_h(n),ch_t(n);
    for(uint32_t i=0;i<n;++i){ ch_h[i]=i; ch_t[i]=i; }
    size_t links=0;
    for(const Cand& c:cands){
        if(nxt[c.a]!=NONE||prv[c.b]!=NONE||ch_h[c.a]==c.b) continue;
        nxt[c.a]=c.b; prv[c.b]=c.a; ovl[c.a]=c.ov;
        uint32_t h=ch_h[c.a],t=ch_t[c.b]; ch_t[h]=t; ch_h[t]=h; ++links;
    }
    lap("greedy link");

    // ── emit chains; singletons held back for mapping ────────────────────────
    std::string pg; pg.reserve((size_t)n*40);
    std::vector<uint32_t> leftovers;
    uint32_t multi=0;
    for(uint32_t i=0;i<n;++i){
        if(prv[i]!=NONE) continue;
        if(nxt[i]==NONE){ leftovers.push_back(i); continue; }
        ++multi;
        uint32_t cur=i; pg+=reads[cur];
        while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                               pg.append(reads[cur].data()+o,reads[cur].size()-o); }
    }
    fprintf(stderr,"links=%zu chains(multi)=%u leftovers=%zu pg after chains=%zu\n",
            links,multi,leftovers.size(),pg.size());
    lap("emit chains");

    // ── Stage C: pigeonhole mapping into the pg, forward then RC ─────────────
    std::vector<uint8_t> matched(n,0);
    size_t n_matched=0;
    {
        const uint32_t NPARTS=(uint32_t)MAXMM+1;
        std::unordered_map<uint64_t,std::vector<std::pair<uint32_t,uint8_t>>> seeds;
        seeds.reserve(leftovers.size()*NPARTS*2);
        for(uint32_t rid:leftovers){
            const std::string& R=reads[rid];
            if(R.size()<NPARTS*SEED) continue;
            for(uint32_t p=0;p<NPARTS;++p){
                uint64_t k; if(pack(R.data()+p*SEED,k)) seeds[k].push_back({rid,(uint8_t)p});
            }
        }
        auto scan=[&](const std::string& text){
            if(text.size()<SEED) return;
            uint64_t k=0; uint32_t good=0;
            for(size_t p=0;p<text.size();++p){
                int v=b2(text[p]);
                if(v<0){ good=0; k=0; continue; }
                k=((k<<2)|(uint64_t)v)&SEEDMSK; ++good;
                if(good<SEED) continue;
                auto it=seeds.find(k); if(it==seeds.end()) continue;
                const size_t seedStart=p-SEED+1;
                for(auto& pr:it->second){
                    const uint32_t rid=pr.first; if(matched[rid]) continue;
                    const size_t off=(size_t)pr.second*SEED;
                    if(seedStart<off) continue;
                    const size_t st=seedStart-off;
                    const std::string& R=reads[rid];
                    if(st+R.size()>text.size()) continue;
                    int mm=0; const char* t=text.data()+st;
                    for(size_t j=0;j<R.size();++j) if(t[j]!=R[j]){ if(++mm>MAXMM) break; }
                    if(mm<=MAXMM){ matched[rid]=1; ++n_matched; }
                }
            }
        };
        scan(pg);
        rc_inplace(pg); scan(pg); rc_inplace(pg);
    }
    lap("pigeonhole mapping");

    size_t appended=0;
    for(uint32_t rid:leftovers) if(!matched[rid]){ pg+=reads[rid]; ++appended; }
    fprintf(stderr,"leftovers=%zu mapped=%zu appended=%zu\n",
            leftovers.size(),n_matched,appended);
    printf("PG_LEN %zu\n",pg.size());
    return 0;
}
