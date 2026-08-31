// Greedy SCS + pigeonhole seed-and-extend mapping, written from the algorithm.
//
// Stage A  greedy exact-overlap chaining              (Tarhio-Ukkonen)
// Stage B  emit a pseudogenome from multi-read chains
// Stage C  map every leftover read INTO that pg, tolerating mismatches, forward
//          and reverse-complement; only genuinely unmappable reads are appended
//
// Stage C is the piece the earlier prototypes lacked, and it is textbook: the
// q-gram / pigeonhole lemma. A read of length L carrying at most m mismatches,
// cut into m+1 equal parts, must have at least one part matching EXACTLY. So
// index the parts, stream the pg once looking for exact part hits, and verify
// each candidate placement over the whole read. Seeds are 32 bases so a part
// packs into one uint64 and the scan rolls in O(1) per position.
//
// Reverse strand is handled the way PgRC does it -- reverse-complement the pg
// in place and run the same scan again -- rather than by indexing both strands.
// One index, one text, half the memory of a two-view index.
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

static const uint32_t NONE = UINT32_MAX;
static const uint32_t SEEDLEN = 32;                 // one 2-bit-packed uint64

static inline uint64_t fnv(const char* p, uint32_t n){
    uint64_t h=1469598103934665603ULL;
    for(uint32_t i=0;i<n;++i){ h^=(uint8_t)p[i]; h*=1099511628211ULL; } return h; }
static inline int b2(char c){
    switch(c){ case 'A':return 0; case 'C':return 1; case 'G':return 2; case 'T':return 3; }
    return -1; }
static void rc_inplace(std::string& s){
    size_t i=0, j=s.size()?s.size()-1:0;
    auto comp=[](char c){ switch(c){case 'A':return 'T';case 'C':return 'G';
                                    case 'G':return 'C';case 'T':return 'A';} return c; };
    while(i<j){ char a=comp(s[i]), b=comp(s[j]); s[i]=b; s[j]=a; ++i; --j; }
    if(i==j) s[i]=comp(s[i]);
}

int main(int argc,char** argv){
    if(argc<2){ fprintf(stderr,"usage: scs3 <in.fq> [maxmm] [minov]\n"); return 1; }
    const int      MAXMM = argc>2?atoi(argv[2]):3;
    const uint32_t MINOV = argc>3?(uint32_t)atoi(argv[3]):30;
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
    // Reverse-complement view of every read. A shotgun library samples both
    // strands, so chaining only forward-suffix to forward-prefix assembles two
    // half-coverage datasets that can never join. Each read may enter a chain
    // in either orientation, and is consumed once whichever it enters in.
    std::vector<std::string> rcv(n);
    for(uint32_t i=0;i<n;++i){ rcv[i]=reads[i]; rc_inplace(rcv[i]); }
    auto view=[&](uint32_t i,uint8_t v)->const std::string&{ return v?rcv[i]:reads[i]; };
    std::vector<uint8_t> vw(n,0);   // orientation each read occupies in its chain
    fprintf(stderr,"reads in=%zu N-filtered=%zu unique=%u maxlen=%u\n",n_in,n_filt,n,Lmax);
    lap("load+filter+dedup");

    // ── Stage A: greedy exact-overlap chaining, longest overlap first ────────
    std::vector<uint32_t> nxt(n,NONE),prv(n,NONE),ovl(n,0),ch_h(n),ch_t(n);
    for(uint32_t i=0;i<n;++i){ ch_h[i]=i; ch_t[i]=i; }
    {
        std::vector<uint32_t> tails(n),heads(n);
        for(uint32_t i=0;i<n;++i){ tails[i]=i; heads[i]=i; }
        std::unordered_map<uint64_t,std::vector<uint32_t>> idx;
        for(uint32_t L=Lmax-1;L>=MINOV;--L){
            heads.erase(std::remove_if(heads.begin(),heads.end(),[&](uint32_t i){
                return prv[i]!=NONE||reads[i].size()<L; }),heads.end());
            tails.erase(std::remove_if(tails.begin(),tails.end(),[&](uint32_t i){
                return nxt[i]!=NONE||reads[i].size()<L; }),tails.end());
            if(heads.empty()||tails.empty()) break;
            // index the L-prefix of BOTH views of every still-open read
            idx.clear(); idx.reserve(heads.size()*4);
            for(uint32_t b:heads){
                idx[fnv(reads[b].data(),L)].push_back(b<<1);        // forward view
                idx[fnv(rcv[b].data(),L)].push_back((b<<1)|1);      // RC view
            }
            for(uint32_t a:tails){
                if(nxt[a]!=NONE) continue;
                const std::string& A=view(a,vw[a]);                 // a's own orientation
                auto it=idx.find(fnv(A.data()+(A.size()-L),L)); if(it==idx.end()) continue;
                for(uint32_t code:it->second){
                    const uint32_t b=code>>1; const uint8_t vb=(uint8_t)(code&1);
                    if(b==a||prv[b]!=NONE||ch_h[a]==b) continue;
                    const std::string& B=view(b,vb);
                    if(memcmp(A.data()+(A.size()-L),B.data(),L)!=0) continue;
                    nxt[a]=b; prv[b]=a; ovl[a]=L; vw[b]=vb;
                    uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h;
                    break;
                }
            }
        }
    }
    lap("greedy overlap sweep");

    // ── Stage B: pg from multi-read chains; singletons held back ─────────────
    std::string pg; pg.reserve((size_t)n*40);
    std::vector<uint32_t> leftovers;
    uint32_t multi=0;
    for(uint32_t i=0;i<n;++i){
        if(prv[i]!=NONE) continue;
        if(nxt[i]==NONE){ leftovers.push_back(i); continue; }   // chain of one
        ++multi;
        uint32_t cur=i; pg+=view(cur,vw[cur]);
        while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                               const std::string& V=view(cur,vw[cur]);
                               pg.append(V.data()+o,V.size()-o); }
    }
    fprintf(stderr,"chains(multi)=%u  leftover singletons=%zu  pg after stage B=%zu\n",
            multi,leftovers.size(),pg.size());
    lap("emit multi-read chains");

    // ── Stage C: pigeonhole map leftovers into pg (forward, then RC) ─────────
    std::vector<uint8_t> matched(n,0);
    size_t n_matched=0, n_false=0;
    {
        // index every leftover read's m+1 seed parts
        std::unordered_map<uint64_t,std::vector<std::pair<uint32_t,uint8_t>>> seeds;
        seeds.reserve(leftovers.size()*(size_t)(MAXMM+1)*2);
        const uint32_t NPARTS=(uint32_t)MAXMM+1;
        for(uint32_t rid:leftovers){
            const std::string& R=reads[rid];
            if(R.size()<NPARTS*SEEDLEN) continue;
            for(uint32_t p=0;p<NPARTS;++p){
                bool ok=true; uint64_t k=0;
                for(uint32_t j=0;j<SEEDLEN;++j){
                    int v=b2(R[p*SEEDLEN+j]); if(v<0){ ok=false; break; }
                    k=(k<<2)|(uint64_t)v;
                }
                if(ok) seeds[k].push_back({rid,(uint8_t)p});
            }
        }
        fprintf(stderr,"seed index: %zu keys\n",seeds.size());

        auto scan=[&](const std::string& text,bool rc){
            if(text.size()<SEEDLEN) return;
            const uint64_t MASK=(SEEDLEN==32)?~0ULL:((1ULL<<(2*SEEDLEN))-1);
            uint64_t k=0; uint32_t good=0;
            for(size_t p=0;p<text.size();++p){
                int v=b2(text[p]);
                if(v<0){ good=0; k=0; continue; }
                k=((k<<2)|(uint64_t)v)&MASK; ++good;
                if(good<SEEDLEN) continue;
                auto it=seeds.find(k); if(it==seeds.end()) continue;
                const size_t seedStart=p-SEEDLEN+1;
                for(auto& pr:it->second){
                    const uint32_t rid=pr.first; if(matched[rid]) continue;
                    const size_t off=(size_t)pr.second*SEEDLEN;
                    if(seedStart<off) continue;
                    const size_t st=seedStart-off;
                    const std::string& R=reads[rid];
                    if(st+R.size()>text.size()) continue;
                    int mm=0; const char* t=text.data()+st;
                    for(size_t j=0;j<R.size();++j) if(t[j]!=R[j]){ if(++mm>MAXMM) break; }
                    if(mm<=MAXMM){ matched[rid]=1; ++n_matched; }
                    else ++n_false;
                }
            }
        };
        scan(pg,false);
        fprintf(stderr,"  forward pass: matched=%zu\n",n_matched);
        rc_inplace(pg);                       // same text, other strand
        scan(pg,true);
        rc_inplace(pg);                       // restore
        fprintf(stderr,"  after RC pass: matched=%zu (false verifications=%zu)\n",n_matched,n_false);
    }
    lap("pigeonhole mapping");

    // ── append only what genuinely could not be placed ───────────────────────
    size_t appended=0;
    for(uint32_t rid:leftovers) if(!matched[rid]){ pg+=reads[rid]; ++appended; }
    lap("append unmapped");

    fprintf(stderr,"leftovers=%zu mapped=%zu appended=%zu\n",
            leftovers.size(),n_matched,appended);
    printf("PG_LEN %zu\n",pg.size());
    return 0;
}
