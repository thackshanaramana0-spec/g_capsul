// Greedy SCS with MINIMIZER-SEEDED, MISMATCH-TOLERANT overlap detection.
//
// The point of the experiment: PgRC's HQ/LQ division exists because its overlap
// test is exact, and one sequencing error destroys an exact overlap — so reads
// carrying errors must be quarantined by quality score and assembled separately.
// That is a workaround for a brittle detector, and it costs twice: the LQ reads
// leave the main assembly, and the split depends on quality scores that modern
// platforms bin down to four levels.
//
// Here the detector itself tolerates errors. Minimizers (Roberts et al. 2004
// winnowing; Li 2016 minimap) sample a read sparsely, so one or two wrong bases
// leave most of the sample intact and the true partner is still found. Every
// read then participates in one assembly, and no quality score is consulted.
//
// Greedy is applied properly: score every read's best partner first, then link
// in globally descending overlap order — not the per-length sweep.
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
static inline uint64_t mix64(uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33;
                                          x*=0xc4ceb9fe1a85ec53ULL; x^=x>>33; return x; }
static inline uint64_t fnv(const char* p, uint32_t n){
    uint64_t h=1469598103934665603ULL;
    for(uint32_t i=0;i<n;++i){ h^=(uint8_t)p[i]; h*=1099511628211ULL; } return h; }
static inline int b2(char c){ switch(c){case 'A':return 0;case 'C':return 1;
                                        case 'G':return 2;case 'T':return 3;} return -1; }

int main(int argc, char** argv){
    if(argc<2){ fprintf(stderr,"usage: scs2 <in.fq> [k] [w] [maxmm] [minov]\n"); return 1; }
    const uint32_t K     = argc>2?(uint32_t)atoi(argv[2]):21;
    const uint32_t W     = argc>3?(uint32_t)atoi(argv[3]):10;
    const int      MAXMM = argc>4?atoi(argv[4]):3;
    const uint32_t MINOV = argc>5?(uint32_t)atoi(argv[5]):30;
    auto T0=std::chrono::steady_clock::now();
    auto lap=[&](const char*w){ auto t=std::chrono::steady_clock::now();
        fprintf(stderr,"  %-24s %6.2f s\n",w,std::chrono::duration<double>(t-T0).count()); T0=t; };

    // ── load, drop N reads, dedup (identical to the exact-greedy baseline) ───
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
    fprintf(stderr,"reads in=%zu N-filtered=%zu unique=%u  (k=%u w=%u maxmm=%d)\n",
            n_in,n_filt,n,K,W,MAXMM);
    lap("load+filter+dedup");

    // ── minimizer index: sparse sample, so a wrong base costs only its own seeds
    struct Post { uint32_t rid; uint16_t pos; };
    std::unordered_map<uint64_t,std::vector<Post>> idx; idx.reserve(1u<<22);
    for(uint32_t r=0;r<n;++r){
        const std::string& S=reads[r]; if(S.size()<K) continue;
        const uint32_t nk=(uint32_t)S.size()-K+1;
        uint64_t last=UINT64_MAX; uint32_t lastp=UINT32_MAX;
        for(uint32_t i=0;i+W<=nk;++i){                       // window of W k-mers
            uint64_t best=UINT64_MAX; uint32_t bp=0;
            for(uint32_t j=i;j<i+W;++j){
                uint64_t h=mix64(fnv(S.data()+j,K));
                if(h<best){ best=h; bp=j; }
            }
            if(best!=last||bp!=lastp){ idx[best].push_back({r,(uint16_t)bp}); last=best; lastp=bp; }
        }
    }
    size_t posts=0; for(auto& kv:idx) posts+=kv.second.size();
    fprintf(stderr,"minimizer index: %zu keys, %zu postings\n",idx.size(),posts);
    lap("minimizer index");

    // ── best tolerant suffix-prefix partner for every read ───────────────────
    struct Best { uint32_t rid; uint32_t ov; };
    std::vector<Best> best(n,{NONE,0});
    std::unordered_map<uint32_t,uint32_t> vote;                 // partner -> implied overlap
    const size_t BUCKET_CAP=64;                                 // ignore repeat-heavy seeds
    for(uint32_t a=0;a<n;++a){
        const std::string& A=reads[a]; if(A.size()<K) continue;
        vote.clear();
        const uint32_t nk=(uint32_t)A.size()-K+1;
        for(uint32_t i=0;i+W<=nk;++i){
            uint64_t bh=UINT64_MAX; uint32_t bp=0;
            for(uint32_t j=i;j<i+W;++j){ uint64_t h=mix64(fnv(A.data()+j,K)); if(h<bh){bh=h;bp=j;} }
            auto it=idx.find(bh); if(it==idx.end()||it->second.size()>BUCKET_CAP) continue;
            for(const Post& p:it->second){
                if(p.rid==a) continue;
                if(bp<p.pos) continue;                          // B must start inside A
                uint32_t start=bp-p.pos;                        // B's start in A's frame
                if(start==0) continue;                          // identical placement
                if(start>=A.size()) continue;
                uint32_t ov=(uint32_t)A.size()-start;           // suffix-prefix overlap
                if(ov<MINOV||ov>reads[p.rid].size()) continue;
                auto& v=vote[p.rid]; if(ov>v) v=ov;
            }
        }
        for(auto& kv:vote){
            uint32_t b=kv.first, ov=kv.second;
            if(ov<=best[a].ov) continue;
            const std::string& B=reads[b];
            int mm=0; const char* pa=A.data()+(A.size()-ov); const char* pb=B.data();
            for(uint32_t i=0;i<ov;++i) if(pa[i]!=pb[i]){ if(++mm>MAXMM) break; }
            if(mm<=MAXMM) best[a]={b,ov};                       // tolerant verification
        }
    }
    lap("tolerant overlap search");

    // ── proper greedy: link in globally descending overlap order ─────────────
    std::vector<uint32_t> order; order.reserve(n);
    for(uint32_t i=0;i<n;++i) if(best[i].rid!=NONE) order.push_back(i);
    std::sort(order.begin(),order.end(),
              [&](uint32_t x,uint32_t y){ return best[x].ov>best[y].ov; });
    std::vector<uint32_t> nxt(n,NONE),prv(n,NONE),ovl(n,0),ch_h(n),ch_t(n);
    for(uint32_t i=0;i<n;++i){ ch_h[i]=i; ch_t[i]=i; }
    size_t links=0;
    for(uint32_t a:order){
        uint32_t b=best[a].rid;
        if(nxt[a]!=NONE||prv[b]!=NONE||ch_h[a]==b) continue;    // O(1) cycle avoidance
        nxt[a]=b; prv[b]=a; ovl[a]=best[a].ov;
        uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h; ++links;
    }
    lap("greedy link");

    // ── emit; mismatches are absorbed by taking the predecessor's bases ──────
    std::string pg; pg.reserve((size_t)n*40); uint32_t chains=0;
    for(uint32_t i=0;i<n;++i){
        if(prv[i]!=NONE) continue;
        ++chains; uint32_t cur=i; pg+=reads[cur];
        while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                               pg.append(reads[cur].data()+o,reads[cur].size()-o); }
    }
    lap("emit");
    fprintf(stderr,"links=%zu chains=%u\n",links,chains);
    printf("PG_LEN %zu\n",pg.size());
    return 0;
}
