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
#include <thread>
#include <unistd.h>

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
    const uint32_t MINOV = argc>3?(uint32_t)atoi(argv[3]):40;
    // Sweep seed width. 32 bases fills a uint64 exactly, which is why it was
    // fixed there -- but it also silently floors the shortest overlap the sweep
    // can ever see, and the measured trend (MINOV 50 -> 40 -> 32 giving 29.88M
    // -> 28.91M -> 28.28M) runs straight into that floor. A narrower seed costs
    // specificity (more candidates to verify) but lets round 2 pair reads that
    // overlap by less than 32, which after the division are still reads known
    // to tile well.
    const uint32_t SW = argc>4?(uint32_t)atoi(argv[4]):32;
    const uint64_t SWMASK = (SW>=32)?~0ULL:((1ULL<<(2*SW))-1);
    auto packSW=[&](const char* p,uint64_t& out)->bool{
        uint64_t k=0;
        for(uint32_t i=0;i<SW;++i){ int v=b2(p[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
        out=k; return true; };   // >= SEED
    auto T0=std::chrono::steady_clock::now();
    auto rss_mb=[]()->size_t{ FILE* f=fopen("/proc/self/statm","r"); if(!f) return 0;
        long pg_=0,res=0; if(fscanf(f,"%ld %ld",&pg_,&res)!=2) res=0; fclose(f);
        return (size_t)((double)res*(double)sysconf(_SC_PAGESIZE)/1048576.0); };
    auto lap=[&](const char* w){ auto t=std::chrono::steady_clock::now();
        fprintf(stderr,"  %-26s %6.2f s   rss=%zu MB\n",w,
                std::chrono::duration<double>(t-T0).count(),rss_mb()); T0=t; };

    // ── load, drop N-reads, collapse exact duplicates ────────────────────────
    // ── STAGE 21: reads stored 2 bits/base ───────────────────────────────────
    // Profiled, peak RSS hit 224 MB at load, before any index existed -- 96% of
    // PgRC2's entire 232 MB peak, spent purely on holding the reads. 851,275
    // std::strings pay a 32-byte header plus allocator rounding on every 150
    // bases, ~192 B for 150 B of sequence.
    //
    // PgRC2 does not do that: SymbolsPackingFacility keeps reads packed and
    // compares packed against unpacked in place (countSequenceMismatchesVsUnpacked,
    // compareSequenceWithUnpacked). At 4 bases/byte the same reads are ~34 MB.
    //
    // Safe here because N-containing reads are already filtered out above, so
    // the alphabet really is 2-bit. Each read starts on a word boundary: it
    // wastes under a word per read and makes a prefix load aligned.
    std::vector<uint64_t> rpk;            // 32 bases per uint64, base 0 in the top bits
    std::vector<uint64_t> woff;           // starting WORD of each read
    std::vector<uint8_t>  rlen;           // length in bases
    // STAGE 22 (their stage 6, OrderInfo): to restore the original file order we
    // need, for every ORIGINAL read, where its sequence sits in the pg. Dedup
    // collapses duplicates, so the original->unique map has to be kept too.
    std::vector<uint32_t> orig2uid;       // original read (N-filtered out) -> unique id
    size_t n_in=0,n_filt=0;
    {
        std::ifstream f(argv[1]); std::string a,b,c,d;
        std::unordered_map<uint64_t,std::vector<uint32_t>> seen; seen.reserve(1u<<21);
        std::vector<uint64_t> tmpw; char ubuf[256];
        auto packstr=[&](const char* q,uint32_t L,std::vector<uint64_t>& out){
            out.clear();
            for(uint32_t i=0;i<L;i+=32){
                uint64_t v=0;
                for(uint32_t j=0;j<32;++j){
                    const uint32_t idx=i+j;
                    const uint64_t c=(idx<L)?(uint64_t)b2(q[idx]):0ULL;
                    v=(v<<2)|c;
                }
                out.push_back(v);
            }
        };
        while(std::getline(f,a)&&std::getline(f,b)&&std::getline(f,c)&&std::getline(f,d)){
            ++n_in;
            if(b.find('N')!=std::string::npos){ ++n_filt; continue; }
            if(b.size()>255) continue;                       // uint8 length field
            auto& bk=seen[fnv(b.data(),(uint32_t)b.size())]; bool dup=false;
            for(uint32_t id:bk){
                if(rlen[id]!=(uint8_t)b.size()) continue;
                const uint64_t* w=&rpk[woff[id]];
                const uint32_t nw=((uint32_t)rlen[id]+31)/32;
                for(uint32_t t=0;t<nw;++t){
                    for(uint32_t j=0;j<32;++j){
                        const uint32_t idx=t*32+j;
                        ubuf[idx]="ACGT"[(w[t]>>(2*(31-j)))&3ULL];
                    }
                }
                if(memcmp(ubuf,b.data(),b.size())==0){ dup=true; orig2uid.push_back(id); break; }
            }
            if(!dup){
                orig2uid.push_back((uint32_t)rlen.size());
                bk.push_back((uint32_t)rlen.size());
                packstr(b.data(),(uint32_t)b.size(),tmpw);
                woff.push_back(rpk.size());
                rpk.insert(rpk.end(),tmpw.begin(),tmpw.end());
                rlen.push_back((uint8_t)b.size());
            }
        }
        rpk.shrink_to_fit(); woff.shrink_to_fit(); rlen.shrink_to_fit();
    }
    const uint32_t n=(uint32_t)rlen.size();
    uint32_t Lmax=0; for(uint32_t i=0;i<n;++i) Lmax=std::max(Lmax,(uint32_t)rlen[i]);

    // 32 bases starting at an arbitrary ABSOLUTE base index, base 0 in top bits.
    auto w32=[&](uint64_t ab)->uint64_t{
        const size_t wi=ab>>5; const uint32_t sh=(uint32_t)(ab&31)*2;
        const uint64_t hi=rpk[wi];
        if(!sh) return hi;
        const uint64_t lo=(wi+1<rpk.size())?rpk[wi+1]:0ULL;
        return (hi<<sh)|(lo>>(64-sh));
    };
    // read i's bases [o, o+SWb) as an integer, SWb <= 32
    auto rseed=[&](uint32_t i,uint32_t o,uint32_t SWb)->uint64_t{
        return w32(woff[i]*32ULL+o)>>(2*(32-SWb));
    };
    // compare read a's [off,off+L) against read b's [0,L)
    auto rcmp=[&](uint32_t a,uint32_t off,uint32_t b,uint32_t L)->bool{
        uint64_t pa=woff[a]*32ULL+off, pb=woff[b]*32ULL;
        while(L>=32){ if(w32(pa)!=w32(pb)) return false; pa+=32; pb+=32; L-=32; }
        if(L){ const uint64_t m=~0ULL<<(2*(32-L)); if((w32(pa)^w32(pb))&m) return false; }
        return true;
    };
    // unpack read i into buf (caller guarantees >= rlen[i] bytes)
    auto runp=[&](uint32_t i,char* buf){
        const uint32_t L=rlen[i]; const uint64_t base=woff[i]*32ULL;
        for(uint32_t j=0;j<L;++j) buf[j]="ACGT"[(w32(base+j)>>62)&3ULL];
    };
    auto rappend=[&](std::string& dst,uint32_t i,uint32_t from){
        char b_[256]; runp(i,b_); dst.append(b_+from,rlen[i]-from);
    };
    fprintf(stderr,"reads in=%zu N-filtered=%zu unique=%u maxlen=%u\n",n_in,n_filt,n,Lmax);
    lap("load+filter+dedup");

    // ── index every read by the 32 bases at its start ────────────────────────
    std::unordered_map<uint64_t,std::vector<uint32_t>> pref;
    pref.reserve(n*2);
    for(uint32_t i=0;i<n;++i){
        if(rlen[i]<SW) continue;
        pref[rseed(i,0,SW)].push_back(i);
    }
    lap("prefix seed index");

    // ── longest verified suffix-prefix partner for each read, one pass ───────
    // ── TRUE descending-length sweep ────────────────────────────────────────
    // Storing one best partner per read is what broke this: sorted by overlap
    // and committed, a read whose single candidate was already consumed links
    // to nothing, and every unlinked read starts a chain that costs a full read
    // length. Measured, that was 78,327 chains x 150 b = 11.7 MB of pure chain
    // heads, against PgRC2's ~5,100 chains on the same read count -- the whole
    // gap, and nothing to do with tolerance or strandedness.
    //
    // The sweep instead revisits EVERY still-open end at every overlap length,
    // longest first, so a read whose best partner is taken still gets its
    // next-best one length down. Naively that re-reads L bytes per read per
    // length (O(n*L^2), the 40 s version). Here each open tail carries a rolling
    // 32-base seed: its suffix offset grows by exactly one as L drops by one, so
    // the seed advances with a single shift-or, and the L-byte memcmp runs only
    // after an exact 64-bit seed hit.
    std::vector<uint64_t> ppos;           // pg position of each unique read (stage 22)
    std::vector<uint32_t> nxt(n,NONE),prv(n,NONE),ovl(n,0),ch_h(n),ch_t(n);
    std::vector<uint64_t> seed(n,0);
    std::vector<uint8_t>  ok(n,0);
    std::vector<uint32_t> tails;
    size_t links=0, probes=0;
    // `admit` selects which reads take part; everything else is reset per round.
    std::vector<uint8_t> admit(n,1);
    // Round 1 only LABELS reads (which ones sit interior to a chain), so its
    // links cost nothing -- a link there buys classification, not bytes. It
    // therefore runs at the widest overlap the seed allows, while round 2, which
    // actually emits sequence, keeps the caller's floor. Running both at the
    // same conservative floor was leaving reads unlabelled that would have
    // classified fine: interior fraction 66% against PgRC2's 83%, which starves
    // the main pseudogenome and leaves the mapping stage with less to hit.
    uint32_t sweep_minov=MINOV;
    auto sweep=[&](){
        links=0; probes=0;
        std::fill(nxt.begin(),nxt.end(),NONE); std::fill(prv.begin(),prv.end(),NONE);
        std::fill(ovl.begin(),ovl.end(),0);    std::fill(ok.begin(),ok.end(),0);
        for(uint32_t i=0;i<n;++i){ ch_h[i]=i; ch_t[i]=i; }
        tails.clear(); tails.reserve(n);
        for(uint32_t i=0;i<n;++i) if(admit[i]&&rlen[i]>=sweep_minov) tails.push_back(i);
    // STAGE 19. The sweep splits cleanly into a search and a commit.
    //
    // Everything expensive here -- the seed roll, the hash probe, and the
    // L-byte memcmp against every candidate -- reads only state that is fixed
    // for the whole of one L level: reads, pref, admit, and each tail's own
    // seed/ok slot. The only conditions that change WITHIN a level are
    // `prv[b]!=NONE` (b already claimed) and `ch_h[a]==b` (would close a cycle),
    // and both are a couple of array loads.
    //
    // So the search runs in parallel and the commit stays serial in the original
    // tail order. Because the parallel phase records EVERY candidate that passes
    // the static tests, in bucket order, the serial phase can apply the dynamic
    // tests to that same list and take the same first survivor the fully serial
    // loop would have taken. This is not an approximation of the old behaviour;
    // it is the old behaviour with the memcmp moved off the critical path.
    const uint32_t CCAP=8;               // candidates retained per tail per level
    std::vector<uint32_t> cand; std::vector<uint8_t> ccnt;
    unsigned NT=std::thread::hardware_concurrency(); if(!NT) NT=1;
    for(uint32_t L=Lmax-1; L>=sweep_minov && L>=SW; --L){
        // compact to tails still open and still long enough at this L
        size_t w=0;
        for(uint32_t a:tails){
            if(nxt[a]!=NONE) continue;                       // already extended
            if(rlen[a]<L) continue;
            const uint32_t off=(uint32_t)rlen[a]-L;          // grows by 1 per length
            if(off+SW>rlen[a]) continue;
            tails[w++]=a;                                    // still open: keep
        }
        tails.resize(w);
        if(tails.empty()) break;

        // ── parallel: roll each seed, probe, keep candidates passing the
        //    static tests, in bucket order ──────────────────────────────────
        cand.assign(w*CCAP,NONE); ccnt.assign(w,0);
        {
            const unsigned T=(w<4096)?1u:NT;
            std::vector<std::thread> th; th.reserve(T);
            std::vector<size_t> pr(T,0);
            const size_t span=(w+T-1)/T;
            for(unsigned t=0;t<T;++t){
                const size_t lo=(size_t)t*span, hi=std::min(w,lo+span);
                if(lo>=hi) break;
                th.emplace_back([&,t,lo,hi]{
                    for(size_t i=lo;i<hi;++i){
                        const uint32_t a=tails[i];
                        const uint32_t off=(uint32_t)rlen[a]-L;
                        // packed reads are always ACGT, so the seed always exists
                        seed[a]=rseed(a,off,SW); ok[a]=1;
                        ++pr[t];
                        auto it=pref.find(seed[a]); if(it==pref.end()) continue;
                        uint8_t c=0;
                        for(uint32_t b:it->second){
                            if(b==a) continue;
                            if(!admit[b]) continue;          // excluded reads are leftovers,
                                                             // never chain members -- without
                                                             // this they get emitted twice
                            if(rlen[b]<L) continue;
                            if(!rcmp(a,off,b,L)) continue;
                            cand[i*CCAP+c]=b;
                            if(++c==CCAP) break;
                        }
                        ccnt[i]=c;
                    }
                });
            }
            for(auto& x:th) x.join();
            for(unsigned t=0;t<T;++t) probes+=pr[t];
        }

        // ── serial: same order, same first survivor ───────────────────────
        // A retained list that filled to CCAP may have been truncated, so if
        // every entry in a full list turns out to be taken, that tail is
        // re-probed serially. Without this the cap silently drops the 9th-and-
        // later candidate and the run loses links -- measured at CCAP=8, 14 of
        // 673,334 links and 1,288 bytes of literal. Rare enough to cost nothing,
        // and it makes the result identical to the fully serial sweep rather
        // than merely close to it.
        for(size_t i=0;i<w;++i){
            const uint32_t a=tails[i];
            bool done=false;
            for(uint8_t c=0;c<ccnt[i];++c){
                const uint32_t b=cand[i*CCAP+c];
                if(prv[b]!=NONE||ch_h[a]==b) continue;       // taken, or would cycle
                nxt[a]=b; prv[b]=a; ovl[a]=L;
                uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h; ++links;
                done=true; break;
            }
            if(done||ccnt[i]<CCAP||!ok[a]) continue;         // list was complete
            const uint32_t off=(uint32_t)rlen[a]-L;
            auto it=pref.find(seed[a]); if(it==pref.end()) continue;
            uint8_t seen=0;
            for(uint32_t b:it->second){
                if(b==a) continue;
                if(!admit[b]) continue;
                if(rlen[b]<L) continue;
                if(!rcmp(a,off,b,L)) continue;
                if(seen++<CCAP) continue;                    // already tried above
                if(prv[b]!=NONE||ch_h[a]==b) continue;
                nxt[a]=b; prv[b]=a; ovl[a]=L;
                uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h; ++links;
                break;
            }
        }
    }
    };

    // ── round 1: label, do not build ────────────────────────────────────────
    // The first sweep exists only to find which reads sit INTERIOR to a chain --
    // overlapped on both sides. Those are the reads that tile the genome well.
    // Reads that end up at a chain end, or unlinked, overlap poorly; keeping
    // them in the pseudogenome is what drags the mean overlap down (measured
    // 119.9 against PgRC2's 133.7, where ideal tiling at this coverage is ~136).
    // They cost far less mapped into the finished sequence than assembled into
    // it, so round 2 rebuilds without them.
    // SW was doing two unrelated jobs: the width of the rolling seed, and the
    // minimum overlap round 1 will classify on. They do not have to be equal --
    // the seed is only a prefix filter, and the L-base memcmp behind it is what
    // decides. Measured, raising SW improves the main pg monotonically and
    // worsens the survivor pg monotonically (32 -> main 11,384,786 / survivor
    // 1,567,625; 16 -> 11,632,777 / 1,294,148), and the best half of each beats
    // PgRC2's total. Splitting the parameter is what makes both reachable at
    // once. Same shape as stage 17's split of acceptance from seed geometry.
    const uint32_t R1MINOV=(argc>7)?(uint32_t)atoi(argv[7]):SW;
    sweep_minov=(R1MINOV>=SW)?R1MINOV:SW;   // classifier floor, seed stays SW
    sweep();
    size_t both_sides=0;
    {
        std::vector<uint8_t> keep(n,0);
        for(uint32_t i=0;i<n;++i)
            if(nxt[i]!=NONE && prv[i]!=NONE){ keep[i]=1; ++both_sides; }
        admit.swap(keep);
    }
    fprintf(stderr,"round1: links=%zu both-sides-overlapped=%zu (%.1f%%)\n",
            links,both_sides,100.0*(double)both_sides/(double)n);
    lap("round 1 (division)");

    // ── round 2: build over the well-overlapped reads only ──────────────────
    sweep_minov=MINOV;                // builder: caller's floor
    sweep();
    fprintf(stderr,"round2: probes=%zu links=%zu\n",probes,links);
    lap("round 2 (assembly)");

    // ── emit chains; singletons held back for mapping ────────────────────────
    std::string pg; pg.reserve((size_t)n*40);
    ppos.assign(n,UINT64_MAX);
    std::vector<uint32_t> leftovers;
    uint32_t multi=0;
    for(uint32_t i=0;i<n;++i){
        if(!admit[i]){ leftovers.push_back(i); continue; }   // excluded in round 1
        if(prv[i]!=NONE) continue;
        if(nxt[i]==NONE){ leftovers.push_back(i); continue; }
        ++multi;
        uint32_t cur=i; ppos[cur]=pg.size(); rappend(pg,cur,0);
        while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                               ppos[cur]=pg.size()-o; rappend(pg,cur,o); }
    }
    fprintf(stderr,"links=%zu chains(multi)=%u leftovers=%zu pg after chains=%zu\n",
            links,multi,leftovers.size(),pg.size());
    lap("emit chains");

    // The prefix index is not consulted anywhere in the mapping stage -- that
    // stage builds and queries its own seed index -- but it was staying live
    // across it, so the two large indexes coexisted at exactly the moment peak
    // RSS is set. It is next needed by the survivor sweep, which runs over ~13k
    // reads rather than 851k, so it is cheaper to drop it here and rebuild a
    // small one there than to carry the full one through.
    { std::unordered_map<uint64_t,std::vector<uint32_t>>().swap(pref); }

    // ── Stage C: pigeonhole mapping into the pg, forward then RC ─────────────
    std::vector<uint8_t> matched(n,0);
    size_t n_matched=0;
    {
        // Part width must SCALE with the tolerance, not stay pinned at 32.
        // Pigeonhole needs MAXMM+1 parts to fit inside one read, so a fixed
        // 32-base part caps tolerance at 150/32-1 = 3 mismatches -- ask for more
        // and the parts no longer fit, no seed is indexed, and mapping silently
        // does nothing. PgRC2 makes the part length the parameter and derives
        // the tolerance from it (targetMismatches = readLength/partLength - 1).
        // Same here: parts are as wide as they can be, capped at 32 because that
        // is what a uint64 holds exactly at 2 bits per base.
        // STAGE 17. The above ties acceptance to the seed geometry: NPARTS =
        // MAXMM+1 parts means a read is rejected past 3 mismatches, because that
        // is all the pigeonhole lemma can certify with parts this wide.
        //
        // PgRC2 does not tie them together. ReadsMatchers.cpp:700 reads
        //     uint8_t maxMismatches = readLength / minCharsPerMismatch;
        // with minCharsPerMismatch = 3 at CODER_LEVEL_NORMAL, so a 150-base read
        // is accepted at up to FIFTY mismatches. The seed (their
        // readsExactMatchingChars = 38) only proposes candidates; acceptance is a
        // separate, far looser test. Their log shows what that buys:
        //     290,521 LQ reads -> 11,298 survivors   (96.1% placed)
        // against ours at 74.8%. The survivor pseudogenome that holds 94% of the
        // remaining gap is made almost entirely of reads this limit rejected.
        //
        // Their seed is 38 bases; ours stays 32 because that is what a uint64
        // holds exactly at 2 bits per base. A shorter seed is more sensitive, not
        // less, so this is not a handicap -- it proposes more candidates, and
        // acceptance is what decides.
        //
        // Not copied: they re-scan and keep the BEST match per read (the
        // "better-matches" column in their log, ReadsMatchers.cpp:315-328). That
        // only shrinks their mismatch stream, which this progression does not
        // measure, so first-acceptable wins here. Noted rather than silently
        // skipped.
        // STAGE 20. Seeds were indexed at fixed multiples of SEEDW (0,32,64,96),
        // which has two defects at a 150-base read: a read whose only clean
        // 32-mer starts at, say, offset 17 is invisible, and bases 128-149 are
        // in no seed at all. PgRC2 does not do this -- at CODER_LEVEL_NORMAL the
        // matcher is mode 'c', CopMEMReadsApproxMatcher, which finds maximal
        // exact matches at ANY offset (mapReadsIntoPg, ReadsMatchers.cpp:715-735).
        // Their fixed 38 is a minimum MEM length, not a position.
        //
        // Note their pigeonhole guarantee is only
        //     targetMismatches = readLength/exactMatchingChars - 1 = 150/38-1 = 2
        // while maxMismatches is 50, so past 2 mismatches their seed is a
        // heuristic filter too. The difference is purely that their filter looks
        // everywhere and ours looked in four places.
        //
        // Sliding by SEEDSTRIDE recovers most of that for a proportional index
        // cost and no extra scan cost -- the pg is still swept once, one hash
        // probe per position, only the candidate lists get longer.
        const uint32_t MAXMAP=(uint32_t)(Lmax/3);        // ReadsMatchers.cpp:700
        // Sensitivity is SEEDW + SEEDSTRIDE - 1: the shortest exact read/pg
        // stretch guaranteed to be found. copMEM's K is 28 with k1*k2 = 10, so
        // theirs is 28+10-1 = 37 against our 32+8-1 = 39. A SHORTER seed with a
        // WIDER stride therefore improves sensitivity and shrinks the index at
        // the same time -- unlike shortening the stride, which improves
        // sensitivity and inflates the index (measured: stride 4 mapped 262 more
        // reads for 150 MB).
        uint32_t SEEDW=(argc>6)?(uint32_t)atoi(argv[6]):32;   // <= 32 (uint64 at 2 bits/base)
        if(SEEDW>32) SEEDW=32; if(SEEDW<8) SEEDW=8;
        const uint32_t SEEDSTRIDE=(argc>5)?(uint32_t)atoi(argv[5]):8;
        const uint32_t NPARTS=(Lmax>=SEEDW)?((Lmax-SEEDW)/SEEDSTRIDE+1):1;
        if(SEEDW<8) SEEDW=8;                     // below this a seed is noise
        const uint64_t SW_MASK=(SEEDW>=32)?~0ULL:((1ULL<<(2*SEEDW))-1);
        fprintf(stderr,"  sensitivity floor: %u bases (SEEDW %u + stride %u - 1)\n",SEEDW+SEEDSTRIDE-1,SEEDW,SEEDSTRIDE);
        auto packW=[&](const char* p,uint64_t& out)->bool{
            uint64_t k=0;
            for(uint32_t i=0;i<SEEDW;++i){ int v=b2(p[i]); if(v<0) return false;
                                          k=(k<<2)|(uint64_t)v; }
            out=k; return true; };
        fprintf(stderr,"  mapping: %u seeds x %u bases, stride %u, accepting up to %u mismatches\n",
                NPARTS,SEEDW,SEEDSTRIDE,MAXMAP);
        // Flat CSR, not unordered_map<uint64,vector<...>>. At stride 8 this index
        // holds ~4M entries over ~3.5M distinct keys, and a map would put a node
        // header plus a separately-allocated vector behind every one of them --
        // measured, that alone took peak RSS from 389 MB to 673 MB. Sorting
        // (key,rid,part) once and pointing an open-addressing table at each key's
        // first occurrence answers the same queries from three flat arrays.
        std::vector<uint64_t> mkey; std::vector<uint32_t> mrid; std::vector<uint8_t> mpart;
        {
            struct E { uint64_t k; uint32_t rid; uint8_t part; };
            std::vector<E> tmp; tmp.reserve(leftovers.size()*NPARTS/2+16);
            for(uint32_t rid:leftovers){
                if(rlen[rid]<SEEDW) continue;
                for(uint32_t p=0;p<NPARTS;++p){
                    const uint32_t off=p*SEEDSTRIDE;
                    if(off+SEEDW>rlen[rid]) break;
                    tmp.push_back({rseed(rid,off,SEEDW),rid,(uint8_t)p});
                }
            }
            std::sort(tmp.begin(),tmp.end(),
                      [](const E&a,const E&b){ return a.k<b.k; });
            mkey.resize(tmp.size()); mrid.resize(tmp.size()); mpart.resize(tmp.size());
            for(size_t i=0;i<tmp.size();++i){ mkey[i]=tmp[i].k; mrid[i]=tmp[i].rid; mpart[i]=tmp[i].part; }
        }
        size_t msize=1; while(msize < mkey.size()*2+1) msize<<=1;
        const uint64_t MMASK=msize-1;
        std::vector<uint32_t> mtab(msize,UINT32_MAX);
        auto mmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
        for(size_t i=0;i<mkey.size();++i){
            if(i && mkey[i]==mkey[i-1]) continue;
            size_t h=mmix(mkey[i])&MMASK;
            while(mtab[h]!=UINT32_MAX) h=(h+1)&MMASK;
            mtab[h]=(uint32_t)i;
        }
        auto mfind=[&](uint64_t k)->uint32_t{
            size_t h=mmix(k)&MMASK;
            while(mtab[h]!=UINT32_MAX){ if(mkey[mtab[h]]==k) return mtab[h]; h=(h+1)&MMASK; }
            return UINT32_MAX;
        };
        fprintf(stderr,"  seed index: %zu entries (CSR)\n",mkey.size());
        // The pg sweep is read-only against `seeds` and `reads`; the only shared
        // write is the matched flag. Threads take disjoint slices and collect
        // hits locally, so the merge is a union -- no atomics, no false sharing,
        // and the result cannot depend on thread count. Reading `matched` inside
        // the loop stays an optimisation only: a stale read costs a redundant
        // verify, never a wrong answer.
        auto scan=[&](const std::string& text){
            if(text.size()<SEEDW) return;
            unsigned T=std::thread::hardware_concurrency(); if(!T) T=1;
            if(text.size()<(1u<<20)) T=1;
            std::vector<std::vector<std::pair<uint32_t,uint32_t>>> hit(T);
            std::vector<std::thread> th; th.reserve(T);
            const size_t chunk=(text.size()+T-1)/T;
            for(unsigned t=0;t<T;++t){
                const size_t lo=(size_t)t*chunk;
                size_t hi=std::min(text.size(),lo+chunk);
                if(lo>=hi) break;
                th.emplace_back([&,t,lo,hi]{
                    // start SEEDW-1 early so seeds spanning the boundary are not lost
                    const size_t from=(lo>=(size_t)(SEEDW-1))?lo-(SEEDW-1):0;
                    uint64_t k=0; uint32_t good=0;
                    for(size_t p=from;p<hi;++p){
                        int v=b2(text[p]);
                        if(v<0){ good=0; k=0; continue; }
                        k=((k<<2)|(uint64_t)v)&SW_MASK; ++good;
                        if(good<SEEDW) continue;
                        const size_t seedStart=p-SEEDW+1;
                        if(seedStart<lo) continue;      // owned by the previous slice
                        const uint32_t ix=mfind(k); if(ix==UINT32_MAX) continue;
                        for(uint32_t q=ix;q<mkey.size()&&mkey[q]==k;++q){
                            const uint32_t rid=mrid[q]; if(matched[rid]) continue;
                            const size_t off=(size_t)mpart[q]*SEEDSTRIDE;
                            if(seedStart<off) continue;
                            const size_t st=seedStart-off;
                            const uint32_t RL=rlen[rid];
                            if(st+RL>text.size()) continue;
                            char rb[256]; runp(rid,rb);
                            uint32_t mm=0; const char* t2=text.data()+st;
                            for(uint32_t j=0;j<RL;++j) if(t2[j]!=rb[j]){ if(++mm>MAXMAP) break; }
                            if(mm<=MAXMAP) hit[t].push_back({rid,(uint32_t)st});
                        }
                    }
                });
            }
            for(auto& x:th) x.join();
            for(auto& v:hit) for(auto& pr:v)
                if(!matched[pr.first]){ matched[pr.first]=1; ppos[pr.first]=pr.second; ++n_matched; }
        };
        scan(pg);
        rc_inplace(pg); scan(pg); rc_inplace(pg);
    }
    lap("pigeonhole mapping");

    // Survivors get their own pseudogenome rather than being appended whole.
    // They are the reads that neither tiled well enough to assemble in round 2
    // nor matched anything already built -- but they still overlap EACH OTHER,
    // so a second sweep over just this set recovers most of their length.
    // Appending them raw costs a full read apiece (measured: 43,075 x 150 =
    // 6.46 MB); PgRC2 assembles its 11,298 survivors into 1,504,035 bytes.
    size_t appended=0, second_pg=0;
    const size_t main_pg_end = pg.size();
    {
        std::vector<uint8_t> keep(n,0);
        for(uint32_t rid:leftovers) if(!matched[rid]){ keep[rid]=1; ++appended; }
        admit.swap(keep);
        // Rebuild the prefix index over the survivors only (see the free above).
        pref.clear(); pref.reserve(appended*2+16);
        for(uint32_t i=0;i<n;++i){
            if(!admit[i]||rlen[i]<SW) continue;
            pref[rseed(i,0,SW)].push_back(i);
        }
        sweep();
        const size_t before=pg.size();
        for(uint32_t i=0;i<n;++i){
            if(!admit[i]||prv[i]!=NONE) continue;            // not a chain head here
            uint32_t cur=i; ppos[cur]=pg.size(); rappend(pg,cur,0);
            while(nxt[cur]!=NONE){ uint32_t o=ovl[cur]; cur=nxt[cur];
                                   ppos[cur]=pg.size()-o; rappend(pg,cur,o); }
        }
        second_pg=pg.size()-before;
        fprintf(stderr,"second pg: %zu reads -> %zu B (raw would be %zu B)\n",
                appended,second_pg,appended*(size_t)Lmax);
    }
    fprintf(stderr,"leftovers=%zu mapped=%zu appended=%zu\n",
            leftovers.size(),n_matched,appended);

    // ── STAGE 22: the order permutation (their stage 6) ─────────────────────
    // To restore the original file order the archive must carry, for every
    // original read, where it sits in emission order. PgRC2 stores the inverse
    // permutation as a raw uint32 array and LZMAs it
    // (SeparatedPseudoGenomePersistence.cpp:226-233) -- no permutation-specific
    // coding at all in single-end mode.
    //
    // Measured, that costs them 7,063,459 - 4,898,620 = 2,164,839 B on this
    // input. The floor for a RANDOM permutation of 999,340 elements is
    // log2(n!) = 2,309,000 B -- so their output is 144 KB BELOW the random
    // floor, which means the permutation is not random and LZMA is finding real
    // structure in it. Any coder that assumes uniformity (Lehmer / factorial
    // base) would therefore do WORSE than what they already have, not better.
    // That is why this stage dumps the permutation for analysis instead of
    // going straight to a clever coder.
    //
    // Note positions from the RC pass are in reverse-complement coordinates.
    // That is fine here: only the ORDER matters for this measurement, and a
    // real encoder would store the strand flag alongside (ARCS already does).
    if(getenv("DUMP_PERM")){
        std::vector<uint32_t> uidOrder; uidOrder.reserve(n);
        for(uint32_t i=0;i<n;++i) if(ppos[i]!=UINT64_MAX) uidOrder.push_back(i);
        std::sort(uidOrder.begin(),uidOrder.end(),
                  [&](uint32_t a,uint32_t b){ return ppos[a]!=ppos[b]?ppos[a]<ppos[b]:a<b; });
        std::vector<uint32_t> rank(n,UINT32_MAX);
        for(uint32_t i=0;i<uidOrder.size();++i) rank[uidOrder[i]]=i;
        // emission order = unique reads in pg order, each followed by every
        // original read that collapsed onto it, in ascending original index
        std::vector<uint32_t> cnt(uidOrder.size()+1,0);
        size_t placed=0;
        for(uint32_t o=0;o<orig2uid.size();++o){
            const uint32_t u=orig2uid[o];
            if(u<n&&rank[u]!=UINT32_MAX){ ++cnt[rank[u]]; ++placed; }
        }
        std::vector<uint32_t> start(uidOrder.size()+1,0);
        for(size_t i=0;i<uidOrder.size();++i) start[i+1]=start[i]+cnt[i];
        std::vector<uint32_t> fill=start, rev(orig2uid.size(),UINT32_MAX);
        for(uint32_t o=0;o<orig2uid.size();++o){
            const uint32_t u=orig2uid[o];
            if(u<n&&rank[u]!=UINT32_MAX) rev[o]=fill[rank[u]]++;
        }
        FILE* f=fopen("perm.u32","wb");
        for(uint32_t o=0;o<rev.size();++o) if(rev[o]!=UINT32_MAX) fwrite(&rev[o],4,1,f);
        fclose(f);
        fprintf(stderr,"[ORDER] originals=%zu placed=%zu unique-in-pg=%zu -> perm.u32\n",
                orig2uid.size(),placed,uidOrder.size());
    }

    // Nothing below this point reads the reads, the prefix index, or any of the
    // sweep's per-read arrays -- the MEM stage works purely on `pg`. Holding
    // them through it was costing more than the stage itself uses: 851,275
    // std::strings at ~150 bases plus a 32-byte header apiece is ~168 MB, and
    // the prefix index another ~80 MB. Released before the index is built, so
    // the two allocations never coexist.
    { std::vector<uint64_t>().swap(rpk); std::vector<uint64_t>().swap(woff);
      std::vector<uint8_t>().swap(rlen); }
    { std::unordered_map<uint64_t,std::vector<uint32_t>>().swap(pref); }
    { std::vector<uint32_t>().swap(nxt);  std::vector<uint32_t>().swap(prv);
      std::vector<uint32_t>().swap(ovl);  std::vector<uint32_t>().swap(ch_h);
      std::vector<uint32_t>().swap(ch_t); std::vector<uint32_t>().swap(tails);
      std::vector<uint64_t>().swap(seed); std::vector<uint8_t>().swap(ok);
      std::vector<uint8_t>().swap(admit); std::vector<uint8_t>().swap(matched); }

    // ── MEM matching between and within pseudogenomes ────────────────────────
    // PgRC2's stage 7 (SimplePgMatcher::matchPgsInPg) and the single largest
    // thing this progression was missing. They build ONE matcher over the HQ pg
    // and then, reverse-complement matching enabled throughout:
    //   - match the LQ pg into it        (SimplePgMatcher.cpp:187)
    //   - match the N pg into it         (:194)
    //   - match the HQ pg into ITSELF    (:202)
    // replacing each maximal exact match with a mark plus (offset,length).
    //
    // Measured from their own stdout on yeast_sub.fq:
    //     HQ 21,104,500 -> 11,463,320   (45.8% removed)
    //     LQ  1,504,035 ->  1,230,006   (18.5%)
    //     N       2,100 ->      1,849   (12.1%)
    // Their real final literal is 12,695,175. The 22,610,635 this progression
    // has been calling "their 22.60M" is their pg BEFORE this stage runs.
    //
    // Their parameters at CODER_LEVEL_NORMAL (pgrc-params.h:138-145): minimum
    // match length 45, reverse complement on.
    //
    // NOT copied: they use CopMEM with sparse sampling, which their own paper
    // concedes can miss matches. This reuses the packed-32-base seed already in
    // this file and samples the SOURCE every STEP = 45-32+1 = 14. That is exact,
    // not heuristic: a match of length >= 45 spans >= 14 consecutive candidate
    // seed starts, which must contain a multiple of 14, so no qualifying match
    // can be missed.
    size_t rem_main=0, rem_second=0, nm_main=0, nm_second=0;
    {
        const size_t MINMEM = 45;
        const size_t STEP   = MINMEM - SEED + 1;   // 14
        const size_t MAXCAND = 64;                 // bound work in repetitive regions

        // Flat CSR index instead of unordered_map<uint64_t,vector<uint32_t>>.
        // That map allocated a separate heap block per distinct key -- ~1.5M of
        // them here -- and paid a node header on each. Sorting (key,pos) pairs
        // once and pointing an open-addressing table at the first occurrence of
        // each key gives the same answers from two flat arrays: all occurrences
        // of a key are adjacent after the sort, so a lookup is one probe plus a
        // forward scan. This is the same exact-size-CSR fix already carried into
        // ARCS from this progression.
        std::vector<uint64_t> skey; std::vector<uint32_t> spos;
        {
            std::vector<std::pair<uint64_t,uint32_t>> tmp;
            tmp.reserve(main_pg_end/STEP+16);
            for(size_t p=0;p+SEED<=main_pg_end;p+=STEP){
                uint64_t k; if(pack(pg.data()+p,k)) tmp.push_back({k,(uint32_t)p});
            }
            std::sort(tmp.begin(),tmp.end());
            skey.resize(tmp.size()); spos.resize(tmp.size());
            for(size_t i=0;i<tmp.size();++i){ skey[i]=tmp[i].first; spos[i]=tmp[i].second; }
        }
        size_t tsize=1; while(tsize < skey.size()*2+1) tsize<<=1;
        const uint64_t TMASK=tsize-1;
        std::vector<uint32_t> htab(tsize,UINT32_MAX);
        auto hmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
        for(size_t i=0;i<skey.size();++i){
            if(i && skey[i]==skey[i-1]) continue;       // first occurrence only
            size_t h=hmix(skey[i])&TMASK;
            while(htab[h]!=UINT32_MAX) h=(h+1)&TMASK;
            htab[h]=(uint32_t)i;
        }
        auto lookup=[&](uint64_t k)->uint32_t{
            size_t h=hmix(k)&TMASK;
            while(htab[h]!=UINT32_MAX){ if(skey[htab[h]]==k) return htab[h]; h=(h+1)&TMASK; }
            return UINT32_MAX;
        };

        enum Mode { CROSS, SELF_FWD, SELF_RC };
        // Greedy left-to-right parse over ONE slice of the destination. Matches
        // are collected rather than written straight into the shared consumed
        // bitmap, so threads never touch the same bytes and the merge below is
        // an ordinary union.
        auto parse_range=[&](const char* Q,size_t qlen,Mode mode,size_t lo,size_t hi,
                             std::vector<std::pair<uint32_t,uint32_t>>& out){
            size_t qp=lo;
            while(qp<hi && qp+MINMEM<=qlen){
                uint64_t k;
                if(!pack(Q+qp,k)){ ++qp; continue; }
                const uint32_t idx=lookup(k);
                if(idx==UINT32_MAX){ ++qp; continue; }
                size_t best=0, tried=0;
                for(uint32_t i=idx;i<skey.size()&&skey[i]==k;++i){
                    if(++tried>MAXCAND) break;
                    const size_t s=spos[i];
                    // Cap usable length BEFORE extending, never after -- in a
                    // self-match the seed at qp also occurs at qp itself and an
                    // uncapped extension runs to end of text.
                    size_t capL=SIZE_MAX;
                    if(mode==SELF_FWD){ if(s>=qp) continue; capL=qp-s; }
                    else if(mode==SELF_RC){ if(s>=qlen-qp) continue; capL=(qlen-qp-s)/2; }
                    if(capL<MINMEM) continue;
                    size_t L=0;
                    while(L<capL && qp+L<qlen && s+L<main_pg_end && Q[qp+L]==pg[s+L]) ++L;
                    if(L>best) best=L;
                }
                if(best>=MINMEM){ out.push_back({(uint32_t)qp,(uint32_t)best}); qp+=best; }
                else ++qp;
            }
        };
        // Split the destination across threads. The greedy parse is serial by
        // nature (it skips ahead by each accepted match), so what is parallelised
        // is the search, exactly as ARCS's repeat_elim does: every query reads
        // only the shared, already-built index. Slicing can only differ from a
        // single serial parse at the T-1 slice boundaries, bounded by ~45 bytes
        // each -- under 0.005% of a 21 MB pseudogenome.
        auto run=[&](const char* Q,size_t qlen,Mode mode,std::vector<uint8_t>& consumed)->size_t{
            unsigned T=std::thread::hardware_concurrency(); if(!T) T=1;
            if(qlen < (1u<<20)) T=1;
            std::vector<std::vector<std::pair<uint32_t,uint32_t>>> res(T);
            std::vector<std::thread> th;
            const size_t chunk=(qlen+T-1)/T;
            for(unsigned t=0;t<T;++t){
                const size_t lo=(size_t)t*chunk, hi=std::min(qlen,lo+chunk);
                if(lo>=hi) break;
                th.emplace_back([&,t,lo,hi]{ parse_range(Q,qlen,mode,lo,hi,res[t]); });
            }
            for(auto& x:th) x.join();
            size_t nm=0;
            for(auto& v:res){
                nm+=v.size();
                for(auto& m:v)
                    for(uint32_t j=0;j<m.second && (size_t)m.first+j<qlen;++j)
                        consumed[m.first+j]=1;
            }
            return nm;
        };

        // Second pg (their LQ+N) against the main pg, forward then reverse complement.
        {
            const size_t qlen=pg.size()-main_pg_end;
            if(qlen>=MINMEM){
                std::string Q(pg,main_pg_end,qlen);
                std::vector<uint8_t> c(qlen,0), cr(qlen,0);
                nm_second =run(Q.data(),qlen,CROSS,c);
                std::string R=Q; rc_inplace(R);
                nm_second+=run(R.data(),qlen,CROSS,cr);
                for(size_t i=0;i<qlen;++i) if(cr[i]) c[qlen-1-i]=1;
                for(size_t i=0;i<qlen;++i) if(c[i]) ++rem_second;
            }
        }
        // Main pg against itself, forward then reverse complement.
        {
            const size_t qlen=main_pg_end;
            std::vector<uint8_t> c(qlen,0), cr(qlen,0);
            nm_main =run(pg.data(),qlen,SELF_FWD,c);
            std::string R(pg,0,qlen); rc_inplace(R);
            nm_main+=run(R.data(),qlen,SELF_RC,cr);
            for(size_t i=0;i<qlen;++i) if(cr[i]) c[qlen-1-i]=1;
            for(size_t i=0;i<qlen;++i) if(c[i]) ++rem_main;
        }
    }
    lap("pg MEM matching");

    const size_t lit_main   = main_pg_end-rem_main;
    const size_t lit_second = (pg.size()-main_pg_end)-rem_second;
    fprintf(stderr,"MEM main   : %zu -> %zu (removed %zu, %.1f%%, %zu matches)\n",
            main_pg_end,lit_main,rem_main,rem_main*100.0/(main_pg_end?main_pg_end:1),nm_main);
    fprintf(stderr,"MEM second : %zu -> %zu (removed %zu, %.1f%%, %zu matches)\n",
            pg.size()-main_pg_end,lit_second,rem_second,
            rem_second*100.0/((pg.size()-main_pg_end)?(pg.size()-main_pg_end):1),nm_second);
    fprintf(stderr,"PgRC2 for reference: main 21104500 -> 11463320 (45.8%%), "
                   "lq+n 1506135 -> 1231855, final literal 12695175\n");
    printf("PG_LEN %zu\n",pg.size());
    printf("PG_LITERAL %zu\n",lit_main+lit_second);
    return 0;
}
