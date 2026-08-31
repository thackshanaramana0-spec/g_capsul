// ============================================================================
// STAGE 90 -- CRITICAL FIX: silent data loss on reads >255 bases.
//
// Found while properly SCOPING generalization (per direct user instruction:
// "do not [win on] one dataset... first define scope... how many types of
// ds and variation") -- testing M. tuberculosis (ERR552797, one of the 10
// locked accessions, never touched before this) instead of re-testing
// files already measured. ERR552797 uses 300bp reads.
//
// The bug: `rlen` was `std::vector<uint8_t>` (stage 87 line 152), and the
// loader silently `continue`d past any read >255 bases (stage 87 line 177)
// -- NOT counted in n_filt like the N-filter above it, just dropped with no
// record at all. On a 100,000-read real slice of ERR552797, 85,622 reads
// (299-301bp) were silently discarded before dedup even ran; the loader
// then reported "unique=8189" out of the 8,189 reads that happened to be
// <=255bp, which looked like a plausible (if high) dedup rate until checked
// directly: `sort -u` on the raw sequences shows 99,279 of 100,000 are
// genuinely distinct. This was never a dedup problem -- it was near-total
// silent data loss, invisible on every file actually tested this session
// (yeast 150bp, E. coli 150bp, P. falciparum 100bp, SARS-CoV-2 ~100-221bp)
// because none of them happened to exceed 255bp. Exactly the failure mode
// proper scoping exists to catch: an architectural assumption baked in
// early (uint8 length field), never exercised against the real locked
// scope's actual read-length range until this file was tried.
//
// Fix: widen the length field to uint16_t (max 65535, far beyond any real
// short-read length -- MiSeq 2x300 tops out at 300bp) and the two matching
// fixed-size unpack buffers (`ubuf`/`b_`, stage 87 lines 161/231) from 256
// to 1024 bytes so they can't overflow once longer reads are actually kept.
// The size-cutoff filter moves from 255 to 1023 to match. Every other use
// of uint8_t in this file (readMM's per-read mismatch count, 255 as an
// "unplaced" sentinel) is unrelated to read length and is left untouched --
// mismatch counts realistically never approach that range.
// ============================================================================
//
// ============================================================================
// STAGE 87 -- CRITICAL FIX: a real, reproducible race condition in the
// combined pipeline's task-level overlap (stages 76/82/86), found while
// investigating something else entirely (a RAM measurement wrapper
// perturbed timing enough to expose it) and confirmed deliberately before
// accepting it, per direct user instruction to think and plan before
// changing anything.
//
// The bug: `fopen("literal.txt","wb")` (and the same pattern for
// perm.u32) makes the file VISIBLE to a concurrent reader immediately,
// but its content is not COMPLETE until fclose(). The pipeline's overlap
// mechanism polls file EXISTENCE (`[ -f literal.txt ]`), not completion,
// to decide when to launch seqpar/permcoder early. For literal.txt
// specifically, the file stays open across the ENTIRE MEM-matching
// compute phase (both self-match passes), not just a write loop --
// measured directly with steady_clock timestamps on real data: the gap
// between fopen and fclose is ~0.66 SECONDS, not a microsecond fluke.
//
// Reproduced the actual failure, deterministically, using the real
// pipeline's exact polling logic (10ms sleep, existence check) against
// BOTH this file's instrumented build AND the plain, unmodified
// production `best` binary already used for every earlier measurement
// this session: seqpar read literal.txt as completely empty (bases=0,
// coded=0 B), 3 times out of 3 test runs. Every earlier "beat SPRING"
// speed/RAM measurement in this project won this same race by timing
// luck specific to running inside the FULL combined pipeline (where
// names/quality also compete for CPU) -- it was never actually a proven
// safe mechanism, just one that happened not to fail under that specific
// load pattern.
//
// The fix is the standard, established one for exactly this class of bug:
// write to a temporary filename, then rename() onto the final name only
// after fclose() succeeds. POSIX rename() is atomic on the same
// filesystem, so a concurrent reader can only ever observe "not there
// yet" or "fully there" -- structurally impossible to race, not merely
// less likely to race. Applied to both literal.txt and perm.u32 (the two
// files the polling loop watches; mem_triples.bin is only read by
// refcoder AFTER the whole assembly process exits, via `wait $PID_ASM`,
// so it was never exposed to this specific race).
//
// Verified: byte-identical output vs the unfixed binary (md5sum match on
// literal.txt/perm.u32/mem_triples.bin/mem_gaps.bin/mem_lens.bin). Real
// stress test, exact same methodology that reproduced the bug: 5 runs
// with the fix, 5/5 correct (bases=12,506,313, coded=3,007,051 B every
// time) -- versus 3/3 failures on the unfixed binary using the identical
// test.
// ============================================================================
//
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
#include <omp.h>
#include <unistd.h>
#if defined(__GLIBC__)
#include <malloc.h>
#endif

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
    std::vector<uint16_t> rlen;           // length in bases (uint16: 300bp reads are real, uint8 silently dropped them)
    // STAGE 22 (their stage 6, OrderInfo): to restore the original file order we
    // need, for every ORIGINAL read, where its sequence sits in the pg. Dedup
    // collapses duplicates, so the original->unique map has to be kept too.
    std::vector<uint32_t> orig2uid;       // original read (N-filtered out) -> unique id
    size_t n_in=0,n_filt=0;
    {
        std::ifstream f(argv[1]); std::string a,b,c,d;
        std::unordered_map<uint64_t,std::vector<uint32_t>> seen; seen.reserve(1u<<21);
        std::vector<uint64_t> tmpw; char ubuf[1024];
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
            if(b.size()>1023){ ++n_filt; continue; }         // uint16 length field, matches ubuf/b_ buffer size
            auto& bk=seen[fnv(b.data(),(uint32_t)b.size())]; bool dup=false;
            for(uint32_t id:bk){
                if(rlen[id]!=(uint16_t)b.size()) continue;
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
                rlen.push_back((uint16_t)b.size());
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
        char b_[1024]; runp(i,b_); dst.append(b_+from,rlen[i]-from);
    };

    // ── STAGE 42: prefix-containment removal, generalising dedup ─────────────
    // Real variable-length data comes from 3' quality trimming, so a trimmed
    // read is a strict PREFIX of its untrimmed twin. Exact-hash dedup cannot see
    // that, and measured on 30% of reads trimmed to 130-149 bases the cost is
    // severe: unique reads 851,275 -> 918,598, both-side-overlapped 69.8% ->
    // 40.0%, survivors 12,065 -> 57,230, literal +24%.
    //
    // Myers' string graph (2005) and SGA (2012) both drop contained reads before
    // building, because a contained read carries nothing its container does not.
    // The important property here is that containment SUBSUMES dedup: at fixed
    // length one read contains another only when they are equal, so this
    // degenerates to exactly the existing behaviour and the fixed-length results
    // are unchanged -- which is checked, not assumed.
    //
    // Sorting lexicographically puts a prefix immediately before its extensions,
    // so one linear scan finds every prefix containment. A contained read is then
    // handled exactly like a duplicate: its originals are redirected to the
    // container (offset 0, so no extra bookkeeping) and it leaves the read set.
    std::vector<uint8_t> contained(n,0);
    size_t n_contained=0;
    // At constant read length one read can only contain another by equalling it,
    // which dedup has already handled -- so the whole pass is guaranteed to find
    // nothing and is pure cost. Measured, it was 0.31 s inside the load stage
    // (0.56 -> 0.87 s), which is 8% of the entire run spent proving a negative.
    bool varlen=false;
    for(uint32_t i=1;i<n;++i) if(rlen[i]!=rlen[0]){ varlen=true; break; }
    if(varlen){
        std::vector<uint32_t> ord(n);
        for(uint32_t i=0;i<n;++i) ord[i]=i;
        auto pfxcmp=[&](uint32_t x,uint32_t y)->int{     // compare min(len) bases
            uint32_t L=(rlen[x]<rlen[y])?rlen[x]:rlen[y];
            uint64_t px=woff[x]*32ULL, py=woff[y]*32ULL;
            while(L>=32){ const uint64_t a=w32(px),b=w32(py);
                          if(a!=b) return a<b?-1:1; px+=32; py+=32; L-=32; }
            if(L){ const uint64_t m=~0ULL<<(2*(32-L));
                   const uint64_t a=w32(px)&m,b=w32(py)&m;
                   if(a!=b) return a<b?-1:1; }
            return 0;
        };
        std::sort(ord.begin(),ord.end(),[&](uint32_t x,uint32_t y){
            const int c=pfxcmp(x,y);
            return c ? c<0 : rlen[x]<rlen[y];        // shorter first among prefixes
        });
        // a read is contained when the next entry shares all of its bases
        std::vector<uint32_t> alias(n,UINT32_MAX);
        for(size_t i=0;i+1<ord.size();++i){
            const uint32_t a=ord[i], b=ord[i+1];
            if(rlen[a]<=rlen[b] && pfxcmp(a,b)==0 && rlen[a]<rlen[b]){
                contained[a]=1; alias[a]=b; ++n_contained;
            }
        }
        // a container may itself be contained, so follow the chain to the end
        for(uint32_t i=0;i<n;++i){
            if(!contained[i]) continue;
            uint32_t t=alias[i]; int guard=0;
            while(t!=UINT32_MAX && contained[t] && ++guard<64) t=alias[t];
            alias[i]=t;
        }
        for(uint32_t& u:orig2uid) if(u<n && contained[u] && alias[u]!=UINT32_MAX) u=alias[u];
        if(n_contained)
            fprintf(stderr,"contained reads folded into their container: %zu\n",n_contained);
    }
    // ord and alias are n x 4 bytes each and go out of scope above, but glibc
    // keeps the arena -- measured, the containment pass cost 12 MB of peak
    // (193,736 -> 205,844 KB) purely by being retained. Same lesson as stage 41.
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
    fprintf(stderr,"reads in=%zu N-filtered=%zu unique=%u maxlen=%u\n",n_in,n_filt,n,Lmax);
    lap("load+filter+dedup");

    // `admit` selects which reads take part; hoisted above the index because
    // buildPref filters on it when rebuilding for the survivor pass.
    std::vector<uint8_t> admit(n,1);
    for(uint32_t i=0;i<n;++i) if(contained[i]) admit[i]=0;   // stage 42

    // ── STAGE 43: the prefix index, flat ─────────────────────────────────────
    // This was still an unordered_map<uint64,vector<uint32>> -- the exact
    // structure stages 18 and 40 replaced elsewhere for large wins, left in the
    // one place it is hit hardest: rounds 1 and 2 probe it 20.9M times and cost
    // 2.21 s of a 4.7 s run. 851k distinct keys each carry a node header plus a
    // separately allocated vector, ~95 MB, and every probe is a pointer chase.
    //
    // Flat instead: (key32 << 32) | rid in one sorted array, with an
    // open-addressing table pointing at each key's first entry. Sorting the whole
    // word sorts by key because the key is in the high bits, and rid ascends
    // within a key -- so candidates are visited in the same order the map's
    // insertion-ordered vector gave, and the output is unchanged. ~15 MB.
    //
    // The 32-bit key is exact at the default SW 16 and partial above it, which is
    // safe: a collision only proposes a candidate, and rcmp verifies every
    // candidate against the full L bases anyway.
    std::vector<uint64_t> pent;
    std::vector<uint32_t> ptab;
    size_t pmask=0;
    auto pmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
    // STAGE 97 -- real fix, read directly from SPRING's actual source
    // (bitset_util.cpp: bbhashdict::remove()), not guessed. SPRING's
    // dictionary physically removes a read from its bucket the moment it's
    // matched (binary search for it, shift-compact the rest left), so
    // buckets shrink continuously as reordering progresses. Ours never did
    // this -- `pent` was built once and reused unchanged across every one
    // of round1/round2's ~70-100 L-levels, so a repeat-driven bucket on a
    // genome like T. cacao gets walked at its FULL original size on every
    // level, even for reads that were linked (taken) many levels ago and are
    // already skipped via the `prv[b]!=NONE` check -- the check prevents a
    // WRONG match, but does nothing to shrink the cost of finding it. Three
    // earlier ideas (stage 93 cap, stage 95 stride, stage 96 longer key) were
    // all tested directly and disproved -- this is the first one grounded in
    // what a real, comparable tool's source code actually does differently.
    //
    // `liveLen[s]` (s = a bucket's start index) tracks each bucket's CURRENT
    // length, mutated by `removeRead` as reads get taken -- `pent`/`ptab`
    // themselves stay fixed (correct: only entries WITHIN a bucket's own
    // span move on removal, never across bucket boundaries, since removal is
    // a same-bucket compaction found via that read's own key).
    std::vector<uint32_t> liveLen;
    // STAGE 98 -- stage 97's removal was correct but O(current bucket size)
    // per removal (find-by-scan + shift-compact, mirroring SPRING's own
    // remove() exactly) -- measured to be a net loss on well-behaved genomes
    // (C. elegans 42.56s -> 45.97s, +8%) since that cost is paid on EVERY
    // link, and normal genomes' buckets are already tiny (E. coli avg scan
    // 0.12) -- nothing expensive was being saved, pure overhead added. A
    // REMOVE_THRESH gate (skip removal below bucket size 64) recovered some
    // of that but still cost +4% on C. elegans -- not yet a clean win even
    // on the easy case.
    //
    // Real fix: nothing downstream depends on WITHIN-bucket order (that was
    // only ever an artifact of sorting pent by the full packed key, never a
    // correctness requirement -- SPRING keeps order because ITS remove()
    // uses lower_bound to locate a read, which needs sorted order; ours can
    // locate a read in O(1) instead via a direct reverse index, so it never
    // needs order preservation in the first place). `posOf[read_id]` maps
    // straight to that read's current slot in `pent` -- removal becomes:
    // look up its slot (O(1)), swap it with the bucket's current last live
    // slot (O(1)), fix up the swapped read's posOf entry (O(1)), shrink
    // liveLen (O(1)). No threshold needed: removal is now cheap enough to
    // do unconditionally, every link, restoring full shrinking benefit with
    // none of the O(bucket) cost that made it a net loss before.
    std::vector<uint32_t> posOf;
    auto buildPref=[&](bool useAdmit){
        pent.clear(); pent.reserve(n);
        for(uint32_t i=0;i<n;++i){
            if(useAdmit && !admit[i]) continue;
            if(rlen[i]<SW) continue;
            pent.push_back(((uint64_t)(uint32_t)rseed(i,0,SW)<<32)|(uint64_t)i);
        }
        std::sort(pent.begin(),pent.end());
        size_t sz=1; while(sz<pent.size()*2+1) sz<<=1;
        pmask=sz-1; ptab.assign(sz,UINT32_MAX);
        liveLen.assign(pent.size(),0);
        posOf.assign(n,UINT32_MAX);
        size_t i=0;
        while(i<pent.size()){
            const uint32_t k=(uint32_t)(pent[i]>>32);
            size_t j=i+1;
            while(j<pent.size() && (uint32_t)(pent[j]>>32)==k) ++j;
            liveLen[i]=(uint32_t)(j-i);
            for(size_t q=i;q<j;++q) posOf[(uint32_t)(pent[q]&0xFFFFFFFFULL)]=(uint32_t)q;
            size_t h=pmix(k)&pmask;
            while(ptab[h]!=UINT32_MAX) h=(h+1)&pmask;
            ptab[h]=(uint32_t)i;
            i=j;
        }
    };
    auto pfind=[&](uint64_t key)->uint32_t{
        const uint32_t k=(uint32_t)key;
        size_t h=pmix(k)&pmask;
        while(ptab[h]!=UINT32_MAX){
            if((uint32_t)(pent[ptab[h]]>>32)==k) return ptab[h];
            h=(h+1)&pmask;
        }
        return UINT32_MAX;
    };
    // Serial-only (called from the level's commit phase, never the parallel
    // probe phase). O(1): no scan, no shift.
    auto removeRead=[&](uint32_t b){
        if(rlen[b]<SW) return;                // never indexed -- nothing to remove
        const uint32_t q=posOf[b];
        if(q==UINT32_MAX) return;             // already removed, or never indexed
        const uint32_t bk=(uint32_t)rseed(b,0,SW);
        const uint32_t s=pfind(bk);
        if(s==UINT32_MAX) return;
        const uint32_t last=s+liveLen[s]-1;
        if(q!=last){
            const uint32_t moved=(uint32_t)(pent[last]&0xFFFFFFFFULL);
            pent[q]=pent[last];
            posOf[moved]=q;
        }
        posOf[b]=UINT32_MAX;
        --liveLen[s];
    };
    buildPref(false);
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
    size_t mm_total=0, mm_reads=0; std::vector<size_t> mm_hist(256,0);
    // STAGE 28. Best match, not first acceptable.
    // Stage 27 measured what first-acceptable costs: 10.74 mismatches per placed
    // read against PgRC2's 1.39. Every one of those has to be stored, so at
    // their coded rate (265,900 B for 1,369,413 mismatches = 1.55 bits each)
    // first-acceptable was quietly spending ~260 KB more than they do -- over
    // half the order-coding win, in a stream this progression never counted.
    //
    // PgRC2 re-scans and keeps the best placement (ReadsMatchers.cpp:315-328).
    // The mechanism that makes it cheap is that the current best doubles as the
    // early-exit bound: a candidate is abandoned the moment it is worse than
    // what the read already has, so better placements also mean less work.
    std::vector<uint8_t> readMM;          // 255 = unplaced
    std::vector<uint32_t> nxt(n,NONE),prv(n,NONE),ovl(n,0),ch_h(n),ch_t(n);
    std::vector<uint64_t> seed(n,0);
    std::vector<uint8_t>  ok(n,0);
    std::vector<uint32_t> tails;
    size_t links=0, probes=0;
    // `admit` selects which reads take part; everything else is reset per round.
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
        //
        // STAGE 45. This block runs once per overlap level -- 134 levels in
        // round 1, 110 in round 2 -- so anything paid per level is paid 244
        // times. Three such costs were here, none of them the actual search:
        //
        //   cand.assign(w*CCAP,NONE) memset 27 MB at the widest level and
        //   several MB at most others, to write a sentinel that is never read:
        //   the serial phase reads cand[i*CCAP+c] only for c < ccnt[i]. Only
        //   ccnt has to be cleared, because a tail whose seed misses the index
        //   `continue`s before assigning it. So the buffer is grown once to the
        //   high-water mark and never cleared again.
        //
        //   12 std::threads created and joined per level -- 2,928 spawns for a
        //   run whose whole budget is under four seconds. OpenMP keeps one pool
        //   alive across all of them, and schedule(static) over [0,w) hands each
        //   thread the same contiguous span the manual split did, so every tail
        //   is searched by exactly the work it was before.
        //
        //   ++pr[t] on a vector<size_t> puts 12 counters in two cache lines and
        //   increments them 20.9M times from 12 cores: every increment steals
        //   the line back. A local counter added once at the end is the same
        //   number with none of the traffic.
        //
        // The output is unchanged by construction -- the search is per-tail
        // independent and the serial commit below still walks i in order.
        if(cand.size()<w*CCAP) cand.resize(w*CCAP);
        ccnt.assign(w,0);
        {
            const unsigned T=(w<4096)?1u:NT;
            size_t pr_total=0;
            #pragma omp parallel for schedule(static) num_threads(T) reduction(+:pr_total)
            for(long long ii=0;ii<(long long)w;++ii){
                const size_t i=(size_t)ii;
                const uint32_t a=tails[i];
                const uint32_t off=(uint32_t)rlen[a]-L;
                // packed reads are always ACGT, so the seed always exists
                seed[a]=rseed(a,off,SW); ok[a]=1;
                ++pr_total;
                const uint32_t pix=pfind(seed[a]);
                if(pix==UINT32_MAX) continue;
                uint8_t c=0;
                // STAGE 97: bounded by liveLen, not the original static
                // bucket end -- reads removeRead()'d out by earlier levels'
                // commits are already physically gone from this span, so
                // this scan only ever walks currently-live candidates.
                const uint32_t qend=pix+liveLen[pix];
                for(uint32_t q=pix;q<qend;++q){
                    const uint32_t b=(uint32_t)(pent[q]&0xFFFFFFFFULL);
                    if(b==a) continue;
                    if(!admit[b]) continue;                  // excluded reads are leftovers,
                                                             // never chain members -- without
                                                             // this they get emitted twice
                    if(rlen[b]<L) continue;
                    if(!rcmp(a,off,b,L)) continue;
                    cand[i*CCAP+c]=b;
                    if(++c==CCAP) break;
                }
                ccnt[i]=c;
            }
            probes+=pr_total;
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
                nxt[a]=b; prv[b]=a; ovl[a]=L; removeRead(b);
                uint32_t h=ch_h[a],t=ch_t[b]; ch_t[h]=t; ch_h[t]=h; ++links;
                done=true; break;
            }
            if(done||ccnt[i]<CCAP||!ok[a]) continue;         // list was complete
            const uint32_t off=(uint32_t)rlen[a]-L;
            const uint32_t pix=pfind(seed[a]);
            if(pix==UINT32_MAX) continue;
            const uint32_t pk=(uint32_t)seed[a];
            uint8_t seen=0;
            for(uint32_t q=pix;q<pent.size()&&(uint32_t)(pent[q]>>32)==pk;++q){
                const uint32_t b=(uint32_t)(pent[q]&0xFFFFFFFFULL);
                if(b==a) continue;
                if(!admit[b]) continue;
                if(rlen[b]<L) continue;
                if(!rcmp(a,off,b,L)) continue;
                if(seen++<CCAP) continue;                    // already tried above
                if(prv[b]!=NONE||ch_h[a]==b) continue;
                nxt[a]=b; prv[b]=a; ovl[a]=L; removeRead(b);
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
    // STAGE 97: round1's removeRead() calls shrank `pent`/`liveLen` as it
    // matched reads -- round2 runs its OWN independent sweep (fresh nxt/prv,
    // different admit set) and must not inherit round1's removals, or it
    // would silently lose reads round1 happened to link but round2 might
    // validly want too. Rebuild fresh before round2 starts.
    buildPref(false);
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
        if(contained[i]) continue;                 // stage 42: its container carries it
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
    // STAGE 41: return the freed arena to the OS.
    //
    // Three targeted reductions in a row measured ZERO -- ppos uint64->uint32
    // (-3.4 MB of array), the packed mapping index (13 B -> 8 B/entry and a
    // 30 MB temp removed), and the parallel loader. Each shrank something real
    // and none moved peak RSS, which is the signal: the peak is not the live
    // set. Live structures at mapping total ~155 MB against a measured 232 MB.
    //
    // The gap is glibc holding freed arenas. ARCS hit exactly this and says so
    // in src/vodbg_pg.cpp: "malloc_trim is what makes the release visible: glibc
    // keeps freed arenas by default, so an earlier measurement recovered only
    // 52 MB of a 312 MB free until the arena was returned explicitly." The
    // sweep's per-level candidate array and the 851k-entry prefix map are both
    // freed before mapping and both stay resident without this.
    { std::vector<uint64_t>().swap(pent); std::vector<uint32_t>().swap(ptab); }
#if defined(__GLIBC__)
    malloc_trim(0);
#endif

    // ── Stage C: pigeonhole mapping into the pg, forward then RC ─────────────
    std::vector<uint8_t> matched(n,0);
    readMM.assign(n,255);
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
        const uint32_t K2=(argc>8)?(uint32_t)atoi(argv[8]):1;   // query-side stride
        const uint8_t MMGOOD=(uint8_t)(getenv("MMGOOD")?atoi(getenv("MMGOOD")):0);
        const uint32_t NPARTS=(Lmax>=SEEDW)?((Lmax-SEEDW)/SEEDSTRIDE+1):1;
        if(SEEDW<8) SEEDW=8;                     // below this a seed is noise
        const uint64_t SW_MASK=(SEEDW>=32)?~0ULL:((1ULL<<(2*SEEDW))-1);
        fprintf(stderr,"  sensitivity floor: %u bases (SEEDW %u + k1 %u * k2 %u - 1)\n",SEEDW+SEEDSTRIDE*K2-1,SEEDW,SEEDSTRIDE,K2);
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
        // STAGE 40: one packed array instead of three plus a temp.
        //
        // The modelled orientation flip (index the pg, stream reads) is
        // unaffordable: at 2.5x fewer operations it needs 141 MB of index and
        // peak RSS is already 233 MB against PgRC2's 232 MB. At EQUAL index size
        // the flip is 80% worse. So their mapping advantage is not orientation
        // but a denser index -- copMEM at their parameters samples the pg every
        // 5, giving ~4.3M entries against our 1.85M inside the same footprint.
        // The way to afford density is a smaller entry, which is this.
        //
        // Was: mkey(8) + mrid(4) + mpart(1) = 13 B/entry, built through a
        // vector<E> of 16 B/entry that coexists with all three during the copy --
        // 29 B/entry at the build moment, ~54 MB here.
        //
        // Now: (key32 << 32) | (rid << 3) | part in ONE uint64, sorted in place.
        // Sorting the whole word sorts by key, because the key occupies the high
        // bits -- so the run-detection and first-occurrence logic are unchanged.
        // 8 B/entry, no temp, no copy: ~15 MB.
        //
        // The key is truncated to 32 bits, which is exact for SEEDW <= 16 (the
        // default) and a partial key above it. Partial is still correct: a
        // colliding key only produces a candidate, and every candidate is
        // verified against the text anyway. It costs a few wasted verifies, never
        // a wrong placement.
        std::vector<uint64_t> ment;
        {
            // Exactly the upper bound, so push_back never doubles and no
            // shrink_to_fit copy is needed. Reserving half instead made the
            // vector reallocate mid-fill and then copy again on shrink, which
            // cost more peak than the packing saved: mapping RSS went 232 -> 240
            // MB on the first attempt.
            ment.reserve((size_t)leftovers.size()*NPARTS);
            for(uint32_t rid:leftovers){
                if(rlen[rid]<SEEDW) continue;
                for(uint32_t p=0;p<NPARTS;++p){
                    const uint32_t off=p*SEEDSTRIDE;
                    if(off+SEEDW>rlen[rid]) break;
                    const uint64_t k32=(uint64_t)(uint32_t)rseed(rid,off,SEEDW);
                    ment.push_back((k32<<32)|((uint64_t)rid<<3)|(uint64_t)p);
                }
            }
            std::sort(ment.begin(),ment.end());
        }
        auto MKEY =[&](size_t i)->uint32_t{ return (uint32_t)(ment[i]>>32); };
        auto MRID =[&](size_t i)->uint32_t{ return (uint32_t)((ment[i]>>3)&0x1FFFFFFFULL); };
        auto MPART=[&](size_t i)->uint32_t{ return (uint32_t)(ment[i]&7ULL); };
        size_t msize=1; while(msize < ment.size()*2+1) msize<<=1;
        const uint64_t MMASK=msize-1;
        std::vector<uint32_t> mtab(msize,UINT32_MAX);
        auto mmix=[](uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; return x; };
        for(size_t i=0;i<ment.size();++i){
            if(i && MKEY(i)==MKEY(i-1)) continue;
            size_t h=mmix(MKEY(i))&MMASK;
            while(mtab[h]!=UINT32_MAX) h=(h+1)&MMASK;
            mtab[h]=(uint32_t)i;
        }
        // STAGE 26. The mapping stage is probe-bound, not verify-bound and not
        // probe-COUNT-bound. 43M probes (21.6M pg positions x 2 directions) hit
        // an index of ~24 MB of CSR arrays plus a 16 MB open-addressing table --
        // past L3, so essentially every probe is a memory access. That is also
        // why stage 24's coprime sampling lost: it halved the probes but doubled
        // the table, and the table was the problem.
        //
        // Almost every probe MISSES: most pg positions carry no read seed. A one
        // bit per slot presence filter is 2 MB at this size, stays in cache, and
        // rejects those misses without ever touching the 40 MB structure. Only a
        // survivor pays for the real lookup. False positives cost one wasted
        // lookup and cannot change the result.
        std::vector<uint64_t> filt((size_t)1<<18,0);      // 2^24 bits = 2 MB
        const uint64_t FMASK=((uint64_t)1<<24)-1;
        for(size_t i=0;i<ment.size();++i){
            if(i && MKEY(i)==MKEY(i-1)) continue;
            const uint64_t h=mmix(MKEY(i))&FMASK;
            filt[h>>6]|=1ULL<<(h&63);
        }
        auto mmaybe=[&](uint64_t k)->bool{
            const uint64_t h=mmix(k)&FMASK;
            return (filt[h>>6]>>(h&63))&1ULL;
        };
        auto mfind=[&](uint64_t k)->uint32_t{
            size_t h=mmix(k)&MMASK;
            while(mtab[h]!=UINT32_MAX){ if(MKEY(mtab[h])==(uint32_t)k) return mtab[h]; h=(h+1)&MMASK; }
            return UINT32_MAX;
        };
        fprintf(stderr,"  seed index: %zu entries (packed, %.1f MB)\n",
                ment.size(),ment.size()*8.0/1048576.0);
        // The pg sweep is read-only against `seeds` and `reads`; the only shared
        // write is the matched flag. Threads take disjoint slices and collect
        // hits locally, so the merge is a union -- no atomics, no false sharing,
        // and the result cannot depend on thread count. Reading `matched` inside
        // the loop stays an optimisation only: a stale read costs a redundant
        // verify, never a wrong answer.
        auto scan=[&](const std::string& text){
            if(text.size()<SEEDW) return;
            // STAGE 25. The verify loop was the whole of the 2.6x mapping
            // deficit. Stage 21 packed the reads but then, for every candidate,
            // unpacked all 150 bases into a stack buffer and compared them
            // byte-by-byte -- ~300 operations per candidate, run ~26M times.
            //
            // PgRC2 does not unpack (SymbolsPackingFacility.cpp:344-362): it
            // packs the pattern a word at a time and compares whole packed
            // words, only falling back to per-symbol work when a word differs.
            // Packing the text ONCE per direction goes further -- both sides are
            // then packed, so a 150-base comparison is 5 word loads, 5 XORs and
            // 5 popcounts instead of 300 byte operations.
            std::vector<uint64_t> tpk((text.size()+31)/32+2,0);
            for(size_t i=0;i<text.size();i+=32){
                uint64_t v=0;
                for(uint32_t j=0;j<32;++j){
                    const size_t idx=i+j;
                    v=(v<<2)|(uint64_t)((idx<text.size())?std::max(0,b2(text[idx])):0);
                }
                tpk[i>>5]=v;
            }
            auto tw32=[&](uint64_t ab)->uint64_t{
                const size_t wi=ab>>5; const uint32_t sh=(uint32_t)(ab&31)*2;
                const uint64_t hi=(wi<tpk.size())?tpk[wi]:0ULL;
                if(!sh) return hi;
                const uint64_t lo=(wi+1<tpk.size())?tpk[wi+1]:0ULL;
                return (hi<<sh)|(lo>>(64-sh));
            };
            unsigned T=std::thread::hardware_concurrency(); if(!T) T=1;
            if(text.size()<(1u<<20)) T=1;
            std::vector<std::vector<std::pair<uint32_t,uint32_t>>> hit(T);
            std::vector<std::vector<uint32_t>> hmm(T);
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
                        // STAGE 24: copMEM's lemma applies to the QUERY side too.
                        // Stage 20 took it for sensitivity (sampling the reads)
                        // but still probed every one of the pg's 21.6M positions,
                        // i.e. k2 = 1. Their k1*k2 <= L-K+1 with BOTH sides
                        // sampled is what makes their matching 0.41 s against our
                        // 1.06 s. Sampling the query every K2 and the reads every
                        // SEEDSTRIDE keeps the same floor as long as the product
                        // is unchanged, so K2=2 with stride halved is free
                        // sensitivity-wise and halves the probes.
                        if(K2>1 && (seedStart%K2)) continue;
                        if(!mmaybe(k)) continue;              // 2 MB filter, stays cached
                        const uint32_t ix=mfind(k); if(ix==UINT32_MAX) continue;
                        for(uint32_t q=ix;q<ment.size()&&MKEY(q)==(uint32_t)k;++q){
                            const uint32_t rid=MRID(q);
                            const size_t off=(size_t)MPART(q)*SEEDSTRIDE;
                            if(seedStart<off) continue;
                            const size_t st=seedStart-off;
                            const uint8_t cur=readMM[rid];
                            // Rescanning every already-placed read is what the
                            // best-match pass costs. A read already placed at or
                            // below MMGOOD is left alone: the remaining gain per
                            // read is small and the work per read is not.
                            if(cur<=MMGOOD) continue;
                            const uint32_t lim=(cur==255)?MAXMAP:(uint32_t)cur-1;
                            const uint32_t RL=rlen[rid];
                            if(st+RL>text.size()) continue;
                            // packed vs packed: 2 bits per base, so a XOR leaves a
                            // non-zero 2-bit field exactly where the bases differ.
                            // Folding the two bits together and masking to the low
                            // bit of each pair turns that into a popcount.
                            uint32_t mm=0;
                            const uint64_t rbase=woff[rid]*32ULL;
                            for(uint32_t j=0;j<RL;j+=32){
                                const uint32_t take=(RL-j<32)?(RL-j):32;
                                uint64_t d=w32(rbase+j)^tw32((uint64_t)st+j);
                                if(take<32) d&=~0ULL<<(2*(32-take));
                                d=(d|(d>>1))&0x5555555555555555ULL;
                                mm+=(uint32_t)__builtin_popcountll(d);
                                if(mm>lim) break;
                            }
                            if(mm<=lim){ hit[t].push_back({rid,(uint32_t)st}); hmm[t].push_back(mm); }
                        }
                    }
                });
            }
            for(auto& x:th) x.join();
            // Keep the fewest-mismatch placement per read. Threads may each
            // have found a different one; the merge is serial and authoritative,
            // so the readMM they raced on was only ever a bound hint.
            for(size_t t=0;t<hit.size();++t)
                for(size_t q=0;q<hit[t].size();++q){
                    const auto& pr=hit[t][q]; const uint32_t m=hmm[t][q];
                    if(m<readMM[pr.first]){
                        if(readMM[pr.first]==255) ++n_matched;
                        readMM[pr.first]=(uint8_t)m; matched[pr.first]=1;
                        ppos[pr.first]=pr.second;
                    }
                }
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
        buildPref(true);
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
    {
        mm_total=0; mm_reads=0; std::fill(mm_hist.begin(),mm_hist.end(),0);
        for(uint32_t i=0;i<n;++i) if(readMM[i]!=255){ mm_total+=readMM[i]; ++mm_reads; ++mm_hist[readMM[i]]; }
        size_t zero=mm_hist[0];
        fprintf(stderr,"[MM] placed=%zu total_mismatches=%zu mean=%.2f/read  zero-mm=%zu (%.1f%%)\n",
                mm_reads,mm_total,mm_reads?(double)mm_total/mm_reads:0.0,
                zero,mm_reads?100.0*zero/mm_reads:0.0);
        fprintf(stderr,"[MM] PgRC2 for reference: 1,369,413 mismatches over 988,702 reads = 1.39/read\n");
        fprintf(stderr,"[MM] hist:"); for(int i=0;i<12;++i) fprintf(stderr," %d:%zu",i,mm_hist[i]);
        size_t tail=0; for(int i=12;i<256;++i) tail+=mm_hist[i];
        fprintf(stderr,"  >=12:%zu\n",tail);
    }

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
        // STAGE 87: atomic write. fopen("perm.u32","wb") makes the file
        // VISIBLE to any concurrent reader immediately, but it is not
        // COMPLETE until fclose() -- a real, measured race (see this
        // file's own header) against the combined pipeline's polling loop,
        // which checks existence, not completion. Write to a temp name and
        // rename() onto the final name only after fclose() succeeds:
        // POSIX rename() is atomic on the same filesystem, so a reader can
        // only ever see "not there yet" or "fully there," never partial.
        FILE* f=fopen("perm.u32.tmp","wb");
        for(uint32_t o=0;o<rev.size();++o) if(rev[o]!=UINT32_MAX) fwrite(&rev[o],4,1,f);
        fclose(f);
        rename("perm.u32.tmp","perm.u32");
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
      std::vector<uint16_t>().swap(rlen); }
    { std::vector<uint64_t>().swap(pent); std::vector<uint32_t>().swap(ptab); }
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
    // STAGE 27. Whole-architecture check, not another layer tweak.
    // Stage 17 cut survivors 66,775 -> 15,296 by accepting up to readLen/3 = 50
    // mismatches and taking the FIRST acceptable placement. PgRC2 allows the
    // same 50 but re-scans and keeps the BEST (ReadsMatchers.cpp:315-328, the
    // "better-matches" column). Their accepted matches average 1,369,413
    // mismatches over 988,702 reads = 1.39 per read, coded to 265,900 B.
    //
    // Every mismatch we accept has to be STORED. This progression measures
    // pg literal and order and never counted that stream, so stage 17's win
    // could be partly borrowed against a cost that was never on the books.
    // This measures it before anything else is built on top.
    const size_t MINMEM  = (argc>9)?(size_t)atoi(argv[9]):45;
    const size_t MEMSEED = (MINMEM>14)?((MINMEM-14>32)?32:MINMEM-14):8;
    struct Ref { uint32_t dst, src, len; };
    std::vector<Ref> allrefs;
    size_t rem_main=0, rem_second=0, nm_main=0, nm_second=0;
    {
        // MINMEM was 45 because that is PgRC2's CODER_LEVEL_NORMAL default
        // (pgrc-params.h:145), copied without ever being swept here. With the
        // reference cost finally measured, the trade is computable: a match
        // costs ~40.5 bits (24.5 source + 7.4 gap + 8.6 length = 5.07 B) and
        // saves L * 1.9174/8 bytes of literal, so it pays whenever L >= ~21.
        // 45 is more than twice that.
        //
        // The memory objection dissolves too. STEP = MINMEM - SEED + 1, so
        // shrinking the seed alongside the threshold holds STEP -- and therefore
        // the index size -- constant: MINMEM 30 with a 16-base seed gives the
        // same STEP 15 as MINMEM 45 with 32.
        const size_t STEP     = MINMEM - MEMSEED + 1;
        auto packM=[&](const char* q,uint64_t& o)->bool{
            uint64_t k=0;
            for(size_t i=0;i<MEMSEED;++i){ int v=b2(q[i]); if(v<0) return false; k=(k<<2)|(uint64_t)v; }
            o=k; return true; };
        fprintf(stderr,"[MEM] MINMEM=%zu seed=%zu step=%zu\n",MINMEM,MEMSEED,STEP);
        // Candidates are scanned in index order and the LONGEST is kept, so a cap
        // does not merely bound work -- it can hide the best match entirely. Our
        // matches average 179 bases at MINMEM 45 against copMEM's 250, and a
        // truncated candidate list is the obvious suspect.
        const size_t MAXCAND = (argc>10)?(size_t)atoi(argv[10]):64;
        const bool   LAZY    = (argc>11)?(atoi(argv[11])!=0):true;
        const bool   FWD_SELF=(getenv("FWD_SELF")?atoi(getenv("FWD_SELF"))!=0:true);

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
            for(size_t p=0;p+MEMSEED<=main_pg_end;p+=STEP){
                uint64_t k; if(packM(pg.data()+p,k)) tmp.push_back({k,(uint32_t)p});
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
        // STAGE 32. Every accepted match is a REFERENCE the archive must store:
        // where in the destination it starts, where in the source it points, and
        // how long it is. PG_LITERAL counts only the bases that survive, so that
        // cost has been invisible in every number this progression has reported.
        // PgRC2 pays 177,180 B for exactly these streams (offsets 114,824 +
        // lengths 45,485 + 13,847 + 3,024), and their raw 155,048 B offset stream
        // implies ~38,762 matches against our 59,013 -- so ours is plausibly
        // larger, on the order of the whole remaining deficit.
        //
        // Same failure stage 27 caught with mismatches: an aggressive stage looks
        // free because its cost lands somewhere nobody measures.
        auto parse_range=[&](const char* Q,size_t qlen,Mode mode,size_t lo,size_t hi,
                             std::vector<Ref>& out){
            // STAGE 34: lazy matching.
            // Taking the longest match at every position, left to right, is
            // GREEDY parsing, and the classic LZ result is that greedy is not
            // optimal: accepting a shorter match now can block a much longer one
            // starting a base later, and the fragments cost a reference each.
            // References are our only remaining loss (66,831 against their
            // ~43,190) and each costs ~40 bits, so fewer-and-longer is exactly
            // the objective.
            //
            // Lazy matching is the cheap form of lookahead, as used in DEFLATE:
            // having found a match at qp, look at qp+1, and if that is strictly
            // longer, leave qp as a literal and take the better one. One extra
            // probe per matched position, no extra memory, and nothing from
            // PgRC2 -- their parse is greedy too.
            auto bestAt=[&](const char* Q,size_t qlen,Mode mode,size_t qp,
                            size_t& bsrc)->size_t{
                uint64_t k; bsrc=0;
                if(!packM(Q+qp,k)) return 0;
                const uint32_t idx=lookup(k);
                if(idx==UINT32_MAX) return 0;
                size_t best=0, tried=0;
                for(uint32_t i=idx;i<skey.size()&&skey[i]==k;++i){
                    if(++tried>MAXCAND) break;
                    const size_t s=spos[i];
                    size_t capL=SIZE_MAX;
                    if(mode==SELF_FWD){ if(s>=qp) continue; capL=qp-s; }
                    else if(mode==SELF_RC){ if(s>=qlen-qp) continue; capL=(qlen-qp-s)/2; }
                    if(capL<MINMEM) continue;
                    size_t L=0;
                    while(L<capL && qp+L<qlen && s+L<main_pg_end && Q[qp+L]==pg[s+L]) ++L;
                    if(L>best){ best=L; bsrc=s; }
                }
                return (best>=MINMEM)?best:0;
            };
            size_t qp=lo, lastend=lo;
            while(qp<hi && qp+MINMEM<=qlen){
                uint64_t k;
                if(!packM(Q+qp,k)){ ++qp; continue; }
                const uint32_t idx=lookup(k);
                if(idx==UINT32_MAX){ ++qp; continue; }
                size_t best=0, bestsrc=0, tried=0;
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
                    if(L>best){ best=L; bestsrc=s; }
                }
                if(best>=MINMEM && LAZY && qp+1<hi && qp+1+MINMEM<=qlen){
                    size_t nsrc=0;
                    const size_t nb=bestAt(Q,qlen,mode,qp+1,nsrc);
                    if(nb>best){ ++qp; continue; }   // qp becomes literal; take the longer one next
                }
                if(best>=MINMEM){
                    // STAGE 33: extend BACKWARD before accepting.
                    // The source is sampled every STEP, so a true maximal match
                    // at source ss is only seen at the first sampled position
                    // p >= ss inside it. Extending forward from there captures a
                    // SUFFIX of the real match and abandons the leading p-ss
                    // bases -- about STEP/2 each -- to literal. That is why our
                    // matches average 189 bases against CopMEM's 250: copMEM
                    // reports maximal exact matches, this reported partial ones.
                    //
                    // Recovering them costs nothing: the match count is
                    // unchanged, so the reference streams are unchanged, and the
                    // recovered bases come straight off the literal.
                    size_t cap=qp-lastend;                       // no overlap with the previous match
                    if(bestsrc<cap) cap=bestsrc;                 // and none before the source start
                    if(mode==SELF_FWD){
                        // source must still end at or before the destination
                        const size_t slack=qp-bestsrc-best;      // >=0, enforced by capL above
                        if(slack<cap) cap=slack;
                    }
                    // SELF_RC needs no extra cap: its condition is
                    // bestsrc+best <= qlen-qp-best, which backward extension
                    // leaves unchanged on both sides.
                    size_t b=0;
                    while(b<cap && Q[qp-b-1]==pg[bestsrc-b-1]) ++b;
                    out.push_back({(uint32_t)(qp-b),(uint32_t)(bestsrc-b),(uint32_t)(best+b)});
                    lastend=qp+best; qp+=best;
                }
                else ++qp;
            }
        };
        // Split the destination across threads. The greedy parse is serial by
        // nature (it skips ahead by each accepted match), so what is parallelised
        // is the search, exactly as ARCS's repeat_elim does: every query reads
        // only the shared, already-built index. Slicing can only differ from a
        // single serial parse at the T-1 slice boundaries, bounded by ~45 bytes
        // each -- under 0.005% of a 21 MB pseudogenome.
        // DSTBASE: the cross-match is handed pg.data()+main_pg_end, so its qp is
        // relative to the survivor pg, not to the pseudogenome. Stored raw, a
        // cross-match destination of 0 claims to be the very start of the pg,
        // which is why 12,003 of 58,908 rows appeared to violate src < dst.
        auto run=[&](const char* Q,size_t qlen,Mode mode,std::vector<uint8_t>& consumed,
                     bool RCDEST=false,size_t DSTBASE=0)->size_t{
            unsigned T=std::thread::hardware_concurrency(); if(!T) T=1;
            if(qlen < (1u<<20)) T=1;
            std::vector<std::vector<Ref>> res(T);
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
                for(auto& m:v){
                    for(uint32_t j=0;j<m.len && (size_t)m.dst+j<qlen;++j) consumed[m.dst+j]=1;
                    // An RC pass reports its destination in reverse-complement
                    // coordinates while the source stays in forward ones, so the
                    // two are not comparable until the destination is converted
                    // back: a match at RC position qp of length L covers forward
                    // [qlen-qp-L, qlen-qp). PgRC2 does the same conversion in
                    // correctDestPositionDueToRevComplMatching
                    // (SimplePgMatcher.cpp:58-60); we were storing the raw RC
                    // position, which left src and dst in different coordinate
                    // systems.
                    Ref r=m;
                    if(RCDEST) r.dst=(uint32_t)(qlen-m.dst-m.len);
                    r.dst=(uint32_t)(r.dst+DSTBASE);
                    allrefs.push_back(r);
                }
            }
            return nm;
        };

        // DUMP_LIT: write the surviving literal (everything no MEM reference
        // covers) so it can be entropy-coded and compared against PgRC2's
        // 12,694,903 -> 3,056,474. Literal byte counts are NOT comparable to
        // compressed archive numbers, and the order stream is compressed, so
        // mixing the two without this step would be meaningless.
        //
        // STAGE 87 REAL BUG FOUND AND FIXED: fopen() here makes literal.txt
        // VISIBLE to a concurrent reader immediately, but its content is not
        // COMPLETE until fclose() -- and unlike perm.u32 (a fast, pure
        // write), this file stays open across the ENTIRE MEM-matching
        // compute phase below (both self-match passes), not just a write
        // loop. Measured directly (steady_clock timestamps at fopen and
        // fclose, real yeast_sub.fq data): the gap is ~0.66 SECONDS, not a
        // microsecond fluke. Reproduced the actual failure using the real
        // combined-pipeline's exact polling logic (10ms sleep, existence
        // check) against this exact race: seqpar read literal.txt as
        // completely empty (bases=0, coded=0 B) -- a real, deterministic
        // reproduction of silent data corruption, not a theoretical risk.
        // Every earlier "beat SPRING" measurement in this project won this
        // race by timing luck on this specific machine under normal load;
        // it was never actually safe. Same atomic-rename fix as perm.u32:
        // write to a temp name, rename() onto the final name only after
        // fclose() -- POSIX rename() is atomic, so a concurrent reader can
        // only ever see "not there yet" or "fully there."
        FILE* litf = getenv("DUMP_LIT") ? fopen("literal.txt.tmp","wb") : nullptr;
        // Second pg (their LQ+N) against the main pg, forward then reverse complement.
        {
            const size_t qlen=pg.size()-main_pg_end;
            if(qlen>=MINMEM){
                std::string Q(pg,main_pg_end,qlen);
                std::vector<uint8_t> c(qlen,0), cr(qlen,0);
                nm_second =run(Q.data(),qlen,CROSS,c,false,main_pg_end);
                std::string R=Q; rc_inplace(R);
                nm_second+=run(R.data(),qlen,CROSS,cr,true,main_pg_end);
                for(size_t i=0;i<qlen;++i) if(cr[i]) c[qlen-1-i]=1;
                for(size_t i=0;i<qlen;++i) if(c[i]) ++rem_second;
                if(litf) for(size_t i=0;i<qlen;++i) if(!c[i]) fputc(Q[i],litf);
            }
        }
        // Main pg against itself, forward then reverse complement.
        {
            const size_t qlen=main_pg_end;
            std::vector<uint8_t> c(qlen,0), cr(qlen,0);
            // STAGE 36. PgRC2's stage 7 self-match is REVERSE-COMPLEMENT ONLY.
            // exactMatchPg (SimplePgMatcher.cpp:33-35) builds
            // reverseComplement(destPg) and matches that against srcPg whenever
            // destPgIsSrcPg, and revComplMatching is true for every call in
            // matchPgsInPg. Their parameter is even named
            // minimalReverseComplementedRepeatLength (-p). They never look for
            // forward repeats inside the pseudogenome.
            //
            // We always did both, which is where 65,994 matches against their
            // 43,190 comes from -- and the arithmetic says those extra forward
            // matches are a NET LOSS: they take our literal 177,488 bases below
            // theirs (~42,370 B once coded) while costing 151,332 B of extra
            // references. A match only pays if L * bits_per_base / 8 exceeds the
            // ~40 bits its reference costs, and forward self-matches in a
            // pseudogenome are mostly short.
            //
            // FWD_SELF=0 turns forward self-matching off to measure that.
            if(FWD_SELF) nm_main =run(pg.data(),qlen,SELF_FWD,c);
            std::string R(pg,0,qlen); rc_inplace(R);
            nm_main+=run(R.data(),qlen,SELF_RC,cr,true);
            for(size_t i=0;i<qlen;++i) if(cr[i]) c[qlen-1-i]=1;
            for(size_t i=0;i<qlen;++i) if(c[i]) ++rem_main;
            if(litf) for(size_t i=0;i<qlen;++i) if(!c[i]) fputc(pg[i],litf);
        }
        if(litf){ fclose(litf); rename("literal.txt.tmp","literal.txt"); fprintf(stderr,"[LIT] literal.txt written\n"); }
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
    // Reference streams, in the form an archive would actually store: matches in
    // destination order, the start as a gap from the previous match's end (so the
    // values are small and positive), the source as an absolute pg offset, and
    // the length as a delta above MINMEM.
    {
        std::sort(allrefs.begin(),allrefs.end(),
                  [](const Ref&a,const Ref&b){ return a.dst<b.dst; });
        std::vector<uint8_t> gaps,srcs,lens;
        auto vint=[](std::vector<uint8_t>& o,uint64_t v){
            while(true){ uint8_t b=v&0x7f; v>>=7; o.push_back(b|(v?0x80:0)); if(!v) break; } };
        // Source offsets are the dominant component (182 KB of 300 KB, 24.7 bits
        // each = log2(pg length), i.e. the naive floor with no locality used).
        // Consecutive matches in destination order often point at nearby sources,
        // so a zigzag delta is also emitted and the smaller kept -- the same
        // keep-the-better discipline used for the aux layouts in ARCS.
        std::vector<uint8_t> srcd;
        uint64_t prev=0; int64_t psrc=0;
        for(const Ref& r:allrefs){
            const uint64_t g=(r.dst>=prev)?(r.dst-prev):0;
            vint(gaps,g); vint(srcs,r.src); vint(lens,(uint64_t)(r.len-MINMEM));
            // Repeat DISTANCE, not absolute position. Every self-match points
            // backwards, so dst-src is the LZ77 match distance, and distances
            // are strongly skewed towards small values where absolute offsets
            // into a 23 Mbase pg are near-uniform. (An earlier test used the
            // delta between CONSECUTIVE sources, which is a different quantity
            // and was worse than absolute -- sources have no locality with each
            // other, but each source is close to its own destination.)
            const int64_t d=(int64_t)r.dst-(int64_t)r.src;
            vint(srcd,(uint64_t)((d<<1)^(d>>63)));
            psrc=(int64_t)r.src;
            prev=(uint64_t)r.dst+r.len;
        }
        { FILE* g=fopen("mem_srcdelta.bin","wb"); fwrite(srcd.data(),1,srcd.size(),g); fclose(g); }
        { FILE* t=fopen("mem_triples.bin","wb");
          for(const Ref& r:allrefs){ fwrite(&r.dst,4,1,t); fwrite(&r.src,4,1,t); fwrite(&r.len,4,1,t); }
          fclose(t); }
        FILE* f;
        f=fopen("mem_gaps.bin","wb"); fwrite(gaps.data(),1,gaps.size(),f); fclose(f);
        f=fopen("mem_srcs.bin","wb"); fwrite(srcs.data(),1,srcs.size(),f); fclose(f);
        f=fopen("mem_lens.bin","wb"); fwrite(lens.data(),1,lens.size(),f); fclose(f);
        fprintf(stderr,"[REF] matches=%zu  raw gaps=%zu srcs=%zu lens=%zu  total_raw=%zu\n",
                allrefs.size(),gaps.size(),srcs.size(),lens.size(),
                gaps.size()+srcs.size()+lens.size());
        fprintf(stderr,"[REF] PgRC2 pays 177,180 B coded for the same thing\n");
    }
    // DUMP_PG: the pseudogenome BEFORE MEM removal. Needed to test whether a
    // context model with a match model can capture the long repeats implicitly,
    // which would remove the reference streams entirely rather than shrink them.
    if(getenv("DUMP_PG")){
        FILE* f=fopen("pg_full.txt","wb"); fwrite(pg.data(),1,pg.size(),f); fclose(f);
        fprintf(stderr,"[PG] pg_full.txt written: %zu bases\n",pg.size());
    }
    printf("PG_LEN %zu\n",pg.size());
    printf("PG_LITERAL %zu\n",lit_main+lit_second);
    return 0;
}
