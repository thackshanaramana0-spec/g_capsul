// STAGE 72 -- bounded-memory streaming pipeline, names coder.
//
// Investigated Genozip's REAL dispatcher (cloned divonlan/genozip,
// src/dispatcher.c, read in full -- Genozip publishes source specifically
// so users can verify their own decompression, reading it for
// understanding is the same legitimate use already made of SPRING's
// licensed source all session). Found something genuinely different from
// both SPRING and our own stage 69/70/71: Genozip's parallelism is NOT
// "split the whole file into blocks up front, dispatch by index" (that's
// what SPRING's OpenMP `block_num += num_thr` does, and what stage 69-71
// do too). It's a single main-thread priority loop, every iteration:
//   1. a VBlock is ready and a compute thread is free -> dispatch it
//   2. else a finished VBlock is waiting -> JOIN it (out-of-order allowed,
//      frees its memory immediately, doesn't wait for a specific index)
//   3. else a pool slot is free and input isn't exhausted -> read/generate
//      the NEXT VBlock from the input stream
//   4. else -> sleep briefly
// The real, concrete consequence: memory is bounded by (pool_size x
// VBlock size), NOT by file size -- Genozip never holds the whole input in
// RAM. Stage 69-71 all read the ENTIRE file into a std::vector<std::string>
// before any threading starts -- a real scalability gap for this project
// specifically, since the actual benchmark includes a ~15 GB dataset
// (T. cacao, SRR870667) that would sit fully materialized in memory before
// compression even begins.
//
// Genozip's own priority-loop is hand-rolled coordination logic, likely
// needed because IT also has strict output-ordering requirements for other
// file formats. Our case doesn't have that constraint: each compressed
// block is self-contained (block index already labels it, same as stage
// 69's atomic-counter dispatch), so out-of-order JOIN doesn't need special
// handling here -- the standard, correct concurrency primitive for the
// SAME memory-bounding property (read-ahead limited, backpressure once
// full) is a bounded blocking queue: a reader thread pushes fixed-size
// chunks read directly from the file (blocking once the queue is full),
// worker threads pop and compress. Simpler than reimplementing Genozip's
// raw priority loop, same real property.
//
// Uses stage 69's EXACT IdModels/compress_id/decompress_id (fresh model
// per block, SPRING-style) -- this stage isolates the MEMORY ARCHITECTURE
// as the only variable being tested, not mixed with stage 70's dictionary
// idea.
//
// REAL, MEASURED RESULT (names_1m.txt, 69 MB, 1M names, 12 threads,
// block=20,000-100,000): a genuine memory-vs-throughput tradeoff, not a
// free win.
//   peak RSS:  stage 69 (load-all)      400 MB
//              stage 72 (bounded queue) 216 MB   (-46%, real, /usr/bin/time -v)
//   time:      stage 69 read+encode     ~0.16s (estimated: 0.059s file read,
//                                        measured separately, + 0.099s
//                                        encode-only, measured after read)
//              stage 72 read+encode     ~0.21s (measured directly -- reader
//                                        thread and worker threads overlap,
//                                        so read time is INSIDE this number,
//                                        unlike stage 69 where it's hidden
//                                        before the timer starts)
// First timing comparison attempt directly compared stage 69's "encode
// only" (excludes file read, timer starts after full load) against stage
// 72's "encode only" (INCLUDES file read, since reading and encoding
// genuinely overlap) -- an unfair, apples-to-oranges comparison, caught by
// separately timing raw file-read (0.059s) before drawing a conclusion.
// Even after correcting for that, stage 72 is genuinely ~30-50% slower on
// THIS file size -- real mutex/condvar synchronization overhead on every
// chunk push/pop, not a measurement artifact. Expected, disclosed cost:
// worth paying specifically when the whole file doesn't comfortably fit in
// RAM to begin with (the real motivating case: T. cacao, ~15 GB, in the
// actual ARCS benchmark) -- not a win on a file this small, where stage
// 69's simpler design was never going to run out of memory anyway.
//
//   g++ -O3 -march=native -pthread -o namebq 72_names_boundedqueue.cpp

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <unordered_map>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <algorithm>

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
struct Model {
    std::vector<uint32_t> f; uint32_t N;
    Model(uint32_t n=256):f(n,1),N(n){}
    void enc(RangeEnc& rc,uint32_t v){
        if(v>=N) v=N-1;
        uint32_t tot=0; for(uint32_t i=0;i<N;++i) tot+=f[i];
        uint32_t lo=0; for(uint32_t i=0;i<v;++i) lo+=f[i];
        rc.encode(lo,lo+f[v],tot);
        f[v]+=24; if(tot+24>60000){ for(uint32_t i=0;i<N;++i) f[i]=(f[i]>>1)|1; }
    }
    uint32_t dec(RangeDec& rc){
        uint32_t tot=0; for(uint32_t i=0;i<N;++i) tot+=f[i];
        uint32_t target=rc.getFreq(tot);
        uint32_t v=0; uint32_t lo=0; while(lo+f[v]<=target){ lo+=f[v]; ++v; }
        rc.decodeUpdate(lo,lo+f[v]);
        f[v]+=24; if(tot+24>60000){ for(uint32_t i=0;i<N;++i) f[i]=(f[i]>>1)|1; }
        return v;
    }
};

enum TokType { ID_ALPHA, ID_DIGIT, ID_CHAR, ID_MATCH, ID_ZEROS, ID_DELTA, ID_END, ID_ZDELTA };
static const uint32_t MAXTOK=1024;

struct ValueDict {
    std::unordered_map<uint32_t,uint32_t> val2sym;
    std::vector<uint32_t> sym2val;
    std::vector<uint32_t> freq;
    uint32_t total;
    ValueDict():freq(1,1),total(1){}
    uint32_t lookup(uint32_t v) const { auto it=val2sym.find(v); return it==val2sym.end()?0:it->second; }
    void encSym(RangeEnc& rc,uint32_t sym){
        uint32_t lo=0; for(uint32_t i=0;i<sym;++i) lo+=freq[i];
        rc.encode(lo,lo+freq[sym],total);
        freq[sym]+=16; total+=16;
        if(total>60000){ total=0; for(auto& f:freq){ f=(f>>1)|1; total+=f; } }
    }
    uint32_t decSym(RangeDec& rc){
        uint32_t target=rc.getFreq(total);
        uint32_t sym=0; uint32_t lo=0; while(lo+freq[sym]<=target){ lo+=freq[sym]; ++sym; }
        rc.decodeUpdate(lo,lo+freq[sym]);
        freq[sym]+=16; total+=16;
        if(total>60000){ total=0; for(auto& f:freq){ f=(f>>1)|1; total+=f; } }
        return sym;
    }
    void registerNew(uint32_t v){
        uint32_t sym=(uint32_t)freq.size();
        val2sym[v]=sym; sym2val.push_back(v); freq.push_back(1); total+=1;
    }
};

struct IdModels {
    std::vector<Model> token_type, alpha_len, alpha_value, chars, zero_run, delta;
    std::vector<Model> integer;
    std::vector<Model> zdelta_hi, zdelta_lo;
    std::vector<ValueDict> valdict;
    std::vector<uint32_t> hit, seen;
    std::vector<uint32_t> dict_hit, dict_seen;
    IdModels():
        token_type(MAXTOK,Model(8)), alpha_len(MAXTOK,Model(256)),
        alpha_value(MAXTOK,Model(128)), chars(MAXTOK,Model(128)),
        zero_run(MAXTOK,Model(256)), delta(MAXTOK,Model(256)),
        integer(MAXTOK*4,Model(256)),
        zdelta_hi(MAXTOK,Model(256)), zdelta_lo(MAXTOK,Model(256)),
        hit(MAXTOK,0), seen(MAXTOK,0), valdict(MAXTOK),
        dict_hit(MAXTOK,0), dict_seen(MAXTOK,0) {}
};

void compress_id(RangeEnc& enc, IdModels& m, const char* id, std::string& prev_id,
                  std::array<uint32_t,MAXTOK>& prev_tok_ptr){
    uint32_t token_len=0, match_len=0, token_ctr=0, i=0;
    const char* id_ptr=id;
    const char* prevbuf = prev_id.c_str();
    const size_t prevlen = prev_id.size();
    auto prevc=[&](uint32_t off)->char{ return off<prevlen ? prevbuf[off] : 0; };

    while(*id_ptr){
        token_len=0;
        match_len = (*id_ptr==prevc(prev_tok_ptr[token_ctr]+token_len)); token_len=1;
        const char* id_ptr_tok=id_ptr+1;
        if(isalpha((unsigned char)*id_ptr)){
            while(isalpha((unsigned char)*id_ptr_tok)){
                match_len += (*id_ptr_tok==prevc(prev_tok_ptr[token_ctr]+token_len));
                ++token_len; ++id_ptr_tok;
            }
            if(match_len==token_len && !isalpha((unsigned char)prevc(prev_tok_ptr[token_ctr]+token_len))){
                m.token_type[token_ctr].enc(enc,ID_MATCH);
            } else {
                m.token_type[token_ctr].enc(enc,ID_ALPHA);
                m.alpha_len[token_ctr].enc(enc,token_len);
                for(uint32_t k=0;k<token_len;++k) m.alpha_value[token_ctr].enc(enc,(unsigned char)id_ptr[k]&0x7f);
            }
        } else if(*id_ptr=='0'){
            while(*id_ptr_tok=='0'){
                match_len += ('0'==prevc(prev_tok_ptr[token_ctr]+token_len));
                ++token_len; ++id_ptr_tok;
            }
            if(match_len==token_len && prevc(prev_tok_ptr[token_ctr]+token_len)!='0'){
                m.token_type[token_ctr].enc(enc,ID_MATCH);
            } else {
                m.token_type[token_ctr].enc(enc,ID_ZEROS);
                m.zero_run[token_ctr].enc(enc,token_len);
            }
        } else if(isdigit((unsigned char)*id_ptr)){
            uint32_t digit_value=(uint32_t)(*id_ptr-'0');
            bool prev_is_digit=true; uint32_t prev_digit=0;
            if(!prev_id.empty()){
                char pc = prevc(prev_tok_ptr[token_ctr]+token_len-1);
                if(isdigit((unsigned char)pc) && pc!='0') prev_digit=(uint32_t)(pc-'0');
                else prev_is_digit=false;
            } else prev_is_digit=false;
            if(prev_is_digit){
                uint32_t tmp=1;
                while(isdigit((unsigned char)prevc(prev_tok_ptr[token_ctr]+tmp)) && prev_digit<(1u<<28)){
                    prev_digit = prev_digit*10 + (uint32_t)(prevc(prev_tok_ptr[token_ctr]+tmp)-'0');
                    ++tmp;
                }
            }
            while(isdigit((unsigned char)*id_ptr_tok) && digit_value<(1u<<28)){
                digit_value = digit_value*10 + (uint32_t)(*id_ptr_tok-'0');
                match_len += (*id_ptr_tok==prevc(prev_tok_ptr[token_ctr]+token_len));
                ++token_len; ++id_ptr_tok;
            }
            int64_t delta = (int64_t)digit_value - (int64_t)prev_digit;
            const bool can_delta  = prev_is_digit && delta>0 && delta<256;
            const bool can_zdelta = prev_is_digit && delta>=-32768 && delta<=32767;
            if(prev_is_digit && match_len==token_len && !isdigit((unsigned char)prevc(prev_tok_ptr[token_ctr]+token_len))){
                m.token_type[token_ctr].enc(enc,ID_MATCH);
            } else {
                if(prev_is_digit){
                    ++m.seen[token_ctr];
                    if(delta>=-2048 && delta<=2048) ++m.hit[token_ctr];
                }
                const uint32_t seen=m.seen[token_ctr], hit=m.hit[token_ctr];
                const bool trust_wide = (seen>=20) && (hit*10 >= seen*3);
                if(can_delta){
                    m.token_type[token_ctr].enc(enc,ID_DELTA);
                    m.delta[token_ctr].enc(enc,(uint32_t)delta);
                } else if(can_zdelta && trust_wide){
                    uint32_t z=(uint32_t)((delta<<1) ^ (delta>>63));
                    m.token_type[token_ctr].enc(enc,ID_ZDELTA);
                    m.zdelta_hi[token_ctr].enc(enc,(z>>8)&0xff);
                    m.zdelta_lo[token_ctr].enc(enc,z&0xff);
                } else {
                    m.token_type[token_ctr].enc(enc,ID_DIGIT);
                    ValueDict& vd = m.valdict[token_ctr];
                    uint32_t sym = vd.lookup(digit_value);
                    const uint32_t dseen=m.dict_seen[token_ctr], dhit=m.dict_hit[token_ctr];
                    ++m.dict_seen[token_ctr]; if(sym!=0) ++m.dict_hit[token_ctr];
                    const bool trust_dict = (dseen<30) || (dhit*2 >= dseen);
                    if(trust_dict){
                        vd.encSym(enc,sym);
                        if(sym==0){
                            m.integer[token_ctr*4+0].enc(enc,(digit_value>>0)&0xff);
                            m.integer[token_ctr*4+1].enc(enc,(digit_value>>8)&0xff);
                            m.integer[token_ctr*4+2].enc(enc,(digit_value>>16)&0xff);
                            m.integer[token_ctr*4+3].enc(enc,(digit_value>>24)&0xff);
                            vd.registerNew(digit_value);
                        }
                    } else {
                        m.integer[token_ctr*4+0].enc(enc,(digit_value>>0)&0xff);
                        m.integer[token_ctr*4+1].enc(enc,(digit_value>>8)&0xff);
                        m.integer[token_ctr*4+2].enc(enc,(digit_value>>16)&0xff);
                        m.integer[token_ctr*4+3].enc(enc,(digit_value>>24)&0xff);
                        if(sym==0) vd.registerNew(digit_value);
                    }
                }
            }
        } else {
            if(match_len==token_len){
                m.token_type[token_ctr].enc(enc,ID_MATCH);
            } else {
                m.token_type[token_ctr].enc(enc,ID_CHAR);
                m.chars[token_ctr].enc(enc,(unsigned char)*id_ptr&0x7f);
            }
        }
        prev_tok_ptr[token_ctr]=i;
        i+=token_len; id_ptr=id_ptr_tok; token_ctr++;
        if(token_ctr>=MAXTOK-1) break;
    }
    prev_id.assign(id);
    m.token_type[token_ctr].enc(enc,ID_END);
    for(uint32_t k=token_ctr+1;k<MAXTOK;++k) prev_tok_ptr[k]=0;
}

std::string decompress_id(RangeDec& dec, IdModels& m, std::string& prev_id,
                           std::array<uint32_t,MAXTOK>& prev_tok_ptr,
                           std::array<uint32_t,MAXTOK>& prev_tok_len){
    std::string id; id.reserve(64);
    uint32_t token_ctr=0, i=0;
    const char* prevbuf = prev_id.c_str();
    const size_t prevlen = prev_id.size();
    auto prevc=[&](uint32_t off)->char{ return off<prevlen ? prevbuf[off] : 0; };
    for(;;){
        uint32_t tok = m.token_type[token_ctr].dec(dec);
        if(tok==ID_END) break;
        uint32_t token_len=0;
        if(tok==ID_MATCH){
            uint32_t off=prev_tok_ptr[token_ctr], len=prev_tok_len[token_ctr];
            for(uint32_t k=0;k<len;++k) id.push_back(prevc(off+k));
            token_len=len;
        } else if(tok==ID_ALPHA){
            token_len = m.alpha_len[token_ctr].dec(dec);
            for(uint32_t k=0;k<token_len;++k) id.push_back((char)m.alpha_value[token_ctr].dec(dec));
        } else if(tok==ID_DIGIT){
            ValueDict& vd = m.valdict[token_ctr];
            const uint32_t dseen=m.dict_seen[token_ctr], dhit=m.dict_hit[token_ctr];
            const bool trust_dict = (dseen<30) || (dhit*2 >= dseen);
            uint32_t v, sym;
            if(trust_dict){
                sym = vd.decSym(dec);
                if(sym==0){
                    v = 0;
                    v |= m.integer[token_ctr*4+0].dec(dec)<<0;
                    v |= m.integer[token_ctr*4+1].dec(dec)<<8;
                    v |= m.integer[token_ctr*4+2].dec(dec)<<16;
                    v |= m.integer[token_ctr*4+3].dec(dec)<<24;
                    vd.registerNew(v);
                } else {
                    v = vd.sym2val[sym-1];
                }
            } else {
                v = 0;
                v |= m.integer[token_ctr*4+0].dec(dec)<<0;
                v |= m.integer[token_ctr*4+1].dec(dec)<<8;
                v |= m.integer[token_ctr*4+2].dec(dec)<<16;
                v |= m.integer[token_ctr*4+3].dec(dec)<<24;
                sym = vd.lookup(v);
                if(sym==0) vd.registerNew(v);
            }
            ++m.dict_seen[token_ctr]; if(sym!=0) ++m.dict_hit[token_ctr];
            char buf[16]; int L=snprintf(buf,sizeof(buf),"%u",v);
            id.append(buf,(size_t)L); token_len=(uint32_t)L;
        } else if(tok==ID_DELTA){
            uint32_t delta = m.delta[token_ctr].dec(dec);
            uint32_t off=prev_tok_ptr[token_ctr], len=prev_tok_len[token_ctr];
            std::string prevtokstr; for(uint32_t k=0;k<len;++k) prevtokstr.push_back(prevc(off+k));
            uint32_t prevval=(uint32_t)atoi(prevtokstr.c_str());
            uint32_t v=prevval+delta;
            char buf[16]; int L=snprintf(buf,sizeof(buf),"%u",v);
            id.append(buf,(size_t)L); token_len=(uint32_t)L;
        } else if(tok==ID_ZDELTA){
            uint32_t hi=m.zdelta_hi[token_ctr].dec(dec);
            uint32_t lo=m.zdelta_lo[token_ctr].dec(dec);
            uint32_t z=(hi<<8)|lo;
            int64_t delta = (int64_t)(z>>1) ^ -(int64_t)(z&1);
            uint32_t off=prev_tok_ptr[token_ctr], len=prev_tok_len[token_ctr];
            std::string prevtokstr; for(uint32_t k=0;k<len;++k) prevtokstr.push_back(prevc(off+k));
            int64_t prevval=atoll(prevtokstr.c_str());
            int64_t v=prevval+delta;
            char buf[16]; int L=snprintf(buf,sizeof(buf),"%lld",(long long)v);
            id.append(buf,(size_t)L); token_len=(uint32_t)L;
        } else if(tok==ID_ZEROS){
            token_len = m.zero_run[token_ctr].dec(dec);
            for(uint32_t k=0;k<token_len;++k) id.push_back('0');
        } else if(tok==ID_CHAR){
            id.push_back((char)m.chars[token_ctr].dec(dec));
            token_len=1;
        }
        prev_tok_ptr[token_ctr]=i; prev_tok_len[token_ctr]=token_len;
        i+=token_len; ++token_ctr;
        if(token_ctr>=MAXTOK-1) break;
    }
    prev_id=id;
    if(token_ctr==0) for(uint32_t k=0;k<MAXTOK;++k) prev_tok_ptr[k]=0;
    return id;
}

// ---- STAGE 72: bounded blocking queue, streaming reader/writer ----

struct Chunk {
    size_t block_idx;
    std::vector<std::string> lines; // BLOCK read lines, NOT the whole file
    bool is_end=false;
};

// STAGE 85 REAL FIX: identical bug to stage 73/84 (quality), same
// BoundedQueue design copied from the same pattern -- fixed the same way,
// verified the same way. See stage 84's header for the full diagnostic
// story (real /proc/self/smaps_rollup instrumentation found QCAP was
// never actually the bound; a worker held its popped chunk's raw data
// live, uncounted against QCAP, for its whole processing time). Real fix:
// a chunk's slot only frees when a worker calls release() after fully
// finishing it, not merely once popped -- bounds TOTAL in-flight work
// (queued + being-processed) to QCAP regardless of NTHREADS, block count,
// or file size.
template<typename T>
struct BoundedQueue {
    std::deque<T> q;
    size_t cap;
    std::mutex mu;
    std::condition_variable cv_push, cv_pop;
    bool closed=false;
    size_t inflight=0;   // queued + currently held by a worker (not yet released)
    BoundedQueue(size_t c):cap(c){}
    void push(T v){
        std::unique_lock<std::mutex> lk(mu);
        cv_push.wait(lk,[&]{ return inflight<cap; }); // BACKPRESSURE: bounds TOTAL in-flight, not just queue depth
        ++inflight;
        q.push_back(std::move(v));
        cv_pop.notify_one();
    }
    bool pop(T& out){
        std::unique_lock<std::mutex> lk(mu);
        cv_pop.wait(lk,[&]{ return !q.empty() || closed; });
        if(q.empty() && closed) return false;
        out=std::move(q.front()); q.pop_front();
        // Deliberately do NOT decrement inflight or notify cv_push here --
        // the chunk's data is still live in the caller's hands until
        // release() is called.
        return true;
    }
    void release(){
        std::unique_lock<std::mutex> lk(mu);
        --inflight;
        cv_push.notify_one();
    }
    void close(){ std::unique_lock<std::mutex> lk(mu); closed=true; cv_pop.notify_all(); }
};

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s names.txt [block_size=100000] [nthreads=nproc] [queue_cap=nthreads]\n",argv[0]); return 1; }
    const size_t BLOCK=(argc>2)?(size_t)atoll(argv[2]):100000;
    const unsigned NTHREADS=(argc>3)?(unsigned)atoi(argv[3]):std::thread::hardware_concurrency();
    // Default changed from 2*NTHREADS to 1*NTHREADS now that QCAP actually
    // bounds real memory -- see stage 84's header for why.
    const size_t QCAP=(argc>4)?(size_t)atoll(argv[4]):(size_t)NTHREADS;

    fprintf(stderr,"block=%zu threads=%u queue_cap=%zu (bounded -- reader never gets more than QCAP blocks ahead)\n",BLOCK,NTHREADS,QCAP);

    BoundedQueue<Chunk> workQ(QCAP);
    std::vector<std::vector<uint8_t>> blockOut; // grows as needed; final size known after reader finishes
    std::mutex outMu;

    auto t0=std::chrono::steady_clock::now();

    // READER thread: streams the file line-by-line, never materializes the
    // whole file -- this is the real property being tested (bounded RAM).
    std::thread reader([&](){
        FILE* f=fopen(argv[1],"r"); if(!f){perror("open");exit(1);}
        char buf[2048];
        size_t block_idx=0;
        Chunk cur; cur.block_idx=0; cur.lines.reserve(BLOCK);
        while(fgets(buf,sizeof(buf),f)){
            size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
            cur.lines.emplace_back(buf,L);
            if(cur.lines.size()>=BLOCK){
                workQ.push(std::move(cur));
                ++block_idx;
                cur=Chunk(); cur.block_idx=block_idx; cur.lines.reserve(BLOCK);
            }
        }
        if(!cur.lines.empty()){ workQ.push(std::move(cur)); ++block_idx; }
        fclose(f);
        workQ.close();
        { std::lock_guard<std::mutex> lk(outMu); blockOut.resize(block_idx); }
    });

    // WORKER threads: pop a chunk, compress it (fresh model, stage 69
    // style), store result indexed by block_idx.
    std::vector<std::thread> workers;
    size_t total_names=0; std::mutex cntMu;
    for(unsigned t=0;t<NTHREADS;++t){
        workers.emplace_back([&](){
            Chunk c;
            while(workQ.pop(c)){
                IdModels m;
                RangeEnc enc; enc.out.reserve(c.lines.size()*4);
                std::string prev_id; std::array<uint32_t,MAXTOK> prev_ptr{}; prev_ptr.fill(0);
                for(auto& nm: c.lines) compress_id(enc, m, nm.c_str(), prev_id, prev_ptr);
                enc.flush();
                { std::lock_guard<std::mutex> lk(outMu);
                  if(blockOut.size()<=c.block_idx) blockOut.resize(c.block_idx+1);
                  blockOut[c.block_idx]=std::move(enc.out); }
                { std::lock_guard<std::mutex> lk(cntMu); total_names+=c.lines.size(); }
                std::vector<std::string>().swap(c.lines);   // free the raw chunk data NOW, not at next pop()
                workQ.release();   // the real fix: only now does this chunk's slot free up
            }
        });
    }
    reader.join();
    for(auto& w:workers) w.join();

    auto t1=std::chrono::steady_clock::now();
    size_t bytes=0; for(auto& b:blockOut) bytes+=b.size();
    fprintf(stderr,"ENCODE ONLY (bounded streaming, %u threads): %.3f s\n", NTHREADS,
            std::chrono::duration<double>(t1-t0).count());

    // Round trip: re-stream the original file exactly
    // like the reader did, decode each block with the matching count, and
    // compare directly -- avoids needing to also encode counts into the
    // stream just for this test harness.
    bool ok=true;
    {
        FILE* f=fopen(argv[1],"r");
        char buf[2048];
        std::vector<std::string> cur; cur.reserve(BLOCK);
        size_t b=0;
        auto verify_block=[&](std::vector<std::string>& lines, size_t bidx)->bool{
            if(bidx>=blockOut.size()) return false;
            IdModels m;
            RangeDec dec; dec.init(blockOut[bidx].data(),blockOut[bidx].size());
            std::string dprev; std::array<uint32_t,MAXTOK> dprev_ptr{}, dprev_len{};
            dprev_ptr.fill(0); dprev_len.fill(0);
            for(auto& want: lines){
                std::string got = decompress_id(dec, m, dprev, dprev_ptr, dprev_len);
                if(got!=want){ fprintf(stderr,"MISMATCH block %zu\n",bidx); return false; }
            }
            return true;
        };
        while(fgets(buf,sizeof(buf),f)){
            size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
            cur.emplace_back(buf,L);
            if(cur.size()>=BLOCK){
                if(!verify_block(cur,b)){ ok=false; break; }
                ++b; cur.clear();
            }
        }
        if(ok && !cur.empty()) ok = verify_block(cur,b);
        fclose(f);
    }

    printf("block=%zu  threads=%u  queue_cap=%zu  names=%zu  coded=%zu B  %.3f B/name  round trip: %s\n",
           BLOCK,NTHREADS,QCAP,total_names,bytes,bytes/(double)total_names, ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
