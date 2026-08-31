// STAGE 69 -- threaded, block-chunked names coder. SPRING's real
// parallelism design, confirmed by reading reorder_compress_quality_id.cpp:
// fixed-size `num_reads_per_block`, OpenMP thread pool, each block gets an
// INDEPENDENT, FRESH model -- no cross-block context. Genozip's real design
// (confirmed via its published architecture, Bioinformatics 2021) is
// smarter: threads CLONE the current shared dictionary, encode against the
// clone, then MERGE new observations back -- avoiding the fresh-model
// ratio cost SPRING pays at every block boundary. That is stage 70 (harder:
// merge order must stay deterministic between encoder/decoder or the
// bitstream breaks). This stage is the simpler, first, real baseline:
// independent blocks, no shared state -- matching SPRING's confirmed real
// approach exactly.
//
// compress_id/decompress_id/IdModels/ValueDict/Model/RangeEnc/RangeDec are
// copied VERBATIM from stage 66 (66_dict_gate.cpp) -- parallelism here is a
// call-site/framing change, not a new algorithm. Each block resets
// prev_id="" and prev_tok_ptr={0} (identical to how stage 66's main()
// already initializes the very first ID) -- a block boundary is just
// "another first ID", no new code path, no re-derivation from memory.
//
//   g++ -O3 -march=native -pthread -o namepar 69_names_blockpar.cpp

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
#include <atomic>

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
    double cost(uint32_t v) const {
        if(v>=N) v=N-1;
        uint32_t tot=0; for(uint32_t i=0;i<N;++i) tot+=f[i];
        double p=(double)f[v]/(double)tot;
        return -std::log2(p);
    }
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

uint32_t compute_num_digits(uint32_t x){
    if(x<10) return 1; if(x<100) return 2; if(x<1000) return 3; if(x<10000) return 4;
    if(x<100000) return 5; if(x<1000000) return 6; if(x<10000000) return 7;
    if(x<100000000) return 8; if(x<1000000000) return 9; return 10;
}

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

// ---- STAGE 69: block-parallel wrapper (new code, not from stage 66) ----

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s names.txt [block_size=100000] [nthreads=nproc]\n",argv[0]); return 1; }
    const size_t BLOCK = (argc>2)?(size_t)atoll(argv[2]):100000;
    const unsigned NTHREADS = (argc>3)?(unsigned)atoi(argv[3]):std::thread::hardware_concurrency();

    std::vector<std::string> names;
    { FILE* f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
      char buf[2048];
      while(fgets(buf,sizeof(buf),f)){ size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
          names.emplace_back(buf,L); }
      fclose(f); }
    const size_t n=names.size();
    const size_t nblocks = (n+BLOCK-1)/BLOCK;
    fprintf(stderr,"names=%zu  block=%zu  blocks=%zu  threads=%u\n",n,BLOCK,nblocks,NTHREADS);

    std::vector<std::vector<uint8_t>> blockOut(nblocks);
    auto t0=std::chrono::steady_clock::now();
    {
        std::atomic<size_t> next{0};
        auto worker=[&](){
            size_t b;
            while((b=next.fetch_add(1))<nblocks){
                size_t s=b*BLOCK, e=std::min(n,s+BLOCK);
                IdModels m;
                RangeEnc enc; enc.out.reserve((e-s)*4);
                std::string prev_id; std::array<uint32_t,MAXTOK> prev_ptr{}; prev_ptr.fill(0);
                for(size_t r=s;r<e;++r) compress_id(enc, m, names[r].c_str(), prev_id, prev_ptr);
                enc.flush();
                blockOut[b]=std::move(enc.out);
            }
        };
        std::vector<std::thread> pool;
        for(unsigned t=0;t<NTHREADS;++t) pool.emplace_back(worker);
        for(auto& th:pool) th.join();
    }
    auto t1=std::chrono::steady_clock::now();
    size_t bytes=0; for(auto& b:blockOut) bytes+=b.size();
    fprintf(stderr,"ENCODE ONLY (parallel, %u threads): %.3f s\n", NTHREADS,
            std::chrono::duration<double>(t1-t0).count());

    // round trip: decode every block back, in parallel, and verify
    std::vector<std::vector<std::string>> blockNames(nblocks);
    std::atomic<bool> ok{true};
    {
        std::atomic<size_t> next{0};
        auto worker=[&](){
            size_t b;
            while((b=next.fetch_add(1))<nblocks){
                size_t s=b*BLOCK, e=std::min(n,s+BLOCK);
                IdModels m;
                RangeDec dec; dec.init(blockOut[b].data(),blockOut[b].size());
                std::string dprev; std::array<uint32_t,MAXTOK> dprev_ptr{}, dprev_len{};
                dprev_ptr.fill(0); dprev_len.fill(0);
                std::vector<std::string> got; got.reserve(e-s);
                for(size_t r=s;r<e;++r) got.push_back(decompress_id(dec, m, dprev, dprev_ptr, dprev_len));
                blockNames[b]=std::move(got);
            }
        };
        std::vector<std::thread> pool;
        for(unsigned t=0;t<NTHREADS;++t) pool.emplace_back(worker);
        for(auto& th:pool) th.join();
    }
    for(size_t b=0;b<nblocks && ok;++b){
        size_t s=b*BLOCK, e=std::min(n,s+BLOCK);
        for(size_t r=s;r<e;++r) if(blockNames[b][r-s]!=names[r]){ ok=false; fprintf(stderr,"MISMATCH block %zu read %zu\n",b,r); break; }
    }

    printf("names=%zu  block=%zu  threads=%u  coded=%zu B  %.3f B/name  round trip: %s\n",
           n,BLOCK,NTHREADS,bytes,bytes/(double)n, ok.load()?"VERIFIED":"FAILED");
    return ok.load()?0:1;
}
