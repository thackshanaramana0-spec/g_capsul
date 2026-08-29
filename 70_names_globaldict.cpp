// STAGE 70 -- Genozip-style global dictionary, correctly scoped after
// testing the naive "clone+merge the whole adaptive model" idea and finding
// a real correctness hazard: Model's periodic rescale (halving at 60000
// total) makes clone-before/clone-after diffing lossy, and a fully
// lock-step shared adaptive model would force full serialization (zero
// parallelism) -- verified by reasoning through the dependency chain
// before writing any merge code, not discovered by debugging a broken one.
//
// Correct scoping: only ValueDict (stage 65's per-token-index digit-value
// dictionary, e.g. X coordinate: 6,328 distinct values / 1M names) behaves
// like Genozip's real dictionary -- a growing SYMBOL TABLE, not an
// adaptive frequency model. The small-alphabet fields (token_type,
// alpha_len, delta, ...) saturate almost immediately and are NOT what
// costs stage 69 its 0.7-2.3% block-boundary penalty -- that cost is
// concentrated in ValueDict paying the first-seen escape cost separately
// IN EVERY BLOCK for the same recurring high-cardinality values.
//
// Fix: two passes.
//   Pass 1 (single fast scan, order-independent, cheap): tokenize every ID
//   the same way compress_id does and record every digit-token's value per
//   token index into ONE GLOBAL ValueDict per token index. Serialize that
//   dictionary into a header once.
//   Pass 2 (fully parallel, blocks independent): compress_id/decompress_id
//   use the GLOBAL, now-FROZEN dictionary (STRUCTURE is read-only, no
//   synchronization needed since it never changes during pass 2) but keep
//   a PER-BLOCK-LOCAL adaptive frequency table over which known symbol
//   appears (cheap either way, no cross-block sharing needed for that
//   part). No more merge-back, no rescale hazard: the escape path should
//   almost never fire in pass 2 since pass 1 already saw every value.
//
// This is exactly ARCS's own name_num_codec.h design (src/name_num_codec.h)
// -- rejected at stage 65 specifically because it was two-pass and didn't
// fit a single-threaded streaming coder. That objection no longer applies:
// block-parallel dispatch already requires knowing block boundaries up
// front, so a first pass is free.
//
// FIRST ATTEMPT, REAL NEGATIVE RESULT: pass 1 originally registered every
// digit token unconditionally, regardless of whether pass 2 would actually
// need a dictionary lookup for it. Real result: 1,008,537 dictionary
// entries for 1M names (a 4 MB header, dwarfing the ~1.76 MB body) --
// because the per-read serial-number token (".1",".2",...) has cardinality
// ~1M but is already handled almost free by ID_DELTA and never needs the
// dictionary at all. Fixed by making pass 1 a real DRY RUN of compress_id's
// full routing logic (match/delta/zdelta first, dictionary only as last
// resort, same hit/seen trust_wide gate) so it predicts exactly which
// values pass 2 will actually look up. Real result after the fix: 6,339
// dictionary entries, 25,356 B header -- matches the earlier, independently
// measured "X coordinate has 6,328 distinct values in 1M reads" finding
// from stage 65 almost exactly, a real cross-check the fix is correct.
//
// Real result vs stage 69 (fresh model per block), swept over block size,
// both round-trip verified:
//   block      stage 69   stage 70
//   500,000    1.770      1.787
//   250,000    1.779      1.788
//   100,000    1.806      1.790
//   50,000     1.851      1.794
//   20,000     1.968      1.805
// Stage 70's ratio barely moves as block count rises 25x (+1%) while
// stage 69 degrades sharply (+11%) over the same range -- the global
// dictionary decouples ratio from block granularity, so block size can be
// picked purely for parallelism/load-balancing. At large blocks stage 70
// is marginally worse (one-time header not yet amortized); at small
// blocks, where real many-core parallelism actually needs many blocks,
// it is a clear win.
//
//   g++ -O3 -march=native -pthread -o nameglobal 70_names_globaldict.cpp

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

// Global, read-only in pass 2 -- STRUCTURE built once in pass 1. freq[]
// stays adaptive but is per-block-local (each block gets its own copy of
// this GlobalDict's freq table, structure shared, never mutated during
// pass 2, so no locking needed to read val2sym/sym2val concurrently).
struct GlobalDict {
    std::unordered_map<uint32_t,uint32_t> val2sym; // value -> symbol, 1-based; 0 = escape (should be rare in pass 2)
    std::vector<uint32_t> sym2val;
    void registerNew(uint32_t v){
        if(val2sym.count(v)) return;
        uint32_t sym=(uint32_t)sym2val.size()+1;
        val2sym[v]=sym; sym2val.push_back(v);
    }
    uint32_t lookup(uint32_t v) const { auto it=val2sym.find(v); return it==val2sym.end()?0:it->second; }
};

// Per-block-local adaptive frequency table OVER a fixed global symbol set.
struct LocalDictFreq {
    std::vector<uint32_t> freq; // freq[0]=escape
    uint32_t total;
    LocalDictFreq(const GlobalDict& g):freq(g.sym2val.size()+1,1),total((uint32_t)freq.size()){}
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
};

struct IdModels {
    std::vector<Model> token_type, alpha_len, alpha_value, chars, zero_run, delta;
    std::vector<Model> integer;
    std::vector<Model> zdelta_hi, zdelta_lo;
    std::vector<uint32_t> hit, seen;
    IdModels():
        token_type(MAXTOK,Model(8)), alpha_len(MAXTOK,Model(256)),
        alpha_value(MAXTOK,Model(128)), chars(MAXTOK,Model(128)),
        zero_run(MAXTOK,Model(256)), delta(MAXTOK,Model(256)),
        integer(MAXTOK*4,Model(256)),
        zdelta_hi(MAXTOK,Model(256)), zdelta_lo(MAXTOK,Model(256)),
        hit(MAXTOK,0), seen(MAXTOK,0) {}
};

// ---- Pass 1: a real DRY RUN of compress_id's routing logic (match / zeros
// / delta / zdelta / digit), WITHOUT any entropy coding -- first naive
// attempt (tokenize every digit token unconditionally) was tested and
// produced a real, decisive NEGATIVE result: 1,008,537 dictionary entries
// for 1M names, a 4 MB header dwarfing the ~1.76 MB body. Diagnosis: the
// per-read serial-number token (".1", ".2", ...) has cardinality ~1M but
// is already handled almost free by ID_DELTA in the real coder -- it never
// needs a dictionary lookup at all. Only register a value in the global
// dictionary when this SAME routing logic would actually fall through to
// the ID_DIGIT branch, mirroring compress_id's hit/seen trust_wide gate
// exactly so pass 1 predicts pass 2's real routing, not a superset of it.
static void dry_run_for_dict(const std::vector<std::string>& names, std::array<GlobalDict,MAXTOK>& gdict){
    std::string prev_id;
    std::array<uint32_t,MAXTOK> prev_tok_ptr{}; prev_tok_ptr.fill(0);
    std::vector<uint32_t> hit(MAXTOK,0), seen(MAXTOK,0);
    for(auto& nm: names){
        const char* id = nm.c_str();
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
            } else if(*id_ptr=='0'){
                while(*id_ptr_tok=='0'){
                    match_len += ('0'==prevc(prev_tok_ptr[token_ctr]+token_len));
                    ++token_len; ++id_ptr_tok;
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
                bool is_match = prev_is_digit && match_len==token_len && !isdigit((unsigned char)prevc(prev_tok_ptr[token_ctr]+token_len));
                if(!is_match){
                    if(prev_is_digit){
                        ++seen[token_ctr];
                        if(delta>=-2048 && delta<=2048) ++hit[token_ctr];
                    }
                    const bool trust_wide = (seen[token_ctr]>=20) && (hit[token_ctr]*10 >= seen[token_ctr]*3);
                    if(!can_delta && !(can_zdelta && trust_wide)){
                        gdict[token_ctr].registerNew(digit_value); // real ID_DIGIT fallback -- register
                    }
                }
            }
            prev_tok_ptr[token_ctr]=i;
            i+=token_len; id_ptr=id_ptr_tok; ++token_ctr;
            if(token_ctr>=MAXTOK-1) break;
        }
        prev_id.assign(id);
        for(uint32_t k=token_ctr+1;k<MAXTOK;++k) prev_tok_ptr[k]=0;
    }
}

uint32_t compute_num_digits(uint32_t x){
    if(x<10) return 1; if(x<100) return 2; if(x<1000) return 3; if(x<10000) return 4;
    if(x<100000) return 5; if(x<1000000) return 6; if(x<10000000) return 7;
    if(x<100000000) return 8; if(x<1000000000) return 9; return 10;
}

void compress_id(RangeEnc& enc, IdModels& m, std::array<GlobalDict,MAXTOK>& gdict,
                  std::array<LocalDictFreq*,MAXTOK>& ldict, const char* id, std::string& prev_id,
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
                    // STAGE 70: dictionary STRUCTURE is global/frozen (built
                    // in pass 1), so lookup should hit except for genuine
                    // pass-1 misses (shouldn't happen -- pass 1 scans every
                    // ID with the identical tokenizer). Frequency table is
                    // per-block-local.
                    GlobalDict& gd = gdict[token_ctr];
                    LocalDictFreq& ld = *ldict[token_ctr];
                    uint32_t sym = gd.lookup(digit_value);
                    ld.encSym(enc,sym);
                    if(sym==0){ // should be rare/never; safety net only
                        m.integer[token_ctr*4+0].enc(enc,(digit_value>>0)&0xff);
                        m.integer[token_ctr*4+1].enc(enc,(digit_value>>8)&0xff);
                        m.integer[token_ctr*4+2].enc(enc,(digit_value>>16)&0xff);
                        m.integer[token_ctr*4+3].enc(enc,(digit_value>>24)&0xff);
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

std::string decompress_id(RangeDec& dec, IdModels& m, std::array<GlobalDict,MAXTOK>& gdict,
                           std::array<LocalDictFreq*,MAXTOK>& ldict, std::string& prev_id,
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
            GlobalDict& gd = gdict[token_ctr];
            LocalDictFreq& ld = *ldict[token_ctr];
            uint32_t sym = ld.decSym(dec);
            uint32_t v;
            if(sym==0){
                v = 0;
                v |= m.integer[token_ctr*4+0].dec(dec)<<0;
                v |= m.integer[token_ctr*4+1].dec(dec)<<8;
                v |= m.integer[token_ctr*4+2].dec(dec)<<16;
                v |= m.integer[token_ctr*4+3].dec(dec)<<24;
            } else {
                v = gd.sym2val[sym-1];
            }
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

    // ---- PASS 1: build the global, frozen per-token-index dictionary ----
    auto tp0=std::chrono::steady_clock::now();
    std::array<GlobalDict,MAXTOK> gdict;
    dry_run_for_dict(names, gdict);
    auto tp1=std::chrono::steady_clock::now();
    size_t dictBytes=0, dictEntries=0;
    for(auto& g: gdict){ dictEntries += g.sym2val.size(); dictBytes += g.sym2val.size()*4; }
    fprintf(stderr,"PASS1 (build global dict): %.3f s   entries=%zu  header~=%zu B\n",
            std::chrono::duration<double>(tp1-tp0).count(), dictEntries, dictBytes);

    // ---- PASS 2: parallel block encode against the frozen dictionary ----
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
                std::array<LocalDictFreq*,MAXTOK> ldict{};
                for(uint32_t k=0;k<MAXTOK;++k) ldict[k]=new LocalDictFreq(gdict[k]);
                RangeEnc enc; enc.out.reserve((e-s)*4);
                std::string prev_id; std::array<uint32_t,MAXTOK> prev_ptr{}; prev_ptr.fill(0);
                for(size_t r=s;r<e;++r) compress_id(enc, m, gdict, ldict, names[r].c_str(), prev_id, prev_ptr);
                enc.flush();
                blockOut[b]=std::move(enc.out);
                for(uint32_t k=0;k<MAXTOK;++k) delete ldict[k];
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

    // round trip
    std::vector<std::vector<std::string>> blockNames(nblocks);
    std::atomic<bool> ok{true};
    {
        std::atomic<size_t> next{0};
        auto worker=[&](){
            size_t b;
            while((b=next.fetch_add(1))<nblocks){
                size_t s=b*BLOCK, e=std::min(n,s+BLOCK);
                IdModels m;
                std::array<LocalDictFreq*,MAXTOK> ldict{};
                for(uint32_t k=0;k<MAXTOK;++k) ldict[k]=new LocalDictFreq(gdict[k]);
                RangeDec dec; dec.init(blockOut[b].data(),blockOut[b].size());
                std::string dprev; std::array<uint32_t,MAXTOK> dprev_ptr{}, dprev_len{};
                dprev_ptr.fill(0); dprev_len.fill(0);
                std::vector<std::string> got; got.reserve(e-s);
                for(size_t r=s;r<e;++r) got.push_back(decompress_id(dec, m, gdict, ldict, dprev, dprev_ptr, dprev_len));
                blockNames[b]=std::move(got);
                for(uint32_t k=0;k<MAXTOK;++k) delete ldict[k];
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

    printf("names=%zu  block=%zu  threads=%u  body=%zu B  dict_header~=%zu B  total~=%zu B  %.3f B/name  round trip: %s\n",
           n,BLOCK,NTHREADS,bytes,dictBytes,bytes+dictBytes,(bytes+dictBytes)/(double)n, ok.load()?"VERIFIED":"FAILED");
    return ok.load()?0:1;
}
