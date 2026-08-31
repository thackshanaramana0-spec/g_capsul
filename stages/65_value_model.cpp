// STAGE 62 -- read-ID coder, SPRING's real algorithm (id_compression.cpp,
// shubhamchandak94/SPRING, itself crediting Samcomp/Mahoney & Bonfield 2012),
// cloned and read line by line, not from a summary. Genuinely different from
// stage 61's fqzcomp approach, and it is why stage 61 lost to SPRING/Genozip
// by 34-43% on real data (measured, not assumed).
//
// The key difference from fqzcomp: correspondence is tracked PER TOKEN
// INDEX, not per byte position. Walk the ID string classifying runs on the
// fly (alphabetic / zeros / digit / single-char, no delimiter list needed);
// for token k, compare against token k of the PREVIOUS id (via
// prev_tokens_ptr[k], a byte offset into the previous id string). A digit-
// width shift in an EARLIER token (e.g. read counter going from 6 digits to
// 7) never misaligns anything AFTER it, because later tokens are found by
// token INDEX, not by counting bytes forward from a shifted position -- the
// exact problem stage 61's byte-position realignment could only partially
// fix.
//
// Per token, in order: exact match with the same-index previous token
// (near-free once learned) > small positive delta for digit tokens (the
// read counter is almost always +1, extremely cheap once the model learns
// that) > raw value. Every field (token_type, alpha_len, digit bytes, delta,
// literal chars, zero-run length) has its OWN adaptive model PER TOKEN INDEX
// -- 1024 slots, matching SPRING's own MAX_NUMBER_TOKENS_ID.
//
//   g++ -O3 -march=native -o idcoder 62_id_tokenizer.cpp
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
// adaptive model over an N-symbol alphabet
struct Model {
    std::vector<uint32_t> f; uint32_t N;
    Model(uint32_t n=256):f(n,1),N(n){}
    // STAGE 64: real bit-cost estimate under CURRENT model state, without
    // encoding or mutating anything -- the same RDO principle already
    // validated earlier this session for read mapping (estimate real cost
    // from the live adaptive model, pick the cheapest candidate, THEN
    // encode). This is what lets the decision generalize safely: a candidate
    // only wins if it is ACTUALLY cheaper right now, not because it matches
    // a hand-picked rule that happens to fire on noise (X's contamination
    // of the blanket "any signed delta in range" trigger, found in the
    // first general-fix attempt).
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
// STAGE 63: SPRING's own ID_DELTA only fires for a SMALL POSITIVE delta
// (0 < delta < 256) -- verified directly against real data (yeast/human
// Illumina headers): the Y tile coordinate's delta has stdev ~348, mean
// ~0, and is 0.0% negative in this data (a real physical raster-scan
// signature) -- SPRING's narrow window catches only 28.9% of Y deltas,
// falling back to 4 raw bytes (ID_DIGIT) for the rest. X's delta, by
// contrast, has stdev ~18,710 (near the full coordinate range) -- no real
// structure there, confirmed before building anything for it.
// ID_ZDELTA extends the SAME underlying idea (delta from the same-index
// token in the previous read) to the FULL SIGNED range via zigzag + a
// 2-byte adaptive model, instead of a single-byte model restricted to
// (0,255) positive-only. This is not a new mechanism, it is SPRING's own
// mechanism minus the artificial window that was throwing away most of
// the real signal.
static const uint32_t MAXTOK=1024;

// STAGE 65. Measured directly on real data: for a field with NO cross-read
// structure (X coordinate -- 0% match rate, delta entropy 12.376 bits ==
// its own raw-value entropy, confirmed no tile-conditioning helps either),
// the ONLY remaining lever is how efficiently the raw VALUE gets coded.
// SPRING's ID_DIGIT spends 4 independent byte models (up to 32 bits) where
// the true entropy floor is ~12.5 bits (X has 6,328 distinct values in 1M
// reads) -- splitting into independent bytes throws away the fact that
// there are far fewer real DISTINCT VALUES than the raw bit-width implies.
//
// ARCS's own name_num_codec.h solves exactly this for numeric ID columns
// with a static order-0 dictionary -- proven, held-out verified -- but it
// is TWO-PASS (needs the whole column's distinct values before encoding
// one symbol), which does not fit this coder's one-read-at-a-time
// streaming design. Instead: reuse THIS PROJECT'S OWN already-proven
// pattern (stage 58/60's sparse per-key adaptive table, there keyed by pg
// position, here keyed by observed VALUE) -- online, growing, streaming-
// compatible, no two-pass requirement. A classic PPM-style ESCAPE symbol
// (Cleary & Witten) handles values not yet seen: pay the raw 4-byte cost
// only the FIRST time a distinct value appears (6,328 times for X, not
// 1,000,000), then near-entropy-floor cost every time after.
struct ValueDict {
    std::unordered_map<uint32_t,uint32_t> val2sym;   // value -> symbol (1-based; 0 = ESCAPE)
    std::vector<uint32_t> sym2val;
    std::vector<uint32_t> freq;                       // freq[0] = ESCAPE's own frequency
    uint32_t total;
    ValueDict():freq(1,1),total(1){}
    // returns symbol if known, else 0 (ESCAPE) -- caller decides what to do
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
    std::vector<Model> integer; // 4 bytes per token slot: used ONLY for the escape
                                 // (first-appearance) case now, not every occurrence
    std::vector<Model> zdelta_hi, zdelta_lo; // stage 63: wide signed delta, 2 adaptive bytes
    std::vector<ValueDict> valdict; // stage 65: per-token-index adaptive value dictionary
    // STAGE 64 v2: per-token-index track record, NOT model-cost comparison
    // (that had a real chicken-and-egg problem -- whichever model got used
    // first specialized and looked cheap regardless of true structure).
    // hit[k]/seen[k]: of the times token k had a real delta computed, how
    // often did |delta| stay small (<2048, well inside zigzag's efficient
    // range)? This observes the DATA directly, not a co-evolving model
    // cost, so it can't get stuck -- a position with real structure (Y)
    // accumulates a high hit rate fast; a noisy position (X) never does,
    // purely from its own history, no hardcoded index needed.
    std::vector<uint32_t> hit, seen;
    IdModels():
        token_type(MAXTOK,Model(8)), alpha_len(MAXTOK,Model(256)),
        alpha_value(MAXTOK,Model(128)), chars(MAXTOK,Model(128)),
        zero_run(MAXTOK,Model(256)), delta(MAXTOK,Model(256)),
        integer(MAXTOK*4,Model(256)),
        zdelta_hi(MAXTOK,Model(256)), zdelta_lo(MAXTOK,Model(256)),
        hit(MAXTOK,0), seen(MAXTOK,0), valdict(MAXTOK) {}
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
        token_len=0;   // BUG FIX: must be reset before use below, the original
                       // C resets it in its end-of-token trailer; this port
                       // dropped that reset, so the first comparison of every
                       // token after the first used the PREVIOUS token's
                       // leftover token_len as a stale index offset.
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
                // STAGE 64 v2: per-token-index TRACK RECORD, not a live
                // model-cost race (that had a real chicken-and-egg problem,
                // see the comment on IdModels::hit/seen). Update the
                // record from every real delta observed, regardless of
                // which path gets used, then decide from the LEARNED rate.
                if(prev_is_digit){
                    ++m.seen[token_ctr];
                    if(delta>=-2048 && delta<=2048) ++m.hit[token_ctr];
                }
                const uint32_t seen=m.seen[token_ctr], hit=m.hit[token_ctr];
                // Bootstrap period: not enough evidence yet, use SPRING's
                // own safe narrow-only default. Once enough evidence exists,
                // trust the position's own observed hit rate.
                const bool trust_wide = (seen>=20) && (hit*10 >= seen*3); // >=30%
                if(can_delta){
                    m.token_type[token_ctr].enc(enc,ID_DELTA);
                    m.delta[token_ctr].enc(enc,(uint32_t)delta);
                } else if(can_zdelta && trust_wide){
                    uint32_t z=(uint32_t)((delta<<1) ^ (delta>>63));
                    m.token_type[token_ctr].enc(enc,ID_ZDELTA);
                    m.zdelta_hi[token_ctr].enc(enc,(z>>8)&0xff);
                    m.zdelta_lo[token_ctr].enc(enc,z&0xff);
                } else {
                    // STAGE 65: route the "no delta structure" fallback
                    // through the adaptive value dictionary instead of raw
                    // bytes every time. token_type still signals ID_DIGIT
                    // to the decoder (the field TYPE is unchanged); what
                    // changes is HOW the value itself gets coded once we're
                    // in that branch.
                    m.token_type[token_ctr].enc(enc,ID_DIGIT);
                    ValueDict& vd = m.valdict[token_ctr];
                    uint32_t sym = vd.lookup(digit_value);
                    vd.encSym(enc,sym);
                    if(sym==0){   // ESCAPE: never seen this value at this token before
                        m.integer[token_ctr*4+0].enc(enc,(digit_value>>0)&0xff);
                        m.integer[token_ctr*4+1].enc(enc,(digit_value>>8)&0xff);
                        m.integer[token_ctr*4+2].enc(enc,(digit_value>>16)&0xff);
                        m.integer[token_ctr*4+3].enc(enc,(digit_value>>24)&0xff);
                        vd.registerNew(digit_value);
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
        if(token_ctr>=MAXTOK-1) break;   // safety
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
            uint32_t sym = vd.decSym(dec);
            uint32_t v;
            if(sym==0){   // ESCAPE: decode the raw value, then register it
                v = 0;
                v |= m.integer[token_ctr*4+0].dec(dec)<<0;
                v |= m.integer[token_ctr*4+1].dec(dec)<<8;
                v |= m.integer[token_ctr*4+2].dec(dec)<<16;
                v |= m.integer[token_ctr*4+3].dec(dec)<<24;
                vd.registerNew(v);
            } else {
                v = vd.sym2val[sym-1];
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
    if(argc<2){ fprintf(stderr,"usage: %s names.txt\n",argv[0]); return 1; }
    std::vector<std::string> names;
    { FILE* f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
      char buf[2048];
      while(fgets(buf,sizeof(buf),f)){ size_t L=strlen(buf); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
          names.emplace_back(buf,L); }
      fclose(f); }
    const size_t n=names.size();
    fprintf(stderr,"names=%zu\n",n);

    IdModels enc_models;
    RangeEnc enc; enc.out.reserve(n*4);
    std::string prev_id; std::array<uint32_t,MAXTOK> prev_ptr{}; prev_ptr.fill(0);
    for(size_t r=0;r<n;++r) compress_id(enc, enc_models, names[r].c_str(), prev_id, prev_ptr);
    enc.flush();
    const size_t bytes=enc.out.size();

    IdModels dec_models;
    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    std::string dprev; std::array<uint32_t,MAXTOK> dprev_ptr{}, dprev_len{};
    dprev_ptr.fill(0); dprev_len.fill(0);
    bool ok=true;
    for(size_t r=0;r<n;++r){
        std::string got = decompress_id(dec, dec_models, dprev, dprev_ptr, dprev_len);
        if(got!=names[r]){ ok=false; fprintf(stderr,"MISMATCH at %zu:\n got : %s\n want: %s\n",r,got.c_str(),names[r].c_str()); break; }
    }
    printf("names=%zu  coded=%zu B  %.3f B/name  round trip: %s\n",
           n,bytes,bytes/(double)n, ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
