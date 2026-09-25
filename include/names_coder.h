#pragma once
// Names / read-ID coder, shared by the encoder (106_inprocess.cpp) and the
// archive decoder (capsule_decode.cpp) so the two cannot diverge -- the same
// reason seqpar_core.h exists.
//
// Algorithm: stage 86 (= stage 70's global frozen dictionary + stage 85's
// bounded-memory streaming queue, merged). Namespaced `nmc` following this
// project's existing per-coder convention (mmc, refc, pgc), each of which
// carries its own RangeEnc/RangeDec.
//
// FOUR THINGS THIS FIXES relative to stage 86, all of which a mechanical
// copy of that file would have carried into the archive:
//
// 1. Stage 86 streams from a names.txt that someone extracted by hand. Here
//    pass 1 and pass 2 both stream the ID column straight out of the FASTQ,
//    so nothing accumulates. Holding the headers instead (the obvious
//    MemStream route, matching every other stream in this pipeline) would be
//    ~5 GB on T. cacao and would defeat the entire point of merging stage 85
//    in -- the bounded queue would sit behind a fully-materialized input.
//
// 2. Stage 86's "dict_header~=25356 B" is `entries*4`, arithmetic, never
//    serialized. Here the dictionary is really serialized, and returned RAW
//    for the caller to run through best_encode() like every other stream --
//    so its cost is real bytes in the archive, and it gets the coder
//    selector instead of raw uint32s.
//
// 3. Stage 86 has NO block index: its round trip re-reads the original file
//    to learn how many names are in each block. A decoder holding only the
//    archive cannot do that. Block byte-lengths and name-counts are a real
//    stream here and are really counted. This is the exact bug class that
//    has already cost this project twice (mm_pos/mm_count omitted from
//    totals; mem_triples' destination/length/rc dropped as "diagnostics
//    only" leaving an archive that could not decode) -- a cost landing in a
//    stream nobody measures.
//
// 4. Decode streams block-by-block to a FILE* rather than returning a
//    vector<string> of every name, so the decoder keeps the same bounded
//    memory property as the encoder.
//
// The '@' is not stored: it is structurally constant on every FASTQ header
// line and is re-added on output. SPRING and Genozip both do the same.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <vector>
#include <string>
#include <array>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <algorithm>
#include <unordered_map>
#include <cmath>

namespace nmc {

// ---- range coder ------------------------------------------------------------
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
    const uint8_t* p=nullptr; const uint8_t* end=nullptr; uint32_t range=0xFFFFFFFFu, code=0;
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

// ID_SEQLEN: a digit token whose value IS the read's own sequence length.
//
// Every SRA-derived FASTQ header here ends in "length=NNN", and that number is
// the length of the SEQ line two lines below it -- a value the CAPSULE archive
// ALREADY stores, per original read, in the `read_lengths` stream for sequence
// reconstruction. Coding it a second time inside the name is paying twice for
// one fact. Measured cost of that duplication before this token existed:
// ERR5181310 38,059 of 39,416 B (96.6% of the entire names stream) and
// ERR552797 269,005 of 1,341,038 B (20.1%). Genozip pays for the same field in
// its own separate `length` context (25,395 B and 177,459 B on the same two
// files); we can pay nothing, because we already hold the value.
//
// So the token carries NO payload -- just the type symbol, which an adaptive
// model drives to near-zero when it dominates a token index. The decoder
// recovers the digits from the read length it has already decoded.
//
// This is a cross-column reference, the same class of mechanism Genozip uses
// (SNIP_REDIRECTION between contexts), but strictly cheaper here because the
// sequence side needs read_lengths regardless of what names do.
enum TokType { ID_ALPHA, ID_DIGIT, ID_CHAR, ID_MATCH, ID_ZEROS, ID_DELTA, ID_END, ID_ZDELTA, ID_SEQLEN };
// Diagnostic only: isolate our two additions over SPRING's baseline coder.
static const bool NO_ZD   = getenv("NMC_NOZD")!=nullptr;    // drop ID_ZDELTA
static const bool NO_DICT = getenv("NMC_NODICT")!=nullptr;  // drop the value dictionary
static const uint32_t MAXTOK=1024;

struct GlobalDict {
    std::unordered_map<uint32_t,uint32_t> val2sym;   // value -> symbol, 1-based; 0 = escape
    std::vector<uint32_t> sym2val;
    void registerNew(uint32_t v){
        if(val2sym.count(v)) return;
        val2sym[v]=(uint32_t)sym2val.size()+1; sym2val.push_back(v);
    }
    uint32_t lookup(uint32_t v) const { auto it=val2sym.find(v); return it==val2sym.end()?0:it->second; }
};
// Adaptive frequency model over the value dictionary.
//
// TWO REAL DEFECTS FIXED HERE, both only visible on high-cardinality fields
// (the X/Y coordinates of raw Illumina headers), which is exactly where
// Genozip beats us:
//
// 1. The rescale ceiling was a hardcoded 60000 while `total` STARTS at N (all
//    frequencies initialised to 1). So any dictionary with more than ~60k
//    symbols was already over the limit before coding its first symbol, and
//    halved every frequency on every symbol thereafter -- the model could
//    never accumulate a statistic, and the dictionary degenerated to a flat
//    log2(N) code (~16.6 bits/value at N=98k) while still paying dictionary
//    overhead. Between ~3.7k and 60k symbols it was impaired rather than
//    dead. The range coder's ACTUAL constraint is tot <= 2^24 (it computes
//    range/=tot after normalising range >= 2^24), so 60000 was needlessly
//    conservative by two orders of magnitude. The ceiling now scales with the
//    alphabet, giving every dictionary a real adaptation window.
//
// 2. encSym/decSym computed the cumulative frequency by LINEAR SCAN, O(N) per
//    symbol. At N=98k over 500k names that is ~10^11 operations -- this is
//    what looked like a hang on raw-Illumina headers (it was quadratic
//    blow-up, not a deadlock), and in decSym the unbounded
//    `while(lo+freq[sym]<=target)` walk could also run past the array.
//    Replaced with a Fenwick tree (binary indexed tree): O(log N) prefix sum
//    to encode, O(log N) binary-lifting search to decode. Standard,
//    decades-proven for large-alphabet adaptive range coding.
struct LocalDictFreq {
    std::vector<uint32_t> freq;   // per-symbol frequency
    std::vector<uint32_t> bit;    // Fenwick tree, 1-based
    uint32_t N, total, limit, LOG;
    explicit LocalDictFreq(const GlobalDict& g)
        : freq(g.sym2val.size()+1,1), N((uint32_t)freq.size()), total(N) {
        limit = (uint32_t)std::min<uint64_t>(1u<<22,
                    std::max<uint64_t>(60000u, (uint64_t)N*8));
        LOG=1; while((LOG<<1) <= N) LOG<<=1;
        build();
    }
    void build(){
        bit.assign(N+1,0);
        for(uint32_t i=1;i<=N;++i){
            bit[i]+=freq[i-1];
            uint32_t j=i+(i&(uint32_t)(-(int32_t)i));
            if(j<=N) bit[j]+=bit[i];
        }
    }
    inline uint32_t prefix(uint32_t i) const {          // sum of freq[0..i-1]
        uint32_t s=0; while(i){ s+=bit[i]; i-=i&(uint32_t)(-(int32_t)i); } return s;
    }
    inline void addAt(uint32_t i,uint32_t v){           // freq[i] += v
        for(uint32_t x=i+1;x<=N;x+=x&(uint32_t)(-(int32_t)x)) bit[x]+=v;
    }
    inline void bump(uint32_t sym){
        freq[sym]+=16; total+=16; addAt(sym,16);
        if(total>limit){ total=0; for(auto& f:freq){ f=(f>>1)|1; total+=f; } build(); }
    }
    void encSym(RangeEnc& rc,uint32_t sym){
        uint32_t lo=prefix(sym);
        rc.encode(lo,lo+freq[sym],total);
        bump(sym);
    }
    uint32_t decSym(RangeDec& rc){
        uint32_t target=rc.getFreq(total);
        uint32_t idx=0, rem=target;
        for(uint32_t pw=LOG; pw; pw>>=1){
            if(idx+pw<=N && bit[idx+pw]<=rem){ idx+=pw; rem-=bit[idx]; }
        }
        uint32_t sym = idx<N ? idx : N-1;
        uint32_t lo = target-rem;                        // == prefix(sym)
        rc.decodeUpdate(lo,lo+freq[sym]);
        bump(sym);
        return sym;
    }
};
struct IdModels {
    std::vector<Model> token_type, alpha_len, alpha_value, chars, zero_run, delta;
    std::vector<Model> integer;
    std::vector<Model> zdelta_hi, zdelta_lo;
    std::vector<uint32_t> hit, seen;
    IdModels():
        token_type(MAXTOK,Model(9)), alpha_len(MAXTOK,Model(256)),
        alpha_value(MAXTOK,Model(128)), chars(MAXTOK,Model(128)),
        zero_run(MAXTOK,Model(256)), delta(MAXTOK,Model(256)),
        integer(MAXTOK*4,Model(256)),
        zdelta_hi(MAXTOK,Model(256)), zdelta_lo(MAXTOK,Model(256)),
        hit(MAXTOK,0), seen(MAXTOK,0) {}
};

// ---- varint helpers ---------------------------------------------------------
static inline void pv(std::vector<uint8_t>& o, uint64_t v){
    while(v>=0x80){ o.push_back((uint8_t)(v|0x80)); v>>=7; }
    o.push_back((uint8_t)v);
}
static inline uint64_t gv(const uint8_t*& p, const uint8_t* end){
    uint64_t v=0; int s=0;
    while(p<end){ uint8_t b=*p++; v |= (uint64_t)(b&0x7f)<<s; if(!(b&0x80)) break; s+=7; }
    return v;
}

// ---- file-constant token elision -------------------------------------------
// Measured on real SRA headers: 23 of 27 tokens on E. coli, 24 of 28 on
// DRR976266, are byte-identical in EVERY read of the file -- the accession,
// the instrument, run, flowcell, lane, every separator, and "length=NNN".
// They already cost nothing to represent (each resolves to ID_MATCH), but
// each still spends one token_type symbol per read, and 23 x 500,000 of
// those is ~43 KB on a 1.36 MB stream. Genozip stores such a field once:
// its constant QNAME contexts are 44-223 B for 500,000 reads.
//
// So: pass 1 already reads the whole file, and can record for free which
// token indices never vary. Those are written once into the header and
// carry ZERO per-read bits -- no token_type, no payload. Keyed on a
// measured property of the input (invariance), not on a dataset, so a file
// with no constant field simply gets an empty table and a byte-identical
// stream to before.
//
// Only engaged when every read has the SAME token count. That is what makes
// the token walk terminable without an ID_END marker, and it is checked, not
// assumed -- variable-count files (mixed read groups) keep the old path
// exactly.
struct Layout {
    uint32_t ntok=0;                     // highest token count seen
    std::vector<uint8_t> isconst;        // ntok flags
    std::vector<std::string> constant;   // ntok texts (valid where isconst)
    // A token index where EVERY read's value is that read's own sequence
    // length. Hoisted into the header exactly like a constant index: the
    // per-read cost falls from one ID_SEQLEN type symbol to nothing at all.
    std::vector<uint8_t> isseqlen;
    // Line 3 of a FASTQ record. In every SRA-derived file here it is '+'
    // followed by a verbatim repeat of the header -- 15.5% of the raw text,
    // and pure redundancy. Stored as ONE byte for the whole file:
    //   0 = "+" alone, 1 = "+" + the header, 2 = neither (must be stored, not
    //   supported yet -- an archive claiming losslessness must not see this).
    // Same admission rule as every other elision: verified on EVERY read in
    // pass 1, never inferred from the first one.
    uint8_t line3_mode=2;
    bool active() const { return ntok>0 && !isconst.empty(); }
};

// ---- pass 1: dry run predicting compress_id's ID_DIGIT fallback -------------
// Streams the FASTQ ID column; one header materialized at a time.
static void dict_pass(const char* fq_path, std::array<GlobalDict,MAXTOK>& gdict,
                      uint64_t* n_names_out, Layout* lay=nullptr){
    FILE* f=fopen(fq_path,"r");
    if(!f){
        // Not reachable today: encoder.cpp's own pre-flight std::ifstream check
        // on argv[1] already refuses before this runs. Kept as a loud diagnostic
        // rather than a silent no-op anyway, since encode_from_fastq() would
        // otherwise return a normally-shaped Encoded{n_names=0} here with
        // nothing to tell a caller "the file couldn't be opened" apart from
        // "there happened to be zero names" -- the same class of silent-success
        // failure this project's own encoder.cpp SCOPE CHECK block exists to
        // eliminate for the oversize-read/malformed-header cases.
        fprintf(stderr,"WARNING: names_coder::dict_pass could not open '%s' -- "
                       "names stream will be silently empty\n", fq_path);
        return;
    }
    std::string prev_id;
    std::array<uint32_t,MAXTOK> prev_tok_ptr{}; prev_tok_ptr.fill(0);
    std::vector<uint32_t> hit(MAXTOK,0), seen(MAXTOK,0);
    std::vector<char> buf(1u<<16);
    uint64_t lineno=0, nn=0;
    // constancy tracking: first read seeds each index, later reads clear the
    // flag on any disagreement; a differing token COUNT disables the whole
    // mechanism (see Layout).
    std::vector<std::string> ctext(MAXTOK);
    std::vector<uint8_t>     cflag(MAXTOK,0);
    std::vector<uint64_t>    cseen(MAXTOK,0);   // reads that HAVE this index
    std::vector<uint8_t>     sflag(MAXTOK,1);   // every occurrence == seq length
    std::vector<uint64_t>    sseen(MAXTOK,0);
    uint32_t ntok_max=0;
    // Per-token-index statistics for the dictionary decision. Counted only on
    // tokens that actually reach the ID_DIGIT fallback -- the same routing
    // pass 2 will take -- so the estimate prices the real alternative.
    std::vector<std::unordered_map<uint32_t,uint32_t>> vcnt(MAXTOK);
    std::vector<std::array<std::array<uint32_t,256>,4>> bhist(MAXTOK);
    for(auto& h:bhist) for(auto& q:h) q.fill(0);
    std::vector<uint64_t> vtot(MAXTOK,0);
    // One-line lookahead: a header's "length=NNN" is priced against the SEQ
    // line that FOLLOWS it, so the header is held until that length is known.
    // Without this, pass 1 would register the length values in the dictionary
    // and the gate would price a path pass 2 no longer takes (ID_SEQLEN).
    std::string pend; bool havepend=false; uint32_t pend_slen=0;
    std::string lasthdr; bool l3_all_plus=true, l3_all_copy=true; uint64_t l3_seen=0;
    while(fgets(buf.data(),(int)buf.size(),f) || havepend){
        const bool eof_flush = havepend && feof(f);
        if(!eof_flush){
            const uint64_t ln = lineno++;
            if(ln%4==0){
                size_t L0=strlen(buf.data()); while(L0&&(buf[L0-1]=='\n'||buf[L0-1]=='\r')) --L0;
                buf[L0]=0; pend.assign(buf.data(),L0); havepend=true; pend_slen=0;
                continue;
            }
            if(ln%4==1 && havepend){
                size_t L1=strlen(buf.data()); while(L1&&(buf[L1-1]=='\n'||buf[L1-1]=='\r')) --L1;
                pend_slen=(uint32_t)L1;
            } else if(ln%4==2){
                size_t L3=strlen(buf.data()); while(L3&&(buf[L3-1]=='\n'||buf[L3-1]=='\r')) --L3;
                buf[L3]=0; ++l3_seen;
                const char* h = lasthdr.c_str();
                const size_t hn = lasthdr.size();
                if(!(L3==1 && buf[0]=='+')) l3_all_plus=false;
                if(!(L3==hn+1 && buf[0]=='+' && memcmp(buf.data()+1,h,hn)==0)) l3_all_copy=false;
                continue;
            } else continue;
        }
        havepend=false;
        const uint32_t seqlen = pend_slen;
        ++nn;
        // header WITHOUT the '@' -- line 3 is '+' plus that, not '+@...'
        lasthdr = (!pend.empty() && pend[0]=='@') ? pend.substr(1) : pend;
        const char* id = pend.c_str();
        if(*id=='@'){ ++id; }               // '@' is structural, not stored
        const char* id_ptr=id;
        uint32_t token_len=0, match_len=0, token_ctr=0, i=0;
        const char* prevbuf = prev_id.c_str();
        const size_t prevlen = prev_id.size();
        auto prevc=[&](uint32_t off)->char{ return off<prevlen ? prevbuf[off] : 0; };
        while(*id_ptr){
            token_len=0; bool tok_is_seqlen=false;
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
                const bool is_match = prev_is_digit && match_len==token_len &&
                                      !isdigit((unsigned char)prevc(prev_tok_ptr[token_ctr]+token_len));
                // Independent of the coding branch: a repeated length would be
                // coded ID_MATCH, which previously hid the property entirely
                // and stopped the hoist from ever firing.
                tok_is_seqlen = seqlen && digit_value==seqlen
                                && !isdigit((unsigned char)*id_ptr_tok);
                if(!is_match){
                    if(prev_is_digit){
                        ++seen[token_ctr];
                        if(delta>=-2048 && delta<=2048) ++hit[token_ctr];
                    }
                    const bool trust_wide = !NO_ZD && (seen[token_ctr]>=20) && (hit[token_ctr]*10 >= seen[token_ctr]*3);
                    const bool is_seqlen = tok_is_seqlen;
                    if(!can_delta && !(can_zdelta && trust_wide) && !NO_DICT && !is_seqlen){
                        ++vcnt[token_ctr][digit_value];
                        ++vtot[token_ctr];
                        for(int bb=0;bb<4;++bb) ++bhist[token_ctr][bb][(digit_value>>(bb*8))&0xff];
                    }
                }
            }
            if(lay){
                if(cseen[token_ctr]==0){ ctext[token_ctr].assign(id_ptr,token_len); cflag[token_ctr]=1; }
                else if(cflag[token_ctr]){
                    if(ctext[token_ctr].size()!=token_len ||
                       memcmp(ctext[token_ctr].data(),id_ptr,token_len)!=0) cflag[token_ctr]=0;
                }
                ++cseen[token_ctr];
                if(tok_is_seqlen) ++sseen[token_ctr]; else sflag[token_ctr]=0;
            }
            prev_tok_ptr[token_ctr]=i;
            i+=token_len; id_ptr=id_ptr_tok; ++token_ctr;
            if(token_ctr>=MAXTOK-1) break;
        }
        if(lay && token_ctr>ntok_max) ntok_max=token_ctr;
        prev_id.assign(id);
        for(uint32_t k=token_ctr+1;k<MAXTOK;++k) prev_tok_ptr[k]=0;
    }
    fclose(f);
    // ---- dictionary decision, per token index, by measured cost ------------
    // Stage 66 gated its LOCAL dictionary on an observed hit rate; that gate
    // was lost when stage 70's GLOBAL dictionary replaced it, and stage 70
    // was only ever measured on one file. Measured across 14 real datasets
    // the unconditional dictionary is a large LOSS on high-cardinality name
    // fields (+351,569 B on HG002, +125,846 B on ERR552797) and a large win
    // on low-cardinality ones (-113,766 B on DRR976266).
    //
    // So price both routes from pass 1's own histograms and keep the cheaper:
    //   dictionary : N*H(values) + header (4 B/entry, conservative)
    //   raw        : N*sum of the four byte-model order-0 entropies
    // An index that loses simply gets no entries; an empty GlobalDict makes
    // LocalDictFreq a one-symbol alphabet, which the range coder emits in
    // zero bits, so the fallback needs no flag and no format change.
    for(uint32_t k=0;k<MAXTOK;++k){
        const uint64_t N=vtot[k];
        if(!N) continue;
        double Hv=0.0;
        for(auto& e:vcnt[k]){ const double p=(double)e.second/(double)N; Hv-=p*std::log2(p); }
        const double dict_bits = (double)N*Hv + (double)vcnt[k].size()*4.0*8.0;
        double Hb=0.0;
        for(int bb=0;bb<4;++bb){
            double h=0.0;
            for(int v=0;v<256;++v){
                const uint32_t c=bhist[k][bb][v]; if(!c) continue;
                const double p=(double)c/(double)N; h-=p*std::log2(p);
            }
            Hb+=h;
        }
        const double raw_bits = (double)N*Hb;
        if(dict_bits < raw_bits)
            for(auto& e:vcnt[k]) gdict[k].registerNew(e.first);
    }
    if(lay){
        lay->line3_mode = (l3_seen && l3_all_plus) ? 0
                        : (l3_seen && l3_all_copy) ? 1 : 2;
    }
    if(lay && ntok_max>0){
        lay->ntok = ntok_max;
        lay->isconst.assign(ntok_max,0);
        lay->constant.assign(ntok_max,std::string());
        lay->isseqlen.assign(ntok_max,0);
        for(uint32_t k=0;k<ntok_max;++k){
            if(cflag[k] && cseen[k]==nn){          // invariant AND present in every read
                lay->isconst[k]=1; lay->constant[k]=ctext[k];
            }
            // Same admission rule as a constant index: the property must hold
            // in EVERY read, so the token walk stays terminable. A constant
            // index wins if both apply (it needs no read_lengths lookup).
            else if(sflag[k] && sseen[k]==nn && nn>0) lay->isseqlen[k]=1;
        }
    }
    if(n_names_out) *n_names_out = nn;
}

// ---- per-name encode / decode ----------------------------------------------
static void compress_id(RangeEnc& enc, IdModels& m, std::array<GlobalDict,MAXTOK>& gdict,
                        std::array<LocalDictFreq*,MAXTOK>& ldict, const char* id,
                        std::string& prev_id, std::array<uint32_t,MAXTOK>& prev_tok_ptr,
                        const Layout* lay=nullptr, uint32_t seqlen=0){
    const bool ELIDE = lay && lay->active();
    uint32_t token_len=0, match_len=0, token_ctr=0, i=0;
    const char* id_ptr=id;
    const char* prevbuf = prev_id.c_str();
    const size_t prevlen = prev_id.size();
    auto prevc=[&](uint32_t off)->char{ return off<prevlen ? prevbuf[off] : 0; };
    while(*id_ptr){
        token_len=0;
        match_len = (*id_ptr==prevc(prev_tok_ptr[token_ctr]+token_len)); token_len=1;
        const char* id_ptr_tok=id_ptr+1;
        // A token index that is invariant across the whole file is in the
        // header already: measure its length, advance, spend NO bits at all.
        const bool skip = ELIDE && token_ctr<lay->ntok && lay->isconst[token_ctr];
        const bool skipsl = ELIDE && seqlen && token_ctr<lay->isseqlen.size()
                            && lay->isseqlen[token_ctr];
        if(skip){
            token_len=(uint32_t)lay->constant[token_ctr].size();
            id_ptr_tok=id_ptr+token_len;
        } else if(skipsl){
            // Hoisted: this index is the read's own length in every read, so
            // no type symbol and no payload -- just advance over the digits.
            while(isdigit((unsigned char)*id_ptr_tok)) ++id_ptr_tok;
            token_len=(uint32_t)(id_ptr_tok-id_ptr);
        } else if(isalpha((unsigned char)*id_ptr)){
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
            // The read's own length, already in the archive -- emit the type
            // symbol and nothing else. Guarded on the token ending here so the
            // decoder's digit run reconstructs to exactly these characters.
            const bool is_seqlen = seqlen && digit_value==seqlen
                                   && !isdigit((unsigned char)*id_ptr_tok)
                                   && digit_value>=(uint32_t)1;
            if(prev_is_digit && match_len==token_len && !isdigit((unsigned char)prevc(prev_tok_ptr[token_ctr]+token_len))){
                m.token_type[token_ctr].enc(enc,ID_MATCH);
            } else if(is_seqlen){
                m.token_type[token_ctr].enc(enc,ID_SEQLEN);
            } else {
                if(prev_is_digit){
                    ++m.seen[token_ctr];
                    if(delta>=-2048 && delta<=2048) ++m.hit[token_ctr];
                }
                const uint32_t seen=m.seen[token_ctr], hit=m.hit[token_ctr];
                const bool trust_wide = !NO_ZD && (seen>=20) && (hit*10 >= seen*3);
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
                    uint32_t sym = NO_DICT ? 0u : gdict[token_ctr].lookup(digit_value);
                    if(!NO_DICT) ldict[token_ctr]->encSym(enc,sym);
                    if(sym==0){   // pass 1 missed it; safety net
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

static std::string decompress_id(RangeDec& dec, IdModels& m, std::array<GlobalDict,MAXTOK>& gdict,
                                 std::array<LocalDictFreq*,MAXTOK>& ldict, std::string& prev_id,
                                 std::array<uint32_t,MAXTOK>& prev_tok_ptr,
                                 std::array<uint32_t,MAXTOK>& prev_tok_len,
                                 const Layout* lay=nullptr, uint32_t seqlen=0){
    const bool ELIDE = lay && lay->active();
    std::string id; id.reserve(64);
    uint32_t token_ctr=0, i=0;
    const char* prevbuf = prev_id.c_str();
    const size_t prevlen = prev_id.size();
    auto prevc=[&](uint32_t off)->char{ return off<prevlen ? prevbuf[off] : 0; };
    for(;;){
        // Constant index: text is in the header, consume no bits. Safe as a
        // silent append because such an index is present in EVERY read, so it
        // can never be where a shorter header would have stopped.
        if(ELIDE && seqlen && token_ctr<lay->isseqlen.size() && lay->isseqlen[token_ctr]){
            char b2[16]; int L2=snprintf(b2,sizeof(b2),"%u",seqlen);
            id.append(b2,(size_t)L2);
            prev_tok_ptr[token_ctr]=i; prev_tok_len[token_ctr]=(uint32_t)L2;
            i+=(uint32_t)L2; ++token_ctr;
            if(token_ctr>=MAXTOK-1) break;
            continue;
        }
        if(ELIDE && token_ctr<lay->ntok && lay->isconst[token_ctr]){
            const std::string& t=lay->constant[token_ctr];
            id.append(t);
            prev_tok_ptr[token_ctr]=i; prev_tok_len[token_ctr]=(uint32_t)t.size();
            i+=(uint32_t)t.size(); ++token_ctr;
            if(token_ctr>=MAXTOK-1) break;
            continue;
        }
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
            uint32_t sym = NO_DICT ? 0u : ldict[token_ctr]->decSym(dec);
            uint32_t v;
            if(sym==0){
                v = 0;
                v |= m.integer[token_ctr*4+0].dec(dec)<<0;
                v |= m.integer[token_ctr*4+1].dec(dec)<<8;
                v |= m.integer[token_ctr*4+2].dec(dec)<<16;
                v |= m.integer[token_ctr*4+3].dec(dec)<<24;
            } else {
                v = gdict[token_ctr].sym2val[sym-1];
            }
            char buf[16]; int L=snprintf(buf,sizeof(buf),"%u",v);
            id.append(buf,(size_t)L); token_len=(uint32_t)L;
        } else if(tok==ID_SEQLEN){
            // No payload: the value is this read's own sequence length, which
            // the archive already carries in read_lengths.
            char buf[16]; int L=snprintf(buf,sizeof(buf),"%u",seqlen);
            id.append(buf,(size_t)L); token_len=(uint32_t)L;
        } else if(tok==ID_DELTA){
            uint32_t delta = m.delta[token_ctr].dec(dec);
            uint32_t off=prev_tok_ptr[token_ctr], len=prev_tok_len[token_ctr];
            std::string prevtokstr; for(uint32_t k=0;k<len;++k) prevtokstr.push_back(prevc(off+k));
            uint32_t v=(uint32_t)atoi(prevtokstr.c_str())+delta;
            char buf[16]; int L=snprintf(buf,sizeof(buf),"%u",v);
            id.append(buf,(size_t)L); token_len=(uint32_t)L;
        } else if(tok==ID_ZDELTA){
            uint32_t hi=m.zdelta_hi[token_ctr].dec(dec);
            uint32_t lo=m.zdelta_lo[token_ctr].dec(dec);
            uint32_t z=(hi<<8)|lo;
            int64_t delta = (int64_t)(z>>1) ^ -(int64_t)(z&1);
            uint32_t off=prev_tok_ptr[token_ctr], len=prev_tok_len[token_ctr];
            std::string prevtokstr; for(uint32_t k=0;k<len;++k) prevtokstr.push_back(prevc(off+k));
            int64_t v=atoll(prevtokstr.c_str())+delta;
            char buf[24]; int L=snprintf(buf,sizeof(buf),"%lld",(long long)v);
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

// ---- bounded queue (stage 85's, with its in-flight accounting fix) ---------
struct Chunk { size_t block_idx=0; std::vector<std::string> lines;
               std::vector<uint32_t> slen; };   // seq length per read (for ID_SEQLEN)
template<typename T>
struct BoundedQueue {
    std::deque<T> q; size_t cap;
    std::mutex mu; std::condition_variable cv_push, cv_pop;
    bool closed=false; size_t inflight=0;
    BoundedQueue(size_t c):cap(c?c:1){}
    void push(T v){
        std::unique_lock<std::mutex> lk(mu);
        cv_push.wait(lk,[&]{ return inflight<cap; });
        ++inflight; q.push_back(std::move(v)); cv_pop.notify_one();
    }
    bool pop(T& out){
        std::unique_lock<std::mutex> lk(mu);
        cv_pop.wait(lk,[&]{ return !q.empty() || closed; });
        if(q.empty() && closed) return false;
        out=std::move(q.front()); q.pop_front();
        return true;   // slot frees only at release(), not here
    }
    void release(){ std::unique_lock<std::mutex> lk(mu); --inflight; cv_push.notify_one(); }
    void close(){ std::unique_lock<std::mutex> lk(mu); closed=true; cv_pop.notify_all(); }
};

// ---- archive-level entry points --------------------------------------------
// body  : concatenated per-block range-coded payloads (already entropy coded)
// dict  : RAW serialization of the global dictionary -- caller runs it
//         through best_encode(), like every other stream
// index : RAW varints, per block (byte length, name count) -- ditto
struct Encoded {
    std::vector<uint8_t> body, dict, index;
    uint64_t n_names=0, n_blocks=0, block=0;
    uint32_t tokens=0, const_tokens=0;   // reported: how much the header absorbed
    uint8_t line3_mode=2;                // 0='+' 1='+'+header 2=must be stored
};

static std::vector<uint8_t> serialize_dict(const std::array<GlobalDict,MAXTOK>& gdict,
                                          const Layout& lay){
    std::vector<uint8_t> o;
    // layout first: fixed flag, token count, then (index, len, bytes) per
    // file-constant token. Rides in the existing dict stream rather than
    // adding a fourth, so the archive gains no new section.
    pv(o, lay.active()?1:0);
    if(lay.active()){
        pv(o, lay.ntok);
        uint32_t nc=0; for(uint32_t k=0;k<lay.ntok;++k) if(lay.isconst[k]) ++nc;
        pv(o, nc);
        for(uint32_t k=0;k<lay.ntok;++k){
            if(!lay.isconst[k]) continue;
            pv(o,k); pv(o,lay.constant[k].size());
            o.insert(o.end(), lay.constant[k].begin(), lay.constant[k].end());
        }
        pv(o, lay.line3_mode);
        uint32_t ns=0;
        for(uint32_t k=0;k<lay.ntok;++k) if(k<lay.isseqlen.size() && lay.isseqlen[k]) ++ns;
        pv(o, ns);
        for(uint32_t k=0;k<lay.ntok;++k)
            if(k<lay.isseqlen.size() && lay.isseqlen[k]) pv(o,k);
    }
    uint32_t used=0;
    for(uint32_t k=0;k<MAXTOK;++k) if(!gdict[k].sym2val.empty()) ++used;
    pv(o,used);
    for(uint32_t k=0;k<MAXTOK;++k){
        const auto& g=gdict[k];
        if(g.sym2val.empty()) continue;
        pv(o,k); pv(o,g.sym2val.size());
        // Values at one token index are drawn from one field (a coordinate,
        // a lane, ...) so they cluster; sorting makes them a monotone
        // sequence and delta-varint shrinks it well before best_encode even
        // sees it. Symbol identity is the ORIGINAL index, so the permutation
        // has to travel with the values.
        std::vector<std::pair<uint32_t,uint32_t>> sv;   // (value, symbol-1)
        sv.reserve(g.sym2val.size());
        for(size_t s=0;s<g.sym2val.size();++s) sv.emplace_back(g.sym2val[s],(uint32_t)s);
        std::sort(sv.begin(),sv.end());
        uint64_t prev=0;
        for(auto& e:sv){ pv(o,(uint64_t)e.first-prev); prev=e.first; }
        for(auto& e:sv) pv(o,e.second);
    }
    return o;
}
static void deserialize_dict(const std::vector<uint8_t>& in, std::array<GlobalDict,MAXTOK>& gdict,
                             Layout& lay){
    if(in.empty()) return;
    const uint8_t* p=in.data(); const uint8_t* end=p+in.size();
    if(gv(p,end)){
        lay.ntok=(uint32_t)gv(p,end);
        const uint64_t nc=gv(p,end);
        lay.isconst.assign(lay.ntok,0);
        lay.constant.assign(lay.ntok,std::string());
        for(uint64_t c=0;c<nc && p<end;++c){
            const uint64_t k=gv(p,end); const uint64_t len=gv(p,end);
            if(k<lay.ntok && (size_t)(end-p)>=len){
                lay.isconst[k]=1; lay.constant[k].assign((const char*)p,(size_t)len);
            }
            p+=len;
        }
        lay.line3_mode=(uint8_t)gv(p,end);
        lay.isseqlen.assign(lay.ntok,0);
        const uint64_t ns=gv(p,end);
        for(uint64_t c=0;c<ns && p<end;++c){
            const uint64_t k=gv(p,end);
            if(k<lay.ntok) lay.isseqlen[k]=1;
        }
    }
    const uint64_t used=gv(p,end);
    for(uint64_t u=0;u<used && p<end;++u){
        const uint64_t k=gv(p,end); const uint64_t cnt=gv(p,end);
        if(k>=MAXTOK) return;
        std::vector<uint32_t> vals; vals.reserve(cnt);
        uint64_t prev=0;
        for(uint64_t i=0;i<cnt;++i){ prev+=gv(p,end); vals.push_back((uint32_t)prev); }
        std::vector<uint32_t> syms; syms.reserve(cnt);
        for(uint64_t i=0;i<cnt;++i) syms.push_back((uint32_t)gv(p,end));
        GlobalDict& g=gdict[k];
        g.sym2val.assign(cnt,0);
        for(uint64_t i=0;i<cnt;++i){
            if(syms[i]<cnt){ g.sym2val[syms[i]]=vals[i]; g.val2sym[vals[i]]=syms[i]+1; }
        }
    }
}

// Encode the ID column of a FASTQ. Two streaming passes over the file; at no
// point is more than (queue cap x block) names resident.
// BLOCK: names per parallel block. Each block starts with fresh adaptive
// models, so a finer split re-pays the model warm-up more often -- measured,
// not assumed: E. coli 500 K names costs 1,356,520 B at 100 K/block,
// 1,354,790 at 250 K and 1,354,147 as a single block, and DRR976266 moves
// 1,780,510 -> 1,778,873 the same way. 250 K names is ~16 MB of header text,
// which is also where Genozip sits (its own --STATS reports 18.0 MB vblocks),
// and it holds peak RSS to 118 MB. Not fitted to either file: both improve,
// and the single-block end of the sweep shows what the warm-up was costing.
static Encoded encode_from_fastq(const char* fq_path, size_t BLOCK=250000,
                                 unsigned NTHREADS=0, size_t QCAP=0){
    Encoded E; E.block=BLOCK;
    if(!NTHREADS) NTHREADS=std::max(1u,std::thread::hardware_concurrency());
    if(!QCAP) QCAP=NTHREADS;

    std::array<GlobalDict,MAXTOK> gdict;
    Layout lay;
    dict_pass(fq_path, gdict, &E.n_names, &lay);
    E.line3_mode = lay.line3_mode;
    E.dict = serialize_dict(gdict, lay);
    E.const_tokens=0;
    if(lay.active()) for(uint32_t k=0;k<lay.ntok;++k) if(lay.isconst[k]) ++E.const_tokens;
    E.tokens = lay.ntok;

    BoundedQueue<Chunk> workQ(QCAP);
    std::vector<std::vector<uint8_t>> blockOut;
    std::vector<uint64_t> blockCnt;
    std::mutex outMu;

    std::thread reader([&]{
        FILE* f=fopen(fq_path,"r"); if(!f){ workQ.close(); return; }
        std::vector<char> buf(1u<<16);
        uint64_t lineno=0; size_t block_idx=0;
        Chunk cur; cur.block_idx=0; cur.lines.reserve(BLOCK);
        while(fgets(buf.data(),(int)buf.size(),f)){
            const uint64_t ln=lineno++; 
            if(ln%4==1){ // SEQ line: its length is what "length=NNN" duplicates
                size_t L=strlen(buf.data()); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
                if(!cur.slen.empty()) cur.slen.back()=(uint32_t)L;
                continue;
            }
            if(ln%4!=0) continue;
            size_t L=strlen(buf.data()); while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
            const char* s=buf.data(); size_t off=(L&&s[0]=='@')?1:0;
            cur.lines.emplace_back(s+off,L-off); cur.slen.push_back(0);
            if(cur.lines.size()>=BLOCK){
                workQ.push(std::move(cur));
                ++block_idx; cur=Chunk(); cur.block_idx=block_idx; cur.lines.reserve(BLOCK); cur.slen.reserve(BLOCK);
            }
        }
        if(!cur.lines.empty()){ workQ.push(std::move(cur)); ++block_idx; }
        fclose(f);
        workQ.close();
        { std::lock_guard<std::mutex> lk(outMu);
          if(blockOut.size()<block_idx){ blockOut.resize(block_idx); blockCnt.resize(block_idx,0); } }
    });

    std::vector<std::thread> workers;
    for(unsigned t=0;t<NTHREADS;++t) workers.emplace_back([&]{
        Chunk c;
        while(workQ.pop(c)){
            IdModels m;
            std::array<LocalDictFreq*,MAXTOK> ldict{};
            for(uint32_t k=0;k<MAXTOK;++k) ldict[k]=new LocalDictFreq(gdict[k]);
            RangeEnc enc; enc.out.reserve(c.lines.size()*4);
            std::string prev_id; std::array<uint32_t,MAXTOK> pp{}; pp.fill(0);
            for(size_t z=0;z<c.lines.size();++z)
                compress_id(enc,m,gdict,ldict,c.lines[z].c_str(),prev_id,pp,&lay,
                            z<c.slen.size()?c.slen[z]:0u);
            enc.flush();
            { std::lock_guard<std::mutex> lk(outMu);
              if(blockOut.size()<=c.block_idx){ blockOut.resize(c.block_idx+1); blockCnt.resize(c.block_idx+1,0); }
              blockCnt[c.block_idx]=c.lines.size();
              blockOut[c.block_idx]=std::move(enc.out); }
            for(uint32_t k=0;k<MAXTOK;++k) delete ldict[k];
            std::vector<std::string>().swap(c.lines);
            workQ.release();
        }
    });
    reader.join();
    for(auto& w:workers) w.join();

    E.n_blocks=blockOut.size();
    for(size_t b=0;b<blockOut.size();++b){
        pv(E.index, blockOut[b].size());
        pv(E.index, blockCnt[b]);
        E.body.insert(E.body.end(), blockOut[b].begin(), blockOut[b].end());
        std::vector<uint8_t>().swap(blockOut[b]);
    }
    return E;
}

// Decode straight to a FILE*, one block at a time -- the decoder keeps the
// same bounded-memory property as the encoder. `at` prefixes each line with
// '@' when a FASTQ header column is wanted rather than a bare name column.
// `slens` supplies each read's sequence length, in original read order, for
// ID_SEQLEN. The archive already decodes read_lengths before it reaches the
// names column, so this costs the container nothing. Passing nullptr is only
// valid for inputs that contain no ID_SEQLEN token.
static uint64_t decode_to_file(const std::vector<uint8_t>& body,
                               const std::vector<uint8_t>& dictRaw,
                               const std::vector<uint8_t>& indexRaw,
                               FILE* out, bool at=false,
                               const std::vector<uint32_t>* slens=nullptr){
    std::array<GlobalDict,MAXTOK> gdict;
    Layout lay;
    deserialize_dict(dictRaw, gdict, lay);
    const uint8_t* ip=indexRaw.data(); const uint8_t* iend=ip+indexRaw.size();
    size_t boff=0; uint64_t written=0;
    std::string line;
    while(ip<iend){
        const uint64_t blen=gv(ip,iend);
        const uint64_t bcnt=gv(ip,iend);
        if(boff+blen>body.size()) break;
        IdModels m;
        std::array<LocalDictFreq*,MAXTOK> ldict{};
        for(uint32_t k=0;k<MAXTOK;++k) ldict[k]=new LocalDictFreq(gdict[k]);
        RangeDec dec; dec.init(body.data()+boff, (size_t)blen);
        std::string prev; std::array<uint32_t,MAXTOK> pp{}, pl{}; pp.fill(0); pl.fill(0);
        for(uint64_t r=0;r<bcnt;++r){
            const uint32_t sl = (slens && written<slens->size()) ? (*slens)[(size_t)written] : 0u;
            std::string got=decompress_id(dec,m,gdict,ldict,prev,pp,pl,&lay,sl);
            line.clear();
            if(at) line.push_back('@');
            line += got; line.push_back('\n');
            fwrite(line.data(),1,line.size(),out);
            ++written;
        }
        for(uint32_t k=0;k<MAXTOK;++k) delete ldict[k];
        boff += blen;
    }
    return written;
}

} // namespace nmc
