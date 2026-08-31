// STAGE 61 -- read-name/header coder, the first stage of method-c's
// extension PAST PgRC2's scope. Verified directly: PgRC2's own decoder
// (method_c/pgrc/pgrc-decoder.cpp) writes sequence-only output -- no '@'
// line, no quality line, anywhere -- so PgRC2 has no names/quality stages to
// reimplement. SPRING and Genozip DO handle both; this is the next real
// milestone on the way to being comparable with them, not PgRC2 anymore.
//
// Technique: fqzcomp's real algorithm (jkbonfield/fqzcomp, cloned and read
// line by line, not from a summary). Per name: find the length of the
// shared PREFIX with the PREVIOUS name, and the shared SUFFIX; code those
// two lengths plus the total length, each with its own adaptive model
// keyed on the previous read's corresponding value; code only the MIDDLE
// bytes, one at a time, with a context built from the corresponding
// position in the previous name PLUS a realignment counter that resyncs on
// ':' and ' ' delimiters, so a field that changes WIDTH (3-digit tile
// becoming 4-digit) doesn't throw off every byte after it. No explicit
// tokenizer -- the delimiter-triggered realignment does that job implicitly.
//
//   g++ -O3 -march=native -o namecoder 61_name_coder.cpp
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <array>

// ── adaptive order-0 model over a small alphabet, range-coded ──────────────
// (same shape as this project's other coders: RangeEnc/RangeDec, adaptive
// frequency table, periodic rescale)
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

// A context slot: adaptive frequency table over a byte alphabet (0..255),
// but sized down to a practical symbol set to keep tables small and fast to
// adapt -- names are printable ASCII (33..126), so this maps into a 96-entry
// alphabet (space through '~').
static const int ALPHA_LO=32, ALPHA_N=96;   // ' '(32) .. 127 exclusive
struct ByteModel {
    std::vector<uint16_t> f; // ALPHA_N freqs
    ByteModel():f(ALPHA_N,1){}
    void enc(RangeEnc& rc,int sym){
        int idx=sym-ALPHA_LO; if(idx<0) idx=0; if(idx>=ALPHA_N) idx=ALPHA_N-1;
        uint32_t tot=0; for(int i=0;i<ALPHA_N;++i) tot+=f[i];
        uint32_t lo=0; for(int i=0;i<idx;++i) lo+=f[i];
        rc.encode(lo,lo+f[idx],tot);
        f[idx]+=32; if(tot+32>60000){ for(int i=0;i<ALPHA_N;++i) f[i]=(f[i]>>1)|1; }
    }
    int dec(RangeDec& rc){
        uint32_t tot=0; for(int i=0;i<ALPHA_N;++i) tot+=f[i];
        uint32_t target=rc.getFreq(tot);
        int idx=0; uint32_t lo=0; while(lo+f[idx]<=target){ lo+=f[idx]; ++idx; }
        rc.decodeUpdate(lo,lo+f[idx]);
        f[idx]+=32; if(tot+32>60000){ for(int i=0;i<ALPHA_N;++i) f[i]=(f[i]>>1)|1; }
        return idx+ALPHA_LO;
    }
};
// small-integer model (lengths, prefix/suffix counts), alphabet 0..255
struct LenModel {
    std::vector<uint32_t> f;
    LenModel():f(256,1){}
    void enc(RangeEnc& rc,int v){
        if(v<0) v=0; if(v>255) v=255;
        uint32_t tot=0; for(int i=0;i<256;++i) tot+=f[i];
        uint32_t lo=0; for(int i=0;i<v;++i) lo+=f[i];
        rc.encode(lo,lo+f[v],tot);
        f[v]+=32; if(tot+32>60000){ for(int i=0;i<256;++i) f[i]=(f[i]>>1)|1; }
    }
    int dec(RangeDec& rc){
        uint32_t tot=0; for(int i=0;i<256;++i) tot+=f[i];
        uint32_t target=rc.getFreq(tot);
        int v=0; uint32_t lo=0; while(lo+f[v]<=target){ lo+=f[v]; ++v; }
        rc.decodeUpdate(lo,lo+f[v]);
        f[v]+=32; if(tot+32>60000){ for(int i=0;i<256;++i) f[i]=(f[i]>>1)|1; }
        return v;
    }
};

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

    LenModel model_prefix, model_suffix, model_len;
    // STAGE 61 FIX: fqzcomp's real design keys the middle-byte model by a
    // CONTEXT HASH -- the corresponding character in the last name, whether
    // the previous byte matched, and a position counter -- not one shared
    // model. This is what lets it predict "this exact position, after this
    // exact preceding context, is virtually always X" instead of learning
    // only the aggregate character distribution across the whole format.
    // 8192 contexts, matching fqzcomp's own modulus exactly.
    std::vector<ByteModel> model_middle(8192);

    auto run=[&](bool decode, RangeEnc* enc, RangeDec* dec, std::vector<std::string>* outv)->void{
        std::string last; int last_p=0;
        for(size_t r=0;r<n;++r){
            const std::string& cur = decode ? std::string() : names[r]; // encode reads real name
            int len, p_len, s_len;
            std::string decoded;
            if(!decode){
                const std::string& name=names[r];
                len=(int)name.size();
                int i=0; while(i<len && i<(int)last.size() && name[i]==last[i]) ++i;
                p_len=i;
                int a=len-1,b=(int)last.size()-1;
                while(a>=0 && b>=0 && name[a]==last[b]){ --a; --b; }
                s_len=len-1-a;
                if(len - s_len - p_len < 0) s_len = len - p_len;
                model_prefix.enc(*enc,p_len);
                model_suffix.enc(*enc,s_len);
                model_len.enc(*enc,len);
                int len2=len-s_len, lc2=p_len?1:0;
                int i2=p_len, j=p_len, k=0;
                for(; i2<len2; ++i2, ++j, ++k){
                    const int lastc = (j>=0 && j<(int)last.size()) ? (unsigned char)last[j] : 32;
                    const int ctx = (((lastc-32)*2 + lc2 + k*64) % 8192 + 8192) % 8192;
                    model_middle[ctx].enc(*enc, (unsigned char)name[i2]);
                    const bool wasmatch = (j<(int)last.size() && name[i2]==last[j]);
                    if(name[i2]==' ' && j<(int)last.size() && last[j]!=' ') ++j;
                    if(name[i2]!=' ' && j<(int)last.size() && last[j]==' ') --j;
                    if(name[i2]==':' && j<(int)last.size() && last[j]!=':') ++j;
                    if(name[i2]!=':' && j<(int)last.size() && last[j]==':') --j;
                    if(name[i2]==':' || name[i2]==' ') k=((k+3)>>2<<2)-1;   // -1: loop's ++k restores it
                    lc2 = wasmatch;
                }
                last=name; last_p=p_len;
            } else {
                p_len=model_prefix.dec(*dec);
                s_len=model_suffix.dec(*dec);
                len=model_len.dec(*dec);
                decoded.resize(len);
                for(int i=0;i<p_len && i<(int)last.size();++i) decoded[i]=last[i];
                int len2=len-s_len, lc2=p_len?1:0;
                int i2=p_len, j=p_len, k=0;
                for(; i2<len2; ++i2, ++j, ++k){
                    const int lastc = (j>=0 && j<(int)last.size()) ? (unsigned char)last[j] : 32;
                    const int ctx = (((lastc-32)*2 + lc2 + k*64) % 8192 + 8192) % 8192;
                    int c=model_middle[ctx].dec(*dec);
                    decoded[i2]=(char)c;
                    const bool wasmatch = (j<(int)last.size() && c==(unsigned char)last[j]);
                    if(c==' ' && j<(int)last.size() && last[j]!=' ') ++j;
                    if(c!=' ' && j<(int)last.size() && last[j]==' ') --j;
                    if(c==':' && j<(int)last.size() && last[j]!=':') ++j;
                    if(c!=':' && j<(int)last.size() && last[j]==':') --j;
                    if(c==':' || c==' ') k=((k+3)>>2<<2)-1;
                    lc2 = wasmatch;
                }
                // suffix: copy from END of last name
                for(int k=0;k<s_len;++k){
                    int srcpos=(int)last.size()-s_len+k;
                    if(srcpos>=0 && srcpos<(int)last.size()) decoded[len-s_len+k]=last[srcpos];
                }
                last=decoded;
                outv->push_back(decoded);
            }
        }
    };

    RangeEnc enc; enc.out.reserve(n*8);
    run(false,&enc,nullptr,nullptr);
    enc.flush();
    const size_t bytes=enc.out.size();

    // reset models for decode
    model_prefix=LenModel(); model_suffix=LenModel(); model_len=LenModel(); model_middle.assign(8192,ByteModel());
    RangeDec dec; dec.init(enc.out.data(),enc.out.size());
    std::vector<std::string> outv; outv.reserve(n);
    run(true,nullptr,&dec,&outv);

    bool ok = (outv.size()==n);
    if(ok) for(size_t i=0;i<n;++i) if(outv[i]!=names[i]){ ok=false;
        fprintf(stderr,"MISMATCH at %zu:\n  got : %s\n  want: %s\n",i,outv[i].c_str(),names[i].c_str()); break; }

    printf("names=%zu  coded=%zu B  %.3f B/name  round trip: %s\n",
           n,bytes,bytes/(double)n, ok?"VERIFIED":"FAILED");
    return ok?0:1;
}
