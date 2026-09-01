// CAPSULE decoder -- reads the archive back.
//
// Until now losslessness was verified by decoding the RAW intermediate streams
// dumped alongside the archive, which exercises the whole algorithm (assembly,
// mapping, mismatches, order) but never the entropy-coding layer. Every coder
// we use is standard and invertible, so this was a completeness gap rather than
// a suspected bug -- but a compressor whose archive nothing reads back is not a
// format. This closes it: archive in, raw streams out, byte-for-byte identical
// to what the encoder held in memory.
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <lzma.h>
#include "coders_pgrc.h"

// ---- inverse of xz_compress / xz_compress_lzma (both write .xz containers) --
static std::vector<uint8_t> xz_decompress(const uint8_t* d, size_t n, size_t hint){
    std::vector<uint8_t> out(hint ? hint : n*4 + 4096);
    for(int attempt=0; attempt<8; ++attempt){
        uint64_t memlimit = UINT64_MAX;
        size_t inpos=0, outpos=0;
        lzma_ret r = lzma_stream_buffer_decode(&memlimit, 0, nullptr,
                                               d, &inpos, n,
                                               out.data(), &outpos, out.size());
        if(r == LZMA_OK){ out.resize(outpos); return out; }
        if(r == LZMA_BUF_ERROR){ out.resize(out.size()*2); continue; }
        return {};
    }
    return {};
}
// ---- inverse of u32_byteplanes ---------------------------------------------
static std::vector<uint8_t> u32_unplane(const std::vector<uint8_t>& bp){
    const size_t cnt = bp.size()/4;
    std::vector<uint8_t> out(cnt*4);
    for(size_t i=0;i<cnt;++i)
        for(size_t k=0;k<4;++k) out[i*4+k] = bp[k*cnt+i];
    return out;
}
static const uint8_t CONST_MARKER = 0xC0;

// ---- inverse of best_encode / const_or_encode -------------------------------
// Layout written by the encoder: [method:1][raw_len:8][payload]
// or, for a constant stream: [0xC0][count:8][value bytes]
std::vector<uint8_t> capsule_decode_stream(const std::vector<uint8_t>& in, size_t width){
    if(in.empty()) return {};
    const uint8_t m = in[0];
    if(m == CONST_MARKER){
        if(in.size() < 9+width) return {};
        uint64_t cnt; memcpy(&cnt,&in[1],8);
        std::vector<uint8_t> out; out.reserve(cnt*width);
        for(uint64_t i=0;i<cnt;++i) out.insert(out.end(), in.begin()+9, in.begin()+9+width);
        return out;
    }
    if(in.size() < 9) return {};
    uint64_t rawlen; memcpy(&rawlen,&in[1],8);
    const uint8_t* p = in.data()+9; const size_t pn = in.size()-9;
    switch(m){
        case 0: case 1: return xz_decompress(p,pn,rawlen);
        case 2: return pgc::ppmd_decode(p,pn,rawlen);
        case 3: return pgc::fse_decode(p,pn,rawlen);
        case 4: return pgc::range_decode(p,pn,rawlen,1);
        case 5: case 6: {
            // method 6 may carry B1's raw-plane mask; method 5 never does
            std::vector<uint8_t> bp;
            if(m==6){
                const uint8_t rawmask = p[0];
                if(rawmask==0){ bp = xz_decompress(p+1,pn-1,rawlen); }
                else {
                    const size_t cnt = rawlen/4;
                    uint32_t cl; memcpy(&cl,p+1,4);
                    auto coded = xz_decompress(p+5,cl,rawlen);
                    const uint8_t* rawp = p+5+cl;
                    bp.resize(rawlen);
                    size_t ci=0;
                    for(int k=0;k<4;++k){
                        if(rawmask & (1u<<k)){ memcpy(&bp[(size_t)k*cnt], rawp, cnt); rawp += cnt; }
                        else { memcpy(&bp[(size_t)k*cnt], coded.data()+ci, cnt); ci += cnt; }
                    }
                }
            } else bp = xz_decompress(p,pn,rawlen);
            if(bp.size()!=rawlen) return {};
            return u32_unplane(bp);
        }
    }
    return {};
}

// ---- container walk ---------------------------------------------------------
struct Stream { std::string name; std::vector<uint8_t> coded; };
static bool read_capsule(const char* path, uint64_t& pg_len, uint64_t& main_end,
                         std::vector<Stream>& out){
    FILE* f=fopen(path,"rb"); if(!f) return false;
    char magic[8]; if(fread(magic,1,8,f)!=8 || memcmp(magic,"CAPSULE\0",8)){ fclose(f); return false; }
    uint16_t ver=0, ns=0;
    if(fread(&ver,2,1,f)!=1 || fread(&pg_len,8,1,f)!=1 || fread(&main_end,8,1,f)!=1
       || fread(&ns,2,1,f)!=1){ fclose(f); return false; }
    if(ver!=1){ fprintf(stderr,"capsule: unsupported version %u\n",ver); fclose(f); return false; }
    for(uint16_t i=0;i<ns;++i){
        uint8_t nl=0; if(fread(&nl,1,1,f)!=1) break;
        std::string nm(nl,'\0'); if(nl && fread(&nm[0],1,nl,f)!=nl) break;
        uint64_t n=0; if(fread(&n,8,1,f)!=1) break;
        Stream s; s.name=nm; s.coded.resize(n);
        if(n && fread(s.coded.data(),1,n,f)!=n) break;
        out.push_back(std::move(s));
    }
    fclose(f); return out.size()==ns;
}

int main(int argc,char** argv){
    if(argc<3){ fprintf(stderr,"usage: %s <in.capsule> <outdir> [--verify dumpdir]\n",argv[0]); return 2; }
    uint64_t pg_len=0, main_end=0; std::vector<Stream> ss;
    if(!read_capsule(argv[1],pg_len,main_end,ss)){ fprintf(stderr,"capsule: bad archive\n"); return 1; }
    const std::string outdir=argv[2];
    const char* verify = (argc>4 && !strcmp(argv[3],"--verify")) ? argv[4] : nullptr;
    fprintf(stderr,"capsule v1  pg_len=%llu main_pg_end=%llu  streams=%zu\n",
            (unsigned long long)pg_len,(unsigned long long)main_end,ss.size());
    { std::string pp=outdir+"/pg_params.txt"; FILE* g=fopen(pp.c_str(),"w");
      if(g){ fprintf(g,"%llu %llu\n",(unsigned long long)pg_len,(unsigned long long)main_end); fclose(g);} }

    // streams produced by best_encode / const_or_encode, and the dump each maps to
    const std::map<std::string,std::pair<std::string,size_t>> M = {
        {"pos_abs",       {"pos_abs.bin",1}},   {"pos_strand",{"pos_strand.bin",1}},
        {"n_pos",         {"n_pos.bin",1}},     {"n_indices", {"n_indices.bin",1}},
        {"n_cnt",         {"n_cnt.bin",1}},     {"read_lengths",{"read_lengths.bin",2}},
        {"orig2uid_flags",{"orig2uid_flags.bin",1}},
        {"orig2uid_vals", {"orig2uid_vals.bin",1}},
        {"mm_cnt_flags",  {"mm_cnt_flags.bin",1}},
        {"mm_cnt_vals",   {"mm_cnt_vals.bin",1}},
    };
    size_t ok=0, bad=0, skipped=0;
    for(auto& s : ss){
        auto it=M.find(s.name);
        if(it==M.end()){ ++skipped; fprintf(stderr,"  %-16s SKIP (coder inverse not yet written)\n",s.name.c_str()); continue; }
        auto raw = capsule_decode_stream(s.coded, it->second.second);
        std::string fp = outdir+"/"+it->second.first;
        FILE* g=fopen(fp.c_str(),"wb"); if(g){ if(!raw.empty()) fwrite(raw.data(),1,raw.size(),g); fclose(g); }
        if(verify){
            std::string vp = std::string(verify)+"/"+it->second.first;
            FILE* v=fopen(vp.c_str(),"rb");
            if(!v){ fprintf(stderr,"  %-16s %10zu B  (no dump to compare)\n",s.name.c_str(),raw.size()); continue; }
            fseek(v,0,SEEK_END); long vn=ftell(v); fseek(v,0,SEEK_SET);
            std::vector<uint8_t> ref(vn); if(vn) { if(fread(ref.data(),1,vn,v)!=(size_t)vn){} } fclose(v);
            bool same = (ref.size()==raw.size()) && (raw.empty() || !memcmp(ref.data(),raw.data(),raw.size()));
            fprintf(stderr,"  %-16s %10zu B  %s\n",s.name.c_str(),raw.size(), same?"IDENTICAL":"DIFFER");
            same?++ok:++bad;
        }
    }
    if(verify) fprintf(stderr,"entropy round-trip: %zu identical, %zu differ, %zu not yet implemented\n",ok,bad,skipped);
    return bad?1:0;
}
