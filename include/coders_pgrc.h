#pragma once
// Real entropy coders, matching what PgRC2 actually uses per stream.
//
// Until now our "per-stream coder selection" only chose among xz / bzip2 /
// zstd, which was selection in name only -- PgRC2's advantage on the streams
// where it still beat us comes from PPMd7, FSE/Huff0 and range coders with
// per-position models, not from a different general-purpose compressor.
// (See their PropsLibrary.cpp: mismatched symbols use PPMd order 5; mismatch
// counts use range_coder period 1; mismatch positions use a selector over
// range/FSE/PPMd per bucket.)
//
// PPMd7 (LZMA SDK, public domain) and FSE/Huff0 (Yann Collet, BSD) are
// standard third-party libraries that PgRC2 merely bundles -- using them is
// no different from our already using liblzma, and is not a copy of their
// implementation. The range coder below is our own, modelled on the same
// idea their rangecoder/simple_model.h uses.
#include <cstdint>
#include <cstring>
#include <vector>
#include <cstdio>

extern "C" {
#include "thirdparty/ppmd/Ppmd7.h"
#include "thirdparty/fse/fse.h"
#include "thirdparty/fse/huf.h"
#include "thirdparty/ppmd/Alloc.h"
}

namespace pgc {

// ---------------- PPMd7 ----------------
struct OutBuf { IByteOut vt; uint8_t* cur; uint8_t* lim; size_t written=0; bool overflow=false; };
static void ppmd_write(const IByteOut* p, Byte b){
    OutBuf* o = (OutBuf*)((char*)p - offsetof(OutBuf, vt));
    if(o->cur == o->lim){ o->overflow=true; return; }
    *o->cur++ = b; ++o->written;
}
// order/memSize mirror PgRC2's getDefaultCoderProps(PPMD7_CODER, level, order).
static std::vector<uint8_t> ppmd_encode(const uint8_t* src, size_t n,
                                        unsigned order=5, uint32_t memMB=32){
    if(!n) return {};
    CPpmd7 ppmd;
    Ppmd7_Construct(&ppmd);
    uint32_t memSize = memMB << 20;
    if(!Ppmd7_Alloc(&ppmd, memSize, &g_Alloc)) return {};
    std::vector<uint8_t> out(n + n/3 + 256);
    OutBuf ob; ob.vt.Write = ppmd_write;
    ob.cur = out.data()+5; ob.lim = out.data()+out.size();
    out[0]=(uint8_t)order;
    out[1]=(uint8_t)(memSize); out[2]=(uint8_t)(memSize>>8);
    out[3]=(uint8_t)(memSize>>16); out[4]=(uint8_t)(memSize>>24);
    Ppmd7z_Init_RangeEnc(&ppmd);
    Ppmd7_Init(&ppmd, order);
    ppmd.rc.enc.Stream = &ob.vt;
    Ppmd7z_EncodeSymbols(&ppmd, src, src+n);
    Ppmd7z_Flush_RangeEnc(&ppmd);
    Ppmd7_Free(&ppmd, &g_Alloc);
    if(ob.overflow) return {};
    out.resize(5 + ob.written);
    return out;
}

// ---------------- FSE / Huff0 ----------------
static std::vector<uint8_t> fse_encode(const uint8_t* src, size_t n){
    if(!n) return {};
    std::vector<uint8_t> out(FSE_compressBound(n)+16);
    size_t r = FSE_compress(out.data(), out.size(), src, n);
    if(FSE_isError(r) || r==0 || r>=n) return {};   // 0 = not compressible, 1 = RLE
    out.resize(r); return out;
}
static std::vector<uint8_t> huf_encode(const uint8_t* src, size_t n){
    if(!n) return {};
    std::vector<uint8_t> out(HUF_compressBound(n)+16);
    size_t r = HUF_compress(out.data(), out.size(), src, n);
    if(HUF_isError(r) || r==0 || r>=n) return {};
    out.resize(r); return out;
}

// ---------------- adaptive range coder, order-0 with `period` models -------
// PgRC2's RangeCoderCompressTemplate interleaves `period` independent models
// so that column k of a fixed-stride record gets its own statistics. That is
// the transpose idea expressed inside the coder.
struct RC {
    std::vector<uint8_t> out; uint64_t low=0; uint32_t range=0xFFFFFFFFu;
    uint8_t cache=0; uint64_t cacheSize=1;
    void shiftLow(){
        if((uint32_t)(low>>32)!=0 || (uint32_t)low < 0xFF000000u){
            uint8_t t=cache; do{ out.push_back((uint8_t)(t+(uint8_t)(low>>32))); t=0xFF; }while(--cacheSize);
            cache=(uint8_t)((uint32_t)low>>24);
        }
        ++cacheSize; low=(uint64_t)((uint32_t)low<<8);
    }
    void encode(uint32_t lo,uint32_t hi,uint32_t tot){
        range/=tot; low+=(uint64_t)lo*range; range*=(hi-lo);
        while(range<(1u<<24)){ range<<=8; shiftLow(); }
    }
    void flush(){ for(int i=0;i<5;++i) shiftLow(); }
};
struct Model256 {
    uint16_t f[256]; uint32_t tot;
    Model256(){ for(int i=0;i<256;++i) f[i]=1; tot=256; }
    void encode(RC& rc, uint8_t s){
        uint32_t lo=0; for(int i=0;i<s;++i) lo+=f[i];
        rc.encode(lo, lo+f[s], tot);
        f[s]+=32; tot+=32;
        if(tot>60000){ tot=0; for(int i=0;i<256;++i){ f[i]=(uint16_t)((f[i]>>1)|1); tot+=f[i]; } }
    }
};
static std::vector<uint8_t> range_encode(const uint8_t* src, size_t n, unsigned period=1){
    if(!n) return {};
    if(period<1) period=1;
    std::vector<Model256> m(period);
    RC rc; rc.out.reserve(n/2+64);
    for(size_t i=0;i<n;++i) m[i%period].encode(rc, src[i]);
    rc.flush();
    return std::move(rc.out);
}

} // namespace pgc
