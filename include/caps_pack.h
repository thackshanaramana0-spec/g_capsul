#pragma once
// ═══════════════════════════════════════════════════════════════════════════
//  The read and quality encodings the graph caller consumes.
//
//  ONE DEFINITION, TWO CALLERS. The encoder builds these while compressing;
//  capsule_decode builds them from a stored archive. If the two ever disagree
//  by a byte the caller silently sees different data and the archive path
//  produces different variants -- and this project has already lost reads to
//  exactly this class of bug once, when packing mapped non-ACGT to 'A' and
//  added 2,273 k-mers that no F1 check could see. Hence: shared, not copied.
//
//  SEQUENCE  pure ACGT : [0][len_lo][len_hi][2-bit packed, 4 bases/byte]
//            has non-ACGT: [1][the raw bases]
//            The raw fallback exists because 2-bit packing cannot represent N,
//            and mapping N to a real base invents k-mers that were never read.
//
//  QUALITY   one bit per base, set when that base met QMIN. The caller asks
//            only "did this base clear the threshold", so a bitmap answers it
//            in 1/8th the memory of the phred string (233 MB vs 2.27 GB at
//            full chr20).
// ═══════════════════════════════════════════════════════════════════════════
#include <string>
#include <cstddef>

namespace capspack {

// Sequence -> the caller's packed form. Mirrors caps_caller.h::unpack_read.
inline std::string pack_seq(const std::string& b){
    bool pure = true;
    for(size_t i=0;i<b.size();++i){
        const char c=b[i];
        if(c!='A'&&c!='C'&&c!='G'&&c!='T'&&c!='a'&&c!='c'&&c!='g'&&c!='t'){ pure=false; break; }
    }
    if(!pure){
        std::string raw; raw.reserve(1+b.size());
        raw.push_back((char)1);
        raw += b;
        return raw;
    }
    std::string pk(3 + (b.size()+3)/4, '\0');
    pk[0] = (char)0;
    pk[1] = (char)(b.size() & 0xFF);
    pk[2] = (char)((b.size() >> 8) & 0xFF);
    for(size_t i=0;i<b.size();++i){
        int bb;
        switch(b[i]){ case 'A': case 'a': bb=0; break;
                      case 'C': case 'c': bb=1; break;
                      case 'G': case 'g': bb=2; break;
                      default:            bb=3; }
        pk[3 + (i>>2)] = (char)(pk[3 + (i>>2)] | (bb << (2*(i&3))));
    }
    return pk;
}

// Phred string -> one bit per base, set when the base met qmin (phred+33).
inline std::string pack_qual(const std::string& d, int qmin){
    std::string packed((d.size()+7)/8, '\0');
    for(size_t i=0;i<d.size();++i)
        if((int)(d[i]-33) >= qmin) packed[i>>3] |= (char)(1u << (i&7));
    return packed;
}

} // namespace capspack
