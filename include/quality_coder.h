#pragma once
// Quality coder for the CAPSULE archive, wrapping the vendored fqzcomp_qual
// codec (thirdparty/htscodecs, BSD 3-clause, James Bonfield / Genome Research
// Ltd -- see thirdparty/htscodecs/LICENSE.md). Namespaced `qlc`, following the
// project's per-coder convention (mmc, refc, pgc, nmc).
//
// WHY VENDOR RATHER THAN REIMPLEMENT (measured, see docs/REIMPL_NOTES.md):
// our own quality coder (stages 67-92) beats SPRING 7/8 and Genozip 8/8 but
// loses to real fqzcomp on 8/8 by 1.1-4.5%, while fqzcomp is also 2x faster
// single-threaded than ours is on 12 threads at the same RAM. htscodecs is BSD,
// unlike PgRC2 (GPL-3) which is why THAT had to be reimplemented. Every
// cross-column signal we hold and fqzcomp cannot see was tested by held-out
// entropy and refuted: tile +12.20% worse, base call +1.04% worse, is-N 0.00%.
//
// TWO THINGS DELIBERATELY NOT DONE, both recorded so they are not "fixed" later
// by someone assuming they were oversights:
//
// 1. flags[] is left at 0, NOT fed from pos_strand. FQZ_FREVERSE exists because
//    in BAM/CRAM the stored quality of a reverse-strand alignment is already
//    reversed relative to the original read; fqzcomp un-reverses it to restore
//    the position-quality correlation (quality degrades toward the 3' end). In
//    FASTQ the quality is ALREADY in original orientation, and our pos_strand
//    describes pseudogenome placement, not read orientation -- feeding it would
//    reverse strings that were never reversed and DESTROY that correlation.
//    FQZ_FREAD2 likewise applies to paired files, which this single-file path
//    does not produce.
//
// 2. Record lengths are still stored by fqzcomp, not suppressed in favour of
//    the archive's own read_lengths. fqzcomp auto-detects constant length and
//    then stores exactly one length for the whole block (fqz_pick_parameters
//    sets fixed_len, encoder guards on `!pm->fixed_len || state->first_len`),
//    so on fixed-length files this already costs ~nothing. Suppressing it for
//    variable-length files means patching vendored code and keeping encoder and
//    decoder in lockstep forever; that is only worth doing if measured to pay.
//
// Block layout mirrors the names column exactly: blocks are independent, and a
// real index stream carries per-block (compressed length, read count) because a
// decoder holding only the archive cannot recover them otherwise. That is the
// same omission that has already cost this project twice.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

extern "C" {
#include "fqzcomp_qual.h"
}

namespace qlc {

// ---- varint helpers (same encoding as nmc) ---------------------------------
static inline void pv(std::vector<uint8_t>& o, uint64_t v){
    while(v>=0x80){ o.push_back((uint8_t)(v|0x80)); v>>=7; }
    o.push_back((uint8_t)v);
}
static inline uint64_t gv(const uint8_t*& p, const uint8_t* end){
    uint64_t v=0; int s=0;
    while(p<end){ uint8_t b=*p++; v |= (uint64_t)(b&0x7f)<<s; if(!(b&0x80)) break; s+=7; }
    return v;
}

struct Encoded {
    std::vector<uint8_t> body, index;
    uint64_t n_reads=0, n_blocks=0, n_qbytes=0;
};

// fqzcomp's own cap is BLK_SIZE = 300 MB of quality per call. 256 MB keeps
// every dataset in the locked set to a single block (largest quality column
// measured: 145 MB) -- which is also the best ratio, since a block pays its
// model warm-up once -- while still bounding memory on inputs of any size.
static const size_t QBLOCK_BYTES = 256u*1024u*1024u;

// Compress one block, trying every fqzcomp strategy and keeping the smallest.
// The compressed stream is self-describing (fqz_decompress takes no strategy
// argument), so the winner needs no flag of its own. Measured worth 1.5-3.0%
// against fqzcomp's own default on 3 of 8 datasets.
// Quality is handed to fqzcomp as raw values with the block's own minimum
// subtracted, NOT as raw ASCII. fqzcomp sizes its models on max_sym, so ASCII
// (35..74 on real data) makes it carry nearly twice the symbol space of the
// phred range it actually needs, diluting every context. CRAM feeds it raw
// phred for the same reason; subtracting the observed minimum is strictly
// tighter than a fixed -33 and is a pure bijection, so it cannot lose
// information. The offset rides in the index, one varint per block.
static bool encode_block(const std::string& qbuf,
                         const std::vector<uint32_t>& lens,
                         std::vector<uint8_t>& out,
                         uint8_t& qmin_out)
{
    if(lens.empty()) return true;
    uint8_t qmin=255;
    for(unsigned char c : qbuf) if(c<qmin) qmin=c;

    std::vector<uint32_t> flags(lens.size(), 0u);   // see note 1 in the header
    fqz_slice s;
    s.num_records = (int)lens.size();
    s.len   = const_cast<uint32_t*>(lens.data());
    s.flags = flags.data();

    // Two candidate offsets, both bijections: the block's own minimum (the
    // tightest alphabet) and the FASTQ standard 33 (raw phred, what CRAM
    // feeds it). Tighter is NOT reliably smaller -- fqzcomp's qmap/qshift
    // selection is not monotone in max_sym -- so both are tried and the
    // smaller kept, the same "try it, keep it if it measures smaller"
    // discipline used for the stream coders and for MAXMAP. The winning
    // offset is already carried per block in the index, so this costs no
    // extra format.
    uint8_t cands[2] = { qmin, 33 };
    const int ncand = (qmin!=33 && qmin>=33) ? 2 : 1;

    bool got=false; size_t best=0;
    std::string shifted; shifted.resize(qbuf.size());
    for(int ci=0; ci<ncand; ++ci){
        const uint8_t off=cands[ci];
        for(size_t i=0;i<qbuf.size();++i) shifted[i]=(char)((unsigned char)qbuf[i]-off);
        for(int strat=0; strat<4; ++strat){
            size_t osz=0;
            char* c = fqz_compress(4 /*CRAM 4.0 vers*/, &s,
                                   const_cast<char*>(shifted.data()), shifted.size(),
                                   &osz, strat, nullptr);
            if(!c) continue;
            if(!got || osz<best){
                out.assign((uint8_t*)c, (uint8_t*)c+osz);
                best=osz; got=true; qmin_out=off;
            }
            free(c);
        }
    }
    return got;
}

// Streams the quality column straight out of the FASTQ, one block at a time --
// never materializing the whole column, matching the bounded-memory property
// the names coder has.
static Encoded encode_from_fastq(const char* fq_path, size_t BLOCK_BYTES=QBLOCK_BYTES){
    Encoded E;
    FILE* f=fopen(fq_path,"r");
    if(!f) return E;
    std::vector<char> buf(1u<<16);
    std::string qbuf; qbuf.reserve(BLOCK_BYTES+65536);
    std::vector<uint32_t> lens;
    uint64_t lineno=0;

    auto flush=[&](){
        if(lens.empty()) return;
        std::vector<uint8_t> blk; uint8_t qmin=0;
        if(encode_block(qbuf, lens, blk, qmin)){
            pv(E.index, blk.size());
            pv(E.index, lens.size());
            pv(E.index, qmin);
            E.body.insert(E.body.end(), blk.begin(), blk.end());
            ++E.n_blocks;
        }
        E.n_reads  += lens.size();
        E.n_qbytes += qbuf.size();
        qbuf.clear(); lens.clear();
    };

    while(fgets(buf.data(),(int)buf.size(),f)){
        const uint64_t ln = lineno++;
        if(ln%4!=3) continue;                       // quality is line 4
        size_t L=strlen(buf.data());
        while(L&&(buf[L-1]=='\n'||buf[L-1]=='\r')) --L;
        if(qbuf.size()+L > BLOCK_BYTES && !lens.empty()) flush();
        qbuf.append(buf.data(), L);
        lens.push_back((uint32_t)L);
    }
    flush();
    fclose(f);
    return E;
}

// Inverts the above from the archive alone. Record boundaries come from the
// caller's read_lengths (which the archive decodes before it reaches this
// column), so fqz_decompress is asked for no length array of its own.
static uint64_t decode_to_file(const std::vector<uint8_t>& body,
                               const std::vector<uint8_t>& indexRaw,
                               FILE* out,
                               const std::vector<uint32_t>& lengths)
{
    const uint8_t* ip=indexRaw.data(); const uint8_t* iend=ip+indexRaw.size();
    size_t boff=0; uint64_t written=0, li=0;
    while(ip<iend){
        const uint64_t blen=gv(ip,iend);
        const uint64_t bcnt=gv(ip,iend);
        const uint8_t  qmin=(uint8_t)gv(ip,iend);
        if(boff+blen>body.size()) break;
        size_t osz=0;
        char* d = fqz_decompress((char*)body.data()+boff, blen, &osz, nullptr, 0);
        if(!d) break;
        for(size_t i=0;i<osz;++i) d[i]=(char)((unsigned char)d[i]+qmin);
        size_t off=0;
        for(uint64_t r=0; r<bcnt && li<lengths.size(); ++r,++li){
            const uint32_t L=lengths[li];
            if(off+L>osz) break;
            fwrite(d+off,1,L,out); fputc('\n',out);
            off+=L; ++written;
        }
        free(d);
        boff+=blen;
    }
    return written;
}

} // namespace qlc
