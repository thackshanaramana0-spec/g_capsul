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
#include <chrono>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <vector>
#include <map>
#include <lzma.h>
#include <fstream>
#include <sys/stat.h>
#include <omp.h>
#include "caps_pack.h"
#include "caps_caller.h"
#include "coders_pgrc.h"
#include "seqpar_core.h"
#include "coders_inproc.h"
#include "names_coder.h"
#include "quality_coder.h"

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
        case 7: {
            // Chunked: each chunk is itself a complete best_encode output, so
            // this recurses rather than duplicating the per-method logic.
            if(pn < 4) return {};
            uint32_t nc; memcpy(&nc,p,4);
            size_t off = 4;
            if((uint64_t)off + 4ull*nc > pn) return {};
            std::vector<uint32_t> lens(nc);
            for(uint32_t c=0;c<nc;++c){ memcpy(&lens[c], p+off, 4); off+=4; }
            std::vector<uint8_t> all; all.reserve(rawlen);
            for(uint32_t c=0;c<nc;++c){
                if((uint64_t)off + lens[c] > pn) return {};
                std::vector<uint8_t> part(p+off, p+off+lens[c]);
                auto dec = capsule_decode_stream(part, width);
                all.insert(all.end(), dec.begin(), dec.end());
                off += lens[c];
            }
            return all;
        }
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
                         uint32_t& minmem, std::vector<Stream>& out){
    FILE* f=fopen(path,"rb"); if(!f) return false;
    char magic[8]; if(fread(magic,1,8,f)!=8 || memcmp(magic,"CAPSULE\0",8)){ fclose(f); return false; }
    uint16_t ver=0, ns=0;
    if(fread(&ver,2,1,f)!=1 || fread(&pg_len,8,1,f)!=1 || fread(&main_end,8,1,f)!=1){ fclose(f); return false; }
    if(ver!=2){ fprintf(stderr,"capsule: unsupported version %u (expect 2)\n",ver); fclose(f); return false; }
    if(fread(&minmem,4,1,f)!=1 || fread(&ns,2,1,f)!=1){ fclose(f); return false; }
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

// Decode every stream, rebuild the pseudogenome, and write the raw streams the
// read-reconstruction step consumes. The ONLY input is the archive.
//
// mm_sym is coded against ref under an adaptive model, and ref is never stored:
// it is pg[read position + mismatch offset]. So the order is forced -- rebuild
// the pseudogenome first, derive every ref in the encoder's own order (unique
// read ascending, mismatch offset ascending), then decode the observed bases.
static void put_file(const std::string& p, const void* d, size_t n){
    FILE* f=fopen(p.c_str(),"wb"); if(!f) return;
    if(n) fwrite(d,1,n,f); fclose(f);
}
static std::vector<uint64_t> varints(const std::vector<uint8_t>& v){
    std::vector<uint64_t> o; size_t i=0;
    while(i<v.size()){ uint64_t x=0; int sh=0;
        while(i<v.size()){ uint8_t b=v[i++]; x|=(uint64_t)(b&0x7f)<<sh; sh+=7; if(!(b&0x80)) break; }
        o.push_back(x); }
    return o;
}

// ── CLAIM 3: ADDRESSABLE ────────────────────────────────────────────────────
// export / coverage / query are served DIRECTLY from the archive. The work a
// conventional pipeline does at query time -- assembling contigs, or indexing
// and aligning reads to compute depth -- CAPSULE already did at compress time,
// and stored. So these are stream decodes, not computations:
//   export   : literal + mem_triples  -> the pseudogenome         (no assembly)
//   coverage : pos_abs + read_lengths -> per-position depth       (no alignment)
//   query    : pos_abs + pg           -> reads overlapping a range (no full decode)
// Each stops as soon as the streams it needs are decoded, so none pays for the
// full reconstruction.
// IN-MEMORY HANDOFF (out_flat/out_rowoff). The decoder already materialises
// every read into ONE contiguous buffer; writing that 1.86 GB out just so the
// caller can read it straight back costs a full write plus a full re-parse for
// nothing. When these are non-null the buffer is handed over directly and the
// file is not written at all. Every other caller passes nothing and is
// byte-for-byte unaffected.
int capsule_decode_all(const char* arcpath, const std::string& outdir,
                       const std::string& outreads = std::string(),
                       const std::string& mode = std::string(),
                       const std::string& modearg = std::string(),
                       std::vector<uint8_t>* out_flat = nullptr,
                       std::vector<size_t>* out_rowoff = nullptr,
                       std::vector<uint8_t>* out_qflat = nullptr,
                       std::vector<std::string>* out_qbits = nullptr,
                       int qbits_qmin = 20,
                       std::vector<std::string>* out_qtext = nullptr,
                       std::vector<std::string>* out_contigs = nullptr){
    auto _dt0 = std::chrono::steady_clock::now();
    // [DEC-RSS] temporary: the run's PEAK is 13.2 GB while the caller's own
    // RSS is 4.8 GB, so ~8.4 GB is transient in here. Find out where.
    auto _rssmb = []() -> long { FILE* f=fopen("/proc/self/status","r"); if(!f) return -1;
        char l[256]; long kb=-1; while(fgets(l,sizeof l,f)) if(!strncmp(l,"VmRSS:",6)){sscanf(l+6,"%ld",&kb);break;}
        fclose(f); return kb/1024; };
    auto _rsspk = []() -> long { FILE* f=fopen("/proc/self/status","r"); if(!f) return -1;
        char l[256]; long kb=-1; while(fgets(l,sizeof l,f)) if(!strncmp(l,"VmHWM:",6)){sscanf(l+6,"%ld",&kb);break;}
        fclose(f); return kb/1024; };
    auto _rss = [&](const char* where){ fprintf(stderr,"[DEC-RSS] %-26s rss=%ldMB peak=%ldMB\n", where, _rssmb(), _rsspk()); };
    _rss("entry");
    uint64_t PGLEN=0, MAINEND=0; uint32_t MINMEM=0; std::vector<Stream> ss;
    if(!read_capsule(arcpath,PGLEN,MAINEND,MINMEM,ss)){ fprintf(stderr,"bad archive\n"); return 1; }
    std::map<std::string,std::vector<uint8_t>> S;
    for(auto& x:ss) S[x.name]=std::move(x.coded);
    auto has=[&](const char* n){ return S.count(n)>0; };
    auto dec=[&](const char* n,size_t w=1){ return has(n)?capsule_decode_stream(S[n],w):std::vector<uint8_t>(); };

    // ---- literal: 2-bit codes -> ACGT --------------------------------------
    auto litcode = seq_decode_mem(S["literal"].data(), S["literal"].size());
    std::vector<uint8_t> literal(litcode.size());
    { const char M[4]={'A','C','G','T'};
      for(size_t i=0;i<litcode.size();++i) literal[i]=(uint8_t)M[litcode[i]&3]; }

    // ---- references --------------------------------------------------------
    auto gaps = varints(dec("mem_dstgap"));
    auto lraw = varints(dec("mem_len"));
    auto rcb  = dec("mem_rc");
    const size_t NR = gaps.size();
    std::vector<uint32_t> dst(NR), mlen(NR);
    { uint64_t prev=0;
      for(size_t i=0;i<NR;++i){ uint64_t d=prev+gaps[i]; dst[i]=(uint32_t)d;
          mlen[i]=(uint32_t)(lraw[i]+MINMEM); prev=d+mlen[i]; } }
    // mem_self marks references whose source lives in the second region; their
    // sources were coded relative to main_pg_end. Absent stream == none.
    auto selfflags = has("mem_self") ? dec("mem_self") : std::vector<uint8_t>();
    auto src = refc::decode(S["mem_triples"].data(), S["mem_triples"].size(), dst, MAINEND, selfflags);

    // ── CLAIM 3 / coverage — hoisted ABOVE the pseudogenome rebuild ─────────
    // Per-base depth needs only the pseudogenome LENGTH (already in the header)
    // and the per-read placements. It does NOT need pg CONTENT, so decoding the
    // literal stream and replaying every reference is pure waste for this
    // operation. Doing it here skips both.
    if(mode=="coverage"){
        // INDEXING: pos_abs is per-UNIQUE read, read_lengths is per-ORIGINAL
        // read (see 106_inprocess.cpp's "read_lengths is indexed by ORIGINAL
        // read, always" invariant). Indexing both with one counter silently
        // undercounts every DUPLICATE read, since duplicates share a unique's
        // placement but are separate reads for coverage purposes.
        // Measured before this fix on E. coli: 185,346,900 covered bases
        // against the reads' true 232,988,850 -- a 20% undercount.
        // Expand originals through orig2uid, exactly as the read
        // reconstruction path does.
        auto posb2=dec("pos_abs"), lenb2=dec("read_lengths",2);
        auto o2f2=dec("orig2uid_flags"), o2v2=dec("orig2uid_vals");
        std::vector<uint32_t> P(posb2.size()/4);
        memcpy(P.data(),posb2.data(),P.size()*4);
        std::vector<uint16_t> L(lenb2.size()/2);
        memcpy(L.data(),lenb2.data(),L.size()*2);
        const size_t NO2 = L.size();
        std::vector<uint32_t> o2u2;
        if(!o2v2.empty()||!o2f2.empty()){
            size_t k=0;
            for(size_t i=0;i<NO2;++i){
                bool nz = (i>>3)<o2f2.size() ? ((o2f2[i>>3]>>(7-(i&7)))&1) : 0;
                int32_t d=0; if(nz && k+4<=o2v2.size()){ memcpy(&d,&o2v2[k],4); k+=4; }
                o2u2.push_back((uint32_t)d);
            }
            uint32_t exp=0;
            for(auto& v:o2u2){ int32_t d=(int32_t)v; if(d==0){ v=exp; ++exp; } else v=exp-d; }
        } else {
            o2u2.resize(NO2); for(size_t i=0;i<NO2;++i) o2u2[i]=(uint32_t)i;
        }
        std::vector<int32_t> diff(PGLEN+2,0);
        size_t placed=0;
        for(size_t o=0;o<NO2;++o){
            uint32_t u = o<o2u2.size()? o2u2[o] : (uint32_t)o;
            if(u>=P.size()) continue;
            uint64_t a=P[u]; uint16_t l=L[o];          // length from the ORIGINAL read
            if(!l||a>=PGLEN) continue;
            uint64_t b=std::min<uint64_t>(PGLEN,a+l);
            ++diff[a]; --diff[b]; ++placed;
        }
        FILE* f=fopen(outdir.c_str(),"wb");
        if(!f){ fprintf(stderr,"cannot write %s\n",outdir.c_str()); return 1; }
        fprintf(f,"#region\tstart\tend\tdepth\n");
        long cur=0, run=0; uint64_t rs=0;
        for(uint64_t i=0;i<PGLEN;++i){
            cur+=diff[i];
            if(i==0){ run=cur; rs=0; continue; }
            // Force a run break at the main/second region boundary: without it
            // a run spanning MAINEND is emitted as one row labelled by its
            // START, so part of it is attributed to the wrong region.
            if(cur!=run || i==MAINEND){
                fprintf(f,"%s\t%llu\t%llu\t%ld\n", rs<MAINEND?"pg_main":"pg_second",
                        (unsigned long long)rs,(unsigned long long)i,run);
                run=cur; rs=i;
            }
        }
        fprintf(f,"%s\t%llu\t%llu\t%ld\n", rs<MAINEND?"pg_main":"pg_second",
                (unsigned long long)rs,(unsigned long long)PGLEN,run);
        fclose(f);
        fprintf(stderr,"[coverage] %zu placements over %llu bp -> %s (no pg rebuild)\n",
                placed,(unsigned long long)PGLEN,outdir.c_str());
        return 0;
    }

    // ---- rebuild the pseudogenome -----------------------------------------
    // Extension mismatches: a reference may carry positions where the plain
    // copy is WRONG relative to the true content. Two passes, because the
    // "ref" side of the mismatch coder (mmc) is "what the plain copy places",
    // which is only known once that copy has actually happened -- exactly
    // mirroring the encoder, which reads it from `pg` at final emission time.
    // Pass 1: copy everything as before, AND read off ref_bytes right after
    // each ref's own copy (before any override). Pass 2: decode obs against
    // those ref_bytes in one call (the coder is adaptive and order-dependent),
    // then apply each override.
    std::vector<uint8_t> pg(PGLEN,0);
    auto extmmcnt = dec("mem_extmm_cnt");
    auto extmmposraw = has("mem_extmm_pos") ? dec("mem_extmm_pos") : std::vector<uint8_t>();
    auto extmmpos = varints(extmmposraw);
    { size_t li=0, pos=0, mmi=0, applied=0;
      const char CB[256]={0}; (void)CB;
      auto comp=[](uint8_t b)->uint8_t{ return b=='A'?'T':b=='C'?'G':b=='G'?'C':'A'; };
      // A later match's SOURCE can fall inside an EARLIER match's DESTINATION
      // -- exactly what happens through a tandem repeat, where these
      // mismatches concentrate. A batch decode-then-apply (compute every ref
      // byte, decode all obs in one call, apply afterwards) reads the earlier
      // match's un-corrected byte in that case: verified directly on
      // Drosophila SRR40104920, ref i=7 (dst=15789): encoder ref=G obs=C,
      // batch-decoder ref=C obs=T -- same entry, wrong ref because an earlier
      // ref's own override at that source position hadn't been applied yet.
      // So each override is applied immediately, before the NEXT ref's copy
      // can read it, and the range coder is driven one symbol at a time
      // through the SAME adaptive state `decode()` uses (mmc::StreamDecoder).
      mmc::StreamDecoder msd;
      if(has("mem_extmm_obs")) msd.init(S["mem_extmm_obs"].data(), S["mem_extmm_obs"].size());
      for(size_t i=0;i<NR;++i){
          const size_t d=dst[i], L=mlen[i];
          if(d>pos){ memcpy(&pg[pos],&literal[li],d-pos); li+=d-pos; }
          if(rcb.size()>i && rcb[i]){ for(size_t k=0;k<L;++k) pg[d+k]=comp(pg[src[i]+L-1-k]); }
          else memcpy(&pg[d],&pg[src[i]],L);
          const uint8_t cnt = (i<extmmcnt.size())?extmmcnt[i]:0;
          for(uint8_t k=0;k<cnt && mmi<extmmpos.size();++k,++mmi){
              const uint32_t off=(uint32_t)extmmpos[mmi];
              if(off>=L) continue;
              const uint8_t refb = pg[d+off];        // what the copy just placed, NOW, fully up to date
              pg[d+off] = msd.next(refb);
              ++applied;
          }
          pos=d+L;
      }
      if(PGLEN>pos){ memcpy(&pg[pos],&literal[li],PGLEN-pos); li+=PGLEN-pos; }
      if(li!=literal.size()){ fprintf(stderr,"literal not fully consumed: %zu vs %zu\n",li,literal.size()); return 1; }
      if(applied) fprintf(stderr,"  extension mismatches applied: %zu\n", applied);
    }
    fprintf(stderr,"  pg rebuilt: %llu bytes from %zu refs\n",(unsigned long long)PGLEN,NR);
    _rss("after read_capsule / pg built");
    if(getenv("DUMP_PG")){
        FILE* f=fopen("pg_full_dec.txt","wb"); fwrite(pg.data(),1,pg.size(),f); fclose(f);
    }

    _rss("after quality");
    // ── CONTIGS HANDED BACK IN MEMORY ───────────────────────────────────────
    // The caller used to get its contigs by invoking this whole function a
    // SECOND time in "export" mode: re-open the archive, re-decode the ref and
    // literal streams, re-apply the extension mismatches, rebuild the same
    // ~148 MB pseudogenome, write ~450k FASTA records to disk, and then parse
    // that file back into strings. 5.60 s of export plus 0.97 s of re-parse,
    // for a pseudogenome that is sitting in `pg` right here.
    //
    // Same contigs, character for character: this is the export loop's own
    // span walk, including its skip of empty or out-of-range spans, taking
    // substrings of pg instead of writing 60-column FASTA that the caller
    // would immediately re-concatenate.
    if(out_contigs && has("contig_spans")){
        auto sb = dec("contig_spans");
        size_t p2 = 0;
        auto getv = [&]() -> uint64_t {                 // LEB128
            uint64_t x = 0; int sh = 0;
            while (p2 < sb.size()) { uint8_t b = sb[p2++];
                x |= (uint64_t)(b & 0x7F) << sh;
                if (!(b & 0x80)) break; sh += 7; }
            return x;
        };
        const uint64_t nsp = getv();
        uint64_t prev = 0;
        out_contigs->reserve((size_t)nsp);
        for (uint64_t i = 0; i < nsp; ++i) {
            const uint64_t gap = getv(), len = getv();
            const uint64_t a0 = prev + gap, b0 = a0 + len;
            prev = b0;
            if (b0 > pg.size() || len == 0) continue;
            out_contigs->emplace_back((const char*)pg.data() + a0, (size_t)len);
        }
        fprintf(stderr, "[export] %zu contigs (from contig_spans) handed over in memory\n",
                out_contigs->size());
    }

    // ── CLAIM 3 / export ────────────────────────────────────────────────────
    // The pseudogenome IS the assembly; it was built at compress time. Emit it
    // as FASTA and stop -- no read reconstruction, no assembler.
    if(mode=="export"){
        FILE* f=fopen(outdir.c_str(),"wb");
        if(!f){ fprintf(stderr,"cannot write %s\n",outdir.c_str()); return 1; }
        // main region and second region as separate records
        auto emit=[&](const char* nm,size_t a,size_t b){
            if(b<=a) return;
            fprintf(f,">%s len=%zu\n",nm,b-a);
            for(size_t i=a;i<b;i+=60){ size_t e=std::min(b,i+60);
                fwrite(pg.data()+i,1,e-i,f); fputc('\n',f); }
        };
        // PER-CONTIG EXPORT when the archive carries contig_spans.
        //
        // The two-record form (pg_main / pg_second) is the pseudogenome as one
        // concatenated sequence. That is right for Claim 3's `export`, but the
        // CALLER needs the individual assembled contigs -- `build_substrate`
        // collapses and re-places reads per contig, and handing it 2 giant
        // records is a completely different operation from handing it the
        // ~451k real ones.
        //
        // contig_spans is written by the encoder under CAPS_CALL (see
        // 106_inprocess.cpp): the chain boundaries, which the decoder cannot
        // recompute because it rebuilds pg faithfully but never sees where one
        // chain ended and the next began.
        if (has("contig_spans") && getenv("CAPSULE_EXPORT_CONTIGS")) {
            auto sb = dec("contig_spans");
            size_t p2 = 0;
            auto getv = [&]() -> uint64_t {                 // LEB128
                uint64_t x = 0; int sh = 0;
                while (p2 < sb.size()) { uint8_t b = sb[p2++];
                    x |= (uint64_t)(b & 0x7F) << sh;
                    if (!(b & 0x80)) break; sh += 7; }
                return x;
            };
            const uint64_t nsp = getv();
            uint64_t prev = 0; size_t nw = 0;
            for (uint64_t i = 0; i < nsp; ++i) {
                const uint64_t gap = getv(), len = getv();
                const uint64_t a0 = prev + gap, b0 = a0 + len;
                prev = b0;
                if (b0 > pg.size() || len == 0) continue;
                fprintf(f, ">contig_%llu\n", (unsigned long long)i);
                for (size_t k = (size_t)a0; k < (size_t)b0; k += 60) {
                    size_t e = std::min((size_t)b0, k + 60);
                    fwrite(pg.data() + k, 1, e - k, f); fputc('\n', f);
                }
                ++nw;
            }
            fclose(f);
            fprintf(stderr, "[export] %zu contigs (from contig_spans) -> %s\n", nw, outdir.c_str());
            return 0;
        }
        emit("capsule_pg_main",0,(size_t)MAINEND);
        emit("capsule_pg_second",(size_t)MAINEND,pg.size());
        fclose(f);
        fprintf(stderr,"[export] %zu bp pseudogenome -> %s\n",pg.size(),outdir.c_str());
        return 0;
    }

    // ---- per-read streams --------------------------------------------------
    auto posb=dec("pos_abs"), lenb=dec("read_lengths",2), strb=dec("pos_strand");
    auto o2f=dec("orig2uid_flags"), o2v=dec("orig2uid_vals");
    auto cf=dec("mm_cnt_flags"), cv=dec("mm_cnt_vals"), cflat=dec("mm_cnt");

    // ── CLAIM 3 / coverage and query ────────────────────────────────────────
    // Both need only the per-read PLACEMENTS, which the compressor computed and
    // stored. No alignment, no index build, no read reconstruction.
    if(mode=="coverage" || mode=="query"){
        std::vector<uint32_t> P(posb.size()/4);
        memcpy(P.data(),posb.data(),P.size()*4);
        std::vector<uint16_t> L(lenb.size()/2);
        memcpy(L.data(),lenb.data(),L.size()*2);
        // orig2uid: 1 bit/read "is a duplicate" + sparse alias values
        std::vector<uint32_t> vals; { auto vv=varints(o2v); vals.assign(vv.begin(),vv.end()); }
        const size_t NORIG = L.size();
        auto uid_of=[&](size_t o)->uint32_t{
            if(o>=NORIG) return UINT32_MAX;
            size_t byte=o>>3, bit=o&7;
            bool dup = (byte<o2f.size()) && ((o2f[byte]>>bit)&1);
            if(!dup) return (uint32_t)o < (uint32_t)P.size() ? (uint32_t)o : UINT32_MAX;
            return UINT32_MAX;   // duplicate: same placement as its representative
        };
        if(mode=="coverage"){
            std::vector<uint32_t> depth(pg.size()+1,0);
            size_t placed=0;
            for(size_t u=0;u<P.size();++u){
                uint64_t a=P[u]; uint16_t l = u<L.size()?L[u]:0;
                if(!l || a==UINT32_MAX || a>=pg.size()) continue;
                uint64_t b=std::min<uint64_t>(pg.size(),a+l);
                ++depth[a]; if(b<depth.size()) --depth[b];       // difference array
                ++placed;
            }
            FILE* f=fopen(outdir.c_str(),"wb");
            if(!f){ fprintf(stderr,"cannot write %s\n",outdir.c_str()); return 1; }
            fprintf(f,"#region\tstart\tend\tdepth\n");
            long run=0; uint64_t runstart=0; long cur=0;
            for(size_t i=0;i<pg.size();++i){
                cur+=(long)depth[i];
                if(i==0){ run=cur; runstart=0; continue; }
                if(cur!=run){
                    const char* rg = runstart<MAINEND ? "pg_main":"pg_second";
                    fprintf(f,"%s\t%llu\t%zu\t%ld\n",rg,(unsigned long long)runstart,i,run);
                    run=cur; runstart=i;
                }
            }
            { const char* rg = runstart<MAINEND ? "pg_main":"pg_second";
              fprintf(f,"%s\t%llu\t%zu\t%ld\n",rg,(unsigned long long)runstart,pg.size(),run); }
            fclose(f);
            fprintf(stderr,"[coverage] %zu placements over %zu bp -> %s\n",placed,pg.size(),outdir.c_str());
            return 0;
        }
        // query: "START-END" -> every read overlapping that pseudogenome range
        uint64_t qa=0,qb=0;
        { const char* d=strchr(modearg.c_str(),'-');
          if(!d){ fprintf(stderr,"query needs START-END\n"); return 2; }
          qa=strtoull(modearg.c_str(),nullptr,10); qb=strtoull(d+1,nullptr,10); }
        FILE* f=fopen(outdir.c_str(),"wb");
        if(!f){ fprintf(stderr,"cannot write %s\n",outdir.c_str()); return 1; }
        // Emit one record per UNIQUE read. A duplicate read has the same
        // placement AND the same sequence, so emitting it again would return
        // byte-identical records; callers wanting original multiplicity should
        // use the coverage output, which does count every original read.
        size_t n=0;
        for(size_t u=0;u<P.size();++u){
            uint64_t a=P[u]; uint16_t l=u<L.size()?L[u]:0;
            if(!l||a==UINT32_MAX||a>=pg.size()) continue;
            uint64_t b=a+l;
            if(b<=qa||a>=qb) continue;                       // no overlap
            uint64_t e=std::min<uint64_t>(pg.size(),b);
            fprintf(f,">r%zu pos=%llu len=%llu\n",u,(unsigned long long)a,(unsigned long long)(e-a));
            fwrite(pg.data()+a,1,e-a,f); fputc('\n',f); ++n;
        }
        fclose(f);
        fprintf(stderr,"[query] %zu reads overlap %llu-%llu -> %s\n",n,
                (unsigned long long)qa,(unsigned long long)qb,outdir.c_str());
        return 0;
    }
    std::vector<uint32_t> positions(posb.size()/4);
    memcpy(positions.data(),posb.data(),positions.size()*4);
    std::vector<uint16_t> lengths(lenb.size()/2);
    memcpy(lengths.data(),lenb.data(),lengths.size()*2);
    const size_t NU=positions.size();
    std::vector<uint8_t> strand(NU,0);
    for(size_t i=0;i<NU;++i){ size_t B=i>>3,b=i&7; if(B<strb.size()) strand[i]=(strb[B]>>(7-b))&1; }
    // orig2uid: delta coded
    std::vector<uint32_t> o2u;
    { std::vector<uint8_t> raw;
      if(!o2v.empty()||!o2f.empty()){
          size_t k=0; const size_t NO=lengths.size();
          for(size_t i=0;i<NO;++i){
              bool nz = (i>>3)<o2f.size() ? ((o2f[i>>3]>>(7-(i&7)))&1) : 0;
              int32_t d=0; if(nz && k+4<=o2v.size()){ memcpy(&d,&o2v[k],4); k+=4; }
              o2u.push_back((uint32_t)d);
          }
      } }
    { uint32_t exp=0; for(auto& v:o2u){ int32_t d=(int32_t)v; if(d==0){ v=exp; ++exp; } else v=exp-d; } }
    // mm counts
    std::vector<uint16_t> mmcount;
    if(!cf.empty()||!cv.empty()) mmcount = mmcnt_join(cf,cv,NU);
    else { mmcount.resize(cflat.size()/2); memcpy(mmcount.data(),cflat.data(),mmcount.size()*2); }
    // mm positions (flat or bucketed)
    std::vector<uint8_t> mmpos;
    { auto& raw=S["mm_pos"];
      if(!raw.empty()){
          const uint8_t form=raw[0];
          std::vector<uint8_t> body(raw.begin()+1, raw.end());
          mmpos = form ? mmpos_decode_buckets(body.data(),body.size(),mmcount)
                       : capsule_decode_stream(body,1);
      } }
    // Positions are one byte per mismatch while Lmax<=256, and varint-coded
    // deltas above that -- a byte cannot hold a position past 255, which is
    // what used to make long-read archives lossy. Both forms are widened to
    // uint32 here so the two consumers below index one array either way.
    std::vector<uint32_t> mmpos32;
    {
        uint16_t Lm=0; for(auto L:lengths) if(L>Lm) Lm=L;
        if(Lm<=256){ mmpos32.assign(mmpos.begin(), mmpos.end()); }
        else {
            const uint8_t* p=mmpos.data(); const uint8_t* e=p+mmpos.size();
            while(p<e){ uint32_t v=0; int sh=0;
                        while(p<e){ uint8_t b=*p++; v |= (uint32_t)(b&0x7f)<<sh;
                                    if(!(b&0x80)) break; sh+=7; }
                        mmpos32.push_back(v); }
        }
    }

    // ---- derive ref, then decode obs --------------------------------------
    uint16_t Lmax=0; for(auto L:lengths) if(L>Lmax) Lmax=L;
    const bool MMDELTA = (Lmax<=256);
    std::vector<uint16_t> rlenU(NU,0);
    // The encoder indexes mismatches by the UNIQUE read's own length. Under
    // containment a shorter read is absorbed into a longer one, so several
    // ORIGINAL reads of DIFFERENT lengths can share a unique id -- and taking
    // whichever happened to come first gave the wrong length, which shifts the
    // reverse-complement index (q+RL-1-j) and so hands the adaptive mismatch
    // decoder the wrong `ref` byte. From there the model diverges and the
    // stream decodes to wrong symbols: 3,712 of 100,000 reads on C. jejuni,
    // most of them wrong in exactly one base. Containment only arises with
    // variable-length input, which is why fixed-length datasets never showed
    // it. The unique read is the LONGEST of its group, so take the max.
    for(size_t o=0;o<o2u.size() && o<lengths.size();++o)
        if(o2u[o]<NU && lengths[o]>rlenU[o2u[o]]) rlenU[o2u[o]]=lengths[o];
    std::vector<uint8_t> refs; refs.reserve(mmpos.size());
    { size_t off=0;
      auto comp=[](uint8_t b)->uint8_t{ return b=='A'?'T':b=='C'?'G':b=='G'?'C':'A'; };
      for(size_t u=0;u<NU;++u){
          const uint16_t cnt = u<mmcount.size()?mmcount[u]:0;
          if(!cnt) continue;
          const int64_t RL=rlenU[u]; const bool rc=strand[u];
          // pos_abs already holds the coordinate the mismatch loop needs: the
          // encoder converts RC positions once, when writing the stream
          // (q = PL - ppos - rl), so no further conversion belongs here.
          // Verified directly against the encoder's mm_ref: for the first
          // mismatch-carrying read (u=0, rc=1, pos=1,538,107, j=137, ref='T')
          // only pg[pos + RL-1-j] = pg[1,538,120] gives 'T'. Converting with
          // main_pg_end or the full pg length both give the wrong base.
          const int64_t q = (int64_t)positions[u];
          uint32_t prevj=0;
          for(uint16_t m=0;m<cnt;++m){
              uint32_t j = off+m<mmpos32.size()?mmpos32[off+m]:0;
              j=prevj+j; prevj=j;
              int64_t idx = rc ? q+RL-1-(int64_t)j : q+(int64_t)j;
              refs.push_back((idx>=0 && idx<(int64_t)PGLEN) ? pg[idx] : 'A');
          }
          off+=cnt;
      } }
    _rss("before mm_sym decode");
    auto obs = mmc::decode(S["mm_sym"].data(), S["mm_sym"].size(), refs);
    _rss("after mm_sym decode");

    // ---- reconstruct the reads (was decode_105.py) -------------------------
    // That script was 80% of decompression wall clock and set the memory peak
    // at 1563 MB, above our own compress peak, because it built one Python
    // string per read for millions of reads. Everything it needs is already
    // decoded here.
    if(!outreads.empty()){
        const size_t NO = o2u.size();
        std::vector<size_t> mmoff(NU+1,0);
        for(size_t u=0;u<NU;++u) mmoff[u+1]=mmoff[u]+(u<mmcount.size()?mmcount[u]:0);
        auto comp=[](uint8_t b)->uint8_t{ return b=='A'?'T':b=='C'?'G':b=='G'?'C':'A'; };

        // One flat buffer: every read plus its newline. Read lengths are known,
        // so the layout is computed up front and filled in parallel.
        std::vector<size_t> rowoff(NO+1,0);
        for(size_t o=0;o<NO;++o) rowoff[o+1]=rowoff[o]+(o<lengths.size()?lengths[o]:0)+1;
        _rss("before flat alloc");
        std::vector<uint8_t> flat(rowoff[NO], '\n');
        _rss("after flat alloc");

        std::atomic<size_t> nxo{0};
        unsigned T=std::max(1u,std::thread::hardware_concurrency());
        std::vector<std::thread> th;
        for(unsigned t=0;t<T;++t) th.emplace_back([&]{
            for(;;){
                size_t lo=nxo.fetch_add(4096); if(lo>=NO) break;
                size_t hi=std::min(NO,lo+4096);
                for(size_t o=lo;o<hi;++o){
                    const uint32_t u=o2u[o];
                    const uint32_t L=(o<lengths.size())?lengths[o]:0;
                    if(u>=NU || !L) continue;
                    const uint64_t pp=positions[u];
                    const bool rc=strand[u];
                    uint8_t* dst=&flat[rowoff[o]];
                    // pg[pp .. pp+RLu) holds the reverse complement of the
                    // UNIQUE read, so a contained original -- a PREFIX in read
                    // space -- is a SUFFIX in pg space. Indexing with the
                    // original's own L instead of the unique's RLu therefore
                    // read the wrong end for every contained reverse-strand
                    // read: 95 of 100,000 on C. jejuni, all of them RC and all
                    // shorter than their unique. Containment only occurs with
                    // variable-length input, which is why fixed-length data
                    // never showed it. Forward reads are prefixes in both
                    // spaces and are unaffected.
                    const uint32_t RLu = (u<rlenU.size() && rlenU[u]>L) ? rlenU[u] : L;
                    if(!rc){ for(uint32_t k=0;k<L;++k) dst[k]=(pp+k<PGLEN)?pg[pp+k]:'A'; }
                    else   { for(uint32_t k=0;k<L;++k) dst[k]=(pp+RLu-1-k<PGLEN)?comp(pg[pp+RLu-1-k]):'A'; }
                    const uint16_t cnt=(u<mmcount.size())?mmcount[u]:0;
                    if(!cnt) continue;
                    size_t off=mmoff[u]; uint32_t prevj=0;
                    for(uint16_t m=0;m<cnt;++m){
                        if(off+m>=mmpos32.size()||off+m>=obs.size()) break;
                        uint32_t j = off+m<mmpos32.size()?mmpos32[off+m]:0;
                        j=prevj+j; prevj=j;
                        // The old `j==255 -> skip` guard is gone with the clamp
                        // that created it: 255 is now an ordinary delta value.
                        if(j>=L) continue;             // mismatch past this read
                        const uint8_t ob=obs[off+m];
                        dst[j]= rc ? comp(ob) : ob;
                    }
                }
            }
        });
        for(auto& x:th) x.join();

        // N restoration: N-reads went through the same pipeline with each N
        // replaced by 'A', so only the characters need putting back.
        { auto ni=dec("n_indices"), nc=dec("n_cnt"), np=dec("n_pos");
          const size_t NI=ni.size()/4; size_t k=0;
          for(size_t r=0;r<NI && r<nc.size();++r){
              uint32_t oi; memcpy(&oi,&ni[r*4],4);
              const uint8_t c=nc[r];
              if(oi<NO) for(uint8_t m=0;m<c && k+m<np.size();++m){
                  const size_t j=np[k+m];
                  if(rowoff[oi]+j < rowoff[oi+1]-1) flat[rowoff[oi]+j]='N';
              }
              k+=c;
          } }

        if(out_flat && out_rowoff){
            out_flat->swap(flat);          // move, no copy
            *out_rowoff = rowoff;
            fprintf(stderr,"  reads handed over in memory: %zu (no file written)\n", NO);
            fprintf(stderr,"  [dec-timing] read reconstruction %7.2fs\n",
                    std::chrono::duration<double>(std::chrono::steady_clock::now()-_dt0).count());
            _dt0 = std::chrono::steady_clock::now();
        } else {
            FILE* of=fopen(outreads.c_str(),"wb");
            if(of){ fwrite(flat.data(),1,flat.size(),of); fclose(of); }
            fprintf(stderr,"  reads written: %zu\n", NO);
        }
    }

    // ---- names / read-ID column (Phase 2) ----------------------------------
    // Present only when the archive was written with CAPS_NAMES=1; absent
    // archives decode exactly as before. Streams block-by-block straight to
    // the output file so the decoder keeps the encoder's bounded-memory
    // property instead of materializing every name.
    if(has("names_body")){
        const std::string npath = outreads.empty() ? (outdir+"/names.txt") : (outreads+".names");
        FILE* nf=fopen(npath.c_str(),"wb");
        if(nf){
            auto ndict = dec("names_dict");
            auto nindex= dec("names_index");
            // ID_SEQLEN codes "length=NNN" in the header as a reference to the
            // read's own length, which we have already decoded above into
            // `lengths` (per ORIGINAL read, the same order names are in). Hand
            // it over: without this the names column cannot be reconstructed.
            std::vector<uint32_t> slens(lengths.begin(), lengths.end());
            const uint64_t nw = nmc::decode_to_file(S["names_body"], ndict, nindex, nf, false, &slens);
            fclose(nf);
            fprintf(stderr,"  names written: %llu -> %s\n",(unsigned long long)nw,npath.c_str());
        }
    }

    // ---- quality column (Phase 3) ------------------------------------------
    // Present only when the archive was written with CAPS_QUAL=1. Record
    // boundaries come from `lengths` (per ORIGINAL read, already decoded
    // above) rather than from anything fqzcomp stores, so the two columns
    // cannot disagree about where a read ends.
    // Straight to the caller's bitmaps, in parallel, no text and no file. The
    // consumer reduces quality to one bit per base anyway, so decoding 1.86 GB
    // of text and writing it out only to re-read it was the wrong
    // representation as well as the wrong number of cores (1 of 12).
    _rss("before quality");
    if(has("qual_body") && out_qtext && !getenv("CAPS_SKIP_QUAL")){
        // Full-caller route: quality as TEXT, in memory, in parallel. Same
        // characters decode_to_file produced; no 590 MB write and re-read, and
        // 12 cores instead of 1.
        _dt0 = std::chrono::steady_clock::now();
        auto qindex = dec("qual_index");
        std::vector<uint32_t> qlens(lengths.begin(), lengths.end());
        const uint64_t qw = qlc::decode_to_strings(S["qual_body"], qindex, qlens, *out_qtext);
        fprintf(stderr,"  quality -> text in parallel: %llu\n",(unsigned long long)qw);
        fprintf(stderr,"  [dec-timing] quality decode     %7.2fs\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now()-_dt0).count());
    } else if(has("qual_body") && out_qbits && !getenv("CAPS_SKIP_QUAL")){
        _dt0 = std::chrono::steady_clock::now();
        auto qindex = dec("qual_index");
        std::vector<uint32_t> qlens(lengths.begin(), lengths.end());
        const uint64_t qw = qlc::decode_to_bitmaps(S["qual_body"], qindex, qlens,
                                                   qbits_qmin, *out_qbits);
        fprintf(stderr,"  quality -> bitmaps in parallel: %llu\n",(unsigned long long)qw);
        fprintf(stderr,"  [dec-timing] quality decode     %7.2fs\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now()-_dt0).count());
    } else if(has("qual_body") && !getenv("CAPS_SKIP_QUAL")){
        _dt0 = std::chrono::steady_clock::now();
        const std::string qpath = outreads.empty() ? (outdir+"/qual.txt") : (outreads+".qual");
        FILE* qf=fopen(qpath.c_str(),"wb");
        if(qf){
            auto qindex = dec("qual_index");
            std::vector<uint32_t> qlens(lengths.begin(), lengths.end());
            const uint64_t qw = qlc::decode_to_file(S["qual_body"], qindex, qf, qlens);
            fclose(qf);
            fprintf(stderr,"  quality written: %llu -> %s\n",(unsigned long long)qw,qpath.c_str());
            fprintf(stderr,"  [dec-timing] quality decode     %7.2fs\n",
                    std::chrono::duration<double>(std::chrono::steady_clock::now()-_dt0).count());
        }
    }

    // ---- emit what the read-reconstruction step reads ----------------------
    std::string O=outdir;
    put_file(O+"/literal.txt", literal.data(), literal.size());
    { std::vector<uint8_t> tri(NR*13);
      for(size_t i=0;i<NR;++i){ uint8_t* p=&tri[i*13];
          memcpy(p,&dst[i],4); memcpy(p+4,&src[i],4); memcpy(p+8,&mlen[i],4);
          p[12]=(rcb.size()>i)?rcb[i]:0; }
      put_file(O+"/mem_triples.bin",tri.data(),tri.size()); }
    put_file(O+"/pos_abs.bin",posb.data(),posb.size());
    put_file(O+"/pos_strand.bin",strb.data(),strb.size());
    put_file(O+"/read_lengths.bin",lenb.data(),lenb.size());
    { std::vector<int32_t> d(o2u.size()); uint32_t exp=0;
      for(size_t i=0;i<o2u.size();++i){ d[i]= (o2u[i]==exp)?0:(int32_t)(exp-o2u[i]); if(o2u[i]==exp) ++exp; }
      put_file(O+"/orig2uid.bin",d.data(),d.size()*4); }
    put_file(O+"/mm_count_per_read.bin",mmcount.data(),mmcount.size()*2);
    put_file(O+"/mm_ref.bin",refs.data(),refs.size());
    put_file(O+"/mm_obs.bin",obs.data(),obs.size());
    put_file(O+"/mm_pos.bin",mmpos.data(),mmpos.size());
    { auto a=dec("n_pos"),b=dec("n_indices"),c=dec("n_cnt");
      put_file(O+"/n_pos.bin",a.data(),a.size());
      put_file(O+"/n_indices.bin",b.data(),b.size());
      put_file(O+"/n_cnt.bin",c.data(),c.size()); }
    { FILE* f=fopen((O+"/pg_params.txt").c_str(),"w");
      if(f){ fprintf(f,"%llu %llu\n",(unsigned long long)PGLEN,(unsigned long long)MAINEND); fclose(f);} }
    fprintf(stderr,"  streams written to %s\n",O.c_str());
    return 0;
}

// ═══════════════════════════════════════════════════════════════════════════
//  CALL VARIANTS FROM A STORED ARCHIVE
//
//  Until now Claim 2 called at COMPRESS time, from the encoder's in-memory
//  state, so the honest sentence was "compress once and get calls as a
//  byproduct" -- not "open an archive and call variants". This closes that.
//
//  It reuses the decoder's existing modes rather than reaching into them:
//    1. normal decode  -> reads.seq and reads.seq.qual, in ORIGINAL read order
//                         (the same order the caller indexes by)
//    2. export         -> the pseudogenome, which supplies cd.contigs
//    3. pack both with capspack (the SAME definition the encoder uses -- see
//       include/caps_pack.h; two copies of that format is how the two paths
//       would silently disagree)
//    4. run_variant_call
//
//  CAPS_DBG / CAPS_DBG_ONLY are set here because the graph caller keys its
//  input format off them: with both set it expects packed reads and a quality
//  BITMAP, which is exactly what capspack produces. Setting them in the
//  environment and forgetting to pack, or packing and forgetting to set them,
//  are both silent wrong-data failures.
// Decode a single named stream straight from an archive path. The main
// decoder's `dec()` is a lambda over state local to capsule_decode_all, and the
// call path needs a few streams (pos_abs, pos_strand, contig_spans, orig2uid)
// without re-running the whole decode.
static std::vector<uint8_t> arc_stream(const std::string& path, const char* name, size_t w = 1)
{
    uint64_t pg_len = 0, main_end = 0; uint32_t minmem = 0;
    std::vector<Stream> ss;
    if (!read_capsule(path.c_str(), pg_len, main_end, minmem, ss)) return {};
    for (auto& st : ss)
        if (st.name == name) return capsule_decode_stream(st.coded, w);
    return {};
}

static int capsule_call_from_archive(const std::string& in, const std::string& out_vcf,
                                     const std::string& workdir)
{
    const std::string wd = workdir.empty() ? (in + ".calltmp") : workdir;
    mkdir(wd.c_str(), 0755);
    const std::string rp = wd + "/reads.seq";
    const std::string cf = wd + "/contigs.fa";

    // CAPS_CALL_INDELS=1 runs the FULL caller (SNV pileup + indel_pass) from a
    // stored archive instead of the graph-only path. Two reasons:
    //  * indels are unreachable from the archive otherwise -- DBG_ONLY skips
    //    indel_pass entirely, so Option 0 cannot emit them at all;
    //  * it is the fast iteration harness for indel work. The encoder + SNV
    //    stages ahead of indel_pass cost ~650 s of every full run; entering
    //    from the archive skips the encoder and starts the caller in ~18 s.
    // Default is unchanged (graph-only), so Option 0's measured numbers stand.
    const bool WANT_INDELS = std::getenv("CAPS_CALL_INDELS") != nullptr;
    if (!WANT_INDELS) {
        setenv("CAPS_DBG", "1", 1);
        setenv("CAPS_DBG_ONLY", "1", 1);
    } else {
        // CAREFUL: `SEQ_PACKED` in caps_caller.h is
        // `getenv("CAPS_DBG_ONLY") && getenv("CAPS_DBG")`, so clearing
        // DBG_ONLY also switches the caller to expecting PLAIN TEXT reads.
        // Packing them anyway would be read as text -- silent corruption, not
        // an error. So in this mode the reads are handed over UNPACKED, which
        // is what the full caller expects.
        unsetenv("CAPS_DBG");
        unsetenv("CAPS_DBG_ONLY");
    }
    // Phase timing. The 4 steps below were reported as bare progress lines, so
    // "decode+export+pack = 71.95 s" was only knowable by subtracting the
    // caller's own total from the wall clock -- which says nothing about WHICH
    // of them to attack.
    auto _t0 = std::chrono::steady_clock::now();
    auto _lap = [&](const char* what){
        auto now = std::chrono::steady_clock::now();
        fprintf(stderr, "[call-timing] %-22s %7.2fs\n", what,
                std::chrono::duration<double>(now - _t0).count());
        _t0 = now;
    };
    const int QMIN = getenv("CAPS_DBG_MINQ") ? atoi(getenv("CAPS_DBG_MINQ")) : 20;

    fprintf(stderr, "[call] 1/4 decoding reads and quality from the archive\n");
    // Reads come back IN MEMORY (no 1.86 GB write, no re-parse). Quality still
    // goes to a file: qlc::decode_to_file writes straight to a FILE*, and
    // routing it through memory would cost a 1.86 GB copy on top of an already
    // 5.4 GB peak -- so it is mmap'd below instead, which is page-cache backed.
    // CAPS_CALL_NOQUAL=1: skip the quality column entirely. The coherence pass
    // uses quality only to skip low-Q bases when counting per-path read support
    // (`qs` bitmap); with no quality it counts every read, which is what kc's
    // own min1/min2 already do. Decoding quality costs 35.25 s of a 147 s run
    // -- 24% -- so whether that filter is load-bearing is worth one measurement.
    if(getenv("CAPS_CALL_NOQUAL")) setenv("CAPS_SKIP_QUAL","1",1);
    std::vector<uint8_t> rflat; std::vector<size_t> rowoff;
    std::vector<std::string> qbits;
    std::vector<std::string> qtext_mem;
    std::vector<std::string> contigs_mem;
    // In WANT_INDELS mode the full caller reverses `quals[oi]` PER BASE, so it
    // needs quality as TEXT, not as the Q>=QMIN bitmap the graph path uses.
    // Requesting bitmaps here would be silently wrong rather than an error.
    if(capsule_decode_all(in.c_str(), wd, rp, std::string(), std::string(),
                          &rflat, &rowoff, nullptr,
                          WANT_INDELS ? nullptr : &qbits, QMIN,
                          WANT_INDELS ? &qtext_mem : nullptr,
                          WANT_INDELS ? &contigs_mem : nullptr) != 0){
        fprintf(stderr,"[call] decode failed\n"); return 1; }

    _lap("1 decode reads+qual");
    // In indel mode the caller needs the INDIVIDUAL assembled contigs, not the
    // two concatenated pseudogenome records -- build_substrate collapses and
    // re-places reads per contig, so 2 giant records is a different operation.
    if (WANT_INDELS) setenv("CAPSULE_EXPORT_CONTIGS", "1", 1);
    // Step 1 hands the contigs back in memory when the archive carries
    // contig_spans, so the second full archive pass runs only as a fallback.
    if (contigs_mem.empty()) {
        fprintf(stderr, "[call] 2/4 exporting the retained pseudogenome for contigs\n");
        if(capsule_decode_all(in.c_str(), cf, std::string(), "export", std::string()) != 0)
            fprintf(stderr,"[call] export failed -- ploidy gate will see no contigs\n");
    } else {
        fprintf(stderr, "[call] 2/4 contigs came back with step 1 -- no second archive pass\n");
    }

    _lap("2 export pseudogenome");
    fprintf(stderr, "[call] 3/4 packing reads and quality (capspack, shared with the encoder)\n");
    std::vector<std::string> seqs, quals;
    {
        // Reads: straight out of the handed-over buffer. rowoff[o+1]-rowoff[o]
        // includes the trailing newline the decoder wrote, hence the -1.
        // WANT_INDELS: quality went to a file as text (no bitmaps requested).
        std::vector<std::string> qtext;
        if (WANT_INDELS) {
            // Handed over in memory by decode_to_strings; the file round trip
            // it replaced is gone.
            qtext.swap(qtext_mem);
            fprintf(stderr, "  [call] quality as text: %zu records (in memory)\n", qtext.size());
        }
        const bool haveq = WANT_INDELS ? !qtext.empty() : !qbits.empty();
        const size_t NO = rowoff.empty() ? 0 : rowoff.size() - 1;
        // Quality already arrived as packed bitmaps (parallel, straight from
        // the archive), so this loop only packs sequence. Sized up front so
        // every slot is written exactly once and no two threads share one.
        seqs.resize(NO);
        if(WANT_INDELS && qtext.size() >= NO) quals.swap(qtext);
        else if(!WANT_INDELS && haveq && qbits.size() >= NO) quals.swap(qbits);
        else { quals.assign(NO, std::string());
               if(haveq) fprintf(stderr,"  [call] quality count %zu < reads %zu\n", qbits.size(), NO); }
        #pragma omp parallel for schedule(static)
        for(long long o = 0; o < (long long)NO; ++o){
            const size_t b = rowoff[o];
            size_t len = (rowoff[o+1] > b) ? (rowoff[o+1] - b - 1) : 0;
            if(b + len > rflat.size()) len = (b < rflat.size()) ? (rflat.size() - b) : 0;
            seqs[o] = WANT_INDELS
                    ? std::string((const char*)rflat.data() + b, len)   // plain text
                    : capspack::pack_seq((const char*)rflat.data() + b, len);
        }
        std::vector<uint8_t>().swap(rflat);      // release 1.86 GB before calling
        std::vector<size_t>().swap(rowoff);
        std::vector<std::string>().swap(qbits);
    }
    _lap("3 read back + pack");
    capscall::CallData cd;
    if(!contigs_mem.empty()){
        cd.contigs.swap(contigs_mem);
    } else {   // contigs from the exported pseudogenome; the ploidy gate samples these
        std::ifstream fc(cf);
        std::string line, cur;
        while(std::getline(fc, line)){
            if(!line.empty() && line[0]=='>'){ if(!cur.empty()){ cd.contigs.push_back(cur); cur.clear(); } }
            else { while(!line.empty() && (line.back()=='\n'||line.back()=='\r')) line.pop_back(); cur += line; }
        }
        if(!cur.empty()) cd.contigs.push_back(cur);
    }
    cd.valid = true;
    fprintf(stderr, "[call]     %zu contigs from the archive's pseudogenome\n", cd.contigs.size());

    // ── READ PLACEMENTS FROM THE ARCHIVE (what makes indels possible) ───────
    // The full caller's build_substrate needs each read's (contig id, offset).
    // DBG_ONLY does not -- bubbles read only k-mers -- which is why SNVs work
    // from the archive today and indels do not: with empty placement arrays
    // build_substrate would place NO reads and the indel pass would emit
    // nothing, silently.
    //
    // The archive has the pieces: pos_abs (pg position per UNIQUE read),
    // pos_strand, and contig_spans (pg -> contig boundaries). Two traps, both
    // already documented in this file and both silent if got wrong:
    //   * pos_abs is per-UNIQUE read; the caller indexes per-ORIGINAL read.
    //     They differ by duplicates, and conflating them undercounted coverage
    //     by 20% on E. coli. orig2uid is the expansion.
    //   * read_clip has NO archive source. The encoder derives it during
    //     placement, not from a stream. It is left zero here, which is
    //     CORRECT ONLY IF the caller treats 0 as "no left overhang" -- stated
    //     here so the first VCF comparison against the FASTQ path tells us
    //     whether that assumption holds rather than us assuming it does.
    if (WANT_INDELS) {
        // w=1, matching BOTH proven consumers (the coverage path at ~line 231
        // and the read path at ~line 401): they decode pos_abs as raw bytes and
        // then reinterpret them as uint32_t. Passing w=4 tells the stream
        // decoder a different element width and yields garbage, not an error.
        auto pb  = arc_stream(in, "pos_abs");
        auto sbv = arc_stream(in, "pos_strand");
        auto spb = arc_stream(in, "contig_spans");
        auto ofl = arc_stream(in, "orig2uid_flags");
        auto ovl = arc_stream(in, "orig2uid_vals");
        if (pb.empty() || spb.empty()) {
            fprintf(stderr, "[call] ARCHIVE LACKS %s -- indels need it; "
                            "re-compress with CAPS_CALL=1 so contig_spans is written\n",
                    spb.empty() ? "contig_spans" : "pos_abs");
            return 1;
        }
        // spans -> a sorted table of contig starts
        std::vector<uint64_t> cstart, cend;
        { size_t p2 = 0;
          auto getv=[&]()->uint64_t{ uint64_t x=0; int sh=0;
              while(p2<spb.size()){ uint8_t b=spb[p2++]; x |= (uint64_t)(b&0x7F)<<sh;
                                    if(!(b&0x80)) break; sh+=7; } return x; };
          const uint64_t nsp=getv(); uint64_t prev=0;
          cstart.reserve(nsp); cend.reserve(nsp);
          for(uint64_t i=0;i<nsp;++i){ const uint64_t g=getv(), l=getv();
              const uint64_t a0=prev+g; cstart.push_back(a0); cend.push_back(a0+l); prev=a0+l; } }
        std::vector<uint32_t> P(pb.size()/4);
        memcpy(P.data(), pb.data(), P.size()*4);
        const size_t NU2 = P.size();
        std::vector<uint8_t> st(NU2, 0);
        for(size_t i=0;i<NU2;++i){ size_t B=i>>3,b=i&7; if(B<sbv.size()) st[i]=(sbv[B]>>(7-b))&1; }
        // expand unique -> original exactly as the read path does
        // EXACTLY the expansion the coverage path uses (this file, ~line 239).
        // My first version had the flag sense INVERTED and treated the stored
        // value as an absolute id; it is a DELTA. Flag set -> back-reference to
        // an earlier unique at (exp - d); flag clear -> a new unique at exp.
        // Getting this wrong mis-maps every duplicate read silently -- the same
        // class of error that undercounted coverage by 20% on E. coli.
        std::vector<uint32_t> o2u;
        const size_t NO2 = seqs.size();
        if(!ovl.empty()||!ofl.empty()){
            size_t k=0; o2u.reserve(NO2);
            for(size_t i=0;i<NO2;++i){
                bool nz = (i>>3)<ofl.size() ? ((ofl[i>>3]>>(7-(i&7)))&1) : 0;
                int32_t d=0; if(nz && k+4<=ovl.size()){ memcpy(&d,&ovl[k],4); k+=4; }
                o2u.push_back((uint32_t)d);
            }
            uint32_t exp=0;
            for(auto& v:o2u){ int32_t d=(int32_t)v; if(d==0){ v=exp; ++exp; } else v=(uint32_t)(exp-d); }
        } else { o2u.resize(NO2); for(size_t i=0;i<NO2;++i) o2u[i]=(uint32_t)i; }
        cd.read_cid.assign(NO2, UINT32_MAX);
        cd.read_pos.assign(NO2, 0);
        cd.read_rc.assign(NO2, 0);
        cd.read_clip.assign(NO2, 0);
        size_t placed=0;
        for(size_t o=0;o<NO2;++o){
            const uint32_t u=o2u[o];
            if(u>=NU2) continue;
            const uint64_t gp=P[u];
            auto it=std::upper_bound(cstart.begin(), cstart.end(), gp);
            if(it==cstart.begin()) continue;
            const size_t ci=(size_t)(it-cstart.begin()-1);
            if(gp>=cend[ci]) continue;                 // in a gap between contigs
            cd.read_cid[o]=(uint32_t)ci;
            cd.read_pos[o]=(uint32_t)(gp-cstart[ci]);
            cd.read_rc[o]=st[u];
            ++placed;
        }
        fprintf(stderr, "[call]     %zu/%zu read placements rebuilt from the archive "
                        "(%zu contig spans)\n", placed, NO2, cstart.size());
    }

    _lap("3b parse contigs");
    fprintf(stderr, "[call] 4/4 calling variants\n");
    const int n = capscall::run_variant_call(seqs, quals, cd, out_vcf);
    _lap("4 run_variant_call");
    fprintf(stderr, "[call] %d records -> %s\n", n, out_vcf.c_str());
    return n >= 0 ? 0 : 1;
}

#ifndef CAPSULE_NO_MAIN
int main(int argc,char** argv){
    // Claim 2 from a STORED archive: capsule_decode call <in.capsule> <out.vcf> [workdir]
    if(argc>=4 && !strcmp(argv[1],"call"))
        return capsule_call_from_archive(argv[2], argv[3], argc>4?argv[4]:std::string());
    // Claim 3 modes:  capsule_decode export|coverage|query <in.capsule> <out> [range]
    if(argc>=4 && (!strcmp(argv[1],"export")||!strcmp(argv[1],"coverage")||!strcmp(argv[1],"query")))
        return capsule_decode_all(argv[2], argv[3], std::string(), argv[1],
                                  argc>4?argv[4]:std::string());
    if(argc<3){ fprintf(stderr,"usage: capsule_decode <in.capsule> <outdir> [reads.out]\n"
                               "       capsule_decode export   <in.capsule> <out.fa>\n"
                               "       capsule_decode coverage <in.capsule> <out.tsv>\n"
                               "       capsule_decode query    <in.capsule> <out.fa> <START-END>\n"
                               "       capsule_decode call     <in.capsule> <out.vcf> [workdir]\n"); return 2; }
    return capsule_decode_all(argv[1], argv[2], argc>3?argv[3]:std::string());
}
#endif  // CAPSULE_NO_MAIN
