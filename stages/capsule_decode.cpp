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
#include "seqpar_core.h"
#include "coders_inproc.h"

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

int capsule_decode_all(const char* arcpath, const std::string& outdir,
                       const std::string& outreads = std::string()){
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
    if(getenv("DUMP_PG")){
        FILE* f=fopen("pg_full_dec.txt","wb"); fwrite(pg.data(),1,pg.size(),f); fclose(f);
    }

    // ---- per-read streams --------------------------------------------------
    auto posb=dec("pos_abs"), lenb=dec("read_lengths",2), strb=dec("pos_strand");
    auto o2f=dec("orig2uid_flags"), o2v=dec("orig2uid_vals");
    auto cf=dec("mm_cnt_flags"), cv=dec("mm_cnt_vals"), cflat=dec("mm_cnt");
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

    // ---- derive ref, then decode obs --------------------------------------
    uint16_t Lmax=0; for(auto L:lengths) if(L>Lmax) Lmax=L;
    const bool MMDELTA = (Lmax<=256);
    std::vector<uint16_t> rlenU(NU,0);
    for(size_t o=0;o<o2u.size() && o<lengths.size();++o) if(o2u[o]<NU && !rlenU[o2u[o]]) rlenU[o2u[o]]=lengths[o];
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
              uint32_t j = off+m<mmpos.size()?mmpos[off+m]:0;
              if(MMDELTA){ j=prevj+j; prevj=j; }
              int64_t idx = rc ? q+RL-1-(int64_t)j : q+(int64_t)j;
              refs.push_back((idx>=0 && idx<(int64_t)PGLEN) ? pg[idx] : 'A');
          }
          off+=cnt;
      } }
    auto obs = mmc::decode(S["mm_sym"].data(), S["mm_sym"].size(), refs);

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
        std::vector<uint8_t> flat(rowoff[NO], '\n');

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
                    if(!rc){ for(uint32_t k=0;k<L;++k) dst[k]=(pp+k<PGLEN)?pg[pp+k]:'A'; }
                    else   { for(uint32_t k=0;k<L;++k) dst[k]=(pp+L-1-k<PGLEN)?comp(pg[pp+L-1-k]):'A'; }
                    const uint16_t cnt=(u<mmcount.size())?mmcount[u]:0;
                    if(!cnt) continue;
                    size_t off=mmoff[u]; uint32_t prevj=0;
                    for(uint16_t m=0;m<cnt;++m){
                        if(off+m>=mmpos.size()||off+m>=obs.size()) break;
                        uint32_t j=mmpos[off+m];
                        if(MMDELTA){ j=prevj+j; prevj=j; }
                        else if(j==255) continue;      // capped, >255 bp reads
                        if(j>=L) continue;             // container mismatch past this read
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

        FILE* of=fopen(outreads.c_str(),"wb");
        if(of){ fwrite(flat.data(),1,flat.size(),of); fclose(of); }
        fprintf(stderr,"  reads written: %zu\n", NO);
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

#ifndef CAPSULE_NO_MAIN
int main(int argc,char** argv){
    if(argc<3){ fprintf(stderr,"usage: capsule_decode <in.capsule> <outdir> [reads.out]\n"); return 2; }
    return capsule_decode_all(argv[1], argv[2], argc>3?argv[3]:std::string());
}
#endif  // CAPSULE_NO_MAIN
