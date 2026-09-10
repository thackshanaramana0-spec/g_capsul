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

// Rebuild the full per-read position array from the region-split streams.
// Old archives carry pos_abs alone; new ones carry main positions in pos_abs,
// second-region zigzag-varint deltas in pos_sec, and a region bitmap prefixed
// with the exact count in pos_region. Returns the raw uint32 byte vector the
// rest of the decoder already expects, so nothing downstream changes.
static std::vector<uint8_t> caps_join_positions(const std::vector<uint8_t>& mainb,
                                                const std::vector<uint8_t>& secb,
                                                const std::vector<uint8_t>& regb,
                                                uint64_t MAINEND)
{
    if (regb.size() < 4) return mainb;                 // old format: pos_abs is complete
    uint32_t np = 0; memcpy(&np, regb.data(), 4);
    std::vector<uint8_t> out((size_t)np * 4);
    const uint32_t* mp = reinterpret_cast<const uint32_t*>(mainb.data());
    const size_t nm = mainb.size() / 4;
    size_t mi = 0, sp = 0; int64_t prev = (int64_t)MAINEND;
    for (uint32_t i = 0; i < np; ++i) {
        const bool is_sec = (regb[4 + (i >> 3)] >> (i & 7)) & 1u;
        uint32_t v = 0;
        if (is_sec) {
            uint64_t z = 0; int sh = 0;
            while (sp < secb.size()) { const uint8_t b = secb[sp++];
                z |= (uint64_t)(b & 0x7f) << sh; if (!(b & 0x80)) break; sh += 7; }
            const int64_t d = (int64_t)((z >> 1) ^ (~(z & 1) + 1));
            prev += d; v = (uint32_t)prev;
        } else if (mi < nm) {
            v = mp[mi++];
        }
        memcpy(out.data() + (size_t)i * 4, &v, 4);
    }
    return out;
}

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
    // ── QUALITY DECODED CONCURRENTLY WITH THE SEQUENCE ──────────────────────
    //
    // MEASURED: the literal (sequence) decode runs on 4 threads because the
    // archive carries 4 chunks, and the quality decode on 3 because
    // QBLOCK_BYTES yields 3 blocks. Both counts are fixed by the archive and
    // widening either changes the format, which Claim 1's locked sizes forbid.
    // But the two were also run BACK TO BACK -- 5.1 s then 5.2 s -- with never
    // more than 4 of 12 cores busy.
    //
    // They are independent: quality needs only the per-read lengths, its own
    // small stream, not anything the sequence path produces. So it starts here
    // and is joined where it used to be decoded. Same bytes, same order; the
    // two costs overlap instead of adding.
    std::thread qthread; bool qual_async = false; uint64_t qw_async = 0;
    double qsecs_async = 0.0;
    if (has("qual_body") && (out_qtext || out_qbits) && !getenv("CAPS_SKIP_QUAL")) {
        auto lenb_e = dec("read_lengths", 2);
        std::vector<uint16_t> L16(lenb_e.size()/2);
        if (!L16.empty()) memcpy(L16.data(), lenb_e.data(), L16.size()*2);
        std::vector<uint32_t> qlens_e(L16.begin(), L16.end());
        auto qindex_e = dec("qual_index");
        if (!qlens_e.empty()) {
            const std::vector<uint8_t>& qbody = S.find("qual_body")->second;
            qual_async = true;
            qthread = std::thread([&, qlens_e, qindex_e]() {
                auto t0 = std::chrono::steady_clock::now();
                if (out_qtext) qw_async = qlc::decode_to_strings(qbody, qindex_e, qlens_e, *out_qtext);
                else           qw_async = qlc::decode_to_bitmaps(qbody, qindex_e, qlens_e, qbits_qmin, *out_qbits);
                qsecs_async = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
            });
        }
    }

    // ── COVERAGE EXITS BEFORE THE ASSEMBLY LAYER ────────────────────────────
    //
    // MEASURED: coverage was the SLOWEST of the three Claim 3 modes (0.26 s vs
    // export 0.18 s and query 0.20 s) despite doing the least work -- it needs
    // no pseudogenome and peaks at 4 MB against their 36 MB. The cause was
    // placement: its early exit sat AFTER seq_decode_mem (the context-mixing
    // literal decode, the most expensive stream in the archive), after
    // refc::decode, and after mem_dstgap/mem_len/mem_rc/mem_self -- the entire
    // assembly layer, none of which it reads.
    //
    // It needs only pos_abs, read_lengths, orig2uid, and the scalars PGLEN and
    // MAINEND from the header. Moving the exit above the assembly layer is what
    // makes the claim's own architectural statement -- "coverage is served
    // without reconstructing the assembly" -- true of the code as well as the
    // prose.
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
        posb2 = caps_join_positions(posb2, dec("pos_sec"), dec("pos_region"), MAINEND);
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
    // QUERY BUDGET INSTRUMENT (CAPS_QTIME=1). The point of query is that its
    // cost should track the ANSWER, not the archive. To know whether that is
    // true we need the split, not the total.
    const bool QTIME = getenv("CAPS_QTIME")!=nullptr;
    auto _qt0 = std::chrono::steady_clock::now();
    auto QLAP=[&](const char* what){ if(!QTIME) return;
        auto n=std::chrono::steady_clock::now();
        fprintf(stderr,"  [qtime] %-28s %6.3f s\n", what,
                std::chrono::duration<double>(n-_qt0).count()); _qt0=n; };
    QLAP("archive read + stream map");
    // ── QUERY FROM A SIDECAR INDEX — no assembly layer decoded at all ──────
    // Measured (CAPS_QTIME): 98.6% of a query is spent materialising the whole
    // pseudogenome (literal decode 55%, reference replay 29%, reference stream
    // decode 10%) to hand back a slice. The emit itself is 1.3%.
    //
    // The archive is NOT changed to fix this. Claim 1 is the stronger claim and
    // must not be spent buying Claim 3 a constant. Instead the pseudogenome is
    // written once to an OPTIONAL sidecar, 2-bit packed, exactly the pattern
    // BAM/.bai and CRAM/.crai already establish in this field. Absent the
    // sidecar everything behaves as before; present, query skips the assembly
    // layer entirely and reads bases straight out of it.
    //   build:  capsule_decode index <in.capsule> <out.qidx>
    //   use:    automatic if <in.capsule>.qidx exists, or set CAPS_QIDX
    if(mode=="query"){
        std::string qidx = getenv("CAPS_QIDX") ? getenv("CAPS_QIDX") : (std::string(arcpath) + ".qidx");
        FILE* qf = fopen(qidx.c_str(), "rb");
        if(qf){
            uint64_t magic=0, plen=0, mend=0;
            if(fread(&magic,8,1,qf)==1 && magic==0x5158494451494451ULL &&
               fread(&plen,8,1,qf)==1 && fread(&mend,8,1,qf)==1 && plen==PGLEN){
                std::vector<uint8_t> packed((plen+3)/4);
                if(fread(packed.data(),1,packed.size(),qf)==packed.size()){
                    QLAP("sidecar: pg load");
                    auto base_at=[&](uint64_t i)->char{
                        static const char M[4]={'A','C','G','T'};
                        return M[(packed[i>>2] >> (2*(3-(i&3)))) & 3]; };
                    uint64_t np=0,nl=0; std::vector<uint32_t> P; std::vector<uint16_t> L;
                    if(fread(&np,8,1,qf)==1 && fread(&nl,8,1,qf)==1){
                        P.resize(np); L.resize(nl);
                        if(fread(P.data(),1,np*4,qf)!=np*4 || fread(L.data(),1,nl*2,qf)!=nl*2){
                            fprintf(stderr,"[query] sidecar truncated -- rebuild it\n"); fclose(qf); return 1; }
                    } else { fprintf(stderr,"[query] sidecar has no placements -- rebuild it\n"); fclose(qf); return 1; }
                    // per-read deviations, if this sidecar carries them
                    std::unordered_map<uint32_t,std::pair<uint32_t,uint32_t>> mmix; // u -> (start,count)
                    std::vector<uint32_t> mmoffs; std::vector<uint8_t> mmsyms;
                    { uint64_t nid=0,nmm=0;
                      if(fread(&nid,8,1,qf)==1 && fread(&nmm,8,1,qf)==1 && nid && nmm){
                          std::vector<uint32_t> ids(nid*2);
                          mmoffs.resize(nmm); mmsyms.resize(nmm);
                          if(fread(ids.data(),1,ids.size()*4,qf)==ids.size()*4 &&
                             fread(mmoffs.data(),1,nmm*4,qf)==nmm*4 &&
                             fread(mmsyms.data(),1,nmm,qf)==nmm){
                              uint32_t at=0;
                              for(uint64_t i=0;i<nid;++i){ mmix[ids[i*2]]={at,ids[i*2+1]}; at+=ids[i*2+1]; }
                              fprintf(stderr,"[query] sidecar carries deviations for %llu reads\n",
                                      (unsigned long long)nid);
                          } else { mmix.clear(); mmoffs.clear(); mmsyms.clear(); }
                      } }
                    fclose(qf);
                    QLAP("  placement streams decode");
                    std::vector<std::pair<uint64_t,uint64_t>> rr;
                    bool sq = !modearg.empty() &&
                              modearg.find_first_not_of("ACGTNacgtn")==std::string::npos;
                    if(sq){
                        std::string q=modearg; for(auto&c:q) c=(char)toupper((unsigned char)c);
                        std::string rc(q.rbegin(),q.rend());
                        for(auto&c:rc) c=(c=='A'?'T':c=='T'?'A':c=='C'?'G':c=='G'?'C':c);
                        std::string hay(plen,'\0');
                        for(uint64_t i=0;i<plen;++i) hay[i]=base_at(i);
                        for(const std::string* pat : { &q, &rc }){
                            if(pat==&rc && rc==q) break;
                            for(size_t at=hay.find(*pat); at!=std::string::npos; at=hay.find(*pat,at+1))
                                rr.push_back({(uint64_t)at,(uint64_t)(at+pat->size())});
                        }
                        std::sort(rr.begin(),rr.end());
                        fprintf(stderr,"[query] sequence of %zu bp -> %zu occurrence(s)\n",q.size(),rr.size());
                    } else {
                        const char* d=strchr(modearg.c_str(),'-');
                        if(!d){ fprintf(stderr,"query needs START-END or a DNA sequence\n"); return 2; }
                        rr.push_back({strtoull(modearg.c_str(),nullptr,10), strtoull(d+1,nullptr,10)});
                    }
                    QLAP("  locate ranges");
                    FILE* of=fopen(outdir.c_str(),"wb");
                    if(!of){ fprintf(stderr,"cannot write %s\n",outdir.c_str()); return 1; }
                    size_t n=0; std::string obuf; obuf.reserve(1024);
                    setvbuf(of,nullptr,_IOFBF,1<<20);
                    for(size_t u=0;u<P.size();++u){
                        uint64_t aa=P[u]; uint16_t l=u<L.size()?L[u]:0;
                        if(!l||aa==UINT32_MAX||aa>=plen) continue;
                        uint64_t b=aa+l; bool hit=false;
                        for(const auto& r:rr) if(b>r.first && aa<r.second){ hit=true; break; }
                        if(!hit) continue;
                        uint64_t e=std::min<uint64_t>(plen,b);
                        char hdr[96];
                        int hl=snprintf(hdr,sizeof hdr,">r%zu pos=%llu len=%llu\n",u,
                                        (unsigned long long)aa,(unsigned long long)(e-aa));
                        obuf.resize(0); obuf.append(hdr,hl);
                        const size_t base = obuf.size();
                        for(uint64_t i=aa;i<e;++i) obuf.push_back(base_at(i));
                        // apply this read's own deviations, so what is emitted
                        // is the READ and not the consensus under it
                        auto it = mmix.find((uint32_t)u);
                        if(it != mmix.end()){
                            for(uint32_t k=0;k<it->second.second;++k){
                                const uint32_t o = mmoffs[it->second.first+k];
                                if(base+o < obuf.size()) obuf[base+o] = (char)mmsyms[it->second.first+k];
                            }
                        }
                        obuf.push_back('\n');
                        fwrite(obuf.data(),1,obuf.size(),of); ++n;
                    }
                    fclose(of);
                    QLAP("emit answer from sidecar");
                    fprintf(stderr,"[query] %zu reads over %zu range(s) -> %s (via sidecar, no assembly decode)\n",
                            n,rr.size(),outdir.c_str());
                    return 0;
                }
            }
            fclose(qf);
            fprintf(stderr,"[query] sidecar %s unusable -- falling back to full rebuild\n",qidx.c_str());
        }
    }

    // ---- literal: 2-bit codes -> ACGT --------------------------------------
    auto litcode = seq_decode_mem(S["literal"].data(), S["literal"].size());
    std::vector<uint8_t> literal(litcode.size());
    { const char M[4]={'A','C','G','T'};
      for(size_t i=0;i<litcode.size();++i) literal[i]=(uint8_t)M[litcode[i]&3]; }
    QLAP("literal decode (seq_decode_mem)");

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
    QLAP("reference streams decode");
    // CAPS_DUMP_REFS=1 -> the reference list, so the DEPENDENCY CLOSURE of a
    // pg window can be measured offline. The question windowed rebuild turns
    // on is not whether references reach far back (they do, 99.7% over 100kb)
    // but how much of pg a window transitively depends on.
    if(getenv("CAPS_DUMP_REFS")){
        FILE* rf=fopen(getenv("CAPS_DUMP_REFS"),"wb");
        if(rf){ fprintf(rf,"dst\tsrc\tlen\n");
            for(size_t i=0;i<NR;++i)
                fprintf(rf,"%u\t%u\t%u\n",dst[i],(unsigned)src[i],mlen[i]);
            fclose(rf);
            fprintf(stderr,"[refs] %zu references dumped\n",NR); }
    }

    // ── CLAIM 3 / coverage — hoisted ABOVE the pseudogenome rebuild ─────────
    // Per-base depth needs only the pseudogenome LENGTH (already in the header)
    // and the per-read placements. It does NOT need pg CONTENT, so decoding the
    // literal stream and replaying every reference is pure waste for this
    // operation. Doing it here skips both.


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
    if(out_contigs && !has("contig_spans")){
        // no per-contig spans in this archive: nothing to hand over, the
        // caller falls back to the export pass exactly as before.
    } else if(out_contigs && !getenv("CAPSULE_EXPORT_CONTIGS")){
        // TWO-RECORD FORM, for the graph (CAPS_DBG_ONLY) path.
        //
        // That path does not want the individual contigs -- it wants the
        // pseudogenome as the two concatenated records `export` writes, which
        // is what its ploidy gate samples. Handing it 451k separate contigs
        // instead would be a different input, not a faster route to the same
        // one. So this mirrors emit("capsule_pg_main", 0, MAINEND) and
        // emit("capsule_pg_second", MAINEND, pg.size()) exactly, including
        // their `if (b <= a) return` skip, minus the FASTA round trip.
        if ((size_t)MAINEND > 0)
            out_contigs->emplace_back((const char*)pg.data(), (size_t)MAINEND);
        if (pg.size() > (size_t)MAINEND)
            out_contigs->emplace_back((const char*)pg.data() + MAINEND, pg.size() - (size_t)MAINEND);
        fprintf(stderr, "[export] %zu pseudogenome records handed over in memory\n",
                out_contigs->size());
    } else if(out_contigs){
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
    // ── build the sidecar index: the pseudogenome, 2-bit packed ────────────
    // One full rebuild, written once, so every later query skips the assembly
    // layer. Costs an optional file (~9% of the archive here); costs the
    // ARCHIVE nothing, which is the point.
    // CAPS_PILEUP=1 makes `index` continue past the pg dump and also emit the
    // VARIANT-SITE TABLE (see below). Off by default so `index` stays fast.
    // `index` returns early with pg+placements only. With CAPS_PILEUP=1 it
    // continues past the mismatch decode so the sidecar can also carry each
    // read's own deviations -- which is what makes `query` emit the READ
    // rather than the consensus beneath it.
    const bool WANT_PILEUP = (mode=="index") && getenv("CAPS_PILEUP")!=nullptr;
    if(mode=="index"){
        FILE* f=fopen(outdir.c_str(),"wb");
        if(!f){ fprintf(stderr,"cannot write %s\n",outdir.c_str()); return 1; }
        const uint64_t magic=0x5158494451494451ULL, plen=pg.size(), mend=MAINEND;
        fwrite(&magic,8,1,f); fwrite(&plen,8,1,f); fwrite(&mend,8,1,f);
        std::vector<uint8_t> packed((plen+3)/4, 0);
        for(uint64_t i=0;i<plen;++i){
            unsigned c = pg[i]=='C'?1u : pg[i]=='G'?2u : pg[i]=='T'?3u : 0u;
            packed[i>>2] |= (uint8_t)(c << (2*(3-(i&3))));
        }
        fwrite(packed.data(),1,packed.size(),f);
        // Placements go in too. Without them the sidecar still had to decode
        // pos_abs + read_lengths for ALL reads to answer a query about a few
        // thousand -- measured 0.073 s of a 0.10 s query, i.e. the same
        // decode-everything-to-use-a-slice disease one level down.
        { auto pb=dec("pos_abs"), lb=dec("read_lengths",2);
          if(has("pos_sec")||has("pos_region"))
              pb = caps_join_positions(pb, has("pos_sec")?dec("pos_sec"):std::vector<uint8_t>(),
                                       has("pos_region")?dec("pos_region"):std::vector<uint8_t>(), MAINEND);
          const uint64_t np=pb.size()/4, nl=lb.size()/2;
          fwrite(&np,8,1,f); fwrite(&nl,8,1,f);
          fwrite(pb.data(),1,np*4,f); fwrite(lb.data(),1,nl*2,f);
          fprintf(stderr,"[index] %llu bp + %llu placements -> %zu B sidecar -> %s\n",
                  (unsigned long long)plen,(unsigned long long)np,
                  packed.size()+40+np*4+nl*2, outdir.c_str()); }

        fclose(f);
        // Without CAPS_PILEUP the sidecar is complete here. With it, fall
        // through so the deviations can be appended once the mismatch streams
        // have been decoded further down.
        if(!WANT_PILEUP) return 0;
    }
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
    posb = caps_join_positions(posb, dec("pos_sec"), dec("pos_region"), MAINEND);
    auto o2f=dec("orig2uid_flags"), o2v=dec("orig2uid_vals");
    auto cf=dec("mm_cnt_flags"), cv=dec("mm_cnt_vals"), cflat=dec("mm_cnt");

    // ── CLAIM 3 / coverage and query ────────────────────────────────────────
    // Both need only the per-read PLACEMENTS, which the compressor computed and
    // stored. No alignment, no index build, no read reconstruction.
    QLAP("pg rebuild (replay references)");
    if(mode=="coverage" || mode=="query"){
        std::vector<uint32_t> P(posb.size()/4);
        memcpy(P.data(),posb.data(),P.size()*4);
        std::vector<uint16_t> L(lenb.size()/2);
        memcpy(L.data(),lenb.data(),L.size()*2);
        // orig2uid: 1 bit/read "is a duplicate" + sparse alias values
        std::vector<uint32_t> vals; { auto vv=varints(o2v); vals.assign(vv.begin(),vv.end()); }
        const size_t NORIG = L.size();
        QLAP("placement streams unpack");
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
        // query accepts EITHER "START-END" (pseudogenome offsets) or a literal
        // DNA SEQUENCE.
        //
        // WHY THE SEQUENCE FORM EXISTS. Offsets are pseudogenome offsets, not
        // chromosome positions, so "give me BRCA1" was unanswerable without
        // bringing back a reference -- the one dependency this archive exists
        // to avoid. Asking by SEQUENCE removes the coordinate system from the
        // interface entirely: the user supplies the gene/primer/probe they
        // already have, we locate it in the pseudogenome, and return the reads
        // sitting there. Still no reference, and now biologically meaningful.
        //
        // It is nearly free: pg is already fully decoded and in memory for the
        // offset form, so this adds a scan over a buffer we are holding anyway.
        // Both strands are searched, because a read's placement carries no
        // guarantee about which strand the query was written on.
        std::vector<std::pair<uint64_t,uint64_t>> ranges;
        bool seqmode = !modearg.empty() &&
                       modearg.find_first_not_of("ACGTNacgtn") == std::string::npos;
        if(seqmode){
            std::string q=modearg; for(auto& c:q) c=(char)toupper((unsigned char)c);
            std::string rc(q.rbegin(), q.rend());
            for(auto& c:rc) c = (c=='A'?'T':c=='T'?'A':c=='C'?'G':c=='G'?'C':c);
            const std::string hay((const char*)pg.data(), pg.size());
            for(const std::string* pat : { &q, &rc }){
                if(pat==&rc && rc==q) break;            // palindrome: do not double-count
                for(size_t at=hay.find(*pat); at!=std::string::npos; at=hay.find(*pat,at+1))
                    ranges.push_back({(uint64_t)at,(uint64_t)(at+pat->size())});
            }
            // LONG PROBES FAIL EXACT MATCH, and silently returning nothing is
            // the wrong answer. Reads carry sequencing errors and the
            // pseudogenome is a CONSENSUS, so the chance a long query matches
            // exactly falls off fast: measured on E. coli, a 150 bp read and a
            // 100 bp chunk of it both give 0 occurrences, while 60/40/30 bp
            // chunks of the SAME read all hit. A user querying with a whole
            // read or a gene should get its locus, not an empty file.
            //
            // So on a miss, fall back to sliding a 40 bp window along the
            // query and taking the first window that hits. 40 is not tuned --
            // it is the size at which a probe is long enough to be unique in a
            // human-scale pseudogenome (4^40 >> 3e9) and short enough to have
            // a good chance of being error-free.
            if(ranges.empty() && q.size() > 60){
                // UNION every window that hits, never the first one.
                //
                // Taking the first hit was wrong and measurably so. Windows of
                // one 150 bp read all agree on the LOCUS (every hitting window
                // resolved pg 5,628,008) but return different READ SETS --
                // Jaccard 0.41 between offsets 0 and 60 -- because each window
                // asks about a different 40 bp span. First-hit therefore
                // answers "reads around bases 0-40 of your query", not "reads
                // around your query", and which span you got depended on where
                // the sequencing error happened to fall.
                //
                // Unioning the hitting windows spans the whole query, which is
                // the question actually asked. Windows carrying an error simply
                // contribute nothing.
                const size_t W = 40; size_t hitw = 0, tried = 0;
                for(size_t off=0; off+W<=q.size(); off+=W/2){
                    ++tried;
                    std::string sub=q.substr(off,W), rsub(sub.rbegin(),sub.rend());
                    for(auto& c:rsub) c=(c=='A'?'T':c=='T'?'A':c=='C'?'G':c=='G'?'C':c);
                    const size_t before = ranges.size();
                    for(const std::string* pat : { &sub, &rsub }){
                        if(pat==&rsub && rsub==sub) break;
                        for(size_t at=hay.find(*pat); at!=std::string::npos; at=hay.find(*pat,at+1))
                            ranges.push_back({(uint64_t)at,(uint64_t)(at+pat->size())});
                    }
                    if(ranges.size()>before) ++hitw;
                }
                if(hitw)
                    fprintf(stderr,"[query] %zu bp query had no exact match; %zu of %zu "
                                   "%zu bp windows matched, union taken\n",
                            q.size(), hitw, tried, W);
            }
            std::sort(ranges.begin(),ranges.end());
            fprintf(stderr,"[query] sequence of %zu bp -> %zu occurrence(s) in the pseudogenome\n",
                    q.size(), ranges.size());
            // A probe inside a repeat family resolves to many loci. That is the
            // correct answer -- the sequence really does occur there -- and it
            // is cheap: the worst case measured on HG002 (a 12 bp probe, 2,029
            // loci, 22,371 reads) cost 3.9 MB of output and 15.7 s, with RAM
            // flat because it is dominated by loading the pseudogenome. So this
            // is NOT capped: silently truncating would turn a correct answer
            // into an arbitrary one. It is flagged instead, because a caller
            // could otherwise read a repeat hit as a specific locus.
            if(ranges.size() > 32)
                fprintf(stderr,"[query] NOTE: %zu loci -- this sequence is repetitive, so the "
                               "reads returned span many places. A longer probe is more "
                               "specific (measured on HG002: 12bp->2029 loci, 40bp->76, "
                               "100bp->2).\n", ranges.size());
            if(ranges.empty()){
                fprintf(stderr,"[query] not found -- no reads emitted. A long query may "
                               "carry a sequencing error; try a 30-60 bp sub-sequence.\n");
                FILE* fe=fopen(outdir.c_str(),"wb"); if(fe) fclose(fe);
                return 0;
            }
        } else {
            uint64_t qa=0,qb=0;
            const char* d=strchr(modearg.c_str(),'-');
            if(!d){ fprintf(stderr,"query needs START-END or a DNA sequence\n"); return 2; }
            qa=strtoull(modearg.c_str(),nullptr,10); qb=strtoull(d+1,nullptr,10);
            ranges.push_back({qa,qb});
        }
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
            bool hit=false;
            for(const auto& r : ranges) if(b>r.first && a<r.second){ hit=true; break; }
            if(!hit) continue;
            uint64_t e=std::min<uint64_t>(pg.size(),b);
            fprintf(f,">r%zu pos=%llu len=%llu\n",u,(unsigned long long)a,(unsigned long long)(e-a));
            fwrite(pg.data()+a,1,e-a,f); fputc('\n',f); ++n;
        }
        fclose(f);
        QLAP("scan all reads + emit answer");
        fprintf(stderr,"[query] %zu reads over %zu range(s) -> %s\n",n,ranges.size(),outdir.c_str());
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
    // ── VARIANT-SITE TABLE — the projection compression discards ────────────
    // The compressor must compute, for every read, how it deviates from the
    // consensus; that IS how it compresses. Doing so computes the complete
    // pileup. It then stores only the READ-KEYED projection, because that is
    // the only one decompression walks, and codes it with an adaptive model in
    // read order -- which makes read #k's deviations unreachable without
    // decoding all k-1 before it. The POSITION-KEYED projection, the pileup,
    // is therefore not merely absent but structurally unrecoverable, even
    // though the encoder held it.
    //
    // Measured on HG002 chr20: 7,466,871 deviations over 6,430,543 positions,
    // of which only 15,515 (position,base) pairs recur >=5x. 98.6% of the
    // stream is sequencing noise; 1.4% is variant signal. Retaining just the
    // recurrent part, delta-coded, costs 0.0606% of the archive at >=3x.
    if(WANT_PILEUP){
        // Append each read's deviations to the sidecar written above.
        // ── PER-READ DEVIATIONS ────────────────────────────────────────────
        // Without these, `query` emits the CONSENSUS at each read's position,
        // so every returned read agrees with every other -- measured, 121
        // overlapping pairs with zero mismatching bases. That is not a pileup,
        // it is a consensus repeated N times, and it is why variants encoded as
        // per-read mismatches (rather than as separate haplotype contigs) were
        // invisible to a locus query.
        //
        // They cannot be reached from the archive on demand: mm_sym is coded
        // with an ADAPTIVE model in READ order, so read k's deviations require
        // decoding all k-1 before it. That is a property of the coding, not an
        // oversight -- decompression only ever walks reads in order. The
        // sidecar breaks the ordering dependency by decoding once.
        //
        // Stored CSR-style and ONLY for reads that have deviations, so the cost
        // scales with error count and not with read count.
        { FILE* f=fopen(outdir.c_str(),"ab");
          if(!f){ fprintf(stderr,"cannot append to %s\\n",outdir.c_str()); return 1; }
          std::vector<uint32_t> ids; std::vector<uint32_t> offs; std::vector<uint8_t> syms;
          size_t run=0;
          for(size_t u=0; u<NU; ++u){
              const uint16_t c=(u<mmcount.size())?mmcount[u]:0;
              if(!c){ continue; }
              uint32_t prevj=0; size_t wrote=0;
              for(uint16_t m=0;m<c;++m){
                  if(run+m>=mmpos32.size()||run+m>=obs.size()) break;
                  uint32_t j=prevj+mmpos32[run+m]; prevj=j;
                  offs.push_back(j); syms.push_back(obs[run+m]); ++wrote;
              }
              if(wrote){ ids.push_back((uint32_t)u); ids.push_back((uint32_t)wrote); }
              run += c;
          }
          const uint64_t nid=ids.size()/2, nmm=offs.size();
          fwrite(&nid,8,1,f); fwrite(&nmm,8,1,f);
          fwrite(ids.data(),1,ids.size()*4,f);
          fwrite(offs.data(),1,offs.size()*4,f);
          fwrite(syms.data(),1,syms.size(),f);
          fprintf(stderr,"[index] + %llu reads carrying %llu deviations (%zu B)\n",
                  (unsigned long long)nid,(unsigned long long)nmm,
                  (size_t)(16+ids.size()*4+offs.size()*4+syms.size())); fclose(f); }
        std::vector<size_t> vmoff(NU+1,0);
        for(size_t u=0;u<NU;++u) vmoff[u+1]=vmoff[u]+(u<mmcount.size()?mmcount[u]:0);
        std::vector<uint32_t> vpos_(positions.size());
        memcpy(vpos_.data(),posb.data(),std::min(posb.size(),vpos_.size()*4));
        std::map<std::pair<uint32_t,uint8_t>,uint32_t> tally;
        for(size_t u=0;u<NU && u<vpos_.size();++u){
            const uint16_t c=(u<mmcount.size())?mmcount[u]:0;
            if(!c) continue;
            const uint32_t a0=vpos_[u];
            if(a0==UINT32_MAX) continue;
            size_t off=vmoff[u]; uint32_t prevj=0;
            for(uint16_t m=0;m<c;++m){
                if(off+m>=mmpos32.size()||off+m>=obs.size()) break;
                uint32_t j=prevj+mmpos32[off+m]; prevj=j;
                tally[{a0+j,obs[off+m]}]++;
            }
        }
        const uint32_t KMIN = getenv("CAPS_PILEUP_MIN")?(uint32_t)atoi(getenv("CAPS_PILEUP_MIN")):3;
        FILE* vf=fopen((outdir+".sites").c_str(),"wb");
        if(vf){
            uint64_t magic=0x5345544953504756ULL, nsite=0;
            for(auto& kv:tally) if(kv.second>=KMIN) ++nsite;
            fwrite(&magic,8,1,vf); fwrite(&nsite,8,1,vf);
            uint32_t prev=0;
            for(auto& kv:tally){
                if(kv.second<KMIN) continue;
                uint32_t d=kv.first.first-prev; prev=kv.first.first;
                while(d>=128){ uint8_t b=(uint8_t)(0x80|(d&0x7F)); fwrite(&b,1,1,vf); d>>=7; }
                { uint8_t b=(uint8_t)d; fwrite(&b,1,1,vf); }
                uint8_t base=kv.first.second;
                uint8_t cnt=(uint8_t)std::min<uint32_t>(255,kv.second);
                fwrite(&base,1,1,vf); fwrite(&cnt,1,1,vf);
            }
            long sz=ftell(vf); fclose(vf);
            fprintf(stderr,"[index] variant sites >=%ux: %llu -> %ld B (%.4f%% of archive)\n",
                    KMIN,(unsigned long long)nsite,sz,100.0*sz/(double)PGLEN);
        }
        return 0;
    }
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
    if(qual_async){
        qthread.join();
        fprintf(stderr,"  quality -> %s, overlapped with the sequence decode: %llu\n",
                out_qtext ? "text" : "bitmaps", (unsigned long long)qw_async);
        fprintf(stderr,"  [dec-timing] quality decode     %7.2fs (overlapped)\n", qsecs_async);
    } else if(has("qual_body") && out_qtext && !getenv("CAPS_SKIP_QUAL")){
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
    // Set BEFORE step 1 now, because step 1 is what hands the contigs back and
    // it needs to know which form this mode wants.
    if (WANT_INDELS) setenv("CAPSULE_EXPORT_CONTIGS", "1", 1);
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
                          &contigs_mem) != 0){
        fprintf(stderr,"[call] decode failed\n"); return 1; }

    _lap("1 decode reads+qual");
    // In indel mode the caller needs the INDIVIDUAL assembled contigs, not the
    // two concatenated pseudogenome records -- build_substrate collapses and
    // re-places reads per contig, so 2 giant records is a different operation.
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
        // READ THE ARCHIVE ONCE, NOT FIVE TIMES.
        //
        // arc_stream() calls read_capsule(), which opens and parses the WHOLE
        // container, and then returns a single stream. Five calls meant five
        // full passes over a 592 MB archive -- measured as 4.07 s in the
        // "3b parse contigs" lap, a lap that otherwise does nothing but swap a
        // vector now that the contigs arrive in memory.
        //
        // Same streams, same decode, same order; the file is just read once.
        // pos_sec/pos_region are REQUIRED here: pos_abs alone carries only the
        // MAIN-region positions since the region split, so reading it without
        // joining silently drops every second-region read's placement. This was
        // measured, not assumed -- the end-to-end benchmark reported
        // 371,009/460,501 placements, and the missing 89,492 is exactly the
        // second-region read count.
        std::vector<uint8_t> pb, psec, preg, sbv, spb, ofl, ovl;
        uint64_t _pl = 0, _me = 0;          // pg length / main-region end, needed by the join below
        {
            uint32_t _mm = 0;
            std::vector<Stream> _ss;
            if (read_capsule(in.c_str(), _pl, _me, _mm, _ss)) {
                for (auto& st : _ss) {
                    if      (st.name == "pos_abs")        pb  = capsule_decode_stream(st.coded, 1);
                    else if (st.name == "pos_sec")        psec = capsule_decode_stream(st.coded, 1);
                    else if (st.name == "pos_region")     preg = capsule_decode_stream(st.coded, 1);
                    else if (st.name == "pos_strand")     sbv = capsule_decode_stream(st.coded, 1);
                    else if (st.name == "contig_spans")   spb = capsule_decode_stream(st.coded, 1);
                    else if (st.name == "orig2uid_flags") ofl = capsule_decode_stream(st.coded, 1);
                    else if (st.name == "orig2uid_vals")  ovl = capsule_decode_stream(st.coded, 1);
                }
            }
        }
        if (pb.empty() || spb.empty()) {
            fprintf(stderr, "[call] ARCHIVE LACKS %s -- indels need it; "
                            "re-compress with CAPS_CALL=1 so contig_spans is written\n",
                    spb.empty() ? "contig_spans" : "pos_abs");
            return 1;
        }
        // Restore the full position array: pos_abs holds only main-region
        // positions since the split, and reading it alone silently drops every
        // second-region read.
        pb = caps_join_positions(pb, psec, preg, _me);
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
#ifndef CAPS_VERSION
#define CAPS_VERSION "1.0.0"
#endif

int main(int argc,char** argv){
    // --version must be answered BEFORE any argument parsing, or it is
    // treated as an input filename.
    if(argc>=2 && (!strcmp(argv[1],"--version")||!strcmp(argv[1],"-V"))){
        printf("g_capsul decoder %s (archive format v2)\n", CAPS_VERSION);
        return 0; }
    // Claim 2 from a STORED archive: capsule_decode call <in.capsule> <out.vcf> [workdir]
    if(argc>=4 && !strcmp(argv[1],"call"))
        return capsule_call_from_archive(argv[2], argv[3], argc>4?argv[4]:std::string());
    // Claim 3 modes:  capsule_decode export|coverage|query <in.capsule> <out> [range]
    if(argc>=4 && (!strcmp(argv[1],"export")||!strcmp(argv[1],"coverage")||!strcmp(argv[1],"query")||!strcmp(argv[1],"index")))
        return capsule_decode_all(argv[2], argv[3], std::string(), argv[1],
                                  argc>4?argv[4]:std::string());
    if(argc<3){ fprintf(stderr,"usage: capsule_decode <in.capsule> <outdir> [reads.out]\n"
                               "       capsule_decode export   <in.capsule> <out.fa>\n"
                               "       capsule_decode index    <in.capsule> <out.qidx>\n"
                               "       capsule_decode coverage <in.capsule> <out.tsv>\n"
                               "       capsule_decode query    <in.capsule> <out.fa> <START-END>\n"
                               "       capsule_decode call     <in.capsule> <out.vcf> [workdir]\n"); return 2; }
    return capsule_decode_all(argv[1], argv[2], argc>3?argv[3]:std::string());
}
#endif  // CAPSULE_NO_MAIN
