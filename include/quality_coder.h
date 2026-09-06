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
#include <atomic>
#include <algorithm>

extern "C" {
#include <omp.h>
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
// ── ONE POOL FOR (BLOCK x TRIAL), NOT A POOL PER BLOCK ──────────────────────
// encode_block parallelises its own trials, so with 8 trials on a 12-core box
// four cores sat idle for the whole quality stage and blocks were still
// processed one at a time. Splitting the block into prep / trial / select lets
// the caller schedule trials from SEVERAL blocks in one flat parallel loop,
// which fills the machine.
//
// Output is unchanged: the same blocks, the same trial set per block, the same
// fixed-index scan for the smallest result, and results appended in block
// order -- never completion order.
struct BlockJob {
    std::vector<std::string>          shifted;   // one per offset candidate
    std::vector<uint32_t>             lens;
    uint8_t                           cands[2] = {0,0};
    int                               nstrat = 4, nc = 1, ntrial = 4;
    std::vector<std::vector<uint8_t>> res;
    std::vector<char>                 ok;
    size_t                            qbytes = 0;
};

static void prep_block(const std::string& qbuf, const std::vector<uint32_t>& lens,
                       BlockJob& J){
    J.lens = lens; J.qbytes = qbuf.size();
    uint8_t qmin = 255;
    for(size_t i=0;i<qbuf.size();++i){ const uint8_t c=(uint8_t)qbuf[i]; if(c<qmin) qmin=c; }
    if(qbuf.empty()) qmin = 33;
    J.cands[0]=qmin; J.cands[1]=33;
    J.nc = (qmin!=33 && qmin>=33) ? 2 : 1;
    J.nstrat = 4;
    if(const char* e=getenv("CAPS_QSTRAT")){ int v=atoi(e); if(v>=1&&v<=4) J.nstrat=v; }
    if(const char* e=getenv("CAPS_QOFF")){ int v=atoi(e); if(v>=1&&v<=J.nc) J.nc=v; }
    J.ntrial = J.nc * J.nstrat;
    J.shifted.assign(J.nc, std::string());
    for(int ci=0; ci<J.nc; ++ci){
        J.shifted[ci].resize(qbuf.size());
        const uint8_t off=J.cands[ci];
        for(size_t i=0;i<qbuf.size();++i)
            J.shifted[ci][i]=(char)((unsigned char)qbuf[i]-off);
    }
    J.res.assign(J.ntrial, std::vector<uint8_t>());
    J.ok.assign(J.ntrial, 0);
}

static void run_trial(BlockJob& J, int t){
    const int ci = t / J.nstrat, strat = t % J.nstrat;
    // fqz_compress MUTATES the slice (flags at fqzcomp_qual.c:655, len at :790),
    // and copying the struct still shares those arrays -- eight threads
    // scribbling on each other's parameters grew the archive 1.55% once. Each
    // trial therefore owns its len/flags storage.
    std::vector<uint32_t> mylen(J.lens.begin(), J.lens.end());
    std::vector<uint32_t> myflags(J.lens.size(), 0u);
    fqz_slice st;
    st.num_records = (int)J.lens.size();
    st.len   = mylen.data();
    st.flags = myflags.data();
    size_t osz=0;
    char* c = fqz_compress(4, &st, const_cast<char*>(J.shifted[ci].data()),
                           J.shifted[ci].size(), &osz, strat, nullptr);
    if(c){ J.res[t].assign((uint8_t*)c,(uint8_t*)c+osz); J.ok[t]=1; free(c); }
}

static bool select_block(BlockJob& J, std::vector<uint8_t>& out, uint8_t& qmin_out){
    bool got=false; size_t best=0;
    for(int t=0; t<J.ntrial; ++t){
        if(!J.ok[t]) continue;
        if(!got || J.res[t].size() < best){
            out = J.res[t]; best = J.res[t].size(); got = true;
            qmin_out = J.cands[t/J.nstrat];
        }
    }
    return got;
}

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

    // ── THE 8 TRIALS RUN CONCURRENTLY ────────────────────────────────────
    // This tries every (offset, strategy) pair and keeps the smallest, which
    // is 2 x 4 = 8 FULL compressions of the same block. That is why this coder
    // measured 16.6 MB/s where fqzcomp alone runs an order of magnitude faster
    // -- the throughput was never the coder, it was doing the coder's work
    // eight times.
    //
    // The eight trials are completely independent: each reads a shifted copy
    // of the block and writes its own buffer. Running them concurrently is
    // therefore FREE SPEED WITH NO SIZE COST -- the same eight candidates are
    // evaluated and the same one wins.
    //
    // DETERMINISM. The serial loop kept the first STRICTLY smaller result, so
    // ties resolved to the lowest (ci, strat). Results are collected into a
    // fixed-index array and scanned in that same order, so the selection is
    // identical regardless of completion order -- verified byte-identical.
    // CAPS_QSTRAT / CAPS_QOFF limit the trial grid, for measuring what the 8
    // trials actually buy. Unset, both are the full grid and the output is
    // bit-identical to before. Strategy order matters: the scan below keeps the
    // first strictly-smaller result, so restricting to the first k strategies
    // evaluates a PREFIX of the same ordered candidate list.
    // ── ONE STRATEGY, NOT FOUR. MEASURED 2026-09-05 ─────────────────────
    // The four fqzcomp strategies were being tried in full on every block and
    // buy NOTHING. Measured by varying the trial grid and comparing archives:
    //
    //   E. coli   2 offsets x 4 strategies  68,677,977   2 x 2  68,677,977  (same)
    //             1 offset  x 4 strategies  68,771,735   1 x 2  68,771,735  (same)
    //             1 offset  x 1 strategy    68,771,735               (same)
    //   HG002     2 x 4  574,014,786        2 x 2  574,014,786       (same)
    //
    // So strategy contributes 0.000% on both, while the OFFSET dimension is
    // worth 0.137% and is kept. That turns 8 full compressions of every 256 MB
    // block into 2 -- this coder's own comment already identified the 8x as why
    // it ran at 16.6 MB/s when fqzcomp alone is an order of magnitude faster.
    //
    // CAPS_QSTRAT=4 restores the full strategy search.
    // REVERTED TO 4 AFTER MEASUREMENT. Dropping to one strategy looked free on
    // E. coli and HG002 (both 0.000%) and is NOT: on S. acidocaldarius it costs
    // +476 B (+0.003%), isolated exactly --
    //     reference            15,159,145
    //     new code, QSTRAT=4   15,159,145  IDENTICAL  (mext + pipelining clean)
    //     new code, QSTRAT=1   15,159,621  DIFFERS
    // Two datasets are not evidence for a default that touches 82% of the
    // archive. CAPS_QSTRAT=1 remains available for anyone accepting that cost.
    int nstrat = 4;
    if(const char* e=getenv("CAPS_QSTRAT")){ int v=atoi(e); if(v>=1&&v<=4) nstrat=v; }
    int nc = ncand;
    if(const char* e=getenv("CAPS_QOFF")){ int v=atoi(e); if(v>=1&&v<=ncand) nc=v; }
    const int ntrial = nc * nstrat;
    std::vector<std::string> shifted(ncand);
    for(int ci=0; ci<ncand; ++ci){
        shifted[ci].resize(qbuf.size());
        const uint8_t off=cands[ci];
        for(size_t i=0;i<qbuf.size();++i)
            shifted[ci][i]=(char)((unsigned char)qbuf[i]-off);
    }
    std::vector<std::vector<uint8_t>> res(ntrial);
    std::vector<char> ok(ntrial, 0);
    #pragma omp parallel for schedule(dynamic,1)
    for(int t=0; t<ntrial; ++t){
        const int ci = t / nstrat, strat = t % nstrat;
        // fqz_compress MUTATES the slice -- it writes s->flags[rec]
        // (fqzcomp_qual.c:655) and s->len[i] (:790). Copying the struct alone
        // is not enough because the copy still points at the SAME len/flags
        // arrays, so eight threads were scribbling over each other's
        // parameters. Measured cost of that bug: the archive grew 1.55% and
        // stopped matching the serial output. Each trial therefore gets its
        // own len and flags storage.
        std::vector<uint32_t> mylen(lens.begin(), lens.end());
        std::vector<uint32_t> myflags(lens.size(), 0u);
        fqz_slice st;
        st.num_records = (int)lens.size();
        st.len   = mylen.data();
        st.flags = myflags.data();
        size_t osz=0;
        char* c = fqz_compress(4 /*CRAM 4.0 vers*/, &st,
                               const_cast<char*>(shifted[ci].data()),
                               shifted[ci].size(), &osz, strat, nullptr);
        if(c){ res[t].assign((uint8_t*)c,(uint8_t*)c+osz); ok[t]=1; free(c); }
    }
    bool got=false; size_t best=0;
    for(int t=0; t<ntrial; ++t){
        if(!ok[t]) continue;
        if(!got || res[t].size() < best){
            out = res[t]; best = res[t].size(); got = true;
            qmin_out = cands[t/nstrat];
        }
    }
    return got;
}

// Streams the quality column straight out of the FASTQ, one block at a time --
// never materializing the whole column, matching the bounded-memory property
// the names coder has.
static Encoded encode_from_fastq(const char* fq_path, size_t BLOCK_BYTES=QBLOCK_BYTES){
    // Block size is overridable for measurement: it trades compression (each
    // block resets the fqzcomp context) against parallelism (blocks are
    // independent). Default is unchanged.
    if(const char* e = getenv("CAPS_QBLOCK_MB")){ long v=atol(e); if(v>0) BLOCK_BYTES=(size_t)v*1024u*1024u; }
    Encoded E;
    FILE* f=fopen(fq_path,"r");
    if(!f) return E;
    std::vector<char> buf(1u<<16);
    std::string qbuf; qbuf.reserve(BLOCK_BYTES+65536);
    std::vector<uint32_t> lens;
    uint64_t lineno=0;

    // ── PIPELINE BLOCKS, NOT JUST TRIALS ────────────────────────────────────
    // Blocks were encoded strictly one at a time while the trials inside a
    // block ran in parallel. With the strategy search removed there are only 2
    // trials per block, so 2 of 12 cores would work and the wall time would not
    // move at all -- dropping the trials FREES cores, it does not by itself
    // save time. Batching blocks is what converts that into speed.
    //
    // Output is unchanged: the same blocks are cut at the same boundaries, each
    // is encoded by the same encode_block, and results are appended in BLOCK
    // ORDER after the batch completes -- not in completion order.
    //
    // Batch size is bounded by memory, not by core count: a block in flight
    // costs its own bytes plus the two shifted copies encode_block makes, so
    // ~3x BLOCK_BYTES each. Four 256 MB blocks is ~3 GB, which sits under the
    // run's existing peak (measured 4.2 GB, reached later during mapping), so
    // this buys parallelism without moving peak RSS.
    // NB MUST TRACK THE TRIAL COUNT, or this makes things WORSE. encode_block
    // parallelises its trials; batching blocks puts that inside another
    // parallel region, and with nested parallelism off (the default) the inner
    // loop collapses to serial. At 8 trials and NB=4 that turns 17 blocks x 1
    // parallel trial-round into 5 batches x 8 SERIAL trials -- a 2.4x
    // regression. So the batch is sized to keep total parallelism at ~one team:
    // NB = threads / trials-per-block, which is 1 at the default 8 trials
    // (exactly the previous behaviour) and 6 when CAPS_QSTRAT=1 makes it 2.
    size_t NB;
    {
        unsigned thr = (unsigned)omp_get_max_threads(); if(!thr) thr=1;
        int nstrat_est = 4;
        if(const char* e=getenv("CAPS_QSTRAT")){ int v=atoi(e); if(v>=1&&v<=4) nstrat_est=v; }
        const unsigned trials = (unsigned)(2*nstrat_est);        // 2 offsets, worst case
        // ONE BLOCK. MEASURED AND REVERTED 2026-09-06. Filling the idle cores
        // with a second block looked obviously right and LOSES once quality is
        // overlapped with the assembly (which it now is): the worker already
        // shares the machine with round 1, so a second block only doubles its
        // footprint and oversubscribes. HG002, same archive every time:
        //     overlap, 1 block   279.37 s   peak 4870 MB
        //     overlap, 2 blocks  282.73 s   peak 6190 MB
        //                        280.54 s   peak 6271 MB
        // 1-3 s slower for +1.3 GB. CAPS_QBATCH>1 re-enables it for measurement.
        (void)thr; (void)trials;
        NB = 1;
    }
    if(const char* e=getenv("CAPS_QBATCH")){ long v=atol(e); if(v>0) NB=(size_t)v; }
    {
        size_t avail_mb = 0;
        if(FILE* mi = fopen("/proc/meminfo","r")){
            char k[64]; unsigned long v; char u[16];
            while(fscanf(mi,"%63s %lu %15s",k,&v,u)>=2)
                if(!strcmp(k,"MemAvailable:")){ avail_mb=v/1024; break; }
            fclose(mi);
        }
        const size_t per_mb = (BLOCK_BYTES*3)/(1024*1024) + 1;
        if(avail_mb){ size_t fit = (avail_mb/4)/per_mb; if(fit<1) fit=1; if(NB>fit) NB=fit; }
        if(NB<1) NB=1;
    }
    std::vector<std::string>            pend_q;
    std::vector<std::vector<uint32_t>>  pend_l;

    auto drain=[&](){
        const size_t nb = pend_q.size();
        if(!nb) return;
        std::vector<std::vector<uint8_t>> blk(nb);
        std::vector<uint8_t> qmn(nb, 0);
        std::vector<char>    okb(nb, 0);
        std::vector<BlockJob> jobs(nb);
        for(size_t b=0;b<nb;++b) prep_block(pend_q[b], pend_l[b], jobs[b]);
        std::vector<std::pair<size_t,int>> tasks;          // (block, trial)
        for(size_t b=0;b<nb;++b)
            for(int t=0;t<jobs[b].ntrial;++t) tasks.push_back({b,t});
        #pragma omp parallel for schedule(dynamic,1)
        for(long long ti=0; ti<(long long)tasks.size(); ++ti)
            run_trial(jobs[tasks[(size_t)ti].first], tasks[(size_t)ti].second);
        for(size_t b=0;b<nb;++b){
            uint8_t qm=0;
            if(select_block(jobs[b], blk[b], qm)){ okb[b]=1; qmn[b]=qm; }
        }
        for(size_t b=0;b<nb;++b){                    // strict block order
            if(okb[b]){
                pv(E.index, blk[b].size());
                pv(E.index, pend_l[b].size());
                pv(E.index, qmn[b]);
                E.body.insert(E.body.end(), blk[b].begin(), blk[b].end());
                ++E.n_blocks;
            }
            E.n_reads  += pend_l[b].size();
            E.n_qbytes += pend_q[b].size();
        }
        pend_q.clear(); pend_l.clear();
    };

    auto flush=[&](){
        if(lens.empty()) return;
        pend_q.emplace_back(std::move(qbuf));
        pend_l.emplace_back(std::move(lens));
        qbuf.clear(); lens.clear();
        qbuf.reserve(BLOCK_BYTES+65536);
        if(pend_q.size() >= NB) drain();
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
    drain();
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

// ── PARALLEL DECODE STRAIGHT TO THE CALLER'S BITMAP ───────────────────────
// decode_to_file above is a SERIAL loop over blocks that are entirely
// independent: each carries its own length, count and qmin, is fqz_decompress-ed
// on its own, and nothing crosses a block boundary. It measured 35.25 s of a
// 147 s call-from-archive run -- 24% -- on one core of twelve.
//
// It is also decoding to the wrong representation for this consumer. The caller
// reduces every quality string to ONE BIT per base ("is this base >= QMIN"), so
// materialising 1.86 GB of text, writing it to disk and reading it back to
// extract those bits is pure overhead. This produces the packed bitmaps
// directly, in parallel, in the same layout capspack::pack_qual yields.
//
// Two passes: the index is walked once to record each block's byte offset and
// first read index (cheap -- it is a few thousand varints), then the blocks are
// decompressed concurrently. Every block writes a disjoint range of `out`, so
// no synchronisation is needed and the result is order-independent.
static uint64_t decode_to_bitmaps(const std::vector<uint8_t>& body,
                                  const std::vector<uint8_t>& indexRaw,
                                  const std::vector<uint32_t>& lengths,
                                  int qmin_want,
                                  std::vector<std::string>& out)
{
    struct Blk { size_t boff, blen, bcnt, first; uint8_t qmin; };
    std::vector<Blk> blks;
    {
        const uint8_t* ip=indexRaw.data(); const uint8_t* iend=ip+indexRaw.size();
        size_t boff=0; uint64_t li=0;
        while(ip<iend){
            const uint64_t blen=gv(ip,iend);
            const uint64_t bcnt=gv(ip,iend);
            const uint8_t  qmin=(uint8_t)gv(ip,iend);
            if(boff+blen>body.size()) break;
            blks.push_back({boff,(size_t)blen,(size_t)bcnt,(size_t)li,qmin});
            boff+=blen; li+=bcnt;
            if(li>=lengths.size()+1) break;
        }
    }
    if(out.size()<lengths.size()) out.resize(lengths.size());
    std::atomic<uint64_t> written{0};
    #pragma omp parallel for schedule(dynamic,1)
    for(long long b=0;b<(long long)blks.size();++b){
        const Blk& B=blks[(size_t)b];
        size_t osz=0;
        char* d=fqz_decompress((char*)body.data()+B.boff,B.blen,&osz,nullptr,0);
        if(!d) continue;
        size_t off=0; uint64_t w=0;
        for(size_t r=0;r<B.bcnt && (B.first+r)<lengths.size();++r){
            const uint32_t L=lengths[B.first+r];
            if(off+L>osz) break;
            std::string packed((L+7)/8,'\0');
            for(uint32_t i=0;i<L;++i){
                const int q=(int)((unsigned char)d[off+i]+B.qmin)-33;
                if(q>=qmin_want) packed[i>>3]|=(char)(1u<<(i&7));
            }
            out[B.first+r].swap(packed);
            off+=L; ++w;
        }
        free(d);
        written+=w;
    }
    return written.load();
}

} // namespace qlc
