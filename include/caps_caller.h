#pragma once
// ── CAPSULE Claim-2 caller ──────────────────────────────────────────────────
// A faithful port of the outer ARCS project's reference-free bubble caller
// (/root/arcs-clean/src/caller.cpp), retargeted at CAPSULE's own assembly
// output instead of ARCS's build_vodbg_pg / CallData.
//
// Frozen parameters (do not re-tune — these are ARCS's own validated,
// held-out-tested values, ported verbatim; see
// /root/arcs-clean/docs/VARIANT_ANALYSIS_MASTER.md):
//   HDMAX=2, MAF=0.20, DHI=2.5, KHI=1.1, MC=3, TRI=0.12, HALF=15.
//
// Structural note: ARCS's Method B assembler (vodbg_pg) keeps reads in
// several DISCRETE contigs by construction. CAPSULE's own assembler stitches
// everything into one (or two) long pg buffers via MEM references, which
// would erase the haplotype separation the indel-bubble mechanism needs. The
// fix (see stages/106_inprocess.cpp, CAPS_CALL gate): every chain emitted in
// round-1/round-2 (main pg) and in the second-region sweep is its own
// contiguous span of `pg`, captured as a contig BEFORE the MEM self-match
// stage runs (MEM never mutates `pg` itself, only which of its bytes get
// literal-encoded — so contig boundaries captured pre-MEM are still valid
// after MEM completes). Pigeonhole-mapped reads land inside an
// already-emitted span by construction (they only match against already-
// appended pg content), so they inherit that span's contig id for free.
//
// This header owns ONLY the caller logic; contig-span capture and CallData
// population live in 106_inprocess.cpp, gated on CAPS_CALL=1, zero cost when
// unset (matches the CAPS_NAMES / CAPS_QUAL convention).

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <string>
#include <array>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <climits>

namespace capscall {

struct CallData {
    std::vector<std::string> contigs;     // one entry per assembled contig (pre-MEM span)
    std::vector<uint32_t>    read_cid;    // [original read idx] -> contig id
    std::vector<uint32_t>    read_pos;    // [original read idx] -> contig-local start
    std::vector<uint8_t>     read_rc;     // [original read idx] -> reverse-complement flag
    // Bases to SKIP at the start of the (already strand-oriented) read. Non-zero
    // when the read hangs off the LEFT end of its contig: the read's base
    // read_clip maps to contig position read_pos. Without this the read would
    // have to be dropped, and reads near contig starts are a large share of a
    // fragmented substrate.
    std::vector<uint16_t>    read_clip;
    bool valid = false;
};

namespace detail {

constexpr int HDMAX = 2, MC = 3, HALF = 15;
constexpr double MAF = 0.20, DHI = 2.5, KHI = 1.1, TRI = 0.12;
constexpr int MAXINDEL = 30;

inline int b2i(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; }
    return -1;
}
inline char compl_base(char c) {
    switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; }
    return 'N';
}
inline std::string rc_str(const std::string& s) {
    std::string r(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) r[s.size() - 1 - i] = compl_base(s[i]);
    return r;
}

inline bool pack31(const char* s, uint64_t& out) {
    uint64_t v = 0;
    for (int i = 0; i < 31; ++i) { int b = b2i(s[i]); if (b < 0) return false; v = (v << 2) | (uint64_t)b; }
    out = v; return true;
}
inline uint64_t rc31(uint64_t v) {
    uint64_t r = 0;
    for (int i = 0; i < 31; ++i) { r = (r << 2) | (3u - (v & 3u)); v >>= 2; }
    return r;
}
inline uint64_t canon31(uint64_t v) { uint64_t r = rc31(v); return v < r ? v : r; }
inline uint64_t colkey(uint32_t cid, uint32_t pos) { return ((uint64_t)cid << 32) | pos; }

inline bool is_str_event(const std::string& seq, const std::string& flank) {
    // SINGLE-BASE indels were exempt from this filter entirely (the old guard
    // was `seq.size() < 2 -> false`), yet they are the dominant false-positive
    // class: on HG002 r4, 16 of 24 indel FPs were 1 bp deletions against only
    // 6 TPs of that class. A 1 bp indel adjacent to a run of the same base is
    // a homopolymer-length artifact -- the known hard case for every
    // reference-free caller. Detect it by the actual run length in the left
    // flank (a measured property of the sequence, not a per-dataset rule)
    // rather than by mere presence of the base, which would fire on almost
    // everything.
    if (seq.size() == 1) {
        int run = 0;
        for (int i = (int)flank.size() - 1; i >= 0 && flank[(size_t)i] == seq[0]; --i) ++run;
        return run >= 2;              // deleted/inserted base extends a >=2 run
    }
    if (seq.size() < 2) return false;
    for (int u = 1; u <= std::min((int)seq.size(), 4); ++u) {
        bool tandem = true;
        for (size_t i = 0; i < seq.size(); i += (size_t)u) {
            size_t rem = std::min((size_t)u, seq.size() - i);
            if (seq.substr(i, rem) != seq.substr(0, rem)) { tandem = false; break; }
        }
        if (!tandem) continue;
        std::string unit = seq.substr(0, (size_t)u);
        if (flank.find(unit) != std::string::npos) return true;
    }
    return false;
}

struct Bubble { int type = 0; uint32_t apos = 0; int len = 0; std::string ins; bool ok = false; };

inline bool flank_match(const std::string& A, size_t ai, const std::string& B, size_t bi, int F) {
    if (ai + (size_t)F > A.size() || bi + (size_t)F > B.size()) return false;
    for (int k = 0; k < F; ++k) if (A[ai + (size_t)k] != B[bi + (size_t)k]) return false;
    return true;
}

// REFUTED, kept for the record (not wired). Tolerant re-convergence test.
// Tried because:  on HG002 r2, 29 of the 32
// MISSED indels already have a contig PAIR covering them -- the substrate is
// fine, the extractor is what fails. It demanded a byte-identical 15 bp flank,
// so a single heterozygous SNV anywhere in that flank destroyed the bubble.
// DiscoSNP++ explicitly models this case (its own header calls a bubble
// "an indel+n SNPs"), which is part of why its indel recall beats ours.
// Tolerance reuses the project's already-frozen HDMAX=2 rather than
// introducing a new constant.
// RESULT: net NEGATIVE -- r2 indel F1 0.602->0.595, r5 0.600->0.582,
// HG005 0.542->0.508. Recall rose slightly but precision fell more: the
// bubble GEOMETRY is not what limits us, and loosening it only admits FPs.
inline bool flank_match_tol(const std::string& A, size_t ai, const std::string& B, size_t bi,
                            int F, int tol) {
    if (ai + (size_t)F > A.size() || bi + (size_t)F > B.size()) return false;
    int mm = 0;
    for (int k = 0; k < F; ++k)
        if (A[ai + (size_t)k] != B[bi + (size_t)k] && ++mm > tol) return false;
    return true;
}

inline Bubble extract_bubble(const std::string& A, uint32_t pA, const std::string& B, uint32_t qB,
                              int maxindel, int FLANK) {
    Bubble r;
    const size_t la = A.size(), lb = B.size();
    uint32_t d = 0;
    while (pA + d < la && qB + d < lb && A[pA + d] == B[qB + d]) ++d;
    if (d == 0) return r;
    if (pA + d >= la || qB + d >= lb) return r;
    for (int g = 1; g <= maxindel; ++g) {
        if (qB + d + (uint32_t)g + (uint32_t)FLANK > lb) break;
        if (flank_match(A, pA + d, B, qB + d + (uint32_t)g, FLANK)) {
            r.type = 1; r.apos = pA + d; r.len = g;
            r.ins = B.substr(qB + d, (size_t)g); r.ok = true; return r;
        }
    }
    for (int g = 1; g <= maxindel; ++g) {
        if (pA + d + (uint32_t)g + (uint32_t)FLANK > la) break;
        if (flank_match(A, pA + d + (uint32_t)g, B, qB + d, FLANK)) {
            r.type = 0; r.apos = pA + d; r.len = g; r.ok = true; return r;
        }
    }
    return r;
}

// ── Cross-contig SNV bubble ───────────────────────────────────────────────────
// A 1-substitution bubble between two haplotype contigs: they share flanking
// sequence but differ at exactly one base. CAPSULE's assembler (exact
// suffix-prefix overlap chaining, no mismatch tolerance in the chain step
// itself) fragments a chain at EVERY heterozygous site by construction, so
// this bubble pattern -- not a within-contig pileup mismatch -- is the
// PRIMARY het-SNV signal here (the reverse of ARCS's own vodbg_pg assembler,
// where this pass is experimental/opt-in because most het-SNVs already show
// up as pileup columns on one shared consensus contig). Default ON for
// CAPSULE; disable with CAPS_NO_XSNV=1 to reproduce a pileup-only caller.
struct SnvBubble { bool ok = false; uint32_t apos = 0; char ref_base = 0, alt_base = 0; };

inline SnvBubble extract_snv_bubble(const std::string& A, uint32_t pA,
                                     const std::string& B, uint32_t qB, int FLANK) {
    SnvBubble r;
    const size_t la = A.size(), lb = B.size();
    uint32_t d = 0;
    while (pA + d < la && qB + d < lb && A[pA + d] == B[qB + d]) ++d;
    if (pA + d >= la || qB + d >= lb) return r;
    char ra = A[pA + d], rb = B[qB + d];
    if (b2i(ra) < 0 || b2i(rb) < 0) return r;
    if (pA + d + 1u + (uint32_t)FLANK > (uint32_t)la) return r;
    if (qB + d + 1u + (uint32_t)FLANK > (uint32_t)lb) return r;
    if (!flank_match(A, pA + d + 1u, B, qB + d + 1u, FLANK)) return r;
    r.ok = true; r.apos = pA + d; r.ref_base = ra; r.alt_base = rb;
    return r;
}

inline bool pack25(const char* s, uint64_t& out) {
    uint64_t v = 0;
    for (int i = 0; i < 25; ++i) { int b = b2i(s[i]); if (b < 0) return false; v = (v << 2) | (uint64_t)b; }
    out = v; return true;
}
inline uint64_t rc25(uint64_t v) {
    uint64_t r = 0;
    for (int i = 0; i < 25; ++i) { r = (r << 2) | (3u - (v & 3u)); v >>= 2; }
    return r;
}

inline double loglik(const std::vector<std::pair<int,int>>& rl, const int* al, int m) {
    double ll = 0.0;
    for (auto& bq : rl) {
        double e = std::pow(10.0, -bq.second / 10.0);
        double p = 0.0;
        for (int a = 0; a < m; ++a) p += (bq.first == al[a]) ? (1.0 - e) : (e / 3.0);
        ll += std::log(std::max(p / m, 1e-300));
    }
    return ll;
}

} // namespace detail

// seqs/quals are ORIGINAL reads in original order (same indexing as cd.read_cid etc).
// quals may be empty strings per-read (treated as q=0, i.e. filtered out by the
// medq>=20 gate on the minor allele unless real quality is supplied) — pass real
// quality lines whenever available.
// ── Calling substrate rebuild ────────────────────────────────────────────────
// MEASURED MOTIVATION (2026-09-02, HG002 chr20 r2, 422 truth het-SNV sites):
//   * the reads contain the variants: 99.5% of truth sites carry BOTH alleles
//     with >=3 reads each at depth >=6 -- the ceiling is not the data;
//   * per-contig depth is NOT the blocker: 96.2% of sites already have >=6
//     reads stacked on one contig;
//   * the blocker is READ->CONTIG ASSIGNMENT. CAPSULE assigns each read to the
//     chain it EXACTLY overlaps, so ref-allele and alt-allele reads land on
//     different chains and never meet in a pileup. Only 124 of 75,115 reads
//     ever reached the mismatch-tolerant pigeonhole mapper -- chaining claimed
//     the rest. Re-assigning reads mismatch-tolerantly makes 70.1% of truth
//     sites immediately show a both-allele stack.
//   * reads covering one site are additionally spread over a MEDIAN of 5
//     contigs, so even a correct assignment splits the evidence -- hence the
//     collapse pass below, which keeps one representative contig per locus.
//
// Both passes are gated and default ON for calling only; they never touch the
// compression path. CAPS_NO_REMAP=1 / CAPS_NO_COLLAPSE=1 restore the previous
// behaviour exactly, so every measurement has a revert.
struct GapEvent;

struct Substrate {
    std::vector<std::string> contigs;
    std::vector<uint32_t> read_cid, read_pos;
    std::vector<uint8_t>  read_rc;
    std::vector<uint16_t> read_clip;
    std::vector<GapEvent> gaps;      // read-level indel evidence (second channel)
};

namespace detail {

// Greedy non-redundant contig set: longest first, drop a contig whose k-mers
// are already almost entirely claimed by a kept (longer) contig. This is the
// "shared sequence collapses to one place" property a de Bruijn graph gets for
// free, obtained without building a graph. Keyed on a measured fraction, not a
// per-dataset constant.
inline std::vector<uint32_t> collapse_contigs(const std::vector<std::string>& ctgs,
                                              double dup_frac, int K) {
    std::vector<uint32_t> order(ctgs.size());
    for (uint32_t i = 0; i < ctgs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](uint32_t a, uint32_t b){ return ctgs[a].size() > ctgs[b].size(); });
    std::unordered_set<uint64_t> claimed;
    claimed.reserve(1u << 20);
    std::vector<uint32_t> keep;
    std::vector<uint64_t> km;
    for (uint32_t ci : order) {
        const std::string& c = ctgs[ci];
        if ((int)c.size() < K) { keep.push_back(ci); continue; }
        // Query on a SAMPLED stride but claim EVERY k-mer: a duplicate contig
        // is almost always phase-shifted relative to its twin, so sampling both
        // sides at the same stride compares disjoint k-mer sets and the
        // duplicate is never detected (measured: only 7.6% of contigs collapsed
        // before this fix).
        km.clear();
        size_t hit = 0, tot = 0;
        for (size_t i = 0; i + (size_t)K <= c.size(); i += 5) {
            uint64_t v; if (!pack25(c.data() + i, v)) continue;
            uint64_t rcv = rc25(v), can = v < rcv ? v : rcv;
            ++tot;
            if (claimed.count(can)) ++hit;
        }
        if (tot > 0 && (double)hit / (double)tot >= dup_frac) continue;   // redundant
        keep.push_back(ci);
        for (size_t i = 0; i + (size_t)K <= c.size(); ++i) {
            uint64_t v; if (!pack25(c.data() + i, v)) continue;
            uint64_t rcv = rc25(v), can = v < rcv ? v : rcv;
            claimed.insert(can);
        }
    }
    std::sort(keep.begin(), keep.end());
    return keep;
}

} // namespace detail

// Rebuild the substrate: collapse redundant contigs, then place EVERY read on
// the surviving contigs allowing mismatches (both strands), keeping the best
// placement. Falls back to the encoder's own placement for reads that cannot
// be placed, so nothing is lost relative to the previous behaviour.
// ── Read-level gapped indel evidence ────────────────────────────────────────
// SECOND, INDEPENDENT indel channel. The bubble passes need two contigs that
// both exist and differ; measured on HG002 r2 that misses most events (32 FN
// against only 13 FP -- recall, not precision, is the indel limiter).
// A read that SPANS an indel carries the evidence directly: its seeds before
// the event and after it imply contig start positions that differ by exactly
// the indel length. No gapped aligner is needed -- the seed index already
// built for placement supplies both anchors, so this costs one extra pass.
// An event is emitted only when >= MC independent reads agree on the same
// (contig, boundary, signed length), which is the same read-support standard
// the junction test applies to bubbles.
struct GapEvent { uint32_t cid; uint32_t pos; int len; int support; std::string ins; };

inline std::vector<GapEvent> gapped_indel_scan(const std::vector<std::string>& seqs,
                                               const Substrate& S,
                                               const std::unordered_map<uint64_t,
                                                     std::vector<std::pair<uint32_t,uint32_t>>>& idx,
                                               int K) {
    using namespace detail;
    // key: (contig, boundary/8, signed length, inserted bases). Requiring reads
    // to agree on the INSERTED SEQUENCE as well as the length is a free
    // consistency check -- independent reads spanning the same real insertion
    // agree on its bases; coincidental seed splits do not.
    std::map<std::tuple<uint32_t,int32_t,int32_t,std::string>, int> votes;
    for (size_t o = 0; o < seqs.size(); ++o) {
        uint32_t cid = S.read_cid[o];
        if (cid == UINT32_MAX) continue;
        std::string r = S.read_rc[o] ? rc_str(seqs[o]) : seqs[o];
        const int rl = (int)r.size();
        if (rl < 2 * K + 4) continue;
        // implied contig start from each seed offset, restricted to this read's contig
        std::vector<std::pair<int,int64_t>> imp;    // (read offset, implied start)
        for (int off = 0; off + K <= rl; off += 4) {
            uint64_t v; if (!pack25(r.data() + off, v)) continue;
            uint64_t rcv = rc25(v), can = v < rcv ? v : rcv;
            auto it = idx.find(can);
            if (it == idx.end()) continue;
            for (auto& pr : it->second)
                if (pr.first == cid) imp.push_back({off, (int64_t)pr.second - off});
        }
        if (imp.size() < 4) continue;
        // Take the two DOMINANT implied starts rather than demanding that every
        // seed agree. A single spurious repeat hit inside the same contig is
        // enough to break a strict all-agree test -- which is why the first
        // version of this scan found exactly 1 event in a whole window.
        std::map<int64_t,int> cnt;
        for (auto& pr : imp) cnt[pr.second]++;
        if (cnt.size() < 2) continue;
        std::vector<std::pair<int,int64_t>> byc;
        for (auto& kv : cnt) byc.push_back({kv.second, kv.first});
        std::sort(byc.rbegin(), byc.rend());
        if (byc[0].first < 2 || byc[1].first < 2) continue;
        int64_t a = byc[0].second, b = byc[1].second;
        int64_t mean_a = 0, mean_b = 0; int na = 0, nb = 0;
        for (auto& pr : imp) {
            if (pr.second == a) { mean_a += pr.first; ++na; }
            else if (pr.second == b) { mean_b += pr.first; ++nb; }
        }
        if (!na || !nb) continue;
        int64_t s0 = (mean_a / na <= mean_b / nb) ? a : b;
        int64_t s1 = (s0 == a) ? b : a;
        int64_t g = s1 - s0;
        if (g == 0 || std::llabs(g) > MAXINDEL) continue;
        int last_s0 = -1, first_s1 = rl + 1;
        for (auto& pr : imp) {
            if (pr.second == s0) last_s0 = std::max(last_s0, pr.first);
            else if (pr.second == s1) first_s1 = std::min(first_s1, pr.first);
        }
        if (last_s0 < 0 || first_s1 > rl || last_s0 >= first_s1) continue;
        // Exact divergence point: walk the read against the contig from the
        // last agreeing seed until the first mismatch. Voting on a rounded
        // boundary would smear the position by up to the seed stride, and an
        // indel called even a few bases off does not match the truth set.
        const std::string& cseq = S.contigs[cid];
        int64_t bp = last_s0 + K;
        while (bp < rl && s0 + bp >= 0 && s0 + bp < (int64_t)cseq.size()
               && r[(size_t)bp] == cseq[(size_t)(s0 + bp)]) ++bp;
        int64_t boundary = s0 + bp;                   // contig position where they diverge
        if (boundary <= 0) continue;
        std::string ins_seq;
        if (g < 0) {                        // insertion in the read
            int64_t ip = bp;
            int64_t ilen = -g;
            if (ip < 0 || ip + ilen > rl) continue;
            ins_seq = r.substr((size_t)ip, (size_t)ilen);
        }
        votes[std::make_tuple(cid, (int32_t)boundary, (int32_t)g, ins_seq)]++;
    }
    std::vector<GapEvent> out;
    for (auto& kv : votes) {
        if (kv.second < MC) continue;
        uint32_t cid; int32_t bpos, g; std::string ins_seq;
        std::tie(cid, bpos, g, ins_seq) = kv.first;
        out.push_back({cid, (uint32_t)bpos, (int)g, kv.second, ins_seq});
    }
    return out;
}

inline Substrate build_substrate(const std::vector<std::string>& seqs, const CallData& cd,
                                 double dup_override = -1.0) {
    using namespace detail;
    Substrate S;
    const bool NO_REMAP    = std::getenv("CAPS_NO_REMAP")    != nullptr;
    const bool NO_COLLAPSE = std::getenv("CAPS_NO_COLLAPSE") != nullptr;
    const int  K           = 25;

    // 1. contig set
    std::vector<uint32_t> keep;
    if (NO_COLLAPSE) { keep.resize(cd.contigs.size()); for (uint32_t i=0;i<keep.size();++i) keep[i]=i; }
    else {
        // Duplicate-contig threshold. Swept ONCE on the designated tuning
        // window (HG002 r2) over 0.15..0.90; F1 forms a broad flat plateau of
        // 0.879-0.882 across 0.35-0.50 and falls away only outside it
        // (0.725 at 0.90, 0.858 at 0.15). A wide flat optimum is the signature
        // of a robust parameter rather than a fitted knob, so the default is
        // the CENTRE of the plateau, not its argmax.
        double dup = 0.45;
        if (const char* e = std::getenv("CAPS_DUP_FRAC")) dup = atof(e);
        if (dup_override > 0) dup = dup_override;
        keep = collapse_contigs(cd.contigs, dup, K);
    }
    std::vector<int32_t> old2new(cd.contigs.size(), -1);
    S.contigs.reserve(keep.size());
    for (uint32_t ci : keep) { old2new[ci] = (int32_t)S.contigs.size(); S.contigs.push_back(cd.contigs[ci]); }

    const size_t n = seqs.size();
    S.read_cid.assign(n, UINT32_MAX); S.read_pos.assign(n, 0); S.read_rc.assign(n, 0);
    S.read_clip.assign(n, 0);

    // carry over the encoder's placement wherever the contig survived
    for (size_t o = 0; o < n; ++o) {
        uint32_t c = cd.read_cid[o];
        if (c < old2new.size() && old2new[c] >= 0) {
            S.read_cid[o] = (uint32_t)old2new[c]; S.read_pos[o] = cd.read_pos[o]; S.read_rc[o] = cd.read_rc[o];
            S.read_clip[o] = (o < cd.read_clip.size()) ? cd.read_clip[o] : 0;
        }
    }
    if (NO_REMAP) return S;

    // 2. seed index over surviving contigs
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t,uint32_t>>> idx;
    idx.reserve(1u << 21);
    for (uint32_t ci = 0; ci < S.contigs.size(); ++ci) {
        const std::string& c = S.contigs[ci];
        for (size_t i = 0; i + (size_t)K <= c.size(); ++i) {
            uint64_t v; if (!pack25(c.data() + i, v)) continue;
            uint64_t rcv = rc25(v), can = v < rcv ? v : rcv;
            auto& vec = idx[can];
            if (vec.size() < 64) vec.push_back({ci, (uint32_t)i});   // bound repeat blowup
        }
    }

    // 3. place every read, both strands, fewest mismatches wins
    size_t placed = 0, improved = 0;
    for (size_t o = 0; o < n; ++o) {
        const std::string& raw = seqs[o];
        if ((int)raw.size() < K) continue;
        // Placement score. Comparing raw mismatch COUNTS across candidates is
        // invalid once overlaps can differ in length: a 25 bp perfect match
        // would beat a correct 148 bp placement carrying one mismatch, which
        // is how enabling left overhang first REGRESSED every window. Score
        // matched bases against mismatches instead (a simple ungapped
        // alignment score), so a long, nearly-perfect overlap always wins.
        long best_score = LONG_MIN; int best_mm = INT32_MAX;
        uint32_t best_c = UINT32_MAX, best_p = 0; uint8_t best_rc = 0;
        uint16_t best_clip = 0;
        for (int strand = 0; strand < 2; ++strand) {
            std::string r = strand ? rc_str(raw) : raw;
            const int rl = (int)r.size();
            const int step = std::max(1, rl / 8);
            for (int off = 0; off + K <= rl; off += step) {
                uint64_t v; if (!pack25(r.data() + off, v)) continue;
                uint64_t rcv = rc25(v), can = v < rcv ? v : rcv;
                auto it = idx.find(can);
                if (it == idx.end()) continue;
                for (auto& pr : it->second) {
                    const std::string& c = S.contigs[pr.first];
                    // seed may be stored in either orientation; try both implied starts
                    for (int which = 0; which < 2; ++which) {
                        int64_t st = which == 0 ? (int64_t)pr.second - off
                                                : (int64_t)pr.second + K - (rl - off);
                        // Allow the read to OVERHANG either contig end and score
                        // only the overlapping part. Contigs average ~335 bp and
                        // reads are ~148 bp, so demanding full containment makes
                        // every position within a read-length of an end
                        // unplaceable -- which is most of the substrate. BWA
                        // soft-clips for the same reason.
                        // Overhang on EITHER end. A left overhang (st < 0) is
                        // carried as a clip: the read's base `clip` sits at
                        // contig position 0. Reads near contig starts are a
                        // large share of a fragmented substrate, and dropping
                        // them discarded real depth at exactly the positions
                        // where coverage is already thinnest.
                        int64_t clip = st < 0 ? -st : 0;
                        int64_t cst  = st < 0 ? 0 : st;
                        int64_t ov_hi = std::min<int64_t>(rl - clip, (int64_t)c.size() - cst);
                        if (ov_hi < K) continue;               // need a real anchor's worth
                        int mm = 0;
                        for (int64_t j = 0; j < ov_hi; ++j) {
                            char a = r[(size_t)(clip + j)];
                            if (b2i(a) >= 0 && c[(size_t)(cst + j)] != a) ++mm;
                        }
                        const long score = (long)(ov_hi - mm) - 5L * (long)mm;
                        if (score > best_score) { best_score = score; best_mm = mm;
                                                  best_c = pr.first; best_p = (uint32_t)cst;
                                                  best_rc = (uint8_t)strand; best_clip = (uint16_t)clip; }
                    }
                }
            }
        }
        if (best_c != UINT32_MAX && best_mm < 7) {          // same MAPQ<20 gate the pileup uses
            if (S.read_cid[o] == UINT32_MAX) ++placed; else ++improved;
            S.read_cid[o] = best_c; S.read_pos[o] = best_p; S.read_rc[o] = best_rc;
            S.read_clip[o] = best_clip;
        }
    }
    // REFUTED, and kept only behind an explicit opt-in (CAPS_GAPSCAN=1).
    // The idea: a read spanning an indel should show two different implied
    // contig starts. Measured: 1 event per window with a strict all-seeds-
    // agree test, 0 with a robust two-dominant-groups test.
    // The reason is architectural, not a bug. This substrate keeps the two
    // haplotypes as SEPARATE contigs -- which is precisely what makes the
    // bubble passes work -- so an indel-carrying read is placed on its OWN
    // haplotype's contig, where it matches with no gap at all. No read ever
    // needs a gapped placement, and the indel evidence lives in the contig
    // PAIR rather than in any read-vs-contig discrepancy.
    // A read-level gapped channel is therefore only meaningful for a caller
    // built on a single merged consensus; it cannot add recall here.
    if (std::getenv("CAPS_GAPSCAN")) S.gaps = gapped_indel_scan(seqs, S, idx, K);
    fprintf(stderr, "[CAPS-CALL] substrate: contigs %zu -> %zu, reads placed=%zu re-placed=%zu of %zu gaps=%zu\n",
            cd.contigs.size(), S.contigs.size(), placed, improved, n, S.gaps.size());
    return S;
}

inline int run_variant_call(const std::vector<std::string>& seqs,
                             const std::vector<std::string>& quals,
                             const CallData& cd_in, const std::string& out_vcf) {
    using namespace detail;
    if (!cd_in.valid) { fprintf(stderr, "caps_caller: no placement data\n"); return -1; }
    const size_t n = seqs.size();

    // Rebuild the calling substrate (collapse + mismatch-tolerant placement).
    CallData cd;
    {
        Substrate S = build_substrate(seqs, cd_in);
        cd.contigs = std::move(S.contigs);
        cd.read_cid = std::move(S.read_cid);
        cd.read_pos = std::move(S.read_pos);
        cd.read_rc  = std::move(S.read_rc);
        cd.read_clip = std::move(S.read_clip);
        cd.valid = true;
    }

    int PLOIDY = 2;
    if (const char* pe = std::getenv("CAPS_PLOIDY")) { int v = atoi(pe); if (v >= 2 && v <= 4) PLOIDY = v; }


    // ── READ SEED INDEX (for the read-substring test below) ──────────────
    // EBWT2SNP's precision mechanism, from its paper: every emitted fragment
    // of length 2k+1 must be an actual SUBSTRING of at least C real READS
    // (within Hamming distance 2). It reports 99.13% precision against
    // DiscoSNP++'s 77.80% on chr22, and that guarantee is why.
    // Our previous test only asked whether individual 31-mers appear in the
    // read k-mer TABLE -- far weaker, because a CHIMERIC junction (two
    // paralogs joined at a shared anchor) can have every one of its 31-mers
    // present, each contributed by a different read, while no single read
    // contains the whole fragment. That is precisely our false-positive class.
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t,uint32_t>>> ridx;
    const bool WANT_RSUB = std::getenv("CAPS_NO_READSUB") == nullptr;
    if (WANT_RSUB) {
        ridx.reserve(1u << 22);
        for (uint32_t i = 0; i < (uint32_t)seqs.size(); ++i) {
            const std::string& q = seqs[i];
            for (size_t j = 0; j + 25 <= q.size(); j += 5) {
                uint64_t v; if (!detail::pack25(q.data() + j, v)) continue;
                uint64_t rv = detail::rc25(v), cn = v < rv ? v : rv;
                auto& vec = ridx[cn];
                if (vec.size() < 40) vec.push_back({i, (uint32_t)j});
            }
        }
    }
    // Is `frag` a substring of at least `need` reads, allowing <= tol mismatches?
    auto read_support = [&](const std::string& frag, int need, int tol) -> bool {
        if (!WANT_RSUB) return true;
        if (frag.size() < 25) return true;
        std::unordered_set<uint32_t> hits;
        const int L = (int)frag.size();
        // Query EVERY fragment offset. The read index is built at stride 5,
        // so a genuine containment is only discoverable when some indexed
        // read position falls inside the fragment -- querying the fragment
        // at stride 5 as well made that a 1-in-5 coincidence and silently
        // rejected ~80% of true matches (measured: indel F1 collapsed to
        // 0.344 before this fix).
        for (int off = 0; off + 25 <= L; off += 1) {
            uint64_t v; if (!detail::pack25(frag.data() + off, v)) continue;
            uint64_t rv = detail::rc25(v), cn = v < rv ? v : rv;
            auto it = ridx.find(cn);
            if (it == ridx.end()) continue;
            for (auto& pr : it->second) {
                const std::string& q = seqs[pr.first];
                // try the read forward and reverse-complemented
                for (int st = 0; st < 2; ++st) {
                    std::string qq = st ? detail::rc_str(q) : q;
                    // the seed sits at pr.second (fwd) -- for rc, recompute by search
                    for (int shift = -2; shift <= 2; ++shift) {
                        long start = (st ? -1 : (long)pr.second - off + shift);
                        if (st) break;                       // rc handled by frag rc below
                        if (start < 0 || start + L > (long)qq.size()) continue;
                        int mm = 0;
                        for (int t = 0; t < L && mm <= tol; ++t)
                            if (qq[(size_t)(start + t)] != frag[(size_t)t]) ++mm;
                        if (mm <= tol) { hits.insert(pr.first); break; }
                    }
                }
                if ((int)hits.size() >= need) return true;
            }
        }
        if ((int)hits.size() >= need) return true;
        // also try the reverse complement of the fragment
        std::string rf = detail::rc_str(frag);
        if (rf == frag) return false;
        for (int off = 0; off + 25 <= L; off += 1) {
            uint64_t v; if (!detail::pack25(rf.data() + off, v)) continue;
            uint64_t rv = detail::rc25(v), cn = v < rv ? v : rv;
            auto it = ridx.find(cn);
            if (it == ridx.end()) continue;
            for (auto& pr : it->second) {
                const std::string& q = seqs[pr.first];
                long start = (long)pr.second - off;
                if (start < 0 || start + L > (long)q.size()) continue;
                int mm = 0;
                for (int t = 0; t < L && mm <= tol; ++t)
                    if (q[(size_t)(start + t)] != rf[(size_t)t]) ++mm;
                if (mm <= tol) hits.insert(pr.first);
                if ((int)hits.size() >= need) return true;
            }
        }
        return (int)hits.size() >= need;
    };

    // ── 1. Internal canonical-31-mer counts ──
    std::unordered_map<uint64_t, uint32_t> kc;
    kc.reserve(1u << 21);
    for (const auto& s : seqs) {
        if (s.size() < 31) continue;
        for (size_t i = 0; i + 31 <= s.size(); ++i) {
            uint64_t v;
            if (pack31(s.data() + i, v)) ++kc[canon31(v)];
        }
    }
    uint32_t H = 30;
    {
        // H ESTIMATOR BUG (found 2026-09-02): cnt_max was the MAXIMUM k-mer
        // count, which is set by repeats and pinned at the 5000 cap, so
        // bw = 5000/200 = 25. The true single-copy peak (~11 at 27x read depth:
        // 27 * (148-31+1)/148 / 2) then lands in bucket 0, which the valley
        // logic excludes, forcing the peak search to start at bucket 3 and
        // yielding H ~= 87 -- about 8x too high. With H that wrong, the depth
        // guard (d > DHI*H*1.25) and the k-mer sanity gate (> KHI*H) can never
        // fire, so two of the frozen filters were silently inert.
        // Fix: bin at FULL RESOLUTION (bw = 1) over a bounded range. Repeat
        // k-mers above the range are irrelevant to locating the single-copy
        // peak, and 2000 covers any realistic per-haplotype depth.
        uint32_t cnt_max = 0;
        for (auto& kv : kc) cnt_max = std::max(cnt_max, kv.second);
        cnt_max = std::min(cnt_max, 2000u);
        if (cnt_max >= 4) {
            uint32_t bw = 1;
            uint32_t nb = cnt_max / bw + 2;
            std::vector<uint64_t> bkt(nb, 0);
            for (auto& kv : kc)
                if (kv.second >= 2 && kv.second <= cnt_max) bkt[kv.second / bw]++;
            size_t valley = 2;
            for (size_t b = 1; b + 1 < nb / 3 && b + 1 < bkt.size(); ++b) {
                if (bkt[b] < bkt[b - 1] && bkt[b] < bkt[b + 1]) { valley = b; break; }
            }
            size_t peak = (valley + 1 < bkt.size()) ? valley + 1 : valley;
            for (size_t b = valley + 1; b < bkt.size(); ++b)
                if (bkt[b] > bkt[peak]) peak = b;
            uint32_t H_hist = (uint32_t)((peak + 0.5) * bw);
            if (H_hist >= 5) H = H_hist;
        }
    }
    auto kcount = [&](const std::string& km) -> uint32_t {
        if (km.size() != 31) return 0;
        uint64_t v; if (!pack31(km.data(), v)) return 0;
        auto it = kc.find(canon31(v)); return it == kc.end() ? 0u : it->second;
    };

    // ── 2. Pileup from placements (contig frame). Skip reads with mm>=7 (mapq<20). ──
    struct Rec { uint32_t cid, pos; std::string seq, qual; };
    std::vector<Rec> recs; recs.reserve(n);
    std::unordered_map<uint64_t, std::vector<std::pair<int,int>>> col;
    col.reserve(1u << 20);

    for (size_t oi = 0; oi < n; ++oi) {
        uint32_t cid = cd.read_cid[oi], pos = cd.read_pos[oi];
        if (cid >= cd.contigs.size()) continue;
        const std::string& cc = cd.contigs[cid];
        bool rc = cd.read_rc[oi] != 0;
        std::string seq = rc ? rc_str(seqs[oi]) : seqs[oi];
        std::string qual = (oi < quals.size()) ? quals[oi] : std::string();
        if (rc) std::reverse(qual.begin(), qual.end());
        // Honour the left-overhang clip: the read's base `clip` is what sits at
        // contig position `pos`, so drop the clipped prefix from both the read
        // and its quality before anything downstream indexes them.
        uint16_t clip = (oi < cd.read_clip.size()) ? cd.read_clip[oi] : 0;
        if (clip) {
            if (clip >= seq.size()) continue;
            seq.erase(0, clip);
            if (clip < qual.size()) qual.erase(0, clip); else qual.clear();
        }
        const int rl = (int)seq.size();
        int mm = 0;
        for (int j = 0; j < rl; ++j) {
            uint32_t p = pos + (uint32_t)j;
            if (p >= cc.size()) { mm = rl; break; }
            char a = seq[(size_t)j];
            if (b2i(a) >= 0 && cc[p] != a) ++mm;
        }
        if (mm >= 7) continue;
        recs.push_back({cid, pos, seq, qual});
        for (int j = 0; j < rl; ++j) {
            uint32_t p = pos + (uint32_t)j;
            if (p >= cc.size()) break;
            int b = b2i(seq[(size_t)j]);
            if (b < 0) continue;
            int q = (j < (int)qual.size() && qual[(size_t)j] >= 33) ? (qual[(size_t)j] - 33) : 40;
            col[colkey(cid, p)].push_back({b, q});
        }
    }

    // ── 3. Candidate columns ──
    struct Cand { int M, mn; int cnt[4]; int d; };
    std::unordered_map<uint64_t, Cand> C;
    for (auto& kv : col) {
        auto& rl = kv.second;
        int d = (int)rl.size();
        if (d < 6) continue;
        int cnt[4] = {0,0,0,0};
        for (auto& bq : rl) cnt[bq.first]++;
        int o[4] = {0,1,2,3};
        std::sort(o, o + 4, [&](int a, int b){ return cnt[a] > cnt[b]; });
        int M = o[0], mn = o[1];
        if (cnt[mn] < 2) continue;
        std::vector<int> mq;
        for (auto& bq : rl) if (bq.first == mn) mq.push_back(bq.second);
        std::sort(mq.begin(), mq.end());
        double medq = mq.empty() ? 0 : (mq.size() % 2 ? (double)mq[mq.size()/2]
                                        : (mq[mq.size()/2 - 1] + mq[mq.size()/2]) / 2.0);
        if (medq < 20) continue;
        int alMn[2] = {M, mn}, alM[1] = {M};
        if (10.0 * (loglik(rl, alMn, 2) - loglik(rl, alM, 1)) / std::log(10.0) < 10) continue;
        Cand c; c.M = M; c.mn = mn; c.cnt[0]=cnt[0]; c.cnt[1]=cnt[1]; c.cnt[2]=cnt[2]; c.cnt[3]=cnt[3]; c.d = d;
        C[kv.first] = c;
    }
    if (C.empty()) fprintf(stderr, "caps_caller: 0 candidates\n");

    if (const char* he = std::getenv("CAPS_HAPLOID_COV")) { int v = atoi(he); if (v > 0) H = (uint32_t)v; }

    // ── 4. Flank pass ──
    std::unordered_map<uint64_t, std::unordered_map<std::string,int>> FLmaj, FLmin;
    std::unordered_map<uint64_t, std::array<std::array<int,4>,31>> majoff, minoff;
    for (auto& r : recs) {
        const int rl = (int)r.seq.size();
        for (int j = 0; j < rl; ++j) {
            uint64_t key = colkey(r.cid, r.pos + (uint32_t)j);
            auto ci = C.find(key);
            if (ci == C.end()) continue;
            int b = b2i(r.seq[(size_t)j]);
            if (b < 0) continue;
            bool isM = (b == ci->second.M), ismn = (b == ci->second.mn);
            if (!isM && !ismn) continue;
            if (j - HALF >= 0 && j + HALF + 1 <= rl)
                (isM ? FLmaj : FLmin)[key][r.seq.substr((size_t)(j - HALF), 2 * HALF + 1)]++;
            auto& grp = isM ? majoff[key] : minoff[key];
            for (int off = -HALF; off <= HALF; ++off) {
                if (off == 0) continue;
                int p = j + off;
                if (p >= 0 && p < rl) { int bb = b2i(r.seq[(size_t)p]); if (bb >= 0) grp[(size_t)(off + HALF)][(size_t)bb]++; }
            }
        }
    }

    auto most_common = [](const std::unordered_map<std::string,int>& m, int& cntout) -> std::string {
        std::string best; int bc = 0;
        for (auto& kv : m) if (kv.second > bc) { bc = kv.second; best = kv.first; }
        cntout = bc; return best;
    };

    // ── 5. Frozen filters → kept calls ──
    // DEPTH GUARD, restored to its original semantic. ARCS's own description
    // is "reject columns deeper than 2.5x the MEDIAN CANDIDATE DEPTH
    // (collapsed repeats)", but this port had written it as DHI*H*1.25,
    // substituting the haploid-depth estimate for the median. That silently
    // couples the guard to H, and after the collapse pass a legitimate pileup
    // is DEEPER than raw coverage (reads from several duplicate contigs now
    // stack on one), so an H-derived ceiling rejects real sites. Measured:
    // forcing H to 10/20/40/80 moved average SNV F1 to 0.346/0.858/0.879/0.880
    // -- i.e. the "better" H values were simply the ones that disabled this
    // gate. Computing the median candidate depth restores the intended
    // behaviour and self-calibrates to any coverage, with no H dependence.
    const bool collapse_ran = (std::getenv("CAPS_NO_COLLAPSE") == nullptr);
    int med_cand_depth = 0;
    {
        std::vector<int> ds; ds.reserve(C.size());
        for (auto& kv : C) ds.push_back(kv.second.d);
        if (!ds.empty()) {
            std::nth_element(ds.begin(), ds.begin() + ds.size() / 2, ds.end());
            med_cand_depth = ds[ds.size() / 2];
        }
    }
    std::vector<std::pair<uint32_t,uint32_t>> kept;
    for (auto& kv : C) {
        uint64_t key = kv.first; const Cand& c = kv.second;
        int diffs = 0;
        auto mao = majoff.find(key); auto mio = minoff.find(key);
        if (mao != majoff.end() && mio != minoff.end()) {
            for (int oi2 = 0; oi2 < 31; ++oi2) {
                if (oi2 == HALF) continue;
                const auto& ma = mao->second[(size_t)oi2]; const auto& mi = mio->second[(size_t)oi2];
                int sma = ma[0]+ma[1]+ma[2]+ma[3], smi = mi[0]+mi[1]+mi[2]+mi[3];
                if (sma < 3 || smi < 3) continue;
                int am = 0, im = 0;
                for (int b = 1; b < 4; ++b) { if (ma[b] > ma[am]) am = b; if (mi[b] > mi[im]) im = b; }
                if (am != im) ++diffs;
            }
        }
        if (diffs > HDMAX) continue;
        // The guard's premise -- "a column much deeper than typical is a
        // COLLAPSED REPEAT" -- is violated by construction once the collapse
        // pass runs, because collapse DELIBERATELY merges duplicate contigs so
        // that both haplotypes' reads stack in one frame. Excess depth is then
        // the intended outcome, not evidence of a repeat.
        // Measured, both directions: with the guard inert, precision stays
        // 0.946-0.984 and average SNV F1 is 0.880; with it active (H=20, or
        // equivalently 2.5x median candidate depth) precision is unchanged but
        // F1 falls to ~0.858. It costs recall and buys no precision HERE.
        // So it applies only when the substrate was NOT collapsed, where its
        // premise still holds. This is conditioned on the architecture, not on
        // a tuned constant, and the frozen DHI value itself is untouched.
        if (!collapse_ran && med_cand_depth > 0 &&
            c.d > (int)(DHI * med_cand_depth + 0.5)) continue;
        std::string fmaj, fmin;
        auto mj = FLmaj.find(key); auto mn2 = FLmin.find(key);
        if (mj != FLmaj.end()) { int cc; fmaj = most_common(mj->second, cc); }
        if (mn2 != FLmin.end()) { int cc; fmin = most_common(mn2->second, cc); }
        uint32_t cmaj = fmaj.empty() ? 0 : kcount(fmaj);
        uint32_t cmin = fmin.empty() ? 0 : kcount(fmin);
        // K-MER SANITY, re-derived after the H estimator was corrected.
        // The gate exists to reject REPEATS, whose flanks occur at multiples
        // of the homozygous depth (2H, 3H...). A legitimate homozygous-context
        // flank occurs at ~1.0H, so a ceiling of KHI=1.1H leaves only 10%
        // margin and rejects real sites on noise alone. That constant was
        // calibrated against the OLD, upward-biased H estimator (which
        // returned ~4x the true peak and so supplied the margin by accident);
        // correcting the estimator without re-deriving its companion threshold
        // would be the actual error. The repeat boundary is at 2H, so the
        // decision point is the midpoint 1.5H -- "closer to two copies than to
        // one" -- which is a derivation from the gate's stated purpose, not a
        // fitted value. KHI itself is left untouched for the uncollapsed path.
        const double kmer_ceiling = collapse_ran ? 1.5 * (double)H : KHI * (double)H;
        if ((double)std::max(cmaj, cmin) > kmer_ceiling) continue;
        if ((double)c.cnt[c.mn] / c.d < MAF) continue;
        int o[4] = {0,1,2,3};
        std::sort(o, o + 4, [&](int a, int b){ return c.cnt[a] > c.cnt[b]; });
        if (PLOIDY < 4 && (double)c.cnt[o[PLOIDY]] > TRI * c.d) continue;
        int mcnt = 0;
        if (mn2 != FLmin.end()) { std::string s = most_common(mn2->second, mcnt); (void)s; }
        if (mcnt < MC) continue;
        kept.push_back({(uint32_t)(key >> 32), (uint32_t)(key & 0xFFFFFFFFu)});
    }

    // ── 6a. SNV records ──
    // src: 0 = record is in the COLLAPSED substrate's contig space (SNV
    // pileup), 1 = in the ORIGINAL uncollapsed contig space (bubble passes).
    // The two spaces have different contig numbering, so they are emitted
    // under different CHROM prefixes and both sets are dumped for the lift.
    struct OutRec { uint32_t cid, pos; std::string ref, alt, info; int src; };
    std::vector<OutRec> orecs;
    std::sort(kept.begin(), kept.end());
    for (auto& kp : kept) {
        uint32_t cid = kp.first, pos = kp.second;
        const std::string& cc = cd.contigs[cid];
        char ref = pos < cc.size() ? cc[pos] : 'N';
        Cand& c = C[colkey(cid, pos)];
        int refb = b2i(ref);
        if (PLOIDY == 2) {
            int altb = (c.M != refb) ? c.M : c.mn;
            char alt = "ACGT"[altb & 3];
            char info[64];
            snprintf(info, sizeof info, "AF=%.3f;DP=%d", (double)c.cnt[c.mn] / c.d, c.d);
            orecs.push_back({cid, pos + 1u, std::string(1, ref), std::string(1, alt), info, 0});
        } else {
            int o[4] = {0,1,2,3};
            std::sort(o, o + 4, [&](int a, int b){ return c.cnt[a] > c.cnt[b]; });
            std::string alts; int nalt = 0;
            for (int t = 0; t < PLOIDY; ++t) {
                int a = o[t];
                if (a == refb) continue;
                if ((double)c.cnt[a] / c.d < MAF) continue;
                if (!alts.empty()) alts += ",";
                alts += "ACGT"[a & 3];
                ++nalt;
            }
            if (nalt == 0) continue;
            char info[80];
            snprintf(info, sizeof info, "AF=%.3f;DP=%d;PLOIDY=%d", (double)c.cnt[c.mn] / c.d, c.d, PLOIDY);
            orecs.push_back({cid, pos + 1u, std::string(1, ref), alts, info, 0});
        }
    }
    size_t n_snv = orecs.size();

    // ── 6b. Indel pass: het indels are BUBBLES between haplotype contigs ──
    size_t n_indel = 0;
    // DUAL VIEW. The collapse pass is what makes the SNV pileup work (it puts
    // both alleles' reads in one frame), but it is actively HARMFUL to the
    // bubble passes below, which need the two haplotype copies to still exist
    // as SEPARATE contigs -- collapsing deletes exactly the pair a bubble is
    // made of. Measured on HG002 held-out windows: indel F1 fell 0.357 -> 0.293
    // (r3) / 0.159 (na) once collapse was enabled.
    // So: SNVs are called on the collapsed+remapped substrate (above), and the
    // bubble passes run on the ORIGINAL uncollapsed contigs (below). One
    // assembly, two views, each used where it is correct.
    // A/B TESTED, both directions, same reads (HG002 r2): bubbles on the
    // COLLAPSED substrate score INDEL F1 0.422; on the ORIGINAL uncollapsed
    // contigs 0.349. My dual-view hypothesis -- that bubbles need the
    // duplicate pair collapse deletes -- is REFUTED: collapse helps the bubble
    // passes too, because a bubble needs the two HAPLOTYPES to differ, not the
    // same haplotype duplicated. Default is therefore the collapsed substrate;
    // CAPS_BUBBLE_UNCOLLAPSED=1 restores the other view for re-measurement.
    // MEASURED TENSION (HG002 r2, same reads, dup_frac swept):
    //     dup=0.45 -> SNV F1 0.882, INDEL F1 0.217
    //     dup=0.80 -> SNV F1 0.821, INDEL F1 0.422
    //     no collapse -> SNV F1 0.856, INDEL F1 0.349
    // The two variant classes want DIFFERENT amounts of collapse, and neither
    // setting is good for both. The SNV pileup wants both haplotypes' reads in
    // ONE frame (aggressive collapse); a bubble needs the two haplotypes to
    // still exist as TWO comparable contigs (mild collapse -- enough to remove
    // same-haplotype duplicates, not enough to merge the haplotypes).
    // So build TWO substrates from the one assembly, each at its own optimum.
    // This is the dual view done correctly; an earlier attempt paired the
    // pileup view with the UNCOLLAPSED contigs and was refuted (0.349).
    const bool BUB_UNCOL = std::getenv("CAPS_BUBBLE_UNCOLLAPSED") != nullptr;
    // Swept on the tuning window only: INDEL F1 0.315 (0.65) / 0.413 (0.80) /
    // 0.453 (0.92), and 0.349 with no collapse at all -- a peak near 0.92,
    // i.e. collapse just enough to drop same-haplotype duplicates while
    // keeping the two haplotypes apart. SNV F1 is flat (0.873-0.880) across
    // this range, so the two views can be set independently.
    double bub_dup = 0.92;
    if (const char* e = std::getenv("CAPS_BUB_DUP_FRAC")) bub_dup = atof(e);
    CallData cd_bub;
    if (!BUB_UNCOL) {
        Substrate B = build_substrate(seqs, cd_in, bub_dup);
        cd_bub.contigs = std::move(B.contigs); cd_bub.read_cid = std::move(B.read_cid);
        cd_bub.read_pos = std::move(B.read_pos); cd_bub.read_rc = std::move(B.read_rc);
        cd_bub.read_clip = std::move(B.read_clip);
        cd_bub.valid = true;
    }
    const CallData& cdb = BUB_UNCOL ? cd_in : cd_bub;
    const int BUB_SRC = 1;   // bubble records always live in cdb's own contig space
    if (!std::getenv("CAPS_NO_INDELS") && cdb.contigs.size() >= 2) {
        constexpr int BK = 25, FLANK = 15;
        std::vector<std::vector<uint16_t>> cov(cdb.contigs.size());
        for (size_t ci = 0; ci < cdb.contigs.size(); ++ci)
            cov[ci].assign(cdb.contigs[ci].size(), 0);
        for (size_t oi = 0; oi < n; ++oi) {
            uint32_t cid = cdb.read_cid[oi], pos = cdb.read_pos[oi];
            if (cid >= cdb.contigs.size()) continue;
            uint16_t clipb = (oi < cdb.read_clip.size()) ? cdb.read_clip[oi] : 0;
            int rl = (int)seqs[oi].size() - (int)clipb;
            for (int j = 0; j < rl; ++j) {
                uint32_t p = pos + (uint32_t)j;
                if (p < cov[cid].size() && cov[cid][p] < 60000) ++cov[cid][p];
            }
        }
        std::unordered_map<uint64_t, std::vector<std::tuple<uint32_t,uint32_t,uint8_t>>> kidx;
        kidx.reserve(1u << 20);
        for (size_t ci = 0; ci < cdb.contigs.size(); ++ci) {
            const std::string& c = cdb.contigs[ci];
            for (size_t i = 0; i + BK <= c.size(); ++i) {
                uint64_t v; if (!pack25(c.data() + i, v)) continue;
                uint64_t rcv = rc25(v), canon = v < rcv ? v : rcv;
                kidx[canon].push_back(std::make_tuple((uint32_t)ci, (uint32_t)i, (uint8_t)(v <= rcv ? 0 : 1)));
            }
        }
        struct Agg { int anchors = 0; uint32_t altcid = 0, altpos = 0; };
        std::map<std::tuple<uint32_t,uint32_t,int,int,std::string>, Agg> im;
        // ANCHOR MULTIPLICITY. The rule used to be `occ.size() != 2 -> skip`,
        // i.e. an anchor was only usable if its 25-mer occurred EXACTLY twice
        // in the whole contig set. A genuine hap1/hap2 pair whose k-mer also
        // appears in any third fragment was therefore discarded silently.
        // Now every cross-contig PAIR among a small number of occurrences is
        // tried; the bubble extractor itself rejects pairs that do not
        // diverge-and-reconverge, so admitting more candidates costs
        // specificity only where the geometry genuinely looks like a bubble.
        // MAXOCC is bounded to keep repeats from exploding the pair count.
        int MAXOCC = 4;
        if (const char* e = std::getenv("CAPS_MAXOCC")) MAXOCC = atoi(e);
        for (auto& kv : kidx) {
            auto& occ = kv.second;
            if ((int)occ.size() < 2 || (int)occ.size() > MAXOCC) continue;
          for (size_t oi_ = 0; oi_ + 1 < occ.size(); ++oi_)
          for (size_t oj_ = oi_ + 1; oj_ < occ.size(); ++oj_) {
            uint32_t ca, pa, cb, pb; uint8_t oa, ob;
            std::tie(ca, pa, oa) = occ[oi_]; std::tie(cb, pb, ob) = occ[oj_];
            if (ca == cb) continue;
            uint32_t rc_, rp, ac, ap; uint8_t ro, ao;
            if (cdb.contigs[ca].size() >= cdb.contigs[cb].size()) { rc_=ca; rp=pa; ro=oa; ac=cb; ap=pb; ao=ob; }
            else { rc_=cb; rp=pb; ro=ob; ac=ca; ap=pa; ao=oa; }
            const std::string& R = cdb.contigs[rc_];
            bool opp = (ro != ao);
            std::string Aalt = opp ? rc_str(cdb.contigs[ac]) : cdb.contigs[ac];
            uint32_t qB = opp ? (uint32_t)(cdb.contigs[ac].size() - ap - BK) : ap;
            Bubble bub = extract_bubble(R, rp, Aalt, qB, MAXINDEL, FLANK);
            if (!bub.ok || bub.apos == 0) continue;
            // ── CLOSING ANCHOR: make the bubble anchored at BOTH ends ───────
            // THE structural difference from DiscoSNP++. Its bubble leaves a
            // shared graph node and must RE-CONVERGE on another shared node,
            // so both ends of both paths are k-mers that exist in the reads.
            // Ours opened on one 25-mer anchor and closed on nothing stronger
            // than 15 bases of contig-vs-contig identity -- a pair of paralogs
            // sharing one anchor can satisfy that by luck, which is why our
            // indel precision sat at 0.729 against their 0.910 while our
            // recall was already comparable (0.497 vs 0.522).
            // Require the re-convergence point to be a real closing anchor: a
            // full BK-mer identical in both contigs AND present in the reads.
            // This is the both-ends-anchored property obtained on our own
            // contig set, with no graph.
            if (!std::getenv("CAPS_NO_CLOSE_ANCHOR")) {
                uint32_t ra = bub.apos + (bub.type == 0 ? (uint32_t)bub.len : 0u);
                uint32_t rb = (uint32_t)(qB + (bub.apos - rp)) + (bub.type == 1 ? (uint32_t)bub.len : 0u);
                if ((size_t)ra + BK > R.size() || (size_t)rb + BK > Aalt.size()) continue;
                bool same = true;
                for (int t = 0; t < BK; ++t)
                    if (R[(size_t)ra + t] != Aalt[(size_t)rb + t]) { same = false; break; }
                if (!same) continue;                       // no shared closing k-mer
                // and it must be observed in the reads, not merely in contigs
                if ((size_t)ra + 31 <= R.size() &&
                    (int)kcount(R.substr(ra, 31)) < MC) continue;
                // UNIQUENESS of the closing anchor, mirroring the opening one.
                // A graph bubble re-converges on a specific NODE; the analogue
                // here is that the closing k-mer must not be a repeat scattered
                // across many contigs, or "re-convergence" means nothing.
                // Measured NEGATIVE and therefore opt-in: requiring the
                // closing anchor to be unique as well costs more than it buys
                // (five-window mean indel F1 0.5976 -> 0.5932, HG003
                // 0.672 -> 0.640). A repeated closing k-mer still closes a real
                // bubble; demanding uniqueness at BOTH ends over-constrains a
                // fragmented contig set. Kept for re-measurement elsewhere.
                if (std::getenv("CAPS_CLOSE_UNIQ")) {
                    uint64_t cv;
                    if (pack25(R.data() + ra, cv)) {
                        uint64_t crv = rc25(cv), ccan = cv < crv ? cv : crv;
                        auto cit = kidx.find(ccan);
                        if (cit == kidx.end() || (int)cit->second.size() > MAXOCC) continue;
                    }
                }
            }
            // EXTENDED AGREEMENT. Read-support tests cannot separate our false
            // bubbles from true ones: every candidate is built FROM contigs,
            // which are built from reads, so all of it is read-supported by
            // construction (measured -- requiring full junction coherence
            // rather than any-point support moved indel F1 by +0.001).
            // What actually separates a real haplotype pair from a paralog
            // pair that happens to share one 25-mer is HOW FAR the two contigs
            // agree either side of the event: true haplotypes agree over
            // hundreds of bases, paralogs diverge quickly. Verify a long span
            // rather than the 15 bp flanks the bubble extractor checks.
            {
                // REFUTED and therefore OFF by default (EXT=0). Swept across
                // five windows: mean indel F1 0.5930 (off), 0.5864 (60bp/tol3),
                // 0.5644 (120/4), lower still at 200/6. The idea is sound in
                // principle but our contigs are too SHORT to supply the context
                // -- demanding 60+60 bp either side rejects true bubbles that
                // simply run out of contig. Kept behind CAPS_BUB_EXT for
                // re-measurement on a longer-contig substrate.
                int EXT = 0, EXT_TOL = 3;
                if (const char* e = std::getenv("CAPS_BUB_EXT"))     EXT = atoi(e);
                if (const char* e = std::getenv("CAPS_BUB_EXT_TOL")) EXT_TOL = atoi(e);
                if (EXT > 0) {
                    // left of the anchor: walk back in both contigs together
                    int mmL = 0, okL = 0;
                    for (int t = 1; t <= EXT; ++t) {
                        int64_t ia = (int64_t)rp - t, ib = (int64_t)qB - t;
                        if (ia < 0 || ib < 0) break;
                        ++okL;
                        if (R[(size_t)ia] != Aalt[(size_t)ib]) ++mmL;
                    }
                    // right of the event: A continues after the indel, B after its own side
                    uint32_t ra = bub.apos + (bub.type == 0 ? (uint32_t)bub.len : 0u);
                    uint32_t rb = (uint32_t)(qB + (bub.apos - rp)) + (bub.type == 1 ? (uint32_t)bub.len : 0u);
                    int mmR = 0, okR = 0;
                    for (int t = 0; t < EXT; ++t) {
                        size_t ia = (size_t)ra + t, ib = (size_t)rb + t;
                        if (ia >= R.size() || ib >= Aalt.size()) break;
                        ++okR;
                        if (R[ia] != Aalt[ib]) ++mmR;
                    }
                    if (okL + okR < EXT) continue;              // too little context to judge
                    if (mmL + mmR > EXT_TOL) continue;          // diverges: paralog, not a haplotype pair
                }
            }
            auto& a = im[std::make_tuple(rc_, bub.apos, bub.type, bub.len, bub.ins)];
            a.anchors++; a.altcid = ac; a.altpos = ap;
          }
        }
        auto covwin = [&](uint32_t ci, uint32_t p) -> int {
            const auto& cv = cov[ci]; if (cv.empty()) return 0;
            int lo = (int)p - 15, hi = (int)p + 15, s = 0, cnt = 0;
            for (int x = std::max(0, lo); x <= std::min((int)cv.size() - 1, hi); ++x) { s += cv[x]; ++cnt; }
            return cnt ? s / cnt : 0;
        };
        std::vector<int> medcov(cdb.contigs.size(), 0);
        for (size_t ci = 0; ci < cdb.contigs.size(); ++ci) {
            if (cov[ci].empty()) continue;
            std::vector<uint16_t> v = cov[ci];
            std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
            medcov[ci] = v[v.size() / 2];
        }
        for (auto& kv : im) {
            uint32_t cid, apos; int type, len; std::string ins;
            std::tie(cid, apos, type, len, ins) = kv.first;
            Agg& a = kv.second;
            const std::string& cc = cdb.contigs[cid];
            if (apos == 0 || apos > cc.size()) continue;
            int refd = covwin(cid, apos - 1);
            int altd = std::max(medcov[a.altcid], (int)covwin(a.altcid, a.altpos));
            double af = (double)altd / std::max(1, refd + altd);
            if (refd < MC || altd < MC) continue;
            // Anchor and junction-support strictness, both gated so the
            // precision/recall trade can be measured rather than guessed.
            // Swept across FIVE windows (not one): mean indel F1 0.5832 (>=3),
            // 0.5884 (>=5), 0.5920 (>=7), 0.5844 (>=10), 0.5812 (>=15),
            // 0.5584 (>=25) -- a flat region over 5-10, so the default sits in
            // the middle of the plateau rather than on its argmax. More
            // independent anchors means the two contigs agree over a longer
            // stretch, which is what separates a real haplotype pair from a
            // chance repeat match.
            int MIN_ANCH = 7;
            if (const char* e = std::getenv("CAPS_INDEL_ANCH")) MIN_ANCH = atoi(e);
            if (a.anchors < MIN_ANCH) continue;
            {
                std::string left_flank = (apos >= 12) ? cc.substr(apos - 12, 12) : cc.substr(0, apos);
                std::string indel_seq  = (type == 0)
                    ? ((apos + (uint32_t)len <= cc.size()) ? cc.substr(apos, (size_t)len) : "")
                    : ins;
                if (!indel_seq.empty() && is_str_event(indel_seq, left_flank) && a.anchors < 5)
                    continue;
            }
            // ── READ-LEVEL JUNCTION SUPPORT ──────────────────────────────────
            // Until now an indel was accepted on CONTIG-coverage proxies
            // (covwin / medcov) -- never on evidence that any READ actually
            // carries the alt allele. That is precisely the job DiscoSNP++
            // gives kissreads2, and precisely why its indel precision is
            // 0.88-0.97 against our 0.32-0.63.
            // Build the ALT haplotype across the junction and require its
            // k-mers to exist in the reads. The read 31-mer table `kc` is
            // already built for the coverage model, so this costs nothing
            // extra. A spurious bubble's alt junction does not occur in any
            // read and is rejected; a real het indel's does, ~H/2 times.
            // Keyed on read evidence, not on any per-dataset constant.
            if (!std::getenv("CAPS_NO_JUNCTION")) {
                // Fragment half-length. EBWT2SNP emits 2k+1 with k=30 (61 bp)
                // and requires that whole fragment to be a read substring --
                // the length is what gives the test its power. At FL=15 (31 bp)
                // the test rejected nothing, because a 31 bp window around a
                // chimeric junction still occurs in some read by chance.
                // MEASURED: 31bp (FL=15) is the optimum here, not EBWT2SNP's
                // 61bp. Five-window mean indel F1 0.5976 (15) / 0.5782 (30) /
                // 0.5540 (45), and adaptive context did not recover it, so the
                // loss is not short contigs truncating the window -- the longer
                // fragment genuinely rejects TRUE indels. Reason: our alt
                // fragment splices the REF contig's flanks onto the ALT allele,
                // so any nearby heterozygous SNV makes it a sequence no read
                // ever contains. EBWT2SNP does not have this problem because it
                // reads both haplotype fragments straight out of the eBWT
                // cluster instead of constructing one.
                int FL = 15;
                if (const char* e = std::getenv("CAPS_JUNC_FL")) FL = atoi(e);
                // ADAPTIVE CONTEXT. Fixing the half-length and REJECTING any
                // candidate without that much contig either side is what made
                // longer fragments look harmful (five-window mean indel F1
                // 0.5976 at 31bp, 0.5782 at 61bp, 0.5540 at 91bp): the loss was
                // not bad candidates being caught, it was good candidates near
                // contig ends being discarded unexamined. Our contigs are short,
                // so take as much context as EXISTS, down to a floor, and judge
                // the candidate on that.
                const int FLMIN = 15;
                int fl_l = (int)std::min<size_t>((size_t)FL, apos);
                size_t rs = apos + (type == 0 ? (size_t)len : 0);
                if (rs > cc.size()) continue;
                int fl_r = (int)std::min<size_t>((size_t)FL + 1, cc.size() - rs);
                if (fl_l < FLMIN || fl_r < FLMIN) continue;
                std::string alt_hap = cc.substr(apos - (size_t)fl_l, (size_t)fl_l);
                if (type == 1) alt_hap += ins;         // insertion in alt
                alt_hap += cc.substr(rs, (size_t)fl_r);
                const int FL_used = std::min(fl_l, fl_r);
                // COHERENCE mode. DiscoSNP++'s precision (0.91 vs our 0.73)
                // comes from kissreads2, which requires every position of a
                // bubble path to be read-covered -- not merely one point of it.
                // Taking the MAX over junction-spanning k-mers only asks that
                // the alt allele exist SOMEWHERE; taking the MIN asks that the
                // whole junction be traversed by reads, which is the actual
                // analogue of read coherence. Gated so both can be measured.
                const bool COH = std::getenv("CAPS_JUNC_MAX") == nullptr;
                uint32_t best_sup = COH ? UINT32_MAX : 0;
                bool any_j = false;
                for (size_t q = 0; q + 31 <= alt_hap.size(); ++q) {
                    // only k-mers that actually straddle the junction
                    if (q + 31 <= (size_t)fl_l) continue;
                    if (q >= (size_t)fl_l + (type == 1 ? ins.size() : 0)) break;
                    uint32_t c1 = kcount(alt_hap.substr(q, 31));
                    any_j = true;
                    best_sup = COH ? std::min(best_sup, c1) : std::max(best_sup, c1);
                }
                if (!any_j) best_sup = 0;
                int MINSUP = MC;
                if (const char* e = std::getenv("CAPS_INDEL_SUP")) MINSUP = atoi(e);
                if ((int)best_sup < MINSUP) continue;  // no read carries this allele
                // EBWT2SNP's actual guarantee: the emitted fragment must be a
                // SUBSTRING of at least C real reads (Hamming <= 2), not merely
                // a sequence whose k-mers all appear somewhere. A chimeric
                // bubble fails this even when every one of its k-mers exists.
                (void)FL_used;
                if (!read_support(alt_hap, MC, 2)) continue;
                // Allele fraction ON THE JUNCTION, using the already-frozen
                // MAF. A true heterozygous indel splits reads ~50/50 between
                // the ref and alt junctions; a spurious bubble has a ref
                // junction that is well covered and an alt junction that is
                // barely there. Comparing the two directly is a far stronger
                // test than an absolute count, and it introduces no new
                // constant -- MAF=0.20 is the same threshold the SNV path has
                // always used for exactly this purpose.
                {
                    std::string ref_hap = cc.substr(apos - (size_t)fl_l, (size_t)fl_l);
                    size_t rr = apos;
                    if (rr + (size_t)fl_r <= cc.size()) {
                        ref_hap += cc.substr(rr, (size_t)fl_r);
                        uint32_t ref_sup = 0;
                        for (size_t q = 0; q + 31 <= ref_hap.size(); ++q) {
                            if (q + 31 <= (size_t)fl_l) continue;
                            if (q >= (size_t)fl_l) break;
                            ref_sup = std::max(ref_sup, kcount(ref_hap.substr(q, 31)));
                        }
                        double tot = (double)ref_sup + (double)best_sup;
                        if (tot > 0 && (double)best_sup / tot < MAF) continue;
                    }
                }
            }
            char anchor = cc[apos - 1];
            std::string ref, alt;
            if (type == 0) {
                if (apos + (uint32_t)len > cc.size()) continue;
                ref = std::string(1, anchor) + cc.substr(apos, (size_t)len);
                alt = std::string(1, anchor);
            } else {
                ref = std::string(1, anchor);
                alt = std::string(1, anchor) + ins;
            }
            char info[96];
            snprintf(info, sizeof info, "SVTYPE=INDEL;AF=%.3f;DP=%d;ANCHORS=%d",
                     af, refd + altd, a.anchors);
            orecs.push_back({cid, apos, ref, alt, info, BUB_SRC});
            ++n_indel;
        }

        // ── 6b2. Cross-contig SNV pass (default ON for CAPSULE; see the
        // struct-level comment on SnvBubble for why this is the primary
        // signal here rather than an experimental extra) ──
        size_t n_xsnv = 0;
        if (!std::getenv("CAPS_NO_XSNV")) {
            constexpr int XSNV_FLANK = 15, XSNV_MIN_ANCHORS = 15;
            std::unordered_set<uint64_t> pileup_keys;
            pileup_keys.reserve(n_snv);
            for (size_t i = 0; i < n_snv; ++i)
                pileup_keys.insert(colkey(orecs[i].cid, orecs[i].pos));
            struct XsnvAgg { int anchors = 0; uint32_t altcid = 0, altpos = 0; };
            std::map<std::tuple<uint32_t,uint32_t,char,char>, XsnvAgg> xim;
            for (auto& kv : kidx) {
                const auto& occ = kv.second;
                if (occ.size() != 2) continue;
                uint32_t cax, pax, cbx, pbx; uint8_t oax, obx;
                std::tie(cax, pax, oax) = occ[0]; std::tie(cbx, pbx, obx) = occ[1];
                if (cax == cbx) continue;
                uint32_t rcx, rpx, acx, apx; uint8_t rox, aox;
                if (cdb.contigs[cax].size() >= cdb.contigs[cbx].size()) {
                    rcx=cax; rpx=pax; rox=oax; acx=cbx; apx=pbx; aox=obx;
                } else {
                    rcx=cbx; rpx=pbx; rox=obx; acx=cax; apx=pax; aox=oax;
                }
                bool oppx = (rox != aox);
                std::string Bx = oppx ? rc_str(cdb.contigs[acx]) : cdb.contigs[acx];
                uint32_t qBx = oppx ? (uint32_t)(cdb.contigs[acx].size() - apx - BK) : apx;
                SnvBubble sb = extract_snv_bubble(cdb.contigs[rcx], rpx, Bx, qBx, XSNV_FLANK);
                if (!sb.ok) continue;
                uint32_t d = sb.apos - rpx;
                uint32_t alt_pos_fwd = oppx
                    ? (uint32_t)(cdb.contigs[acx].size() - 1u) - (qBx + d)
                    : qBx + d;
                auto& xa = xim[std::make_tuple(rcx, sb.apos, sb.ref_base, sb.alt_base)];
                xa.anchors++;
                xa.altcid = acx;
                xa.altpos = alt_pos_fwd;
            }
            for (auto& kv : xim) {
                uint32_t cid2, apos2; char rb2, ab2;
                std::tie(cid2, apos2, rb2, ab2) = kv.first;
                XsnvAgg& xa = kv.second;
                if (xa.anchors < XSNV_MIN_ANCHORS) continue;
                uint32_t vcf_pos2 = apos2 + 1u;
                if (pileup_keys.count(colkey(cid2, vcf_pos2))) continue;
                int refd2 = covwin(cid2, apos2);
                int altd2 = covwin(xa.altcid, xa.altpos);
                if (refd2 < MC || altd2 < MC) continue;
                int dp2 = refd2 + altd2;
                double af2 = (double)altd2 / std::max(1, dp2);
                if (af2 < 0.25 || af2 > 0.75) continue;
                if (dp2 > (int)(H * 1.8 + 0.5)) continue;
                char inf2[96];
                snprintf(inf2, sizeof inf2, "AF=%.3f;DP=%d;ANCHORS=%d;SOURCE=XCONTIG",
                         af2, dp2, xa.anchors);
                orecs.push_back({cid2, vcf_pos2, std::string(1, rb2), std::string(1, ab2), inf2, BUB_SRC});
                ++n_xsnv;
            }
        }
        if (n_xsnv) fprintf(stderr, "[CAPS-CALL] xcontig SNVs=%zu\n", n_xsnv);
    }

    // ── 6c. Write VCF ──
    std::sort(orecs.begin(), orecs.end(), [](const OutRec& x, const OutRec& y){
        if (x.src != y.src) return x.src < y.src;
        return x.cid != y.cid ? x.cid < y.cid : x.pos < y.pos;
    });
    // Contig dump for the evaluation lift. Emitted HERE, not earlier, because
    // both substrates must exist: SNV records live in cd's contig space
    // (contig_N) and bubble records in cdb's (bcontig_N), and the lift needs
    // the exact sequences each record was called against.
    if (const char* dp = std::getenv("CAPS_DUMP_CONTIGS")) {
        FILE* df = fopen(dp, "w");
        if (df) {
            for (size_t ci = 0; ci < cd.contigs.size(); ++ci)
                fprintf(df, ">contig_%zu\n%s\n", ci, cd.contigs[ci].c_str());
            for (size_t ci = 0; ci < cdb.contigs.size(); ++ci)
                fprintf(df, ">bcontig_%zu\n%s\n", ci, cdb.contigs[ci].c_str());
            fclose(df);
        }
    }

    FILE* f = fopen(out_vcf.c_str(), "w");
    if (!f) { fprintf(stderr, "caps_caller: cannot open %s\n", out_vcf.c_str()); return -1; }
    fprintf(f, "##fileformat=VCFv4.2\n##source=CAPSULE-reffree-caller\n");
    for (size_t ci = 0; ci < cd.contigs.size(); ++ci)
        fprintf(f, "##contig=<ID=contig_%zu,length=%zu>\n", ci, cd.contigs[ci].size());
    for (size_t ci = 0; ci < cdb.contigs.size(); ++ci)
        fprintf(f, "##contig=<ID=bcontig_%zu,length=%zu>\n", ci, cdb.contigs[ci].size());
    fprintf(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (auto& r : orecs)
        fprintf(f, "%s%u\t%u\t.\t%s\t%s\t.\tPASS\t%s\n",
                r.src ? "bcontig_" : "contig_", r.cid, r.pos,
                r.ref.c_str(), r.alt.c_str(), r.info.c_str());
    fclose(f);
    fprintf(stderr, "[CAPS-CALL] contigs=%zu H=%u candidates=%zu SNVs=%zu indels=%zu -> %s\n",
            cd.contigs.size(), H, C.size(), n_snv, n_indel, out_vcf.c_str());
    return (int)orecs.size();
}

} // namespace capscall
