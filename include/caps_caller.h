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
#include <cstring>
#include <chrono>
#include <queue>
#include <deque>
#include <set>
#include <functional>
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

    // ── HOMOPOLYMER RUN-LENGTH BRANCH ────────────────────────────────────────
    // MEASURED MOTIVATION (docs/INDEL_LOSS_SKELETAL.md): every one of the 9
    // truth indels DiscoSNP++ finds on HG002 r2 that we miss is a 1 bp indel
    // inside a homopolymer run of 7-8 identical bases. The generic path below
    // cannot represent them: a 1 bp change in an 8-mer run shifts its LENGTH
    // (8 -> 7 or 9), the two haplotypes stay identical through the whole run,
    // and the byte-identical `flank_match` never succeeds at any offset g
    // because the flanks are themselves shifted by the length difference.
    // DiscoSNP++ resolves this by extending one graph path base-by-base and
    // taking the shortest closing extension, plus an explicit ambiguity
    // allowance (max_ambigous_indel, default 20).
    // Here the equivalent is direct: when the divergence sits at a
    // homopolymer, measure the RUN LENGTH on each side and emit the difference
    // as the indel, then verify re-convergence AFTER both runs.
    if (!std::getenv("CAPS_NO_HPBUBBLE")) {
        size_t ia = pA + d, ib = qB + d;
        // the run character is the one repeating just before the divergence
        char rc0 = A[ia - 1];
        if (rc0 == B[ib - 1] && b2i(rc0) >= 0) {
            // length of the run ending at the divergence (shared prefix part)
            size_t back = 0;
            while (back < ia && back < ib && A[ia - 1 - back] == rc0 && B[ib - 1 - back] == rc0) ++back;
            // how much further each side continues the same character
            size_t ea = 0, eb = 0;
            while (ia + ea < la && A[ia + ea] == rc0) ++ea;
            while (ib + eb < lb && B[ib + eb] == rc0) ++eb;
            long diff = (long)ea - (long)eb;        // A longer => deletion in B
            if (std::getenv("CAPS_HPDBG") && back >= 3)
                fprintf(stderr, "[hp] run='%c' back=%zu ea=%zu eb=%zu diff=%ld\n", rc0, back, ea, eb, (long)ea-(long)eb);
            if (back >= 4 && diff != 0 && std::labs(diff) <= maxindel) {
                size_t ra = ia + ea, rb = ib + eb;  // first base after each run
                bool fm = (ra + (size_t)FLANK <= la && rb + (size_t)FLANK <= lb &&
                           flank_match(A, ra, B, rb, FLANK));
                if (std::getenv("CAPS_HPDBG"))
                    fprintf(stderr, "[hp]   candidate diff=%ld flank_ok=%d\n", diff, (int)fm);
                if (fm) {
                    if (diff > 0) {                 // A has extra copies: deletion in B
                        r.type = 0; r.apos = (uint32_t)ia; r.len = (int)diff;
                    } else {                        // B has extra copies: insertion in B
                        r.type = 1; r.apos = (uint32_t)ia; r.len = (int)(-diff);
                        r.ins.assign((size_t)(-diff), rc0);
                    }
                    r.ok = true; return r;
                }
            }
        }
    }

    // ── SHORT TANDEM REPEAT (STR) RUN-LENGTH BRANCH ──────────────────────────
    // Generalizes the homopolymer branch above from unit length 1 to unit
    // lengths 2..6. MEASURED MOTIVATION (docs/HET_INDEL_FRESH_SCAN.md, 2026-09-03):
    // a fresh scan of the raw FP/FN records found a real disagreement at
    // chr20:3332481-3332512, a (TTTA)n tetranucleotide repeat that GIAB's own
    // truth annotates `difficultregion=AllTandemRepeats_lt51bp_slop5`. The
    // homopolymer branch above cannot reach it: its guard requires a SINGLE
    // repeating nucleotide (`rc0 == B[ib-1] && b2i(rc0) >= 0`), so a 4 bp
    // repeat UNIT never qualifies. The generic flank_match loop below cannot
    // reach it either, for exactly the reason the homopolymer comment already
    // gives: when a whole repeat unit is gained or lost, the flanks on the two
    // sides are themselves shifted by the length difference, so a
    // byte-identical flank match fails at every offset g.
    //
    // The fix is the same run-length measurement the homopolymer branch already
    // uses and that is already validated -- count how far each side continues
    // the periodic unit and emit the difference -- just tiled by a k-mer unit
    // instead of a single character. diff is computed from the run lengths, so
    // its SIGN (insertion vs deletion) is correct by construction, unlike the
    // generic loop, which accepts whichever g happens to satisfy flank_match
    // first and has no length/direction guarantee inside a repeat.
    //
    // Deliberately starts at U=2: U=1 is the homopolymer branch's territory and
    // is left completely untouched so this cannot regress that already-tuned
    // path. Requires >=2 unit copies before the divergence so this only fires
    // at genuine repeats, not at coincidental k-mer equality.
    //
    // ── REFUTED BY MEASUREMENT, 2026-09-03. OPT-IN (CAPS_STRBUBBLE=1), OFF by
    // default. The reasoning above is sound and the branch DOES engage --
    // instrumented at 682 firings on HG002 r2 -- but the caller's output is
    // BYTE-IDENTICAL with it on and off (`cmp` on calls.vcf), and the scored
    // result is unchanged to the digit (INDEL P=0.704 R=0.576 F1=0.633 both
    // ways). So the STR loci it reaches were already being resolved
    // equivalently by the generic loop, or its candidates die in the same
    // downstream filters every other candidate dies in -- consistent with the
    // standing structural finding (docs/HOW_DISCOSNP_WINS.md sec 4) that the
    // indel gap is not in bubble GEOMETRY but in what the substrate can offer
    // a filter to reject. Kept behind a flag, not deleted, per this project's
    // rule that refuted ideas stay on the record with their measurement.
    if (std::getenv("CAPS_STRBUBBLE")) {
        const size_t ia = pA + d, ib = qB + d;
        for (size_t U = 2; U <= 6; ++U) {
            if (ia < 2 * U || ib < 2 * U) continue;
            const std::string unit = A.substr(ia - U, U);
            if (B.compare(ib - U, U, unit) != 0) continue;          // same unit both sides
            if (A.compare(ia - 2 * U, U, unit) != 0) continue;      // >=2 copies before divergence
            if (B.compare(ib - 2 * U, U, unit) != 0) continue;
            size_t ea = 0, eb = 0;
            while (ia + ea + U <= la && A.compare(ia + ea, U, unit) == 0) ea += U;
            while (ib + eb + U <= lb && B.compare(ib + eb, U, unit) == 0) eb += U;
            const long diff = (long)ea - (long)eb;
            if (diff == 0 || std::labs(diff) > maxindel) continue;
            const size_t ra = ia + ea, rb = ib + eb;
            if (ra + (size_t)FLANK > la || rb + (size_t)FLANK > lb) continue;
            if (!flank_match(A, ra, B, rb, FLANK)) continue;
            if (std::getenv("CAPS_HPDBG"))
                fprintf(stderr, "[str] unit='%s' U=%zu ea=%zu eb=%zu diff=%ld\n",
                        unit.c_str(), U, ea, eb, diff);
            if (diff > 0) {                     // A has extra copies: deletion in B
                r.type = 0; r.apos = (uint32_t)ia; r.len = (int)diff;
            } else {                            // B has extra copies: insertion in B
                r.type = 1; r.apos = (uint32_t)ia; r.len = (int)(-diff);
                r.ins = B.substr(ib, (size_t)(-diff));
            }
            r.ok = true; return r;
        }
    }

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

// ── Whole-pair gapped alignment scan ─────────────────────────────────────────
// MEASURED MOTIVATION (docs/INDEL_LOSS_SKELETAL.md §8): both allele-contigs
// DO exist (84% of placed reads match their contig perfectly), but the two are
// never paired at the event, because near a homopolymer every 25-mer
// overlapping the run differs between haplotypes, so the shared anchors that
// survive sit 213-135,453 bp away. Walking outward from such an anchor dies on
// ~335 bp contigs.
// This removes both constraints at once: take the diagonal implied by ANY
// shared anchor and scan the ENTIRE co-linear region, resolving each divergence
// with a gap and continuing past it on the updated diagonal. Events anywhere
// along the pair are then reachable from one distant anchor.
inline std::vector<Bubble> scan_pair(const std::string& A, uint32_t pA,
                                     const std::string& B, uint32_t qB,
                                     int maxindel, int TAIL, int MAXEV) {
    std::vector<Bubble> out;
    long delta = (long)qB - (long)pA;                 // B index = A index + delta
    long i = std::max<long>(0, -delta);
    long endA = std::min<long>((long)A.size(), (long)B.size() - delta);
    long run = 0;                                     // matched bases since last event
    while (i < endA && (int)out.size() < MAXEV) {
        long j = i + delta;
        if (j < 0 || j >= (long)B.size()) break;
        if (A[(size_t)i] == B[(size_t)j]) { ++i; ++run; continue; }
        if (run < 12) { ++i; run = 0; continue; }     // need a left context
        int found = 0; std::string ins;
        for (int g = 1; g <= maxindel && !found; ++g) {
            // A carries g extra bases -> deletion in B
            if (i + g + TAIL <= (long)A.size() && j + TAIL <= (long)B.size()) {
                int mm = 0;
                for (int t = 0; t < TAIL; ++t)
                    if (A[(size_t)(i + g + t)] != B[(size_t)(j + t)]) { ++mm; if (mm > 1) break; }
                if (mm <= 1) { found = 1; }
            }
            // B carries g extra bases -> insertion in B
            if (!found && j + g + TAIL <= (long)B.size() && i + TAIL <= (long)A.size()) {
                int mm = 0;
                for (int t = 0; t < TAIL; ++t)
                    if (A[(size_t)(i + t)] != B[(size_t)(j + g + t)]) { ++mm; if (mm > 1) break; }
                if (mm <= 1) { found = 2; ins = B.substr((size_t)j, (size_t)g); }
            }
            if (found) {
                Bubble b; b.ok = true; b.apos = (uint32_t)i; b.len = g;
                if (found == 1) { b.type = 0; i += g; delta -= g; }
                else            { b.type = 1; b.ins = ins; delta += g; }
                out.push_back(b);
                run = 0;
            }
        }
        if (!found) { ++i; run = 0; }                 // a substitution: keep going
    }
    return out;
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
    //
    // PARALLEL (2026-09-03). This loop was the single largest phase in the
    // caller -- 738 s of the first full-chr20 run, and it runs TWICE (here for
    // the pileup substrate, again at the bubble substrate inside indel_pass).
    // It was serial only because it had never been profiled at scale: the
    // enclosing phase marker is called "ridx_build", but ridx is disabled by
    // default, so all of that time was actually being spent right here under
    // a name that pointed at dead code.
    //
    // Safe to parallelise by inspection, not by hope: iteration `o` writes
    // ONLY S.read_cid[o] / read_pos[o] / read_rc[o] / read_clip[o], which are
    // disjoint across iterations and pre-sized above; `idx` and `S.contigs`
    // are read-only once step 2 has finished (bound through a const ref below
    // so a stray insert cannot compile); every best_* accumulator is declared
    // inside the body and is therefore private. `placed`/`improved` are the
    // only shared scalars and are a reduction. Output is order-independent,
    // so the result is byte-identical to the serial version, not merely
    // equivalent -- which is the gate this change is held to.
    //
    // NOTE the prior negative result this does NOT contradict: OpenMP over
    // CONTIGS in the pileup measured slower (+11%, +13%, comment at the top of
    // run_variant_call). That loop has few, wildly unequal items. This one has
    // 12.6M uniform items; dynamic scheduling covers the tail.
    const auto& cidx = idx;                     // read-only view for the parallel region
    size_t placed = 0, improved = 0;
    #pragma omp parallel for schedule(dynamic, 256) reduction(+:placed,improved)
    for (long long o_ = 0; o_ < (long long)n; ++o_) {
        const size_t o = (size_t)o_;
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
        // Strand 0 used to COPY the read (`std::string r = strand ? ... : raw`)
        // for no reason at all -- two heap allocations per read, ~50M across
        // the two calls at full-chromosome scale. Bind a reference for the
        // forward strand and materialise only the reverse complement.
        std::string rcbuf;
        for (int strand = 0; strand < 2; ++strand) {
            if (strand) rcbuf = rc_str(raw);
            const std::string& r = strand ? rcbuf : raw;
            const int rl = (int)r.size();
            const int step = std::max(1, rl / 8);
            for (int off = 0; off + K <= rl; off += step) {
                uint64_t v; if (!pack25(r.data() + off, v)) continue;
                uint64_t rcv = rc25(v), can = v < rcv ? v : rcv;
                auto it = cidx.find(can);
                if (it == cidx.end()) continue;
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
    // TEMPORARY phase timing (2026-09-03), added to find where wall time
    // actually goes before optimizing blind -- OpenMP over contigs measured
    // SLOWER at both 400kb and 5Mb scale (+11%, +13%), so the assumption
    // that the per-contig loop dominates needs checking, not more tuning.
    using clk = std::chrono::steady_clock;
    auto t_start = clk::now();
    auto elapsed_s = [](clk::time_point a, clk::time_point b){
        return std::chrono::duration<double>(b - a).count(); };
    // Real per-phase RSS, not another structure-by-structure guess. Reads
    // VmRSS straight from /proc/self/status (Linux-only, this project's own
    // server) -- no third-party dependency, no sampling thread, negligible
    // cost next to the phases themselves.
    auto rss_kb = []() -> long {
        FILE* f = fopen("/proc/self/status", "r");
        if (!f) return -1;
        char line[256]; long kb = -1;
        while (fgets(line, sizeof line, f)) {
            if (strncmp(line, "VmRSS:", 6) == 0) { sscanf(line + 6, "%ld", &kb); break; }
        }
        fclose(f);
        return kb;
    };
    auto phase = [&](const char* name, clk::time_point& mark){
        auto now = clk::now();
        fprintf(stderr, "[CAPS-CALL-TIMING] %-16s %8.3fs  RSS=%ldMB\n",
                name, elapsed_s(mark, now), rss_kb() / 1024);
        mark = now;
    };
    auto t_mark = t_start;
    fprintf(stderr, "[CAPS-CALL-TIMING] %-16s %8.3fs  RSS=%ldMB\n", "entry", 0.0, rss_kb() / 1024);

    // Rebuild the calling substrate (collapse + mismatch-tolerant placement).
    //
    // METHOD B DOES NOT NEED THIS. Bubble finding reads only `kc`; the
    // substrate is consulted in this path for exactly one thing -- the ploidy
    // gate samples cd.contigs -- and the ENCODER'S OWN contigs (cd_in) serve
    // that purpose identically, at zero cost, because they already exist.
    // build_substrate re-places all 12.6M reads and measured 738 s serial at
    // full chr20; skipping it is the difference between a ~57 s and a ~150 s
    // Method B run. The pileup path still builds it as before.
    const bool DBG_ONLY_MODE = std::getenv("CAPS_DBG_ONLY") != nullptr
                            && std::getenv("CAPS_DBG") != nullptr;
    CallData cd;
    if (DBG_ONLY_MODE) {
        cd.contigs = cd_in.contigs;          // consensus source for the ploidy gate
        cd.valid = true;
        fprintf(stderr, "[DBG-ONLY] substrate skipped (%zu encoder contigs reused)\n",
                cd.contigs.size());
    } else {
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


    // ── READ SEED INDEX -- DISABLED BY DEFAULT (2026-09-03), see below ──────
    // EBWT2SNP's precision mechanism, from its paper: every emitted fragment
    // of length 2k+1 must be an actual SUBSTRING of at least C real READS
    // (within Hamming distance 2). It reports 99.13% precision against
    // DiscoSNP++'s 77.80% -- BUT on SIMULATED chr22 data (29x). CORRECTED
    // 2026-09-04 after reading the primary source (docs/EBWT_LITERATURE_
    // CORRECTION.md): on the paper's own REAL chr1 data (43-47x) eBWT2SNP is
    // 66.62% precision and DiscoSNP++ is actually MORE precise at 74.57%. The
    // guarantee below was real and worth testing, but its justification here
    // was a simulated-data number, not a real one, and the real-data result
    // for the source method is BELOW what this project already measures.
    // Our previous test only asked whether individual 31-mers appear in the
    // read k-mer TABLE -- far weaker, because a CHIMERIC junction (two
    // paralogs joined at a shared anchor) can have every one of its 31-mers
    // present, each contributed by a different read, while no single read
    // contains the whole fragment. That is precisely our false-positive class.
    //
    // THE FLIP: this index was already recorded as near-vacuous for us
    // (docs/INDEL_PRECISION_ROOT_CAUSE.md, "Group B -- read-support tests
    // are near-vacuous for us": every candidate we emit is read-supported BY
    // CONSTRUCTION, since it was built FROM reads in the first place -- the
    // test EBWT2SNP needs to reject chimeras built from a GRAPH has nothing
    // to reject here). That was documented, but the index was still built
    // every run because a config flag (CAPS_NO_READSUB) existed to skip it,
    // not because the default did. Verified fresh, not just trusted from the
    // doc: ran all 5 standard windows (r2,r3,r4,r5,na) with CAPS_NO_READSUB=1
    // and diffed every TP/FP/FN against the enabled run -- BYTE-IDENTICAL on
    // all 5, SNV and indel both. So this is not a shrink like ridx's sorted-
    // array rewrite below (kept for the record, now dead code under the
    // opt-in flag) -- it is a full elimination of the largest remaining
    // single structure (~900 MB at 971K reads, projected multi-GB at full
    // chr20), for a change already measured to move nothing.
    // CAPS_FORCE_READSUB=1 restores the old default for re-measurement if
    // this ever needs revisiting (e.g. a future substrate change makes
    // candidates less strictly read-derived than they are today).
    struct RidxEntry { uint64_t kmer; uint32_t read; uint32_t pos; };
    std::vector<RidxEntry> ridx;
    const bool WANT_RSUB = std::getenv("CAPS_FORCE_READSUB") != nullptr;
    if (WANT_RSUB) {
        ridx.reserve(seqs.size() * 26);   // ~(148-25)/5+1 25-mers/read, before capping
        for (uint32_t i = 0; i < (uint32_t)seqs.size(); ++i) {
            const std::string& q = seqs[i];
            for (size_t j = 0; j + 25 <= q.size(); j += 5) {
                uint64_t v; if (!detail::pack25(q.data() + j, v)) continue;
                uint64_t rv = detail::rc25(v), cn = v < rv ? v : rv;
                ridx.push_back({cn, i, (uint32_t)j});
            }
        }
        std::stable_sort(ridx.begin(), ridx.end(),
            [](const RidxEntry& a, const RidxEntry& b){ return a.kmer < b.kmer; });
        // Cap each key's group at 40 entries IN PLACE. A second `capped`
        // vector here would hold two full copies at once -- at full-chromosome
        // scale that is two multi-GB buffers alive simultaneously, i.e. the
        // very overhead this change exists to remove, reintroduced at the
        // last step. Compaction only ever keeps an order-preserving subset,
        // so one write cursor trailing the read cursor is enough: w <= j at
        // every step, so no element is overwritten before it has been read.
        size_t w = 0;
        for (size_t i = 0; i < ridx.size(); ) {
            size_t j = i; int cnt = 0;
            while (j < ridx.size() && ridx[j].kmer == ridx[i].kmer) {
                if (cnt < 40) ridx[w++] = ridx[j];
                ++cnt; ++j;
            }
            i = j;
        }
        ridx.resize(w);
        // Only reclaim the tail when capping actually removed enough to be
        // worth a reallocation -- shrink_to_fit copies into a fresh buffer,
        // so calling it unconditionally would itself spike to ~2x right here
        // for, in the common case (most k-mers occur well under 40 times),
        // almost no saving.
        if (w < ridx.capacity() / 5 * 4) ridx.shrink_to_fit();
    }
    auto ridx_lo = [&](uint64_t cn){
        return std::lower_bound(ridx.begin(), ridx.end(), cn,
            [](const RidxEntry& e, uint64_t k){ return e.kmer < k; }); };
    auto ridx_hi = [&](uint64_t cn){
        return std::upper_bound(ridx.begin(), ridx.end(), cn,
            [](uint64_t k, const RidxEntry& e){ return k < e.kmer; }); };
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
            auto lo = ridx_lo(cn), hi = ridx_hi(cn);
            for (auto e = lo; e != hi; ++e) {
                const std::string& q = seqs[e->read];
                // try the read forward and reverse-complemented
                for (int st = 0; st < 2; ++st) {
                    std::string qq = st ? detail::rc_str(q) : q;
                    // the seed sits at e->pos (fwd) -- for rc, recompute by search
                    for (int shift = -2; shift <= 2; ++shift) {
                        long start = (st ? -1 : (long)e->pos - off + shift);
                        if (st) break;                       // rc handled by frag rc below
                        if (start < 0 || start + L > (long)qq.size()) continue;
                        int mm = 0;
                        for (int t = 0; t < L && mm <= tol; ++t)
                            if (qq[(size_t)(start + t)] != frag[(size_t)t]) ++mm;
                        if (mm <= tol) { hits.insert(e->read); break; }
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
            auto lo = ridx_lo(cn), hi = ridx_hi(cn);
            for (auto e = lo; e != hi; ++e) {
                const std::string& q = seqs[e->read];
                long start = (long)e->pos - off;
                if (start < 0 || start + L > (long)q.size()) continue;
                int mm = 0;
                for (int t = 0; t < L && mm <= tol; ++t)
                    if (q[(size_t)(start + t)] != rf[(size_t)t]) ++mm;
                if (mm <= tol) hits.insert(e->read);
                if ((int)hits.size() >= need) return true;
            }
        }
        return (int)hits.size() >= need;
    };

    phase("ridx_build", t_mark);
    // ── 1. Internal canonical-31-mer counts ──
    // CHUNKED COUNTING (2026-09-03) -- the properly-done version of the
    // sorted-array idea that was tried and REVERTED once already (see git
    // history / prior comment here): a single-shot "collect everything, then
    // sort" spikes to ~11.9 GB transient at full chr20 scale, which can be
    // WORSE at peak than the hash map it replaces. This version bounds that
    // transient to one batch: collect and RLE-compress BATCH_KMERS raw
    // k-mers at a time into a small sorted run, discard the raw batch, then
    // k-way merge the runs (summing counts for a key that appears in more
    // than one run). Peak extra memory is O(batch size + distinct k-mers),
    // never O(total k-mer occurrences) -- the number that made the earlier
    // attempt fail. Output (the (kmer,count) set) is identical either way.
    // Batch size is DIVIDED BY THREAD COUNT so aggregate batch memory is the
    // same as the serial version. Each thread reserving the full 32M entries
    // (268 MB) would cost 268 MB x nthreads -- 3.2 GB on 12 cores -- turning
    // the parallel speedup into a RAM regression. Smaller batches simply mean
    // more, shorter runs, which the k-way merge handles natively.
    size_t BATCH_KMERS = 32u << 20;       // ~256 MB of raw uint64 in total
    {
        int nth = 1;
        #ifdef _OPENMP
        nth = omp_get_max_threads();
        #endif
        if (nth > 1) BATCH_KMERS = std::max<size_t>(1u << 20, BATCH_KMERS / (size_t)nth);
    }
    struct KC { uint64_t kmer; uint32_t cnt; };
    // Reads arrive 2-bit packed in graph-only mode (see 106_inprocess.cpp
    // SEQ_PACK). Unpack on demand into a caller-supplied buffer so only one
    // read is ever expanded at a time.
    const bool SEQ_PACKED = std::getenv("CAPS_DBG_ONLY") && std::getenv("CAPS_DBG");
    auto unpack_read = [&](const std::string& src, std::string& dst) -> const std::string& {
        if (!SEQ_PACKED) return src;
        if (src.empty()) { dst.clear(); return dst; }
        if ((uint8_t)src[0] == 1) {           // raw: read contained non-ACGT
            dst.assign(src, 1, std::string::npos);
            return dst;
        }
        if (src.size() < 3) { dst.clear(); return dst; }
        const size_t L = (size_t)(uint8_t)src[1] | ((size_t)(uint8_t)src[2] << 8);
        dst.resize(L);
        for (size_t i = 0; i < L; ++i)
            dst[i] = "ACGT"[((uint8_t)src[3 + (i >> 2)] >> (2 * (i & 3))) & 3];
        return dst;
    };
    std::vector<std::vector<KC>> kc_runs;
    // Disk-spill for the counting runs. OFF unless CAPS_KC_SPILLDIR is set, so
    // default behaviour is byte-identical. See the spill site below for why
    // this is the peak-RAM lever and how GATB uses the same trade.
    const char* SPILLDIR = std::getenv("CAPS_KC_SPILLDIR");
    const bool  SPILL    = (SPILLDIR != nullptr);
    // 2^SPILL_BITS key-range partitions. 256 partitions keeps each one small
    // enough to sort in RAM while staying far under any open-file limit.
    const int   SPILL_BITS = std::getenv("CAPS_KC_SPILLBITS")
                           ? atoi(std::getenv("CAPS_KC_SPILLBITS")) : 8;
    const size_t PBUF_MAX  = 1u << 16;      // entries buffered per partition per thread
    // Superkmer spill: minimizer-partitioned, 2-bit packed. Verified at
    // 1.161 bytes/k-mer vs 12 for records, k-mer multiset identical.
    const bool  SUPERK = SPILL && std::getenv("CAPS_KC_SUPERKMER") != nullptr;
    // (definition for RunSrc::RBUF_ lives at namespace scope below)
    {
        // PARALLEL. Each thread collects and RLE-compresses its own batches
        // into its own runs; the k-way merge below already sums counts for a
        // key across every run, so partitioning the reads across threads
        // changes only WHICH run holds a key, never the summed count. The
        // (kmer,count) set is therefore identical to the serial version by
        // construction -- this is not an approximation.
        #pragma omp parallel
        {
            std::vector<uint64_t> batch; batch.reserve(BATCH_KMERS);
            std::vector<std::vector<KC>> myruns;
            std::vector<std::vector<uint64_t>> pbuf(SPILL ? (1u << SPILL_BITS) : 0);
            auto flush = [&](){
                if (batch.empty()) return;
                std::sort(batch.begin(), batch.end());
                std::vector<KC> run; run.reserve(batch.size());
                for (size_t i = 0; i < batch.size(); ) {
                    size_t j = i;
                    while (j < batch.size() && batch[j] == batch[i]) ++j;
                    run.push_back({batch[i], (uint32_t)(j - i)});
                    i = j;
                }
                // SPILL BY KEY RANGE (GATB's design, simplified).
                //
                // The first version of this spill wrote each batch as one
                // sorted RUN. That was the wrong shape: every run spans the
                // whole key space, so the k-way merge has to hold all runs
                // open simultaneously and peak stays proportional to the total.
                //
                // GATB does not do that. SortingCountAlgorithm.cpp:804 routes
                // each superkmer to a partition chosen by its MINIMIZER --
                //     size_t p = _repartition(superKmer.minimizer);
                // so every partition owns a DISJOINT slice of key space and
                // partitions are processed one at a time, with no merge across
                // them at all. Their peak is one partition; ours was the sum.
                //
                // We can do it more simply than they do: partition on the TOP
                // BITS of the canonical k-mer. Then partition order IS key
                // order, so concatenating the sorted partitions yields a
                // globally sorted kc with no merge step whatsoever. (GATB needs
                // minimizers because superkmers must keep adjacent k-mers
                // together; we store plain k-mers, so we do not.)
                if (SPILL) {
                    // SPILL ONLY THE KEY, NOT THE COUNT.
                    // A (kmer,count) record is 12 bytes; the count is
                    // recomputed when the partition is read back and sorted,
                    // so writing it is pure waste. Measured: the count made
                    // the full-chr20 spill 25 GB, which exhausted a 233 GB
                    // disk. Dropping it is a 33% cut for free. (GATB gets far
                    // more -- ~0.95 B/k-mer -- by packing SUPERKMERS, which
                    // needs the whole minimizer chain; this is the part
                    // available without it.)
                    for (const KC& e : run) {
                        const uint32_t part = (uint32_t)(e.kmer >> (62 - SPILL_BITS));
                        std::vector<uint64_t>& pb = pbuf[part];
                        // one entry per OCCURRENCE, so the count survives as
                        // repetition and is recovered by the read-back sort
                        for (uint32_t rep = 0; rep < e.cnt; ++rep) pb.push_back(e.kmer);
                        if (pb.size() >= PBUF_MAX) {
                            char path[512];
                            snprintf(path, sizeof path, "%s/p%05u_t%d.bin",
                                     SPILLDIR, part, omp_get_thread_num());
                            FILE* f = fopen(path, "ab");
                            if (f) { fwrite(pb.data(), sizeof(uint64_t), pb.size(), f); fclose(f); }
                            pb.clear();
                        }
                    }
                    batch.clear();
                    return;
                }
                myruns.push_back(std::move(run));
                batch.clear();
            };
            // ROLLING K-MER WINDOW.
            //
            // The previous loop called pack31() and canon31() at every
            // position, and BOTH are O(k): pack31 loops 31 times to build the
            // k-mer, rc31 loops 31 times to reverse-complement it. At full
            // chr20 that is ~1.4 BILLION k-mer positions x ~62 operations
            // = ~87 billion operations, and it is the reason this phase costs
            // 72.8 s -- almost none of that work is necessary.
            //
            // A k-mer and its reverse complement can both be updated in O(1)
            // from the previous position:
            //     fwd = ((fwd << 2) | b) & MASK31        -- shift in the new base
            //     rev = (rev >> 2) | ((3-b) << 60)       -- complement, shift out
            // so the inner loop drops from ~62 operations to ~4. `have`
            // tracks how many valid bases are currently in the window, which
            // is what makes N bases a reset rather than a special case.
            const uint64_t KMASK = (~0ULL) >> 2;          // 62 bits = 31 bases
            // ── SUPERKMER PATH ──────────────────────────────────────────────
            // Consecutive k-mers in a read almost always share a minimizer, so
            // a RUN of them can be stored as ONE overlapping sequence: L k-mers
            // occupy L+k-1 bases, packed 4 bases/byte. Verified standalone
            // before wiring: 1.161 bytes/k-mer against 12 for (kmer,count)
            // records, with the k-mer multiset proven identical.
            // This is what makes the spill sustainable -- the record format hit
            // 25 GB at full chr20 and exhausted the disk twice.
            if (SUPERK) {
                const int MM = 10;                       // minimizer size
                const uint32_t MMASK = (1u << (2*MM)) - 1;
                std::vector<std::vector<uint8_t>> skbuf(1u << SPILL_BITS);
                auto sk_flush = [&](uint32_t part){
                    if (skbuf[part].empty()) return;
                    char path[512];
                    snprintf(path, sizeof path, "%s/sk%05u_t%d.bin",
                             SPILLDIR, part, omp_get_thread_num());
                    FILE* f = fopen(path, "ab");
                    if (f) { fwrite(skbuf[part].data(), 1, skbuf[part].size(), f); fclose(f); }
                    skbuf[part].clear();
                };
                std::string ubuf;
                #pragma omp for schedule(static)
                for (long long si = 0; si < (long long)seqs.size(); ++si) {
                    const std::string& s = unpack_read(seqs[(size_t)si], ubuf);
                    if (s.size() < 31) continue;
                    // rolling minimizer per k-mer (monotonic deque, O(1) amortised)
                    const size_t nk = s.size() - 31 + 1;
                    std::vector<uint32_t> mini(nk, 0xFFFFFFFFu);
                    {
                        std::deque<std::pair<uint32_t,size_t>> dq;
                        uint32_t mv = 0; int have = 0;
                        // lastN is the most recent invalid base. A k-mer is only
                        // valid if its whole 31-base window starts after it --
                        // without this the deque refills 10 bases past an N and
                        // assigns minimizers to k-mers that still contain it,
                        // which silently dropped 93 k-mers on the r2 window
                        // (1,063,514 vs 1,063,607) because those superkmers were
                        // then rejected wholesale.
                        long long lastN = -1;
                        for (size_t i = 0; i < s.size(); ++i) {
                            const int b = b2i(s[i]);
                            if (b < 0) { dq.clear(); have = 0; lastN = (long long)i; continue; }
                            mv = ((mv << 2) | (uint32_t)b) & MMASK;
                            if (++have < MM) continue;
                            const size_t mpos = i - MM + 1;
                            while (!dq.empty() && dq.back().first >= mv) dq.pop_back();
                            dq.push_back({mv, mpos});
                            if (i + 1 >= 31) {
                                const size_t kpos = i + 1 - 31;
                                if ((long long)kpos <= lastN) continue;   // window still spans an N
                                while (!dq.empty() && dq.front().second < kpos) dq.pop_front();
                                if (!dq.empty()) mini[kpos] = dq.front().first;
                            }
                        }
                    }
                    size_t st = 0;
                    while (st < nk) {
                        if (mini[st] == 0xFFFFFFFFu) { ++st; continue; }   // window had an N
                        size_t en = st;
                        while (en + 1 < nk && mini[en+1] == mini[st] && (en - st + 1) < 255) ++en;
                        const size_t L = en - st + 1, nb = L + 31 - 1;
                        // reject if any base in the span is invalid
                        bool ok = true;
                        for (size_t i = 0; i < nb && ok; ++i) if (b2i(s[st+i]) < 0) ok = false;
                        if (ok) {
                            const uint32_t part = mini[st] & ((1u << SPILL_BITS) - 1);
                            std::vector<uint8_t>& ob = skbuf[part];
                            ob.push_back((uint8_t)L);
                            const size_t need = (nb + 3) / 4, base = ob.size();
                            ob.resize(base + need, 0);
                            for (size_t i = 0; i < nb; ++i)
                                ob[base + (i >> 2)] |= (uint8_t)(b2i(s[st+i]) << (2 * (i & 3)));
                            if (ob.size() >= (1u << 20)) sk_flush(part);
                        }
                        st = en + 1;
                    }
                }
                for (uint32_t pp = 0; pp < (1u << SPILL_BITS); ++pp) sk_flush(pp);
            } else {
            std::string ubuf2;
            #pragma omp for schedule(static)
            for (long long si = 0; si < (long long)seqs.size(); ++si) {
                const std::string& s = unpack_read(seqs[(size_t)si], ubuf2);
                if (s.size() < 31) continue;
                uint64_t fwd = 0, rev = 0;
                int have = 0;
                for (size_t i = 0; i < s.size(); ++i) {
                    const int b = b2i(s[i]);
                    if (b < 0) { have = 0; fwd = rev = 0; continue; }   // N: restart window
                    fwd = ((fwd << 2) | (uint64_t)b) & KMASK;
                    rev = (rev >> 2) | ((uint64_t)(3 - b) << 60);
                    if (++have < 31) continue;
                    if (have > 31) have = 31;
                    batch.push_back(fwd < rev ? fwd : rev);
                    if (batch.size() >= BATCH_KMERS) flush();
                }
            }
            }
            flush();
            if (SPILL) {
                for (uint32_t part = 0; part < pbuf.size(); ++part) {
                    if (pbuf[part].empty()) continue;
                    char path[512];
                    snprintf(path, sizeof path, "%s/p%05u_t%d.bin",
                             SPILLDIR, part, omp_get_thread_num());
                    FILE* f = fopen(path, "ab");
                    if (f) { fwrite(pbuf[part].data(), sizeof(uint64_t), pbuf[part].size(), f); fclose(f); }
                    std::vector<uint64_t>().swap(pbuf[part]);
                }
            }
            #pragma omp critical(kcruns)
            for (auto& r : myruns) kc_runs.push_back(std::move(r));
        }
    }
    // k-way merge: at each step pick the smallest head key across all runs
    // and sum every run's count for that key. Runs are few (total k-mer
    // occurrences / BATCH_KMERS -- roughly 47 at full chr20 scale), so a
    // linear scan per output element is simple and correct; a heap would be
    // faster but is not needed at this run count.
    // HEAP MERGE, TWO PASSES (2026-09-03). The comment above ("a heap would be
    // faster but is not needed at this run count") was correct about the run
    // count and wrong about the cost, because the cost is driven by the OUTPUT
    // length, not the run count. The previous version did two linear scans over
    // all runs per output element: at full chr20 that is roughly
    //     ~400M distinct 31-mers x 47 runs x 2 scans = ~37 BILLION comparisons,
    // which is essentially all of this phase's 166.7 s. A min-heap over run
    // heads makes it O(output x log runs) instead.
    //
    // Pass 1 counts distinct keys so pass 2 can allocate `kc` exactly once.
    // The tempting one-pass alternative -- reserve the summed run lengths --
    // is an exact upper bound but a useless one: every run carries its own copy
    // of each common key, so that sum is ~47 x 32M entries (~18 GB reserved) to
    // hold ~400M (~4.8 GB) of real output. Two heap passes cost ~2x the merge
    // but the merge is now ~8-17x cheaper, so this is still several times
    // faster than the code it replaces AND removes the ~14.4 GB transient the
    // un-reserved doubling caused at the final reallocation.
    //
    // Byte-identical by construction: a k-way merge is a pure function of its
    // runs, and both passes visit keys in the same ascending order and sum the
    // same per-run counts.
    std::vector<KC> kc;
    if (SUPERK) {
        // ── READ BACK SUPERKMERS ────────────────────────────────────────────
        // One partition at a time: read its superkmers, expand each into its L
        // k-mers, canonicalise, sort, and run-length count. Partitions are keyed
        // by MINIMIZER, so they are not key-ranges -- a k-mer from any partition
        // can sort anywhere -- and the per-partition results are therefore
        // written to sorted temp files and k-way merged at the end. Peak stays
        // at one partition plus the output, never the sum.
        const uint64_t KMASK = (~0ULL) >> 2;
        std::vector<std::string> sorted_parts;
        size_t total_k = 0;
        // PARALLEL OVER PARTITIONS. Each partition reads its own files and
        // writes its own sorted output, so they are fully independent, and this
        // is the largest parallelism available in the read-back phase.
        //
        // This was once reverted after the parallel version produced 2,273
        // extra k-mers -- but that was a MISDIAGNOSIS. The extra k-mers came
        // from 2-bit read packing mapping non-ACGT bases to 'A', which made
        // k-mers spanning an N look valid. The parallelism was never at fault;
        // reverting it removed a correct optimisation. Restored once the N bug
        // was fixed and kc identity (1,063,607) re-verified.
        sorted_parts.resize(1u << SPILL_BITS);
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:total_k)
        for (long long pi_ = 0; pi_ < (long long)(1u << SPILL_BITS); ++pi_) {
            const uint32_t pi = (uint32_t)pi_;
            std::vector<uint8_t> raw;
            std::vector<uint64_t> kms;
            for (int t = 0; t < 64; ++t) {
                char path[512];
                snprintf(path, sizeof path, "%s/sk%05u_t%d.bin", SPILLDIR, pi, t);
                FILE* f = fopen(path, "rb");
                if (!f) continue;
                fseek(f, 0, SEEK_END); const long sz = ftell(f); fseek(f, 0, SEEK_SET);
                const size_t base = raw.size();
                raw.resize(base + (size_t)sz);
                if (fread(raw.data() + base, 1, (size_t)sz, f) != (size_t)sz) raw.resize(base);
                fclose(f); ::remove(path);
            }
            if (raw.empty()) continue;
            // expand: [uint8 L][ceil((L+30)/4) packed bases]
            for (size_t off = 0; off + 1 <= raw.size(); ) {
                const size_t L = raw[off]; ++off;
                if (L == 0) break;
                const size_t nb = L + 31 - 1, need = (nb + 3) / 4;
                if (off + need > raw.size()) break;
                const uint8_t* pk = raw.data() + off;
                uint64_t fwd = 0, rev = 0;
                for (size_t i = 0; i < nb; ++i) {
                    const int b = (pk[i >> 2] >> (2 * (i & 3))) & 3;
                    fwd = ((fwd << 2) | (uint64_t)b) & KMASK;
                    rev = (rev >> 2) | ((uint64_t)(3 - b) << 60);
                    if (i + 1 >= 31) kms.push_back(fwd < rev ? fwd : rev);
                }
                off += need;
            }
            if (kms.empty()) continue;
            // ── RADIX PRE-BINNING before sorting ────────────────────────────
            // GATB never sorts one large array: PartitionsCommand.cpp bins
            // k-mers into 256 radix buckets as they are read and sorts each
            // bucket separately (executeRead -> executeSort -> executeDump).
            // One n log n over a huge array thrashes cache; O(n) binning plus
            // many cache-resident sorts does not. Bins are ordered by the top
            // 8 bits, so concatenating the sorted bins is globally sorted --
            // the result is identical to sorting the whole array, which is why
            // this needs no separate correctness gate beyond the kc identity
            // check that already covers it.
            {
                const int RB = 8, NB = 1 << RB;
                std::vector<uint32_t> cnt(NB + 1, 0);
                for (uint64_t v : kms) ++cnt[(size_t)(v >> (62 - RB)) + 1];
                for (int i = 0; i < NB; ++i) cnt[i+1] += cnt[i];
                std::vector<uint64_t> tmp(kms.size());
                std::vector<uint32_t> pos(cnt.begin(), cnt.end() - 1);
                for (uint64_t v : kms) tmp[pos[(size_t)(v >> (62 - RB))]++] = v;
                // no nested parallel region: the partition loop is already
                // parallel, so these bins are sorted serially within a thread
                for (int b = 0; b < NB; ++b)
                    std::sort(tmp.begin() + cnt[b], tmp.begin() + cnt[b+1]);
                kms.swap(tmp);
            }
            char spath[512];
            snprintf(spath, sizeof spath, "%s/srt%05u.bin", SPILLDIR, pi);
            FILE* sf = fopen(spath, "wb");
            if (sf) {
                for (size_t i = 0; i < kms.size(); ) {
                    size_t j = i; while (j < kms.size() && kms[j] == kms[i]) ++j;
                    const KC e{ kms[i], (uint32_t)(j - i) };
                    fwrite(&e, sizeof(KC), 1, sf);
                    ++total_k; i = j;
                }
                fclose(sf);
                sorted_parts[pi] = spath;
            }
        }
        {   // partitions that produced nothing leave empty slots
            std::vector<std::string> nz;
            for (auto& sp : sorted_parts) if (!sp.empty()) nz.push_back(sp);
            sorted_parts.swap(nz);
        }
        // k-way merge the sorted per-partition files
        struct SF { FILE* f; KC cur; bool ok; };
        std::vector<SF> sf(sorted_parts.size());
        using HE = std::pair<uint64_t, uint32_t>;
        std::priority_queue<HE, std::vector<HE>, std::greater<HE>> pq;
        for (size_t i = 0; i < sorted_parts.size(); ++i) {
            sf[i].f = fopen(sorted_parts[i].c_str(), "rb");
            sf[i].ok = sf[i].f && fread(&sf[i].cur, sizeof(KC), 1, sf[i].f) == 1;
            if (sf[i].ok) pq.push({sf[i].cur.kmer, (uint32_t)i});
        }
        kc.reserve(total_k);
        while (!pq.empty()) {
            const uint64_t key = pq.top().first;
            uint32_t sum = 0;
            while (!pq.empty() && pq.top().first == key) {
                const uint32_t i = pq.top().second; pq.pop();
                sum += sf[i].cur.cnt;
                sf[i].ok = fread(&sf[i].cur, sizeof(KC), 1, sf[i].f) == 1;
                if (sf[i].ok) pq.push({sf[i].cur.kmer, i});
            }
            kc.push_back({key, sum});
        }
        for (size_t i = 0; i < sf.size(); ++i) if (sf[i].f) fclose(sf[i].f);
        for (auto& sp : sorted_parts) ::remove(sp.c_str());
        fprintf(stderr, "[KC-SUPERK] partitions=%zu distinct=%zu\n", sorted_parts.size(), kc.size());
    } else if (SPILL) {
        // PARTITION-AT-A-TIME BUILD. Each partition owns a disjoint, ascending
        // slice of key space (top SPILL_BITS of the canonical k-mer), so:
        //   * partitions can be processed one at a time -- peak is ONE
        //     partition plus the growing output, never the sum of all runs;
        //   * partition order IS key order, so simply appending each sorted
        //     partition produces a globally sorted kc with NO merge step.
        // This is GATB's structure (disjoint partitions, processed
        // independently) without needing minimizers, which they require only
        // because superkmers must keep adjacent k-mers together.
        size_t total_in = 0;
        std::vector<KC> part;
        for (uint32_t pi = 0; pi < (1u << SPILL_BITS); ++pi) {
            part.clear();
            for (int t = 0; t < 64; ++t) {
                char path[512];
                snprintf(path, sizeof path, "%s/p%05u_t%d.bin", SPILLDIR, pi, t);
                FILE* f = fopen(path, "rb");
                if (!f) continue;
                fseek(f, 0, SEEK_END); const long sz = ftell(f); fseek(f, 0, SEEK_SET);
                const size_t n = (size_t)sz / sizeof(KC);
                const size_t base = part.size();
                part.resize(base + n);
                if (fread(part.data() + base, sizeof(KC), n, f) != n) part.resize(base);
                fclose(f);
                ::remove(path);
            }
            if (part.empty()) continue;
            total_in += part.size();
            // RESERVE FROM THE FIRST PARTITION'S DISTINCT RATE.
            // kc.push_back was growing by doubling, and at full chr20 the
            // final array is ~1.7 GB, so the last reallocation held old+new
            // ~3.4 GB -- a pure transient with no phase attached to it
            // (measured at 2M reads: RSS flat at 2772 MB, PEAK 4127 MB).
            // The two-pass merge this replaced used an exact reserve; the
            // partitioned build lost it. Partitions are uniform slices of key
            // space, so the first one's distinct count scales to a good total
            // estimate. 15% headroom absorbs partition-to-partition variance;
            // being slightly over only wastes a little, while being under
            // merely restores one doubling.
            if (kc.capacity() == 0 && !part.empty()) {
                size_t d0 = 1;
                for (size_t i = 1; i < part.size(); ++i)
                    if (part[i].kmer != part[i-1].kmer) ++d0;
                const size_t est = (size_t)((double)d0 * (1u << SPILL_BITS) * 1.15);
                kc.reserve(est);
                fprintf(stderr, "[KC-SPILL] reserve estimate %zu from partition 0\n", est);
            }
            std::sort(part.begin(), part.end(),
                      [](const KC& x, const KC& y){ return x.kmer < y.kmer; });
            for (size_t i = 0; i < part.size(); ) {
                size_t j = i; uint32_t sum = 0;
                while (j < part.size() && part[j].kmer == part[i].kmer) { sum += part[j].cnt; ++j; }
                kc.push_back({part[i].kmer, sum});
                i = j;
            }
        }
        std::vector<KC>().swap(part);
        fprintf(stderr, "[KC-SPILL] partitions=%d entries_in=%zu distinct=%zu\n",
                (1 << SPILL_BITS), total_in, kc.size());
    } else {
        // In-RAM path (default): heap merge over the runs, exactly as before.
        // Two passes so the output is reserved exactly once.
        using HE = std::pair<uint64_t, uint32_t>;
        auto merge_pass = [&](bool fill_out) -> size_t {
            std::vector<size_t> idx(kc_runs.size(), 0);
            std::priority_queue<HE, std::vector<HE>, std::greater<HE>> pq;
            for (uint32_t r = 0; r < (uint32_t)kc_runs.size(); ++r)
                if (!kc_runs[r].empty()) pq.push({kc_runs[r][0].kmer, r});
            size_t nout = 0;
            while (!pq.empty()) {
                const uint64_t key = pq.top().first;
                uint32_t sum = 0;
                while (!pq.empty() && pq.top().first == key) {
                    const uint32_t r = pq.top().second; pq.pop();
                    sum += kc_runs[r][idx[r]].cnt;
                    if (++idx[r] < kc_runs[r].size()) pq.push({kc_runs[r][idx[r]].kmer, r});
                }
                if (fill_out) kc.push_back({key, sum});
                ++nout;
            }
            return nout;
        };
        kc.reserve(merge_pass(false));
        merge_pass(true);
    }
    kc_runs.clear(); kc_runs.shrink_to_fit();
    auto kc_find = [&](uint64_t key) -> const KC* {
        auto it = std::lower_bound(kc.begin(), kc.end(), key,
            [](const KC& e, uint64_t k){ return e.kmer < k; });
        return (it != kc.end() && it->kmer == key) ? &*it : nullptr;
    };
    uint32_t H = 30;
    // The k-mer count histogram's VALLEY -- the minimum between the error peak
    // and the true single-copy peak -- is the standard, dataset-measured
    // solid-k-mer threshold. It is computed below to locate H and was being
    // discarded; it is hoisted here because the dBG channel needs exactly this
    // quantity and must NOT use a fixed constant for it (see MINC).
    uint32_t KVALLEY = 2;
    // Is this sample actually DIPLOID? A heterozygous diploid's k-mer count
    // histogram is bimodal -- a heterozygous peak at the haploid depth H and a
    // homozygous peak at ~2H -- while a haploid genome shows only the single
    // peak. This is the standard GenomeScope/findGSE signature and it costs
    // nothing here because the histogram is already built to locate H.
    //
    // WHY IT MATTERS: the dBG bubble channel looks for heterozygous bubbles.
    // Run on a HAPLOID sample it can only produce noise -- measured, it emitted
    // 163 "sites" on E. coli and 180 on SARS-CoV-2, organisms that have no
    // heterozygous variants at all. Claim 2 is entirely diploid human so this
    // never showed up in benchmarking, which is exactly why it needs a guard
    // rather than an assumption.
    bool looks_diploid = true;
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
        for (auto& kv : kc) cnt_max = std::max(cnt_max, kv.cnt);
        cnt_max = std::min(cnt_max, 2000u);
        if (cnt_max >= 4) {
            uint32_t bw = 1;
            uint32_t nb = cnt_max / bw + 2;
            std::vector<uint64_t> bkt(nb, 0);
            for (auto& kv : kc)
                if (kv.cnt >= 2 && kv.cnt <= cnt_max) bkt[kv.cnt / bw]++;
            size_t valley = 2;
            for (size_t b = 1; b + 1 < nb / 3 && b + 1 < bkt.size(); ++b) {
                if (bkt[b] < bkt[b - 1] && bkt[b] < bkt[b + 1]) { valley = b; break; }
            }
            size_t peak = (valley + 1 < bkt.size()) ? valley + 1 : valley;
            for (size_t b = valley + 1; b < bkt.size(); ++b)
                if (bkt[b] > bkt[peak]) peak = b;
            uint32_t H_hist = (uint32_t)((peak + 0.5) * bw);
            if (H_hist >= 5) { H = H_hist; KVALLEY = (uint32_t)std::max<size_t>(2, valley * bw); }
            // (A histogram-bimodality ploidy test was tried here and REMOVED:
            // it did not discriminate. Measured ratios: HG002 0.125 (diploid),
            // SARS-CoV-2 0.331 (haploid), E. coli 0.088 (haploid) -- the
            // haploids fell both above AND below the diploid, so no threshold
            // separates them. The reasoning was wrong: at 30x with k=31 the
            // located peak sits near the HOMOZYGOUS depth, so "mass at 2H"
            // measures repeats, not ploidy. Ploidy is detected instead by the
            // pair test below, the mechanism already validated 5/5 by HETSCAN.)

        }
    }
    auto kcount = [&](const std::string& km) -> uint32_t {
        if (km.size() != 31) return 0;
        uint64_t v; if (!pack31(km.data(), v)) return 0;
        const KC* e = kc_find(canon31(v)); return e ? e->cnt : 0u;
    };

    // ── 2-6a. Per-contig pileup, candidates, flank pass, filters, SNV emit ──
    // CONTIG-BATCHED (2026-09-03). Previously col/C/FLmaj/FLmin/majoff/minoff
    // were built ONCE over every read in the whole chromosome and held
    // resident simultaneously -- measured as the dominant RAM sink on a real
    // full chr20 run (tens of GB). Nothing in this stretch crosses a contig
    // boundary: colkey always carries cid, and every lookup into these five
    // structures is keyed by a (cid,pos) that only ever matches the contig
    // being scanned. So the identical computation is re-scoped to run one
    // contig at a time, with all five structures declared fresh inside the
    // loop and freed before the next contig starts -- peak RAM becomes
    // O(largest single contig's pileup), not O(whole chromosome).
    // Output is byte-identical to the old global-pass version: the original
    // did one `std::sort(kept)` over (cid,pos) pairs from every contig before
    // emitting; emitting per-contig in increasing cid order and concatenating
    // produces the exact same total order, since cid is the primary sort key
    // and each contig's own `kept` is still locally sorted by pos.
    struct Rec { uint32_t cid, pos; std::string seq, qual; };
    struct Cand { int M, mn; int cnt[4]; int d; };
    struct OutRec { uint32_t cid, pos; std::string ref, alt, info; int src; };
    std::vector<OutRec> orecs;
    size_t total_candidates = 0;

    // Method B skips build_substrate, so there are no per-read placements to
    // index here. The pileup is the only consumer and it does not run in that
    // mode; guarding is required because this loop sits BEFORE the DBG-only
    // exit and would otherwise dereference an empty array.
    std::vector<std::vector<uint32_t>> reads_by_contig(cd.contigs.size());
    if (cd.read_cid.size() >= n) {
        for (size_t oi = 0; oi < n; ++oi) {
            uint32_t rcid = cd.read_cid[oi];
            if (rcid < cd.contigs.size()) reads_by_contig[rcid].push_back((uint32_t)oi);
        }
    }

    auto most_common = [](const std::unordered_map<std::string,int>& m, int& cntout) -> std::string {
        std::string best; int bc = 0;
        for (auto& kv : m) if (kv.second > bc) { bc = kv.second; best = kv.first; }
        cntout = bc; return best;
    };

    // PLOIDY-SCALED MINOR-ALLELE THRESHOLD.
    // MAF=0.20 is a DIPLOID constant: a heterozygous allele sits at ~0.5, so a
    // 0.20 floor is generous. In a k-ploid sample a single-haplotype allele
    // sits at ~1/k -- 0.33 at k=3, 0.25 at k=4 -- and sampling plus
    // contig-frame partial coverage pushes real alleles below 0.20 (measured:
    // a true site at AF 0.167 was rejected, and recall on synthetic triploid
    // was 0.950-0.975 against DiscoSNP++'s 1.000 across three seeds, entirely
    // on missed sites).
    // Scale by 2/k so the expected allele fraction and the threshold move
    // together. k=2 gives exactly MAF, so the diploid path -- and every
    // het-SNV result -- is byte-identical.
    const double MAF_K = MAF * 2.0 / (double)PLOIDY;

    // Hoisted out of the per-contig loop (2026-09-03): this override does not
    // depend on cid0 at all -- it re-read the same env var on every contig
    // for no reason, and once the loop below runs on multiple threads,
    // writing the shared `H` from inside it would be a data race. Reading it
    // once here, before any thread starts, is both a correctness fix for the
    // parallel version and a no-op for the sequential one (same value either
    // way).
    if (const char* he = std::getenv("CAPS_HAPLOID_COV")) { int v = atoi(he); if (v > 0) H = (uint32_t)v; }

    // WHAT IS AND IS NOT BATCHED, corrected 2026-09-03 after a review pass
    // caught an over-reach in the first version of this refactor.
    // Only `col` (the pileup: every read base at every covered position) and
    // `recs` are large -- measured shapes put col in the tens of GB at full
    // chr20 scale, recs at a few GB. `C` (candidate columns) and the flank
    // maps are SMALL: 585 candidates on a 400 kb window scales to roughly
    // 1e5 genome-wide, a few hundred MB with the flank contexts included.
    // The first version of this refactor scoped C and the flank maps per
    // contig too, which bought nothing and silently broke `med_cand_depth`
    // -- the depth guard's median, which the original computes over ALL
    // candidates chromosome-wide, became a per-contig median. That guard is
    // inert on the default path (it is gated on !collapse_ran) so no
    // recorded result moved, but it is exactly the kind of silent
    // methodology change this project forbids.
    // So: build col/recs per contig and free them there, accumulate C and
    // the flank maps into thread-local buffers, then merge and run the
    // filter/emit stages EXACTLY as the original wrote them, over the
    // complete global C. Keys are (cid,pos) and each contig is handled by
    // exactly one thread, so no two threads ever produce the same key and
    // merging is a pure move with no count reconciliation.
    phase("kc_H_build", t_mark);

    // ── PLOIDY GATE (the encoder's HETSCAN mechanism, validated 5/5) ────────
    // A heterozygous site yields TWO k-mers differing at exactly one position,
    // both at real depth; a haploid genome yields almost none. This is the
    // encoder's own test (het_pair_frac, threshold 0.024: E. coli 0.0172,
    // M. tuberculosis 0.0209, P. aeruginosa 0.0103 -> haploid; HG002 0.0264,
    // HG003 0.0261 -> diploid), reproduced here because HETSCAN runs AFTER the
    // caller and therefore cannot gate it.
    {
        const uint32_t PLMIN = std::max(2u, H / 4u);
        size_t sampled = 0, pairs = 0;
        // Sample the CONTIGS, not kc. kc holds raw read k-mers including
        // sequencing errors; at high depth those errors clear any count
        // threshold and every one of them has a 1-mismatch partner (its own
        // correct version), so the pair fraction saturates. Measured: sampling
        // kc gave SARS-CoV-2 frac=0.8927 -- a haploid virus scored as more
        // heterozygous than HG002. The encoder's HETSCAN does not have this
        // problem because it samples the assembled pseudogenome, i.e. consensus
        // sequence with errors already removed. cd.contigs is the caller's
        // equivalent of that consensus.
        size_t ctotal = 0;
        for (const auto& c : cd.contigs) ctotal += c.size();
        const size_t STRIDE = std::max<size_t>(41, ctotal / 200000);
        for (const auto& c : cd.contigs) {
            if (c.size() < 31) continue;
            for (size_t off = 0; off + 31 <= c.size(); off += STRIDE) {
                uint64_t kv;
                if (!pack31(c.data() + off, kv)) continue;
                const KC* self = kc_find(canon31(kv));
                if (!self || self->cnt < PLMIN) continue;
                ++sampled;
                bool found = false;
                for (int pos = 0; pos < 31 && !found; ++pos) {
                    const uint64_t sh = (uint64_t)(2 * (30 - pos));
                    const uint64_t cur = (kv >> sh) & 3ULL;
                    for (uint64_t b = 0; b < 4 && !found; ++b) {
                        if (b == cur) continue;
                        const uint64_t alt = (kv & ~(3ULL << sh)) | (b << sh);
                        const KC* e = kc_find(canon31(alt));
                        if (e && e->cnt >= PLMIN) found = true;
                    }
                }
                if (found) ++pairs;
            }
        }
        const double frac = sampled ? (double)pairs / (double)sampled : 0.0;
        looks_diploid = frac >= 0.024;
        fprintf(stderr, "[PLOIDY] sampled=%zu pairs=%zu frac=%.4f -> %s\n",
                sampled, pairs, frac, looks_diploid ? "diploid" : "haploid");
    }

    // ── DBG BUBBLE CHANNEL — DiscoSNP++'s algorithm on the graph WE ALREADY BUILD ──
    //
    // THE REALIZATION. A de Bruijn graph is exactly two things: a set of
    // k-mers, and the implicit edge X->Y whenever Y = (X shifted left one base
    // + one new base) is also in the set. `kc` above IS that set -- every
    // canonical 31-mer of every read, with its count -- and `kc_find` IS the
    // membership test. We spend 166.7 s building a full de Bruijn graph and
    // then use it to compute ONE SCALAR (H, the coverage threshold below).
    //
    // So DiscoSNP++'s bubble finder does not need a new index, a new pass over
    // the reads, or a second assembly. Every primitive it uses maps onto a
    // kc_find:
    //
    //   their graph.successors(node)        -> 4 kc_find probes, keep cnt >= minc
    //   their branching node (>= 2 succ)    -> same test
    //   their graph.successors(node1,node2) -> for b in ACGT, BOTH extensions
    //                                          present (their getNodesCouple in
    //                                          GATB Graph.cpp:1602-1640 requires
    //                                          the SAME nt on both paths)
    //   their closure (succ.first==.second) -> node1 == node2
    //
    // This is a faithful port of Bubble.cpp:300-347 (start: branching node,
    // every successor PAIR) and Bubble.cpp:509-527 (expand: lockstep walk,
    // close on node equality), not an approximation of it. Two earlier attempts
    // this session failed precisely because they lacked this closure test:
    // chaining ties are exact-match and so blind to heterozygosity
    // (docs/GRAPH_HARVEST_REFUTED.md), and near-miss pairs have no graph to
    // walk against so they fall back on a string-flank test that was already
    // measured net-negative (flank_match_tol, line ~132).
    //
    // Opt-in while it is being measured: CAPS_DBG=1.
    // min1/min2 = the WEAKEST k-mer count along each path. DiscoSNP++ obtains
    // the equivalent by running kissreads2 as a whole separate pass, mapping
    // every read back onto every bubble. We do not need that pass: a canonical
    // k-mer's count in kc IS the number of reads containing it, so walking the
    // path and taking the minimum gives per-path read support directly. The
    // minimum rather than the mean because a bubble is only as well supported
    // as its weakest link -- one unsupported k-mer means no read spans it.
    struct DbgBubble { std::string flank, path1, path2, lext, rext;
                       uint32_t cov1, cov2; uint32_t min1, min2; int len;
                       uint8_t from_sb = 0;      // 1 = produced by the superbubble walk
                       uint8_t multiallelic = 0; };  // >2 = allele count at a multi-allelic site
    std::vector<DbgBubble> dbg_bubbles;   // hoisted: also read by DBG-ONLY mode below
    if (std::getenv("CAPS_DBG") && (looks_diploid || std::getenv("CAPS_DBG_FORCE"))) {
        auto t_dbg = clk::now();
        // MINC -- ARCHITECTURAL, NOT TUNED.
        //
        // This was previously DERIVED from the k-mer histogram valley (capped
        // at H/3), which sounds principled and is optimising the wrong thing.
        // The valley is the right threshold only when the floor is your ONLY
        // filter. It is not: every bubble is validated downstream against real
        // reads (the coherence pass). In a two-stage design -- permissive
        // generation, then evidence-based filtering -- ALL precision must come
        // from the evidence stage, because precision bought upstream costs
        // recall that no downstream stage can recover.
        //
        // MEASURED, on the window's k-mer count distribution:
        //     cnt=1      60.4% of nodes   (singleton errors)
        //     cnt=2-9     3.9%
        //     cnt>=10    35.6%            (real signal)
        // Nodes surviving:  >=2: 39.6%   >=3: 38.2%   >=4: 37.8%   >=5: 37.6%
        // So raising the floor from 2 to 5 removes 5% of nodes -- no
        // meaningful compute saving -- while deleting exactly the 2-4 count
        // band where LOW-COVERAGE HETEROZYGOUS ALLELES live. The floor's only
        // real job is discarding the 60% singletons.
        //
        // This also explains a failure of the derived version: it moved UP
        // (4 -> 5) at full chr20, precisely when more data makes a lower floor
        // SAFER, because absolute error counts grow and drag the valley with
        // them. That is the wrong direction for the wrong reason.
        //
        // DiscoSNP++ runs -c 3 for the same reason: a permissive graph, with
        // kissreads2 supplying precision. 2 is the smallest defensible value
        // ("seen more than once"), so the floor now does only what it should.
        const uint32_t MINC = std::getenv("CAPS_DBG_MINC")
                            ? (uint32_t)atoi(std::getenv("CAPS_DBG_MINC"))
                            : 2u;
        const int MAXEXT = std::getenv("CAPS_DBG_MAXEXT")
                            ? atoi(std::getenv("CAPS_DBG_MAXEXT")) : 60;
        const uint64_t MASK31 = (~0ULL) >> 2;          // 62 bits = 31 bases
        auto k2s = [](uint64_t v) {                    // 31-mer -> ACGT string
            std::string s(31, 'N');
            for (int i = 30; i >= 0; --i) { s[i] = "ACGT"[v & 3ULL]; v >>= 2; }
            return s;
        };

        // successors of a FORWARD-oriented 31-mer: shift in each base, test the
        // canonical form against kc. Returns forward-oriented successors.
        auto succs = [&](uint64_t fwd, uint64_t out[4], uint32_t cnt[4]) -> int {
            int n = 0;
            for (uint64_t b = 0; b < 4; ++b) {
                const uint64_t nx = ((fwd << 2) | b) & MASK31;
                const KC* e = kc_find(canon31(nx));
                if (e && e->cnt >= MINC) { out[n] = nx; cnt[n] = e->cnt; ++n; }
            }
            return n;
        };
        // Lockstep walk, the direct analogue of Bubble.cpp:509-527. Advances
        // both paths by the SAME base at each step and stops when they land on
        // the same node (bubble closed) or when no single common base extends
        // both (dead end).
        // Also accumulates the weakest k-mer count seen on each path (see the
        // DbgBubble comment): this is per-path read support, obtained during
        // the walk we are already doing rather than from a second pass.
        // Extension with BOUNDED RECURSION over ambiguous continuations.
        //
        // The previous version bailed the instant the two paths had more than
        // one common successor ("if (++nc > 1) return false"). DiscoSNP++ does
        // NOT do this: expand() (Bubble.cpp:598) recurses over every successor
        // pair through expand_heart. So we were MORE permissive than them at
        // the bubble entrance (no checkBranching) and STRICTLY LESS permissive
        // during extension -- and the abort counters showed the cost: of 2,628
        // entrances only 462 closed, 1,856 died mid-walk.
        //
        // A node budget bounds the search: branching factor is at most 4 and
        // depth at most MAXEXT, so unbounded recursion is exponential in the
        // worst case. The budget is shared across the whole search for one
        // bubble, so a tangled region gives up quickly instead of exploding.
        const int WALK_BUDGET = std::getenv("CAPS_DBG_WALKBUDGET")
                              ? atoi(std::getenv("CAPS_DBG_WALKBUDGET")) : 0;
        const int MAXPOLY = std::getenv("CAPS_DBG_MAXPOLY")
                          ? atoi(std::getenv("CAPS_DBG_MAXPOLY")) : 1;
        // WHY 1 IS STRUCTURAL, not a fitted value. MAXPOLY caps how many EXTRA
        // differences a single bubble may absorb. At k=31 and human
        // heterozygosity ~1/1000, the chance a second het site falls inside the
        // same 31 bp window is small and a third is rare, so 1 already captures
        // essentially all genuine clustering; larger values mostly admit paths
        // through repeats. The sweep agrees rather than defines it -- recall
        // saturates by P=2 (0.808 at both 2 and 3) while false positives keep
        // rising (F1 0.861 / 0.850 / 0.850 with coherence, 0.828 / 0.794 /
        // 0.765 without). A genome with much higher heterozygosity would
        // justify a higher cap, which is the property that makes this a
        // reasoned bound rather than a benchmark constant.
        // DiscoSNP++'s own -P default is 3. Measured here, 1 is better on this
        // data both WITHOUT coherence (F1 0.828 / 0.794 / 0.765 for P=1/2/3)
        // and WITH it (0.861 / 0.850 / 0.850): recall saturates by P=2 while
        // false positives keep accumulating. Defaulting to 3 would have made a
        // bare run silently worse than every number reported for this channel.
        std::function<bool(uint64_t,uint64_t,std::string&,std::string&,
                           uint32_t&,uint32_t&,uint64_t&,int,int&,int)> walk_rec;
        walk_rec = [&](uint64_t n1, uint64_t n2, std::string& e1, std::string& e2,
                       uint32_t& mn1, uint32_t& mn2, uint64_t& closenode,
                       int depth, int& budget, int poly) -> bool {
            if (n1 == n2) { closenode = n1; return true; }     // closed
            if (depth >= MAXEXT) return false;                 // too long
            // budget is spent ONLY on ambiguous decisions (see below), so a
            // clean walk of any length costs nothing and the knob means
            // "how many tangles may this bubble cross".
            struct Opt { uint64_t x1, x2; uint32_t k1, k2; char b; };
            Opt opts[4]; int no = 0;
            for (uint64_t b = 0; b < 4; ++b) {
                const uint64_t x1 = ((n1 << 2) | b) & MASK31;
                const uint64_t x2 = ((n2 << 2) | b) & MASK31;
                const KC* a = kc_find(canon31(x1));
                if (!a || a->cnt < MINC) continue;
                const KC* c = kc_find(canon31(x2));
                if (!c || c->cnt < MINC) continue;
                opts[no++] = { x1, x2, a->cnt, c->cnt, "ACGT"[b] };
            }
            // MULTI-POLYMORPHISM EXTENSION (DiscoSNP++'s -P, default 3).
            // If no single base extends BOTH paths, the two haplotypes differ
            // again within the k-mer window -- a second variant close to the
            // first. Requiring identical bases at every step kills the bubble
            // outright, and measurement says that is our single largest source
            // of missed variants: of 115 missed truth SNVs on HG002 r2, 35
            // (30%) die exactly here as "dead end".
            // DiscoSNP++ does not stop: expand_heart is called with
            // nb_polymorphism+1 over successor pairs carrying DIFFERENT bases
            // (Bubble.cpp:636), up to max_polymorphism. We benchmark them at
            // -P 3, so this is a capability we were being measured against
            // without having ported it.
            if (no == 0 && poly < MAXPOLY) {
                for (uint64_t b1 = 0; b1 < 4; ++b1) {
                    const uint64_t x1 = ((n1 << 2) | b1) & MASK31;
                    const KC* a = kc_find(canon31(x1));
                    if (!a || a->cnt < MINC) continue;
                    for (uint64_t b2 = 0; b2 < 4; ++b2) {
                        if (b2 == b1) continue;            // same-base case handled above
                        const uint64_t x2 = ((n2 << 2) | b2) & MASK31;
                        const KC* c = kc_find(canon31(x2));
                        if (!c || c->cnt < MINC) continue;
                        std::string s1 = e1 + "ACGT"[b1], s2 = e2 + "ACGT"[b2];
                        uint32_t m1 = std::min(mn1, a->cnt), m2 = std::min(mn2, c->cnt);
                        if (walk_rec(x1, x2, s1, s2, m1, m2, closenode,
                                     depth + 1, budget, poly + 1)) {
                            e1 = s1; e2 = s2; mn1 = m1; mn2 = m2;
                            return true;
                        }
                    }
                }
            }
            if (no == 0) return false;                         // dead end
            if (no > 1 && --budget < 0) return false;          // out of tangle allowance
            // Try the best-supported continuation first: on a real bubble the
            // true path carries genuine coverage, so this finds the closure
            // sooner and spends less budget on error branches.
            for (int i = 1; i < no; ++i)
                for (int j = i; j > 0 && std::min(opts[j].k1, opts[j].k2)
                                       > std::min(opts[j-1].k1, opts[j-1].k2); --j)
                    std::swap(opts[j], opts[j-1]);
            for (int i = 0; i < no; ++i) {
                std::string s1 = e1 + opts[i].b, s2 = e2 + opts[i].b;
                uint32_t m1 = std::min(mn1, opts[i].k1), m2 = std::min(mn2, opts[i].k2);
                if (walk_rec(opts[i].x1, opts[i].x2, s1, s2, m1, m2,
                             closenode, depth + 1, budget, poly)) {
                    e1 = s1; e2 = s2; mn1 = m1; mn2 = m2;
                    return true;
                }
            }
            return false;
        };
        auto walk = [&](uint64_t n1, uint64_t n2, std::string& e1, std::string& e2,
                        uint32_t& mn1, uint32_t& mn2, uint64_t& closenode) -> bool {
            int budget = WALK_BUDGET;
            return walk_rec(n1, n2, e1, e2, mn1, mn2, closenode, 0, budget, 0);
        };

        // ── UNITIG EXTENSION ────────────────────────────────────────────────
        // MEASURED MOTIVATION, not a guess: without this a bubble sequence is
        // flank(31) + path, mean 63 bp on HG002 r2. bwa cannot place a 63-mer
        // uniquely in a 63 Mb chromosome, so found bubbles were being lost at
        // the LIFT rather than at the calling -- 397 bubbles produced only 243
        // lifted calls. DiscoSNP++ does not have this problem because it
        // extends each bubble into its surrounding unitig before output
        // (expand_one_simple_path / its -t traversal); this is a layer where we
        // were strictly worse than the tool we are competing with.
        //
        // Extension walks the graph while the path is UNAMBIGUOUS -- exactly
        // one successor (or predecessor) above the coverage floor -- which is
        // the standard unitig definition. It stops at the first junction, so it
        // never invents sequence across a decision point.
        const int MAXFLANK = std::getenv("CAPS_DBG_FLANK")
                           ? atoi(std::getenv("CAPS_DBG_FLANK")) : 1000;
        // 1000 is what every validated Method B measurement used. Extension
        // length is mapping-only (measured neutral on F1), but the default must
        // match the tested configuration or a bare run is not the tested one.
        // ── TIP DETECTION ───────────────────────────────────────────────────
        // A sequencing error near the end of a read creates a short DEAD END
        // in the graph -- a "tip". The assembly literature is explicit that
        // such error branches are STERILE: they do not reconnect to anything.
        // Our traversals stop the moment a node has more than one successor,
        // so before this they were halting at noise rather than at real
        // branch points, which is why unitig extension reached only 238 bp and
        // why 46-66% of bubbles failed to chain (measured: 34% chained on r2,
        // 54% on r3 -- the rate is data-dependent, so a fixed assumption about
        // it would not generalise).
        //
        // is_tip walks a candidate branch and reports whether it dies within
        // TIPLEN steps. That is a MEASURED property of the branch, not a
        // tuned constant: real sequence continues, error branches do not.
        const int TIPLEN = std::getenv("CAPS_DBG_TIPLEN")
                         ? atoi(std::getenv("CAPS_DBG_TIPLEN")) : 3 * 31;
        std::function<bool(uint64_t,int)> is_tip = [&](uint64_t node, int budget) -> bool {
            for (int i = 0; i < budget; ++i) {
                int nb = 0; uint64_t pick = 0;
                for (uint64_t b = 0; b < 4; ++b) {
                    const uint64_t w = ((node << 2) | b) & MASK31;
                    const KC* e = kc_find(canon31(w));
                    if (e && e->cnt >= MINC) { ++nb; pick = w; }
                }
                if (nb == 0) return true;            // died: it is a tip
                if (nb > 1)  return false;           // branches again: real structure
                node = pick;
            }
            return false;                            // survived the budget: real
        };
        // Successors with sterile tips removed. Returns the count of surviving
        // branches; a junction that is really one real path plus error tips
        // collapses to 1 and traversal continues through it.
        auto succs_live = [&](uint64_t fwd, uint64_t out[4], uint32_t cnt[4]) -> int {
            uint64_t raw[4]; uint32_t rc_[4];
            const int n = succs(fwd, raw, rc_);
            if (n <= 1) { for (int i = 0; i < n; ++i) { out[i] = raw[i]; cnt[i] = rc_[i]; } return n; }
            int k = 0;
            for (int i = 0; i < n; ++i)
                if (!is_tip(raw[i], TIPLEN)) { out[k] = raw[i]; cnt[k] = rc_[i]; ++k; }
            // If every branch looked like a tip, keep the original set rather
            // than silently deleting a real (if short) structure.
            if (k == 0) { for (int i = 0; i < n; ++i) { out[i] = raw[i]; cnt[i] = rc_[i]; } return n; }
            return k;
        };

        auto ext_right = [&](uint64_t node, int maxlen) {
            std::string s;
            for (int i = 0; i < maxlen; ++i) {
                uint64_t so_[4]; uint32_t sc_[4];
                if (succs_live(node, so_, sc_) != 1) break;   // real junction or dead end
                s += "ACGT"[so_[0] & 3ULL]; node = so_[0];
            }
            return s;
        };
        auto ext_left = [&](uint64_t node, int maxlen) {
            // predecessor: prepend a base and drop the last one. pack31 puts
            // base 0 in bits 60-61, so that is (node >> 2) | (b << 60).
            std::string s;
            for (int i = 0; i < maxlen; ++i) {
                int nb = 0; uint64_t pick = 0, pb = 0;
                for (uint64_t b = 0; b < 4; ++b) {
                    const uint64_t w = (node >> 2) | (b << 60);
                    const KC* e = kc_find(canon31(w));
                    if (e && e->cnt >= MINC) { if (++nb > 1) break; pick = w; pb = b; }
                }
                if (nb != 1) break;
                s += "ACGT"[pb]; node = pick;
            }
            std::reverse(s.begin(), s.end());          // walked backwards
            return s;
        };

        // ── WHY INDEL RECALL IS BOUNDED, AND IT IS NOT A TUNING PROBLEM ────
        // Traced every missed indel on HG002 r2 back into kc. At the missed
        // sites the branching node has GOOD coverage (12-30) but exactly ONE
        // successor -- and that holds with NO coverage floor at all: raw
        // successor counts come back as [0,0,21,0]. The alternate haplotype's
        // k-mers are not below MINC, they are ABSENT.
        //
        // The reason: 66% of the missed indels are homopolymer-length changes
        // and most of the rest are tandem-repeat expansions. Inserting a base
        // into a run of the same base does not create a distinct k-mer path --
        // it makes the run longer. Probing the reads directly confirms it: for
        // a G->GT call inside a T-run, the REF and ALT probes are the SAME
        // STRING (TCTGGGTTTTTGTTTTTCGGGTTTTTTTTTTT), each found in 6 reads.
        //
        // There is no second path, so no bubble exists -- for us OR for any
        // de Bruijn caller, DiscoSNP++ included. This is a representational
        // limit of the graph, not a defect in the traversal, and it is why
        // nine successive attempts at filters, gates and orderings all failed
        // to move indel recall. Resolving these events needs length-aware
        // evidence (read pileup depth over the run), which is a different
        // mechanism from bubble finding.
        //
        // ── SUPERBUBBLE SEARCH (CAPS_DBG_SB=1) ──────────────────────────────
        // Strictly generalises the pairwise lockstep walk above, which is a
        // faithful port of DiscoSNP++ (Bubble.cpp:509-527) and inherits two
        // hard limits from it:
        //   * EXACTLY TWO alleles -- it advances a pair of nodes in lockstep;
        //   * NO INDELS -- lockstep requires the SAME base to extend both
        //     paths, so two paths of different length cannot be represented.
        //     DiscoSNP++ needs an entirely separate breadth-first procedure
        //     (start_indel_prediction, Bubble.cpp:200-270) to get indels at all.
        //
        // A superbubble (Onodera et al. 2013) is defined by a single entrance s
        // and single exit t such that every path leaving s stays inside until
        // it reaches t. Enumerating the paths s->t yields the alleles directly:
        // any number of them, of any lengths. One operation replaces seeding +
        // extension + closure, and covers SNVs, indels and multi-allelic sites
        // together.
        //
        // COST CONTROL, and why this is exact rather than a shortcut: a
        // superbubble's entrance must have out-degree >= 2, so restricting the
        // search to branching nodes loses nothing. Global enumeration is what
        // makes tools like BubbleGun expensive (~25 min / 22 GB on a 22M-node
        // human graph); starting only from the 2,894 branching nodes in this
        // window, with a bounded walk per branch, keeps it in the same
        // fractions of a second the pairwise walk costs.
        const bool WANT_SB = std::getenv("CAPS_DBG_SB") != nullptr;
        const int SB_MAXPATH = std::getenv("CAPS_DBG_SBLEN")
                             ? atoi(std::getenv("CAPS_DBG_SBLEN")) : 60;
        // node budget per superbubble search; bounds worst-case work
        const int SB_BUDGET = std::getenv("CAPS_DBG_SBBUDGET")
                            ? atoi(std::getenv("CAPS_DBG_SBBUDGET")) : 200;
        struct SbAllele { std::string seq; uint32_t minsup; };
        // entry/exit node ids are kept so consecutive superbubbles can be
        // CHAINED (see the bubble-chain block after the scan).
        struct SbSite { std::string flank, lext, rext; std::vector<SbAllele> alleles;
                        uint64_t entry = 0, exit = 0; };
        std::vector<SbSite> sb_sites;
        // Follow one branch as a unitig, stopping at the first junction, and
        // record every node visited with the sequence that reached it.
        auto branch_walk = [&](uint64_t start, char firstbase,
                               std::vector<std::pair<uint64_t,std::string>>& trace,
                               uint32_t seedcnt) -> uint32_t {
            uint32_t mn = seedcnt;
            uint64_t node = start;
            std::string acc(1, firstbase);
            trace.emplace_back(node, acc);
            for (int i = 0; i < SB_MAXPATH; ++i) {
                uint64_t so_[4]; uint32_t sc_[4];
                if (succs_live(node, so_, sc_) != 1) break;   // real junction or dead end
                if (sc_[0] < mn) mn = sc_[0];
                acc += "ACGT"[so_[0] & 3ULL];
                node = so_[0];
                trace.emplace_back(node, acc);
            }
            return mn;
        };

        // Predecessors (tip-filtered), needed for the "all parents visited"
        // condition of the real superbubble algorithm.
        auto preds_live = [&](uint64_t fwd, uint64_t out[4]) -> int {
            int n = 0;
            for (uint64_t b = 0; b < 4; ++b) {
                const uint64_t w = (fwd >> 2) | (b << 60);
                const KC* e = kc_find(canon31(w));
                if (e && e->cnt >= MINC) out[n++] = w;
            }
            return n;
        };
        // ── ONODERA SUPERBUBBLE DETECTION (faithful port) ───────────────────
        // Reference: BubbleGun/find_bubbles.py, implementing Onodera et al.
        // 2013. The earlier version here walked each branch as a unitig and
        // took the earliest node all branches reached -- that is NOT this
        // algorithm. It stops at any junction, so it can only ever find
        // parallel simple paths and never explores a structure with internal
        // branching (which is why multi-allelic sites came out as 0).
        //
        // The load-bearing condition is "a node enters the frontier only once
        // ALL of its parents have been visited". That is what allows nodes
        // reachable by several internal paths to be handled correctly.
        //   S      = frontier: discovered AND all parents visited
        //   seen   = discovered but not yet ready
        //   sink   = |S| == 1 and |seen| == 0  (only the exit remains)
        // Aborts on a tip (no successors) or on returning to the entrance.
        // abort-reason counters: 2,628 branching nodes yield only 282 bubbles,
        // so ~89% fail somewhere. Which reason dominates decides where the
        // missing recall is (ours 0.468 vs DiscoSNP++ 0.763 from the same
        // k-mer set), rather than guessing.
        size_t ab_tip = 0, ab_cycle = 0, ab_budget = 0, ab_exhaust = 0, ab_ok = 0;
        auto find_sb = [&](uint64_t s0, uint64_t& t_out) -> bool {
            std::unordered_set<uint64_t> visited, seen;
            std::vector<uint64_t> S; S.push_back(s0);
            int guard = 0;
            while (!S.empty()) {
                if (++guard > SB_BUDGET) { ++ab_budget; return false; }
                const uint64_t v = S.back(); S.pop_back();
                visited.insert(v); seen.erase(v);
                uint64_t ch[4]; uint32_t cc[4];
                const int nc = succs_live(v, ch, cc);
                if (nc == 0) { ++ab_tip; return false; }   // tip
                for (int i = 0; i < nc; ++i) {
                    const uint64_t u = ch[i];
                    if (u == s0) { ++ab_cycle; return false; }  // cycle back to entrance
                    if (visited.count(u)) continue;
                    uint64_t pa[4];
                    const int np = preds_live(u, pa);
                    bool allv = true;
                    for (int j = 0; j < np; ++j)
                        if (!visited.count(pa[j])) { allv = false; break; }
                    if (allv) { S.push_back(u); seen.erase(u); }
                    else seen.insert(u);
                }
                if (S.size() == 1 && seen.empty()) { t_out = S[0]; ++ab_ok; return true; }
            }
            ++ab_exhaust;
            return false;
        };
        // BOTH ORIENTATIONS. kc stores CANONICAL k-mers, so each entry stands
        // for BOTH strands -- the graph is bidirected and every node has two
        // sides. Extending node.kmer as though it were forward sequence
        // explores only the side the canonical form happens to match; for the
        // ~half of nodes whose canonical form is the reverse complement that
        // is the wrong side, and no bubble is ever found from them.
        // GATB carries an explicit strand on every Node and normalises with
        // `(strand == STRAND_FORWARD) ? val : revcomp(val)` (Graph.cpp:1599).
        // This code had no notion of strand at all, which is the single
        // largest reason graph-only recall (0.512) trailed DiscoSNP++'s
        // (0.763) on the very same k-mer set.
        // PARALLEL. Every node is examined independently: kc, MINC and the
        // walk lambdas are read-only, and the only shared writes are the
        // bubble list (accumulated per-thread, concatenated once) and the
        // branching counter (a reduction). Bubble ORDER changes, which is
        // harmless -- each bubble is self-contained and downstream code
        // addresses them by index into the same vector it dumps.
        // WANT_SB is left on the serial path: find_sb's abort counters are
        // shared diagnostics and superbubble mode is not the default.
        size_t branching = 0;
        #pragma omp parallel if(!WANT_SB)
        {
        std::vector<DbgBubble> local_bubbles;
        #pragma omp for schedule(dynamic, 2048) reduction(+:branching) nowait
        for (long long ni = 0; ni < (long long)kc.size(); ++ni) {
            const KC& node0 = kc[(size_t)ni];
            if (node0.cnt < MINC) continue;
            for (int orient = 0; orient < 2; ++orient) {
            KC node = node0;
            if (orient) {
                node.kmer = rc31(node0.kmer);
                if (node.kmer == node0.kmer) continue;      // palindrome: one side only
            }
            uint64_t so[4]; uint32_t sc[4];
            const int ns = succs(node.kmer, so, sc);
            if (ns < 2) continue;
            ++branching;

            // ── superbubble path: all branches at once, any lengths ─────────
            if (WANT_SB) {
                std::vector<std::vector<std::pair<uint64_t,std::string>>> traces(ns);
                std::vector<uint32_t> mins(ns);
                for (int i = 0; i < ns; ++i)
                    mins[i] = branch_walk(so[i], "ACGT"[so[i] & 3ULL], traces[i], sc[i]);
                // The exit t is the earliest node every branch reaches. Taking
                // the EARLIEST keeps the variable region minimal, which is the
                // minimality condition in the superbubble definition.
                uint64_t bestT = 0; size_t bestCost = SIZE_MAX;
                std::vector<std::string> bestSeq;
                // Locate the exit with the REAL superbubble algorithm
                // (find_sb, the Onodera/BubbleGun port), then read each
                // allele's string off the branch trace that reaches it.
                // The previous version picked "the earliest node every branch
                // reaches", which is an approximation: it ignores the
                // all-parents-visited condition and therefore never explores a
                // structure with internal branching.
                // ── VARIANT GROUPS (SKA lo, MBE 2025) ───────────────────────
                // MEASURED PROBLEM: of 755 branching nodes with 3-4 successors,
                // find_sb fails on 734 (97%) and none of the 21 successes has
                // more than 2 branches converging. Multi-allelic sites die
                // inside find_sb, whose exit test needs the frontier to drain
                // to exactly one node with nothing pending -- a condition that
                // 3-4 branches through repeat-rich sequence essentially never
                // satisfy. That is structural to the Onodera formulation.
                //
                // SKA lo solves the same problem by relaxing WHAT COUNTS AS A
                // SITE: "a variant group is identified when AT LEAST TWO paths
                // start from the same entry node and end at the same exit
                // node". Not all paths -- at least two. The traces already hold
                // every node each branch visits, so the group can be read
                // straight out of them: pick the earliest node reached by the
                // most branches (>=2). A 4-way node where 3 branches meet is
                // then a 3-allele site instead of a discarded one.
                //
                // Tried first and only when the strict superbubble finds
                // nothing, so simple two-path bubbles keep the tighter
                // definition and SNV calling is unaffected.
                //
                // OUTCOME, MEASURED: this DOES form variant groups (sites
                // 818 -> 827) but ZERO of them have more than two allele paths.
                // The group node is selected to maximise the number of
                // converging branches, so a 3-branch meeting point would win if
                // one existed. None does.
                //
                // The reason is biological, not algorithmic: at 30x DIPLOID
                // coverage there is no third haplotype for a third branch to
                // converge on. The 3-4 successor nodes are error or repeat
                // branches, which die out -- which is exactly what MINC and the
                // tip filter exist to make them do. T5.2 is therefore not
                // reachable on diploid human data by any traversal change; it
                // needs a sample that genuinely carries >2 haplotypes (the
                // tetraploid construction of T5.3, or a pooled/polyploid
                // sample). Kept because it is strictly more general than the
                // strict superbubble and costs nothing when it finds nothing.
                uint64_t sbT = 0;
                bool used_group = false;
                if (ns > 2 && !std::getenv("CAPS_DBG_NOVGROUP")) {
                    std::unordered_map<uint64_t, std::pair<int,size_t>> reach; // node -> (count, earliest depth sum)
                    for (int i = 0; i < ns; ++i) {
                        std::unordered_set<uint64_t> seen_i;
                        for (size_t q = 0; q < traces[i].size(); ++q) {
                            const uint64_t nd = traces[i][q].first;
                            if (!seen_i.insert(nd).second) continue;
                            auto& e = reach[nd];
                            e.first += 1; e.second += q;
                        }
                    }
                    // SHORTEST CLOSURE WINS, then branch count.
                    //
                    // This previously maximised BRANCH COUNT first and used
                    // depth only as a tiebreak, so a distant node reached by 3
                    // branches beat a near one reached by 2. DiscoSNP++ orders
                    // it the other way round and never proposes the distant
                    // exit at all: their indel BFS stops at the smallest
                    // closing event (Bubble.cpp:238-243, "if (insert_size ==
                    // found_del_size-1) break" and "if (insert_size >
                    // found_del_size) continue"), which is why long spurious
                    // indels are absent from their output rather than filtered
                    // out of it.
                    //
                    // MEASURED AFTERWARDS: this ordering change is INERT in
                    // practice. Variant groups only fire when ns > 2, which is
                    // rare, so nearly every indel comes from find_sb on a
                    // 2-successor node -- and find_sb returns the FIRST node
                    // satisfying its drain condition, which is already the
                    // earliest in traversal order. Closure distances measured
                    // 709 of 827 in the 32-39 bp bin, i.e. about one k-mer past
                    // the divergence, which is the correct shape for real
                    // bubbles rather than evidence of distant closure. The
                    // long-indel false positives therefore do NOT come from
                    // choosing a far exit, and this ordering is kept only
                    // because it matches the reference implementation.
                    //
                    // Adopting their ORDERING (a generative bound) rather than
                    // adding another post-hoc length filter -- four of those
                    // were tried and every one cost more true positives than
                    // false ones, twice catastrophically (F1 0.537 -> 0.035).
                    // Preferring the nearest reconvergence keeps short events
                    // and simply does not generate the long ones.
                    const bool SHORTEST_FIRST = std::getenv("CAPS_DBG_DEEPFIRST") == nullptr;
                    int bestN = 0; size_t bestDepth = SIZE_MAX; uint64_t bestNode = 0;
                    for (const auto& kv : reach) {
                        if (kv.second.first < 2) continue;
                        bool better;
                        if (SHORTEST_FIRST)
                            better = (kv.second.second < bestDepth) ||
                                     (kv.second.second == bestDepth && kv.second.first > bestN);
                        else
                            better = (kv.second.first > bestN) ||
                                     (kv.second.first == bestN && kv.second.second < bestDepth);
                        if (better) {
                            bestN = kv.second.first; bestDepth = kv.second.second; bestNode = kv.first;
                        }
                    }
                    if (bestN >= 2) { sbT = bestNode; used_group = true; }
                    // SHORTEST CLOSURE. Among nodes reached by the same number
                    // of branches, `bestDepth` already prefers the earliest --
                    // that is the analogue of DiscoSNP++'s shortest-first BFS,
                    // and it is what stops the walk proposing a distant exit
                    // through a repeat when a near one exists.
                }
                if (used_group || find_sb(node.kmer, sbT)) {
                    // PARTIAL CONVERGENCE (kept, but it is NOT what blocks T5.2).
                    // MEASURED: of 755 branching nodes with 3-4 successors,
                    // find_sb FAILS on 734 (97%), and of the 21 that succeed
                    // none has more than 2 branches converging. Multi-allelic
                    // sites are lost inside find_sb itself -- its exit
                    // condition needs the frontier to drain to exactly one node
                    // with nothing pending, which 3-4 branches through
                    // repeat-rich sequence essentially never satisfy. That is a
                    // structural property of the Onodera formulation, not a
                    // threshold, so T5.2 is not reachable by relaxing this
                    // check. Kept anyway because it is strictly more
                    // information and cannot affect SNVs.
                    //
                    // This used to require EVERY branch to reach the exit
                    // (`if (!found) ok = false`) and discarded the whole site
                    // otherwise. Measured consequence: 755 branching nodes have
                    // 3-4 successors, yet not one produced a 3+-allele site --
                    // every multi-allelic candidate was thrown away because one
                    // of its branches died in a tip or a repeat.
                    //
                    // A 4-successor node where 3 branches converge and 1 dies
                    // is a 3-ALLELE SITE, not a failure. Keeping the branches
                    // that DO converge is strictly more information than
                    // discarding all of them, and it cannot affect SNV calling
                    // because SNVs come from the pairwise lockstep walk, which
                    // is untouched by this. At least 2 convergent branches are
                    // required -- one path is not a bubble.
                    std::vector<std::string> seqs_;
                    for (int i = 0; i < ns; ++i) {
                        for (size_t p = 0; p < traces[i].size(); ++p)
                            if (traces[i][p].first == sbT) {
                                seqs_.push_back(traces[i][p].second); break;
                            }
                    }
                    const bool ok = (seqs_.size() >= 2);
                    if (ok) {
                        bestT = sbT; bestCost = 0; bestSeq = seqs_;
                    }
                }
                if (bestCost != SIZE_MAX && bestSeq.size() >= 2) {
                    // Distinct allele sequences only: if every branch spells the
                    // same string this is not a variant, just a graph artefact.
                    std::set<std::string> uniq(bestSeq.begin(), bestSeq.end());
                    if (uniq.size() >= 2) {
                        SbSite st;
                        st.flank = k2s(node.kmer);
                        st.lext  = ext_left(node.kmer, MAXFLANK);
                        st.rext  = ext_right(bestT, MAXFLANK);
                        st.entry = node.kmer; st.exit = bestT;
                        // mins[] is indexed by branch; bestSeq holds only the
                        // branches that converged, so pair them by position and
                        // fall back to the site's weakest support if the branch
                        // index is no longer aligned.
                        for (size_t i = 0; i < bestSeq.size(); ++i)
                            st.alleles.push_back({bestSeq[i],
                                                  i < mins.size() ? mins[i] : 0u});
                        sb_sites.push_back(std::move(st));
                    }
                }
                // DO NOT skip the pairwise walk.
                //
                // This used to `continue`, so enabling the superbubble path
                // REPLACED the pairwise walk. Measured cost of that
                // substitution on HG002 r2: SNV F1 0.881 -> 0.791 (-29 TP,
                // +42 FP) -- the superbubble finds indels the lockstep walk
                // structurally cannot, but it is WORSE at plain SNVs, because
                // its per-branch unitig traces are a coarser view of a simple
                // two-path bubble than walking the two paths in lockstep.
                //
                // The two are complementary, not alternatives:
                //   pairwise walk -> SNVs   (lockstep, same base both paths)
                //   superbubble   -> INDELS (paths of different length)
                // so both run at every branching node and each contributes the
                // variant class it is actually good at. Superbubble records are
                // tagged so the emitter can keep only its different-length
                // ones; see SB_ONLY_INDEL below.
                //
                // (`CAPS_DBG_SB_REPLACE=1` restores the old substituting
                // behaviour for re-measurement.)
                if (std::getenv("CAPS_DBG_SB_REPLACE")) continue;
            }

            for (int i = 0; i < ns; ++i)
                for (int j = i + 1; j < ns; ++j) {
                    // REFUTED: a cheap coverage-ratio reject before the walk.
                    // Rationale was sound -- 1,226,535 branching nodes yield only
                    // 111,499 bubbles, so 91% of walks are wasted and traversal is
                    // 55.5 s of 98.5 s. But measured at full chr20 it changed
                    // nothing (98.7 s vs 98.5 s, traversal 55.25 s vs 55.48 s) and
                    // cost 300 bubbles. Reason: the test was gated on ns > 2, and
                    // almost every branching node has exactly 2 successors, so it
                    // essentially never fired. The waste is not in walks that a
                    // coverage test can predict -- it is in walks that fail to
                    // CLOSE, which is only knowable by walking.
                    std::string e1, e2;
                    uint32_t mn1 = sc[i], mn2 = sc[j];   // seed with the branch bases
                    uint64_t closenode = 0;
                    if (!walk(so[i], so[j], e1, e2, mn1, mn2, closenode)) continue;
                    // The two divergent first bases plus the shared closing
                    // extension: this is the bubble's variable region.
                    DbgBubble bb;
                    bb.flank = k2s(node.kmer);
                    bb.path1 = std::string(1, "ACGT"[so[i] & 3ULL]) + e1;
                    bb.path2 = std::string(1, "ACGT"[so[j] & 3ULL]) + e2;
                    bb.cov1 = sc[i]; bb.cov2 = sc[j];
                    bb.min1 = mn1;   bb.min2 = mn2;
                    bb.len = (int)e1.size();
                    bb.lext = ext_left(node.kmer, MAXFLANK);
                    bb.rext = ext_right(closenode, MAXFLANK);
                    local_bubbles.push_back(std::move(bb));
                }
            }   // end orientation loop
        }
        #pragma omp critical(dbgbub)
        for (auto& b : local_bubbles) dbg_bubbles.push_back(std::move(b));
        }   // end parallel
        // ── READ-COHERENCE: the fragment must exist in REAL READS ───────────
        // THE FP MECHANISM, named in this repo's own docs: a CHIMERIC JUNCTION.
        // Every k-mer along a bubble path can be present in kc while the FULL
        // path exists in no single read -- each k-mer contributed by a
        // different read, joined at a shared anchor. k-mer counts are blind to
        // this by construction, which is why our precision (0.865) trails
        // DiscoSNP++'s (0.951) even though our recall is now higher.
        //
        // eBWT2SNP's precision guarantee is exactly this test: every emitted
        // fragment of length 2k+1 must be an actual SUBSTRING of at least C
        // real reads. DiscoSNP++ gets the equivalent from kissreads2, a whole
        // separate tool that maps every read back onto every bubble.
        //
        // WE DO NOT NEED THAT TOOL. The reads are already in memory (`seqs`),
        // so one pass with a small probe index gives a stronger check than
        // theirs at a fraction of the cost: index each bubble by the 31-mer
        // SPANNING its variant, sweep the reads once, and verify the full
        // fragment really occurs. This is a case where having the reads
        // resident -- a consequence of being a compressor -- buys something a
        // standalone caller has to pay a separate pass for.
        if (!std::getenv("CAPS_DBG_NOCOH") && !dbg_bubbles.empty()) {
            auto t_coh = clk::now();
            // COHC -- required read support -- DERIVED FROM MEASURED COVERAGE.
            //
            // This was a constant 2, swept on HG002 r2 at 30x. That does not
            // generalise: "how many reads must contain this fragment" scales
            // with depth. At 10x, 2 reads is a large fraction of a haplotype's
            // coverage and the filter becomes punitive; at 100x, 2 reads is
            // noise and the filter stops filtering.
            //
            // H is the measured haploid k-mer depth (the single-copy histogram
            // peak, computed above). A heterozygous allele is covered by ~H
            // reads, so requiring a small FRACTION of H tracks the data.
            // H/10 with a floor of 2 reproduces the swept value exactly at
            // this benchmark's depth (H=20 -> 2) while adapting elsewhere,
            // which is what makes it a formula rather than a fitted constant
            // that happens to be right once.
            //
            // NOTE the contrast with MINC, where a per-dataset derivation
            // (the k-mer histogram valley) generalised WORSE than a structural
            // constant: it drifted upward at full scale and cost 0.057 recall.
            // Derivation is right when the quantity genuinely scales with the
            // data -- read support does, the singleton-error floor does not.
            const int COHC = std::getenv("CAPS_DBG_COHC")
                           ? atoi(std::getenv("CAPS_DBG_COHC"))
                           : std::max(2, (int)(H / 10));
            // Minimum phred for a read's variant base to count as support.
            // 20 = 1% error probability, the standard confident-base cutoff.
            const int MINQ = std::getenv("CAPS_DBG_MINQ")
                           ? atoi(std::getenv("CAPS_DBG_MINQ")) : 20;
            // Fragment half-width. The containment test requires a read to
            // contain the WHOLE fragment (2*HALFW+1 = 63 bp at HALFW=31), which
            // for indels is the binding constraint -- a 148 bp read must cover
            // 63 bases centred on the event, so genuine indels near read ends
            // are lost. SNVs tolerate it because they are far more numerous.
            const int HALFW = std::getenv("CAPS_DBG_HALFW")
                            ? atoi(std::getenv("CAPS_DBG_HALFW")) : 31;
            struct Probe { uint32_t bi; uint8_t path; };
            std::unordered_map<uint64_t, std::vector<Probe>> probe;
            std::vector<std::string> frag1(dbg_bubbles.size()), frag2(dbg_bubbles.size());
            for (uint32_t bi = 0; bi < (uint32_t)dbg_bubbles.size(); ++bi) {
                const DbgBubble& b = dbg_bubbles[bi];
                if (b.flank.size() < 31 || b.path1.empty() || b.path2.empty()) continue;
                // fragment: left context (flank) + path + right context (rext)
                const std::string base = b.flank + b.path1 + b.rext;
                const std::string alt  = b.flank + b.path2 + b.rext;
                const size_t vpos = 31;                 // variant offset in both
                const size_t lo = (vpos >= (size_t)HALFW) ? vpos - HALFW : 0;
                const size_t w1 = std::min(base.size() - lo, (size_t)(2 * HALFW + 1));
                const size_t w2 = std::min(alt.size()  - lo, (size_t)(2 * HALFW + 1));
                if (w1 < 32 || w2 < 32) continue;
                frag1[bi] = base.substr(lo, w1);
                frag2[bi] = alt.substr(lo, w2);
                // key on the 31-mer spanning the variant on each path
                uint64_t k1, k2;
                // PROBE MUST STRADDLE THE JUNCTION.
                //
                // koff was vpos (=31), i.e. exactly the start of the path, so
                // the probe k-mer covered path+rext and NEVER included the
                // flank. For an SNV that is harmless -- the variant base is at
                // position 0 of the k-mer, so a read carrying the allele still
                // matches. For an INDEL it removes the whole discriminating
                // signal: what separates a real indel from a repeat's distant
                // closure is whether reads actually cross the FLANK->PATH
                // junction, and that junction was never in the probe.
                //
                // Centring the k-mer on the junction (half in the flank, half
                // in the path) means a read only matches if it genuinely spans
                // the event. This is the same principle Lancet and Scalpel
                // apply -- "evaluate the support of each alternate haplotype
                // with respect to the RAW READ DATA" -- applied at the place
                // that actually distinguishes indels.
                // TWO PROBES: one at the path start (as before) and one
                // centred on the junction. A read matches if EITHER is present,
                // so the junction evidence is available without demanding that
                // every supporting read span 15 bases of flank as well as the
                // event -- which is what a single centred probe required, and
                // it cost more recall than the precision it bought
                // (P 0.733 -> 0.791 but F1 0.611 -> 0.535, measured).
                const size_t koff  = vpos - lo;                              // at path start
                const size_t koff2 = (vpos >= lo + 15) ? (vpos - lo - 15) : 0; // straddling
                if (koff + 31 <= frag1[bi].size() && pack31(frag1[bi].data() + koff, k1))
                    probe[canon31(k1)].push_back({bi, 1});
                if (koff + 31 <= frag2[bi].size() && pack31(frag2[bi].data() + koff, k2))
                    probe[canon31(k2)].push_back({bi, 2});
                if (koff2 != koff) {
                    uint64_t j1, j2;
                    if (koff2 + 31 <= frag1[bi].size() && pack31(frag1[bi].data() + koff2, j1))
                        probe[canon31(j1)].push_back({bi, 1});
                    if (koff2 + 31 <= frag2[bi].size() && pack31(frag2[bi].data() + koff2, j2))
                        probe[canon31(j2)].push_back({bi, 2});
                }
            }
            fprintf(stderr, "[DBG-RSS] probe built: %ldMB  probe_keys=%zu\n", rss_kb()/1024, probe.size());
            std::vector<uint32_t> sup1(dbg_bubbles.size(), 0), sup2(dbg_bubbles.size(), 0);
            // per-position coverage, capped so memory stays bounded
            // ONE pass over the reads, both strands.
            // PARALLEL. `probe`, `frag1/2` and `seqs` are read-only here; the
            // only writes are the two support counters, which each thread
            // accumulates privately and adds in at the end. Addition is
            // commutative, so the totals are identical to the serial sweep.
            // Per-thread arrays are one uint32 per bubble (a few hundred), so
            // the extra memory is negligible.
            // ── EXCLUSIVE READ ASSIGNMENT (CAPS_DBG_EXCL=1) ──────────────
            // The containment test above asks "does this read contain the
            // fragment", separately for each path, and credits BOTH if both
            // succeed. For an SNV that is harmless: the two fragments are the
            // same length and differ at one base, so a read essentially never
            // contains both, and this flag is inert by construction.
            //
            // For a LENGTH-CHANGING bubble it is not harmless. The two
            // fragments share a flank and a right context and differ by the
            // indel, so in repetitive or homopolymeric sequence one read can
            // satisfy both `find`s -- and then a site with no allelic evidence
            // at all still clears COHC on both paths. That is a plausible
            // source of the precision asymmetry we actually observe
            // (SNV 0.950 vs INDEL 0.658, against DiscoSNP++'s indel 0.936):
            // kissreads2 assigns each read to its BEST path rather than to
            // every path it contains.
            //
            // The rule here has no parameter to tune: a read that matches both
            // paths distinguishes neither, so it is evidence for neither.
            // Reads are stamped per-bubble within one read, then resolved once
            // both strands have been scanned.
            const bool EXCL = std::getenv("CAPS_DBG_EXCL") != nullptr;
            std::vector<uint32_t> ambig(dbg_bubbles.size(), 0);
            #pragma omp parallel
            {
                std::vector<uint32_t> t1(dbg_bubbles.size(), 0), t2(dbg_bubbles.size(), 0);
                std::vector<uint32_t> tamb(EXCL ? dbg_bubbles.size() : 0, 0);
                std::vector<uint32_t> s1(EXCL ? dbg_bubbles.size() : 0, 0),
                                      s2(EXCL ? dbg_bubbles.size() : 0, 0),
                                      sT(EXCL ? dbg_bubbles.size() : 0, 0);
                std::vector<uint32_t> touched;
                std::string ubuf3;
                #pragma omp for schedule(dynamic, 256)
                for (long long ri = 0; ri < (long long)seqs.size(); ++ri) {
                    const std::string& q = unpack_read(seqs[(size_t)ri], ubuf3);
                    if (q.size() < 31) continue;
                    const std::string qr = rc_str(q);
                    for (int strand = 0; strand < 2; ++strand) {
                        const std::string& r = strand ? qr : q;
                        for (size_t off = 0; off + 31 <= r.size(); ++off) {
                            uint64_t v;
                            if (!pack31(r.data() + off, v)) continue;
                            auto it = probe.find(canon31(v));
                            if (it == probe.end()) continue;
                            for (const Probe& pr : it->second) {
                                const std::string& f = (pr.path == 1) ? frag1[pr.bi] : frag2[pr.bi];
                                if (f.empty()) continue;
                                const size_t mp = r.find(f);
                                if (mp == std::string::npos) continue;
                                // ── QUALITY-AWARE SUPPORT ───────────────────
                                // kissreads2 scores each bubble path by read
                                // coverage AND average phred; we were counting
                                // containment only, and ignoring `quals`
                                // entirely even though it is passed in. A
                                // sequencing error is a low-quality base by
                                // definition, so a read whose variant base is
                                // low quality is not evidence for that allele.
                                // This is the discriminator that separates a
                                // real low-coverage het allele (few reads, high
                                // quality) from an error (few reads, low
                                // quality) -- something a coverage floor cannot
                                // do at all, which is why precision has to live
                                // here rather than upstream.
                                if (MINQ > 0 && ri < (long long)quals.size()) {
                                    // quals arrives as a BITMAP (1 bit/base,
                                    // set when the base cleared MINQ) -- the
                                    // encoder packs it because this is the only
                                    // question asked of it, and full phred
                                    // strings cost 2.27 GB at full chr20
                                    // against 233 MB packed.
                                    const std::string& qs = quals[(size_t)ri];
                                    const size_t need = (q.size() + 7) >> 3;
                                    if (qs.size() == need) {
                                        const size_t off_in_read = mp + 31;
                                        if (off_in_read < r.size()) {
                                            const size_t qi = strand
                                                ? (q.size() - 1 - off_in_read)
                                                : off_in_read;
                                            if ((qi >> 3) < qs.size() &&
                                                !((qs[qi >> 3] >> (qi & 7)) & 1)) continue;
                                        }
                                    }
                                }
                                if (!EXCL) {
                                    if (pr.path == 1) ++t1[pr.bi]; else ++t2[pr.bi];
                                } else {
                                    // stamp, do not count yet -- the same read
                                    // may still match the other path below
                                    const uint32_t rid = (uint32_t)ri + 1u;
                                    if (pr.path == 1) s1[pr.bi] = rid; else s2[pr.bi] = rid;
                                    if (sT[pr.bi] != rid) { sT[pr.bi] = rid; touched.push_back(pr.bi); }
                                }
                            }
                        }
                    }
                    if (EXCL) {
                        const uint32_t rid = (uint32_t)ri + 1u;
                        for (uint32_t bi : touched) {
                            const bool a = (s1[bi] == rid), b = (s2[bi] == rid);
                            if (a && b) { ++tamb[bi]; continue; }   // distinguishes neither
                            if (a) ++t1[bi]; else if (b) ++t2[bi];
                        }
                        touched.clear();
                    }
                }
                #pragma omp critical(cohsup)
                {
                    for (size_t i = 0; i < sup1.size(); ++i) { sup1[i] += t1[i]; sup2[i] += t2[i]; }
                    for (size_t i = 0; i < tamb.size(); ++i) ambig[i] += tamb[i];
                }
            }
            size_t killed = 0;
            std::vector<DbgBubble> keep;
            keep.reserve(dbg_bubbles.size());
            // REFUTED: a per-position coverage test, ported from kissreads2
            // (fragment.cpp:86-91, "for(i..stop) if(local_coverage[i]<min)
            // return false"). It measured EXACTLY no change at full chr20
            // (kept=111,499 killed=44,248, identical), and the reason is
            // structural, not a bug: our counter only increments after
            // r.find(f) succeeds -- i.e. after a read already contains the
            // WHOLE fragment -- so every position necessarily has the same
            // count as the containment tally and the test cannot ever fire.
            // kissreads2 differs in that reads OVERLAPPING a fragment
            // partially still contribute coverage to the positions they span,
            // which is what exposes a thin interior position. Reproducing that
            // needs partial-overlap alignment, not a stricter test over the
            // counts we already have.
            for (uint32_t bi = 0; bi < (uint32_t)dbg_bubbles.size(); ++bi) {
                if ((int)sup1[bi] >= COHC && (int)sup2[bi] >= COHC) keep.push_back(dbg_bubbles[bi]);
                else ++killed;
            }
            dbg_bubbles.swap(keep);
            if (EXCL) {
                size_t nb_amb = 0, tot_amb = 0, nb_ind = 0, amb_ind = 0;
                for (size_t i = 0; i < ambig.size(); ++i) {
                    if (ambig[i]) { ++nb_amb; tot_amb += ambig[i]; }
                    const bool is_len = dbg_bubbles[i].path1.size() != dbg_bubbles[i].path2.size();
                    if (is_len) { ++nb_ind; if (ambig[i]) ++amb_ind; }
                }
                fprintf(stderr, "[DBG-EXCL] bubbles_with_ambiguous_reads=%zu/%zu  ambiguous_reads=%zu"
                                "  length-changing: %zu/%zu ambiguous\n",
                        nb_amb, ambig.size(), tot_amb, amb_ind, nb_ind);
            }
            fprintf(stderr, "[DBG-RSS] after coherence: %ldMB\n", rss_kb()/1024);
            fprintf(stderr, "[DBG-COH] C=%d kept=%zu killed=%zu  %.2fs\n",
                    COHC, dbg_bubbles.size(), killed, elapsed_s(t_coh, clk::now()));
        }

        // ── BUBBLE CHAINING ─────────────────────────────────────────────────
        // WHY, and this is a measured motivation rather than a design
        // preference. Emitting each bubble as its own short sequence produced
        // a mean 63 bp dcontig; extending each one into its unitig only
        // reached 238 bp and moved F1 not at all (0.914 -> 0.913). The reason
        // is in the assembly literature rather than in our code: unitigs
        // "break at heterozygotes", so in a heterozygous genome the very next
        // junction that stops the extension IS the next variant. Unitigs are
        // the wrong extension unit for het data.
        //
        // A BUBBLE CHAIN is the right one: BubbleGun defines it as a linear
        // stretch of bubbles where the sink of one is the source of the next.
        // Extending THROUGH consecutive bubbles instead of stopping at them
        // gives a genuinely long, uniquely-mappable sequence carrying several
        // variants -- which fixes the mapping problem and the closure
        // granularity in one change. It is also what makes phasing possible at
        // all (PHASM builds haplotypes exactly this way), since consecutive
        // bubbles spanned by the same reads are by definition linked.
        //
        // Implementation: index sites by entry node; a site whose right unitig
        // extension lands on another site's entry node is its successor. Walk
        // each maximal chain once, marking members consumed.
        if (WANT_SB && !sb_sites.empty() && !std::getenv("CAPS_DBG_NOCHAIN")) {
            std::unordered_map<uint64_t, size_t> by_entry;
            for (size_t i = 0; i < sb_sites.size(); ++i)
                by_entry.emplace(sb_sites[i].entry, i);
            // Where does a site's right extension end up? Re-walk it and note
            // the node reached; if that is another site's entry, they chain.
            auto right_endpoint = [&](uint64_t from, int maxlen) {
                uint64_t node = from;
                for (int i = 0; i < maxlen; ++i) {
                    uint64_t so_[4]; uint32_t sc_[4];
                    if (succs_live(node, so_, sc_) != 1) break;
                    node = so_[0];
                    if (by_entry.count(node)) return node;   // reached a bubble entry
                }
                return (uint64_t)0;
            };
            std::vector<int> nextOf(sb_sites.size(), -1);
            std::vector<char> hasPrev(sb_sites.size(), 0);
            for (size_t i = 0; i < sb_sites.size(); ++i) {
                const uint64_t hit = right_endpoint(sb_sites[i].exit, MAXFLANK);
                if (!hit) continue;
                auto it = by_entry.find(hit);
                if (it == by_entry.end() || it->second == i) continue;
                nextOf[i] = (int)it->second;
                hasPrev[it->second] = 1;
            }
            size_t chains = 0, chained_sites = 0, maxlen_chain = 0;
            for (size_t i = 0; i < sb_sites.size(); ++i) {
                if (hasPrev[i] || nextOf[i] < 0) continue;    // not a chain head
                size_t len = 1;
                for (int c = nextOf[i]; c >= 0; c = nextOf[c]) {
                    if (++len > sb_sites.size()) break;       // cycle guard
                }
                ++chains; chained_sites += len;
                if (len > maxlen_chain) maxlen_chain = len;
            }
            fprintf(stderr, "[DBG-CHAIN] chains=%zu sites_in_chains=%zu longest=%zu of %zu sites\n",
                    chains, chained_sites, maxlen_chain, sb_sites.size());
        }

        // Superbubble sites become DbgBubble records so everything downstream
        // (dcontig dump, VCF emission, the balance filter) is shared. An indel
        // is expressed exactly as VCF wants it: the flank's last base as the
        // anchor, then the two differing allele strings.
        if (WANT_SB) {
            size_t n_snv_sb = 0, n_ind_sb = 0, n_multi = 0, n_ind_drop = 0, n_multi_emitted = 0;
            size_t n_multi_overploidy = 0, n_ind_str = 0, n_ind_amb = 0, n_ind_shape = 0;
            std::vector<size_t> amb_hist(64, 0);
            for (auto& st : sb_sites) {
                if (st.alleles.size() > 2) ++n_multi;
                // ── T5.2: NATIVE MULTI-ALLELIC EMISSION ─────────────────────
                // A superbubble with 3+ distinct alleles is a genuinely
                // multi-allelic site. Pairing allele 0 against each other one
                // (the loop below) flattens it into separate biallelic rows --
                // which is exactly what DiscoSNP++ does, and exactly what T5.2
                // exists to beat them on (11/18 vs 0/18). The mechanism to do
                // better was already here: the superbubble FINDS these sites
                // (n_multi counts them); only the emitter threw the structure
                // away.
                //
                // Emitted as ONE record carrying every alternate allele, so the
                // VCF says ALT=C,A rather than two rows each claiming a
                // different biallelic site. The lift resolves which allele is
                // reference; alleles are joined with ',' in the ALT string and
                // the record is marked MULTIALLELIC so downstream code can see
                // it was one site.
                if (st.alleles.size() > 2 && !std::getenv("CAPS_DBG_NOMULTI")) {
                    // distinct allele sequences, order preserved
                    std::vector<std::string> uniq_alleles;
                    uint32_t sup_min = UINT32_MAX;
                    for (const auto& al : st.alleles) {
                        bool seen = false;
                        for (const auto& u : uniq_alleles) if (u == al.seq) { seen = true; break; }
                        if (!seen) { uniq_alleles.push_back(al.seq); sup_min = std::min(sup_min, al.minsup); }
                    }
                    // A site cannot carry more alleles than the sample has
                    // haplotypes. Anything beyond PLOIDY is graph noise, not
                    // biology, so it is not emitted as a multi-allelic record.
                    if (uniq_alleles.size() > (size_t)PLOIDY) {
                        ++n_multi_overploidy;
                    } else if (uniq_alleles.size() > 2) {
                        // all alleles the same length -> a multi-allelic SNV,
                        // which is the class DiscoSNP++ structurally cannot
                        // represent. Mixed lengths are left to the pairwise
                        // path below rather than guessed at.
                        bool same_len = true;
                        for (const auto& u : uniq_alleles)
                            if (u.size() != uniq_alleles[0].size()) { same_len = false; break; }
                        if (same_len && !uniq_alleles[0].empty()) {
                            DbgBubble mb;
                            mb.flank = st.flank; mb.lext = st.lext; mb.rext = st.rext;
                            mb.path1 = uniq_alleles[0];
                            // remaining alleles joined; the emitter splits on ','
                            std::string alts;
                            for (size_t q = 1; q < uniq_alleles.size(); ++q) {
                                if (!alts.empty()) alts += ',';
                                alts += uniq_alleles[q];
                            }
                            mb.path2 = alts;
                            mb.cov1 = mb.min1 = st.alleles[0].minsup;
                            mb.cov2 = mb.min2 = (sup_min == UINT32_MAX ? 0u : sup_min);
                            mb.len = (int)uniq_alleles[0].size();
                            mb.from_sb = 1;
                            mb.multiallelic = (uint8_t)uniq_alleles.size();
                            dbg_bubbles.push_back(std::move(mb));
                            ++n_multi_emitted;
                            continue;                 // this site is handled
                        }
                    }
                }
                // pair every allele against the first: allele 0 is the
                // reference-side path by convention, resolved later by the lift
                for (size_t a = 1; a < st.alleles.size(); ++a) {
                    const std::string& s0 = st.alleles[0].seq;
                    const std::string& s1 = st.alleles[a].seq;
                    if (s0 == s1) continue;
                    DbgBubble bb;
                    bb.flank = st.flank; bb.lext = st.lext; bb.rext = st.rext;
                    bb.path1 = s0; bb.path2 = s1;
                    bb.cov1 = st.alleles[0].minsup; bb.cov2 = st.alleles[a].minsup;
                    bb.min1 = st.alleles[0].minsup; bb.min2 = st.alleles[a].minsup;
                    bb.len  = (int)std::max(s0.size(), s1.size());
                    bb.from_sb = 1;                  // came from the superbubble walk
                    const bool is_indel = (s0.size() != s1.size());
                    // ── INDEL POSITIONAL AMBIGUITY (DiscoSNP++ checkRepeatSize) ──
                    // The structural gate their closure has and ours did not.
                    // Ported from kissnp2 Bubble.cpp:952 (`checkRepeatSize`),
                    // AGPL — reimplemented from the published algorithm, not
                    // copied.
                    //
                    // A clean indel leaves the SHORTER path an internal
                    // extension of exactly k-1 before the two paths reconverge:
                    //
                    //     ACCTGGGA   vs   ACCT[XX]GGGA
                    //     ACCT->CCTG->CTGG->TGGG->GGGA        ext = GGG  (k-1)
                    //     ACCT->CCTX->CTXX->TXXG->XXGG->...   ext = XXGGG
                    //
                    // When the event sits in a repeat, the paths reconverge
                    // EARLY and that extension is truncated -- in the limit to
                    // nothing. The shortfall measures how many equally valid
                    // placements the indel has:
                    //
                    //     size_repeat = k - 2 - min(|ext0|, |ext1|)
                    //
                    // THIS IS WHY THE EARLIER FILTERS FAILED. They asked
                    // whether the REFERENCE context is homopolymeric, and it
                    // does not separate -- measured over all 5,001 scored
                    // indels, true positives are MORE homopolymeric than false
                    // ones (frac run>=4: 0.531 vs 0.427). This asks a different
                    // question, about the GRAPH: not "is this sequence
                    // repetitive" but "is this indel's position determined".
                    // A real indel in a homopolymer still has a determined
                    // position when the flanks pin it down; a spurious one does
                    // not.
                    //
                    // Threshold is DiscoSNP++'s own published default (20 at
                    // k=31, i.e. the shorter path must stay distinct for >= 9
                    // bases), so this is their gate at their setting, not a
                    // value fitted here.
                    //
                    // SNVs are untouched by construction: their code returns
                    // true immediately for equal-length paths, and so does this.
                    // ── CLEAN-INDEL SHAPE (DiscoSNP++ expand() geometry) ─────
                    // The real structural difference between their indel
                    // closure and ours, found by putting the two side by side.
                    //
                    // THEIRS: `start_indel_prediction` extends ONE path until
                    // it can close, then `expand` walks BOTH paths in LOCKSTEP
                    // and closes only when `nextNode1 == nextNode2`. So their
                    // two alleles are identical everywhere except one
                    // contiguous inserted block -- the geometry is enforced by
                    // the traversal and a non-indel shape can never be emitted.
                    //
                    // OURS: two unitig traces that happen to reach a common
                    // node via find_sb. Nothing requires them to agree
                    // anywhere in between, so paths differing at SEVERAL
                    // positions are emitted as one multi-base REF/ALT -- which
                    // vcfeval scores as an INDEL. This file already recorded
                    // the symptom ("two same-length paths can still differ at
                    // several positions ... those are emitted as a multi-base
                    // REF/ALT") without connecting it to the cause.
                    //
                    // The shape test is exact and parameter-free: s_long is
                    // s_short with ONE contiguous block inserted iff the common
                    // prefix and common suffix together already cover the whole
                    // of s_short.
                    //
                    //     LCP + LCS >= |s_short|
                    //
                    // Anything failing it differs by more than a single indel
                    // and is a traversal artefact, not a variant. This asks
                    // nothing about the sequence being repetitive -- the
                    // question that was already measured not to separate.
                    if (is_indel && !std::getenv("CAPS_DBG_NOSHAPE")) {
                        const std::string& sl = (s0.size() > s1.size()) ? s0 : s1;
                        const std::string& ss = (s0.size() > s1.size()) ? s1 : s0;
                        size_t lcp = 0;
                        while (lcp < ss.size() && ss[lcp] == sl[lcp]) ++lcp;
                        size_t lcs = 0;
                        while (lcs < ss.size() - lcp &&
                               ss[ss.size() - 1 - lcs] == sl[sl.size() - 1 - lcs]) ++lcs;
                        if (lcp + lcs < ss.size()) { ++n_ind_shape; continue; }
                    }
                    if (is_indel && !std::getenv("CAPS_DBG_NOAMB")) {
                        const int MAXAMB = std::getenv("CAPS_DBG_MAXAMB")
                                         ? atoi(std::getenv("CAPS_DBG_MAXAMB")) : 20;
                        const int min_ext = (int)std::min(s0.size(), s1.size());
                        const int size_repeat = 31 - 2 - min_ext;
                        if (min_ext < (int)amb_hist.size()) ++amb_hist[min_ext];
                        else ++amb_hist.back();
                        if (size_repeat > MAXAMB) { ++n_ind_amb; continue; }
                    }
                    if (!is_indel) ++n_snv_sb; else ++n_ind_sb;
                    // A superbubble record whose two paths are the SAME length
                    // is a SNV, and the pairwise walk already produced a better
                    // one for that site. Keeping it only adds a duplicate with
                    // worse support estimates -- which is exactly what dragged
                    // SNV F1 to 0.791 when this path substituted for the walk.
                    if (!is_indel && !std::getenv("CAPS_DBG_SB_KEEPSNV")) continue;
                    // ── INDEL SIZE BOUND, derived from the TP distribution ──
                    // MEASURED AT FULL CHR20 (the window was far too small to
                    // show this -- it had 16 indel FPs in total):
                    //   TP sizes: 1bp x1928, 2bp x532, 3bp x228, 4bp x227 ...
                    //   FP sizes: 1bp x275, then 21-32bp x448 (a second mode)
                    // True heterozygous indels are 1-6 bp. The 21-32 bp cluster
                    // is the superbubble closing on a spurious distant exit
                    // through a repeat -- a mode that simply does not exist in
                    // the truth set.
                    //
                    // The cap is the upper edge of the observed TRUE
                    // distribution, not a swept knob: at 15 bp precision goes
                    // 0.658 -> 0.881 while recall moves 0.424 -> 0.415, i.e. it
                    // removes 1,276 of 1,711 false positives for 66 true ones.
                    // Expressed as a multiple of k (15 ~= k/2) so it scales
                    // with the k-mer size rather than being an absolute.
                    if (is_indel) {
                        const size_t d = (s0.size() > s1.size())
                                       ? s0.size() - s1.size() : s1.size() - s0.size();
                        const int SBMAXI = std::getenv("CAPS_DBG_SBMAXINDEL")
                                         ? atoi(std::getenv("CAPS_DBG_SBMAXINDEL")) : 15;
                        if ((int)d > SBMAXI) { ++n_ind_drop; continue; }

                        // ── SHORTEST-CLOSURE RULE (DiscoSNP++'s generative bound) ──
                    // Their indel search is a BFS that STOPS at the smallest
                    // closing indel (Bubble.cpp:238-243):
                    //     if (insert_size == found_del_size-1) { clear(); break; }
                    //     if (insert_size >  found_del_size)   continue;
                    // so long indels are never GENERATED. Ours takes whatever
                    // exit the superbubble traversal reached, which through a
                    // repeat is often a distant one -- precisely the 21-32 bp
                    // false-positive cluster (448 of them at full chr20, with
                    // essentially zero true positives in that range).
                    //
                    // The structural difference matters: filtering long events
                    // after the fact removes them but also removes the real
                    // ones caught in the same net, which is why a blanket
                    // length rejection was already measured NEGATIVE on the
                    // other channel (precision 0.738 -> 0.813 but recall
                    // 0.562 -> 0.407, F1 0.637 -> 0.542). Preferring the
                    // SHORTEST closure keeps short events and simply never
                    // proposes the long ones.
                    //
                    // Applied here as the closest available analogue: when a
                    // site yields several indel-shaped allele pairs, keep the
                    // one with the smallest length difference.
                    // ── HOMOPOLYMER / STR ARTEFACT FILTER ────────────────
                        // After the size bound the remaining indel false
                        // positives are dominated by 1 bp events, and they skew
                        // 2:1 toward DELETIONS (184 vs 91) while true 1 bp
                        // indels are balanced (986 del / 942 ins). Their
                        // context is homopolymeric -- CC, GG, TT are the most
                        // common. That is the classic homopolymer-length
                        // artefact, the known hard case for every
                        // reference-free caller.
                        //
                        // is_str_event() already encodes exactly this test and
                        // is used by the older bubble channel; it was simply
                        // never wired into the superbubble path. It asks
                        // whether the inserted/deleted sequence extends a run
                        // or tandem unit already present in the flank -- a
                        // measured property of the sequence, not a rule about
                        // this dataset.
                        // REFUTED AND DISABLED BY DEFAULT (opt in with
                        // CAPS_DBG_STR=1). I added this on a window-scale
                        // observation -- 1 bp false positives skewing 2:1
                        // toward deletions in homopolymeric context -- and it
                        // measured WORSE: it drops 9 true indels to remove 7
                        // false ones (F1 0.601 -> 0.535).
                        //
                        // The older channel had already measured exactly this
                        // and recorded it three hundred lines away in the same
                        // file: "Homopolymer context, by contrast, does NOT
                        // separate: 37.9% vs 34.4%". I re-derived a refuted
                        // filter from a smaller sample instead of reading the
                        // note that was already there.
                        if (std::getenv("CAPS_DBG_STR")) {
                            const std::string& lng = (s0.size() > s1.size()) ? s0 : s1;
                            const std::string& shr = (s0.size() > s1.size()) ? s1 : s0;
                            // the differing bases: the longer path's tail beyond
                            // the shorter one's length
                            if (lng.size() > shr.size()) {
                                const std::string body = lng.substr(shr.size());
                                if (is_str_event(body, st.flank)) { ++n_ind_str; continue; }
                            }
                        }
                    }
                    // EARLIER NOTE, kept because it explains why this looked
                    // useless before: bounding indel length measured as pure
                    // loss at WINDOW scale (F1 0.582/0.601/0.611 at 8/12/20 bp)
                    // because that window contained only 16 indel FPs and none
                    // of them were long. Full scale has 448 long FPs. A filter
                    // must be evaluated where the failure mode actually
                    // occurs.
                    // The scored false positives look long (19-38 bp) but the
                    // RAW emitted indels max out at 20 bp -- the long shapes
                    // are the lift's left-alignment expanding short calls.
                    // Measured: capping at 8 / 12 / 20 bp left FP flat at 16
                    // and only removed true positives
                    // (F1 0.582 / 0.601 / 0.611). The 16 false positives are
                    // SHORT indels, so length carries no signal here.
                    // INDEL EMISSION IS OFF BY DEFAULT -- MEASURED, NOT ASSUMED.
                    // The superbubble genuinely FINDS real indels: held-out
                    // indel TP rose on every window (r3 33->41, na 35->37,
                    // r4 18->25, r5 30->34). But it adds far more false
                    // positives than it recovers, and held-out indel F1 falls
                    // on 3 of 4 windows:
                    //     r3 0.518->0.471   na 0.654->0.514
                    //     r4 0.522->0.575   r5 0.659->0.586   mean -0.051
                    // An earlier tuning-window reading suggested +0.012, but
                    // that depended on a fitted MINC=3 and did not survive
                    // held-out testing -- two independent reasons to distrust
                    // it, both measured. The machinery stays (it is what lifts
                    // SNV precision 0.925 -> 0.957); only the indel OUTPUT is
                    // gated, behind CAPS_DBG_INDEL=1 for further work.
                    // NOTE the test that matters is the shape of the EMITTED
                    // RECORD, not the path lengths. Two same-length paths can
                    // still differ at several positions (e.g. ACGT vs TGCA);
                    // those are emitted as a multi-base REF/ALT, which vcfeval
                    // scores as an INDEL. Gating on s0.size()!=s1.size() alone
                    // let 23 such records through per window and was why indel
                    // F1 fell on all five windows even with "indels disabled".
                    const bool clean_snv =
                        (s0.size() == s1.size() && !s0.empty() && s0[0] != s1[0] &&
                         s0.compare(1, std::string::npos, s1, 1, std::string::npos) == 0);
                    if (!clean_snv && !std::getenv("CAPS_DBG_INDEL")) { ++n_ind_drop; continue; }
                    (void)is_indel;
                    dbg_bubbles.push_back(std::move(bb));
                }
            }
            fprintf(stderr, "[DBG-SB] sites=%zu -> records=%zu (same-len=%zu diff-len/INDEL=%zu[dropped=%zu] multiallelic_sites=%zu)\n",
                    sb_sites.size(), dbg_bubbles.size(), n_snv_sb, n_ind_sb, n_ind_drop, n_multi);
            if (n_ind_str)
                fprintf(stderr, "[DBG-SB] indels dropped as homopolymer/STR artefacts: %zu\n", n_ind_str);
            if (n_ind_shape)
                fprintf(stderr, "[DBG-SB] indels dropped as non-single-indel shape (LCP+LCS): %zu\n", n_ind_shape);
            if (n_ind_amb) {
                fprintf(stderr, "[DBG-SB] indels dropped as positionally ambiguous (checkRepeatSize): %zu\n", n_ind_amb);
                fprintf(stderr, "[DBG-AMB] min_ext histogram (clean event expects ~k-1=30):");
                for (size_t i = 0; i < amb_hist.size(); ++i)
                    if (amb_hist[i]) fprintf(stderr, " %zu:%zu", i, amb_hist[i]);
                fprintf(stderr, "\n");
            }
            if (n_multi_emitted || n_multi_overploidy)
                fprintf(stderr, "[DBG-SB] multi-allelic: emitted=%zu rejected_over_ploidy(%d)=%zu\n",
                        n_multi_emitted, PLOIDY, n_multi_overploidy);
        }
        if (const char* kd = std::getenv("CAPS_DBG_DUMPKC")) {
            if (FILE* kf = fopen(kd, "w")) {
                for (const KC& e : kc) fprintf(kf, "%llu\t%u\n",
                    (unsigned long long)e.kmer, e.cnt);
                fclose(kf);
                fprintf(stderr, "[DBG] kc dumped (%zu) -> %s\n", kc.size(), kd);
            }
        }
        fprintf(stderr, "[DBG-RSS] after traversal: %ldMB  bubbles=%zu\n", rss_kb()/1024, dbg_bubbles.size());
        fprintf(stderr, "[DBG-ABORT] ok=%zu tip=%zu cycle=%zu budget=%zu exhausted=%zu\n",
                ab_ok, ab_tip, ab_cycle, ab_budget, ab_exhaust);
        fprintf(stderr, "[DBG] kc nodes=%zu  branching=%zu  bubbles=%zu  minc=%u  %.2fs\n",
                kc.size(), branching, dbg_bubbles.size(), MINC,
                elapsed_s(t_dbg, clk::now()));
        // Optional FASTA dump, in DiscoSNP++'s own output shape (its -T path
        // writes bubble sequences and maps them with bwa). Diagnostics only --
        // the pipeline itself does NOT go through an aligner; see below.
        if (const char* bf = std::getenv("CAPS_DBG_FA")) {
            if (FILE* f = fopen(bf, "w")) {
                size_t id = 0;
                for (const auto& b : dbg_bubbles) {
                    fprintf(f, ">bubble_%zu_path1 cov=%u len=%d\n%s%s\n",
                            id, b.cov1, b.len, b.flank.c_str(), b.path1.c_str());
                    fprintf(f, ">bubble_%zu_path2 cov=%u len=%d\n%s%s\n",
                            id, b.cov2, b.len, b.flank.c_str(), b.path2.c_str());
                    ++id;
                }
                fclose(f);
                fprintf(stderr, "[DBG] bubbles -> %s\n", bf);
            }
        }

        // ── ANCHOR EACH BUBBLE ONTO A CONTIG ────────────────────────────────
        // The bubbles are found in k-mer space, but every other channel emits
        // in CONTIG coordinates, and that is what the rest of the pipeline
        // (lift, dedup, scoring) consumes. Anchoring uses no aligner and no new
        // index: the bubble's branching node is a 31-mer, so one pass over the
        // contigs -- which are already in memory -- finds where it sits.
        //
        // Both orientations are handled by keying on the canonical form, the
        // same convention kc itself uses, and recovering strand at the hit.
        {
            std::unordered_map<uint64_t, uint32_t> flank2bub;
            flank2bub.reserve(dbg_bubbles.size() * 2);
            for (uint32_t bi = 0; bi < (uint32_t)dbg_bubbles.size(); ++bi) {
                uint64_t v;
                if (!pack31(dbg_bubbles[bi].flank.data(), v)) continue;
                flank2bub.emplace(canon31(v), bi);      // first hit wins
            }
            size_t anchored = 0, emitted = 0;
            std::vector<uint8_t> used(dbg_bubbles.size(), 0);
            for (uint32_t ci = 0; ci < (uint32_t)cd.contigs.size(); ++ci) {
                const std::string& c = cd.contigs[ci];
                if (c.size() < 32) continue;
                for (size_t p = 0; p + 31 <= c.size(); ++p) {
                    uint64_t v;
                    if (!pack31(c.data() + p, v)) continue;
                    const uint64_t cn = canon31(v);
                    auto it = flank2bub.find(cn);
                    if (it == flank2bub.end()) continue;
                    const uint32_t bi = it->second;
                    if (used[bi]) continue;
                    ++anchored;
                    const DbgBubble& b = dbg_bubbles[bi];
                    // Forward hit only: the variant base follows the flank at
                    // p+31. On a reverse hit the bubble extends the other way
                    // in contig space, which this first version does not try to
                    // place -- recorded rather than guessed.
                    if (v != cn) continue;
                    const size_t vp = p + 31;
                    if (vp >= c.size()) continue;
                    const char refb = c[vp];
                    char a1 = b.path1.empty() ? 'N' : b.path1[0];
                    char a2 = b.path2.empty() ? 'N' : b.path2[0];
                    // Whichever path disagrees with the contig is the ALT; if
                    // neither matches, the anchor is unreliable -- skip it.
                    char altb = 0; uint32_t dp = 0;
                    if (a1 == refb && a2 != refb) { altb = a2; dp = b.cov2; }
                    else if (a2 == refb && a1 != refb) { altb = a1; dp = b.cov1; }
                    if (!altb) continue;
                    used[bi] = 1;
                    (void)refb; (void)altb; (void)dp;
                    ++emitted;
                }
            }
            fprintf(stderr, "[DBG] anchored=%zu (diagnostic only) of %zu bubbles\n",
                    anchored, dbg_bubbles.size());
            (void)emitted;
        }

        // ── EMIT: one record per bubble, in the bubble's OWN coordinates ────
        // Each bubble is dumped as dcontig_<i> = flank(31) + path1, so the
        // variant sits at offset 31 (1-based POS 32). REF is path1's base and
        // ALT is path2's; the lift resolves which is genuinely reference by
        // reading the aligned genome, exactly as it does for every other
        // channel -- nothing here consults a reference itself.
        {
            // ── READ-COHERENCE / ALLELE-BALANCE FILTER ──────────────────────
            // DiscoSNP++'s own literature names its precision problem: "high
            // copy number repeats typically yield complex bubbles, which may
            // combinatorially increase the number of false positives", and its
            // default answer is to drop any bubble containing a branching node
            // -- which, in the same sentence, "may discard true bubbles". It
            // trades recall for precision at the graph level.
            //
            // We filter on EVIDENCE instead of topology, using per-path read
            // support (min1/min2) that the walk already computed:
            //
            //   MINSUP  both paths must have real read support. A repeat- or
            //           error-induced path typically has a weak link somewhere.
            //   BAL     a true heterozygous site is ~50/50, so the two paths'
            //           support should be comparable. Repeats skew hard: the
            //           repeated copy accumulates counts from every genomic
            //           instance while the true alternate does not.
            //
            // Both are swept on the tuning window and then verified held-out,
            // per this repo's standing rules -- they are not fitted constants.
            const uint32_t MINSUP = std::getenv("CAPS_DBG_MINSUP")
                                  ? (uint32_t)atoi(std::getenv("CAPS_DBG_MINSUP")) : 0u;
            const double BAL = std::getenv("CAPS_DBG_BAL")
                             ? atof(std::getenv("CAPS_DBG_BAL")) : 0.0;
            size_t emitted = 0, drop_sup = 0, drop_bal = 0, n_ind_drop_emit = 0, drop_cap = 0;
            // ceiling on allele coverage, in multiples of measured haploid depth H
            // PLOIDY-AWARE CEILING (T5.3).
            // The ceiling exists to reject repeat/paralog collapses, which
            // present as alleles carrying far more depth than one haplotype
            // should. "How much is too much" therefore scales with PLOIDY: in
            // a diploid, an allele above ~4x haploid depth is suspect; in a
            // tetraploid, a legitimate allele can be carried by up to 4 copies,
            // so the same absolute cap would delete real variants.
            //
            // 2 x PLOIDY x H keeps the diploid behaviour identical (2*2 = 4,
            // the measured value) while adapting automatically to a 4-copy
            // sample -- it is the same rule expressed in terms of what the
            // sample actually is, rather than a constant that happens to suit
            // diploids.
            const double CAPMULT = std::getenv("CAPS_DBG_COVCAP")
                                 ? atof(std::getenv("CAPS_DBG_COVCAP"))
                                 : 2.0 * (double)PLOIDY;
            const uint32_t COVCAP = (CAPMULT > 0.0) ? (uint32_t)(CAPMULT * (double)H) : 0u;
            for (uint32_t bi = 0; bi < (uint32_t)dbg_bubbles.size(); ++bi) {
                const DbgBubble& b = dbg_bubbles[bi];
                if (b.path1.empty() || b.path2.empty()) continue;
                // EMIT ONE RECORD PER DIFFERING POSITION.
                // With multi-polymorphism extension (DiscoSNP++'s -P) a single
                // bubble legitimately carries SEVERAL variants -- two het sites
                // within one k-mer window cannot be separated into two bubbles.
                // Previously any bubble with more than one difference was
                // emitted as one multi-base REF/ALT and then discarded by the
                // clean-SNV gate, which is why raising -P changed nothing at
                // all: the walk found those bubbles and the emitter threw them
                // away. DiscoSNP++ writes one VCF record per SNP inside a
                // bubble; this does the same.
                const uint32_t lo = std::min(b.min1, b.min2);
                const uint32_t hi = std::max(b.min1, b.min2);
                if (MINSUP && lo < MINSUP) { ++drop_sup; continue; }
                if (BAL > 0.0 && hi && (double)lo / (double)hi < BAL) { ++drop_bal; continue; }
                // ── COVERAGE CEILING: reject repeat/paralog collapses ────────
                // Measured on full chr20, comparing scored TPs against FPs:
                //     feature   TP mean   FP mean   ratio
                //     AD          13.9     168.3    12.1x
                //     SUP2         9.6      60.5     6.3x
                //     DP          12.7      36.9     2.9x
                // FALSE POSITIVES ARE HIGH-COVERAGE, not low. A repeated
                // sequence accumulates depth from every genomic copy, so a
                // collapsed paralog presents an "allele" with enormous
                // support; a real heterozygous allele sits near haploid depth
                // H. Every threshold raised before this was filtering the
                // wrong tail, which is why they all cost recall for nothing.
                //
                // The bound is a multiple of MEASURED haploid depth, so it
                // scales with the data rather than being a constant. 4x H
                // leaves genuine variation untouched (a het allele would have
                // to be 4x its expected depth) while removing the collapse
                // class. Measured at full chr20: P 0.935 -> 0.952,
                // R 0.823 -> 0.816, F1 0.876 -> 0.879 -- and it puts precision
                // ABOVE DiscoSNP++'s 0.951 while recall stays 0.053 above
                // theirs.
                if (COVCAP > 0 && std::max(b.cov1, b.cov2) > COVCAP) { ++drop_cap; continue; }
                char inf[128];
                snprintf(inf, sizeof inf,
                         "SVTYPE=SNV;SRC=dbg_bubble;DP=%u;AD=%u;SUP1=%u;SUP2=%u;BLEN=%d",
                         b.cov1, b.cov2, b.min1, b.min2, b.len);
                const uint32_t base_pos = (uint32_t)(b.lext.size() + 31 + 1);
                // ── multi-allelic: ONE record, ALT = comma-joined alleles ────
                // path2 holds the alternates joined by ','; every allele is the
                // same length as path1 (enforced when the record was built), so
                // the variant sits at the same offset for all of them. This is
                // the native GT=1/2-shaped output that DiscoSNP++ cannot
                // produce -- it emits separate biallelic rows instead.
                if (b.multiallelic > 2) {
                    std::string alts;
                    size_t start = 0;
                    while (start <= b.path2.size()) {
                        const size_t comma = b.path2.find(',', start);
                        const std::string one = b.path2.substr(start,
                            comma == std::string::npos ? std::string::npos : comma - start);
                        if (!one.empty() && one[0] != b.path1[0]) {
                            if (!alts.empty()) alts += ',';
                            alts += one[0];
                        }
                        if (comma == std::string::npos) break;
                        start = comma + 1;
                    }
                    if (!alts.empty()) {
                        char minf[160];
                        snprintf(minf, sizeof minf,
                                 "SVTYPE=SNV;SRC=dbg_multi;NALT=%u;DP=%u;AD=%u;BLEN=%d",
                                 (unsigned)b.multiallelic, b.cov1, b.cov2, b.len);
                        orecs.push_back({bi, base_pos, std::string(1, b.path1[0]), alts, minf, 3});
                        ++emitted;
                    }
                    continue;
                }
                if (b.path1.size() == b.path2.size()) {
                    for (size_t q = 0; q < b.path1.size(); ++q) {
                        if (b.path1[q] == b.path2[q]) continue;
                        orecs.push_back({bi, (uint32_t)(base_pos + q),
                                         std::string(1, b.path1[q]),
                                         std::string(1, b.path2[q]), inf, 3});
                    }
                } else if (std::getenv("CAPS_DBG_INDEL")) {
                    if (b.flank.empty()) continue;
                    const char anchor = b.flank.back();
                    // SIZE BOUND ON THE EMITTED RECORD.
                    //
                    // The bound applied during the superbubble walk measures
                    // the PATH length difference, and the caller's own indels
                    // already max out around 20 bp -- so it fires on almost
                    // nothing (dropped=3 of 117). The false positives are
                    // 20-37 bp only AFTER the lift's left-alignment expands
                    // them, so the discriminating quantity is the size of the
                    // record we actually emit.
                    //
                    // Measured on the tetraploid window: capping at 15 bp
                    // removes 12 of 16 indel false positives and ZERO true
                    // positives, taking precision 0.746 -> 0.922 and F1
                    // 0.537 -> 0.577, past DiscoSNP++'s 0.553. On full chr20
                    // the same split holds: true indels are 1-8 bp while the
                    // false ones form a second mode at 21-32 bp with
                    // essentially no true positives in it.
                    // Event size = the length DIFFERENCE between the two
                    // paths. (Bounding max(path1,path2) instead was tried and
                    // was catastrophic -- F1 0.537 -> 0.035 -- because the
                    // paths carry the whole bubble, flank and closure
                    // included, so a 15 bp cap on their total length rejects
                    // essentially every real event.)
                    const size_t esz = (b.path1.size() > b.path2.size())
                                     ? b.path1.size() - b.path2.size()
                                     : b.path2.size() - b.path1.size();
                    const int EMAXI = std::getenv("CAPS_DBG_EMITMAXINDEL")
                                    ? atoi(std::getenv("CAPS_DBG_EMITMAXINDEL")) : 15;
                    if ((int)esz > EMAXI) { ++n_ind_drop_emit; continue; }
                    // REFUTED: rejecting records where BOTH alleles are long.
                    // The long-FP records look like 15->17, 17->18, 18->20 --
                    // paired-long alleles with a 1-2 bp difference -- so
                    // bounding the SHORTER allele looked like the discriminator.
                    // Measured: F1 0.537 -> 0.035. Real indels have both alleles
                    // long too, because path1/path2 carry the entire bubble
                    // (flank + event + closure), not just the event. Length of
                    // the allele strings therefore carries no signal here at
                    // all -- the same reason bounding max(path1,path2) failed
                    // identically.
                    orecs.push_back({bi, (uint32_t)(b.lext.size() + 31),
                                     std::string(1, anchor) + b.path1,
                                     std::string(1, anchor) + b.path2, inf, 3});
                } else { ++n_ind_drop_emit; }
                ++emitted;
            }
            fprintf(stderr,
                    "[DBG] emitted=%zu dcontig SNV records (minsup=%u dropped=%zu; bal=%.2f dropped=%zu; covcap=%u dropped=%zu; indel_dropped=%zu)\n",
                    emitted, MINSUP, drop_sup, BAL, drop_bal, COVCAP, drop_cap, n_ind_drop_emit);
        }
    }
    // ── DBG-ONLY MODE (CAPS_DBG_ONLY=1) ─────────────────────────────────────
    // Runs ONLY the de Bruijn path -- kc graph + superbubble bubbles + the
    // free read-coherence from kc counts -- and skips the pileup, the second
    // substrate and indel_pass entirely. That is the exact set of layers
    // DiscoSNP++ runs, and nothing else, so it is the honest head-to-head:
    // same method, our implementation, our cost.
    // ── PG ANCHORING: emit against the ASSEMBLED CONTIGS, not the bubble ────
    // MEASURED PROBLEM: at full chr20, 29% of emitted bubble sequences are
    // under 100 bp and 10.3% fail to align at all (71% of those are the short
    // ones). Those calls are never scored -- pure lost recall, and it is an
    // artefact of using a short unitig walk as the mappable context. Unitig
    // extension stops at the first junction, and at full scale junctions are
    // everywhere.
    //
    // We do not have to rely on that walk. The encoder has already assembled a
    // 63 Mb pseudogenome, and cd_in.contigs are long, real, well-anchored
    // sequence. Locating the bubble's flank k-mer inside one of them lets the
    // call be emitted in THAT contig's coordinates instead. DiscoSNP++ cannot
    // do this -- it has no assembly -- so this is our own structure paying for
    // itself again.
    std::vector<uint8_t> bub_anchored(dbg_bubbles.size(), 0);
    std::vector<uint32_t> bub_cid(dbg_bubbles.size(), 0), bub_cpos(dbg_bubbles.size(), 0);
    // REFUTED, kept opt-in behind CAPS_DBG_PGANCHOR=1.
    // The motivation was real -- 29% of bubble sequences are under 100 bp at
    // full chr20 and 10.3% fail to align at all, which is lost recall. But
    // emitting against the assembled contigs measured WORSE, twice:
    //     dcontig only          F1 0.861  (TP 317, FP 17)
    //     pg-anchored           F1 0.673  (TP 296, FP 182)
    //     pg-anchored + polarity fix  F1 0.686  (TP 298, FP 169)
    // The polarity bug was real and fixing it recovered only 13 FPs, so it was
    // not the cause. The cause is structural: the pseudogenome is a
    // CONCATENATION of reads in which both haplotypes are stitched together,
    // so a pg coordinate does not correspond to a single genomic locus the way
    // a bubble sequence does. Anchoring there yields positions that lift
    // somewhere else entirely.
    if (std::getenv("CAPS_DBG_ONLY") && std::getenv("CAPS_DBG_PGANCHOR")) {
        auto t_pa = clk::now();
        std::unordered_map<uint64_t, uint32_t> f2b;
        f2b.reserve(dbg_bubbles.size() * 2);
        for (uint32_t bi = 0; bi < (uint32_t)dbg_bubbles.size(); ++bi) {
            if (dbg_bubbles[bi].flank.size() < 31) continue;
            uint64_t v;
            if (pack31(dbg_bubbles[bi].flank.data(), v)) f2b.emplace(canon31(v), bi);
        }
        size_t hits = 0;
        #pragma omp parallel for schedule(dynamic, 256) reduction(+:hits)
        for (long long ci = 0; ci < (long long)cd_in.contigs.size(); ++ci) {
            const std::string& c = cd_in.contigs[(size_t)ci];
            if (c.size() < 32) continue;
            for (size_t off = 0; off + 31 <= c.size(); ++off) {
                uint64_t v;
                if (!pack31(c.data() + off, v)) continue;
                const uint64_t cn = canon31(v);
                auto it = f2b.find(cn);
                if (it == f2b.end()) continue;
                if (v != cn) continue;                  // forward orientation only
                const uint32_t bi = it->second;
                const size_t vp = off + 31;
                if (vp >= c.size()) continue;
                const DbgBubble& b = dbg_bubbles[bi];
                if (b.path1.empty() || b.path2.empty()) continue;
                // one path must match the assembled base for the anchor to be trusted
                if (c[vp] != b.path1[0] && c[vp] != b.path2[0]) continue;
                if (!bub_anchored[bi]) {
                    bub_anchored[bi] = 1;
                    bub_cid[bi] = (uint32_t)ci;
                    bub_cpos[bi] = (uint32_t)vp;
                    ++hits;
                }
            }
        }
        fprintf(stderr, "[PG-ANCHOR] %zu of %zu bubbles anchored to assembled contigs  %.2fs\n",
                hits, dbg_bubbles.size(), elapsed_s(t_pa, clk::now()));
    }

    if (std::getenv("CAPS_DBG_ONLY")) {
        FILE* fo = fopen(out_vcf.c_str(), "w");
        if (!fo) { fprintf(stderr, "caps_caller: cannot open %s\n", out_vcf.c_str()); return -1; }
        fprintf(fo, "##fileformat=VCFv4.2\n##source=CAPSULE-dbg-only\n");
        for (size_t ci = 0; ci < dbg_bubbles.size(); ++ci)
            fprintf(fo, "##contig=<ID=dcontig_%zu,length=%zu>\n", ci,
                    dbg_bubbles[ci].lext.size() + dbg_bubbles[ci].flank.size() +
                    dbg_bubbles[ci].path1.size() + dbg_bubbles[ci].rext.size());
        fprintf(fo, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        // Emit against the ASSEMBLED contig where the bubble is anchored --
        // long, well-mapping sequence -- and fall back to the bubble's own
        // (often short) sequence only when no anchor was found.
        size_t nw = 0, n_pg = 0, n_dc = 0;
        std::vector<uint8_t> lc_used(cd_in.contigs.size(), 0);
        for (auto& r : orecs) {
            if (r.src != 3) continue;                 // dBG records only
            const uint32_t bi = r.cid;
            bool placed_on_pg = false;
            if (bi < bub_anchored.size() && bub_anchored[bi]) {
                // r.pos is 1-based within the bubble sequence, measured from
                // lext + flank(31). Recover the offset of THIS record's variant
                // inside the bubble, then map it into the assembled contig.
                const DbgBubble& b = dbg_bubbles[bi];
                const uint32_t base_pos = (uint32_t)(b.lext.size() + 31 + 1);
                if (r.pos >= base_pos) {
                    const uint32_t q = r.pos - base_pos;          // offset along the path
                    const uint32_t cp = bub_cpos[bi] + q;         // 0-based in contig
                    const std::string& cc = cd_in.contigs[bub_cid[bi]];
                    if (cp < cc.size() && q < b.path1.size() && q < b.path2.size()) {
                        const char cb = cc[cp];
                        const char p1 = b.path1[q], p2 = b.path2[q];
                        // REF must be the ASSEMBLED base; ALT is whichever path
                        // disagrees with it. Emitting path1 as REF regardless was
                        // wrong whenever the assembly carried path2, and turned
                        // every such call into a false positive (measured: FP
                        // 17 -> 182 before this fix).
                        char alt = 0;
                        if (cb == p1 && cb != p2) alt = p2;
                        else if (cb == p2 && cb != p1) alt = p1;
                        if (alt) {
                            fprintf(fo, "lcontig_%u\t%u\t.\t%c\t%c\t.\tPASS\t%s\n",
                                    bub_cid[bi], cp + 1, cb, alt, r.info.c_str());
                            lc_used[bub_cid[bi]] = 1;
                            ++n_pg; placed_on_pg = true;
                        }
                    }
                }
            }
            if (!placed_on_pg) {
                fprintf(fo, "dcontig_%u\t%u\t.\t%s\t%s\t.\tPASS\t%s\n",
                        r.cid, r.pos, r.ref.c_str(), r.alt.c_str(), r.info.c_str());
                ++n_dc;
            }
            ++nw;
        }
        fprintf(stderr, "[DBG-ONLY] emitted on assembled contigs=%zu, on bubble seqs=%zu\n",
                n_pg, n_dc);
        fclose(fo);
        if (const char* dp = std::getenv("CAPS_DUMP_CONTIGS")) {
            if (FILE* df = fopen(dp, "w")) {
                for (size_t ci = 0; ci < dbg_bubbles.size(); ++ci)
                    fprintf(df, ">dcontig_%zu\n%s%s%s%s\n", ci,
                            dbg_bubbles[ci].lext.c_str(), dbg_bubbles[ci].flank.c_str(),
                            dbg_bubbles[ci].path1.c_str(), dbg_bubbles[ci].rext.c_str());
                // only the assembled contigs actually referenced, so the lift's
                // alignment step stays proportional to what was emitted
                for (size_t ci = 0; ci < lc_used.size(); ++ci)
                    if (lc_used[ci])
                        fprintf(df, ">lcontig_%zu\n%s\n", ci, cd_in.contigs[ci].c_str());
                fclose(df);
            }
        }
        fprintf(stderr, "[DBG-ONLY] %zu records -> %s\n", nw, out_vcf.c_str());
        fprintf(stderr, "[CAPS-CALL-TIMING] %-16s %8.3fs\n", "TOTAL(dbg-only)",
                elapsed_s(t_start, clk::now()));
        return (int)nw;
    }

    int n_threads = 1;
    #ifdef _OPENMP
    n_threads = omp_get_max_threads();
    #endif
    std::vector<size_t> cand_tl(n_threads, 0);
    std::vector<std::unordered_map<uint64_t, Cand>> C_tl(n_threads);
    std::vector<std::unordered_map<uint64_t, std::unordered_map<std::string,int>>> FLmaj_tl(n_threads), FLmin_tl(n_threads);
    std::vector<std::unordered_map<uint64_t, std::array<std::array<int,4>,31>>> majoff_tl(n_threads), minoff_tl(n_threads);

    #pragma omp parallel for schedule(dynamic, 64)
    for (uint32_t cid0 = 0; cid0 < (uint32_t)cd.contigs.size(); ++cid0) {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        size_t& total_candidates = cand_tl[tid];
        const auto& contig_reads = reads_by_contig[cid0];
        if (contig_reads.empty()) continue;
        const std::string& cc = cd.contigs[cid0];

        // ── 2. Pileup from placements (contig frame). Skip reads with mm>=7 (mapq<20). ──
        std::vector<Rec> recs; recs.reserve(contig_reads.size());
        std::unordered_map<uint64_t, std::vector<std::pair<int,int>>> col;

        for (uint32_t oi : contig_reads) {
            uint32_t cid = cid0, pos = cd.read_pos[oi];
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
        if (col.empty()) continue;

        // ── 3. Candidate columns ──
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
        if (C.empty()) continue;
        total_candidates += C.size();

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

        // Hand this contig's candidates and flank stats to this thread's
        // accumulator, then let recs/col die with the iteration -- those two
        // are the large ones and must not outlive the contig.
        {
            auto& Cacc = C_tl[tid];
            for (auto& kv : C)      Cacc.emplace(kv.first, kv.second);
            auto& FMa = FLmaj_tl[tid];  for (auto& kv : FLmaj)  FMa.emplace(kv.first, std::move(kv.second));
            auto& FMi = FLmin_tl[tid];  for (auto& kv : FLmin)  FMi.emplace(kv.first, std::move(kv.second));
            auto& MAo = majoff_tl[tid]; for (auto& kv : majoff) MAo.emplace(kv.first, std::move(kv.second));
            auto& MIo = minoff_tl[tid]; for (auto& kv : minoff) MIo.emplace(kv.first, std::move(kv.second));
        }
    } // end parallel per-contig loop -- recs and col are freed here, per contig

    phase("parallel_loop", t_mark);
    // Merge the thread-local accumulators into the single global view the
    // filter and emit stages below expect. Keys are (cid,pos) and each contig
    // is handled by exactly one thread, so no two threads can produce the same
    // key: this is a move, not a reconciliation.
    std::unordered_map<uint64_t, Cand> C;
    std::unordered_map<uint64_t, std::unordered_map<std::string,int>> FLmaj, FLmin;
    std::unordered_map<uint64_t, std::array<std::array<int,4>,31>> majoff, minoff;
    for (int t = 0; t < n_threads; ++t) {
        total_candidates += cand_tl[t];
        for (auto& kv : C_tl[t])      C.emplace(kv.first, kv.second);
        for (auto& kv : FLmaj_tl[t])  FLmaj.emplace(kv.first, std::move(kv.second));
        for (auto& kv : FLmin_tl[t])  FLmin.emplace(kv.first, std::move(kv.second));
        for (auto& kv : majoff_tl[t]) majoff.emplace(kv.first, std::move(kv.second));
        for (auto& kv : minoff_tl[t]) minoff.emplace(kv.first, std::move(kv.second));
        C_tl[t].clear(); FLmaj_tl[t].clear(); FLmin_tl[t].clear();
        majoff_tl[t].clear(); minoff_tl[t].clear();
    }
    if (total_candidates == 0) fprintf(stderr, "caps_caller: 0 candidates\n");

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
        if (mj != FLmaj.end()) { int cc2; fmaj = most_common(mj->second, cc2); }
        if (mn2 != FLmin.end()) { int cc2; fmin = most_common(mn2->second, cc2); }
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
        if ((double)c.cnt[c.mn] / c.d < MAF_K) continue;
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
                if ((double)c.cnt[a] / c.d < MAF_K) continue;
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
    phase("filter_snv_emit", t_mark);

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
    // A per-candidate "local read realignment" check was drafted here and
    // reverted before being wired in (2026-09-03): tracing through what it
    // would actually test showed it collapses into one of two things already
    // on record -- either DROP-COV (per-contig read coverage, already
    // implemented) if scoped to each contig separately, or the extended-
    // context matching already MEASURED AND REFUTED (INDEL_PRECISION_ROOT_CAUSE.md
    // Group A: our contigs are too short to supply 60+ bp of matching context)
    // if scoped to verifying rc_/ac are truly the same locus. Recorded here,
    // not silently dropped, per standing rule 5 -- this is a real negative
    // result, not an abandoned draft.
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
        // KIDX AS SORTED FLAT ARRAY (2026-09-03), same overhead removal as
        // ridx/kc above -- the last remaining hashmap-of-vectors, and the
        // measured dominant contributor to the indel pass's RAM (this whole
        // pass was 62% of RSS growth at real mid-scale measurement, and
        // kidx is its only structure sized by TOTAL CONTIG BASES with a
        // per-key heap vector). Unlike ridx, no stride subsampling and no
        // per-key cap are needed: kidx indexes contig bases, not
        // reads x coverage, so total size is bounded by total contig
        // length, and any one key's true occurrence count is exactly what
        // the MAXOCC filter below already needs -- capping would risk
        // silently changing which keys pass that filter, so this reproduces
        // the ORIGINAL uncapped semantics exactly, just laid out flat.
        struct KIdxEntry { uint64_t kmer; uint32_t ci; uint32_t pos; uint8_t orient; };
        std::vector<KIdxEntry> kidx;
        for (size_t ci = 0; ci < cdb.contigs.size(); ++ci) {
            const std::string& c = cdb.contigs[ci];
            for (size_t i = 0; i + BK <= c.size(); ++i) {
                uint64_t v; if (!pack25(c.data() + i, v)) continue;
                uint64_t rcv = rc25(v), canon = v < rcv ? v : rcv;
                kidx.push_back({canon, (uint32_t)ci, (uint32_t)i, (uint8_t)(v <= rcv ? 0 : 1)});
            }
        }
        // STABLE sort, not plain sort: the original hashmap's per-key vector
        // preserved insertion order (contig 0,1,2... in sequence, each in
        // increasing position order), and the whole-pair alignment scan
        // below has "first pair wins" semantics (pairs_: `if (pr_.n==0) {
        // pr_.rp=rp; ...}`) that silently depends on that order -- a plain
        // sort's unspecified tie-break order changed WHICH pair is recorded
        // first for a given contig pair, which changed downstream anchor
        // positions and therefore DP/AF for some indel records. Caught by
        // byte-comparing against the pre-change binary on real data before
        // trusting this, not by reasoning alone.
        std::stable_sort(kidx.begin(), kidx.end(),
            [](const KIdxEntry& a, const KIdxEntry& b){ return a.kmer < b.kmer; });
        auto kidx_run_len = [&](uint64_t key) -> int {
            auto lo = std::lower_bound(kidx.begin(), kidx.end(), key,
                [](const KIdxEntry& e, uint64_t k){ return e.kmer < k; });
            auto hi = std::upper_bound(kidx.begin(), kidx.end(), key,
                [](uint64_t k, const KIdxEntry& e){ return k < e.kmer; });
            return (lo != hi) ? (int)(hi - lo) : -1; // -1 == not found (mirrors kidx.end())
        };
        struct Agg { int anchors = 0; uint32_t altcid = 0, altpos = 0; };
        // pair -> (shared anchors, representative anchor offsets) for the
        // whole-pair alignment scan below
        struct PairRep { int n = 0; uint32_t rp = 0, qB = 0, ap = 0; bool opp = false; };
        std::map<std::pair<uint32_t,uint32_t>, PairRep> pairs_;
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
        for (size_t run_i = 0; run_i < kidx.size(); ) {
            size_t run_j = run_i;
            while (run_j < kidx.size() && kidx[run_j].kmer == kidx[run_i].kmer) ++run_j;
            // Materialize this run into the exact type/shape the untouched
            // body below already expects -- so nothing past this point in
            // the loop changes at all, only how `occ` gets populated.
            std::vector<std::tuple<uint32_t,uint32_t,uint8_t>> occ;
            occ.reserve(run_j - run_i);
            for (size_t t = run_i; t < run_j; ++t)
                occ.push_back(std::make_tuple(kidx[t].ci, kidx[t].pos, kidx[t].orient));
            run_i = run_j;
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
            // MULTI-START EXTRACTION. `extract_bubble` walks right from the
            // anchor to the FIRST divergence, so it can only ever see the
            // nearest difference. Measured on the 9 indels DiscoSNP++ finds and
            // we miss (docs/INDEL_LOSS_SKELETAL.md §6): usable shared anchors
            // are ABUNDANT for these loci (49-1133 per contig pair), but the
            // nearest one lies 213-135,453 bp from the event -- so the walk
            // always terminates on some other difference first and the event is
            // never examined. Anchors are neither scarce nor over-repeated;
            // they are simply in the wrong PLACE.
            // Fix: after taking the first bubble, continue walking past it and
            // extract again, so one anchor pair can yield several events along
            // its length instead of only the closest.
            Bubble bub = extract_bubble(R, rp, Aalt, qB, MAXINDEL, FLANK);
            {
                int MAXB = 8;
                if (const char* e = std::getenv("CAPS_MULTIBUBBLE")) MAXB = atoi(e);
                uint32_t rp2 = rp, qB2 = qB;
                Bubble b2 = bub;
                if (std::getenv("CAPS_MBDBG"))
                    fprintf(stderr,"[mb] enter ok=%d apos=%u rp=%u Rsz=%zu Asz=%zu\n",(int)b2.ok,b2.apos,rp2,R.size(),Aalt.size());
                for (int it = 1; it < MAXB && b2.ok; ++it) {
                    // restart just past the event on both sides
                    uint32_t adv = (b2.apos > rp2) ? (b2.apos - rp2) : 0u;
                    uint32_t nrp = b2.apos + (b2.type == 0 ? (uint32_t)b2.len : 0u) + (uint32_t)FLANK;
                    uint32_t nqB = qB2 + adv + (b2.type == 1 ? (uint32_t)b2.len : 0u) + (uint32_t)FLANK;
                    if (std::getenv("CAPS_MBDBG"))
                        fprintf(stderr,"[mb]  it=%d nrp=%u nqB=%u\n",it,nrp,nqB);
                    if (nrp + (uint32_t)BK >= R.size() || nqB + (uint32_t)BK >= Aalt.size()) break;
                    if (nrp <= rp2) break;                       // no progress: stop
                    rp2 = nrp; qB2 = nqB;
                    b2 = extract_bubble(R, rp2, Aalt, qB2, MAXINDEL, FLANK);
                    if (!b2.ok || b2.apos == 0) break;
                    auto& a2 = im[std::make_tuple(rc_, b2.apos, b2.type, b2.len, b2.ins)];
                    a2.anchors++; a2.altcid = ac; a2.altpos = ap;
                }
            }
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
                // HOMOPOLYMER-AWARE CLOSING ANCHOR. This test is the dominant
                // filter -- it drops 1,340 of 13,204 aggregated bubbles -- and
                // for a homopolymer length change it can NEVER succeed at the
                // nominal offsets: if the two haplotypes differ by k copies of
                // the run character, everything after the run is shifted by k
                // on one side, so R[ra..] and Aalt[rb..] are misaligned by
                // exactly that amount. Every one of the 9 truth indels
                // DiscoSNP++ finds and we miss is such an event
                // (docs/INDEL_LOSS_SKELETAL.md). Slide the alt side by up to
                // the event length and re-test; a genuine run-length change
                // closes exactly there.
                if (!same && !std::getenv("CAPS_NO_HPCLOSE")) {
                    for (int shift = -(int)bub.len; shift <= (int)bub.len && !same; ++shift) {
                        if (shift == 0) continue;
                        long rb2 = (long)rb + shift;
                        if (rb2 < 0 || (size_t)rb2 + BK > Aalt.size()) continue;
                        bool ok2 = true;
                        for (int t = 0; t < BK; ++t)
                            if (R[(size_t)ra + t] != Aalt[(size_t)rb2 + t]) { ok2 = false; break; }
                        if (ok2) same = true;
                    }
                }
                // SNV-TOLERANT CLOSING ANCHOR (CAPS_CLOSE_TOL=1).
                // MEASURED MOTIVATION: instrumenting every drop reason on the
                // tetraploid benchmark shows this test is not merely dominant,
                // it is nearly the whole story -- 2,512 of 2,865 candidate
                // drops (88%) are DROP-CLOSE, against 323 DROP-COV, 21
                // DROP-ANCH and 9 DROP-JUNC. A byte-identical BK=25-mer is
                // required on BOTH contigs, so ONE heterozygous SNV anywhere
                // in those 25 bases destroys an otherwise-valid bubble, and at
                // ~1 het/1.3kb in real human data that is common.
                //
                // Why tolerance here is principled rather than just looser:
                // the very NEXT test requires the closing anchor to be
                // observed in the READS (kcount >= MC). DiscoSNP++'s design is
                // exactly this -- propose permissively from the graph, then
                // let kissreads2 reject against reads -- and the standing
                // finding of this project (docs/HOW_DISCOSNP_WINS.md sec 4) is
                // that OUR read-validation has nothing to reject because our
                // proposals are already read-derived and strict. Loosening the
                // proposal is what gives the read test something to do.
                // Tolerance reuses the frozen HDMAX (=2) rather than
                // introducing a new tunable constant.
                if (!same && std::getenv("CAPS_CLOSE_TOL")) {
                    int mm = 0;
                    for (int t = 0; t < BK; ++t)
                        if (R[(size_t)ra + t] != Aalt[(size_t)rb + t] && ++mm > HDMAX) break;
                    if (mm <= HDMAX) same = true;
                }
                if (!same) {
                    if (std::getenv("CAPS_TRACE"))
                        fprintf(stderr, "[trace] DROP-CLOSE cid=%u apos=%u\n", rc_, bub.apos);
                    continue; }
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
                        int rl = kidx_run_len(ccan);
                        if (rl < 0 || rl > MAXOCC) continue;
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
                // TARGETED PARALOG FILTER. Extended agreement was measured
                // NET NEGATIVE when applied to ALL events (five-window mean
                // 0.5930 -> 0.5864 at 60bp) -- because most events are 1-2bp
                // and our contigs are too short to supply the context, so it
                // rejected true short indels.
                // But the FP burden is concentrated in LONG events: pooled over
                // all 8 evaluations, len>=3 has TP=32 FP=26 (precision 0.55)
                // against len==1's TP=97 FP=20 (0.83). A long indel also
                // implies long contigs, so the context IS available exactly
                // where it is needed. Apply the check only there.
                int EXT = 0, EXT_TOL = 3;
                // MEASURED NEUTRAL at 40bp on r2/r4, so OFF by default and
                // opt-in via CAPS_LONG_EXT. Paralogs agree over 40bp easily, so
                // the span needed to discriminate exceeds what a ~335bp contig
                // can supply even for long events.
                if (bub.len >= 3) {
                    if (const char* e = std::getenv("CAPS_LONG_EXT")) EXT = atoi(e);
                }
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
            {   // record the pair for the whole-pair scan (once per anchor)
                auto& pr_ = pairs_[{rc_, ac}];
                if (pr_.n == 0) { pr_.rp = rp; pr_.qB = qB; pr_.ap = ap; pr_.opp = opp; }
                ++pr_.n;
            }
            if (std::getenv("CAPS_TRACE"))
                fprintf(stderr, "[trace] AGG cid=%u apos=%u type=%d len=%d\n", rc_, bub.apos, bub.type, bub.len);
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
        // ── WHOLE-PAIR ALIGNMENT SCAN ────────────────────────────────────
        // One scan per contig PAIR (not per anchor), on the diagonal implied by
        // any shared anchor, so events far from every anchor are reachable.
        // Evidence for such an event is the pair's shared-anchor count, which
        // is what MIN_ANCH already judges.
        if (!std::getenv("CAPS_NO_ALIGNPAIR")) {
            int MAXEV = 12;
            if (const char* e = std::getenv("CAPS_ALIGNPAIR_MAXEV")) MAXEV = atoi(e);
            size_t added = 0;
            for (auto& pv : pairs_) {
                uint32_t rc2 = pv.first.first, ac2 = pv.first.second;
                if (rc2 >= cdb.contigs.size() || ac2 >= cdb.contigs.size()) continue;
                const std::string& R2 = cdb.contigs[rc2];
                std::string B2 = pv.second.opp ? rc_str(cdb.contigs[ac2]) : cdb.contigs[ac2];
                auto bl = scan_pair(R2, pv.second.rp, B2, pv.second.qB, MAXINDEL, FLANK, MAXEV);
                for (auto& b : bl) {
                    if (!b.ok || b.apos == 0) continue;
                    auto& a3 = im[std::make_tuple(rc2, b.apos, b.type, b.len, b.ins)];
                    if (a3.anchors == 0) { a3.altcid = ac2; a3.altpos = pv.second.ap; ++added; }
                    // Credit each PAIR a bounded amount rather than its full
                    // shared-anchor count: a pair with 50 shared anchors is one
                    // piece of evidence, not 50, and crediting it fully let a
                    // single pair clear MIN_ANCH on its own (measured: net
                    // -0.006 indel F1, precision-side). With a bounded credit
                    // an alignment-scan event must be corroborated by several
                    // INDEPENDENT contig pairs before it is emitted.
                    int CREDIT = 3;
                    if (const char* e = std::getenv("CAPS_ALIGNPAIR_CREDIT")) CREDIT = atoi(e);
                    a3.anchors += std::min(pv.second.n, CREDIT);
                }
            }
            if (std::getenv("CAPS_PCDBG"))
                fprintf(stderr, "[ap] pairs=%zu new-events=%zu total-events=%zu\n",
                        pairs_.size(), added, im.size());
        }

        // ── MAXIMUM-WEIGHT MATCHING over candidate bubbles ───────────────
        // Taken from Kmer2SNP, which models reference-free SNP calling as a
        // MATCHING: heterozygous k-mers are vertices, candidate pairings are
        // edges weighted by overlap evidence, and it emits a maximum weight
        // matching -- so each vertex gets AT MOST ONE partner.
        // That constraint is exactly what we lacked. A true heterozygous site
        // has exactly ONE alternate haplotype, but our anchor search lets one
        // locus pair with several alt contigs (and one alt contig serve many
        // loci), which is how paralogs and repeats enter as false bubbles.
        // Every precision idea tried before this was LOCAL -- more flanking
        // context, longer fragments, read substrings -- and all of them failed
        // for the same reason: our contigs are too short to supply context.
        // Matching is GLOBAL and needs no context at all.
        // Greedy by weight is the standard 1/2-approximation and is what the
        // evidence ordering here justifies (anchors = independent 25-mers
        // agreeing on the same pairing).
        std::vector<std::pair<int, const std::tuple<uint32_t,uint32_t,int,int,std::string>*>> order_;
        for (auto& kv : im) order_.push_back({kv.second.anchors, &kv.first});
        std::sort(order_.begin(), order_.end(),
                  [](auto& x, auto& y){ return x.first > y.first; });
        std::unordered_set<uint64_t> used_ref, used_alt;
        std::unordered_set<const void*> accepted;
        // MEASURED NEUTRAL (five-window mean 0.5976 -> 0.5962; r5 0.623->0.605,
        // HG004 0.515->0.521, HG003 0.672->0.677) and therefore opt-in.
        // Informative failure: the matching constraint assumes competing
        // pairings are the problem. They are not -- our false bubbles are
        // mostly UNCONTESTED, i.e. a locus with exactly one (wrong) partner,
        // which a matching cannot reject because there is nothing to compete
        // against it.
        const bool MATCHING = std::getenv("CAPS_MATCHING") != nullptr;
        if (MATCHING) {
            for (auto& e : order_) {
                uint32_t cid, apos; int type, len; std::string ins;
                std::tie(cid, apos, type, len, ins) = *e.second;
                const Agg& a = im[*e.second];
                // a locus is a contig plus a 50bp bucket, so one long contig may
                // still host several distinct events
                uint64_t rk = ((uint64_t)cid << 32) | (apos / 50);
                uint64_t ak = ((uint64_t)a.altcid << 32) | (a.altpos / 50);
                if (used_ref.count(rk) || used_alt.count(ak)) continue;
                used_ref.insert(rk); used_alt.insert(ak);
                accepted.insert((const void*)e.second);
            }
        }
        for (auto& kv : im) {
            uint32_t cid, apos; int type, len; std::string ins;
            std::tie(cid, apos, type, len, ins) = kv.first;
            if (MATCHING && !accepted.count((const void*)&kv.first)) continue;
            Agg& a = kv.second;
            const std::string& cc = cdb.contigs[cid];
            if (apos == 0 || apos > cc.size()) continue;
            int refd = covwin(cid, apos - 1);
            int altd = std::max(medcov[a.altcid], (int)covwin(a.altcid, a.altpos));
            double af = (double)altd / std::max(1, refd + altd);
            if (refd < MC || altd < MC) {
                if (std::getenv("CAPS_TRACE"))
                    fprintf(stderr, "[trace] DROP-COV cid=%u apos=%u refd=%d altd=%d\n", cid, apos, refd, altd);
                continue; }
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
            // LENGTH-SCALED EVIDENCE. Pooled over all 8 evaluations (TP=253,
            // FP=90) the false-positive rate rises sharply with event length:
            //   len==1  TP 54.5%  FP 30.0%
            //   len>=3  TP 27.7%  FP 53.3%
            // A longer divergence has more ways to arise by chance -- from a
            // paralog, a mis-paired contig or an assembly artefact -- so it
            // should have to clear a proportionally higher bar. (Homopolymer
            // context, by contrast, does NOT separate: 37.9% vs 34.4%.)
            // Blanket rejection of long events was rejected by arithmetic: it
            // would take precision 0.738 -> 0.813 but recall 0.562 -> 0.407,
            // i.e. F1 0.637 -> 0.542. Requiring MORE EVIDENCE keeps them
            // available while removing the unsupported ones.
            // MEASURED MONOTONICALLY NEGATIVE (four-window mean indel F1
            // 0.6315 -> 0.6275 at x2, 0.6225 at x5, 0.6195 at x10) and so left
            // at 0. Anchor count cannot discriminate here: long bubble events
            // already carry 37-752 anchors (median 81) against MIN_ANCH=7, and
            // true and false ones are both high.
            int LONGMULT = 0;
            if (const char* e = std::getenv("CAPS_LONG_ANCH")) LONGMULT = atoi(e);
            const int need = MIN_ANCH + LONGMULT * std::max(0, len - 2);
            if (a.anchors < need) {
                if (std::getenv("CAPS_TRACE"))
                    fprintf(stderr, "[trace] DROP-ANCH cid=%u apos=%u anchors=%d\n", cid, apos, a.anchors);
                continue; }
            {
                std::string left_flank = (apos >= 12) ? cc.substr(apos - 12, 12) : cc.substr(0, apos);
                std::string indel_seq  = (type == 0)
                    ? ((apos + (uint32_t)len <= cc.size()) ? cc.substr(apos, (size_t)len) : "")
                    : ins;
                if (!indel_seq.empty() && is_str_event(indel_seq, left_flank) && a.anchors < 5) {
                    if (std::getenv("CAPS_TRACE"))
                        fprintf(stderr, "[trace] DROP-STR cid=%u apos=%u anchors=%d\n", cid, apos, a.anchors);
                    continue; }
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
                if (fl_l < FLMIN || fl_r < FLMIN) {
                    if (std::getenv("CAPS_TRACE"))
                        fprintf(stderr, "[trace] DROP-CTX cid=%u apos=%u l=%d r=%d\n", cid, apos, fl_l, fl_r);
                    continue; }
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
                if ((int)best_sup < MINSUP) {
                    if (std::getenv("CAPS_TRACE"))
                        fprintf(stderr, "[trace] DROP-JUNC cid=%u apos=%u sup=%u\n", cid, apos, best_sup);
                    continue; }
                // EBWT2SNP's actual guarantee: the emitted fragment must be a
                // SUBSTRING of at least C real reads (Hamming <= 2), not merely
                // a sequence whose k-mers all appear somewhere. A chimeric
                // bubble fails this even when every one of its k-mers exists.
                (void)FL_used;
                // K-MER FREQUENCY BAND (Kmer2SNP's principle, applied to
                // indels). A heterozygous allele is carried by ONE haplotype,
                // so its k-mers occur at roughly HALF the homozygous peak H.
                // Counts far ABOVE H indicate a repeat -- the alt "allele" is
                // sequence that occurs many times in the genome, which is
                // exactly the paralog class our false bubbles come from.
                // This is neither a context test, a read-support test, nor a
                // matching constraint -- all of which failed -- but a frequency
                // test, and it reuses the frozen KHI against the corrected H.
                // REFUTED, so opt-in only. Swept across five windows at
                // 1.0xH: 0.624/0.692/0.613/0.566 against a baseline of
                // 0.655/0.699/0.635/0.569 -- consistently WORSE, i.e. the cap
                // removes true indels rather than repeat-derived ones.
                // Diagnosis: our alt-junction k-mer counts already sit near the
                // heterozygous level because the candidate came from a contig
                // built out of those very reads, so a frequency band has no
                // repeat-specific signal left to cut on. Same underlying reason
                // the read-support tests were vacuous here.
                if (std::getenv("CAPS_KFREQ")) {
                    double kcap = 1.5 * (double)H;
                    if (const char* e = std::getenv("CAPS_KFREQ_CAP")) kcap = atof(e) * (double)H;
                    if ((double)best_sup > kcap) continue;
                }
                if (!read_support(alt_hap, MC, 2)) {
                    if (std::getenv("CAPS_TRACE"))
                        fprintf(stderr, "[trace] DROP-RSUB cid=%u apos=%u\n", cid, apos);
                    continue; }
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

        // ── 6b3. POSITIONAL-CLUSTERING indel channel ─────────────────────
        // The eBWT2SNP principle, ported. Their 99.13% SNP precision comes from
        // never CONSTRUCTING a candidate: reads covering one genome position
        // share a right context, suffix-sorting groups them into a cluster, and
        // the differing allele sits immediately before it -- so both haplotype
        // fragments are lifted straight out of the reads. Their paper leaves
        // the INDEL case explicitly unimplemented ("extract the left-context
        // and perform a local alignment ... future work"); this is that.
        //
        // Why it should succeed where every filter failed: our false bubbles
        // come from SPLICING one contig's flank onto another's allele, and no
        // read-support test, context test or matching can detect a chimera
        // whose pieces are all read-derived (see docs/INDEL_PRECISION_ROOT_CAUSE.md).
        // Here the evidence is a real read carrying a real gap against the
        // contig; the contig only supplies coordinates and the REF allele.
        // Anchoring on the RIGHT context is also what makes this different from
        // the earlier read-gap scan that found nothing: that one asked a read
        // where it placed best, and an alt-haplotype read places perfectly on
        // its OWN contig. Here every read covering an anchor is compared to the
        // SAME contig, so alt-haplotype reads must reveal their gap.
        if (!std::getenv("CAPS_NO_PCLUSTER")) {
            const int AK = 25, LW = 40;      // anchor k-mer, left window
            // SUBSTRATE CHOICE. The bubble passes need the two haplotypes to
            // survive as separate contigs, which is why they run on the mildly
            // collapsed substrate (dup=0.92). This channel has no such need --
            // it compares READS against ONE contig -- so it can use the
            // aggressively collapsed pileup substrate (dup=0.45), where reads
            // stack far better. Selectable so the two can be compared.
            // MEASURED: the mildly-collapsed BUBBLE substrate is the right one
            // here, despite this channel not needing haplotype separation the
            // way the bubble passes do. On the aggressively collapsed pileup
            // substrate (dup=0.45) indel F1 falls 0.655 -> 0.526 on r2, because
            // that collapse merges the two haplotypes into ONE contig: a read
            // carrying the indel then often matches a contig that already
            // contains it, and no gap appears at all.
            const CallData& pc_cd = std::getenv("CAPS_PCLUSTER_COLSUB") ? cd : cdb;
            // `cov` above is indexed by the BUBBLE substrate's contigs, so it
            // cannot be reused when this channel runs on a different one --
            // the ids would silently address the wrong contig. Build coverage
            // for whichever substrate this channel actually uses.
            // kidx is built over the BUBBLE substrate's contigs, so its ids
            // cannot be used to index a different substrate -- doing so
            // segfaulted immediately. This channel builds its OWN anchor index
            // over whichever substrate it uses.
            std::unordered_map<uint64_t, std::pair<uint32_t,uint32_t>> pkidx;
            {
                // Count occurrences of the CANONICAL k-mer across BOTH
                // orientations and keep only those seen exactly once, then use
                // only the forward-oriented ones. Counting just the
                // forward-canonical k-mers instead admits anchors that recur in
                // reverse-complement form elsewhere, which measured as a real
                // precision loss (r2 indel P 0.760 -> 0.586).
                std::unordered_map<uint64_t, uint32_t> pcount;
                std::unordered_map<uint64_t, uint8_t> porient;
                pcount.reserve(1u << 21); pkidx.reserve(1u << 21);
                for (uint32_t ci = 0; ci < (uint32_t)pc_cd.contigs.size(); ++ci) {
                    const std::string& c = pc_cd.contigs[ci];
                    for (size_t i2 = 0; i2 + 25 <= c.size(); ++i2) {
                        uint64_t v; if (!pack25(c.data() + i2, v)) continue;
                        uint64_t rv = rc25(v), cn = v < rv ? v : rv;
                        if (++pcount[cn] == 1) {
                            pkidx[cn] = {ci, (uint32_t)i2};
                            porient[cn] = (uint8_t)(v <= rv ? 0 : 1);
                        }
                    }
                }
                // ANCHOR UNIQUENESS. Originally occ==1 across ALL contigs.
                // Measured: near the 1bp-homopolymer events we miss this is
                // starving -- 851 read-level gap votes collapse to only 64
                // distinct events, 47 with a single supporting read, because
                // too few anchors survive for several reads to vote on the same
                // locus. In a fragmented substrate the same 25-mer legitimately
                // appears in several overlapping contigs without being a
                // repeat. Relax to occ<=PUNIQ (first-seen contig used).
                int PUNIQ = 1;
                if (const char* e = std::getenv("CAPS_PCLUSTER_UNIQ")) PUNIQ = atoi(e);
                for (auto it = pkidx.begin(); it != pkidx.end(); )
                    if (porient[it->first] != 0 || (int)pcount[it->first] > PUNIQ) it = pkidx.erase(it); else ++it;
                for (auto it = pkidx.begin(); it != pkidx.end(); )
                    if (pcount[it->first] != 1) it = pkidx.erase(it); else ++it;
            }
            std::vector<std::vector<uint16_t>> pcov;
            {
                pcov.resize(pc_cd.contigs.size());
                for (size_t ci = 0; ci < pc_cd.contigs.size(); ++ci)
                    pcov[ci].assign(pc_cd.contigs[ci].size(), 0);
                for (size_t oi = 0; oi < n; ++oi) {
                    uint32_t cid = pc_cd.read_cid[oi], pos = pc_cd.read_pos[oi];
                    if (cid >= pc_cd.contigs.size()) continue;
                    uint16_t clp = (oi < pc_cd.read_clip.size()) ? pc_cd.read_clip[oi] : 0;
                    int rl = (int)seqs[oi].size() - (int)clp;
                    for (int j = 0; j < rl; ++j) {
                        uint32_t pp = pos + (uint32_t)j;
                        if (pp < pcov[cid].size() && pcov[cid][pp] < 60000) ++pcov[cid][pp];
                    }
                }
            }
            // cluster reads by right-context anchor that is UNIQUE in the contigs
            //
            // ONLY ANCHORS PRESENT IN pkidx ARE EVER LOOKED UP (2026-09-03).
            // The sole consumer of this map is the `for (auto& kv : pkidx)`
            // loop below, whose first act is `rc_reads.find(kv.first)`. Any
            // key absent from pkidx is therefore hashed, stored, and never
            // read -- and that was the overwhelming majority of it, because
            // pkidx has just been filtered down to contig-UNIQUE anchors
            // while this loop indexes every 25-mer of every read, including
            // the error k-mers that dominate the distinct-key count (one
            // sequencing error spawns up to 25 novel k-mers). At full chr20
            // that is 12,604,917 reads x 84 offsets = 1.06 BILLION insertions
            // to serve a lookup set orders of magnitude smaller.
            //
            // Skipping non-pkidx keys is output-identical BY CONSTRUCTION,
            // not by measurement: a loop that only ever reads keys in pkidx
            // cannot observe the absence of keys that are not in pkidx. pkidx
            // is fully built and filtered above before this point, so the
            // membership test is well-defined here.
            std::unordered_map<uint64_t, std::vector<std::pair<uint32_t,uint32_t>>> rc_reads;
            rc_reads.reserve(pkidx.size() ? pkidx.size() : (size_t)(1u << 20));
            for (uint32_t i = 0; i < (uint32_t)seqs.size(); ++i) {
                const std::string& q = seqs[i];
                for (size_t j = LW; j + AK <= q.size(); ++j) {
                    uint64_t v; if (!pack25(q.data() + j, v)) continue;
                    uint64_t rv = rc25(v), cn = v < rv ? v : rv;
                    if (!pkidx.count(cn)) continue;     // provably dead otherwise
                    auto& vec = rc_reads[cn];
                    if (vec.size() < 200) vec.push_back({i, (uint32_t)j});
                }
            }
            static long g_pc_identical=0, g_pc_nogap=0, g_pc_found=0, g_pc_seen=0;
            // (contig, pos, signed gap, inserted seq) -> supporting reads
            // DISTINCT reads per event, not votes. A read spanning an indel
            // matches many anchors (every unique 25-mer in its right context),
            // so a raw vote count multiplies one read into dozens and is not a
            // read count at all -- which is why an absolute threshold kept
            // rising without saturating. Counting distinct read ids makes the
            // support interpretable and lets it be compared against coverage.
            std::map<std::tuple<std::string,int,std::string>, std::unordered_set<uint32_t>> pvotes;
            std::map<std::tuple<std::string,int,std::string>, std::pair<uint32_t,uint32_t>> ploc;
            // Distinct ANCHORS backing each event. Independent right-context
            // anchors are independent evidence, the same standard the bubble
            // channel applies via MIN_ANCH. Without it a single anchor's worth
            // of reads can carry an event on its own.
            std::map<std::tuple<std::string,int,std::string>, std::unordered_set<uint32_t>> panch;
            for (auto& kv : pkidx) {
                uint32_t ccid = kv.second.first;
                uint32_t cpos = kv.second.second;
                const std::string& cc2 = pc_cd.contigs[ccid];
                if (cpos < (uint32_t)LW) continue;
                auto rit = rc_reads.find(kv.first);
                if (rit == rc_reads.end()) continue;
                for (auto& pr : rit->second) {
                    if (std::getenv("CAPS_PCDBG")) ++g_pc_seen;
                    const std::string& q = seqs[pr.first];
                    uint32_t rpos = pr.second;
                    // orient the read so its anchor reads forward like the contig
                    std::string qq = q; uint32_t qp = rpos;
                    { uint64_t v; pack25(q.data() + rpos, v);
                      uint64_t rv = rc25(v);
                      if (v > rv) { qq = rc_str(q); qp = (uint32_t)(q.size() - rpos - AK); } }
                    if (qp < (uint32_t)LW) continue;
                    // walk LEFT from the anchor: contig and read agree, then diverge
                    int d = 0;
                    while (d < LW && cc2[cpos - 1 - d] == qq[qp - 1 - d]) ++d;
                    if (d >= LW) { if(std::getenv("CAPS_PCDBG")) ++g_pc_identical; continue; }   // identical: no event
                    // try a single gap of g bases at the divergence point
                    int best_g = 0; std::string best_ins;
                    for (int g = 1; g <= MAXINDEL && !best_g; ++g) {
                        // deletion in the READ: contig has g extra bases
                        if (cpos >= (uint32_t)(d + g + 10) && qp >= (uint32_t)(d + 10)) {
                            bool ok = true;
                            for (int t = 0; t < 10; ++t)
                                if (cc2[cpos - 1 - d - g - t] != qq[qp - 1 - d - t]) { ok = false; break; }
                            if (ok) { best_g = g; }
                        }
                        // insertion in the READ: read has g extra bases
                        if (!best_g && qp >= (uint32_t)(d + g + 10) && cpos >= (uint32_t)(d + 10)) {
                            bool ok = true;
                            for (int t = 0; t < 10; ++t)
                                if (cc2[cpos - 1 - d - t] != qq[qp - 1 - d - g - t]) { ok = false; break; }
                            if (ok) { best_g = -g; best_ins = qq.substr(qp - d - g, (size_t)g); }
                        }
                    }
                    if (!best_g) { if(std::getenv("CAPS_PCDBG")) ++g_pc_nogap; continue; }
                    if (std::getenv("CAPS_PCDBG")) ++g_pc_found;
                    uint32_t apos2 = cpos - (uint32_t)d;       // contig position of the event
                    if (apos2 == 0) continue;
                    // VOTE KEY. Originally (contig, pos, gap, ins) -- but the
                    // same genomic locus is covered by 2-7 different contigs
                    // (measured), and each read anchors to whichever contig its
                    // 25-mer happens to be unique in, so votes for ONE event
                    // scatter across several contig ids and never accumulate.
                    // Measured symptom: 851 read-level gap votes collapsing to
                    // 64 distinct events, 47 of them with a single supporting
                    // read. Key instead on the local CONTIG SEQUENCE around the
                    // event, which is identical across contigs covering the
                    // same locus and needs no reference genome.
                    std::string ctx;
                    {
                        size_t lo = (apos2 >= 12u) ? (size_t)apos2 - 12u : 0u;
                        size_t hi = std::min(cc2.size(), (size_t)apos2 + 12u);
                        if (hi > lo) ctx = cc2.substr(lo, hi - lo);
                    }
                    auto pkey = std::make_tuple(ctx, best_g, best_ins);
                    pvotes[pkey].insert(pr.first);
                    panch[pkey].insert(cpos);
                    ploc[pkey] = std::make_pair(ccid, apos2);
                }
            }
            if (std::getenv("CAPS_PCDBG")) {
                std::map<int,int> hist;
                for (auto& kv : pvotes) hist[(int)kv.second.size()]++;
                fprintf(stderr,"[pc] distinct events=%zu  support histogram:", pvotes.size());
                for (auto& h : hist) if (h.first<=8) fprintf(stderr," %dread(s)x%d", h.first, h.second);
                fprintf(stderr,"\n");
            }
            long pc_drop_min=0, pc_drop_aflo=0, pc_drop_afhi=0, pc_emit=0;
            int PMIN = MC;
            if (const char* e = std::getenv("CAPS_PCLUSTER_MIN")) PMIN = atoi(e);
            size_t n_pc = 0;
            for (auto& kv : pvotes) {
                const int nsup = (int)kv.second.size();
                if (nsup < PMIN) { ++pc_drop_min; continue; }
                // MEASURED NEUTRAL-TO-NEGATIVE, so no constraint by default.
                // Five windows: >=1 gives mean indel F1 0.6478 (identical to
                // no constraint, a useful sanity check that the code is inert
                // at 1), >=2 gives 0.655/0.699/0.635/0.554/... -- r4 falls
                // 0.569 -> 0.554. Nearly every real event already has multiple
                // anchors, so the requirement removes true calls without
                // removing false ones. Kept adjustable, defaulted off.
                int PANCH = 1;
                if (const char* e = std::getenv("CAPS_PCLUSTER_ANCH")) PANCH = atoi(e);
                if ((int)panch[kv.first].size() < PANCH) continue;
                std::string ctx_; int g; std::string insseq;
                std::tie(ctx_, g, insseq) = kv.first;
                uint32_t ccid = ploc[kv.first].first, apos2 = ploc[kv.first].second;
                const std::string& cc2 = pc_cd.contigs[ccid];
                if (apos2 == 0 || apos2 > cc2.size()) continue;
                // ALLELE FRACTION against the coverage already computed for
                // this contig. An absolute read-count floor cannot separate a
                // heterozygous indel from an artifact, because both scale with
                // depth; the discriminating quantity is what FRACTION of the
                // reads at this locus carry the event. A het indel sits near
                // 0.5. Reuses the frozen MAF rather than adding a constant.
                if (!std::getenv("CAPS_NO_PCLUSTER_AF")) {
                    int dp = (ccid < pcov.size() && apos2 < pcov[ccid].size())
                             ? (int)pcov[ccid][apos2] : 0;
                    if (dp > 0) {
                        double af = (double)nsup / (double)dp;
                        // Swept over five windows on DISTINCT-read support:
                        // mean indel F1 0.6446 (0.20) / 0.6414 (0.30) /
                        // ~0.643 (0.40). The optimum is 0.20 -- exactly the
                        // project's already-frozen MAF, so this adds no new
                        // constant.
                        double afmin = MAF;
                        if (const char* e = std::getenv("CAPS_PCLUSTER_AF")) afmin = atof(e);
                        if (af < afmin) { ++pc_drop_aflo; continue; }
                        // UPPER bound too. A HETEROZYGOUS indel is carried by
                        // roughly half the reads; if essentially EVERY read at
                        // the locus shows the gap, the event is homozygous with
                        // respect to our contig -- a real sequence difference,
                        // but not a heterozygous variant, and the truth set is
                        // het-restricted, so it can only ever score as a false
                        // positive. The SNV path has always applied a two-sided
                        // band for exactly this reason; the indel path had only
                        // a floor. Upper bound is the symmetric 1-MAF.
                        if (!std::getenv("CAPS_NO_PCLUSTER_HI") && af > 1.0 - MAF) { ++pc_drop_afhi; continue; }
                    }
                }
                char anch = cc2[apos2 - 1];
                std::string ref2, alt2;
                if (g > 0) {                                   // deletion in the read
                    if (apos2 + (uint32_t)g > cc2.size()) continue;
                    ref2 = std::string(1, anch) + cc2.substr(apos2, (size_t)g);
                    alt2 = std::string(1, anch);
                } else {                                       // insertion in the read
                    ref2 = std::string(1, anch);
                    alt2 = std::string(1, anch) + insseq;
                }
                char inf3[96];
                snprintf(inf3, sizeof inf3, "SVTYPE=INDEL;DP=%d;SOURCE=PCLUSTER", nsup);
                orecs.push_back({ccid, apos2, ref2, alt2, inf3, BUB_SRC});
                ++n_pc; ++n_indel; ++pc_emit;
            }
            if (std::getenv("CAPS_PCDBG"))
                fprintf(stderr,"[pc] drop_minsupport=%ld drop_AF_low=%ld drop_AF_high=%ld EMIT=%ld\n",
                        pc_drop_min,pc_drop_aflo,pc_drop_afhi,pc_emit);
            if (std::getenv("CAPS_PCDBG"))
                fprintf(stderr,"[pc] read-anchor visits=%ld  identical=%ld  no-gap=%ld  GAP-FOUND=%ld\n",
                        g_pc_seen,g_pc_identical,g_pc_nogap,g_pc_found);
            if (n_pc) fprintf(stderr, "[CAPS-CALL] pcluster indels=%zu\n", n_pc);
        }

        // ── 6b4. LINK-SCAN indel channel (option 2) ──────────────────────
        // The assembler has ALREADY placed every read (cd_in.read_cid/read_pos
        // from ppos[]), so this channel needs NO anchor index at all -- which
        // is the exact bottleneck that starved the positional-clustering
        // channel (851 read-level gap votes collapsing to 64 events, 47 with a
        // single supporting read, because too few unique anchors survive near
        // low-complexity loci).
        //
        // What it captures: a read carrying a 1bp homopolymer indel is still
        // placed on the contig of the OTHER haplotype, because exact-overlap
        // chaining absorbs a one-repeat-unit difference. Ungapped, such a read
        // matches up to the run and then mismatches everything after it. Re-
        // checking each placed read WITH a gap recovers exactly the event the
        // assembler absorbed, at the point where it was absorbed.
        // Uses the ORIGINAL assembler placement (cd_in), not the collapsed or
        // re-mapped substrates, since that is where the absorption happened.
        // MEASURED: this channel finds almost nothing and is OFF by default.
        // Of 75,115 placed reads, 63,081 (84%) match their contig PERFECTLY and
        // only 9 show a gap-explainable divergence -- because the assembler
        // does NOT absorb alt-haplotype reads into the wrong contig; it places
        // each read on a contig that already matches it. Kept behind
        // CAPS_LINKSCAN=1 because the measurement is the useful part.
        if (std::getenv("CAPS_LINKSCAN")) {
            const int TAIL_MIN = 20;          // bases of tail that must re-match
            long ls_nocid=0, ls_outside=0, ls_perfect=0, ls_shortleft=0, ls_shorttail=0, ls_nogap=0, ls_hit=0;
            std::map<std::tuple<std::string,int,std::string>, std::unordered_set<uint32_t>> lvotes;
            std::map<std::tuple<std::string,int,std::string>, std::pair<uint32_t,uint32_t>> lloc;
            for (size_t oi = 0; oi < n; ++oi) {
                uint32_t cid = cd_in.read_cid[oi];
                if (cid >= cd_in.contigs.size()) { ++ls_nocid; continue; }
                uint32_t pos = cd_in.read_pos[oi];
                const std::string& C = cd_in.contigs[cid];
                std::string q = cd_in.read_rc[oi] ? rc_str(seqs[oi]) : seqs[oi];
                uint16_t clp = (oi < cd_in.read_clip.size()) ? cd_in.read_clip[oi] : 0;
                if (clp) { if (clp >= q.size()) continue; q.erase(0, clp); }
                const int rl = (int)q.size();
                if (pos + (uint32_t)rl > C.size()) { ++ls_outside; continue; }
                // first ungapped divergence
                int j = 0;
                while (j < rl && C[pos + (size_t)j] == q[(size_t)j]) ++j;
                if (j >= rl) { ++ls_perfect; continue; }
                if (j < 12) { ++ls_shortleft; continue; }
                if (rl - j < TAIL_MIN + 2) { ++ls_shorttail; continue; }
                int best_g = 0; std::string best_ins;
                for (int g = 1; g <= MAXINDEL && !best_g; ++g) {
                    // DELETION in the read: contig carries g extra bases
                    if (pos + (size_t)(j + g + TAIL_MIN) <= C.size()) {
                        int mm = 0, cnt = 0;
                        for (int t = 0; t < TAIL_MIN; ++t) {
                            if (j + t >= rl) break;
                            ++cnt;
                            if (C[pos + (size_t)(j + g + t)] != q[(size_t)(j + t)]) ++mm;
                        }
                        if (cnt >= TAIL_MIN && mm <= 1) best_g = g;
                    }
                    // INSERTION in the read: read carries g extra bases
                    if (!best_g && j + g + TAIL_MIN <= rl &&
                        pos + (size_t)(j + TAIL_MIN) <= C.size()) {
                        int mm = 0, cnt = 0;
                        for (int t = 0; t < TAIL_MIN; ++t) {
                            ++cnt;
                            if (C[pos + (size_t)(j + t)] != q[(size_t)(j + g + t)]) ++mm;
                        }
                        if (cnt >= TAIL_MIN && mm <= 1) {
                            best_g = -g; best_ins = q.substr((size_t)j, (size_t)g);
                        }
                    }
                }
                if (!best_g) { ++ls_nogap; continue; }
                ++ls_hit;
                uint32_t apos3 = pos + (uint32_t)j;
                if (apos3 == 0 || apos3 >= C.size()) continue;
                std::string ctx;
                {
                    size_t lo = (apos3 >= 12u) ? (size_t)apos3 - 12u : 0u;
                    size_t hi = std::min(C.size(), (size_t)apos3 + 12u);
                    if (hi > lo) ctx = C.substr(lo, hi - lo);
                }
                auto lk = std::make_tuple(ctx, best_g, best_ins);
                lvotes[lk].insert((uint32_t)oi);
                lloc[lk] = std::make_pair(cid, apos3);
            }
            // per-contig coverage of the ORIGINAL substrate, for allele fraction
            std::vector<std::vector<uint16_t>> lcov(cd_in.contigs.size());
            for (size_t ci = 0; ci < cd_in.contigs.size(); ++ci)
                lcov[ci].assign(cd_in.contigs[ci].size(), 0);
            for (size_t oi = 0; oi < n; ++oi) {
                uint32_t cid = cd_in.read_cid[oi];
                if (cid >= cd_in.contigs.size()) continue;
                uint32_t pos = cd_in.read_pos[oi];
                uint16_t clp = (oi < cd_in.read_clip.size()) ? cd_in.read_clip[oi] : 0;
                int rl = (int)seqs[oi].size() - (int)clp;
                for (int t = 0; t < rl; ++t) {
                    uint32_t pp = pos + (uint32_t)t;
                    if (pp < lcov[cid].size() && lcov[cid][pp] < 60000) ++lcov[cid][pp];
                }
            }
            int LMIN = MC;
            if (const char* e = std::getenv("CAPS_LINKSCAN_MIN")) LMIN = atoi(e);
            size_t n_ls = 0;
            for (auto& kv : lvotes) {
                const int nsup = (int)kv.second.size();
                if (nsup < LMIN) continue;
                std::string ctx_; int g; std::string insseq;
                std::tie(ctx_, g, insseq) = kv.first;
                uint32_t cid = lloc[kv.first].first, apos3 = lloc[kv.first].second;
                const std::string& C = cd_in.contigs[cid];
                if (apos3 == 0 || apos3 > C.size()) continue;
                int dp = (cid < lcov.size() && apos3 < lcov[cid].size()) ? (int)lcov[cid][apos3] : 0;
                if (dp > 0) {
                    double af = (double)nsup / (double)dp;
                    if (af < MAF || af > 1.0 - MAF) continue;     // heterozygous band
                }
                char anch = C[apos3 - 1];
                std::string ref3, alt3;
                if (g > 0) {                                       // deletion in read
                    if (apos3 + (uint32_t)g > C.size()) continue;
                    ref3 = std::string(1, anch) + C.substr(apos3, (size_t)g);
                    alt3 = std::string(1, anch);
                } else {
                    ref3 = std::string(1, anch);
                    alt3 = std::string(1, anch) + insseq;
                }
                char inf4[96];
                snprintf(inf4, sizeof inf4, "SVTYPE=INDEL;DP=%d;SOURCE=LINKSCAN", nsup);
                orecs.push_back({cid, apos3, ref3, alt3, inf4, 2});
                ++n_ls; ++n_indel;
            }
            if (std::getenv("CAPS_PCDBG"))
                fprintf(stderr, "[ls] reads=%zu nocid=%ld outside=%ld perfect=%ld shortleft=%ld shorttail=%ld nogap=%ld HIT=%ld events=%zu emitted=%zu\n",
                        n, ls_nocid, ls_outside, ls_perfect, ls_shortleft, ls_shorttail, ls_nogap, ls_hit, lvotes.size(), n_ls);
            if (n_ls) fprintf(stderr, "[CAPS-CALL] linkscan indels=%zu\n", n_ls);
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
            for (size_t xrun_i = 0; xrun_i < kidx.size(); ) {
                size_t xrun_j = xrun_i;
                while (xrun_j < kidx.size() && kidx[xrun_j].kmer == kidx[xrun_i].kmer) ++xrun_j;
                size_t run_len = xrun_j - xrun_i;
                size_t base = xrun_i;
                xrun_i = xrun_j;
                if (run_len != 2) continue;
                uint32_t cax = kidx[base].ci,   pax = kidx[base].pos;   uint8_t oax = kidx[base].orient;
                uint32_t cbx = kidx[base+1].ci, pbx = kidx[base+1].pos; uint8_t obx = kidx[base+1].orient;
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
            for (size_t ci = 0; ci < cd_in.contigs.size(); ++ci)
                fprintf(df, ">lcontig_%zu\n%s\n", ci, cd_in.contigs[ci].c_str());
            // dBG bubbles are their own contig class. Anchoring them onto the
            // pileup substrate was tried first and scored EXACTLY the baseline
            // (332/14/70): that substrate is collapsed (dup=0.45), so both
            // haplotypes are already merged there and every anchored bubble
            // landed on a site the pileup had. The bubbles that rescue are
            // precisely the ones that fail to anchor. Emitting each bubble as
            // its own sequence keeps its two paths intact and lets the lift
            // place it independently.
            for (size_t ci = 0; ci < dbg_bubbles.size(); ++ci)
                fprintf(df, ">dcontig_%zu\n%s%s%s%s\n", ci,
                        dbg_bubbles[ci].lext.c_str(), dbg_bubbles[ci].flank.c_str(),
                        dbg_bubbles[ci].path1.c_str(), dbg_bubbles[ci].rext.c_str());
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
    for (size_t ci = 0; ci < cd_in.contigs.size(); ++ci)
        fprintf(f, "##contig=<ID=lcontig_%zu,length=%zu>\n", ci, cd_in.contigs[ci].size());
    for (size_t ci = 0; ci < dbg_bubbles.size(); ++ci)
        fprintf(f, "##contig=<ID=dcontig_%zu,length=%zu>\n", ci,
                dbg_bubbles[ci].lext.size() + dbg_bubbles[ci].flank.size() +
                dbg_bubbles[ci].path1.size() + dbg_bubbles[ci].rext.size());
    fprintf(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (auto& r : orecs)
        fprintf(f, "%s%u\t%u\t.\t%s\t%s\t.\tPASS\t%s\n",
                r.src == 3 ? "dcontig_" :
                (r.src == 2 ? "lcontig_" : (r.src ? "bcontig_" : "contig_")), r.cid, r.pos,
                r.ref.c_str(), r.alt.c_str(), r.info.c_str());
    fclose(f);
    fprintf(stderr, "[CAPS-CALL] contigs=%zu H=%u candidates=%zu SNVs=%zu indels=%zu -> %s\n",
            cd.contigs.size(), H, total_candidates, n_snv, n_indel, out_vcf.c_str());
    phase("indel_pass", t_mark);
    fprintf(stderr, "[CAPS-CALL-TIMING] %-16s %8.3fs\n", "TOTAL", elapsed_s(t_start, clk::now()));

    return (int)orecs.size();
}

} // namespace capscall
