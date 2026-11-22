#pragma once
// ── CAPSULE Claim-2 caller ───────────────────────────────────────────────────
// Reference-free heterozygous SNV and indel caller that works on the reads of a
// G_CAPSUL archive, with no reference genome. It is a port of the ARCS project's
// bubble caller (a separate codebase, not in this repository), retargeted at
// CAPSULE's own assembly output instead of ARCS's build_vodbg_pg / CallData.
//
// ENTRY POINT. capscall::run_variant_call(seqs, quals, cd, out_vcf) writes a VCF
// to out_vcf and returns the number of records written, or -1 on failure (cd.valid
// false, or the output file cannot be opened). `seqs` and `quals` are the original
// reads in file order. `cd` (struct CallData) says where each read sits: its contig,
// its offset in that contig, its strand, its left clip and its global pg offset.
// Two callers exist:
//   * `capsule_decode call` (src/decoder.cpp, capsule_call_from_archive) builds
//     CallData from an archive compressed with CAPS_SPANS=1 or CAPS_CALL=1 (which
//     store the contig spans). By default it runs the GRAPH-ONLY path: it sets
//     CAPS_DBG=1 and CAPS_DBG_ONLY=1, packs the reads with capspack (include/caps_pack.h)
//     and passes a quality bitmap. With CAPS_CALL_INDELS=1 it unsets both and runs the
//     FULL caller on plain-text reads. The packed or plain input format depends on
//     those two variables (SEQ_PACKED below), so they must agree with how the reads
//     were prepared.
//   * the encoder, when CAPS_CALL=1, calls it at the end of compression, with the
//     in-memory reads and no DBG variables set unless the user sets them.
//
// PIPELINE inside run_variant_call (the banner comments below carry these names):
//   substrate   build_substrate: drop contigs that are near-duplicates of a longer
//               one (collapse_contigs, off with CAPS_NO_COLLAPSE), then re-place the
//               reads on the surviving contigs (off with CAPS_NO_REMAP). SNVs are
//               called on this collapsed substrate. Indel bubbles use it too by default.
//   1. kc       canonical 31-mer counts of all reads. Spills to disk past a memory
//               ceiling (CAPS_MAXRAM_MB), see the LAYER 2 / LAYER 3 comments.
//   2-6a.       per-contig pileup, SNV candidates, flank pass, filters, SNV records.
//   DBG channel de Bruijn bubble channel (CAPS_DBG=1, and only on a sample that looks
//               diploid or with CAPS_DBG_FORCE), DiscoSNP++'s algorithm run on the graph
//               already built. With CAPS_DBG_ONLY=1 it is the ONLY channel: the pileup, the
//               second substrate and the indel pass are skipped (see DBG-ONLY MODE).
//   5. filters  the frozen filters (parameters below) turn candidates into kept calls.
//   6b. indels  het indels are BUBBLES between two haplotype contigs. Channels: the
//               contig-pair bubble scan (with closing anchor, junction support and
//               maximum-weight matching), 6b3 positional clustering (off with
//               CAPS_NO_PCLUSTER), 6b4 link-scan (CAPS_LINKSCAN), 6b2 cross-contig SNVs
//               (off with CAPS_NO_XSNV).
//   6c. VCF     sort and write. Records from the SNV pass and from the bubble pass use
//               different contig spaces (contig_N and bcontig_N), see CAPS_DUMP_CONTIGS.
// The switches are listed with their defaults in docs/configuration_reference.md.
// The stages are all inside one function, so a reader should follow the banner
// comments, not function boundaries.
//
// Frozen parameters (do not re-tune: these are ARCS's own validated, held-out-tested
// values, ported verbatim):
//   HDMAX=2, MAF=0.20, DHI=2.5, KHI=1.1, MC=3, TRI=0.12, HALF=15.
// They are described where they are used. The long-tail tuning and diagnostic
// switches (CAPS_DBG_*, CAPS_PCLUSTER_*, CAPS_TRACE and so on) are experiments that
// do not change the default result. The measured outcomes of those experiments are
// quoted next to the code that was tried.
//
// WHY THE CONTIGS COME FROM THE ENCODER. ARCS's assembler keeps reads in several
// DISCRETE contigs by construction. CAPSULE's own assembler stitches everything
// into one (or two) long pg buffers via MEM references, which would erase the
// haplotype separation the indel-bubble mechanism needs. The fix (src/encoder.cpp,
// CAPS_CALL / CAPS_SPANS gate): every chain emitted in round 1 and round 2 (main
// pg) and in the second-region sweep is its own contiguous span of `pg`, captured as
// a contig BEFORE the MEM self-match stage runs. MEM never mutates `pg` itself, only
// which of its bytes get literal-encoded, so contig boundaries captured pre-MEM are
// still valid after MEM completes. Pigeonhole-mapped reads land inside an
// already-emitted span by construction (they only match against already-appended pg
// content), so they inherit that span's contig id for free.
//
// This header owns ONLY the caller logic. Contig-span capture and CallData
// population live in src/encoder.cpp (gated, zero cost when unset, like CAPS_NAMES
// and CAPS_QUAL) and in src/decoder.cpp (`call`).
//
// MEASUREMENT HISTORY. The comments below state each design decision and the
// measurement behind it. The full measurements and rejected approaches are in
// docs/extras/ENGINEERING_NOTES_CALLER.md; the reported results are in docs/claim2/.
// Measurement labels such as "HG002 r2" are the chr20 tuning and held-out windows.
// "DiscoSNP++" source files that are named (Bubble.cpp, Graph.cpp and so on) are the
// third-party DiscoSNP++ and gatb-core code.
//
// The numeric results of this caller are produced by the code below. This header has
// no unit tests of its own. scripts/claim2_scripts/test_claim2.sh is the fast synthetic
// check and the benchmark scripts are the real-data check.

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <atomic>
#include <sys/mman.h>
#ifdef _OPENMP
#include <parallel/algorithm>
#endif
#include <string>
#include <array>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cerrno>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
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
    // GLOBAL pseudogenome offset of each original read (UINT64_MAX = unplaced).
    // contig id + contig-local offset cannot be compared ACROSS contigs, and
    // the substrate is fragmented enough that one genomic locus routinely spans
    // several -- measured: allele placements were "disjoint" for 369 of 486
    // bubbles at contig granularity, which turned out to carry no signal
    // (precision 0.952 disjoint vs 0.916 shared, against 0.944 overall). A
    // single global coordinate is what makes "were these reads near each other"
    // answerable at all. The encoder already has it as ppos[u].
    std::vector<uint64_t>    read_ppos;
    bool valid = false;
};

namespace detail {

constexpr int HDMAX = 2, MC = 3, HALF = 15;
constexpr double MAF = 0.20, DHI = 2.5, KHI = 1.1, TRI = 0.12;
constexpr int MAXINDEL = 30;

// Base helpers. b2i maps A,C,G,T to 0..3 and anything else (N, lowercase) to -1.
// compl_base complements a base and returns 'N' for a non-ACGT character. rc_str is
// the reverse complement of a string. pack31 packs 31 bases into the low 62 bits
// (first base in the high bits) and returns false if any base is not ACGT.
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
// BIT-PARALLEL, not a 31-iteration loop.
//
// canon31 calls this on all four candidate successors of every solid node, and
// the traversal calls canon31 twice per node on top of that -- so an O(k) loop
// here is ~124 iterations per node before any lookup happens.
//
// Complementing a 2-bit code is a bitwise NOT (3-b == ~b & 3), and reversing
// the 2-bit groups is the standard halving shuffle, so the whole thing is ~10
// operations regardless of k. The final shift re-aligns the 62-bit word.
// Checked against the loop over 8,000,004 values including edge cases: zero
// mismatches.
inline uint64_t rc31(uint64_t v) {
    v = ~v;
    v = ((v >> 2)  & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
    v = ((v >> 4)  & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
    v = ((v >> 8)  & 0x00FF00FF00FF00FFULL) | ((v & 0x00FF00FF00FF00FFULL) << 8);
    v = ((v >> 16) & 0x0000FFFF0000FFFFULL) | ((v & 0x0000FFFF0000FFFFULL) << 16);
    v = (v >> 32) | (v << 32);
    return (v >> 2) & ((~0ULL) >> 2);
}
// Sum of the indel_pass section timers, checked against the stage total.
static double g_ipsum = 0.0;
// Last lap timestamp, published so the segment AFTER the indel block -- the
// cross-contig SNV pass, which is default-ON and was never instrumented -- can
// be timed from outside the block where the _iplap lambda lives.
static std::chrono::steady_clock::time_point g_ip_last;
static std::chrono::steady_clock::time_point g_ip_start;
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

// One indel read off a pair of contigs A and B that diverge. type 0: A has bases that
// B lacks (a deletion in B), len bases long. type 1: B has bases that A lacks (an
// insertion in B), with the inserted bases in `ins`. apos is the position in A where
// they diverge. ok is false when no bubble was found.
struct Bubble { int type = 0; uint32_t apos = 0; int len = 0; std::string ins; bool ok = false; };

inline bool flank_match(const std::string& A, size_t ai, const std::string& B, size_t bi, int F) {
    if (ai + (size_t)F > A.size() || bi + (size_t)F > B.size()) return false;
    for (int k = 0; k < F; ++k) if (A[ai + (size_t)k] != B[bi + (size_t)k]) return false;
    return true;
}

// Not called. Re-convergence test that tolerates up to `tol` mismatches in the
// flank, so that a heterozygous SNV next to an indel does not break the bubble.
// It lowered indel F1 on every chr20 test window (precision fell more than recall
// rose), so extract_bubble keeps the exact flank_match. Measurements: docs/extras/ENGINEERING_NOTES_CALLER.md.
inline bool flank_match_tol(const std::string& A, size_t ai, const std::string& B, size_t bi,
                            int F, int tol) {
    if (ai + (size_t)F > A.size() || bi + (size_t)F > B.size()) return false;
    int mm = 0;
    for (int k = 0; k < F; ++k)
        if (A[ai + (size_t)k] != B[bi + (size_t)k] && ++mm > tol) return false;
    return true;
}

// Looks for an indel between contig A (from position pA) and contig B (from qB), the
// positions where the comparison starts. It first extends the exact match to the
// point d where they differ. A bubble is accepted only if the two sequences then re-converge
// exactly, over FLANK identical bases, after an indel of at most `maxindel` bases. Tried in
// order: the homopolymer run-length branch, the short-tandem-repeat branch (opt-in,
// CAPS_STRBUBBLE), then the generic shift search, insertion before deletion. Returns ok=false
// when nothing closes.
inline Bubble extract_bubble(const std::string& A, uint32_t pA, const std::string& B, uint32_t qB,
                              int maxindel, int FLANK) {
    Bubble r;
    const size_t la = A.size(), lb = B.size();
    uint32_t d = 0;
    while (pA + d < la && qB + d < lb && A[pA + d] == B[qB + d]) ++d;
    if (d == 0) return r;
    if (pA + d >= la || qB + d >= lb) return r;

    // ── HOMOPOLYMER RUN-LENGTH BRANCH ────────────────────────────────────────
    // A 1 bp indel inside a homopolymer changes the run LENGTH, so the two
    // haplotypes stay identical through the run and their flanks are offset by the
    // length difference. No shift g in the generic search below then gives an exact
    // flank match. Here the run length is measured on each side, the difference is
    // emitted as the indel, and re-convergence is checked after both runs. (DiscoSNP++
    // reaches the same events by extending a graph path base by base.) These events
    // were the whole of the 1 bp indels missed on the HG002 tuning window; see
    // docs/extras/ENGINEERING_NOTES_CALLER.md.
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
    // Generalizes the homopolymer branch above from unit length 1 to repeat units of
    // 2 to 6 bases: count how far each side continues the unit and emit the length
    // difference, so the sign (insertion or deletion) is correct by construction.
    // U=1 is left to the homopolymer branch. At least two copies of the unit must
    // precede the divergence, so the branch fires only inside a repeat.
    //
    // OPT-IN (CAPS_STRBUBBLE=1), OFF by default: the branch fires, but the caller's
    // VCF was byte-identical with it on and off on the HG002 tuning window, so the
    // loci it reaches are already resolved by the generic search. See
    // docs/extras/ENGINEERING_NOTES_CALLER.md.
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
// Both allele contigs usually exist, but near a homopolymer every 25-mer that
// overlaps the run differs between haplotypes, so the nearest shared anchor can be
// hundreds of bases away and an outward walk from it ends on short contigs. This
// takes the diagonal implied by ANY shared anchor and scans the whole co-linear
// region, resolving each divergence with a gap and continuing on the updated
// diagonal, so an event anywhere along the pair is reachable from one anchor.
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

// ── ROLLING CANONICAL 25-MER ────────────────────────────────────────────────
// Packing each 25-mer from scratch costs 25 steps per strand per position, although
// consecutive windows differ by one base. Roll25 keeps both strands and updates them
// in O(1) per base.
//
// A sliding window updates both strands in O(1):
//   forward  drop the top base, shift in the new one at the bottom
//   reverse  drop the bottom base, shift in the new complement at the top
// which is exactly rc25's output for the shifted window.
//
// N-HANDLING IS EXACT, not approximate. pack25 fails if ANY of the 25 bases is
// non-ACGT, so `run` counts consecutive valid bases ending at the current one
// and ok() requires 25 -- that is the same predicate, position for position.
struct Roll25 {
    static constexpr uint64_t M50 = (1ULL << 50) - 1;
    uint64_t f = 0, r = 0; int run = 0;
    inline void push(char ch) {
        const int b = b2i(ch);
        if (b < 0) { run = 0; f = 0; r = 0; return; }
        f = ((f << 2) | (uint64_t)b) & M50;
        r = (r >> 2) | ((uint64_t)(3 - b) << 48);
        if (run < 25) ++run;
    }
    inline bool ok() const { return run >= 25; }
    inline uint64_t canon() const { return f < r ? f : r; }
};

// Log-likelihood of the observed (base, quality) pairs `rl` under a pileup whose
// alleles are al[0..m). Each observation is explained by an allele chosen uniformly: it
// matches with probability 1-e and mismatches with e/3, where e = 10^(-Q/10).
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
// Why the substrate is rebuilt. The encoder assigns each read to the chain it
// overlaps EXACTLY, so reads carrying the reference and the alternate allele land on
// different chains and never meet in one pileup, and reads covering one site are
// spread over several redundant contigs. The rebuild therefore (1) collapses
// redundant contigs to one representative per locus and (2) re-places every read
// with mismatches allowed. On the HG002 tuning window the reads themselves carried
// both alleles at almost every truth site; the assignment was the limit. Figures in
// docs/extras/ENGINEERING_NOTES_CALLER.md.
//
// Both passes run for calling only and never touch the compression path.
// CAPS_NO_REMAP=1 / CAPS_NO_COLLAPSE=1 disable them.
struct GapEvent;

// The contigs and read placements the caller works on after build_substrate has
// collapsed redundant contigs and re-placed the reads. Same fields as CallData
// (placement per original read), plus the optional read-level gapped evidence.
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
// ── FLAT OPEN-ADDRESSING SET ────────────────────────────────────────────────
// The collapse walk is almost entirely hash insertion. A flat open-addressing table
// (one contiguous array, one cache line per probe, no per-element allocation) runs
// the same insert/query mix about 16x faster than std::unordered_set. Semantics are
// identical: incremental insert and query, so the accept/reject decision for each
// contig is unchanged. Measurements: docs/extras/ENGINEERING_NOTES_CALLER.md.
// Ask for 2 MB pages on a big array that will be probed randomly.
//
// ORDER MATTERS: madvise only affects pages not yet faulted in, so each call must sit
// immediately after the allocation and before the first write. The log line names
// the arrays that received the hint.
inline long caps_hwm_mb() {
    FILE* f = fopen("/proc/self/status", "r"); if (!f) return -1;
    char l[256]; long kb = -1;
    while (fgets(l, sizeof l, f)) if (!strncmp(l, "VmHWM:", 6)) { sscanf(l + 6, "%ld", &kb); break; }
    fclose(f); return kb / 1024;
}
inline long caps_rss_mb() {
    FILE* f = fopen("/proc/self/status", "r"); if (!f) return -1;
    char l[256]; long kb = -1;
    while (fgets(l, sizeof l, f)) if (!strncmp(l, "VmRSS:", 6)) { sscanf(l + 6, "%ld", &kb); break; }
    fclose(f); return kb / 1024;
}

inline void hp_hint(void* base, size_t len, const char* what) {
    if (!base || len < (2u << 20)) return;
    uintptr_t a0 = ((uintptr_t)base + (2u<<20) - 1) & ~(uintptr_t)((2u<<20) - 1);
    size_t off = a0 - (uintptr_t)base;
    if (len <= off) return;
    size_t l2 = (len - off) & ~(size_t)((2u<<20) - 1);
    if (!l2) return;
#ifdef __linux__
    if (madvise((void*)a0, l2, MADV_HUGEPAGE) == 0)
        fprintf(stderr, "[HP] MADV_HUGEPAGE on %.2f GB of %s\n", (double)l2/1073741824.0, what);
    else
        fprintf(stderr, "[HP] madvise failed (errno %d) -- %s stays on 4K pages\n", errno, what);
#else
    (void)a0; (void)l2; (void)what;
#endif
}

struct FlatKmerSet {
    std::vector<uint64_t> t;
    size_t mask = 0, live = 0, lim = 0;
    static constexpr uint64_t EMPTY = 0xFFFFFFFFFFFFFFFFULL;   // no 50-bit k-mer is all ones
    // SIZE FROM WHAT GOES IN, NOT FROM A LOOSE UPPER BOUND.
    //
    // The table starts small and doubles past 50% load. A hash set's contents do not
    // depend on its capacity, so which contigs survive the collapse is unchanged;
    // only the memory the probes range over shrinks. (A table sized from total k-mer
    // positions ran at about 8% load over several GB.)
    void alloc(size_t b) {
        t.assign(b, EMPTY); mask = b - 1; lim = b / 2;
    }
    void init(size_t cap) {
        // Start at the largest power of two not exceeding `cap` (total k-mer positions,
        // an upper bound on distinct k-mers). This lands near 45% load and avoids a
        // grow-and-rehash step with both tables live. Only the starting size changes,
        // never the contents; grow() remains as the safety net.
        size_t want = 1024;
        while (want * 2 <= cap) want <<= 1;
        alloc(want);
        live = 0;
        hp_hint((void*)t.data(), t.size() * sizeof(uint64_t), "collapse claimed set");
    }
    static inline uint64_t mix(uint64_t x) {
        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; return x;
    }
    inline bool count(uint64_t k) const {
        size_t i = mix(k) & mask;
        for (;;) { const uint64_t v = t[i]; if (v == EMPTY) return false;
                   if (v == k) return true; i = (i + 1) & mask; }
    }
    inline void put(uint64_t k) {                    // no growth check: used by grow()
        size_t i = mix(k) & mask;
        for (;;) { const uint64_t v = t[i]; if (v == k) return;
                   if (v == EMPTY) { t[i] = k; return; } i = (i + 1) & mask; }
    }
    void grow() {
        std::vector<uint64_t> old; old.swap(t);
        alloc(old.size() * 2);
        hp_hint((void*)t.data(), t.size() * sizeof(uint64_t), "collapse claimed set (grown)");
        for (uint64_t v : old) if (v != EMPTY) put(v);
    }
    inline void insert(uint64_t k) {
        size_t i = mix(k) & mask;
        for (;;) { const uint64_t v = t[i]; if (v == k) return;
                   if (v == EMPTY) { t[i] = k; if (++live > lim) grow(); return; }
                   i = (i + 1) & mask; }
    }
};

inline std::vector<uint32_t> collapse_contigs(const std::vector<std::string>& ctgs,
                                              double dup_frac, int K) {
    std::vector<uint32_t> order(ctgs.size());
    for (uint32_t i = 0; i < ctgs.size(); ++i) order[i] = i;
    // KNOWN NONDETERMINISM, left as is: equal-length contigs are ordered by
    // std::sort's unspecified order. Adding an index tie-break changes which contigs
    // survive the collapse and therefore the downstream work and calls, so it is a
    // behaviour change that would need its own F1 measurement, not a cleanup.
    std::sort(order.begin(), order.end(),
              [&](uint32_t a, uint32_t b){ return ctgs[a].size() > ctgs[b].size(); });
    // Sized from the actual input instead of a fixed 1M: total k-mer positions
    // is a hard upper bound on distinct claims, and rehashing a set that grows
    // to tens of millions from a 1M seed is pure copying. Capped so a huge
    // input does not over-allocate.
    FlatKmerSet claimed;
    {
        size_t tot_k = 0;
        for (const auto& c : ctgs) if ((int)c.size() >= K) tot_k += c.size() - (size_t)K + 1;
        claimed.init(tot_k ? tot_k : 1024);
    }
    std::vector<uint32_t> keep;
    std::vector<uint64_t> km;
    std::vector<uint8_t> kmq;   // was this k-mer at a stride-5 position?
    for (uint32_t ci : order) {
        const std::string& c = ctgs[ci];
        if ((int)c.size() < K) { keep.push_back(ci); continue; }
        // Query on a SAMPLED stride but claim EVERY k-mer: a duplicate contig
        // is almost always phase-shifted relative to its twin, so sampling both
        // sides at the same stride compares disjoint k-mer sets and the
        // duplicate is never detected (measured: only 7.6% of contigs collapsed
        // with same-stride sampling).
        // ONE walk per contig collects every canonical k-mer (all are claimed) and
        // flags the ones at stride-5 contig positions (those are queried). The source
        // position is recorded per entry because positions where pack25 fails are
        // skipped, so "every 5th collected k-mer" would not be the same set. `hit`,
        // `tot` and the accept/reject decision are those of a separate query walk
        // followed by a claim walk.
        km.clear(); kmq.clear();
        {   // rolling window; identical output to the pack25/rc25 form above
            Roll25 rw; const size_t cn2 = c.size();
            for (size_t b2 = 0; b2 < cn2; ++b2) {
                rw.push(c[b2]);
                if (b2 + 1 < (size_t)K) continue;
                if (!rw.ok()) continue;
                const size_t i = b2 + 1 - (size_t)K;
                km.push_back(rw.canon());
                kmq.push_back((uint8_t)(i % 5 == 0));
            }
        }
        size_t hit = 0, tot = 0;
        for (size_t w = 0; w < km.size(); ++w)
            if (kmq[w]) { ++tot; if (claimed.count(km[w])) ++hit; }
        if (tot > 0 && (double)hit / (double)tot >= dup_frac) continue;   // redundant
        keep.push_back(ci);
        for (uint64_t kk : km) claimed.insert(kk);
    }
    fprintf(stderr, "[COLLAPSE] claimed set: %zu distinct in %zu slots (%.0f MB, load %.0f%%)\n",
            claimed.live, claimed.t.size(), (double)claimed.t.size() * 8 / 1048576.0,
            100.0 * (double)claimed.live / (double)claimed.t.size());
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
        // Take the two DOMINANT implied starts rather than requiring every seed to
        // agree: one spurious repeat hit inside the same contig would otherwise
        // break the test.
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
    // [BS-SPLIT] stderr timers for the collapse, index and placement steps. They
    // only print and do not affect results.
    using _bsclk = std::chrono::steady_clock;
    auto _bst0 = _bsclk::now(); auto _bsl = _bst0;
    auto _bslap = [&](const char* nm){ auto now=_bsclk::now();
        fprintf(stderr, "[BS-SPLIT] %-18s %7.2fs  rss=%ldMB peak=%ldMB\n", nm,
                std::chrono::duration<double>(now-_bsl).count(),
                caps_rss_mb(), caps_hwm_mb()); _bsl=now; };

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
    _bslap("collapse");
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

    // 2. seed index over surviving contigs.
    //
    // PACKED SORTED ARRAY, not a hash map: one binary search and a bounded scan per
    // lookup instead of pointer chasing, on the caller's largest phase.
    //
    // INVARIANT: each key keeps at most 64 hits, and they must be the FIRST 64 in
    // (contig, position) order, the order in which the hits are produced. The
    // comparator (mix, key, ci, pos) states that order explicitly, so an unstable
    // in-place sort inside each bucket gives the same sequence a stable sort would,
    // and the in-place 64-cap keeps exactly those 64. One array is used throughout,
    // with no sort temporary or second copy (this step set the run's peak memory).
    struct SIdx { uint64_t key; uint32_t ci, pos; };
    std::vector<SIdx> sidx;
    static constexpr int BDIR = 24;          // directory bits (unchanged)
    static constexpr int SBB  = 12;          // scatter buckets: a prefix of BDIR,
                                             // so a directory group never straddles two
    std::vector<uint32_t> sdir;
    {
        const size_t NC = S.contigs.size();
        const size_t NSB = (size_t)1 << SBB;
        int T = 1;
        #ifdef _OPENMP
        T = omp_get_max_threads();
        #endif
        // 1. per-contig k-mer counts -> a balanced contig range per thread.
        // The ranges are EXPLICIT rather than an OpenMP schedule because the
        // histogram pass and the scatter pass must see the same contig->thread
        // map; a dynamic schedule does not promise that.
        std::vector<size_t> ccnt(NC + 1, 0);
        #pragma omp parallel for schedule(dynamic, 256)
        for (long long ci = 0; ci < (long long)NC; ++ci) {
            const std::string& c = S.contigs[(size_t)ci];
            size_t cnt = 0; Roll25 rw; const size_t cn2 = c.size();
            for (size_t b2 = 0; b2 < cn2; ++b2) {
                rw.push(c[b2]);
                if (b2 + 1 >= (size_t)K && rw.ok()) ++cnt;
            }
            ccnt[(size_t)ci + 1] = cnt;
        }
        for (size_t i = 0; i < NC; ++i) ccnt[i + 1] += ccnt[i];
        const size_t total = ccnt[NC];
        std::vector<size_t> tlo((size_t)T + 1, NC);
        tlo[0] = 0;
        for (int t = 1; t < T; ++t) {
            const size_t want = total * (size_t)t / (size_t)T;
            size_t at = (size_t)(std::lower_bound(ccnt.begin(), ccnt.end(), want) - ccnt.begin());
            if (at > NC) at = NC;
            if (at < tlo[(size_t)t - 1]) at = tlo[(size_t)t - 1];
            tlo[(size_t)t] = at;
        }
        tlo[(size_t)T] = NC;
        // 2. per-thread bucket histogram over its own contig range
        std::vector<std::vector<uint32_t>> hist((size_t)T, std::vector<uint32_t>(NSB, 0));
        #pragma omp parallel for schedule(static, 1)
        for (int t = 0; t < T; ++t) {
            std::vector<uint32_t>& h = hist[(size_t)t];
            for (size_t ci = tlo[(size_t)t]; ci < tlo[(size_t)t + 1]; ++ci) {
                const std::string& c = S.contigs[ci];
                Roll25 rw; const size_t cn2 = c.size();
                for (size_t b2 = 0; b2 < cn2; ++b2) {
                    rw.push(c[b2]);
                    if (b2 + 1 < (size_t)K || !rw.ok()) continue;
                    ++h[(size_t)(FlatKmerSet::mix(rw.canon()) >> (64 - SBB))];
                }
            }
        }
        // 3. bucket-major, thread-minor write offsets
        std::vector<size_t> boff(NSB + 1, 0);
        for (size_t b = 0; b < NSB; ++b) {
            size_t sum = 0;
            for (int t = 0; t < T; ++t) sum += hist[(size_t)t][b];
            boff[b + 1] = sum;
        }
        for (size_t b = 0; b < NSB; ++b) boff[b + 1] += boff[b];
        std::vector<std::vector<size_t>> wpos((size_t)T, std::vector<size_t>(NSB, 0));
        for (size_t b = 0; b < NSB; ++b) {
            size_t cur = boff[b];
            for (int t = 0; t < T; ++t) { wpos[(size_t)t][b] = cur; cur += hist[(size_t)t][b]; }
        }
        sidx.resize(total);
        hp_hint((void*)sidx.data(), sidx.size() * sizeof(SIdx), "seed index");
        // 4. scatter
        #pragma omp parallel for schedule(static, 1)
        for (int t = 0; t < T; ++t) {
            std::vector<size_t>& w = wpos[(size_t)t];
            for (size_t ci = tlo[(size_t)t]; ci < tlo[(size_t)t + 1]; ++ci) {
                const std::string& c = S.contigs[ci];
                Roll25 rw; const size_t cn2 = c.size();
                for (size_t b2 = 0; b2 < cn2; ++b2) {
                    rw.push(c[b2]);
                    if (b2 + 1 < (size_t)K || !rw.ok()) continue;
                    const uint64_t kk = rw.canon();
                    sidx[w[(size_t)(FlatKmerSet::mix(kk) >> (64 - SBB))]++] =
                        SIdx{ kk, (uint32_t)ci, (uint32_t)(b2 + 1 - (size_t)K) };
                }
            }
        }
        _bslap("idx: flat fill");
        // 5. sort inside each bucket, in place
        #pragma omp parallel for schedule(dynamic, 8)
        for (long long b = 0; b < (long long)NSB; ++b) {
            std::sort(sidx.begin() + (long)boff[(size_t)b],
                      sidx.begin() + (long)boff[(size_t)b + 1],
                      [](const SIdx& x, const SIdx& y) {
                          const uint64_t hx = FlatKmerSet::mix(x.key), hy = FlatKmerSet::mix(y.key);
                          if (hx != hy)      return hx < hy;
                          if (x.key != y.key) return x.key < y.key;
                          if (x.ci  != y.ci)  return x.ci  < y.ci;
                          return x.pos < y.pos;
                      });
        }
        _bslap("idx: sort");
        // 6. 64-cap, compacting in place inside each bucket (a key never
        //    straddles buckets: equal keys have equal mix), then one
        //    left-moving pass in ascending bucket order to close the gaps.
        std::vector<size_t> kept(NSB, 0);
        #pragma omp parallel for schedule(dynamic, 8)
        for (long long b = 0; b < (long long)NSB; ++b) {
            const size_t lo0 = boff[(size_t)b], hi0 = boff[(size_t)b + 1];
            size_t w = lo0, i = lo0;
            while (i < hi0) {
                size_t j = i; uint32_t k0 = 0;
                while (j < hi0 && sidx[j].key == sidx[i].key) {
                    if (k0 < 64) { sidx[w++] = sidx[j]; ++k0; }
                    ++j;
                }
                i = j;
            }
            kept[(size_t)b] = w - lo0;
        }
        size_t out = 0;
        for (size_t b = 0; b < NSB; ++b) {
            if (kept[b] && out != boff[b])
                memmove(&sidx[out], &sidx[boff[b]], kept[b] * sizeof(SIdx));
            out += kept[b];
        }
        sidx.resize(out);
        // 7. directory over the final array
        sdir.assign(((size_t)1 << BDIR) + 1, 0);
        for (size_t i = 0; i < sidx.size(); ++i)
            ++sdir[(size_t)(FlatKmerSet::mix(sidx[i].key) >> (64 - BDIR)) + 1];
        for (size_t i = 0; i < ((size_t)1 << BDIR); ++i) sdir[i + 1] += sdir[i];
    }
    _bslap("idx: cap+emit");
    // lookup(k): [first,last) range in sidx matching canonical key k, or (0,0).
    auto lookup_range = [&](uint64_t k) -> std::pair<size_t,size_t> {
        const size_t h = (size_t)(FlatKmerSet::mix(k) >> (64 - BDIR));
        size_t i = sdir[h]; const size_t e = sdir[h + 1];
        while (i < e && sidx[i].key != k) ++i;
        if (i == e) return {0, 0};
        size_t j = i + 1;
        while (j < e && sidx[j].key == k) ++j;
        return {i, j};
    };
    // idx (a hash map) serves only gapped_indel_scan, which is opt-in (CAPS_GAPSCAN),
    // so it is built only then.
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t,uint32_t>>> idx;
    if (std::getenv("CAPS_GAPSCAN")) {
        idx.reserve(1u << 21);
        for (uint32_t ci = 0; ci < S.contigs.size(); ++ci) {
            const std::string& c = S.contigs[ci];
            for (size_t i = 0; i + (size_t)K <= c.size(); ++i) {
                uint64_t v; if (!pack25(c.data() + i, v)) continue;
                uint64_t rcv = rc25(v), can = v < rcv ? v : rcv;
                auto& vec = idx[can];
                if (vec.size() < 64) vec.push_back({ci, (uint32_t)i});
            }
        }
    }

    // 3. place every read, both strands, fewest mismatches wins
    //
    // PARALLEL. This placement loop is the caller's largest phase and runs twice
    // (here for the pileup substrate, again for the bubble substrate inside
    // indel_pass). It is safe to parallelise: iteration `o` writes ONLY
    // S.read_cid[o] / read_pos[o] / read_rc[o] / read_clip[o], which are disjoint
    // across iterations and pre-sized above; `idx` and `S.contigs` are read-only
    // once step 2 has finished (bound through a const ref below so a stray insert
    // cannot compile); every best_* accumulator is declared inside the body and is
    // private. `placed`/`improved` are the only shared scalars and are a reduction.
    // The output does not depend on order, so it is byte-identical to a serial
    // run. (OpenMP over CONTIGS in the pileup is slower, because that loop has few,
    // very unequal items; this one has millions of uniform items and dynamic
    // scheduling covers the tail.)
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
        // ── STOP AT A PROVABLY OPTIMAL PLACEMENT ────────────────────────────
        //
        // score = ov_hi - 6*mm, and ov_hi <= rl, so rl is the MAXIMUM any
        // candidate can score -- reached only by a full-length, zero-mismatch,
        // zero-clip placement. Once one is in hand no later candidate can beat
        // it, and the update below is a STRICT `>`, so none could displace it
        // even if evaluated. Abandoning the rest of the seed offsets AND the
        // other strand is therefore exact, not a heuristic cutoff.
        //
        // This is the cheap way to cut the candidate count. The obvious way --
        // deduplicating repeated (contig, start) proposals -- was implemented
        // and MEASURED: it skipped 723 MILLION repeats at full chr20 and made
        // placement slightly SLOWER, because branch-and-bound already kills a
        // bad candidate in 19-30 bytes of sequential compare on a cache line
        // that is already loaded, while the dedup probe is a fresh random miss.
        // Trading one cache miss for another buys nothing; not evaluating the
        // candidate at all is what helps.
        bool perfect = false;
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
                const auto rng = lookup_range(can);
                for (size_t ri = rng.first; ri < rng.second; ++ri) {
                    const std::pair<uint32_t,uint32_t> pr{sidx[ri].ci, sidx[ri].pos};
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
                        // BRANCH-AND-BOUND, exact.
                        //
                        // The final score is ov_hi - 6*mm and mm only ever
                        // grows, so ov_hi - 6*mm_so_far is a valid UPPER BOUND
                        // on what this candidate can still reach. Once that
                        // bound is <= best_score the candidate provably cannot
                        // satisfy the strict `score > best_score` test below,
                        // so abandoning it changes nothing that is kept.
                        //
                        // This is a pruning rule, not a heuristic cutoff: it
                        // is exact for ANY arrival order of candidates, so the
                        // surviving best is bit-identical to the exhaustive
                        // scan. Its force comes from the 5x mismatch penalty
                        // already in the score -- once a perfect ~148 bp
                        // placement is in hand, a rival of the same length is
                        // dead at its FIRST mismatch.
                        // WORD-WISE SKIP over matching stretches.
                        //
                        // The scan compared one byte at a time through a b2i()
                        // call. Where eight bytes are EQUAL they contribute
                        // zero mismatches no matter what they are -- including
                        // N, since the test is `read base differs from contig
                        // base`, which equal bytes never satisfy. So an equal
                        // 8-byte word can be skipped outright, and only a word
                        // containing a difference needs the per-byte path.
                        //
                        // Mismatches are still found strictly left to right, so
                        // the sequence of mm values, and therefore the exact
                        // point the bound below fires, is unchanged. Bounds are
                        // safe by construction: j + 8 <= ov_hi <= rl - clip and
                        // <= c.size() - cst, so both loads stay inside their
                        // strings.
                        //
                        // This matters because branch-and-bound already kills
                        // bad candidates in a few bytes; what is left is the
                        // long, nearly-perfect alignments, which are exactly
                        // the ones made of matching words.
                        int mm = 0; int64_t j = 0; bool dead = false;
                        const char* rp = r.data() + (size_t)clip;
                        const char* cp = c.data() + (size_t)cst;
                        while (j + 8 <= ov_hi) {
                            uint64_t xw, yw;
                            memcpy(&xw, rp + j, 8); memcpy(&yw, cp + j, 8);
                            if (xw == yw) { j += 8; continue; }
                            const int64_t e8 = j + 8;
                            for (; j < e8; ++j) {
                                const char a = rp[j];
                                if (b2i(a) >= 0 && cp[j] != a) {
                                    ++mm;
                                    if ((long)ov_hi - 6L * (long)mm <= best_score) { dead = true; break; }
                                }
                            }
                            if (dead) break;
                        }
                        if (!dead) for (; j < ov_hi; ++j) {
                            const char a = rp[j];
                            if (b2i(a) >= 0 && cp[j] != a) {
                                ++mm;
                                if ((long)ov_hi - 6L * (long)mm <= best_score) { dead = true; break; }
                            }
                        }
                        if (dead) continue;               // provably cannot win
                        const long score = (long)(ov_hi - mm) - 5L * (long)mm;
                        if (score > best_score) { best_score = score; best_mm = mm;
                                                  best_c = pr.first; best_p = (uint32_t)cst;
                                                  best_rc = (uint8_t)strand; best_clip = (uint16_t)clip;
                                                  if (score == (long)rl) { perfect = true; break; } }
                    }
                    if (perfect) break;
                }
                if (perfect) break;
            }
            if (perfect) break;
        }
        if (best_c != UINT32_MAX && best_mm < 7) {          // same MAPQ<20 gate the pileup uses
            if (S.read_cid[o] == UINT32_MAX) ++placed; else ++improved;
            S.read_cid[o] = best_c; S.read_pos[o] = best_p; S.read_rc[o] = best_rc;
            S.read_clip[o] = best_clip;
        }
    }
    _bslap("place reads");
    // OPT-IN (CAPS_GAPSCAN=1), not used by any reported configuration. This substrate
    // keeps the two haplotypes as separate contigs, so a read carrying an indel is
    // placed on its own haplotype's contig with no gap; the indel evidence is in the
    // contig PAIR, not in a read-versus-contig gap. A read-level gapped channel can only
    // add recall to a caller built on one merged consensus. It found 0-1 events per
    // test window here.
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
    // Per-phase wall time and RSS, printed to stderr. Diagnostic only; results do not
    // depend on it.
    using clk = std::chrono::steady_clock;
    auto t_start = clk::now();
    auto elapsed_s = [](clk::time_point a, clk::time_point b){
        return std::chrono::duration<double>(b - a).count(); };
    // Reads VmRSS from /proc/self/status (Linux only; elsewhere it reports -1).
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
        // ── READ PLACEMENTS, CARRIED THROUGH ─────────────────────────────────
        // cd_in holds, for every original read, its contig and offset, derived by the
        // encoder from the pseudogenome position it computed in order to compress.
        // k-mer counts cannot reconstruct this provenance. Linked de Bruijn graphs
        // (McCortex, Turner et al. 2018) and colored graphs (LueVari) use the same
        // idea of read provenance to resolve repeats; here compression has already
        // computed it, so carrying it costs three vector copies.
        cd.read_cid  = cd_in.read_cid;
        cd.read_pos  = cd_in.read_pos;
        cd.read_rc   = cd_in.read_rc;
        cd.read_clip = cd_in.read_clip;
        cd.read_ppos = cd_in.read_ppos;
        cd.valid = true;
        fprintf(stderr, "[PLACE] carried %zu read placements from the encoder\n",
                cd.read_cid.size());
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


    // ── READ SEED INDEX -- DISABLED BY DEFAULT ─────────────────────────────
    // Read-substring support test after eBWT2SNP: an emitted fragment must occur
    // as a substring of at least C reads (Hamming distance <= 2), which rejects
    // chimeric junctions whose k-mers each occur in some read while no single read
    // contains the whole fragment. Here every candidate is already built from reads,
    // so the test has nothing to reject: on the five chr20 tuning windows the calls
    // were byte-identical with and without it, and the index is the largest single
    // structure. It is therefore off by default. CAPS_FORCE_READSUB=1 builds it.
    // Literature figures and the measurement: docs/extras/ENGINEERING_NOTES_CALLER.md.
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
        // rejected ~80% of true matches (measured: indel F1 fell to
        // 0.344 with stride-5 queries).
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
    // CHUNKED COUNTING. Collecting every k-mer occurrence and sorting once would
    // need memory proportional to total occurrences (about 12 GB at full chr20).
    // Instead BATCH_KMERS occurrences at a time are sorted and run-length counted
    // into a run, the raw batch is discarded, and the runs are k-way merged with
    // counts summed per key. Peak extra memory is O(batch + distinct k-mers) and the
    // (kmer, count) set is the same. The batch size is divided by the thread count so
    // total batch memory does not grow with threads.
    size_t BATCH_KMERS = 32u << 20;       // ~256 MB of raw uint64 in total
    {
        int nth = 1;
        #ifdef _OPENMP
        nth = omp_get_max_threads();
        #endif
        if (nth > 1) BATCH_KMERS = std::max<size_t>(1u << 20, BATCH_KMERS / (size_t)nth);
    }
    // 12 BYTES, NOT 16.
    //
    // {uint64 kmer; uint32 cnt} is 16 bytes because the compiler pads it to
    // the alignment of its widest member -- 4 bytes of nothing, on the single
    // largest array in the run. At full chr20 kc holds 140,719,632 entries, so
    // the padding alone is 537 MB.
    //
    // Packing is safe here: the k-mer occupies 62 bits so nothing is truncated,
    // the spill files are written and read by the same binary in one run so the
    // on-disk record size stays self-consistent, and every offset into those
    // files is already computed as `index * sizeof(KC)`. x86 unaligned loads
    // cost nothing measurable next to the cache misses this removes.
    #pragma pack(push,1)
    struct KC { uint64_t kmer; uint32_t cnt; };
    #pragma pack(pop)
    static_assert(sizeof(KC) == 12, "kc record must stay packed at 12 B");
    // [KC-SPLIT] stderr timers for the spill, partition sort and merge steps of
    // kc_H_build. Diagnostic only.
    using _kcclk = std::chrono::steady_clock;
    auto _kcl = _kcclk::now();
    auto _kclap = [&](const char* nm){ auto now=_kcclk::now();
        fprintf(stderr, "[KC-SPLIT] %-20s %7.2fs  rss=%ldMB peak=%ldMB\n", nm,
                std::chrono::duration<double>(now-_kcl).count(),
                caps_rss_mb(), caps_hwm_mb()); _kcl=now; };
    // ── HUGE PAGES FOR kc ───────────────────────────────────────
    // kc is a multi-GB array that the bubble traversal binary-searches, so that
    // walk is TLB-bound. With transparent_hugepage=[madvise] a region gets huge
    // pages only if it asks. The hint must follow reserve() and precede the fill
    // (madvise does not affect pages already faulted in), and all three reserve()
    // sites (in-RAM, spill/superkmer, k-way merge) carry it.
    // Dropping cnt==1 nodes from kc does not shrink it at full chr20 scale: the
    // superkmer counting path has already summed duplicate keys, so almost no
    // singletons remain (the singleton share quoted in the MINC comment is a
    // window-scale figure). Measure the cnt==1 fraction at the scale being run
    // before revisiting this.
    auto kc_hp_hint = [](std::vector<KC>& v){
        void* base = (void*)v.data();
        size_t len = v.capacity() * sizeof(KC);
        if (!base || len < (2u << 20)) return;
        uintptr_t a0 = ((uintptr_t)base + (2u<<20) - 1) & ~(uintptr_t)((2u<<20) - 1);
        size_t off = a0 - (uintptr_t)base;
        if (len <= off) return;
        size_t l2 = (len - off) & ~(size_t)((2u<<20) - 1);
        if (!l2) return;
#ifdef __linux__
        if (madvise((void*)a0, l2, MADV_HUGEPAGE) == 0)
            fprintf(stderr, "[KC-HP] MADV_HUGEPAGE on %.2f GB of kc\n", (double)l2/1073741824.0);
        else
            fprintf(stderr, "[KC-HP] madvise failed (errno %d) -- kc stays on 4K pages\n", errno);
#else
        (void)a0; (void)l2;
#endif
    };
    // Reads arrive 2-bit packed in graph-only mode (see src/encoder.cpp
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
    // ── LAYER 3: DECLARED MEMORY CEILING ──────────────────────────────────
    // As in GATB and KMC (`-max-memory`), the counter declares a memory budget,
    // estimates the in-RAM counting cost from the input, and turns the disk spill on
    // by itself when the estimate exceeds the budget. Without this a full-chr20 run
    // started without the spill variable counted entirely in RAM at about 4x the
    // memory. Spilling changes where k-mers are counted, never which ones exist:
    // the merge sums counts per key either way, which the k-mer-count gate checks.
    auto mem_available_mb = []() -> size_t {
        FILE* f = fopen("/proc/meminfo", "r");
        if (!f) return 0;
        char key[64]; unsigned long val; char unit[16]; size_t avail = 0;
        while (fscanf(f, "%63s %lu %15s", key, &val, unit) >= 2) {
            if (!strcmp(key, "MemAvailable:")) { avail = val / 1024; break; }
        }
        fclose(f);
        return avail;
    };
    const size_t MEM_AVAIL_MB = mem_available_mb();
    // ── AN ABSOLUTE DEFAULT BUDGET, NOT A FRACTION OF WHAT IS FREE ──────────
    // A budget of 60% of MemAvailable disabled the spill on large machines (the
    // projection fit, so counting ran in RAM at several times the memory). The
    // default is therefore an absolute 5000 MB, capped by system memory, which is
    // GATB's rule (ConfigurationAlgorithm.cpp:334-343):
    //
    //     if (_max_memory == 0) _max_memory = 5000;                  // MB
    //     if (_max_memory > (system_mem*2)/3) _max_memory = (system_mem*2)/3;
    //
    // An explicit CAPS_MAXRAM_MB overrides it. This changes no output, for the
    // reason given above.
    size_t mem_ceil_default = 5000;
    if (MEM_AVAIL_MB && mem_ceil_default > (MEM_AVAIL_MB * 2) / 3)
        mem_ceil_default = (MEM_AVAIL_MB * 2) / 3;
    const size_t MEM_CEIL_MB  = std::getenv("CAPS_MAXRAM_MB")
                              ? (size_t)atoll(std::getenv("CAPS_MAXRAM_MB"))
                              : mem_ceil_default;
    // Projected in-RAM counting cost: every k-mer of every read becomes one
    // 16-byte KC record in the run buffers before RLE, which is what actually
    // sets the peak. Reads and lengths are already resident and are not part
    // of the delta being budgeted here.
    // Length is read from the packed header without unpacking the bases.
    auto packed_len = [&](const std::string& src) -> size_t {
        if (!SEQ_PACKED) return src.size();
        if (src.empty()) return 0;
        if ((uint8_t)src[0] == 1) return src.size() - 1;      // raw (had non-ACGT)
        if (src.size() < 3) return 0;
        return (size_t)(uint8_t)src[1] | ((size_t)(uint8_t)src[2] << 8);
    };
    size_t proj_kmers = 0;
    for (size_t i = 0; i < seqs.size(); ++i) {
        const size_t L = packed_len(seqs[i]);
        if (L >= 31) proj_kmers += L - 30;
    }
    const size_t PROJ_MB = (proj_kmers * sizeof(KC)) / (1024 * 1024);

    const char* SPILLDIR_ENV = std::getenv("CAPS_KC_SPILLDIR");
    // Auto-enable when the projection exceeds the ceiling and a target
    // directory is derivable. Explicit settings always win over the estimate.
    static std::string spill_auto;
    const char* SPILLDIR = SPILLDIR_ENV;
    bool spill_forced_off = std::getenv("CAPS_KC_NOSPILL") != nullptr;
    if (!SPILLDIR && !spill_forced_off && MEM_CEIL_MB && PROJ_MB > MEM_CEIL_MB) {
        const char* base = std::getenv("TMPDIR");
        if (!base) base = "/tmp";
        spill_auto = std::string(base) + "/caps_kc_spill_" + std::to_string((long)getpid());
        if (mkdir(spill_auto.c_str(), 0700) == 0 || errno == EEXIST)
            SPILLDIR = spill_auto.c_str();
    }
    fprintf(stderr, "[KC-BUDGET] avail=%zuMB ceiling=%zuMB projected=%zuMB -> %s\n",
            MEM_AVAIL_MB, MEM_CEIL_MB, PROJ_MB,
            SPILLDIR ? (SPILLDIR_ENV ? "spill (explicit)" : "spill (auto)") : "in-RAM");
    const bool  SPILL    = (SPILLDIR != nullptr);
    // 2^SPILL_BITS key-range partitions. 256 partitions keeps each one small
    // enough to sort in RAM while staying far under any open-file limit.
    // Range-checked like CAPS_PLOIDY above: SPILL_BITS feeds `1u<<SPILL_BITS`
    // (partition table size) and `62-SPILL_BITS` (a k-mer bit-shift) a few
    // lines below and throughout this block; without the check one bad
    // value would give an undefined-behaviour shift (>=62) or an
    // out-of-memory partition table (SPILL_BITS in the high 20s+).
    int SPILL_BITS = 8;
    if (const char* sb = std::getenv("CAPS_KC_SPILLBITS")) {
        int v = atoi(sb); if (v >= 1 && v <= 20) SPILL_BITS = v;
    }
    const size_t PBUF_MAX  = 1u << 16;      // entries buffered per partition per thread
    // Superkmer spill: minimizer-partitioned, 2-bit packed. Verified at
    // 1.161 bytes/k-mer vs 12 for records, k-mer multiset identical.
    // SUPERKMER IS THE DEFAULT WHENEVER THE COUNTER SPILLS, FOR CORRECTNESS.
    // The key-only spill path loses part of the k-mer set (HG002 tuning window,
    // spill forced: 741,350 kc nodes against 1,063,607 in RAM, SNV F1 0.000 against
    // 0.886). The superkmer path reproduces the in-RAM multiset exactly.
    // CAPS_KC_NOSUPERKMER=1 selects the key-only path and exists only to reproduce
    // that measurement; do not use it for results.
    const bool  SUPERK = SPILL && std::getenv("CAPS_KC_NOSUPERKMER") == nullptr;

    // ── LAYER 2: FREQUENCY-RANKED MINIMIZERS ────────────────────────────────
    // Frequency-ranked minimizers, as in KMC2 and GATB: order m-mers by sampled
    // frequency (ties broken by the m-mer, so the order is deterministic) instead of
    // lexicographically.
    //
    // THIS CANNOT CHANGE THE K-MER SET: a superkmer is a container, and the merge
    // sums counts per key whatever the grouping.
    //
    // OPT-IN (CAPS_KC_FREQMIN=1), OFF by default. It does not reduce spill volume:
    // on the HG002 tuning window it wrote 14% more bytes per k-mer, because rare
    // minimizers change more often and give shorter runs. KMC2 and GATB use it for
    // partition balance, which would matter only if the merge became memory-bound
    // per partition. Measurements: docs/extras/ENGINEERING_NOTES_CALLER.md.
    const bool FREQMIN = SUPERK && std::getenv("CAPS_KC_FREQMIN") != nullptr;
    std::vector<uint32_t> mm_freq;
    size_t sk_bytes_total = 0;   // superkmer spill volume, for bytes/k-mer
    if (FREQMIN) {
        const int MM_ = 10;
        const uint32_t MMASK_ = (1u << (2*MM_)) - 1;
        mm_freq.assign((size_t)1 << (2*MM_), 0u);
        const size_t STRIDE = std::getenv("CAPS_KC_FREQSTRIDE")
                            ? (size_t)atoll(std::getenv("CAPS_KC_FREQSTRIDE")) : 32;
        std::string fbuf;
        for (size_t si = 0; si < seqs.size(); si += STRIDE) {
            const std::string& s = unpack_read(seqs[si], fbuf);
            uint32_t mv = 0; int have = 0;
            for (size_t i = 0; i < s.size(); ++i) {
                const int b = b2i(s[i]);
                if (b < 0) { have = 0; continue; }
                mv = ((mv << 2) | (uint32_t)b) & MMASK_;
                if (++have >= MM_) ++mm_freq[mv];
            }
        }
        fprintf(stderr, "[KC-FREQMIN] sampled every %zu reads, %zu-mer ranking table built\n",
                STRIDE, (size_t)MM_);
    }
    // (definition for RunSrc::RBUF_ lives at namespace scope below)
    {
        // PARALLEL. Each thread collects and RLE-compresses its own batches
        // into its own runs; the k-way merge below already sums counts for a
        // key across every run, so partitioning the reads across threads
        // changes only WHICH run holds a key, never the summed count. The
        // (kmer,count) set is therefore identical to the serial version by
        // construction -- this is not an approximation.
        // A failed spill write (disk full, quota, permissions) is recorded in a flag
        // rather than aborting inside the parallel region, and checked once after the
        // region joins, so a short spill is an error instead of silently dropped
        // k-mer occurrences.
        std::atomic<bool> spill_write_ok{true};
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
                // Partition on the TOP BITS of the canonical k-mer, so every partition owns
                // a disjoint key range and partitions are processed one at a time (GATB
                // partitions by minimizer for the same reason, SortingCountAlgorithm.cpp:804).
                // Partition order is key order, so concatenating the sorted partitions gives a
                // globally sorted kc with no merge. (GATB needs minimizers because superkmers
                // must keep adjacent k-mers together; kc stores plain k-mers.)
                if (SPILL) {
                    // SPILL ONLY THE KEY: the count is recovered when the partition is read back
                    // and sorted (one entry per occurrence), so writing it would waste a third
                    // of the spill.
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
                            else { spill_write_ok.store(false, std::memory_order_relaxed); }
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
            // Packing and reverse-complementing each k-mer from scratch is O(k) per
            // position, at about 1.4 billion positions for full chr20.
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
            // occupy L+k-1 bases, packed 4 bases/byte.
            // The format was verified at 1.161 bytes per k-mer, against 12 for
            // (kmer, count) records, with an identical k-mer multiset.
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
                    else { spill_write_ok.store(false, std::memory_order_relaxed); }
                    #pragma omp atomic
                    sk_bytes_total += skbuf[part].size();
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
                        // Ordering key: frequency first, m-mer as a
                        // deterministic tie-break. Falls back to the plain
                        // m-mer when frequency ranking is off.
                        auto mkey = [&](uint32_t m) -> uint64_t {
                            return FREQMIN ? (((uint64_t)mm_freq[m] << 20) | m)
                                           : (uint64_t)m;
                        };
                        std::deque<std::pair<uint64_t,size_t>> dq;
                        uint32_t mv = 0; int have = 0;
                        // lastN is the most recent invalid base. A k-mer is only
                        // valid if its whole 31-base window starts after it --
                        // without this the deque refills 10 bases past an N and
                        // assigns minimizers to k-mers that still contain it,
                        // which would drop k-mers next to an N (their superkmers are
                        // rejected whole).
                        long long lastN = -1;
                        for (size_t i = 0; i < s.size(); ++i) {
                            const int b = b2i(s[i]);
                            if (b < 0) { dq.clear(); have = 0; lastN = (long long)i; continue; }
                            mv = ((mv << 2) | (uint32_t)b) & MMASK;
                            if (++have < MM) continue;
                            const size_t mpos = i - MM + 1;
                            const uint64_t kv = mkey(mv);
                            while (!dq.empty() && dq.back().first >= kv) dq.pop_back();
                            dq.push_back({kv, mpos});
                            if (i + 1 >= 31) {
                                const size_t kpos = i + 1 - 31;
                                if ((long long)kpos <= lastN) continue;   // window still spans an N
                                while (!dq.empty() && dq.front().second < kpos) dq.pop_front();
                                if (!dq.empty())
                                    mini[kpos] = FREQMIN
                                               ? (uint32_t)(dq.front().first & 0xFFFFFu)
                                               : (uint32_t)dq.front().first;
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
                    else { spill_write_ok.store(false, std::memory_order_relaxed); }
                    std::vector<uint64_t>().swap(pbuf[part]);
                }
            }
            #pragma omp critical(kcruns)
            for (auto& r : myruns) kc_runs.push_back(std::move(r));
        }
        if (!spill_write_ok.load(std::memory_order_relaxed)) {
            fprintf(stderr, "caps_caller: FATAL: a k-mer spill file write failed "
                            "(disk full, quota, or permissions under %s) -- "
                            "refusing to call variants against an incomplete "
                            "k-mer set rather than silently under-counting\n",
                            SPILLDIR);
            return -1;
        }
    }
    // k-way merge: at each step pick the smallest head key across all runs
    // and sum every run's count for that key. Runs are few (total k-mer
    // occurrences / BATCH_KMERS -- roughly 47 at full chr20 scale), so a
    // linear scan per output element is simple and correct; a heap would be
    // faster but is not needed at this run count.
    // HEAP MERGE, TWO PASSES. The cost is driven by the output length, not the run
    // count, so a linear scan of all run heads per output key is too slow at full
    // chr20; a min-heap over run heads makes the merge O(output x log runs).
    // Pass 1 counts distinct keys so that pass 2 allocates kc exactly once (reserving
    // the summed run lengths would over-reserve several-fold, because every run holds
    // its own copy of each common key). Both passes visit keys in the same ascending
    // order and sum the same per-run counts, so the output is unchanged.
    std::vector<KC> kc;
    _kclap("spill write");
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
        // The parallel read-back reproduces the serial k-mer set exactly (kc node
        // count checked on the HG002 tuning window). An earlier discrepancy was caused
        // by 2-bit packing mapping non-ACGT bases to 'A', not by the parallelism.
        sorted_parts.resize(1u << SPILL_BITS);
        // ── PARALLEL MERGE BY KEY RANGE ─────────────────────────────────────
        // A single serial priority-queue merge over all partition files was most of
        // kc_H_build's time.
        //
        // Partitions are keyed by minimizer, so they are not key-ranges and a
        // global ordering did have to come from somewhere. But each partition
        // file is written in ascending key order, so the records for any given
        // top-8-bit BIN occupy one contiguous slice of it. Recording those 256
        // slice boundaries per partition makes bin b's inputs addressable
        // directly, and bins partition the key space -- so the 256 bins can be
        // merged independently, in parallel, and concatenated in bin order to
        // get exactly the globally sorted array the serial merge produced.
        static constexpr int MB = 256, MSHIFT = 54;   // top 8 bits of a 62-bit k-mer
        std::vector<std::vector<uint32_t>> bin_off((size_t)1u << SPILL_BITS);
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
                // Records go out in ascending key order, so recording the
                // record index at which each top-8-bit BIN starts costs one
                // comparison per record and turns the global merge below into
                // 256 independent ones. (The bins are the same ones the radix
                // pre-binning above already sorted by, so they are contiguous.)
                std::vector<uint32_t>& bo = bin_off[pi];
                bo.assign(MB + 1, 0);
                std::vector<KC> outbuf; outbuf.reserve(kms.size() / 4 + 16);
                for (size_t i = 0; i < kms.size(); ) {
                    size_t j = i; while (j < kms.size() && kms[j] == kms[i]) ++j;
                    outbuf.push_back(KC{ kms[i], (uint32_t)(j - i) });
                    ++total_k; i = j;
                }
                for (const KC& e : outbuf) ++bo[(size_t)(e.kmer >> MSHIFT) + 1];
                for (int bq = 0; bq < MB; ++bq) bo[bq + 1] += bo[bq];
                // one bulk write instead of a 12-byte fwrite per record
                fwrite(outbuf.data(), sizeof(KC), outbuf.size(), sf);
                fclose(sf);
                sorted_parts[pi] = spath;
            }
        }
        _kclap("partition expand+sort");
        std::vector<uint32_t> live;                 // partitions that produced records
        for (uint32_t i = 0; i < sorted_parts.size(); ++i)
            if (!sorted_parts[i].empty()) live.push_back(i);
        // Merge each key-space bin independently. A bin's inputs are the
        // corresponding slice of every partition file; each slice is already
        // ascending, and after concatenating them one sort restores order --
        // the slices are ~n/256 records, so this sorts cache-resident pieces
        // rather than streaming 100M records through a heap.
        // TWO PASSES, EXACT ALLOCATION.
        //
        // On the graph-only path (the published Claim 2 configuration) kc sets the peak
        // memory, so the merged result must never be resident twice, and allocating kc
        // at the pre-dedup bound would hold the excess for the rest of the run (resize()
        // down does not release capacity).
        //
        // So: pass A sorts each bin and counts its distinct keys, retaining
        // nothing; that gives the EXACT final size, so kc is allocated once at
        // exactly the right size; pass B redoes the read and sort and writes
        // the deduped run straight into kc's slice. The price is one extra
        // read+sort of each bin (the bins are ~n/256 and cache-resident, and
        // the whole merge is ~1 s at full chr20); the gain is that the merged
        // array is never resident twice.
        //
        // Identical contents and order either way: bins partition the key
        // space and each is sorted by k-mer.
        auto load_bin = [&](int b2, std::vector<KC>& acc) {
            acc.clear();
            for (uint32_t pi : live) {
                const std::vector<uint32_t>& bo = bin_off[pi];
                if (bo.size() != (size_t)MB + 1) continue;
                const uint32_t lo = bo[(size_t)b2], hi = bo[(size_t)b2 + 1];
                if (hi <= lo) continue;
                FILE* f = fopen(sorted_parts[pi].c_str(), "rb");
                if (!f) continue;
                const size_t base = acc.size();
                acc.resize(base + (size_t)(hi - lo));
                if (fseek(f, (long)lo * (long)sizeof(KC), SEEK_SET) != 0 ||
                    fread(acc.data() + base, sizeof(KC), (size_t)(hi - lo), f) != (size_t)(hi - lo))
                    acc.resize(base);
                fclose(f);
            }
            if (!acc.empty())
                std::sort(acc.begin(), acc.end(),
                          [](const KC& x, const KC& y){ return x.kmer < y.kmer; });
        };
        std::vector<size_t> bkept((size_t)MB, 0);
        #pragma omp parallel
        {
            std::vector<KC> acc;
            #pragma omp for schedule(dynamic, 1)
            for (int b2 = 0; b2 < MB; ++b2) {
                load_bin(b2, acc);
                size_t d = 0;
                for (size_t i = 0; i < acc.size(); ) {
                    size_t j = i;
                    while (j < acc.size() && acc[j].kmer == acc[i].kmer) ++j;
                    ++d; i = j;
                }
                bkept[(size_t)b2] = d;
            }
        }
        std::vector<size_t> boff2((size_t)MB + 1, 0);
        for (int b2 = 0; b2 < MB; ++b2) boff2[(size_t)b2 + 1] = boff2[(size_t)b2] + bkept[(size_t)b2];
        // reserve -> hint -> resize, in that order: the hint must land on
        // memory that has not been faulted in yet, and resize() zero-fills.
        kc.reserve(boff2[(size_t)MB]);
        kc_hp_hint(kc);
        kc.resize(boff2[(size_t)MB]);
        #pragma omp parallel
        {
            std::vector<KC> acc;
            #pragma omp for schedule(dynamic, 1)
            for (int b2 = 0; b2 < MB; ++b2) {
                if (!bkept[(size_t)b2]) continue;
                load_bin(b2, acc);
                size_t o = boff2[(size_t)b2];
                for (size_t i = 0; i < acc.size(); ) {
                    size_t j = i; uint32_t sum = 0;
                    const uint64_t kk = acc[i].kmer;
                    while (j < acc.size() && acc[j].kmer == kk) { sum += acc[j].cnt; ++j; }
                    kc[o].kmer = kk; kc[o].cnt = sum; ++o;
                    i = j;
                }
            }
        }
        for (auto& sp : sorted_parts) if (!sp.empty()) ::remove(sp.c_str());
        _kclap("k-way merge");
        fprintf(stderr, "[KC-SUPERK] partitions=%zu distinct=%zu spill_bytes=%zu bytes/kmer=%.3f freqmin=%d\n",
                sorted_parts.size(), kc.size(), sk_bytes_total,
                kc.empty() ? 0.0 : (double)sk_bytes_total / (double)kc.size(), (int)FREQMIN);
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
            // Growing kc by doubling would hold old and new arrays together at the
            // last reallocation (several GB at full chr20).
            // Partitions are uniform slices of key
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
                kc_hp_hint(kc);
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
                    if (++idx[r] < kc_runs[r].size()) pq.push({(uint64_t)kc_runs[r][idx[r]].kmer, r});
                }
                if (fill_out) kc.push_back({key, sum});
                ++nout;
            }
            return nout;
        };
        kc.reserve(merge_pass(false));
        kc_hp_hint(kc);
        merge_pass(true);
    }
    kc_runs.clear(); kc_runs.shrink_to_fit();
    // ── DIRECTORY OVER kc ───────────────────────────────────────────────────
    //
    // kc_find was a lower_bound over the WHOLE array. At full chr20 that is
    // 140,719,632 entries, so ~27 dependent probes -- and the graph traversal
    // calls it for every successor of every solid node, twice per orientation,
    // over the whole array. It is the same shape already fixed twice in this
    // file (the seed index and the pcluster anchor map): a sorted array probed
    // at random is a latency problem, not a comparison problem.
    //
    // kc is sorted by RAW k-mer and is final here (nothing mutates it after
    // this point), so the top KDB bits of a key name a contiguous slice.
    // Recording each slice's start turns the search into one directory probe
    // plus a lower_bound over ~8 entries: about 4 touches instead of 27.
    //
    // The directory is keyed on the raw k-mer, NOT on a hash, because the array
    // is in raw order -- so buckets are as uneven as the k-mer distribution.
    // That is fine for a range narrowing: an over-full bucket just means a
    // slightly longer bounded search, never a wrong answer. 2^24 entries x 4 B
    // = 67 MB, independent of input size.
    //
    // EXACT: this narrows the search range and changes nothing else. A key is
    // found iff it was found before.
    static constexpr int KDB = 24;
    std::vector<uint32_t> kdir;
    {
        const int sh = 62 - KDB;
        kdir.assign(((size_t)1 << KDB) + 1, 0);
        for (size_t i = 0; i < kc.size(); ++i) ++kdir[(size_t)(kc[i].kmer >> sh) + 1];
        for (size_t i = 0; i < ((size_t)1 << KDB); ++i) kdir[i + 1] += kdir[i];
        fprintf(stderr, "[KC-DIR] %d-bit directory over %zu k-mers (%.1f MB, %.1f entries/bucket)\n",
                KDB, kc.size(), (double)kdir.size() * 4 / 1048576.0,
                kc.empty() ? 0.0 : (double)kc.size() / (double)((size_t)1 << KDB));
    }
    auto kc_find = [&](uint64_t key) -> const KC* {
        const size_t h = (size_t)(key >> (62 - KDB));
        const size_t lo = kdir[h], hi = kdir[h + 1];
        if (lo == hi) return nullptr;
        auto it = std::lower_bound(kc.begin() + (long)lo, kc.begin() + (long)hi, key,
            [](const KC& e, uint64_t k){ return e.kmer < k; });
        return (it != kc.begin() + (long)hi && it->kmer == key) ? &*it : nullptr;
    };
    uint32_t H = 30;
    // The k-mer count histogram's VALLEY -- the minimum between the error peak
    // and the true single-copy peak -- is the standard, dataset-measured
    // solid-k-mer threshold. It is computed here, with H, because the dBG channel
    // needs this quantity and must not use a fixed constant for it (see MINC).
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
        // H ESTIMATION bins the count histogram at FULL RESOLUTION (bw = 1) over a
        // bounded range. Bin width derived from the maximum count (set by repeats) would
        // put the single-copy peak in bucket 0, which the valley logic excludes, and
        // overestimate H several-fold; the depth guard (d > DHI*H*1.25) and the k-mer
        // sanity gate (> KHI*H) would then never fire. Counts above the range are repeats
        // and do not affect the single-copy peak; 2000 covers any realistic
        // per-haplotype depth.
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
            // (Histogram bimodality does not identify ploidy here: at 30x with k=31 the
            // located peak sits near the homozygous depth, so mass at 2H measures repeats.
            // The measured ratios for haploid samples fell both above and below HG002.
            // Ploidy is detected by the k-mer pair test below.)

        }
    }
    auto kcount = [&](const std::string& km) -> uint32_t {
        if (km.size() != 31) return 0;
        uint64_t v; if (!pack31(km.data(), v)) return 0;
        const KC* e = kc_find(canon31(v)); return e ? e->cnt : 0u;
    };
    // Same as kcount(std::string) without building a 31-character string per probe;
    // the junction scans in the indel emit path call it once per offset.
    auto kcount_at = [&](const char* p31, size_t avail) -> uint32_t {
        if (avail < 31) return 0;
        uint64_t v; if (!pack31(p31, v)) return 0;
        const KC* e = kc_find(canon31(v)); return e ? e->cnt : 0u;
    };

    // ── 2-6a. Per-contig pileup, candidates, flank pass, filters, SNV emit ──
    // CONTIG-BATCHED. Nothing in this stretch crosses a contig boundary: every key
    // into the pileup structures carries the contig id. So the pileup is built one
    // contig at a time and freed before the next, and peak memory is that of the
    // largest single contig's pileup rather than the whole chromosome's (tens of GB
    // at full chr20). Output order is unchanged: contigs are emitted in increasing
    // cid, cid is the primary sort key, and each contig's `kept` is sorted by pos.
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

    // Read once, before any thread starts: the override does not depend on the
    // contig, and writing the shared `H` inside the parallel loop would be a data
    // race.
    if (const char* he = std::getenv("CAPS_HAPLOID_COV")) { int v = atoi(he); if (v > 0) H = (uint32_t)v; }

    // WHAT IS AND IS NOT BATCHED. Only `col` (every read base at every covered
    // position) and `recs` are large, so only they are built and freed per contig.
    // The candidate columns `C` and the flank maps are small and are accumulated into
    // thread-local buffers, then merged, and the filter and emit stages run over the
    // complete global C exactly as in the unbatched form. This matters for
    // `med_cand_depth`, the depth guard's median, which is defined over ALL candidates
    // chromosome-wide and must not become a per-contig median. Keys are (cid,pos) and
    // each contig is handled by one thread, so the merge is a move with no count
    // reconciliation.
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

    // ── DBG BUBBLE CHANNEL — DiscoSNP++'s algorithm on the existing k-mer graph ──
    //
    // A de Bruijn graph is a set of k-mers plus the implicit edge X->Y whenever Y
    // (X shifted left one base plus one new base) is also in the set. `kc` is that
    // set (every canonical 31-mer of the reads, with its count) and `kc_find` is the
    // membership test, so DiscoSNP++'s bubble finder needs no new index, no new pass
    // over the reads and no second assembly. Each primitive maps onto kc_find:
    //
    //   their graph.successors(node)        -> 4 kc_find probes, keep cnt >= minc
    //   their branching node (>= 2 succ)    -> same test
    //   their graph.successors(node1,node2) -> for b in ACGT, BOTH extensions
    //                                          present (their getNodesCouple in
    //                                          GATB Graph.cpp:1602-1640 requires
    //                                          the SAME nt on both paths)
    //   their closure (succ.first==.second) -> node1 == node2
    //
    // The traversal reimplements DiscoSNP++ Bubble.cpp:300-347 (start: a branching
    // node, every successor PAIR) and Bubble.cpp:509-527 (expand: lockstep walk,
    // close on node equality). The closure test is the essential part: approaches
    // without it (chaining ties, which are exact-match; string-flank tests on
    // near-miss pairs) cannot see heterozygosity or lowered F1.
    //
    // Enabled by CAPS_DBG=1, which `capsule_decode call` sets.
    // min1/min2 = the WEAKEST k-mer count along each path. DiscoSNP++ obtains
    // the equivalent by running kissreads2 as a whole separate pass, mapping
    // every read back onto every bubble. Here that pass is not needed: a canonical
    // k-mer's count in kc is the number of reads containing it, so the minimum along
    // the path is per-path read support. The minimum rather than the mean, because
    // one unsupported k-mer means no read spans that part of the path.
    struct DbgBubble { std::string flank, path1, path2, lext, rext;
                       uint32_t cov1, cov2; uint32_t min1, min2; int len;
                       uint8_t from_sb = 0;      // 1 = produced by the superbubble walk
                       uint8_t multiallelic = 0; };  // >2 = allele count at a multi-allelic site
    std::vector<DbgBubble> dbg_bubbles;   // hoisted: also read by DBG-ONLY mode below
    if (std::getenv("CAPS_DBG") && (looks_diploid || std::getenv("CAPS_DBG_FORCE"))) {
        auto t_dbg = clk::now();
        // MINC -- ARCHITECTURAL, NOT TUNED.
        //
        // The floor is not the precision filter: every bubble is validated downstream
        // against the reads (the coherence pass), and in a two-stage design (permissive
        // generation, then evidence-based filtering) precision bought upstream costs recall
        // that no later stage can recover. On the tuning window, 60.4% of graph nodes have
        // count 1 (singleton errors), 3.9% have 2-9 and 35.6% have >= 10. Raising the
        // floor from 2 to 5 removes only 5% of nodes while deleting the 2-4 band where
        // low-coverage heterozygous alleles are. A floor derived from the histogram valley
        // also rises with depth, the wrong direction. MINC = 2 ("seen more than once")
        // discards the singletons and nothing else. DiscoSNP++ uses a permissive -c 3 for
        // the same reason, with kissreads2 supplying precision.
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
        // FOUR SEARCHES, INTERLEAVED, NOT SEQUENTIAL.
        //
        // The four successors are known before any lookup starts, and each
        // lookup is an independent ~27-hop binary search over a 2.1 GB table.
        // Run back to back, each hop's cache miss is waited out in full; run in
        // lockstep with a prefetch issued for all four midpoints before any of
        // them is read, the four latencies OVERLAP.
        //
        // This helps only because the loop is latency-bound rather than bandwidth-bound:
        // a model of the access pattern and the observed traversal time agree on that,
        // and a standalone benchmark on a table of the same shape ran 1.52x faster
        // interleaved.
        //
        // The loop count is fixed rather than data-dependent so all four lanes
        // step together; a lane whose range has collapsed simply stops
        // narrowing. Results are identical to four lower_bound calls by
        // construction: same invariant, same midpoint rule.
        auto succs = [&](uint64_t fwd, uint64_t out[4], uint32_t cnt[4]) -> int {
            uint64_t key[4]; uint64_t nxs[4];
            for (uint64_t b = 0; b < 4; ++b) {
                nxs[b] = ((fwd << 2) | b) & MASK31;
                key[b] = canon31(nxs[b]);
            }
            const KC* base = kc.data();
            const size_t N = kc.size();          // still the bound for the final check
            // START FROM THE DIRECTORY SLICE, NOT THE WHOLE ARRAY.
            //
            // This 4-way interleaved search is the traversal's successor lookup (kc_find is
            // not on this path).
            //
            // kdir already maps the top KDB bits of a k-mer to its contiguous
            // slice, so each key starts bounded to ~8 entries: about 3 steps
            // instead of 27. Exact -- the search is unchanged, only its initial
            // bounds are tightened to the slice the key must lie in.
            size_t lo[4], hi[4];
            for (int j = 0; j < 4; ++j) {
                const size_t hh = (size_t)(key[j] >> (62 - KDB));
                lo[j] = kdir[hh]; hi[j] = kdir[hh + 1];
            }
            for (int step = 0; step < 64; ++step) {
                bool any = false;
                size_t mid[4];
                for (int j = 0; j < 4; ++j) {
                    if (lo[j] < hi[j]) {
                        mid[j] = lo[j] + ((hi[j] - lo[j]) >> 1);
                        __builtin_prefetch(base + mid[j], 0, 1);
                        any = true;
                    }
                }
                if (!any) break;
                for (int j = 0; j < 4; ++j) {
                    if (lo[j] >= hi[j]) continue;
                    if (base[mid[j]].kmer < key[j]) lo[j] = mid[j] + 1;
                    else                            hi[j] = mid[j];
                }
            }
            int n = 0;
            for (int j = 0; j < 4; ++j) {
                if (lo[j] < N && base[lo[j]].kmer == key[j] && base[lo[j]].cnt >= MINC) {
                    out[n] = nxs[j]; cnt[n] = base[lo[j]].cnt; ++n;
                }
            }
            return n;
        };
        // Lockstep walk, the direct analogue of Bubble.cpp:509-527. Advances
        // both paths by the SAME base at each step and stops when they land on
        // the same node (bubble closed) or when no single common base extends
        // both (dead end).
        // Also accumulates the weakest k-mer count seen on each path (see the
        // DbgBubble comment), so per-path read support comes from this walk rather than
        // from a second pass.
        // Extension with BOUNDED RECURSION over ambiguous continuations.
        //
        // When the two paths have more than one common successor, the search recurses over
        // every successor pair, as DiscoSNP++ expand() does (Bubble.cpp:598, through
        // expand_heart). Stopping at the first ambiguity lost most bubble entrances
        // mid-walk.
        //
        // A node budget bounds the search: branching factor is at most 4 and
        // depth at most MAXEXT, so unbounded recursion is exponential in the
        // worst case. The budget is shared across the whole search for one
        // bubble, so a tangled region gives up quickly instead of exploding.
        const int WALK_BUDGET = std::getenv("CAPS_DBG_WALKBUDGET")
                              ? atoi(std::getenv("CAPS_DBG_WALKBUDGET")) : 0;
        const int MAXPOLY = std::getenv("CAPS_DBG_MAXPOLY")
                          ? atoi(std::getenv("CAPS_DBG_MAXPOLY")) : 1;
        // MAXPOLY = 1, MEASURED. SNV F1 on the tuning window (r2) and four held-out chr20
        // windows:
        //     window      MAXPOLY=1   MAXPOLY=2   MAXPOLY=3 (DiscoSNP++ -P default)
        //     r2 (tune)     0.886       0.881       0.881
        //     r3            0.892       0.886       0.879
        //     na            0.958       0.958       0.956
        //     r4            0.894       0.895       0.897
        //     r5            0.937       0.923       0.916
        //     held-out mean 0.9203      0.9155      0.9120
        // A higher cap adds recall (+6 TP on r2) but about twice as many false positives,
        // because a second polymorphism inside one walk is more often two error branches
        // than two real heterozygous sites. The structural reason: at k=31 and human
        // heterozygosity of about 1/1000, a second het site inside one 31 bp window is
        // rare, so larger values mostly admit paths through repeats. A genome with much
        // higher heterozygosity would justify a higher cap.
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
            // outright; on the HG002 tuning window 30% of missed truth SNVs ended here.
            // DiscoSNP++ continues: expand_heart is called with nb_polymorphism+1 over
            // successor pairs carrying DIFFERENT bases (Bubble.cpp:636), up to
            // max_polymorphism.
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
        // Without extension a bubble sequence is its 31 bp flank plus the path (mean
        // 63 bp on the tuning window), which bwa cannot place uniquely in chr20, so
        // bubbles were lost at the LIFT rather than at calling. DiscoSNP++ likewise
        // extends each bubble into its surrounding unitig before output
        // (expand_one_simple_path, its -t traversal).
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
        // At the missed indel sites on the tuning window the branching node has good
        // coverage (12-30) but exactly ONE successor even with no coverage floor: the
        // alternate haplotype's k-mers are absent, not below MINC. Most of these events
        // are homopolymer-length changes and the rest mostly tandem-repeat expansions.
        // Inserting a base into a run of that base makes the run longer without creating
        // a distinct k-mer path, so no bubble exists for any de Bruijn caller, DiscoSNP++
        // included. This is a representational limit of the graph, not a traversal
        // defect, and filters, gates or orderings cannot recover these events. They need
        // length-aware evidence such as read depth over the run.
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
        // 2013. A unitig walk that takes the earliest node all branches reach can only
        // find parallel simple paths; this is the full algorithm.
        //
        // The load-bearing condition is "a node enters the frontier only once
        // ALL of its parents have been visited". That is what allows nodes
        // reachable by several internal paths to be handled correctly.
        //   S      = frontier: discovered AND all parents visited
        //   seen   = discovered but not yet ready
        //   sink   = |S| == 1 and |seen| == 0  (only the exit remains)
        // Aborts on a tip (no successors) or on returning to the entrance.
        // The ab_* abort-reason counters are diagnostics.
        size_t ab_tip = 0, ab_cycle = 0, ab_budget = 0, ab_exhaust = 0, ab_ok = 0;
        // Independent indel discovery (their start_indel_prediction).
        // OPT-IN (CAPS_DBG_IBFS=1), off by default and not used by any reported result;
        // unmeasured at scale. IBFS_MAXINDEL mirrors DiscoSNP++'s D=100.
        const bool WANT_IBFS = std::getenv("CAPS_DBG_IBFS") != nullptr;
        const int IBFS_MAXINDEL = std::getenv("CAPS_DBG_IBFSMAX")
                                ? atoi(std::getenv("CAPS_DBG_IBFSMAX")) : 100;
        size_t n_ibfs_found = 0;
        auto find_sb = [&](uint64_t s0, uint64_t& t_out) -> bool {
            std::unordered_set<uint64_t> visited, seen;
            std::vector<uint64_t> S; S.push_back(s0);
            int guard = 0;
            while (!S.empty()) {
                if (++guard > SB_BUDGET) { _Pragma("omp atomic") ++ab_budget; return false; }
                const uint64_t v = S.back(); S.pop_back();
                visited.insert(v); seen.erase(v);
                uint64_t ch[4]; uint32_t cc[4];
                const int nc = succs_live(v, ch, cc);
                if (nc == 0) { _Pragma("omp atomic") ++ab_tip; return false; }   // tip
                for (int i = 0; i < nc; ++i) {
                    const uint64_t u = ch[i];
                    if (u == s0) { _Pragma("omp atomic") ++ab_cycle; return false; }  // cycle back to entrance
                    if (visited.count(u)) continue;
                    uint64_t pa[4];
                    const int np = preds_live(u, pa);
                    bool allv = true;
                    for (int j = 0; j < np; ++j)
                        if (!visited.count(pa[j])) { allv = false; break; }
                    if (allv) { S.push_back(u); seen.erase(u); }
                    else seen.insert(u);
                }
                if (S.size() == 1 && seen.empty()) { t_out = S[0]; _Pragma("omp atomic") ++ab_ok; return true; }
            }
            _Pragma("omp atomic") ++ab_exhaust;
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
        // Extending both sides is required for recall: without it graph-only recall was
        // 0.512 against DiscoSNP++'s 0.763 on the same k-mer set.
        // PARALLEL. Every node is examined independently: kc, MINC and the
        // walk lambdas are read-only, and the only shared writes are the
        // bubble list (accumulated per-thread, concatenated once) and the
        // branching counter (a reduction). Bubble ORDER changes, which is
        // harmless -- each bubble is self-contained and downstream code
        // addresses them by index into the same vector it dumps.
        // WANT_SB is left on the serial path: find_sb's abort counters are
        // shared diagnostics and superbubble mode is not the default.
        size_t branching = 0;
        // ── PARALLEL SUPERBUBBLE TRAVERSAL ─────────────────────────────────
        // Traversal is the caller's longest stage, so it runs in parallel in superbubble
        // mode too. Shared state:
        //   * dbg_bubbles  -- per-thread (`local_bubbles`), merged once
        //   * sb_sites     -- per-thread (`local_sb`), merged once
        //   * ab_* counters-- diagnostics only, atomic
        // Nothing else in the loop writes shared state; kc, MINC and the walk lambdas are
        // read-only. Bubble and site ORDER depend on thread scheduling, which is harmless:
        // every record is self-contained and the VCF is sorted downstream. The bubble
        // COUNT and F1 do not change.
        #pragma omp parallel
        {
        std::vector<DbgBubble> local_bubbles;
        std::vector<SbSite> local_sb;
        #pragma omp for schedule(dynamic, 2048) reduction(+:branching,n_ibfs_found) nowait
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
                // The exit is located with the superbubble algorithm (find_sb, the
                // Onodera/BubbleGun port), which respects the all-parents-visited condition,
                // and each allele's string is read off the branch trace that reaches it.
                // ── VARIANT GROUPS (SKA lo, MBE 2025) ───────────────────────
                // find_sb's exit test needs the frontier to drain to exactly one node with
                // nothing pending, which 3-4 branches through repeat-rich sequence almost never
                // satisfy, so multi-allelic sites fail there. SKA lo relaxes what counts as a
                // site: "a variant group is identified when AT LEAST TWO paths start from the
                // same entry node and end at the same exit node". The traces already hold every
                // node each branch visits, so the group is the earliest node reached by the most
                // branches (>= 2). This is tried only when the strict superbubble finds nothing,
                // so two-path bubbles and SNV calling are unaffected.
                //
                // On 30x diploid human data no group has more than two allele paths, because no
                // third haplotype exists to converge on; the extra successors are error or repeat
                // branches. Multi-allelic recovery therefore needs a sample with more than two
                // haplotypes (the tetraploid construction, or a pooled sample).
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
                    // Among candidate exits the nearest reconvergence is preferred, as in DiscoSNP++,
                    // whose indel BFS stops at the smallest closing event (Bubble.cpp:238-243), so
                    // long spurious events are not generated rather than filtered afterwards. In
                    // practice the ordering rarely matters: groups fire only when ns > 2, and
                    // find_sb already returns the earliest exit.
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
                    // PARTIAL CONVERGENCE. A 4-successor node where 3 branches converge and 1 dies
                    // in a tip or repeat is a 3-allele site, so the branches that do converge are
                    // kept instead of discarding the site. At least 2 convergent branches are
                    // required (one path is not a bubble). SNV calling is unaffected, because SNVs
                    // come from the pairwise lockstep walk. Multi-allelic loss on diploid data is
                    // inside find_sb (see above), not here.
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
                        local_sb.push_back(std::move(st));
                    }
                }
                // DO NOT skip the pairwise walk.
                //
                // The pairwise walk and the superbubble path are complementary:
                //   pairwise walk -> SNVs   (lockstep, same base both paths)
                //   superbubble   -> INDELS (paths of different length)
                // Both run at every branching node. Replacing the pairwise walk with the
                // superbubble path lowered SNV F1 (0.881 -> 0.791 on the tuning window), because
                // per-branch unitig traces are a coarser view of a two-path bubble. Superbubble
                // records are tagged so the emitter keeps only the different-length ones (see
                // SB_ONLY_INDEL). CAPS_DBG_SB_REPLACE=1 reproduces the substituting variant for
                // re-measurement only.
                if (std::getenv("CAPS_DBG_SB_REPLACE")) continue;
            }

            for (int i = 0; i < ns; ++i)
                for (int j = i + 1; j < ns; ++j) {
                    // No cheap pre-walk reject is applied: most failed walks fail to CLOSE, which
                    // only walking reveals, and a coverage-ratio pre-test changed neither run time
                    // nor results at full chr20.
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

            // ── INDEL DISCOVERY (DiscoSNP++ start_indel_prediction) ─────────
            // Reimplemented from the published DiscoSNP++ algorithm (kissnp2
            // Bubble.cpp:200-285). At every branching node and every successor PAIR,
            // DiscoSNP++ runs its SNP search and its indel search INDEPENDENTLY
            // (Bubble.cpp:310-346): the indel search needs neither the SNP search nor a
            // superbubble to succeed. Indels taken only from find_sb would be attempted at
            // the small share of nodes where find_sb succeeds (15.9% at full chr20).
            //
            // THE ALGORITHM. Extend ONE path while the other waits, until the extended path's
            // last base equals the base that would let the two reconverge, then close in
            // LOCKSTEP with walk(), which is DiscoSNP++'s expand(). With their default
            // `authorised_branching=0` the queue never holds more than one live node (the
            // search stops at the first junction), so this is a unitig extension, not a search.
            //
            // The shape test is applied inline, as for superbubble records: lockstep closure
            // guarantees the tail agrees, but the inserted block must be ONE contiguous event.
            if (WANT_IBFS) {
                const uint64_t begin_[2] = { 0, 0 };
                (void)begin_;
                for (int i = 0; i < ns; ++i)
                for (int j = i + 1; j < ns; ++j) {
                    const uint64_t bg[2] = { so[i], so[j] };
                    const uint32_t bc[2] = { sc[i], sc[j] };
                    int found_del = IBFS_MAXINDEL + 10;
                    for (int ext_id = 0; ext_id < 2; ++ext_id) {
                        const uint64_t other = bg[1 - ext_id];
                        const char end_ins = "ACGT"[other & 3ULL];
                        uint64_t cur = bg[ext_id];
                        std::string ins;                 // the inserted block
                        for (int step = 0; step < IBFS_MAXINDEL; ++step) {
                            uint64_t ch[4]; uint32_t cc[4];
                            const int nc = succs_live(cur, ch, cc);
                            // authorised_branching == 0: abandon at a junction
                            if (nc != 1) break;
                            cur = ch[0];
                            ins += "ACGT"[cur & 3ULL];
                            const int insert_size = (int)ins.size();
                            if (insert_size >= found_del) break;
                            if ("ACGT"[cur & 3ULL] != end_ins) continue;
                            // candidate closure: lockstep from here
                            std::string e1, e2;
                            uint32_t m1 = bc[ext_id], m2 = bc[1 - ext_id];
                            uint64_t closenode = 0;
                            if (!walk(cur, other, e1, e2, m1, m2, closenode)) continue;
                            std::string p_ext =
                                std::string(1, "ACGT"[bg[ext_id] & 3ULL]) + ins + e1;
                            std::string p_oth =
                                std::string(1, "ACGT"[other & 3ULL]) + e2;
                            if (p_ext.size() == p_oth.size()) continue;   // not an indel
                            // one contiguous event only (LCP + LCS test)
                            {
                                const std::string& sl =
                                    (p_ext.size() > p_oth.size()) ? p_ext : p_oth;
                                const std::string& ss =
                                    (p_ext.size() > p_oth.size()) ? p_oth : p_ext;
                                size_t lcp = 0;
                                while (lcp < ss.size() && ss[lcp] == sl[lcp]) ++lcp;
                                size_t lcs = 0;
                                while (lcs < ss.size() - lcp &&
                                       ss[ss.size()-1-lcs] == sl[sl.size()-1-lcs]) ++lcs;
                                if (lcp + lcs < ss.size()) continue;
                            }
                            DbgBubble bb2;
                            bb2.flank = k2s(node.kmer);
                            bb2.path1 = (ext_id == 0) ? p_ext : p_oth;
                            bb2.path2 = (ext_id == 0) ? p_oth : p_ext;
                            bb2.cov1 = sc[i];  bb2.cov2 = sc[j];
                            bb2.min1 = (ext_id == 0) ? m1 : m2;
                            bb2.min2 = (ext_id == 0) ? m2 : m1;
                            bb2.len  = (int)std::max(p_ext.size(), p_oth.size());
                            bb2.lext = ext_left(node.kmer, MAXFLANK);
                            bb2.rext = ext_right(closenode, MAXFLANK);
                            bb2.from_sb = 1;             // indel-shaped record
                            local_bubbles.push_back(std::move(bb2));
                            ++n_ibfs_found;
                            found_del = insert_size;     // prefer the shortest
                            break;
                        }
                    }
                }
            }
            }   // end orientation loop
        }
        #pragma omp critical(dbgbub)
        {
            for (auto& b : local_bubbles) dbg_bubbles.push_back(std::move(b));
            for (auto& t : local_sb)      sb_sites.push_back(std::move(t));
        }
        }   // end parallel
        // ── READ-COHERENCE: the fragment must exist in REAL READS ───────────
        // The false-positive mechanism this targets is a CHIMERIC JUNCTION: every k-mer
        // along a bubble path can be present in kc while the FULL path occurs in no single
        // read (each k-mer contributed by a different read, joined at a shared anchor).
        // k-mer counts cannot see this. The test, as in eBWT2SNP, is that the fragment
        // must be a substring of real reads. DiscoSNP++ obtains the equivalent from
        // kissreads2, a separate pass mapping every read onto every bubble. Here the
        // reads are already in memory (`seqs`), so bubbles are indexed by the 31-mer
        // spanning the variant, the reads are swept once, and the full fragment is
        // verified.
        if (!std::getenv("CAPS_DBG_NOCOH") && !dbg_bubbles.empty()) {
            auto t_coh = clk::now();
            // COHC -- required read support -- DERIVED FROM MEASURED COVERAGE.
            //
            // The number of reads that must contain a fragment scales with depth, so COHC is a
            // fraction of H, the measured haploid k-mer depth (a heterozygous allele is
            // covered by about H reads): H/10 with a floor of 2. At the benchmark depth (H=20)
            // this equals 2, the value swept on the tuning window. Unlike the error floor
            // MINC, which does not scale with depth, read support does.
            const int COHC = std::getenv("CAPS_DBG_COHC")
                           ? atoi(std::getenv("CAPS_DBG_COHC"))
                           : std::max(2, (int)(H / 10));
            // Minimum phred for a read's variant base to count as support.
            // 20 = 1% error probability, the standard confident-base cutoff.
            const int MINQ = std::getenv("CAPS_DBG_MINQ")
                           ? atoi(std::getenv("CAPS_DBG_MINQ")) : 20;
            // FRAGMENT HALF-WIDTH. A supporting read must contain the WHOLE fragment
            // (2*HALFW+1 bases), so HALFW sets how much unique context anchors each allele.
            // More context raises precision (repeat and paralog collapses stop finding
            // support) but loses true variants near read ends, so the useful ceiling depends
            // on read length. halfw_auto below is 0.30 x the median read length (floor 31,
            // because the probe must straddle a 31-mer junction): 44 on 148 bp reads.
            int halfw_auto = 31;
            {
                std::vector<size_t> lens;
                lens.reserve(std::min<size_t>(seqs.size(), 4096));
                const size_t stride = std::max<size_t>(1, seqs.size() / 4096);
                std::string lbuf;
                for (size_t i = 0; i < seqs.size(); i += stride)
                    lens.push_back(unpack_read(seqs[i], lbuf).size());
                if (!lens.empty()) {
                    std::nth_element(lens.begin(), lens.begin() + lens.size()/2, lens.end());
                    const size_t med = lens[lens.size()/2];
                    halfw_auto = (int)std::max<size_t>(31, (size_t)(0.30 * (double)med));
                }
            }
            // The derived value is NOT the default: it raised F1 on the tuning window but not
            // on the four held-out chr20 windows (mean change -0.002), so HALFW stays 31
            // (63 bp), the value of every reported measurement. CAPS_DBG_HALFW opts in.
            // Measurements: docs/extras/ENGINEERING_NOTES_CALLER.md.
            (void)halfw_auto;
            const int HALFW = std::getenv("CAPS_DBG_HALFW")
                            ? atoi(std::getenv("CAPS_DBG_HALFW")) : 31;
            fprintf(stderr, "[HALFW] derived=%d (fragment %d bp) from measured read length\n",
                    HALFW, 2*HALFW+1);
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
                // A probe that starts at the path never includes the flank. For an SNV that is
                // harmless (the variant base is position 0 of the k-mer), but for an INDEL the
                // discriminating signal is whether reads cross the FLANK->PATH junction. As in
                // Lancet and Scalpel, support is evaluated against the raw reads at the place that
                // distinguishes the event.
                // TWO PROBES: one at the path start and one centred on the junction. A read
                // matches if EITHER is present. A single centred probe alone required every
                // supporting read to span 15 flank bases as well as the event and cost more
                // recall than the precision it bought (measured F1 0.611 -> 0.535).
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
            // ── PLACEMENT DISPERSION (instrumentation) ──
            // For every bubble, record which CONTIGS its supporting reads were
            // placed on by the compressor. A genuine heterozygous site draws
            // both alleles from reads at ONE locus. A repeat/paralog collapse
            // -- the dominant false-positive mechanism here, measured at 12x
            // the allele depth of true sites -- draws support from reads the
            // compressor placed at DIFFERENT loci, because they came from
            // different genomic copies.
            //
            // k-mer counts cannot see this distinction by construction: the
            // collapsed copies are identical in k-mer space, which is why they
            // collapsed. Placement is the only signal that separates them, and
            // it is exactly what a linked/read-coloured graph is built to
            // supply.
            //
            // This pass only COUNTS; it rejects nothing unless PLACEMODE (a diagnostic) is
            // set.
            const bool PLACE_ON = !cd.read_ppos.empty();
            // PER-PATH locus sets. A distinct-contig COUNT would measure
            // fragmentation, not collapse: the encoder's contigs are fragmented
            // (451,578 at full chr20), so one locus legitimately spans many of them.
            // What separates a real
            // het site from a collapse is whether the TWO ALLELES come from the
            // SAME contigs: one locus with two haplotypes shares its placements,
            // two collapsed copies do not. That comparison is immune to
            // fragmentation because both paths fragment identically.
            // Global pseudogenome span of the supporting reads. A real
            // heterozygous site is ONE locus, so every supporting read -- both
            // alleles -- was placed within about a read length of the same
            // pseudogenome offset. A repeat/paralog collapse pools reads the
            // compressor placed at genuinely different offsets, because they
            // came from different copies. That is a distance, and distance is
            // only meaningful in a single coordinate system.
            std::vector<uint64_t> pmin(PLACE_ON ? dbg_bubbles.size() : 0, UINT64_MAX);
            std::vector<uint64_t> pmax(PLACE_ON ? dbg_bubbles.size() : 0, 0);
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
            // satisfy both `find`s, and a site with no allelic evidence then
            // clears COHC on both paths. kissreads2 assigns each read to its BEST
            // path instead of every path it contains. This is a candidate cause of the
            // lower indel than SNV precision.
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
                                // coverage AND average phred. A sequencing
                                // error is a low-quality base by definition,
                                // so a read whose variant base is low quality
                                // is not evidence for that allele. This
                                // separates a real low-coverage het allele
                                // (few reads, high quality) from an error (few
                                // reads, low quality), which a coverage floor
                                // cannot do.
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
                                if (PLACE_ON && ri < (long long)cd.read_ppos.size()) {
                                    const uint64_t gp = cd.read_ppos[(size_t)ri];
                                    if (gp != UINT64_MAX) {
                                        #pragma omp critical(bloci)
                                        { if (gp < pmin[pr.bi]) pmin[pr.bi] = gp;
                                          if (gp > pmax[pr.bi]) pmax[pr.bi] = gp; }
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
            // No per-position coverage test (kissreads2 fragment.cpp:86-91) is applied: the
            // counter here increments only after a read contains the WHOLE fragment, so every
            // position has the containment count and such a test cannot fire. kissreads2
            // also counts partially overlapping reads; reproducing that needs partial-overlap
            // alignment.
            // DIAGNOSTIC SPLIT, not a shipped filter: 1 = keep only bubbles
            // whose two alleles share a placement contig, 2 = keep only those
            // whose placements are disjoint. Running both and comparing TP/FP
            // measures whether placement discriminates at all, which is the
            // question that decides whether this signal is worth a filter.
            const int PLACEMODE = std::getenv("CAPS_DBG_PLACEMODE")
                                ? atoi(std::getenv("CAPS_DBG_PLACEMODE")) : 0;
            for (uint32_t bi = 0; bi < (uint32_t)dbg_bubbles.size(); ++bi) {
                if (PLACEMODE && PLACE_ON) {
                    // PLACEMODE = maximum allowed pseudogenome span, in bp.
                    // Split the bubbles by locus tightness and score each half.
                    const bool has = (pmin[bi] != UINT64_MAX);
                    const uint64_t sp = has ? (pmax[bi] - pmin[bi]) : UINT64_MAX;
                    if (PLACEMODE > 0  && !(has && sp <= (uint64_t)PLACEMODE))  { ++killed; continue; }
                    if (PLACEMODE < 0 && !(has && sp >  (uint64_t)(-PLACEMODE))) { ++killed; continue; }
                }
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
            if (PLACE_ON) {
                size_t hb[7] = {0,0,0,0,0,0,0};
                for (uint32_t bi = 0; bi < (uint32_t)dbg_bubbles.size(); ++bi) {
                    if (pmin[bi] == UINT64_MAX) { ++hb[6]; continue; }
                    const uint64_t sp = pmax[bi] - pmin[bi];
                    ++hb[sp <= 200 ? 0 : sp <= 1000 ? 1 : sp <= 10000 ? 2
                        : sp <= 100000 ? 3 : sp <= 1000000 ? 4 : 5];
                }
                fprintf(stderr, "[PLACE-SPAN] pseudogenome span of supporting reads over %zu bubbles: "
                                "<=200:%zu <=1k:%zu <=10k:%zu <=100k:%zu <=1M:%zu >1M:%zu none:%zu\n",
                        dbg_bubbles.size(), hb[0], hb[1], hb[2], hb[3], hb[4], hb[5], hb[6]);
            }
            fprintf(stderr, "[DBG-RSS] after coherence: %ldMB\n", rss_kb()/1024);
            fprintf(stderr, "[DBG-COH] C=%d kept=%zu killed=%zu  %.2fs\n",
                    COHC, dbg_bubbles.size(), killed, elapsed_s(t_coh, clk::now()));
        }

        // ── BUBBLE CHAINING ─────────────────────────────────────────────────
        // Unitigs "break at heterozygotes": in a heterozygous genome the junction that
        // stops a unitig extension is the next variant, so extending each bubble into its
        // unitig barely lengthens it (mean 63 -> 238 bp on the tuning window, F1
        // unchanged). A BUBBLE CHAIN (BubbleGun: a linear stretch of bubbles where the
        // sink of one is the source of the next) extends THROUGH consecutive bubbles and
        // gives a long, uniquely mappable sequence carrying several variants. It is also
        // the basis of phasing (PHASM), since consecutive bubbles spanned by the same
        // reads are linked.
        //
        // Implementation: index sites by entry node; a site whose right unitig extension
        // lands on another site's entry node is its successor. Walk each maximal chain
        // once, marking members consumed.
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
                    // Reimplemented from the published DiscoSNP++ algorithm (kissnp2
                    // Bubble.cpp:952, `checkRepeatSize`, AGPL); no code is copied.
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
                    // This is a question about the GRAPH ("is this indel's position determined"),
                    // not about the sequence ("is the reference context repetitive"). The latter does
                    // not separate: over all 5,001 scored indels, true positives are MORE
                    // homopolymeric than false ones (fraction with a run >= 4: 0.531 against 0.427).
                    // A real indel in a homopolymer still has a determined position when the flanks
                    // pin it down; a spurious one does not.
                    //
                    // Threshold is DiscoSNP++'s own published default (20 at
                    // k=31, i.e. the shorter path must stay distinct for >= 9
                    // bases), so this is their gate at their setting, not a
                    // value fitted here.
                    //
                    // SNVs are untouched by construction: their code returns
                    // true immediately for equal-length paths, and so does this.
                    // ── CLEAN-INDEL SHAPE (DiscoSNP++ expand() geometry) ─────
                    // DiscoSNP++ extends ONE path until it can close and then walks BOTH paths in
                    // LOCKSTEP, closing only when the next nodes are equal, so its two alleles differ
                    // in exactly one contiguous block by construction. Two unitig traces that reach a
                    // common node via find_sb carry no such guarantee: paths differing at SEVERAL
                    // positions would be emitted as one multi-base REF/ALT, which vcfeval scores as an
                    // INDEL. The shape test is exact and parameter-free: s_long is s_short with ONE
                    // contiguous block inserted iff the common prefix and common suffix together
                    // cover the whole of s_short:
                    //
                    //     LCP + LCS >= |s_short|
                    //
                    // Anything failing it differs by more than a single indel and is a traversal
                    // artefact, not a variant.
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
                    // DiscoSNP++'s indel search is a BFS that STOPS at the smallest closing indel
                    // (Bubble.cpp:238-243):
                    //     if (insert_size == found_del_size-1) { clear(); break; }
                    //     if (insert_size >  found_del_size)   continue;
                    // so long indels are never generated. A superbubble exit reached through a repeat
                    // is often a distant one, which produces the 21-32 bp false-positive mode. A
                    // blanket length rejection after the fact also removes real events (measured on
                    // the other channel: F1 0.637 -> 0.542), whereas preferring the shortest closure
                    // does not propose the long ones. The closest analogue here: when a site yields
                    // several indel-shaped allele pairs, keep the one with the smallest length
                    // difference.
                    // ── HOMOPOLYMER / STR ARTEFACT FILTER: OPT-IN (CAPS_DBG_STR=1) ──
                    // Rejects an indel whose inserted or deleted sequence extends a run or tandem
                    // unit already present in the flank (is_str_event). Off by default: it removed
                    // more true indels than false ones (F1 0.601 -> 0.535), consistent with the
                    // earlier finding that homopolymer context does not separate true from false
                    // indels here.
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
                    // At window scale a length cap only removed true positives (the window's few
                    // indel false positives were short); the long false-positive mode appears only at
                    // full chr20, so the cap is evaluated there.
                    // GRAPH-CHANNEL INDEL EMISSION IS OFF BY DEFAULT (CAPS_DBG_INDEL=1 enables it).
                    // The superbubble path finds real indels (held-out indel TP rose on every window)
                    // but adds more false positives than it recovers: held-out indel F1 fell on 3 of
                    // 4 windows (mean -0.051). The machinery stays because it raises SNV precision
                    // (0.925 -> 0.957); only its indel OUTPUT is gated. The reported indel results
                    // come from the full caller (indel_pass), selected in `capsule_decode call` by
                    // CAPS_CALL_INDELS=1.
                    // The gate tests the shape of the EMITTED RECORD, not the path lengths: two
                    // same-length paths can differ at several positions (e.g. ACGT vs TGCA) and would
                    // be emitted as a multi-base REF/ALT, which vcfeval scores as an INDEL.
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
        fprintf(stderr, "[DBG-IBFS] independent indel discovery: %zu records (enabled=%d, maxindel=%d)\n",
                n_ibfs_found, (int)WANT_IBFS, IBFS_MAXINDEL);
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
                    // in contig space; such bubbles are not placed (counted, not
                    // guessed).
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
            // This filters on EVIDENCE instead of topology, using per-path read
            // support (min1/min2) that the walk already computed:
            //
            //   MINSUP  both paths must have real read support. A repeat- or
            //           error-induced path typically has a weak link somewhere.
            //   BAL     a true heterozygous site is ~50/50, so the two paths'
            //           support should be comparable. Repeats skew hard: the
            //           repeated copy accumulates counts from every genomic
            //           instance while the true alternate does not.
            //
            // Both are swept on the tuning window and then verified on held-out
            // windows.
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
                // DiscoSNP++ writes one VCF record per SNP inside a bubble, and so does this;
                // emitting one multi-base REF/ALT would be discarded by the clean-SNV gate.
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
                // H. Low-coverage thresholds therefore filter the wrong tail.
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
                    // The bound during the superbubble walk measures the PATH length difference and
                    // rarely fires; the false positives become 20-37 bp only after the lift's
                    // left-alignment, so the bound is applied to the emitted record. On the
                    // tetraploid window a 15 bp cap removed 12 of 16 indel false positives and no true
                    // positives (F1 0.537 -> 0.577); at full chr20 true indels are 1-8 bp and the false
                    // ones form a second mode at 21-32 bp.
                    // Event size is the length DIFFERENCE between the two paths. Bounding the path
                    // lengths themselves rejects almost every real event (F1 0.537 -> 0.035), because
                    // the paths carry the whole bubble, flank and closure included.
                    const size_t esz = (b.path1.size() > b.path2.size())
                                     ? b.path1.size() - b.path2.size()
                                     : b.path2.size() - b.path1.size();
                    const int EMAXI = std::getenv("CAPS_DBG_EMITMAXINDEL")
                                    ? atoi(std::getenv("CAPS_DBG_EMITMAXINDEL")) : 15;
                    if ((int)esz > EMAXI) { ++n_ind_drop_emit; continue; }
                    // Allele-string length is not used as a filter: real indels also have long
                    // alleles (path1/path2 carry flank, event and closure), and bounding the shorter
                    // allele gave the same collapse in F1 as bounding the path lengths.
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
    // DiscoSNP++ runs, and nothing else, so it is the like-for-like comparison:
    // same method, this implementation, this implementation's cost.
    // ── PG ANCHORING: emit against the ASSEMBLED CONTIGS, not the bubble ────
    // MEASURED PROBLEM: at full chr20, 29% of emitted bubble sequences are
    // under 100 bp and 10.3% fail to align at all (71% of those are the short
    // ones). Those calls are never scored -- pure lost recall, and it is an
    // artefact of using a short unitig walk as the mappable context. Unitig
    // extension stops at the first junction, and at full scale junctions are
    // everywhere.
    //
    // The alternative is to locate the bubble's flank k-mer in one of the encoder's
    // long assembled contigs (cd_in.contigs) and emit the call in that contig's
    // coordinates, which DiscoSNP++ cannot do without an assembly.
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
                        // REF must be the ASSEMBLED base; ALT is whichever path disagrees with it.
                        // Taking path1 as REF regardless is wrong whenever the assembly carries path2.
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
    // DEPTH GUARD: reject columns deeper than 2.5x the MEDIAN CANDIDATE DEPTH
    // (collapsed repeats), as in the ARCS caller this pileup follows. It uses the
    // median of the candidates, not a ceiling derived from H, because after the
    // collapse pass a legitimate pileup is deeper than raw coverage (reads from
    // several duplicate contigs stack on one); an H-derived ceiling rejects real
    // sites. The median self-calibrates to any coverage.
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
        // K-MER SANITY. The gate rejects REPEATS, whose flanks occur at multiples of the
        // homozygous depth (2H, 3H...). A legitimate flank occurs at about 1.0H, so a
        // ceiling of 1.1H would reject real sites on noise. The repeat boundary is 2H, so
        // the decision point is the midpoint 1.5H ("closer to two copies than to one"),
        // derived from the gate's purpose. KHI is kept for the uncollapsed path.
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
    //
    // TWO SUBSTRATES FROM ONE ASSEMBLY. SNVs and indels need different amounts of
    // contig collapse. The SNV pileup wants both haplotypes' reads in ONE frame
    // (aggressive collapse, dup_frac 0.45). A bubble needs the two haplotypes to
    // remain TWO comparable contigs, with only same-haplotype duplicates removed (mild
    // collapse, bub_dup 0.92). On the tuning window (HG002 r2):
    //     dup=0.45    -> SNV F1 0.882, INDEL F1 0.217
    //     dup=0.80    -> SNV F1 0.821, INDEL F1 0.422
    //     no collapse -> SNV F1 0.856, INDEL F1 0.349
    // so the SNV pileup and the bubble passes each run on their own substrate.
    // CAPS_BUBBLE_UNCOLLAPSED=1 runs the bubble passes on the uncollapsed contigs
    // instead (re-measurement only).
    // indel_pass timing starts here, so the second build_substrate is included.
    g_ip_start = clk::now();
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
    g_ipsum += elapsed_s(g_ip_start, clk::now());
    fprintf(stderr, "[INDEL-PROF] %-24s %8.2fs\n", "2nd build_substrate",
            elapsed_s(g_ip_start, clk::now()));
    g_ip_last = clk::now();
    const CallData& cdb = BUB_UNCOL ? cd_in : cd_bub;
    const int BUB_SRC = 1;   // bubble records always live in cdb's own contig space
    // No per-candidate local read realignment is applied: scoped per contig it
    // reduces to DROP-COV (per-contig read coverage, implemented), and scoped to
    // cross-locus verification it needs 60+ bp of matching context that these
    // contigs are too short to supply.
    if (!std::getenv("CAPS_NO_INDELS") && cdb.contigs.size() >= 2) {
        constexpr int BK = 25, FLANK = 15;
        // PARTITIONED BY CONTIG, NOT BY READ.
        //
        // This is 4M reads x ~148 bases = ~590M increments, serial. Reads
        // cannot simply be split across threads because two reads on the same
        // contig would race on the same counters. Splitting by `cid % T`
        // instead gives each thread a disjoint set of CONTIGS, so no two
        // threads ever touch the same counter, and each thread still visits
        // reads in increasing oi order -- so the increment sequence at every
        // position is exactly the serial one, including the 60000 clamp's
        // order dependence. The cost is that every thread scans the 4M read
        // headers (cheap) to do 1/T of the ~590M increments.
        std::vector<std::vector<uint16_t>> cov(cdb.contigs.size());
        #pragma omp parallel for schedule(dynamic, 256)
        for (long long ci = 0; ci < (long long)cdb.contigs.size(); ++ci)
            cov[(size_t)ci].assign(cdb.contigs[(size_t)ci].size(), 0);
        {
            int cov_nt = 1;
            #ifdef _OPENMP
            cov_nt = omp_get_max_threads();
            #endif
            #pragma omp parallel for schedule(static, 1)
            for (int t = 0; t < cov_nt; ++t) {
                for (size_t oi = 0; oi < n; ++oi) {
                    uint32_t cid = cdb.read_cid[oi], pos = cdb.read_pos[oi];
                    if (cid >= cdb.contigs.size()) continue;
                    if ((int)(cid % (uint32_t)cov_nt) != t) continue;
                    uint16_t clipb = (oi < cdb.read_clip.size()) ? cdb.read_clip[oi] : 0;
                    int rl = (int)seqs[oi].size() - (int)clipb;
                    for (int j = 0; j < rl; ++j) {
                        uint32_t p = pos + (uint32_t)j;
                        if (p < cov[cid].size() && cov[cid][p] < 60000) ++cov[cid][p];
                    }
                }
            }
        }
        // KIDX AS A SORTED FLAT ARRAY (one record per contig k-mer position), not a
        // hash map of vectors; it is the indel pass's largest structure. No stride
        // subsampling or per-key cap: each key's true occurrence count is what the
        // MAXOCC filter below needs, so capping could change which keys pass.
        // 16-BYTE RECORDS: {uint64 kmer; uint32 ci, pos; uint8 orient} would pad to
        // 24 B. A canonical 25-mer is 50 bits, so the orientation bit rides in bit 63.
        struct KIdxEntry {
            uint64_t kf;                       // bit 63 = orient, bits 49..0 = k-mer
            uint32_t ci, pos;
            enum : uint64_t { ORI = 1ULL << 63 };
            inline uint64_t kmer()   const { return kf & ~(uint64_t)ORI; }
            inline uint8_t  orient() const { return (uint8_t)(kf >> 63); }
        };
        // BUILT IN PLACE, like the seed index, so no sort temporary doubles it. Entries
        // are bucketed by the top bits of the RAW key (not of mix(key): this index is
        // binary-searched by key and must stay in key order), and each bucket is sorted
        // in place with an explicit (key, ci, pos) comparator. The fill is contig-major,
        // position-minor, so this is the order a stable sort would give. With no per-key
        // cap, buckets stay contiguous and no compaction is needed.
        std::vector<KIdxEntry> kidx;
        {
            const size_t NC = cdb.contigs.size();
            static constexpr int KBB = 12;
            const size_t NKB = (size_t)1 << KBB;
            int T = 1;
            #ifdef _OPENMP
            T = omp_get_max_threads();
            #endif
            std::vector<size_t> ccnt(NC + 1, 0);
            #pragma omp parallel for schedule(dynamic, 256)
            for (long long ci = 0; ci < (long long)NC; ++ci) {
                const std::string& c = cdb.contigs[(size_t)ci];
                size_t cnt = 0; Roll25 rw; const size_t cn2 = c.size();
                for (size_t b2 = 0; b2 < cn2; ++b2) {
                    rw.push(c[b2]);
                    if (b2 + 1 >= (size_t)BK && rw.ok()) ++cnt;
                }
                ccnt[(size_t)ci + 1] = cnt;
            }
            for (size_t i = 0; i < NC; ++i) ccnt[i + 1] += ccnt[i];
            const size_t total = ccnt[NC];
            std::vector<size_t> tlo((size_t)T + 1, NC);
            tlo[0] = 0;
            for (int t = 1; t < T; ++t) {
                const size_t want = total * (size_t)t / (size_t)T;
                size_t at = (size_t)(std::lower_bound(ccnt.begin(), ccnt.end(), want) - ccnt.begin());
                if (at > NC) at = NC;
                if (at < tlo[(size_t)t - 1]) at = tlo[(size_t)t - 1];
                tlo[(size_t)t] = at;
            }
            tlo[(size_t)T] = NC;
            std::vector<std::vector<uint32_t>> hist((size_t)T, std::vector<uint32_t>(NKB, 0));
            #pragma omp parallel for schedule(static, 1)
            for (int t = 0; t < T; ++t) {
                std::vector<uint32_t>& h = hist[(size_t)t];
                for (size_t ci = tlo[(size_t)t]; ci < tlo[(size_t)t + 1]; ++ci) {
                    const std::string& c = cdb.contigs[ci];
                    Roll25 rw; const size_t cn2 = c.size();
                    for (size_t b2 = 0; b2 < cn2; ++b2) {
                        rw.push(c[b2]);
                        if (b2 + 1 < (size_t)BK || !rw.ok()) continue;
                        ++h[(size_t)(rw.canon() >> (50 - KBB))];
                    }
                }
            }
            std::vector<size_t> boff(NKB + 1, 0);
            for (size_t b = 0; b < NKB; ++b) {
                size_t sum = 0;
                for (int t = 0; t < T; ++t) sum += hist[(size_t)t][b];
                boff[b + 1] = sum;
            }
            for (size_t b = 0; b < NKB; ++b) boff[b + 1] += boff[b];
            std::vector<std::vector<size_t>> wpos((size_t)T, std::vector<size_t>(NKB, 0));
            for (size_t b = 0; b < NKB; ++b) {
                size_t cur = boff[b];
                for (int t = 0; t < T; ++t) { wpos[(size_t)t][b] = cur; cur += hist[(size_t)t][b]; }
            }
            kidx.resize(total);
            hp_hint((void*)kidx.data(), kidx.size() * sizeof(KIdxEntry), "kidx");
            #pragma omp parallel for schedule(static, 1)
            for (int t = 0; t < T; ++t) {
                std::vector<size_t>& w = wpos[(size_t)t];
                for (size_t ci = tlo[(size_t)t]; ci < tlo[(size_t)t + 1]; ++ci) {
                    const std::string& c = cdb.contigs[ci];
                    Roll25 rw; const size_t cn2 = c.size();
                    for (size_t b2 = 0; b2 < cn2; ++b2) {
                        rw.push(c[b2]);
                        if (b2 + 1 < (size_t)BK || !rw.ok()) continue;
                        const uint64_t kk = rw.canon();
                        kidx[w[(size_t)(kk >> (50 - KBB))]++] =
                            KIdxEntry{ kk | (rw.f <= rw.r ? 0ULL : (uint64_t)KIdxEntry::ORI),
                                       (uint32_t)ci, (uint32_t)(b2 + 1 - (size_t)BK) };
                    }
                }
            }
            #pragma omp parallel for schedule(dynamic, 8)
            for (long long b = 0; b < (long long)NKB; ++b) {
                std::sort(kidx.begin() + (long)boff[(size_t)b],
                          kidx.begin() + (long)boff[(size_t)b + 1],
                          [](const KIdxEntry& x, const KIdxEntry& y) {
                              const uint64_t kx = x.kmer(), ky = y.kmer();
                              if (kx != ky)       return kx < ky;
                              if (x.ci != y.ci)   return x.ci < y.ci;
                              return x.pos < y.pos;
                          });
            }
        }
        auto kidx_run_len = [&](uint64_t key) -> int {
            auto lo = std::lower_bound(kidx.begin(), kidx.end(), key,
                [](const KIdxEntry& e, uint64_t k){ return e.kmer() < k; });
            auto hi = std::upper_bound(kidx.begin(), kidx.end(), key,
                [](uint64_t k, const KIdxEntry& e){ return k < e.kmer(); });
            return (lo != hi) ? (int)(hi - lo) : -1; // -1 == not found (mirrors kidx.end())
        };
        struct Agg { int anchors = 0; uint32_t altcid = 0, altpos = 0; };
        // pair -> (shared anchors, representative anchor offsets) for the
        // whole-pair alignment scan below
        struct PairRep { int n = 0; uint32_t rp = 0, qB = 0, ap = 0; bool opp = false; };
        // ANCHOR MULTIPLICITY. Every cross-contig PAIR among a small number of
        // occurrences is tried, not only anchors that occur exactly twice: a real
        // hap1/hap2 pair whose k-mer also occurs in a third fragment must not be lost. The
        // bubble extractor rejects pairs that do not diverge and reconverge, so extra
        // candidates cost specificity only where the geometry looks like a bubble.
        // MAXOCC is bounded to keep repeats from exploding the pair count --
        // that bound only holds if CAPS_MAXOCC itself is range-checked: the
        // nested pair scan below is O(occ^2), so an unbounded override
        // defeats the exact guard this comment describes. Same convention as
        // CAPS_PLOIDY: an out-of-range override is ignored, default kept.
        int MAXOCC = 4;
        if (const char* e = std::getenv("CAPS_MAXOCC")) {
            int v = atoi(e); if (v >= 1 && v <= 64) MAXOCC = v;
        }
        // ── PARALLELISE THE ANCHOR SCAN ─────────────────────────────────────
        // This loop and its nested O(occ^2) pair scan are the bulk of
        // indel_pass, which is 65% of the full CAPS_CALL run (1250 s of ~1900)
        // and had NOT ONE `#pragma omp` in ~1100 lines -- single-threaded on a
        // 12-core box, immediately after `parallel_loop`, which parallelises
        // the same kind of per-contig work.
        //
        // Each kidx RUN (one canonical 25-mer's occurrences) is independent:
        // the run body reads cdb.contigs/cov and writes only two shared maps.
        // But both are ORDER-DEPENDENT, so a bare `#pragma omp for` would
        // silently change results -- exactly the failure the stable_sort note
        // above records (a plain sort moved anchor positions and DP/AF, caught
        // only by byte-comparison):
        //
        //   pairs_[{rc_,ac}]  FIRST write wins   (`if (pr_.n == 0) {...}`)
        //   im[key]           anchors++ accumulates, but altcid/altpos are
        //                     overwritten, so the LAST write wins
        //
        // So: find the run boundaries first, let each thread accumulate into
        // its OWN maps, then merge strictly in ascending run order. Replaying
        // the merge in that order reproduces first-wins and last-wins exactly,
        // for all THREE `im` write sites (the multi-start extraction at the
        // b2 site, the main aggregate, and the whole-pair scan below).
        //
        // Verification is by byte-comparing the VCF against the serial build,
        // never by inspection.
        // NOTE: g_ip_start is set at the TOP of the indel block so the laps
        // below cover the same span phase("indel_pass") measures. Without it
        // the second build_substrate, which runs before the first lap, is counted.
        auto _ip0 = g_ip_last;   // continue from the 2nd build_substrate lap
        auto _iplap = [&](const char* what){
            const double d = elapsed_s(_ip0, clk::now());
            g_ipsum += d;
            g_ip_last = clk::now();
            fprintf(stderr, "[INDEL-PROF] %-24s %8.2fs  rss=%ldMB peak=%ldMB\n", what, d,
                    caps_rss_mb(), caps_hwm_mb());
            _ip0 = clk::now();
        };
        // KEEP ONLY THE RUNS THAT DO WORK.
        //
        // Most runs are no-ops (the body returns at once when occ.size() < 2 or
        // > MAXOCC), so threading over ALL runs spends its time on scheduling.
        //
        // The pair scan is O(occ^2), so cost is concentrated in the few runs
        // with many occurrences. Filtering first turns tens of millions of
        // trivial tasks into a small set of substantial ones, and lets the
        // schedule be guided by actual cost.
        std::vector<std::pair<size_t,size_t>> runs;
        {
            size_t total = 0, kept = 0;
            for (size_t r0 = 0; r0 < kidx.size(); ) {
                size_t r1 = r0;
                while (r1 < kidx.size() && kidx[r1].kmer() == kidx[r0].kmer()) ++r1;
                ++total;
                const size_t occn = r1 - r0;
                if ((int)occn >= 2 && (int)occn <= MAXOCC) { runs.push_back({r0, r1}); ++kept; }
                r0 = r1;
            }
            fprintf(stderr, "[INDEL-PAR] %zu runs total, %zu carry work (%.2f%%)\n",
                    total, kept, total ? 100.0*kept/total : 0.0);
            // Largest first: with O(occ^2) bodies, a big run scheduled last
            // leaves 11 threads idle waiting for it.
            // Sorting largest-first changes the ORDER runs are visited, and
            // the two shared maps are order-sensitive (pairs_ first-wins, im
            // last-wins). So each run keeps its ORIGINAL index and the merge
            // replays those semantics explicitly; relying on thread-visit order would not
            // match the serial result, because schedule(dynamic) does not preserve order.
            std::sort(runs.begin(), runs.end(),
                      [](const std::pair<size_t,size_t>& a, const std::pair<size_t,size_t>& b){
                          return (a.second - a.first) > (b.second - b.first); });
        }
        const bool IND_PAR = std::getenv("CAPS_INDEL_SERIAL") == nullptr;
        int n_thr = 1;
        #ifdef _OPENMP
        if (IND_PAR) n_thr = omp_get_max_threads();
        #endif
        // (payload, original-run-index) so the merge can pick min/max
        struct PairRepR { PairRep v; size_t run = (size_t)-1; };
        struct AggR     { Agg v;     size_t run = 0; };
        std::vector<std::map<std::pair<uint32_t,uint32_t>, PairRepR>> tp((size_t)n_thr);
        std::vector<std::map<std::tuple<uint32_t,uint32_t,int,int,std::string>, AggR>> ti((size_t)n_thr);
        std::vector<size_t> run_order(runs.size());
        {   // map each (sorted) slot back to its position in kidx order
            std::vector<std::pair<size_t,size_t>> tmp(runs.size());
            for (size_t i = 0; i < runs.size(); ++i) tmp[i] = {runs[i].first, i};
            std::sort(tmp.begin(), tmp.end());
            for (size_t k = 0; k < tmp.size(); ++k) run_order[tmp[k].second] = k;
        }
        _iplap("indel setup (cov+kidx+runs)");
        fprintf(stderr, "[INDEL-PAR] %zu working runs over %d threads\n", runs.size(), n_thr);
        #pragma omp parallel for schedule(dynamic, 1) if(IND_PAR)
        for (long long ri = 0; ri < (long long)runs.size(); ++ri) {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            auto& tpairs = tp[(size_t)tid];
            auto& tim    = ti[(size_t)tid];
            const size_t myrun = run_order[(size_t)ri];   // position in ORIGINAL kidx order
            size_t run_i = runs[(size_t)ri].first, run_j = runs[(size_t)ri].second;
            {
            // Materialize this run into the exact type/shape the untouched
            // body below already expects -- so nothing past this point in
            // the loop changes at all, only how `occ` gets populated.
            std::vector<std::tuple<uint32_t,uint32_t,uint8_t>> occ;
            occ.reserve(run_j - run_i);
            for (size_t t = run_i; t < run_j; ++t)
                occ.push_back(std::make_tuple(kidx[t].ci, kidx[t].pos, kidx[t].orient()));
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
            // nearest difference. For the missed indels on the tuning window, usable shared
            // anchors were abundant but far from the event (213-135,453 bp), so the walk ended
            // on another difference first. After taking a bubble the walk continues past it
            // and extracts again, so one anchor pair can yield several events along its
            // length.
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
                    auto& a2r = tim[std::make_tuple(rc_, b2.apos, b2.type, b2.len, b2.ins)];
                    a2r.v.anchors++;
                    if (a2r.v.anchors == 1 || myrun >= a2r.run) { a2r.v.altcid = ac; a2r.v.altpos = ap; a2r.run = myrun; }
                }
            }
            if (!bub.ok || bub.apos == 0) continue;
            // ── CLOSING ANCHOR: make the bubble anchored at BOTH ends ───────
            // DiscoSNP++'s bubble leaves a shared graph node and must RE-CONVERGE on another
            // shared node, so both ends of both paths are k-mers present in the reads. An
            // opening 25-mer anchor plus 15 bases of contig identity at the close can be
            // satisfied by chance by two paralogs, which lowers indel precision. So the
            // re-convergence point must be a real closing anchor: a full BK-mer identical in
            // both contigs AND present in the reads. This is the both-ends-anchored property,
            // obtained on the contig set without a graph.
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
                // exactly that amount (all of the 1 bp indels missed on the tuning window are
                // such events). The alt side is slid by up to the event length and re-tested; a
                // real run-length change closes exactly there.
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
                // On the tetraploid benchmark 88% of candidate drops were at the closing-anchor
                // test: a byte-identical BK=25-mer is required on BOTH contigs, so one
                // heterozygous SNV in those 25 bases (common at about 1 het per 1.3 kb in human
                // data) destroys an otherwise valid bubble. Tolerance is justified because the next
                // test requires the closing anchor to be observed in the READS (kcount >= MC):
                // propose permissively, then validate against reads, as DiscoSNP++ does with
                // kissreads2. Tolerance reuses the frozen HDMAX (=2) instead of a new constant.
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
                // OPT-IN: requiring a unique closing anchor lowered indel F1 (five-window mean
                // 0.5976 -> 0.5932); a repeated closing k-mer still closes a real bubble.
                if (std::getenv("CAPS_CLOSE_UNIQ")) {
                    uint64_t cv;
                    if (pack25(R.data() + ra, cv)) {
                        uint64_t crv = rc25(cv), ccan = cv < crv ? cv : crv;
                        int rl = kidx_run_len(ccan);
                        if (rl < 0 || rl > MAXOCC) continue;
                    }
                }
            }
            // EXTENDED AGREEMENT. Read-support tests cannot separate false bubbles here,
            // because every candidate is built from contigs built from reads. What separates
            // a real haplotype pair from a paralog pair sharing one 25-mer is how far the two
            // contigs agree on either side of the event: haplotypes agree over hundreds of
            // bases, paralogs diverge quickly. This verifies a long span rather than the
            // 15 bp flanks the bubble extractor checks.
            {
                // OFF by default (EXT=0, CAPS_BUB_EXT to enable): the contigs are too short to
                // supply 60+ bp either side, so it rejected true bubbles (five-window mean indel
                // F1 0.5930 -> 0.5864 at 60 bp).
                // TARGETED PARALOG FILTER: extended agreement applied only to long events. The
                // false-positive burden is concentrated at len >= 3 (pooled over 8 evaluations:
                // precision 0.55 against 0.83 at len 1), and a long indel implies longer contigs,
                // so the context is available where it is needed.
                int EXT = 0, EXT_TOL = 3;
                // OFF by default (CAPS_LONG_EXT): neutral at 40 bp, because paralogs agree over
                // 40 bp easily and longer spans exceed the contig lengths.
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
                auto& prr = tpairs[{rc_, ac}];
                if (prr.v.n == 0 || myrun < prr.run) {
                    if (prr.v.n == 0) { prr.v.rp = rp; prr.v.qB = qB; prr.v.ap = ap; prr.v.opp = opp; prr.run = myrun; }
                    else { prr.v.rp = rp; prr.v.qB = qB; prr.v.ap = ap; prr.v.opp = opp; prr.run = myrun; }
                }
                ++prr.v.n;
            }
            if (std::getenv("CAPS_TRACE"))
                fprintf(stderr, "[trace] AGG cid=%u apos=%u type=%d len=%d\n", rc_, bub.apos, bub.type, bub.len);
            auto& ar = tim[std::make_tuple(rc_, bub.apos, bub.type, bub.len, bub.ins)];
            ar.v.anchors++;
            if (ar.v.anchors == 1 || myrun >= ar.run) { ar.v.altcid = ac; ar.v.altpos = ap; ar.run = myrun; }
          }
          }
        }
        // ── DETERMINISTIC MERGE ─────────────────────────────────────────────
        // Threads accumulated into their own maps in whatever order the
        // dynamic schedule handed them runs, so merging must RESTORE the
        // serial order's outcome, not just combine counts.
        //
        // LIMIT OF THIS MERGE. Visiting threads in ascending id order
        // does NOT reproduce serial run order: with `schedule(dynamic)` a
        // higher-numbered thread may have processed an EARLIER run. So:
        //
        //   anchors / n   sums are order-free and are therefore EXACT.
        //   pairs_ rep    (rp/qB/ap/opp) and im's altcid/altpos are picked by
        //                 thread-visit order, which may differ from serial.
        //
        // Those fields feed the whole-pair alignment scan and the alt-locus
        // annotation, so a difference here CAN move a record. This is exactly
        // the class of change the stable_sort note above says must be settled
        // by byte-comparison, not by reasoning -- so the VCF is compared
        // against a CAPS_INDEL_SERIAL=1 build on real data, and if it differs
        // the merge must be upgraded to carry a per-key run index and pick
        // min (first-wins) / max (last-wins) explicitly. Do not assume.
        std::map<std::pair<uint32_t,uint32_t>, PairRep> pairs_;
        std::map<std::tuple<uint32_t,uint32_t,int,int,std::string>, Agg> im;
        std::map<std::pair<uint32_t,uint32_t>, size_t> pr_run;
        std::map<std::tuple<uint32_t,uint32_t,int,int,std::string>, size_t> im_run;
        for (int t = 0; t < n_thr; ++t) {
            for (auto& kv : tp[(size_t)t]) {
                auto it = pairs_.find(kv.first);
                if (it == pairs_.end()) { pairs_.emplace(kv.first, kv.second.v); pr_run[kv.first] = kv.second.run; }
                else {
                    it->second.n += kv.second.v.n;               // counts are order-free
                    if (kv.second.run < pr_run[kv.first]) {      // FIRST run wins the rep
                        const int keep = it->second.n;
                        it->second = kv.second.v; it->second.n = keep;
                        pr_run[kv.first] = kv.second.run;
                    }
                }
            }
            for (auto& kv : ti[(size_t)t]) {
                auto it = im.find(kv.first);
                if (it == im.end()) { im.emplace(kv.first, kv.second.v); im_run[kv.first] = kv.second.run; }
                else {
                    it->second.anchors += kv.second.v.anchors;   // order-free
                    if (kv.second.run >= im_run[kv.first]) {     // LAST run wins altcid/altpos
                        it->second.altcid = kv.second.v.altcid;
                        it->second.altpos = kv.second.v.altpos;
                        im_run[kv.first] = kv.second.run;
                    }
                }
            }
        }
        { std::vector<std::map<std::pair<uint32_t,uint32_t>, PairRepR>>().swap(tp);
          std::vector<std::map<std::tuple<uint32_t,uint32_t,int,int,std::string>, AggR>>().swap(ti); }
        _iplap("anchor scan + merge");
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
        _iplap("medcov");
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
        // A true heterozygous site has exactly ONE alternate haplotype, but the anchor
        // search lets one locus pair with several alt contigs (and one alt contig serve
        // many loci), which is how paralogs and repeats can enter as false bubbles.
        // Matching is GLOBAL and needs no flanking context, which these short contigs
        // cannot supply.
        // Greedy by weight is the standard 1/2-approximation and is what the
        // evidence ordering here justifies (anchors = independent 25-mers
        // agreeing on the same pairing).
        std::vector<std::pair<int, const std::tuple<uint32_t,uint32_t,int,int,std::string>*>> order_;
        for (auto& kv : im) order_.push_back({kv.second.anchors, &kv.first});
        std::sort(order_.begin(), order_.end(),
                  [](auto& x, auto& y){ return x.first > y.first; });
        std::unordered_set<uint64_t> used_ref, used_alt;
        _iplap("scan_pair over pairs_");
        std::unordered_set<const void*> accepted;
        // MEASURED NEUTRAL (five-window mean 0.5976 -> 0.5962; r5 0.623->0.605,
        // HG004 0.515->0.521, HG003 0.672->0.677) and therefore opt-in.
        // The matching constraint assumes competing pairings are the problem, but false
        // bubbles here are mostly UNCONTESTED (a locus with exactly one wrong partner),
        // which a matching cannot reject.
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
        _iplap("matching");
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
            // Contig-coverage proxies (covwin / medcov) do not show that any READ carries the
            // alt allele. DiscoSNP++ gives that job to kissreads2.
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
                // fragment rejects TRUE indels: the alt fragment splices the REF contig's flanks
                // onto the ALT allele, so any nearby heterozygous SNV makes it a sequence no read
                // contains. eBWT2SNP avoids this by reading both haplotype fragments out of the
                // eBWT cluster instead of constructing one.
                int FL = 15;
                if (const char* e = std::getenv("CAPS_JUNC_FL")) FL = atoi(e);
                // ADAPTIVE CONTEXT. Fixing the half-length and REJECTING any
                // candidate without that much contig either side is what made
                // longer fragments look harmful (five-window mean indel F1
                // 0.5976 at 31bp, 0.5782 at 61bp, 0.5540 at 91bp): the loss was
                // good candidates near contig ends being discarded unexamined. So the context is
                // as much as exists, down to a floor, and the candidate is judged on that.
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
                // COHERENCE mode. DiscoSNP++'s indel precision (0.91 against 0.73 here)
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
                    uint32_t c1 = kcount_at(alt_hap.data() + q, alt_hap.size() - q);
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
                // exactly the paralog class bubble false positives come from.
                // It reuses the frozen KHI against H.
                // OPT-IN only: across five windows at 1.0xH it lowered indel F1 on every window
                // (it removed true indels), because alt-junction k-mer counts already sit near the
                // heterozygous level when the candidate comes from a contig built from those
                // reads.
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
                            ref_sup = std::max(ref_sup, kcount_at(ref_hap.data() + q, ref_hap.size() - q));
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
        // After eBWT2SNP: reads covering one genome position share a right context, so
        // grouping reads by a right-context anchor puts the differing allele immediately
        // before it, and both haplotype fragments come from the reads rather than from a
        // constructed candidate. The eBWT2SNP paper leaves the indel case unimplemented
        // ("extract the left-context and perform a local alignment ... future work"); this
        // channel implements it. The evidence is a real read carrying a real gap against
        // the contig; the contig supplies only coordinates and the REF allele. Bubble
        // false positives come from splicing one contig's flank onto another's allele,
        // which read-support tests cannot detect when every piece is read-derived. Every
        // read covering an anchor is compared with the SAME contig (not placed where it
        // fits best), so alt-haplotype reads reveal their gap.
        _iplap("im emit loop");
        if (!std::getenv("CAPS_NO_PCLUSTER")) {
            const int AK = 25, LW = 40;      // anchor k-mer, left window
            // SUBSTRATE CHOICE: the mildly collapsed bubble substrate (dup=0.92). On the
            // aggressively collapsed pileup substrate (dup=0.45, CAPS_PCLUSTER_COLSUB) the two
            // haplotypes merge into one contig, a read carrying the indel often matches a
            // contig that already contains it, and no gap appears (indel F1 0.655 -> 0.526 on
            // the tuning window).
            const CallData& pc_cd = std::getenv("CAPS_PCLUSTER_COLSUB") ? cd : cdb;
            // `cov` above is indexed by the BUBBLE substrate's contigs, so it
            // cannot be reused when this channel runs on a different one --
            // the ids would silently address the wrong contig. Build coverage
            // for whichever substrate this channel actually uses.
            // kidx is built over the BUBBLE substrate's contigs, so its ids
            // cannot index a different substrate. This channel builds its OWN
            // anchor index over whichever substrate it uses.
            // FLAT SORTED MAP, NOT A HASH MAP. `pv` is sorted by key and the selection walk
            // emits entries in that order, so the result is inherently sorted; ordered append
            // plus membership is what the semantics need.
            struct FlatAnchorMap {
                std::vector<uint64_t> k;
                std::vector<std::pair<uint32_t,uint32_t>> v;
                inline void append(uint64_t key, uint32_t ci, uint32_t pos) {
                    k.push_back(key); v.emplace_back(ci, pos);
                }
                inline size_t size() const { return k.size(); }
                inline bool count(uint64_t key) const {
                    return std::binary_search(k.begin(), k.end(), key);
                }
                inline size_t index_of(uint64_t key) const {
                    auto it = std::lower_bound(k.begin(), k.end(), key);
                    if (it == k.end() || *it != key) return (size_t)-1;
                    return (size_t)(it - k.begin());
                }
                // ── ONE-SHOT HASH DIRECTORY ─────────────────────────────────
                // The indel-pass read scan probes this map about 85 times per read. This resolves
                // a query in TWO memory touches whether it hits or misses: one probe into a
                // 2^23-entry offset table, then one cache line of a bucket array of (key, index)
                // pairs (about 2.6 entries per bucket). Buckets come from mix(key), because raw
                // canonical k-mers are strongly non-uniform. `k` is NOT reordered: the anchor
                // loop's iteration order over pkidx is load-bearing.
                enum { PDB = 23 };   // local class: enum, not a static member
                std::vector<uint32_t> dir;
                std::vector<std::pair<uint64_t,uint32_t>> buck;
                void build_dir() {
                    const size_t NB = (size_t)1 << PDB;
                    dir.assign(NB + 1, 0);
                    for (uint64_t key : k)
                        ++dir[(size_t)(FlatKmerSet::mix(key) >> (64 - PDB)) + 1];
                    for (size_t i = 0; i < NB; ++i) dir[i + 1] += dir[i];
                    buck.resize(k.size());
                    std::vector<uint32_t> w(dir.begin(), dir.end() - 1);
                    for (size_t i = 0; i < k.size(); ++i)
                        buck[w[(size_t)(FlatKmerSet::mix(k[i]) >> (64 - PDB))]++] =
                            { k[i], (uint32_t)i };
                }
                inline size_t lookup(uint64_t key) const {
                    const size_t h = (size_t)(FlatKmerSet::mix(key) >> (64 - PDB));
                    for (uint32_t i = dir[h], e = dir[h + 1]; i < e; ++i)
                        if (buck[i].first == key) return (size_t)buck[i].second;
                    return (size_t)-1;
                }
                inline const std::pair<uint32_t,uint32_t>* find_ptr(uint64_t key) const {
                    auto it = std::lower_bound(k.begin(), k.end(), key);
                    if (it == k.end() || *it != key) return nullptr;
                    return &v[(size_t)(it - k.begin())];
                }
                void reserve(size_t n) { k.reserve(n); v.reserve(n); }
            } pkidx;
            {
                // ── FLAT SORTED ARRAY, NOT THREE HASH MAPS ──────────────────
                // Computed by sorting rather than with hash maps (this block dominated indel_pass
                // time). Semantics:
                //   * count occurrences of the CANONICAL 25-mer across BOTH orientations
                //     (counting forward-canonical only lowered indel precision);
                //   * the FIRST occurrence supplies (contig, pos, orient); contigs are visited
                //     in increasing ci, so "first in sort order" IS "first seen";
                //   * survivors are the keys with count 1 that are forward-oriented. Both erase
                //     conditions are applied, so PUNIQ > 1 keeps its meaning (at the default
                //     PUNIQ=1 the second subsumes the first).
                // The build is count -> prefix-sum -> parallel fill by index, as for kidx, so entry
                // order matches the serial push order.
                // 16 BYTES, NOT 24.
                //
                // {uint64 k; uint32 ci, pos; uint8 fwd} would pad to 24 B on an array with one
                // entry per contig k-mer position. A canonical 25-mer occupies 50 bits, so the
                // orientation flag rides in bit 63 and the record is two 8-byte words. The sort
                // masks the flag off, so the ordering and the stability the run grouping depends
                // on are unchanged.
                struct PE {
                    uint64_t kf;                       // bit 63 = fwd, bits 49..0 = key
                    uint32_t ci, pos;
                    enum : uint64_t { FWD = 1ULL << 63 };   // local class: enum, not a static member
                    inline uint64_t key() const { return kf & ~(uint64_t)FWD; }
                    inline uint8_t  fwd() const { return (uint8_t)(kf >> 63); }
                };
                // pv is derived from kidx instead of being rebuilt. Both walk the same contigs
                // (pc_cd IS cdb unless CAPS_PCLUSTER_COLSUB substitutes the collapsed SNV
                // substrate), take 25-mers with the same Roll25 and ok() predicate, pack
                // canon | (f<=r ? 0 : bit63), store the same (ci,pos) and are stable-sorted by the
                // masked key, so they are the same array. kidx is still live here, so a rebuild
                // would duplicate it at the run's memory peak. The derivation is exact: same fill
                // order, same stable sort, so each run's first-seen entry (the one the occ==1 test
                // keeps) is the same entry. CAPS_PV_REBUILD=1 rebuilds it for A/B checking.
                const bool PV_FROM_KIDX = (&pc_cd == &cdb) && !std::getenv("CAPS_PV_REBUILD");
                std::vector<PE> pv;
                const size_t PNC = pc_cd.contigs.size();
                if (!PV_FROM_KIDX) {
                    std::vector<size_t> off(PNC + 1, 0);
                    #pragma omp parallel for schedule(dynamic, 256)
                    for (long long ci = 0; ci < (long long)PNC; ++ci) {
                        const std::string& c = pc_cd.contigs[(size_t)ci];
                        size_t cnt = 0;
                        Roll25 rw; const size_t cn2 = c.size();
                        for (size_t b2 = 0; b2 < cn2; ++b2) {
                            rw.push(c[b2]);
                            if (b2 + 1 >= 25 && rw.ok()) ++cnt;
                        }
                        off[(size_t)ci + 1] = cnt;
                    }
                    for (size_t i2 = 0; i2 < PNC; ++i2) off[i2 + 1] += off[i2];
                    pv.resize(off[PNC]);
                    #pragma omp parallel for schedule(dynamic, 256)
                    for (long long ci = 0; ci < (long long)PNC; ++ci) {
                        const std::string& c = pc_cd.contigs[(size_t)ci];
                        size_t w = off[(size_t)ci];
                        Roll25 rw; const size_t cn2 = c.size();
                        for (size_t b2 = 0; b2 < cn2; ++b2) {
                            rw.push(c[b2]);
                            if (b2 + 1 < 25 || !rw.ok()) continue;
                            const size_t i2 = b2 + 1 - 25;
                            pv[w++] = { rw.canon() | (rw.f <= rw.r ? 0ULL : (uint64_t)PE::FWD),
                                        (uint32_t)ci, (uint32_t)i2 };
                        }
                    }
                }
                int PUNIQ = 1;
                if (const char* e = std::getenv("CAPS_PCLUSTER_UNIQ")) PUNIQ = atoi(e);
                _iplap("pcluster: build flat array");
                if (PV_FROM_KIDX) {
                    // identical grouping, reading kidx in place
                    pkidx.reserve(kidx.size() / 8 + 1);
                    for (size_t a0 = 0; a0 < kidx.size(); ) {
                        size_t b0 = a0;
                        while (b0 < kidx.size() && kidx[b0].kmer() == kidx[a0].kmer()) ++b0;
                        const size_t occ = b0 - a0;
                        if (kidx[a0].orient() == 0 && (int)occ <= PUNIQ && occ == 1)
                            pkidx.append(kidx[a0].kmer(), kidx[a0].ci, kidx[a0].pos);
                        a0 = b0;
                    }
                    fprintf(stderr, "[PCLUSTER] %zu 25-mers -> %zu unique forward anchors (from kidx)\n",
                            kidx.size(), pkidx.size());
                } else {
                // STABLE sort: ties must keep (ci, pos) ascending so the first
                // entry of a run is the first-seen occurrence.
                #if defined(_OPENMP) && defined(_GLIBCXX_PARALLEL_ALGO_H)
                __gnu_parallel::stable_sort(pv.begin(), pv.end(),
                    [](const PE& a, const PE& b){ return a.key() < b.key(); });
                #else
                std::stable_sort(pv.begin(), pv.end(),
                    [](const PE& a, const PE& b){ return a.key() < b.key(); });
                #endif
                pkidx.reserve(pv.size() / 8 + 1);
                for (size_t a0 = 0; a0 < pv.size(); ) {
                    size_t b0 = a0;
                    while (b0 < pv.size() && pv[b0].key() == pv[a0].key()) ++b0;
                    const size_t occ = b0 - a0;
                    // first-seen entry decides contig/pos/orientation
                    if (pv[a0].fwd() == 0 && (int)occ <= PUNIQ && occ == 1)
                        pkidx.append(pv[a0].key(), pv[a0].ci, pv[a0].pos);
                    a0 = b0;
                }
                fprintf(stderr, "[PCLUSTER] %zu 25-mers -> %zu unique forward anchors\n",
                        pv.size(), pkidx.size());
                }
                { std::vector<PE>().swap(pv); }   // free immediately either way
                // SELF-CHECK (CAPS_PCLUSTER_VERIFY=1): replay the ORIGINAL
                // hash-map selection on the same data and diff. The synthetic
                // test said the predicates agree; this says whether they agree
                // ON THIS INPUT, which is the only question that matters.
                if (std::getenv("CAPS_PCLUSTER_VERIFY")) {
                    std::unordered_map<uint64_t,uint32_t> ocount;
                    std::unordered_map<uint64_t,uint8_t> oorient;
                    std::unordered_map<uint64_t,std::pair<uint32_t,uint32_t>> ok2;
                    ocount.reserve(1u<<21); ok2.reserve(1u<<21);
                    for (uint32_t ci = 0; ci < (uint32_t)pc_cd.contigs.size(); ++ci) {
                        const std::string& c = pc_cd.contigs[ci];
                        for (size_t i2 = 0; i2 + 25 <= c.size(); ++i2) {
                            uint64_t v; if (!pack25(c.data() + i2, v)) continue;
                            uint64_t rv = rc25(v), cn = v < rv ? v : rv;
                            if (++ocount[cn] == 1) { ok2[cn] = {ci,(uint32_t)i2};
                                                     oorient[cn] = (uint8_t)(v <= rv ? 0 : 1); }
                        }
                    }
                    for (auto it=ok2.begin(); it!=ok2.end(); )
                        if (oorient[it->first]!=0 || (int)ocount[it->first]>PUNIQ) it=ok2.erase(it); else ++it;
                    for (auto it=ok2.begin(); it!=ok2.end(); )
                        if (ocount[it->first]!=1) it=ok2.erase(it); else ++it;
                    size_t onlyA=0, onlyB=0, valdiff=0;
                    for (auto& kv : ok2) { const auto* pv2 = pkidx.find_ptr(kv.first);
                        if (!pv2) ++onlyA; else if (*pv2 != kv.second) ++valdiff; }
                    for (size_t ki = 0; ki < pkidx.k.size(); ++ki)
                        if (!ok2.count(pkidx.k[ki])) ++onlyB;
                    fprintf(stderr, "[PCLUSTER-VERIFY] original=%zu flat=%zu | only_orig=%zu "
                                    "only_flat=%zu value_diff=%zu\n",
                            ok2.size(), pkidx.size(), onlyA, onlyB, valdiff);
                }
                _iplap("pcluster: sort + select");
            }
            // pcov aliases cov when pc_cd IS cdb (the default): the coverage build would be
            // identical (same contigs, read_cid/read_pos/read_clip and 60000 clamp). It is
            // rebuilt only for the opt-in CAPS_PCLUSTER_COLSUB substrate.
            std::vector<std::vector<uint16_t>> pcov_own;
            if (&pc_cd != &cdb) {
                pcov_own.resize(pc_cd.contigs.size());
                for (size_t ci = 0; ci < pc_cd.contigs.size(); ++ci)
                    pcov_own[ci].assign(pc_cd.contigs[ci].size(), 0);
                for (size_t oi = 0; oi < n; ++oi) {
                    uint32_t cid = pc_cd.read_cid[oi], pos = pc_cd.read_pos[oi];
                    if (cid >= pc_cd.contigs.size()) continue;
                    uint16_t clp = (oi < pc_cd.read_clip.size()) ? pc_cd.read_clip[oi] : 0;
                    int rl = (int)seqs[oi].size() - (int)clp;
                    for (int j = 0; j < rl; ++j) {
                        uint32_t pp = pos + (uint32_t)j;
                        if (pp < pcov_own[cid].size() && pcov_own[cid][pp] < 60000) ++pcov_own[cid][pp];
                    }
                }
            }
            const std::vector<std::vector<uint16_t>>& pcov =
                (&pc_cd == &cdb) ? cov : pcov_own;
            _iplap("  pe: pcov build");
            // cluster reads by right-context anchor that is UNIQUE in the contigs
            //
            // Only anchors present in pkidx are stored. The sole consumer is the loop over
            // pkidx below, which looks up rc_reads.find(kv.first), so a key absent from pkidx
            // is never read; most read 25-mers (including error k-mers) are such keys. This
            // is output-identical by construction, and pkidx is complete before this point.
            // ── BITSET PRE-FILTER ON THE ANCHOR PROBE ───────────────────────
            // Each read probes about 85 k-mers against pkidx, and almost all probes miss.
            //
            // A Bloom-style bitset answers "definitely absent" without touching
            // the map. It can produce FALSE POSITIVES but never false negatives,
            // so a survivor is still checked against pkidx exactly as before --
            // the filter can only skip lookups that would have missed. Output is
            // therefore identical BY CONSTRUCTION, not by measurement.
            //
            // Sized at 16 bits per anchor (~43 MB for 21.6M anchors) with two
            // independent hashes: ~0.4% false-positive rate, so ~99.6% of the
            // ~340M probes end at a single cache line instead of a hash walk.
            pkidx.build_dir();
            fprintf(stderr, "[PCLUSTER] anchor directory %zu buckets + %zu entries (%.1f MB) for %zu anchors\n",
                    pkidx.dir.size() - 1, pkidx.buck.size(),
                    (double)(pkidx.dir.size() * 4 + pkidx.buck.size() * 12) / 1048576.0,
                    pkidx.size());
            _iplap("  pe: anchor directory");
            // ── PARALLEL READ SCAN, FUSED WITH THE ANCHOR WORK ──────────────────
            // Reads are independent, so each thread scans its own reads and does the indel
            // comparison inline as each anchor hit is found, with per-thread accumulators
            // merged afterwards. No per-anchor grouping is built. That grouping existed only
            // for a 200-hits-per-anchor cap, which never binds on contig-unique 25-mers at 30x
            // (longest row measured: 70 hits). With the cap inert nothing downstream depends
            // on hit order: pvotes and panch are set inserts and ploc keeps the smallest
            // (contig,pos), so all three merge by union or minimum.
            static long g_pc_identical=0, g_pc_nogap=0, g_pc_found=0, g_pc_seen=0;
            // (contig, pos, signed gap, inserted seq) -> supporting reads
            // DISTINCT reads per event, not votes. A read spanning an indel
            // matches many anchors (every unique 25-mer in its right context),
            // so a raw vote count multiplies one read into dozens and is not a
            // read count at all.
            std::map<std::tuple<std::string,int,std::string>, std::unordered_set<uint32_t>> pvotes;
            std::map<std::tuple<std::string,int,std::string>, std::pair<uint32_t,uint32_t>> ploc;
            // Distinct ANCHORS backing each event -- independent right-context
            // anchors are independent evidence, the same standard the bubble
            // channel applies via MIN_ANCH.
            std::map<std::tuple<std::string,int,std::string>, std::unordered_set<uint32_t>> panch;
            int pc_nt = 1;
            #ifdef _OPENMP
            pc_nt = omp_get_max_threads();
            #endif
            std::vector<std::map<std::tuple<std::string,int,std::string>, std::unordered_set<uint32_t>>> tvotes((size_t)pc_nt), tanch((size_t)pc_nt);
            std::vector<std::map<std::tuple<std::string,int,std::string>, std::pair<uint32_t,uint32_t>>> tloc((size_t)pc_nt);
            // getenv is read once here, not per hit. The counters it guards are
            // debug-only and racy, so they stay behind it.
            const bool PCDBG = std::getenv("CAPS_PCDBG") != nullptr;
            #pragma omp parallel for schedule(static)
            for (long long i = 0; i < (long long)seqs.size(); ++i) {
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif
                auto& pvotes_t = tvotes[(size_t)tid];
                auto& panch_t  = tanch[(size_t)tid];
                auto& ploc_t   = tloc[(size_t)tid];
                const std::string& q = seqs[(size_t)i];
                // Reverse complement of this read, built at most ONCE per read.
                std::string qq_rc; bool rc_built = false;
                Roll25 rw; const size_t qn = q.size();
                for (size_t b3 = LW; b3 < qn; ++b3) {
                    rw.push(q[b3]);
                    if (b3 + 1 < (size_t)LW + (size_t)AK) continue;
                    if (!rw.ok()) continue;
                    const size_t j = b3 + 1 - (size_t)AK;
                    const size_t ai = pkidx.lookup(rw.canon());
                    if (ai == (size_t)-1) continue;   // provably dead otherwise
                    const uint32_t ccid = pkidx.v[ai].first;
                    const uint32_t cpos = pkidx.v[ai].second;
                    if (cpos < (uint32_t)LW) continue;
                    const std::string& cc2 = pc_cd.contigs[ccid];
                    if (PCDBG) ++g_pc_seen;
                    // Orientation comes straight from the rolling window:
                    // rw.f/rw.r ARE the forward k-mer at j and its reverse
                    // complement (pack25 + rc25 at this position).
                    const std::string* qqp = &q;
                    uint32_t qp = (uint32_t)j;
                    if (rw.f > rw.r) {
                        if (!rc_built) { qq_rc = rc_str(q); rc_built = true; }
                        qqp = &qq_rc;
                        qp = (uint32_t)(q.size() - j - AK);
                    }
                    const std::string& qq = *qqp;
                    if (qp < (uint32_t)LW) continue;
                    int d = 0;
                    while (d < LW && cc2[cpos - 1 - d] == qq[qp - 1 - d]) ++d;
                    if (d >= LW) { if (PCDBG) ++g_pc_identical; continue; }   // identical: no event
                    int best_g = 0; std::string best_ins;
                    for (int g = 1; g <= MAXINDEL && !best_g; ++g) {
                        if (cpos >= (uint32_t)(d + g + 10) && qp >= (uint32_t)(d + 10)) {
                            bool ok = true;
                            for (int t = 0; t < 10; ++t)
                                if (cc2[cpos - 1 - d - g - t] != qq[qp - 1 - d - t]) { ok = false; break; }
                            if (ok) { best_g = g; }
                        }
                        if (!best_g && qp >= (uint32_t)(d + g + 10) && cpos >= (uint32_t)(d + 10)) {
                            bool ok = true;
                            for (int t = 0; t < 10; ++t)
                                if (cc2[cpos - 1 - d - t] != qq[qp - 1 - d - g - t]) { ok = false; break; }
                            if (ok) { best_g = -g; best_ins = qq.substr(qp - d - g, (size_t)g); }
                        }
                    }
                    if (!best_g) { if (PCDBG) ++g_pc_nogap; continue; }
                    if (PCDBG) ++g_pc_found;
                    const uint32_t apos2 = cpos - (uint32_t)d;   // contig position of the event
                    if (apos2 == 0) continue;
                    std::string ctx;
                    {
                        size_t lo = (apos2 >= 12u) ? (size_t)apos2 - 12u : 0u;
                        size_t hi = std::min(cc2.size(), (size_t)apos2 + 12u);
                        if (hi > lo) ctx = cc2.substr(lo, hi - lo);
                    }
                    auto pkey = std::make_tuple(ctx, best_g, best_ins);
                    pvotes_t[pkey].insert((uint32_t)i);
                    panch_t[pkey].insert(cpos);
                    {
                        auto lit = ploc_t.find(pkey);
                        const std::pair<uint32_t,uint32_t> cand(ccid, apos2);
                        if (lit == ploc_t.end()) ploc_t.emplace(pkey, cand);
                        else if (cand < lit->second) lit->second = cand;
                    }
                }
            }
            _iplap("  pe: read scan + anchor work");
            // MERGE. pvotes/panch are set unions (order-free). ploc keeps the
            // smallest (contig,pos), so the merged value is the minimum over
            // all threads -- identical to the serial minimum regardless of how
            // anchors were distributed.
            for (int t = 0; t < pc_nt; ++t) {
                for (auto& kv : tvotes[(size_t)t]) pvotes[kv.first].insert(kv.second.begin(), kv.second.end());
                for (auto& kv : tanch[(size_t)t])  panch[kv.first].insert(kv.second.begin(), kv.second.end());
                for (auto& kv : tloc[(size_t)t]) {
                    auto it = ploc.find(kv.first);
                    if (it == ploc.end()) ploc.emplace(kv.first, kv.second);
                    else if (kv.second < it->second) it->second = kv.second;
                }
            }
            { decltype(tvotes)().swap(tvotes); decltype(tanch)().swap(tanch); decltype(tloc)().swap(tloc);
              }
            if (std::getenv("CAPS_PCDBG")) {
                std::map<int,int> hist;
                for (auto& kv : pvotes) hist[(int)kv.second.size()]++;
                fprintf(stderr,"[pc] distinct events=%zu  support histogram:", pvotes.size());
                for (auto& h : hist) if (h.first<=8) fprintf(stderr," %dread(s)x%d", h.first, h.second);
                fprintf(stderr,"\n");
            }
            long pc_drop_min=0, pc_drop_aflo=0, pc_drop_afhi=0, pc_emit=0;
            _iplap("pcluster: pre-emit");
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
                        // respect to the contig -- a real sequence difference,
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
            _iplap("pcluster: scan + emit");
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
            _iplap("pcluster: teardown");
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
            _iplap("lvotes: read scan");
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
            _iplap("lvotes: coverage");
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

        _iplap("pcluster: emit loop");
        _iplap("XSNV: second kidx walk");
        // ── 6b2. Cross-contig SNV pass (default ON for CAPSULE; see the
        // struct-level comment on SnvBubble for why this is the primary
        // signal here rather than an experimental extra) ──
        size_t n_xsnv = 0;
        _iplap("gap before XSNV");
        if (!std::getenv("CAPS_NO_XSNV")) {
            constexpr int XSNV_FLANK = 15, XSNV_MIN_ANCHORS = 15;
            std::unordered_set<uint64_t> pileup_keys;
            pileup_keys.reserve(n_snv);
            for (size_t i = 0; i < n_snv; ++i)
                pileup_keys.insert(colkey(orecs[i].cid, orecs[i].pos));
            struct XsnvAgg { int anchors = 0; uint32_t altcid = 0, altpos = 0; };
            // PARALLEL OVER RUNS, FILTERED FIRST. This walks all ~65M kidx
            // entries and the body bails immediately unless a run has EXACTLY
            // two occurrences, so the overwhelming majority are no-ops. That is
            // the same failure mode that made a first attempt at the anchor
            // scan deliver 2.9%: threading 86M trivial tasks costs more in
            // scheduling than the work is worth. Collect the runs that carry
            // work, then thread only those.
            //
            // `xim` accumulates anchors (order-free) but overwrites
            // altcid/altpos (last writer). Each thread keeps its own map and
            // they merge afterwards, taking the entry from the LOWEST original
            // run index so the result does not depend on scheduling.
            std::vector<std::pair<size_t,size_t>> xruns;
            for (size_t r0 = 0; r0 < kidx.size(); ) {
                size_t r1 = r0;
                while (r1 < kidx.size() && kidx[r1].kmer() == kidx[r0].kmer()) ++r1;
                if (r1 - r0 == 2) xruns.push_back({r0, r1});
                r0 = r1;
            }
            int x_nt = 1;
            #ifdef _OPENMP
            x_nt = omp_get_max_threads();
            #endif
            std::vector<std::map<std::tuple<uint32_t,uint32_t,char,char>, std::pair<XsnvAgg,size_t>>> txim((size_t)x_nt);
            fprintf(stderr, "[XSNV] %zu runs of exactly 2 over %d threads\n", xruns.size(), x_nt);
            #pragma omp parallel for schedule(dynamic, 256)
            for (long long xi = 0; xi < (long long)xruns.size(); ++xi) {
                int xtid = 0;
                #ifdef _OPENMP
                xtid = omp_get_thread_num();
                #endif
                auto& xim = txim[(size_t)xtid];
                const size_t base = xruns[(size_t)xi].first;
                const size_t myx  = (size_t)xi;
                uint32_t cax = kidx[base].ci,   pax = kidx[base].pos;   uint8_t oax = kidx[base].orient();
                uint32_t cbx = kidx[base+1].ci, pbx = kidx[base+1].pos; uint8_t obx = kidx[base+1].orient();
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
                auto& xe = xim[std::make_tuple(rcx, sb.apos, sb.ref_base, sb.alt_base)];
                xe.first.anchors++;
                if (xe.first.anchors == 1 || myx < xe.second) {
                    xe.first.altcid = acx; xe.first.altpos = alt_pos_fwd; xe.second = myx;
                }
            }
            // merge: anchors sum, altcid/altpos from the lowest run index
            std::map<std::tuple<uint32_t,uint32_t,char,char>, XsnvAgg> xim;
            {
                std::map<std::tuple<uint32_t,uint32_t,char,char>, size_t> xrun_of;
                for (int t = 0; t < x_nt; ++t)
                    for (auto& kv : txim[(size_t)t]) {
                        auto it = xim.find(kv.first);
                        if (it == xim.end()) { xim.emplace(kv.first, kv.second.first); xrun_of[kv.first] = kv.second.second; }
                        else {
                            it->second.anchors += kv.second.first.anchors;
                            if (kv.second.second < xrun_of[kv.first]) {
                                it->second.altcid = kv.second.first.altcid;
                                it->second.altpos = kv.second.first.altpos;
                                xrun_of[kv.first] = kv.second.second;
                            }
                        }
                    }
                decltype(txim)().swap(txim);
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
        // Benchmark-only (CAPS_DUMP_CONTIGS, so the lift can map calls to the genome).
        // About a million contigs are written; a large stdio buffer (setvbuf) keeps
        // this from distorting the timing of the run.
        static std::vector<char> dumpbuf;
        if (df) { dumpbuf.resize(1u << 22); setvbuf(df, dumpbuf.data(), _IOFBF, dumpbuf.size()); }
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
    fprintf(stderr, "[INDEL-PROF] %-24s %8.2fs\n", "6b2 cross-contig SNV",
            elapsed_s(g_ip_last, clk::now()));
    g_ipsum += elapsed_s(g_ip_last, clk::now());
    fprintf(stderr, "[CAPS-CALL] contigs=%zu H=%u candidates=%zu SNVs=%zu indels=%zu -> %s\n",
            cd.contigs.size(), H, total_candidates, n_snv, n_indel, out_vcf.c_str());
    // RECONCILIATION. Six timers once summed to 75.6 s inside a stage that ran
    // 1000+ s, two numbers that contradict each other. A stage total that does not
    // match its parts means the parts are mis-placed, so the check is printed on
    // every run.
    fprintf(stderr, "[INDEL-PROF] === accounted %.2fs of stage; if this differs "
                    "from the indel_pass total below, the timers are WRONG ===\n",
            g_ipsum);
    phase("indel_pass", t_mark);
    fprintf(stderr, "[CAPS-CALL-TIMING] %-16s %8.3fs\n", "TOTAL", elapsed_s(t_start, clk::now()));

    return (int)orecs.size();
}

} // namespace capscall
