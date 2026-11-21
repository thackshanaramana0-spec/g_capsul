// caps_fmi.h -- FM-index over the reads the pseudogenome cannot represent.
//
// WHY THIS EXISTS
// ---------------
// `query` answers "which reads CONTAIN this string" by searching the
// pseudogenome (pg). That search is complete for every read whose emitted
// sequence is a pure pg substring -- a read with no deviations and no N is
// emitted byte-for-byte as pg[pos, pos+len), so any string it contains also
// occurs in pg inside the read's span, and the pg search finds it.
//
// The only reads the pg search can miss are the others: reads carrying
// deviations from the consensus, and reads carrying N (pg never holds an N).
// Compression has already separated the reads into these two groups, so this
// index covers exactly the second group and nothing else. Together with the pg
// search the answer is complete BY CONSTRUCTION -- no tolerance, no seed length
// and no probe-length precondition is involved -- and the index grows with the
// number of reads compression could not absorb, not with the number of reads.
//
// It is built AFTER compression, from the archive alone, and is a separate
// file. The archive itself is not changed.
//
// WHAT IS STORED
// --------------
// Text   T = r_0 # r_1 # ... # r_{m-1} $   over {$ < # < N < A < C < G < T}.
// BWT    2 bits per row for A/C/G/T. The rare non-ACGT rows ($, #, N) are kept
//        as a sorted exception list, so the bulk stays at 2 bits per base.
// Occ    rank checkpoints every 512 rows (0.5 bit per row).
// Locate a read is identified by walking LF back to its start: the row whose
//        BWT symbol is '#' or '$' is the read's first suffix, and a per-'#'-row
//        table names the read. Walks are bounded by the read length, and for
//        long reads by an optional sample every D bases inside a read, so no
//        suffix-array samples are needed.
//
// Construction is linear time: the suffix array is built with SA-IS
// (Nong, Zhang and Chan, 2009) and discarded once the BWT is taken.
//
// VERBATIM READS
// --------------
// A read the archive reproduces exactly but that `query` cannot emit as
// itself (several originals sharing one unique placement with different
// lengths, N patterns or deviations) is stored verbatim in the same file and
// indexed like any other. Its identity carries VERBATIM_FLAG, so a caller
// tells it from a placement-keyed read and emits the stored sequence. The
// decoder stores it in pseudogenome orientation with its own pseudogenome
// start, exactly as every other record `query` emits.
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>

namespace capsfmi {

// ---------------------------------------------------------------------------
// Bit-packed integer vector. Every array of row numbers or read numbers is
// stored at exactly the width its largest value needs.
// ---------------------------------------------------------------------------
class IntVec {
public:
    static unsigned bits_for(uint64_t maxv) {
        unsigned w = 1;
        while (w < 64 && (maxv >> w)) ++w;
        return w;
    }
    void init(size_t n, unsigned width) {
        n_ = n; w_ = width ? width : 1;
        words_.assign((n_ * w_ + 63) / 64 + 1, 0);
    }
    size_t size() const { return n_; }
    void set(size_t i, uint64_t v) {
        const uint64_t mask = (w_ == 64) ? ~0ull : ((1ull << w_) - 1);
        v &= mask;
        const size_t b = i * w_, k = b >> 6; const unsigned o = (unsigned)(b & 63);
        words_[k] = (words_[k] & ~(mask << o)) | (v << o);
        if (o + w_ > 64) {
            const unsigned hi = o + w_ - 64;
            words_[k + 1] = (words_[k + 1] & ~((1ull << hi) - 1)) | (v >> (64 - o));
        }
    }
    uint64_t get(size_t i) const {
        const uint64_t mask = (w_ == 64) ? ~0ull : ((1ull << w_) - 1);
        const size_t b = i * w_, k = b >> 6; const unsigned o = (unsigned)(b & 63);
        uint64_t x = words_[k] >> o;
        if (o + w_ > 64) x |= words_[k + 1] << (64 - o);
        return x & mask;
    }
    // number of elements < x, for a vector sorted ascending
    size_t rank_lt(uint64_t x) const {
        size_t lo = 0, hi = n_;
        while (lo < hi) { const size_t mid = lo + (hi - lo) / 2; if (get(mid) < x) lo = mid + 1; else hi = mid; }
        return lo;
    }
    size_t bytes() const { return words_.size() * 8 + 16; }
    bool write(FILE* f) const {
        const uint64_t n = n_, w = w_, nw = words_.size();
        return fwrite(&n, 8, 1, f) == 1 && fwrite(&w, 8, 1, f) == 1 && fwrite(&nw, 8, 1, f) == 1 &&
               (nw == 0 || fwrite(words_.data(), 8, nw, f) == nw);
    }
    bool read(FILE* f) {
        uint64_t n = 0, w = 0, nw = 0;
        if (fread(&n, 8, 1, f) != 1 || fread(&w, 8, 1, f) != 1 || fread(&nw, 8, 1, f) != 1) return false;
        if (w < 1 || w > 64 || nw != (n * w + 63) / 64 + 1) return false;
        n_ = (size_t)n; w_ = (unsigned)w; words_.resize((size_t)nw);
        return nw == 0 || fread(words_.data(), 8, nw, f) == nw;
    }
private:
    size_t n_ = 0; unsigned w_ = 1; std::vector<uint64_t> words_;
};

// ---------------------------------------------------------------------------
// SA-IS suffix array construction (Nong, Zhang and Chan, 2009).
// s[0..n) over alphabet [0,K), n >= 1. Preconditions the caller must meet:
// s[n-1] is the unique smallest symbol (the text ends in SYM_END), every symbol
// is < K, and I can represent n. Writes the suffix array to SA[0..n). Linear
// time. Besides SA, each recursion level holds one bit per position and a
// bucket table of K entries, and the recursion is on a problem at most half the size.
// ---------------------------------------------------------------------------
namespace detail {

template <typename C, typename I>
void sais(const C* s, I* SA, I n, I K) {
    const I EMPTY = std::numeric_limits<I>::max();
    if (n == 1) { SA[0] = 0; return; }
    std::vector<bool> t((size_t)n);          // true = S-type
    t[(size_t)n - 1] = true;
    t[(size_t)n - 2] = false;                 // s[n-1] is the unique minimum
    for (int64_t i = (int64_t)n - 3; i >= 0; --i)
        t[(size_t)i] = s[i] < s[i + 1] || (s[i] == s[i + 1] && t[(size_t)i + 1]);
    auto lms = [&](I i) -> bool { return i > 0 && i != EMPTY && t[(size_t)i] && !t[(size_t)i - 1]; };

    std::vector<I> bkt((size_t)K);
    auto buckets = [&](bool end) {
        std::fill(bkt.begin(), bkt.end(), (I)0);
        for (I i = 0; i < n; ++i) ++bkt[(size_t)s[i]];
        I sum = 0;
        for (size_t k = 0; k < (size_t)K; ++k) { const I c = bkt[k]; sum += c; bkt[k] = end ? sum : (I)(sum - c); }
    };
    auto induce = [&]() {
        buckets(false);
        for (I i = 0; i < n; ++i) {
            const I p = SA[i];
            if (p != EMPTY && p > 0 && !t[(size_t)p - 1]) SA[bkt[(size_t)s[p - 1]]++] = p - 1;
        }
        buckets(true);
        for (int64_t i = (int64_t)n - 1; i >= 0; --i) {
            const I p = SA[i];
            if (p != EMPTY && p > 0 && t[(size_t)p - 1]) SA[--bkt[(size_t)s[p - 1]]] = p - 1;
        }
    };

    // Stage 1: sort the LMS substrings.
    buckets(true);
    std::fill(SA, SA + n, EMPTY);
    for (I i = 1; i < n; ++i) if (lms(i)) SA[--bkt[(size_t)s[i]]] = i;
    induce();

    I n1 = 0;
    for (I i = 0; i < n; ++i) if (lms(SA[i])) SA[n1++] = SA[i];
    std::fill(SA + n1, SA + n, EMPTY);
    I name = 0, prev = EMPTY;
    for (I i = 0; i < n1; ++i) {
        const I pos = SA[i];
        bool diff = false;
        for (I d = 0; d < n; ++d) {
            if (prev == EMPTY || s[pos + d] != s[prev + d] || t[(size_t)(pos + d)] != t[(size_t)(prev + d)]) { diff = true; break; }
            if (d > 0 && (lms(pos + d) || lms(prev + d))) break;
        }
        if (diff) { ++name; prev = pos; }
        SA[n1 + pos / 2] = name - 1;
    }
    for (int64_t i = (int64_t)n - 1, j = (int64_t)n - 1; i >= (int64_t)n1; --i)
        if (SA[i] != EMPTY) SA[j--] = SA[i];

    // Stage 2: sort the reduced string, recursively if names are not unique.
    I* s1 = SA + n - n1;
    if (name < n1) sais<I, I>(s1, SA, n1, name);
    else for (I i = 0; i < n1; ++i) SA[s1[i]] = i;

    // Stage 3: induce the full order from the sorted LMS suffixes.
    buckets(true);
    for (I i = 1, j = 0; i < n; ++i) if (lms(i)) s1[j++] = i;
    for (I i = 0; i < n1; ++i) SA[i] = s1[SA[i]];
    std::fill(SA + n1, SA + n, EMPTY);
    for (int64_t i = (int64_t)n1 - 1; i >= 0; --i) {
        const I p = SA[i]; SA[i] = EMPTY;
        SA[--bkt[(size_t)s[p]]] = p;
    }
    induce();
}

inline unsigned popc(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return (unsigned)__builtin_popcountll(x);
#else
    x = x - ((x >> 1) & 0x5555555555555555ull);
    x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
    return (unsigned)((((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Full) * 0x0101010101010101ull) >> 56);
#endif
}
// one bit (at the low bit of each 2-bit lane) per symbol equal to `code`
inline uint64_t lanes_equal(uint64_t w, unsigned code) {
    const uint64_t x = w ^ (0x5555555555555555ull * code);
    return ~(x | (x >> 1)) & 0x5555555555555555ull;
}

} // namespace detail

// Symbol codes in the indexed text. The order is the sort order.
enum : uint8_t { SYM_END = 0, SYM_SEP = 1, SYM_N = 2, SYM_A = 3, SYM_C = 4, SYM_G = 5, SYM_T = 6, SIGMA = 7 };

inline int sym_of(char ch) {
    switch (ch) {
        case 'A': case 'a': return SYM_A;
        case 'C': case 'c': return SYM_C;
        case 'G': case 'g': return SYM_G;
        case 'T': case 't': return SYM_T;
        case 'N': case 'n': return SYM_N;
        default: return -1;
    }
}

// FM-index over a set of reads. Typical use:
//   FMIndex fm;  fm.build(reads, ids, D, err);  fm.save(path, err);     // at index time
//   FMIndex fm;  fm.load(path, err);  fm.find(probe, hits, err);        // at query time
// find() reports the identity (ids[i] given to build) of every indexed read that
// contains the pattern, one entry per occurrence, so a read containing the
// pattern twice is reported twice. Callers deduplicate.
//
// Thread safety: after build() or load() the object is never modified, and
// find() is const and uses no shared scratch, so concurrent find() calls are safe.
//
// Memory: about 2 bits per indexed base for the BWT, plus 0.5 bit per base for
// rank checkpoints, plus the read identities and the (usually empty) locate samples.
class FMIndex {
public:
    static constexpr uint64_t MAGIC = 0x31494D4653504143ull;   // "CAPSFMI1" little-endian
    static constexpr uint32_t VERSION = 2;
    static constexpr uint32_t VERBATIM_FLAG = 0x80000000u;     // identity bit: stored sequence, not a placement
    static constexpr unsigned BLOCK_SHIFT = 9;                 // checkpoint every 512 rows

    // Build over reads[i] (alphabet ACGTN, any other byte is indexed as N) and
    // remember ids[i] as the identity a hit on reads[i] reports. Empty reads
    // contain nothing and are skipped. sample_interval D > 0 bounds every
    // locate walk to D steps; reads no longer than D add no samples at all.
    bool build(const std::vector<std::string>& reads, const std::vector<uint32_t>& ids,
               uint32_t sample_interval, std::string& err) {
        if (reads.size() != ids.size()) { err = "reads/ids size mismatch"; return false; }
        auto vs = std::move(vb_seq_); auto vp = std::move(vb_place_);
        clear();
        vb_seq_ = std::move(vs); vb_place_ = std::move(vp);
        D_ = sample_interval;
        std::vector<size_t> keep;
        uint64_t bases = 0;
        for (size_t i = 0; i < reads.size(); ++i)
            if (!reads[i].empty()) { keep.push_back(i); bases += reads[i].size(); if (reads[i].size() > maxlen_) maxlen_ = (uint32_t)std::min<size_t>(reads[i].size(), UINT32_MAX); }
        m_ = keep.size();
        if (m_ == 0) { n_ = 0; bwt_.assign(1, 0); ckpt_.assign(4, 0); return true; }
        n_ = bases + m_;

        std::vector<uint8_t> T((size_t)n_);
        std::vector<uint64_t> starts(m_);
        uint32_t maxid = 0;
        for (size_t k = 0; k < m_; ++k) maxid = std::max(maxid, ids[keep[k]]);
        docid_.init(m_, IntVec::bits_for(maxid));
        size_t p = 0;
        for (size_t k = 0; k < m_; ++k) {
            const std::string& r = reads[keep[k]];
            starts[k] = p;
            docid_.set(k, ids[keep[k]]);
            for (char ch : r) { const int c = sym_of(ch); T[p++] = (uint8_t)(c < 0 ? SYM_N : c); }
            T[p++] = (k + 1 == m_) ? SYM_END : SYM_SEP;
        }

        if (n_ < (uint64_t)std::numeric_limits<uint32_t>::max()) {
            std::vector<uint32_t> SA((size_t)n_);
            detail::sais<uint8_t, uint32_t>(T.data(), SA.data(), (uint32_t)n_, SIGMA);
            from_sa(T, SA, starts);
        } else {
            std::vector<uint64_t> SA((size_t)n_);
            detail::sais<uint8_t, uint64_t>(T.data(), SA.data(), n_, SIGMA);
            from_sa(T, SA, starts);
        }
        return true;
    }

    // Identities of every indexed read that contains `pat` exactly, on the
    // strand given. Callers search the reverse complement separately. A pattern
    // character outside A, C, G, T, N (either case) matches nothing, and an N in
    // the pattern matches an N in a read. An empty pattern returns no hits.
    bool find(const std::string& pat, std::vector<uint32_t>& out, std::string& err) const {
        if (n_ == 0 || pat.empty()) return true;
        uint64_t sp = 0, ep = n_;
        for (size_t k = pat.size(); k-- > 0;) {
            const int c = sym_of(pat[k]);
            if (c < 0) return true;                       // symbol never indexed: no read contains it
            sp = C_[c] + occ((unsigned)c, sp);
            ep = C_[c] + occ((unsigned)c, ep);
            if (sp >= ep) return true;
        }
        for (uint64_t row = sp; row < ep; ++row) {
            uint64_t doc = 0;
            if (!locate_doc(row, doc, err)) return false;
            out.push_back((uint32_t)docid_.get((size_t)doc));
        }
        return true;
    }

    // File layout, host byte order (little-endian on every supported platform):
    //   magic u64, version u32 (2), n u64 (text length), m u64 (indexed reads),
    //   maxlen u32 (longest indexed read), D u32 (locate sample interval),
    //   dollar_row u64 (row whose BWT symbol is SYM_END), C[0..SIGMA] u64 x 8,
    //   nb u64 then nb u64 words of 2-bit BWT codes (32 rows per word),
    //   nc u64 then nc u64 rank counters (four per 512-row block),
    //   seven integer vectors, each [length u64, bit width u64, word count u64, words]:
    //     exception rows, exception symbols, '#' rank at each exception,
    //     read of each '#' row, sampled rows, read of each sampled row, caller identities,
    //   nv u64 then nv verbatim reads, each [position u32, length u32, bytes].
    bool save(const std::string& path, std::string& err) const {
        FILE* f = fopen(path.c_str(), "wb");
        if (!f) { err = "cannot write " + path; return false; }
        bool ok = fwrite(&MAGIC, 8, 1, f) == 1 && fwrite(&VERSION, 4, 1, f) == 1 &&
                  fwrite(&n_, 8, 1, f) == 1 && fwrite(&m_, 8, 1, f) == 1 &&
                  fwrite(&maxlen_, 4, 1, f) == 1 && fwrite(&D_, 4, 1, f) == 1 &&
                  fwrite(&dollar_row_, 8, 1, f) == 1 && fwrite(C_, 8, SIGMA + 1, f) == SIGMA + 1;
        const uint64_t nb = bwt_.size(), nc = ckpt_.size();
        ok = ok && fwrite(&nb, 8, 1, f) == 1 && (nb == 0 || fwrite(bwt_.data(), 8, nb, f) == nb);
        ok = ok && fwrite(&nc, 8, 1, f) == 1 && (nc == 0 || fwrite(ckpt_.data(), 8, nc, f) == nc);
        ok = ok && exc_rows_.write(f) && exc_type_.write(f) && sep_before_.write(f) &&
             sep_doc_.write(f) && smp_rows_.write(f) && smp_doc_.write(f) && docid_.write(f);
        const uint64_t nv = vb_seq_.size();
        ok = ok && fwrite(&nv, 8, 1, f) == 1;
        for (size_t i = 0; ok && i < vb_seq_.size(); ++i) {
            const uint32_t len = (uint32_t)vb_seq_[i].size();
            ok = fwrite(&vb_place_[i], 4, 1, f) == 1 && fwrite(&len, 4, 1, f) == 1 &&
                 (len == 0 || fwrite(vb_seq_[i].data(), 1, len, f) == len);
        }
        ok = (fclose(f) == 0) && ok;
        if (!ok) err = "write failed for " + path;
        return ok;
    }

    bool load(const std::string& path, std::string& err) {
        clear();
        FILE* f = fopen(path.c_str(), "rb");
        if (!f) { err = "cannot open " + path; return false; }
        uint64_t magic = 0; uint32_t ver = 0;
        bool ok = fread(&magic, 8, 1, f) == 1 && magic == MAGIC && fread(&ver, 4, 1, f) == 1 && ver == VERSION;
        if (!ok) { fclose(f); err = path + ": not a CAPS FM-index (or unsupported version)"; return false; }
        ok = fread(&n_, 8, 1, f) == 1 && fread(&m_, 8, 1, f) == 1 && fread(&maxlen_, 4, 1, f) == 1 &&
             fread(&D_, 4, 1, f) == 1 && fread(&dollar_row_, 8, 1, f) == 1 && fread(C_, 8, SIGMA + 1, f) == SIGMA + 1;
        uint64_t nb = 0, nc = 0;
        ok = ok && fread(&nb, 8, 1, f) == 1 && nb == n_ / 32 + 1;
        if (ok) { bwt_.resize((size_t)nb); ok = fread(bwt_.data(), 8, nb, f) == nb; }
        ok = ok && fread(&nc, 8, 1, f) == 1 && nc == ((n_ >> BLOCK_SHIFT) + 1) * 4;
        if (ok) { ckpt_.resize((size_t)nc); ok = fread(ckpt_.data(), 8, nc, f) == nc; }
        ok = ok && exc_rows_.read(f) && exc_type_.read(f) && sep_before_.read(f) &&
             sep_doc_.read(f) && smp_rows_.read(f) && smp_doc_.read(f) && docid_.read(f);
        ok = ok && docid_.size() == m_ && exc_type_.size() == exc_rows_.size() &&
             sep_before_.size() == exc_rows_.size() + 1 && smp_doc_.size() == smp_rows_.size();
        uint64_t nv = 0;
        ok = ok && fread(&nv, 8, 1, f) == 1 && nv <= (1ull << 31);
        for (uint64_t i = 0; ok && i < nv; ++i) {
            uint32_t place = 0, len = 0;
            ok = fread(&place, 4, 1, f) == 1 && fread(&len, 4, 1, f) == 1 && len <= (1u << 28);
            if (!ok) break;
            std::string s(len, 'A');
            ok = len == 0 || fread(&s[0], 1, len, f) == len;
            if (ok) { vb_seq_.push_back(std::move(s)); vb_place_.push_back(place); }
        }
        fclose(f);
        if (!ok) { clear(); err = path + ": truncated or inconsistent FM-index"; return false; }
        return true;
    }

    // Verbatim sequences, addressed by (identity & ~VERBATIM_FLAG).
    uint32_t add_verbatim(const std::string& seq, uint32_t place) {
        vb_seq_.push_back(seq); vb_place_.push_back(place);
        return (uint32_t)(vb_seq_.size() - 1) | VERBATIM_FLAG;
    }
    size_t num_verbatim() const { return vb_seq_.size(); }
    const std::string& verbatim_seq(uint32_t id) const { return vb_seq_[id & ~VERBATIM_FLAG]; }
    uint32_t verbatim_place(uint32_t id) const { return vb_place_[id & ~VERBATIM_FLAG]; }

    uint64_t text_length() const { return n_; }
    uint64_t num_reads() const { return m_; }
    uint64_t num_samples() const { return smp_rows_.size(); }
    size_t bytes() const {
        size_t vb = 0; for (const auto& s : vb_seq_) vb += s.size() + 8;
        return vb + 64 + 8 * (SIGMA + 1) + bwt_.size() * 8 + ckpt_.size() * 8 + exc_rows_.bytes() + exc_type_.bytes() +
               sep_before_.bytes() + sep_doc_.bytes() + smp_rows_.bytes() + smp_doc_.bytes() + docid_.bytes();
    }

private:
    uint64_t n_ = 0, m_ = 0, dollar_row_ = 0;
    uint32_t maxlen_ = 0, D_ = 0;
    uint64_t C_[SIGMA + 1] = {0};
    std::vector<uint64_t> bwt_;      // 2-bit A/C/G/T codes, exceptions stored as 0
    std::vector<uint64_t> ckpt_;     // per block, raw counts of codes 0..3 before the block
    IntVec exc_rows_;                // sorted rows whose BWT symbol is $, # or N
    IntVec exc_type_;                // that symbol (SYM_END, SYM_SEP or SYM_N)
    IntVec sep_before_;              // number of '#' rows among exc_rows_[0..k)
    IntVec sep_doc_;                 // read (index into docid_) whose first suffix sits at the k-th '#' row
    IntVec smp_rows_, smp_doc_;      // sorted sampled rows inside long reads, and their read
    IntVec docid_;                   // read index -> caller identity
    std::vector<std::string> vb_seq_;     // verbatim sequences
    std::vector<uint32_t> vb_place_;      // caller-defined position of each verbatim read (decoder: its pg start)

    void clear() {
        n_ = m_ = dollar_row_ = 0; maxlen_ = D_ = 0;
        vb_seq_.clear(); vb_place_.clear();
        std::fill(C_, C_ + SIGMA + 1, 0);
        bwt_.clear(); ckpt_.clear();
        exc_rows_.init(0, 1); exc_type_.init(0, 2); sep_before_.init(1, 1);
        sep_doc_.init(0, 1); smp_rows_.init(0, 1); smp_doc_.init(0, 1); docid_.init(0, 1);
    }

    template <typename I>
    void from_sa(const std::vector<uint8_t>& T, const std::vector<I>& SA, const std::vector<uint64_t>& starts) {
        const size_t n = (size_t)n_;
        auto read_of = [&](uint64_t pos) -> size_t {      // read containing text position pos
            return (size_t)(std::upper_bound(starts.begin(), starts.end(), pos) - starts.begin()) - 1;
        };
        bwt_.assign(n / 32 + 1, 0);
        uint64_t cnt[SIGMA] = {0};
        std::vector<uint64_t> erow; std::vector<uint8_t> etype; std::vector<uint64_t> sdoc;
        std::vector<uint64_t> srow, sdocs;
        const bool sampling = D_ > 0 && maxlen_ > D_;
        for (size_t i = 0; i < n; ++i) {
            const uint64_t p = (uint64_t)SA[i];
            const uint8_t c = p == 0 ? T[n - 1] : T[(size_t)p - 1];
            ++cnt[c];
            if (c >= SYM_A) {
                bwt_[i >> 5] |= (uint64_t)(c - SYM_A) << (2 * (i & 31));
            } else {
                erow.push_back(i); etype.push_back(c);
                if (c == SYM_END) dollar_row_ = i;
                if (c == SYM_SEP) sdoc.push_back(read_of(p));
            }
            if (sampling && T[(size_t)p] >= SYM_N && c != SYM_SEP && c != SYM_END) {
                const size_t r = read_of(p);
                const uint64_t off = p - starts[r];
                if (off % D_ == 0) { srow.push_back(i); sdocs.push_back(r); }
            }
        }
        // C[c] = number of text symbols smaller than c
        C_[0] = 0;
        for (unsigned c = 1; c <= SIGMA; ++c) C_[c] = C_[c - 1] + cnt[c - 1];

        const size_t nblk = (n >> BLOCK_SHIFT) + 1;
        ckpt_.assign(nblk * 4, 0);
        uint64_t run[4] = {0, 0, 0, 0};
        for (size_t b = 0; b < nblk; ++b) {
            for (unsigned c = 0; c < 4; ++c) ckpt_[b * 4 + c] = run[c];
            const size_t lo = b << BLOCK_SHIFT, hi = std::min(n, lo + ((size_t)1 << BLOCK_SHIFT));
            for (size_t i = lo; i < hi; ++i) ++run[(bwt_[i >> 5] >> (2 * (i & 31))) & 3];
        }

        const unsigned wrow = IntVec::bits_for(n_);
        exc_rows_.init(erow.size(), wrow);
        exc_type_.init(etype.size(), 2);
        sep_before_.init(erow.size() + 1, IntVec::bits_for(erow.size()));
        uint64_t seps = 0;
        for (size_t k = 0; k < erow.size(); ++k) {
            exc_rows_.set(k, erow[k]); exc_type_.set(k, etype[k]);
            sep_before_.set(k, seps);
            if (etype[k] == SYM_SEP) ++seps;
        }
        sep_before_.set(erow.size(), seps);
        const unsigned wdoc = IntVec::bits_for(m_ ? m_ - 1 : 0);
        sep_doc_.init(sdoc.size(), wdoc);
        for (size_t k = 0; k < sdoc.size(); ++k) sep_doc_.set(k, sdoc[k]);
        smp_rows_.init(srow.size(), wrow);
        smp_doc_.init(sdocs.size(), wdoc);
        for (size_t k = 0; k < srow.size(); ++k) { smp_rows_.set(k, srow[k]); smp_doc_.set(k, sdocs[k]); }
    }

    // raw count of 2-bit code `code` in BWT rows [0, i), exceptions counted as code 0
    uint64_t occ_raw(unsigned code, uint64_t i) const {
        const size_t blk = (size_t)(i >> BLOCK_SHIFT);
        uint64_t c = ckpt_[blk * 4 + code];
        size_t w = blk << (BLOCK_SHIFT - 5);
        const size_t wend = (size_t)(i >> 5);
        for (; w < wend; ++w) c += detail::popc(detail::lanes_equal(bwt_[w], code));
        const unsigned r = (unsigned)(i & 31);
        if (r) c += detail::popc(detail::lanes_equal(bwt_[wend], code) & ((1ull << (2 * r)) - 1));
        return c;
    }
    // occurrences of text symbol `sym` in BWT rows [0, i)
    uint64_t occ(unsigned sym, uint64_t i) const {
        if (sym >= SYM_A) {
            uint64_t c = occ_raw(sym - SYM_A, i);
            if (sym == SYM_A) c -= exc_rows_.rank_lt(i);
            return c;
        }
        const size_t k = exc_rows_.rank_lt(i);
        const uint64_t seps = sep_before_.get(k);
        const uint64_t ends = (m_ && dollar_row_ < i) ? 1 : 0;
        if (sym == SYM_SEP) return seps;
        if (sym == SYM_END) return ends;
        return (uint64_t)k - seps - ends;                  // SYM_N
    }
    // BWT symbol at `row`; for an exception row also its index in exc_rows_
    unsigned sym_at(uint64_t row, size_t& exc_index) const {
        const unsigned code = (unsigned)((bwt_[(size_t)(row >> 5)] >> (2 * (row & 31))) & 3);
        if (code) return SYM_A + code;
        const size_t k = exc_rows_.rank_lt(row);
        if (k < exc_rows_.size() && exc_rows_.get(k) == row) { exc_index = k; return (unsigned)exc_type_.get(k); }
        return SYM_A;
    }
    bool locate_doc(uint64_t row, uint64_t& doc, std::string& err) const {
        const uint64_t bound = (uint64_t)maxlen_ + 1;
        for (uint64_t steps = 0; steps <= bound; ++steps) {
            size_t k = 0;
            const unsigned s = sym_at(row, k);
            if (s == SYM_END) { doc = 0; return true; }                     // first suffix of read 0
            if (s == SYM_SEP) { doc = sep_doc_.get((size_t)sep_before_.get(k)); return true; }
            if (smp_rows_.size()) {
                const size_t j = smp_rows_.rank_lt(row);
                if (j < smp_rows_.size() && smp_rows_.get(j) == row) { doc = smp_doc_.get(j); return true; }
            }
            row = C_[s] + occ(s, row);
        }
        err = "FM-index locate walk exceeded the longest read -- index is corrupt";
        return false;
    }
};

} // namespace capsfmi
