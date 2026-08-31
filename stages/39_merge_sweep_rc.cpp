// PgRC2's actual overlap method: no index at all.
//
// Every prototype here so far, and every backend in src/, answers "which reads
// overlap" by BUILDING A STRUCTURE over the read text -- a suffix array (11
// bytes/char), a BWT with checkpoints (4.77), a rope (0.4), a seed hash. The
// structures differ by 20x among themselves and all of them are the wrong shape,
// because they are random-access indexes serving a sequential-access algorithm.
// Greedy growth consumes the longest overlap first and never revisits; it does
// not need "the longest overlap for every pair at every length, on demand".
//
// PgRC2 has no suffix array anywhere -- grep divsufsort/libsais/sais over their
// tree returns nothing. GreedySwipingPackedOverlapPseudoGenomeGenerator instead
// sweeps the overlap length downward and, at each length, MERGE-JOINS two sorted
// index lists: reads ordered by prefix against reads ordered by suffix-at-offset.
// Their working set is sortedReadsIdxs + sortedSuffixIdxs + a 256-entry symbol
// queue: about 8 bytes per READ, not per character, and it shrinks as reads are
// consumed (their readsLeft).
//
//   ARCS   SA + LCP + PLCP    11 bytes/char   2,730 MB on 257 Mchar
//   PgRC2  two index arrays    8 bytes/read       ~7 MB on 851k reads
//
// This prototype measures whether that shape actually holds, and whether the
// assembly it produces is as good. Two things are checked: peak RSS, and
// pseudogenome length against the suffix-array path's 16,347,983.
//
// Reads are packed at 2 bits/base so suffix comparison is word-at-a-time rather
// than byte-at-a-time -- that is why their sweep costs 1.5 s where an unpacked
// version of the same loop is O(n*L^2) byte reads and takes 40 s.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <chrono>

static const uint32_t NONE = UINT32_MAX;

static size_t rss_mb() {
    FILE* f = fopen("/proc/self/statm", "r"); long s = 0, r = 0;
    if (f) { if (fscanf(f, "%ld %ld", &s, &r) != 2) r = 0; fclose(f); }
    return (size_t)r * 4096 / 1048576;
}
static size_t peak_mb() {
    FILE* f = fopen("/proc/self/status", "r"); char l[256]; size_t v = 0;
    if (f) { while (fgets(l, sizeof l, f)) if (!strncmp(l, "VmHWM:", 6)) { sscanf(l + 6, "%zu", &v); break; } fclose(f); }
    return v / 1024;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: 39_merge_sweep_rc <in.fq> [minov]\n"); return 1; }
    const uint32_t MINOV = argc > 2 ? (uint32_t)atoi(argv[2]) : 24;
    auto T0 = std::chrono::steady_clock::now();
    auto lap = [&](const char* w) {
        auto t = std::chrono::steady_clock::now();
        fprintf(stderr, "  %-26s %6.2f s   rss %zu MB\n", w,
                std::chrono::duration<double>(t - T0).count(), rss_mb());
        T0 = t;
    };

    // ── load, drop N-reads, dedup ────────────────────────────────────────────
    std::vector<std::string> reads;
    {
        std::ifstream f(argv[1]); std::string a, b, c, d;
        while (std::getline(f, a) && std::getline(f, b) && std::getline(f, c) && std::getline(f, d))
            if (b.find('N') == std::string::npos) reads.push_back(b);
        std::sort(reads.begin(), reads.end());
        reads.erase(std::unique(reads.begin(), reads.end()), reads.end());
    }
    const uint32_t n = (uint32_t)reads.size();
    const uint32_t L = n ? (uint32_t)reads[0].size() : 0;
    for (const auto& r : reads) if (r.size() != L) { fprintf(stderr, "variable read length; this prototype assumes fixed\n"); return 1; }
    fprintf(stderr, "reads=%u  len=%u  minov=%u\n", n, L, MINOV);
    lap("load+dedup");

    // ── pack 2 bits/base, one uint64 per 32 bases ────────────────────────────
    // Both views, laid out as entry = 2*rid + view -- the same interleaving
    // build_apsp_candidates uses. PgRC2 assembles forward-only and recovers the
    // reverse-complement redundancy afterwards by matching the pseudogenome
    // against its own RC; this variant instead makes RC visible DURING the
    // sweep, to separate "the 3x longer pseudogenome is the missing RC" from
    // "the 3x is the greedy itself".
    const uint32_t W = (L + 31) / 32;
    const uint32_t m = n * 2;
    std::vector<uint64_t> pk((size_t)m * W, 0);
    for (uint32_t i = 0; i < n; ++i) {
        const std::string& s = reads[i];
        for (uint32_t j = 0; j < L; ++j) {
            const char ch = s[j];
            const uint64_t v = (ch == 'A') ? 0 : (ch == 'C') ? 1 : (ch == 'G') ? 2 : 3;
            pk[(size_t)(2 * i) * W + j / 32] |= v << (62 - 2 * (j % 32));
            // reverse complement: base j of the RC is complement of base L-1-j
            const uint64_t rv = 3ULL - v;
            const uint32_t rj = L - 1 - j;
            pk[(size_t)(2 * i + 1) * W + rj / 32] |= rv << (62 - 2 * (rj % 32));
        }
    }
    { std::vector<std::string> tmp; tmp.swap(reads); }        // packed form is enough from here
    lap("pack 2bit");

    // base at position j of read i
    auto base_at = [&](uint32_t i, uint32_t j) -> uint32_t {
        return (uint32_t)((pk[(size_t)i * W + j / 32] >> (62 - 2 * (j % 32))) & 3ULL);
    };
    // compare read a's substring [ao, ao+len) against read b's [bo, bo+len)
    auto cmp_sub = [&](uint32_t a, uint32_t ao, uint32_t b, uint32_t bo, uint32_t len) -> int {
        for (uint32_t k = 0; k < len; ++k) {
            const uint32_t x = base_at(a, ao + k), y = base_at(b, bo + k);
            if (x != y) return x < y ? -1 : 1;
        }
        return 0;
    };

    // ── sortedPrefix: reads ordered by their prefix ──────────────────────────
    std::vector<uint32_t> byPrefix(m);
    for (uint32_t i = 0; i < m; ++i) byPrefix[i] = i;
    std::sort(byPrefix.begin(), byPrefix.end(),
              [&](uint32_t x, uint32_t y) { return cmp_sub(x, 0, y, 0, L) < 0; });
    lap("sort by prefix");

    // ── descending sweep with a merge-join at each overlap length ────────────
    std::vector<uint32_t> nxt(m, NONE), prv(m, NONE), ovl(m, 0), head(m), tail(m);
    for (uint32_t i = 0; i < m; ++i) { head[i] = i; tail[i] = i; }
    // A read and its RC are the same read, so a read may appear in the
    // pseudogenome in exactly ONE view -- but it must still be extendable on
    // both sides within that view. Marking the read "used" at link time blocks
    // the second extension too, which cost 64k links and left the pseudogenome
    // twice as long as forward-only. Commit the VIEW instead.
    std::vector<uint8_t> view_of(n, 2);            // 2 = not yet committed

    std::vector<uint32_t> bySuffix;                 // reused every iteration
    bySuffix.reserve(n);
    size_t links = 0;

    for (uint32_t ov = L - 1; ov >= MINOV; --ov) {
        const uint32_t off = L - ov;                // suffix starts here
        // reads that can still take a successor
        bySuffix.clear();
        for (uint32_t i = 0; i < m; ++i) if (nxt[i] == NONE) bySuffix.push_back(i);
        if (bySuffix.empty()) break;
        std::sort(bySuffix.begin(), bySuffix.end(),
                  [&](uint32_t x, uint32_t y) { return cmp_sub(x, off, y, off, ov) < 0; });

        // Merge-join: suffix(A, off, ov) == prefix(B, 0, ov).
        //
        // Both sides carry duplicate keys, so this cannot be a plain two-pointer
        // walk. Advancing past a B that turned out to be unusable (already has a
        // predecessor, or would close a cycle) would hide it from every later A
        // with the same key -- which is not a small effect: it cut links from
        // ~140k to 9,982 and left the pseudogenome barely compressed. PgRC2
        // saves curPreIt and rescans the equal range for exactly this reason.
        size_t ia = 0, ib = 0;
        while (ia < bySuffix.size() && ib < byPrefix.size()) {
            const uint32_t A = bySuffix[ia];
            // skip prefixes that sort before this suffix (ib is monotone overall)
            while (ib < byPrefix.size() && cmp_sub(A, off, byPrefix[ib], 0, ov) > 0) ++ib;
            if (ib >= byPrefix.size()) break;
            // scan the whole equal range from a saved start
            for (size_t j = ib; j < byPrefix.size(); ++j) {
                const uint32_t B = byPrefix[j];
                if (cmp_sub(A, off, B, 0, ov) != 0) break;
                if (A == B || (A >> 1) == (B >> 1)) continue;      // self or own RC
                if (nxt[A] != NONE || prv[B] != NONE || head[A] == B) continue;
                const uint8_t va = view_of[A >> 1], vb = view_of[B >> 1];
                if (va != 2 && va != (A & 1)) continue;            // committed to the other view
                if (vb != 2 && vb != (B & 1)) continue;
                nxt[A] = B; prv[B] = A; ovl[A] = ov;
                view_of[A >> 1] = (uint8_t)(A & 1);
                view_of[B >> 1] = (uint8_t)(B & 1);
                const uint32_t h = head[A], t = tail[B];
                tail[h] = t; head[t] = h;
                ++links;
                break;
            }
            ++ia;
        }
        if (ov == MINOV) break;
    }
    lap("descending merge sweep");
    fprintf(stderr, "links=%zu\n", links);

    // ── emit pseudogenome ────────────────────────────────────────────────────
    size_t pg_len = 0;
    std::vector<uint8_t> emitted(n, 0);
    // Real chains first. An entry with no predecessor AND no successor is not a
    // chain -- it is usually the unused view of a read that is already placed
    // via its other view, and emitting those in index order before the real
    // chains double-counted most of the reads (pg 112 MB instead of ~48 MB).
    for (uint32_t i = 0; i < m; ++i) {
        if (prv[i] != NONE || nxt[i] == NONE) continue;
        if (emitted[i >> 1]) continue;
        uint32_t cur = i;
        pg_len += L; emitted[cur >> 1] = 1;
        while (nxt[cur] != NONE) { pg_len += L - ovl[cur]; cur = nxt[cur]; emitted[cur >> 1] = 1; }
    }
    // reads that never entered any chain still have to be stored
    for (uint32_t r = 0; r < n; ++r) if (!emitted[r]) pg_len += L;
    lap("emit");

    fprintf(stderr, "\nPG_LEN %zu   (suffix-array path: 16,347,983)\n", pg_len);
    fprintf(stderr, "PEAK   %zu MB  (suffix-array path: 3,602 MB process / 2,730 MB structures)\n", peak_mb());
    printf("PG_LEN %zu PEAK_MB %zu\n", pg_len, peak_mb());
    return 0;
}
