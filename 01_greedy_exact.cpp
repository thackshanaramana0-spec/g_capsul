// Textbook greedy shortest-common-superstring assembler.
//
// Written from the algorithm (Tarhio & Ukkonen 1988; Blum et al. 1994), not from
// any existing implementation: take the pair with the largest overlap, merge,
// repeat. Realised the standard way — sweep the overlap length L downwards, and
// at each L pair every still-open tail with a still-open head sharing that
// L-length overlap.
//
// Cycles are AVOIDED rather than repaired: a link is refused when it would close
// the chain onto itself. Each chain keeps a head pointer at its tail and a tail
// pointer at its head, so the check and the update are both O(1) and no cycle is
// ever created. (The alternative — build cycles then break each at its weakest
// link — is what PgRC2 does, and its implementation shadows the running minimum
// inside the search loop, so it always breaks at the scan's entry point instead
// of the smallest overlap.)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <unordered_map>
#include <algorithm>

static const uint32_t NONE = UINT32_MAX;

static inline uint64_t fnv(const char* p, uint32_t n) {
    uint64_t h = 1469598103934665603ULL;
    for (uint32_t i = 0; i < n; ++i) { h ^= (uint8_t)p[i]; h *= 1099511628211ULL; }
    return h;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: scs <in.fq> [minOverlap]\n"); return 1; }
    const uint32_t MIN_OV = argc > 2 ? (uint32_t)atoi(argv[2]) : 30;
    auto T0 = std::chrono::steady_clock::now();
    auto lap = [&](const char* w){
        auto t = std::chrono::steady_clock::now();
        fprintf(stderr, "  %-22s %6.2f s\n", w, std::chrono::duration<double>(t-T0).count());
        T0 = t;
    };

    // ── load, drop reads containing N, collapse exact duplicates ─────────────
    std::vector<std::string> reads;
    std::vector<uint32_t>    mult;      // how many input reads each entry stands for
    size_t n_in = 0, n_filtered = 0;
    {
        std::ifstream f(argv[1]);
        std::string l1, l2, l3, l4;
        std::unordered_map<uint64_t, std::vector<uint32_t>> seen;
        seen.reserve(1u << 21);
        while (std::getline(f, l1) && std::getline(f, l2) &&
               std::getline(f, l3) && std::getline(f, l4)) {
            ++n_in;
            if (l2.find('N') != std::string::npos) { ++n_filtered; continue; }
            uint64_t h = fnv(l2.data(), (uint32_t)l2.size());
            auto& b = seen[h];
            bool dup = false;
            for (uint32_t id : b)
                if (reads[id] == l2) { ++mult[id]; dup = true; break; }
            if (!dup) { b.push_back((uint32_t)reads.size());
                        reads.push_back(l2); mult.push_back(1); }
        }
    }
    const uint32_t n = (uint32_t)reads.size();
    uint32_t Lmax = 0;
    for (auto& r : reads) Lmax = std::max(Lmax, (uint32_t)r.size());
    fprintf(stderr, "reads in=%zu  N-filtered=%zu  unique=%u  maxlen=%u\n",
            n_in, n_filtered, n, Lmax);
    lap("load+filter+dedup");

    // ── greedy: sweep overlap length downwards ───────────────────────────────
    std::vector<uint32_t> nxt(n, NONE), prv(n, NONE), ovl(n, 0);
    std::vector<uint32_t> chain_head(n), chain_tail(n);
    for (uint32_t i = 0; i < n; ++i) { chain_head[i] = i; chain_tail[i] = i; }

    std::vector<uint32_t> tails(n), heads(n);
    for (uint32_t i = 0; i < n; ++i) { tails[i] = i; heads[i] = i; }

    std::unordered_map<uint64_t, std::vector<uint32_t>> idx;
    size_t links = 0;
    for (uint32_t L = Lmax - 1; L >= MIN_OV; --L) {
        // compact: keep only reads still open on the relevant side and long enough
        heads.erase(std::remove_if(heads.begin(), heads.end(), [&](uint32_t i){
            return prv[i] != NONE || reads[i].size() < L; }), heads.end());
        tails.erase(std::remove_if(tails.begin(), tails.end(), [&](uint32_t i){
            return nxt[i] != NONE || reads[i].size() < L; }), tails.end());
        if (heads.empty() || tails.empty()) break;

        idx.clear();
        idx.reserve(heads.size() * 2);
        for (uint32_t b : heads) idx[fnv(reads[b].data(), L)].push_back(b);

        for (uint32_t a : tails) {
            if (nxt[a] != NONE) continue;
            const std::string& A = reads[a];
            auto it = idx.find(fnv(A.data() + (A.size() - L), L));
            if (it == idx.end()) continue;
            for (uint32_t b : it->second) {
                if (b == a || prv[b] != NONE) continue;
                if (chain_head[a] == b) continue;                 // would close a cycle
                if (memcmp(A.data() + (A.size() - L), reads[b].data(), L) != 0) continue;
                nxt[a] = b; prv[b] = a; ovl[a] = L;
                uint32_t h = chain_head[a], t = chain_tail[b];
                chain_tail[h] = t; chain_head[t] = h;
                ++links;
                break;
            }
        }
    }
    lap("greedy overlap sweep");

    // ── emit ─────────────────────────────────────────────────────────────────
    std::string pg;
    pg.reserve((size_t)n * Lmax / 3);
    uint32_t n_chains = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (prv[i] != NONE) continue;                              // not a chain head
        ++n_chains;
        uint32_t cur = i;
        pg += reads[cur];
        while (nxt[cur] != NONE) {
            uint32_t o = ovl[cur]; cur = nxt[cur];
            pg.append(reads[cur].data() + o, reads[cur].size() - o);
        }
    }
    lap("emit");
    fprintf(stderr, "links=%zu  chains=%u\n", links, n_chains);
    printf("PG_LEN %zu\n", pg.size());
    if (argc > 3) { std::ofstream o(argv[3], std::ios::binary); o.write(pg.data(), pg.size()); }
    return 0;
}
