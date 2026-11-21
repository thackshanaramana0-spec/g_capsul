// fmi_selftest -- checks include/caps_fmi.h against brute force.
//
// 1. SA-IS against a naive suffix sort, on random and adversarial texts
//    (single-symbol runs, periodic strings, tiny alphabets, length 1 and 2).
// 2. FM-index read retrieval against a direct scan of every read, for
//    random reads, repetitive reads, reads with N, reads longer than the
//    locate sample interval, and patterns that are substrings, random
//    strings, strings with N and strings with symbols that are never indexed.
// 3. Save / load round trip gives identical answers.
//
// Exit status is 0 only if every comparison agrees.
#include "caps_fmi.h"
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++failures; fprintf(stderr, "FAIL: " __VA_ARGS__); fprintf(stderr, "\n"); } } while (0)

static void test_sais(std::mt19937_64& rng) {
    size_t cases = 0;
    auto run = [&](std::vector<uint8_t> s) {
        // s holds symbols >= 1, terminator 0 appended
        s.push_back(0);
        const uint32_t n = (uint32_t)s.size();
        uint32_t K = 0; for (auto c : s) K = std::max<uint32_t>(K, c);
        std::vector<uint32_t> sa(n);
        capsfmi::detail::sais<uint8_t, uint32_t>(s.data(), sa.data(), n, K + 1);
        std::vector<uint32_t> ref(n);
        for (uint32_t i = 0; i < n; ++i) ref[i] = i;
        std::sort(ref.begin(), ref.end(), [&](uint32_t a, uint32_t b) {
            return std::lexicographical_compare(s.begin() + a, s.end(), s.begin() + b, s.end());
        });
        CHECK(sa == ref, "SA-IS differs from naive sort (n=%u)", n);
        ++cases;
    };
    run({});
    run({1});
    run({1, 1});
    for (int len : {1, 2, 3, 7, 64, 500}) {
        run(std::vector<uint8_t>(len, 3));                       // single-symbol run
        std::vector<uint8_t> per; for (int i = 0; i < len; ++i) per.push_back((uint8_t)(1 + i % 3));
        run(per);                                                // periodic
    }
    for (int it = 0; it < 400; ++it) {
        const int len = 1 + (int)(rng() % 600);
        const int K = 1 + (int)(rng() % 6);
        std::vector<uint8_t> s(len);
        for (auto& c : s) c = (uint8_t)(1 + rng() % K);
        run(s);
    }
    fprintf(stderr, "[selftest] SA-IS: %zu texts compared with a naive sort\n", cases);
}

static std::string rc(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (auto& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : c;
    return r;
}

static void check_reads(const std::vector<std::string>& reads, uint32_t D, std::mt19937_64& rng, size_t& queries) {
    std::vector<uint32_t> ids(reads.size());
    for (size_t i = 0; i < ids.size(); ++i) ids[i] = (uint32_t)(1000 + 7 * i);   // non-trivial identities
    capsfmi::FMIndex fm, fm2;
    std::string err;
    // Every third read is registered as a verbatim read: its identity must come
    // back flagged, and its stored sequence must survive save / load unchanged.
    std::vector<std::string> vseq;
    for (size_t i = 0; i < reads.size(); ++i)
        if (i % 3 == 0 && !reads[i].empty()) {
            ids[i] = fm.add_verbatim(reads[i], (uint32_t)(i + 5));
            vseq.push_back(reads[i]);
        }
    CHECK(fm.build(reads, ids, D, err), "build failed: %s", err.c_str());
    const std::string path = "fmi_selftest.tmp";
    CHECK(fm.save(path, err), "save failed: %s", err.c_str());
    CHECK(fm2.load(path, err), "load failed: %s", err.c_str());
    std::remove(path.c_str());
    CHECK(fm2.num_verbatim() == vseq.size(), "verbatim count %zu != %zu", fm2.num_verbatim(), vseq.size());
    for (size_t i = 0, k = 0; i < reads.size(); ++i)
        if (i % 3 == 0 && !reads[i].empty()) {
            CHECK((ids[i] & capsfmi::FMIndex::VERBATIM_FLAG) != 0, "verbatim flag missing");
            CHECK(fm2.verbatim_seq(ids[i]) == vseq[k], "verbatim sequence changed by save/load");
            CHECK(fm2.verbatim_place(ids[i]) == (uint32_t)(i + 5), "verbatim placement changed by save/load");
            ++k;
        }

    std::vector<std::string> pats;
    for (int k = 0; k < 120; ++k) {
        if (!reads.empty() && k % 3 != 2) {                      // substring of a read
            const std::string& r = reads[rng() % reads.size()];
            if (r.empty()) continue;
            const size_t a = rng() % r.size(), len = 1 + rng() % std::min<size_t>(r.size() - a, 60);
            pats.push_back(r.substr(a, len));
        } else {                                                  // random string
            const char* al = (k % 6 == 5) ? "ACGTN" : "ACGT";
            std::string p(1 + rng() % 12, 'A');
            for (auto& c : p) c = al[rng() % strlen(al)];
            pats.push_back(p);
        }
    }
    pats.push_back("ACGX");                                       // symbol never indexed
    pats.push_back("N");
    for (const std::string& p : pats) {
        for (const std::string& q : {p, rc(p)}) {
            std::set<uint32_t> truth;
            for (size_t i = 0; i < reads.size(); ++i) if (reads[i].find(q) != std::string::npos) truth.insert(ids[i]);
            for (const capsfmi::FMIndex* f : {&fm, &fm2}) {
                std::vector<uint32_t> got;
                CHECK(f->find(q, got, err), "find failed: %s", err.c_str());
                std::set<uint32_t> gs(got.begin(), got.end());
                CHECK(gs == truth, "FM-index answer differs from scan for pattern %s (%zu vs %zu reads)",
                      q.c_str(), gs.size(), truth.size());
            }
            ++queries;
        }
    }
}

static void test_fmi(std::mt19937_64& rng) {
    size_t queries = 0, sets = 0;
    auto rand_read = [&](size_t len, const char* al) {
        std::string s(len, 'A');
        for (auto& c : s) c = al[rng() % strlen(al)];
        return s;
    };
    check_reads({}, 256, rng, queries); ++sets;                   // empty index
    check_reads({"A"}, 256, rng, queries); ++sets;
    check_reads({"", "ACGT", ""}, 256, rng, queries); ++sets;     // empty reads are skipped
    for (int it = 0; it < 60; ++it) {
        std::vector<std::string> reads;
        const int m = 1 + (int)(rng() % 80);
        const int kind = it % 5;
        for (int i = 0; i < m; ++i) {
            if (kind == 0) reads.push_back(rand_read(1 + rng() % 160, "ACGT"));
            else if (kind == 1) reads.push_back(rand_read(1 + rng() % 160, "ACGTN"));
            else if (kind == 2) reads.push_back(rand_read(1 + rng() % 40, "AC"));          // highly repetitive
            else if (kind == 3) reads.push_back(std::string(1 + rng() % 50, "ACGT"[rng() % 4]));
            else reads.push_back(rand_read(300 + rng() % 1500, "ACGT"));                   // longer than D
        }
        // shared sequence across reads, as overlapping reads have
        if (m > 3) { const std::string core = rand_read(30, "ACGT"); for (int i = 0; i < 3; ++i) reads[i] += core; }
        check_reads(reads, (kind == 4) ? 64 : 256, rng, queries);
        ++sets;
    }
    fprintf(stderr, "[selftest] FM-index: %zu read sets, %zu pattern queries compared with a scan\n", sets, queries);
}

int main() {
    std::mt19937_64 rng(20261120);
    test_sais(rng);
    test_fmi(rng);
    if (failures) { fprintf(stderr, "[selftest] %d FAILURES\n", failures); return 1; }
    fprintf(stderr, "[selftest] all checks passed\n");
    return 0;
}
