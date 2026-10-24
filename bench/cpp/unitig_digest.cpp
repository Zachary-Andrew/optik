// unitig_digest - an order- and orientation-independent digest of a unitig set, for comparing graphs too large for the Python validators (a human genome has 35M unitigs and 3.5G bases).
//
// The input is the S lines of a GFA, or a FASTA with one record per unitig at any line width.
// Each unitig is turned to canonical orientation, the smaller of the sequence and its reverse complement, and hashed with 64-bit FNV-1a followed by a finalizer.
// The digest holds the count, the total length, and the sum and xor of the hashes.
// Two partitions of the same k-mer set into the same maximal unitigs have equal digests, while any differing unitig changes the sum and the xor.
//
// Usage: unitig_digest FILE   prints unitigs, bases, sum and xor
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

static uint64_t hash_canon(const std::string& s) {
    struct Comp {
        char t[256] = {};
        constexpr Comp() {
            constexpr char kFwd[] = "ACGTNacgtn", kRev[] = "TGCANTGCAN";
            for (int i = 0; i < 10; ++i) t[static_cast<unsigned char>(kFwd[i])] = kRev[i];
        }
    };
    static constexpr Comp ct{};
    const char* comp = ct.t;
    const std::size_t n = s.size();
    // compare the sequence with its reverse complement position by position
    int dir = 0;
    for (std::size_t i = 0; i < n && dir == 0; ++i) {
        const char f = s[i] & ~0x20, r = comp[static_cast<unsigned char>(s[n - 1 - i])];
        if (f != r) dir = f < r ? 1 : -1;
    }
    uint64_t h = 0xcbf29ce484222325ULL;
    for (std::size_t i = 0; i < n; ++i) {
        const char c = dir >= 0 ? (s[i] & ~0x20) : comp[static_cast<unsigned char>(s[n - 1 - i])];
        h = (h ^ static_cast<unsigned char>(c)) * 0x100000001b3ULL;
    }
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL; h ^= h >> 33;
    return h;
}

int main(int argc, char** argv) {
    if (argc != 2) { std::fprintf(stderr, "usage: unitig_digest FILE\n"); return 2; }
    std::FILE* f = std::fopen(argv[1], "r");
    if (!f) { std::perror(argv[1]); return 2; }
    uint64_t cnt = 0, bases = 0, sum = 0, xr = 0;
    std::string cur;
    bool have = false;
    auto done = [&] {
        if (!have) return;
        const uint64_t h = hash_canon(cur);
        ++cnt; bases += cur.size(); sum += h; xr ^= h;
        cur.clear(); have = false;
    };
    char* line = nullptr;
    std::size_t cap = 0;
    ssize_t len;
    while ((len = getline(&line, &cap, f)) > 0) {
        while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
        if (line[0] == 'S' && line[1] == '\t') {             // GFA segment
            done();
            const char* p = std::strchr(line + 2, '\t');
            if (!p) continue;
            const char* q = std::strchr(p + 1, '\t');
            cur.assign(p + 1, q ? q : line + len);
            have = true;
            done();
        } else if (line[0] == '>') {                          // FASTA record
            done();
            have = true;
        } else if (have && line[0] != 'H' && line[0] != 'L') {
            cur.append(line, static_cast<std::size_t>(len));
        }
    }
    done();
    std::free(line);
    std::fclose(f);
    std::printf("%llu %llu %016llx %016llx\n", (unsigned long long)cnt, (unsigned long long)bases,
                (unsigned long long)sum, (unsigned long long)xr);
    return 0;
}
