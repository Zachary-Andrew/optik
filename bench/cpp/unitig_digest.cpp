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

/// Number of distinct byte values, the size of the complement table.
constexpr int kByteValues = 256;
/// Clearing this bit upper-cases an ASCII letter.
constexpr char kAsciiLowerBit = 0x20;
/// 64-bit FNV-1a offset basis and prime.
constexpr uint64_t kFnvOffsetBasis = 0xcbf29ce484222325ULL;
constexpr uint64_t kFnvPrime = 0x100000001b3ULL;

/// MurmurHash3's fmix64 finalizer, so that digests of similar unitig sets differ in every bit.
/// This tool builds standalone, so it keeps its own copy of mix64() from include/wyhash.hpp.
static uint64_t fmix64(uint64_t value) {
    // Shift and multipliers chosen by Austin Appleby for MurmurHash3 by search for the best avalanche.
    constexpr unsigned shift = 33;
    constexpr uint64_t first_multiplier = 0xff51afd7ed558ccdULL;
    constexpr uint64_t second_multiplier = 0xc4ceb9fe1a85ec53ULL;

    value ^= value >> shift;
    value *= first_multiplier;
    value ^= value >> shift;
    value *= second_multiplier;
    value ^= value >> shift;
    return value;
}

static uint64_t hash_canon(const std::string& seq) {
    struct Comp {
        char table[kByteValues] = {};
        constexpr Comp() {
            constexpr char kFwd[] = "ACGTNacgtn", kRev[] = "TGCANTGCAN";
            for (int i = 0; i < int(sizeof(kFwd)) - 1; ++i) table[static_cast<unsigned char>(kFwd[i])] = kRev[i];
        }
    };
    static constexpr Comp comp_table{};
    const char* comp = comp_table.table;
    const std::size_t len = seq.size();
    // compare the sequence with its reverse complement position by position
    int dir = 0;
    for (std::size_t i = 0; i < len && dir == 0; ++i) {
        const char fwd_base = seq[i] & ~kAsciiLowerBit, rev_base = comp[static_cast<unsigned char>(seq[len - 1 - i])];
        if (fwd_base != rev_base) dir = fwd_base < rev_base ? 1 : -1;
    }
    uint64_t hash = kFnvOffsetBasis;
    for (std::size_t i = 0; i < len; ++i) {
        const char base = dir >= 0 ? (seq[i] & ~kAsciiLowerBit) : comp[static_cast<unsigned char>(seq[len - 1 - i])];
        hash = (hash ^ static_cast<unsigned char>(base)) * kFnvPrime;
    }
    return fmix64(hash);
}

int main(int argc, char** argv) {
    if (argc != 2) { std::fprintf(stderr, "usage: unitig_digest FILE\n"); return 2; }
    std::FILE* file = std::fopen(argv[1], "r");
    if (!file) { std::perror(argv[1]); return 2; }
    uint64_t cnt = 0, bases = 0, sum = 0, xor_sum = 0;
    std::string cur;
    bool have = false;
    auto done = [&] {
        if (!have) return;
        const uint64_t hash = hash_canon(cur);
        ++cnt; bases += cur.size(); sum += hash; xor_sum ^= hash;
        cur.clear(); have = false;
    };
    char* line = nullptr;
    std::size_t cap = 0;
    ssize_t len;
    while ((len = getline(&line, &cap, file)) > 0) {
        while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
        if (line[0] == 'S' && line[1] == '\t') {             // GFA segment
            done();
            const char* seq_start = std::strchr(line + 2, '\t');
            if (!seq_start) continue;
            const char* seq_end = std::strchr(seq_start + 1, '\t');
            cur.assign(seq_start + 1, seq_end ? seq_end : line + len);
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
    std::fclose(file);
    std::printf("%llu %llu %016llx %016llx\n", (unsigned long long)cnt, (unsigned long long)bases,
                (unsigned long long)sum, (unsigned long long)xor_sum);
    return 0;
}
