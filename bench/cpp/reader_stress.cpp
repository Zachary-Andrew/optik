// reader_stress.cpp - corner cases of the FASTA/FASTQ reader (fasta_reader.hpp), namely whitespace and CRLF line ends, blank lines, '@' quality lines, multi-line FASTQ, gzip (by suffix, by magic bytes, multi-member, quoted paths), and truncated gzip, which must throw instead of yielding a prefix of the input.
// base_encode must also read RNA U as T.
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>
#include <unistd.h>
#include "fasta_reader.hpp"
#include "kmer.hpp"
using namespace tpoptoa;

static int g_failures = 0;
#define CHECK(cond, ...)                                                    \
    do { if (!(cond)) { ++g_failures; std::fprintf(stderr, "FAIL: " __VA_ARGS__); \
         std::fprintf(stderr, "\n"); } } while (0)

static std::string g_dir;
static std::string put(const std::string& name, const std::string& data) {
    const std::string p = g_dir + "/" + name;
    std::FILE* f = std::fopen(p.c_str(), "wb");
    std::fwrite(data.data(), 1, data.size(), f);
    std::fclose(f);
    return p;
}
static std::string slurp(const std::string& p) {
    std::string s; char buf[65536]; std::size_t n;
    std::FILE* f = std::fopen(p.c_str(), "rb");
    while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) s.append(buf, n);
    std::fclose(f);
    return s;
}
static std::string gz(const std::string& name, const std::string& data) {
    const std::string p = put(name, data);
    const std::string cmd = "gzip -nf '" + p + "'";
    if (std::system(cmd.c_str()) != 0) return "";
    return p + ".gz";
}
using Recs = std::vector<std::pair<std::string, std::string>>;
static Recs read_all(const std::string& p, bool* threw = nullptr) {
    Recs r;
    try {
        iterate_sequences(p.c_str(), [&](const SeqRecord& x) { r.push_back({x.name, x.seq}); });
        if (threw) *threw = false;
    } catch (const std::exception&) {
        if (threw) *threw = true;
    }
    return r;
}

int main() {
    char tmpl[] = "/tmp/optik_reader_XXXXXX";
    g_dir = mkdtemp(tmpl);
    const Recs want = {{"a", "ACGTACGTAA"}, {"b", "GGGTTTCCCA"}};

    CHECK(read_all(put("ws.fa", ">a desc\t\nACGTA \t\nCGTAA\t\n>b\nGGGTTTCCCA  \n")) == want,
          "trailing spaces/tabs must not be sequence");
    CHECK(read_all(put("crlf.fa", ">a\r\nACGTA\r\nCGTAA\r\n>b\r\nGGGTTTCCCA\r\n")) == want, "CRLF");
    CHECK(read_all(put("blank.fa", ">a\nACGTA\n\nCGTAA\n\n>b\nGGGTTTCCCA")) == want,
          "blank lines / no final newline");
    CHECK(read_all(put("q.fq", "@a\nACGTA\nCGTAA\n+\n@IIII\nIIIII\n@b\nGGGTTTCCCA\n+\n>IIIIIIIII\n")) == want,
          "multi-line FASTQ with '@'/'>' leading quality lines");

    const std::string big = [] {
        std::string s;
        for (int i = 0; i < 20000; ++i) s += ">r" + std::to_string(i) + "\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        return s;
    }();
    const std::string g1 = gz("big.fa", big);
    CHECK(!g1.empty(), "gzip available");
    if (!g1.empty()) {
        bool threw = true;
        CHECK(read_all(g1, &threw).size() == 20000 && !threw, "gzip by suffix");
        const std::string raw = slurp(g1);
        CHECK(read_all(put("noext.fa", raw)).size() == 20000, "gzip detected by magic bytes");
        const std::string mm = put("mm.fa.gz", raw + raw);
        CHECK(read_all(mm).size() == 40000, "multi-member gzip");
        const std::string q = put("it's here.fa.gz", raw);
        CHECK(read_all(q, &threw).size() == 20000 && !threw, "path with a quote and a space");
        for (double frac : {0.1, 0.5, 0.95}) {
            const std::string t = put("trunc.fa.gz", raw.substr(0, std::size_t(raw.size() * frac)));
            read_all(t, &threw);
            CHECK(threw, "truncated gzip (%.0f%%) must throw", frac * 100);
        }
        bool missing_threw = false;
        read_all(g_dir + "/nonexist.fa.gz", &missing_threw);
        CHECK(missing_threw, "missing .gz must throw");
    }

    CHECK(base_encode('U') == base_encode('T') && base_encode('u') == base_encode('t'),
          "U encodes as T");
    CHECK(base_encode('N') == 0xFF && base_encode(' ') == 0xFF, "N / space invalid");

    std::error_code ec;
    std::filesystem::remove_all(g_dir, ec);  // a leftover scratch directory is not a test failure
    std::printf(g_failures ? "reader_stress: %d FAILURE(S)\n" : "reader_stress: ALL PASS\n", g_failures);
    return g_failures != 0;
}
