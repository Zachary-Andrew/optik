// reader_stress.cpp - corner cases of the FASTA/FASTQ reader (fasta_reader.hpp).
// It feeds whitespace and CRLF line ends, blank lines, '@' quality lines and multi-line FASTQ.
// Gzip input is opened by suffix, by magic bytes, as several members and through quoted paths.
// Truncated gzip must throw instead of yielding a prefix of the input.
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
    const std::string path = g_dir + "/" + name;
    std::FILE* file = std::fopen(path.c_str(), "wb");
    std::fwrite(data.data(), 1, data.size(), file);
    std::fclose(file);
    return path;
}
static std::string slurp(const std::string& path) {
    std::string data; char buf[65536]; std::size_t got;
    std::FILE* file = std::fopen(path.c_str(), "rb");
    while ((got = std::fread(buf, 1, sizeof buf, file)) > 0) data.append(buf, got);
    std::fclose(file);
    return data;
}
static std::string write_gz(const std::string& name, const std::string& data) {
    const std::string path = put(name, data);
    const std::string cmd = "gzip -nf '" + path + "'";
    if (std::system(cmd.c_str()) != 0) return "";
    return path + ".gz";
}
using Recs = std::vector<std::pair<std::string, std::string>>;
static Recs read_all(const std::string& path, bool* threw = nullptr) {
    Recs recs;
    try {
        iterate_sequences(path.c_str(), [&](const SeqRecord& rec) { recs.push_back({rec.name, rec.seq}); });
        if (threw) *threw = false;
    } catch (const std::exception&) {
        if (threw) *threw = true;
    }
    return recs;
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
        std::string text;
        for (int i = 0; i < 20000; ++i) text += ">r" + std::to_string(i) + "\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        return text;
    }();
    const std::string gz_path = write_gz("big.fa", big);
    CHECK(!gz_path.empty(), "gzip available");
    if (!gz_path.empty()) {
        bool threw = true;
        CHECK(read_all(gz_path, &threw).size() == 20000 && !threw, "gzip by suffix");
        const std::string raw = slurp(gz_path);
        CHECK(read_all(put("noext.fa", raw)).size() == 20000, "gzip detected by magic bytes");
        const std::string multi_member = put("mm.fa.gz", raw + raw);
        CHECK(read_all(multi_member).size() == 40000, "multi-member gzip");
        const std::string quoted = put("it's here.fa.gz", raw);
        CHECK(read_all(quoted, &threw).size() == 20000 && !threw, "path with a quote and a space");
        for (double frac : {0.1, 0.5, 0.95}) {
            const std::string truncated = put("trunc.fa.gz", raw.substr(0, std::size_t(raw.size() * frac)));
            read_all(truncated, &threw);
            CHECK(threw, "truncated gzip (%.0f%%) must throw", frac * 100);
        }
        bool missing_threw = false;
        read_all(g_dir + "/nonexist.fa.gz", &missing_threw);
        CHECK(missing_threw, "missing .gz must throw");
    }

    CHECK(base_encode('U') == base_encode('T') && base_encode('u') == base_encode('t'),
          "U encodes as T");
    CHECK(base_encode('N') == 0xFF && base_encode(' ') == 0xFF, "N / space invalid");

    std::error_code err_code;
    std::filesystem::remove_all(g_dir, err_code);  // a leftover scratch directory is not a test failure
    std::printf(g_failures ? "reader_stress: %d FAILURE(S)\n" : "reader_stress: ALL PASS\n", g_failures);
    return g_failures != 0;
}
