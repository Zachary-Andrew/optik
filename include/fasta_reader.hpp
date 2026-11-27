#pragma once

/**
 * @file
 * @brief Streaming FASTA/FASTQ reader, plain or gzip.
 *
 * Records are parsed line by line with getline(3), so memory stays at one line buffer plus the current record.
 * Inputs can be plain text, gzip (detected by extension or magic bytes), or "-" for standard input, which is read as plain text.
 *
 * Gzip is inflated in the calling thread, so each reader costs exactly one thread of the caller's -t budget.
 * The gzip member layout and magic bytes come from RFC 1952.
 * FASTQ follows Cock et al. 2010, where a record's quality section has as many lines as its sequence section.
 * Quality lines can start with '@' (Phred 31) or '>' (Phred 29), which is why they are counted rather than recognised by their first byte.
 */

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <string_view>
#include <stdexcept>
#include <vector>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>

namespace tpoptoa {

inline constexpr std::size_t kGzInputBufferBytes = std::size_t(1) << 20;
/// zlib counts output space in uInt, so a single inflate() call never asks for more than this.
inline constexpr std::size_t kGzMaxInflateChunk = std::size_t(1) << 30;
/// inflateInit2() takes the window size (15 = 32 KiB) plus a flag (16) that selects the gzip wrapper over zlib's.
inline constexpr int kDeflateWindowBits = 15;
inline constexpr int kGzipWrapperFlag = 16;
inline constexpr int kGzipWindowBits = kDeflateWindowBits + kGzipWrapperFlag;
/// Every gzip member starts with these two bytes (RFC 1952, 2.3.1).
inline constexpr unsigned char kGzipId1 = 0x1f, kGzipId2 = 0x8b;
inline constexpr std::size_t kReaderStdioBufferBytes = std::size_t(1) << 20;
inline constexpr int kNoFd = -1;
inline constexpr char kGzSuffix[] = ".gz", kGzSuffixUpper[] = ".GZ";
inline constexpr std::size_t kGzSuffixLen = sizeof(kGzSuffix) - 1;
inline constexpr ssize_t kGetlineEnd = -1;     // getline() at end of file or on error

/// A gzip file exposed as a FILE* through glibc's fopencookie, so the line parser treats compressed and plain input alike.
///
/// We drive inflate directly instead of using gzread(), because zlib 1.3.2's gzread reports a truncated stream as a clean end of file.
/// The same inputs `gzip -dc` rejects are errors here: truncation inside a member, corrupt data, trailing garbage after a member and an empty file.
/// Each of them makes the read return -1, which sets the stream's error flag that is checked after the last line.
/// Concatenated members, like those pigz or `cat a.gz b.gz` produce, are read one after another.
struct GzCookie {
    int file_desc = kNoFd;
    z_stream stream{};
    std::vector<unsigned char> in_buf = std::vector<unsigned char>(kGzInputBufferBytes);
    bool in_member = false;     // a member has started but its end hasn't been seen
    bool any = false;           // at least one member started
    bool done = false;
    bool bad = false;
};

/// fopencookie read hook.
/// It returns the bytes produced, 0 at a clean end and -1 on any error.
inline ssize_t gz_cookie_read(void* cookie, char* buf, size_t max_bytes) {
    auto* gz_state = static_cast<GzCookie*>(cookie);

    if (gz_state->bad) return -1;

    if (gz_state->done || max_bytes == 0) return 0;

    gz_state->stream.next_out = reinterpret_cast<Bytef*>(buf);
    gz_state->stream.avail_out = static_cast<uInt>(max_bytes < kGzMaxInflateChunk ? max_bytes : kGzMaxInflateChunk);
    const uInt want = gz_state->stream.avail_out;

    while (gz_state->stream.avail_out == want) {
        if (gz_state->stream.avail_in == 0) {
            const ssize_t bytes_read = ::read(gz_state->file_desc, gz_state->in_buf.data(), gz_state->in_buf.size());

            if (bytes_read < 0) { gz_state->bad = true; return -1; }

            if (bytes_read == 0) {
                // Running out of input is only clean between members.
                if (gz_state->in_member || !gz_state->any) { gz_state->bad = true; return -1; }

                gz_state->done = true;
                break;                  // return what we have, and 0 on the next call
            }

            gz_state->stream.next_in = gz_state->in_buf.data();
            gz_state->stream.avail_in = static_cast<uInt>(bytes_read);
        }

        if (!gz_state->in_member) {
            // Anything after a member has to be another member, so check its magic bytes.
            if (gz_state->stream.next_in[0] != kGzipId1 || (gz_state->stream.avail_in > 1 && gz_state->stream.next_in[1] != kGzipId2)) {
                gz_state->bad = true;
                return -1;
            }

            if (inflateReset(&gz_state->stream) != Z_OK) { gz_state->bad = true; return -1; }

            gz_state->in_member = gz_state->any = true;
        }

        const int status = inflate(&gz_state->stream, Z_NO_FLUSH);

        if (status == Z_STREAM_END) gz_state->in_member = false;
        else if (status != Z_OK && !(status == Z_BUF_ERROR && gz_state->stream.avail_in == 0)) { gz_state->bad = true; return -1; }
    }

    return static_cast<ssize_t>(want - gz_state->stream.avail_out);
}

/// fopencookie close hook, which returns EOF after a stream error so that fclose() reports it.
inline int gz_cookie_close(void* cookie) {
    auto* gz_state = static_cast<GzCookie*>(cookie);

    inflateEnd(&gz_state->stream);
    ::close(gz_state->file_desc);
    // Reporting inflate failures at close turns a truncated .gz into a read error.
    const bool bad = gz_state->bad;

    delete gz_state;
    return bad ? EOF : 0;
}

inline FILE* gz_fopen(const char* path) {
    const int file_desc = ::open(path, O_RDONLY | O_CLOEXEC);

    if (file_desc < 0) return nullptr;

    auto* gz_state = new GzCookie;

    gz_state->file_desc = file_desc;

    if (inflateInit2(&gz_state->stream, kGzipWindowBits) != Z_OK) { ::close(file_desc); delete gz_state; return nullptr; }

    // fopencookie gives the zlib stream a FILE* interface, letting getline() read .gz input through the same path as plain text.
    cookie_io_functions_t hooks{gz_cookie_read, nullptr, nullptr, gz_cookie_close};
    FILE* file = fopencookie(gz_state, "r", hooks);

    if (!file) { gz_cookie_close(gz_state); return nullptr; }

    setvbuf(file, nullptr, _IOFBF, kReaderStdioBufferBytes);
    return file;
}

struct SeqRecord {
    std::string name;
    std::string seq;
    std::string qual;   // empty for FASTA
};

/// Calls @p callback with every record of @p path in file order.
///
/// @p path may be "-" for standard input (never decompressed), a gzip file, or plain text.
/// Gzip is recognised by a .gz/.GZ name or, for regular files, by the magic bytes.
///
/// A record's name is its header up to the first space or tab, and an empty header becomes `seq_<n>` counting from 0.
/// Trailing CR, space and tab are stripped from each line, while empty lines are skipped entirely.
/// Quality lines end up in SeqRecord::qual, which the k-mer tools ignore but the mapper writes to SAM QUAL.
///
/// @p keep_empty also reports records with no sequence, because the mapper owes one SAM record per read.
/// The index builders, by contrast, skip such records.
///
/// @p max_seq truncates long records to their first @p max_seq bases and reads past the rest.
/// Sampling the start of an input this way avoids holding a whole chromosome in memory.
///
/// Unopenable files, directories and truncated or corrupt gzip all throw std::runtime_error.
/// Lines before the first header are ignored, with a warning on stderr.
template <typename Callback>
void iterate_sequences(const char* path, Callback&& callback, bool keep_empty = false,
                       std::size_t max_seq = SIZE_MAX) {
    FILE* file = nullptr;
    bool gz_in = false;
    std::size_t path_len = std::strlen(path);

    if (path[0] == '-' && path[1] == '\0') {
        file = stdin;
    } else {
        // Checking the magic bytes too means a gzipped file without ".gz" isn't parsed as binary garbage.
        bool is_gzip = path_len >= kGzSuffixLen && (std::strcmp(path + path_len - kGzSuffixLen, kGzSuffix) == 0 ||
                                    std::strcmp(path + path_len - kGzSuffixLen, kGzSuffixUpper) == 0);

        // Probing a pipe would consume its first bytes, so only regular files get checked.
        struct stat pst;

        if (::stat(path, &pst) != 0)
            throw std::runtime_error(std::string("cannot open file: ") + path);

        if (S_ISREG(pst.st_mode)) {
            FILE* probe = std::fopen(path, "rb");

            if (!probe) throw std::runtime_error(std::string("cannot open file: ") + path);

            constexpr std::size_t kGzipIdBytes = 2;
            unsigned char magic[kGzipIdBytes] = {0, 0};

            if (std::fread(magic, 1, kGzipIdBytes, probe) == kGzipIdBytes && magic[0] == kGzipId1 && magic[1] == kGzipId2) is_gzip = true;

            std::fclose(probe);
        }

        if (is_gzip) {
            file = gz_fopen(path);
            gz_in = true;

            if (!file)
                throw std::runtime_error(std::string("cannot open gzipped file: ") + path);
        } else {
            file = std::fopen(path, "r");
        }
    }

    if (!file)
        throw std::runtime_error(std::string("cannot open file: ") + path);

    {
        // On Linux a directory opens fine and then reads as empty, so it has to be rejected explicitly.
        struct stat info;

        if (::fstat(::fileno(file), &info) == 0 && S_ISDIR(info.st_mode)) {
            if (file != stdin) std::fclose(file);

            throw std::runtime_error(std::string("is a directory, not a sequence file: ") + path);
        }
    }

    char*   line_buf  = nullptr;
    size_t  buf_cap   = 0;
    ssize_t line_len  = 0;

    SeqRecord rec;
    // Bare headers (">" or "@" alone) are named seq_<n>, matching how the benchmark readers name them.
    std::size_t record_index = 0;
    // A FASTQ record is @id, sequence lines, '+', then exactly as many quality lines as sequence lines.
    // Those lines are taken as quality whatever their first byte, because Phred+33 scores 31 and 29 print as '@' and '>'.
    bool have_record    = false;
    bool fastq_record   = false;
    std::size_t seq_lines = 0;
    std::size_t qual_left = 0;    // quality lines still expected
    std::size_t orphan_lines = 0; // text before the first header, reported once at the end
    bool capped = false;          // record hit max_seq and was already delivered

    auto flush_record = [&]() {
        // A capped record went out when it reached max_seq.
        if (have_record && !capped && (keep_empty || !rec.seq.empty()))
            callback(rec);

        capped = false;
        rec.name.clear();
        rec.seq.clear();
        rec.qual.clear();
        have_record  = false;
        fastq_record = false;
        seq_lines    = 0;
        qual_left    = 0;
    };

    auto close_input = [&]() -> int {
        if (file == stdin) return 0;

        const int status = std::fclose(file);

        file = nullptr;
        return status;
    };

    try {
    while ((line_len = getline(&line_buf, &buf_cap, file)) != kGetlineEnd) {
        // Stripping CR, spaces and tabs lets CRLF files and padded lines parse alike.
        while (line_len > 0 &&
               (line_buf[line_len - 1] == '\n' || line_buf[line_len - 1] == '\r' ||
                line_buf[line_len - 1] == ' '  || line_buf[line_len - 1] == '\t'))
            line_buf[--line_len] = '\0';

        if (line_len == 0) continue;

        char first = line_buf[0];

        // The quality-block test runs before the header tests because a Phred+33 quality line can begin with @ or +.
        if (qual_left > 0) {      // quality lines are kept verbatim
            --qual_left;
            rec.qual.append(line_buf, static_cast<std::size_t>(line_len));
            continue;
        }

        if (first == '>' || first == '@') {
            flush_record();
            have_record  = true;
            fastq_record = (first == '@');
            std::string_view header(line_buf + 1, static_cast<std::size_t>(line_len - 1));
            auto space = header.find_first_of(" \t");

            rec.name = (space == std::string_view::npos) ? std::string(header)
                                                       : std::string(header.substr(0, space));

            if (rec.name.empty())
                rec.name = "seq_" + std::to_string(record_index);

            ++record_index;
            continue;
        }

        if (first == '+' && have_record) {
            // Multi-line FASTQ has one quality line per sequence line.
            qual_left = seq_lines;
            continue;
        }

        if (have_record) {
            if (fastq_record) ++seq_lines;

            if (capped) continue;

            rec.seq.append(line_buf, static_cast<std::size_t>(line_len));

            // A record longer than max_seq is emitted truncated, and skip_rest discards its remaining lines.
            if (rec.seq.size() >= max_seq) {
                rec.seq.resize(max_seq);
                capped = true;
                callback(rec);
            }
        } else {
            ++orphan_lines;
        }
    }

    // A truncated or corrupt .gz delivers a prefix of the data before it fails.
    // So the read error is checked before the last record goes out, turning a partial input into an error instead of a silently short result.
    if (gz_in) {
        const bool bad = std::ferror(file) != 0;
        const int status = close_input();

        if (bad || status != 0)
            throw std::runtime_error(std::string(
                "gzip decompression failed (truncated or corrupt file): ") + path);
    }

    flush_record();

    // Without this warning, a file that starts mid-record (split, or truncated with its header lost) would drop those bases silently.
    if (orphan_lines)
        std::fprintf(stderr, "[reader] warning: %s: %zu line(s) before the first '>'/'@' header "
                     "were ignored (not FASTA/FASTQ, or the file starts mid-record)\n",
                     path, orphan_lines);
    } catch (...) {
        std::free(line_buf);

        if (file) close_input();

        throw;
    }

    std::free(line_buf);

    if (file) close_input();
}

/// Splits one argument into input paths.
/// The CLI joins repeated -i values (paired-end R1 and R2, say) with '\n', and a lone path comes back unchanged.
inline std::vector<std::string> split_inputs(const char* paths) {
    std::vector<std::string> out;
    std::string cur;

    for (const char* ch = paths; *ch; ++ch) {
        if (*ch == '\n') { out.push_back(cur); cur.clear(); }
        else cur.push_back(*ch);
    }

    out.push_back(cur);
    return out;
}

template <typename Callback>
void iterate_inputs(const char* paths, Callback&& callback, bool keep_empty = false) {
    for (const std::string& path : split_inputs(paths)) iterate_sequences(path.c_str(), callback, keep_empty);
}

}