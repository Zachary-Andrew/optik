#pragma once

/*
 * fasta_reader.hpp - streaming FASTA/FASTQ reader, plain or gzip.
 *
 * Records are parsed line by line with getline(3), so memory is one line buffer plus the current record.
 * Inputs may be plain text, gzip (by extension or by magic bytes), or "-" for standard input (read as plain text).
 *
 * Gzip input is inflated in the calling thread through zlib, so a reader is exactly one thread of the caller's -t budget.
 * The format sources are RFC 1952 for the member layout and magic bytes and zlib for inflate.
 * FASTQ follows Cock et al. 2010, where a record's quality section has as many lines as its sequence section, and a quality line may begin with '@', Phred 31, as it may with '>', Phred 29, so quality lines are counted rather than recognised by their first byte.
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

/// Compressed bytes read from the file per read(2).
inline constexpr std::size_t kGzInputBufferBytes = std::size_t(1) << 20;
/// Largest output request handed to one inflate() call, since zlib counts in uInt.
inline constexpr std::size_t kGzMaxInflateChunk = std::size_t(1) << 30;
/// inflateInit2() window bits, 15 for the 32 KiB deflate window plus 16, which makes zlib expect the gzip wrapper instead of the zlib one.
inline constexpr int kGzipWindowBits = 15 + 16;
/// The two magic bytes that start every gzip member (RFC 1952, 2.3.1).
inline constexpr unsigned char kGzipId1 = 0x1f, kGzipId2 = 0x8b;
/// stdio buffer of the FILE* that wraps the decoder.
inline constexpr std::size_t kReaderStdioBufferBytes = std::size_t(1) << 20;

/// GzCookie - state of a gzip file opened as a FILE* (glibc fopencookie over zlib's inflate), so the line parser reads compressed and plain input alike.
///
/// inflate is driven directly rather than through gzread(), because gzread in zlib 1.3.2 reports a stream cut short as a clean end of file.
/// As for `gzip -dc`, the following are errors, each a read of -1 that sets the stream's error flag checked after the last line.
/// They are input truncated inside a member, corrupt data, trailing bytes after a member that do not start another one, and an empty file.
/// Concatenated members, as pigz or cat a.gz b.gz produce, are read in turn.
struct GzCookie {
    int fd = -1;
    z_stream zs{};
    std::vector<unsigned char> in = std::vector<unsigned char>(kGzInputBufferBytes);
    bool in_member = false;     // inside a member whose end is not yet seen
    bool any = false;           // a member was started
    bool done = false;
    bool bad = false;
};

/// gz_cookie_read() - fopencookie read hook that inflates up to @n bytes into
/// @buf.
/// Returns the bytes produced, 0 at a clean end, -1 on any error.
inline ssize_t gz_cookie_read(void* cookie, char* buf, size_t n) {
    auto* gz = static_cast<GzCookie*>(cookie);
    if (gz->bad) return -1;
    if (gz->done || n == 0) return 0;
    gz->zs.next_out = reinterpret_cast<Bytef*>(buf);
    gz->zs.avail_out = static_cast<uInt>(n < kGzMaxInflateChunk ? n : kGzMaxInflateChunk);
    const uInt want = gz->zs.avail_out;
    while (gz->zs.avail_out == want) {
        if (gz->zs.avail_in == 0) {
            const ssize_t r = ::read(gz->fd, gz->in.data(), gz->in.size());
            if (r < 0) { gz->bad = true; return -1; }
            if (r == 0) {
                // end of the compressed file, which is clean only between members
                if (gz->in_member || !gz->any) { gz->bad = true; return -1; }
                gz->done = true;
                break;                  // what this call decoded, then 0 next time
            }
            gz->zs.next_in = gz->in.data();
            gz->zs.avail_in = static_cast<uInt>(r);
        }
        if (!gz->in_member) {
            // a new member must start with the gzip magic bytes
            if (gz->zs.next_in[0] != kGzipId1 || (gz->zs.avail_in > 1 && gz->zs.next_in[1] != kGzipId2)) {
                gz->bad = true;
                return -1;
            }
            if (inflateReset(&gz->zs) != Z_OK) { gz->bad = true; return -1; }
            gz->in_member = gz->any = true;
        }
        const int status = inflate(&gz->zs, Z_NO_FLUSH);
        if (status == Z_STREAM_END) gz->in_member = false;
        else if (status != Z_OK && !(status == Z_BUF_ERROR && gz->zs.avail_in == 0)) { gz->bad = true; return -1; }
    }
    return static_cast<ssize_t>(want - gz->zs.avail_out);
}

/// gz_cookie_close() - fopencookie close hook, returning EOF if the stream had an error so that fclose() reports it.
inline int gz_cookie_close(void* cookie) {
    auto* gz = static_cast<GzCookie*>(cookie);
    inflateEnd(&gz->zs);
    ::close(gz->fd);
    const bool bad = gz->bad;
    delete gz;
    return bad ? EOF : 0;
}

/// gz_fopen() - open a gzip file for reading as a FILE*, or nullptr.
inline FILE* gz_fopen(const char* path) {
    const int fd = ::open(path, O_RDONLY | O_CLOEXEC);
    if (fd < 0) return nullptr;
    auto* gz = new GzCookie;
    gz->fd = fd;
    if (inflateInit2(&gz->zs, kGzipWindowBits) != Z_OK) { ::close(fd); delete gz; return nullptr; }
    cookie_io_functions_t io{gz_cookie_read, nullptr, nullptr, gz_cookie_close};
    FILE* f = fopencookie(gz, "r", io);
    if (!f) { gz_cookie_close(gz); return nullptr; }
    setvbuf(f, nullptr, _IOFBF, kReaderStdioBufferBytes);
    return f;
}

struct SeqRecord {
    std::string name;
    std::string seq;
    std::string qual;   // FASTQ quality string, empty for FASTA
};

/// iterate_sequences() - call @cb with every record of the FASTA/FASTQ file
/// @path, in file order.
///
/// @path is "-" for standard input, which is not decompressed, or a gzip file, whose name ends in .gz or .GZ or which is a regular file starting with the gzip magic bytes), or a plain text file.
///
/// A record's name is its header up to the first space or tab.
/// An empty header is named "seq_<n>", with n counting records from 0.
/// Trailing CR, space and tab are stripped from every line, and empty lines are skipped.
/// FASTQ quality lines are kept in SeqRecord::qual, which the k-mer tools ignore and the mapper writes to SAM QUAL.
///
/// @keep_empty: also report records with an empty sequence, since the mapper must write one SAM record per input read, while the index builders skip them.
///
/// @max_seq: a record reaching @max_seq bases is delivered with its first
/// @max_seq bases, and the rest of it is read past without being kept (a sample of the input's start never holds a whole chromosome).
///
/// Throws std::runtime_error if the file cannot be opened, is a directory, or is a truncated or corrupt gzip file.
/// Lines before the first header are ignored with a warning on stderr.
template <typename Callback>
void iterate_sequences(const char* path, Callback&& cb, bool keep_empty = false,
                       std::size_t max_seq = SIZE_MAX) {
    FILE* fh = nullptr;
    bool gz_in = false;
    std::size_t path_len = std::strlen(path);
    if (path[0] == '-' && path[1] == '\0') {
        fh = stdin;
    } else {
        // gzip by extension or by magic bytes, so a gzipped file without ".gz" is not parsed as binary garbage
        bool gz = path_len >= 3 && (std::strcmp(path + path_len - 3, ".gz") == 0 ||
                                    std::strcmp(path + path_len - 3, ".GZ") == 0);
        // regular files only, since probing a pipe would consume its first bytes
        struct stat pst;
        if (::stat(path, &pst) != 0)
            throw std::runtime_error(std::string("cannot open file: ") + path);
        if (S_ISREG(pst.st_mode)) {
            FILE* probe = std::fopen(path, "rb");
            if (!probe) throw std::runtime_error(std::string("cannot open file: ") + path);
            unsigned char magic[2] = {0, 0};
            if (std::fread(magic, 1, 2, probe) == 2 && magic[0] == kGzipId1 && magic[1] == kGzipId2) gz = true;
            std::fclose(probe);
        }
        if (gz) {
            fh = gz_fopen(path);
            gz_in = true;
            if (!fh)
                throw std::runtime_error(std::string("cannot open gzipped file: ") + path);
        } else {
            fh = std::fopen(path, "r");
        }
    }
    if (!fh)
        throw std::runtime_error(std::string("cannot open file: ") + path);
    {
        // a directory opens fine on Linux and then reads as empty, so reject it
        struct stat st;
        if (::fstat(::fileno(fh), &st) == 0 && S_ISDIR(st.st_mode)) {
            if (fh != stdin) std::fclose(fh);
            throw std::runtime_error(std::string("is a directory, not a sequence file: ") + path);
        }
    }

    char*   line_buf  = nullptr;
    size_t  buf_cap   = 0;
    ssize_t line_len  = 0;

    SeqRecord rec;
    // Record count for the seq_<n> names of bare headers, ">" or "@" with no name, which the benchmark readers name the same way.
    std::size_t record_index = 0;
    // A FASTQ record is @id / sequence lines / + / quality lines, with as many quality lines as sequence lines.
    // After the '+' exactly that many lines are taken as quality whatever their first byte, since Phred+33 scores 31 and 29 are the characters '@' and '>', which would otherwise be read as record headers.
    bool have_record    = false;
    bool fastq_record   = false;
    std::size_t seq_lines = 0;    // sequence lines in the current section
    std::size_t qual_left = 0;    // quality lines still to skip
    std::size_t orphan_lines = 0; // lines before the first header (ignored, reported)
    bool capped = false;          // the record reached max_seq and was delivered

    auto flush_record = [&]() {
        if (have_record && !capped && (keep_empty || !rec.seq.empty()))
            cb(rec);
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
        if (fh == stdin) return 0;
        const int st = std::fclose(fh);
        fh = nullptr;
        return st;
    };
    try {
    while ((line_len = getline(&line_buf, &buf_cap, fh)) != -1) {
        while (line_len > 0 &&
               (line_buf[line_len - 1] == '\n' || line_buf[line_len - 1] == '\r' ||
                line_buf[line_len - 1] == ' '  || line_buf[line_len - 1] == '\t'))
            line_buf[--line_len] = '\0';

        if (line_len == 0) continue;

        char first = line_buf[0];

        if (qual_left > 0) {      // a quality line, kept verbatim
            --qual_left;
            rec.qual.append(line_buf, static_cast<std::size_t>(line_len));
            continue;
        }

        if (first == '>' || first == '@') {
            flush_record();
            have_record  = true;
            fastq_record = (first == '@');
            std::string_view sv(line_buf + 1, static_cast<std::size_t>(line_len - 1));
            auto space = sv.find_first_of(" \t");
            rec.name = (space == std::string_view::npos) ? std::string(sv)
                                                       : std::string(sv.substr(0, space));
            if (rec.name.empty())
                rec.name = "seq_" + std::to_string(record_index);
            ++record_index;
            continue;
        }

        if (first == '+' && have_record) {
            qual_left = seq_lines;  // mirror the sequence line count
            continue;
        }

        if (have_record) {
            if (fastq_record) ++seq_lines;
            if (capped) continue;
            rec.seq.append(line_buf, static_cast<std::size_t>(line_len));
            if (rec.seq.size() >= max_seq) {
                rec.seq.resize(max_seq);
                capped = true;
                cb(rec);
            }
        } else {
            ++orphan_lines;       // text before the first '>' / '@' header
        }
    }

    // A truncated or corrupt .gz delivers a prefix of the data and then a A read error is checked before delivering the last record, so a partial input is an error rather than a silently short result.
    if (gz_in) {
        const bool bad = std::ferror(fh) != 0;
        const int st = close_input();
        if (bad || st != 0)
            throw std::runtime_error(std::string(
                "gzip decompression failed (truncated or corrupt file): ") + path);
    }
    flush_record();
    // a file that starts mid-record (a split or truncated file, a header lost) would otherwise lose those bases without a word
    if (orphan_lines)
        std::fprintf(stderr, "[reader] warning: %s: %zu line(s) before the first '>'/'@' header "
                     "were ignored (not FASTA/FASTQ, or the file starts mid-record)\n",
                     path, orphan_lines);
    } catch (...) {
        std::free(line_buf);
        if (fh) close_input();
        throw;
    }

    std::free(line_buf);
    if (fh) close_input();
}

/// split_inputs() - the input files of one argument, as paths joined by '\n' (the CLI's repeated -i, e.g. paired-end R1 and R2).
/// A single path is itself.
inline std::vector<std::string> split_inputs(const char* paths) {
    std::vector<std::string> out;
    std::string cur;
    for (const char* c = paths; *c; ++c) {
        if (*c == '\n') { out.push_back(cur); cur.clear(); }
        else cur.push_back(*c);
    }
    out.push_back(cur);
    return out;
}

/// iterate_inputs() - iterate_sequences() over every input of split_inputs(@paths), in order.
template <typename Callback>
void iterate_inputs(const char* paths, Callback&& cb, bool keep_empty = false) {
    for (const std::string& p : split_inputs(paths)) iterate_sequences(p.c_str(), cb, keep_empty);
}

}