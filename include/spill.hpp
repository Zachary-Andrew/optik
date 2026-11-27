#pragma once

/**
 * @file
 * @brief Disk spill for the k-mer counter's partition buffers.
 *
 * KMC 2 (Deorowicz et al. 2015) writes super-k-mers to disk bins by minimizer and counts one bin at a time.
 * We keep the bins in memory while they fit and write them out only once the memory budget is exceeded.
 * Afterwards the bins are read back one at a time, as in KMC.
 *
 * @par Memory budget
 * It is the user's --max-mem, or else half of what the process may still take.
 * "What it may take" is the smaller of MemAvailable and the tightest cgroup v2 headroom on the path to the root.
 * Checking cgroups is what makes systemd-run -p MemoryMax=... scopes and containers work.
 *
 * @par Spill files
 * Each worker owns one file, opened with O_TMPFILE (or created and immediately unlinked), so a crash leaves nothing behind.
 * Workers append whole partition buffers with pwritev and index the written extents per bin.
 * Reading a bin back then takes only a few large preads.
 * Short writes, ENOSPC and EIO all throw.
 * After a bin's last read, its extents are punched out (FALLOC_FL_PUNCH_HOLE), so disk use shrinks as counting proceeds.
 */

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/statfs.h>
#include <sys/statvfs.h>
#include <sys/uio.h>
#include <unistd.h>

namespace tpoptoa {

/// Process-wide counter limits, set by the CLI before a build (ok_set_count_limits() in src/ffi.cpp).
struct CountLimits {
    uint64_t max_mem = 0;      // bytes, 0 = half of what's available
    std::string tmp_dir;       // --tmp-dir, else the output's directory, empty = "."
    int graph_mode = 0;        // --graph: 0 auto, 1 partitioned, 2 table
};
inline CountLimits& count_limits() { static CountLimits limits; return limits; }

namespace spill {

/// Used when neither /proc/meminfo nor any cgroup reports a limit.
inline constexpr uint64_t kFallbackAvailableBytes = uint64_t(8) << 30;

/// statfs f_type values of RAM-backed filesystems (linux/magic.h).
inline constexpr long kTmpfsMagicNumber = 0x01021994;
inline constexpr long kRamfsMagicNumber = 0x858458f6;
/// The cgroup v2 line of /proc/self/cgroup reads "0::<path>".
inline constexpr char kCgroupV2Prefix[] = "0::";
inline constexpr std::size_t kCgroupV2PrefixLen = sizeof(kCgroupV2Prefix) - 1;
inline constexpr mode_t kSpillFileMode = 0600;   // owner read/write only
inline constexpr int kNoFd = -1;
/// /proc/meminfo labels its values "kB", but they are binary kilobytes.
inline constexpr uint64_t kBytesPerMeminfoKb = 1024;

/// Reads the first token of @p path as an integer.
/// A missing or unparsable file yields UINT64_MAX, and so does "max", which is how cgroups spell "no limit".
inline uint64_t read_u64_file(const std::string& path) {
    std::ifstream file(path);
    std::string token;

    if (!(file >> token) || token == "max") return UINT64_MAX;

    try { return std::stoull(token); } catch (...) { return UINT64_MAX; }
}

/// MemAvailable, clamped by every cgroup v2 limit between this process's cgroup and the root.
inline uint64_t available_bytes() {
    uint64_t avail = UINT64_MAX;

    {
        std::ifstream meminfo("/proc/meminfo");
        std::string key; uint64_t kib; std::string unit;

        while (meminfo >> key >> kib) {
            std::getline(meminfo, unit);

            if (key == "MemAvailable:") { avail = kib * kBytesPerMeminfoKb; break; }
        }
    }
    // The cgroup v2 line of /proc/self/cgroup names this process's cgroup.
    std::ifstream cgroup_file("/proc/self/cgroup");
    std::string line, cgroup_path;

    while (std::getline(cgroup_file, line))
        if (line.rfind(kCgroupV2Prefix, 0) == 0) { cgroup_path = line.substr(kCgroupV2PrefixLen); break; }

    if (!cgroup_path.empty()) {
        // Any ancestor's memory.max can bind, so the walk climbs to the root and keeps the smallest room left.
        for (std::string path = cgroup_path;;) {
            const std::string dir = "/sys/fs/cgroup" + (path == "/" ? std::string() : path);
            const uint64_t limit = read_u64_file(dir + "/memory.max");

            if (limit != UINT64_MAX) {
                const uint64_t used = read_u64_file(dir + "/memory.current");
                const uint64_t room = (used != UINT64_MAX && used < limit) ? limit - used : 0;

                if (room < avail) avail = room;
            }

            if (path == "/" || path.empty()) break;

            const auto slash = path.find_last_of('/');

            path = slash == 0 || slash == std::string::npos ? "/" : path.substr(0, slash);
        }
    }

    return avail == UINT64_MAX ? kFallbackAvailableBytes : avail;
}

inline uint64_t resident_bytes() {
    long pages = 0, resident = 0;

    if (std::FILE* statm = std::fopen("/proc/self/statm", "r")) {
        constexpr int kStatmFields = 2;   // total and resident pages

        if (std::fscanf(statm, "%ld %ld", &pages, &resident) != kStatmFields) resident = 0;

        std::fclose(statm);
    }

    return static_cast<uint64_t>(resident) * static_cast<uint64_t>(sysconf(_SC_PAGESIZE));
}

/// Returns --max-mem when set, or half of available_bytes() otherwise.
/// @p explicit_limit reports which of the two applied.
inline uint64_t memory_budget(bool& explicit_limit) {
    const uint64_t max_mem = count_limits().max_mem;

    explicit_limit = max_mem != 0;
    constexpr uint64_t kAutoBudgetShareDen = 2;   // half of what is available

    return max_mem ? max_mem : available_bytes() / kAutoBudgetShareDen;
}

/// Space available to this user on @p dir's filesystem, or 0 when unknown.
inline uint64_t free_bytes(const std::string& dir) {
    struct statvfs fs_stats;

    if (statvfs(dir.c_str(), &fs_stats) != 0) return 0;

    return static_cast<uint64_t>(fs_stats.f_bavail) * static_cast<uint64_t>(fs_stats.f_frsize);
}

inline std::string spill_dir() {
    const std::string& dir = count_limits().tmp_dir;

    return dir.empty() ? "." : dir;
}

/// Spilling onto tmpfs or ramfs saves no memory, so callers check this first.
inline bool ram_backed(const std::string& dir) {
    struct statfs fs_info;

    if (statfs(dir.c_str(), &fs_info) != 0) return false;

    return fs_info.f_type == kTmpfsMagicNumber || fs_info.f_type == static_cast<decltype(fs_info.f_type)>(kRamfsMagicNumber);
}

/// One worker's anonymous spill file, supporting appends and positional reads.
class File {
public:
    File() = default;
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    ~File() { if (file_desc_ >= 0) ::close(file_desc_); }

    void open(const std::string& dir) {
        if (file_desc_ >= 0) return;
#ifdef O_TMPFILE
        file_desc_ = ::open(dir.c_str(), O_TMPFILE | O_RDWR | O_CLOEXEC, kSpillFileMode);

#endif
        if (file_desc_ < 0) {   // no O_TMPFILE support, so create and unlink by hand
            std::string tmpl = dir + "/optik_spill_XXXXXX";
            std::vector<char> path_buf(tmpl.begin(), tmpl.end()); path_buf.push_back('\0');

            file_desc_ = ::mkstemp(path_buf.data());

            if (file_desc_ >= 0) ::unlink(path_buf.data());
        }

        if (file_desc_ < 0)
            throw std::runtime_error("cannot create a spill file in '" + dir + "': " +
                                     std::strerror(errno) + " (set --tmp-dir)");
    }
    bool is_open() const noexcept { return file_desc_ >= 0; }
    uint64_t size() const noexcept { return end_; }

    /// Appends every buffer in @p iov, retrying short writes and EINTR, and consumes @p iov in the process.
    void append(std::vector<struct iovec>& iov) {
        std::size_t next_iov = 0;

        while (next_iov < iov.size()) {
            // pwritev accepts at most IOV_MAX buffers per call.
            const int cnt = static_cast<int>(std::min<std::size_t>(iov.size() - next_iov, IOV_MAX));
            const ssize_t written = ::pwritev(file_desc_, iov.data() + next_iov, cnt, static_cast<off_t>(end_));

            if (written < 0) {
                if (errno == EINTR) continue;

                throw std::runtime_error(std::string("spill write failed: ") + std::strerror(errno) +
                                         (errno == ENOSPC ? " (the spill directory's disk is full: set --tmp-dir "
                                                            "to a larger disk, or raise --max-mem so less is spilled)"
                                                          : ""));
            }

            // A write of zero bytes would loop forever.
            if (written == 0) throw std::runtime_error("spill write failed: no progress");

            end_ += static_cast<uint64_t>(written);
            std::size_t left = static_cast<std::size_t>(written);

            // pwritev may stop partway through the list. Fully written iovecs are skipped here and the partial one is advanced below.
            while (next_iov < iov.size() && left >= iov[next_iov].iov_len) { left -= iov[next_iov].iov_len; ++next_iov; }

            if (left) {
                iov[next_iov].iov_base = static_cast<char*>(iov[next_iov].iov_base) + left;
                iov[next_iov].iov_len -= left;
            }
        }
    }
    void read(uint64_t off, uint64_t len, uint8_t* dst) const {
        while (len) {
            const ssize_t got = ::pread(file_desc_, dst, len, static_cast<off_t>(off));

            if (got < 0) {
                if (errno == EINTR) continue;

                throw std::runtime_error(std::string("spill read failed: ") + std::strerror(errno));
            }

            // pread returning 0 before len bytes means the spill file shrank after it was written.
            if (got == 0) throw std::runtime_error("spill read failed: file shorter than written");

            dst += got; off += static_cast<uint64_t>(got); len -= static_cast<uint64_t>(got);
        }
    }
    /// Releases blocks that won't be read again, leaving the file size unchanged.
    /// Without hole-punching support the call fails quietly, and the blocks are freed at close instead.
    void punch(uint64_t off, uint64_t len) const {
        if (file_desc_ >= 0 && len)
            (void)::fallocate(file_desc_, FALLOC_FL_PUNCH_HOLE | FALLOC_FL_KEEP_SIZE,
                              static_cast<off_t>(off), static_cast<off_t>(len));
    }
    void will_need(uint64_t off, uint64_t len) const {
        if (file_desc_ >= 0) ::posix_fadvise(file_desc_, static_cast<off_t>(off), static_cast<off_t>(len), POSIX_FADV_WILLNEED);
    }
private:
    int file_desc_ = kNoFd;
    uint64_t end_ = 0;
};

/// Records that many threads distribute over buckets, for trivially copyable Rec.
///
/// Each thread buffers in memory up to its byte budget, and beyond that it appends to its own spill file.
/// load(b) then gathers bucket b from every thread's memory and file.
template <typename Rec>
class BucketSpill {
public:
    BucketSpill(unsigned nthreads, std::size_t nbuckets, uint64_t per_thread_bytes, std::string dir)
        : nbuckets_(nbuckets), budget_(per_thread_bytes), dir_(std::move(dir)), threads_(nthreads) {
        for (auto& thread_bufs : threads_) { thread_bufs.buf.resize(nbuckets_); thread_bufs.extents.resize(nbuckets_); }
    }
    /// May flush the thread's buffers to disk when they outgrow the budget.
    void add(unsigned tid, std::size_t bucket, const Rec& rec) {
        ThreadBuffers& thread_bufs = threads_[tid];
        auto& bucket_buf = thread_bufs.buf[bucket];
        const std::size_t cap0 = bucket_buf.capacity();

        bucket_buf.push_back(rec);

        if (bucket_buf.capacity() != cap0) thread_bufs.held += (bucket_buf.capacity() - cap0) * sizeof(Rec);

        if (thread_bufs.held > budget_) flush(thread_bufs);
    }
    bool spilled() const noexcept {
        for (const ThreadBuffers& thread_bufs : threads_) if (thread_bufs.file.is_open()) return true;

        return false;
    }
    /// Moves every record of @p bucket into @p out and releases the bucket's memory.
    void load(std::size_t bucket, std::vector<Rec>& out) {
        out.clear();
        std::size_t total = 0;

        for (ThreadBuffers& thread_bufs : threads_) {
            total += thread_bufs.buf[bucket].size();

            for (const auto& extent : thread_bufs.extents[bucket]) total += extent.second / sizeof(Rec);
        }

        out.resize(total);
        std::size_t filled = 0;

        for (ThreadBuffers& thread_bufs : threads_) {
            for (const auto& extent : thread_bufs.extents[bucket]) {
                thread_bufs.file.read(extent.first, extent.second, reinterpret_cast<uint8_t*>(out.data() + filled));
                filled += extent.second / sizeof(Rec);
            }

            std::copy(thread_bufs.buf[bucket].begin(), thread_bufs.buf[bucket].end(), out.begin() + static_cast<std::ptrdiff_t>(filled));
            filled += thread_bufs.buf[bucket].size();
            std::vector<Rec>().swap(thread_bufs.buf[bucket]);
        }
    }
    /// Counts records both in memory and on disk.
    std::size_t size(std::size_t bucket) const {
        std::size_t total = 0;

        for (const ThreadBuffers& thread_bufs : threads_) {
            total += thread_bufs.buf[bucket].size();

            for (const auto& extent : thread_bufs.extents[bucket]) total += extent.second / sizeof(Rec);
        }

        return total;
    }
    /// Asks the kernel to read the bucket's on-disk parts ahead of time.
    void will_need(std::size_t bucket) const {
        for (const ThreadBuffers& thread_bufs : threads_) for (const auto& extent : thread_bufs.extents[bucket]) thread_bufs.file.will_need(extent.first, extent.second);
    }

private:
    struct ThreadBuffers {
        std::vector<std::vector<Rec>> buf;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> extents;   // (offset, bytes) per bucket
        uint64_t held = 0;                                                // capacity in bytes, not size
        File file;
    };
    /// Writes all of a thread's buffers to its spill file and frees them.
    void flush(ThreadBuffers& thread_bufs) {
        thread_bufs.file.open(dir_);
        std::vector<struct iovec> iov;
        uint64_t off = thread_bufs.file.size();

        for (std::size_t bucket = 0; bucket < nbuckets_; ++bucket) {
            auto& bucket_buf = thread_bufs.buf[bucket];

            if (bucket_buf.empty()) continue;

            const uint64_t len = bucket_buf.size() * sizeof(Rec);

            // Every non-empty bucket becomes one extent that load() reads back later.
            iov.push_back({bucket_buf.data(), len});
            thread_bufs.extents[bucket].push_back({off, len});
            off += len;
        }

        thread_bufs.file.append(iov);

        for (auto& bucket_buf : thread_bufs.buf) std::vector<Rec>().swap(bucket_buf);

        thread_bufs.held = 0;
    }
    std::size_t nbuckets_;
    uint64_t budget_;
    std::string dir_;
    std::vector<ThreadBuffers> threads_;
};

}
}
