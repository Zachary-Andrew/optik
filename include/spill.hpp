#pragma once

/*
 * spill.hpp - disk spill for the k-mer counter's partition buffers.
 *
 * KMC 2 (Deorowicz et al. 2015) writes super-k-mers to disk bins by minimizer and then counts one bin at a time.
 * optik keeps its bins in memory while they fit and only writes them to disk past the memory budget.
 * The bins are then read back one at a time, as in KMC.
 *
 * Memory budget.
 * The user's --max-mem, or else half of what the process may still take, i.e. the smaller of MemAvailable and the tightest cgroup v2 room (memory.max - memory.current) from this process's cgroup up to the root, which covers a systemd-run -p MemoryMax=... scope or a container.
 *
 * Spill files.
 * One file per worker thread, opened with O_TMPFILE (or created and unlinked at once), so a crash leaves nothing behind.
 * A worker appends whole partition buffers with pwritev, and the extents it wrote are indexed per bin, a run of consecutive partitions, so a bin is read back with a few large preads.
 * Short writes, ENOSPC and EIO throw.
 * A bin read for the last time is punched out of its file (FALLOC_FL_PUNCH_HOLE), so disk use falls as the counting pass advances.
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

/// CountLimits - process-wide counter limits, set by the CLI before a build (ok_set_count_limits() in src/ffi.cpp).
struct CountLimits {
    uint64_t max_mem = 0;      // bytes, where 0 means half of what is available
    std::string tmp_dir;       // spill directory (the CLI's --tmp-dir, else the output's), where empty means "."
    int graph_mode = 0;        // --graph mode, 0 auto, 1 partitioned, 2 table
};
inline CountLimits& count_limits() { static CountLimits l; return l; }

namespace spill {

/// Budget assumed when neither /proc/meminfo nor a cgroup gives a limit.
inline constexpr uint64_t kFallbackAvailableBytes = uint64_t(8) << 30;

/// statfs f_type values of RAM-backed filesystems (linux/magic.h).
inline constexpr long kTmpfsMagicNumber = 0x01021994;
inline constexpr long kRamfsMagicNumber = 0x858458f6;

/// read_u64_file() - the first token of @path as an integer, or UINT64_MAX when the file is missing, unparsable or says "max", a cgroup without a limit.
inline uint64_t read_u64_file(const std::string& path) {
    std::ifstream f(path);
    std::string s;
    if (!(f >> s) || s == "max") return UINT64_MAX;
    try { return std::stoull(s); } catch (...) { return UINT64_MAX; }
}

/// available_bytes() - bytes this process may still allocate, i.e. MemAvailable bounded by every cgroup v2 limit from its own cgroup up to the root.
inline uint64_t available_bytes() {
    uint64_t avail = UINT64_MAX;
    {
        std::ifstream f("/proc/meminfo");
        std::string key; uint64_t kib; std::string unit;
        while (f >> key >> kib) {
            std::getline(f, unit);
            if (key == "MemAvailable:") { avail = kib * 1024; break; }
        }
    }
    // the cgroup v2 entry of /proc/self/cgroup is "0::<path>"
    std::ifstream cgroup_file("/proc/self/cgroup");
    std::string line, cgroup_path;
    while (std::getline(cgroup_file, line))
        if (line.rfind("0::", 0) == 0) { cgroup_path = line.substr(3); break; }
    if (!cgroup_path.empty()) {
        for (std::string p = cgroup_path;;) {
            const std::string dir = "/sys/fs/cgroup" + (p == "/" ? std::string() : p);
            const uint64_t limit = read_u64_file(dir + "/memory.max");
            if (limit != UINT64_MAX) {
                const uint64_t used = read_u64_file(dir + "/memory.current");
                const uint64_t room = (used != UINT64_MAX && used < limit) ? limit - used : 0;
                if (room < avail) avail = room;
            }
            if (p == "/" || p.empty()) break;
            const auto slash = p.find_last_of('/');
            p = slash == 0 || slash == std::string::npos ? "/" : p.substr(0, slash);
        }
    }
    return avail == UINT64_MAX ? kFallbackAvailableBytes : avail;
}

/// resident_bytes() - resident set of this process, in bytes.
inline uint64_t resident_bytes() {
    long pages = 0, resident = 0;
    if (std::FILE* f = std::fopen("/proc/self/statm", "r")) {
        if (std::fscanf(f, "%ld %ld", &pages, &resident) != 2) resident = 0;
        std::fclose(f);
    }
    return static_cast<uint64_t>(resident) * static_cast<uint64_t>(sysconf(_SC_PAGESIZE));
}

/// memory_budget() - the memory the counter plans with, --max-mem or half of available_bytes().
/// @explicit_limit says which.
inline uint64_t memory_budget(bool& explicit_limit) {
    const uint64_t m = count_limits().max_mem;
    explicit_limit = m != 0;
    return m ? m : available_bytes() / 2;
}

/// free_bytes() - bytes free to this user on the filesystem holding @dir (0 if unknown).
inline uint64_t free_bytes(const std::string& dir) {
    struct statvfs v;
    if (statvfs(dir.c_str(), &v) != 0) return 0;
    return static_cast<uint64_t>(v.f_bavail) * static_cast<uint64_t>(v.f_frsize);
}

/// spill_dir() - the --tmp-dir setting, or "." when unset.
inline std::string spill_dir() {
    const std::string& d = count_limits().tmp_dir;
    return d.empty() ? "." : d;
}

/// ram_backed() - true if @dir is on tmpfs or ramfs, where spilling saves no memory.
inline bool ram_backed(const std::string& dir) {
    struct statfs s;
    if (statfs(dir.c_str(), &s) != 0) return false;
    return s.f_type == kTmpfsMagicNumber || s.f_type == static_cast<decltype(s.f_type)>(kRamfsMagicNumber);
}

/// File - one worker's anonymous spill file, with appends and positional reads.
class File {
public:
    File() = default;
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    ~File() { if (fd_ >= 0) ::close(fd_); }

    void open(const std::string& dir) {
        if (fd_ >= 0) return;
#ifdef O_TMPFILE
        fd_ = ::open(dir.c_str(), O_TMPFILE | O_RDWR | O_CLOEXEC, 0600);
#endif
        if (fd_ < 0) {   // filesystems without O_TMPFILE get a file created, then unlinked
            std::string tmpl = dir + "/optik_spill_XXXXXX";
            std::vector<char> b(tmpl.begin(), tmpl.end()); b.push_back('\0');
            fd_ = ::mkstemp(b.data());
            if (fd_ >= 0) ::unlink(b.data());
        }
        if (fd_ < 0)
            throw std::runtime_error("cannot create a spill file in '" + dir + "': " +
                                     std::strerror(errno) + " (set --tmp-dir)");
    }
    bool is_open() const noexcept { return fd_ >= 0; }
    uint64_t size() const noexcept { return end_; }

    /// append() - write every buffer of @iov at the end of the file, retrying short writes and EINTR.
    /// @iov is consumed.
    void append(std::vector<struct iovec>& iov) {
        std::size_t i = 0;
        while (i < iov.size()) {
            const int cnt = static_cast<int>(std::min<std::size_t>(iov.size() - i, IOV_MAX));
            const ssize_t w = ::pwritev(fd_, iov.data() + i, cnt, static_cast<off_t>(end_));
            if (w < 0) {
                if (errno == EINTR) continue;
                throw std::runtime_error(std::string("spill write failed: ") + std::strerror(errno) +
                                         (errno == ENOSPC ? " (the spill directory's disk is full: set --tmp-dir "
                                                            "to a larger disk, or raise --max-mem so less is spilled)"
                                                          : ""));
            }
            if (w == 0) throw std::runtime_error("spill write failed: no progress");
            end_ += static_cast<uint64_t>(w);
            std::size_t left = static_cast<std::size_t>(w);
            while (i < iov.size() && left >= iov[i].iov_len) { left -= iov[i].iov_len; ++i; }
            if (left) {
                iov[i].iov_base = static_cast<char*>(iov[i].iov_base) + left;
                iov[i].iov_len -= left;
            }
        }
    }
    /// read() - read exactly @len bytes at @off into @dst.
    void read(uint64_t off, uint64_t len, uint8_t* dst) const {
        while (len) {
            const ssize_t r = ::pread(fd_, dst, len, static_cast<off_t>(off));
            if (r < 0) {
                if (errno == EINTR) continue;
                throw std::runtime_error(std::string("spill read failed: ") + std::strerror(errno));
            }
            if (r == 0) throw std::runtime_error("spill read failed: file shorter than written");
            dst += r; off += static_cast<uint64_t>(r); len -= static_cast<uint64_t>(r);
        }
    }
    /// punch() - free the blocks of [@off, @off + @len), which will not be read again.
    /// The file keeps its size.
    /// On a filesystem without hole punching the call fails silently and the blocks are freed at close.
    void punch(uint64_t off, uint64_t len) const {
        if (fd_ >= 0 && len)
            (void)::fallocate(fd_, FALLOC_FL_PUNCH_HOLE | FALLOC_FL_KEEP_SIZE,
                              static_cast<off_t>(off), static_cast<off_t>(len));
    }
    /// will_need() - ask the kernel to read [@off, @off + @len) ahead.
    void will_need(uint64_t off, uint64_t len) const {
        if (fd_ >= 0) ::posix_fadvise(fd_, static_cast<off_t>(off), static_cast<off_t>(len), POSIX_FADV_WILLNEED);
    }
private:
    int fd_ = -1;
    uint64_t end_ = 0;
};

/// BucketSpill - records distributed over buckets by many threads.
///
/// Each thread keeps its records in memory up to a per-thread byte budget and appends them to its own spill file past it.
/// load(b) gathers bucket b from every thread's memory and file.
/// Rec must be trivially copyable.
template <typename Rec>
class BucketSpill {
public:
    BucketSpill(unsigned nthreads, std::size_t nbuckets, uint64_t per_thread_bytes, std::string dir)
        : nbuckets_(nbuckets), budget_(per_thread_bytes), dir_(std::move(dir)), threads_(nthreads) {
        for (auto& t : threads_) { t.buf.resize(nbuckets_); t.extents.resize(nbuckets_); }
    }
    /// add() - thread @tid adds record @r to bucket @b, which may write the thread's buffers to disk.
    void add(unsigned tid, std::size_t b, const Rec& r) {
        ThreadBuffers& t = threads_[tid];
        auto& v = t.buf[b];
        const std::size_t cap0 = v.capacity();
        v.push_back(r);
        if (v.capacity() != cap0) t.held += (v.capacity() - cap0) * sizeof(Rec);
        if (t.held > budget_) flush(t);
    }
    bool spilled() const noexcept {
        for (const ThreadBuffers& t : threads_) if (t.file.is_open()) return true;
        return false;
    }
    /// load() - every record of bucket @b into @out, releasing the bucket's memory.
    void load(std::size_t b, std::vector<Rec>& out) {
        out.clear();
        std::size_t n = 0;
        for (ThreadBuffers& t : threads_) {
            n += t.buf[b].size();
            for (const auto& e : t.extents[b]) n += e.second / sizeof(Rec);
        }
        out.resize(n);
        std::size_t at = 0;
        for (ThreadBuffers& t : threads_) {
            for (const auto& e : t.extents[b]) {
                t.file.read(e.first, e.second, reinterpret_cast<uint8_t*>(out.data() + at));
                at += e.second / sizeof(Rec);
            }
            std::copy(t.buf[b].begin(), t.buf[b].end(), out.begin() + static_cast<std::ptrdiff_t>(at));
            at += t.buf[b].size();
            std::vector<Rec>().swap(t.buf[b]);
        }
    }
    /// size() - records in bucket @b, in memory and on disk.
    std::size_t size(std::size_t b) const {
        std::size_t n = 0;
        for (const ThreadBuffers& t : threads_) {
            n += t.buf[b].size();
            for (const auto& e : t.extents[b]) n += e.second / sizeof(Rec);
        }
        return n;
    }
    /// will_need() - ask the kernel to read bucket @b's disk parts ahead.
    void will_need(std::size_t b) const {
        for (const ThreadBuffers& t : threads_) for (const auto& e : t.extents[b]) t.file.will_need(e.first, e.second);
    }

private:
    struct ThreadBuffers {
        std::vector<std::vector<Rec>> buf;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> extents;   // (offset, bytes) per bucket
        uint64_t held = 0;                                                // bytes of buf capacity
        File file;
    };
    /// flush() - append all of @t's buffers to its spill file and free them.
    void flush(ThreadBuffers& t) {
        t.file.open(dir_);
        std::vector<struct iovec> iov;
        uint64_t off = t.file.size();
        for (std::size_t b = 0; b < nbuckets_; ++b) {
            auto& v = t.buf[b];
            if (v.empty()) continue;
            const uint64_t len = v.size() * sizeof(Rec);
            iov.push_back({v.data(), len});
            t.extents[b].push_back({off, len});
            off += len;
        }
        t.file.append(iov);
        for (auto& v : t.buf) std::vector<Rec>().swap(v);
        t.held = 0;
    }
    std::size_t nbuckets_;
    uint64_t budget_;
    std::string dir_;
    std::vector<ThreadBuffers> threads_;
};

}
}
