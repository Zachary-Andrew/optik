#pragma once

/// Unprivileged nanoBench-style microbenchmarks, without sudo or a kernel module.
///
/// The method is nanoBench's (Abel and Reineke, ISPASS 2020).
/// The code under test runs as a sequence of n operations with hardware performance counters read around it, and the difference between the lengths n and 2n cancels the cost of reading the counters and of the loop scaffolding.
/// This is repeated, the warm-up runs are discarded, and the median and minimum per operation are reported.
///
/// nanoBench programs the counters by writing model-specific registers and enables user rdpmc through sysfs, both of which need root.
/// Here the counters come from perf_event_open on the calling thread with exclude_kernel and exclude_hv set, which the common default kernel.perf_event_paranoid <= 2 allows to any user.
///   - Events are opened as groups led by cycles, pinned, with at most MAX_GROUP per group, so they are counted together and never multiplexed.
///     A group whose time_running < time_enabled (the kernel could not schedule it on the counters) is rejected, not scaled.
///   - The thread is pinned to its current CPU (sched_setaffinity, allowed to any user), so one core's counters see the whole run.
///   - Counts are read with one read() of the group before and one after.
///
/// Calibration kernels (dependent and independent add chains, an empty operation, and L1 and DRAM pointer chases) check that the subtraction removes the overhead and that the counters report what the code does.

#include <algorithm>
#include <cstdlib>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include <linux/perf_event.h>
#include <sched.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>

namespace nb {

struct Event { const char* name; uint32_t type; uint64_t config; };

inline uint64_t cache_cfg(uint64_t id, uint64_t op, uint64_t res) { return id | (op << 8) | (res << 16); }

/// Events grouped so that each group fits the core's counters (AMD Zen has 6 general-purpose counters, one of which the NMI watchdog may hold).
inline std::vector<std::vector<Event>> default_groups() {
    return {
        {{"cycles", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
         {"instructions", PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS},
         {"branch_misses", PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_MISSES},
         {"l1d_misses", PERF_TYPE_HW_CACHE, cache_cfg(PERF_COUNT_HW_CACHE_L1D, PERF_COUNT_HW_CACHE_OP_READ,
                                                    PERF_COUNT_HW_CACHE_RESULT_MISS)}},
        {{"cycles", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
         {"llc_misses", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES},
         {"dtlb_misses", PERF_TYPE_HW_CACHE, cache_cfg(PERF_COUNT_HW_CACHE_DTLB, PERF_COUNT_HW_CACHE_OP_READ,
                                                     PERF_COUNT_HW_CACHE_RESULT_MISS)}},
    };
}

class Group {
public:
    explicit Group(const std::vector<Event>& evs) : evs_(evs) {
        for (std::size_t i = 0; i < evs.size(); ++i) {
            perf_event_attr a;
            std::memset(&a, 0, sizeof(a));
            a.size = sizeof(a);
            a.type = evs[i].type;
            a.config = evs[i].config;
            a.disabled = i == 0;
            a.exclude_kernel = 1;
            a.exclude_hv = 1;
            a.pinned = i == 0;
            a.read_format = PERF_FORMAT_GROUP | PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
            const int fd = static_cast<int>(syscall(SYS_perf_event_open, &a, 0, -1, i ? fds_[0] : -1, 0));
            if (fd < 0) { err_ = std::string(evs[i].name) + ": " + std::strerror(errno); close_all(); return; }
            fds_.push_back(fd);
        }
        ioctl(fds_[0], PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
        ioctl(fds_[0], PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
    }
    ~Group() { close_all(); }
    Group(const Group&) = delete;
    Group& operator=(const Group&) = delete;
    bool ok() const { return err_.empty() && !fds_.empty(); }
    const std::string& error() const { return err_; }
    const std::vector<Event>& events() const { return evs_; }
    /// Current counts, or false if the group was not on the counters the whole time it was enabled, since such counts would need scaling, which is not accepted.
    bool read(std::vector<uint64_t>& v, bool& full) const {
        std::vector<uint64_t> buf(3 + evs_.size());
        if (::read(fds_[0], buf.data(), buf.size() * 8) != static_cast<ssize_t>(buf.size() * 8)) return false;
        full = buf[1] == buf[2];
        v.assign(buf.begin() + 3, buf.end());
        return true;
    }
private:
    void close_all() { for (int fd : fds_) close(fd); fds_.clear(); }
    std::vector<Event> evs_;
    std::vector<int> fds_;
    std::string err_;
};

/// Pin the calling thread to the CPU it runs on, which any user may do.
inline int pin_here() {
    const int cpu = sched_getcpu();
    cpu_set_t s; CPU_ZERO(&s); CPU_SET(cpu, &s);
    sched_setaffinity(0, sizeof(s), &s);
    return cpu;
}

struct Stat { double median = 0, min = 0, q1 = 0, q3 = 0; };
struct Result {
    std::string name;
    std::vector<std::pair<std::string, Stat>> per_op;   // event -> per-operation
    std::string error;
};

inline Stat summarize(std::vector<double> v) {
    if (v.empty()) return {};
    std::sort(v.begin(), v.end());
    return {v[v.size() / 2], v.front(), v[v.size() / 4], v[(3 * v.size()) / 4]};
}

/// Repetitions and warm-ups of every measure() call, from NB_RUNS and NB_WARMUP in the environment (bench/micro_bench.py --runs and --warmup), else 31 and 5.
inline int env_int(const char* k, int dflt) {
    const char* v = std::getenv(k);
    return v && *v ? std::max(1, std::atoi(v)) : dflt;
}

/// Measure `run(n)`, which performs n operations, per operation.
/// For every event group it takes `warmup` + `reps` measurements of the difference between n and 2n operations, divided by n, and discards the warm-ups.
inline Result measure(const std::string& name, const std::function<void(uint64_t)>& run, uint64_t n,
                      int reps = -1, int warmup = -1) {
    if (reps < 0) reps = env_int("NB_RUNS", 31);
    if (warmup < 0) warmup = std::getenv("NB_WARMUP") ? std::max(0, std::atoi(std::getenv("NB_WARMUP"))) : 5;
    Result r{name, {}, {}};
    std::vector<std::pair<std::string, std::vector<double>>> acc;
    for (const auto& evs : default_groups()) {
        Group g(evs);
        if (!g.ok()) { r.error = g.error(); return r; }
        std::vector<std::vector<double>> per(evs.size());
        for (int it = 0; it < warmup + reps; ++it) {
            std::vector<uint64_t> a, b, c, d;
            bool f1 = true, f2 = true, f3 = true, f4 = true;
            g.read(a, f1); run(n); g.read(b, f2);
            g.read(c, f3); run(2 * n); g.read(d, f4);
            if (!(f1 && f2 && f3 && f4)) { r.error = "counters multiplexed (group did not fit)"; return r; }
            if (it < warmup) continue;
            for (std::size_t e = 0; e < evs.size(); ++e) {
                const double one = double(b[e] - a[e]), two = double(d[e] - c[e]);
                per[e].push_back((two - one) / double(n));
            }
        }
        for (std::size_t e = 0; e < evs.size(); ++e) {
            if (std::string(evs[e].name) == "cycles" && !acc.empty()) continue;   // counted once
            acc.push_back({evs[e].name, per[e]});
        }
    }
    for (auto& [ev, v] : acc) r.per_op.push_back({ev, summarize(v)});
    return r;
}

inline std::string to_json(const Result& r) {
    std::string s = "{\"kernel\":\"" + r.name + "\"";
    if (!r.error.empty()) return s + ",\"error\":\"" + r.error + "\"}";
    for (const auto& [ev, st] : r.per_op) {
        char b[256];
        std::snprintf(b, sizeof b, ",\"%s\":%.4f,\"%s_min\":%.4f,\"%s_q1\":%.4f,\"%s_q3\":%.4f", ev.c_str(),
                      st.median, ev.c_str(), st.min, ev.c_str(), st.q1, ev.c_str(), st.q3);
        s += b;
    }
    return s + "}";
}

/// Keep a value alive without letting the compiler drop the work.
template <typename T> inline void keep(const T& v) { asm volatile("" : : "r,m"(v) : "memory"); }

}
