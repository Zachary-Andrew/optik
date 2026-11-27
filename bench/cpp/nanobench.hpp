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

namespace nanobench {
/// PERF_TYPE_HW_CACHE config layout: cache id, then the operation at kCacheOpShift and the result at kCacheResultShift.
constexpr unsigned kCacheOpShift = 8;
constexpr unsigned kCacheResultShift = 16;
/// perf_event_open arguments: any CPU, and no group leader for the first event.
constexpr int kAnyCpu = -1;
constexpr int kNoGroupLeader = -1;
/// A PERF_FORMAT_GROUP read starts with the event count, the enabled time and the running time.
constexpr std::size_t kGroupReadHeaderWords = 3;
/// measure(): reps and warmup of -1 take NB_RUNS / NB_WARMUP, defaulting to kDefaultRuns and kDefaultWarmup.
constexpr int kFromEnv = -1;
constexpr int kDefaultRuns = 31;
constexpr int kDefaultWarmup = 5;
/// Buffer for one event's JSON fields.
constexpr std::size_t kJsonFieldBufLen = 256;


struct Event { const char* name; uint32_t type; uint64_t config; };

inline uint64_t cache_cfg(uint64_t cache_id, uint64_t cache_op, uint64_t res) { return cache_id | (cache_op << kCacheOpShift) | (res << kCacheResultShift); }

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
            perf_event_attr attr;
            std::memset(&attr, 0, sizeof(attr));
            attr.size = sizeof(attr);
            attr.type = evs[i].type;
            attr.config = evs[i].config;
            attr.disabled = i == 0;
            attr.exclude_kernel = 1;
            attr.exclude_hv = 1;
            attr.pinned = i == 0;
            attr.read_format = PERF_FORMAT_GROUP | PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
            const int desc = static_cast<int>(syscall(SYS_perf_event_open, &attr, 0, kAnyCpu, i ? fds_[0] : kNoGroupLeader, 0));
            if (desc < 0) { err_ = std::string(evs[i].name) + ": " + std::strerror(errno); close_all(); return; }
            fds_.push_back(desc);
        }
        ioctl(fds_[0], PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
        ioctl(fds_[0], PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
    }
    ~Group() { close_all(); }
    Group(const Group&) = delete;
    Group& operator=(const Group&) = delete;
    bool usable() const { return err_.empty() && !fds_.empty(); }
    const std::string& error() const { return err_; }
    const std::vector<Event>& events() const { return evs_; }
    /// Current counts, or false if the group was not on the counters the whole time it was enabled, since such counts would need scaling, which is not accepted.
    bool read(std::vector<uint64_t>& counts, bool& full) const {
        std::vector<uint64_t> buf(kGroupReadHeaderWords + evs_.size());
        if (::read(fds_[0], buf.data(), buf.size() * sizeof(uint64_t)) != static_cast<ssize_t>(buf.size() * sizeof(uint64_t))) return false;
        full = buf[1] == buf[2];
        counts.assign(buf.begin() + kGroupReadHeaderWords, buf.end());
        return true;
    }
private:
    void close_all() { for (int desc : fds_) close(desc); fds_.clear(); }
    std::vector<Event> evs_;
    std::vector<int> fds_;
    std::string err_;
};

/// Pin the calling thread to the CPU it runs on, which any user may do.
inline int pin_here() {
    const int cpu = sched_getcpu();
    cpu_set_t cpus; CPU_ZERO(&cpus); CPU_SET(cpu, &cpus);
    sched_setaffinity(0, sizeof(cpus), &cpus);
    return cpu;
}

struct Stat { double median = 0, min = 0, first_quartile = 0, third_quartile = 0; };
struct Result {
    std::string name;
    std::vector<std::pair<std::string, Stat>> per_op;   // event -> per-operation
    std::string error;
};

inline Stat summarize(std::vector<double> samples) {
    if (samples.empty()) return {};
    std::sort(samples.begin(), samples.end());
    return {samples[samples.size() / 2], samples.front(), samples[samples.size() / 4], samples[(3 * samples.size()) / 4]};
}

/// Repetitions and warm-ups of every measure() call, from NB_RUNS and NB_WARMUP in the environment (bench/micro_bench.py --runs and --warmup), else 31 and 5.
inline int env_int(const char* key, int dflt) {
    const char* value = std::getenv(key);
    return value && *value ? std::max(1, std::atoi(value)) : dflt;
}

/// Measure `run(n)`, which performs n operations, per operation.
/// For every event group it takes `warmup` + `reps` measurements of the difference between n and 2n operations, divided by n, and discards the warm-ups.
inline Result measure(const std::string& name, const std::function<void(uint64_t)>& run, uint64_t count,
                      int reps = kFromEnv, int warmup = kFromEnv) {
    if (reps < 0) reps = env_int("NB_RUNS", kDefaultRuns);
    if (warmup < 0) warmup = std::getenv("NB_WARMUP") ? std::max(0, std::atoi(std::getenv("NB_WARMUP"))) : kDefaultWarmup;
    Result result{name, {}, {}};
    std::vector<std::pair<std::string, std::vector<double>>> acc;
    for (const auto& evs : default_groups()) {
        Group group(evs);
        if (!group.usable()) { result.error = group.error(); return result; }
        std::vector<std::vector<double>> per(evs.size());
        for (int it = 0; it < warmup + reps; ++it) {
            std::vector<uint64_t> before1, after1, before2, after2;
            bool full1 = true, full2 = true, full3 = true, full4 = true;
            group.read(before1, full1); run(count); group.read(after1, full2);
            group.read(before2, full3); run(2 * count); group.read(after2, full4);
            if (!(full1 && full2 && full3 && full4)) { result.error = "counters multiplexed (group did not fit)"; return result; }
            if (it < warmup) continue;
            for (std::size_t e = 0; e < evs.size(); ++e) {
                const double one = double(after1[e] - before1[e]), two = double(after2[e] - before2[e]);
                per[e].push_back((two - one) / double(count));
            }
        }
        for (std::size_t e = 0; e < evs.size(); ++e) {
            if (std::string(evs[e].name) == "cycles" && !acc.empty()) continue;   // counted once
            acc.push_back({evs[e].name, per[e]});
        }
    }
    for (auto& [event_name, samples] : acc) result.per_op.push_back({event_name, summarize(samples)});
    return result;
}

inline std::string to_json(const Result& result) {
    std::string json = "{\"kernel\":\"" + result.name + "\"";
    if (!result.error.empty()) return json + ",\"error\":\"" + result.error + "\"}";
    for (const auto& [event_name, stat] : result.per_op) {
        char buf[kJsonFieldBufLen];
        std::snprintf(buf, sizeof buf, ",\"%s\":%.4f,\"%s_min\":%.4f,\"%s_q1\":%.4f,\"%s_q3\":%.4f", event_name.c_str(),
                      stat.median, event_name.c_str(), stat.min, event_name.c_str(), stat.first_quartile, event_name.c_str(), stat.third_quartile);
        json += buf;
    }
    return json + "}";
}

/// Keep a value alive without letting the compiler drop the work.
template <typename T> inline void keep(const T& value) { asm volatile("" : : "r,m"(value) : "memory"); }

}
