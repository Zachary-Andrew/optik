#pragma once

/*
 * delta_oa.hpp - concurrent open-addressing set of 64-bit keys.
 *
 * DeltaOA holds the canonical k-mers added to, or erased from, a static SPSS core, forming the dynamic layer of DynSpss in dyn_spss.hpp.
 *
 * Layout.
 * The table is an array of 64-byte lines,
 *
 *     u32 version | u8 count | u8 flags | u16 pad | 7 x u64 keys
 *
 * so one probe reads one cache line and compares its 7 keys in registers.
 *
 * Placement.
 * Level 0 has L lines and level i has L / 2^i.
 * A key hashes to two lines on every level.
 * An insert goes to the emptier of its two level-0 lines, to level 1 only when both are full, and so on down.
 * When an insert passes over a full pair it sets kOverflowFlag on both lines, and a lookup descends past a level only if one of its two lines carries the flag, so a negative lookup usually reads 2 lines.
 *
 * The geometrically shrinking sub-arrays are those of Farach-Colton, Krapivin and Kuszmaul, "Optimal bounds for open addressing without reordering" (FOCS 2024).
 * In this adaptation the sub-arrays are made of cache lines rather than slots, and the insertion policy is the two-choice rule above, not the paper's.
 * As in the paper, a key never moves after insertion.
 * That last property is what keeps the lock-free readers simple.
 *
 * Concurrency.
 *   Readers never block.  Each line is a seqlock, so a reader reads the
 *   version, retrying while it is odd, reads the keys, and retries if the
 *   version changed.
 *
 *   Writers lock the key's two level-0 lines in address order (an odd
 *   version means locked).  Every insert or erase of a key takes those same
 *   two locks first, so writers of one key are serialised (no duplicate
 *   inserts) while writers of different keys run in parallel.  Deeper lines
 *   are locked afterwards, in (level, address) order.
 *
 *   Resizing.  Writers hold the writer mutex shared, and grow() takes it
 *   exclusively, copies every key into a table twice as large and publishes
 *   it with a release store.  The old table is freed after a grace period in
 *   the style of the Linux kernel's SRCU, where readers register in one of
 *   two per-slot generation counters and grow() flips the generation and
 *   waits for the old one to drain.  Readers never wait on writers.
 *
 * The keys kEmpty (~0) and kTomb (~0 - 1) are reserved, and neither is a valid canonical k-mer.
 */

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <thread>
#include <vector>

#include "wyhash.hpp"

namespace tpoptoa {

class DeltaOA {
public:
    static constexpr uint64_t kEmpty = ~uint64_t(0);
    static constexpr uint64_t kTomb  = ~uint64_t(0) - 1;
    static constexpr unsigned kSlots = 7;               // keys per 64-byte line

    explicit DeltaOA(std::size_t expected = 1024) {
        tab_.store(new Table(lines_for(expected)), std::memory_order_release);
    }
    ~DeltaOA() { delete tab_.load(std::memory_order_acquire); }
    DeltaOA(const DeltaOA&) = delete;
    DeltaOA& operator=(const DeltaOA&) = delete;

    std::size_t size() const noexcept { return size_.load(std::memory_order_relaxed); }
    bool empty() const noexcept { return size() == 0; }
    std::size_t bytes() const noexcept {
        const Table* t = tab_.load(std::memory_order_acquire);
        return t->nlines_total * sizeof(Line);
    }

    /// contains() - lock-free membership test.
    bool contains(uint64_t key) const noexcept {
        ReadGuard guard(*this);
        return find_in(guard.t, key);
    }

    /// insert() - add @key, returning true if it was absent.
    bool insert(uint64_t key) {
        for (;;) {
            {
                std::shared_lock<std::shared_mutex> writer_lock(writer_mutex_);
                Table* t = tab_.load(std::memory_order_acquire);
                int r = insert_in(t, key);
                if (r >= 0) {
                    if (r) size_.fetch_add(1, std::memory_order_relaxed);
                    return r == 1;
                }
            }
            grow();                                   // every level is full, so resize
        }
    }

    /// erase() - remove @key, leaving a tombstone, and return true if it was present.
    bool erase(uint64_t key) {
        std::shared_lock<std::shared_mutex> writer_lock(writer_mutex_);
        Table* t = tab_.load(std::memory_order_acquire);
        const std::size_t a = t->line_of(0, key, 0), b = t->line_of(0, key, 1);
        LockPair level0_locks(t->lines[std::min(a, b)], t->lines[std::max(a, b)], a == b);
        for (unsigned lvl = 0; lvl < t->levels; ++lvl) {
            const std::size_t x = t->line_of(lvl, key, 0), y = t->line_of(lvl, key, 1);
            for (std::size_t li : {x, y}) {
                Line& L = t->lines[li];
                const bool already_locked = (lvl == 0);
                if (!already_locked) lock(L);
                for (unsigned s = 0; s < kSlots; ++s) {
                    if (L.keys[s].load(std::memory_order_relaxed) == key) {
                        L.keys[s].store(kTomb, std::memory_order_relaxed);
                        if (!already_locked) unlock(L);
                        size_.fetch_sub(1, std::memory_order_relaxed);
                        return true;
                    }
                }
                if (!already_locked) unlock(L);
            }
            if (!(t->lines[x].flags.load(std::memory_order_relaxed) & kOverflowFlag) &&
                !(t->lines[y].flags.load(std::memory_order_relaxed) & kOverflowFlag)) break;
        }
        return false;
    }

    /// for_each() - call @f with every key.
    /// Must not run concurrently with writers.
    template <typename F> void for_each(F&& f) const {
        const Table* t = tab_.load(std::memory_order_acquire);
        for (std::size_t i = 0; i < t->nlines_total; ++i)
            for (unsigned s = 0; s < kSlots; ++s) {
                const uint64_t k = t->lines[i].keys[s].load(std::memory_order_relaxed);
                if (k != kEmpty && k != kTomb) f(k);
            }
    }

private:
    /// Line::flags bit set when an insert found both lines of this pair full and went one level deeper.
    static constexpr uint8_t kOverflowFlag = 1;

    /// Smallest level-0 size, in lines.
    static constexpr std::size_t kMinLevel0Lines = 16;

    /// Odd multiplier (2^64 / golden ratio) that gives each (level, choice) pair its own seed.
    static constexpr uint64_t kGoldenRatio64 = 0x9E3779B97F4A7C15ULL;

    struct alignas(64) Line {
        std::atomic<uint32_t> version{0};           // seqlock version, odd while locked
        std::atomic<uint8_t> count{0};              // slots ever used, tombstones included
        std::atomic<uint8_t> flags{0};              // kOverflowFlag
        uint16_t pad = 0;
        std::atomic<uint64_t> keys[kSlots];
        Line() { for (auto& k : keys) k.store(kEmpty, std::memory_order_relaxed); }
    };
    static_assert(sizeof(Line) == 64, "a line must be one cache line");

    struct Table {
        std::size_t l0 = 0;                          // lines at level 0, a power of two
        unsigned levels = 0;
        std::size_t nlines_total = 0;
        std::vector<std::size_t> base;               // first line of each level
        std::unique_ptr<Line[]> lines;
        explicit Table(std::size_t level0_lines) : l0(level0_lines) {
            std::size_t n = l0, total = 0;
            while (n >= 1) { base.push_back(total); total += n; ++levels; if (n == 1) break; n >>= 1; }
            nlines_total = total;
            lines.reset(new Line[total]);
        }
        /// line_of() - index of line @choice, 0 or 1, of @key on level @lvl.
        std::size_t line_of(unsigned lvl, uint64_t key, unsigned choice) const noexcept {
            const std::size_t n = l0 >> lvl;
            const uint64_t h = mix64(key ^ (kGoldenRatio64 * (2 * lvl + choice + 1)), 0);
            return base[lvl] + static_cast<std::size_t>((static_cast<__uint128_t>(h) * n) >> 64);
        }
    };

    static void lock(Line& L) noexcept {
        for (;;) {
            uint32_t v = L.version.load(std::memory_order_relaxed);
            if (!(v & 1) && L.version.compare_exchange_weak(v, v + 1, std::memory_order_acquire))
                return;
            std::this_thread::yield();
        }
    }
    static void unlock(Line& L) noexcept {
        L.version.fetch_add(1, std::memory_order_release);
    }

    /// LockPair - scoped lock of two lines, given in address order, with @same set when both choices hash to one line.
    struct LockPair {
        Line& a; Line& b; bool same;
        LockPair(Line& x, Line& y, bool s) : a(x), b(y), same(s) { lock(a); if (!same) lock(b); }
        ~LockPair() { if (!same) unlock(b); unlock(a); }
    };

    /// line_has() - seqlock read of one line, telling whether it holds @key and, through @overflow, whether its kOverflowFlag is set.
    static bool line_has(const Line& L, uint64_t key, bool& overflow) noexcept {
        for (;;) {
            const uint32_t v1 = L.version.load(std::memory_order_acquire);
            if (v1 & 1) { std::this_thread::yield(); continue; }
            bool hit = false;
            for (unsigned s = 0; s < kSlots; ++s)
                hit |= L.keys[s].load(std::memory_order_relaxed) == key;
            const bool overflowed = L.flags.load(std::memory_order_relaxed) & kOverflowFlag;
            std::atomic_thread_fence(std::memory_order_acquire);
            if (L.version.load(std::memory_order_relaxed) == v1) { overflow = overflowed; return hit; }
        }
    }

    static bool find_in(const Table* t, uint64_t key) noexcept {
        for (unsigned lvl = 0; lvl < t->levels; ++lvl) {
            bool overflow0 = false, overflow1 = false;
            if (line_has(t->lines[t->line_of(lvl, key, 0)], key, overflow0)) return true;
            if (line_has(t->lines[t->line_of(lvl, key, 1)], key, overflow1)) return true;
            if (!overflow0 && !overflow1) return false;
        }
        return false;
    }

    /// insert_in() - returns 1 if inserted, 0 if already present, -1 if every level is full.
    int insert_in(Table* t, uint64_t key) {
        const std::size_t a = t->line_of(0, key, 0), b = t->line_of(0, key, 1);
        LockPair level0_locks(t->lines[std::min(a, b)], t->lines[std::max(a, b)], a == b);
        // presence check on every level the key could occupy (writers of the same key are serialised by the level-0 locks)
        {
            for (unsigned lvl = 0; lvl < t->levels; ++lvl) {
                const std::size_t x = t->line_of(lvl, key, 0), y = t->line_of(lvl, key, 1);
                for (std::size_t li : {x, y})
                    for (unsigned s = 0; s < kSlots; ++s)
                        if (t->lines[li].keys[s].load(std::memory_order_relaxed) == key) return 0;
                if (!(t->lines[x].flags.load(std::memory_order_relaxed) & kOverflowFlag) &&
                    !(t->lines[y].flags.load(std::memory_order_relaxed) & kOverflowFlag)) break;
            }
        }
        for (unsigned lvl = 0; lvl < t->levels; ++lvl) {
            const std::size_t x = t->line_of(lvl, key, 0), y = t->line_of(lvl, key, 1);
            Line& X = t->lines[x];
            Line& Y = t->lines[y];
            const bool already_locked = (lvl == 0);
            if (!already_locked) { lock(x < y ? X : Y); if (x != y) lock(x < y ? Y : X); }
            // first free or tombstoned slot of the line
            auto place = [&](Line& L) -> bool {
                for (unsigned s = 0; s < kSlots; ++s) {
                    const uint64_t k = L.keys[s].load(std::memory_order_relaxed);
                    if (k == kEmpty || k == kTomb) {
                        L.keys[s].store(key, std::memory_order_relaxed);
                        if (k == kEmpty) L.count.fetch_add(1, std::memory_order_relaxed);
                        return true;
                    }
                }
                return false;
            };
            const bool x_first = X.count.load(std::memory_order_relaxed) <= Y.count.load(std::memory_order_relaxed);
            bool ok = x_first ? (place(X) || place(Y)) : (place(Y) || place(X));
            if (!ok) {                                   // both full, so mark and descend
                X.flags.fetch_or(kOverflowFlag, std::memory_order_relaxed);
                Y.flags.fetch_or(kOverflowFlag, std::memory_order_relaxed);
            }
            if (!already_locked) { if (x != y) unlock(x < y ? Y : X); unlock(x < y ? X : Y); }
            if (ok) return 1;
        }
        return -1;
    }

    // Reader registration for the resize grace period.
    // Each thread uses one of kReaderSlots counter pairs, chosen by arrival order, with one counter per generation parity.
    static constexpr unsigned kReaderSlots = 64;
    struct alignas(64) ReaderSlot { std::atomic<uint32_t> n[2]; ReaderSlot() { n[0] = 0; n[1] = 0; } };
    static unsigned reader_slot() noexcept {
        static std::atomic<unsigned> next{0};
        thread_local unsigned s = next.fetch_add(1, std::memory_order_relaxed) % kReaderSlots;
        return s;
    }

    /// ReadGuard - registers the calling reader in the current generation and pins the table it then loads until destruction.
    struct ReadGuard {
        const DeltaOA& d; unsigned slot, gen; const Table* t;
        explicit ReadGuard(const DeltaOA& owner) : d(owner), slot(reader_slot()) {
            for (;;) {
                gen = d.generation_.load(std::memory_order_acquire) & 1;
                d.reader_slots_[slot].n[gen].fetch_add(1, std::memory_order_seq_cst);
                if ((d.generation_.load(std::memory_order_seq_cst) & 1) == gen) break;
                d.reader_slots_[slot].n[gen].fetch_sub(1, std::memory_order_release);
            }
            t = d.tab_.load(std::memory_order_acquire);
        }
        ~ReadGuard() { d.reader_slots_[slot].n[gen].fetch_sub(1, std::memory_order_release); }
    };

    /// grow() - replace the table by one with twice the level-0 lines (more if the keys still do not fit), then free the old table once no reader can hold it.
    void grow() {
        std::unique_lock<std::shared_mutex> writer_lock(writer_mutex_);
        Table* old = tab_.load(std::memory_order_acquire);
        // Another writer may have grown the table already.
        // Growing twice costs memory only, and the caller retries its insert either way.
        std::size_t l0 = old->l0 * 2;
        for (;;) {
            Table* bigger = new Table(l0);
            bool ok = true;
            for (std::size_t i = 0; i < old->nlines_total && ok; ++i)
                for (unsigned s = 0; s < kSlots && ok; ++s) {
                    const uint64_t k = old->lines[i].keys[s].load(std::memory_order_relaxed);
                    if (k != kEmpty && k != kTomb) ok = insert_in(bigger, k) == 1;
                }
            if (ok) {
                tab_.store(bigger, std::memory_order_release);
                // Grace period, flipping the generation and waiting for old-generation readers
                const uint32_t g = generation_.fetch_add(1, std::memory_order_seq_cst) & 1;
                for (unsigned i = 0; i < kReaderSlots; ++i)
                    while (reader_slots_[i].n[g].load(std::memory_order_acquire) != 0) std::this_thread::yield();
                delete old;
                return;
            }
            delete bigger;
            l0 *= 2;
        }
    }

    /// lines_for() - level-0 line count for @expected keys, a power of two (at least kMinLevel0Lines) giving about kSlots / 2 keys per line.
    static std::size_t lines_for(std::size_t expected) {
        std::size_t need = expected / (kSlots / 2) + 1, l0 = kMinLevel0Lines;
        while (l0 < need) l0 <<= 1;
        return l0;
    }

    std::atomic<Table*> tab_{nullptr};
    std::atomic<std::size_t> size_{0};
    std::shared_mutex writer_mutex_;
    mutable std::atomic<uint32_t> generation_{0};
    mutable ReaderSlot reader_slots_[kReaderSlots];
};

}
