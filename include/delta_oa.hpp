#pragma once

/**
 * @file
 * @brief Concurrent open-addressing set of 64-bit keys.
 *
 * DeltaOA holds the canonical k-mers added to, or erased from, a static SPSS core.
 * In other words, it is the editable layer of DynSpss (dyn_spss.hpp).
 *
 * @par Layout
 * The table is an array of 64-byte lines:
 *
 *     u32 version | u8 count | u8 flags | u16 pad | 7 x u64 keys
 *
 * One probe therefore reads one cache line and compares all 7 keys in registers.
 *
 * @par Placement
 * Level 0 has L lines, and level i has L / 2^i.
 * On every level a key hashes to two candidate lines.
 * An insert takes the emptier of its two level-0 lines and only moves to level 1 when both are full, and so on down.
 * Passing over a full pair sets kOverflowFlag on both lines.
 * Lookups descend only past a flagged pair, so a negative lookup usually reads just 2 lines.
 *
 * The geometrically shrinking levels come from Farach-Colton, Krapivin and Kuszmaul, "Optimal bounds for open addressing without reordering" (FOCS 2024).
 * Our version builds the levels from cache lines instead of slots and uses the two-choice rule above rather than the paper's policy.
 *
 * @par Concurrency
 * Readers never block, because each line doubles as a seqlock.
 * A reader loads the version, retrying while it is odd, then reads the keys and retries if the version moved in the meantime.
 *
 * Writers lock the key's two level-0 lines in address order, and an odd version marks a held lock.
 * Every insert or erase of a key takes those same two locks first.
 * Writers of one key are therefore serialised and can't insert it twice, while writers of different keys proceed in parallel.
 * Deeper lines are locked afterwards, in (level, address) order.
 *
 * @par Resizing
 * Writers hold the writer mutex shared, whereas grow() takes it exclusively.
 * grow() then copies every key into a table twice as large and publishes it with a release store.
 * Readers may still hold the old table, so it is freed only after a grace period in the style of the Linux kernel's SRCU.
 * Each reader registers in one of two per-slot generation counters, and grow() flips the generation and then waits for the old counters to drain.
 * At no point does a reader wait on a writer.
 *
 * The keys kEmpty (~0) and kTomb (~0 - 1) are reserved, which is safe because neither is a valid canonical k-mer.
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
    static constexpr unsigned kSlots = 7;               // keys per line
    static constexpr std::size_t kDefaultExpected = 1024;

    explicit DeltaOA(std::size_t expected = kDefaultExpected) {
        tab_.store(new Table(lines_for(expected)), std::memory_order_release);
    }
    ~DeltaOA() { delete tab_.load(std::memory_order_acquire); }
    DeltaOA(const DeltaOA&) = delete;
    DeltaOA& operator=(const DeltaOA&) = delete;

    std::size_t size() const noexcept { return size_.load(std::memory_order_relaxed); }
    bool empty() const noexcept { return size() == 0; }

    std::size_t bytes() const noexcept {
        const Table* table = tab_.load(std::memory_order_acquire);

        return table->nlines_total * sizeof(Line);
    }

    bool contains(uint64_t key) const noexcept {
        ReadGuard guard(*this);

        return find_in(guard.table, key);
    }

    /// Returns true if @p key was absent.
    bool insert(uint64_t key) {
        for (;;) {
            {
                std::shared_lock<std::shared_mutex> writer_lock(writer_mutex_);
                Table* table = tab_.load(std::memory_order_acquire);
                const InsertResult result = insert_in(table, key);

                if (result != InsertResult::kTableFull) {
                    if (result == InsertResult::kInserted) size_.fetch_add(1, std::memory_order_relaxed);

                    return result == InsertResult::kInserted;
                }
            }
            // grow() needs the writer mutex to itself. It runs once the shared lock is gone, and the insert retries on the new table.
            grow();
        }
    }

    /// Leaves a tombstone and returns true if @p key was present.
    bool erase(uint64_t key) {
        std::shared_lock<std::shared_mutex> writer_lock(writer_mutex_);
        Table* table = tab_.load(std::memory_order_acquire);
        const std::size_t home_first = table->line_of(0, key, 0), home_second = table->line_of(0, key, 1);
        // Every insert of this key takes these two level-0 locks first, which rules out a race with the erase.
        LockPair level0_locks(table->lines[std::min(home_first, home_second)], table->lines[std::max(home_first, home_second)], home_first == home_second);

        for (unsigned lvl = 0; lvl < table->levels; ++lvl) {
            const std::size_t first_idx = table->line_of(lvl, key, 0), second_idx = table->line_of(lvl, key, 1);

            for (std::size_t line_idx : {first_idx, second_idx}) {
                Line& line = table->lines[line_idx];
                // With both level-0 locks already held by erase(), deeper lines are locked one at a time as the search descends.
                const bool already_locked = (lvl == 0);

                if (!already_locked) lock(line);

                for (unsigned s = 0; s < kSlots; ++s) {
                    if (line.keys[s].load(std::memory_order_relaxed) == key) {
                        // A tombstone lets place() reuse the slot without bumping count again.
                        line.keys[s].store(kTomb, std::memory_order_relaxed);

                        if (!already_locked) unlock(line);

                        size_.fetch_sub(1, std::memory_order_relaxed);
                        return true;
                    }
                }

                if (!already_locked) unlock(line);
            }

            // The key can sit deeper only if some insert found this pair full and set an overflow flag on the way down.
            if (!(table->lines[first_idx].flags.load(std::memory_order_relaxed) & kOverflowFlag) &&
                !(table->lines[second_idx].flags.load(std::memory_order_relaxed) & kOverflowFlag)) break;
        }

        return false;
    }

    /// Not safe while writers are active.
    template <typename F> void for_each(F&& func) const {
        const Table* table = tab_.load(std::memory_order_acquire);

        for (std::size_t i = 0; i < table->nlines_total; ++i)
            for (unsigned s = 0; s < kSlots; ++s) {
                const uint64_t key = table->lines[i].keys[s].load(std::memory_order_relaxed);

                if (key != kEmpty && key != kTomb) func(key);
            }
    }

private:
    /// Set on both lines of a pair once an insert found them full and went a level deeper.
    static constexpr uint8_t kOverflowFlag = 1;
    /// An odd seqlock version means a writer holds the line.
    static constexpr uint32_t kLockedBit = 1;
    /// Each key gets two candidate lines per level.
    static constexpr unsigned kChoicesPerLevel = 2;
    static constexpr std::size_t kMinLevel0Lines = 16;
    static constexpr std::size_t kGrowthFactor = 2;
    /// Sizing aims for half-full lines, which leaves the two-choice rule room to balance.
    static constexpr unsigned kTargetKeysPerLine = kSlots / 2;
    static constexpr unsigned kHashBits = 64;

    /// A Line and a ReaderSlot each fill one x86-64 cache line, so no two threads share one.
    static constexpr std::size_t kCacheLineBytes = 64;

    enum class InsertResult { kPresent, kInserted, kTableFull };

    struct alignas(kCacheLineBytes) Line {
        std::atomic<uint32_t> version{0};           // seqlock version, odd while locked
        std::atomic<uint8_t> count{0};              // slots ever used, tombstones included
        std::atomic<uint8_t> flags{0};
        uint16_t pad = 0;
        std::atomic<uint64_t> keys[kSlots];

        Line() { for (auto& slot_key : keys) slot_key.store(kEmpty, std::memory_order_relaxed); }
    };
    static_assert(sizeof(Line) == kCacheLineBytes, "a line must be one cache line");

    struct Table {
        std::size_t level0_size = 0;                          // in lines, always a power of two
        unsigned levels = 0;
        std::size_t nlines_total = 0;
        std::vector<std::size_t> base;               // first line of each level
        std::unique_ptr<Line[]> lines;

        explicit Table(std::size_t level0_lines) : level0_size(level0_lines) {
            std::size_t level_lines = level0_size, total = 0;

            // Level sizes halve down to a single line, and all deeper levels together hold at most as many lines as level 0.
            while (level_lines >= 1) { base.push_back(total); total += level_lines; ++levels; if (level_lines == 1) break; level_lines >>= 1; }

            nlines_total = total;
            lines.reset(new Line[total]);
        }
        std::size_t line_of(unsigned lvl, uint64_t key, unsigned choice) const noexcept {
            const std::size_t level_lines = level0_size >> lvl;
            // Each (level, choice) pair XORs in a different odd multiple of the golden ratio, which acts as an independent seed.
            const uint64_t hash = mix64(key ^ (kGoldenRatio64 * (kChoicesPerLevel * lvl + choice + 1)), 0);

            return base[lvl] + static_cast<std::size_t>((static_cast<__uint128_t>(hash) * level_lines) >> kHashBits);
        }
    };

    static void lock(Line& line) noexcept {
        for (;;) {
            uint32_t version = line.version.load(std::memory_order_relaxed);

            // An odd version both marks the lock and tells seqlock readers to retry.
            if (!(version & kLockedBit) && line.version.compare_exchange_weak(version, version + 1, std::memory_order_acquire))
                return;

            std::this_thread::yield();
        }
    }
    static void unlock(Line& line) noexcept {
        // The version ends even but changed, which sends any overlapping reader back to retry.
        line.version.fetch_add(1, std::memory_order_release);
    }

    /// Scoped lock over two lines given in address order.
    /// @p same covers the case where both choices hash to one line.
    struct LockPair {
        Line& first; Line& second; bool same;

        LockPair(Line& first_line, Line& second_line, bool same_line) : first(first_line), second(second_line), same(same_line) { lock(first); if (!same) lock(second); }
        ~LockPair() { if (!same) unlock(second); unlock(first); }
    };

    /// Seqlock read that also reports, through @p overflow, whether the line has kOverflowFlag set.
    static bool line_has(const Line& line, uint64_t key, bool& overflow) noexcept {
        for (;;) {
            const uint32_t version_before = line.version.load(std::memory_order_acquire);

            // An odd version means a writer is mid-update, and reading now could return a half-written key.
            if (version_before & kLockedBit) { std::this_thread::yield(); continue; }

            bool hit = false;

            // With the line in cache, comparing all seven slots beats an early exit.
            for (unsigned s = 0; s < kSlots; ++s)
                hit |= line.keys[s].load(std::memory_order_relaxed) == key;

            const bool overflowed = line.flags.load(std::memory_order_relaxed) & kOverflowFlag;

            // The fence keeps the key loads ahead of the second version load.
            std::atomic_thread_fence(std::memory_order_acquire);

            // If the version is unchanged no writer touched the line during the read, otherwise the lookup starts over.
            if (line.version.load(std::memory_order_relaxed) == version_before) { overflow = overflowed; return hit; }
        }
    }

    static bool find_in(const Table* table, uint64_t key) noexcept {
        for (unsigned lvl = 0; lvl < table->levels; ++lvl) {
            bool overflow0 = false, overflow1 = false;

            if (line_has(table->lines[table->line_of(lvl, key, 0)], key, overflow0)) return true;

            if (line_has(table->lines[table->line_of(lvl, key, 1)], key, overflow1)) return true;

            // No flag on either line means no insert ever went past this level.
            if (!overflow0 && !overflow1) return false;
        }

        return false;
    }

    InsertResult insert_in(Table* table, uint64_t key) {
        const std::size_t home_first = table->line_of(0, key, 0), home_second = table->line_of(0, key, 1);
        LockPair level0_locks(table->lines[std::min(home_first, home_second)], table->lines[std::max(home_first, home_second)], home_first == home_second);

        // Checking every level for the key is safe here, because the level-0 locks serialise all writers of this key.
        {
            for (unsigned lvl = 0; lvl < table->levels; ++lvl) {
                const std::size_t first_idx = table->line_of(lvl, key, 0), second_idx = table->line_of(lvl, key, 1);

                for (std::size_t line_idx : {first_idx, second_idx})
                    for (unsigned s = 0; s < kSlots; ++s)
                        if (table->lines[line_idx].keys[s].load(std::memory_order_relaxed) == key) return InsertResult::kPresent;

                if (!(table->lines[first_idx].flags.load(std::memory_order_relaxed) & kOverflowFlag) &&
                    !(table->lines[second_idx].flags.load(std::memory_order_relaxed) & kOverflowFlag)) break;
            }
        }
        for (unsigned lvl = 0; lvl < table->levels; ++lvl) {
            const std::size_t first_idx = table->line_of(lvl, key, 0), second_idx = table->line_of(lvl, key, 1);
            Line& first_line = table->lines[first_idx];
            Line& second_line = table->lines[second_idx];
            const bool already_locked = (lvl == 0);

            // Every writer locks deeper lines in address order, which rules out deadlock.
            if (!already_locked) { lock(first_idx < second_idx ? first_line : second_line); if (first_idx != second_idx) lock(first_idx < second_idx ? second_line : first_line); }

            // Reuse the first free or tombstoned slot.
            auto place = [&](Line& line) -> bool {
                for (unsigned s = 0; s < kSlots; ++s) {
                    const uint64_t slot_key = line.keys[s].load(std::memory_order_relaxed);

                    if (slot_key == kEmpty || slot_key == kTomb) {
                        line.keys[s].store(key, std::memory_order_relaxed);

                        if (slot_key == kEmpty) line.count.fetch_add(1, std::memory_order_relaxed);

                        return true;
                    }
                }

                return false;
            };
            // Placing the key in the emptier of its two lines evens out the fill, and deeper levels are rarely needed as a result.
            const bool x_first = first_line.count.load(std::memory_order_relaxed) <= second_line.count.load(std::memory_order_relaxed);
            bool placed = x_first ? (place(first_line) || place(second_line)) : (place(second_line) || place(first_line));

            // With both lines full, the overflow flags set here tell later lookups and erases to follow the key down a level.
            if (!placed) {
                first_line.flags.fetch_or(kOverflowFlag, std::memory_order_relaxed);
                second_line.flags.fetch_or(kOverflowFlag, std::memory_order_relaxed);
            }

            if (!already_locked) { if (first_idx != second_idx) unlock(first_idx < second_idx ? second_line : first_line); unlock(first_idx < second_idx ? first_line : second_line); }

            if (placed) return InsertResult::kInserted;
        }

        // kTableFull tells insert() that every level on this key's path was full and grow() has to run.
        return InsertResult::kTableFull;
    }

    // Readers register here for the resize grace period.
    // Each thread picks one of kReaderSlots counter pairs by arrival order, and each pair has one counter per generation parity.
    static constexpr unsigned kReaderSlots = 64;
    static constexpr unsigned kGenerationParities = 2;

    struct alignas(kCacheLineBytes) ReaderSlot { std::atomic<uint32_t> counts[kGenerationParities]; ReaderSlot() { counts[0] = 0; counts[1] = 0; } };
    static unsigned reader_slot() noexcept {
        static std::atomic<unsigned> next{0};
        // Reader slots are handed out round-robin, which spreads concurrent readers across separate counter cache lines.
        thread_local unsigned slot = next.fetch_add(1, std::memory_order_relaxed) % kReaderSlots;

        return slot;
    }

    /// Registers the reader in the current generation, which pins the table it loads until the guard is destroyed.
    struct ReadGuard {
        const DeltaOA& delta; unsigned slot, gen; const Table* table;

        explicit ReadGuard(const DeltaOA& owner) : delta(owner), slot(reader_slot()) {
            for (;;) {
                gen = delta.generation_.load(std::memory_order_acquire) % kGenerationParities;
                // Register first, then confirm the generation has not moved.
                delta.reader_slots_[slot].counts[gen].fetch_add(1, std::memory_order_seq_cst);

                if ((delta.generation_.load(std::memory_order_seq_cst) % kGenerationParities) == gen) break;

                // grow() may have flipped the generation meanwhile and be waiting on the other parity. Undo and try again.
                delta.reader_slots_[slot].counts[gen].fetch_sub(1, std::memory_order_release);
            }

            // The table pointer is loaded only after registering. grow() waits for registered readers, so it cannot free the table under this one.
            table = delta.tab_.load(std::memory_order_acquire);
        }
        ~ReadGuard() { delta.reader_slots_[slot].counts[gen].fetch_sub(1, std::memory_order_release); }
    };

    /// Doubles the level-0 lines, or more if the keys still don't fit.
    /// The old table is freed only once no reader can still hold it.
    void grow() {
        std::unique_lock<std::shared_mutex> writer_lock(writer_mutex_);
        Table* old = tab_.load(std::memory_order_acquire);
        // Another writer may have grown the table already.
        // Growing twice only costs memory, though, and the caller retries its insert either way.
        std::size_t level0_size = old->level0_size * kGrowthFactor;

        for (;;) {
            Table* bigger = new Table(level0_size);
            bool all_fit = true;

            // Only live keys are reinserted, which also clears out the tombstones left by erase().
            for (std::size_t i = 0; i < old->nlines_total && all_fit; ++i)
                for (unsigned s = 0; s < kSlots && all_fit; ++s) {
                    const uint64_t key = old->lines[i].keys[s].load(std::memory_order_relaxed);

                    if (key != kEmpty && key != kTomb) all_fit = insert_in(bigger, key) == InsertResult::kInserted;
                }

            if (all_fit) {
                tab_.store(bigger, std::memory_order_release);
                // For the grace period, flip the generation and then wait for the old generation's readers to leave.
                const uint32_t old_gen = generation_.fetch_add(1, std::memory_order_seq_cst) % kGenerationParities;

                for (unsigned i = 0; i < kReaderSlots; ++i)
                    while (reader_slots_[i].counts[old_gen].load(std::memory_order_acquire) != 0) std::this_thread::yield();

                delete old;
                return;
            }

            // A key overflowed every level even at this size, so the loop doubles the capacity again.
            delete bigger;
            level0_size *= kGrowthFactor;
        }
    }

    /// Picks a power-of-two level-0 size, at least kMinLevel0Lines, that holds @p expected keys at kTargetKeysPerLine.
    static std::size_t lines_for(std::size_t expected) {
        std::size_t need = expected / kTargetKeysPerLine + 1, level0_size = kMinLevel0Lines;

        while (level0_size < need) level0_size <<= 1;

        return level0_size;
    }

    std::atomic<Table*> tab_{nullptr};
    std::atomic<std::size_t> size_{0};
    std::shared_mutex writer_mutex_;
    mutable std::atomic<uint32_t> generation_{0};
    mutable ReaderSlot reader_slots_[kReaderSlots];
};

}
