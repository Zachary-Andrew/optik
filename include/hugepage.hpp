#pragma once

/**
 * @file
 * @brief Transparent huge pages for large, randomly accessed arrays.
 *
 * Slotin's Algorithms for Modern Hardware explains the background in its "Memory paging" chapter.
 * Once an array on 4 KiB pages outgrows the L2 TLB's reach (roughly 6-12 MiB), most random accesses pay for a page walk.
 * With 2 MiB pages, by contrast, a few hundred TLB entries cover the same array.
 *
 * HugeAlloc rounds big allocations up to whole huge pages, aligns them to 2 MiB and advises MADV_HUGEPAGE.
 * The advice matters under THP "madvise" mode, which many distributions ship by default.
 * Under "always", the alignment alone does the job.
 * For A/B timing, OPTIK_THP=0 flips the advice to MADV_NOHUGEPAGE.
 */

#include <cstddef>
#include <cstdlib>
#include <new>
#include <vector>

#if defined(__linux__)
#include <sys/mman.h>
#endif

namespace tpoptoa {

constexpr std::size_t kHugePageBytes = std::size_t(1) << 21;

inline bool thp_enabled() noexcept {
    static const bool enabled = [] {
        const char* env = std::getenv("OPTIK_THP");

        return !(env && env[0] == '0');
    }();

    return enabled;
}

/// For memory we didn't allocate ourselves, such as a mapped index file.
inline void advise_huge(void* addr, std::size_t len) noexcept {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    if (!addr || len < kHugePageBytes) return;
    ::madvise(addr, len, thp_enabled() ? MADV_HUGEPAGE : MADV_NOHUGEPAGE);
#else
    (void)addr; (void)len;
#endif
}

/// Uninitialised array pinned to 4 KiB pages, for trivially copyable T.
/// It suits a big buffer that fills in scattered pieces while other memory is being freed, because RSS then tracks only what's written.
/// Under THP "always", on the other hand, touching one byte per 2 MiB region would make the entire buffer resident at once.
template <typename T>
class SmallPageArray {
public:
    static constexpr int kNoFd = -1;   // anonymous mappings take fd -1

    explicit SmallPageArray(std::size_t count) : count_(count) {
        bytes_ = count * sizeof(T);

        if (!bytes_) return;
#if defined(__linux__)
        void* mem = ::mmap(nullptr, bytes_, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, kNoFd, 0);

        if (mem == MAP_FAILED) throw std::bad_alloc();
#if defined(MADV_NOHUGEPAGE)
        // MADV_NOHUGEPAGE stops a sparsely touched array from committing a whole 2 MiB page for each scattered write.
        ::madvise(mem, bytes_, MADV_NOHUGEPAGE);
#endif
        data_ = static_cast<T*>(mem);
#else
        data_ = static_cast<T*>(::operator new(bytes_));
#endif
    }
    ~SmallPageArray() {
        if (!data_) return;
#if defined(__linux__)
        ::munmap(data_, bytes_);
#else
        ::operator delete(data_);
#endif
    }
    SmallPageArray(const SmallPageArray&) = delete;
    SmallPageArray& operator=(const SmallPageArray&) = delete;
    T* data() noexcept { return data_; }
    std::size_t size() const noexcept { return count_; }

private:
    T* data_ = nullptr;
    std::size_t count_ = 0, bytes_ = 0;
};

/// Drop-in std::allocator that puts anything a huge page or larger on aligned, huge-page-advised memory.
template <typename T>
struct HugeAlloc {
    using value_type = T;
    HugeAlloc() noexcept = default;
    template <typename U> HugeAlloc(const HugeAlloc<U>&) noexcept {}

    // Cache-line-aligned blocks and similar need the aligned operator new even when small.
    static constexpr bool kOverAligned = alignof(T) > __STDCPP_DEFAULT_NEW_ALIGNMENT__;

    T* allocate(std::size_t count) {
        const std::size_t bytes = count * sizeof(T);

        // Small blocks gain nothing from huge pages and use the ordinary heap.
        if (bytes < kHugePageBytes) {
            if constexpr (kOverAligned)
                return static_cast<T*>(::operator new(bytes, std::align_val_t(alignof(T))));

            return static_cast<T*>(::operator new(bytes));
        }

        // Rounding to whole huge pages lets the kernel back the block entirely with them.
        const std::size_t len = (bytes + kHugePageBytes - 1) & ~(kHugePageBytes - 1);
        void* mem = std::aligned_alloc(kHugePageBytes, len);

        if (!mem) throw std::bad_alloc();

        advise_huge(mem, len);
        return static_cast<T*>(mem);
    }
    void deallocate(T* ptr, std::size_t count) noexcept {
        // The same size test as allocate() picks the matching release.
        if (count * sizeof(T) >= kHugePageBytes) std::free(ptr);
        else if constexpr (kOverAligned) ::operator delete(ptr, std::align_val_t(alignof(T)));
        else ::operator delete(ptr);
    }
    template <typename U> bool operator==(const HugeAlloc<U>&) const noexcept { return true; }
    template <typename U> bool operator!=(const HugeAlloc<U>&) const noexcept { return false; }
};

template <typename T> using HugeVec = std::vector<T, HugeAlloc<T>>;

}
