#pragma once

/*
 * hugepage.hpp - transparent huge pages for large, randomly accessed arrays.
 *
 * The source is Slotin, Algorithms for Modern Hardware, "Memory paging".
 * An array on 4 KiB pages that is larger than the second-level TLB reach of about 6-12 MiB pays a page walk on most random accesses.
 * On 2 MiB pages a few hundred TLB entries cover it.
 *
 * As applied here, HugeAlloc rounds allocations of at least one huge page up to whole huge pages, aligns them to 2 MiB and marks them MADV_HUGEPAGE, so the kernel backs them with huge pages even in THP "madvise" mode, the default on many distributions.
 * In "always" mode only the alignment matters.
 * Smaller allocations use operator new.
 * Setting OPTIK_THP=0 marks the ranges MADV_NOHUGEPAGE instead, for A/B measurements.
 */

#include <cstddef>
#include <cstdlib>
#include <new>
#include <vector>

#if defined(__linux__)
#include <sys/mman.h>
#endif

namespace tpoptoa {

/// Size of an x86-64 transparent huge page, 2 MiB.
constexpr std::size_t kHugePageBytes = std::size_t(1) << 21;

/// thp_enabled() - false when OPTIK_THP starts with '0', read once.
inline bool thp_enabled() noexcept {
    static const bool on = [] {
        const char* e = std::getenv("OPTIK_THP");
        return !(e && e[0] == '0');
    }();
    return on;
}

/// advise_huge() - mark an existing page-aligned range (an index mmap, for example) for huge pages, or against them when OPTIK_THP=0.
/// Ranges shorter than one huge page are left alone.
inline void advise_huge(void* p, std::size_t len) noexcept {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    if (!p || len < kHugePageBytes) return;
    ::madvise(p, len, thp_enabled() ? MADV_HUGEPAGE : MADV_NOHUGEPAGE);
#else
    (void)p; (void)len;
#endif
}

/// SmallPageArray - uninitialised array kept on 4 KiB pages (MADV_NOHUGEPAGE).
///
/// For a large buffer filled in scattered pieces while other memory is being released, whose resident size then grows only with what has been written.
/// In THP "always" mode one byte written in each 2 MiB region would make the whole buffer resident at once.
/// T must be trivially copyable.
template <typename T>
class SmallPageArray {
public:
    explicit SmallPageArray(std::size_t n) : n_(n) {
        bytes_ = n * sizeof(T);
        if (!bytes_) return;
#if defined(__linux__)
        void* p = ::mmap(nullptr, bytes_, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (p == MAP_FAILED) throw std::bad_alloc();
#if defined(MADV_NOHUGEPAGE)
        ::madvise(p, bytes_, MADV_NOHUGEPAGE);
#endif
        p_ = static_cast<T*>(p);
#else
        p_ = static_cast<T*>(::operator new(bytes_));
#endif
    }
    ~SmallPageArray() {
        if (!p_) return;
#if defined(__linux__)
        ::munmap(p_, bytes_);
#else
        ::operator delete(p_);
#endif
    }
    SmallPageArray(const SmallPageArray&) = delete;
    SmallPageArray& operator=(const SmallPageArray&) = delete;
    T* data() noexcept { return p_; }
    std::size_t size() const noexcept { return n_; }
private:
    T* p_ = nullptr;
    std::size_t n_ = 0, bytes_ = 0;
};

/// HugeAlloc - std::allocator replacement that puts allocations of one huge page or more on 2 MiB-aligned, huge-page-advised memory.
template <typename T>
struct HugeAlloc {
    using value_type = T;
    HugeAlloc() noexcept = default;
    template <typename U> HugeAlloc(const HugeAlloc<U>&) noexcept {}

    // Types aligned past what plain operator new guarantees, such as cache-line blocks, get the aligned operator new for small allocations too.
    static constexpr bool kOverAligned = alignof(T) > __STDCPP_DEFAULT_NEW_ALIGNMENT__;

    T* allocate(std::size_t n) {
        const std::size_t bytes = n * sizeof(T);
        if (bytes < kHugePageBytes) {
            if constexpr (kOverAligned)
                return static_cast<T*>(::operator new(bytes, std::align_val_t(alignof(T))));
            return static_cast<T*>(::operator new(bytes));
        }
        const std::size_t len = (bytes + kHugePageBytes - 1) & ~(kHugePageBytes - 1);
        void* p = std::aligned_alloc(kHugePageBytes, len);
        if (!p) throw std::bad_alloc();
        advise_huge(p, len);
        return static_cast<T*>(p);
    }
    void deallocate(T* p, std::size_t n) noexcept {
        if (n * sizeof(T) >= kHugePageBytes) std::free(p);
        else if constexpr (kOverAligned) ::operator delete(p, std::align_val_t(alignof(T)));
        else ::operator delete(p);
    }
    template <typename U> bool operator==(const HugeAlloc<U>&) const noexcept { return true; }
    template <typename U> bool operator!=(const HugeAlloc<U>&) const noexcept { return false; }
};

template <typename T> using HugeVec = std::vector<T, HugeAlloc<T>>;

}
