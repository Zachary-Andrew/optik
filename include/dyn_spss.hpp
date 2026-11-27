#pragma once

/**
 * @file
 * @brief An exact k-mer set that accepts edits, built from a static SpssCore plus two edit buffers.
 *
 * The layout follows BufBOSS (Alanko et al. 2021), which buffers edits beside a BOSS graph.
 * Here, though, the static part is our SPSS index and the buffers are DeltaOA hash sets (delta_oa.hpp):
 *
 *   added_    k-mers inserted that the core doesn't have
 *   deleted_  core k-mers that were erased
 *
 * Membership is (in core && !in deleted_) || (!in core && in added_).
 * Both buffers hold canonical k-mers, so a k-mer repeated in the original input can't confuse the bookkeeping.
 *
 * compact() rebuilds the core with the edits folded in and then clears the buffers.
 * Nothing else may touch the object while it runs.
 */

#include <memory>
#include <string>
#include <vector>

#include "dbg_graph.hpp"
#include "delta_oa.hpp"
#include "spss_core.hpp"

namespace tpoptoa {

class DynSpss {
public:
    /// Starting capacity of each edit buffer, which DeltaOA grows as needed.
    static constexpr std::size_t kEditSetInitialKeys = 1024;

    explicit DynSpss(std::unique_ptr<SpssCore> core) : core_(std::move(core)) {}

    const SpssCore& core() const noexcept { return *core_; }
    unsigned kmer_len() const noexcept { return core_->kmer_len(); }
    std::size_t added() const noexcept { return added_->size(); }
    std::size_t deleted() const noexcept { return deleted_->size(); }

    /// Accepts either orientation and returns true if @p kmer wasn't already present.
    bool insert(uint64_t kmer) {
        const uint64_t canonical = canon(kmer);

        if (core_->contains(kmer)) return deleted_->erase(canonical);

        return added_->insert(canonical);
    }

    /// True if @p kmer was present.
    bool erase(uint64_t kmer) {
        const uint64_t canonical = canon(kmer);

        if (core_->contains(kmer)) return deleted_->insert(canonical);

        return added_->erase(canonical);
    }

    /// Rebuilds the core from the current set's maximal unitigs, keeping k and the minimizer length.
    /// If the SPSS build fails, however, nothing changes and the result is false.
    bool compact(unsigned nthreads = 0) {
        const unsigned kmer_len = core_->kmer_len(), minimizer_len = core_->minimizer_len();
        // DbgGraph expects (k-mer, count) pairs, and a count of 1 is all we need here.
        std::vector<std::pair<uint64_t, uint32_t>> kmer_counts;

        kmer_counts.reserve(core_->stats().kmers + added_->size());
        core_->for_each_kmer([&](uint64_t kmer) {
            if (deleted_->empty() || !deleted_->contains(kmer)) kmer_counts.push_back({kmer, 1});
        });
        added_->for_each([&](uint64_t kmer) { kmer_counts.push_back({kmer, 1}); });
        // An added k-mer may already be in the core, so the merged list is deduplicated.
        std::sort(kmer_counts.begin(), kmer_counts.end());
        kmer_counts.erase(std::unique(kmer_counts.begin(), kmer_counts.end(),
                                      [](const auto& lhs, const auto& rhs) { return lhs.first == rhs.first; }),
                          kmer_counts.end());
        std::vector<std::string> seqs;

        if (!kmer_counts.empty()) {
            DbgGraph graph(kmer_len);

            graph.build_core(std::move(kmer_counts));
            graph.unitigs([&](const uint8_t* unitig, uint32_t len, uint64_t) {
                seqs.emplace_back(reinterpret_cast<const char*>(unitig), len);
            }, nthreads);
        }

        auto fresh = std::make_unique<SpssCore>();

        // The new core replaces the old one only once it is fully built, which is why a failed build leaves the index unchanged.
        if (!fresh->build(seqs, kmer_len, minimizer_len, nthreads, &seqs)) return false;

        core_ = std::move(fresh);
        added_ = std::make_unique<DeltaOA>(kEditSetInitialKeys);
        deleted_ = std::make_unique<DeltaOA>(kEditSetInitialKeys);
        return true;
    }

private:
    uint64_t canon(uint64_t kmer) const noexcept {
        const uint64_t rev_comp = SpssCore::revcomp_packed(kmer, core_->kmer_len());

        return kmer < rev_comp ? kmer : rev_comp;
    }

    std::unique_ptr<SpssCore> core_;
    std::unique_ptr<DeltaOA> added_ = std::make_unique<DeltaOA>(kEditSetInitialKeys);
    std::unique_ptr<DeltaOA> deleted_ = std::make_unique<DeltaOA>(kEditSetInitialKeys);
};

}
