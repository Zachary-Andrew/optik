#pragma once

/*
 * dyn_spss.hpp - dynamic exact k-mer set, a static SpssCore plus two buffered edit sets.
 *
 * The source is the static-core-plus-buffered-updates design of BufBOSS (Alanko et al. 2021).
 * BufBOSS buffers edits next to a BOSS de Bruijn graph.
 * Here the static part is optik's SPSS index and the buffers are the DeltaOA hash sets of delta_oa.hpp.
 *
 *   added_    canonical k-mers inserted that the static core lacks
 *   deleted_  canonical k-mers of the static core that were erased
 *
 * A k-mer is present if the core has it and it is not in deleted_, or the core lacks it and it is in added_.
 * Edits are keyed by canonical value, so they stay correct when the static input repeats a k-mer.
 *
 * compact() folds the edits into a fresh static core (k-mers -> maximal unitigs -> SpssCore) and empties both sets.
 * It must not run concurrently with any other operation on the object.
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
    /// Initial expected size of each edit set, which DeltaOA grows past.
    static constexpr std::size_t kEditSetInitialKeys = 1024;

    explicit DynSpss(std::unique_ptr<SpssCore> core) : core_(std::move(core)) {}

    const SpssCore& core() const noexcept { return *core_; }
    unsigned k() const noexcept { return core_->k(); }
    std::size_t added() const noexcept { return added_->size(); }
    std::size_t deleted() const noexcept { return deleted_->size(); }

    /// insert() - add k-mer @x in either orientation, reporting true when it was absent.
    bool insert(uint64_t x) {
        const uint64_t c = canon(x);
        if (core_->contains(x)) return deleted_->erase(c);
        return added_->insert(c);
    }

    /// erase() - remove k-mer @x, true when it had been present.
    bool erase(uint64_t x) {
        const uint64_t c = canon(x);
        if (core_->contains(x)) return deleted_->insert(c);
        return added_->erase(c);
    }

    /// compact() - rebuild the static core from the current set (maximal unitigs of the node-centric graph over the k-mers) and clear the edit sets.
    /// Keeps the core's k and minimizer length.
    /// Returns false, with the object unchanged, if the SPSS build fails.
    bool compact(unsigned nthreads = 0) {
        const unsigned K = core_->k(), M = core_->m();
        // DbgGraph takes counted k-mers, so each canonical k-mer gets a count of 1
        std::vector<std::pair<uint64_t, uint32_t>> kmer_counts;
        kmer_counts.reserve(core_->stats().kmers + added_->size());
        core_->for_each_kmer([&](uint64_t c) {
            if (deleted_->empty() || !deleted_->contains(c)) kmer_counts.push_back({c, 1});
        });
        added_->for_each([&](uint64_t c) { kmer_counts.push_back({c, 1}); });
        std::sort(kmer_counts.begin(), kmer_counts.end());
        kmer_counts.erase(std::unique(kmer_counts.begin(), kmer_counts.end(),
                                      [](const auto& a, const auto& b) { return a.first == b.first; }),
                          kmer_counts.end());
        std::vector<std::string> seqs;
        if (!kmer_counts.empty()) {
            DbgGraph g(K);
            g.build_core(std::move(kmer_counts));
            g.unitigs([&](const uint8_t* s, uint32_t len, uint64_t) {
                seqs.emplace_back(reinterpret_cast<const char*>(s), len);
            }, nthreads);
        }
        auto fresh = std::make_unique<SpssCore>();
        if (!fresh->build(seqs, K, M, nthreads, &seqs)) return false;
        core_ = std::move(fresh);
        added_ = std::make_unique<DeltaOA>(kEditSetInitialKeys);
        deleted_ = std::make_unique<DeltaOA>(kEditSetInitialKeys);
        return true;
    }

private:
    uint64_t canon(uint64_t x) const noexcept {
        const uint64_t r = SpssCore::rc(x, core_->k());
        return x < r ? x : r;
    }

    std::unique_ptr<SpssCore> core_;
    std::unique_ptr<DeltaOA> added_ = std::make_unique<DeltaOA>(kEditSetInitialKeys);
    std::unique_ptr<DeltaOA> deleted_ = std::make_unique<DeltaOA>(kEditSetInitialKeys);
};

}
