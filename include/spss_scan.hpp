#pragma once

/**
 * @file
 * @brief SpssCore's streaming k-mer lookups.
 *
 * scan() and scan_reference() look up every k-mer of a sequence, extending along the string match whenever possible.
 * Lookups that do need the index run in interleaved lanes, so the MPHF, record and string misses of different k-mers overlap.
 * spss_core.hpp includes this file at the end, so callers never include it directly.
 */

#include "spss_core.hpp"

namespace tpoptoa {

inline void SpssCore::scan(const char* seq, std::size_t len, ScanStats& out) const noexcept {
    scan_many(&seq, &len, 1, out);
}

inline void SpssCore::scan_reference(const char* seq, std::size_t len, ScanStats& out) const noexcept {
    scan_one(seq, len, out);
}

inline void SpssCore::flush_group(GroupBuf& group, ScanStats& out) const noexcept {
    if (group.count == 0) return;

    const unsigned count = group.count;

    group.count = 0;
    out.valid += count;
    out.searches += count;

    if (stats_.buckets == 0) return;

    MinInfo min_infos[kGroup];
    uint64_t buckets[kGroup], entries[kGroup];

    // Each pass over the group resolves one level and prefetches the next, so the misses of all its k-mers are in flight together.
    // The levels are the MPHF pilot, the bucket record, then the string window or the overflow list.
    for (unsigned i = 0; i < count; ++i) {
        min_infos[i] = minimizer(group.kmers[i]);
        mphf_.prefetch(min_infos[i].key);
    }

    for (unsigned i = 0; i < count; ++i) {
        buckets[i] = mphf_(min_infos[i].key);
        prefetch_read(reinterpret_cast<const uint8_t*>(rec_view_) + bits_to_bytes(buckets[i] * rec_bits_));
    }

    for (unsigned i = 0; i < count; ++i) {
        entries[i] = get_bits(rec_view_, buckets[i] * rec_bits_, rec_bits_);

        // A record with its low bit clear points into the string, while a set bit sends the lookup to the overflow list.
        if (!(entries[i] & 1)) {
            const uint64_t mer_pos = entry_mer_pos(entries[i]);

            prefetch_read(str_view_ + (mer_pos >= kmer_len_ ? mer_pos - kmer_len_ : 0) / kBasesPerStrWord);
            prefetch_read(str_view_ + mer_pos / kBasesPerStrWord + 1);
        } else {
            prefetch_read(reinterpret_cast<const uint8_t*>(ovf_view_) + bits_to_bytes(entry_ovf_offset(entries[i])));
        }
    }

    // By now every load has arrived, so the probes run from cache.
    for (unsigned i = 0; i < count; ++i) {
        const bool hit = entry_has_key(entries[i], min_infos[i].key)
                         && probe_entry(group.kmers[i], min_infos[i], entries[i], nullptr) != kNotFound;

        out.positive += hit;
    }
}

inline void SpssCore::prefetch_window(uint64_t mer_pos) const noexcept {
    prefetch_read(str_view_ + (mer_pos >= kmer_len_ ? mer_pos - kmer_len_ : 0) / kBasesPerStrWord);
    prefetch_read(str_view_ + (mer_pos + kmer_len_) / kBasesPerStrWord);
}

inline void SpssCore::lane_reset(Lane& lane, const char* seq, std::size_t len) const noexcept {
    lane.seq = seq; lane.len = len; lane.pos = 0;
    lane.fwd = lane.fwd_mer = lane.rev_mer = 0; lane.valid = 0; lane.midx = lane.hashed = 0; lane.have_min = false;
    lane.cur = kNotFound; lane.cached_key = ~uint64_t(0); lane.cached_ok = false; lane.stage = kStageNone;
}

inline bool SpssCore::lane_run(Lane& lane, ScanStats& out) const noexcept {
    const unsigned window = kmer_len_ - minimizer_len_ + 1;
    const unsigned mshift = kBitsPerBase * (minimizer_len_ - 1);

    while (lane.pos < lane.len) {
        const uint8_t base = base_encode(lane.seq[lane.pos++]);

        // An N restarts the k-mer and ends any run of extensions.
        if (base == kInvalidBase) { lane.valid = 0; lane.fwd = 0; lane.midx = lane.hashed = 0; lane.have_min = false; lane.cur = kNotFound; continue; }

        // Rolling both strands together makes the canonical form a single compare.
        lane.fwd = ((lane.fwd << kBitsPerBase) | base) & kmask_;
        lane.fwd_mer = ((lane.fwd_mer << kBitsPerBase) | base) & mmask_;
        lane.rev_mer = (lane.rev_mer >> kBitsPerBase) | (uint64_t(base ^ kComplementMask) << mshift);
        ++lane.valid;

        // The ring holds the window's canonical m-mers until they are hashed.
        if (lane.valid >= minimizer_len_) { lane.rcn[lane.midx & kRingMask] = lane.fwd_mer < lane.rev_mer ? lane.fwd_mer : lane.rev_mer; ++lane.midx; }

        if (lane.valid < kmer_len_) continue;

        ++out.valid;

        // After a hit, the next k-mer of a sequence covered by the SPSS is usually the neighbouring one in the string.
        // So comparing one shifted k-mer, in the orientation of the last hit, often answers the lookup with no index access at all.
        if (lane.cur != kNotFound) {
            const bool room = lane.cur_rev ? (lane.cur > lane.seq_lo) : (lane.cur + 1 + kmer_len_ <= lane.seq_hi);
            const uint64_t next_off = lane.cur_rev ? lane.cur - 1 : lane.cur + 1;

            if (room && get_mer(next_off, kmer_len_) == (lane.cur_rev ? revcomp_packed(lane.fwd, kmer_len_) : lane.fwd)) {
                lane.cur = next_off; ++out.positive; ++out.extensions;
                continue;
            }
        }

        const uint64_t first = lane.midx - window;

        // The window minimizer is kept incrementally, and it only has to be recomputed when the old minimum slides out of the window or a gap reset the lane.
        if (lane.hashed < first) { lane.hashed = first; lane.have_min = false; }

        if (lane.have_min && lane.best_lo < first) lane.have_min = false;

        if (!lane.have_min) {
            for (; lane.hashed < lane.midx; ++lane.hashed) lane.ring_hashes[lane.hashed & kRingMask] = mhash(lane.rcn[lane.hashed & kRingMask]);

            // best_lo and best_hi are the first and last positions of the minimum hash, because a tied minimizer can sit at two places in the k-mer.
            uint64_t best = ~uint64_t(0), best_lo = first, best_hi = first;

            for (uint64_t t = first; t < lane.midx; ++t) {
                const uint64_t hash = lane.ring_hashes[t & kRingMask];
                const bool is_less = hash < best, is_le = hash <= best;

                best = is_less ? hash : best;
                best_lo = is_less ? t : best_lo;
                best_hi = is_le ? t : best_hi;
            }

            lane.best = best; lane.best_lo = best_lo; lane.best_hi = best_hi;
            lane.have_min = true;
        } else {
            for (; lane.hashed < lane.midx; ++lane.hashed) {
                const uint64_t hash = mhash(lane.rcn[lane.hashed & kRingMask]);

                lane.ring_hashes[lane.hashed & kRingMask] = hash;

                if (hash < lane.best) { lane.best = hash; lane.best_lo = lane.best_hi = lane.hashed; }
                else if (hash == lane.best) lane.best_hi = lane.hashed;
            }
        }

        lane.min_info.key = lane.rcn[lane.best_lo & kRingMask];
        lane.min_info.pos[0] = static_cast<uint8_t>(lane.best_lo - first);
        lane.min_info.pos[1] = static_cast<uint8_t>(lane.best_hi - first);
        lane.min_info.npos = (lane.best_hi != lane.best_lo) ? kTiedMinPositions : 1;
        lane.kmer = lane.fwd;
        ++out.searches;

        // Consecutive k-mers often share a minimizer, so the bucket record from the last lookup is cached per lane.
        if (lane.min_info.key != lane.cached_key) {                // a new bucket, so stage it
            lane.fresh = true;
            mphf_.prefetch(lane.min_info.key);
            lane.stage = kStagePilot;
            return true;
        }

        if (!lane.cached_ok) { lane.cur = kNotFound; continue; }   // minimizer absent

        lane.fresh = false;                                // bucket known, so only the string remains
        lane.entry = lane.cached_entry;

        if (!entry_overflows(lane.entry)) { prefetch_window(entry_mer_pos(lane.entry)); lane.stage = kStageString; }
        else {
            prefetch_read(reinterpret_cast<const uint8_t*>(ovf_view_) + bits_to_bytes(entry_ovf_offset(lane.entry)));
            lane.stage = kStageOverflow;
        }

        return true;
    }

    return false;
}

inline bool SpssCore::lane_stage(Lane& lane, ScanStats& out) const noexcept {
    // Each stage issues one load and returns, so the scan advances other lanes while that load is in flight.
    switch (lane.stage) {
    case kStagePilot: {
        const uint64_t bkt = mphf_(lane.min_info.key);

        lane.entry = bkt;                                  // bucket id until stage 2
        prefetch_read(reinterpret_cast<const uint8_t*>(rec_view_) + bits_to_bytes(bkt * rec_bits_));
        lane.stage = kStageRecord;
        return true;
    }
    case kStageRecord: {
        lane.entry = get_bits(rec_view_, lane.entry * rec_bits_, rec_bits_);

        if (!entry_overflows(lane.entry)) { prefetch_window(entry_mer_pos(lane.entry)); lane.stage = kStageString; }
        else {
            prefetch_read(reinterpret_cast<const uint8_t*>(ovf_view_) + bits_to_bytes(entry_ovf_offset(lane.entry)));
            lane.stage = kStageOverflow;
        }

        return true;
    }
    case kStageOverflow: {                                           // first overflow record
        // Prefetching more records' windows measured no faster, since a probe usually stops early and extra prefetches crowd the fill buffers.
        unsigned count_width;

        ovf_count(entry_ovf_offset(lane.entry), count_width);
        prefetch_window(get_bits(ovf_view_, entry_ovf_offset(lane.entry) + count_width, lo_bits_ + kNearBits) >> kNearBits);
        lane.stage = kStageString;
        return true;
    }
    default: {                                          // stage 4 compares and probes
        lane.stage = kStageNone;

        if (lane.fresh) {
            lane.cached_key = lane.min_info.key;
            lane.cached_entry = lane.entry;
            lane.cached_ok = entry_has_key(lane.entry, lane.min_info.key);

            if (!lane.cached_ok) { lane.cur = kNotFound; return false; }
        }

        bool rev = false;
        const uint64_t offset = probe_entry(lane.kmer, lane.min_info, lane.entry, &rev);

        if (offset == kNotFound) { lane.cur = kNotFound; return false; }

        ++out.positive;
        lane.cur = offset; lane.cur_rev = rev;
        seq_bounds(offset, lane.seq_lo, lane.seq_hi);
        return false;
    }
    }
}

inline void SpssCore::scan_one(const char* seq, std::size_t len, ScanStats& out) const noexcept {
    if (stats_.buckets == 0) {                    // empty index, which counts but never hits
        for (std::size_t i = 0, run_len = 0; i < len; ++i) {
            if (base_encode(seq[i]) == kInvalidBase) { run_len = 0; continue; }

            if (++run_len >= kmer_len_) { ++out.valid; ++out.searches; }
        }

        return;
    }

    const unsigned window = kmer_len_ - minimizer_len_ + 1;           // m-mers per k-mer (<= 32)
    const unsigned mshift = kBitsPerBase * (minimizer_len_ - 1);
    uint64_t rcn[kMaxWindowMmers], ring_hashes[kMaxWindowMmers];   // rings of canonical m-mers and hashes
    uint64_t fwd = 0, fwd_mer = 0, rev_mer = 0;
    std::size_t valid = 0;
    uint64_t midx = 0, hashed = 0;             // m-mers seen and hashed in the segment
    // Incremental window minimum over hashed m-mers at absolute indices, with the leftmost and rightmost positions of the minimum hash.
    uint64_t best = ~uint64_t(0), best_lo = 0, best_hi = 0;
    bool have_min = false;
    uint64_t cur = kNotFound; bool cur_rev = false; uint64_t seq_lo = 0, seq_hi = 0;
    uint64_t cached_key = ~uint64_t(0), cached_entry = 0;
    bool cached_ok = false;

    for (std::size_t i = 0; i < len; ++i) {
        const uint8_t base = base_encode(seq[i]);

        // An N restarts the k-mer and ends any run of extensions.
        if (base == kInvalidBase) { valid = 0; fwd = 0; midx = hashed = 0; have_min = false; cur = kNotFound; continue; }

        fwd = ((fwd << kBitsPerBase) | base) & kmask_;
        fwd_mer = ((fwd_mer << kBitsPerBase) | base) & mmask_;
        rev_mer = (rev_mer >> kBitsPerBase) | (uint64_t(base ^ kComplementMask) << mshift);
        ++valid;

        if (valid >= minimizer_len_) { rcn[midx & kRingMask] = fwd_mer < rev_mer ? fwd_mer : rev_mer; ++midx; }

        if (valid < kmer_len_) continue;

        ++out.valid;

        // Within a stored string the next query k-mer is usually the next stored one. A shifted compare then replaces the MPHF lookup.
        if (cur != kNotFound) {
            const bool room = cur_rev ? (cur > seq_lo) : (cur + 1 + kmer_len_ <= seq_hi);
            const uint64_t next_off = cur_rev ? cur - 1 : cur + 1;

            if (room && get_mer(next_off, kmer_len_) == (cur_rev ? revcomp_packed(fwd, kmer_len_) : fwd)) {
                cur = next_off; ++out.positive; ++out.extensions;
                continue;
            }
        }

        // A full search hashes the window's not-yet-hashed m-mers and updates the window minimum incrementally, rescanning (branch-free) only when the minimum left the window or hashing restarted after a run.
        const uint64_t first = midx - window;

        if (hashed < first) { hashed = first; have_min = false; }

        if (have_min && best_lo < first) have_min = false;

        if (!have_min) {
            for (; hashed < midx; ++hashed) ring_hashes[hashed & kRingMask] = mhash(rcn[hashed & kRingMask]);

            best = ~uint64_t(0); best_lo = best_hi = first;

            for (uint64_t t = first; t < midx; ++t) {
                const uint64_t hash = ring_hashes[t & kRingMask];
                const bool is_less = hash < best, is_le = hash <= best;

                best = is_less ? hash : best;
                best_lo = is_less ? t : best_lo;
                best_hi = is_le ? t : best_hi;
            }

            have_min = true;
        } else {
            for (; hashed < midx; ++hashed) {
                const uint64_t hash = mhash(rcn[hashed & kRingMask]);

                ring_hashes[hashed & kRingMask] = hash;

                if (hash < best) { best = hash; best_lo = best_hi = hashed; }
                else if (hash == best) best_hi = hashed;
            }
        }

        MinInfo min_info;

        // Positions are relative to the window start, as probe() expects.
        min_info.key = rcn[best_lo & kRingMask];
        min_info.pos[0] = static_cast<uint8_t>(best_lo - first);
        min_info.pos[1] = static_cast<uint8_t>(best_hi - first);
        min_info.npos = (best_hi != best_lo) ? kTiedMinPositions : 1;
        ++out.searches;

        // The bucket entry is refetched only when the minimizer changes, because neighbouring k-mers usually share one.
        if (min_info.key != cached_key) {
            cached_key = min_info.key;
            cached_entry = bucket_entry(min_info.key);
            cached_ok = entry_has_key(cached_entry, min_info.key);
        }

        if (!cached_ok) {                                  // minimizer absent
            cur = kNotFound; continue;
        }

        bool rev = false;
        const uint64_t offset = probe_entry(fwd, min_info, cached_entry, &rev);

        if (offset == kNotFound) { cur = kNotFound; continue; }

        ++out.positive;
        cur = offset; cur_rev = rev;
        seq_bounds(offset, seq_lo, seq_hi);
    }
}

}
