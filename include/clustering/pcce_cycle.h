/**
 * @file pcce_cycle.h
 * @brief Detection of a cycle in the deterministic pCCE partition iteration.
 */
#ifndef CCEX_PCCE_CYCLE_H
#define CCEX_PCCE_CYCLE_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

// Exact, bounded-memory detection for a deterministic iteration. Checkpoints
// are spaced in doubling windows; no growing history or hash-only match is used.
// Observe completed iterations only, after testing the original convergence rule.
class PcceCycleTracker {
    size_t center_bytes_, assignment_bytes_;
    std::vector<unsigned char> saved_;
    uint64_t window_=1, distance_=0;
    bool initialized_=false;

    void remember(const void* centers, const void* assignments) {
        std::memcpy(saved_.data(),centers,center_bytes_);
        std::memcpy(saved_.data()+center_bytes_,assignments,assignment_bytes_);
    }
public:
    PcceCycleTracker(size_t center_bytes,size_t assignment_bytes)
        : center_bytes_(center_bytes), assignment_bytes_(assignment_bytes),
          saved_(center_bytes+assignment_bytes) {}

    void reset() { window_=1; distance_=0; initialized_=false; }

    // Returns the exact period, or zero if none has been proven. Production
    // max_iter <= INT_MAX bounds distance; uint64 windows avoid doubling overflow.
    uint64_t observe(const void* centers,const void* assignments) {
        if(!initialized_) {
            remember(centers,assignments); initialized_=true; return 0;
        }
        ++distance_;
        if(!std::memcmp(saved_.data(),centers,center_bytes_) &&
           !std::memcmp(saved_.data()+center_bytes_,assignments,assignment_bytes_))
            return distance_;
        if(distance_==window_) {
            remember(centers,assignments); window_*=2; distance_=0;
        }
        return 0;
    }
};

#endif
