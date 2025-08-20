#pragma once
#include <cstddef>
#include <ostream>
#include <vector>
#include <algorithm>

#include "../pipeline/shared.hpp"

namespace knotergy {

struct BasePair {
    size_t i;
    size_t j;
    std::vector<BasePair> children;

    BasePair() = default;
    BasePair(size_t left_index, size_t right_index) : i{left_index}, j{right_index} {}
    bool is_stack(BasePair child) const { return i + 1 == child.i && j - 1 == child.j; }
};

// Visual representation of the band:
// (((((((((((((((....[.......)))))))..))))))))....]
// ^             ^            ^               ^
// left_border   left_inner   right_inner     right_border

// non-pseudoknots don't have bands
class Band {
   public:
    Band(size_t lb, size_t li, size_t ri, size_t rb, const std::vector<size_t>& pairings, const std::vector<size_t>& cr_pairings)
        : left_border_{lb}, left_inner_{li}, right_inner_{ri}, right_border_{rb} {
        
        if (std::max({pairings[lb], pairings[li], pairings[ri], pairings[rb]}) == NULL_INDEX){
            THROW_ERROR("One or more indices are not base-pairs");
        }

        if (pairings[lb] != rb || pairings[li] != ri){
            THROW_ERROR("Incorrect pairings in Band");
        }

        // find all base pairs
        base_pairs_.emplace_back(left_border_, pairings[left_border_]);
        for (size_t idx = left_border_ + 1; idx <= left_inner_; ++idx) {
            
            // skip closed regions
            if (cr_pairings[idx] != NULL_INDEX){
                base_pairs_.back().children.emplace_back(idx, cr_pairings[idx]);
                idx = cr_pairings[idx];
                continue;
            }
            
            // add base pair
            size_t paired = pairings[idx];
            if (paired >= right_inner_ && paired <= right_border_) {
                base_pairs_.emplace_back(idx, paired);
            }
        }
    }

    bool contains(size_t idx) const {
        return (idx >= left_border_ && idx <= left_inner_) ||
               (idx >= right_inner_ && idx <= right_border_);
    }
    bool contains(size_t idx, size_t idx2) const { return contains(idx) && contains(idx2); }

    bool nests(size_t idx) const { return (idx > left_inner_ && idx < right_inner_); }

    bool nests(size_t idx, size_t idx2) const { return nests(idx) && nests(idx2); }

    // === Read-only accessors ===
    size_t left_border() const { return left_border_; }
    size_t left_inner() const { return left_inner_; }
    size_t right_inner() const { return right_inner_; }
    size_t right_border() const { return right_border_; }

    const std::vector<BasePair>& base_pairs() const { return base_pairs_; }

   private:
    size_t left_border_;   // i
    size_t left_inner_;    // i`
    size_t right_inner_;   // j`
    size_t right_border_;  // j

    std::vector<BasePair> base_pairs_;
};

// === Operator overload for printing ===
inline std::ostream& operator<<(std::ostream& os, const Band& band) {
    os << "Band(" << band.left_border() << ", " << band.left_inner() << ", " << band.right_inner()
       << ", " << band.right_border() << ")";
    return os;
}

}  // namespace knotergy
