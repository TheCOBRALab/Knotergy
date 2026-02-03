#pragma once
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>

#include "../io/common.hpp"
#include "BasePair.hpp"

namespace knotergy {

// Visual representation of a band:
// (((((..(...)..((((((((((....[.......)))))))..))))))))....]
// ^                       ^            ^               ^
// left_border             left_inner   right_inner     right_border

/**
 * @brief Represents a pseudoknot band in an RNA secondary structure.
 *
 * A band is the structural feature of a pseudoknot where base pairs cross each other.
 * It consists of four key positions:
 * - left_border: Leftmost position of the band
 * - left_inner: Rightmost position of left arm
 * - right_inner: Leftmost position of right arm
 * - right_border: Rightmost position of the band
 *
 * The band contains all base pairs that participate in the pseudoknot crossing.
 * Non-pseudoknotted structures do not have bands.
 */
class Band {
   public:

    /**
     * @brief Construct a Band from boundary positions. 
     * 
     * A simple constructor without validation or base pair extraction.
     *
     * @param lb Left border position.
     * @param li Left inner position.
     * @param ri Right inner position.
     * @param rb Right border position.
     */
    Band(size_t lb, size_t li, size_t ri, size_t rb)
        : left_border_{lb}, left_inner_{li}, right_inner_{ri}, right_border_{rb} {}

    /**
     * @brief Construct a Band from boundary positions and pairing information.
     *
     * Validates the band structure and extracts all base pairs that participate
     * in the band, including nested closed regions.
     *
     * @param lb Left border position.
     * @param li Left inner position.
     * @param ri Right inner position.
     * @param rb Right border position.
     * @param pairings Base-pair index mapping for the structure.
     * @param cr_pairings Closed region pairing indices.
     * @throws DetailedException if band structure is invalid.
     */
    Band(size_t lb, size_t li, size_t ri, size_t rb, const std::vector<size_t>& pairings,
         const std::vector<size_t>& cr_pairings)
        : left_border_{lb}, left_inner_{li}, right_inner_{ri}, right_border_{rb} {
        populate_base_pairs(pairings, cr_pairings);
    }

    void populate_base_pairs(const std::vector<size_t>& pairings, const std::vector<size_t>& cr_pairings) {
        if (std::max({pairings[left_border_], pairings[left_inner_], pairings[right_inner_], pairings[right_border_]}) == NULL_INDEX) {
            THROW_ERROR("One or more indices are not base-pairs");
        }

        if (pairings[left_border_] != right_border_ || pairings[left_inner_] != right_inner_) {
            THROW_ERROR("Incorrect pairings in Band");
        }

        // find all base pairs
        base_pairs_.emplace_back(left_border_, pairings[left_border_]);
        for (size_t idx = left_border_ + 1; idx <= left_inner_; ++idx) {
            // skip closed regions
            if (cr_pairings[idx] != NULL_INDEX) {
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

        // find remaining child base pairs on right side
        size_t current_bp_idx = 0;
        size_t next_swap_idx = NULL_INDEX;

        for (size_t idx = right_border_ - 1; idx > right_inner_; --idx) {
            BasePair& current_bp = base_pairs_[current_bp_idx];
            if ((current_bp_idx + 1 < base_pairs_.size()) && next_swap_idx > right_inner_) {
                next_swap_idx = base_pairs_[current_bp_idx + 1].j;
            }

            if (idx == next_swap_idx) {
                ++current_bp_idx;
                continue;
            }

            // skip closed regions
            if (cr_pairings[idx] != NULL_INDEX) {
                idx = cr_pairings[idx];
                current_bp.children.emplace_back(idx, cr_pairings[idx]);
                if (idx == 0) break;
                continue;
            }
        }
    }

    /**
     * @brief Check if an index is within the band region.
     *
     * @param idx Sequence index to check.
     * @return True if the index is in the left or right arm of the band.
     */
    bool contains(size_t idx) const {
        return (idx >= left_border_ && idx <= left_inner_) ||
               (idx >= right_inner_ && idx <= right_border_);
    }

    /**
     * @brief Check if two indices are both within the band region.
     *
     * @param idx First sequence index.
     * @param idx2 Second sequence index.
     * @return True if both indices are in the band.
     */
    bool contains(size_t idx, size_t idx2) const { return contains(idx) && contains(idx2); }

    /**
     * @brief Check if an index is nested within the band (between inner positions).
     *
     * @param idx Sequence index to check.
     * @return True if the index is between left_inner and right_inner.
     */
    bool nests(size_t idx) const { return (idx > left_inner_ && idx < right_inner_); }

    /**
     * @brief Check if two indices are both nested within the band.
     *
     * @param idx First sequence index.
     * @param idx2 Second sequence index.
     * @return True if both indices are nested within the band.
     */
    bool nests(size_t idx, size_t idx2) const { return nests(idx) && nests(idx2); }

    // === Read-only accessors ===

    /// @return Left border position of the band.
    size_t left_border() const { return left_border_; }

    /// @return Left inner position of the band.
    size_t left_inner() const { return left_inner_; }

    /// @return Right inner position of the band.
    size_t right_inner() const { return right_inner_; }

    /// @return Right border position of the band.
    size_t right_border() const { return right_border_; }

    /// @return Vector of base pairs that participate in this band.
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
