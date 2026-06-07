#pragma once
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>

#include "../io/common.hpp"
#include "BasePair.hpp"

namespace knotergy {
/**
 * @brief Represents a pseudoknot band in an RNA secondary structure.
 *
 * Stores the four key positions that define the band.
 * Stores every base pair that participates in the band
 * each base pair also stores any nested closed regions that are contained within it.
 *
 * Visual representation of a band:
 * (((((..(...)..((((((((((....[.......)))))))..))))))))....]
 * ^                       ^            ^               ^
 * left_border             left_inner   right_inner     right_border
 *
 * The band contains all base pairs that participate in the pseudoknot crossing.
 * Non-pseudoknotted structures do not have bands.
 */
class Band {
   public:
    // ------------- Constructors -------------
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
     * @param pair_table Base-pair index mapping for the structure.
     * @param cr_pair_table Closed region pairing indices.
     * @throws DetailedException if band structure is invalid.
     */
    Band(size_t lb, size_t li, size_t ri, size_t rb, const std::vector<size_t>& pair_table,
         const std::vector<size_t>& cr_pair_table)
        : left_border_{lb}, left_inner_{li}, right_inner_{ri}, right_border_{rb} {
        populate_base_pairs(pair_table, cr_pair_table);
    }

    // ------------------- Public methods -------------------

    /**
     * @brief Finds all base pairs that participate in the band
     *
     * Validates the band structure and extracts all base pairs that participate in the band
     * and each base pair's nested closed region.
     *
     * @param pair_table Base-pair index mapping for the structure.
     * @param cr_pair_table Closed region pairing indices.
     *
     * @throws DetailedException if band structure is invalid or if indices are not base pairs.
     */
    void populate_base_pairs(const std::vector<size_t>& pair_table,
                             const std::vector<size_t>& cr_pair_table) {
        // -------------- Validate band structure -------------
        if (std::max({pair_table[left_border_], pair_table[left_inner_], pair_table[right_inner_],
                      pair_table[right_border_]}) == NULL_INDEX) {
            THROW_ERROR("One or more indices are not base-pairs");
        }

        if (pair_table[left_border_] != right_border_ || pair_table[left_inner_] != right_inner_) {
            THROW_ERROR("Incorrect pair_table in Band");
        }

        // -------------- Extract base pairs -------------
        base_pairs_.reserve(std::min(right_border_ - right_inner_, left_inner_ - left_border_) + 1);

        // find all consecutive base pairs on left side and their nested closed regions
        base_pairs_.emplace_back(left_border_, pair_table[left_border_]);
        for (size_t idx = left_border_ + 1; idx <= left_inner_; ++idx) {
            // skip closed regions and add them as children to the current base pair
            if (cr_pair_table[idx] != NULL_INDEX) {
                base_pairs_.back().children.emplace_back(idx, cr_pair_table[idx]);
                idx = cr_pair_table[idx];
                ++number_of_children;
                continue;
            }

            // add base pair if it's part of the band
            size_t paired = pair_table[idx];
            if (paired >= right_inner_ && paired <= right_border_) {
                base_pairs_.emplace_back(idx, paired);
            }
        }

        // If nested closed regions on the right is impossible, we can skip the rest of the process
        if ((right_border_ - right_inner_) <= 2) {
            return;
        }

        // Find all nested closed regions on the right side of the band
        // And add them as children to the correct base pair
        size_t current_bp_idx = 0;
        size_t next_bp_right_border = base_pairs_.size() > 1 ? base_pairs_[1].j : NULL_INDEX;
        for (size_t idx = right_border_ - 1; idx > right_inner_; --idx) {
            BasePair& current_bp = base_pairs_[current_bp_idx];

            // When we reach the next pair in the band, we swap to the next base pair
            if (idx == next_bp_right_border) {
                ++current_bp_idx;
                // If the next base pair is valid, we update the next_bp_right_border
                if ((base_pairs_.size() > current_bp_idx + 1) &&
                    next_bp_right_border > right_inner_) {
                    next_bp_right_border = base_pairs_[current_bp_idx + 1].j;
                }
                continue;
            }

            // skip closed regions and add them as children to the current base pair
            if (cr_pair_table[idx] != NULL_INDEX) {
                idx = cr_pair_table[idx];
                current_bp.children.emplace_back(idx, cr_pair_table[idx]);
                ++number_of_children;
                if (idx == 0) break;
                continue;
            }
        }
    }

    /**
     * @brief Check if an index is within a band's borders (either in the left or right arm).
     *
     * @param idx Index to check.
     * @return True if the index is in the left or right arm of the band.
     */
    [[nodiscard]] bool contains(size_t idx) const {
        return (idx >= left_border_ && idx <= left_inner_) ||
               (idx >= right_inner_ && idx <= right_border_);
    }

    /**
     * @brief Check if two indices are both within the borders of the band.
     *
     * @param idx First index.
     * @param idx2 Second index.
     * @return True if both indices are in the band.
     */
    [[nodiscard]] bool contains(size_t idx, size_t idx2) const {
        return contains(idx) && contains(idx2);
    }

    /**
     * @brief Check if an index is nested within the band (between inner positions).
     *
     * Nests means not within the borders, but the area between the inner positions of the band
     *
     * @param idx Index to check.
     * @return True if the index is between left_inner and right_inner.
     */
    [[nodiscard]] bool nests(size_t idx) const { return (idx > left_inner_ && idx < right_inner_); }

    /**
     * @brief Check if two indices are both nested within the band.
     *
     * @param idx First index.
     * @param idx2 Second index.
     * @return True if both indices are nested within the band.
     */
    [[nodiscard]] bool nests(size_t idx, size_t idx2) const { return nests(idx) && nests(idx2); }

    // === Read-only accessors ===

    /// @return Left border position of the band.
    [[nodiscard]] size_t left_border() const { return left_border_; }

    /// @return Left inner position of the band.
    [[nodiscard]] size_t left_inner() const { return left_inner_; }

    /// @return Right inner position of the band.
    [[nodiscard]] size_t right_inner() const { return right_inner_; }

    /// @return Right border position of the band.
    [[nodiscard]] size_t right_border() const { return right_border_; }

    /// @return Number of children within the band.
    [[nodiscard]] int get_number_of_children() const { return number_of_children; }

    /// @return Vector of base pairs that participate in this band.
    [[nodiscard]] const std::vector<BasePair>& base_pairs() const { return base_pairs_; }

   private:
    size_t left_border_;   // i
    size_t left_inner_;    // i`
    size_t right_inner_;   // j`
    size_t right_border_;  // j

    std::vector<BasePair> base_pairs_;
    int number_of_children = 0;
};

// === Operator overload for printing ===
inline std::ostream& operator<<(std::ostream& os, const Band& band) {
    os << "Band(" << band.left_border() << ", " << band.left_inner() << ", " << band.right_inner()
       << ", " << band.right_border() << ")";
    for (const auto& base_pair : band.base_pairs()) {
        os << "    BasePair(" << base_pair.i << ", " << base_pair.j << ")\n";
    }
    return os;
}

}  // namespace knotergy
