#pragma once
#include "preprocessing/ClosedRegion.hpp"
#include "utils/common.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>

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
struct BandBounds {
    size_t left_border;
    size_t left_inner;
    size_t right_inner;
    size_t right_border;
};

/**
 * @brief Construct a base pair from two indices and its children.
 *
 * @param left_index 5' position of the base pair.
 * @param right_index 3' position of the base pair.
 * @param child_regions Nested base pairs (closed regions) within this pair.
 */
struct BasePair {
    BasePair(size_t left_index, size_t right_index, std::vector<ClosedRegion> child_regions = {})
        : i{left_index}, j{right_index}, children{std::move(child_regions)} {}
    size_t i;
    size_t j;
    std::vector<ClosedRegion> children;

    [[nodiscard]] bool is_stack(const BasePair& child) const {
        return i + 1 == child.i && j - 1 == child.j;
    }
};

struct Band {
    Band(size_t lb, size_t li, size_t ri, size_t rb, std::vector<BasePair> base_pairs,
         int number_of_children = 0)
        : left_border_{lb},
          left_inner_{li},
          right_inner_{ri},
          right_border_{rb},
          base_pairs_{std::move(base_pairs)},
          number_of_children_{number_of_children} {};

    Band(BandBounds bounds, std::vector<BasePair> base_pairs, int number_of_children = 0)
        : Band(bounds.left_border, bounds.left_inner, bounds.right_inner, bounds.right_border,
               std::move(base_pairs), number_of_children) {}

    [[nodiscard]] bool contains(size_t idx) const {
        return (idx >= left_border_ && idx <= left_inner_) ||
               (idx >= right_inner_ && idx <= right_border_);
    }
    [[nodiscard]] bool contains(size_t idx, size_t idx2) const {
        return contains(idx) && contains(idx2);
    }
    [[nodiscard]] bool nests(size_t idx) const { return (idx > left_inner_ && idx < right_inner_); }

    [[nodiscard]] bool nests(size_t idx, size_t idx2) const { return nests(idx) && nests(idx2); }

    [[nodiscard]] size_t left_border() const { return left_border_; }

    [[nodiscard]] size_t left_inner() const { return left_inner_; }

    [[nodiscard]] size_t right_inner() const { return right_inner_; }

    [[nodiscard]] size_t right_border() const { return right_border_; }

    [[nodiscard]] int get_number_of_children() const { return number_of_children_; }

    [[nodiscard]] const std::vector<BasePair>& base_pairs() const { return base_pairs_; }

   private:
    size_t left_border_;
    size_t left_inner_;
    size_t right_inner_;
    size_t right_border_;

    std::vector<BasePair> base_pairs_;
    int number_of_children_ = 0;
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
