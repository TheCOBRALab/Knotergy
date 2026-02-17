#pragma once

#include <iostream>
#include <vector>
#include "../preprocessing/ClosedRegion.hpp"

namespace knotergy {

/**
 * @brief Represents a base pair in an RNA secondary structure.
 *
 * A base pair consists of two indices (i, j) where i < j, representing
 * positions in the RNA sequence that are paired together.
 * 
 * These base pairs can also have nested closed regions as children
 * 
 * Example:
 * (.....(.......).......([.....)].......)
 * ^     ^       ^       ^       ^       ^
 * i     nested_cr       nested_cr       j   
 * 
 * 
 * @param i 5' position of the base pair.
 * @param j 3' position of the base pair.
 * @param children Nested base pairs (closed regions) within this pair.
 */
struct BasePair {
    size_t i;  ///< 5' position of the base pair.
    size_t j;  ///< 3' position of the base pair.
    std::vector<ClosedRegion> children;  ///< Nested base pairs (closed regions) within this pair.

    BasePair() = default;

    /**
     * @brief Construct a base pair from two indices.
     *
     * @param left_index 5' position of the base pair.
     * @param right_index 3' position of the base pair.
     */
    BasePair(size_t left_index, size_t right_index) : i{left_index}, j{right_index} {}

    /**
     * @brief Construct a base pair from two indices and its children.
     * 
     * @param left_index 5' position of the base pair.
     * @param right_index 3' position of the base pair.
     * @param child_regions Nested base pairs (closed regions) within this pair.
     */
    BasePair(size_t left_index, size_t right_index, std::vector<ClosedRegion> child_regions)
        : i{left_index}, j{right_index}, children{std::move(child_regions)} {}

    /**
     * @brief Check if this base pair forms a stack with a child base pair.
     *
     * A stack occurs when the child base pair is immediately adjacent:
     * (i, j) stacks with (i+1, j-1).
     *
     * @param child The potential child base pair.
     * @return True if this forms a stack with the child.
     */
    bool is_stack(BasePair child) const { return i + 1 == child.i && j - 1 == child.j; }
};

// === Operator overload for printing ===
inline std::ostream& operator<<(std::ostream& os, const BasePair& bp) {
    os << "(" << bp.i << ", " << bp.j << ")";
    
    // Display children base pairs if any
    if (!bp.children.empty()) {
        os << "\n  Children: ";
        for (const ClosedRegion& child : bp.children) {
            os << child << " ";
        }
    }
    return os;
}
}