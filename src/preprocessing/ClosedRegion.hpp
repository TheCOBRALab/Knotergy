#pragma once

#include "preprocessing/RNAEntry.hpp"

#include <array>
#include <ostream>

namespace knotergy {

/**
 * @brief Half-open interval that delimits a closed region in an RNA structure.
 *
 * A *closed region* is a contiguous range [begin, end] in which every opening bracket
 * has its matching closing bracket *within the same range* (i.e., no pairing escapes
 * the region). Nested regions are allowed.
 *
 * Example (four closed regions labeled 1..4):
 *
 *   ..(...)...([...)]...(.(...).)
 *     ^   ^   ^     ^   ^ ^   ^ ^
 *     1   1   2     2   3 4   4 3
 *
 * Using 0-based indices, the regions are:
 *   cr1(2, 6), cr2(10, 16), cr3(20, 28), cr4(22, 26)
 *
 * Notes:
 * - Indices are 0-based.
 * - The end index is inclusive.
 * - A closed region contains at least one base pair.
 * - A closed region may be a single base pair, e.g., "( )".
 * - Nested closed regions are possible, e.g., "( ( ) )".
 *
 * @see RNAProcessor::compute_closed_regions() for region extraction logic.
 */
struct ClosedRegion {
    size_t begin{};
    size_t end{};

    ClosedRegion() = default;
    ClosedRegion(size_t b, size_t e) : begin{b}, end{e} {}

    [[nodiscard]] bool operator==(const ClosedRegion& rhs) const {
        return begin == rhs.begin && end == rhs.end;
    }

    // Should never be needed, but it's to safeguard against unsorted closed regions breaking the
    // loop factory algorithm.
    [[nodiscard]] bool operator<(const ClosedRegion& rhs) const {
        return std::tie(begin, end) < std::tie(rhs.begin, rhs.end);
    }
};

/**
 * @brief Stream insertion for ClosedRegion.
 *
 * Prints in the form: "ClosedRegion(<begin>, <end>)"
 *
 * Example:
 *   ClosedRegion cr{2, 6};
 *   std::cout << cr << std::endl;  // ClosedRegion(2, 6)
 */
inline std::ostream& operator<<(std::ostream& os, const ClosedRegion& region) {
    os << "ClosedRegion(" << region.begin << ", " << region.end << ")";
    return os;
}

}  // namespace knotergy
