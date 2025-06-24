#pragma once

#include <array>

#include "RNAEntry.hpp"

namespace compute_energy {
struct ClosedRegion {
    size_t begin{};
    size_t end{};
    bool pseudoknotted = false;

    ClosedRegion() = default;
    ClosedRegion(size_t b, size_t e) : begin(b), end(e) {}

    bool operator==(const ClosedRegion& rhs) const { return begin == rhs.begin && end == rhs.end; }
    bool operator<(const ClosedRegion& rhs) const {
        if (end != rhs.end) {
            return end < rhs.end;
        }
        return begin < rhs.begin;
    }
};

// Lets you print out the Region (overloading the << operator )
inline std::ostream& operator<<(std::ostream& os, const ClosedRegion& region) {
    os << "ClosedRegion(" << region.begin << ", " << region.end << ")";
    return os;
}
}  // namespace compute_energy