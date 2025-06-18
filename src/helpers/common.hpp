#pragma once

#include <cstddef>
#include <iostream>

namespace compute_energy {

#define THROW_ERROR(msg) throw DetailedException((msg), __FILE__, __LINE__, __func__)

// Max size of size_t
constexpr size_t NULL_INDEX = static_cast<size_t>(-1);

struct Region {
    size_t begin{};
    size_t end{};

    Region() = default;
    Region(size_t b, size_t e) : begin(b), end(e) {}

    bool operator==(const Region& rhs) const { return begin == rhs.begin && end == rhs.end; }
};

// Lets you print out the Region (overloading the << operator )
inline std::ostream& operator<<(std::ostream& os, const Region& region) {
    os << "Region(" << region.begin << ", " << region.end << ")";
    return os;
}

enum class LoopType { stackloop, hairpin, interior, multi, external, pseudo };

/******************************************
//possible location status for the loops
*******************************************/
enum class PseudoNestedType {
    none,
    inBand,
    unBand,
    inMulti

};
}  // namespace compute_energy