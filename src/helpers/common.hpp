#pragma once

#include <cstddef>

namespace compute_energy {

#define THROW_ERROR(msg) throw DetailedException((msg), __FILE__, __LINE__, __func__)

// Max size of size_t
constexpr size_t NULL_INDEX = static_cast<size_t>(-1);

struct Pair {
    size_t start = NULL_INDEX;
    size_t end = NULL_INDEX;
    bool pseudo = false;

    Pair() = default;
    Pair(size_t s, size_t e, bool p = false) : start(s), end(e), pseudo(p) {}
};

enum class LoopType{
      stackloop, hairpin,	interior,	multi,	external,	pseudo
};


/******************************************
//possible location status for the loops
*******************************************/
enum class PseudoNestedType{
	none, inBand, unBand, inMulti

};
}  // namespace compute_energy