#pragma once

#include <cstddef>

namespace compute_energy {

#define THROW_ERROR(msg) throw DetailedException((msg), __FILE__, __LINE__, __func__)

// Max size of size_t
constexpr size_t NULL_INDEX = static_cast<size_t>(-1);

struct Region {
    size_t begin = NULL_INDEX;
    size_t end = NULL_INDEX;
};
}  // namespace compute_energy