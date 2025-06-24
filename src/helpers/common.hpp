#pragma once

#include <cstddef>
#include <iostream>

namespace compute_energy {

#define THROW_ERROR(msg) throw DetailedException((msg), __FILE__, __LINE__, __func__)

// Max size of size_t
constexpr size_t NULL_INDEX = static_cast<size_t>(-1);

}  // namespace compute_energy