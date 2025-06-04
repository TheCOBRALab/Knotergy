#pragma once

#include <cstddef>

namespace compute_energy{
    constexpr size_t null_index = static_cast<size_t>(-1);

    struct Region {
    int begin = -1;
    int end = -1;
};
}