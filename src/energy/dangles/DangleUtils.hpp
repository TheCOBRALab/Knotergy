#pragma once

#include <ViennaRNA/params/basic.hpp>

namespace knotergy {
[[nodiscard]] inline constexpr int add_or_inf(int a, int b) {
    if (a >= INF || b >= INF) {
        return INF;
    }
    return a + b;
}
}  // namespace knotergy