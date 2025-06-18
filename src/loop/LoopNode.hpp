#pragma once

#include <cstddef>

#include "../helpers/common.hpp"

namespace compute_energy {
class LoopNode {
    LoopNode(size_t start, size_t end, LoopType loop_type, PseudoNestedType pseudo_nested_type)
        : i(start), j(end), type(loop_type), pseudo_type(pseudo_nested_type) {};

    size_t i;
    size_t j;
    LoopType type;
    PseudoNestedType pseudo_type = PseudoNestedType::none;
};
}  // namespace compute_energy