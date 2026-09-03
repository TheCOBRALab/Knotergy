#pragma once

#include "loop_tree/LoopTypes.hpp"

#include <vector>

namespace knotergy {

struct PKLoopBreakdown {
    LoopType loop_type = LoopType::Unknown;
    double energy = 0;
};

struct PKEnergyBreakdown {
    bool is_inf = false;

    double init_penalty = 0;
    LoopType parent_loop_type = LoopType::Unknown;
    PseudoNestedType pk_nested_type = PseudoNestedType::None;

    double band_penalty = 0;
    int number_of_bands = 0;

    double unpaired_penalty = 0;
    int unpaired_count = 0;

    double cr_penalty = 0;
    int number_of_outsideband_children = 0;

    double total_loop_energy = 0;
    std::vector<PKLoopBreakdown> loop_breakdowns;  // pair of loop type and penalty

    void reserve(size_t n) { loop_breakdowns.reserve(n); }

    double get_total_energy() const {
        return init_penalty + band_penalty + unpaired_penalty + cr_penalty + total_loop_energy;
    }
};

}  // namespace knotergy