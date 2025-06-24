#pragma once

#include <algorithm>
#include <memory>

#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {
enum class LoopType { Stack, Hairpin, Internal, Multi, External, Pseudoknot };
enum class PseudoNestedType { None, InsideBand, OutsideBand, InsideMultiloop };

struct LoopNode {
   public:
    LoopNode(ClosedRegion cr) : begin{cr.begin}, end{cr.end} {}
    LoopNode() : begin{NULL_INDEX}, end{NULL_INDEX} {}

    size_t begin;
    size_t end;

    LoopType loop_type;
    PseudoNestedType pseudo_type = PseudoNestedType::None;
    int number_of_unpaired_bases = 0;
    int number_of_children_inside_band = 0;
    int number_of_children_outside_band = 0;
    int number_of_unpaired_bases_in_children_outside_band = 0;

    std::shared_ptr<LoopNode> parent;
    std::vector<std::shared_ptr<LoopNode>> children;
};
}  // namespace compute_energy