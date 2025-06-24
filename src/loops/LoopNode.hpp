#pragma once

#include <algorithm>
#include <memory>

#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {
enum class LoopType { Stack, Hairpin, Internal, Multi, External, Pseudoknot };

struct LoopNode {
   public:
    LoopNode(ClosedRegion cr) : begin{cr.begin}, end{cr.end} {}
    LoopNode() : begin{NULL_INDEX}, end{NULL_INDEX} {}
    size_t begin;
    size_t end;
    int num_of_unpaired_bases;
    LoopType loop_type;
    std::shared_ptr<LoopNode> parent;
    std::vector<std::shared_ptr<LoopNode>> children;
};
}  // namespace compute_energy