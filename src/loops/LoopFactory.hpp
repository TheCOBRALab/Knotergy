#pragma once

#include "../rna_regions/RNAEntry.hpp"
#include "LoopNode.hpp"

namespace compute_energy {
class LoopFactory {
   public:
    LoopFactory(const RNAEntry& entry);
    void print_tree() const;
    void print_tree(const std::shared_ptr<LoopNode>& node, size_t depth) const;

   private:
    std::vector<ClosedRegion> closed_regions_;
    const RNAEntry& entry_;
    std::vector<LoopNode> loop_tree;
    std::shared_ptr<LoopNode> root_;
    std::vector<std::shared_ptr<LoopNode>> nodes_;
    std::shared_ptr<LoopNode> root_node_;

    LoopType get_loop_type(const LoopNode& node);
    void PseudoNestedCheck(const LoopNode& node);

    void build_tree();
};

}  // namespace compute_energy