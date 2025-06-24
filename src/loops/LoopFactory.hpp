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
    std::shared_ptr<LoopNode> root_node_;
    size_t structure_length_;

    LoopType get_loop_type(const LoopNode& node);
    void pseudo_nested_check(std::shared_ptr<LoopNode> node);
    std::vector<ClosedRegion> closed_region_bucket_sort(std::vector<ClosedRegion>& closed_regions,
                                                        size_t structure_length);

    void build_tree();
};

}  // namespace compute_energy