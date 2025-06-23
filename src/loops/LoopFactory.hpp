#pragma once

#include "../rna_regions/RNAEntry.hpp"
#include "LoopNode.hpp"

namespace compute_energy {
    class LoopFactory {
    public:
        LoopFactory(const RNAEntry& entry);
        const std::shared_ptr<LoopNode>& root()  const { return root_; }
        const std::vector<std::shared_ptr<LoopNode>>& nodes() const { return nodes_; }
        void print_tree() const;
        void print_tree(const std::shared_ptr<LoopNode>& node, int depth) const ;
        
    private:
        std::vector<ClosedRegion> closed_regions_;
        const RNAEntry& entry_;
        std::vector<LoopNode> loop_tree;
        std::shared_ptr<LoopNode> root_;
        std::vector<std::shared_ptr<LoopNode>> nodes_;
        LoopType get_loop_type(const LoopNode& node);
        std::shared_ptr<LoopNode> root_node;

        void build_tree();
    };

}