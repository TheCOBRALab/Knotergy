#pragma once
#include "../loop_tree/LoopNode.hpp"
#include "PseudoknotFunctions.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
class ComputeEnergy {
   public:
    ComputeEnergy(std::shared_ptr<LoopNode> root_node, const std::string& sequence, bool round = false)
        : root_node_{root_node}, sequence_{sequence}, round_{round} {
        process_tree(*root_node_);
    };

    // Add methods to compute energy, etc.
    float getEnergy() const { return energy_; };

   private:
    ViennaFunctions vienna;
    PseudoknotFunctions pseudo;
    std::shared_ptr<LoopNode> root_node_;
    const std::string& sequence_;
    float energy_ = 0.0f;
    bool round_ = false;
    void process_tree(const LoopNode& root_node);
    float process_node(const LoopNode& node);
};
}  // namespace knotergy