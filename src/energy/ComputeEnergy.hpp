#pragma once
#include "../loop_tree/LoopNode.hpp"
#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "PseudoknotFunctions.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
class ComputeEnergy {
   public:
    ComputeEnergy(std::shared_ptr<LoopNode> root_node, const std::string& sequence,
                  ProcessedRNAEntry processed_rna, int dangle = 2, bool round = false, bool verbose = false)
        : vienna(dangle), pseudo(dangle), root_node_{root_node}, sequence_{sequence}, processed_rna_{processed_rna}, round_{round} {
        process_tree(*root_node_, verbose);
    };

    // Add methods to compute energy, etc.
    float getEnergy() const { return energy_; };

   private:
    ViennaFunctions vienna;
    PseudoknotFunctions pseudo;
    std::shared_ptr<LoopNode> root_node_;
    const std::string& sequence_;
    ProcessedRNAEntry processed_rna_;
    float energy_ = 0.0f;
    bool round_ = false;
    void process_tree(LoopNode& root_node, bool verbose = false);
    float process_node(LoopNode& node);
};

}  // namespace knotergy