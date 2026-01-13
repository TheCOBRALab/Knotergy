#pragma once
#include "../loop_tree/LoopNode.hpp"
#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "../preprocessing/RNAProcessor.hpp"
#include "PseudoknotFunctions.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
class ComputeEnergy {
   public:
    ComputeEnergy(std::shared_ptr<LoopNode> root_node, ProcessedRNAEntry processed_rna, std::vector<modified_base_params> mod_params = {}, bool round = false, bool verbose = false)
        : root_node_{root_node}, processed_rna_{processed_rna}, mod_params_{mod_params}, sequence_{processed_rna.get_sequence()}, round_{round} {
        process_tree(*root_node_, verbose);
    };

    // Add methods to compute energy, etc.
    float getEnergy() const { return energy_; };

   private:
    PseudoknotFunctions pseudo;
    std::shared_ptr<LoopNode> root_node_;
    ProcessedRNAEntry processed_rna_;
    std::vector<modified_base_params> mod_params_;
    const std::string& sequence_;
    float energy_ = 0.0f;
    bool round_ = false;
    void process_tree(LoopNode& root_node, bool verbose = false);
    float process_node(LoopNode& node);
};

}  // namespace knotergy