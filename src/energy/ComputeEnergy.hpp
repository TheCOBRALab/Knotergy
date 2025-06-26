#pragma once
#include "../loops/LoopNode.hpp"

namespace compute_energy {
class ComputeEnergy {
   public:
    ComputeEnergy(std::shared_ptr<LoopNode> root_node) : root_node_(root_node) {
        energy_ = calculate_energy(root_node_);
    };

    // Add methods to compute energy, etc.
    float getEnergy() const { return energy_; };

   private:
    std::shared_ptr<LoopNode> root_node_;
    float energy_;
    float calculate_energy(std::shared_ptr<LoopNode> node) const;
    float ComputeEnergy::process_node(std::shared_ptr<LoopNode> node) const;
};
}  // namespace compute_energy