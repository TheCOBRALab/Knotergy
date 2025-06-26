#include "ComputeEnergy.hpp"

namespace compute_energy {

float ComputeEnergy::process_node(std::shared_ptr<LoopNode> node) const {
    float node_energy = 0.0f;
    LoopType type = node->loop_type;

    switch (type) {
        case LoopType::Stack:
            break;
        case LoopType::Hairpin:
            break;
        case LoopType::Internal:
            break;
        case LoopType::Multi:
            break;
        case LoopType::Pseudoknot:
            break;
        case LoopType::External:
            break;
    }
    return node_energy;
}

}  // namespace compute_energy