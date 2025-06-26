#include "ComputeEnergy.hpp"

#include "ViennaFunctions.hpp"

namespace knotergy {

ViennaFunctions vienna;

void ComputeEnergy::process_tree(const LoopNode& root_node) {
    std::vector<std::shared_ptr<LoopNode>> children = root_node.children;
    for (std::shared_ptr<LoopNode> child : children) {
        energy_ += process_node(*child);
        process_tree(*child);
    }
}

float ComputeEnergy::process_node(const LoopNode& node) const {
    float node_energy = 0.0f;
    LoopType type = node.loop_type;

    switch (type) {
        case LoopType::Stack:
            break;
        case LoopType::Hairpin:
            node_energy += vienna.hairpin_energy(node.begin, node.end, sequence_);
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

}  // namespace knotergy