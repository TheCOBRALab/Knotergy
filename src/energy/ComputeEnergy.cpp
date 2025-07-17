#include "ComputeEnergy.hpp"

namespace knotergy {

void ComputeEnergy::process_tree(const LoopNode& node) {
    energy_ += process_node(node);
    for (std::shared_ptr<LoopNode> child : node.children) {
        process_tree(*child);
    }
}

float ComputeEnergy::process_node(const LoopNode& node) {
    double node_energy = 0.0;
    LoopType type = node.loop_type;

    switch (type) {
        case LoopType::Stack:
            node_energy += vienna.stack_energy(node.begin, node.end, node.children[0]->begin,
                                               node.children[0]->end, sequence_);
            std::cout << "Stack (" << node.begin << ", " << node.end << "): " << node_energy
                      << std::endl;
            break;
        case LoopType::Hairpin:
            node_energy += vienna.hairpin_energy(node.begin, node.end, sequence_);
            std::cout << "Hairpin (" << node.begin << ", " << node.end << "): " << node_energy
                      << std::endl;
            break;
        case LoopType::Internal:
            node_energy += vienna.internal_loop_energy(
                node.begin, node.end, node.children[0]->begin, node.children[0]->end, sequence_);
            std::cout << "Internal (" << node.begin << ", " << node.end << "): " << node_energy
                      << std::endl;
            break;
        case LoopType::Multibranch:
            node_energy += vienna.multibranch_energy(node, sequence_);
            std::cout << "Multibranch (" << node.begin << ", " << node.end << "): " << node_energy
                      << std::endl;
            break;
        case LoopType::Pseudoknot:
            node_energy += pseudo.pseudoknot_energy(node, sequence_);
            break;
        case LoopType::External:
            node_energy += vienna.external_energy(node.children, sequence_);
            std::cout << "External (" << node.begin << ", " << node.end << "): " << node_energy
                      << std::endl;
            break;
    }
    return static_cast<float>(node_energy) / 100.0f;
}

}  // namespace knotergy