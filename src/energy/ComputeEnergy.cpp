#include "ComputeEnergy.hpp"

namespace knotergy {

void ComputeEnergy::process_tree(LoopNode& node, bool verbose) {
    energy_ += process_node(node);
    if (verbose) std::cout << node.energy_breakdown(sequence_.size());
    for (std::shared_ptr<LoopNode> child : node.children) {
        process_tree(*child, verbose);
    }
    
}

float ComputeEnergy::process_node(LoopNode& node) {
    double node_energy = 0.0;
    LoopType type = node.loop_type;

    switch (type) {
        case LoopType::Stack:
            node_energy += vienna.stack_energy(node.begin, node.end, node.children[0]->begin,
                                               node.children[0]->end, sequence_);
            break;
        case LoopType::Hairpin:
            node_energy += vienna.hairpin_energy(node.begin, node.end, sequence_);
            break;
        case LoopType::Internal:
            node_energy += vienna.internal_loop_energy(
                node.begin, node.end, node.children[0]->begin, node.children[0]->end, sequence_);
            break;
        case LoopType::Multibranch:
            node_energy += vienna.multibranch_energy(node, sequence_);
            break;
        case LoopType::Pseudoknot:
            node_energy += pseudo.pseudoknot_energy(node, sequence_, processed_rna_, round_);
            break;
        case LoopType::External:
            node_energy += vienna.external_energy(node.children, sequence_);
            break;
    }
    node.energy = node_energy;
    return static_cast<float>(node_energy) / 100.0f;
}

}  // namespace knotergy