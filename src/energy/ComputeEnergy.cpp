#include "ComputeEnergy.hpp"

#include "./ModifiedBasesFunctions.hpp"
namespace knotergy {

void ComputeEnergy::process_tree(LoopNode& node, bool verbose) {
    energy_ += process_node(node);
    if (verbose) std::cout << node.energy_breakdown(sequence_.size());
    for (std::shared_ptr<LoopNode> child : node.children) {
        process_tree(*child, verbose);
    }
}

float ComputeEnergy::process_node(LoopNode& node) {
    double node_energy = INF;
    LoopType type = node.loop_type;
    const std::vector<std::string_view>& mod_sequence = processed_rna_.get_modified_sequence();

    switch (type) {
        case LoopType::Stack: {
            // Check for modified base stacking energy
            if (processed_rna_.has_modified_bases()) {
                node_energy = ModifiedBasesFunctions::find_mod_stack_energy(node.begin, node.end, node.children[0]->begin, node.children[0]->end, sequence_, mod_sequence, mod_params_);
                break;
            } 
            node_energy = ViennaFunctions::stack_energy(node.begin, node.end, node.children[0]->begin, node.children[0]->end, sequence_);
            break;
            case LoopType::Hairpin:
                node_energy = ViennaFunctions::hairpin_energy(node.begin, node.end, sequence_);
                break;
            case LoopType::Internal:
                node_energy = ViennaFunctions::internal_loop_energy(
                    node.begin, node.end, node.children[0]->begin, node.children[0]->end,
                    sequence_);
                break;
            case LoopType::Multibranch:
                node_energy = ViennaFunctions::multibranch_energy(node, sequence_);
                break;
            case LoopType::Pseudoknot:
                node_energy = PseudoknotFunctions::pseudoknot_energy(node, sequence_, processed_rna_, round_);
                break;
            case LoopType::External:
                if (processed_rna_.has_modified_bases()) {
                    node_energy = ModifiedBasesFunctions::find_mod_external_energy(node.children, sequence_, mod_sequence, mod_params_);
                    break;
                }
                node_energy = ViennaFunctions::external_energy(node.children, sequence_);
                break;
            default:
                THROW_ERROR("Unknown loop type encountered during energy computation: " + std::to_string(static_cast<int>(type)));
        }
    }
    node.energy = node_energy;
    return static_cast<float>(node_energy) / 100.0f;
}

}  // namespace knotergy