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

void ComputeEnergy::process_modified_tree(LoopNode& node, bool verbose) {
    energy_ += process_modified_node(node);
    if (verbose) std::cout << node.energy_breakdown(sequence_.size());
    for (std::shared_ptr<LoopNode> child : node.children) {
        process_modified_tree(*child, verbose);
    }
}

float ComputeEnergy::process_node(LoopNode& node) {
    double node_energy;
    LoopType loop_type = node.loop_type;
    bool is_inf = false;

    switch (loop_type) {
        case LoopType::Stack:
            node_energy =
                ViennaFunctions::stack_energy(node.begin, node.end, node.children[0]->begin,
                                              node.children[0]->end, sequence_, vp_);
            break;
        case LoopType::Hairpin:
            node_energy =
                ViennaFunctions::hairpin_energy(node.begin, node.end, sequence_, is_inf, vp_);
            break;
        case LoopType::Internal:
            node_energy =
                ViennaFunctions::internal_loop_energy(node.begin, node.end, node.children[0]->begin,
                                                      node.children[0]->end, sequence_, vp_);
            break;
        case LoopType::Multibranch:
            node_energy = ViennaFunctions::multibranch_energy(node, processed_rna_, vp_);
            break;
        case LoopType::Pseudoknot:
            node_energy = PseudoknotFunctions::pseudoknot_energy(node, processed_rna_, vp_, mp_,
                                                                 pkp_, is_inf, round_);
            break;
        case LoopType::External:
            node_energy = ViennaFunctions::external_energy(node.children, processed_rna_, vp_);
            break;
        default:
            THROW_ERROR("Unknown loop type encountered during energy computation: " +
                        std::to_string(static_cast<int>(loop_type)));
    }

    node.energy = node_energy;
    node.is_inf = is_inf;
    infinite_energy_flag_ = infinite_energy_flag_ || is_inf;
    return static_cast<float>(node_energy) / 100.0f;
}

float ComputeEnergy::process_modified_node(LoopNode& node) {
    double node_energy;
    LoopType loop_type = node.loop_type;
    const std::vector<std::string_view>& mod_sequence = processed_rna_.get_modified_sequence();
    bool is_inf = false;

    switch (loop_type) {
        case LoopType::Stack:
            node_energy = ModifiedBasesFunctions::find_mod_stack_energy(
                node.begin, node.end, node.children[0]->begin, node.children[0]->end, sequence_,
                mod_sequence, vp_, mp_);
            break;
        case LoopType::Hairpin:
            node_energy =
                ViennaFunctions::hairpin_energy(node.begin, node.end, sequence_, is_inf, vp_);
            break;
        case LoopType::Internal:
            node_energy =
                ViennaFunctions::internal_loop_energy(node.begin, node.end, node.children[0]->begin,
                                                      node.children[0]->end, sequence_, vp_);
            break;
        case LoopType::Multibranch:
            node_energy = ModifiedBasesFunctions::find_mod_multiloop_energy(node, processed_rna_,
                                                                            mod_sequence, vp_, mp_);
            break;
        case LoopType::Pseudoknot:
            node_energy = PseudoknotFunctions::pseudoknot_energy(node, processed_rna_, vp_, mp_,
                                                                 pkp_, is_inf, round_);
            break;
        case LoopType::External:
            node_energy = ModifiedBasesFunctions::find_mod_external_energy(
                node.children, processed_rna_, mod_sequence, vp_, mp_);
            break;
        default:
            THROW_ERROR("Unknown loop type encountered during modified energy computation: " +
                        std::to_string(static_cast<int>(loop_type)));
    }

    node.energy = node_energy;
    node.is_inf = is_inf;
    infinite_energy_flag_ = infinite_energy_flag_ || is_inf;
    return static_cast<float>(node_energy) / 100.0f;
}

}  // namespace knotergy