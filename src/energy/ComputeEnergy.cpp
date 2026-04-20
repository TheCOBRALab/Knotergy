#include "ComputeEnergy.hpp"

#include "./ModifiedBasesFunctions.hpp"
namespace knotergy {

// A stack is used instead of recursion to avoid stack overflows on deeply nested structures.
void ComputeEnergy::process_tree(LoopNode& root, bool verbose) {
    std::vector<LoopNode*> stack;
    stack.push_back(&root);

    while (!stack.empty()) {
        LoopNode* node = stack.back();
        stack.pop_back();

        energy_ += process_node(*node);
        if (verbose) {
            std::cout << node->energy_breakdown(sequence_.size());
        }

        for (auto it = node->children.rbegin(); it != node->children.rend(); ++it) {
            stack.push_back(it->get());
        }
    }
}

double ComputeEnergy::process_node(LoopNode& node) {
    double node_energy = 0.0;
    bool is_inf = false;
    const bool has_mod = processed_rna_.has_modified_bases();
    const std::vector<std::string_view>& mod_sequence = processed_rna_.get_modified_sequence();

    switch (node.loop_type) {
        case LoopType::Stack:
            if (has_mod) {
                node_energy = ModifiedBasesFunctions::find_mod_stack_energy(
                    node.begin, node.end, node.children[0]->begin, node.children[0]->end, sequence_,
                    mod_sequence, vp_, mp_);
            } else {
                node_energy =
                    ViennaFunctions::stack_energy(node.begin, node.end, node.children[0]->begin,
                                                  node.children[0]->end, sequence_, vp_);
            }
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
            if (has_mod) {
                node_energy = ModifiedBasesFunctions::find_mod_multiloop_energy(
                    node, processed_rna_, mod_sequence, vp_, mp_);
            } else {
                node_energy = ViennaFunctions::multibranch_energy(node, processed_rna_, vp_);
            }
            break;

        case LoopType::Pseudoknot:
            node_energy = PseudoknotFunctions::pseudoknot_energy(node, processed_rna_, vp_, mp_,
                                                                 pkp_, is_inf, round_);
            break;

        case LoopType::External:
            if (has_mod) {
                node_energy = ModifiedBasesFunctions::find_mod_external_energy(
                    node.children, processed_rna_, mod_sequence, vp_, mp_);
            } else {
                node_energy = ViennaFunctions::external_energy(node.children, processed_rna_, vp_);
            }
            break;

        default:
            THROW_ERROR("Unknown loop type encountered during energy computation: " +
                        std::to_string(static_cast<int>(node.loop_type)));
    }

    node.energy = node_energy;
    node.is_inf = is_inf;
    infinite_energy_flag_ = infinite_energy_flag_ || is_inf;
    return node_energy / 100.0;
}

}  // namespace knotergy