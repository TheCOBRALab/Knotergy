/** @file
 *  @brief Implementation of the Dangle1 class.
 *
 * Dangle 1 is computed by pre-computing all dangle possibilities for each child loop node,
 * and then walking the children in 5'->3' order to find the optimal configuration.
 *
 * The optimal walk is computed using dynamic programming, where the state is defined by whether the
 * last stem used its 3' neighbor as a dangle or not. If there is a gap of at least 2 unpaired
 * nucleotides between two stems, the 5' dangle of the next stem is free to use, and the state is
 * reset to allow both possibilities.
 */

#include "Dangle1.hpp"

#include "energy/dangles/DangleUtils.hpp"
#include "energy/vienna/ViennaUtils.hpp"

#include <ViennaRNA/eval/exterior.hpp>
#include <ViennaRNA/eval/multibranch.hpp>
#include <ViennaRNA/params/basic.hpp>

#include <algorithm>

namespace viennarna = thermorna::viennarna;

namespace knotergy {
namespace {

// Lowest energy of the stems walked so far, split by whether the last stem used the
// nucleotide after it as its 3' dangle.
struct WalkState {
    int free;   // the nucleotide after the last stem is still available
    int taken;  // the last stem used it
};

// Minimum gap between two stems to not share a nucleotide.
constexpr std::size_t no_share_gap = 2;

// Adds the next stem to the walk. `gap` is the number of unpaired nucleotides between
// the previous stem and this one.
WalkState add_stem(WalkState state, const DangleSet& e, std::size_t gap) {
    // If the gap is at least 2, the 5' dangle of the next stem is free to use.
    if (gap >= no_share_gap) {
        state = {std::min(state.free, state.taken), INF};
    }
    return {
        // The stem does not use its 3' neighbor.
        std::min({add_or_inf(state.free, e.no_dangle), add_or_inf(state.free, e.left_dangle),
                  add_or_inf(state.taken, e.no_dangle)}),
        // The stem uses its 3' neighbor.
        std::min({add_or_inf(state.free, e.right_dangle), add_or_inf(state.free, e.both_dangle),
                  add_or_inf(state.taken, e.right_dangle)}),
    };
}

// Walks the children of a loop in 5'->3' order, starting from `state`. `first_gap` is
// the number of unpaired nucleotides before the first child.
WalkState walk_children(const std::vector<LoopNode*>& children,
                        const std::vector<DangleSet>& dangle_energies, WalkState state = {0, INF},
                        std::size_t first_gap = no_share_gap) {
    for (std::size_t k = 0; k < children.size(); ++k) {
        const std::size_t gap = k == 0 ? first_gap : children[k]->begin - children[k - 1]->end - 1;
        state = add_stem(state, dangle_energies[k], gap);
    }
    return state;
}

}  // namespace

int Dangle1::get_external_dangle_1(const std::vector<LoopNode*>& children,
                                   const std::vector<DangleSet>& dangle_energies) {
    // Nothing comes before the first stem, so its 5' neighbor is free.
    const WalkState state = walk_children(children, dangle_energies);
    return std::min(state.free, state.taken);
}

int Dangle1::get_external_dangle_1(const std::vector<LoopNode*>& children,
                                   const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    std::vector<DangleSet> dangle_energies = populate_children_dangle_energies(children, pRNA, vp);
    return get_external_dangle_1(children, dangle_energies);
}

int Dangle1::get_multibranch_dangle_1(const LoopNode& node,
                                      const std::vector<DangleSet>& dangle_energies,
                                      const DangleSet& closing) {
    const std::vector<LoopNode*>& children = node.children;
    if (children.empty()) {
        return closing.best();
    }

    // The gap between the closing pair and the first child
    const std::size_t front_gap = children.front()->begin - node.begin - 1;
    // The gap between the last child and the closing pair
    const std::size_t back_gap = node.end - children.back()->end - 1;

    // Assumes the closing pair uses the j-1 nucleotide for its dangle.
    const WalkState leaves = walk_children(children, dangle_energies,
                                           {closing.no_dangle, closing.left_dangle}, front_gap);

    // assumes the closing pair does not use the j-1 nucleotide for its dangle.
    const WalkState uses = walk_children(children, dangle_energies,
                                         {closing.right_dangle, closing.both_dangle}, front_gap);
    const int uses_energy = back_gap == 1 ? uses.free : std::min(uses.free, uses.taken);

    return std::min({leaves.free, leaves.taken, uses_energy});
}

int Dangle1::get_multibranch_dangle_1(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                      vrna_md_param& vp) {
    const bool is_external = false;
    std::vector<DangleSet> dangle_energies =
        populate_children_dangle_energies(node.children, pRNA, vp, is_external);
    DangleSet closing = get_ml_closing_dangle_energy(node, pRNA, vp);

    return get_multibranch_dangle_1(node, dangle_energies, closing);
}

DangleSet Dangle1::get_ml_closing_dangle_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                                vrna_md_param& vp) {
    const std::string& sequence = pRNA.get_sequence();
    std::size_t pi = node.begin;
    std::size_t pj = node.end;

    auto [n5d, n3d] = ViennaUtils::encode_inner_dangles(pi, pj, pRNA, vp.md);
    unsigned int pair_type = ViennaUtils::reverse_pair_type(sequence[pi], sequence[pj], vp.md);

    // closing pair dangles
    DangleSet ml_dangle{
        viennarna::vrna_E_multibranch_stem(pair_type, -1, -1, vp.p),   // No dangle
        viennarna::vrna_E_multibranch_stem(pair_type, -1, n5d, vp.p),  // Left dangle
        viennarna::vrna_E_multibranch_stem(pair_type, n3d, -1, vp.p),  // Right dangle
        viennarna::vrna_E_multibranch_stem(pair_type, n3d, n5d, vp.p)  // Both dangles
    };

    return ml_dangle;
}

DangleSet Dangle1::get_child_dangle_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                           vrna_md_param& vp, bool is_external) {
    const std::string& sequence = pRNA.get_sequence();
    std::size_t ci = node.begin;
    std::size_t cj = node.end;

    auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(ci, cj, pRNA, vp.md);
    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj], vp.md);

    auto vrna_E_stem =
        is_external ? viennarna::vrna_E_exterior_stem : viennarna::vrna_E_multibranch_stem;
    // child stem dangles
    DangleSet ml_dangle{
        vrna_E_stem(pair_type, -1, -1, vp.p),   // No dangle
        vrna_E_stem(pair_type, n5d, -1, vp.p),  // Left dangle
        vrna_E_stem(pair_type, -1, n3d, vp.p),  // Right dangle
        vrna_E_stem(pair_type, n5d, n3d, vp.p)  // Both dangles
    };

    return ml_dangle;
}

// Precompute dangle energies for all children in the loop
std::vector<DangleSet> Dangle1::populate_children_dangle_energies(
    const std::vector<LoopNode*>& children, const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
    bool is_external) {
    std::vector<DangleSet> dangle_energies;
    dangle_energies.reserve(children.size());
    for (const LoopNode* child : children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            // basepair penalty for pseudoknot stems must be added for multibranch loops.
            int b = is_external ? 0 : vp.p->MLintern[1];
            dangle_energies.push_back(DangleSet{b, INF, INF, INF});
            continue;
        }
        dangle_energies.push_back(get_child_dangle_energy(*child, pRNA, vp, is_external));
    }

    return dangle_energies;
}

}  // namespace knotergy
