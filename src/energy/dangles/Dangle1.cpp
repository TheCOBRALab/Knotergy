#include "Dangle1.hpp"

#include "energy/vienna/ViennaUtils.hpp"

extern "C" {
#include <ViennaRNA/eval/exterior.h>
#include <ViennaRNA/eval/hairpin.h>
#include <ViennaRNA/eval/internal.h>
#include <ViennaRNA/eval/multibranch.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/sequences/alphabet.h>
#include <ViennaRNA/utils/basic.h>
}

namespace knotergy {
int Dangle1::get_external_dangle_1(const std::vector<std::unique_ptr<LoopNode>>& children,
                                   const std::vector<DangleSet>& dangle_energies) {
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);
    return process_chains(dangle_chains, dangle_energies);
}

int Dangle1::get_external_dangle_1(const std::vector<std::unique_ptr<LoopNode>>& children,
                                   const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    std::vector<DangleSet> dangle_energies = populate_children_dangle_energies(children, pRNA, vp);
    return get_external_dangle_1(children, dangle_energies);
}

int Dangle1::get_multibranch_dangle_1(const LoopNode& node, std::vector<DangleSet> dangle_energies,
                                      DangleSet closing) {
    const std::vector<std::unique_ptr<LoopNode>>& children = node.children;
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);

    return process_ml_chains(dangle_chains, children, dangle_energies, node, closing);
}

int Dangle1::get_multibranch_dangle_1(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                      vrna_md_param& vp) {
    bool is_external = false;
    const std::vector<std::unique_ptr<LoopNode>>& children = node.children;
    std::vector<DangleSet> dangle_energies =
        populate_children_dangle_energies(children, pRNA, vp, is_external);
    DangleSet closing = get_ml_closing_dangle_energy(node, pRNA, vp);

    return get_multibranch_dangle_1(node, dangle_energies, closing);
}

DangleSet Dangle1::get_ml_closing_dangle_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                                vrna_md_param& vp) {
    const std::string& sequence = pRNA.get_sequence();
    size_t pi = node.begin;
    size_t pj = node.end;

    auto [n5d, n3d] = ViennaUtils::encode_inner_dangles(pi, pj, pRNA, vp.md);
    unsigned int pair_type = ViennaUtils::reverse_pair_type(sequence[pi], sequence[pj], vp.md);

    // closing pair dangles
    DangleSet ml_dangle{
        vrna_E_multibranch_stem(pair_type, -1, -1, vp.p),   // No dangle
        vrna_E_multibranch_stem(pair_type, -1, n5d, vp.p),  // Left dangle
        vrna_E_multibranch_stem(pair_type, n3d, -1, vp.p),  // Right dangle
        vrna_E_multibranch_stem(pair_type, n3d, n5d, vp.p)  // Both dangles
    };

    return ml_dangle;
}

DangleSet Dangle1::get_child_dangle_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                           vrna_md_param& vp, bool is_external) {
    const std::string& sequence = pRNA.get_sequence();
    size_t ci = node.begin;
    size_t cj = node.end;

    auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(ci, cj, pRNA, vp.md);
    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj], vp.md);

    auto vrna_E_stem = is_external ? vrna_E_exterior_stem : vrna_E_multibranch_stem;
    // closing pair dangles
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
    const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
    vrna_md_param& vp, bool is_external) {
    std::vector<DangleSet> dangle_energies;
    dangle_energies.reserve(children.size());
    for (const auto& child : children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            dangle_energies.push_back(DangleSet{0, 0, 0, 0});
            continue;
        }
        dangle_energies.push_back(get_child_dangle_energy(*child, pRNA, vp, is_external));
    }

    return dangle_energies;
}

// Identify chains of children that share dangles
std::vector<std::vector<size_t>> Dangle1::get_dangle_chains(
    const std::vector<std::unique_ptr<LoopNode>>& children) {
    std::vector<std::vector<size_t>> dangle_chains;
    // Stores child indices in each chain. Initialize with first child.
    dangle_chains.reserve(children.size());
    if (!children.empty()) dangle_chains.push_back({0});

    // Iterate through children to identify chains
    for (size_t i = 1; i < children.size(); ++i) {
        const std::unique_ptr<LoopNode>& child = children[i];
        const std::unique_ptr<LoopNode>& prev_child = children[i - 1];
        if (child->loop_type == LoopType::Pseudoknot) {
            continue;
        }

        // Check if current child is contiguous with previous child (shares a dangle or is adjacent)
        if (child->begin - prev_child->end <= 2) {
            dangle_chains.back().push_back(i);
        } else {  // start new chain
            dangle_chains.push_back({i});
        }
    }
    // for (const auto& chain : dangle_chains) {
    //     std::cout << "Dangle chain: ";
    //     for (size_t idx : chain) {
    //         std::cout << idx << " ";
    //     }
    //     std::cout << std::endl;
    // }
    return dangle_chains;
}

// Dynamic programming to compute optimal dangle energies for a single chain of children
int Dangle1::process_chain(const std::vector<size_t>& chain,
                           const std::vector<DangleSet>& dangle_energies,
                           bool disable_last_right_dangle, std::array<int, 2> init,
                           DangleSet closing) {
    std::array<int, 2> prev = init;  // default {0, INF}
    for (size_t idx : chain) {
        const DangleSet& energies = dangle_energies[idx];
        std::array<int, 2> cur = {INF, INF};

        // Check if left or right dangle is possible based on adjacency (no unpaired bases in
        // between)

        const bool disable_right_dangle = (idx == chain.back() && disable_last_right_dangle);

        int RFreeD0 = prev[RightFree] + energies.no_dangle;
        int RFreeDL = prev[RightFree] + energies.left_dangle;
        int RFreeDR = prev[RightFree] + energies.right_dangle;
        int RFreeDB = prev[RightFree] + energies.both_dangle;
        int RTakenD0 = prev[RightTaken] + energies.no_dangle;
        int RTakenDR = prev[RightTaken] + energies.right_dangle;

        // Update touching_right based on previous state and current possibilities
        if (disable_right_dangle) {
            cur[RightFree] = std::min({RFreeD0, RFreeDL, RTakenD0});
        } else {
            cur[RightFree] = std::min({RFreeD0, RFreeDL, RTakenD0});
            cur[RightTaken] = std::min({RFreeDR, RTakenDR, RFreeDB});
        }

        // Sanity check for overflow
        if ((cur[RightFree] > INF) || (cur[RightTaken] > INF)) {
            THROW_ERROR("Dangle energy overflow detected in external loop dangle calculation.");
        }

        prev = cur;
    }
    return std::min(prev[RightFree] + closing.best(), prev[RightTaken] + closing.best_left());
}

// Process multiple chains of children and aggregate their dangle energies
int Dangle1::process_chains(const std::vector<std::vector<size_t>>& dangle_chains,
                            const std::vector<DangleSet>& dangle_energies,
                            bool disable_last_right_dangle, std::array<int, 2> init,
                            DangleSet closing) {
    if (dangle_chains.empty()) {
        return closing.best();
    }

    int total = 0;
    const size_t last = dangle_chains.size() - 1;

    for (size_t i = 0; i < dangle_chains.size(); ++i) {
        total += process_chain(
            dangle_chains[i], dangle_energies, i == last ? disable_last_right_dangle : false,
            i == 0 ? init : std::array<int, 2>{0, INF}, i == last ? closing : DangleSet{});
    }

    return total;
}

// Specialized processing for multibranch loops with closing pair dangles
int Dangle1::process_ml_chains(const std::vector<std::vector<size_t>>& dangle_chains,
                               const std::vector<std::unique_ptr<LoopNode>>& children,
                               const std::vector<DangleSet>& dangle_energies, const LoopNode& node,
                               const DangleSet closing) {
    if (children.empty()) {
        return closing.best();
    }

    const bool front_dangle = children.front()->begin - node.begin <= 2;
    const bool back_dangle = node.end - children.back()->end <= 2;

    if (!front_dangle && !back_dangle) {
        return closing.best() + process_chains(dangle_chains, dangle_energies);
    }

    if (front_dangle && !back_dangle) {
        // ((....)..(....)...) or (.(....)..(....)...)
        // closing pair is touching first child or dangles with first child
        return process_chains(dangle_chains, dangle_energies, false,
                              {closing.best_right(), closing.best()});
    }

    if (!front_dangle && back_dangle) {
        // (...(....)..((....)) or (...(....)..((....).)
        // closing pair is touching or dangles with last child
        return process_chains(dangle_chains, dangle_energies, false, {0, INF}, closing);
    }

    // Both ends are involved:
    // ((....)..(....)) / ((....)..(.....).) / (.(....)..(....)) / (.(....)..(.....).)
    const bool front_contig = contiguous(node.begin, children.front()->begin);
    const bool back_contig = contiguous(node.end, children.back()->end);

    if (front_contig && back_contig) {
        return process_chains(dangle_chains, dangle_energies, false, {closing.no_dangle, INF});
    }

    if (front_contig) {
        const int chain1 =
            process_chains(dangle_chains, dangle_energies, false, {closing.no_dangle, INF});

        const int chain2 =
            process_chains(dangle_chains, dangle_energies, true, {closing.best_right(), INF});

        return std::min(chain1, chain2);
    }

    if (back_contig) {
        return process_chains(dangle_chains, dangle_energies, false,
                              {closing.best_right(), closing.best()});
    }

    // Neither end is contiguous, but both can dangle
    const int closing_does_not_take_back = process_chains(dangle_chains, dangle_energies, false,
                                                          {closing.no_dangle, closing.left_dangle});

    const int closing_takes_back = process_chains(dangle_chains, dangle_energies, true,
                                                  {closing.right_dangle, closing.both_dangle});

    return std::min(closing_does_not_take_back, closing_takes_back);
}

}  // namespace knotergy
