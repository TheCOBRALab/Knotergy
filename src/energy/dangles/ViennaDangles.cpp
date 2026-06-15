#include "ViennaDangles.hpp"

#include "energy/modified_bases/ModStack.hpp"

/**
 * This file implements dangle models 1 and 3 in different styles.
 *
 * D1:
 * Precomputes each child's possible dangle contributions in DangleSet objects,
 * groups neighboring children into chains that can compete for shared dangles,
 * then scores each chain with a small dynamic program. The DP state records
 * whether the previous child already used its right-side dangle.
 *
 * This approach is modular and easy to test, but uses extra vectors and passes
 * compared with a more direct walk.
 *
 * D3:
 * Uses a direct ViennaRNA-style fixed-structure multiloop walk. It iterates
 * around the multiloop stems while carrying two running energies: the committed
 * best energy so far, and a pending coaxial-stacking candidate.
 *
 * This approach is lower-overhead and closer to ViennaRNA's control flow, but
 * is less intuitive and less modular than the d1 implementation.
 *
 *
 * Having separate implementations for d1 and d3 is a good way to show multiple
 * algorithmic approaches to the same problem.
 *
 *
 * NOTE: D1: if there are dangles on both sides, it uses mismatched dangle energies,
 *       D3: if there are dangles on both sides, it uses 5'+ 3' dangle energies.
 *       This difference has a real impact on the final energy and follows ViennaRNA implementation.
 *
 */

namespace knotergy {

// State for if the previous pair in a chain took the right dangle
enum TouchingRight { RightFree = 0, RightTaken = 1 };

// Calculate dangle energies for external loops (dangle type 1) with precomputed dangle energies
// This is used by modified bases to inject custom energies. To inject values,
// first get a DangleSet for each child (check ViennaDangles::populate_children_dangle_energies).
// Then inject modified energies into the DangleSet, then this function will compute the d1 energy
int ViennaDangles::get_external_dangle_1(const std::vector<std::unique_ptr<LoopNode>>& children,
                                         const std::vector<DangleSet>& dangle_energies) {
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);
    return process_chains(dangle_chains, dangle_energies);
}

// Calculate dangle energies for external loops (dangle type 1)
// This is used when you don't have any modified bases and can just compute energies directly
int ViennaDangles::get_external_dangle_1(const std::vector<std::unique_ptr<LoopNode>>& children,
                                         const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    std::vector<DangleSet> dangle_energies = populate_children_dangle_energies(children, pRNA, vp);
    return get_external_dangle_1(children, dangle_energies);
}

// Calculate dangle energies for multibranch loops (dangle type 1)
// This (similar to get_external_dangle_1) is used by modified bases to inject custom energies.
int ViennaDangles::get_multibranch_dangle_1(const LoopNode& node,
                                            std::vector<DangleSet> dangle_energies,
                                            DangleSet closing) {
    const std::vector<std::unique_ptr<LoopNode>>& children = node.children;
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);

    return process_ml_chains(dangle_chains, children, dangle_energies, node, closing);
}

// Calculate dangle energies for multibranch loops (dangle type 1)
// This is used when you don't have any modified bases and can just compute energies directly
int ViennaDangles::get_multibranch_dangle_1(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                            vrna_md_param& vp) {
    bool is_external = false;
    const std::vector<std::unique_ptr<LoopNode>>& children = node.children;
    std::vector<DangleSet> dangle_energies =
        populate_children_dangle_energies(children, pRNA, vp, is_external);
    DangleSet closing = get_ml_closing_dangle_energy(node, pRNA, vp);

    return get_multibranch_dangle_1(node, dangle_energies, closing);
}

DangleSet ViennaDangles::get_ml_closing_dangle_energy(const LoopNode& node,
                                                      const ProcessedRNAEntry& pRNA,
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

DangleSet ViennaDangles::get_child_dangle_energy(const LoopNode& node,
                                                 const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                                 bool is_external) {
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
std::vector<DangleSet> ViennaDangles::populate_children_dangle_energies(
    const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
    vrna_md_param& vp, bool is_external) {
    std::vector<DangleSet> dangle_energies;
    dangle_energies.reserve(children.size());
    for (const auto& child : children) {
        dangle_energies.push_back(get_child_dangle_energy(*child, pRNA, vp, is_external));
    }

    return dangle_energies;
}

// Identify chains of children that share dangles
std::vector<std::vector<size_t>> ViennaDangles::get_dangle_chains(
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
int ViennaDangles::process_chain(const std::vector<size_t>& chain,
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
            cur[RightTaken] = std::min({RFreeDL, RFreeDR, RTakenDR, RFreeDB});
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
int ViennaDangles::process_chains(const std::vector<std::vector<size_t>>& dangle_chains,
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
int ViennaDangles::process_ml_chains(const std::vector<std::vector<size_t>>& dangle_chains,
                                     const std::vector<std::unique_ptr<LoopNode>>& children,
                                     const std::vector<DangleSet>& dangle_energies,
                                     const LoopNode& node, const DangleSet closing) {
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
        return process_chains(dangle_chains, dangle_energies, false, {0, closing.best_left()});
    }

    // Neither end is contiguous, but both can dangle
    const int chain1 =
        process_chains(dangle_chains, dangle_energies, false, {0, closing.best_left()});

    const int chain2 = process_chains(dangle_chains, dangle_energies, true, {0, closing.best()});

    return std::min(chain1, chain2);
}

// ------------------- Dangle 3 -----------------------

int add_or_inf(int a, int b) {
    if (a >= INF || b >= INF) {
        return INF;
    }
    return a + b;
}

// applies to closing pair and first child
int ViennaDangles::compute_initial_ld5_for_d3(const MultiloopStem& stem,
                                              const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    const std::string& sequence = pRNA.get_sequence();
    const std::vector<size_t>& pair_table = pRNA.get_pair_table();

    if (stem.begin == 0) {
        THROW_ERROR(
            "Unexpected multiloop stem at position 0."
            "ML with no closing pair does not make sense.");
    }

    const size_t dangle_pos = stem.begin - 1;
    const int encoding = vrna_nucleotide_encode(sequence[dangle_pos], &vp.md);

    int ld5;
    if (pRNA.has_modified_bases()) {
        ld5 = ModBaseUtils::get_dangle5_mod_energy(stem.begin, stem.end, sequence,
                                                   pRNA.get_modified_sequence(), vp,
                                                   all_mod_params(), false);
    } else {
        ld5 = vp.p->dangle5[stem.type][encoding];
    }

    // begin is end for closing pair
    if (stem.begin < 2) {
        return ld5;
    }

    // begin is end for closing pair
    if (pair_table[stem.begin - 2] == NULL_INDEX) {
        return ld5;
    }

    const size_t closing_5 = stem.begin - 2;
    const size_t closing_3 = pair_table[closing_5];

    if (closing_3 != NULL_INDEX) {
        const unsigned int closing_type = ViennaUtils::get_pair_type(
            sequence[static_cast<size_t>(closing_3)], sequence[closing_5], vp.md);

        const int competing_dangle3 = vp.p->dangle3[closing_type][encoding];

        if (competing_dangle3 < ld5) {
            ld5 = 0;
        }
    }

    return ld5;
}

std::vector<MultiloopStem> ViennaDangles::populate_multiloop_stems(const LoopNode& node,
                                                                   const ProcessedRNAEntry& pRNA,
                                                                   vrna_md_param& vp) {
    std::vector<MultiloopStem> stems;
    stems.reserve(node.children.size() + 1);

    for (const auto& child : node.children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            continue;
        }
        int dangle5 = child->n5d_outer >= 0 ? vp.p->dangle5[child->pair_type][child->n5d_outer] : 0;
        int dangle3 = child->n3d_outer >= 0 ? vp.p->dangle3[child->pair_type][child->n3d_outer] : 0;
        stems.push_back(MultiloopStem{child->begin, child->end, child->end, child->pair_type,
                                      dangle5, dangle3});
    }

    // The closing pair is encountered from its 3' side while walking the multiloop.
    // n3d_inner and n5d_inner are reversed because the closing pair is flipped in orientation
    int dangle5_closing = node.n3d_inner >= 0 ? vp.p->dangle5[node.r_pair_type][node.n3d_inner] : 0;
    int dangle3_closing = node.n5d_inner >= 0 ? vp.p->dangle3[node.r_pair_type][node.n5d_inner] : 0;
    stems.push_back(MultiloopStem{node.end, node.begin, node.begin, node.r_pair_type,
                                  dangle5_closing, dangle3_closing});

    MultiloopStem& closing_stem = stems.back();
    closing_stem.initial_ld5 = compute_initial_ld5_for_d3(closing_stem, pRNA, vp);

    MultiloopStem& start_stem = stems.front();
    start_stem.initial_ld5 = compute_initial_ld5_for_d3(start_stem, pRNA, vp);

    return stems;
}

int walk_multiloop_d3_from_start(const ProcessedRNAEntry& pRNA, size_t start_prev,
                                 const std::vector<MultiloopStem>& stems, vrna_md_param& vp,
                                 const all_mod_params& mp) {
    const std::string& sequence = pRNA.get_sequence();
    const size_t stem_count = stems.size();

    const MultiloopStem& start_stem = stems[start_prev];

    unsigned int prev_type = start_stem.type;

    // This is the index of the last base of the previous stem
    size_t prev_end_idx = start_stem.prev_end;
    size_t current = (start_prev + 1) % stem_count;

    // ld5 is the 5' dangle of the previous stem. (left dangle 5)
    int ld5 = start_stem.initial_ld5;
    int energy = 0;       // True energy of the current walk (including coaxial)
    int cx_energy = INF;  // Energy for coaxial
    int coaxial_ml_base_energy = vp.p->MLintern[1];

    for (size_t step = 0; step < stem_count; ++step) {
        const MultiloopStem& stem = stems[current];
        const MultiloopStem& prev_stem = stems[(current + stem_count - 1) % stem_count];

        const size_t begin = stem.begin;
        const size_t end = stem.end;
        const unsigned int current_type = stem.type;

        int new_cx = INF;  // potential new coaxial energy

        // Unpaired bases between the previous stem and the current stem.
        // This determines the dangle case
        const int unpaired_between =
            (begin > prev_end_idx) ? static_cast<int>(begin - prev_end_idx - 1) : 0;

        const int current_ml = vp.p->MLintern[current_type];  // Base ML energy (no dangles)

        energy = add_or_inf(energy, current_ml);
        cx_energy = add_or_inf(cx_energy, current_ml);  // Coaxial energy of previous stem

        int curr_dang5 = 0;  // 5` dangle of the current child
        int prev_dang3 = 0;  // 3` dangle of the previous child

        if (begin > 0) {
            curr_dang5 =
                stem.dangle5;  // vp.p->dangle5[current_type][vrna_nucleotide_encode(sequence[begin
                               // - 1], &vp.md)];
            curr_dang5 = std::min(curr_dang5, 0);  // don't apply dangle bonus if it's positive
                                                   // (dangles should never be a penalty)
        }

        if (prev_end_idx + 1 < sequence.size()) {
            prev_dang3 = prev_stem.dangle3;
            prev_dang3 = std::min(prev_dang3, 0);  // don't apply dangle bonus if it's positive
                                                   // (dangles should never be a penalty)
        }

        switch (unpaired_between) {
            case 0: {
                // adjacent helices: possible flush coaxial stacking
                int stack_energy;
                if (pRNA.has_modified_bases()) {
                    stack_energy = ModStack::find_mod_stack_energy(
                        prev_stem.end, prev_stem.begin, stem.begin, stem.end, sequence,
                        pRNA.get_modified_sequence(), vp, mp);
                    break;
                } else {
                    stack_energy = vp.p->stack[vp.md.rtype[prev_type]][vp.md.rtype[current_type]];
                }
                new_cx = energy + stack_energy;

                // swaps Base ML energy for Coaxial base energy,
                // and removes the 5' dangle of the previous stem (dangles don't apply in coaxial
                // stacking)
                new_cx += -ld5 - vp.p->MLintern[current_type] - vp.p->MLintern[prev_type] +
                          2 * coaxial_ml_base_energy;

                ld5 = 0;
                energy = std::min(energy, cx_energy);
                break;
            }

            case 1: {
                // one unpaired base between helices: ordinary odd-dangle treatment
                const int dang = std::min(prev_dang3, curr_dang5);

                energy = add_or_inf(energy, dang);
                ld5 = dang - prev_dang3;

                // Since each nucleotide or helix end can participate in only one favorable
                // interaction if the previous stem was coaxially stacked, only the 5' dangle of the
                // current stem can be used https://rna.urmc.rochester.edu/NNDB/turner04/mb/
                if (add_or_inf(cx_energy, curr_dang5) < energy) {
                    energy = add_or_inf(cx_energy, curr_dang5);
                    ld5 = curr_dang5;
                }

                new_cx = INF;
                break;
            }

            default: {
                // many unpaired bases between helices
                energy = add_or_inf(energy, curr_dang5 + prev_dang3);
                energy = std::min(energy, add_or_inf(cx_energy, curr_dang5));

                new_cx = INF;
                ld5 = curr_dang5;
                break;
            }
        }

        prev_type = current_type;
        cx_energy = new_cx;
        prev_end_idx = end;
        current = (current + 1) % stem_count;
    }

    // Match ViennaRNA: don't use cx_energy here, because that would allow
    // the final helix to stack with the first helix in this walk.
    return energy;
}

int ViennaDangles::get_multibranch_dangle_3(const LoopNode& node,
                                            const std::vector<MultiloopStem>& stems,
                                            const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                            const all_mod_params& mp) {
    if (node.children.empty()) {
        return ViennaDangles::get_multibranch_dangle_1(node, pRNA, vp);
    }

    // First walk:
    // start from the multiloop closing pair. This disallows stacking of the
    // last child back into the closing pair at the final edge of this walk.
    int best = walk_multiloop_d3_from_start(pRNA, node.children.size(), stems, vp, mp);

    // Second walk:
    // start from the first child. This disallows stacking of the closing pair
    // back into the first child at the final edge of this walk.
    //
    // ViennaRNA does the same "walk around the loop twice" trick for d3.
    best = std::min(best, walk_multiloop_d3_from_start(pRNA, 0, stems, vp, mp));

    return best;
}

int ViennaDangles::get_multibranch_dangle_3(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                            vrna_md_param& vp, const all_mod_params& mp) {
    std::vector<MultiloopStem> stems = populate_multiloop_stems(node, pRNA, vp);
    return get_multibranch_dangle_3(node, stems, pRNA, vp, mp);
};

}  // namespace knotergy