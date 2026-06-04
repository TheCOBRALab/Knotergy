#include "ViennaDangles.hpp"

namespace knotergy {

// State for if the previous pair in a chain took the right dangle
enum TouchingRight { RightFree = 0, RightTaken = 1 };

// Calculate dangle energies for external loops (dangle type 1) with precomputed dangle energies
int ViennaDangles::get_external_dangle_1(const std::vector<std::unique_ptr<LoopNode>>& children,
                                         const std::vector<DangleSet>& dangle_energies) {
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);
    return process_chains(dangle_chains, dangle_energies);
}

// Calculate dangle energies for external loops (dangle type 1)
int ViennaDangles::get_external_dangle_1(const std::vector<std::unique_ptr<LoopNode>>& children,
                                         const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    std::vector<DangleSet> dangle_energies = populate_children_dangle_energies(children, pRNA, vp);
    return get_external_dangle_1(children, dangle_energies);
}

// Calculate dangle energies for multibranch loops (dangle type 1)
int ViennaDangles::get_multibranch_dangle_1(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                            vrna_md_param& vp) {
    bool is_external = false;
    const std::vector<std::unique_ptr<LoopNode>>& children = node.children;
    std::vector<DangleSet> dangle_energies =
        populate_children_dangle_energies(children, pRNA, vp, is_external);
    DangleSet closing = get_ml_closing_dangle_energy(node, pRNA, vp);
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);

    return process_ml_chains(dangle_chains, children, dangle_energies, node, closing);
}

int ViennaDangles::get_multibranch_dangle_1(const LoopNode& node,
                                            std::vector<DangleSet> dangle_energies,
                                            DangleSet closing) {
    const std::vector<std::unique_ptr<LoopNode>>& children = node.children;
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);

    return process_ml_chains(dangle_chains, children, dangle_energies, node, closing);
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

// Precompute dangle energies for all children in the loop
std::vector<DangleSet> ViennaDangles::populate_children_dangle_energies(
    const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
    vrna_md_param& vp, bool is_external) {
    const std::string& sequence = pRNA.get_sequence();
    std::vector<DangleSet> dangle_energies;
    dangle_energies.reserve(children.size());

    auto vrna_E_stem = is_external ? vrna_E_exterior_stem : vrna_E_multibranch_stem;

    for (const auto& child : children) {
        const size_t ci = child->begin;
        const size_t cj = child->end;

        auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(ci, cj, pRNA, vp.md);
        const unsigned int pair_type =
            ViennaUtils::get_pair_type(sequence[ci], sequence[cj], vp.md);

        dangle_energies.push_back(DangleSet{
            vrna_E_stem(pair_type, -1, -1, vp.p),   // No dangle
            vrna_E_stem(pair_type, n5d, -1, vp.p),  // Left dangle
            vrna_E_stem(pair_type, -1, n3d, vp.p),  // Right dangle
            vrna_E_stem(pair_type, n5d, n3d, vp.p)  // Both dangles
        });
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

int add_or_inf(int a, int b) {
    if (a >= INF || b >= INF) {
        return INF;
    }
    return a + b;
}

std::vector<int> build_pair_table_for_loop(const LoopNode& node, const std::string& sequence) {
    std::vector<int> pair_table(sequence.size(), -1);

    pair_table[node.begin] = static_cast<int>(node.end);
    pair_table[node.end] = static_cast<int>(node.begin);

    for (const auto& child : node.children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            continue;
        }

        pair_table[child->begin] = static_cast<int>(child->end);
        pair_table[child->end] = static_cast<int>(child->begin);
    }

    return pair_table;
}

std::vector<MultiloopStem> build_multiloop_stems(const LoopNode& node, const std::string& sequence,
                                                 vrna_md_param& vp) {
    std::vector<MultiloopStem> stems;
    stems.reserve(node.children.size() + 1);

    for (const auto& child : node.children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            continue;
        }
        stems.push_back(MultiloopStem{child->begin, child->end, child->end,
                                      ViennaUtils::get_pair_type(*child, sequence, vp.md)});
    }

    // The closing pair is encountered from its 3' side while walking the multiloop.
    stems.push_back(MultiloopStem{node.end, node.begin, node.begin,
                                  ViennaUtils::reverse_pair_type(node, sequence, vp.md)});

    return stems;
}

int prime_ld5_for_start_stem(const MultiloopStem& stem, const std::string& sequence,
                             const std::vector<size_t>& pairings, vrna_md_param& vp) {
    int ld5 = 0;

    // ViennaRNA equivalent:
    //
    // if (sn[j - 1] == sn[j]) {
    //   ld5 = P->dangle5[type][s1[j - 1]];
    //   if ((p = pt[j - 2]) && (sn[j - 2] == sn[j - 1]))
    //     if (P->dangle3[pair[s[p]][s[j - 2]]][s1[j - 1]] < ld5)
    //       ld5 = 0;
    // }
    //
    // Here everything is single-strand / 0-based.
    if (stem.p == 0) {
        return 0;
    }

    const size_t dangle_pos = stem.p - 1;
    int encoding = vrna_nucleotide_encode(sequence[dangle_pos], &vp.md);  // Function from ViennaRNA
    ld5 = vp.p->dangle5[stem.type][encoding];

    if (stem.p >= 2) {
        const size_t check_pos = stem.p - 2;
        const size_t partner = pairings[check_pos];

        if (partner != NULL_INDEX) {
            const size_t partner_pos = static_cast<size_t>(partner);
            const unsigned int other_type =
                ViennaUtils::get_pair_type(sequence[partner_pos], sequence[check_pos], vp.md);

            const int competing_dangle3 = vp.p->dangle3[other_type][encoding];

            if (competing_dangle3 < ld5) {
                ld5 = 0;
            }
        }
    }

    return ld5;
}

int walk_multiloop_d3_from_start(const std::vector<MultiloopStem>& stems, size_t start_prev,
                                 const std::string& sequence, const std::vector<size_t>& pairings,
                                 vrna_md_param& vp) {
    const size_t stem_count = stems.size();

    const MultiloopStem& start_stem = stems[start_prev];

    unsigned int type = start_stem.type;
    size_t i1 = start_stem.exit_position;
    size_t current = (start_prev + 1) % stem_count;

    int ld5 = prime_ld5_for_start_stem(start_stem, sequence, pairings, vp);
    int energy = 0;
    int cx_energy = INF;

    for (size_t step = 0; step < stem_count; ++step) {
        const MultiloopStem& stem = stems[current];

        const size_t p = stem.p;
        const size_t q = stem.q;
        const unsigned int tt = stem.type;

        int new_cx = INF;

        const int unpaired_between = (p > i1) ? static_cast<int>(p - i1 - 1) : 0;

        const int current_ml = vp.p->MLintern[tt];

        energy = add_or_inf(energy, current_ml);
        cx_energy = add_or_inf(cx_energy, current_ml);

        int dang5 = 0;
        int dang3 = 0;

        if (p > 0) {
            dang5 = vp.p->dangle5[tt][vrna_nucleotide_encode(sequence[p - 1], &vp.md)];
        }

        if (i1 + 1 < sequence.size()) {
            dang3 = vp.p->dangle3[type][vrna_nucleotide_encode(sequence[i1 + 1], &vp.md)];
        }

        switch (unpaired_between) {
            case 0: {
                // adjacent helices: possible flush coaxial stacking
                new_cx = energy + vp.p->stack[vp.md.rtype[type]][vp.md.rtype[tt]];
                new_cx += -ld5 - vp.p->MLintern[tt] - vp.p->MLintern[type] + 2 * vp.p->MLintern[1];

                ld5 = 0;
                energy = std::min(energy, cx_energy);
                break;
            }

            case 1: {
                // one unpaired base between helices: ordinary odd-dangle treatment
                const int dang = std::min(dang3, dang5);

                energy = add_or_inf(energy, dang);
                ld5 = dang - dang3;

                if (add_or_inf(cx_energy, dang5) < energy) {
                    energy = add_or_inf(cx_energy, dang5);
                    ld5 = dang5;
                }

                new_cx = INF;
                break;
            }

            default: {
                // many unpaired bases between helices
                energy = add_or_inf(energy, dang5 + dang3);
                energy = std::min(energy, add_or_inf(cx_energy, dang5));

                new_cx = INF;
                ld5 = dang5;
                break;
            }
        }

        type = tt;
        cx_energy = new_cx;
        i1 = q;
        current = (current + 1) % stem_count;
    }

    // Match ViennaRNA: don't use cx_energy here, because that would allow
    // the final helix to stack with the first helix in this walk.
    return energy;
}

int ViennaDangles::get_multibranch_dangle_3(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                            vrna_md_param& vp) {
    const std::string& sequence = pRNA.get_sequence();

    if (node.children.empty()) {
        return ViennaDangles::get_multibranch_dangle_1(node, pRNA, vp);
    }

    // If a pseudoknot is mixed into the multiloop children, fall back to the
    // existing d1 machinery. ViennaRNA's ordinary fixed-structure d3 model does
    // not know about Knotergy pseudoknot children.
    for (const auto& child : node.children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            return ViennaDangles::get_multibranch_dangle_1(node, pRNA, vp);
        }
    }

    const std::vector<MultiloopStem> stems = build_multiloop_stems(node, sequence, vp);

    if (stems.empty()) {
        return 0;
    }

    const std::vector<size_t> pairings = pRNA.get_pairings();

    // First walk:
    // start from the multiloop closing pair. This disallows stacking of the
    // last child back into the closing pair at the final edge of this walk.
    int best = walk_multiloop_d3_from_start(stems, stems.size() - 1, sequence, pairings, vp);

    // Second walk:
    // start from the first child. This disallows stacking of the closing pair
    // back into the first child at the final edge of this walk.
    //
    // ViennaRNA does the same "walk around the loop twice" trick for d3.
    if (stems.size() > 1) {
        best = std::min(best, walk_multiloop_d3_from_start(stems, 0, sequence, pairings, vp));
    }

    return best;
}

}  // namespace knotergy