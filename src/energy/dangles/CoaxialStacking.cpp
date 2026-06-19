#include "CoaxialStacking.hpp"

#include "energy/modified_bases/ModStack.hpp"

namespace knotergy {

int CoaxialStacking::add_or_inf(int a, int b) {
    if (a >= INF || b >= INF) {
        return INF;
    }
    return a + b;
}

// applies to closing pair and first child
int CoaxialStacking::compute_initial_ld5_for_d3(const MultiloopStem& stem,
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

std::vector<MultiloopStem> CoaxialStacking::populate_multiloop_stems(const LoopNode& node,
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

int CoaxialStacking::walk_multiloop_d3_from_start(const ProcessedRNAEntry& pRNA, size_t start_prev,
                                                  const std::vector<MultiloopStem>& stems,
                                                  vrna_md_param& vp, const all_mod_params& mp) {
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

int CoaxialStacking::get_multibranch_dangle_3(const LoopNode& node,
                                              const std::vector<MultiloopStem>& stems,
                                              const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                              const all_mod_params& mp) {
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

int CoaxialStacking::get_multibranch_dangle_3(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                              vrna_md_param& vp, const all_mod_params& mp) {
    std::vector<MultiloopStem> stems = populate_multiloop_stems(node, pRNA, vp);
    return get_multibranch_dangle_3(node, stems, pRNA, vp, mp);
};

}  // namespace knotergy