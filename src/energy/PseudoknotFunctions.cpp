#include "PseudoknotFunctions.hpp"

#include <cmath>
#include <iostream>

namespace knotergy {

double PseudoknotFunctions::pseudoknot_energy(const LoopNode& node,
                                              const ProcessedRNAEntry& processed_rna,
                                              vrna_md_param& vp,
                                              const std::vector<modified_base_param>& mp,
                                              const pk_param& pkp, bool& is_inf, bool round) {
    int unpaired = node.exclusive_unpaired_bases_count;

    // remove unpaired bases within bands since they're already included in ViennaRNA's energy
    // calculations for internal loops
    for (const Band& band : node.bands) {
        unpaired -= processed_rna.get_unpaired_count(band.left_border(), band.left_inner());
        unpaired -= processed_rna.get_unpaired_count(band.right_inner(), band.right_border());
    }

    // Previous loop removed ALL unpaired bases within bands, this includes base pairs of children
    // that are within the band. Since the base pairs of all children were already removed in
    // exclusive_unpaired_bases_count, we need to add them back due to double counting. We can
    // identify these base pairs as the children that are within bands (pseudo_type == WithinBand)
    for (const std::shared_ptr<LoopNode>& child : node.children) {
        if (child->pseudo_type == PseudoNestedType::WithinBand) {
            unpaired += child->total_unpaired_bases_count;
        }
    }

    double energy = 0;

    energy += PseudoknotFunctions::init_penalty(node, pkp);
    energy += pkp.band * node.number_of_bands;
    energy += pkp.unpaired_in_pk * unpaired;
    energy += pkp.cr_in_pk * node.number_of_nested_children;
    energy += PseudoknotFunctions::loop_penalties(node, processed_rna, vp, mp, pkp, round, is_inf);

    // Children that are nested within a band are considered to be in a pseudoknotted multiloop, 
    // so we multiply the number of children within bands by the multiloop base pair penalty.
    energy += node.number_of_withinband_children * pkp.pk_mloop_bp;

    // Personal note: I find it dumb that the number of children is what used for base pair penalty
    // If a child is a pseudoknot, it can have multiple base pairs.
    // Like an H-type pseudoknot has 2 bands. Why tf does it only get 1 base pair penalty?
    // But this is how the original HotKnotsV2 implementation did it. HFold and other programs
    // also use the same convention, so I guess we have to do it too for consistency.
    // But if you're reading this, maybe this could be a paper? idk.

    return energy;
}

double PseudoknotFunctions::init_penalty(const LoopNode& node, const knotergy::pk_param& pkp) {
    // initialization penalties
    double energy = 0;
    if (const std::shared_ptr<LoopNode>& parent = node.parent.lock()) {
        switch (parent->loop_type) {
            case (LoopType::External):
                energy += pkp.pk_in_ext;
                break;
            case (LoopType::Multibranch):
                energy += pkp.pk_in_mloop;
                break;
            case (LoopType::Pseudoknot):
                energy += node.pseudo_type == PseudoNestedType::WithinBand ? pkp.pk_in_mloop
                                                                           : pkp.pk_in_pk;
                break;
            default:
                std::cerr << "Warning: Parent of this node is not a pseudoknot, external, or "
                             "multiloop"
                          << node << std::endl;
                break;
        }
    } else {
        THROW_ERROR("Parent node of pseudoknot (" + std::to_string(node.begin) + ", " +
                    std::to_string(node.end) + ") has expired.");
    }
    return energy;
}

double PseudoknotFunctions::loop_penalties(const LoopNode& node,
                                           const ProcessedRNAEntry& processed_rna,
                                           vrna_md_param& vp,
                                           const std::vector<modified_base_param>& mp,
                                           const knotergy::pk_param& pkp, bool round,
                                           bool& is_inf) {
    double energy = 0;

    for (const Band& band : node.bands) {
        const std::vector<BasePair>& bps = band.base_pairs();
        const size_t n = bps.size();

        // Sanity check: left inner border must be less than right inner border
        if (band.left_inner() >= band.right_inner()) {
            THROW_ERROR("Invalid band with borders (" + std::to_string(band.left_border()) + ", " +
                        std::to_string(band.right_border()) + ") in pseudoknot (" +
                        std::to_string(node.begin) + ", " + std::to_string(node.end) +
                        "). Left inner border must be less than right inner border.");
        }

        // check if the band is valid (has at least 3 base pairs to avoid infinite energy)
        size_t size = band.right_inner() - band.left_inner() - 1;

        if (size <= 30 && vp.p->hairpin[size] == INF) {
            std::cout << "Warning: Band with borders (" << band.left_border() << ", "
                      << band.right_border()
                      << ") are too close (usually < 3 base pairs), resulting in infinite energy."
                      << std::endl;
            energy += INF;
            is_inf = true;
            continue;
        }

        // loops through each base pair in band (except last one)
        for (size_t idx = 0; idx + 1 < n; ++idx) {
            const BasePair& bp = bps[idx];
            const BasePair& next_bp = bps[idx + 1];

            if (bp.is_stack(next_bp)) {
                energy += PseudoknotFunctions::pk_stack_energy(bp, next_bp, processed_rna, vp, mp,
                                                               pkp, round);
            } else if (bp.children.empty()) {
                // if no nested structure between two base pairs of a band, it's an internal loop
                energy += PseudoknotFunctions::pk_internal_energy(bp, next_bp, processed_rna, vp,
                                                                  pkp, round);
            } else {
                energy += PseudoknotFunctions::pk_multiloop_energy(bp, next_bp, processed_rna, pkp);
            }
        }
    }

    return energy;
}

double PseudoknotFunctions::pk_stack_energy(const BasePair& bp, const BasePair& next_bp,
                                            const ProcessedRNAEntry& processed_rna,
                                            vrna_md_param& vp,
                                            const std::vector<modified_base_param>& mp,
                                            const knotergy::pk_param& pkp, bool round) {
    const std::string& sequence = processed_rna.get_sequence();

    int stack_energy =
        processed_rna.has_modified_bases()
            ? ModifiedBasesFunctions::find_mod_stack_energy(bp, next_bp, processed_rna, vp, mp)
            : ViennaFunctions::stack_energy(bp, next_bp, sequence, vp);

    double stack_penalty = stack_energy * pkp.pk_stack_x;
    return round ? std::round(stack_penalty) : stack_penalty;
}

double PseudoknotFunctions::pk_internal_energy(const BasePair& bp, const BasePair& next_bp,
                                               const ProcessedRNAEntry& processed_rna,
                                               vrna_md_param& vp, const knotergy::pk_param& pkp,
                                               bool round) {
    const std::string& sequence = processed_rna.get_sequence();
    double internal_penalty =
        ViennaFunctions::internal_loop_energy(bp, next_bp, sequence, vp) * pkp.pk_internal_x;

    return round ? std::round(internal_penalty) : internal_penalty;
}

double PseudoknotFunctions::pk_multiloop_energy(const BasePair& bp, const BasePair& next_bp,
                                                const ProcessedRNAEntry& processed_rna,
                                                const knotergy::pk_param& pkp) {
    double multiloop_penalty = pkp.pk_mloop_init;

    // Since a multiloop is nested between two base pairs, we add 2 * bp_penalty
    // We add the child base pairs at a different part of the energy calculation
    multiloop_penalty += pkp.pk_mloop_bp * 2;

    // Get unpaired bases between the two base pairs of the multiloop
    // then subtract any unpaired bases that are part of children
    int unpaired = processed_rna.get_unpaired_count(bp.i, next_bp.i);
    unpaired += processed_rna.get_unpaired_count(next_bp.j, bp.j);
    for (ClosedRegion nested_cr : bp.children) {
        unpaired -= processed_rna.get_unpaired_count(nested_cr.begin, nested_cr.end);
    }

    // get unpaired penalty
    int pk_mloop_unpaired_energy = unpaired * pkp.pk_mloop_unpaired;
    multiloop_penalty += pk_mloop_unpaired_energy;

    return multiloop_penalty;
}

}  // namespace knotergy