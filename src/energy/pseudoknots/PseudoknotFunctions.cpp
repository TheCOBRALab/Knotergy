#include "PseudoknotFunctions.hpp"

#include "energy/modified_bases/ModInternal.hpp"
#include "energy/modified_bases/ModStack.hpp"

#include <cmath>
#include <iostream>

namespace knotergy {

double PseudoknotFunctions::pseudoknot_energy(const LoopNode& node,
                                              const ProcessedRNAEntry& processed_rna,
                                              vrna_md_param& vp, const all_mod_params& mp,
                                              const pk_param& pkp, bool& is_inf,
                                              const bool pk_dangles) {
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
    for (const LoopNode* child : node.children) {
        if (child->pseudo_type == PseudoNestedType::WithinBand) {
            unpaired += child->total_unpaired_bases_count;
        }
    }

    double energy = 0;

    energy += init_penalty(node, vp, pkp);
    energy += pkp.band_penalty * static_cast<int>(node.bands.size());
    energy += pkp.unpaired_in_pk * unpaired;
    energy += pkp.cr_in_pk * node.number_of_outsideband_children;
    energy += loop_penalties(node, processed_rna, vp, mp, pkp, is_inf);
    if (pk_dangles) {
        energy += pk_dangling_energy(node, processed_rna, vp, mp);
    }

    return energy;
}

double PseudoknotFunctions::init_penalty(const LoopNode& node, vrna_md_param& vp,
                                         const knotergy::pk_param& pkp) {
    // initialization penalties
    double energy = 0;
    switch (node.parent->loop_type) {
        case (LoopType::External):    energy += pkp.pk_in_ext; break;
        case (LoopType::Multibranch): energy += pkp.pk_in_mloop + vp.p->MLintern[1]; break;
        case (LoopType::Pseudoknot):
            energy +=
                node.pseudo_type == PseudoNestedType::WithinBand ? pkp.pk_in_mloop : pkp.pk_in_pk;
            break;
        default:
            std::cerr << "Warning: Parent of this node is not a pseudoknot, external, or "
                         "multiloop"
                      << node << std::endl;
            break;
    }
    return energy;
}

double PseudoknotFunctions::loop_penalties(const LoopNode& node,
                                           const ProcessedRNAEntry& processed_rna,
                                           vrna_md_param& vp, const all_mod_params& mp,
                                           const knotergy::pk_param& pkp, bool& is_inf) {
    double energy = 0;

    for (const Band& band : node.bands) {
        const std::vector<PKBasePair>& bps = band.base_pairs();
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
            const PKBasePair& bp = bps[idx];
            const PKBasePair& next_bp = bps[idx + 1];

            if (bp.is_stack(next_bp)) {
                energy += pk_stack_energy(bp, next_bp, processed_rna, vp, pkp, mp);
            } else if (bp.children.empty()) {
                // if no nested structure between two base pairs of a band, it's an internal loop
                energy += pk_internal_energy(bp, next_bp, processed_rna, vp, pkp, mp);
            } else {
                energy += pk_multiloop_energy(bp, next_bp, processed_rna, pkp);
            }
        }
    }

    return energy;
}

double PseudoknotFunctions::pk_dangling_energy(const LoopNode& node,
                                               const ProcessedRNAEntry& processed_rna,
                                               vrna_md_param& vp, const all_mod_params& mp) {
    double energy = 0;
    const std::string& sequence = processed_rna.get_sequence();
    const std::vector<std::string_view>& mod_sequence = processed_rna.get_modified_sequence();

    for (const Band& band : node.bands) {
        size_t i = band.left_border();
        size_t j = band.right_border();

        unsigned int pair_type = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
        auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(i, j, processed_rna, vp.md);

        n5d = vp.md.dangles != 0 ? n5d : -1;
        n3d = vp.md.dangles != 0 ? n3d : -1;

        if (n5d > 0) {
            int n5d_unmod_energy = vp.p->dangle5[pair_type][n5d];
            if (processed_rna.has_modified_bases()) {
                int mod_d5 = ModBaseUtils::get_dangle5_mod_energy(i, j, mod_sequence, mp);
                energy += mod_d5 != NULL_ENERGY ? mod_d5 : n5d_unmod_energy;
            } else {
                energy += n5d_unmod_energy;
            }
        }

        if (n3d > 0) {
            int n3d_unmod_energy = vp.p->dangle3[pair_type][n3d];
            if (processed_rna.has_modified_bases()) {
                int mod_d3 = ModBaseUtils::get_dangle3_mod_energy(i, j, mod_sequence, mp);
                energy += mod_d3 != NULL_ENERGY ? mod_d3 : n3d_unmod_energy;
            } else {
                energy += n3d_unmod_energy;
            }
        }
    }

    return energy;
}

double PseudoknotFunctions::pk_stack_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                            const ProcessedRNAEntry& processed_rna,
                                            vrna_md_param& vp, const knotergy::pk_param& pkp,
                                            const all_mod_params& mp) {
    const std::string& sequence = processed_rna.get_sequence();

    int stack_energy = processed_rna.has_modified_bases()
                           ? ModStack::find_mod_stack_energy(bp, next_bp, processed_rna, vp, mp)
                           : ViennaFunctions::stack_energy(bp, next_bp, sequence, vp);

    double stack_pk_energy = stack_energy * pkp.pk_stack_x;
    return round_energy(stack_pk_energy, pkp.round);
}

double PseudoknotFunctions::pk_internal_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                               const ProcessedRNAEntry& processed_rna,
                                               vrna_md_param& vp, const knotergy::pk_param& pkp,
                                               const all_mod_params& mp) {
    const std::string& sequence = processed_rna.get_sequence();
    int internal_energy =
        processed_rna.has_modified_bases()
            ? ModInternal::find_mod_internal_energy(bp, next_bp, processed_rna, vp, mp)
            : ViennaFunctions::internal_loop_energy(bp, next_bp, sequence, vp);
    double internal_pk_energy = internal_energy * pkp.pk_internal_x;

    return round_energy(internal_pk_energy, pkp.round);
}

double PseudoknotFunctions::pk_multiloop_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                                const ProcessedRNAEntry& processed_rna,
                                                const knotergy::pk_param& pkp) {
    double multiloop_penalty = pkp.pk_mloop_init;

    // Since a multiloop is nested between two base pairs, we add 2 * bp_penalty
    // plus the number of children * bp_penalty for each child nested within the multiloop

    // Personal note: I find it weird that the number of children is what used
    // for base pair penalty. If a child is a pseudoknot, it can have multiple base pairs.
    // Like an H-type pseudoknot has 2 bands. Why does it only get 1 base pair penalty?
    // But this is how the original HotKnotsV2 implementation did it. HFold and other programs
    // also use the same convention, so we have to do it too for consistency.
    // But if you're reading this, maybe this could be a paper? idk.
    multiloop_penalty +=
        pkp.pk_mloop_bp * 2 + static_cast<int>(bp.children.size()) * pkp.pk_mloop_bp;

    // Get unpaired bases between the two base pairs of the multiloop
    // then subtract any unpaired bases that are part of children
    int unpaired = processed_rna.get_unpaired_count(bp.i, next_bp.i);
    unpaired += processed_rna.get_unpaired_count(next_bp.j, bp.j);
    for (ClosedRegion nested_cr : bp.children) {
        unpaired -= processed_rna.get_unpaired_count(nested_cr.begin, nested_cr.end);
    }

    // get unpaired penalty
    multiloop_penalty += unpaired * pkp.pk_mloop_unpaired;

    return multiloop_penalty;
}

double PseudoknotFunctions::round_energy(double energy, RoundMethod round) {
    switch (round) {
        case RoundMethod::None:           return energy;              // no rounding
        case RoundMethod::Bankers:        return std::rint(energy);   // banker's rounding
        case RoundMethod::RoundToNearest: return std::round(energy);  // round to nearest integer
        case RoundMethod::RoundDown:      return std::floor(energy);  // round down
        case RoundMethod::RoundUp:        return std::ceil(energy);   // round up
        case RoundMethod::Truncate:       return std::trunc(energy);        // truncate
        default:
            THROW_ERROR("Invalid round value: " + std::to_string(static_cast<int>(round)) +
                        ". Valid values are 0 (None), 1 (Bankers), 2 (RoundToNearest), 3 "
                        "(RoundDown), 4 (RoundUp), 5 (Truncate).");
    }
}
}  // namespace knotergy