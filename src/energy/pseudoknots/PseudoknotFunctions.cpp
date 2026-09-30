#include "PseudoknotFunctions.hpp"

#include "energy/modified_bases/ModInternal.hpp"
#include "energy/modified_bases/ModStack.hpp"

#include <cmath>
#include <iostream>

namespace knotergy {

double PseudoknotFunctions::pseudoknot_energy(LoopNode& node,
                                              const ProcessedRNAEntry& processed_rna,
                                              vrna_md_param& vp, const all_mod_params& mp,
                                              const pk_param& pkp, bool& is_inf,
                                              const bool pk_dangles) {
    // Where most of the energy calculations are done
    // Populate the breakdown with all the energies and information for reporting
    double energy = get_pk_energy(node, processed_rna, vp, mp, pkp, is_inf);

    // WIP, will add dangles to get_pk_energy once implementation is complete
    if (pk_dangles) {
        energy += pk_dangling_energy(node, processed_rna, vp, mp);
    }

    return energy;
}

double PseudoknotFunctions::get_pk_energy(LoopNode& node, const ProcessedRNAEntry& processed_rna,
                                          vrna_md_param& vp, const all_mod_params& mp,
                                          const pk_param& pkp, bool& is_inf) {
    // Calculate the number of unpaired bases outside of bands

    // Store breakdown information (excluding energies)
    node.unpaired_outside_bands_count = get_unpaired_outside_of_bands(node, processed_rna);
    int number_of_bands = static_cast<int>(node.bands.size());

    // Calculate the total energy of the pseudoknot
    double init_penalty = get_init_penalty(node, pkp);
    int band_penalty = pkp.band_penalty * number_of_bands;
    int unpaired_penalty = pkp.unpaired_in_pk * node.unpaired_outside_bands_count;
    int cr_penalty = pkp.cr_in_pk * node.number_of_outsideband_children;

    node.pk_level_energy = init_penalty + band_penalty + unpaired_penalty + cr_penalty;

    // Compute loop-specific energies for each band in the pseudoknot
    double loop_energies = compute_loop_energies(node, processed_rna, vp, mp, pkp, is_inf);

    return round_energy(node.pk_level_energy + loop_energies, pkp.round);
}

int PseudoknotFunctions::get_unpaired_outside_of_bands(const LoopNode& node,
                                                       const ProcessedRNAEntry& processed_rna) {
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

    return unpaired;
}

double PseudoknotFunctions::get_init_penalty(const LoopNode& node, const knotergy::pk_param& pkp) {
    // initialization penalties
    double energy = 0;
    switch (node.parent->loop_type) {
        case (LoopType::External):    energy += pkp.pk_in_ext; break;
        case (LoopType::Multibranch): energy += pkp.pk_in_mloop; break;
        case (LoopType::Pseudoknot):
            energy +=
                node.pseudo_type == PseudoNestedType::WithinBand ? pkp.pk_in_mloop : pkp.pk_in_pk;
            break;
        default:
            std::cerr << WARNING
                      << "Parent of this node is not a pseudoknot, external, or "
                         "multiloop"
                      << node << std::endl;
            break;
    }
    return energy;
}

double PseudoknotFunctions::compute_loop_energies(LoopNode& node,
                                                  const ProcessedRNAEntry& processed_rna,
                                                  vrna_md_param& vp, const all_mod_params& mp,
                                                  const knotergy::pk_param& pkp, bool& is_inf) {
    double energy = 0;

    for (Band& band : node.bands) {
        // Sanity check: left inner border must be less than right inner border
        if (band.left_inner() >= band.right_inner()) {
            THROW_ERROR("Invalid band with borders (" + std::to_string(band.left_border()) + ", " +
                        std::to_string(band.right_border()) + ") in pseudoknot (" +
                        std::to_string(node.begin) + ", " + std::to_string(node.end) +
                        "). Left inner border must be less than right inner border.");
        }

        std::vector<PKBasePair>& bps = band.base_pairs();

        // loops through each base pair in band (except last one)
        for (std::size_t idx = 0; idx + 1 < bps.size(); ++idx) {
            PKBasePair& bp = bps[idx];
            PKBasePair& next_bp = bps[idx + 1];
            switch (bp.loop_type) {
                case LoopType::Stack:
                    bp.energy = pk_stack_energy(bp, next_bp, processed_rna, vp, pkp, mp);
                    break;
                case LoopType::Internal:
                    bp.energy = pk_internal_energy(bp, next_bp, processed_rna, vp, pkp, mp);
                    break;
                case LoopType::Multibranch:
                    bp.energy = pk_multiloop_energy(bp, next_bp, processed_rna, pkp);
                    break;
                default:
                    THROW_ERROR("Invalid loop type for base pair (" + std::to_string(bp.i) + ", " +
                                std::to_string(bp.j) + ") in pseudoknot (" +
                                std::to_string(node.begin) + ", " + std::to_string(node.end) +
                                "). Loop type must be Stack, Internal, or Multibranch.");
            }
            energy += bp.energy;
        }

        // Handle the innermost base pair of the band (last base pair)
        PKBasePair& innermost_bp = bps.back();
        innermost_bp.energy = pk_innermost_energy(innermost_bp, vp, is_inf);
        energy += innermost_bp.energy;

        // Sanity check: the innermost base pair must have loop type InnermostBP
        if (innermost_bp.loop_type != LoopType::InnermostBP) {
            THROW_ERROR("Invalid loop type for innermost base pair (" +
                        std::to_string(innermost_bp.i) + ", " + std::to_string(innermost_bp.j) +
                        ") in pseudoknot (" + std::to_string(node.begin) + ", " +
                        std::to_string(node.end) + "). Loop type must be InnermostBP.");
        }
    }

    return energy;
}

double PseudoknotFunctions::pk_dangling_energy(const LoopNode& node,
                                               const ProcessedRNAEntry& processed_rna,
                                               vrna_md_param& vp, const all_mod_params& mp) {
    double energy = 0;

    for (const Band& band : node.bands) {
        std::size_t i = band.left_border();
        std::size_t j = band.right_border();

        if (i == node.begin || j == node.end) {
            // If the band is at the edge of the loop, skip it. This will be handled by the
            // multiloop or external loop energy calculation.
            continue;
        }

        auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(i, j, processed_rna, vp.md);

        if (n5d && n3d) {
            energy += ModBaseUtils::get_mismatch_energy(i, j, processed_rna, vp, mp);
        } else if (n5d) {
            energy += ModBaseUtils::get_dangle5_energy(i, j, processed_rna, vp, mp);
        } else if (n3d) {
            energy += ModBaseUtils::get_dangle3_energy(i, j, processed_rna, vp, mp);
        }
    }

    return energy;
}

// Checks if innermost base pair is infinite energy or not
double PseudoknotFunctions::pk_innermost_energy(const PKBasePair& bp, vrna_md_param& vp,
                                                bool& is_inf) {
    // check if the band is valid (has at least 3 base pairs to avoid infinite energy)
    std::size_t size = bp.j - bp.i - 1;

    if (size <= 30 && vp.p->hairpin[size] == INF) {
        std::cout << WARNING << "Band with borders (" << bp.i + 1 << ", " << bp.j + 1
                  << ") are too close (usually < 3 base pairs), resulting in infinite energy."
                  << std::endl;
        is_inf = true;
        return INF;
    }
    return 0;
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