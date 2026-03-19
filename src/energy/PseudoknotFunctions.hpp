#pragma once

#include <iostream>

#include "../io/PseudoknotParams.hpp"
#include "../io/ViennaParams.hpp"
#include "../loop_tree/LoopNode.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
/**
 * @brief Energy calculation functions for pseudoknotted RNA structures.
 *
 * This class implements the pseudoknot energy model, which includes penalties for
 * pseudoknot initialization, bands, unpaired bases, and nested structures within
 * pseudoknots.
 */
class PseudoknotFunctions {
   public:
    /**
     * @brief Calculate the total energy of a pseudoknot loop.
     *
     * Computes initialization penalties, band penalties, unpaired base penalties,
     * and loop-specific energies for a pseudoknotted structure.
     *
     * @param node The loop node representing the pseudoknot.
     * @param sequence The RNA nucleotide sequence.
     * @param processed_rna The processed RNA entry with structural information.
     * @param round Whether to round energy values (default: false).
     * @param is_inf If the energy is infinite (distance between base pairs < 3)
     * @return Total pseudoknot energy in centicalories.
     */
    static double pseudoknot_energy(const LoopNode& node, const std::string& sequence,
                                    const ProcessedRNAEntry& processed_rna,  
                                    vrna_md_param& vp, const knotergy::pk_param& pkp,
                                    bool& is_inf, bool round = false) {
                                        
        int unpaired = node.exclusive_unpaired_bases_count;
        
        // remove unpaired bases within bands since they're already included in ViennaRNA's energy calculations for internal loops
        for (Band band : node.bands) {
            unpaired -= processed_rna.get_unpaired_count(band.left_border(), band.left_inner());
            unpaired -= processed_rna.get_unpaired_count(band.right_inner(), band.right_border());
        }

        // Previous loop removed ALL unpaired bases within bands, this includes base pairs of children that are within the band.
        // Since the base pairs of all children were already removed in exclusive_unpaired_bases_count, we need to add them back
        // due to double counting. We can identify these base pairs as the children that are within bands (pseudo_type == WithinBand)
        for (std::shared_ptr<LoopNode> child : node.children) {
            if (child->pseudo_type == PseudoNestedType::WithinBand) {
                unpaired += child->total_unpaired_bases_count;
            }
        }

        double energy = 0;

        energy += PseudoknotFunctions::init_penalty(node, pkp);
        energy += pkp.band * node.number_of_bands;
        energy += pkp.unpaired_in_pk * unpaired;
        energy += pkp.cr_in_pk * node.number_of_nested_children;
        energy += PseudoknotFunctions::loop_penalties(node, sequence, processed_rna, vp, pkp, round, is_inf);

        // Children that are within a band are nested inside of a pseudoknotted multiloop
        // Add the base pair penalty for each child that is within a band
        for (std::shared_ptr<LoopNode> c : node.children) {
            if (c->pseudo_type == PseudoNestedType::WithinBand) {
                energy += pkp.pk_mloop_bp;
            }
        }
        return energy;
    }

   private:
    /**
     * @brief Calculate initialization penalty for a pseudoknot based on its parent loop type.
     *
     * Different penalties apply depending on whether the pseudoknot is in an external loop,
     * multibranch loop, or nested within another pseudoknot.
     *
     * @param node The pseudoknot loop node.
     * @return Initialization penalty in centicalories.
     */
    [[nodiscard]] static double init_penalty(const LoopNode& node, const knotergy::pk_param& pkp) {
        // initialization penalties
        double energy = 0;
        if (std::shared_ptr<LoopNode> parent = node.parent.lock()) {
            switch (parent->loop_type) {
                case (LoopType::External):
                    energy += pkp.pk_in_ext;
                    break;
                case (LoopType::Multibranch):
                    energy += pkp.pk_in_mloop;
                    break;
                case (LoopType::Pseudoknot):
                    energy += node.pseudo_type == PseudoNestedType::WithinBand
                                  ? pkp.pk_in_mloop
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

    /**
     * @brief Calculate loop-specific energy penalties in a pseudoknot.
     *
     * Iterates through each band in the pseudoknot and calculates stacking, internal loop, and multiloop
     * energies for each base pair in the band. 
     *
     * @param node The pseudoknot loop node.
     * @param sequence The RNA nucleotide sequence.
     * @param processed_rna The processed RNA entry with structural information.
     * @param round Whether to round energy values.
     * @return Total loop penalties in centicalories.
     */
    [[nodiscard]] static double loop_penalties(const LoopNode& node, const std::string& sequence,
                                               const ProcessedRNAEntry& processed_rna, vrna_md_param& vp,
                                               const knotergy::pk_param& pkp, bool round, bool& is_inf) {
        double energy = 0;

        for (const Band& band : node.bands) {
            const std::vector<BasePair>& bps = band.base_pairs();
            const size_t n = bps.size();

            // check if the band is valid (has at least 3 base pairs to avoid infinite energy)
            if (band.right_inner() - band.left_inner() < 4) {
                std::cout << "Warning: Band with borders (" << band.left_border() << ", " << band.right_border()
                << ") has less than 3 unpaired bases between its inner borders, resulting in infinite energy." << std::endl;
                energy = INF;
                is_inf = true;
                continue;
            }

            // loops through each base pair in band (except last one)
            for (size_t idx = 0; idx + 1 < n; ++idx) {
                const BasePair& bp = bps[idx];
                const BasePair& next_bp = bps[idx + 1];

                if (bp.is_stack(next_bp)) {
                    energy += PseudoknotFunctions::pk_stack_energy(bp, next_bp, sequence, vp, pkp, round);
                } else if (bp.children.empty()) {
                    // if no nested structure between two base pairs of a band, it's an internal loop
                    energy += PseudoknotFunctions::pk_internal_energy(bp, next_bp, sequence, vp, pkp, round);
                } else {
                    energy += PseudoknotFunctions::pk_multiloop_energy(bp, next_bp, processed_rna, pkp);
                }
            }
        }
        return energy;
    }

    /**
     * @brief Calculate stacking energy for base pairs in a pseudoknot band.
     *
     * Applies the pseudoknot stacking multiplier to the standard ViennaRNA stacking energy.
     *
     * @param bp The outer base pair.
     * @param next_bp The inner (stacked) base pair.
     * @param sequence The RNA nucleotide sequence.
     * @param round Whether to round the energy value.
     * @return Stacking energy with pseudoknot multiplier in centicalories.
     */
    [[nodiscard]] static double pk_stack_energy(const BasePair& bp, const BasePair& next_bp,
                                                const std::string& sequence,
                                                vrna_md_param& vp,
                                                const knotergy::pk_param& pkp,
                                                const bool& round) {
        double stack_penalty = ViennaFunctions::stack_energy(bp, next_bp, sequence, vp) *
                                                                                pkp.pk_stack_x;
        return round ? std::round(stack_penalty) : stack_penalty;
    }

    /**
     * @brief Calculate internal loop energy for base pairs in a pseudoknot band.
     *
     * Applies the pseudoknot internal loop multiplier to the standard ViennaRNA internal loop energy.
     *
     * @param bp The outer base pair.
     * @param next_bp The inner base pair.
     * @param sequence The RNA nucleotide sequence.
     * @param round Whether to round the energy value.
     * @return Internal loop energy with pseudoknot multiplier in centicalories.
     */
    [[nodiscard]] static double pk_internal_energy(const BasePair& bp, const BasePair& next_bp,
                                                   const std::string& sequence, 
                                                   vrna_md_param& vp,
                                                   const knotergy::pk_param& pkp,
                                                   const bool& round) {
        double internal_penalty = ViennaFunctions::internal_loop_energy(bp, next_bp, sequence, vp) *
                                  pkp.pk_internal_x;

        return round ? std::round(internal_penalty) : internal_penalty;
    }

    /**
     * @brief Calculate multiloop energy nested between base pairs in a pseudoknot band.
     *
     * Computes pseudoknot-specific multiloop initialization, base pair, and unpaired penalties.
     *
     * @param bp The outer base pair.
     * @param next_bp The inner base pair.
     * @param processed_rna The processed RNA entry with structural information.
     * @return Multiloop energy in centicalories.
     */
    [[nodiscard]] static double pk_multiloop_energy(const BasePair& bp, const BasePair& next_bp,
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
};

}  // namespace knotergy
