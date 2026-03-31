#pragma once

#include "../io/PseudoknotParams.hpp"
#include "../io/ViennaParams.hpp"
#include "../loop_tree/LoopNode.hpp"
#include "./ModifiedBasesFunctions.hpp"
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
     * @param processed_rna The processed RNA entry with structural information.
     * @param vp ViennaRNA model parameters.
     * @param mp Modified base parameters.
     * @param pkp Pseudoknot parameters.
     * @param is_inf If the energy is infinite (distance between base pairs < 3)
     * @param round Whether to round energy values (default: false).
     * @return Total pseudoknot energy in centicalories.
     */
    static double pseudoknot_energy(const LoopNode& node, const ProcessedRNAEntry& processed_rna,
                                    vrna_md_param& vp, const std::vector<modified_base_param>& mp,
                                    const pk_param& pkp, bool& is_inf, bool round = false);

   private:
    /**
     * @brief Calculate initialization penalty for a pseudoknot based on its parent loop type.
     *
     * Different penalties apply depending on whether the pseudoknot is in an external loop,
     * multibranch loop, or nested within another pseudoknot.
     *
     * @param node The pseudoknot loop node.
     * @param pkp Pseudoknot parameters.
     * @return Initialization penalty in centicalories.
     */
    [[nodiscard]] static double init_penalty(const LoopNode& node, const knotergy::pk_param& pkp);

    /**
     * @brief Calculate loop-specific energy penalties in a pseudoknot.
     *
     * Iterates through each band in the pseudoknot and calculates stacking, internal loop, and
     * multiloop energies for each base pair in the band.
     *
     * @param node The pseudoknot loop node.
     * @param processed_rna The processed RNA entry with structural information.
     * @param vp ViennaRNA model parameters.
     * @param mp Modified base parameters.
     * @param pkp Pseudoknot parameters.
     * @param round Whether to round energy values.
     * @param is_inf Whether an infinite energy condition was encountered.
     * @return Total loop penalties in centicalories.
     */
    [[nodiscard]] static double loop_penalties(const LoopNode& node,
                                               const ProcessedRNAEntry& processed_rna,
                                               vrna_md_param& vp,
                                               const std::vector<modified_base_param>& mp,
                                               const knotergy::pk_param& pkp, bool round,
                                               bool& is_inf);

    /**
     * @brief Calculate stacking energy for base pairs in a pseudoknot band.
     *
     * Applies the pseudoknot stacking multiplier to the standard ViennaRNA stacking energy.
     *
     * @param bp The outer base pair.
     * @param next_bp The inner (stacked) base pair.
     * @param processed_rna The processed RNA entry with structural information.
     * @param vp ViennaRNA model parameters.
     * @param mp Modified base parameters.
     * @param pkp Pseudoknot parameters.
     * @param round Whether to round the energy value.
     * @return Stacking energy with pseudoknot multiplier in centicalories.
     */
    [[nodiscard]] static double pk_stack_energy(const BasePair& bp, const BasePair& next_bp,
                                                const ProcessedRNAEntry& processed_rna,
                                                vrna_md_param& vp,
                                                const std::vector<modified_base_param>& mp,
                                                const knotergy::pk_param& pkp, bool round);

    /**
     * @brief Calculate internal loop energy for base pairs in a pseudoknot band.
     *
     * Applies the pseudoknot internal loop multiplier to the standard ViennaRNA internal loop
     * energy.
     *
     * @param bp The outer base pair.
     * @param next_bp The inner base pair.
     * @param processed_rna The processed RNA entry with structural information.
     * @param vp ViennaRNA model parameters.
     * @param pkp Pseudoknot parameters.
     * @param round Whether to round the energy value.
     * @return Internal loop energy with pseudoknot multiplier in centicalories.
     */
    [[nodiscard]] static double pk_internal_energy(const BasePair& bp, const BasePair& next_bp,
                                                   const ProcessedRNAEntry& processed_rna,
                                                   vrna_md_param& vp, const knotergy::pk_param& pkp,
                                                   bool round);

    /**
     * @brief Calculate multiloop energy nested between base pairs in a pseudoknot band.
     *
     * Computes pseudoknot-specific multiloop initialization, base pair, and unpaired penalties.
     *
     * @param bp The outer base pair.
     * @param next_bp The inner base pair.
     * @param processed_rna The processed RNA entry with structural information.
     * @param pkp Pseudoknot parameters.
     * @return Multiloop energy in centicalories.
     */
    [[nodiscard]] static double pk_multiloop_energy(const BasePair& bp, const BasePair& next_bp,
                                                    const ProcessedRNAEntry& processed_rna,
                                                    const knotergy::pk_param& pkp);
};

}  // namespace knotergy