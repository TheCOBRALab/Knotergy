#pragma once

#include "energy/vienna/ViennaFunctions.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"

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
     * @param breakdown Reference to a PKEnergyBreakdown object to store detailed energy components.
     * @param pk_dangles Whether to include pseudoknot dangle energy contributions (default: false).
     * @return Total pseudoknot energy in centicalories.
     */
    [[nodiscard]] static double pseudoknot_energy(const LoopNode& node,
                                                  const ProcessedRNAEntry& processed_rna,
                                                  vrna_md_param& vp, const all_mod_params& mp,
                                                  const pk_param& pkp, PKEnergyBreakdown& breakdown,
                                                  const bool pk_dangles = false);

   private:
    static void populate_pk_energy_breakdown(const LoopNode& node,
                                             const ProcessedRNAEntry& processed_rna,
                                             vrna_md_param& vp, const all_mod_params& mp,
                                             const pk_param& pkp, PKEnergyBreakdown& breakdown);

    /**
     * @brief Calculate initialization penalty for a pseudoknot based on its parent loop type.
     *
     * Different penalties apply depending on whether the pseudoknot is in an external loop,
     * multibranch loop, or nested within another pseudoknot.
     *
     * @param node The pseudoknot loop node.
     * @param vp ViennaRNA model parameters.
     * @param pkp Pseudoknot parameters.
     * @return Initialization penalty in centicalories.
     */
    [[nodiscard]] static double init_penalty(const LoopNode& node, const knotergy::pk_param& pkp);

    /**
     * @brief Calculate the number of unpaired bases in a pseudoknot, excluding those within bands.
     *
     * This function counts unpaired bases in the pseudoknot loop while excluding those that
     * are part of pseudoknot bands, as they are already accounted for in ViennaRNA's energy
     * calculations for internal loops. And multiloops spanning bands are handled separately
     * in the pseudoknot energy model.
     *
     * @param node The pseudoknot loop node.
     * @param processed_rna The processed RNA entry with structural information.
     * @return Count of unpaired bases outside of bands.
     */
    [[nodiscard]] static int get_unpaired_outside_of_bands(const LoopNode& node,
                                                           const ProcessedRNAEntry& processed_rna);

    /**
     * @brief Calculate loop-specific energies in a pseudoknot.
     *
     * Iterates through each band in the pseudoknot and calculates stacking, internal loop, and
     * multiloop energies for each base pair in the band.
     *
     * @param node The pseudoknot loop node.
     * @param processed_rna The processed RNA entry with structural information.
     * @param vp ViennaRNA model parameters.
     * @param mp Modified base parameters.
     * @param pkp Pseudoknot parameters.
     * @param breakdown Reference to a PKEnergyBreakdown object to store detailed energy components.
     * @return Total loop energies in centicalories.
     */
    static void loop_energies(const LoopNode& node, const ProcessedRNAEntry& processed_rna,
                              vrna_md_param& vp, const all_mod_params& mp,
                              const knotergy::pk_param& pkp, PKEnergyBreakdown& breakdown);

    /**
     * @brief Checks if the innermost base pair has infinite energy.
     *
     * This function checks if the innermost base pair of a pseudoknot band is too close
     * (less than 3 base pairs apart)
     *
     * @param bp The innermost base pair of the pseudoknot band.
     * @param vp ViennaRNA model parameters.
     * @param is_inf Reference to a boolean that will be set to true if the energy is infinite.
     * @return The energy of the innermost base pair (0 if valid, INF if invalid).
     */
    [[nodiscard]] static double pk_innermost_energy(const PKBasePair& bp, vrna_md_param& vp,
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
     * @return Stacking energy with pseudoknot multiplier in centicalories.
     */
    [[nodiscard]] static double pk_stack_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                                const ProcessedRNAEntry& processed_rna,
                                                vrna_md_param& vp, const knotergy::pk_param& pkp,
                                                const all_mod_params& mp);

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
     * @return Internal loop energy with pseudoknot multiplier in centicalories.
     */
    [[nodiscard]] static double pk_internal_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                                   const ProcessedRNAEntry& processed_rna,
                                                   vrna_md_param& vp, const knotergy::pk_param& pkp,
                                                   const all_mod_params& mp);

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
    [[nodiscard]] static double pk_multiloop_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                                    const ProcessedRNAEntry& processed_rna,
                                                    const knotergy::pk_param& pkp);

    [[nodiscard]] static double round_energy(double energy, RoundMethod round);

    /**
     * @brief Calculate dangling end energy for a loop node in a pseudoknot.
     *
     * @param node The loop node for which to calculate dangling end energy.
     * @param processed_rna The processed RNA entry with structural information.
     * @param vp ViennaRNA model parameters.
     * @param mp Modified base parameters.
     * @param pkp Pseudoknot parameters.
     * @param is_inf A reference to a boolean indicating if the energy is infinite.
     * @return Dangling end energy in centicalories.
     */
    [[nodiscard]] static double pk_dangling_energy(const LoopNode& node,
                                                   const ProcessedRNAEntry& processed_rna,
                                                   vrna_md_param& vp, const all_mod_params& mp);
};

}  // namespace knotergy