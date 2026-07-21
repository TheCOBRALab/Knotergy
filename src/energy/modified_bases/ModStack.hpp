#pragma once

#include "energy/modified_bases/ModBaseUtils.hpp"
#include "energy/vienna/ViennaFunctions.hpp"

namespace knotergy {
class ModStack {
   public:
    /**
     * @brief Calculate stacking energy for base pairs with modified bases.
     *
     * @param i 5' position of outer base pair.
     * @param j 3' position of outer base pair.
     * @param ci 5' position of inner base pair.
     * @param cj 3' position of inner base pair.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @param mp Vector of modified base parameters.
     * @return Stacking energy in centicalories, accounting for modified bases.
     */
    // Gets the modified energy of a stack
    static int find_mod_stack_energy(size_t i, size_t j, size_t ci, size_t cj,
                                     const std::string&                   sequence,
                                     const std::vector<std::string_view>& mod_sequence,
                                     vrna_md_param& vp, const all_mod_params& mp) {
        // Get unmodified energy first to use as fallback if no modified nucleotides are found
        int unmod_energy = ViennaFunctions::stack_energy(i, j, ci, cj, sequence, vp);

        // Find all modified bases at the inner edge of the stack (i, j, i+1, j-1)
        std::vector<std::string_view> unique_mod_bases =
            ModBaseUtils::unique_modified_bases_at_inner_edge(i, j, mod_sequence);
        if (unique_mod_bases.empty()) return unmod_energy;

        // Used to look up stacking energies in modified base parameters
        std::string l_key = ModBaseUtils::join_string_views({i, ci, j, cj}, mod_sequence);

        // Get mod base energy correction (returns original energy if no modifications found)
        int e = ModBaseUtils::get_mod_energy(l_key, unique_mod_bases, mp, unmod_energy,
                                             ModLookup::Stacking);

        // Due to the symmetrical nature of stacks, if the key was not found, we check the reverse
        // order (ci, i, cj, j) for modified bases
        if (e == unmod_energy) {
            std::string r_key = ModBaseUtils::join_string_views({cj, j, ci, i}, mod_sequence);
            e = ModBaseUtils::get_mod_energy(r_key, unique_mod_bases, mp, unmod_energy,
                                             ModLookup::Stacking);
        }

        return e;
    }

    /**
     * @brief Overload of find_mod_stack_energy that takes BasePair objects and a ProcessedRNAEntry.
     *
     * This is a convenience function that extracts the necessary information from the
     * BasePair objects and ProcessedRNAEntry to call the main find_mod_stack_energy function.
     *
     * @param bp The outer base pair.
     * @param next_bp The inner base pair that stacks with the outer base pair.
     * @param processed_rna The processed RNA entry containing sequence and modifications.
     * @param vp ViennaRNA model and parameters.
     * @param mp Vector of modified base parameters.
     * @return Stacking energy in centicalories, accounting for modified bases.
     *
     */
    [[nodiscard]] static int find_mod_stack_energy(const BasePair& bp, const BasePair& next_bp,
                                                   const ProcessedRNAEntry& processed_rna,
                                                   vrna_md_param& vp, const all_mod_params& mp) {
        const std::string&                   sequence     = processed_rna.get_sequence();
        const std::vector<std::string_view>& mod_sequence = processed_rna.get_modified_sequence();
        return find_mod_stack_energy(bp.i, bp.j, next_bp.i, next_bp.j, sequence, mod_sequence, vp,
                                     mp);
    }
};
}  // namespace knotergy