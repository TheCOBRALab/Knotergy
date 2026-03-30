#pragma once

#include <utility>

#include "../preprocessing/RNAProcessor.hpp"
#include "./ViennaDangles.hpp"
#include "./ViennaFunctions.hpp"

namespace knotergy {

/**
 * @brief Enumeration of modified base energy lookup types.
 *
 * Specifies which type of energy parameter to look up for modified bases.
 */
enum class ModLookup { Stacking, Terminal, Mismatch, Dangle5, Dangle3 };

// Stores the differences in energy contributions due to modified bases
struct ModDiffs {
   ModDiffs(int terminal_diff, int mismatch_diff, int n5d_diff, int n3d_diff):
         terminalAU{terminal_diff},
         mismatch{mismatch_diff},
         n5d{n5d_diff},
         n3d{n3d_diff} {}
    const int terminalAU;
    const int mismatch;
    const int n5d;
    const int n3d;
};

class ModifiedBasesFunctions {
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
     * @param mod_params Vector of modified base parameters.
     * @return Stacking energy in centicalories, accounting for modified bases.
     */
    static int find_mod_stack_energy(const size_t& i, const size_t& j, const size_t& ci,
                                        const size_t& cj, const std::string& sequence,
                                        const std::vector<std::string_view>& mod_sequence,
                                        vrna_md_param& vp,
                                        const std::vector<modified_base_param>& mod_params);

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
   static int find_mod_stack_energy(const BasePair& bp, const BasePair& next_bp,
                                    const ProcessedRNAEntry& processed_rna,
                                    vrna_md_param& vp, 
                                    const std::vector<modified_base_param>& mp) {
      const std::string& sequence = processed_rna.get_sequence();
      const std::vector<std::string_view>& mod_sequence = processed_rna.get_modified_sequence();
      return find_mod_stack_energy(bp.i, bp.j, next_bp.i, next_bp.j, sequence, mod_sequence, vp, mp);
    }

   /**
     * @brief Calculate multiloop energy with modified bases.
     *
     * @param node The multiloop node for which to calculate energy.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @param mod_params Vector of modified base parameters.
     * @return Multiloop energy in centicalories, accounting for modified bases.
     */
    static int find_mod_multiloop_energy(const LoopNode& node, const std::string& sequence,
                                           const std::vector<std::string_view>& mod_sequence,
                                           vrna_md_param& vp,
                                           const std::vector<modified_base_param>& mod_params);

    /**
     * @brief Calculate external loop energy with modified bases.
     *
     * @param children Vector of child loop nodes in the external loop.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @param mod_params Vector of modified base parameters.
     * @return External loop energy in centicalories, accounting for modified bases.
     */
    static int find_mod_external_energy(const std::vector<std::shared_ptr<LoopNode>>& children, 
                                        const std::string& sequence,
                                        const std::vector<std::string_view>& mod_sequence,
                                        vrna_md_param& vp,
                                        const std::vector<modified_base_param>& mod_params);


   
   /**
    * @brief Update the energy of a loop based on modified bases.
    * 
    * Supports updating the energy value for 3 things
    * 1) Children of an external loop (if is_external = true  and is_closing = false)
    * 2) Closing pair of a multiloop  (if is_external = false and is_closing = true)
    * 3) Children of a multiloop      (if is_external = false and is_closing = false)
    * 
    * @param node The loop node to update energy for.
    * @param sequence The unmodified RNA nucleotide sequence.
    * @param mod_sequence The modified RNA sequence (grapheme views).
    * @param vp ViennaRNA model and parameters.
    * @param mod_params Vector of modified base parameters.
    * @param current_set The DangleSet of the current loop (used if dangles == 1).
    * @param is_external Whether this loop is an external loop (default false).
    * @param is_closing Whether this loop is a closing pair (default false).
    * @return Energy difference in centicalories to apply to the loop's energy due to modified bases.
    * 
    */
   static int update_energy(const LoopNode& node, const std::string& sequence, 
                            const std::vector<std::string_view>& mod_sequence,
                            vrna_md_param& vp,const std::vector<modified_base_param>& mod_params, 
                            DangleSet& current_set, bool is_external = false, 
                            bool is_closing = false);

   private:
    /**
     * @brief Find unique modified bases at specified sequence positions.
     *
     * @param indices Vector of sequence indices to check.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @return Vector of unique modified base string views found at those positions.
     */
    static std::vector<std::string_view> unique_modified_bases_at_indices(
                   std::vector<size_t> indices, const std::vector<std::string_view>& mod_sequence);

    /**
     * @brief Join string views at specified indices into a single string.
     *
     * @param indices Vector of sequence indices.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @return Concatenated string of bases at the specified indices.
     */
    static std::string join_string_views(std::vector<size_t> indices,
                                         const std::vector<std::string_view>& mod_sequence);

    /**
     * @brief Get modified energy or fall back to unmodified energy.
     *
     * Looks up the modified energy parameter for a given key and modified bases.
     * Returns the unmodified energy if no modified parameter is found.
     *
     * @param key The parameter key to look up.
     * @param modified Vector of modified base identifiers.
     * @param mod_params Vector of modified base parameters.
     * @param unmod_energy The unmodified energy to use as fallback.
     * @param lookup_type Type of energy lookup (Stacking, Terminal, etc.).
     * @return Energy value in centicalories.
     */
    static int get_mod_energy(const std::string& key,
                                 const std::vector<std::string_view>& modified,
                                 const std::vector<modified_base_param>& mod_params,
                                 int unmod_energy, ModLookup lookup_type);

   /**
    * @brief Modify a DangleSet's energies based on modified base energy
    * 
    * @param original_set The original DangleSet to modify.
    * @param diffs The ModDiffs containing the energy differences to apply.
    */
    static void modify_dangle_set( DangleSet& original_set, ModDiffs diffs);
   
   /**
    * @brief Get the energy difference between modified and unmodified for a given key and lookup type.
     *
     * @param key The parameter key to look up.
     * @param unique_mod_bases Vector of unique modified base identifiers relevant to this energy calculation.
     * @param mod_params Vector of modified base parameters.
     * @param unmod_energy The unmodified energy to use as fallback.
     * @param lookup_type Type of energy lookup (Stacking, Terminal, etc.).
     * @return Energy difference in centicalories (mod_energy - unmod_energy).
    */
    static int get_mod_energy_difference(const std::string& key,
                                          const std::vector<std::string_view>& modified,
                                          const std::vector<modified_base_param>& mod_params,
                                          int unmod_energy, ModLookup lookup_type);
   
   /**
     * @brief Generates a ModDiffs object which stores the energy differences for modified bases
     *  
     * Calculates the energy differences for terminal AU penalty, mismatch, 5' dangle,
     * and 3' dangle energies due to modified bases in a loop.
     * 
     * @param node The loop node for which to calculate energy differences.
     * @param n5d Encoded 5' dangle nucleotide (or -1 if no dangle).
     * @param n3d Encoded 3' dangle nucleotide (or -1 if no dangle).
     * @param type The pair type of the closing pair of the loop.
     * @param unique_mod_bases Vector of unique modified base identifiers relevant to this loop.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @param vp ViennaRNA model and parameters.
     * @param mod_params Vector of modified base parameters.
     * @param is_external Whether this loop is an external loop (affects which energy parameters to use).
     * @return ModDiffs object containing the energy differences for terminal AU penalty, mismatch, 5' dangle, and 3' dangle energies.
     * 
     */
   static ModDiffs get_mod_dangle_energy_diffs(const LoopNode& node, const int n5d,
                                             const int n3d, const unsigned int type,
                                             const std::vector<std::string_view>& unique_mod_bases,
                                             const std::vector<std::string_view>& mod_sequence,
                                             vrna_md_param& vp,
                                             const std::vector<modified_base_param>& mod_params, 
                                             bool is_external) ;
};

}  // namespace knotergy